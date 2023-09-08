/******************************************************************************
 * Copyright 1986, 2016 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// CUDA kernel for CUDA IPC volume editing

#include <cstdio>
#include <string.h>
#include <cuda_runtime.h>


#include "cuda_ipc_volume_editing_kernel.cuh"

namespace cuda {

//
// Helper functions
//

//----------------------------------------------------------------------
// Terminates with an error message if a CUDA error is detected
__host__ void check(cudaError_t error, const char* msg = 0)
{
    if (error != cudaSuccess)
    {
        int current_gpu = -1;
        cudaGetDevice(&current_gpu);

        printf("CUDA error on GPU %d: '%s'", current_gpu, cudaGetErrorString(error));
        if (msg)
        {
            printf(", debug message: '%s'", msg);
        }
        printf("\n");

        exit(1);
    }
}

//----------------------------------------------------------------------
// 3D surface to access the volume brick (stored as a CUDA array on the device), may include
// boundary voxels
surface<void, cudaSurfaceType3D> volume_brick_surf;

//----------------------------------------------------------------------
// Writes a voxel to the 3D surface
__device__ __forceinline__ void write_voxel(

    int x, int y, int z,
    unsigned char value)
{
    surf3Dwrite(
        value, volume_brick_surf,
        x * sizeof(unsigned char), y, z);
}

//----------------------------------------------------------------------
// Reads a voxel from the 3D surface
__device__ __forceinline__ unsigned char read_voxel(
    int x, int y, int z)
{
    unsigned char value;
    surf3Dread(
        &value, volume_brick_surf,
        x * sizeof(unsigned char), y, z);
    return value;
}

//----------------------------------------------------------------------
//
// Compute kernel
//
__global__ void set_voxel_values_kernel(
    int3          brick_size,
    unsigned char *values)
{
   // Threads iterate over z
   const int z = blockIdx.x * blockDim.x + threadIdx.x;
   if (z >= brick_size.z)
       return;

   // For each thread iterate over y and x
   values = values + z*brick_size.x*brick_size.y;
   for (int y = 0; y < brick_size.y; ++y)
   {
       for (int x = 0; x < brick_size.x; ++x, values++)
       {
           write_voxel(x, y, z, *values);
       }
   }
}

//----------------------------------------------------------------------
#define USE_COPY_KERNEL     0   
// USE_COPY_KERNEL == 0,  Copy brick directly by CUDA cudaMemcpy3D(). Faster (default)
// USE_COPY_KERNEL == 1,  Copy brick through a kernel. Slower but more flexible in case of conversions.

void set_voxel_values(
    cudaArray*                  dst_brick,
    const int                   dst_brick_size[3],
    const cudaIpcMemHandle_t&   src_brick_ipc_handle,
    unsigned char*              d_src_brick,
    const int                   src_brick_size[3],
    size_t                      src_mem_offset,
    int                         device_id)
{
    cudaError_t error;
    
    error = cudaSetDevice(device_id);
    check(error, "setting device");
        
    // device_id == 0 is running on the same thread so no need to get the brick address by its CUDA IPC handle.
    unsigned char *d_brick = NULL;
    if(device_id != 0) 
    {
        error = cudaIpcOpenMemHandle((void **)&d_brick, src_brick_ipc_handle, cudaIpcMemLazyEnablePeerAccess);
        check(error, "opening memory handle");
    }
    else
    {
        d_brick = d_src_brick;
    }
    
    // d_brick += src_mem_offset;
    
#if USE_COPY_KERNEL == 1

    // Brick size components are flipped in device memory
    const int3 brick_size = make_int3(dst_brick_size[2], dst_brick_size[1], dst_brick_size[0]);
    
    // Bind the surface to the CUDA array
    error = cudaBindSurfaceToArray(volume_brick_surf, dst_brick);
    check(error, "bind surface to array");
    
    // Set up kernel launch configuration: Let the threads iterate over z
    const unsigned int block_dim = 32;

    unsigned int block_count = brick_size.z / block_dim;
    if (brick_size.z % block_dim > 0)
    {
        block_count++;
    }

    dim3 grid(block_count, 1, 1);
    dim3 block(block_dim, 1, 1);

    // Launch the compute kernel
    set_voxel_values_kernel<<<grid, block, 0>>>(brick_size, d_brick + src_mem_offset);

    // Check for CUDA errors
    check(cudaPeekAtLastError(), "kernel launch arguments");
    // Syncronizes device to get errors immediately (prevents asynchronous execution, should not be
    // used in production code)
    check(cudaDeviceSynchronize(), "kernel execution");

#else
    
    uint3 src_extent = make_uint3(src_brick_size[2], src_brick_size[1], src_brick_size[0]);
    uint3 src_offset;
    src_offset.x = (src_brick_size[2] - dst_brick_size[2])/2;
    src_offset.y = (src_brick_size[1] - dst_brick_size[1])/2;
    src_offset.z = (src_brick_size[0] - dst_brick_size[0])/2;
    uint3 dst_offset = make_uint3(0, 0, 0);
    uint3 copy_extend = make_uint3(dst_brick_size[2], dst_brick_size[1], dst_brick_size[0]);
    
    error = copy_device_to_device(dst_brick, d_brick + src_mem_offset, sizeof(unsigned char),
        src_extent, src_offset, dst_offset, copy_extend, 0);
    
    check(error, "copying brick to array");
    
#endif // USE_COPY_KERNEL

    if(device_id != 0) 
    {
        // error = cudaSetDevice(device_id);
        // check(error, "setting device before closing memory handle+");

        error = cudaIpcCloseMemHandle(d_brick);
        check(error, "closing memory handle");
    }
    
}

//----------------------------------------------------------------------
__host__ void alloc_IPC_brick(
    cudaIpcMemHandle_t&     ipc_handle,
    unsigned char**         d_voxels,
    const unsigned char*    h_voxels, 
    int                     nb_voxels,
    int                     device_id)
{
    cudaError_t error;
    
    error = cudaSetDevice(device_id);
    check(error, "setting device");
    
    error = cudaMalloc(&(*d_voxels), sizeof(unsigned char)*nb_voxels);    
    check(error, "vector allocation");
    
    error = cudaMemcpy(*d_voxels, h_voxels, sizeof(unsigned char)*nb_voxels, cudaMemcpyHostToDevice);  
    check(error, "copying vector from host to device");

    error = cudaIpcGetMemHandle(&ipc_handle, *d_voxels);
    check(error, "getting IPC handle");
}

//----------------------------------------------------------------------
__host__ void get_IPC_brick(
    unsigned char*      h_voxels, 
    int                 nb_voxels, 
    cudaIpcMemHandle_t& ipc_handle, 
    int                 device_id)
{
    unsigned char *d_vector;
    cudaError_t error;
    
    error = cudaSetDevice(device_id);
    check(error, "setting device");
    
    error = cudaIpcOpenMemHandle((void **)&d_vector, ipc_handle, cudaIpcMemLazyEnablePeerAccess);
    check(error, "opening memory handle");
        
    error = cudaMemcpy(h_voxels, d_vector, sizeof(unsigned char)*nb_voxels, cudaMemcpyDeviceToHost);  
    check(error, "copying vector from device to host");
    
}

//----------------------------------------------------------------------
__host__ void free_IPC_brick(
    cudaIpcMemHandle_t& ipc_handle, 
    int                 device_id)
{
    unsigned char *d_vector;
    cudaError_t error;
    
    error = cudaSetDevice(device_id);
    check(error, "setting device");
    
    error = cudaIpcOpenMemHandle((void**)&d_vector, ipc_handle, cudaIpcMemLazyEnablePeerAccess);
    check(error, "opening memory handle");
    
    error = cudaIpcCloseMemHandle((void *) d_vector);
    check(error, "closing memory handle");
    
    error = cudaFree(d_vector);
    check(error, "vector deallocation");
}

//----------------------------------------------------------------------
//
// Update brick kernel
//
__global__ void update_brick_kernel(unsigned char* d_voxels, int num_lines, int line_size)
{
    int idx = (blockIdx.x*blockDim.x + threadIdx.x);
    if(idx < num_lines)
    {
        idx *= line_size;
        for(int i=idx; i<idx+line_size; i++)
        {
            d_voxels[i] = (d_voxels[i] + 1)%255;
        }
    }
}

//----------------------------------------------------------------------
__host__ void update_brick(
    unsigned char*  d_voxels,
    int             brick_size[3],
    int             device_id)
{
    cudaError_t error;
    
    error = cudaSetDevice(device_id);
    check(error, "setting device");
    
    // Set up kernel launch configuration: Let the threads iterate over a z-line
    const int BLOCK_SIZE = 32;
    int num_lines = brick_size[0]*brick_size[1];
    int num_blocks = num_lines/BLOCK_SIZE + ((num_lines%BLOCK_SIZE > 0) ? 1 : 0);

    update_brick_kernel<<<dim3(num_blocks), dim3(BLOCK_SIZE), 0>>>(d_voxels, num_lines, brick_size[2]);
    
    // Check for CUDA errors
    check(cudaPeekAtLastError(), "update brick kernel launch arguments");
    // Syncronizes device to get errors immediately (prevents asynchronous execution, should not be
    // used in production code)
    check(cudaDeviceSynchronize(), "update brick kernel execution");
    
}

//----------------------------------------------------------------------
__host__ cudaError_t copy_device_to_device(
    cudaArray*                              pdst,
    const void*                             psrc,
    size_t                                fmt_size,
    uint3&  src_extent,
    uint3&  src_offset,
    uint3&  dst_offset,
    uint3&  copy_extend,
    const cudaStream_t cuda_stream)
{
    cudaMemcpy3DParms cp_param = {0};

    const size_t            data_start_off =   static_cast<size_t>(src_offset.x)
                                               + static_cast<size_t>(src_offset.y) * src_extent.x
                                               + static_cast<size_t>(src_offset.z) * src_extent.x * src_extent.y;
    const unsigned char* data_ptr       = (const unsigned char*)(psrc);
    const unsigned char* data_begin_ptr = data_ptr + data_start_off;
    cp_param.srcPtr =
        make_cudaPitchedPtr(
            const_cast<unsigned char*>(data_begin_ptr),
            src_extent.x * fmt_size,
            src_extent.x, src_extent.y);

    cp_param.dstArray = pdst;
    cp_param.dstPos   = make_cudaPos(dst_offset.x, dst_offset.y, dst_offset.z);
    cp_param.extent   = make_cudaExtent(copy_extend.x, copy_extend.y, copy_extend.z);

    cp_param.kind     = cudaMemcpyDeviceToDevice;

    cudaError cu_err  = cudaSuccess;
    if (cuda_stream != 0) 
        cu_err = cudaMemcpy3DAsync(&cp_param, cuda_stream);
    else 
        cu_err = cudaMemcpy3D(&cp_param);
    
    return cu_err;
}

//----------------------------------------------------------------------
__host__ void enum_cuda_devices(
    char    *buff, 
    int     buff_len)
{
    char tmp_buff[buff_len];
    int i, dev_count;
    
    cudaGetDeviceCount(&dev_count);
    sprintf(tmp_buff, "Number of devices: %d", dev_count);
    strncat(buff, tmp_buff, buff_len);
    
    for(i=0; i<dev_count; i++)
    {
        cudaDeviceProp dev_prop;
        
        cudaGetDeviceProperties(&dev_prop, i);
        sprintf(tmp_buff, " %d:%s", i, dev_prop.name);
        strncat(buff, tmp_buff, buff_len);
        
    }
}

//----------------------------------------------------------------------
__host__ int get_nb_gpus()
{
    int dev_count;

    cudaError_t error;
    
    error = cudaGetDeviceCount(&dev_count);
    check(error, "getting number of GPUs");
    
    return dev_count;
}

__host__ void set_flags(int device_id)
{
    cudaError_t error;
    error = cudaSetDeviceFlags(cudaDeviceLmemResizeToMax);
    check(error, "setting device flags");
}

} // cuda
