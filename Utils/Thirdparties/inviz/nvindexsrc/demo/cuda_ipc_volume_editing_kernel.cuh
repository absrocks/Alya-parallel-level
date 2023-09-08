/******************************************************************************
 * Copyright 1986, 2016 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// CUDA kernel for volume editing

#ifndef NVIDIA_INDEX_CUDA_IPC_VOLUME_EDITING_KERNEL_CUH
#define NVIDIA_INDEX_CUDA_IPC_VOLUME_EDITING_KERNEL_CUH

namespace cuda
{
//----------------------------------------------------------------------
/// Copy a brick already allocated in device memory (by the application)
/// to the internal IndeX cudaArray representation by using CUDA IPC
///
/// \param[out] dst_brick       The destination IndeX's volume brick pointer
/// \param[in]  dst_brick_size  The size (sx,sy,sz) of the destination brick
/// \param[in]  src_ipc_handle  The CUDA IPC handle of the source brick
/// \param[in]  d_src_brick     The device memory pointer of the source brick
///                             Used when CUDA IPC handle is not available
/// \param[in]  src_brick_size  The size (sx,sy,sz) of the source brick
/// \param[in]  src_mem_offset  Memory offset for the source brick pointer
/// \param[in]  device_id       The device id of the GPU that owns the brick
///
void set_voxel_values(
    cudaArray*                  dst_brick,
    const int                   dst_brick_size[3],
    const cudaIpcMemHandle_t&   src_brick_ipc_handle,
    unsigned char*              d_src_brick,
    const int                   src_brick_size[3],
    size_t                      src_mem_offset,
    int                         device_id);

//----------------------------------------------------------------------
/// Upload(allocate) a host brick on device memory and return its CUDA IPC handle
///
/// \param[out] ipc_handle  The CUDA IPC handle of the destination brick
/// \param[out] d_voxels    Device pointer of the destination brick
/// \param[in] h_voxels     Host pointer of the source brick.
/// \param[in] nb_voxels    The total number of voxels for the brick.
/// \param[in] device_id    The device id of the GPU that owns the brick
///
void alloc_IPC_brick(
    cudaIpcMemHandle_t&     ipc_handle, 
    unsigned char**         d_voxels, 
    const unsigned char*    h_voxels, 
    int                     nb_voxels,
    int                     device_id);

//----------------------------------------------------------------------
/// Download a pre allocated brick to host array by using its CUDA IPC handle
///
/// \param[out] h_voxels    Host pointer of the destination brick.
/// \param[in] nb_voxels    The total number of voxels for the brick.
/// \param[in] ipc_handle   The CUDA IPC handle of the source brick
/// \param[in] device_id    The device id of the GPU that owns the brick
///
void get_IPC_brick(
    unsigned char*      h_voxels, 
    int                 nb_voxels, 
    cudaIpcMemHandle_t& ipc_handle, 
    int                 device_id);

//----------------------------------------------------------------------
/// Free device memory associated to the brick
///
/// \param[in] ipc_handle   The CUDA IPC handle of the source brick
/// \param[in] device_id    The device id of the GPU that owns the brick
///
void free_IPC_brick(
    cudaIpcMemHandle_t& ipc_handle, 
    int                 device_id);
    
//----------------------------------------------------------------------
/// Sample brick update function. It changes the content of a given brick
/// For this example the functions increases all voxels values in 1.
///
/// \param[in] d_voxels     Device pointer of the source brick
/// \param[in] brick_size   The brick size (sx, sy, sz)
/// \param[in] device_id    The device id of the GPU that owns the brick
///
void update_brick(
    unsigned char*  d_voxels,
    int             brick_size[3],
    int             device_id);

//----------------------------------------------------------------------
/// Utility functions
    
//----------------------------------------------------------------------
/// Device to device memory copy from some source to a CudaArray
cudaError_t copy_device_to_device(
    cudaArray*                              pdst,
    const void*                             psrc,
    size_t                                  fmt_size,
    uint3&  src_extent,
    uint3&  src_offset,
    uint3&  dst_offset,
    uint3&  copy_extend,
    const cudaStream_t cuda_stream);
    
//----------------------------------------------------------------------
/// Build a string with all device information for the current host
void enum_cuda_devices(
    char    *buff, 
    int     buff_len);

//----------------------------------------------------------------------
/// Get the number of GPU available for the current host
int get_nb_gpus();

void set_flags(int device_id);
    
} // cuda


#endif // NVIDIA_INDEX_CUDA_IPC_VOLUME_EDITING_KERNEL_CUH
