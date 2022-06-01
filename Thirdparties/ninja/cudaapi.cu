#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusparse_v2.h>
#include "cudakernel.h"

typedef size_t devptr_t;

extern "C" void gpumallochost_(devptr_t *ptr,int *sz)
{
  void *x;
  cudaError_t err = cudaMallocHost(&x,*sz);
  if(err != cudaSuccess)
    {
      printf("CUDA malloc host Error %s\n",cudaGetErrorString(err));
      exit(1);
    }
  *ptr = (devptr_t)x;
  return;  
}
extern "C" void gpufreehost_(devptr_t *ptr)
{
  cudaError_t err = cudaFreeHost((void*)(*ptr));
  if(err != cudaSuccess)
    {
      printf("CUDA host free failed Error %s\n",cudaGetErrorString(err));
      exit(1);
    }
  return;  
}
extern "C" void gpuallstreamsync_()
{
  cudaError_t err = cudaThreadSynchronize();
  if(err != cudaSuccess)
    {
      printf("CUDA event ThreadSync failed Error %s\n",cudaGetErrorString(err));
      exit(1);
    }
  return;  
}

extern "C" void gpugetsms_(int *nsms)
{
  int nDevices;  
  cudaDeviceProp prop;
  cudaGetDeviceCount(&nDevices);
  if(nDevices == 0)
    {
      printf("NUMERO DE GPU DEVICE EN SERIO:: %d",nDevices);exit(1);
    }
  cudaGetDeviceProperties(&prop, 0);
  
  *nsms = prop.multiProcessorCount;
  return;
}
extern "C" void gpueventcreate_(devptr_t* x)
{
  cudaEvent_t eve;
  cudaError_t err = cudaEventCreate(&eve);
  if(err != cudaSuccess)
    {
      printf("CUDA event create failed Error %s\n",cudaGetErrorString(err));
      exit(1);
    }
  *x = (devptr_t)(eve);
  return;
}

extern "C" void gpueventdestroy_(devptr_t* x)
{
  cudaError_t err = cudaEventDestroy((cudaEvent_t)(*x));
  if(err != cudaSuccess)
    {
      printf("CUDA event destroy failed Error %s\n",cudaGetErrorString(err));
      exit(1);
    }
  return;
}
extern "C" void gpueventrecord_(devptr_t* x)
{
  cudaError_t err = cudaEventRecord((cudaEvent_t)(*x));
  if(err != cudaSuccess)
    {
      printf("CUDA event Record failed Error %s\n",cudaGetErrorString(err));
      exit(1);
    }
  return;
}
extern "C"  void gpueventelapsedtime_(double *ti,devptr_t *start,devptr_t* stop)
{
  float out;
  cudaError_t err = cudaEventElapsedTime(&out,(cudaEvent_t)(*start),(cudaEvent_t)(*stop));
  if(err != cudaSuccess)
    {
      printf("CUDA event elapsed time failed Error %s\n",cudaGetErrorString(err));
      exit(1);
    }
  *ti = (double)out;
  return;
}

extern "C" void gpureset_()
{
  cudaError_t stat;
  stat =   cudaDeviceReset();
  if(stat != cudaSuccess)
    {
      printf("FAiled to reset GPU  %s\n",cudaGetErrorString(stat));exit(1);
    }
  return;

}
extern "C" void gpuset_(int *id)
{
  int nDevices;  
  int ID = *id;
  if (ID == -1) 
    {
      ID = 0;
    }
  else
    {
      cudaGetDeviceCount(&nDevices);
      if(nDevices == 0)
	{
	  printf("NUMERO DE GPU DEVICE:: %d ,EN SERIO\n",nDevices);exit(1);
	}
      else
	{
	  ID = *id % nDevices;
	}
    }

  if(cudaSetDevice(ID) != cudaSuccess)
    {
      printf("Could not set GPU device id:: %d for MPI rank %d\n",ID,*id );exit(1);
    }

  return;
}

extern "C" void gpusetmany_(int *id, int *rcpugpu)
{
  int nDevices;  
  int ID = *id;
  if (ID == -1) 
    {
      ID = 0;
    }
  else
    {
      cudaGetDeviceCount(&nDevices);
      if(nDevices == 0)
	{
	  printf("NUMERO DE GPU DEVICE:: %d ,EN SERIO\n",nDevices);exit(1);
	}
      else
	{
	    ID = *id / *rcpugpu;
            ID = ID % nDevices;
	}
    }

  if(cudaSetDevice(ID) != cudaSuccess)
    {
      printf("Could not set GPU device id:: %d for MPI rank %d\n",ID,*id );exit(1);
    }

  return;
}



extern "C" void getgpualignment_(int *alignm)
{
  int nDevices;  
  cudaDeviceProp prop;
  cudaGetDeviceCount(&nDevices);
  if(nDevices == 0)
    {
      printf("NUMERO DE GPU DEVICE EN SERIO:: %d",nDevices);exit(1);
    }
  cudaGetDeviceProperties(&prop, 0);
  *alignm = prop.textureAlignment;
}
extern "C" void gpumemgetinfo_(devptr_t *free,devptr_t *total)
{
  size_t fr,tot;
  cudaMemGetInfo(&fr,&tot);
  *free = fr;
  *total = tot;
}

extern "C" void gpumemset_(devptr_t* mem,int *val,size_t* count)
{
  cudaError_t err;
  err = cudaMemset((void*)(*mem),*val,*count);
    if(err != cudaSuccess)
      {
	printf("CUDA memset error %s",cudaGetErrorString(err));exit(1);
      }
    return;
}

extern "C" void gpumalloc_(devptr_t *pt,size_t *N,devptr_t *caller)
{
  void *x;
  cudaError_t err;
  err = cudaMalloc(&x,*N);
  if(err != cudaSuccess)
    {
      printf("CUDA malloc error %s called by %s",cudaGetErrorString(err),(char*)caller);exit(1);
    }
  *pt = (devptr_t)x;
  return;
}

extern "C" void memcpytogpuapi_ (devptr_t *dst,devptr_t *src,int *N,devptr_t *caller)
{
  cudaError_t err;
  err = cudaMemcpy((void*)(*dst),(void*)src,*N,cudaMemcpyHostToDevice);
  if(err != cudaSuccess)
    {
      printf("CUDA memcpy  error to GPU %s called by %s",cudaGetErrorString(err),(char*)caller);exit(1);
    }
  return;
}
extern "C" void memcpyfromgpuapi_ (devptr_t *dst,devptr_t *src,int *N,devptr_t* caller)
{
  cudaError_t err;
  err = cudaMemcpy((void*)dst,(void*)*src,*N,cudaMemcpyDeviceToHost);
    if(err != cudaSuccess)
      {
	printf("CUDA memcpy error from GPU %s called by %s",cudaGetErrorString(err),(char*)caller);exit(1);
      }
    return;
}

extern "C" void memcpytogpuapi2_ (devptr_t *dst,devptr_t *src,size_t *N,devptr_t *caller)
{
  cudaError_t err;
  err = cudaMemcpy((void*)(*dst),(void*)src,*N,cudaMemcpyHostToDevice);
  if(err != cudaSuccess)
    {
      printf("CUDA memcpy  error to GPU %s called by %s",cudaGetErrorString(err),(char*)caller);exit(1);
    }
  return;
}
extern "C" void memcpyfromgpuapi2_ (devptr_t *dst,devptr_t *src,size_t *N,devptr_t* caller)
{
  cudaError_t err;
  err = cudaMemcpy((void*)dst,(void*)*src,*N,cudaMemcpyDeviceToHost);
    if(err != cudaSuccess)
      {
	printf("CUDA memcpy error from GPU %s called by %s",cudaGetErrorString(err),(char*)caller);exit(1);
      }
    return;
}


extern "C" void  gpustreamcreate_(devptr_t *str)
{
  cudaError_t stat;
  cudaStream_t s;
  stat = cudaStreamCreate(&s);
  if(stat != cudaSuccess)
    {
      printf("FAiled create cuda stream\n");exit(1);
    }
  *str  = (devptr_t)s;
  return;
}

extern "C" void memcpyingpuasyncapi_ (devptr_t *dst,devptr_t *src,int *N,devptr_t *str,devptr_t *caller)
{
  cudaError_t err;
  err = cudaMemcpyAsync((void*)(*dst),(void*)(*src),*N,cudaMemcpyDeviceToDevice,(cudaStream_t) (*str));
    if(err != cudaSuccess)
      {
	printf("CUDA memcpy  device to device error %s called by %s",cudaGetErrorString(err),(char*)(caller));exit(1);
      }
    return;
}

extern "C" void memcpytogpuasyncapi_ (devptr_t *dst,devptr_t *src,int *N,devptr_t *str,devptr_t *caller)
{
  cudaError_t err;
  err = cudaMemcpyAsync((void*)(*dst),(void*)src,*N,cudaMemcpyHostToDevice,(cudaStream_t) (*str));
    if(err != cudaSuccess)
      {
	printf("CUDA memcpy  error %s called by %s",cudaGetErrorString(err),(char*)(caller));exit(1);
      }
    return;
}
extern "C" void memcpyfromgpuasyncapi_ (devptr_t *dst,devptr_t *src,int *N,devptr_t *str,devptr_t *caller)
{
  cudaError_t err;
  err = cudaMemcpyAsync((void*)dst,(void*)*src,*N,cudaMemcpyDeviceToHost,(cudaStream_t) (*str));
    if(err != cudaSuccess)
      {
	printf("CUDA memcpy error %s called by %s",cudaGetErrorString(err),(char*)(caller));exit(1);
      }
    return;
}

extern "C" void memcpytogpuasyncapi2_ (devptr_t *dst,devptr_t *src,size_t *N,devptr_t *str,devptr_t *caller)
{
  cudaError_t err;
  err = cudaMemcpyAsync((void*)(*dst),(void*)src,*N,cudaMemcpyHostToDevice,(cudaStream_t) (*str));
    if(err != cudaSuccess)
      {
	printf("CUDA memcpy  error %s called by %s",cudaGetErrorString(err),(char*)(caller));exit(1);
      }
    return;
}
extern "C" void memcpyfromgpuasyncapi2_ (devptr_t *dst,devptr_t *src,size_t *N,devptr_t *str,devptr_t *caller)
{
  cudaError_t err;
  err = cudaMemcpyAsync((void*)dst,(void*)*src,*N,cudaMemcpyDeviceToHost,(cudaStream_t) (*str));
    if(err != cudaSuccess)
      {
	printf("CUDA memcpy error %s called by %s",cudaGetErrorString(err),(char*)(caller));exit(1);
      }
    return;
}

extern "C" void memcpytogpuasyncapi3_ (devptr_t *dst,devptr_t *src,size_t *N,devptr_t *str,devptr_t *caller)
{
  cudaError_t err;
  err = cudaMemcpyAsync((void*)(*dst),(void*)(*src),*N,cudaMemcpyHostToDevice,(cudaStream_t) (*str));
    if(err != cudaSuccess)
      {
	printf("CUDA memcpy  error %s called by %s",cudaGetErrorString(err),(char*)(caller));exit(1);
      }
    return;
}
extern "C" void memcpyfromgpuasyncapi3_ (devptr_t *dst,devptr_t *src,size_t *N,devptr_t *str,devptr_t *caller)
{
  cudaError_t err;
  err = cudaMemcpyAsync((void*)(*dst),(void*)(*src),*N,cudaMemcpyDeviceToHost,(cudaStream_t) (*str));
    if(err != cudaSuccess)
      {
	printf("CUDA memcpy error %s called by %s",cudaGetErrorString(err),(char*)(caller));exit(1);
      }
    return;
}


extern "C" void  gpustreamsync_(devptr_t *str)
{
  cudaError_t stat;
  stat = cudaStreamSynchronize((cudaStream_t)(*str));
  if(stat != cudaSuccess)
    {
      printf("FAiled to sync cuda stream error %s\n",cudaGetErrorString(stat));exit(1);
    }
  return;
}

extern "C" void  gpusync_()
{
  cudaError_t stat;
  stat = cudaDeviceSynchronize();
  if(stat != cudaSuccess)
    {
      printf("FAiled to sync GPU  %s\n",cudaGetErrorString(stat));exit(1);
    }
  return;
}

extern "C" void int2cptr_(devptr_t *pt,void **out)
{
  *out = (void*)(*pt);
  return;
}


extern "C" void gpuhostregister_(void *ptr,size_t *sz)
{

  cudaError_t stat;
  stat = cudaHostRegister(ptr,*sz,cudaHostRegisterDefault);
  if(stat != cudaSuccess)
    {
      printf("FAiled to pin host memory %s\n",cudaGetErrorString(stat));exit(1);
    }
  return;
}

extern "C" void gpuhostunregister_(void *ptr)
{
  cudaError_t stat;
  stat = cudaHostUnregister(ptr);
  if(stat != cudaSuccess)
    {
      printf("FAiled to unpin host memory %s\n",cudaGetErrorString(stat));exit(1);
    }
  return;
}

extern "C" void  gpulasterror_()
{
  cudaError_t stat;
  stat = cudaGetLastError();
  if(stat != cudaSuccess)
    {
      printf("Last Error GPU  %s\n",cudaGetErrorString(stat));exit(1);
    }
  return;
}


