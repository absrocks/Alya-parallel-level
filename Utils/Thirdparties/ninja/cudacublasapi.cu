#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include "cudakernel.h"

typedef size_t devptr_t;

extern "C" void cudaddot_(devptr_t *han,const int *n, const devptr_t *devPtrx, const int *incx,const devptr_t *devPtry, const int *incy,double *out)
{
  cublasStatus_t stat;
  double result;
  double *x = (double *)(*devPtrx);
  double *y = (double *)(*devPtry);
  stat = cublasDdot ((cublasHandle_t)*han,*n, x, *incx, y, *incy,&result);
  if(stat != CUBLAS_STATUS_SUCCESS)
    {
      printf("Failed DDOT %s \n",cublasGetError(stat));exit(1);
    }
  *out = result;
  return;
}
extern "C" void cudadaxpy_(devptr_t *han,const int *n,const double *alpha,const devptr_t *devPtrx, const int *incx,const devptr_t *devPtry, const int *incy)
{
  double co = *alpha;
  cublasStatus_t stat;
  double *x = (double *)(*devPtrx);
  double *y = (double *)(*devPtry);
  stat = cublasDaxpy ((cublasHandle_t)(*han),*n,&co, x, *incx, y, *incy);
  if(stat != CUBLAS_STATUS_SUCCESS)
    {
      printf("Failed AXPY\n");exit(1);
    }
  return;
}
extern "C" void cudadcopy_(devptr_t *han,const int *n, const devptr_t *devPtrx, const int *incx,const devptr_t *devPtry, const int *incy)
{
  cublasStatus_t stat;
  double *x = (double *)(*devPtrx);
  double *y = (double *)(*devPtry);
  stat = cublasDcopy((cublasHandle_t)*han,*n, x, *incx, y, *incy);
  if(stat != CUBLAS_STATUS_SUCCESS)
    {
      printf("Failed DCOPY\n");exit(1);
    }
  return;
}
extern "C" void cudadscal_(devptr_t *han,const int *n,double *alpha,devptr_t *devPtry, const int *incy)
{
  double al = *alpha;
  cublasStatus_t stat;
  double *y = (double *)(*devPtry);
  stat = cublasDscal ((cublasHandle_t)*han,*n,&al, y, *incy);
  if(stat != CUBLAS_STATUS_SUCCESS)
    {
      printf("Failed DSCAl\n");exit(1);
    }
  return;
}
extern "C" void cudadnrm2_(devptr_t *han,const int *n,devptr_t *devPtry, const int *incy,double *out)
{
  cublasStatus_t stat;
  double *y = (double *)(*devPtry);
  double result;
  
  stat = cublasDnrm2 ((cublasHandle_t)*han,*n,y, *incy,&result);
  if(stat != CUBLAS_STATUS_SUCCESS)
    {
      printf("Failed DNRM2\n");exit(1);
    }
  *out = result;
  return;
}
extern "C" void  cublashandle_(devptr_t *han)
{
  cublasStatus_t stat;
  cublasHandle_t hand=0;
  stat = cublasCreate(&hand);
  if(stat != CUBLAS_STATUS_SUCCESS)
    {
      printf("FAiled handle CUBLAS\n");exit(1);
    }
  *han  = (devptr_t)hand;
  return;
}

extern "C" void  cublashandledestroy_(devptr_t *han)
{
  cublasStatus_t stat;
  stat = cublasDestroy((cublasHandle_t)(*han));
  if(stat != CUBLAS_STATUS_SUCCESS)
    {
      printf("FAiled handle destroy CUBLAS\n");exit(1);
    }
  return;
}

extern "C" void cublasattachstream_(devptr_t *han,devptr_t *str)
{
  cublasStatus_t stat;
  stat = cublasSetStream((cublasHandle_t)(*han),(cudaStream_t)(*str));
  if(stat != CUBLAS_STATUS_SUCCESS)
    {
      printf("FAiled handle CUBLAS attach stream\n");exit(1);
    }
}

extern "C" void cudasdot_(devptr_t *han,const int *n, const devptr_t *devPtrx, const int *incx,const devptr_t *devPtry, const int *incy,float *out,devptr_t *caller)
{
  cublasStatus_t stat;
  float result;
  float *x = (float *)(*devPtrx);
  float *y = (float *)(*devPtry);
  stat = cublasSdot ((cublasHandle_t)*han,*n, x, *incx, y, *incy,&result);
  if(stat != CUBLAS_STATUS_SUCCESS)
    {
      printf("Failed DDOT %s called by %s\n",cublasGetError(stat),(char*)caller);exit(1);
    }
  *out = result;
  return;
}
extern "C" void cudasaxpy_(devptr_t *han,const int *n,const float *alpha,const devptr_t *devPtrx, const int *incx,const devptr_t *devPtry, const int *incy)
{
  float co = *alpha;
  cublasStatus_t stat;
  float *x = (float *)(*devPtrx);
  float *y = (float *)(*devPtry);
  stat = cublasSaxpy ((cublasHandle_t)(*han),*n,&co, x, *incx, y, *incy);
  if(stat != CUBLAS_STATUS_SUCCESS)
    {
      printf("Failed AXPY\n");exit(1);
    }
  return;
}
extern "C" void cudascopy_(devptr_t *han,const int *n, const devptr_t *devPtrx, const int *incx,const devptr_t *devPtry, const int *incy)
{
  cublasStatus_t stat;
  float *x = (float *)(*devPtrx);
  float *y = (float *)(*devPtry);
  stat = cublasScopy((cublasHandle_t)*han,*n, x, *incx, y, *incy);
  if(stat != CUBLAS_STATUS_SUCCESS)
    {
      printf("Failed DCOPY\n");exit(1);
    }
  return;
}
extern "C" void cudasscal_(devptr_t *han,const int *n,float *alpha,devptr_t *devPtry, const int *incy)
{
  float al = *alpha;
  cublasStatus_t stat;
  float *y = (float *)(*devPtry);
  stat = cublasSscal ((cublasHandle_t)*han,*n,&al, y, *incy);
  if(stat != CUBLAS_STATUS_SUCCESS)
    {
      printf("Failed DSCAl\n");exit(1);
    }
  return;
}
extern "C" void cudasnrm2_(devptr_t *han,const int *n,devptr_t *devPtry, const int *incy,float *out)
{
  cublasStatus_t stat;
  float *y = (float *)(*devPtry);
  float result;
  
  stat = cublasSnrm2 ((cublasHandle_t)*han,*n,y, *incy,&result);
  if(stat != CUBLAS_STATUS_SUCCESS)
    {
      printf("Failed DNRM2\n");exit(1);
    }
  *out = result;
  return;
}
