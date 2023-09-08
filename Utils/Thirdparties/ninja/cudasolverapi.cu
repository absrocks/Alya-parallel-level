#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusparse_v2.h>
#include <cusolverSp.h>
#include <cusolverSp_LOWLEVEL_PREVIEW.h>
#include "cudakernel.h"

typedef size_t devptr_t;

extern "C" void cusolverspdcsrlsvchol_(devptr_t* han,int *m,int *nnz,devptr_t* desc,devptr_t* val,devptr_t* row,devptr_t* col,devptr_t* rhs,double *tol,int *reod,devptr_t *unk)
{
  cusolverStatus_t stat;
  int sing;
  stat = cusolverSpDcsrlsvchol((cusolverSpHandle_t) (*han), 
			       *m, 
			       *nnz, 
			       (cusparseMatDescr_t) (*desc), 
			       (double*) (*val), 
			       (int*) (*row), 
			       (int*) (*col), 
			       (double*) (*rhs), 
			       *tol, 
			       *reod, 
			       (double*)(*unk), 
			       &sing);
  if(stat != CUSOLVER_STATUS_SUCCESS)
    {
      printf("FAiled CUSOLVER cholesky direct %s\n",cusolverGetError(stat));exit(1);
    }
  if( sing != -1)
    {
      printf("Singularity observerd..value of singularity ---> %d\n",sing);
    }
  return;
}


extern "C" void cusolversphandle_(devptr_t *han)
{
  cusolverStatus_t stat;
  cusolverSpHandle_t hand=0;
  stat = cusolverSpCreate(&hand);
  if(stat != CUSOLVER_STATUS_SUCCESS)
    {
      printf("FAiled handle CUSOLVER %s\n",cusolverGetError(stat));exit(1);
    }
  *han  = (devptr_t)hand;
  return;
}

extern "C" void cusolverspattachstream_(devptr_t *han,devptr_t *str)
{
  cusolverStatus_t stat;
  stat = cusolverSpSetStream((cusolverSpHandle_t)(*han),(cudaStream_t)(*str));
  if(stat != CUSOLVER_STATUS_SUCCESS)
    {
      printf("FAiled handle CUSOLVER attach stream %s\n",cusolverGetError((cusolverStatus_t)stat));exit(1);
    }
  return;
}

extern "C" void cusolverspcreatecsrqrinfo_(devptr_t *info)
{
  csrqrInfo_t infoc=0;
  cusolverStatus_t stat;
  stat = cusolverSpCreateCsrqrInfo(&infoc);
  if(stat != CUSOLVER_STATUS_SUCCESS)
    {
      printf("FAiled cusolver QR csr info object %s\n",cusolverGetError((cusolverStatus_t)stat));exit(1);
    }
  *info = (devptr_t)(infoc);
  return;
}

extern "C" void cusolverspxcsrqranalysis_(devptr_t *hand,int *m,int *n,int *nnza,devptr_t *desc,devptr_t *r,devptr_t *c,devptr_t* info)
{
  cusolverStatus_t stat;
  stat = cusolverSpXcsrqrAnalysis((cusolverSpHandle_t)(*hand),
				  *m,
				  *n,
				  *nnza,
				  (cusparseMatDescr_t)(*desc),
				  (int*)(*r),
				  (int*)(*c),
				  (csrqrInfo_t)(*info));
				  
  if(stat != CUSOLVER_STATUS_SUCCESS)
    {
      printf("FAiled cusolver QR csr analysis %s\n",cusolverGetError((cusolverStatus_t)stat));exit(1);
    }
  return;
}

extern "C" void cusolverspdcsrqrbufferinfo_(devptr_t *hand,int *m,int *n,int *nnza,devptr_t *desc,devptr_t *mat,devptr_t *r,devptr_t *c,devptr_t* info,devptr_t* datain,devptr_t* workin)
{
  cusolverStatus_t stat;
  size_t sz1,sz2;
  stat = cusolverSpDcsrqrBufferInfo((cusolverSpHandle_t)(*hand), 
				    *m,
				    *n, 
				    *nnza, 
				    (cusparseMatDescr_t)(*desc),
				    (double*)(*mat),
				    (int*)(*r), 
				    (int*)(*c), 
				    (csrqrInfo_t)(*info),
				    &sz1,
				    &sz2);
  
  if(stat != CUSOLVER_STATUS_SUCCESS)
    {
      printf("FAiled cusolver QR csr buffer info %s\n",cusolverGetError((cusolverStatus_t)stat));exit(1);
    }
  
  *datain = sz1;
  *workin = sz2;
  return;
}

extern "C" void cusolverspdcsrqrfactor_(devptr_t *hand,int *m,int *n,int *nnza,devptr_t* info,devptr_t* pbuff)
{
  cusolverStatus_t stat;
  double *nnull = NULL;
  stat = cusolverSpDcsrqrFactor((cusolverSpHandle_t)(*hand), 
				*m,
				*n, 
				*nnza, 
				nnull,
				nnull,
				(csrqrInfo_t)(*info),
				(void*)(*pbuff));
  if(stat != CUSOLVER_STATUS_SUCCESS)
    {
      printf("FAiled cusolver QR csr factor info %s\n",cusolverGetError((cusolverStatus_t)stat));exit(1);
    }
  return;
}

extern "C" void cusolverspdcsrqrsetup_(devptr_t *hand,int *m,int *n,int *nnza,devptr_t *desc,devptr_t *mat,devptr_t *r,devptr_t *c,double *mu,devptr_t* info)
{
  cusolverStatus_t stat;
  stat = cusolverSpDcsrqrSetup((cusolverSpHandle_t)(*hand), 
				   *m,
				   *n, 
				   *nnza, 
				   (cusparseMatDescr_t)(*desc),
				   (double*)(*mat),
				   (int*)(*r), 
				   (int*)(*c), 
				   *mu,
				   (csrqrInfo_t)(*info));

  if(stat != CUSOLVER_STATUS_SUCCESS)
    {
      printf("FAiled cusolver QR csr setup %s\n",cusolverGetError((cusolverStatus_t)stat));exit(1);
    }
  return;
}

extern "C" void cusolverspdcsrqrsolve_(devptr_t *hand,int *m,int *n,devptr_t *b,devptr_t *x,devptr_t* info,devptr_t* pbuff)
{
  cusolverStatus_t stat;
  stat = cusolverSpDcsrqrSolve((cusolverSpHandle_t)(*hand), 
				   *m,
				   *n, 
				   (double*)(*b),
				   (double*)(*x), 
				   (csrqrInfo_t)(*info),
				   (void*)(*pbuff));

  if(stat != CUSOLVER_STATUS_SUCCESS)
    {
      printf("FAiled cusolver QR csr solve %s\n",cusolverGetError((cusolverStatus_t)stat));exit(1);
    }
  return;
}

extern "C" void cusolverspdcsrlsvqr_(devptr_t *hand,int *m,int *nnza,devptr_t *desc,devptr_t *mat,devptr_t *row,devptr_t *col,devptr_t *b,double *tol,int *reorder,devptr_t *x)
{
  cusolverStatus_t stat;
  int sing;
  stat = cusolverSpDcsrlsvqr((cusolverSpHandle_t) (*hand), 
			     *m, 
			     *nnza, 
			     (cusparseMatDescr_t) (*desc), 
			     (double*) (*mat), 
			     (int*) (*row), 
			     (int*) (*col), 
			     (double*) (*b), 
			     *tol, 
			     *reorder, 
			     (double*)(*x), 
			     &sing);
  if(stat != CUSOLVER_STATUS_SUCCESS)
    {
      printf("FAiled CUSOLVER QR solver direct %s the number is %d\n",cusolverGetError((cusolverStatus_t)stat),(int)stat);exit(1);
    }
  if( sing != -1)
    {
      printf("Singularity observerd..value of singularity ---> %d\n",sing);
    }
  return;
}

extern "C" void cusolverspxcsrsymrcm_(devptr_t *handle, int *n, int *nnz,devptr_t *desc,int *csrrow,int *csrcol, int *p)
{
  cusolverStatus_t stat;
  stat = cusolverSpXcsrsymrcmHost((cusolverSpHandle_t)(*handle), 
				  *n, 
				  *nnz,
				  (cusparseMatDescr_t)(*desc),
				  csrrow,
				  csrcol,
				  p);
  
  if(stat != CUSOLVER_STATUS_SUCCESS)
    {
      printf("FAiled CUSOLVER computing permutations \n",cusolverGetError((cusolverStatus_t)stat));exit(1);
    }
}

extern "C" void cusolverspxcsrperm_buffersize_(devptr_t *handle, int *m, int *n, int *nnz, devptr_t *desc, int *csrrow, int *csrcol, int *p,int *q, size_t *bsize)
{
  cusolverStatus_t stat;
  size_t sz;
  stat = cusolverSpXcsrperm_bufferSizeHost((cusolverSpHandle_t)(*handle), 
					   *m,
					   *n,
					   *nnz, (cusparseMatDescr_t)(*desc), 
					   csrrow, 
					   csrcol,
					   p,
					   q,
					   &sz);

  *bsize = sz;
  if(stat != CUSOLVER_STATUS_SUCCESS)
    {
      printf("FAiled CUSOLVER buffer size in permutations \n",cusolverGetError((cusolverStatus_t)stat));exit(1);
    }
}

extern "C" void cusolverspxcsrperm_(devptr_t *handle,int *m,int *n,int *nnz,devptr_t *desc,int *csrrow,int *csrcol,int *p,int *q,int *map,void *buff)
{
  cusolverStatus_t stat;
  stat = cusolverSpXcsrpermHost((cusolverSpHandle_t)(*handle),
				*m,
				*n,
				*nnz,
				(cusparseMatDescr_t)(*desc),
				csrrow,
				csrcol, 
				p,
				q, 
				map,
				buff);
  
  if(stat != CUSOLVER_STATUS_SUCCESS)
    {
      printf("FAiled CUSOLVER buffer size in permutations \n",cusolverGetError((cusolverStatus_t)stat));exit(1);
    }
}

extern "C" void cusolverspdestroycsrqrinfo_(devptr_t *info)
{
  cusolverStatus_t stat;
  stat = cusolverSpDestroyCsrqrInfo((csrqrInfo_t)(*info));
  
  if(stat != CUSOLVER_STATUS_SUCCESS)
    {
      printf("FAiled CUSOLVER QR destroy info%s\n",cusolverGetError(stat));exit(1);
    }
  return;
}
