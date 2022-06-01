#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusparse_v2.h>
#include "cudakernel.h"

typedef size_t devptr_t;

extern "C" void cusparsehandle_(devptr_t *handle)
{
  cusparseHandle_t thandle = 0; 
  cusparseStatus_t stat = cusparseCreate(&thandle);
  if(stat != CUSPARSE_STATUS_SUCCESS)
    {
      printf("CUDA CUSPARSE create handle error ");exit(1);
    }
  *handle  = (devptr_t)thandle;
  return;
}

extern "C" void cusparse_mat_descr_create_(devptr_t *descrA)
{
  cusparseMatDescr_t tdescrA = 0;
  cusparseStatus_t stat = cusparseCreateMatDescr(&tdescrA);
  if(stat != CUSPARSE_STATUS_SUCCESS)
    {
      printf("CUDA matrix desc create error ");exit(1);
    }
  *descrA  = (devptr_t)tdescrA;
  return;
}

extern "C" void  cusparse_set_mat_base_(devptr_t *descrA, int *base)
{
  cusparseStatus_t stat = cusparseSetMatIndexBase((cusparseMatDescr_t)(*descrA),(cusparseIndexBase_t)(*base));
  if(stat != CUSPARSE_STATUS_SUCCESS)
    {
      printf("CUDA mat base set  error ");exit(1);
    }
  return;
}
extern "C" void  cusparse_set_mat_type_(devptr_t *descrA, int *base)
{
  cusparseStatus_t stat = cusparseSetMatType((cusparseMatDescr_t)(*descrA), (cusparseMatrixType_t)(*base));
  if(stat != CUSPARSE_STATUS_SUCCESS)
    {
      printf("CUDA mat base type  error ");exit(1);
    }
  return;
}
extern "C" void cudadbsrmv_(devptr_t *handle,int *dir,int *transA,int *mb,int *nb,int *nnzb,double *alpha,devptr_t *descrA,devptr_t *bsrValA,devptr_t *bsrRowPtrA,devptr_t *bsrColIndA,int *blockdim,devptr_t *x,double *beta, devptr_t *y)
{
  const double al = *alpha;
  const double bl = *beta;
  
  double *rhs = (double *)(*x);
  double *unk = (double *)(*y);
  double *A = (double *)(*bsrValA);
  int *r = (int *)(*bsrRowPtrA);
  int *c = (int *)(*bsrColIndA);
  cusparseStatus_t stat;
  stat = cusparseDbsrmv((cusparseHandle_t)(*handle),
			(cusparseDirection_t)(*dir),
			(cusparseOperation_t)(*transA),
			*mb, 
			*nb, 
			*nnzb,
			&al,
			(cusparseMatDescr_t)(*descrA), 
			A, 
			r, 
			c, 
			*blockdim,
			rhs,
			&bl, 
			unk);
  if(stat != CUSPARSE_STATUS_SUCCESS)
    {
      printf("CUDA mat DBSRMV  error ");exit(1);
    }
  return;
}

extern "C" void cusparseattachstream_(devptr_t *han,devptr_t *str)
{
  cusparseStatus_t stat;
  stat = cusparseSetStream((cusparseHandle_t)(*han),(cudaStream_t)(*str));
  if(stat != CUSPARSE_STATUS_SUCCESS)
    {
      printf("FAiled handle sparse attach stream\n");exit(1);
    }
}


extern "C" void cudadbsr2csr_(devptr_t* han,int *dir,int *mb,int *nb,devptr_t* descbsr,devptr_t* valbsr,devptr_t* rowbsr,devptr_t* colbsr,int *blockdim,devptr_t* desccsr,devptr_t* valcsr,devptr_t* rowcsr,devptr_t* colcsr)
{
  double *A = (double *)(*valbsr);
  int *ra = (int *)(*rowbsr);
  int *ca = (int *)(*colbsr);
  double *C = (double *)(*valcsr);
  int *rc = (int *)(*rowcsr);
  int *cc = (int *)(*colcsr);
  
  cusparseStatus_t stat;
  stat = cusparseDbsr2csr((cusparseHandle_t)(*han),
			(cusparseDirection_t)(*dir),
			*mb, 
			*nb, 
			(cusparseMatDescr_t)(*descbsr), 
			A, 
			ra, 
			ca, 
			*blockdim,
			(cusparseMatDescr_t)(*desccsr), 
			C, 
			rc, 
			cc);
			
  if(stat != CUSPARSE_STATUS_SUCCESS)
    {
      printf("CUDA mat format convert error %s\n",cublasGetError((cublasStatus_t)stat));exit(1);
    }
  return;  

}


extern "C" void cudadcsrmv_(devptr_t *handle,int *transA,int *mb,int *nb,int *nnzb,double *alpha,devptr_t *descrA,devptr_t *bsrValA,devptr_t *bsrRowPtrA,devptr_t *bsrColIndA,devptr_t *x,double *beta, devptr_t *y)
{
  const double al = *alpha;
  const double bl = *beta;
  
  double *rhs = (double *)(*x);
  double *unk = (double *)(*y);
  double *A = (double *)(*bsrValA);
  int *r = (int *)(*bsrRowPtrA);
  int *c = (int *)(*bsrColIndA);
  cusparseStatus_t stat;
  stat = cusparseDcsrmv((cusparseHandle_t)(*handle),
			(cusparseOperation_t)(*transA),
			*mb, 
			*nb, 
			*nnzb,
			&al,
			(cusparseMatDescr_t)(*descrA), 
			A, 
			r, 
			c, 
			rhs,
			&bl, 
			unk);
  if(stat != CUSPARSE_STATUS_SUCCESS)
    {
      printf("CUDA mat CSRMV  error ");exit(1);
    }
  return;
}

extern "C" void cusparse_mat_hyb_create_(devptr_t *mat)
{
 cusparseStatus_t stat;
 cusparseHybMat_t hybmat=0;
 stat = cusparseCreateHybMat(&hybmat);
 if(stat != CUSPARSE_STATUS_SUCCESS)
   {
     printf("CUDA mat hyb mat create  error ");exit(1);
   }
 *mat  = (devptr_t)hybmat;
 return;
}

extern "C" void cudadcsr2hyb_(devptr_t* han,int *mb,int *nb,devptr_t* desccsr,devptr_t* valcsr,devptr_t* rowcsr,devptr_t* colcsr,devptr_t* hybmat,int *width,int *part)
{
  double *C = (double *)(*valcsr);
  int *rc = (int *)(*rowcsr);
  int *cc = (int *)(*colcsr);
  
  cusparseStatus_t stat;
  stat = cusparseDcsr2hyb((cusparseHandle_t)(*han),
			  *mb, 
			  *nb, 
			  (cusparseMatDescr_t)(*desccsr), 
			  C, 
			  rc, 
			  cc,
			  (cusparseHybMat_t)(*hybmat),
			  *width,
			  (cusparseHybPartition_t)(*part));
  
  if(stat != CUSPARSE_STATUS_SUCCESS)
    {
      printf("CUDA mat format convert CSR2HYB error %s\n",cublasGetError((cublasStatus_t)stat));exit(1);
    }
  return;  

}

extern "C" void cudadhybmv_(devptr_t *handle,int *transA,double *alpha,devptr_t *descrA,devptr_t *hybmat,devptr_t *x,double *beta, devptr_t *y)
{
  const double al = *alpha;
  const double bl = *beta;
  
  double *rhs = (double *)(*x);
  double *unk = (double *)(*y);
  cusparseStatus_t stat;
  stat = cusparseDhybmv((cusparseHandle_t)(*handle),
			(cusparseOperation_t)(*transA),
			&al,
			(cusparseMatDescr_t)(*descrA), 
			(cusparseHybMat_t)(*hybmat),
			rhs,
			&bl, 
			unk);
  if(stat != CUSPARSE_STATUS_SUCCESS)
    {
      printf("CUDA mat HYBMV  error error %s\n",cublasGetError((cublasStatus_t)stat));exit(1);
    }
  return;
}
