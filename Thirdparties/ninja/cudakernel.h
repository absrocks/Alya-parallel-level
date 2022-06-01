#ifndef _CUDA_KERNELS_
#define _CUDA_KERNELS_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cublas_v2.h>
#include <stdint.h>
#include <cusolverSp.h>

#define MIN(a, b) ((a < b) ? a : b)
#define MAX(a, b) ((a > b) ? a : b)

__global__ void initboundaries(int off,int V,int bdim,int *perm,double *send,double *recv,double *bdr);
__global__ void sumboundaries(int off,int V,int bdim,int *perm,double *send,double *recv,double *bdr);
__global__ void genbdrrowidxker(int *x,int *y,int len,int val);
__global__ void dmulbyelemker(double *x,double *y,double *out,int N);
__global__ void dmultiplybyelemker(double *x,double *y,double *out,int N);
__global__ void linearcombo(double *x,double *v,double *y,int m,int n);
__global__ void linearcomboprecon(double *x,double *v,double *y,double *pre,int m,int n);
__global__ void linearcombostatic(double *x,double *v,double *y,int m,int n);
__global__ void dbsrmvkernel(double *matrix,int *row,int *col,double *rhs,double *unk,int m,int n,int dim);
__global__ void csrdeflkernel(double *matrix,int *row,int *col,int NGRP,int ND,int V,double *grpmatrix,int *iagro,int *jagro,int *lgrp);
__global__ void wtveckernel(double *small,double *big,int *lgrp,int npoi1,int npoi2,int npoi3,int V,int smallsz);
__global__ void wveckernel(double *big,double *small,int *lgrp,int ND,int V);

const char* cublasGetError(cublasStatus_t error);
const char* cusolverGetError(cusolverStatus_t error);

__global__ void csr64DynamicVectorker(uint32_t THREADS_PER_VECTOR,uint32_t* __restrict cudaRowCounter, uint32_t _cudaNumRows,const uint32_t* __restrict rowOffsets, const uint32_t* __restrict colIndexValues,const double* __restrict numericalValues, const double* __restrict vectorX, double* vectorY);


__global__ void csr64DynamicWarpker(uint32_t THREADS_PER_VECTOR,uint32_t* __restrict cudaRowCounter, uint32_t _cudaNumRows,const uint32_t* __restrict rowOffsets, const uint32_t* __restrict colIndexValues,const double* __restrict numericalValues, const double* __restrict vectorX, double* vectorY);

__global__ void tobasezeroker(int *row,int *col,int *rowcsr,int* colcsr,int ND,int MD,int NNZB);
__global__ void diagonalpregen(int N,double *matrix,int *row,int *col,double *diag);

__global__ void applypermutationkernel(double *matrix,double *tempmatrix,int *map,int NNZB);
__global__ void applyrevpermutationkernel(double *matrix,double *tempmatrix,int *map,int NNZB);

__global__ void find_bandwidth_kernel(int *array,int *colarray,int *max,int n);
#endif
