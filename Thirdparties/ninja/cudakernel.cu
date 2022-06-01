#include <stdio.h>
#include <stdlib.h>
#include "cudakernel.h"
#include <stdint.h>

const char* cublasGetError(cublasStatus_t error)
{
  switch (error)
    {
    case CUBLAS_STATUS_SUCCESS:
      return "CUBLAS_STATUS_SUCCESS";

    case CUBLAS_STATUS_NOT_INITIALIZED:
      return "CUBLAS_STATUS_NOT_INITIALIZED";

    case CUBLAS_STATUS_ALLOC_FAILED:
      return "CUBLAS_STATUS_ALLOC_FAILED";

    case CUBLAS_STATUS_INVALID_VALUE:
      return "CUBLAS_STATUS_INVALID_VALUE";

    case CUBLAS_STATUS_ARCH_MISMATCH:
      return "CUBLAS_STATUS_ARCH_MISMATCH";

    case CUBLAS_STATUS_MAPPING_ERROR:
      return "CUBLAS_STATUS_MAPPING_ERROR";

    case CUBLAS_STATUS_EXECUTION_FAILED:
      return "CUBLAS_STATUS_EXECUTION_FAILED";

    case CUBLAS_STATUS_INTERNAL_ERROR:
      return "CUBLAS_STATUS_INTERNAL_ERROR";
    }

  return "<unknown>";
}

const char* cusolverGetError(cusolverStatus_t error)
{
  switch (error)
    {
    case CUSOLVER_STATUS_SUCCESS:
      return "CUSOLVER_STATUS_SUCCESS";

    case CUSOLVER_STATUS_NOT_INITIALIZED:
      return "CUSOLVER_STATUS_NOT_INITIALIZED";

    case CUSOLVER_STATUS_ALLOC_FAILED:
      return "CUSOLVER_STATUS_ALLOC_FAILED";

    case CUSOLVER_STATUS_INVALID_VALUE:
      return "CUSOLVER_STATUS_INVALID_VALUE";

    case CUSOLVER_STATUS_ARCH_MISMATCH:
      return "CUSOLVER_STATUS_ARCH_MISMATCH";

    case CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
      return "CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED";

    case CUSOLVER_STATUS_EXECUTION_FAILED:
      return "CUSOLVER_STATUS_EXECUTION_FAILED";

    case CUSOLVER_STATUS_INTERNAL_ERROR:
      return "CUSOLVER_STATUS_INTERNAL_ERROR";
    }

  return "<unknown>";
}

__global__ void genbdrrowidxker(int *x,int *y,int len,int val)
{
  int tx = threadIdx.x + blockIdx.x*blockDim.x;
  if(tx<len)
    {
      y[tx]=x[tx]-val;
    }
  return;
}

__global__ void dmulbyelemker(double *x,double *y,double *out,int N)
{
  int tx = threadIdx.x + blockIdx.x*blockDim.x;
  if(tx<N)
    {
      out[tx]=y[tx]/x[tx];
    }
  return;
}
__global__ void dmultiplybyelemker(double *x,double *y,double *out,int N)
{
  int tx = threadIdx.x + blockIdx.x*blockDim.x;
  if(tx<N)
    {
      out[tx]=y[tx]*x[tx];
    }
  return;
}

__global__ void linearcombo(double *x,double *v,double *y,int m,int n)
{
  extern __shared__ double temp[];
  int tx = threadIdx.x + blockIdx.x*blockDim.x;
  double out;
  if(tx<m)
    {
      out = x[tx];
      for(int i=threadIdx.x; i < n; i += blockDim.x) 
	{
	  temp[i] = y[ i ];
	}
      __syncthreads();
      for(int i=0;i<n;i++)
	{
	  out += v[tx + i*m]*temp[i];
	}
      x[tx] = out;
    }
  return;
}

__global__ void linearcomboprecon(double *x,double *v,double *y,double *pre,int m,int n)
{
  extern __shared__ double temp[];
  int tx = threadIdx.x + blockIdx.x*blockDim.x;
  double out;
  if(tx<m)
    {
      out = x[tx];
      for(int i=threadIdx.x; i < n; i += blockDim.x) 
	{
	  temp[i] = y[ i ];
	}
      __syncthreads();
      for(int i=0;i<n;i++)
	{
	  out += v[tx + i*m]*temp[i];
	}
      x[tx] = out;
    }
  return;
}

__global__ void linearcombostatic(double *x,double *v,double *y,int m,int n)
{
  int tx = threadIdx.x + blockIdx.x*blockDim.x;
  double out;
  if(tx<m)
    {
      out = 0;
      for(int i=0;i<n;i++)
	{
	  out += v[tx + i*m]*y[i];
	}
      x[tx] = x[tx] + out;
    }
  return;
}

__global__ void dbsrmvkernel(double *matrix,int *row,int *col,double *rhs,double *unk,int m,int n,int dim)
{
  int tx = threadIdx.x + blockIdx.x*blockDim.x;
  int blkst;
  int blkend;
  int blkoff;
  int colid;
  double temp;
  double raux;
  if (tx < m)
    {
      blkst = row[tx] - 1;
      blkend = row[tx+1] - 1;
      blkoff = blkst*dim*dim;
      for(int k=blkst;k<blkend;k++)
	{
	  colid = col[k]-1;
	  for(int i=0;i<dim;i++)
	    {
	      temp=0;
	      raux = rhs[colid*dim + i];
	      for(int j=0;j<dim;j++)
		{
		  temp += matrix[blkoff + j*dim + i]*raux;
		}
	      unk[colid*dim + i] = temp;
	    }
	}
    }
  return;
}
__global__ void csrdeflkernel(double *matrix,int *row,int *col,int NGRP,int ND,int V,double *grpmatrix,int *iagro,int *jagro,int *lgrp)
{
  int tx = threadIdx.x + blockIdx.x*blockDim.x;
  int igr,jgr,ty;
  int tmp1;
  if(tx < ND)
    {
      if(lgrp[tx] >= 0)
	{
	  igr = lgrp[tx]-1;
	  for(int i = row[tx]-1;i<row[tx+1]-1;i++)
	    {
	      ty = col[i] -1;
	      if( lgrp[ty] >= 0 )
		{
		  jgr = lgrp[ty]-1;
		  tmp1 = iagro[igr]-1;
		  while ( jagro[tmp1]-1 != jgr )
		    {
		      tmp1++;
		    }
		  for(int j=0;j<V*V;j++)
		    {
		      grpmatrix[tmp1*V*V +  j ] += matrix[ i*V*V +j ];
		    }
		}
	    }
	}
    }
  return;
}

__device__ __forceinline__ double myatomicAdd(double* address, double val)
{
  unsigned long long int* address_as_ull = (unsigned long long int*) address;
  unsigned long long int old = *address_as_ull, assumed;
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
}

__global__ void wtveckernel(double *small,double *big,int *lgrp,int npoi1,int npoi2,int npoi3,int V,int smallsz)
{
  extern __shared__ double temp[];
  int tx = threadIdx.x + blockIdx.x*blockDim.x;
  int lgp;
  for(int i=threadIdx.x; i < smallsz; i += blockDim.x) 
    {
      temp[i] = 0;
    }
  __syncthreads();
  if(tx < npoi1*V || ( tx >= ((npoi2-1)*V) && tx < npoi3*V) )
    {
      lgp = lgrp[(int)(tx/V)] - 1;
      if(lgp >= 0)
	{
	  //temp[lgp*V + (tx%V)] += big[tx];
	  myatomicAdd( &temp[lgp*V+(tx%V)] , big[tx]);
	}
    }
  __syncthreads();
  for(int i=threadIdx.x; i < smallsz; i += blockDim.x) 
    {
      myatomicAdd( &small[i] , temp[i]);
    }
}
__global__ void wveckernel(double *big,double *small,int *lgrp,int ND,int V)
{
  int tx = threadIdx.x + blockIdx.x*blockDim.x;
  int lgp;
  if(tx < ND)
    {
      lgp = lgrp[tx] - 1;
      if(lgp >= 0)
	{
	  for(int j=0;j<V;j++)
	    {
	      big[tx*V + j] = small[lgp*V +j]; 
	    }
	}
      else
	{
	  for(int j=0;j<V;j++)
	    {
	      big[tx*V + j] = 0.0; 
	    }
	}
    }
}

__global__ void initboundaries(int off,int v,int bdim,int *perm,double *send,double *recv,double *bdrv)
{

 int tx = threadIdx.x + blockIdx.x*blockDim.x;	
if(v == 1)
{
 if(tx < bdim)
 {
 int ipoin = perm[tx];
 send[tx] = bdrv[ ipoin-off-1 ];
recv[tx] = 0.0;
 }
}

else
{ 
  if(tx < v*bdim)
 {
 int ipoin = perm[(int)(tx/v)];
 send[tx] = bdrv[ (ipoin-off-1)*v + tx%v ];
recv[tx] = 0.0;
 }
}
}

__global__ void sumboundaries(int off,int v,int bdim,int *perm,double *send,double *recv,double *bdrv)
{

 int tx = threadIdx.x + blockIdx.x*blockDim.x;	
 if(v == 1)
{
 if(tx < bdim)
 {
 int ipoin = perm[tx];
 myatomicAdd( &bdrv[ipoin-off-1] , recv[tx]);	
 // bdrv[ ipoin-off-1 ] += recv[tx];
 }
}
else
{
  if(tx < v*bdim)
 {
 int ipoin = perm[(int)(tx/v)];
 myatomicAdd( &bdrv[(ipoin-off-1)*v+(tx%v)] , recv[tx]);	
 //bdrv[ (ipoin-off-1)*v + tx%v ] += recv[tx];
 }
}
}

/*device functions*/
template < typename T>
__device__ inline T shfl_down_64bits(T var, int32_t srcLane,
				     int32_t width) {

  int2 a = *reinterpret_cast<int2*>(&var);

  /*exchange the data*/
  a.x = __shfl_down_sync(a.x, srcLane, width);
  a.y = __shfl_down_sync(a.y, srcLane, width);

  return *reinterpret_cast<T*>(&a);
}


__device__ inline double DOUBLE_VECTOR_GET (const double* __restrict vectorX, uint32_t index){
  return vectorX[index];
}


/*64-bit functions*/
__global__ void csr64DynamicVectorker(uint32_t THREADS_PER_VECTOR,uint32_t* __restrict cudaRowCounter, uint32_t _cudaNumRows,const uint32_t* __restrict rowOffsets, const uint32_t* __restrict colIndexValues,const double* __restrict numericalValues, const double* __restrict vectorX, double* vectorY)
{
  uint32_t i;
  double sum;
  uint32_t row;
  uint32_t rowStart, rowEnd;
  const uint32_t laneId = threadIdx.x % THREADS_PER_VECTOR; /*lane index in the vector*/
  const uint32_t vectorId = threadIdx.x / THREADS_PER_VECTOR; /*vector index in the block*/
  
  extern __shared__ volatile uint32_t space[];
  
  /*get the row index*/
  if (laneId == 0) {
    row = atomicAdd(cudaRowCounter, 1);
  }
  /*broadcast the value to other lanes from lane 0*/
  row = __shfl_sync(0xffffffff,row, 0, THREADS_PER_VECTOR);
  
  /*check the row range*/
  while (row < _cudaNumRows) {
    
    /*use two threads to fetch the row offset*/
    if (laneId < 2) {
      space[vectorId*2 + laneId] = rowOffsets[row + laneId];
    }
    rowStart = space[vectorId*2 + 0 ];
    rowEnd = space[vectorId*2  + 1 ];
    
    /*there are non-zero elements in the current row*/
    sum = 0;
    /*compute dot product*/
    if (THREADS_PER_VECTOR == 32) {
      
      /*ensure aligned memory access*/
      i = rowStart - (rowStart & (THREADS_PER_VECTOR - 1)) + laneId;
      
      /*process the unaligned part*/
      if (i >= rowStart && i < rowEnd) {
	sum += numericalValues[i] * DOUBLE_VECTOR_GET(vectorX, colIndexValues[i]);
      }
      
      /*process the aligned part*/
      for (i += THREADS_PER_VECTOR; i < rowEnd; i += THREADS_PER_VECTOR) {
	sum += numericalValues[i] * DOUBLE_VECTOR_GET(vectorX, colIndexValues[i]);
      }
    } else {
      /*regardless of the global memory access alignment*/
      for (i = rowStart + laneId; i < rowEnd; i +=
	     THREADS_PER_VECTOR) {
	sum += numericalValues[i] * DOUBLE_VECTOR_GET(vectorX, colIndexValues[i]);
      }
    }
    /*intra-vector reduction*/
    for (i = THREADS_PER_VECTOR >> 1; i > 0; i >>= 1) {
      sum += shfl_down_64bits<double>(sum, i, THREADS_PER_VECTOR);
    }
    
    /*save the results and get a new row*/
    if (laneId == 0) {
      /*save the results*/
      vectorY[row] = sum;
      
      /*get a new row index*/
      row = atomicAdd(cudaRowCounter, 1);
    }
    row = __shfl_sync(0xffffffff,row, 0, THREADS_PER_VECTOR);
  }/*while*/
}

__global__ void tobasezeroker(int *row,int *col,int *rowcsr,int* colcsr,int ND,int MD,int NNZB)
{

  int tx = threadIdx.x + blockIdx.x*blockDim.x;
  if( tx < (ND+1) )
    {
      rowcsr[tx] = row[tx]-1;
    }
  if( tx < NNZB )
    {
      colcsr[tx] = col[tx]-1;
    }

}

__global__ void csr64DynamicWarpker(uint32_t THREADS_PER_VECTOR, uint32_t* __restrict cudaRowCounter, uint32_t _cudaNumRows,const uint32_t* __restrict rowOffsets, const uint32_t* __restrict colIndexValues,const double* __restrict numericalValues, const double* __restrict vectorX, double* vectorY)
  
{
  uint32_t i;
  double sum;
  uint32_t row;
  uint32_t rowStart, rowEnd;
  const uint32_t laneId = threadIdx.x % THREADS_PER_VECTOR; /*lane index in the vector*/
  const uint32_t vectorId = threadIdx.x / THREADS_PER_VECTOR; /*vector index in the thread block*/
  const uint32_t warpLaneId = threadIdx.x & 31;	/*lane index in the warp*/
  const uint32_t warpVectorId = warpLaneId / THREADS_PER_VECTOR;	/*vector index in the warp*/
  
  extern __shared__ volatile uint32_t space[];

  /*get the row index*/
  if (warpLaneId == 0) {
    row = atomicAdd(cudaRowCounter, 32 / THREADS_PER_VECTOR);
  }
  /*broadcast the value to other threads in the same warp*/
  row = __shfl_sync(0xffffffff,row, 0) + warpVectorId;
  
  /*check the row range*/
  while (row < _cudaNumRows) {
    
    /*use two threads to fetch the row offset*/
    if (laneId < 2) {
      space[vectorId*2 + laneId] = rowOffsets[row + laneId];
    }
    rowStart = space[vectorId*2 + 0 ];
    rowEnd = space[vectorId*2  + 1 ];
    
    /*there are non-zero elements in the current row*/
    sum = 0;
    /*compute dot product*/
    if (THREADS_PER_VECTOR == 32) {
      
      /*ensure aligned memory access*/
      i = rowStart - (rowStart & (THREADS_PER_VECTOR - 1)) + laneId;
      
      /*process the unaligned part*/
      if (i >= rowStart && i < rowEnd) {
	sum += numericalValues[i] * DOUBLE_VECTOR_GET(vectorX, colIndexValues[i]);
      }
      
      /*process the aligned part*/
      for (i += THREADS_PER_VECTOR; i < rowEnd; i += THREADS_PER_VECTOR) {
	sum += numericalValues[i] * DOUBLE_VECTOR_GET(vectorX, colIndexValues[i]);
      }
    } else {
      /*regardless of the global memory access alignment*/
      for (i = rowStart + laneId; i < rowEnd; i +=
	     THREADS_PER_VECTOR) {
	sum += numericalValues[i] * DOUBLE_VECTOR_GET(vectorX, colIndexValues[i]);
      }
    }
    
    /*intra-vector reduction*/
    for (i = THREADS_PER_VECTOR >> 1; i > 0; i >>= 1) {
      sum += shfl_down_64bits<double>(sum, i, THREADS_PER_VECTOR);
    }
    
    /*save the results and get a new row*/
    if (laneId == 0) {
      /*save the results*/
      vectorY[row] = sum;
    }
    
    /*get a new row index*/
    if(warpLaneId == 0){
      row = atomicAdd(cudaRowCounter, 32 / THREADS_PER_VECTOR);
    }
    /*broadcast the value to other threads in the same warp*/
    row = __shfl_sync(0xffffffff,row, 0) + warpVectorId;
    
  }/*while*/
}

__global__ void diagonalpregen(int N,double *matrix,int *row,int *col,double *diag)
{
 int tx = threadIdx.x + blockIdx.x*blockDim.x;	
 int rst,rend;
 if(tx < N){

 rst = row[tx];
 rend = row[tx+1];
 for(int i=rst;i<rend;i++){
 if(col[i] == tx){
 diag[tx] = matrix[i];
 break;
}}}

}

__global__ void applypermutationkernel(double *matrix,double *tempmatrix,int *map,int NNZB)
{

 int tx = threadIdx.x + blockIdx.x*blockDim.x;	
 if(tx < NNZB)
   {
     matrix[tx] = tempmatrix[map[tx]];
   }
 /* __syncthreads();
 if(tx < NNZB)
   {
     matrix[tx] = tempmatrix[tx];
     }*/
}

__global__ void applyrevpermutationkernel(double *matrix,double *tempmatrix,int *map,int NNZB)
{

 int tx = threadIdx.x + blockIdx.x*blockDim.x;	
 if(tx < NNZB)
   {
     matrix[map[tx]] = tempmatrix[tx];
   }
 /* __syncthreads();
 if(tx < NNZB)
   {
     matrix[tx] = tempmatrix[tx];
     }*/
}

__global__ void find_bandwidth_kernel(int *array, int *colarray,int *max,int n)
{
  unsigned int index = threadIdx.x + blockIdx.x*blockDim.x;
  unsigned int stride = gridDim.x*blockDim.x;
  unsigned int offset = 0;

  extern __shared__ int cache[];
  int temp = INT_MIN;
  
  if(index == 0){
    *max = temp;
  }

  while(index + offset < n){
    temp = MAX(temp, colarray[ array[index+offset+1] - 1] - colarray[ array[index + offset] ]);
    offset += stride;
  }

  cache[threadIdx.x] = temp;

  __syncthreads();


  // reduction
  unsigned int i = blockDim.x/2;
  while(i != 0){
    if(threadIdx.x < i){
      cache[threadIdx.x] = MAX(cache[threadIdx.x], cache[threadIdx.x + i]);
    }
    __syncthreads();
    i /= 2;
  }

  if(threadIdx.x == 0){
    atomicMax(max, cache[0]);
  }
}
