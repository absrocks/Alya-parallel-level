#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <cublas_v2.h>
#include "cudakernel.h"
#define threads 256
#define devptr_t size_t

extern "C" void applypermutation_(devptr_t* mat,devptr_t* tempmat,devptr_t* map,int *NNZB)
{
  cudaError_t stat;
  int N = *NNZB;
  applypermutationkernel<<<((int)(N/threads) + 1),threads>>>((double*)(*mat),(double*)(*tempmat),(int*)(*map),N);
  stat = cudaGetLastError();
  if(stat != cudaSuccess)
    {
      printf("FAiled launch apply permutation kernel %s\n",cudaGetErrorString(stat));exit(1);
    }

  return;

}

extern "C" void applyrevpermutation_(devptr_t* mat,devptr_t* tempmat,devptr_t* map,int *NNZB)
{
  cudaError_t stat;
  int N = *NNZB;
  applyrevpermutationkernel<<<((int)(N/threads) + 1),threads>>>((double*)(*mat),(double*)(*tempmat),(int*)(*map),N);
  stat = cudaGetLastError();
  if(stat != cudaSuccess)
    {
      printf("FAiled launch apply permutation kernel %s\n",cudaGetErrorString(stat));exit(1);
    }

  return;

}
extern "C" void genbdrrowidx_(devptr_t *ia,devptr_t *iaa,int *npoi,int *ND,int *val,devptr_t *str)
{
  cudaError_t stat;
  int N = *ND;
  int off = *npoi;
  int *x = (int *)(*ia + off*4);
  int *y = (int *)(*iaa);
  genbdrrowidxker<<<(int)((N-off+1)/threads)+1,threads>>>(x,y,N-off+1,*val);
  stat = cudaGetLastError();
  if(stat != cudaSuccess)
    {
      printf("FAiled launch BDR rowidx kernel %s\n",cudaGetErrorString(stat));exit(1);
    }

  return;
}

extern "C" void ddiagonalgen_(int *N,devptr_t *mat,devptr_t *ia,devptr_t *ja,devptr_t *diag,devptr_t *str)
{
  cudaError_t stat;
  diagonalpregen<<<(int)((*N)/threads)+1,threads>>>(*N,(double*)(*mat),(int*)(*ia),(int*)(*ja),(double*)(*diag));
  stat = cudaGetLastError();
  if(stat != cudaSuccess)
    {
      printf("FAiled launch diag genkernel %s\n",cudaGetErrorString(stat));exit(1);
    }

  return;
}

extern "C" void ddivbyelem_(devptr_t *xx,devptr_t *yy,devptr_t *out,int *len,devptr_t *str)
{
  cudaError_t stat;
  double *y = (double *)(*yy);
  double *x = (double *)(*xx);
  double *z = (double *)(*out);
  int N = *len;
  if( *str == 0)
    {
      dmulbyelemker<<<(int)(N/threads)+1,threads>>>(x,y,z,N);
    }
  else
    {
      dmulbyelemker<<<(int)(N/threads)+1,threads,0,(cudaStream_t)(*str)>>>(x,y,z,N);
    }
  stat = cudaGetLastError();
  if(stat != cudaSuccess)
    {
      printf("FAiled launch ddiv by element kernel %s\n",cudaGetErrorString(stat));exit(1);
    }
  return;
}

extern "C" void dmultiplybyelem_(devptr_t *xx,devptr_t *yy,devptr_t *out,int *len,devptr_t *str)
{
  cudaError_t stat;
  double *y = (double *)(*yy);
  double *x = (double *)(*xx);
  double *z = (double *)(*out);
  int N = *len;
  if( *str == 0)
    {
      dmultiplybyelemker<<<(int)(N/threads)+1,threads>>>(x,y,z,N);
    }
  else
    {
      dmultiplybyelemker<<<(int)(N/threads)+1,threads,0,(cudaStream_t)(*str)>>>(x,y,z,N);
    }
  stat = cudaGetLastError();
  if(stat != cudaSuccess)
    {
      printf("FAiled launch dmultiply by element kernel %s\n",cudaGetErrorString(stat));exit(1);
    }
  return;
}


extern "C" void linearcombogpu_(devptr_t *X,devptr_t *V,devptr_t *Y,int *M,int *N,devptr_t *str)
{
  cudaError_t stat;
  double *y = (double *)(*Y);
  double *x = (double *)(*X);
  double *v = (double *)(*V);
  int n = *N;
  int m = *M;
  if( *str == 0)
    {
      linearcombo<<<(int)(m/threads)+1,threads,n*sizeof(double)>>>(x,v,y,m,n);
      //linearcombostatic<<<(int)(m/threads)+1,threads>>>(x,v,y,m,n);
    }
  else
    {
      linearcombo<<<(int)(m/threads)+1,threads,n*sizeof(double),(cudaStream_t)(*str)>>>(x,v,y,m,n);
    }
  stat = cudaGetLastError();
  if(stat != cudaSuccess)
    {
      printf("FAiled launch linear combo kernel %s\n",cudaGetErrorString(stat));
      exit(1);
    }
  return;
}


extern "C" void dbsrmv_(int *mb,int *nb,int *nnzb,devptr_t *bsrValA,devptr_t *bsrRowPtrA,devptr_t *bsrColIndA,int *blockdim,devptr_t *x,devptr_t *y,devptr_t *str)
{
  double *rhs = (double *)(*x);
  double *unk = (double *)(*y);
  double *A = (double *)(*bsrValA);
  int *r = (int *)(*bsrRowPtrA);
  int *c = (int *)(*bsrColIndA);
  cudaError_t stat;
  //  if( *str == 0)
  //  {
  dbsrmvkernel<<<(int)(*mb/threads)+1,threads>>>(A,r,c,rhs,unk,*mb,*nb,*blockdim);
  //  }
  // else
  // {
  //  linearcombo<<<(int)(m/threads)+1,threads,0,(cudaStream_t)(*str)>>>(x,v,y,m,n);
  // }
  stat = cudaGetLastError();
  if(stat != cudaSuccess)
    {
      printf("FAiled launch linear combo kernel %s\n",cudaGetErrorString(stat));
      exit(1);
    }
  return;
}

extern "C" void csrdefl_(devptr_t *mat,devptr_t *row,devptr_t *col,int *NGRP,int *ND,int *V,devptr_t *grpmatrix,devptr_t *iagro,devptr_t* jagro,devptr_t *LGRP,devptr_t *str)
{
  double *A = (double*)(*mat);
  int *R = (int*)(*row);
  int *C = (int*)(*col);
  int *lgrp = (int*)(*LGRP);
  int *ig = (int*)(*iagro);
  int *jg = (int*)(*jagro);
  double *GPM = (double*)(*grpmatrix);
  cudaError_t stat;
  csrdeflkernel<<<(int)(*ND/threads)+1,threads>>>(A,R,C,*NGRP,*ND,*V,GPM,ig,jg,lgrp);
  stat = cudaGetLastError();
  if(stat != cudaSuccess)
    {
      printf("FAiled launch CSR DEFL kernel %s\n",cudaGetErrorString(stat));
      exit(1);
    }
  return;
}

extern "C" void cudawtvec_(devptr_t *SMALL,devptr_t *BIG,devptr_t* LGRP,int *npoi1,int *npoi2,int *npoi3,int *V,int *NGRP,devptr_t* str)
{
  double *s = (double*)(*SMALL);
  double *b = (double*)(*BIG);
  int *lgrp = (int*)(*LGRP);
  cudaError_t stat;
  int launchsz =0;
  if( *npoi2 <= *npoi3 )
    {
      launchsz = *npoi3;
      launchsz *= *V;
    }
  else
    {
      launchsz = *npoi1;
      launchsz *= *V;
    }
  wtveckernel<<<(int)(launchsz/threads)+1,threads,(*NGRP)*sizeof(double)>>>(s,b,lgrp,*npoi1,*npoi2,*npoi3,*V,*NGRP);
  stat = cudaGetLastError();
  if(stat != cudaSuccess)
    {
      printf("FAiled launch WTVEC kernel %s\n",cudaGetErrorString(stat));
      exit(1);
    }
  return;
}

extern "C" void cudawvec_(devptr_t *BIG,devptr_t *SMALL,devptr_t* LGRP,int *ND,int *V,devptr_t* str)
{
  double *s = (double*)(*SMALL);
  double *b = (double*)(*BIG);
  int *lgrp = (int*)(*LGRP);
  cudaError_t stat;
  wveckernel<<<(int)(*ND/threads)+1,threads>>>(b,s,lgrp,*ND,*V);
  stat = cudaGetLastError();
  if(stat != cudaSuccess)
    {
      printf("FAiled launch WVEC kernel %s\n",cudaGetErrorString(stat));
      exit(1);
    }
  return;
}

extern "C" void initialize_boundaries_(int *npoi1,int *v,int *bdim,devptr_t *perm,devptr_t *send,devptr_t *recv,devptr_t *bdrv)
{
  int V = *v;
  int dim = *bdim;
  cudaError_t stat;
  initboundaries<<<(int)((V*dim)/threads)+1,threads>>>(*npoi1,V,dim,(int*)(*perm),(double*)(*send),(double*)(*recv),(double*)(*bdrv));
  stat = cudaGetLastError();
  if(stat != cudaSuccess)
    {
      printf("FAiled launch init boundaries kernel %s\n",cudaGetErrorString(stat));
      exit(1);
    }
  return;
}

extern "C" void permute_sum_(int *npoi1,int *v,int *bdim,devptr_t *perm,devptr_t *send,devptr_t *recv,devptr_t *bdrv)
{
  int V = *v;
  int dim = *bdim;
  cudaError_t stat;
  sumboundaries<<<(int)((V*dim)/threads)+1,threads>>>(*npoi1,V,dim,(int*)(*perm),(double*)(*send),(double*)(*recv),(double*)(*bdrv));
  stat = cudaGetLastError();
  if(stat != cudaSuccess)
    {
      printf("FAiled launch init boundaries kernel %s\n",cudaGetErrorString(stat));
      exit(1);
    }
  return;
}

extern "C" void lightspmvvector_(size_t *rc,int *ND,size_t *row,size_t *col,size_t* val,size_t* x,size_t *y,devptr_t *str,int *blocklnch)
{
  cudaError_t stat;
  if(*str == 0)
    {
      csr64DynamicVectorker<<<(*blocklnch),1024,128*2*sizeof(uint32_t)>>>(32,(uint32_t*)(*rc),(uint32_t)(*ND),(uint32_t*)(*row),(uint32_t*)(*col),(double*)(*val),(double*)(*x),(double*)(*y));
    }
  else
    {
      csr64DynamicVectorker<<<(*blocklnch),1024,128*2*sizeof(uint32_t),(cudaStream_t)(*str)>>>(32,(uint32_t*)(*rc),(uint32_t)(*ND),(uint32_t*)(*row),(uint32_t*)(*col),(double*)(*val),(double*)(*x),(double*)(*y));
    }
  stat = cudaGetLastError();
  if(stat != cudaSuccess)
    {
      printf("FAiled launch light SPMV kernel %s\n",cudaGetErrorString(stat));
      exit(1);
    }
  return;
}

extern "C" void lightspmvwarp_(size_t *rc,int *ND,size_t *row,size_t *col,size_t* val,size_t* x,size_t *y,devptr_t *str,int *blocklnch)
{
  cudaError_t stat;
  if(*str == 0)
    {
      csr64DynamicWarpker<<<(*blocklnch),1024,128*2*sizeof(uint32_t)>>>(32,(uint32_t*)(*rc),(uint32_t)(*ND),(uint32_t*)(*row),(uint32_t*)(*col),(double*)(*val),(double*)(*x),(double*)(*y));
    }
  else
    {
      csr64DynamicWarpker<<<(*blocklnch),1024,128*2*sizeof(uint32_t),(cudaStream_t)(*str)>>>(32,(uint32_t*)(*rc),(uint32_t)(*ND),(uint32_t*)(*row),(uint32_t*)(*col),(double*)(*val),(double*)(*x),(double*)(*y));
    }  
  stat = cudaGetLastError();
  if(stat != cudaSuccess)
    {
      printf("FAiled launch light SPMV kernel %s\n",cudaGetErrorString(stat));
      exit(1);
    }
  return;
}

extern "C" void tobasezero_(size_t *row,size_t *col,size_t *rowcsr,size_t *colcsr,int *ND,int *MD,int *NNZB,devptr_t *str)
{
  cudaError_t stat;
  int big;
  
  big = MAX(*NNZB,*ND);
  if(*str == 0)
    {  
      tobasezeroker<<<((int)(big/threads)+1),threads>>>((int*)(*row),(int*)(*col),(int*)(*rowcsr),(int*)(*colcsr),*ND,*MD,*NNZB);
    }
  else
    {
      tobasezeroker<<<((int)(big/threads)+1),threads,0,(cudaStream_t)(*str)>>>((int*)(*row),(int*)(*col),(int*)(*rowcsr),(int*)(*colcsr),*ND,*MD,*NNZB);
    }
      stat = cudaGetLastError();
      if(stat != cudaSuccess)
    {
      printf("FAiled launch light SPMV kernel %s\n",cudaGetErrorString(stat));
      exit(1);
    }
  return;
}

extern "C" void find_bandwidth_(devptr_t *rowoff,devptr_t *colarr,devptr_t *max, int *ND)
{
  cudaError_t stat;
  find_bandwidth_kernel<<<32,threads,threads*sizeof(int)>>>((int*)(*rowoff),(int*)(*colarr),(int*)(*max),*ND);
  stat = cudaGetLastError();
  if(stat != cudaSuccess)
    {
      printf("FAiled launch bandwidth subtract kernel %s\n",cudaGetErrorString(stat));
      exit(1);
    }
  return;
}
