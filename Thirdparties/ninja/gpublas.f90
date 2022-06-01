MODULE GPUBLAS

  interface gpudot
     module procedure gpudoubledot,gpusingledot
  end interface gpudot

  interface gpuaxpy
     module procedure gpudoubleaxpy,gpusingleaxpy
  end interface gpuaxpy
  
  interface gpuscal
     module procedure gpudoublescal,gpusinglescal
  end interface gpuscal
  
  interface gpucopy
     module procedure gpudoublecopy,gpusinglecopy
  end interface gpucopy
  
  interface gpunrm2
     module procedure gpudoublenrm2,gpusinglenrm2
  end interface gpunrm2

contains

  subroutine gpudoubledot(handle,N,V,x,y,ans)
    implicit none
    real*8    :: ans
    integer*4 :: N,V,NV
    integer*8 :: handle
    integer*8 :: x,y
    NV = N*V
    call cudaddot(handle,NV,x,1,y,1,ans)
    
  end subroutine gpudoubledot

  subroutine gpusingledot(handle,N,V,x,y,ans)
    implicit none
    real*4    :: ans
    integer*4 :: N,V,NV
    integer*8 :: handle
    integer*8 :: x,y
    NV = N*V
    call cudaddot(handle,NV,x,1,y,1,ans)
    
  end subroutine gpusingledot

  subroutine gpudoubleaxpy(handle,N,V,alpha,x,y)
    implicit none
    real*8    :: alpha
    integer*4 :: N,V,NV
    integer*8 :: handle
    integer*8 :: x,y
    
    NV = N*V
    call CUDADAXPY(handle,NV,alpha,x,1,y,1)

  end subroutine gpudoubleaxpy

  subroutine gpusingleaxpy(handle,N,V,alpha,x,y)
    implicit none
    real*4    :: alpha
    integer*4 :: N,V,NV
    integer*8 :: handle
    integer*8 :: x,y
    
    NV = N*V
    call CUDASAXPY(handle,NV,alpha,x,1,y,1)

  end subroutine gpusingleaxpy

  subroutine gpudoublescal(handle,N,V,alpha,x)
    implicit none
    real*8    :: alpha
    integer*4 :: N,V,NV
    integer*8 :: handle
    integer*8 :: x
    
    NV = N*V
    call CUDADSCAL(handle,NV,alpha,x,1)

  end subroutine gpudoublescal

  subroutine gpusinglescal(handle,N,V,alpha,x)
    implicit none
    real*4    :: alpha
    integer*4 :: N,V,NV
    integer*8 :: handle
    integer*8 :: x
    
    NV = N*V
    call CUDASSCAL(handle,NV,alpha,x,1)

  end subroutine gpusinglescal

  subroutine gpudoublecopy(handle,N,V,x,y,dummr)
    implicit none
    integer*4 :: N,V,NV
    integer*8 :: handle
    integer*8 :: x,y
    real*8    :: dummr
    
    NV = N*V
    call CUDADCOPY(handle,NV,x,1,y,1)
    
  end subroutine gpudoublecopy

  subroutine gpusinglecopy(handle,N,V,x,y,dummr)
    implicit none
    integer*4 :: N,V,NV
    integer*8 :: handle
    integer*8 :: x,y
    real*4    :: dummr

    NV = N*V
    call CUDASCOPY(handle,NV,x,1,y,1)
    
  end subroutine gpusinglecopy

  subroutine gpudoublenrm2(handle,N,V,x,out)
    implicit none
    real*8    :: out
    integer*4 :: N,V,NV
    integer*8 :: handle
    integer*8 :: x
    
    NV = N*V
    call CUDADNRM2(handle,NV,x,1,out)
    
  end subroutine gpudoublenrm2
  
  subroutine gpusinglenrm2(handle,N,V,x,out)
    implicit none
    real*4    :: out
    integer*4 :: N,V,NV
    integer*8 :: handle
    integer*8 :: x
    
    NV = N*V
    call CUDASNRM2(handle,NV,x,1,out)
    
  end subroutine gpusinglenrm2

end MODULE GPUBLAS
