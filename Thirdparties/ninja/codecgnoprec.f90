subroutine GPUCGNOPREC(ND,V,maxiter,amatr,ia,ja,bb,xx)
  use gpuvarcg
  use gpumanager
  use def_master, only : kfl_paral
  use def_solver, only : resin,resfi,resi1,resi2,solve_sol
  implicit none
  integer*4 :: ND,V,maxiter
  integer*4 :: ia(*),ja(*)
  real*8    :: amatr(*),bb(*),xx(*)
  integer*8 :: tsize
  integer*4 :: NNZB,N,NV,NDOT
  real*8    :: a,b,na,gammaold,r1,gamma,delta,deno
  integer*4 :: NVOWN,off,exlen
  real*8    :: alp,bet,nalp
  real*8    :: invb,rnorm,eps
  integer*8 :: free,total,usepergpu,error,emptystr

  if(INOTMASTER) then
     exlen = (ND-npoi1) * V
     off   = (npoi1) * V 
     alp   =  1.0_rp
     bet   =  0.0_rp
     nalp  = -1*alp
     NNZB  = ia(ND+1) -1
     N     = NNZB * V * V
     NV    = ND * V
     if (npoi3 .ne. -1) then
        NDOT  = npoi3 * V
     else
        NDOT = ND * V
     end if
  end if

  emptystr = 0
  call gpumeminitiate(ND,V)        
  if( INOTMASTER ) then
     if (ti == 1) then
        call CUBLASHANDLE(handlebl)
        call CUSPARSEHANDLE(handlesp)
        ti = ti + 1
     end if
     tsize = N
     tsize = tsize*8
     call GPUMEMASSIGN(d_val,tsize)
!     call CPUMEMREGISTER(amatr,tsize)

     tsize = NV
     tsize = tsize*8
     call GPUMEMASSIGN(d_x,tsize)     
     call GPUMEMASSIGN(d_r,tsize)
     call GPUMEMASSIGN(d_p,tsize)
     call GPUMEMASSIGN(d_u,tsize)
     call GPUMEMASSIGN(d_s,tsize)
!     call CPUMEMREGISTER(xx,tsize)
!     call CPUMEMREGISTER(bb,tsize)

     tsize = ND+1
     tsize = tsize*4
     call GPUMEMASSIGN(d_row,tsize)
!     call CPUMEMREGISTER(ia,tsize)

     tsize = NNZB
     tsize = tsize*4
     call GPUMEMASSIGN(d_col,NNZB*4)
!     call CPUMEMREGISTER(ja,tsize)

     tsize = N
     tsize = tsize*8

     call MEMCPYTOGPU(d_val,amatr,tsize,"CG matrix") ! Matrix copy
     call MEMCPYTOGPU(d_col,ja,NNZB*4  ,"CG col")    ! Graph IA opy
     call MEMCPYTOGPU(d_row,ia,(ND+1)*4,"CG row")    ! Graph JA copy

     tsize=NV
     tsize = tsize*8
     call MEMCPYTOGPU(d_x,xx,tsize,"CG guess")       ! Guess copy
     call MEMCPYTOGPU(d_r,bb,tsize,"CG rhs")         ! RHS copy
     call PREPARECSR(handlesp,ND,ND,NNZB,V,d_val,d_row,d_col,emptystr)
     !
     ! calculate deno and invb
     !
!@   call PREPAREPRECON(amatr,ia,ja,d_val,d_row,d_col,ND,V) !@
!@   call GPUPRECON(d_r,d_s,NV,emptystr)       |@
     call CUDADCOPY(handlebl,NV,d_r,1,d_s,1)   !@
   
     call CUDADDOT(handlebl,NDOT,d_r,1,d_s,1,deno)
     call CUDADDOT(handlebl,NDOT,d_r,1,d_r,1,invb)
  end if

  call PAR_SUM(deno,'IN MY CODE')
  call PAR_SUM(invb,'IN MY CODE')

  deno = sqrt(deno)
  invb = 1.0_rp/sqrt(invb)

  if( INOTMASTER ) then
     !
     ! Matrix multiply
     !

     call GPUDMV(handlesp,ND,ND,NNZB,V,alp,&
          d_val,d_row,d_col,d_x,bet,d_s,emptystr)

     !start conjugate gradient

     call CUDADAXPY(handlebl,NV,nalp,d_s,1,d_r,1)
     call CUDADCOPY(handlebl,NV,d_r,1,d_p,1)   !@
   ! call GPUPRECON(d_r,d_p,NV,emptystr)       !@


     call CUDADCOPY(handlebl,NV,d_p,1,d_u,1)
     call CUDADDOT(handlebl,NDOT,d_r,1,d_u,1,gamma)
     call CUDADDOT(handlebl,NDOT,d_r,1,d_r,1,rnorm)
  end if

  call PAR_SUM(gamma,'IN MY CODE')
  call PAR_SUM(rnorm,'IN MY CODE')

  rnorm = sqrt(rnorm)
  r1 = sqrt(gamma)
  solve_sol(1) % iters = 1
  !-------------------ADAPTIVE RESIDUAL--------------------

  solve_sol(1) % resi2 = rnorm * invb     !normalized residual
  solve_sol(1) % resin = r1 / deno     !preconditioned residual
  resin = solve_sol(1) % resin
  if( solve_sol(1) % kfl_adres == 1 ) then
     eps = max( solve_sol(1) % resin * solve_sol(1) % adres , solve_sol(1) % solmi )
  else
     eps = solve_sol(1) % solco
  end if
  tolsta = eps*deno

  do while (r1 > tolsta .and.   solve_sol(1) % iters <= maxiter)
     if( INOTMASTER ) then

        call GPUDMV(handlesp,ND,ND,NNZB,V,alp,&
             d_val,d_row,d_col,d_p,bet,d_s,emptystr)
        !
        ! dot = d_Ax . d_p
        !   
        call CUDADDOT(handlebl,NDOT,d_s,1,d_p,1,delta)
     end if

     call PAR_SUM(delta,'IN MY CODE')
     a  = gamma/delta
     gammaold = gamma
     if( INOTMASTER ) then
        !    PRINT *,dot,a
        ! 
        ! d_x = d_x + a * d_p 
        !
        call CUDADAXPY(handlebl,NV,a,d_p,1,d_x,1)
        !
        ! d_r = d_r + na * d_s
        !
        na = -1.0_rp * a
        call CUDADAXPY(handlebl,NV,na,d_s,1,d_r,1)
        !
        ! ui+1 = M-1 ri+1
        !
          call CUDADCOPY(handlebl,NV,d_r,1,d_u,1)   !@
       !   call GPUPRECON(d_r,d_u,NV,emptystr)      !@


        call CUDADDOT(handlebl,NDOT,d_r,1,d_u,1,gamma)
     end if
     call PAR_SUM(gamma,'IN MY CODE')     
     b = gamma/gammaold
     if( INOTMASTER ) then
        call CUDADSCAL(handlebl,NV,b,d_p,1)
        call CUDADAXPY(handlebl,NV,alp,d_u,1,d_p,1)
     end if

     resi2 = resi1

     r1 = sqrt(gamma)
     resi1 = r1/deno

     solve_sol(1) % iters = solve_sol(1) % iters + 1
  enddo
  !------SUbmitting outputs--------------
  resfi = resi1
  solve_sol(1) % iters =   solve_sol(1) % iters - 1 

  if(INOTMASTER ) then
     call MEMCPYFROMGPU(xx,d_x,tsize,"CG unk") ! Unknown copied back to CPU
!     call CPUMEMUNREGISTER(amatr)
!     call CPUMEMUNREGISTER(ia)
!     call CPUMEMUNREGISTER(ja)
!     call CPUMEMUNREGISTER(bb)
!     call CPUMEMUNREGISTER(xx)
  end if
  call destroyprecon()
  call gpumemreset()


end subroutine GPUCGNOPREC
