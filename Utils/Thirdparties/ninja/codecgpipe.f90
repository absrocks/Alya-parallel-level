subroutine GPUPIPECG(ND,V,maxiter,amatr,ia,ja,bb,xx)
  use gpuvarpipecg
  use gpumanager
  use def_master, only : kfl_paral
  use def_solver, only : resin,resfi,resi1,resi2,solve_sol
  implicit none
  integer*4 :: ND,V,maxiter
  integer*4 :: ia(*),ja(*)
  real*8 :: amatr(*),bb(*),xx(*)
  integer*4 :: NNZB,N,NV,NDOT
  integer*8 :: tsize,emptystr
  real*8 :: a,b,na,gammapre,r1,gamma,delta,deno,alpha,alphapre,beta,nalpha
  integer*4 :: NVOWN,off,exlen
  real*8 ::alp,bet,nalp
  real*8 :: invb,rnorm,eps

  if(INOTMASTER) then
     exlen = (ND-npoi1) * V
     off   =  (npoi1) * V 
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

     if (tiiii == 1) then
        call GPUSTREAMCREATE(str1)
        call GPUSTREAMCREATE(str2)
        call CUBLASHANDLE(handlebl1)
        call CUBLASHANDLE(handlebl2)
        call CUSPARSEHANDLE(handlesp)
        call CUBLASATTACHSTREAM(handlebl1,str1)
        call CUBLASATTACHSTREAM(handlebl2,str2)
        call CUSPARSEATTACHSTREAM(handlesp,str1)
        tiiii = tiiii + 1
     end if

     tsize = N
     tsize = tsize*8     
     call GPUMEMASSIGN(d_val,tsize)
     call CPUMEMREGISTER(amatr,tsize)
     

     tsize = ND+1
     tsize = tsize*4
     call GPUMEMASSIGN(d_row,tsize)
     call CPUMEMREGISTER(ia,tsize)
     
     tsize = NNZB
     tsize = tsize*4
     call GPUMEMASSIGN(d_col,NNZB*4)
     call CPUMEMREGISTER(ja,tsize)
     
     tsize = NV
     tsize = tsize*8  
     call GPUMEMASSIGN(d_x,tsize)
     call GPUMEMASSIGN(d_r,tsize)
     call GPUMEMASSIGN(d_p,tsize)
     call GPUMEMASSIGN(d_w,tsize)
     call GPUMEMASSIGN(d_z,tsize)
     call GPUMEMASSIGN(d_q,tsize)
     call GPUMEMASSIGN(d_u,tsize)
     call GPUMEMASSIGN(d_m,tsize)
     call GPUMEMASSIGN(d_n,tsize)
     call GPUMEMASSIGN(d_s,tsize)
     call CPUMEMREGISTER(xx,tsize)
     call CPUMEMREGISTER(bb,tsize)


     tsize = N
     tsize = tsize*8     
     call MEMCPYTOGPUASYNC(d_val,amatr,tsize,str2,"pipecg val")
     call MEMCPYTOGPUASYNC(d_col,ja,NNZB*4,str1,"pipecg col")
     call MEMCPYTOGPUASYNC(d_row,ia,(ND+1)*4,str1,"pipecg row")
     tsize = NV
     tsize = tsize*8     
     call MEMCPYTOGPUASYNC(d_x,xx,tsize,str2,"pipecg xx")
     call MEMCPYTOGPUASYNC(d_r,bb,tsize,str2,"pipecg rr")
     
     !
     ! calculate deno and invb
     !
     
     call PREPARECSR(handlesp,ND,ND,NNZB,V,d_val,d_row,d_col,str2)

     call PREPAREPRECON(amatr,ia,ja,d_val,d_row,d_col,ND,V)
     call GPUPRECON(d_r,d_u,NV,str1)

     call CUDADDOT(handlebl1,NDOT,d_r,1,d_u,1,deno)
     call CUDADDOT(handlebl2,NDOT,d_r,1,d_r,1,invb)
  end if
  call GPUSYNC()
  call PAR_SUM(deno,'IN MY CODE')
  call PAR_SUM(invb,'IN MY CODE')

  deno = sqrt(deno)
  invb = 1/sqrt(invb)
  
  if( INOTMASTER ) then
     !
     ! Matrix multiply
     !

     call GPUDMV(handlesp,ND,ND,NNZB,V,alp,&
          d_val,d_row,d_col,d_x,bet,d_u,str1)
     !start conjugate gradient

     call CUDADAXPY(handlebl1,NV,nalp,d_u,1,d_r,1)
     call CUDADDOT(handlebl2,NDOT,d_r,1,d_r,1,rnorm)
     call GPUPRECON(d_r,d_u,NV,str1)
     call GPUDMV(handlesp,ND,ND,NNZB,V,alp,&
          d_val,d_row,d_col,d_u,bet,d_w,str1)
     call CUDADDOT(handlebl2,NDOT,d_r,1,d_u,1,gamma)
     
  end if
  call GPUSYNC()
  call PAR_SUM(gamma,'IN MY CODE')
  call PAR_SUM(rnorm,'IN MY CODE')
  
  rnorm = sqrt(rnorm)
  r1 = sqrt(gamma)
  solve_sol(1) % iters = 0
  !-------------------ADAPTIVE RESIDUAL--------------------
  !solve_sol(1) % reni2 = rnorm            !non-normalized residual
  solve_sol(1) % resi2 = rnorm * invb     !normalized residual
  solve_sol(1) % resin = r1 / deno     !preconditioned residual
  resin = solve_sol(1) % resin
  if( solve_sol(1) % kfl_adres == 1 ) then
     eps = max( solve_sol(1) % resin * solve_sol(1) % adres , solve_sol(1) % solmi )
  else
     eps = solve_sol(1) % solco
  end if
  tolsta = eps*deno

  do while (r1 > tolsta .and. solve_sol(1) % iters <= maxiter)
     if( INOTMASTER ) then

        call GPUPRECON(d_w,d_m,NV,str2)
        call CUDADDOT(handlebl1,NDOT,d_r,1,d_u,1,gamma)
        call CUDADDOT(handlebl1,NDOT,d_w,1,d_u,1,delta)

        call GPUDMV(handlesp,ND,ND,NNZB,V,alp,&
             d_val,d_row,d_col,d_m,bet,d_n,str2)
     
     end if
     call GPUSYNC()
     call PAR_SUM(gamma,'IN MY CODE')
     call PAR_SUM(delta,'IN MY CODE')
     
     if (solve_sol(1) % iters > 1) then
        beta = gamma / gammapre
        alpha = gamma/(delta-((beta*gamma)/alphapre))
        nalpha = -1*alpha
     else
        beta=0
        alpha=gamma/delta
        nalpha = -1*alpha
     end if

     if( INOTMASTER ) then
        call CUDADSCAL(handlebl1,NV,beta,d_z,1)
        call CUDADAXPY(handlebl1,NV,alp,d_n,1,d_z,1)

        call CUDADSCAL(handlebl2,NV,beta,d_q,1)
        call CUDADAXPY(handlebl2,NV,alp,d_m,1,d_q,1)

        call CUDADSCAL(handlebl1,NV,beta,d_s,1)
        call CUDADAXPY(handlebl1,NV,alp,d_w,1,d_s,1)

        call CUDADSCAL(handlebl2,NV,beta,d_p,1)
        call CUDADAXPY(handlebl2,NV,alp,d_u,1,d_p,1)

        call CUDADAXPY(handlebl1,NV,nalpha,d_s,1,d_r,1)
        call CUDADAXPY(handlebl2,NV,alpha,d_p,1,d_x,1)
        call CUDADAXPY(handlebl2,NV,nalpha,d_q,1,d_u,1)
        call CUDADAXPY(handlebl1,NV,nalpha,d_z,1,d_w,1)
     end if
     
     alphapre=alpha
     gammapre=gamma
     resi2 = resi1
     
     r1 = sqrt(gamma)
     resi1 = r1/deno

     solve_sol(1) % iters = solve_sol(1) % iters + 1
     call GPUSYNC()
  enddo
  !------SUbmitting outputs--------------
  resfi = resi1
  solve_sol(1) % iters = solve_sol(1) % iters -1 
  if(INOTMASTER ) call MEMCPYFROMGPU(xx,d_x,tsize,"CG unk")

  call destroyprecon()
  if (INOTMASTER) then
     call CPUMEMUNREGISTER(amatr)
     call CPUMEMUNREGISTER(ia)
     call CPUMEMUNREGISTER(ja)
     call CPUMEMUNREGISTER(bb)
     call CPUMEMUNREGISTER(xx)
  end if
  call gpumemreset()
  
end subroutine GPUPIPECG
