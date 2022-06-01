subroutine GPUDEFLCG(ND,V,NGRP,maxiter,amatr,ia,ja,bb,xx)
  use gpuvardeflcg
  use def_solver, only : solve_sol,resin,resfi,resi1,resi2,iters
  use mod_csrdir
  use mod_direct_solver,  only : direct_solver_partialcleaning
  use mod_direct_solver,  only : direct_solver_factorization
  use mod_direct_solver,  only : direct_solver_solution
  use mod_direct_solver,  only : direct_solver_matrix_size
  use mod_direct_solver,  only : direct_allocate_temporary_matrix
  use mod_deflated_cg,    only : matgr2,matgro,wtvect
  use gpumanager
  implicit none
  integer*4 :: ND,V,maxiter
  integer*4 :: ia(*),ja(*),acoarse_size
  real*8    :: amatr(*),bb(*),xx(*)
  integer*4 :: NNZB,N,NV,NDOT
  integer*4 :: plat,info
  integer*8 :: tsize
  real*8    :: a,b,na,gammaold,r1,gamma,delta,deno
  integer*4 :: NVOWN,off,exlen,NGRP
  real*8    :: alp,bet,nalp
  real*8    :: invb,rnorm,eps
  integer*8 :: emptystr
  integer*4,pointer :: lgrp(:)
  
  emptystr = 0
  plat = 1
  if (INOTMASTER) then
     exlen =  (ND-npoi1) * V
     off   =  npoi1 * V 
     alp   =  1.0_rp
     bet   =  0.0_rp
     nalp  = -alp
     NNZB  =  ia(ND+1) -1
     N     =  NNZB * V * V
     NV    =  ND * V
     if (npoi3 .ne. -1) then
        NDOT  = npoi3 * V
     else
        NDOT = ND * V
     end if
  end if
  lgrp => solve_sol(1) % lgrou
  call gpumeminitiate(ND,V)

  if(NGRP /= 0 ) then
     acoarse_size  = direct_solver_matrix_size(solve_sol(1) % direct_solver_Deflation)
     allocate(mu(V*NGRP))
     allocate(murhs(V*NGRP))
     allocate(grpmat(acoarse_size))
     mu = 0.0_rp 
     murhs = 0.0_rp
     grpmat = 0.0_rp
  end if

  if(NGRP /= 0 ) then
     if( solve_sol(1) % kfl_symme == 1 ) then
        call matgro(NGRP,ND,acoarse_size,V,ia,ja,amatr,grpmat)
     else
        call matgr2(NGRP,ND,acoarse_size,V,ia,ja,amatr,grpmat)
     end if
  end if

  if (INOTMASTER)  call direct_solver_factorization(solve_sol(1) % direct_solver_Deflation,grpmat)     

  if( INOTMASTER ) then
     
     if (tii == 1) then
        call GPUSTREAMCREATE(strcopy)
        call GPUSTREAMCREATE(strcomp)
        call CUBLASHANDLE(handlebl)
        call CUSPARSEHANDLE(handlesp)
        call CUBLASATTACHSTREAM(handlebl,strcopy)
        call CUSPARSEATTACHSTREAM(handlesp,strcopy)
        tii = tii + 1
     end if
     tsize = N
     tsize = tsize*8
     call GPUMEMASSIGN(d_val,tsize)
     call CPUMEMREGISTER(amatr,tsize)

     tsize = NV
     tsize = tsize*8
     call GPUMEMASSIGN(d_x,tsize)
     call GPUMEMASSIGN(d_r,tsize)
     call GPUMEMASSIGN(d_p,tsize)
     call GPUMEMASSIGN(d_u,tsize)
     call GPUMEMASSIGN(d_s,tsize)
     call CPUMEMREGISTER(bb,tsize)
     call CPUMEMREGISTER(xx,tsize)

     tsize = NGRP*V
     tsize = tsize*8
     call GPUMEMASSIGN(d_grprhs,tsize)
     call GPUMEMASSIGN(d_grplhs,tsize)
     tsize = ND
     tsize = tsize * 4
     call GPUMEMASSIGN(d_lgrp,tsize)
     call CPUMEMREGISTER(lgrp,tsize)
     tsize = NNZB
     tsize = tsize * 4
     call GPUMEMASSIGN(d_col,tsize)
     call CPUMEMREGISTER(ja,tsize)
     tsize = ND+1
     tsize = tsize * 4
     call GPUMEMASSIGN(d_row,tsize) 
     call CPUMEMREGISTER(ia,tsize)


     tsize = N
     tsize = tsize*8
     call MEMCPYTOGPU(d_val,amatr,tsize,"CG matrix")

     tsize = NNZB
     tsize = tsize*4
     call MEMCPYTOGPU(d_col,ja,tsize,"CG col")
     tsize = ND+1
     tsize = tsize*4
     call MEMCPYTOGPU(d_row,ia,tsize,"CG row")
     tsize = NV
     tsize = tsize*8
     call PREPARECSR(handlesp,ND,ND,NNZB,V,d_val,d_row,d_col,emptystr)

     call MEMCPYTOGPU(d_x,xx,tsize,"CG guess")
     call MEMCPYTOGPU(d_r,bb,tsize,"CG rhs")
     
     call PREPAREPRECON(amatr,ia,ja,d_val,d_row,d_col,ND,V)

     if ( NGRP /=0) then
        tsize = ND
        tsize = tsize * 4
        call MEMCPYTOGPU(d_lgrp,lgrp,tsize,"GRP index")
     end if
     !
     ! calculate deno and invb
     !
     call GPUPRECON(d_r,d_s,NV,emptystr)

     call CUDADDOT(handlebl,NDOT,d_r,1,d_s,1,deno)
     call CUDADDOT(handlebl,NDOT,d_r,1,d_r,1,invb)
  end if
  call PAR_SUM(deno,'IN MY CODE')
  call PAR_SUM(invb,'IN MY CODE')
  deno = sqrt(deno)
  invb = 1.0_rp /sqrt(invb)
  
  if( INOTMASTER ) then          
     !
     ! Matrix multiply
     !
     call GPUDMV(handlesp,ND,ND,NNZB,V,alp,&
          d_val,d_row,d_col,d_x,bet,d_s,strcopy)                                     ! d_s = A d_x

     call CUDADCOPY(handlebl,NV,d_r,1,d_p,1)                                   ! Copy d_p = d_r (d_r is b)
     call CUDADAXPY(handlebl,NV,nalp,d_s,1,d_p,1)                              ! d_p = d_p - d_s
  end if
      
!!!----------do here initial deflation new residual vector----
  if(NGRP /= 0 ) then
     if( INOTMASTER ) then          
        tsize = NV
        tsize = tsize*8
        call MEMCPYFROMGPU(xx,d_p,tsize,"small vec")
     end if
     murhs =0.0_rp
     call wtvect(ND,NGRP,V,murhs,xx)

     if( INOTMASTER ) then
        call direct_solver_solution(solve_sol(1) % direct_solver_Deflation,murhs,mu)
     end if
  end if
  
  if ( INOTMASTER ) then
     if (NGRP /= 0 ) then
        tsize = NGRP*V
        tsize = tsize*8
        call MEMCPYTOGPU(d_grplhs,mu,tsize,"small vec")
        call CUDAWVEC(d_p,d_grplhs,d_lgrp,ND,V,emptystr)
        call CUDADAXPY(handlebl,NV,alp,d_p,1,d_x,1)
        call GPUSYNC()
     end if
     
     call GPUDMV(handlesp,ND,ND,NNZB,V,alp,&
          d_val,d_row,d_col,d_x,bet,d_s,strcopy)
     
     call CUDADAXPY(handlebl,NV,nalp,d_s,1,d_r,1)
     
!!!----------end deflation----
     
     call GPUPRECON(d_r,d_p,NV,emptystr)
     call CUDADCOPY(handlebl,NV,d_p,1,d_u,1)
!!!----preconditioning deflation-------------------------
  end if
  if (NGRP /=0 ) then
     if( INOTMASTER ) then
        call GPUDMV(handlesp,ND,ND,NNZB,V,alp,&
             d_val,d_row,d_col,d_p,bet,d_s,strcopy)
        tsize = NV
        tsize = tsize*8
        call MEMCPYFROMGPU(xx,d_s,tsize,"small vec")
     end if
     murhs = 0.0_rp
     call wtvect(ND,NGRP,V,murhs,xx)

     if( INOTMASTER ) then
        call direct_solver_solution(solve_sol(1) % direct_solver_Deflation,murhs,mu)
     end if
  end if
  
  if ( INOTMASTER ) then
     if (NGRP /= 0 ) then
        tsize = NGRP*V
        tsize = tsize*8
        call MEMCPYTOGPU(d_grplhs,mu,tsize,"small vec")
        
        call CUDAWVEC(d_s,d_grplhs,d_lgrp,ND,V,emptystr)
        
        call CUDADAXPY(handlebl,NV,nalp,d_s,1,d_p,1)
     end if
!!!----end deflation------------------------------------- 
     
     call CUDADDOT(handlebl,NDOT,d_r,1,d_u,1,gamma)
     call CUDADDOT(handlebl,NDOT,d_r,1,d_r,1,rnorm)
  end if
  call PAR_SUM(rnorm,'IN MY CODE')  
  call PAR_SUM(gamma,'IN MY CODE')
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

  do while (r1 > tolsta .and. solve_sol(1) % iters <= maxiter)   
     if( INOTMASTER ) then
        call GPUDMV(handlesp,ND,ND,NNZB,V,alp,&
             d_val,d_row,d_col,d_p,bet,d_s,strcopy)

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
        call GPUPRECON(d_r,d_u,NV,emptystr)

        call CUDADDOT(handlebl,NDOT,d_r,1,d_u,1,gamma)
     end if
     call PAR_SUM(gamma,'IN MY CODE')     
     b = gamma/gammaold
!!!-----------------------Start deflation----------------------------
     
     if ( NGRP /=0 ) then
        if ( INOTMASTER ) then
           call GPUDMV(handlesp,ND,ND,NNZB,V,alp,&
                d_val,d_row,d_col,d_u,bet,d_s,strcopy)

           tsize = NV
           tsize = tsize*8
           call MEMCPYFROMGPU(xx,d_s,tsize,"small vec")
        end if
        murhs =0.0_rp
        call wtvect(ND,NGRP,V,murhs,xx)

        if ( INOTMASTER ) then
           call direct_solver_solution(solve_sol(1) % direct_solver_Deflation,murhs,mu) 
        end if
     end if
     
     if ( INOTMASTER ) then
        if (NGRP /= 0 ) then
           tsize = NGRP*V
           tsize = tsize*8
           call MEMCPYTOGPU(d_grplhs,mu,tsize,"small vec")
           call CUDAWVEC(d_s,d_grplhs,d_lgrp,ND,V,emptystr)
!!!-----------------------End deflation------------------------------ 
        end if
        call CUDADSCAL(handlebl,NV,b,d_p,1)
        call CUDADAXPY(handlebl,NV,alp,d_u,1,d_p,1)
        if (NGRP /=0 ) then
           call CUDADAXPY(handlebl,NV,nalp,d_s,1,d_p,1)     
           call GPUSYNC()
        end if
     end if
     resi2 = resi1
     r1 = sqrt(gamma)
     resi1 = r1 / deno

     solve_sol(1) % iters = solve_sol(1) % iters + 1
  enddo
  !------SUbmitting outputs--------------
  resfi = r1 / deno
  solve_sol(1) % iters = solve_sol(1) % iters -1 

  tsize = NV
  tsize = tsize*8

  if(INOTMASTER ) call MEMCPYFROMGPU(xx,d_x,tsize,"CG unk")



  if(NGRP /= 0) then
     deallocate(mu,murhs,grpmat)
  end if
  call destroyprecon()
  if(INOTMASTER) then
     call CPUMEMUNREGISTER(amatr)
     call CPUMEMUNREGISTER(ia)
     call CPUMEMUNREGISTER(ja)
     call CPUMEMUNREGISTER(bb)
     call CPUMEMUNREGISTER(xx)
     call CPUMEMUNREGISTER(lgrp)
  end if
  call gpumemreset()

end subroutine GPUDEFLCG
