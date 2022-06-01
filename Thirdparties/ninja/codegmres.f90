subroutine GPUGMRES(ND,V,maxiter,restart,amatr,ia,ja,bb,xx)
  use gpuvargmres
  use gpumanager
  use def_kintyp, only : ip,rp
  use def_solver, only : resin,resfi,resi1,resi2,iters,solve_sol
  implicit none
  integer*4 :: ND,V,maxiter,restart,convergence=0,inneriter
  integer*4 :: ia(*),ja(*)
  real*8 :: amatr(*),bb(*),xx(*)
  integer*8 :: tsize,emptystr
  integer*4 :: NNZB,N,NV,exlen,off,NVOWN,NDOT
  real*8 :: bnrm2,alp,bet,nalp,error,temp,rnorm,eps
  integer*4 :: res=1,j=1
  integer*4 :: allocerror
  
  Interface
     subroutine uppertrisolve(H,y,s,N)
       use def_kintyp, only : ip,rp
       implicit none
       integer*4 :: N
       real*8 :: H(:,:),y(:),s(:)
     end subroutine uppertrisolve
  end Interface

  convergence = 0
  emptystr = 0

  if(INOTMASTER) then
     exlen = (ND-npoi1)*V
     off = npoi1 * V
     alp = 1.0_rp
     bet = 0.0_rp
     nalp = -1.0_rp
     NNZB = ia(ND+1) -1
     N = NNZB * V * V
     NV = ND * V
     if (npoi3 .ne. -1) then
        NDOT  = npoi3 * V
     else
        NDOT = ND * V
     end if
  end if
  
  call gpumeminitiate(ND,V)
  allocate(cs(restart))
  cs = 0.0_rp
  allocate(s(restart))
  s = 0.0_rp
  allocate(y(restart))
  y = 0.0_rp
  allocate(sn(restart))
  sn = 0.0_rp
  allocate(H(restart+1,restart))
  H = 0
  
!!!!!!!! INITIALIZE -----------------
  if(INOTMASTER) then
     if (tiii == 1) then
        call GPUSTREAMCREATE(strcopy)
        call GPUSTREAMCREATE(strcomp)
        call CUBLASHANDLE(handlebl)
        call CUSPARSEHANDLE(handlesp)
        call CUBLASATTACHSTREAM(handlebl,strcopy)
        call CUSPARSEATTACHSTREAM(handlesp,strcopy)
        tiii = tiii + 1
     end if

     tsize = N
     tsize = tsize*8
     call GPUMEMASSIGN(d_val,tsize)
     call CPUMEMREGISTER(amatr,tsize)

     tsize = NNZB
     tsize = tsize*4
     call GPUMEMASSIGN(d_col,tsize)
     call CPUMEMREGISTER(ja,tsize)

     tsize = ND + 1
     tsize = tsize*4
     call GPUMEMASSIGN(d_row,tsize)
     call CPUMEMREGISTER(ia,tsize)

     tsize = NV
     tsize = tsize*8
     call GPUMEMASSIGN(d_x,tsize)
     call GPUMEMASSIGN(d_b,tsize)
     call CPUMEMREGISTER(xx,tsize)
     call CPUMEMREGISTER(bb,tsize)           

     call GPUMEMASSIGN(d_r,tsize)
     call GPUMEMASSIGN(d_w,tsize)
     call GPUMEMASSIGN(d_Ax,tsize)
     tsize = NV
     tsize = tsize*(restart+1)
     tsize = tsize*8
     call GPUMEMASSIGN(d_Krylov,tsize)
     tsize = restart
     tsize = tsize*8
     call GPUMEMASSIGN(d_y, tsize)
     
     tsize = N
     tsize = tsize*8
     call MEMCPYTOGPUASYNC(d_val,amatr,tsize,strcopy,"GMRES matrix")
     
     tsize = NNZB
     tsize = tsize*4
     call MEMCPYTOGPUASYNC(d_col,ja,tsize,strcomp,"GMRES col")
     tsize = ND+1
     tsize = tsize*4
     call MEMCPYTOGPUASYNC(d_row,ia,tsize,strcomp,"GMRES row")     
     tsize = NV
     tsize = tsize*8
     call MEMCPYTOGPUASYNC(d_x,xx,tsize,strcomp,"GMRES unk")
     call MEMCPYTOGPUASYNC(d_b,bb,tsize,strcomp,"GMRES rhs")

     call PREPARECSR(handlesp,ND,ND,NNZB,V,d_val,d_row,d_col,emptystr)

!!!!!!!!!!!!INITIAL CALCULATIONS---------------------
     
     call prepareprecon(amatr,ia,ja,d_val,d_row,d_col,ND,V)
     
     call CUDADDOT(handlebl,NDOT,d_b,1,d_b,1,bnrm2)
  end if
  call PAR_SUM(bnrm2,'IN MY CODE')
  bnrm2 = sqrt(bnrm2);

  if(bnrm2 == 0.0_rp) then
     bnrm2 =1.0_rp
  end if
  if ( INOTMASTER ) then
     call CUDADCOPY(handlebl,NV,d_b,1,d_r,1)
     call GPUDMV(handlesp,ND,ND,NNZB,V,alp,&
          d_val,d_row,d_col,d_x,bet,d_Ax,strcopy)

     call CUDADAXPY(handlebl,NV,nalp,d_Ax,1,d_r,1)
     call CUDADDOT(handlebl,NDOT,d_r,1,d_r,1,rnorm)
  end if
  call PAR_SUM(rnorm,'IN MY CODE')

  rnorm = sqrt(rnorm)
  s(1) = rnorm
  temp = 1.0_rp/rnorm
  if ( INOTMASTER ) then
     
     call CUDADSCAL(handlebl,NV,temp,d_r,1)
     call CUDADCOPY(handlebl,NV,d_r,1,d_Krylov,1)
     call GPUPRECONINV(d_x,d_x,NV,strcopy)
  end if
  
  !-------------------ADAPTIVE RESIDUAL--------------------

  solve_sol(1) % resi2 = rnorm  / bnrm2     !normalized residual
  solve_sol(1) % resin = rnorm  / bnrm2     !preconditioned residual

  if( solve_sol(1) % kfl_adres == 1 ) then
     eps = max( solve_sol(1) % resin * solve_sol(1) % adres , solve_sol(1) % solmi )
  else
     eps = solve_sol(1) % solco
  end if
  error = rnorm 
  tolsta = eps*bnrm2
  solve_sol(1) % iters = 1
  
  inneriter = solve_sol(1) % iters
!!!!!---------------Starting iterations-------------------------------
  do while(solve_sol(1) % iters<=maxiter .and. error > tolsta)

     !!-------------KRylov iterations--------------------
     do res=1,restart
        if ( INOTMASTER ) then
           call GPUPRECON(d_Krylov + (res-1)*NV*8,d_Ax,NV,strcopy)
     
           call GPUDMV(handlesp,ND,ND,NNZB,V,alp,&
                d_val,d_row,d_col,d_Ax,bet,d_w,strcopy)
           

        endif
!!!-------Global Re-orthogonalization----------------------------------
        
        do j=1,res
           if ( INOTMASTER ) then 
              call CUDADDOT(handlebl,NDOT,d_w,1,d_Krylov + (j-1)*NV*8 ,1,temp)
           end if

           call PAR_SUM(temp,'IN MY CODE')
           H(j,res) = temp
           if ( INOTMASTER ) then
              temp = -1*temp
              call CUDADAXPY(handlebl,NV,temp,d_Krylov + (j-1)*NV*8,1,d_w,1)
           end if
        end do

        if ( INOTMASTER ) then
           call CUDADDOT(handlebl,NDOT,d_w,1,d_w,1,temp)
        end if
        call PAR_SUM(temp,'IN MY CODE')
        temp = sqrt(temp)           
        H(res+1,res) = temp
        temp = 1.0_rp/temp     
        if ( INOTMASTER ) then
           call CUDADSCAL(handlebl,NV,temp,d_w,1)
           call CUDADCOPY(handlebl,NV,d_w,1,d_Krylov + res*NV*8,1)
        end if
        !!------Givens Rotations------------------------------------
        do j=1,res-1
           temp = cs(j)*H(j,res) + sn(j)*H(j+1,res)
           H(j+1,res) = -sn(j)*H(j,res) + cs(j)*H(j+1,res)
           H(j,res) = temp
        end do
        call GIVENSROT2(cs(res),sn(res),H(res,res),H(res+1,res))
        s(res+1) = -sn(res)*s(res)
        s(res) = cs(res) * s(res)
        H(res,res) = cs(res)*H(res,res) + sn(res)*H(res+1,res)

        !----test convergence-----
        error = abs(s(res+1))
        resi2 = resi1
        resi1 = error / bnrm2

        if(error <= tolsta) then
           convergence = 1
           exit
        end if        
     end do
     
     !!---------upper triangle solve and linear combination------
     if (res == (restart+1) ) res = res - 1

     if ( INOTMASTER ) then
        call UPPERTRISOLVE(H,y,s,res)
        tsize = res*8
        call MEMCPYTOGPU(d_y,y,tsize,"rhs UT sys")
        call LINEARCOMBOGPU(d_x,d_Krylov,d_y,NV,res,strcopy)
!!!----update new vector and residual-------------------

        call GPUPRECON(d_x,d_w,NV,strcopy)
        
        call GPUDMV(handlesp,ND,ND,NNZB,V,alp,&
             d_val,d_row,d_col,d_w,bet,d_Ax,strcopy)
                        
        call CUDADCOPY(handlebl,NV,d_b,1,d_r,1)
        call CUDADAXPY(handlebl,NV,nalp,d_Ax,1,d_r,1)
        
        call CUDADDOT(handlebl,NDOT,d_r,1,d_r,1,temp)     
     end if
     call PAR_SUM(temp,'IN MY CODE')
     temp = sqrt(temp)
     s(1) = temp
     resi2 = resi1
     resi1 = temp / bnrm2
     error = temp
     if ( INOTMASTER ) then
        temp = 1.0_rp/temp
        call CUDADSCAL(handlebl,NV,temp,d_r,1)
        call CUDADCOPY(handlebl,NV,d_r,1,d_Krylov,1)
     end if

     if ( convergence == 1) then
        exit
     end if
     inneriter = inneriter + 1
     solve_sol(1) % iters = (inneriter -1)*restart + res 
  end do

  tsize = NV
  tsize = tsize*8  
  if ( INOTMASTER) then
     call MEMCPYFROMGPU(xx,d_w,tsize,"GMRES unk")
  end if
  !------SUbmitting outputs--------------
  resfi = error / bnrm2
  solve_sol(1) % iters = (inneriter -1)*restart + res 
  call destroyprecon()
  if (INOTMASTER) then
     call CPUMEMUNREGISTER(amatr)
     call CPUMEMUNREGISTER(ia)
     call CPUMEMUNREGISTER(ja)
     call CPUMEMUNREGISTER(bb)
     call CPUMEMUNREGISTER(xx)
  end if
  call gpumemreset()

  deallocate(cs)
  deallocate(s)
  deallocate(y)
  deallocate(sn)
  deallocate(H)
  
end subroutine GPUGMRES

subroutine uppertrisolve(H,y,s,N)
  use def_kintyp, only : ip,rp
  implicit none
  integer*4 :: N,i,j
  real*8 :: H(:,:),y(:),s(:)
  y(N) = s(N) / H(N,N)
  do i=N-1,1,-1
     do j=N,i+1,-1
        s(i) = s(i) - H(i,j)*y(j)
     end do
     y(i) = s(i) / H(i,i)
  end do
end subroutine uppertrisolve

subroutine givensrot2(c,s,a,b)
  use def_kintyp, only : ip,rp
 implicit none
  real*8 :: c,s,a,b
  real*8 :: temp
  temp = a*a + b*b
  temp = sqrt(temp)
  if (temp == 0.0_rp) then
     temp =1.0e-16_rp
  end if
  temp = 1.0_rp / temp
  c = a * temp
  s = b * temp
end subroutine givensrot2
