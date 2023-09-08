subroutine gpuamgxpar(ND,V,maxiter,amatr,ia,ja,bb,xx,partvec)
  use gpuvaramgx
  use gpumanager
  use mod_parall, only : PAR_COMM_MY_CODE_WM4
  use def_master, only : kfl_paral, lninv_loc
  use def_solver, only : solve_sol

  implicit none
  integer*4        :: ND,V,maxiter,ii
  integer*4        :: ia(*),ja(*),partvec(*)
  real*8           :: amatr(*),bb(*),xx(*)
  integer*8        :: tsize
  integer*4        :: NNZB,N,NV
  integer*4        :: NVOWN,off,exlen,iters,idx
  integer*8        :: emptystr
  real*8           :: resid,resin
  
  if(INOTMASTER) then
     exlen = (ND-npoi1) * V
     off   =  (npoi1) * V 
     NNZB  = ia(ND+1) - 1
     N     = NNZB * V * V
     NV    = ND * V
     idx   = 0
  end if
  emptystr = 0

  call gpumeminitiate(ND,V)

  
  if( INOTMASTER ) then

     if(initflag == 1) then

        call AMGX_INIT()
        initflag = initflag + 1
        allocate(ja_global(npoin))
        do ii = 1,npoin
           ja_global(ii) = lninv_loc(ja(ii))
        end do
        
     end if
     
     tsize = N
     tsize = tsize*8
     call GPUMEMASSIGN(d_val,tsize)
     tsize = NV
     tsize = tsize*8
     call GPUMEMASSIGN(d_x,tsize)
     call GPUMEMASSIGN(d_r,tsize)
     tsize = (ND+1)
     tsize = tsize*4
     call GPUMEMASSIGN(d_row,tsize)
     tsize = NNZB
     tsize = tsize*4
     call GPUMEMASSIGN(d_col,tsize)
     
     
     tsize = N
     tsize = tsize*8
     call MEMCPYTOGPU(d_val,amatr,tsize,"AMGX matrix")
     !     call AMGX_MEMORY_PIN(amatr,tsize)
     tsize = NNZB
     tsize = tsize*4
     call MEMCPYTOGPU(d_col,ja,tsize,"AMGX col")
     !     call AMGX_MEMORY_PIN(ja,tsize)
     tsize = (ND+1)
     tsize = tsize*4
     call MEMCPYTOGPU(d_row,ia,tsize,"AMGX row")
     !     call AMGX_MEMORY_PIN(ia,tsize)
     tsize=NV
     tsize = tsize*8
     call MEMCPYTOGPU(d_x,xx,tsize,"AMGX guess")
     call MEMCPYTOGPU(d_r,bb,tsize,"AMGX rhs")
     !     call AMGX_MEMORY_PIN(xx,tsize)
     !     call AMGX_MEMORY_PIN(bb,tsize)
     call tobasezero(d_row,d_col,d_row,d_col,ND,ND,NNZB,emptystr)
  end if
    
  if(INOTMASTER ) then
     
     call AMGX_CREATE_CONFIG(config,solve_sol(1) % conf_file )
     call AMGX_CREATE_RESOURCE(resource,config,PAR_COMM_MY_CODE_WM4,kfl_paral)
     
     call AMGX_CREATE_MATRIX(mathand,resource)
     call AMGX_CREATE_VECTOR(rhshand,resource)
     call AMGX_CREATE_VECTOR(unkhand,resource)
     call AMGX_CREATE_SOLVER(solhand,resource,config)

     call AMGX_UPLOAD_MATRIX(mathand,ND,NNZB,V,V,d_row,d_col,d_val,0)
     call AMGX_BIND_VECTOR(rhshand,mathand)
     call AMGX_BIND_VECTOR(unkhand,mathand)
     
     call AMGX_UPLOAD_VECTOR(rhshand,ND,V,d_r)
     call AMGX_UPLOAD_VECTOR(unkhand,ND,V,d_x)

     call AMGX_SETUP_SOLVER(solhand,mathand)
     call AMGX_SOLVE_SOLVER(solhand,rhshand,unkhand)
     
     call AMGX_DOWNLOAD_VECTOR(unkhand,d_x)

     call AMGX_SOLVER_GET_ITERATIONS_NUMBER(solhand,iters)

     iters = iters - 1
     solve_sol(1) % iters = iters
     
     call AMGX_SOLVER_GET_ITERATION_RESIDUAL(solhand,iters,idx,resid)
     call AMGX_SOLVER_GET_ITERATION_RESIDUAL(solhand,1,idx,resin)
     solve_sol(1) % resin = resin
     solve_sol(1) % resfi = resid

     
     call AMGX_DESTROY_CONFIG(config)
     call AMGX_DESTROY_MATRIX(mathand)
     call AMGX_DESTROY_VECTOR(rhshand)
     call AMGX_DESTROY_VECTOR(unkhand)
     call AMGX_DESTROY_SOLVER(solhand)
     
     call AMGX_DESTROY_RESOURCE(resource)
  end if

  if(INOTMASTER) then
     call MEMCPYFROMGPU(xx,d_x,tsize,"AMGX unk")
  end if
  
  call gpumemreset()
  
end subroutine gpuamgxpar













subroutine GPUAMGX(ND,V,maxiter,amatr,ia,ja,bb,xx)
  use gpuvaramgx
  use gpumanager
  use mod_parall, only : PAR_COMM_MY_CODE_WM4
  use def_master, only : kfl_paral
  use def_solver, only : solve_sol
 
  implicit none
  integer*4        :: ND,V,maxiter
  integer*4        :: ia(*),ja(*)
  real*8           :: amatr(*),bb(*),xx(*)
  integer*8        :: tsize
  integer*4        :: NNZB,N,NV
  integer*4        :: NVOWN,off,exlen,iters,idx
  integer*8        :: emptystr
  real*8           :: resid,resin
  
  if(INOTMASTER) then
     exlen = (ND-npoi1) * V
     off   =  (npoi1) * V 
     NNZB  = ia(ND+1) - 1
     N     = NNZB * V * V
     NV    = ND * V
     idx   = 0
  end if
  emptystr = 0
  
  call gpumeminitiate(ND,V)

  
  if( INOTMASTER ) then

     if(initflag == 1) then
        
        call AMGX_INIT()
        initflag = initflag + 1
        
     end if
     
     tsize = N
     tsize = tsize*8
     call GPUMEMASSIGN(d_val,tsize)
     tsize = NV
     tsize = tsize*8
     call GPUMEMASSIGN(d_x,tsize)
     call GPUMEMASSIGN(d_r,tsize)
     tsize = (ND+1)
     tsize = tsize*4
     call GPUMEMASSIGN(d_row,tsize)
     tsize = NNZB
     tsize = tsize*4
     call GPUMEMASSIGN(d_col,tsize)
     
     
     tsize = N
     tsize = tsize*8
     call MEMCPYTOGPU(d_val,amatr,tsize,"AMGX matrix")
     !     call AMGX_MEMORY_PIN(amatr,tsize)
     tsize = NNZB
     tsize = tsize*4
     call MEMCPYTOGPU(d_col,ja,tsize,"AMGX col")
     !     call AMGX_MEMORY_PIN(ja,tsize)
     tsize = (ND+1)
     tsize = tsize*4
     call MEMCPYTOGPU(d_row,ia,tsize,"AMGX row")
     !     call AMGX_MEMORY_PIN(ia,tsize)
     tsize=NV
     tsize = tsize*8
     call MEMCPYTOGPU(d_x,xx,tsize,"AMGX guess")
     call MEMCPYTOGPU(d_r,bb,tsize,"AMGX rhs")
     !     call AMGX_MEMORY_PIN(xx,tsize)
     !     call AMGX_MEMORY_PIN(bb,tsize)
     call tobasezero(d_row,d_col,d_row,d_col,ND,ND,NNZB,emptystr)
  end if
    
  if(INOTMASTER ) then
     
     call AMGX_CREATE_CONFIG(config,solve_sol(1) % conf_file )
     call AMGX_CREATE_RESOURCE(resource,config,PAR_COMM_MY_CODE_WM4,kfl_paral)
     
     call AMGX_CREATE_MATRIX(mathand,resource)
     call AMGX_CREATE_VECTOR(rhshand,resource)
     call AMGX_CREATE_VECTOR(unkhand,resource)
     call AMGX_CREATE_SOLVER(solhand,resource,config)

     call AMGX_UPLOAD_MATRIX(mathand,ND,NNZB,V,V,d_row,d_col,d_val,0)
     call AMGX_BIND_VECTOR(rhshand,mathand)
     call AMGX_BIND_VECTOR(unkhand,mathand)
     
     call AMGX_UPLOAD_VECTOR(rhshand,ND,V,d_r)
     call AMGX_UPLOAD_VECTOR(unkhand,ND,V,d_x)

     call AMGX_SETUP_SOLVER(solhand,mathand)
     call AMGX_SOLVE_SOLVER(solhand,rhshand,unkhand)
     
     call AMGX_DOWNLOAD_VECTOR(unkhand,d_x)

     call AMGX_SOLVER_GET_ITERATIONS_NUMBER(solhand,iters)

     iters = iters - 1
     solve_sol(1) % iters = iters
     
     call AMGX_SOLVER_GET_ITERATION_RESIDUAL(solhand,iters,idx,resid)
     call AMGX_SOLVER_GET_ITERATION_RESIDUAL(solhand,1,idx,resin)
     solve_sol(1) % resin = resin
     solve_sol(1) % resfi = resid

     
     call AMGX_DESTROY_CONFIG(config)
     call AMGX_DESTROY_MATRIX(mathand)
     call AMGX_DESTROY_VECTOR(rhshand)
     call AMGX_DESTROY_VECTOR(unkhand)
     call AMGX_DESTROY_SOLVER(solhand)
     
     call AMGX_DESTROY_RESOURCE(resource)
  end if

  if(INOTMASTER) then
     call MEMCPYFROMGPU(xx,d_x,tsize,"AMGX unk")
  end if
  
  call gpumemreset()
  
end subroutine GPUAMGX

