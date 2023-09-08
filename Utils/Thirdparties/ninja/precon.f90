subroutine prepareprecon(amatr,ia,ja,d_val,d_ia,d_ja,ND,V)
  use mod_communications, only : PAR_INTERFACE_NODE_EXCHANGE,PAR_INTERFACE_MATRIX_EXCHANGE
  use mod_direct_solver,  only : direct_solver_factorization
  use mod_direct_solver,  only : direct_solver_partialcleaning
  use mod_direct_solver,  only : direct_allocate_temporary_matrix
  use mod_parall,         only : PAR_COMM_MY_CODE_ARRAY
  use mod_matrix
  use def_solver
  use vargpuprecon
  use gpumanager
  use def_master,         only : kfl_paral

  implicit none
  integer*4          :: ND,V,N,ii
  integer*4          :: ia(*),ja(*)
  real*8             :: amatr(*)
  real*8,allocatable :: dia(:)
  integer*8          :: sz,d_val,d_ia,d_ja,str=0
  integer*8          :: d_tempmat
  idprecon = solve_sol(1) % kfl_preco

  if (idprecon == SOL_DIAGONAL ) then


     allocate(dia(ND*V))   

     call GENDIAPRE(amatr,ia,ja,dia,ND,V)
     call PAR_INTERFACE_NODE_EXCHANGE(V,dia,'SUM','IN MY CODE','SYNCHRONOUS')       
     call GENDIAPREUPDATE(amatr,ia,ja,dia,ND,V)
     sz = ND * V
     sz = sz * 8
     call GPUMEMASSIGN(d_dia,sz)        
     call MEMCPYTOGPU(d_dia,dia,sz,"dia prec")
     deallocate(dia)
  else if( idprecon == SOL_RAS ) then     
     !
     ! RAS: Restricted Additive Schwarz
     ! 
     NNZB = ia(ND+1) -1
     N = NNZB * V * V
     nrows = ND
     nvar = V
     permuteflag = 0
     
     call PAR_INTERFACE_MATRIX_EXCHANGE(V,amatr,PAR_COMM_MY_CODE_ARRAY(1))

     if (rasprepflag == 0 ) then
        !call cusolversphandle(handsol)
        call cusparsehandle(handsp)
        !call cusolverspattachstream(handsol,str)
        call CUSPARSE_MAT_DESCR_CREATE(matdescbaseone)
        call CUSPARSE_SET_MAT_TYPE(matdescbaseone,0)
        call CUSPARSE_SET_MAT_BASE(matdescbaseone,1)
        rasprepflag = rasprepflag + 1
     end if

     allocate(rasia((ND+1)*V))
     allocate(rasja(NNZB*V*V))
        
     if (V > 1) then

        sz = NNZB*V*V
        sz = sz*4
        call GPUMEMASSIGN(d_rasja,sz)

        sz = (ND+1)*V
        sz = sz*4
        call GPUMEMASSIGN(d_rasia,sz)
        
        sz = N
        sz = sz*8
        call GPUMEMASSIGN(d_rasmat,sz)  
        call GPUMEMASSIGN(d_tempmat,sz)
        call MEMCPYTOGPU(d_rasmat,amatr,sz,"copy ras mat v greater 1")
        
        call CUDADBSR2CSR(handsp,0,ND,ND,matdescbaseone,d_rasmat,d_ia,&
             d_ja,V,matdescbaseone,d_tempmat,d_rasia,d_rasja)

        sz = (ND+1)*V
        sz = sz*4
        call MEMCPYFROMGPU(rasia,d_rasia,sz,"copy ia v greater 1")
        sz = NNZB*V*V
        sz = sz * 4
        call MEMCPYFROMGPU(rasja,d_rasja,sz,"copy ia v greater 1")

        sz = N
        sz = N * 8
        call MEMCPYFROMGPU(amatr,d_tempmat,sz,"copy mat back v greater 1")
        call GPUMEMUNASSIGN(d_tempmat,sz)
        call GPUMEMUNASSIGN(d_rasmat,sz)

        sz = (ND+1)*V
        sz = sz*4
        call GPUMEMUNASSIGN(d_rasia,sz)
        sz = NNZB*V*V
        sz = sz * 4
        call GPUMEMUNASSIGN(d_rasja,sz)

        NNZB = N
        N = NNZB 
        nrows = ND * V
        nvar = 1
     else
        rasia = ia(1:(ND+1))
        rasja = ja(1:NNZB)
     end if
             
     if ( permuteflag == 1 ) then
        allocate(permvec(nrows))
        allocate(map(NNZB))
        do ii = 1,NNZB
           map(ii)= ii-1
        end do
        !call cusolverspxcsrsymrcm(handsol,nrows,NNZB,matdescbaseone,rasia,rasja,permvec)
        !call cusolverspxcsrperm_buffersize(handsol,nrows,nrows,NNZB,matdescbaseone,rasia,rasja,permvec,permvec,permbuffsz)
        permbuffsz = permbuffsz / 4
        allocate(permbuff(permbuffsz))
        !call cusolverspxcsrperm(handsol,nrows,nrows,NNZB,matdescbaseone,rasia,rasja,permvec,permvec,map,permbuff)
     end if
     
     if ( permuteflag == 1 ) then
        sz = NNZB
        sz = sz * 4
        call GPUMEMASSIGN(d_map,sz)
        call MEMCPYTOGPU(d_map,map,sz,"copy permute map")
        
        sz = nrows
        sz = sz * 4
        call GPUMEMASSIGN(d_permvec,sz)
        call MEMCPYTOGPU(d_permvec,permvec,sz,"copy permute vector")
     end if

     sz = NNZB
     sz = sz*4
     call GPUMEMASSIGN(d_rasja,sz)
     call MEMCPYTOGPU(d_rasja,rasja,sz,"RAS IA")
     sz = nrows+1
     sz = sz*4
     call GPUMEMASSIGN(d_rasia,sz)
     call MEMCPYTOGPU(d_rasia,rasia,sz,"RAS JA")
     
     
     sz = N
     sz = sz*8
     call GPUMEMASSIGN(d_rasmat,sz)
     
     if ( permuteflag == 1 ) then
        call GPUMEMASSIGN(d_tempmat,sz)
        call MEMCPYTOGPU(d_tempmat,amatr,sz,"copy mat RAS")     
        call APPLYPERMUTATION(d_rasmat,d_tempmat,d_map,NNZB)
        call GPUMEMUNASSIGN(d_tempmat,sz)
     else
        call MEMCPYTOGPU(d_rasmat,amatr,sz,"copy mat RAS")     
     end if
     
     !call cusolverspcreatecsrqrinfo(solinfo)
     !call cusolverspxcsrqranalysis(handsol,nrows,nrows,NNZB,matdescbaseone,d_rasia,d_rasja,solinfo)
     !call cusolverspdcsrqrbufferinfo(handsol,nrows,nrows,NNZB,matdescbaseone,d_rasmat,d_rasia,d_rasja,solinfo,datain,workin)
     !call GPUMEMASSIGN(qrbuff,workin)        
     !call cusolverspdcsrqrsetup(handsol,nrows,nrows,NNZB,matdescbaseone,d_rasmat,d_rasia,d_rasja,0,solinfo)
     !call cusolverspdcsrqrfactor(handsol,nrows,nrows,NNZB,solinfo,qrbuff)

     deallocate(rasia,rasja)
     if ( permuteflag == 1 ) then
        deallocate(permvec)
        deallocate(map)
        deallocate(permbuff)
     end if
  else
     
     call runend("No preconditioner found for GPU solvers")
     
  end if
end subroutine prepareprecon
   
subroutine gpuprecon(in,out,NV,str)
  use vargpuprecon
  use def_solver
  use gpumanager, only : gpumemassign,gpumemunassign
  implicit none
  integer*8  :: in,out,str
  integer*4  :: NV  
  integer*8  :: d_tempmat
  integer*8  :: sz
  
  if (idprecon == SOL_DIAGONAL ) then
     
     call DMULTIPLYBYELEM(d_dia,in,out,NV,str)
     
  else if (idprecon == SOL_RAS) then

     sz = NV
     sz = NV*8
     if ( permuteflag == 1 ) then
        call GPUMEMASSIGN(d_tempmat,sz)
        call APPLYPERMUTATION(d_tempmat,in,d_permvec,NV)
        !call cusolverspdcsrqrsolve(handsol,nrows,nrows,d_tempmat,in,solinfo,qrbuff)
        call APPLYREVPERMUTATION(out,in,d_permvec,NV)
        call GPUMEMUNASSIGN(d_tempmat,sz)
     else
        !call cusolverspdcsrqrsolve(handsol,nrows,nrows,in,out,solinfo,qrbuff)
     end if
     
  else
     
     call runend("No preconditioner found for GPU solvers")
  
  end if
  
end subroutine gpuprecon


subroutine gpupreconinv(in,out,NV,str)
  use vargpuprecon
  use def_solver
  implicit none
  
  integer*8  :: in,out,str
  integer*4  :: NV
  real*8               :: alp,bet
  
  if (idprecon == SOL_DIAGONAL ) then

     call DDIVBYELEM(d_dia,in,out,NV,str)

  else if (idprecon == SOL_RAS ) then
     
     bet = 0.0
     alp = 1.0
     
     call cusparseattachstream(handsp,str)
     call CUDADBSRMV(handsp,0,0,nrows,nrows,NNZB,alp,matdescbaseone,&
          d_rasmat,d_rasia,d_rasja,nvar,in,bet,out)
     
  else
     
     call runend("No preconditioner found for GPU solvers")
     
  end if
  
end subroutine gpupreconinv

subroutine destroyprecon()
  use vargpuprecon
  use def_solver
  implicit none

  if (idprecon == SOL_RAS ) then
     !call cusolverspdestroycsrqrinfo(solinfo)
  end if
  
end subroutine destroyprecon

subroutine permute_csr(ia,ja,rasia,rasja,perm,map,ND,NNZB)
  !!ALSO gets converted to base zero along with permuting
  implicit none
  integer*4 :: ND,NNZB
  integer*4 :: ia(ND+1),ja(NNZB),rasia(ND+1),rasja(NNZB)
  integer*4 :: map(NNZB),perm(ND)
  integer*4 :: ii,jj
  
  rasia(1) = 0
  do ii = 1,ND
     rasia(ii+1) = rasia(ii) + ia(perm(ii)+1+1) - ia(perm(ii)+1)
  end do
  
  do ii = 1,NNZB
     rasja(ii) = ja(map(ii)+1)
  end do

  do ii = 1,NNZB
     rasja(ii) = perm(rasja(ii))
  end do
end subroutine permute_csr
