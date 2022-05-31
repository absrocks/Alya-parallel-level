subroutine tem_memall()
  !-----------------------------------------------------------------------
  !****f* temper/tem_memall
  ! NAME 
  !    tem_memall
  ! DESCRIPTION
  !    This routine allocates memory for the arrays needed to solve the
  !    temperature equation
  ! USES
  ! USED BY
  !    tem_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_kermod,     only : kfl_adj_prob, kfl_ndvars_opt
  use def_domain
  use def_solver
  use def_temper
  use mod_memchk
  use mod_memory
  use mod_tem_arrays, only : tem_arrays
  implicit none
  integer(ip) :: ielem,pelty,pnode
  !
  ! Primary variables
  !
  call tem_arrays('ALLOCATE')
  !
  ! variables needed for adjoint
  !     
  if (kfl_adj_prob == 1) then
     !
     ! TEMPER_FORWARD: Temperature known
     !
     call memory_alloca(mem_modul(1:2,modul),'TEMPE_FORW','tem_memall',tempe_forw, nunkn_tem,ncomp_tem)
     if(kfl_ellen_tem==-1) then
        call memory_alloca(mem_modul(1:2,modul),'GRTEM_TEM','tem_memall',grtem_tem,ndime,nunkn_tem)
     end if
     !
     ! resdiff_tem: partial derivatives of R w.r.t design variables 
     !
     call memory_alloca(mem_modul(1:2,modul),'RESDIFF_TEM','tem_memall',resdiff_tem, kfl_ndvars_opt,solve_sol(1) % nzrhs)
     if(postp(1) % npp_stepi(8,0)>0) then
        call memory_alloca(mem_modul(1:2,modul),'TEOLD_TEM','tem_memall',teold_tem,nunkn_tem) 
     end if
     !
     ! Coupling with nastin for coupled adjoint solution
     !
     if(kfl_coupl(ID_TEMPER,ID_NASTIN) == 1 ) then
        call memory_alloca(mem_modul(1:2,modul),'RhsadjNas_tem','tem_memall',RhsadjNas_tem, nelem)
        do ielem=1,nelem
           pelty = ltype(ielem)
           pnode = nnode(pelty)
           call memory_alloca(mem_modul(1:2,modul),'RhsadjNas_tem','tem_memall',RhsadjNas_tem(ielem)%a, ndime,pnode)
        end do
     end if
     !
     ! Coupling with chemic for coupled adjoint solution
     !
     if(kfl_coupl(ID_TEMPER,ID_CHEMIC) == 1 ) then
        call memory_alloca(mem_modul(1:2,modul),'RhsadjChe_tem','tem_memall', RhsadjChe_tem, nelem)
        do ielem=1,nelem
           pelty = ltype(ielem)
           pnode = nnode(pelty)
           call memory_alloca(mem_modul(1:2,modul),'RhsadjChe_tem','tem_memall',RhsadjChe_tem(ielem)%a, pnode,nspec)
        end do
     endif

  endif

end subroutine tem_memall
      
