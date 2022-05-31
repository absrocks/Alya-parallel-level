!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    exm_solmem.f90
!> @date    14/06/2019
!> @author  Mariano Vazquez
!> @brief   Allocate memory
!> @details Allocate memory for solver and variables that should be
!>          reallocated and not redistributed in case of repartitioning
!> @}
!-----------------------------------------------------------------------

subroutine exm_solmem()

  use def_master
  use def_solver
  use def_domain
  use def_exmedi
  use mod_memory
  use mod_exm_sld_eccoupling, only: exm_sld_ecc_allocate_memmory
  implicit none
  integer(ip) :: nzrhs_exm
  !
  ! Residual size
  !
  nzrhs_exm = npoin
  !
  ! Actualize maximum sizes of RHSID
  !
  nzrhs = max(nzrhs,nzrhs_exm)
  !
  ! Solver memory for NVARI variables
  !
  solve_sol => solve(1:1)
  call soldef(4_ip)
  !
  ! Other variables
  !
  call memory_alloca(mem_modul(1:2,modul),'VDIAG_EXM','exm_solmem',vdiag_exm,npoin) 
  call memory_alloca(mem_modul(1:2,modul),'IDIMA_EXM','exm_solmem',idima_exm,npoin) 

  call memory_alloca(mem_modul(1:2,modul),'TICEL_EXM','exm_solmem',ticel_exm,npoin) 
  call memory_alloca(mem_modul(1:2,modul),'JICEL_EXM','exm_solmem',jicel_exm,npoin)

  if (kfl_appli_exm > 0) then
     !
     ! Source fields
     !
     !    Applied current field flag:
     call memory_alloca(mem_modul(1:2,modul),'LAPNO_EXM','exm_solmem',lapno_exm,npoin)

     !    Applied current field:
     call memory_alloca(mem_modul(1:2,modul),'APPFI_EXM','exm_solmem',appfi_exm,npoin)
  end if


  if(kfl_algso_exm==0 .and. kfl_genal_exm > 1) then
     call runend('NOT CODED')
     !allocate(lpexm(npoin),stat=istat)
     !call memchk(zero,istat,memdi,'LPEXM','exm_solmem',lpexm)
     !call mediso(lodof,npoin,lun_solve_exm,lpexm)
  end if

  !
  ! Conductivity, etc.
  !

  call memory_alloca(mem_modul(1:2,modul),'CEDIF_EXM','exm_solmem',cedif_exm,ndime,ndime,npoin) 
  call memory_alloca(mem_modul(1:2,modul),'GRAFI_EXM','exm_solmem',grafi_exm,ndime,ndime,npoin)
  call memory_alloca(mem_modul(1:2,modul),'KGRFI_EXM','exm_solmem',kgrfi_exm,nelem)

  ! 
  ! Electro-mechanical coupling
  !
  call exm_sld_ecc_allocate_memmory(200_ip)

end subroutine exm_solmem
