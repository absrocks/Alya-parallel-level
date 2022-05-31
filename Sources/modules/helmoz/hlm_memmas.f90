subroutine hlm_memmas()

  !-----------------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_memmas.f90
  ! NAME 
  !    hlm_memmas
  ! DESCRIPTION
  !    This routine allocates memory for problem unknowns for master.
  ! USES
  ! USED BY
  !    hlm_output
  !-----------------------------------------------------------------------

  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_solver
  use def_helmoz
  !use mod_memchk
  use mod_memory

  implicit none

  integer(4)  :: istat

  if ( IMASTER ) then

     !Problem unknowns 
     !allocate(smgvp_hlm(ndime,nsite_hlm*nmlsi_hlm),stat=istat)       !Secondary magnetic vector potential
     !call memchk(zero,istat,mem_modul(1:2,modul),'SMGVP_HLM','hlm_memall',smgvp_hlm)
     !allocate(selsp_hlm(nsite_hlm*nmlsi_hlm),      stat=istat)       !Secondary electric scalar potential
     !call memchk(zero,istat,mem_modul(1:2,modul),'SELSP_HLM','hlm_memall',selsp_hlm)
     !
     !allocate(elefi_hlm(ndime,nsite_hlm),stat=istat)       !Vector of electric field intensity
     !call memchk(zero,istat,mem_modul(1:2,modul),'ELEFI_HLM','hlm_memall',elefi_hlm)
     !allocate(magfi_hlm(ndime,nsite_hlm),stat=istat)       !Vector of magnetic field intensity
     !call memchk(zero,istat,mem_modul(1:2,modul),'MAGFI_HLM','hlm_memall',magfi_hlm)    

     call memory_alloca(mem_modul(1:2,modul),'SMGVP_HLM','hlm_memmas',smgvp_hlm,ndime,nsite_hlm*nmlsi_hlm)
     call memory_alloca(mem_modul(1:2,modul),'SELSP_HLM','hlm_memmas',selsp_hlm,nsite_hlm*nmlsi_hlm)
     call memory_alloca(mem_modul(1:2,modul),'ELEFI_HLM','hlm_memmas',elefi_hlm,ndime,nsite_hlm)
     call memory_alloca(mem_modul(1:2,modul),'MAGFI_HLM','hlm_memmas',magfi_hlm,ndime,nsite_hlm)



  end if

end subroutine hlm_memmas
