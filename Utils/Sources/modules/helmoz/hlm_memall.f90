subroutine hlm_memall()

  !-----------------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_memall.f90
  ! NAME 
  !    hlm_memall
  ! DESCRIPTION
  !    This routine allocates memory for arrays needed for the module.
  ! USES
  ! USED BY
  !    hlm_turnon
  !-----------------------------------------------------------------------

  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_solver
  use def_helmoz
  use def_kermod, only : ndivi
  use mod_memory

  implicit none

  integer(4)  :: istat

  if (INOTMASTER) then
     !Problem unknowns 
     !	allocate(smgvp_hlm(ndime,npoin),stat=istat)       !Secondary magnetic vector potential
     !	call memchk(zero,istat,mem_modul(1:2,modul),'SMGVP_HLM','hlm_memall',smgvp_hlm)
     !	allocate(selsp_hlm(npoin),      stat=istat)       !Secondary electric scalar potential
     !	call memchk(zero,istat,mem_modul(1:2,modul),'SELSP_HLM','hlm_memall',selsp_hlm)

     call memory_alloca(mem_modul(1:2,modul),'SMGVP_HLM','hlm_memall',smgvp_hlm,ndime,npoin)
     call memory_alloca(mem_modul(1:2,modul),'SELSP_HLM','hlm_memall',selsp_hlm,npoin)
     call memory_alloca(mem_modul(1:2,modul),'DIFFJ_HLM','hlm_memall',diffj_hlm,npoin)
     call memory_alloca(mem_modul(1:2,modul),'DESIGN_HLM','hlm_memall',design_hlm,npoin)

     if( kfl_edges_hlm == 1 ) then
        call memory_alloca(mem_modul(1:2,modul),'SIGN_EDGES_HLM'  ,'hlm_memall',sign_edges_hlm,  meshe(ndivi) % medge,nelem)
        call memory_alloca(mem_modul(1:2,modul),'LENGTH_EDGES_HLM','hlm_memall',length_edges_hlm,meshe(ndivi) % nedge)
     end if

     !Allocate primary potentials if they haven't been allocated
     if (emmet_hlm > 3_ip .or. ppcod_hlm /= 2) then
        !write(*,*) 'Allocating PMGVP_HLM in MEMALL, I am ', kfl_paral, 'npoin = ', npoin
        !		allocate(pmgvp_hlm(ndime,npoin),stat=istat)       !Primary magnetic vector potential
        !		call memchk(zero,istat,mem_modul(1:2,modul),'PMGVP_HLM','hlm_memall',pmgvp_hlm)
        !		allocate(pelsp_hlm(npoin),      stat=istat)       !Primary electric scalar potential
        !		call memchk(zero,istat,mem_modul(1:2,modul),'PELSP_HLM','hlm_memall',pelsp_hlm)

        call memory_alloca(mem_modul(1:2,modul),'PMGVP_HLM','hlm_memall',pmgvp_hlm,ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'PELSP_HLM','hlm_memall',pelsp_hlm,npoin)


     endif

     if (ISEQUEN) then 
        !		allocate(elefi_hlm(ndime,nsite_hlm),stat=istat)       !Vector of electric field intensity
        !		call memchk(zero,istat,mem_modul(1:2,modul),'ELEFI_HLM','hlm_memall',elefi_hlm)
        !		allocate(magfi_hlm(ndime,nsite_hlm),stat=istat)       !Vector of magnetic field intensity
        !		call memchk(zero,istat,mem_modul(1:2,modul),'MAGFI_HLM','hlm_memall',magfi_hlm)

        call memory_alloca(mem_modul(1:2,modul),'ELEFI_HLM','hlm_memall',elefi_hlm,ndime,nsite_hlm)
        call memory_alloca(mem_modul(1:2,modul),'MAGFI_HLM','hlm_memall',magfi_hlm,ndime,nsite_hlm)

     endif
     !Solver memory
     solve_sol => solve(1:)


     if(kfl_servi(ID_OPTSOL)==1) then
        solad_sol => solad(1:)
     end if
     call soldef(4_ip)
  else
     !allocate(smgvp_hlm(1,1),stat=istat)
     !allocate(selsp_hlm(1),  stat=istat)
     !allocate(elefi_hlm(1,1),stat=istat)
     !allocate(magfi_hlm(1,1),stat=istat)

     !call memory_alloca(mem_modul(1:2,modul),'SMGVP_HLM','hlm_memall',smgvp_hlm,1,1)
     !call memory_alloca(mem_modul(1:2,modul),'SELSP_HLM','hlm_memall',selsp_hlm,1)
     !call memory_alloca(mem_modul(1:2,modul),'ELEFI_HLM','hlm_memall',elefi_hlm,1,1)
     !call memory_alloca(mem_modul(1:2,modul),'MAGFI_HLM','hlm_memall',magfi_hlm,1,1)
     call memory_alloca(mem_modul(1:2,modul),'SMGVP_HLM','hlm_memall',smgvp_hlm,ndime,nsite_hlm*nmlsi_hlm)
     call memory_alloca(mem_modul(1:2,modul),'SELSP_HLM','hlm_memall',selsp_hlm,nsite_hlm*nmlsi_hlm)
     call memory_alloca(mem_modul(1:2,modul),'ELEFI_HLM','hlm_memall',elefi_hlm,ndime,nsite_hlm)
     call memory_alloca(mem_modul(1:2,modul),'MAGFI_HLM','hlm_memall',magfi_hlm,ndime,nsite_hlm)

     call memory_alloca(mem_modul(1:2,modul),'DIFFJ_HLM','hlm_memall',diffj_hlm,1)
     call memory_alloca(mem_modul(1:2,modul),'DESIGN_HLM','hlm_memall',design_hlm,1)

  endif

end subroutine hlm_memall
