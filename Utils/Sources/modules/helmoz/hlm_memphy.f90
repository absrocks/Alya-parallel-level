subroutine hlm_memphy(itask)

  !-----------------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_memphy.f90
  ! NAME
  !    hlm_memphy
  ! DESCRIPTION
  !    This routine allocates memory for the physical problem.
  ! USES
  ! USED BY
  !    hlm_reaphy
  !-----------------------------------------------------------------------

  use def_master
  use def_domain
  use def_helmoz
  !use mod_memchk
  use mod_memory

  implicit none

  integer(ip), intent(in) :: itask

  integer(ip)             :: kpoin,kelem,knode,isour,idime,ipoin
  integer(4)              :: istat

  select case ( itask )

  case ( 1_ip )

	!Properties
!	allocate(perma_hlm(nmate),stat=istat)       !mu   = Magnetic permeability
!	call memchk(zero,istat,mem_modul(1:2,modul),'PERMA_HLM','hlm_memphy',perma_hlm)
!	allocate(epsil_hlm(nmate),stat=istat)       !eps  = Dielectric permittivity
!	call memchk(zero,istat,mem_modul(1:2,modul),'EPSIL_HLM','hlm_memphy',epsil_hlm)
!	allocate(bckco_hlm(ncond_hlm),stat=istat)   !bckco_hlm  = Background electric conductivity tensor
!	call memchk(zero,istat,mem_modul(1:2,modul),'BCKCO_HLM','hlm_memphy',bckco_hlm) 
!	allocate(sigma_hlm(ncond_hlm,nmate),stat=istat)       !sig  = Electric conductivity tensor
!	call memchk(zero,istat,mem_modul(1:2,modul),'SIGMA_HLM','hlm_memphy',sigma_hlm) 
!	allocate(dsigma_hlm(ncond_hlm,nmate),stat=istat)      !dsig = Electric conductivity tensor - Primary electric conductivity tensor
!	call memchk(zero,istat,mem_modul(1:2,modul),'DSIGMA_HLM','hlm_memphy',dsigma_hlm) 
  
        call memory_alloca(mem_modul(1:2,modul),'PERMA_HLM','hlm_memphy',perma_hlm,nmate)
        call memory_alloca(mem_modul(1:2,modul),'EPSIL_HLM','hlm_memphy',epsil_hlm,nmate)
        call memory_alloca(mem_modul(1:2,modul),'BCKCO_HLM','hlm_memphy',bckco_hlm,ncond_hlm)
        call memory_alloca(mem_modul(1:2,modul),'SIGMA_HLM','hlm_memphy',sigma_hlm,ncond_hlm,nmate)
        call memory_alloca(mem_modul(1:2,modul),'DSIGMA_HLM','hlm_memphy',dsigma_hlm,ncond_hlm,nmate)

  case ( 2_ip )

  if (ppcod_hlm == 1) then
write(*,*) 'Allocating PVEPO_HLM, I am ', kfl_paral
	!Primary vector potential in nodes of a reference mesh, of size nz * nr, in the z-r plane 
!	allocate(z_hlm(nz_hlm),stat=istat)
!	call memchk(zero,istat,mem_modul(1:2,modul),'Z_HLM','hlm_memphy',z_hlm)
!	allocate(r_hlm(nr_hlm),stat=istat)
!	call memchk(zero,istat,mem_modul(1:2,modul),'R_HLM','hlm_memphy',r_hlm)
!	allocate(pvepo_hlm(nz_hlm,nr_hlm),stat=istat)
!	call memchk(zero,istat,mem_modul(1:2,modul),'PVEPO_HLM','hlm_memphy',pvepo_hlm)


        call memory_alloca(mem_modul(1:2,modul),'Z_HLM','hlm_memphy',z_hlm,nz_hlm)
        call memory_alloca(mem_modul(1:2,modul),'R_HLM','hlm_memphy',r_hlm,nr_hlm)
        call memory_alloca(mem_modul(1:2,modul),'PVEPO_HLM','hlm_memphy',pvepo_hlm,nz_hlm,nr_hlm)



  endif

  if (ppcod_hlm == 2) then
write(*,*) 'Allocating PMGVP_HLM in MEMPHY, I am ', kfl_paral, 'npoin = ', npoin
	!Primary potentials in the mesh nodes
!	allocate(pmgvp_hlm(ndime,npoin),stat=istat)       !Primary magnetic vector potential
!	call memchk(zero,istat,mem_modul(1:2,modul),'PMGVP_HLM','hlm_memall',pmgvp_hlm)
!	allocate(pelsp_hlm(npoin),      stat=istat)       !Primary electric scalar potential
!	call memchk(zero,istat,mem_modul(1:2,modul),'PELSP_HLM','hlm_memall',pelsp_hlm)

        call memory_alloca(mem_modul(1:2,modul),'PMGVP_HLM','hlm_memphy',pmgvp_hlm,ndime,npoin)
        call memory_alloca(mem_modul(1:2,modul),'PELSP_HLM','hlm_memphy',pelsp_hlm,npoin)


  endif

  case ( 3_ip )

	!Sites
!	allocate(site_hlm(ndime,nsite_hlm),stat=istat)
!	call memchk(zero,istat,mem_modul(1:2,modul),'SITE_HLM','hlm_memphy',site_hlm)

        call memory_alloca(mem_modul(1:2,modul),'SITE_HLM','hlm_memphy',site_hlm,ndime,nsite_hlm)
        call memory_alloca(mem_modul(1:2,modul),'CLSITE1_HLM','hlm_memphy',clsite1_hlm,nsite_hlm*1_ip)
        call memory_alloca(mem_modul(1:2,modul),'CLSITE2_HLM','hlm_memphy',clsite2_hlm,nsite_hlm*1_ip)
        call memory_alloca(mem_modul(1:2,modul),'CLSITE_HLM','hlm_memphy',clsite_hlm,ndime,nsite_hlm*1_ip)

  case ( 4_ip )

	!Closest nodes for MLSI
!	allocate(clnod1_hlm(nsite_hlm*nmlsi_hlm),stat=istat)
!	call memchk(zero,istat,mem_modul(1:2,modul),'CLNOD1_HLM','hlm_memphy',clnod1_hlm)
!	allocate(clnod2_hlm(nsite_hlm*nmlsi_hlm),stat=istat)
!	call memchk(zero,istat,mem_modul(1:2,modul),'CLNOD2_HLM','hlm_memphy',clnod2_hlm)
!	allocate(clcoor_hlm(ndime,nsite_hlm*nmlsi_hlm),stat=istat)
!	call memchk(zero,istat,mem_modul(1:2,modul),'CLCOOR_HLM','hlm_memphy',clcoor_hlm)
        call memory_alloca(mem_modul(1:2,modul),'CLNOD1_HLM','hlm_memphy',clnod1_hlm,nsite_hlm*nmlsi_hlm)
        call memory_alloca(mem_modul(1:2,modul),'CLNOD2_HLM','hlm_memphy',clnod2_hlm,nsite_hlm*nmlsi_hlm)
        call memory_alloca(mem_modul(1:2,modul),'CLCOOR_HLM','hlm_memphy',clcoor_hlm,ndime,nsite_hlm*nmlsi_hlm)

  case ( 5_ip )

	!Cartesian coordinates, current and length of the shots  
	!allocate(xoffsv_hlm(nshot_hlm),stat=istat)
	!call memchk(zero,istat,mem_modul(1:2,modul),'XOFFSV_HLM','hlm_memphy',xoffsv_hlm)
	!allocate(yoffsv_hlm(nshot_hlm),stat=istat)
	!call memchk(zero,istat,mem_modul(1:2,modul),'YOFFSV_HLM','hlm_memphy',yoffsv_hlm)
	!allocate(zoffsv_hlm(nshot_hlm),stat=istat)
	!call memchk(zero,istat,mem_modul(1:2,modul),'ZOFFSV_HLM','hlm_memphy',zoffsv_hlm)
	!allocate(elcurv_hlm(nshot_hlm),stat=istat)
	!call memchk(zero,istat,mem_modul(1:2,modul),'ELCURV_HLM','hlm_memphy',elcurv_hlm)
	!allocate(lengthv_hlm(nshot_hlm),stat=istat)
	!call memchk(zero,istat,mem_modul(1:2,modul),'LENGTHV_HLM','hlm_memphy',lengthv_hlm)
	!allocate(weightv_hlm(nshot_hlm),stat=istat)
	!call memchk(zero,istat,mem_modul(1:2,modul),'WEIGHTV_HLM','hlm_memphy',weightv_hlm)
	!allocate(costfv_hlm(nshot_hlm),stat=istat)
	!call memchk(zero,istat,mem_modul(1:2,modul),'COSTFV_HLM','hlm_memphy',costfv_hlm)

        call memory_alloca(mem_modul(1:2,modul),'XOFFSV_HLM','hlm_memphy',xoffsv_hlm,nshot_hlm)
        call memory_alloca(mem_modul(1:2,modul),'YOFFSV_HLM','hlm_memphy',yoffsv_hlm,nshot_hlm)
        call memory_alloca(mem_modul(1:2,modul),'ZOFFSV_HLM','hlm_memphy',zoffsv_hlm,nshot_hlm)
        call memory_alloca(mem_modul(1:2,modul),'ELCURV_HLM','hlm_memphy',elcurv_hlm,nshot_hlm)
        call memory_alloca(mem_modul(1:2,modul),'LENGTHV_HLM','hlm_memphy',lengthv_hlm,nshot_hlm)
        call memory_alloca(mem_modul(1:2,modul),'WEIGHTV_HLM','hlm_memphy',weightv_hlm,nshot_hlm)
        call memory_alloca(mem_modul(1:2,modul),'COSTFV_HLM','hlm_memphy',costfv_hlm,nshot_hlm)

  case ( 6_ip )

	!allocate(selsp_obs(nsite_hlm*1_ip),stat=istat)
	!call memchk(zero,istat,mem_modul(1:2,modul),'SELSP_OBS','hlm_memphy',selsp_obs)
	!allocate(smgvpX_obs(nsite_hlm*1_ip),stat=istat)
	!call memchk(zero,istat,mem_modul(1:2,modul),'SMGVPX_OBS','hlm_memphy',smgvpX_obs)
	!allocate(smgvpY_obs(nsite_hlm*1_ip),stat=istat)
	!call memchk(zero,istat,mem_modul(1:2,modul),'SMGVPY_OBS','hlm_memphy',smgvpY_obs)
	!allocate(smgvpZ_obs(nsite_hlm*1_ip),stat=istat)
	!call memchk(zero,istat,mem_modul(1:2,modul),'SMGVPZ_OBS','hlm_memphy',smgvpZ_obs)

	!allocate(countobs(nsite_hlm*1_ip),stat=istat)
	!call memchk(zero,istat,mem_modul(1:2,modul),'COUNTOBS','hlm_memphy',countobs)

	!allocate(sum_obs(nsite_hlm*1_ip),stat=istat)
	!call memchk(zero,istat,mem_modul(1:2,modul),'SUM_OBS','hlm_memphy',sum_obs)


        !call memory_alloca(mem_modul(1:2,modul),'SELSP_OBS','hlm_memphy',selsp_obs,nsite_hlm*1_ip)
        !call memory_alloca(mem_modul(1:2,modul),'SMGVPX_OBS','hlm_memphy',smgvpX_obs,nsite_hlm*1_ip)
        !call memory_alloca(mem_modul(1:2,modul),'SMGVPY_OBS','hlm_memphy',smgvpY_obs,nsite_hlm*1_ip)
        !call memory_alloca(mem_modul(1:2,modul),'SMGVPZ_OBS','hlm_memphy',smgvpZ_obs,nsite_hlm*1_ip)
        call memory_alloca(mem_modul(1:2,modul),'COUNTOBS','hlm_memphy',countobs,nsite_hlm*1_ip)
        call memory_alloca(mem_modul(1:2,modul),'SUM_OBS','hlm_memphy',sum_obs,nsite_hlm*1_ip)

        call memory_alloca(mem_modul(1:2,modul),'SELSP_OBS','hlm_memphy',selsp_obs,nshot_hlm,nsite_hlm*1_ip)
        call memory_alloca(mem_modul(1:2,modul),'SMGVPX_OBS','hlm_memphy',smgvpX_obs,nshot_hlm,nsite_hlm*1_ip)
        call memory_alloca(mem_modul(1:2,modul),'SMGVPY_OBS','hlm_memphy',smgvpY_obs,nshot_hlm,nsite_hlm*1_ip)
        call memory_alloca(mem_modul(1:2,modul),'SMGVPZ_OBS','hlm_memphy',smgvpZ_obs,nshot_hlm,nsite_hlm*1_ip)

  case( 7_ip )

        call memory_alloca(mem_modul(1:2,modul),'INCIDENCE_OBS','hlm_memphy',incidence_obs,nshot_hlm,nsite_hlm)
  
  end select

end subroutine hlm_memphy
