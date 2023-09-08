subroutine sld_memphy(itask)
  !-----------------------------------------------------------------------
  !****f* Solidz/sld_memphy
  ! NAME
  !    sld_memphy
  ! DESCRIPTION
  !    Allocate memory for the physical problem
  ! OUTPUT
  ! USES
  ! USED BY
  !    sld_reaphy
  !    sld_sendat
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_solidz
  use mod_memchk
  use mod_memory
  implicit none
  integer(ip), intent(in) :: itask
  integer(4)              :: istat

  select case(itask)

  case(0_ip)

     !----------------------------------------------------------------------
     !
     ! Problem definition
     !
     !----------------------------------------------------------------------

  case(1_ip)

     !----------------------------------------------------------------------
     !
     ! Properties
     !
     !----------------------------------------------------------------------
     !
     ! Density and velocity of sound
     !
     call memory_alloca(mem_modul(1:2,modul),'DENSI_SLD','sld_memphy',densi_sld,ncoef_sld,nmate_sld)
     call memory_alloca(mem_modul(1:2,modul),'VELAS_SLD','sld_memphy',velas_sld,ncoef_sld,nmate_sld)
     !
     ! Rayleigh damping
     !
     allocate(dampi_sld(2,nmate_sld),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'DAMPI_SLD','sld_memphy',dampi_sld)
     allocate(kfl_dampi_sld(nmate_sld),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'KFL_DAMPI_SLD','sld_memphy',kfl_dampi_sld)
     kfl_dampi_sld = 0                          ! Numerical damping
     dampi_sld     = 0.0_rp                     ! Rayleigh damping parameters alpha and beta
     !
     ! Thermal effects
     !
     allocate(deltt_sld(2),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'DELTT_SLD','sld_memphy',deltt_sld)
     deltt_sld = 0.0_rp

     !
     ! Constitutive law: PARCO_SLD, LAWCO_SLD, LAWST_CO, LAWMO_SLD, PARCH_SLD
     !
     call memory_alloca(mem_modul(1:2,modul),'PARCO_SLD','sld_memphy',parco_sld,ncoef_sld,nmate_sld)
     allocate(parcc_sld(ncoef_sld,nmate_sld),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'PARCC_SLD','sld_memphy',parcc_sld)
     call memory_alloca(mem_modul(1:2,modul),'LAWCO_SLD','sld_memphy',lawco_sld,nmate_sld)
     allocate(lawst_sld(nmate_sld),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LAWST_SLD','sld_memphy',lawst_sld)
     allocate(lawch_sld(nmate_sld),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LAWCH_SLD','sld_memphy',lawch_sld)
     allocate(lawpl_sld(nmate_sld),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LAWPL_SLD','sld_memphy',lawpl_sld)
     allocate(lawmo_sld(nmate_sld),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LAWMO_SLD','sld_memphy',lawmo_sld)
     call memory_alloca(mem_modul(1:2,modul),'LAWTA_SLD','sld_memphy',lawta_sld,nmate_sld)
     allocate(lawho_sld(nmate_sld),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LAWHO_SLD','sld_memphy',lawho_sld)
     allocate(modfi_sld(nmate_sld),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'MODFI_SLD','sld_memphy',modfi_sld)
     modfi_sld=-1
     call memory_alloca(mem_modul(1:2,modul),'MODOR_SLD','sld_memphy',modor_sld,2_ip,nmate_sld)
     modor_sld= 0
     allocate(cocof_sld(nmate_sld),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'COCOF_SLD','sld_memphy',cocof_sld)
     cocof_sld=1.0_rp !value for human
     !allocate(prest_sld(nmate_sld),stat=istat)
     !call memchk(zero,istat,mem_modul(1:2,modul),'PREST_SLD','sld_memphy',prest_sld)
     !prest_sld=0.0_rp
     allocate(preti_sld(nmate_sld),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'PRETI_SLD','sld_memphy',preti_sld)
     preti_sld=0.0_rp
     allocate(timec_sld(nmate_sld),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'TIMEC_SLD','sld_memphy',timec_sld)
     timec_sld=0.06_rp  !value for human - in sec
     allocate(hillc_sld(nmate_sld),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'HILLC_SLD','sld_memphy',hillc_sld)
     hillc_sld=3.0_rp  !default value
     allocate(cal50_sld(nmate_sld),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'CAL50_SLD','sld_memphy',cal50_sld)
     cal50_sld=0.5_rp  !default value
     allocate(ortk1_sld(nmate_sld),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'ORTK1_SLD','sld_memphy',ortk1_sld)
     ortk1_sld=0.0_rp  !default value
     allocate(ortk2_sld(nmate_sld),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'ORTK2_SLD','sld_memphy',ortk2_sld)
     ortk2_sld=0.0_rp  !default value
     allocate(trans_sld(2,nmate_sld),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'TRANS_SLD','sld_memphy',trans_sld)
     trans_sld=0.0_rp
     allocate(kfl_coupt_sld(nmate_sld),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'KFL_COUPT_SLD','sld_memphy',kfl_coupt_sld)
     kfl_coupt_sld=0
     allocate(parsp_sld(ncoef_sld,nmate_sld),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'PARSP_SLD','sld_memphy',parsp_sld)
     allocate(parch_sld(ncoef_sld,nmate_sld),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'PARCH_SLD','sld_memphy',parch_sld)
     allocate(parcf_sld(ncoef_sld,nmate_sld),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'PARCF_SLD','sld_memphy',parcf_sld)

     !
     ! Constitutive law :Bridge (added by PIERRE OJO: ADD CHECK MEMORY )
     !
     allocate(gpsl0_sld(ndime,ndime,30),stat=istat)  !30 being the max possible nb of GP
     gpsl0_sld=0.0_rp
     allocate(rstr0_sld(ndime,ndime,30),stat=istat)
     rstr0_sld=0.0_rp
     !
     ! Stiffness matrix for Orthotropic models (undamaged)
     !
     call memory_alloca(mem_modul(1:2,modul),'STIFF0_SLD','sld_memphy',stiff0_sld,ndime*2,ndime*2,nmate_sld)

  case(2_ip)

     !----------------------------------------------------------------------
     !
     ! Parameters
     !
     !----------------------------------------------------------------------
     !
     ! Material axes
     !
     call memory_alloca(mem_modul(1:2,modul),'AXIS1_SLD','sld_memphy',axis1_sld,ndime,nelem)
     call memory_alloca(mem_modul(1:2,modul),'AXIS2_SLD','sld_memphy',axis2_sld,ndime,nelem)
     if (ndime == 3_ip) then
        call memory_alloca(mem_modul(1:2,modul),'AXIS3_SLD','sld_memphy',axis3_sld,ndime,nelem)
     end if

  end select

end subroutine sld_memphy
