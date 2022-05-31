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
  use mod_memory,             only : memory_alloca
  use mod_exm_sld_eccoupling, only : exm_sld_ecc_allocate_memmory
  ! TODO : EXM-SLD-ECC : Clean 
  use mod_exm_sld_eccoupling


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
     call memory_alloca(mem_modul(1:2,modul),'LAWCO_SLD','sld_memphy',lawco_sld,nmate_sld)
     allocate(lawst_sld(nmate_sld),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LAWST_SLD','sld_memphy',lawst_sld)
     allocate(lawpl_sld(nmate_sld),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LAWPL_SLD','sld_memphy',lawpl_sld)
     allocate(lawmo_sld(nmate_sld),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LAWMO_SLD','sld_memphy',lawmo_sld)
     allocate(lawho_sld(nmate_sld),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'LAWHO_SLD','sld_memphy',lawho_sld)
     allocate(modfi_sld(nmate_sld),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'MODFI_SLD','sld_memphy',modfi_sld)
     modfi_sld=-1
     call memory_alloca(mem_modul(1:2,modul),'MODOR_SLD','sld_memphy',modor_sld,2_ip,nmate_sld)
     modor_sld= 0
     
     allocate(parsp_sld(ncoef_sld,nmate_sld),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'PARSP_SLD','sld_memphy',parsp_sld)
     !
     ! Cohesive laws and parameters
     !
     call memory_alloca(mem_modul(1:2,modul),'LAWCH_SLD','sld_memphy',lawch_sld,nmate_sld)
     call memory_alloca(mem_modul(1:2,modul),'PARCH_SLD','sld_memphy',parch_sld,ncoef_sld,nmate_sld)
     !
     ! Contact methods and parameters
     !
     call memory_alloca(mem_modul(1:2,modul),'PARCF_SLD','sld_memphy',parcf_sld,ncoef_sld,nmate_sld)
     !
     ! Tangent calculation methods for material laws
     !
     call memory_alloca(mem_modul(1:2,modul),'LAWTA_SLD','sld_memphy',lawta_sld,nmate_sld)
     !
     ! Stiffness matrix for Orthotropic models (undamaged)
     !
     call memory_alloca(mem_modul(1:2,modul),'STIFF0_SLD','sld_memphy',stiff0_sld,ndime*2,ndime*2,nmate_sld)
     
     ! 
     ! Callocate memmory for the SOLIDZ - EXMEDI coupling
     ! 
     call exm_sld_ecc_allocate_memmory(100_ip)


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
