subroutine chm_memphy(itask)
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_memphy
  ! NAME
  !    chm_memphy
  ! DESCRIPTION
  !    Allocate memory for the physical problem
  ! OUTPUT 
  ! USES
  ! USED BY
  !    chm_reaphy
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_chemic
  use mod_memchk
  use mod_memory,       only : memory_alloca
  use mod_interp_tab,   only : tab_allocate  
  implicit none
  integer(ip), intent(in) :: itask
  integer(4)              :: istat
  integer(ip)             :: ireac,ii,nrow

  select case(itask)

  case(1_ip)
     !
     ! Allocate properties
     !
     allocate(diffu_chm(2,nclas_chm),stat=istat)
     call memchk(zero,istat,mem_modul(1:2,modul),'DIFFU_CHM','chm_memphy',diffu_chm)

     allocate(speci(nclas_chm)) !!**
     !!FER: CANNOT CHECK ALLOCATION OF THE NEWLY DEFINED TYPE, IT IS NOT IN MEMCHK
     !!!FER WHY  call memchk(zero,istat,mem_modul(1:2,modul),'SPECI_CHM','chm_memphy',speci_chm)

     nullify(Le_k)
     call memory_alloca(mem_modul(1:2,modul),'LE_K','chm_memall',Le_k,nclas_chm)
     Le_k = 1.0_rp

  case(2_ip)
     !
     ! Deallocate properties
     !
     call memchk(two,istat,mem_modul(1:2,modul),'diffu_chm','chm_memphy',diffu_chm)
     deallocate(diffu_chm,stat=istat)
     if(istat/=0) call memerr(two,'diffu_chm','chm_memphy',0_ip)

    !!!FER WHY    call memchk(two,istat,mem_modul(1:2,modul),'speci_chm','chm_memphy',speci_chm)
     deallocate(speci,stat=istat) !!**
     if(istat/=0) call memerr(two,'speci','chm_memphy',0_ip)

  case(3_ip)
     !
     ! Allocate mixture fraction vector for CMC model
     !
     nullify(Z_CMC_chm)
     nullify(diff_Z_CMC_chm)
     call memory_alloca(mem_modul(1:2,modul),'MIXT_FRAC_VEC_CMC_CHM','chm_memall',Z_CMC_chm,nZ_CMC_chm)
     call memory_alloca(mem_modul(1:2,modul),'MIXT_FRAC_DIFFERENCE_VEC_CMC_CHM','chm_memall',diff_Z_CMC_chm,nZ_CMC_chm-1_ip)
     Z_CMC_chm(1:nZ_CMC_chm) = 0.0_rp
     diff_Z_CMC_chm(1:nZ_CMC_chm-1) = 0.0_rp
     
   
  case(4_ip)
     !
     ! Allocate memory for reactive scalars (temperature, enthalpy and mass fractions) at boundaries (mixture fraction space) for CMC model
     !
     nullify(T_bc_CMC_chm)
     nullify(react_scal_bc_CMC_chm)
     call memory_alloca(mem_modul(1:2,modul),'TEMPERATURE_BC_CMC_CHM','chm_memall',T_bc_CMC_chm,2_ip)
     call memory_alloca(mem_modul(1:2,modul),'REACTIVE_SCALARS_BC_CMC_CHM','chm_memall',react_scal_bc_CMC_chm,nclas_chm+1_ip,2_ip)
     T_bc_CMC_chm(1:2) = 298.0_rp   ! Temperature
     react_scal_bc_CMC_chm(1:nclas_chm+1,1:2) = 0.0_rp    ! Enthalpy and mass fractions

  
  case(5_ip)
     !
     ! Allocate memory for mixture fraction and segregation factor for AMC model for CMC model
     !
     nullify(Z_AMC_CMC_chm)
     nullify(S_AMC_CMC_chm)
     nullify(Xintegrated_table_AMC_CMC_chm)
     call memory_alloca(mem_modul(1:2,modul),'MIXTURE_FRAC_AMC_MODEL_CMC_CHM','chm_memall',Z_AMC_CMC_chm,nZ_AMC_CMC_chm)
     call memory_alloca(mem_modul(1:2,modul),'SEGREGATION_FACTOR_AMC_MODEL_CMC_CHM','chm_memall',S_AMC_CMC_chm,nS_AMC_CMC_chm)
     call memory_alloca(mem_modul(1:2,modul),'SCAL_DISSIP_RATE_INTEGRATED_PROFILE_CMC_CHM','chm_memall',Xintegrated_table_AMC_CMC_chm,nZ_AMC_CMC_chm,nS_AMC_CMC_chm)

     Z_AMC_CMC_chm(1:nZ_AMC_CMC_chm) = 0.0_rp
     S_AMC_CMC_chm(1:nS_AMC_CMC_chm) = 0.0_rp
     Xintegrated_table_AMC_CMC_chm(1:nZ_AMC_CMC_chm,1:nS_AMC_CMC_chm) = 0.0_rp    


  case(6_ip)
     !
     ! Allocate memory for scalar dissipation rate normalized profile
     !
     nullify(Xnormalized_prof_CMC_chm)
     call memory_alloca(mem_modul(1:2,modul),'SCALAR_DISSIP_RATE_NORMALIZED_PROF_CMC_CHM','chm_memall',Xnormalized_prof_CMC_chm,nZ_CMC_chm)
     Xnormalized_prof_CMC_chm(1:nZ_CMC_chm) = 0.0_rp

  case(7_ip)
     !
     ! Allocate memory for inert and equilibrium profiles and molecular weight.
     ! Note that W_k is only allocated for the master while later in chm_memall
     ! is allocated for the slaves
     !
     nullify(W_k)
     nullify(rscal_inert_CMC_chm)
     nullify(rscal_equil_CMC_chm)
     nullify(temp_inert_CMC_chm)
     nullify(temp_equil_CMC_chm)

     call memory_alloca(mem_modul(1:2,modul),'MOLECULAR_WEIGHT_CMC_CHM','chm_memall',W_k,nclas_chm)
     call memory_alloca(mem_modul(1:2,modul),'INERT_PROFILE_CMC_CHM','chm_memall',rscal_inert_CMC_chm,nZ_CMC_chm,nclas_chm+1)
     call memory_alloca(mem_modul(1:2,modul),'EQUILIBRIUM_PROFILE_CMC_CHM','chm_memall',rscal_equil_CMC_chm,nZ_CMC_chm,nclas_chm+1)
     call memory_alloca(mem_modul(1:2,modul),'INERT_TEMPERATURE_CMC_CHM','chm_memall',temp_inert_CMC_chm,nZ_CMC_chm)
     call memory_alloca(mem_modul(1:2,modul),'EQUILIBRIUM_TEMPERATURE_CMC_CHM','chm_memall',temp_equil_CMC_chm,nZ_CMC_chm)

     W_k(1:nclas_chm)                                = 0.0_rp
     rscal_inert_CMC_chm(1:nZ_CMC_chm,1:nclas_chm+1) = 0.0_rp
     rscal_equil_CMC_chm(1:nZ_CMC_chm,1:nclas_chm+1) = 0.0_rp
     temp_inert_CMC_chm(1:nZ_CMC_chm)                = 0.0_rp
     temp_equil_CMC_chm(1:nZ_CMC_chm)                = 0.0_rp

  end select

end subroutine chm_memphy
 
