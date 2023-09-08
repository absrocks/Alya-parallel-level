subroutine chm_memall()
  !-----------------------------------------------------------------------
  !****f* partis/chm_memall
  ! NAME 
  !    chm_memall
  ! DESCRIPTION
  !    This routine allocates memory for the arrays needed
  ! USES
  ! USED BY
  !    chm_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_domain
  use def_solver
  use def_chemic
  use mod_memchk
  use mod_ADR, only : FULL_OSS
  use mod_ADR, only : A_OSS  
  use mod_ADR, only : AR_OSS 
  use mod_ADR, only : BUBBLE
  use mod_ADR, only : ADR_initialize_type
  use mod_ADR, only : ADR_check_and_compute_data
  use mod_ADR, only : ADR_allocate_projections_bubble_sgs
  implicit none
  integer(ip) :: iclas,ielem,pelty,pgaus
  integer(4)  :: istat
  integer(ip) :: monolithic_dim

  if( kfl_coupl_chm == 2 ) then
     monolithic_dim = nclas_chm
  else
     monolithic_dim = 1
  end if

  if( INOTMASTER ) then
     !
     ! CONCE: concentration
     ! 
     allocate(conce(npoin,nspec_chm,ncomp_chm),stat=istat) !!**
     call memchk(zero,istat,mem_modul(1:2,modul),'CONCE','chm_memall',conce)
     !
     ! COSGS: Subgrid scale concentration 
     !
     if( kfl_sgsti_chm == 1 ) then     
        allocate(cosgs(nelem,nclas_chm),stat=istat)  !!**
        call memchk(zero,istat,mem_modul(1:2,modul),'COSGS','chm_memall',cosgs)
        do iclas = 1,nclas_chm
           do ielem = 1,nelem
              pelty = ltype(ielem)
              pgaus = ngaus(pelty)
              allocate(cosgs(ielem,iclas)%a(pgaus,2),stat=istat) !!**
              call memchk(zero,istat,mem_modul(1:2,modul),'COSGS','chm_memall',cosgs(ielem,iclas)%a)
           end do
        end do
     end if
     !
     ! PROJEC_CHM: Projection
     !
     if( kfl_stabi_chm >= 1 ) then
        allocate(proje_chm(npoin,nclas_chm),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'PROJE_CHM','chm_memall',proje_chm)        
     end if
     !
     ! VELOC_CHM: Advection
     !
     if( kfl_model_chm == 2 .or. ( kfl_advec_chm /= 0 .and. kfl_advec_chm /= -2 ) ) then
        allocate(veloc_chm(ndime,npoin),stat=istat)                                         ! Velocity          (from meteo file or function)
        call memchk(zero,istat,mem_modul(1:2,modul),'VELOC_CHM','chm_memall',veloc_chm)
     end if
     !
     ! METEO model
     !
     if( kfl_model_chm == 2 ) then
  
        if( lawde_chm == -1 ) then
           allocate(densi_chm(npoin),stat=istat)                                            ! Density           (from meteo file)
           call memchk(zero,istat,mem_modul(1:2,modul),'DENSI_CHM','chm_memall',densi_chm)
        end if

        if( lawte_chm == -1 ) then
           allocate(tempe_chm(npoin),stat=istat)                                            ! Temperature       (from meteo file)
           call memchk(zero,istat,mem_modul(1:2,modul),'TEMPE_CHM','chm_memall',tempe_chm)
        end if

        allocate(vfric_chm(npoin),stat=istat)                                               ! Friction velocity (from meteo file)
        call memchk(zero,istat,mem_modul(1:2,modul),'VFRIC_CHM','chm_memall',vfric_chm)
        allocate(hepbl_chm(npoin),stat=istat)                                               ! BL height         (from meteo file)
        call memchk(zero,istat,mem_modul(1:2,modul),'HEPBL_CHM','chm_memall',hepbl_chm)
        allocate(walld_chm(npoin),stat=istat)                                               ! Wall distance     (from meteo file)
        call memchk(zero,istat,mem_modul(1:2,modul),'WALLD_CHM','chm_memall',walld_chm)
        allocate(lmoni_chm(npoin),stat=istat)                                               ! M-O length        (from meteo file)
        call memchk(zero,istat,mem_modul(1:2,modul),'LMONI_CHM','chm_memall',lmoni_chm)
        allocate(tmrat_chm(nclas_chm,npoin),stat=istat)                                     ! Total Mass rate   (from source file)
        call memchk(zero,istat,mem_modul(1:2,modul),'TMRAT_CHM','chm_memall',tmrat_chm)
        allocate(vterm_chm(npoin,nclas_chm),stat=istat)                                     ! Terminal velocity
        call memchk(zero,istat,mem_modul(1:2,modul),'VTERM_CHM','chm_memall',vterm_chm)
        allocate(accum_chm(npoin),stat=istat)                                               ! Accumulation
        call memchk(zero,istat,mem_modul(1:2,modul),'ACCUM_CHM','chm_memall',accum_chm)
     end if
     !
     ! Mechano-biological model
     !
     if( kfl_model_chm == 3 .and. wprob_chm == 'OSTE1' ) then
        allocate(proad_chm(npoin),stat=istat)                                               ! Protein adsorption
        call memchk(zero,istat,mem_modul(1:2,modul),'PROAD_CHM','chm_memall',proad_chm)        
     end if
     !
     ! Combustion
     !
     if( kfl_model_chm == 4 .or. kfl_model_chm == 5) then

        allocate(enthalpy_transport(nelem),stat=istat)  
        call memchk(zero,istat,mem_modul(1:2,modul),'enthalpy_transport','chm_memall',enthalpy_transport)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           allocate(enthalpy_transport(ielem)%a(ndime,pgaus),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'enthalpy_transport','chm_memall',enthalpy_transport(ielem)%a)
        end do
        allocate(div_enthalpy_transport(nelem),stat=istat)  
        call memchk(zero,istat,mem_modul(1:2,modul),'div_enthalpy_transport','chm_memall',div_enthalpy_transport)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           allocate(div_enthalpy_transport(ielem)%a(pgaus),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'div_enthalpy_transport','chm_memall',div_enthalpy_transport(ielem)%a)
        end do
        allocate(chemical_heat(nelem),stat=istat)  
        call memchk(zero,istat,mem_modul(1:2,modul),'chemical_heat','chm_memall',chemical_heat)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           allocate(chemical_heat(ielem)%a(pgaus),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'chemical_heat','chm_memall',chemical_heat(ielem)%a)
        end do
        !
        ! Species fields for radiation for CFI model
        !     
        allocate(rspec_chm(2,npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'rspec_chm','chm_memphy',rspec_chm)

        allocate(radiative_heat(nelem),stat=istat) 
        call memchk(zero,istat,mem_modul(1:2,modul),'radiative_heat','chm_memphy',radiative_heat)
        do ielem = 1,nelem 
           pelty = ltype(ielem)  
           pgaus = ngaus(pelty)      
           allocate(radiative_heat(ielem)%a(pgaus),stat=istat) 
           call memchk(zero,istat,mem_modul(1:2,modul),'radiative_heat','chm_memphy',radiative_heat(ielem)%a)
        end do
        !
        ! Dynamic thickened flame model (DTFLES)
        !

        !
        ! Thickening factor
        !
        allocate(tfles_factor(nelem),stat=istat)  
        call memchk(zero,istat,mem_modul(1:2,modul),'tfles_factor','chm_memall',tfles_factor)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           allocate(tfles_factor(ielem)%a(pgaus),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'tfles_factor','chm_memall',tfles_factor(ielem)%a)
        end do
        !
        ! Flame sensor
        !
        allocate(tfles_sensor(nelem),stat=istat)  
        call memchk(zero,istat,mem_modul(1:2,modul),'tfles_sensor','chm_memall',tfles_sensor)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           allocate(tfles_sensor(ielem)%a(pgaus),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'tfles_sensor','chm_memall',tfles_sensor(ielem)%a)
        end do
        !
        ! Subgrid scale wrinkling 
        !
        allocate(tfles_sgseff(nelem),stat=istat)  
        call memchk(zero,istat,mem_modul(1:2,modul),'tfles_sgseff','chm_memall',tfles_sgseff)
        do ielem = 1,nelem
           pelty = ltype(ielem)
           pgaus = ngaus(pelty)
           allocate(tfles_sgseff(ielem)%a(pgaus),stat=istat)
           call memchk(zero,istat,mem_modul(1:2,modul),'tfles_sgseff','chm_memall',tfles_sgseff(ielem)%a)
        end do

        allocate(visck(npoin,nclas_chm),stat=istat)   !!**
        call memchk(zero,istat,mem_modul(1:2,modul),'VISCK','chm_memall',visck)        
        allocate(condk(npoin,nclas_chm),stat=istat)  !!**
        call memchk(zero,istat,mem_modul(1:2,modul),'CONDK','chm_memall',condk)        
        allocate(sphek(npoin,nclas_chm),stat=istat)  !!**
        call memchk(zero,istat,mem_modul(1:2,modul),'SPHEK','chm_memall',sphek)      
        allocate(sphec(npoin,6,2),stat=istat)  !!**
        call memchk(zero,istat,mem_modul(1:2,modul),'SPHEC','chm_memall',sphec) 
        allocate(entha_chm(npoin,nclas_chm),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'ENTHA_CHM','chm_memall',entha_chm)       
        allocate(massk(npoin,nclas_chm),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'MASSK_CHM','chm_memall',massk)
        allocate(yscale_chm(npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'YSCALE_CHM','chm_memall',yscale_chm)       
        allocate(cvar_chm(npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'CAVR_CHM','chm_memall',cvar_chm)       

        allocate(flsen_chm(npoin),stat=istat)   !!**
        call memchk(zero,istat,mem_modul(1:2,modul),'FLSEN_CHM','chm_memall',flsen_chm)
        allocate(flspe_chm(npoin),stat=istat)   !!**
        call memchk(zero,istat,mem_modul(1:2,modul),'FLSPE_CHM','chm_memall',flspe_chm)
        allocate(flfac_chm(npoin),stat=istat)   !!**
        call memchk(zero,istat,mem_modul(1:2,modul),'FLFAC_CHM','chm_memall',flfac_chm)
        allocate(flsgs_chm(npoin),stat=istat)   !!**
        call memchk(zero,istat,mem_modul(1:2,modul),'FLSGS_CHM','chm_memall',flsgs_chm)
        allocate(flthi_chm(npoin),stat=istat)   !!**
        call memchk(zero,istat,mem_modul(1:2,modul),'FLTHI_CHM','chm_memall',flthi_chm)
        allocate(equiv_chm(npoin),stat=istat)   !!**
        call memchk(zero,istat,mem_modul(1:2,modul),'EQUIV_CHM','chm_memall',equiv_chm)
        allocate(avtem_chm(npoin),stat=istat)   
        call memchk(zero,istat,mem_modul(1:2,modul),'AVTEM_CHM','chm_memall',avtem_chm)
        allocate(avcon_chm(npoin),stat=istat)   
        call memchk(zero,istat,mem_modul(1:2,modul),'AVCON_CHM','chm_memall',avcon_chm)
        allocate(avco2_chm(npoin),stat=istat)   
        call memchk(zero,istat,mem_modul(1:2,modul),'AVCO2_CHM','chm_memall',avco2_chm)
        allocate(avvar_chm(npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'AVVAR_CHM','chm_memall',avvar_chm)
        allocate(avime_chm(npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'AVIME_CHM','chm_memall',avime_chm)
        allocate(avchm_chm(npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'AVCHM_CHM','chm_memall',avchm_chm)
        allocate(avmix_chm(npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'AVMIX_CHM','chm_memall',avmix_chm)
        allocate(avmi2_chm(npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'AVMI2_CHM','chm_memall',avmi2_chm)
        allocate(lescl(npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'LESCL','chm_memall',lescl)
        allocate(encfi(npoin),stat=istat)
        call memchk(zero,istat,mem_modul(1:2,modul),'ENCFI','chm_memall',encfi)

        avtem_chm = 0.0_rp ! Initilization
        avcon_chm = 0.0_rp
        avco2_chm = 0.0_rp
        avvar_chm = 0.0_rp
        avime_chm = 0.0_rp
        avchm_chm = 0.0_rp
        avmix_chm = 0.0_rp
        avmi2_chm = 0.0_rp

        !Kernel shared variables
        allocate(sphea(npoin),stat=istat)   !!**
        call memchk(zero,istat,mem_modul(1:2,modul),'SPHEA','chm_memall',sphea)       
        allocate(wmean(npoin,ncomp_chm),stat=istat)  !!**
        call memchk(zero,istat,mem_modul(1:2,modul),'WMEAN','chm_memall',wmean)       
        allocate(kcond(npoin),stat=istat)    !!**
        call memchk(zero,istat,mem_modul(1:2,modul),'KCOND','chm_memall',kcond)       
        if (.not.associated(visco)) then ! Just checking if not allocated before by nastal, should be done better
           allocate(visco(npoin,2),stat=istat)  !!**
           call memchk(zero,istat,mem_modul(1:2,modul),'VISCO','chm_memall',visco)       
        endif

     end if
     !----------------------------------------------------------------------
     !
     ! Solver
     !
     !----------------------------------------------------------------------
     !
     ! Memory
     !
     if( kfl_coupl_chm == 2 ) then
        solve(1) % ndofn = nspec_chm
     else
        solve(1) % ndofn = 1
     end if
     solve_sol => solve(1:)
     call soldef(4_ip)
     !
     ! Skyline matrix for ODE's
     !
     call chm_skymat()

  else

     allocate(conce(1,nclas_chm,ncomp_chm),stat=istat)   !!**
     call memchk(zero,istat,mem_modul(1:2,modul),'CONCE','chm_memall',conce)

  end if
  !
  ! Residual
  !
  allocate(ripts_chm(nspec_chm*monolithic_dim),stat=istat)
  call memchk(zero,istat,mem_modul(1:2,modul),'RIPTS_CHM','chm_memall',ripts_chm)

  !
  ! Boundary conditions
  !
  solve(1) % bvess     => bvess_chm
  solve(1) % kfl_fixno => kfl_fixno_chm
  
  !
  ! ADR type
  !
  call ADR_initialize_type(ADR_read)
  ADR_read % kfl_time_integration   =  kfl_timei_chm
  ADR_read % kfl_time_step_strategy =  kfl_timco
  ADR_read % kfl_stabilization      =  kfl_stabi_chm
  ADR_read % kfl_shock              =  kfl_shock_chm
  ADR_read % kfl_time_lumped        =  0
  ADR_read % kfl_tau_strategy       =  kfl_taust_chm
  ADR_read % kfl_laplacian          =  0 
  ADR_read % kfl_nonlinear_sgs      =  0                ! Non-linear subscale model is not implemented 
  ADR_read % kfl_time_sgs           =  kfl_sgsti_chm
  ADR_read % kfl_time_bubble        =  kfl_tibub_chm
  ADR_read % kfl_time_scheme        =  kfl_tisch_chm
  ADR_read % kfl_time_order         =  kfl_tiacc_chm    
  ADR_read % kfl_manufactured       =  0                ! Exact solutions not available in chemic
  ADR_read % kfl_length             =  kfl_ellen_chm
  ADR_read % kfl_first_order_sgs    =  1                ! Related to ADR_chm % kfl_time_order
  ADR_read % number_euler_steps     =  neule_chm
  ADR_read % lun_output4            =  int(momod(modul) % lun_outpu,4)
  ADR_read % bemol                  =  bemol_chm
  ADR_read % tau_parameters(1:3)    =  staco_chm(1:3)
  ADR_read % shock                  =  shock_chm

  !
  ! allocate
  !
  allocate(ADR_chm(nclas_chm))

  !
  ! extension to ADR_chm(iclas) 
  !  
  do iclas=1,nclas_chm
     call ADR_initialize_type(ADR_chm(iclas))

     ADR_chm(iclas) = ADR_read
     call ADR_check_and_compute_data(ADR_chm(iclas))
     call ADR_allocate_projections_bubble_sgs(ADR_chm(iclas))
  end do

end subroutine chm_memall

