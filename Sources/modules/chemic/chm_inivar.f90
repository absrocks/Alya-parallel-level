subroutine chm_inivar(itask)
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_inivar
  ! NAME 
  !    chm_inivar
  ! DESCRIPTION
  !    This routine initializes some variables
  ! USES
  ! USED BY
  !    chm_turnon
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_domain
  use def_chemic
  use def_solver
  use def_kintyp
  use mod_ADR,             only : ADR_initialize_type
  use mod_ADR,             only : ADR_check_and_compute_data
  use mod_ADR,             only : ADR_allocate_projections_bubble_sgs
  use mod_arrays,          only : arrays_register
  use mod_chm_rk_explicit, only : chm_rk_explicit_initialization
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ispec,iclas
  integer(ip)             :: nclas_max
  integer(ip), save       :: nbou_set_vars
  type(ADR_typ)           :: ADR_read  ! ADR type for reading

  select case(itask)

  case(0_ip)
     !
     ! Initialize modules
     !
     call chm_rk_explicit_initialization()     
     !
     ! Postprocess
     !
     call arrays_register(  1_ip,(/'CONCE','SCALA','NPOIN','PRIMA'/),conce,ENTITY_POSITION=1_ip,TIME_POSITION=3_ip,EXCLUDE_RST=(/2_ip,3_ip/))
     postp(1) % wopos (1,2)  = 'SOURC'  ! Source terms combustion
     postp(1) % wopos (1,3)  = 'VISCO'  ! Viscosity
     postp(1) % wopos (1,4)  = 'SPHEA'  ! Specific heat
     postp(1) % wopos (1,5)  = 'CONDU'  ! Heat conductivity
     postp(1) % wopos (1,6)  = 'ENTHA'  ! Enthalpy
     postp(1) % wopos (1,7)  = 'DIVEN'  ! Enthalpy transport source term
     postp(1) % wopos (1,8)  = 'CHEMI'  ! Chemical only Heat source term
     postp(1) % wopos (1,9)  = 'SUMCO'  ! Sum of concentration
     postp(1) % wopos (1,10) = 'MOLEC'  ! Molecular weight
     postp(1) % wopos (1,11) = 'TEMPE'  ! Temperature for low-mach Flamelet combustion model
     postp(1) % wopos (1,12) = 'AVY  '  ! Averaged reaction progress variable Yc or C
     postp(1) % wopos (1,13) = 'AVYV '  ! Averaged variance of reaction progress variable Yc or C
     postp(1) % wopos (1,14) = 'AVZV '  ! Averaged variance of mixture fraction Z
     postp(1) % wopos (1,15) = 'AVCHM'  ! Averaged chemical heat
     postp(1) % wopos (1,16) = 'AVZ  '  ! Averaged mixture fraction Z
     postp(1) % wopos (1,17) = 'AVZ2 '  ! Averaged squared of mixture fraction
     postp(1) % wopos (1,18) = 'SPEC1'  ! Species post-processing Flamelet model, radiation
     postp(1) % wopos (1,19) = 'SPEC2'  ! Species post-processing Flamelet model, radiation
     postp(1) % wopos (1,20) = 'RADIA'  ! Radiation source term Flamelet model
     postp(1) % wopos (1,21) = 'ZGRMA'  ! Maximum mixture fraction variance in the direction perpendicular to the flame
     postp(1) % wopos (1,22) = 'PHI  '  ! Weighting parameter for the hybrid model
     postp(1) % wopos (1,23) = 'XYR  '  ! Scalar dissipation rate of Yc (resolved part)
     postp(1) % wopos (1,24) = 'XZR  '  ! Scalar dissipation rate of Z  (resolved part)
     postp(1) % wopos (1,25) = 'XYS  '  ! Scalar dissipation rate of Yc (subgrid part)
     postp(1) % wopos (1,26) = 'XZS  '  ! Scalar dissipation rate of Z  (subgrid part)
     postp(1) % wopos (1,27) = 'AVXYR'  ! Average scalar dissipation rate of Yc (resolved part)
     postp(1) % wopos (1,28) = 'AVXZR'  ! Average scalar dissipation rate of Z  (resolved part)
     postp(1) % wopos (1,29) = 'AVXYS'  ! Average scalar dissipation rate of Yc (subgrid part)
     postp(1) % wopos (1,30) = 'AVXZS'  ! Average scalar dissipation rate of Z  (subgrid part)
     postp(1) % wopos (1,31) = 'AVY2 '  ! Average progress variable squared Yc*Yc or C*C 
     postp(1) % wopos (1,32) = 'AVL  '  ! Average liquid volume fraction phi_L
     postp(1) % wopos (1,33) = 'AVL2 '  ! Average liquid volume fraction squared phi_L*phi_L
     postp(1) % wopos (1,34) = 'AVS  '  ! Average interface surface density Sigma
     postp(1) % wopos (1,35) = 'AVS0 '  ! Average interface surface density Sigma_0
     postp(1) % wopos (1,36) = 'AVD32'  ! Average Sauter mean diameter
     postp(1) % wopos (1,37) = 'AVDEN'  ! Average density
     postp(1) % wopos (1,38) = 'SIGMA'  ! Interface surface density Sigma
     postp(1) % wopos (1,39) = 'SIGM0'  ! Interface surface density Sigma_0
     postp(1) % wopos (1,40) = 'D32  '  ! Sauter mean diameter
     postp(1) % wopos (1,41) = 'GRADY'  ! Gradient species mass fractions
     postp(1) % wopos (1,42) = 'HTRAN'  ! Enthalpy transport by diffusion 
     postp(1) % wopos (1,43) = 'ELEMH'  ! Elemental fraction H 
     postp(1) % wopos (1,44) = 'ELEMO'  ! Elemental fraction O 
     postp(1) % wopos (1,45) = 'ELEMC'  ! Elemental fraction C 
     postp(1) % wopos (1,46) = 'MASSK'  ! Mass source from partis 
     postp(1) % wopos (1,47) = 'ENVIC'  ! Entropy viscosity maximum among species 
     postp(1) % wopos (1,48) = 'SCONC'  ! Scaled control variables
     postp(1) % wopos (1,49) = 'NREAC'  ! Sum of Reactions (finite Rate)
     postp(1) % wopos (1,50) = 'MIXFR'  ! Mixture Fraction (finite rate)
     postp(1) % wopos (1,51) = 'HRR  '  ! Instantaneous heat release
     postp(1) % wopos (1,52) = 'AVHRR'  ! Averaged heat release
     postp(1) % wopos (1,53) = 'AVMSK'  ! Average mass source term from spray
     postp(1) % wopos (1,54) = 'POSTT'  ! Posttab properties
     postp(1) % wopos (1,55) = 'AVPOT'  ! Average posttab properties
     postp(1) % wopos (1,56) = 'MASCN'  ! Conditional mass fractions for CMC model
     postp(1) % wopos (1,57) = 'MASUN'  ! Unconditional mass fractions for CMC model
     postp(1) % wopos (1,58) = 'ENTCN'  ! Conditional enthalpy for CMC model
     postp(1) % wopos (1,59) = 'ENTUN'  ! Unconditional enthalpy for CMC model
     postp(1) % wopos (1,60) = 'TEMCN'  ! Conditional temperature for CMC model
     postp(1) % wopos (1,61) = 'TEMUN'  ! Unconditional temperature for CMC model
     postp(1) % wopos (1,62) = 'SRCCN'  ! Conditional chemical source terms for CMC model
     postp(1) % wopos (1,63) = 'SRCUN'  ! Unconditional chemical source terms for CMC model

     postp(1) % wopos (2,2)  = 'SCALA'
     postp(1) % wopos (2,3)  = 'SCALA'
     postp(1) % wopos (2,4)  = 'SCALA'
     postp(1) % wopos (2,5)  = 'SCALA'
     postp(1) % wopos (2,6)  = 'SCALA'
     postp(1) % wopos (2,7)  = 'SCALA'
     postp(1) % wopos (2,8)  = 'SCALA'
     postp(1) % wopos (2,9)  = 'SCALA'
     postp(1) % wopos (2,10) = 'SCALA'
     postp(1) % wopos (2,11) = 'SCALA'
     postp(1) % wopos (2,12) = 'SCALA'
     postp(1) % wopos (2,13) = 'SCALA'
     postp(1) % wopos (2,14) = 'SCALA'
     postp(1) % wopos (2,15) = 'SCALA'
     postp(1) % wopos (2,16) = 'SCALA'
     postp(1) % wopos (2,17) = 'SCALA'
     postp(1) % wopos (2,18) = 'SCALA'
     postp(1) % wopos (2,19) = 'SCALA'
     postp(1) % wopos (2,20) = 'SCALA'
     postp(1) % wopos (2,21) = 'SCALA'
     postp(1) % wopos (2,22) = 'SCALA'
     postp(1) % wopos (2,23) = 'SCALA'
     postp(1) % wopos (2,24) = 'SCALA'
     postp(1) % wopos (2,25) = 'SCALA'     
     postp(1) % wopos (2,26) = 'SCALA'
     postp(1) % wopos (2,27) = 'SCALA'
     postp(1) % wopos (2,28) = 'SCALA'
     postp(1) % wopos (2,29) = 'SCALA'     
     postp(1) % wopos (2,30) = 'SCALA'
     postp(1) % wopos (2,31) = 'SCALA'
     postp(1) % wopos (2,32) = 'SCALA'
     postp(1) % wopos (2,33) = 'SCALA'
     postp(1) % wopos (2,34) = 'SCALA'
     postp(1) % wopos (2,35) = 'SCALA'
     postp(1) % wopos (2,36) = 'SCALA'
     postp(1) % wopos (2,37) = 'SCALA'
     postp(1) % wopos (2,38) = 'SCALA'
     postp(1) % wopos (2,39) = 'SCALA'
     postp(1) % wopos (2,40) = 'SCALA'
     postp(1) % wopos (2,41) = 'VECTO'
     postp(1) % wopos (2,42) = 'VECTO'
     postp(1) % wopos (2,43) = 'SCALA'
     postp(1) % wopos (2,44) = 'SCALA'
     postp(1) % wopos (2,45) = 'SCALA'
     postp(1) % wopos (2,46) = 'SCALA'
     postp(1) % wopos (2,47) = 'SCALA'
     postp(1) % wopos (2,48) = 'SCALA'
     postp(1) % wopos (2,49) = 'SCALA'
     postp(1) % wopos (2,50) = 'SCALA'
     postp(1) % wopos (2,51) = 'SCALA'
     postp(1) % wopos (2,52) = 'SCALA'
     postp(1) % wopos (2,53) = 'SCALA'
     postp(1) % wopos (2,54) = 'SCALA'
     postp(1) % wopos (2,55) = 'SCALA'
     postp(1) % wopos (2,56) = 'SCALA'
     postp(1) % wopos (2,57) = 'SCALA'
     postp(1) % wopos (2,58) = 'SCALA'
     postp(1) % wopos (2,59) = 'SCALA'
     postp(1) % wopos (2,60) = 'SCALA'
     postp(1) % wopos (2,61) = 'SCALA'
     postp(1) % wopos (2,62) = 'SCALA'
     postp(1) % wopos (2,63) = 'SCALA'

     !
     ! Nodal set variables 
     !
     postp(1) % wonse (1)    = 'CONCE'  

     !
     ! Elemental set variables
     !
     postp(1) % woese (1)    = 'MASS ' !  int_V rho dV
     postp(1) % woese (2)    = 'HEATR' !  int_V omega_h dV, assuming omega_h is in J/(m3 s)
     postp(1) % woese (3:6)  = 'CONCE' !  int_V rho Y_k dV 

     !
     ! Boundary set variables
     !
     postp(1) % wobse (1)    = 'MASS ' ! int_S rho*u.n ds
     postp(1) % wobse (2:9)  = 'CONCE' ! <Yk> = int_S rho*u* Y_k dS / int_S rho*u dS, for k = 1,...,8
     nbou_set_vars = 1                 ! Need to set this for later reentry when we know nspec

     !
     ! Witness variables
     !
     postp(1) % wowit (1:8)  = 'CONCE' ! Species 1 to 8 (update the maximum number of species here)     

     postp(1) % wowit (8+1)  = 'XZR  ' ! Resolved scalar dissipation of mixture fraction
     postp(1) % wowit (8+2)  = 'XZS  ' ! Subgrid scalar dissipation of mixture fraction
     postp(1) % wowit (8+3)  = 'XYR  ' ! Resolved scalar dissipation of progress variable
     postp(1) % wowit (8+4)  = 'XYS  ' ! Subgrid scalar dissipation of progress variable
     postp(1) % wowit (8+5)  = 'D32  ' ! Sauter mean diameter Eulerian atomization 
     postp(1) % wowit (8+6)  = 'SIGMA' ! Liquid surface density S = S_0 + S' 
     postp(1) % wowit (8+7)  = 'SIGM0' ! Mean liquid surface density S_0

     postp(1) % wowit ((8+8):(8+16))  = 'SCONC' ! Scaled control variable 
     !
     ! Solver
     !     
     call soldef(-2_ip)
     solve(1) % wprob       = 'CONCENTRATION'
     solve(1) % kfl_solve   = 1
     
     solve(2) % wprob       = 'CONSISTENT'           ! Consistent matrix
     solve(2) % kfl_solve   = 1
     solve(2) % kfl_iffix   = 2

     !
     ! Others
     !
     cputi_chm = 0.0_rp  ! CPU times
     dtmat_chm = 0.0_rp  ! Matrix-based time step

     !
     ! Nullify pointers
     !
     nullify(kfl_fixno_chm)
     nullify(kfl_fixbo_chm)    
     nullify(kfl_funno_chm)    
     nullify(kfl_funtn_chm)     
     nullify(bvess_chm)            

     ! 
     ! Flamelet combustion model 
     ! 
     nullify(avY_chm)
     nullify(avYv_chm)
     nullify(avZ_chm)
     nullify(avZv_chm)
     nullify(avZ2_chm)
     nullify(avY2_chm)

     nullify(xYr_chm)
     nullify(xYs_chm)
     nullify(xZr_chm)
     nullify(xZs_chm)

     nullify(avxYr_chm)
     nullify(avxYs_chm)
     nullify(avxZr_chm)
     nullify(avxZs_chm)
     nullify(avposttab_chm)

     nullify(zgradmax_chm)
     nullify(phi_chm)

     ! 
     ! ELSA model
     !
     nullify(avL_chm)
     nullify(avL2_chm)
     nullify(avS_chm)
     nullify(avS0_chm)
     nullify(avd32_chm)
     nullify(Sigma_chm)
     nullify(Sigm0_chm)
     nullify(d32_chm)
     nullify(lap_phi_levSet_chm)
     nullify(grad_phi_levSet_chm)
     nullify(grad_phic_levSet_chm)
     
     !
     ! Droplet Identification
     !
     nullify(volume_drop_chm)
     nullify(diameter_drop_chm)
     nullify(centroid_drop_chm)
     nullify(compactness2_drop_chm)
     nullify(volume_cluster_chm)

     ! 
     ! Others
     !
     nullify(avchm_chm)
     nullify(avden_chm)
     nullify(avmsk_chm)
     nullify(dt_rho_chm)
     nullify(dt_chm)
     
     nullify(grad_Yk)
     nullify(grad_T)
     nullify(enthalpy_transport_nodes)
     nullify(coeff_cp_k)
     nullify(W_k)
     nullify(Y_k_n)
     nullify(React_ind)
     nullify(elem_c)
     nullify(elem_h)
     nullify(elem_o)
     nullify(elem_n)
     nullify(src_chm)
     nullify(Corr_chm)
     nullify(mixfr_chm)
     nullify(prog_var_chm)
     nullify(sum_reac_chm)
     nullify(hrr_chm)
     nullify(hrr_avg_chm)

  case(1_ip)

!DMM     if( nbou_set_vars+nclas_chm > nvars ) call runend('CHM_INIVAR')
     kfl_rsta2_chm= .false.                         ! restarting bdf file?

  case(3_ip)

     if(kfl_timei_chm==0) then                      ! Time integration
        dtinv_chm=1.0_rp
     else
        kfl_timei=1
     end if
     momod(modul) % kfl_stead = 0
     kfl_grdif_chm = 0
!!DMM     nspec_chm     = nclas_chm

     !
     ! nspec is the number of species for master
     !
     nspec = nspec_chm

     if( ISLAVE ) then
        do ispec = 1,nspec_chm
           speci(ispec)%name = ''
        enddo
     endif
     !
     ! Time integration definition
     !
     kfl_tiaor_chm = kfl_tiacc_chm                  ! Time accuracy: save original value
     if( kfl_timei_chm == 1 ) then                  

        if( kfl_tisch_chm == 1 ) then
           ncomp_chm = 3                            ! Trapezoidal rule
        else if( kfl_tisch_chm == 2 ) then
           ncomp_chm = 2+kfl_tiacc_chm              ! BDF scheme
        else if( kfl_tisch_chm == 3 ) then
           ncomp_chm = 2+2                          ! AB scheme
        else if( kfl_tisch_chm == 4 ) then
           ncomp_chm = 2+2                          ! RK scheme
        end if

     else
        ncomp_chm = 2     
     end if

     ittot_chm = 0                                  ! Others

     !
     ! To solve equation of state when nastin is not activated
     !
     if ( kfl_model_chm == 1 .and. prthe_chm /= 0 ) then
        prthe(1)      = prthe_chm                    ! Low-Mach: p0
        prthe(2)      = prthe_chm                    ! Low-Mach: p0^{n}
        prthe(3)      = prthe_chm                    ! Low-Mach: p0^{n-1}
        prthe(4)      = prthe_chm                    ! Low-Mach: p0^{0}
     end if

     iclai_chm = 1                                  ! Jacobi starting class
     iclaf_chm = nclas_chm                          ! Jacobi final class
     avtim_chm = 0.0_rp                             ! Accumulated time for time-averaging variables
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
     ADR_read % kfl_nonlinear_sgs      =  0                ! Deprecated
     ADR_read % kfl_time_sgs           =  0                ! Deprecated
     ADR_read % kfl_time_bubble        =  kfl_tibub_chm
     ADR_read % kfl_time_scheme        =  kfl_tisch_chm
     ADR_read % kfl_time_order         =  kfl_tiacc_chm    
     ADR_read % kfl_manufactured       =  0                ! Exact solutions not available in chemic
     ADR_read % kfl_length             =  kfl_ellen_chm
     ADR_read % kfl_first_order_sgs    =  1                ! Related to ADR_chm % kfl_time_order
     ADR_read % number_euler_steps     =  0                ! Deprecated
     ADR_read % lun_output4            =  int(momod(modul) % lun_outpu,4)
     ADR_read % bemol                  =  bemol_chm
     ADR_read % tau_parameters(1:3)    =  staco_chm(1:3)
     ADR_read % shock                  =  shock_chm

     if (kfl_model_chm == 4) then
        nclas_max = nvar_CMC_chm
     else
        nclas_max = nclas_chm
     end if

     !
     ! allocate
     !
     allocate(ADR_chm(nclas_max))
     !
     ! extension to ADR_chm(iclas) 
     !  
     do iclas=1,nclas_max
        call ADR_initialize_type(ADR_chm(iclas))
        ADR_chm(iclas) = ADR_read
        call ADR_check_and_compute_data(ADR_chm(iclas))
     end do
     
  end select

end subroutine chm_inivar
