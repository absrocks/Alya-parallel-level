subroutine chm_reaphy()
  !------------------------------------------------------------------------
  !****f* Chemic/chm_reaphy
  ! NAME 
  !    chm_reaphy
  ! DESCRIPTION
  !    This routine reads the physical problem definition
  !    NREAC_CHM ........................ Number of reactions
  !    LREAC(NREAC_CHM)%L(:) ............ List of reactions
  !    REACT_CHM(NCOEF_CHM,NREAC_CHM) ... Reaction coefficients 
  !  
  ! USES
  ! USED BY
  !    chm_turnon
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_chemic
  use def_domain
  use def_kermod,             only : gasco, kfl_prope, lookup_fw
  use mod_ecoute,             only : ecoute
  use mod_interp_tab,         only : tab_load_file
  use mod_interp_tab,         only : fw_allocate
  use mod_interp_tab,         only : tab_init_fw
  use mod_chm_operations_CMC, only : chm_inert_eq_CMC, chm_AMC_generate_Z_S_vectors_CMC, &
                                     chm_AMC_integrals_CMC, compute_Xnormalized_profile_CMC
#ifdef CANTERA
  use cantera,                only : importPhase, nSpecies, nReactions
#endif
  use iso_fortran_env,        only : iostat_eor,iostat_end
  use mod_messages,           only : livinf
  use mod_chm_sectional_soot_model, only: Vmax_ssm
  use mod_chm_sectional_soot_model, only: RadSoot_ssm 
  use mod_chm_sectional_soot_model, only: RhoC_ssm
  use mod_chm_sectional_soot_model, only: nsect_ssm
  use mod_chm_sectional_soot_model, only: nspec_ssm
  use mod_chm_sectional_soot_model, only: nclas_ssm
  use mod_chm_sectional_soot_model, only: gasCoupling_ssm 

  implicit none
  integer(ip)                   :: iclas,ireac,ispec,kfl_forwa,itab,nsize_spec_name,ipostvar,jpostvar, &
                                   index_zbc,icoun,line_f,state_f
  character(5)                  :: names(maxsp_chm)      ! Names of species 
  character(5)                  :: stoichio_names(3) = ''
  character(80)                 :: eq
  integer(ip)                   :: nbuff, ioerror, ii, jj, order_loc
  character(len=:), allocatable :: buffer, spec_name 
  integer(ip)                   :: stoichio_reaction,kfl_stoic
  real(rp)                      :: sum_Yk   
  real(rp),allocatable          :: Linear_temp_density_coeffs(:,:)
  integer(ip),allocatable       :: Field_ind(:) 
  logical                       :: found
  logical(lg)                   :: give_new_name

  if( INOTSLAVE ) then
     !
     ! Initializations (defaults)
     !
     kfl_model_chm    = 0                                 ! Defect evolution (1) Flamelet model / (2) Finite rate kinetics
     kfl_timei_chm    = 1                                 ! Existence of du/dt
     kfl_advec_chm    = 0                                 ! Existence of (a.grad)u
     kfl_diffu_chm    = 1                                 ! Existence of -div[k. grad(u)], and its options (include density)
     kfl_field_chm(1) = 0                                 ! Initialization by fields (=0 OFF, /=0 ON)
     kfl_field_chm(2) = 0
     kfl_field_chm(3) = 0
     kfl_key_chm      = 0                                 ! Number of Key species for Dynamic Chemistry 
     kfl_z_chm        = 0                                 ! Flag for z
     kfl_tdac_write_chm    = 0                            ! Flag to write reduced table 
     kfl_freq_chm     = 0                                 ! Number of interations between DAC reduction 
     dac_cor_chm      = 0.0_rp                            ! CODAC Corrleation
     bf_fuel_chm      = 0.0_rp                            ! B_F for fuel 
     bo_oxy_chm       = 0.0_rp                            ! B_O for Oxygen  
     kfl_radia_chm    = 0                                 ! radiation model activated?
     kfl_pfa_chm      = -1                                ! Reduction method based on PFA: =-1 NO REDUCTION, =0 EQUILIBRIUM, =1 PFA (static), =2 DAC, =3, CODAC
     dac_crit_chm     = 0.0_rp                            ! Critical PFA Value
     kfl_cotur_chm    = 0                                 ! Turbulence model RANS / LES 
     kfl_premix_chm   = 0                                 ! Combustion model: =0 (NON-PREMIXED) = 1 (PREMIXED)
     kfl_spray_chm    = 0                                 ! Activate Spray model
     kfl_ufpv_chm     = 0                                 ! Unsteady Flamelet Progress Variable model 0: original, 1: steady FPV, 2: omegaYc, 3: Yc_dot
     kfl_heat_loss_chm= 0                                 ! Heat loss included in table    
     kfl_varZ_chm     = 0           
     kfl_varYc_chm    = 0
     kfl_lookg_chm    = 1                                 ! Lookup of table properties, 0: on nodes, 1: on gauss points 
     kfl_soot_chm     = 0                                 ! Activation soot model, 0: OFF, 1: sectional method, -1: validation 
     nsect_ssm        = 0                                 ! Number of sections soot model
     kfl_tab_fw_chm   = -1                                ! Lookup framework
     kfl_post_fw_chm  = -1                                ! Postprocessing lookup framework
     order_loc        = 1_ip
     
     kfl_norma_chm = 0                                    ! Normalize concentrations
     kfl_forwa     = 0                                    ! Forward/backward reaction rate

     prthe_chm     = 0.0_rp                               ! Thermodynamic pressure (if NASTIN is not activated)
     nclas_chm     = 1                                    ! Number of particle classes
     radwt_chm     = 0.0_rp                               ! Wall temperature for radiation model
     sorad_chm     = 0.0_rp                               ! Source radius
     socen_chm     = 0.0_rp                               ! Source center
     nreac_chm     = 0                                    ! Number of reactions
     grnor_chm     = 9.80616_rp                           ! Gravity modulus

     kfl_droplet_id_chm            = 0_ip                 ! Droplet identification flag
     levelSet_threshold_chm        = 0.0_rp               ! Level Set threshold for droplet identification
     droplet_compactness_limit_chm = 1.0_rp               ! Compactness value below which a cluster won't be considered as a droplet
     droplet_max_diameter_chm      = 0.0_rp               ! Max. diameter above which a cluster won't be considered as a droplet
     droplet_h_factor_chm          = 0.0_rp               ! Mesh size factor to define a max. diameter as: droplet_h_factor_chm * h_mesh

     kfl_stoic     = 0_ip                                 ! Stoichiometric mass ratio defined or not
     stoichio_names    = ''                               !
     mechanism_path    = ''                               ! Path to mechanism: can be global path, relative path
     Red_spec          = ''                               ! Reduced Species 

     ! Variables for CMC model
     kfl_weigh_in_eq_CMC_chm = 0_ip                       ! Flag to weigh inert and equilibrium solutions when starting from scratch
     kfl_split_CFD_CMC = 1                                ! Flag to split CFD and CMC into two different executions: 0 CFD and CMC in same execution and 1 if CFD and CMC separated
     kfl_solve_enth_CMC_chm = 0                           ! 0 if enthalpy is not solved and 1 if it is transported
     nZ_CMC_chm = 3                                       ! Number of slices in mixture fraction space. Due to different reasons when using finite rate we take nZ_CMC_chm = 3
     nvar_CMC_chm = 0                                     ! Number of variables to be solved in CMC: nvar_CMC_chm = nclas_chm+1 (species + enthalpy)
     nZ_AMC_CMC_chm = 51                                  ! Number of mixture fractions for scalar dissipation rate integration (AMC model)
     nS_AMC_CMC_chm = 31                                  ! Number of segregation factors for scalar dissipation rate integration (AMC model)
     index_N2 = 0                                         ! Index in the mechanism for N2
     Zs_CMC_chm = 1.0_rp                                  ! Saturation mixture fraction
     Smax_AMC_CMC_chm = 0.3_rp                            ! Maximum segregation factor for AMC model
     S_threshold = 9.0e-3_rp                              ! Threshold segregation factor for integrations


     if (kfl_prope == 0) gasco = 287.17_rp ! Simone no cambies esta linea, Fer

     !
     ! Reach the section
     !
     call ecoute('chm_reaphy')
     do while(words(1)/='PHYSI')
        call ecoute('chm_reaphy')
     end do
     !
     ! Begin to read data
     !
     do while(words(1)/='ENDPH')
        call ecoute('chm_reaphy')

        if( words(1) == 'PROBL' ) then
           !
           ! Problem definition data
           !
           call ecoute('chm_reaphy')

           do while(words(1)/='ENDPR')

              !
              ! Chemistry model:
              !
              if( words(1) == 'MODEL' ) then
                 if( words(2) == 'FINIT' ) then
                    kfl_model_chm   = 3
                    if( exists('LAMIN') ) then
                       kfl_cotur_chm = 0
                    else if ( exists('LES  ') ) then
                       kfl_cotur_chm = -1
                    else if ( exists('RANS ') ) then
                       call runend('CHEMIC REAPHY: RANS model not available for FINITE-RATE chemistry') 
                    endif
                    write(momod(modul) % lun_outpu,*) 'Using finite rate model for calculations'
                    write(momod(modul) % lun_outpu,*)

                 else if (words(2) == 'CMC  ' ) then
                    kfl_model_chm   = 4
                    if( exists('LAMIN') ) then
                        kfl_cotur_chm = 0   ! In this case mixture fraction variance is set to 0 (for PDFs).
                    else if ( exists('LES  ') ) then
                        kfl_cotur_chm = -1
                    else if ( exists('RANS ') ) then
                        call runend('CHEMIC REAPHY: RANS model not available for FINITE-RATE chemistry')
                    endif
                    write(momod(modul) % lun_outpu,*) 'Using CMC model for calculations'
                    write(momod(modul) % lun_outpu,*)

                 else if( words(2) == 'FLAME' ) then
                    kfl_model_chm   = 1
                    !
                    ! Type of combustion problem
                    !
                    nclas_chm                   = 4  ! Non-premixed 
                    kfl_premix_chm = 0_ip            
                    if( .not.exists('NONPR') ) then  ! Premixed
                       nclas_chm                = 2
                       kfl_premix_chm = 1_ip            
                    end if
                    if (exists('NONAD')) kfl_heat_loss_chm = 1

                    !
                    ! Laminar or turbulent calculation
                    !
                    if( exists('LAMIN') ) then
                       kfl_cotur_chm = 0
                    elseif( exists('KEPSI')  ) then
                       kfl_cotur_chm = 1
                    elseif ( exists('KOMEG') ) then
                       kfl_cotur_chm = 2
                    elseif ( exists('LES  ') ) then
                       kfl_cotur_chm = -1
                    else
                       call runend('CHEMIC REAPHY: Turbulence model not valid for flamelet model') 
                    endif

                    !
                    ! Radiation model 
                    !
                    if( exists('RADIA') ) then
                       kfl_radia_chm = 1
                       if ( exists('WALLT') ) then
                          radwt_chm = getrea('WALLT',0.0_rp,'#Wall temperature radiation') 
                          kfl_radia_chm = 2
                       end if
                    end if

                    !
                    ! Unsteady flamelet progress variable model
                    !
                    if( exists('UFPV ') ) then
                       kfl_ufpv_chm = 3   ! default source is dot{Yc}

                       if( exists('OMEGA') ) kfl_ufpv_chm = 2   ! source is omega_Yc
                       if( exists('STEAD') ) kfl_ufpv_chm = 1   ! steady flamelet model
                    end if

                 else if( words(2) == 'SPRAY' ) then
                    kfl_premix_chm = 1_ip
                    kfl_spray_chm  = 1
                    kfl_model_chm  = 1
                    nclas_chm      = 4

                    if( exists('LAMIN') ) then
                       kfl_cotur_chm = 0
                    else if( exists('LES  ') ) then
                       kfl_cotur_chm = -1
                    endif

                 end if

              !
              ! Z variance model for Flamelet models
              !
              else if( words(1) == 'ZVARI' ) then 

                 !
                 ! Transport variance
                 !
                 if( words(2) == 'ZV   ' ) then
                    kfl_varZ_chm = 1
                    if( exists('LEA  ') ) then
                       kfl_varZ_chm = -1
                       call runend('CHEMIC REAPHY: LEA model for Z not implemented') 
                    end if
                 !
                 ! Transport Z*Z
                 !
                 else if( words(2) == 'ZZ  ') then
                    kfl_varZ_chm = 2 
                    if( exists('LEA  ') ) then
                       kfl_varZ_chm = -2
                       call runend('CHEMIC REAPHY: LEA model for Z not implemented') 
                    end if
                 !
                 ! Variance is 0 
                 !
                 else if( words(2) == 'OFF  ') then
                    kfl_varZ_chm = 0
                 end if

              !
              ! Yc variance model for Flamelet models
              !
              else if( words(1) == 'YCVAR' ) then               ! Yc variance model for flamelet combustion model
                 !
                 ! Transport variance
                 !
                 if( words(2) == 'YCV  ' ) then
                    kfl_varYc_chm = 1
                 !
                 ! Transport Yc*Yc
                 !
                 else if( words(2) == 'YCYC ') then
                    kfl_varYc_chm = 2 
                 !
                 ! Variance is 0 
                 !
                 else if( words(2) == 'OFF  ') then
                    kfl_varYc_chm = 0
                 end if

              !
              ! LOOKUP LOCATION 
              !
              else if( words(1) == 'LOOKU' ) then 

                 !
                 ! Transport variance
                 !
                 if( words(2) == 'GAUSS' ) then
                    kfl_lookg_chm = 1
                 elseif ( words(2) == 'NODES' ) then 
                    kfl_lookg_chm = 0
                 endif

              else if( words(1) == 'SPRAY' ) then               ! Problem name 
                 if( words(2) == 'EULER' ) then
                    kfl_spray_chm = 1
                    if (kfl_premix_chm == 0) nclas_chm = nclas_chm + 2

                 else if( words(2) == 'LEVEL') then
                    kfl_spray_chm = 2
                    !call runend('CHEMIC REAPHY: Level set model for ELSA not available yet, be patient')

                 else if( words(2) == 'ELSA ') then
                    kfl_spray_chm = 3
                    call runend('CHEMIC REAPHY: Full ELSA model not available yet, be patient')
                 end if

              else if( words(1) == 'TEMPO' ) then               ! Temporal evolution
                 if( words(2) == 'ON   ' ) then
                    kfl_timei_chm = 1
                 else
                    kfl_timei_chm = 0
                 end if

              else if( words(1) == 'CONVE' ) then               ! Convective term
                 if( words(2) == 'ON   ' ) then
                    kfl_advec_chm = -2
                    if(words(3)=='VELOC' ) then
                       kfl_advec_chm = -2 
                    end if
                 else if( words(2) == 'OFF  ' ) then
                    kfl_advec_chm = 0
                 end if

              else if( words(1) == 'DIFFU' ) then               ! Diffusion
                 if( words(2) == 'ON   ' ) then
                    kfl_diffu_chm = 1
                 else
                    kfl_diffu_chm = 0
                 end if

              else if( words(1) == 'HTRAN' ) then               ! Fluxes for enthalpy transport
                 if( words(2) == 'DETAI' ) then                 ! Detailed flux
                    kfl_htran = 1
                 end if


              else if(words(1).eq.'TURBU') then
                 if(words(2)=='LESMO') then
                    kfl_cotur_chm = -1_ip
                 else if(words(2)=='RANSM' .or. words(2)=='FROMT') then
                    kfl_cotur_chm = 1_ip
                 end if


              ! 
              ! Fuel and Oxidiser definition (used to calculate mixture Fraction) 
              ! 

              else if( words(1) == 'ZSCAL') then
                 bf_fuel_chm  = getrea('FUEL ',-1.0_rp,'#Fuel - Chemic')
                 bo_oxy_chm   = getrea('OXY  ',-1.0_rp,'#Fuel - Chemic')
                 kfl_z_chm    = -1_ip 

                 ! 
                 ! PFA Reduction - Read Non Key Species 
                 ! 
              else if( words(1) == 'REDUC') then 

                 if(words(2)=='PFA  ') then
                    ! 
                    ! Use equilibrium condition for non-Key species
                    !
                    if (exists('EQUIL')) then  
                       kfl_pfa_chm = 0
                       call ecoute('chm_reaphy')
                       backspace(momod(modul) % lun_pdata)
                       nbuff = 1000
                       allocate(character(len=nbuff) :: buffer)
                       read(momod(modul) % lun_pdata,'(a)', advance='no', size=nbuff, iostat=ioerror) buffer
                       if ((ioerror == iostat_eor).or.(ioerror == iostat_end)) ioerror = 0
                       if (ioerror /= 0 ) call runend('CHEMIC REAPHY: Cannot read Reduced Species.')
                       buffer = adjustl(buffer) ! trim leading whitespace
                       buffer = trim(buffer)    ! trim trailing whitespace
                       Red_spec  = buffer
                       nsize_red = len(Red_spec)
                       deallocate(buffer)

                    elseif (exists('STATI')) then
                       kfl_pfa_chm   = 1
                       kfl_key_chm   = getint('KEYSP',0_ip,'#Species starting fields')
                       dac_crit_chm  = getrea('CRITI',0.0_rp,'#Species starting fields')
                       call ecoute('chm_reaphy')
                       backspace(momod(modul) % lun_pdata)
                       nbuff = 1000
                       allocate(character(len=nbuff) :: buffer)
                       read(momod(modul) % lun_pdata,'(a)', advance='no', size=nbuff, iostat=ioerror) buffer
                       if ((ioerror == iostat_eor).or.(ioerror == iostat_end)) ioerror = 0
                       if (ioerror /= 0 ) call runend('CHEMIC REAPHY: Cannot read Reduced Species.')
                       buffer = adjustl(buffer) ! trim leading whitespace
                       buffer = trim(buffer)    ! trim trailing whitespace
                       Red_spec  = buffer
                       nsize_red = len(Red_spec)
                       deallocate(buffer)

                       ! 
                       ! Dynamic Adaptive Chemistry (DAC)
                       !  ("Red_Spec" holds key species)
                       !
                    else if(exists('DAC  ')) then
                       kfl_pfa_chm = 2
                       kfl_key_chm   = getint('KEYSP',0_ip,'#Species starting fields')
                       dac_crit_chm  = getrea('CRITI',0.0_rp,'#Species starting fields')
                       kfl_freq_chm  = getint('FREQU',0_ip,'#Steps between sucssive DAC reductions')
                       call ecoute('chm_reaphy')
                       backspace(momod(modul) % lun_pdata)
                       nbuff = 1000
                       allocate(character(len=nbuff) :: buffer)
                       read(momod(modul) % lun_pdata,'(a)', advance='no', size=nbuff, iostat=ioerror) buffer
                       if ((ioerror == iostat_eor).or.(ioerror == iostat_end)) ioerror = 0
                       if (ioerror /= 0 ) call runend('CHEMIC REAPHY: Cannot read Reduced Species.')
                       buffer = adjustl(buffer) ! trim leading whitespace
                       buffer = trim(buffer)    ! trim trailing whitespace
                       Red_spec  = buffer
                       nsize_red = len(Red_spec)
                       deallocate(buffer)

                       ! 
                       ! Correlated Dynamic Adaptive Chemistry (CODAC)
                       !  ("Red_Spec" holds key species)
                       !
                    else if(exists('CODAC')) then
                       kfl_pfa_chm = 3
                       kfl_key_chm  = getint('KEYSP',0_ip,'#Species starting fields')
                       dac_crit_chm = getrea('CRITI',0.0_rp,'#Species starting fields')
                       dac_cor_chm  = getrea('CORRE',0.0_rp,'#Species starting fields')
                       call ecoute('chm_reaphy')
                       backspace(momod(modul) % lun_pdata)
                       nbuff = 1000
                       allocate(character(len=nbuff) :: buffer)
                       read(momod(modul) % lun_pdata,'(a)', advance='no', size=nbuff, iostat=ioerror) buffer
                       if ((ioerror == iostat_eor).or.(ioerror == iostat_end)) ioerror = 0
                       if (ioerror /= 0 ) call runend('CHEMIC REAPHY: Cannot read Reduced Species.')
                       buffer = adjustl(buffer) ! trim leading whitespace
                       buffer = trim(buffer)    ! trim trailing whitespace
                       Red_spec  = buffer
                       nsize_red = len(Red_spec)
                       deallocate(buffer)

                       ! 
                       ! Tabulated Dynamic Adaptive Chemistry (T-DAC)
                       !  ("Red_Spec" holds defnition of progress variable)
                       !
                    else if(exists('TDAC ')) then
                       kfl_pfa_chm = 4
                       kfl_key_chm  = getint('KEYSP',0_ip,'#Species starting fields')
                       kfl_cont_chm = getint('CONTR',0_ip,'#Controling Varibales TDAC')
                       dac_crit_chm = getrea('CRITI',0.0_rp,'#Species starting fields')
                       if( exists('WRITE') ) kfl_tdac_write_chm = -10_ip
                       if( exists('READ ') ) kfl_tdac_write_chm = 10_ip
                       if( exists('FILL ') ) kfl_tdac_write_chm = 666_ip
                       call ecoute('chm_reaphy')
                       backspace(momod(modul) % lun_pdata)
                       nbuff = 1000
                       allocate(character(len=nbuff) :: buffer)
                       read(momod(modul) % lun_pdata,'(a)', advance='no', size=nbuff, iostat=ioerror) buffer
                       if ((ioerror == iostat_eor).or.(ioerror == iostat_end)) ioerror = 0
                       if (ioerror /= 0 ) call runend('CHEMIC REAPHY: Cannot read Reduced Species.')
                       buffer = adjustl(buffer) ! trim leading whitespace
                       buffer = trim(buffer)    ! trim trailing whitespace
                       Red_spec  = buffer
                       nsize_red = len(Red_spec)
                       deallocate(buffer)
                    end if

                 else if (exists('NONE ') .or. exists('OFF  ')) then
                    !
                    ! Complete Integration of the mechanism 
                    !
                    kfl_pfa_chm = -1_ip

                 else
                    call runend('CHEMIC REAPHY: Reduction activated, but no additional information provided')
                 end if
                 !
                 ! Reactions
                 !
              else if( words(1) == 'MECHA' ) then
                 !
                 ! Cantera
                 !
                 if( words(2) == 'CANTE' ) then

                    !
                    ! Read mechanism path from:   momod(modul) % fil_pdata   with ID: momod(modul) % lun_pdata 
                    !  
                    ii = 0 
                    do while(words(1)/='ENDME')
                       ii = ii + 1
                       !
                       ! Call ecoute to look for the END_REACT, and step back
                       !
                       call ecoute('chm_reaphy')
                       backspace(momod(modul) % lun_pdata)

                       !
                       ! Read next line 
                       !
                       nbuff = 100
                       allocate(character(len=nbuff) :: buffer)
                       read(momod(modul) % lun_pdata,'(a)', advance='no', size=nbuff, iostat=ioerror) buffer   

                       !
                       ! Check read, disregard errors about buffer being larger
                       ! than the line.
                       !
                       if ((ioerror == iostat_eor).or.(ioerror == iostat_end)) ioerror = 0
                       if (ioerror /= 0 ) call runend('CHEMIC REAPHY: Cannot read cantera mechanism path.')

                       !
                       ! Process line
                       !
                       buffer = adjustl(buffer) ! trim leading whitespace
                       buffer = trim(buffer)    ! trim trailing whitespace
                       jj     = 0_ip
                       jj     = index(buffer,'xml')
                       jj     = jj + 2_ip
                       if (mechanism_path=='') then
                          !
                          ! Take first line as path (other lines are ignored
                          ! for now)
                          !
                          mechanism_path  = buffer(1:jj)
                          nsize_mech_name = len(mechanism_path)
                       endif

                       deallocate(buffer)

                    end do
#ifdef CANTERA
                    gas_chm   = importPhase(mechanism_path)
                    nclas_chm = nSpecies(gas_chm)    ! number of species
                    nreac_chm = nReactions(gas_chm)  ! number of reactions
#endif
                 else
                    call runend('CHEMIC REAPHY: Mechanism selection only valid with FINITE RATE model')
                 endif

              end if

              call ecoute('chm_reaphy')
           end do

        else if( words(1) == 'PROPE' ) then
           !
           ! Properties
           !
           call chm_memphy(1_ip)

           call ecoute('chm_reaphy')
           !
           write(momod(modul) % lun_outpu,*)'-----------------'
           write(momod(modul) % lun_outpu,*)' FLOW PROPERTIES '
           write(momod(modul) % lun_outpu,*)'-----------------'
           write(momod(modul) % lun_outpu,*) ''
           !
           do while(words(1)/='ENDPR')

              !-------------------------------------------------------
              !
              ! Properties and constants
              !
              !-------------------------------------------------------

              if ( words(1) == 'GASCO' ) then              ! Universal gas constant
                 gasco = getrea('GASCO',8.3144621_rp,'#Gas constant')

              else if ( words(1) == 'SURFA' ) then              ! Surface tension
                 surf_tension_chm = getrea('SURFA',1.0_rp,'#Surface tension liquid')
                 !
                 ! Variable names for ELSA model
                 !
                 if (kfl_spray_chm /= 0) then
                    if (nclas_chm == 4_ip ) then
                       speci(1_ip)%name = 'NONE'
                       speci(2_ip)%name = 'NONE'
                       speci(3_ip)%name = 'PHI'
                       speci(4_ip)%name = 'SIGMA'
                    elseif (nclas_chm == 6_ip ) then
                       speci(5_ip)%name = 'PHI'
                       speci(6_ip)%name = 'SIGMA'
                    end if
                 end if

              else if ( words(1) == 'THERM' ) then              ! Thermodynamic pressure
                 prthe_chm = getrea('THERM',101325.0_rp,'#Thermodynamic pressure')
                 if (kfl_coupl(ID_CHEMIC,ID_NASTIN) == 0 ) prthe = prthe_chm

              else if( words(1) == 'TURBU' ) then
                 diffu_chm(1,1) = getrea('TURBU',0.9_rp,'#Turbulent Schmidt number')    ! Turbulent Schmidt number
                 write(momod(modul) % lun_outpu,*) 'TURBULENT SCHMIDT NUMBER:',diffu_chm(1,1)
                 write(momod(modul) % lun_outpu,*) ''

              !!DMM else if( words(1) == 'LAWDI' ) then                ! Diffusion law
              !!DMM    if( words(2) == 'CONST' ) then
              !!DMM       if( exists('SPECY') ) then
              !!DMM          ispec = getint('SPECY',1_ip,'#Specy number')                    
              !!DMM          if (exists('VALUE')) diffu_chm(1,ispec) = getrea('VALUE',1.0_rp,'#Diffusion constatnt')
              !!DMM          lawdi_chm(1,ispec) = 0
              !!DMM       else
              !!DMM          lawdi_chm =  0
              !!DMM          diffu_chm = getrea('VALUE',1.0_rp,'#Diffusion constatnt')
              !!DMM       end if
              !!DMM    else if( words(2) == 'PRAND') then
              !!DMM       lawdi_chm = 1
              !!DMM       if(exists('RANS ')) then
              !!DMM          if (exists('TURBU')) diffu_chm(1,1) = getrea('TURBU',0.9_rp,'#Turbulent Schmidt number')
              !!DMM          lawdi_chm = 5
              !!DMM          write(momod(modul) % lun_outpu,*) 'TURBULENT SCHMIDT NUMBER FOR RANS:',diffu_chm(1,1)
              !!DMM          write(momod(modul) % lun_outpu,*) ''
              !!DMM          if (kfl_coupl(ID_CHEMIC,ID_TURBUL) == 0 ) then
              !!DMM             write(momod(modul) % lun_outpu,*) 'CHEMIC: TURBUL is not activated'
              !!DMM             !                            call runend('CHEMIC: For RANS simulation TURBUL must be activated') 
              !!DMM          endif
              !!DMM       elseif ( exists('LES  ') ) then
              !!DMM          if (exists('TURBU')) diffu_chm(1,1) = getrea('TURBU',0.9_rp,'#Turbulent Schmidt number')
              !!DMM          lawdi_chm = 5
              !!DMM          write(momod(modul) % lun_outpu,*) 'TURBULENT SCHMIDT NUMBER FOR LES:',diffu_chm(1,1)
              !!DMM          write(momod(modul) % lun_outpu,*) ''
              !!DMM       endif

              !!DMM    else if( words(2) == 'UNIFO') then
              !!DMM       lawdi_chm = 3
              !!DMM       diffu_chm(1,1) = getrea('PRAND',0.71_rp,'#Uniform Prandtl number')
              !!DMM       diffu_chm(2,1) = getrea('LEWIS',1.0_rp,'#Uniform Lewis number')                    
              !!DMM    else if( words(2) == 'WATER') then
              !!DMM       lawdi_chm = 4
              !!DMM    end if

              else if( words(1) == 'KINET' ) then

                 if ( words(2) == 'CANTE' ) then 
#ifdef CANTERA
                    write(momod(modul) % lun_outpu,*)'Reading mechanism ',mechanism_path
                    write(momod(modul) % lun_outpu,*)''

                    write(momod(modul) % lun_outpu,*)'# Species=',nclas_chm
                    write(momod(modul) % lun_outpu,*)''

                    write(momod(modul) % lun_outpu,*)'# Removed Species=',Red_spec
                    write(momod(modul) % lun_outpu,*)''

                    do iclas=1,nclas_chm
                       call getSpeciesName(gas_chm,iclas,speci(iclas)%name)
                       write(momod(modul) % lun_outpu,*)'Species ',iclas,'=',speci(iclas)%name
                    end do
                    write(momod(modul) % lun_outpu,*)''
                    write(momod(modul) % lun_outpu,*)'Chemical reactions'
                    write(momod(modul) % lun_outpu,*)''
                    write(momod(modul) % lun_outpu,*)'# Reactions=',nreac_chm
                    write(momod(modul) % lun_outpu,*)''

                    do ireac = 1,nreac_chm
                       call getReactionString(gas_chm, ireac,eq)
                       write(momod(modul) % lun_outpu,*)eq
                    end do
                    write(momod(modul) % lun_outpu,*)''
#endif
                 else if ( words(2) == 'USER' ) then 

                 end if

              else if( words(1) == 'TRANS' ) then

                 if( words(2) == 'CANTE' ) then
                    kfl_transport_chm = 1
                 elseif ( words(2) == 'USER ' ) then
                    kfl_transport_chm = 2
                 end if

              !
              ! Specify Lewis numbers
              !
              else if(words(1) == 'LEWIS' ) then  ! User specified Lewis numbers
                 call ecoute('chm_reaphy')
                 do while(words(1) /= 'ENDLE')
                    iclas = 1
                    found = .false.
                    do while( (iclas <= nclas_chm) .and. (.not. found))
                       ! OJO: ¿Y SI LA ESPECIE TIENE MÁS DE 5 LETRAS?
                       if (words(1) == speci(iclas)%name) then
                          Le_k(iclas) = param(1)
                          found = .true.
                       end if
                       iclas = iclas + 1
                    end do
                    if( .not. found ) call runend('CHEMIC REAPHY: Not all the species exist in the mechanism')
                    call ecoute('chm_reaphy')
                 end do

                 write(momod(modul) % lun_outpu,*) 'Lewis number equal to 1 except for species'
                 do iclas = 1, nclas_chm
                    if (Le_k(iclas) /= 1.0_rp) then
                       write(momod(modul) % lun_outpu,*) speci(iclas)%name, ': ', Le_k(iclas)
                    end if
                 end do

              !
              ! Reduction Table
              !
              else if( words(1) == 'REDUC' ) then

                 !
                 ! Flamelet model: read thermochemical database 
                 ! 
                 if( words(2) == 'TABLE' ) then
                    if (kfl_pfa_chm /= 4) call runend('CHEMIC REAPHY: Reduction table only for TDAC')
                    !
                    ! TDAC - Dynamic tabulation scheme 
                    !
                    if ( kfl_tdac_write_chm > 0_ip ) then
                       call tab_load_file(table_coords, table_tab, dont_read_data=.true.)
                       print *, "Not Reading Table Data"
                    else 
                       call tab_load_file(table_coords, table_tab, dont_read_data=.false.)
                       print *, "Reading Table Data"
                    end if
                    table_fw % main_table => table_tab


                    do while(words(1)/='ENDRE')
                       call ecoute('chm_reaphy')
                    end do
                    !
                    ! Print table dimensions
                    !
                    write(momod(modul) % lun_outpu,*)''
                    write(momod(modul) % lun_outpu,*)'---------------------------------------------'
                    write(momod(modul) % lun_outpu,*)'TDAC - REDUCTION:'
                    write(momod(modul) % lun_outpu,*)'---------------------------------------------'
                    write(momod(modul) % lun_outpu,*)''
                    write(momod(modul) % lun_outpu,*)'TABLE DIMENSIONS:'
                    do ii = 1,table_fw % main_table%ndim
                       write(momod(modul) % lun_outpu,'(5X,A5,1X,I5)') table_fw % main_table % coords(ii) % name , table_fw % main_table % coords(ii) % leng
                    enddo
                 end if

              ! Reactions
              !
              else if( words(1) == 'TABFR' ) then
                 !
                 ! Table framework
                 !
                 kfl_tab_fw_chm = getint('TABFR',1_ip,'#Index of table framework')
                 table_fw => lookup_fw( kfl_tab_fw_chm )

              else if( words(1) == 'POSTF' ) then
                 !
                 ! Postprocessing table framework
                 !
                 kfl_post_fw_chm = getint('POSTF',1_ip,'#Index of postprocessing table framework')
                 posttable_fw => lookup_fw( kfl_post_fw_chm  )

              else if( words(1) == 'REACT' ) then

                 !
                 ! Flamelet model: read thermochemical database LEGACY 
                 ! 
                 if( words(2) == 'TABLE' ) then
                    if (kfl_model_chm == 1 .and. kfl_pfa_chm == 4) call runend('CHEMIC REAPHY: Cannot use flamelets and finite rate together DUH !!')
                    if (kfl_model_chm /= 1) call runend('CHEMIC REAPHY: Table should only be used with the flamelet combustion model')

                    !
                    ! Control variable names for Flamlet model
                    !
                    speci(1_ip)%name = 'CMEAN'
                    speci(2_ip)%name = 'CVAR'
                    if (nclas_chm > 2_ip) then
                       speci(3_ip)%name = 'ZMEAN'
                       speci(4_ip)%name = 'ZVAR'
                    endif

                    !
                    ! Control variable names for UFPV model
                    !
                    if (kfl_ufpv_chm > 0) then
                       speci(1_ip)%name = 'CMEAN'
                       speci(2_ip)%name = 'CHIST'
                       speci(3_ip)%name = 'ZMEAN'
                       speci(4_ip)%name = 'ZVAR'

                       !
                       ! NO VARIANCE FOR Yc
                       !
                       kfl_varYc_chm = 0_ip 
                    endif

                    !
                    ! Flamelet model, read tabulated data
                    !
                    if (exists('ORDER')) then
                       order_loc = getint('ORDER',1_ip,'#Order of interpolation.')
                    endif
                    call tab_load_file(table_coords, table_tab, order=order_loc)
                    kfl_tab_fw_chm = 0 
                    allocate(table_fw)
                    call tab_init_fw(table_fw)
                    table_fw % main_table => table_tab
                    call fw_allocate(1_ip,table_fw)

                    !
                    ! Decide scaling method
                    !
                    table_fw % kfl_scale = 0_ip
                    do ii = 1, table_fw % main_table % ndim
                        if (table_fw % main_table % coords(ii) % leng > 1) then 
                           select case (table_fw % main_table % coords(ii) % name)
                           case ('CMEAN','C    ')
                              if (kfl_premix_chm /= 0) then            
                                 table_fw % kfl_scale(ii) = 0
                              else
                                 table_fw % kfl_scale(ii) = 1
                              endif
                           case ('CVAR ')
                              if (abs(kfl_varYc_chm) == 1) then            
                                 table_fw % kfl_scale(ii) = 100 + ii - 1
                              elseif (abs(kfl_varYc_chm) == 2) then            
                                 table_fw % kfl_scale(ii) = -1*(100 + ii - 1)
                              else
                                 table_fw % kfl_scale(ii) = -1
                              endif
                           case ('CHIST')
                              table_fw % kfl_scale(ii) = 1
                           case ('ZMEAN','Z    ')
                              table_fw % kfl_scale(ii) = 0
                           case ('ZVAR ')
                              if (abs(kfl_varZ_chm) == 1) then            
                                 table_fw % kfl_scale(ii) = 100 + ii - 1
                              elseif (abs(kfl_varZ_chm) == 2) then            
                                 table_fw % kfl_scale(ii) = -1*(100 + ii - 1)
                              else
                                 table_fw % kfl_scale(ii) = -1
                              endif
                           case ('IMEAN','I    ')
                              if (kfl_heat_loss_chm /= 0) then
                                  if (kfl_premix_chm /= 0) then            
                                     table_fw % kfl_scale(ii) = 0
                                  else
                                     table_fw % kfl_scale(ii) = 1
                                  endif
                              else
                                 table_fw % kfl_scale(ii) = -1
                              endif
                           case default
                              call runend('chm_reaphy: cannot understand coordinate name: ' // table_fw % main_table % coords(ii) % name)
                           end select
                        else
                           table_fw % kfl_scale(ii) = -1
                        endif
                    enddo
                    !
                    ! Calculate scaling order
                    !
                    call fw_allocate(2_ip,table_fw,0_ip)

                    do while(words(1)/='ENDRE')
                       call ecoute('chm_reaphy')
                    end do

                    !
                    ! Print table dimensions
                    !
                    write(momod(modul) % lun_outpu,*)''
                    write(momod(modul) % lun_outpu,*)'---------------------------------------------'
                    write(momod(modul) % lun_outpu,*)'FLAMELET-MODEL: THERMOCHEMICAL DATABASE '
                    write(momod(modul) % lun_outpu,*)'---------------------------------------------'
                    write(momod(modul) % lun_outpu,*)''
                    if (kfl_premix_chm /= 0) then 
                       write(momod(modul) % lun_outpu,*)'PREMIXED COMBUSTION'
                    else 
                       write(momod(modul) % lun_outpu,*)'NON-PREMIXED COMBUSTION'
                    end if
                    write(momod(modul) % lun_outpu,*)'TABLE DIMENSIONS:'

                    do ii = 1,table_fw % main_table % ndim
                       write(momod(modul) % lun_outpu,'(5X,A5,1X,I5)') table_fw % main_table % coords(ii) % name , table_fw % main_table % coords(ii) % leng
                    enddo

                 end if

                 if (kfl_model_chm /= 3 ) then 
                    do while( words(1) /= 'ENDRE' )
                       call ecoute('chm_reaphy')
                    end do
                 end if

              !
              ! Read Yc scaling table
              !

              else if (words(1)=='MASSF' ) then
                 call tab_load_file(yc_scaling_coords, yc_scaling_tab, order=order_loc)
                 do ii = 1, table_fw % main_table % ndim
                     if  (table_fw % main_table % coords(ii) % name == 'CMEAN' .or. table_fw % main_table % coords(ii) % name == 'C    ') then
                        table_fw % scaling(ii) % tab => yc_scaling_tab
                        call fw_allocate(2_ip,table_fw,ii)
                     endif
                 enddo
                 do while( words(1) /= 'ENDMA' )
                    call ecoute('chm_reaphy')
                 end do
              end if

              !
              ! Read scalar dissipation shape table
              !
              if (words(1)=='CHISC' .and. kfl_ufpv_chm > 0 ) then
                 call tab_load_file(chi_shape_coords, chi_shape_tab, order=order_loc)
                 do ii = 1, table_fw % main_table % ndim
                     if  (table_fw % main_table % coords(ii) % name == 'CHIST' ) then
                        table_fw % scaling(ii) % tab => chi_shape_tab
                        call fw_allocate(2_ip,table_fw,ii)
                     endif
                 enddo
                 do while( words(1) /= 'ENDCH' )
                    call ecoute('chm_reaphy')
                 end do
              end if

              !
              ! Read enthalpy loss table
              !
              if (words(1)=='HEATL' .and. kfl_heat_loss_chm > 0 ) then
                 call tab_load_file(h_scaling_coords, h_scaling_tab, order=order_loc)
                 do ii = 1, table_fw % main_table % ndim
                     if  (table_fw % main_table % coords(ii) % name == 'IMEAN' .or. table_fw % main_table % coords(ii) % name == 'I    ') then
                        table_fw % scaling(ii) % tab => h_scaling_tab
                        call fw_allocate(2_ip,table_fw,ii)
                     endif
                 enddo
                 do while( words(1) /= 'ENDHE' )
                    call ecoute('chm_reaphy')
                 end do
              end if

              !
              ! Read post-processing table
              !
              if (words(1)=='POSTT') then
                 call tab_load_file(posttab_coords, posttab_tab, order=order_loc)
                 kfl_post_fw_chm = 0 

                 allocate(posttable_fw)
                 call tab_init_fw(posttable_fw)
                 posttable_fw % main_table => posttab_tab
                 call fw_allocate(1_ip,posttable_fw)

                 posttable_fw % kfl_scale   = table_fw % kfl_scale 
                 posttable_fw % scaling     = table_fw % scaling
                 posttable_fw % scale_order = table_fw % scale_order
                 do while( words(1) /= 'ENDPO' )
                    call ecoute('chm_reaphy')
                 end do

                 !
                 ! Select names:
                 !
                 do ipostvar=1,posttable_fw % main_table % nvar
                    give_new_name = .false.
                    !
                    ! Check if empty name
                    !
                    if (posttable_fw % main_table % varname(ipostvar) == '') then
                        give_new_name = .true.
                    endif 

                    !
                    ! Check if name is repeated
                    !
                    do jpostvar=1,ipostvar-1
                       if (posttable_fw % main_table % varname(ipostvar) == posttable_fw % main_table % varname(jpostvar)) then
                           give_new_name = .true.
                       endif
                    enddo
                    
                    !
                    ! If necessary give new name
                    !
                    if (give_new_name) then
                        if (ipostvar < 10) then
                            posttable_fw % main_table % varname(ipostvar) = 'POS0'//trim(intost(ipostvar))
                        elseif (ipostvar < 100) then
                            posttable_fw % main_table % varname(ipostvar) = 'POS'//trim(intost(ipostvar))
                        else
                            posttable_fw % main_table % varname(ipostvar) = 'PO'//trim(intost(ipostvar))
                        endif
                    endif 
                 enddo 
              end if

              call ecoute('chm_reaphy')
           end do
        
        !
        ! Reading soot model variables
        !
        else if( words(1) == 'SOOTM' ) then

           call ecoute('chm_reaphy')

           !
           !
           !
           write(momod(modul) % lun_outpu,*)'---------------------'
           write(momod(modul) % lun_outpu,*)' SOOT MODEL VARIBLES '
           write(momod(modul) % lun_outpu,*)'---------------------'
           write(momod(modul) % lun_outpu,*) ''

           !
           ! Soot model reading
           !
           do while(words(1)/='ENDSO')

              if ( words(1) == 'MODEL' ) then
                 if( words(2) == 'SECTI' ) then
                    kfl_soot_chm = 1
                    write(momod(modul) % lun_outpu,*)''
                    write(momod(modul) % lun_outpu,*)'Soot sectional model activated, model =',kfl_soot_chm
                    write(momod(modul) % lun_outpu,*)''
                 else
                    call runend('CHM_REAPHY: Only soot sectional method available, be patient!')
                 end if
              else if ( words(1) == 'SECTI' ) then
                 nsect_ssm = getint('SECTI',1_ip,'#Number of sections soot model')
                 write(momod(modul) % lun_outpu,*)'# Number of sections = ',nsect_ssm
                 write(momod(modul) % lun_outpu,*)''

              else if ( words(1) == 'DENSI' ) then
                 RhoC_ssm = getrea('DENSI',1.0_rp,'#Soot particles density')
                 write(momod(modul) % lun_outpu,*)'# Soot particle density = ',RhoC_ssm,'[kg/m3]'
                 write(momod(modul) % lun_outpu,*)''

              else if ( words(1) == 'VMAXS' ) then
                 Vmax_ssm = getrea('VMAXS',1.0_rp,'#Maximum volume soot particles')
                 write(momod(modul) % lun_outpu,*)'# Maximum volume soot particles = ',Vmax_ssm,'[m3]'
                 write(momod(modul) % lun_outpu,*)''

              else if ( words(1) == 'RADPA' ) then
                 RadSoot_ssm = getrea('RADPA',1.0_rp,'#Radiation parameter')
                 write(momod(modul) % lun_outpu,*)'# Radiation parameter = ',RadSoot_ssm
                 write(momod(modul) % lun_outpu,*)''

              else if ( words(1) == 'VALID' ) then
                 if( words(2) == 'CHEM1' ) then
                    kfl_soot_chm = -1
                    write(momod(modul) % lun_outpu,*)'# Validation test with Chem1d '
                    write(momod(modul) % lun_outpu,*)''
                 endif

              else if ( words(1) == 'GASCO' ) then
                 if( words(2) == 'ON' ) then   
                    gasCoupling_ssm = 1 
                 elseif ( words(2) == 'OFF' ) then
                    gasCoupling_ssm = 0 
                 else
                    call runend('CHM_REAPHY: Coupling soot with gas phase not specified.')
                 endif

                 write(momod(modul) % lun_outpu,*)'# Coupling with gas phase = ',gasCoupling_ssm
                 write(momod(modul) % lun_outpu,*)''

              end if

              call ecoute('chm_reaphy')
           end do

        else if ( words(1) == 'DROPL' ) then
           !
           ! Eulerian droplets identification parameters
           !
           kfl_droplet_id_chm = 1_ip 
           call ecoute('chm_reaphy')
           !
           write(momod(modul) % lun_outpu,*) ''
           write(momod(modul) % lun_outpu,*)'---------------------------------'
           write(momod(modul) % lun_outpu,*)' EULERIAN DROPLET IDENTIFICATION '
           write(momod(modul) % lun_outpu,*)'---------------------------------'
           write(momod(modul) % lun_outpu,*) ''
           !
           do while(words(1)/='ENDDR')

              if ( words(1) == 'LEVEL' ) then   
                 levelSet_threshold_chm = getrea('LEVEL',0.0_rp,'#Droplet identification parameter')
                 write(momod(modul) % lun_outpu,*) 'LEVEL SET THRESHOLD:',levelSet_threshold_chm
              else if ( words(1) == 'COMPA' ) then
                 droplet_compactness_limit_chm = getrea('COMPA',1.0_rp,'#Droplet identification parameter')
                 write(momod(modul) % lun_outpu,*) 'COMPACTNESS LIMIT:',droplet_compactness_limit_chm
              else if ( words(1) == 'MAXDI' ) then
                 droplet_max_diameter_chm = getrea('MAXDI',0.0_rp,'#Droplet identification parameter')
                 write(momod(modul) % lun_outpu,*) 'MAX DIAMETER:',droplet_max_diameter_chm
              else if ( words(1) == 'HFACT' ) then 
                 droplet_h_factor_chm = getrea('HFACT',0.0_rp,'#Droplet identification parameter')
                 write(momod(modul) % lun_outpu,*) 'MESH SIZE FACTOR:',droplet_h_factor_chm
              end if
              call ecoute('chm_reaphy')
           end do
           write(momod(modul) % lun_outpu,*) ''
        !
        ! Read index of starting fields
        !
        else if (words(1)=='FIELD') then
           kfl_field_chm(1) = getint('SPECI',0_ip,'#Species starting fields')
           kfl_field_chm(2) = getint('TEMPE',0_ip,'#Temperature field')
           kfl_field_chm(3) = getint('ALPHA',0_ip,'#Alpha field for CMC model')

        !
        ! Read species names 
        !
        else if (words(1)=='SPECI') then
           kfl_spec_name_chm = getint('SPECI',0_ip,'#Number of species')
           allocate(Field_ind_chm(kfl_spec_name_chm))
           call ecoute('chm_reaphy')
           backspace(momod(modul) % lun_pdata)
           nbuff = 1000
           allocate(character(len=nbuff) :: buffer)
           read(momod(modul) % lun_pdata,'(a)', advance='no', size=nbuff, iostat=ioerror) buffer
           if ((ioerror == iostat_eor).or.(ioerror == iostat_end)) ioerror = 0
           if (ioerror /= 0 ) call runend('CHEMIC REAPHY: Cannot read Reduced Species.')
           buffer = adjustl(buffer) ! trim leading whitespace
           buffer = trim(buffer)    ! trim trailing whitespace
           spec_name  = buffer
           nsize_spec_name = len(spec_name)
#ifdef CANTERA
           call cantera_initialization(1_ip,mechanism_path,nsize_mech_name)
           call cantera_trim(spec_name, nsize_spec_name, Field_ind_chm)
#endif
           deallocate(buffer)


        else if( words(1) == 'CMCMO' ) then
           !
           ! Data for CMC model
           !
           call ecoute('chm_reaphy')

           CMC_model: do while(words(1)/='ENDCM')

              !
              ! CMC model: read mixture fraction vector
              !
              ! Read path and length for mixture fraction vector
              if( words(1) == 'MFVEC') then
                 do icoun = 1,2
                    call ecoute('chm_reaphy' )
                    if ( words(1) == 'MFPAT') then   ! Read path to the file
                       !
                       ! Read next line 
                       !
                       nbuff = 1000
                       allocate(character(len=nbuff) :: buffer)
                       read(momod(modul) % lun_pdata,'(a)', advance='no', size=nbuff, iostat=ioerror) buffer

                       !
                       ! Check read, disregard errors about buffer being larger
                       ! than the line.
                       !
                       if ((ioerror == iostat_eor) .or. (ioerror == iostat_end)) ioerror = 0
                       if (ioerror /= 0 ) call runend('CHEMIC REAPHY: Cannot read mixture fraction vector path.')

                       !
                       ! Process line
                       !
                       buffer = adjustl(buffer) ! trim leading whitespace
                       buffer = trim(buffer)    ! trim trailing whitespace
                       nsize_mf_vec_chm = len(buffer)
                       mf_vec_path_CMC_chm = buffer(1:nsize_mf_vec_chm)
                       deallocate(buffer)
                    else if ( words(1) == 'MFSIZ') then  ! Read size of the vector
                       nZ_CMC_chm = param(1)
                    end if
                 end do

                 if (nsize_mf_vec_chm==0_ip .or. nZ_CMC_chm==1_ip) call runend('CHEMIC REAPHY: bad entries for mixture fraction vector definition')


                 ! Read mixture fraction vector

                 write(momod(modul) % lun_outpu,*) 'Reading mixture fraction vector from ', mf_vec_path_CMC_chm
                 call chm_memphy(3_ip)

                 open(unit=1, file=mf_vec_path_CMC_chm, status='old', action='read', iostat=state_f)
                 if (state_f==0) then
                    do line_f = 1,nZ_CMC_chm
                       read(unit=1, fmt=*) Z_CMC_chm(line_f)
                    end do
                 else
                    call runend('CHEMIC REAPHY: Cannot read mixture fraction vector file')
                 end if
                 close(unit=1, iostat=state_f, status='keep')
                 Zs_CMC_chm = Z_CMC_chm(nZ_CMC_chm)

                 ! Do some checks about the vector consistency (Zs<1, monotonically increasing, etc.)
                 if (Z_CMC_chm(1) /= 0.0_rp)          call runend('CHEMIC REAPHY: mixture fraction value for oxidizer different to 0')
                 if (Z_CMC_chm(nZ_CMC_chm) > 1.0_rp)  call runend('CHEMIC REAPHY: mixture fraction value for fuel not valid')
                 do line_f = 1,nZ_CMC_chm-1
                    if (Z_CMC_chm(line_f+1) <= Z_CMC_chm(line_f))  call runend('CHEMIC REAPHY: mixture fraction vector not strictly monotonically increasing')
                 end do

                 diff_Z_CMC_chm = Z_CMC_chm(2:nZ_CMC_chm) - Z_CMC_chm(1:nZ_CMC_chm-1)

                 write(momod(modul) % lun_outpu,*) 'Values for mixture fraction slices:'
                 do line_f = 1,nZ_CMC_chm
                    write(momod(modul) % lun_outpu,*) line_f, ': ',  Z_CMC_chm(line_f)
                 end do
                 write(momod(modul) % lun_outpu,*)
                 write(momod(modul) % lun_outpu,*) 'Maximum mixture fraction is ', Zs_CMC_chm
                 write(momod(modul) % lun_outpu,*)
                 write(momod(modul) % lun_outpu,*)


              !
              ! Number of mixture fractions and segregation factors for AMC model
              !
              else if( words(1) == 'AMCMF' ) then
                 nZ_AMC_CMC_chm = param(1)

              else if( words(1) == 'AMCSS' ) then
                 nS_AMC_CMC_chm = param(1)

              else if( words(1) == 'AMCSM' ) then
                 Smax_AMC_CMC_chm = param(1)

              else if( words(1) == 'STHRE' ) then
                 S_threshold = param(1)

              !
              ! CMC model: define boundary conditions for species and enthalpy 
              ! at mixt. frac. = 0 and Zs. This information, if appears, is only 
              ! used when starting from scratch.
              !
              else if( words(1) == 'ZBCOX' .or. words(1) == 'ZBCFU' ) then
                 if ( kfl_weigh_in_eq_CMC_chm == 0_ip ) then
                    call chm_memphy(4_ip)  ! It is the first time we enter in this conditional
                    kfl_weigh_in_eq_CMC_chm = 1_ip
                 end if
                 if ( words(1) == 'ZBCOX' ) then
                    index_zbc = 1
                 else
                    index_zbc = 2
                 end if

                 call ecoute('chm_reaphy')
                 if (words(1) == 'TEMPE') then
                    T_bc_CMC_chm(index_zbc) = param(1)
                 else
                    call runend('CHEMIC REAPHY: Temperature not provided at the mixture fraction boundary conditions (CMC model)')
                 end if

                 call ecoute('chm_reaphy')
                 do while(words(1) /= 'ENDZB')
                    iclas = 1
                    found = .false.
                    do while( (iclas <= nclas_chm) .and. (.not. found))
                       ! OJO: ¿Y SI LA ESPECIE TIENE MÁS DE 5 LETRAS?
                       if (words(1) == speci(iclas)%name) then
                          react_scal_bc_CMC_chm(iclas+1,index_zbc) = param(1)
                          found = .true.
                       end if
                       iclas = iclas + 1
                    end do
                    if( .not. found ) call runend('CHEMIC REAPHY: Not all the species exist in the mechanism')
                    call ecoute('chm_reaphy')
                 end do

              else if( words(1) == 'SOLEN' ) then
                 if ( words(2) == 'NO   ') then
                    kfl_solve_enth_CMC_chm = 0_ip
                 else
                    kfl_solve_enth_CMC_chm = 1_ip
                 end if

              else if( words(1) == 'SPLIT' ) then
                 if ( words(2) == 'NO   ') then
                    kfl_split_CFD_CMC = 0_ip
                 else
                    kfl_split_CFD_CMC = 1_ip
                 end if

              end if

              call ecoute('chm_reaphy')

           end do CMC_model

        end if

     end do

     ! Information about enthalpy flux model
     if ( kfl_model_chm  == 3) then
        write(momod(modul) % lun_outpu,*)
        if ( kfl_htran == 1 ) then
           write(momod(modul) % lun_outpu,*) 'Detailed computation for enthalpy flux'
        else
           write(momod(modul) % lun_outpu,*) 'Simplified computation for enthalpy flux'
        end if
        write(momod(modul) % lun_outpu,*)
     end if

     !  
     ! Definition of unknowns for chemic with soot model  
     !
     !   Y*_k = Y_k + Y_s 
     !
     if (kfl_soot_chm /= 0) then
        nspec_chm = nclas_chm
        nclas_chm = nspec_chm + nsect_ssm

        nspec_ssm = nspec_chm
        nclas_ssm = nclas_chm
     else
        nspec_chm = nclas_chm
     end if

     ! Information and actions for CMC model
     if( kfl_model_chm == 4 ) then
        if (nZ_CMC_chm == 3) then
           call runend('CHEMIC REAPHY: mixture fraction vector not provided for CMC model')
        end if

        if (kfl_lookg_chm == 0) then
           call runend('CHEMIC REAPHY: CMC model not compatible with calculation at nodes')
        end if

        if(kfl_solve_enth_CMC_chm == 0) then
           nvar_CMC_chm = nclas_chm ! CMC only solves species
        else
           nvar_CMC_chm = nclas_chm + 1  ! CMC solves species + enthalpy
        end if
      
        ! AMC model for scalar dissipation rate
        call chm_memphy(5)
        call chm_AMC_generate_Z_S_vectors_CMC

        write(momod(modul) % lun_outpu,*)
        write(momod(modul) % lun_outpu,*) 'Mixture fraction (', nZ_AMC_CMC_chm, ' entries) for AMC model:'
        do line_f = 1,nZ_AMC_CMC_chm
           write(momod(modul) % lun_outpu,*) line_f, ': ', Z_AMC_CMC_chm(line_f)
        end do
        write(momod(modul) % lun_outpu,*)

        write(momod(modul) % lun_outpu,*)
        write(momod(modul) % lun_outpu,*) 'Segregation factor (', nS_AMC_CMC_chm, ' entries) for AMC model:'
        do line_f = 1,nS_AMC_CMC_chm
           write(momod(modul) % lun_outpu,*) line_f, ': ', S_AMC_CMC_chm(line_f)
        end do
        write(momod(modul) % lun_outpu,*)
        
        call chm_memphy(6)
        call compute_Xnormalized_profile_CMC   ! Compute scalar dissipation rate normalized profile
        write(momod(modul) % lun_outpu,*)
        write(momod(modul) % lun_outpu,*) 'Normalized scalar dissipation rate profile (', nZ_CMC_chm, ' entries) for AMC model:'
        do line_f = 1,nZ_CMC_chm
           write(momod(modul) % lun_outpu,*) line_f, ': ', Xnormalized_prof_CMC_chm(line_f)
        end do
        write(momod(modul) % lun_outpu,*)
        call chm_AMC_integrals_CMC   ! Compute integrals for AMC model


        ! Compute inert and equilibrium and do some additional checks
        if ( kfl_rstar == 0 .and. kfl_weigh_in_eq_CMC_chm == 1 ) then
           if (kfl_field_chm(3) == 0)  call runend('CHEMIC REAPHY: alpha field not provided')
           if (kfl_field_chm(2) == 0 .and. kfl_solve_enth_CMC_chm /= 0) &
              call runend('CHEMIC REAPHY: temperature not provided but required to solve enthalpy')

           do index_zbc = 1,2
              ! Do some checks of the mass fractions
              react_scal_bc_CMC_chm = abs(react_scal_bc_CMC_chm)
              sum_Yk = 0.0_rp
              do iclas = 2, nclas_chm+1
                 sum_Yk = sum_Yk + react_scal_bc_CMC_chm(iclas,index_zbc)
              end do
              if (sum_Yk /= 0.0_rp)  then  ! If Yk do not sum up 1.0 => normalize
                 react_scal_bc_CMC_chm(2:nclas_chm+1,index_zbc) = react_scal_bc_CMC_chm(2:nclas_chm+1,index_zbc) / sum_Yk
              else
                 call runend('CHEMIC REAPHY: mass fractions at the mixture fractions boundaries not provided')
              end if

              write(momod(modul) % lun_outpu,*) 'Reactive scalar values at the boundary conditions'
              if (index_zbc == 1) then
                 write(momod(modul) % lun_outpu,*) 'OXIDIZER:'
              else
                 write(momod(modul) % lun_outpu,*) 'FUEL:'
              end if
              write(momod(modul) % lun_outpu,*) 'Temperature:', T_bc_CMC_chm(index_zbc), ' K'
              write(momod(modul) % lun_outpu,*) 'Mass fractions'
              do iclas = 1, nclas_chm
                 write(momod(modul) % lun_outpu,*) speci(iclas)%name, ': ', react_scal_bc_CMC_chm(iclas+1,index_zbc)
              end do
              write(momod(modul) % lun_outpu,*)
              write(momod(modul) % lun_outpu,*)
           end do

           call chm_memphy(7)
           call chm_inert_eq_CMC
        end if

        if (kfl_solve_enth_CMC_chm /= 0) then
           iclas = 1
           found = .false.
           do while( (iclas <= nclas_chm) .and. (.not. found))
              if ('N2' == speci(iclas)%name) then
                 index_N2 = iclas
                 found = .true.
              end if
              iclas = iclas + 1
           end do
           if( .not. found )   index_N2 = 0
        end if

     end if

  end if

end subroutine chm_reaphy
