subroutine chm_reaphy()
  !------------------------------------------------------------------------
  !****f* partis/chm_reaphy
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
  use def_kermod, only : gasco, kfl_prope
  use mod_ecoute, only :  ecoute
  use mod_messages, only : livinf

  implicit none
  integer(ip) :: iclas,ireac,icoef,jcoef,ispec,jspec,ipara,kfl_molar,kfl_rhs,kfl_forwa
  character(5):: names(maxsp_chm),corrected_species=''      ! Names of species 
  character(5):: stoichio_names(3) = ''
  integer(ip) :: stoichio_reaction,kfl_stoic,check_cfi,npara
  real(rp)    :: rpara
  real(rp),allocatable    :: Linear_temp_density_coeffs(:,:)

  if( INOTSLAVE ) then
     !
     ! Initializations (defaults)
     !
     kfl_model_chm = 1                                    ! Defect evolution (1) / Meteo (2) / Mechano-biology (3) / Combustion (4) / CFI model (5)
     kfl_timei_chm = 1                                    ! Existence of du/dt
     kfl_advec_chm = 0                                    ! Existence of (a.grad)u
     kfl_diffu_chm = 1                                    ! Existence of -div[k. grad(u)], and its options (include density)
     kfl_sourc_chm = 0                                    ! Existence and type of source term
     kfl_meteo_chm = 0                                    ! No meteo file
     kfl_tfles_chm = 0                                    ! TFLES model activation (=0 OFF, =1 ON)
     kfl_field_chm(1) = 0                                 ! Initialization by fields (=0 OFF, /=0 ON)
     kfl_field_chm(2) = 0
     kfl_radia_chm = 0                                    ! radiation model activated?
     kfl_cotur_chm = 0                                    ! Turbulence model (RANS) for the CFI combustion model
     kfl_tucfi_chm = 3                                    ! Variance for the CFI combustion model: =0 (NO_VARIANCE), =1 (RPV_VARIANCE), =2 (MF_VARIANCE), =3 (ALL_VARIANCE)
     kfl_uncfi_chm = 0                                    ! Activate unscaled model for non-premixed conditions: =0 (scale model -> to disappear), = 1 (unscale model)
     kfl_benme_chm = 0                                    ! Benchmark in meteo problems
     kfl_ansou_chm = -1                                   ! Flag to identify analytic sounding in case 200                          
     kfl_react_chm = 1                                    ! reaction term on LHS
     kfl_corve_chm = 1                                    ! Correction velocity
     kfl_norma_chm = 0                                    ! Normalize concentrations
     kfl_forwa     = 0                                    ! Forward/backward reaction rate
     lawte_chm     = 0                                    ! Law temperature
     lawde_chm     = 0                                    ! Law density
     lawvt_chm     = 0                                    ! Law terminal velocity
     nclas_chm     = 1                                    ! Number of particle classes
     nodes_chm     = 0                                    ! Number of ODE's
     tfles_chm     = 1.0_rp                               ! Thickening factor for TFLES
     flbet_chm     = 10.0_rp                              ! beta for flame sensor in TFLES
     sourc_chm     = 0.0_rp                               ! Parameter for the source term
     radwt_chm     = 0.0_rp                               ! Wall temperature for radiation model
     sorad_chm     = 0.0_rp                               ! Source radius
     socen_chm     = 0.0_rp                               ! Source center
     denma_chm     = 0.0_rp                               ! Material density [g/cm^3]
     radbi_chm     = 0.0_rp                               ! Bi-molecular radius [cm^3]
     temma_chm     = 0.0_rp                               ! Material temperature [K]
     boltz_chm     = 8.617343e-5_rp                       ! Boltzmann constant [eV/K] = 1.3806504e-23 J/K 
     nreac_chm     = 0                                    ! Number of reactions
     ncoef_chm     = 2                                    ! Number of reaction coefficients
     wprob_chm     = 'GENER'                              ! Problem name
     grnor_chm     = 9.80616_rp                           ! Gravity modulus
     kfl_stoic     = 0_ip                                 ! Stoichiometric mass ratio defined or not
     kfl_arreh_chm = 0_ip                                 ! Correct Arrhenius coefficient with equivalence ratio
     kfl_activ_chm = 0_ip                                 ! Are there activity coefficients present?
     check_cfi     = 0_ip                                 ! Is the table option used for reaction when using the CFI combustion model?
     corrected_species = ''                               !
     stoichio_names    = ''                               !

     !
     ! initialization CFI model
     !
     call chm_memphy(6_ip)
     table_cfi(1)%nvcfi = 0          ! CFI model: number of values for each variable in table
     table_cfi(1)%nvcfi = 2          ! CFI model: degree of freedom in tabulation
     table_cfi(1)%nrcfi = 0          ! CFI model: number of rows in table
     table_cfi(1)%nccfi = 17         ! CFI model: number of columns in table
     table_cfi(1)%nfcfi = 15         ! CFI model: numberof material properties, nonpremixed
     table_cfi(1)%fmima(1) = 0.0_rp  ! CFI model: minimum value of mean mixture fraction, nonpremixed
     table_cfi(1)%fmima(2) = 0.0_rp  ! CFI model: maximum value of mean mixture fraction, nonpremixed
     table_cfi(1)%imima(1) = 0.0_rp  ! CFI model: minimum value of enthalpy, nonadiabatic
     table_cfi(1)%imima(2) = 0.0_rp  ! CFI model: maximum value of enthalpy, nonadiabatic
     table_cfi(1)%nclas     = 2
     table_cfi(1)%ndcfi     = 2
     do icoef = 1,15
       table_cfi(1)%inval(1,icoef) = 0.0_rp
       table_cfi(1)%inval(2,icoef) = 0.0_rp  ! CFI model: inlet material properties for nonpremixed cases
     end do

     ! Viscosity amplification in sponge layer
     sponge_chm = 0
     visco_factor_chm = 1.0_rp
     visco_axis = 0.0_rp
     visco_range = 0.0_rp
 
     kfl_lhsas_chm = 0_ip
     !(SM MM)
     cpcoe_chm = 1004.67_rp
     cvcoe_chm = 717.5_rp
     if (kfl_prope == 0) gasco = 287.17_rp ! Simone no cambies esta linea, Fer
     adgam_chm = 1.4_rp
     pbaro_chm = 100000.0_rp

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

              if( words(1) == 'MODEL' ) then                    ! Model
                 if( words(2) == 'DEFEC' ) then
                    kfl_model_chm = 1
                 else if( words(2) == 'METEO' ) then
                    kfl_model_chm = 2
                 else if( words(2) == 'MECHA' ) then
                    kfl_model_chm = 3
                 else if( words(2) == 'COMBU' ) then
                    kfl_model_chm = 4
                    if( words(3) == 'STAGG' ) then
                      kfl_stagg_chm = 1
                    endif
                    ! Ensure 3 Arrhenius coefficients for forward and backward reactions
                    ! plus 4 a b c d coefficients for Troe form
                    ncoef_chm = 10  
                 else if( words(2) == 'CFIMO' ) then
                    kfl_model_chm   = 5
                    nclas_chm       = 4
                    !
                    ! Type of combustion problem
                    !
                    if( exists('NONPR') .and. exists('NONAD') ) then ! nonpremix, nonadiab
                      nclas_chm           = 4
                      table_cfi(1)%nclas     = 4
                      table_cfi(1)%ndcfi     = 5
                    else if( exists('NONPR') .and. .not.exists('NONAD') ) then ! nonpremix, adiab
                      nclas_chm           = 4
                      table_cfi(1)%nclas     = 4
                      table_cfi(1)%ndcfi     = 4
                    else if( .not.exists('NONPR') .and. .not.exists('NONAD') ) then ! premix, adiab
                      nclas_chm           = 2
                      table_cfi(1)%nclas     = 2
                      table_cfi(1)%ndcfi     = 2
                    else if( .not.exists('NONPR') .and. exists('NONAD') ) then ! premix, nonadiab
                      nclas_chm           = 2
                      table_cfi(1)%nclas     = 2
                      table_cfi(1)%ndcfi     = 3
                    end if
                    !
                    ! Activation of variances
                    !
                    if( exists('NOVAR') ) then
                      kfl_tucfi_chm = 0
                    endif
                    if( exists('RPVVA') ) then
                      kfl_tucfi_chm = 1
                    endif
                    if( exists('MFVAR') ) then
                      kfl_tucfi_chm = 2
                    endif
                    if( exists('ALLVA') ) then
                      kfl_tucfi_chm = 3
                    endif
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
                      call runend('CHEMIC REAPHY: Turbulence model not valid for CFI combustion') 
                    endif
                    !
                    ! Unscale model for non-premixed conditions
                    !
                    if( exists('UNSCA') ) then
                      if( .not.exists('NONPR')) call runend('CHEMIC REAPHY: Unscale model only available for non-premixed combustion')
                      kfl_uncfi_chm = 1
                    endif
                    !
                    ! Radiation model 
                    !
                    if( exists('RADIA') ) then
                      kfl_radia_chm = 1
                      table_cfi(1)%nccfi = 19 ! CFI model: number of columns in table
                      table_cfi(1)%nfcfi = 17 ! CFI model: numberof material properties, nonpremixed
                      if ( exists('WALLT') ) then
                         radwt_chm = getrea('WALLT',0.0_rp,'#Wall temperature radiation') 
                         kfl_radia_chm = 2
                      end if
                    end if

                 end if

              else if( words(1) == 'PROBL' ) then               ! Problem name 
                 wprob_chm = words(2)

              else if( words(1) == 'TEMPO' ) then               ! Temporal evolution
                 if( words(2) == 'ON   ' ) then
                    kfl_timei_chm = 1
                 else
                    kfl_timei_chm = 0
                 end if

              else if( words(1) == 'CONVE' ) then               ! Convective term
                 if( words(2) == 'ON   ' ) then
                    kfl_advec_chm = -1
                    if( words(3) == 'FUNCT' ) then
                       kfl_advec_chm = getint('FUNCT',0_ip,'#Velocity function')
                    else if(words(3)=='VELOC' ) then
                       kfl_advec_chm = -2 
                    else if(words(3)=='FILE ' ) then
                       kfl_advec_chm = -1
                    end if
                 else if( words(2) == 'OFF  ' ) then
                    kfl_advec_chm = 0
                 end if
                 if (exists('CORRE')) then
                    kfl_corve_chm = 1
                 else
                    kfl_corve_chm = 0
                 endif

              else if( words(1) == 'DIFFU' ) then               ! Diffusion
                 if( words(2) == 'ON   ' ) then
                    kfl_diffu_chm = 1                           ! Diffusion using mass fractions for combustion
                 else if( words(2) == 'COMPLE' ) then
                    kfl_diffu_chm = 2                           ! Diffusion using mole fractions for combustion
                 else
                    kfl_diffu_chm = 0
                 end if

              else if( words(1) == 'REACT' ) then               ! Reaction on LHS or RHS
                 if( words(2) == 'RHS  ' ) then
                    kfl_react_chm = 0
                 else
                    kfl_react_chm = 1
                 end if

              else if( words(1) == 'NUMBE' ) then
                 nclas_chm=getint('NUMBE',1_ip,'#Number of classes')

              else if( words(1) == 'SPECI' ) then
                 nclas_chm = 0
                 call ecoute('chm_reaphy')
                 do while( words(1) /= 'ENDSP' )
                    nclas_chm=nclas_chm+1
                    if (nclas_chm > maxsp_chm) call runend('CHEMIC REAPHY: Maximum number of species exceeded')
                    names(nclas_chm) = words(1)
                    if (exists('CORRE') .or. exists('NORMA')) corrected_species =  words(1)
                    call ecoute('chm_reaphy')
                 end do

              else if( words(1) == 'PDE  ' ) then
                 nclas_chm=getint('PDE  ',1_ip,'#Number of classes')

              else if( words(1) == 'ODE  ' ) then
                 nodes_chm=getint('ODE  ',0_ip,'#Number of ODEs')

              else if( words(1) == 'ODES ' ) then
                 nodes_chm=getint('ODES ',0_ip,'#Number of ODEs')

              else if( words(1) == 'SOURC' ) then               ! Source term
                 if( words(2) == 'CONST' ) then
                    kfl_sourc_chm =  1
                    sourc_chm = getrea('VALUE',0.0_rp,'#Source term')
                 else if( words(2) == 'SPHER' ) then
                    kfl_sourc_chm =  2
                    sourc_chm     = getrea('VALUE',0.0_rp,'#Source term')
                    sorad_chm     = getrea('RADIU',0.0_rp,'#Source radius')
                    socen_chm(1)  = getrea('X    ',0.0_rp,'#Source center')
                    socen_chm(2)  = getrea('Y    ',0.0_rp,'#Source center')
                    socen_chm(3)  = getrea('Z    ',0.0_rp,'#Source center')
                 else if(words(2)=='DEBUG') then
                    kfl_sourc_chm = -1   
                 else if(words(2)=='FILE ') then
                    if(exists('ASCII')) then
                       kfl_sourc_chm = -2
                    else
                       kfl_sourc_chm = -3                       
                    end if
                 else if(words(2)=='MOIST') then
                    kfl_sourc_chm = -4
                 else if( words(2) == 'REACT' ) then
                    kfl_sourc_chm =  4
                 end if
                 if( exists('LHS  ') ) then ! Left hand side assembly of mass source terms
                    kfl_lhsas_chm = 1_ip  
                 endif

              else if( words(1) == 'METEO' ) then             ! Meteo file
                 if(words(2)=='DEBUG') then
                    kfl_meteo_chm = 1    
                 else if(words(2)=='ASCII') then
                    kfl_meteo_chm = 2
                 else if(words(2)=='NETCD') then
                    kfl_meteo_chm = 3  
                 else if(words(2)=='OFF  ') then
                    kfl_meteo_chm = 0
                 else if(words(2)=='FUNCT') then
                    !
                    ! Initial field of tracers defined by a function (routine) SM
                    !
                    kfl_meteo_chm = -1
                 else
                    kfl_meteo_chm = 4
                 end if

              else if(words(1)=='TERMI') then               ! Terminal velocity
                 if(words(2)=='OFF  ') then
                    lawvt_chm = 0
                 else if(words(2)=='GANSE') then
                    lawvt_chm = 1
                 else if(words(2)=='WILSO') then
                    lawvt_chm = 2
                 else if(words(2)=='DELLI') then
                    lawvt_chm = 3
                 else if(words(2)=='W2PLA') then
                    lawvt_chm = 4
                 end if

              else if (words(1)=='ACTIV') then
                 if (words(2)=='ON   ') then
                    kfl_activ_chm = 1_ip
                 else
                    kfl_activ_chm = 0_ip
                 endif
              
              else if( words(1) == 'COMBU' ) then           ! Tubulent combustion model
                 if(  words(2) == 'TFLES' ) then
                    kfl_tfles_chm = 1
                    if( words(3) == 'FACTO' ) then
                       tfles_chm = getrea('FACTO',1.0_rp,'#Thickening factor')
                    endif
                    if( words(4) == 'BETA ') then
                       flbet_chm = getrea('BETA ',10.0_rp,'#Flame sensor constant')
                    endif
                    write(momod(modul) % lun_outpu,*) 'TFLES ACTIVATED WITH THICKENING FACTOR =',tfles_chm
                    write(momod(modul) % lun_outpu,*) 'Flame sensor controlling transition beta =',flbet_chm
                 end if

              else if(words(1).eq.'TURBU') then
                 if(words(2)=='LESMO') then
                    kfl_cotur_chm = -1_ip
                 else if(words(2)=='RANSM' .or. words(2)=='FROMT') then
                    kfl_cotur_chm = 1_ip
                 end if
              end if 
               
              call ecoute('chm_reaphy')
           end do

        else if( words(1) == 'PROPE' ) then
           !
           ! Properties
           !
           call chm_memphy(1_ip)
           if (kfl_model_chm==4) allocate (Linear_temp_density_coeffs(nclas_chm,8))
           interaction_chm = 0.0_rp                             ! Binary interaction coefficients for UNIQUAC model
     
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

              if( words(1) == 'DENSI' ) then                    ! Material density
                 if( exists('CONST') ) then
                    denma_chm = getrea('VALUE',1.0_rp,'#Material density')
                 else
                    denma_chm = getrea('DENSI',1.0_rp,'#Material density')
                 end if

              else if( words(1) == 'BIMOL' ) then               ! Bi-molecular radius
                 radbi_chm = getrea('BIMOL',0.0_rp,'#Bimolecular radius')

              else if( words(1) == 'TEMPE' ) then               ! Material temperature
                 temma_chm(1:npart_chm) = param(1:npart_chm)

              else if( words(1) == 'BOLTZ' ) then               ! Boltzmann constant
                 boltz_chm = getrea('BOLTZ',0.0_rp,'#Boltzmann constant')

              else if ( words(1) == 'GASCO' ) then               ! Universal gas constant
                 gasco = getrea('GASCO',8.3144621_rp,'#Gas constant')

              else if( words(1) == 'NUMBE' ) then               ! Number of reactions
                 nreac_chm = getint('NUMBE',1_ip,'#Number of reactions')
                 call chm_memphy(3_ip)                          ! Allocate memory for reactions

              else if( words(1) == 'COEFF' ) then               ! Number of coefficients
                 ncoef_chm = getint('COEFF',1_ip,'#Number of coefficients')

              else if( words(1) == 'LAWTE' ) then               ! Temperature law
                 if( words(2) == 'CONST' ) then
                    lawte_chm = 0
                    if (exists('VALUE')) temma_chm = getrea('VALUE',300.0_rp,'#Material density')
                 else if( words(2) == 'FILE ' ) then
                    lawte_chm = -1
                 else if( words(2) == 'TEMPE' ) then
                    lawte_chm = -2
                 else if( words(2) == 'FUNCT' ) then
                    lawte_chm = getint('FUNCT',0_ip,'#Temperature function')
                 end if

              else if( words(1) == 'LAWDE' ) then                ! Density law
                 if( words(2) == 'CONST' ) then
                    lawde_chm =  0
                    if (exists('VALUE')) denma_chm = getrea('VALUE',1.0_rp,'#Material density')
                 else if( words(2) == 'FILE ') then
                    lawde_chm = -1
                 else if( words(2) == 'DENSI') then
                    lawde_chm = -2
                 else if( words(2) == 'LOWMA') then
                    lawde_chm = -3
                 else if( words(2) == 'TEMPE') then
                    lawde_chm = -4
                 else if( words(2) == 'EXTER') then
                    lawde_chm = -10
                 end if

              else if( words(1) == 'SPONG' ) then               ! SPonge layer
                 sponge_chm = 1
                 visco_factor_chm = getrea('FACTO',1.0_rp,'#Viscosity amplification factor')
                 visco_axis(1) = getrea('OX   ',0.0_rp,'#Viscosity amplification axes')
                 visco_axis(2) = getrea('OY   ',0.0_rp,'#Viscosity amplification axes')
                 visco_axis(3) = getrea('OZ   ',0.0_rp,'#Viscosity amplification axes')
                 visco_range(1) = getrea('START',1.0_rp,'#Coordinate to start amplification')
                 visco_range(2) = getrea('END  ',2.0_rp,'#Coordinate to end amplification')

              else if( words(1) == 'LAWDI' ) then                ! Diffusion law
                 if( words(2) == 'CONST' ) then
                    if( exists('SPECY') ) then
                       ispec = getint('SPECY',1_ip,'#Specy number')                    
                       if (exists('VALUE')) diffu_chm(1,ispec) = getrea('VALUE',1.0_rp,'#Diffusion constatnt')
                       lawdi_chm(1,ispec) = 0
                    else
                       lawdi_chm =  0
                       diffu_chm = getrea('VALUE',1.0_rp,'#Diffusion constatnt')
                    end if
                 else if( words(2) == 'PRAND') then
                    lawdi_chm = 1
                    if(exists('RANS ')) then
                        if (exists('TURBU')) diffu_chm(1,1) = getrea('TURBU',0.9_rp,'#Turbulent Schmidt number')
                        lawdi_chm = 5
                        write(momod(modul) % lun_outpu,*) 'TURBULENT SCHMIDT NUMBER FOR RANS:',diffu_chm(1,1)
                        write(momod(modul) % lun_outpu,*) ''
                        if (kfl_coupl(ID_CHEMIC,ID_TURBUL) == 0 ) then
                            write(momod(modul) % lun_outpu,*) 'CHEMIC: TURBUL is not activated'
!                            call runend('CHEMIC: For RANS simulation TURBUL must be activated') 
                        endif
                    elseif ( exists('LES  ') ) then
                        if (exists('TURBU')) diffu_chm(1,1) = getrea('TURBU',0.9_rp,'#Turbulent Schmidt number')
                        lawdi_chm = 5
                        write(momod(modul) % lun_outpu,*) 'TURBULENT SCHMIDT NUMBER FOR LES:',diffu_chm(1,1)
                        write(momod(modul) % lun_outpu,*) ''
                        if (kfl_coupl(ID_CHEMIC,ID_NASTAL) == 0 .and. kfl_coupl(ID_CHEMIC,ID_NASTIN) == 0) then
                            call runend('CHEMIC: For LES simulation NASTAL/NASTIN must be activated') 
                        endif
                    endif
                    if (kfl_tfles_chm >= 1_ip) lawdi_chm = 1 ! DTFLES model

                 else if( words(2) == 'MOLEC') then
                    lawdi_chm = 2
                    call runend('CHEMIC: Molecular diffusion not programmed yet')
                 else if( words(2) == 'UNIFO') then
                    lawdi_chm = 3
                    diffu_chm(1,1) = getrea('PRAND',0.71_rp,'#Uniform Prandtl number')
                    diffu_chm(2,1) = getrea('LEWIS',1.0_rp,'#Uniform Lewis number')                    
                 else if( words(2) == 'WATER') then
                    lawdi_chm = 4
                 end if

              else if (words(1)=='STOIC') then
                 stoichio_reaction=getint('REACT',1_ip,'#Fuel component of the reaction')
                 kfl_stoic = 1_ip

              else if( words(1) == 'CLASS' ) then

                 !-------------------------------------------------------
                 !
                 ! PDE's: Classes
                 !
                 !-------------------------------------------------------

                 if( words(2) == 'DIFFU' ) then
                    call ecoute('chm_reaphy')
                    iclas = 0
                    do while( words(1) /= 'ENDCL' )
                       iclas = iclas + 1
                       if( iclas > nclas_chm ) call runend('TOO MANY CLASS PROPERTIES DEFINED IN DIFFUSION')
                       diffu_chm(1,iclas) = param(1)
                       diffu_chm(2,iclas) = param(2)
                       call ecoute('chm_reaphy')
                    end do
                    if( iclas /= nclas_chm ) call runend('DIFFUSION CLASS PROPERTIES ARE MISSING')

                 else if( words(2) == 'LAWDI' ) then
                    call ecoute('chm_reaphy')
                    iclas = 0
                    do while( words(1) /= 'ENDCL' )
                       iclas = iclas + 1
                       if( iclas > nclas_chm ) call runend('TOO MANY CLASS PROPERTIES DEFINED IN DIFFUSION')
                       !
                       ! In meteo: law for horizontal diffusion
                       !
                       if( words(1) == 'CONST') then
                          lawdi_chm(1,iclas) = 1
                       else if( words(1) == 'RAMS ') then
                          lawdi_chm(1,iclas) = 2
                       else if( words(1) == 'CMAQ ') then
                          lawdi_chm(1,iclas) = 3
                       end if
                       !
                       ! In meteo: law for vertical diffusion
                       !
                       if( words(2) == 'CONST') then
                          lawdi_chm(2,iclas) = 1
                       else if( words(2) == 'SIMIL') then
                          lawdi_chm(2,iclas) = 2
                       end if

                       call ecoute('chm_reaphy')
                    end do
                    if( iclas /= nclas_chm ) call runend('DIFFUSION CLASS PROPERTIES ARE MISSING')

                 else if( words(2) == 'EQUIL' ) then
                    call ecoute('chm_reaphy')
                    iclas = 0
                    do while( words(1) /= 'ENDCL' )
                       iclas = iclas + 1
                       if( iclas > nclas_chm ) call runend('TOO MANY CLASS PROPERTIES DEFINED IN EQUILIBRIUM')
                       equil_chm(1,iclas) = param(1)
                       equil_chm(2,iclas) = param(2)
                       call ecoute('chm_reaphy')
                    end do
                    if( iclas /= nclas_chm ) call runend('EQUILIBRIUM PROPERTIES ARE MISSING')  

                 else if( words(2) == 'PARTI' ) then
                    if( kfl_model_chm /= 2 ) call runend('PARTICLES PROPERTIES ONLY FOR METEO MODEL')
                    call ecoute('pts_reaphy')
                    iclas = 0
                    do while(words(1)/='ENDCL')
                       iclas = iclas +1
                       if( iclas > nclas_chm ) call runend('TOO MANY CLASS PROPERTIES DEFINED')
                       diame_chm(iclas) = param(1)                  ! diameter in mm
                       diame_chm(iclas) = diame_chm(iclas)*1d-3     ! diameter in m
                       rhopa_chm(iclas) = param(2)
                       spher_chm(iclas) = param(3)
                       fract_chm(iclas) = param(4)
                       call ecoute('pts_reaphy')
                    end do
                    if( iclas /= nclas_chm ) call runend('CLASS PROPERTIES ARE MISSING')
                    !
                    !  Calculates shape_chm(nclas_chm) depending on the model
                    !  
                    do iclas = 1,nclas_chm
                       call chm_setpsi(shape_chm(iclas),spher_chm(iclas),diame_chm(iclas),lawvt_chm)
                    end do

                 end if

              else if( words(1) == 'ODES ' ) then

                 !-------------------------------------------------------
                 !
                 ! ODE's
                 !
                 !-------------------------------------------------------
              else if( words(1) == 'SPECI' .and. kfl_model_chm == 4 ) then 
                 !-------------------------------------------------------
                 !
                 ! Species properties (only for combustion)
                 !
                 !-------------------------------------------------------
!
                 write(momod(modul) % lun_outpu,*)'------------------'
                 write(momod(modul) % lun_outpu,*)'SPECIES PROPERTIES'
                 write(momod(modul) % lun_outpu,*)'------------------'
!
                 ispec = 0
                 do while( words(1) /= 'ENDSP' )
                    if( words(1) == 'TYPE ' ) then               
                       do iclas = 1,nclas_chm
                          if (words(2) == names(iclas) ) then ! This type is used in the simulation
                             ispec = ispec + 1
                             if( ispec > nclas_chm ) &
                                  call runend('CHEMIC: More species properties than needed (probably defined twice.)') 
                             speci(ispec) % name = words(2)
                             if (exists('FUEL ')) then
                                 stoichio_names(1)= words(2)
                             else if (exists('OXYGE')) then
                                 stoichio_names(2)=words(2)
                             else if (exists('NITRO')) then
                                 stoichio_names(3)=words(2)
                             endif
                             speci(ispec) % cpcoe = 0.0_rp
                             speci(ispec) % visco = 0.0_rp
                             speci(ispec) % lawvi = 1_ip
                             speci(ispec) % weigh = 1.0_rp
                             speci(ispec) % entha = 0.0_rp
                             speci(ispec) % prand = 0.0_rp
                             speci(ispec) % densi = 1.0_rp
                             speci(ispec) % activ = 1.0_rp
                             speci(ispec) % trang = 0.0_rp
                             if ( corrected_species == names(iclas) )  kfl_norma_chm = ispec
                             kfl_molar=0
!
                             write(momod(modul) % lun_outpu,*) ''
                             write(momod(modul) % lun_outpu,*) 'SPECIES:  ',speci(ispec) % name
!
                             do while( words(1) /= 'ENDTY' )
                                if( words(1) == 'NAME ' ) then               
                                   speci(ispec) % name = words(2)
                                else if( words(1) == 'VISCO' ) then               
                                   speci(ispec) % visco = param(1:2)
                                   write(momod(modul) % lun_outpu,*) 'VISCOSITY',speci(ispec) % visco
                                   if (exists('SUTHE')) then
                                      speci(ispec) % lawvi = 1_ip
                                   else if (exists('LINEA')) then
                                      speci(ispec) % lawvi = 2_ip
                                   endif
                                else if( words(1) == 'ATOMI' .or. words(1) == 'MOLEC' ) then               
                                   speci(ispec) % weigh = param(1)
                                   write(momod(modul) % lun_outpu,*) 'MOLECULAR WEIGHT',speci(ispec) % weigh
                                   if (exists('GRAMS')) speci(ispec) % weigh = speci(ispec) % weigh / 1000.0_rp
                                else if( words(1) == 'DENSI' ) then  
                                   if (exists('TEMPE')) then
                                      ! We expect 8 data points, 4 per each phase A and B 
                                      !  T_1_A densi_1_A T_2_A densi_2_A T_1_B densi_1_B T_2_B densi_2_B
                                      Linear_temp_density_coeffs(ispec,:) = param(1:8)
                                   else ! We expect one density per phase
                                      speci(ispec) % densi(1:2) = param(1:2)
                                      write(momod(modul) % lun_outpu,*) 'DENSITY',speci(ispec) % densi(1:2)
                                   endif

                                else if( words(1) == 'ACTIV' ) then  
                                   speci(ispec)%activ(1)=getrea('R    ',1.0_rp,'#UNIFAC R parameter')
                                   speci(ispec)%activ(1)=getrea('Q    ',1.0_rp,'#UNIQUAC Q parameter')
                                   
                                else if( words(1) == 'ENTHA' ) then               
                                   speci(ispec) % entha = param(1:2)
                                   write(momod(modul) % lun_outpu,*) 'ENTHALPY OF FORMATION',speci(ispec) % entha
                                   if (exists('MOLAR')) speci(ispec) % entha(1) = speci(ispec) % entha(1) / speci(ispec) % weigh
                                else if( words(1) == 'PRAND' ) then               
                                   speci(ispec) % prand = param(1)
                                   write(momod(modul) % lun_outpu,*) 'PRANDTL NUMBER',speci(ispec) % prand
                                else if( words(1) == 'LEWIS' ) then               
                                   speci(ispec) % lewis = param(1)
                                   write(momod(modul) % lun_outpu,*) 'LEWIS NUMBER',speci(ispec) % prand
                                else if( words(1) == 'CPCOE' .or. words(1) == 'SHOMA' ) then               
                                   speci(ispec) % cpcoe(1:8,1) = param(3:10)
                                   speci(ispec) % trang = 0.0_rp
                                   if (exists('MOLAR')) kfl_molar = 1
                                else if( words(1) == 'CPRAN' ) then               
                                   ipara = getint('CPRAN',1_ip,'#Number of temp ranges')
                                   if (ipara > 4) call runend('CHM_REAPHY: max number of Cp temperature ranges exceeded')
                                   if (exists('MOLAR')) kfl_molar = 1
                                   do icoef = 1,ipara
                                      call ecoute('chm_reaphy')
                                      speci(ispec) % cpcoe(1:8,icoef) = param(3:10)
                                      speci(ispec) % trang(icoef:icoef+1) = param(1:2)
                                   enddo
                                endif
                                call ecoute('chm_reaphy')
                             enddo
                             if (kfl_molar ==1) then  ! J/mol K (Molar) Cp coefficients get transformed to J/Kg K
                                do ipara=1,4
                                   do icoef=1,8
                                      speci(ispec) % cpcoe(icoef,ipara) = speci(ispec) % cpcoe(icoef,ipara)/ speci(ispec) % weigh
                                   enddo
                                enddo
                             endif
                          endif
                       enddo
                    endif
                    call ecoute('chm_reaphy')
                 end do
                 if( ispec < nclas_chm ) call runend('CHEMIC: Properties for some species are missing')  

              else if( words(1) == 'REACT' ) then
                 !-------------------------------------------------------
                 !
                 ! Reactions
                 !
                 !-------------------------------------------------------
                 if( words(2) == 'BINDI') then

                    call ecoute('chm_reaphy')
                    ireac = 0
                    do while( words(1) /= 'ENDRE' )
                       ireac = ireac + 1
                       if( ireac > nreac_chm ) call runend('TOO MANY BINDING ENERGIES')
                       do icoef = 1,ncoef_chm
                          react_chm(icoef,ireac) = param(icoef)
                       enddo
                       call ecoute('chm_reaphy')
                    end do
                    if( ireac /= nreac_chm ) call runend('BINDING COEFFICIENTS ARE MISSING')  

                 else if( words(2) == 'CAPTU' ) then

                    call ecoute('chm_reaphy')
                    ireac = 0
                    do while(words(1)/='ENDRE')
                       ireac = ireac + 1
                       if( ireac > nreac_chm ) call runend('TOO MANY REACTIONS DEFINED IN CAPTURE RADIUS')
                       radiu_chm(ireac) = param(1)
                       call ecoute('chm_reaphy')
                    end do
                    if( ireac > nreac_chm ) call runend('CAPTURE RADIUS CLASS PROPERTIES ARE MISSING')  

                 else if( words(2) == 'KINET' .or. words(2) == 'CHEMK' ) then   
                    ! Chemical kinetics similar to CHEMKIN format
                    !
                    !  [n1] R1 + [ [n2] R2 ] + [ [n3] R3 ] + [M] -> [m1] P1 + [ [m2] P2 ] + [ [m3] P3] , [BROKEN_ORDER] [TROE] [LINDEMANN]
                    !  A_f beta_f E_f A_b beta_b E_b a b c d  (last four are for Troe and Lindemann forms, 10 params total)
                    !  ORDER, [FORWARD,BACKWARD]  x1=n1 x2=n2 [x3=n3] y1=m1 [y2=m2] [y3=m3]   (if broken order then the orders go here, default is 0.0)
                    !  EFFICIENCY  X1=e1 X2=e2 Y3=e3    (if M is present, efficiency line changes efficiencies from default 1.0)
                    !  Preferred units are MKS: [A] = m^3/(mol s), [E] = J/mol 
                    !
                    kfl_split_plus = 1 !!Activate the + as a separator, CAREFUL not to use 3.1E+12 for example
                    if (exists('EQUIV')) then
                       kfl_arreh_chm=1
                    endif
                    call ecoute('chm_reaphy')
                    ireac = 0
!
                    write(momod(modul) % lun_outpu,*)''
                    write(momod(modul) % lun_outpu,*)'------------------'
                    write(momod(modul) % lun_outpu,*)'CHEMICAL REACTIONS'
                    write(momod(modul) % lun_outpu,*)'------------------'
                    write(momod(modul) % lun_outpu,*)''
!
                    do while( words(1) /= 'ENDRE' )
                       if (words(1) == 'EFFIC') then
                          do ispec = 1,nclas_chm
                             effic_chm(ispec,ireac) = getrea(speci(ispec)%name,1.0_rp,'#Efficiency of species')
                          enddo
                       else if (words(1) == 'ORDER') then
                          if (exists('BACKW')) then
                             kfl_forwa=2
                          else
                             kfl_forwa=1
                          endif
                          do ispec = 1,nclas_chm
                             order_chm(ispec,ireac,kfl_forwa) = getrea(speci(ispec)%name,0.0_rp,'#Order of species in reaction')
                          enddo
                       else ! Reaction line, parameters must follow right afterwards.
                          ireac = ireac + 1
                          if( ireac > nreac_chm ) call runend('TOO MANY CHEMICAL REACTIONS')
                          do icoef=1,ncoef_chm
                             lreac_chm(ireac)%l(icoef) = 0
                          enddo
                          do ispec = 1,nclas_chm
                             order_chm(ispec,ireac,:) = 0.0_rp
                             stoic_chm(ispec,ireac,:) = 0.0_rp
                          enddo
                          if (exists('LINDE')) then 
                             lreac_chm( ireac ) %l ( 5 ) = 1  ! Lindemann form for the coefficients
                             lreac_chm( ireac ) %l ( 4 ) = 1
                          else if (exists('TROE ')) then
                             lreac_chm( ireac ) %l ( 5 ) = 2  ! Troe form for the coefficients
                             lreac_chm( ireac ) %l ( 4 ) = 1
                          end if
                          if (exists('BROKE')) then
                             lreac_chm( ireac ) %l ( 3 ) = 1
                          endif
                          iclas = 0
                          if (words(1) == '') then ! line starts with a stoichiometric coefficient
                             icoef = 2
                          else
                             icoef = 1 ! Line starts with reactant name so that 1 is assumed
                             param(maxwp) = 1.0_rp ! we impose it on the last parameter because it is usually zero
                          endif                          
                          kfl_rhs = 0
                          do while(icoef <= maxwp )
                             if (icoef > 1) then 
                                jcoef = icoef-1
                             else 
                                jcoef = maxwp ! Periodic reading of parameters
                             endif
                             if (words(icoef)=='->   ' .or. words(icoef)=='<->  ') then
                                kfl_rhs = 1 ! switch to right hand side
                             else if (words(icoef)=='M    ' .or. words(icoef)=='CHAPE') then
                                lreac_chm( ireac ) %l ( 4 ) = 1  ! Third body reaction
                                print *,'--xx WARNING xx-- THIRD BODY REACTIONS NOT TESTED YET'
                             else
                                do ispec = 1,nclas_chm
                                   if (words(icoef) == speci(ispec)%name ) then 
                                      if (param(jcoef)==0.0_rp) param(jcoef)=1.0_rp !Default is 1 mol
                                      if (kfl_rhs == 0) then
                                         lreac_chm(ireac)%l(1) = lreac_chm(ireac)%l(1) + 1
                                         iclas = 5 + lreac_chm(ireac)%l(1)
                                         stoic_chm(ispec,ireac,1) = param(jcoef)
                                      else
                                         lreac_chm(ireac)%l(2) = lreac_chm(ireac)%l(2) + 1
                                         iclas = 5 + lreac_chm(ireac)%l(1) + lreac_chm(ireac)%l(2)
                                         stoic_chm(ispec,ireac,2) = param(jcoef)
                                      endif
                                      lreac_chm( ireac ) %l ( iclas ) = ispec
                                   endif
                                enddo
                             endif
                             icoef = icoef + 1
                          enddo
                          if (lreac_chm( ireac ) %l ( 3 ) /= 1 ) order_chm(:,ireac,:) = stoic_chm(:,ireac,:)
                          ! Finally read the coefficients, they should be 6 for bimolecular reactions and 4 more for Troe forms
                          kfl_split_plus = 0 !! Temporarile deactivate the + as a separator
                          call ecoute('chm_reaphy')
                          do ipara = 1,ncoef_chm
                             react_chm(ipara,ireac) = param(ipara)
                          enddo
                           kfl_split_plus = 1 !! Activate the + as a separator
                       end if
                       call ecoute('chm_reaphy')
                    end do
                    kfl_split_plus = 0 !! Deactivate the + as a separator
                    if( ireac /= nreac_chm ) call runend('CHEMICAL REACTIONS ARE MISSING')  

                    !
                    ! Write out in *.chm.log the chemical kinetic mechanisms
                    !
                    !
                    do ireac=1,nreac_chm
                       !
                       write(momod(modul) % lun_outpu,*)'Reaction ',ireac
                       write(momod(modul) % lun_outpu,*)''
                       !
                       if (kfl_forwa == 1 ) then
                          write(momod(modul) % lun_outpu,*)'Only forward reaction'
                       else if (kfl_forwa == 2 ) then
                          write(momod(modul) % lun_outpu,*)'Forward & backward reactions'
                       else
                          !!call runend('UNDEFINED NUMBER OF REVERSIBLE REACTIONS')
                       endif
                       !   
                       write(momod(modul) % lun_outpu,*)' '
                       do ispec=1,nclas_chm
                          write(momod(modul) % lun_outpu,*)'Efficiency species    ',speci(ispec)%name,' = ',effic_chm(ispec,ireac)
                          write(momod(modul) % lun_outpu,*)'Order species         ',speci(ispec)%name,' = ',order_chm(ispec,ireac,:)
                          write(momod(modul) % lun_outpu,*)'Stoichoimetric coeff. ',speci(ispec)%name,' = ',stoic_chm(ispec,ireac,:)
                          write(momod(modul) % lun_outpu,*)' '
                       enddo
                       !
                       do ipara = 1,ncoef_chm
                          write(momod(modul) % lun_outpu,*)'Arrhenius coefficients ',ipara,' = ',react_chm(ipara,ireac)
                       enddo
                       write(momod(modul) % lun_outpu,*)' '
                    enddo
                    !
                    ! CFI combustion model: read thermochemical database 
                    ! 
                 else if( words(2) == 'TABLE' ) then
                    if (kfl_model_chm /= 5) call runend('CHEMIC REAPHY: Table should only be used with the CFI combustion model')
                    check_cfi = 1_ip

                    speci(1_ip)%name = 'CMEAN'
                    speci(2_ip)%name = 'CVAR'
                    if (nclas_chm > 2_ip) then
                       speci(3_ip)%name = 'FMEAN'
                       speci(4_ip)%name = 'FVAR'
                    endif

                    call ecoute('chm_reaphy')
                    call chm_loatab()

                    do while(words(1)/='ENDRE')
                       call ecoute('chm_reaphy')
                    end do

                    write(momod(modul) % lun_outpu,*)''
                    write(momod(modul) % lun_outpu,*)'---------------------------------------------'
                    write(momod(modul) % lun_outpu,*)'CFI-MODEL: THERMOCHEMICAL DATABASE FROM TABLE'
                    write(momod(modul) % lun_outpu,*)'---------------------------------------------'
                    write(momod(modul) % lun_outpu,*)''
                    if (table_cfi(1)%nclas == 2) then 
                       write(momod(modul) % lun_outpu,*)'PREMIXED COMBUSTION'
                    else 
                       write(momod(modul) % lun_outpu,*)'NON-PREMIXED COMBUSTION'
                    end if
                    write(momod(modul) % lun_outpu,*)''
                    write(momod(modul) % lun_outpu,*)'NUMBER OF VALUES FOR MEAN OF RPV:     ', table_cfi(1)%nvcfi(1)
                    write(momod(modul) % lun_outpu,*)'NUMBER OF VALUES FOR VARIANCE OF RPV: ', table_cfi(1)%nvcfi(2)
                    write(momod(modul) % lun_outpu,*)'NUMBER OF VALUES FOR MEAN OF MF:      ', table_cfi(1)%nvcfi(3)
                    write(momod(modul) % lun_outpu,*)'NUMBER OF VALUES FOR VARIANCE OF MF:  ', table_cfi(1)%nvcfi(4)


                 else if( wprob_chm /= 'GENER' ) then

                    call chm_usrrea()

                 end if

                 do while( words(1) /= 'ENDRE' )
                    call ecoute('chm_reaphy')
                 end do

             else if (words(1)=='INTER') then !Read binary interaction coefficients
                 ! Syntax= SPEC_j=  SPEC_1=x, SPEC_2=y, ... 
                 ! DO NOT PUT SPEC_j here on the rhs
                 call ecoute('chm_reaphy')
                 do while(words(1).ne.'ENDIN')
                    ispec=-1
                    do jspec = 1,nclas_chm
                       if (words(1).eq.speci(jspec)%name) ispec=jspec
                    enddo
                    if (ispec .ne. -1) then
                       do jspec = 1,nclas_chm
                          interaction_chm(ispec,jspec) = getrea(speci(jspec)%name,0.0_rp,'#Binary interaction coefficients')
                       enddo
                    endif
                    call ecoute('chm_reaphy')
                 end do
                 !
                 ! CFI combustion model: mass fractions at equilibrium for scaling RPV for non-premixed conditions
                 !
             else if (words(1)=='MASSF') then

                 call chm_read_mass_equilibrium() 

                 do while( words(1) /= 'ENDMA' )
                    call ecoute('chm_reaphy')
                 end do

             end if
                 call ecoute('chm_reaphy')
           end do

        else if (words(1)=='FIELD') then

           kfl_field_chm(1) = getint('SPECI',0_ip,'#Species starting fields')
           kfl_field_chm(2) = getint('TEMPE',0_ip,'#Temperature field')
           
        else if (words(1)=='METEO') then

           call ecoute('chm_reaphy')
           do while(words(1)/='ENDME')
              if (words(1)=='BENCH') then
                 ! benchmark name
                 if ( words(2) == 'KESSL' .or. &
                      words(2) == 'SQUAL' .or. &
                      words(2) == 'MOIST') then
                    kfl_benme_chm = 200                 ! Squall line with Kessler
                    if ( words(3) == 'SIMPL' .or.  words(3) == 'SUPER') then
                       kfl_benme_chm = 201

                       if ( words(3) == 'SIMPL' .or.  words(3) == 'SUPER') then
                          kfl_benme_chm = 201

                       else if (words(3)=='KLAAS ' .or. words(3)=='KC85 ') then
                          kfl_benme_chm = 203

                       else if (words(3)=='GRABO' .or. words(3)=='GK91 ' .or. words(3)=='G2007') then
                          kfl_benme_chm = 210              !Grabowski JAS 2007
                       end if

                    else if ( words(2) == 'TRACE') then
                       kfl_benme_chm = 100                 ! warm bubble + tracer
           
                    end if
                 else if (words(2)=='WARMB') then
                    if (words(3)=='AHMAD') &
                         kfl_benme_chm = 100             ! warm bubble with tracer
                 end if
                 ! object definition
              else if (words(1)=='ANALY' .or. words(1)=='ANSOU') then
                 kfl_ansou_chm = 1

              end if
              call ecoute('chm_reaphy')
           end do

        end if
     end do

     if( kfl_benme_chm >= 100 ) &
          return
     !
     ! Change laws if a meteo file is used
     !
     if( kfl_model_chm == 2 .and. kfl_meteo_chm >= 2 ) then
        kfl_advec_chm = -1    ! VELOC_CHM:  advection
        lawde_chm     = -1    ! TEMPE_CHM: temperature
        lawte_chm     = -1    ! DENSI_CHM: densi1ty
     end if

     !
     ! Adjust stoichiometric mass ratio s = nuo Wo / nuF WF
     !
     do ispec = 1,nclas_chm
        if (stoichio_names(1) == speci(ispec)%name ) then ! Fuel
           stofu_chm(1) = ispec
        else if (stoichio_names(2) == speci(ispec)%name ) then ! Oxygen
           stofu_chm(2) = ispec
        else if (stoichio_names(3) == speci(ispec)%name ) then ! Oxygen
           stofu_chm(3) = ispec
        endif
     enddo

     if (kfl_stoic == 1_ip ) then
        strat_chm = 1.0_rp
        strat_chm = strat_chm / (speci(stofu_chm(1) )%weigh * stoic_chm(stofu_chm(1),stoichio_reaction,1) )
        strat_chm = strat_chm * (speci(stofu_chm(2))%weigh * stoic_chm(stofu_chm(2),stoichio_reaction,1) )
     endif

     ! Move diffusion coefficents to all the array, for redundancy
     if (kfl_model_chm == 4 .and. lawdi_chm(1,1) == 3_ip) then
        do ispec = 2,nclas_chm
           diffu_chm(1,ispec) = diffu_chm(1,1)
           diffu_chm(2,ispec) = diffu_chm(2,1)
        end do
     endif

     ! Linear dependence of density with temperature needs the a and b coefficients of the law
     if (kfl_model_chm == 4 .and. lawde_chm == -4 ) then
        ! Coeffs are a and b, where rho =a*T+b
        do ispec = 1,nclas_chm           
           !    a =   (rho_2 - rho_1) / (T_2 - T_1)
           speci(ispec)%densi(1) = (Linear_temp_density_coeffs(ispec,4)-Linear_temp_density_coeffs(ispec,2))/ &
                             (Linear_temp_density_coeffs(ispec,3)-Linear_temp_density_coeffs(ispec,1))
           !    b = rho_1 - a  T_1
           speci(ispec)%densi(2) = Linear_temp_density_coeffs(ispec,2) & 
                             - speci(ispec)%densi(1) * Linear_temp_density_coeffs(ispec,1)
           ! and the same for the second phase
           speci(ispec)%densi(3) = (Linear_temp_density_coeffs(ispec,8)-Linear_temp_density_coeffs(ispec,6))/ &
                             (Linear_temp_density_coeffs(ispec,7)-Linear_temp_density_coeffs(ispec,5))
           speci(ispec)%densi(4) = Linear_temp_density_coeffs(ispec,6) & 
                             - speci(ispec)%densi(3) * Linear_temp_density_coeffs(ispec,5)

        enddo
     endif

     if (kfl_model_chm==4) deallocate (Linear_temp_density_coeffs)

     if (kfl_model_chm==4 .and. kfl_norma_chm .ne. 0_ip)  call livinf(-9_ip, & 
          'Will use species '//trim(speci(kfl_norma_chm)%name)//' for correction, #',kfl_norma_chm)

     if (kfl_model_chm == 5 .and. check_cfi /= 1_ip) call runend('CHEMIC REAPHY: reaction,table option needs to be used for the CFI combustion model')

  end if
end subroutine chm_reaphy
