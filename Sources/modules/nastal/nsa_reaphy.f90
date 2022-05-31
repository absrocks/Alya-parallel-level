!------------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_reaphy.f90
!> @date    01/08/2012
!> @author  Mariano Vazquez
!> @brief   Input subroutine for the Physical properties
!> @details Input subroutine for the Physical properties
!> @}
!------------------------------------------------------------------------
subroutine nsa_reaphy
!-----------------------------------------------------------------------
  use      def_parame
  use      def_inpout
  use      def_master
  use      def_domain
  use      mod_memchk
  use      def_kermod
  use      mod_ker_proper 
  use      def_nastal
  use mod_ecoute, only :  ecoute
  use mod_messages, only : livinf

  implicit none
  integer(ip) :: ielem,inode,ielty,pnode,ipoin,kvinu,idime
  real(rp)    :: xsoun,tinfi,axyra,axzra,ayzra,dummr,rgasc
  character(300)           :: messa

  if( INOTSLAVE ) then
     !
     ! Initial and default values
     !
     kfl_stead_nsa = 0                                    ! Steady-state
     kfl_timei_nsa = 1                                    ! Time integration on (ALWAYS in nastal)
     kfl_advec_nsa = 1                                    ! Convection is on
     kfl_visco_nsa = 1                                    ! Viscous term on
     kfl_confi_nsa = 0                                    ! low is not confined
     kfl_locti_nsa = 0                                    ! Global time step (the smallest in the domain)
     kfl_inifi_nsa = 0                                    ! No initial fields given
     kfl_inibu_nsa = 0                                    ! No built-in initial fields given
     kfl_inkee_nsa = 0                                    ! No initial fields given

     kfl_turbu_nsa = kfl_modul(4)                         ! Cross coupling with TURBUL module

     kfl_inico_nsa = 1                                    ! Fluid motion initial condition 

     kfl_zero_initial_velocity_nsa = 0                    ! The problem has non zero initial velocity

     kfl_infun_nsa = 0                                    ! No special initial fields
     kfl_benme_nsa = 0                                    ! METEO: no particular benchmark problem
     kfl_spher_nsa = 1                                    ! METEO: pseudo-3D warm-bubble as a cylinder and not a sphere when kfl_spher_nsa = 0
     kfl_botta_nsa = 0                                    ! METEO: Botta and Klein equilibrium case
     kfl_ansou_nsa = -1                           ! When positive and benme='KESSLER', the sounding is analutic
     kfl_physics_nsa = -1                                  ! METEO: Flag to turn on or off the mirophysics
     kfl_sponge_nsa = -1                                    ! Sponge: by default it is OFF
     kfl_sptyp_nsa = -1                                    ! Sponge type: by default it is OFF
                                                          ! The other option is that of Klemp Lilly 1978 
     kfl_inleb_nsa = 0                                    ! Inlet boundary conditions fixed-fixed
     kfl_dampf_nsa = 0                                    ! Van Driest near-wall damping function for Smagorisnky model
     kfl_mfrco_nsa = 0                                    ! Mass flow rate control activation flag (=0 OFF, =1 ON)
     kfl_foreg_nsa = 0                                    ! Compressible flow
     kfl_relat_nsa = 0                                    ! Non-relativistic flow
     kfl_isent_nsa = 0                                    ! Non-Isentropic flow
     kfl_ncons_nsa = 0                                    ! Flag to use or not Non-conservative set
     kfl_pertu_nsa = 0                                    ! Not-perturbation variables (by default)
     kfl_theta_nsa = 0                                    ! using the usual transport theta equation by default when kfl_thetacons_nsa = 0

     kfl_coupl_nsa = 0                                    ! No coupling 
     kfl_cotur_nsa = 0                                    ! No eddy viscosity closure
     
     lawst_nsa = 1                                        ! Ideal gas
     adgam_nsa = 1.4_rp                                   ! Adiabatic exponent gamma (Cp / Cv)
     cpcoe_nsa = 1.4_rp                                   ! Specific heat at constant volume Cp
     cvcoe_nsa = 1.0_rp                                   ! Specific heat at constant volume Cv
     rgasc_nsa = 0.4_rp
     press_nsa = 0.0_rp
     press_dynamic_nsa = -1.0_rp
     lawvi_nsa = 0                                        ! Constant molecular viscosity
     mfrse_nsa = 1                                        ! Set from which the mass flow rate is calculated (by default = 1)
     kdiff_nsa = 0.0_rp                                   ! Viscosity coefficient of artificial diffusion
     kfl_adiff_nsa= 0                                     ! Flag for artificial diffusion
     kfl_hysta_nsa= 0                                     ! No hydrostatic correction
     kfl_brunt_nsa= 0                                     ! No Brunt frequency (meteo)
     kfl_rearst_nsa=0                                     ! Flag to read local restart file
     
     kfl_refpdt_nsa=0                                     ! Reference p-d-t: compute pressure from t and rho

     kfl_iposi_nsa = 0                                    ! Positional initial fields

     xbubb_nsa = 0.0_rp
     teano_nsa = 0.0_rp

     grnor_nsa     = 0.0_rp                               ! Gravity norm
     gravi_nsa     = 0.0_rp                               ! Gravity vector
     fcons_nsa     = 0.0_rp                               ! Non conservative form default
     fvins_nsa     = 1.0_rp                               ! Divergence form default
     densi_nsa     = 1.0_rp
     mfrgr_nsa     = 1.0_rp                               ! Mass flow rate growth parameter (by default no growth)
     mfrub_nsa     = 1.0_rp                               ! Ub target
     ubpre_nsa     = 0.0_rp                               ! Ub target previous time-step
     tempe_nsa     = 1.0_rp
     entro_nsa     = 1.0_rp
     veloc_nsa     = 0.0_rp
     speed_nsa     = 1.0_rp
     axyin_nsa     = 0.0_rp
     axzin_nsa     = 0.0_rp
     ayzin_nsa     = 0.0_rp
     veloc_nsa(1)  = 1.0_rp
     visco_nsa     = 0.0_rp
     thdif_nsa     = 0.0_rp
     poros_nsa     = 0.0_rp
     dtinv_nsa     = 1.0_rp                               ! 1/dt=1.0
     nmate_nsa     = 1                                    ! Number of materials
     mowei_nsa     = 20.78616_rp                          ! Molecular weight to get R_gas = 0.4 (old default value)
     runiv_nsa     = 8.314621_rp                          ! Universal gas constant Ro [J/K/mol] 

     rmach_nsa     = 0.0_rp
     rreyn_nsa     = 0.0_rp
     prand_nsa     = 0.0_rp
     pratu_nsa     = 0.0_rp
     cppra_nsa     = 0.0_rp
     cpprt_nsa     = 0.0_rp

     trefa_nsa     = 1.0_rp

     tinfi         = -2.0_rp

     axyin_nsa     = 0.0_rp
     axzin_nsa     = 0.0_rp
     spein_nsa     = 1.0_rp

     cleng_nsa     = 1.0_rp

     pbaro_nsa     = 100000_rp
     
     !Sponge and domain parameters from user
     xmin_nsa      = 0.0_rp
     xmax_nsa      = 0.0_rp
     ymin_nsa      = 0.0_rp
     ymax_nsa      = 0.0_rp
     zmin_nsa      = 0.0_rp
     zmax_nsa      = 0.0_rp
     xradi_nsa     = -1.0_rp
     yradi_nsa     = -1.0_rp
     zradi_nsa     = -1.0_rp
     rc_nsa        = -1.0_rp
     xc_nsa        = 0.0_rp
     yc_nsa        = 0.0_rp
     zc_nsa        = 0.0_rp
     dxs_nsa       = 0.0_rp
     dys_nsa       = 0.0_rp
     dzs_nsa       = 0.0_rp
     ampx_nsa      = 0.0_rp
     ampy_nsa      = 0.0_rp
     ampz_nsa      = 0.0_rp
     thetac_nsa    = 0.0_rp
     tracerc_nsa   = 0.0_rp

     kfl_rc_nsa    = 0_ip
     kfl_xc_nsa    = 0_ip
     kfl_yc_nsa    = 0_ip
     kfl_zc_nsa    = 0_ip
     kfl_thetac_nsa= 0_ip

     !------------------------------------------------------------------------
     ! Additional domain variables used only locally with Kessler (SM):
     !------------------------------------------------------------------------
     kfl_adiff_nsa  = 0             ! Flag for artificial diffusion

     istep_nsa      = 0
     nvar_nsa       = 4             ! number of variables (dynamics+tracers)
     nelx_nsa       = 1             ! number of elements in x
     nely_nsa       = 1             ! number of elements in y !used only if ndime > 2. Otherwise z is used for the vertical.
     nelz_nsa       = 1             ! number of elements in z
     nx_nsa         = nelx_nsa + 1  ! number of nodes in x
     ny_nsa         = nely_nsa + 1  ! number of nodes in y !used only if ndime > 2. Otherwise z is used for the vertical.
     nz_nsa         = nelz_nsa + 1  ! number of nodes in z
     ncol_nsa       = nx_nsa        ! number of columns (used for Kessler).
     if(ndime > 2) &
          ncol_nsa  = nx_nsa*ny_nsa ! number of columns (used for Kessler).

     kfl_uniformvelocity_nsa = -999

     if (timef < zensa) timef = 1000000000.0_rp

     !
     ! Reach the section
     !
     rewind(lisda)
     call ecoute('nsa_reaphy')
     do while(words(1)/='PHYSI')
        call ecoute('nsa_reaphy')
     end do
     !
     ! Begin to read data
     !
     !-----------------------------------------------------------------------
     ! ADOC[0,d]> $--tototoooooooo
     ! ADOC[0]> $-----------------------------------
     ! ADOC[0]> $-- Physical properties definition
     ! ADOC[0]> $-----------------------------------
     ! ADOC[0]> PHYSICAL_PROBLEM
     !-----------------------------------------------------------------------
     messa = &
          '        READING PHYSICS...'
     call livinf(0_ip,messa,one)

     do while(words(1)/='ENDPH')
        call ecoute('nsa_reaphy')
        if(words(1)=='PROBL') then
           !-----------------------------------------------------------------------
           ! ADOC[1]> PROBLEM_DEFINITION
           !-----------------------------------------------------------------------
     ! ADOC[1,d]> $--tototooooooooooooooo
     ! ADOC[d]> $--tototooooooooooooooo
     ! ADOC[d]> $--TOTOOOOO
           call ecoute('nsa_reaphy')
           do while(words(1)/='ENDPR')
              if(words(1)=='FORCE') then                         ! Force regime
                 !-----------------------------------------------------------------------
                 ! ADOC[2]> $--- Only when isentropic or perturbation variables (meteo) 
                 ! ADOC[2]> $--- When ISENTROPIC, isentropic constant r must be given  
                 ! ADOC[2]> [ FORCED_REGIME [ISENTROPIC [, VALUE=r] | PERTU] ]
                 !-----------------------------------------------------------------------
                 if(exists('ISENT')) kfl_isent_nsa = 1       ! Isentropic flow
                 if(exists('PERTU')) kfl_pertu_nsa = 1       ! Perturbation variables                      
                 if (kfl_isent_nsa == 1) then
                    entro_nsa = getrea('VALUE',-1.0_rp,'#Isentropic constant A')
                    if (entro_nsa < 0.0_rp) &
                         call runend ('NSA_REAPHY: ISENTROPIC CONSTANT MUST BE GIVEN IN ISENTROPIC MODELS!')
                 end if
              else if(words(1)=='ADVEC') then        
                 ! ADOC[2]> $--- Force advection action
                 ! ADOC[2]> [ ADVECTION [OFF] ]
                 if((words(2)  =='OFF  ') .or. (words(2)  =='INACT')) kfl_advec_nsa=0
              else if(words(1)=='VISCO') then
                 ! ADOC[2]> $--- Force viscosity action
                 ! ADOC[2]> [ VISCOSITY [OFF] ]
                 if((words(2)  =='OFF  ') .or. (words(2)  =='INACT')) then
                    kfl_visco_nsa=0
                    visco_nsa=0.0_rp
                 end if
              else if(words(1)=='GRAVI') then
                 ! ADOC[2]> $--- Gravity 
                 ! ADOC[2]> [ GRAVITY , MODUL= r , XCOMP= r , YCOMP= r , ZCOMP= r ]
                 ! back-compatibility:
 
                 grnor_nsa    = getrea('NORM ',0.0_rp,'#Gravity norm')
                 gravi_nsa(1) = getrea('GX   ',0.0_rp,'#x-component of g')
                 gravi_nsa(2) = getrea('GY   ',0.0_rp,'#y-component of g')
                 gravi_nsa(3) = getrea('GZ   ',0.0_rp,'#z-component of g')
                 ! preferred...
                 grnor_nsa    = getrea('MODUL',0.0_rp,'#Gravity norm')
                 gravi_nsa(1) = getrea('XCOMP',0.0_rp,'#x-component of g')
                 gravi_nsa(2) = getrea('YCOMP',0.0_rp,'#y-component of g')
                 gravi_nsa(3) = getrea('ZCOMP',0.0_rp,'#z-component of g')
                 call vecuni(3_ip,gravi_nsa,dummr)

                 gravm_nsa = 0.0_rp
                 do idime=1,ndime
                    gravm_nsa(idime,ndime+1) = grnor_nsa*gravi_nsa(idime)
                    gravm_nsa(ndime+2,idime) = grnor_nsa*gravi_nsa(idime)
                 end do

              else if(words(1)=='TURBU') then
                 !
                 ! ADOC[2]> TURBULENCE_MODEL:     LES_MODEL [SMAGORINSKY | WALE , PARAMETER=real] | FROM_TURBUL         $ Turbulence model
                 ! ADOC[d]> TURBULENCE_MODEL:
                 ! ADOC[d]> Turbulence Model definition.
                 ! ADOC[d]> For LES_MODEL, a real number is required. For smagorinsky, 
                 ! ADOC[d]> From Turbul options, Nastal uses the turbulent viscosity computed by Turbul.
                 !
                 if(words(2)=='LESMO') then
                    if(exists('SMAGO')) then
                       kfl_cotur_nsa=-1
                       if(exists('DAMPI')) then
                          kfl_dampf_nsa = 1
                       endif
                       turbu_nsa=getrea('PARAM',0.0_rp,'#Coefficient c**2')   ! The typical value is 0.1**2 = 0.01
                    else if(exists('WALE ')) then
                       kfl_cotur_nsa=-2
                       turbu_nsa=getrea('PARAM',0.0_rp,'#Coefficient c**2')   ! The default for WALE is 0.5**2 = 0.25
                    else if(exists('SIGMA')) then
                       kfl_cotur_nsa=-3
                       turbu_nsa=getrea('PARAM',0.0_rp,'#Coefficient c**2')   ! The default for SIGMA is 1.5**2 = 2.25
                    end if
                 else if(words(2)=='RANSD'.or.words(2)=='FROMT') then   ! We are using turmu
                    kfl_cotur_nsa = 1
                 end if

              else if(words(1)=='MASSF') then
                 kfl_mfrco_nsa = 1

                     mfrub_nsa = param(1)
                     mfrgr_nsa = param(2)
                     if (mfrgr_nsa > 2.0_rp) call runend('nsa_reaphy: The factor must be less than 2') 

                     mfrse_nsa = param(3)
 
              else if(words(1)=='TEMPO') then                    
                 if(exists('ON   ') .or. exists('ACTIV')) kfl_timei_nsa = 1 
              else if(words(1)=='INITI') then 
                 ! ADOC[2]> $--- Initial conditions type for the velocity
                 ! ADOC[2]> $---  Fluid at rest
                 ! ADOC[2]> $---  From velocity reference values
                 ! ADOC[2]> $---  From inflow velocity
                 ! ADOC[2]> [ INITIAL_CONDITIONS_TYPE  ATREST | REFERENCE | INFLOW_VELOCITY ]

                 if(words(2)  =='ATRES') kfl_inico_nsa=0         !   0. fluid at rest
                 if(words(2)  =='REFER') kfl_inico_nsa=1         !   1. from reference values
                 if(words(2)  =='INFLO') kfl_inico_nsa=2         !   2. from inflow velocity

                 !! THE FOLLOWING PROPERTIES ARE ALSO INCLUDED IN 
                 !! PROPERTIES DUE TO COMPATIBILITIES WITH NSTINC!!!!!!!!!

              else if(words(1)=='MACHN') then 
                 rmach_nsa = param(1)                       
              else if(words(1)=='REYNO') then               ! 
                 rreyn_nsa = param(1)                       
              else if(words(1)=='PRAND') then               ! 
                 prand_nsa =param(1)
              else if(words(1)=='TUPRA') then               ! 
                 pratu_nsa =param(1)
              else if(words(1)=='THETA') then               ! Theta equation in conservative from
                 if(words(2)  =='YES') &
                      kfl_theta_nsa=1
                  end if
              call ecoute('nsa_reaphy')
           end do
           !-----------------------------------------------------------------------
           ! ADOC[1]> END_PROBLEM_DEFINITION
           !-----------------------------------------------------------------------

        else if(words(1)=='PROPE') then
           !
           ! Properties
           !           
           !-----------------------------------------------------------------------
           ! ADOC[1]> PROPERTIES
           !-----------------------------------------------------------------------
           if(words(2)=='MATER') then
              nmate_nsa=getint('MATER',1_ip,'#Navier-Stokes number of materials')
              if(nmate_nsa>2) &
                   call runend('nsa_reaphy: ONLY 2 MATERIALS ARE POSSIBLE')
              call nsa_memphy
           end if
           
           kvinu = 0   ! when viscosity is given, mu is given

           call ecoute('nsa_reaphy')
           do while(words(1)/='ENDPR')
              if(words(1)=='MACHN') then                    
                 ! ADOC[2]> $--- Reference Mach number (SUPERSEDES Cp)
                 ! ADOC[2]> [ MACH_NUMBER ]
                 rmach_nsa = param(1)                       
              else if(words(1)=='REYNO') then               
                 ! ADOC[2]> $--- Reference Reynolds number (SUPERSEDES mu)
                 ! ADOC[2]> [ REYNOLDS_NUMBER ]
                 rreyn_nsa = param(1)                       
              else if(words(1)=='LENGT') then               ! 
                 ! ADOC[2]> $--- Reference characteristic length
                 ! ADOC[2]> [ LENGTH ]
                 cleng_nsa = param(1)                       
              else if(words(1)=='PRAND') then               ! Prandtl number (SUPERSEDES k)
                 ! ADOC[2]> $--- Prandtl number (SUPERSEDES k)
                 ! ADOC[2]> [ PRANDTL_NUMBER ]
                 prand_nsa =param(1)
              else if(words(1)=='TUPRA') then               ! Turbulent Prandtl number
                 ! ADOC[2]> $--- Turbulent Prandtl number
                 ! ADOC[2]> [ TU_PRANDTL_NUMBER ]
                 pratu_nsa =param(1)
              else if(words(1)=='DENSI') then               ! Density (rho, reference value)
                 densi_nsa=param(1)
              else if(words(1)=='TEMPE') then               ! Temperature (T, reference value)
                 tempe_nsa=param(1)
                 if (exists('DERIV')) kfl_refpdt_nsa = 1    !   Compute T from p and rho
              else if(words(1)=='INFTE') then               ! Temperature (T_inf at inflow, or at infinite)
                 tinfi    =param(1)                         !    it is used just to compute cp
              else if(words(1)=='SPEED') then               ! Module Speed (|u|, reference value)
                 speed_nsa=param(1)
              else if(words(1)=='MOLEC') then               ! Molecular weight
                 mowei_nsa=param(1)
              else if(words(1)=='FLUID') then               ! mks or cgs air under stp, water...
                 if (words(2)=='MKSAI') then
                    visco_nsa = 1.8e-5
                    if (kfl_visco_nsa == 0) visco_nsa = 0.0_rp
                    densi_nsa = 1.2_rp
                    mowei_nsa = 0.0289531_rp 
                    rgasc     = runiv_nsa / mowei_nsa
                 else if (words(2)=='STDAI') then       ! S.M. Standard Air. 
                    visco_nsa = 1.8e-5
                    if (kfl_visco_nsa == 0) visco_nsa = 0.0_rp
                    densi_nsa = 1.2_rp
                    mowei_nsa = 0.0289531_rp 
                    rgasc     = runiv_nsa / mowei_nsa
                 else 
                    call runend('NSA_REAPHY: FLUID MODEL NON-PRELOADED')
                 end if
              else if(words(1)=='STAND') then               ! set velocity to zero everywhere
                 kfl_inico_nsa= 0
                 speed_nsa = 0.0_rp
                 kfl_zero_initial_velocity_nsa = 1
              else if(words(1)=='ATTAC') then               ! Attack angle
                 attxy_nsa = getrea('XYANG',0.0_rp,'#Attack angle X-Y')
                 attxz_nsa = getrea('XZANG',0.0_rp,'#Attack angle X-Z')
                 attyz_nsa = getrea('YZANG',0.0_rp,'#Attack angle Y-Z')
                 axyra = attxy_nsa * acos(-1.0_rp) / 180.0_rp ! value in radians
                 axzra = attxz_nsa * acos(-1.0_rp) / 180.0_rp ! value in radians
                 ayzra = attyz_nsa * acos(-1.0_rp) / 180.0_rp ! value in radians
                 veloc_nsa(1) = speed_nsa * cos(axyra) * cos(axzra)
                 veloc_nsa(2) = speed_nsa * sin(axyra) * cos(axzra) * cos(ayzra)
                 veloc_nsa(3) = speed_nsa * sin(axzra) * sin(ayzra)
              else if(words(1)=='VELOC') then               ! Velocity vect. (u ref value, supersedes speed)
                 veloc_nsa(1:ndime)=param(1:ndime)
                 speed_nsa=sqrt(veloc_nsa(1)*veloc_nsa(1)+veloc_nsa(2)*veloc_nsa(2)+ &
                      veloc_nsa(3)*veloc_nsa(3))
              else if(words(1)=='PRESS') then               ! Pressure (p, reference value)
                 press_nsa=param(1)
              else if(words(1)=='PDYNA') then               ! Dynamic pressure (p, reference value)
                 press_dynamic_nsa=param(1)
              else if(words(1)=='VISCO') then               ! Viscosity (mu, reference value)
                 if (kfl_visco_nsa /= 0) then
                   visco_nsa=param(1)
                 endif
              else if(words(1)=='KINVI') then               ! Kinematic viscosity (nu, reference value)
                 visco_nsa=param(1)
                 kvinu = 1
              else if(words(1)=='DIFFU') then               ! Thermal diffusion (k, reference value)
                 thdif_nsa =param(1)
              else if(words(1)=='PERME') then               ! Porosity (sig)
                 poros_nsa(1:10)=param(1:10)
              else if(words(1)=='RECOV') then               ! Temperature recovery factor
                 trefa_nsa=param(1)

              else if(words(1)=='STATE') then               ! State law definition 

                 do while(words(1)/='ENDST')                 
                    if(words(1)         =='IDEAL') then
                       lawst_nsa = 1
!!                    else if(words(1)    =='ISOEN') then
!!                       lawst_nsa = 2
                    else if(words(1)    =='ISOWR') then
                       lawst_nsa = 3
                    else if(words(1)    =='POTEN') then
                       lawst_nsa = 4
                    else if(words(1)    =='HYDRO') then          ! hydrostatic correction is on
                       kfl_hysta_nsa = 1
                    else if(words(1)    =='CPCOE') then          ! State equation parameters
                       cpcoe_nsa = param(1)
                    else if (words(1)   =='CVCOE') then
                       cvcoe_nsa = param(1)
                    else if (words(1)   =='RGASC') then
                       rgasc_nsa = param(1)
                    else if (words(1)   =='GAMMA') then
                       adgam_nsa = param(1)
                    end if
                    call ecoute('nsa_reaphy')
                 end do

              else if(words(1)=='INLET') then               ! Inlet boundary conditions 
                 if(words(2)  =='REFER') then
                    kfl_inleb_nsa=1                         !   1. from reference values
                 else if(words(2)  =='SETTO') then
                    kfl_inleb_nsa=2                         !   2. from attack angles and speed
                    spein_nsa = getrea('SPEED',1.0_rp,'#Speed at inlets')
                    axyin_nsa = getrea('XYANG',0.0_rp,'#Attack angle X-Y at inlets')
                    axzin_nsa = getrea('XZANG',0.0_rp,'#Attack angle X-Z at inlets')
                    ayzin_nsa = getrea('YZANG',0.0_rp,'#Attack angle Y-Z at inlets')
                 end if
              else if(words(1)=='LAWDE') then               ! Law for rho
                 if(words(2)=='CONST') then
                    lawde_nsa=0
                 end if
              else if(words(1)=='LAWVI') then               ! Law for mu
                 if(words(2)     =='CONST') then
                    lawvi_nsa=0
                 else if(words(2)=='WIDER' .or. words(2)== 'POWER') then            ! Power law
                    lawvi_nsa=1_ip
                 else if(words(2)=='RESTR' .or. words(2)== 'SUTHE') then             ! Sutherland's law
                    lawvi_nsa=2_ip
                 end if
              else if(words(1)=='LAWPO') then               ! Law for sig
                 if(words(2)=='CONST') then
                    lawpo_nsa=0
                 end if
              else if (words(1) == 'BUILT') then
                 if (exists('LAXLI')) kfl_inibu_nsa= 1 
                 call ecoute('nsa_reaphy')
                 llcen_nsa(  1)= getrea('CENTX',1.0_rp,'#Center x')
                 llcen_nsa(  2)= getrea('CENTY',1.0_rp,'#Center y')
                 call ecoute('nsa_reaphy')
                 llrho_nsa(  1)= getrea('DENSI',1.0_rp,'#Density')
                 llvel_nsa(1,1)= getrea('VELOX',1.0_rp,'#Velocity x')
                 llvel_nsa(2,1)= getrea('VELOY',1.0_rp,'#Velocity y')
                 llpre_nsa(  1)= getrea('PRESS',1.0_rp,'#Pressure')
                 call ecoute('nsa_reaphy')
                 llrho_nsa(  2)= getrea('DENSI',1.0_rp,'#Density')
                 llvel_nsa(1,2)= getrea('VELOX',1.0_rp,'#Velocity x')
                 llvel_nsa(2,2)= getrea('VELOY',1.0_rp,'#Velocity y')
                 llpre_nsa(  2)= getrea('PRESS',1.0_rp,'#Pressure')
                 call ecoute('nsa_reaphy')
                 llrho_nsa(  3)= getrea('DENSI',1.0_rp,'#Density')
                 llvel_nsa(1,3)= getrea('VELOX',1.0_rp,'#Velocity x')
                 llvel_nsa(2,3)= getrea('VELOY',1.0_rp,'#Velocity y')
                 llpre_nsa(  3)= getrea('PRESS',1.0_rp,'#Pressure')
                 call ecoute('nsa_reaphy')
                 llrho_nsa(  4)= getrea('DENSI',1.0_rp,'#Density')
                 llvel_nsa(1,4)= getrea('VELOX',1.0_rp,'#Velocity x')
                 llvel_nsa(2,4)= getrea('VELOY',1.0_rp,'#Velocity y')
                 llpre_nsa(  4)= getrea('PRESS',1.0_rp,'#Pressure')

              else if(words(1)=='INITI') then
                 !
                 ! Initial Fields
                 !
                 if(words(2)=='CODES') then ! Initial Fields given by codes
                    call runend('NSA_REAPHY: "CODES" OPTION IN INITIAL_CONDITION IS NOT AVAILABLE')
                 else ! Initial Fields given on nodes
                    call ecoute('nsa_reaphy')
                    do while(words(1)/='ENDIN')
                       if (words(1) == 'VELOC') then
                          kfl_inifi_nsa(1) = 1
                          nfiel_nsa(1) = -getint('FIELD',1_ip,'#Field Number for momentum equation (velocity)')
                       else if (words(1) == 'PRESS') then
                          kfl_inifi_nsa(2) = 2
                          nfiel_nsa(2) = -getint('FIELD',1_ip,'#Field Number for continuity (pressure)')
                       else if (words(1) == 'DENSI') then
                          kfl_inifi_nsa(2) = 1
                          nfiel_nsa(2) = -getint('FIELD',1_ip,'#Field Number for continuity (density)')
                       else if (words(1) == 'TEMPE') then
                          kfl_inifi_nsa(3) = 1
                          nfiel_nsa(3) = -getint('FIELD',1_ip,'#Field Number for energy equation (temperature)')
                       else if (words(1) == 'ENERG') then
                          kfl_inifi_nsa(3) = 2
                          nfiel_nsa(3) = -getint('FIELD',1_ip,'#Field Number for energy equation (energy)')
                       else if (words(1) == 'POSIT') then
                          !
                          ! Positional initial fields
                          !
                          kfl_iposi_nsa = 1
                          if (exists('PRESS')) then
                             iposi_nsa(2)%kflag = 1
                             iposi_nsa(2)%geometry   = 1 ! sphere
                             iposi_nsa(2)%value      = getrea('VALUE',0.0_rp,'#x-center')
                             iposi_nsa(2)%center(1)  = getrea('XCENT',0.0_rp,'#x-center')
                             iposi_nsa(2)%center(2)  = getrea('YCENT',0.0_rp,'#y-center')
                             iposi_nsa(2)%center(3)  = getrea('ZCENT',0.0_rp,'#z-center')
                             iposi_nsa(2)%radius     = getrea('RADIU',0.0_rp,'#z-center')                             
                          end if
                       end if
                       call ecoute('nsa_reaphy')
                    end do
                    if (kfl_inifi_nsa(2) > 0) then
                       if (kfl_inifi_nsa(3) == 0) then
                          ! hay que ver aun...
                          call runend('NSA_REAPHY: PRESS/DENS AND TEMPERATURE/ENERGY MUST BE GIVEN')
                       end if
                    end if
                    if (kfl_inifi_nsa(3) > 0) then
                       if (kfl_inifi_nsa(2) == 0) then
                          ! hay que ver aun...
                          call runend('NSA_REAPHY: PRESS/DENS AND TEMPERATURE/ENERGY MUST BE GIVEN')
                       end if
                    end if
                    if ((kfl_inifi_nsa(1) + kfl_inifi_nsa(2) + kfl_inifi_nsa(3)) == 0) then
                       if (kfl_iposi_nsa == 0) call runend('NSA_REAPHY: NEITHER POSITIONAL NOR EXPLICIT FIELDS WERE GIVEN')
                    end if
                 end if

              else if(words(1)=='REFER') then               ! Reference fields
                 call ecoute('nsa_reaphy')
                 do while(words(1)/='ENDRE')
                    if (words(1) == 'HYDDE') then
                       kfl_inkee_nsa(2) = 1
                       nfiel_nsa(4) = -getint('FIELD',1_ip,'#Field Number for reference pressure')
                    else if (words(1) == 'HYDTE') then
                       kfl_inkee_nsa(3) = 1
                       nfiel_nsa(5) = -getint('FIELD',1_ip,'#Field Number for reference temperature')
                    else if (words(1) == 'HYDPR') then
                       kfl_inkee_nsa(4) = 1
                       nfiel_nsa(6) = -getint('FIELD',1_ip,'#Field Number for reference density')
                    end if
                    call ecoute('nsa_reaphy')
                 end do
                 if (kfl_inkee_nsa(2) == 0) then
                    call runend('NSA_REAPHY: NO REFERENCE DENSITY WAS GIVEN')
                 end if
                 if (kfl_inkee_nsa(3) == 0) then
                    call runend('NSA_REAPHY: NO REFERENCE TEMPERATURE WAS GIVEN')
                 end if
                 if (kfl_inkee_nsa(4) == 0) then
                    call runend('NSA_REAPHY: NO REFERENCE PRESSURE WAS GIVEN')
                 end if
                 
              else if(words(1)=='MATER') then               ! List of materials
                 call ecoute('nsa_reaphy')
                 do while(words(1)/='ENDMA')
                    ielem=int(param(1))
                    if(int(param(2))/=1) lmate_nsa(ielem)=-1
                    call ecoute('nsa_reaphy')
                 end do
              end if
              call ecoute('nsa_reaphy')
           end do
           
           !-----------------------------------------------------------------------
           ! ADOC[1]> END_PROPERTIES
           !-----------------------------------------------------------------------
           
        else if (words(1)=='COUPL') then
           call ecoute('nsa_reaphy')
           do while(words(1)/='ENDCO')
              if (words(1)=='MOIST') then
                 kfl_coupl_nsa = 1
                 lacou_nsa(1) = getrea('LATEN',0.0_rp,'#Reference for vaporization latent heat')                 
              end if
           end do
        else if (words(1)=='METEO') then
           call ecoute('nsa_reaphy')
           do while(words(1)/='ENDME')
              if (words(1)=='INITI') then
                 if (words(2)=='ISOTH') then
                    kfl_infun_nsa = 1            ! isothermal athmosphere (i.e. p barometric law)
                 else if (words(2)=='NEUTR') then
                    kfl_infun_nsa = 2            ! neutral or adiabatic athmosphere
                 else if (words(2)=='BRUNT') then
                    kfl_infun_nsa = 3            ! brunt - vaissala
                    if (exists('CONST')) then
                       brure_nsa = getrea('CONST',1.0_rp,'#Reference Brunt Vaissala frequency') 
                       kfl_brunt_nsa = 1
                    else if (exists('FIELD')) then
                       call runend('NSA_REAPHY: BRUNT FIELD TO BE DONE')
                       kfl_brunt_nsa = 2
                    end if
                 end if
              else if (words(1)=='REFER') then
                 pbaro_nsa= getrea('REFER',100000.0_rp,'#Reference pressure at ground level') 
              else if (words(1)=='ANOMA') then
                 teano_nsa(1:5)= param(1:5)
              else if (words(1)=='BENCH') then
                 ! benchmark name
                 if (words(2)=='WARMB') then
                    kfl_benme_nsa = 1                  ! cosin warm bubble
                    if (words(3)=='ROBER') &
                       kfl_benme_nsa = 11              ! Gaussian warm bubble of Robert 1993
                    if (words(3)=='AHMAD') &
                       kfl_benme_nsa = 12              ! Gaussian warm bubble of Robert 1993
                    
                    if (exists('CYLIN')) &
                         kfl_spher_nsa = 0               !Cylinder instead of sphere in 3D
 
                 else if (words(2)=='DENSI') then
                    kfl_benme_nsa = 2                  ! density current
                 else if (words(2)=='HSLIN') then
                    kfl_benme_nsa = 3
                    
                    if (words(3)=='KL78 ' .or. words(3)=='KLEMP') &! Sponge type
                         kfl_sptyp_nsa = 2
                    if (words(3)=='BOTTA' .or. words(3)=='EQUIL') &! Botta & Klein equilibrium case
                         kfl_botta_nsa = 1
                    
                 else if (words(2)=='HSNLI') then
                    kfl_benme_nsa = 4  

                    if (words(3)=='KL78 ' .or. words(3)=='KLEMP') &! Sponge type
                         kfl_sptyp_nsa = 2! sponge type
                    
                 else if (words(2)=='NHLIN') then
                    kfl_benme_nsa = 5                 ! NON-Hydrostatic linear mountain
                    if (words(3)=='KL78 ' .or. words(3)=='KLEMP') &! Sponge type
                         kfl_sptyp_nsa = 2! sponge type
                    
                 else if (words(2)=='NHNLI') then
                    kfl_benme_nsa = 6                 ! NON-Hydrostatic non-linear mountain
                    if (words(3)=='KL78 ' .or. words(3)=='KLEMP') &! Sponge type
                         kfl_sptyp_nsa = 2! sponge type
                    
                 else if (words(2)=='SCHAR') then
                    kfl_benme_nsa = 7                 ! Schar mountain
                    if (words(3)=='KL78 ' .or. words(3)=='KLEMP') &! Sponge type
                         kfl_sptyp_nsa = 2! sponge type

                 else if (words(2)=='INERT') then
                    kfl_benme_nsa = 15                 ! inertia gravity waves
                    
                 else if (words(2)=='KESSL' .or. words(2)=='MOIST') then
                    kfl_benme_nsa = 200                 ! equilibrium mountain test of Botta and Klein
                    istep_nsa = 0
                    if (words(3)=='SIMPL' .or. words(3)=='SUPER') then
                       kfl_benme_nsa = 201                 ! equilibrium mountain test of Botta and Klein
                       istep_nsa = 0
                    end if

                 else if (words(2)=='KESSL' .or. words(2)=='MOIST') then
                    kfl_benme_nsa = 200                 ! Storm as in Sasa's paper
                    if (words(3)=='KL78 ' .or. words(4)=='KLEMP') &! Sponge type
                            kfl_sptyp_nsa = 2! sponge type
               
                    istep_nsa = 0
                    if (words(3)=='SIMPL' .or. words(3)=='SUPER') then
                       kfl_benme_nsa = 201              ! simple with superaturated cloud from the first step
                       istep_nsa = 0
                    else if (words(3)=='TEST ' .or. words(3)=='HSREF') then
                       kfl_benme_nsa = 202

                       istep_nsa = 0
                    
                    else if (words(3)=='KLAAS ' .or. words(3)=='KC85 ') then
                       kfl_benme_nsa = 203
                       
                       istep_nsa = 0
                       
                    else if (words(3)=='HSMOU' .or. words(3)=='MOUNT') then
                       kfl_benme_nsa = 204              ! Linear HS mountain with moisture

                       if (words(4)=='KL78 ' .or. words(4)=='KLEMP') &! Sponge type
                            kfl_sptyp_nsa = 2! sponge type

                       istep_nsa = 0
                       
                    else if (words(3)=='NHMOU') then
                       kfl_benme_nsa = 205              ! Linear HS mountain with moisture
                       
                       if (words(4)=='KL78 ' .or. words(4)=='KLEMP') &! Sponge type
                            kfl_sptyp_nsa = 2! sponge type
                       
                       istep_nsa = 0

                    else if (words(3)=='GRABO' .or. words(3)=='GK91 ' .or. words(3)=='G2007') then
                       kfl_benme_nsa = 210              !Grabowski JAS 2007
                       istep_nsa = 0
                    end if

                 end if
                 
                 
              else if (words(1)=='ANALY' .or. words(1)=='SNDAN') then
                 kfl_ansou_nsa = 1              ! Analytic sounding

              else if(words(1) == 'SPONG') then          ! Sponge
                 if(words(2) == 'ON   ') then
                    kfl_sponge_nsa = 1
                    kfl_sptyp_nsa = 1
                    if(words(3) == 'KL78 ' .or. words(3) == 'KLEMP') &
                         kfl_sptyp_nsa = 2
                 end if

              else if(words(1) == 'VELOC' .or. words(1) == 'UVELO') then          ! State equation parameters
                 if(words(2) == 'UNIFO') then
                    kfl_uniformvelocity_nsa = 1
                    uvelo_nsa = getrea('UNIFO',0.0_rp,'#Uniform velocity') 
                 end if
                 
              else if(words(1) == 'VELOC' .or. words(1) == 'VVELO') then          ! State equation parameters
                 if(words(2) == 'UNIFO') then
                    kfl_uniformvelocity_nsa = 1
                    vvelo_nsa = getrea('UNIFO',0.0_rp,'#Uniform velocity') 
                 end if

              else if(words(1) == 'VELOC' .or. words(1) == 'WVELO') then          ! State equation parameters
                 if(words(2) == 'UNIFO') then
                    kfl_uniformvelocity_nsa = 1
                    wvelo_nsa = getrea('UNIFO',0.0_rp,'#Uniform velocity') 
                 end if

              else if(words(1) == 'PHYSI' .or. words(1) == 'MICRO') then          ! State equation parameters
                 if(words(2) == 'ON   ') &
                      kfl_physics_nsa = 1

                 
              else if(words(1)    =='DXSPO') then          ! State equation parameters
                 !
                 ! Dx-Sponge given by the user:
                 !
                 dxs_nsa = param(1)

              else if(words(1)    =='DYSPO') then          ! State equation parameters
                 !
                 ! Dy-Sponge given by the user:
                 !
                 dys_nsa = param(1)
                 
              else if(words(1)    =='DZSPO') then          ! State equation parameters
                 !
                 ! Dz-Sponge given by the user:
                 !
                 dzs_nsa = param(1)
              
              else if(words(1)    =='AMPX ') then          ! State equation parameters
                 !
                 ! Sponge amplitude in x
                 !
                 ampx_nsa = param(1)
              
              else if(words(1)    =='AMPY ') then          ! State equation parameters
                 !
                 ! Sponge amplitude in y
                 !
                 ampy_nsa = param(1)
                 
              else if(words(1)    =='AMPZ ') then          ! State equation parameters
                 !
                 ! Sponge amplitude in z
                 !
                 ampz_nsa = param(1)

              else if(words(1)    =='XMIN ') then          ! State equation parameters
                 !
                 ! Domain minimum x:
                 !
                 xmin_nsa = param(1)

              else if(words(1)    =='XMAX ') then          ! State equation parameters
                 !
                 ! Domain max x:
                 !
                 xmax_nsa = param(1)

              else if(words(1)    =='YMIN ') then          ! State equation parameters
                 !
                 ! Domain min y:
                 !
                 ymin_nsa = param(1)

              else if(words(1)    =='YMAX ') then          ! State equation parameters
                 !
                 ! Domain max y:
                 !
                 ymax_nsa = param(1)

              else if(words(1)    =='ZMIN ') then          ! State equation parameters
                 !
                 ! Domain min z:
                 !
                 zmin_nsa = param(1)

              else if(words(1)    =='ZMAX ') then          ! State equation parameters
                 !
                 ! Domain max z:
                 !
                 zmax_nsa = param(1)

              else if(words(1)    =='THETC') then          ! Amplitude of theta perturbation for the Kessler problems
                 !
                 ! Theta perturbation amplitude
                 !
                 thetac_nsa = param(1)

              else if(words(1)    =='XRADI') then          ! X Radius of bubble
                 !
                 ! X radius fo bubble
                 !
                 xradi_nsa = param(1)
                 
              else if(words(1)    =='YRADI') then          ! Y Radius of bubble
                 !
                 ! Y radius fo bubble
                 !
                 yradi_nsa = param(1)

              else if(words(1)    =='ZRADI') then          ! Z Radius of bubble
                 !
                 ! Z radius fo bubble
                 !
                 zradi_nsa = param(1)

              else if(words(1)    =='RC   ') then          ! Radius of bubble
                 !
                 ! rc radius for bubble (in case you want xr=yr=zr for a sphere)
                 !
                 kfl_rc_nsa = 1
                 rc_nsa = param(1)

              else if(words(1)    =='XC   ') then          ! Center x-coord of bubble
                 !
                 ! Xc of bubble
                 !
                 kfl_xc_nsa = 1
                 xc_nsa = param(1)

              else if(words(1)    =='YC   ') then          ! Center y-coord of bubble
                 !
                 ! Yc of bubble
                 !
                 yc_nsa = param(1)

              else if(words(1)    =='ZC   ') then          ! Center z-coord of bubble
                 !
                 ! Zc of bubble
                 !
                 zc_nsa = param(1)

              else if(words(1)    =='SPHER') then          ! Sphere (not cylinder)
                 !
                 ! Spherical perturbation (by default it is a cylinder in 3D, or a circle/ellipses in 2D)
                 !
                 kfl_cylind_nsa = -1
                 
              else if(words(1)    =='TRACC') then          ! Amplitude of tracer perturbation for the Kessler problems
                 !
                 ! Tracer perturbation amplitude
                 !
                 tracerc_nsa = param(1)

              else if (words(1)=='EQUAT') then        ! Equation set
                 if(words(2) == 'NONCO') kfl_ncons_nsa = 1       !Non-conservative set
             
              else if (words(1)=='RESTA') then        ! To read the local (meteo only) restart file
                 if (words(2)=='ON') &
                      kfl_rearst_nsa = 1
                 
              else if (words(1)=='NELZ ' .and. ndime < 3) then
                 nelz_nsa = param(1)
                 nelx_nsa = nelem/nelz_nsa
                 nz_nsa = nelz_nsa + 1_ip
                 nx_nsa = nelx_nsa + 1_ip
                 ncol_nsa = nx_nsa
                 
              else if (words(1)=='NELX ' .and. ndime < 3) then
                 nelx_nsa = param(1)
                 nelz_nsa = nelem/nelx_nsa
                 nz_nsa = nelz_nsa + 1_ip
                 nx_nsa = nelx_nsa + 1_ip
                 ncol_nsa = nx_nsa
                 
              else if (words(1)=='NELY ' .and. ndime < 3) then
                 nely_nsa = param(1)
                 nelz_nsa = nelem/nely_nsa
                 nz_nsa = nelz_nsa + 1_ip
                 ny_nsa = nely_nsa + 1_ip
                 ncol_nsa = nx_nsa

              else if (words(1)=='NELZ ' .and. ndime > 2) then
                 nelz_nsa = param(1)
                 nz_nsa = nelz_nsa + 1_ip

              else if (words(1)=='NELX ' .and. ndime > 2) then
                 nelx_nsa = param(1)
                 nx_nsa = nelx_nsa + 1_ip

              else if (words(1)=='NELY ' .and. ndime > 2) then
                 nely_nsa = param(1)
                 ny_nsa = nely_nsa + 1_ip
                 
                 ! object definition
              else if (words(1)=='YINTE') then   ! for helmholtz
                 xbubb_nsa(2,1)= param(1)                 
              else if (words(1)=='ZINTE') then
                 xbubb_nsa(3,1)= param(1)                 
              else if (words(1)=='XINTE') then
                 xbubb_nsa(1,1)= param(1)                 
              else if (words(1)=='ADIFF' .or. words(1)=='ARTIF') then
                 kfl_adiff_nsa = 1
                 kdiff_nsa = param(1)
              end if
              call ecoute('nsa_reaphy')
           end do
           
        end if
     end do
     !-----------------------------------------------------------------------
     ! ADOC[0]> END_PHYSICAL_PROBLEM
     !-----------------------------------------------------------------------

     !
     ! Checkouts
     !
     if (kvinu == 1) visco_nsa = visco_nsa * densi_nsa
     if (rreyn_nsa > zensa) then
        if (speed_nsa > zensa) then
           visco_nsa = densi_nsa * speed_nsa * cleng_nsa / rreyn_nsa
        else
           visco_nsa = densi_nsa * cleng_nsa / rreyn_nsa
        end if
     end if

     if (tinfi < -1.0_rp) tinfi = tempe_nsa

     if (rmach_nsa > zensa) then
        adgam_nsa = 1.4_rp
        cvcoe_nsa= speed_nsa*speed_nsa / &
             (adgam_nsa*(adgam_nsa-1.0_rp)*rmach_nsa*rmach_nsa*tinfi)
        cpcoe_nsa= cvcoe_nsa*adgam_nsa
        if (cpcoe_nsa < zensa) then
           cpcoe_nsa = 1.4_rp 
           cvcoe_nsa= cpcoe_nsa / adgam_nsa
        end if
        rgasc     = cpcoe_nsa - cvcoe_nsa   
        mowei_nsa = runiv_nsa / rgasc 
        if (kfl_prope /= 0) then
           print*, ' '
           print*, '     ERROR in *.nsa.dat '
           print*, '     You cannot run this problem with kernel properties '
           print*, '     you need to remove the properties from *.ker.dat'
           print*, '     Sorry :('
           print*, ' '
           print*, '     Correct your *.ker.dat file and re-execute the code'
           print*, ' '
           stop
        endif
     else
        ! Only for initialization purposes: prescribe gamma 'adgam_nsa' & molecular_weight 'mowei_nsa'
        adgam_nsa = 1.4_rp
        !
        rgasc     = runiv_nsa / mowei_nsa
        cvcoe_nsa = rgasc / (adgam_nsa - 1.0_rp) 
        cpcoe_nsa = adgam_nsa * cvcoe_nsa
     end if
     
     if (prand_nsa > zensa) then
        thdif_nsa = visco_nsa * cpcoe_nsa / prand_nsa 
     else if (thdif_nsa > zensa) then        
        prand_nsa = visco_nsa * cpcoe_nsa / thdif_nsa
     else if (thdif_nsa < zensa) then        
        prand_nsa = 0.8_rp                                 ! ... which is very reasonable for gases
        thdif_nsa = visco_nsa * cpcoe_nsa / prand_nsa 
     end if
     
     if (prand_nsa > zensa) cppra_nsa = cpcoe_nsa / prand_nsa
     
     if (pratu_nsa > zensa) then
        cpprt_nsa = cpcoe_nsa / pratu_nsa
     else if (visco_nsa > zensa) then        
        pratu_nsa = 1.0_rp                                 ! ... which is very reasonable always
        cpprt_nsa = cpcoe_nsa / pratu_nsa
     end if
     
     ! for diatomic ideal gases
     !        cpcoe_nsa = 1005.095_rp
     !        cvcoe_nsa = 717.925_rp
     !        rgasc_nsa = 287.17_rp
     !        adgam_nsa = 1.4_rp
     
     if(kfl_infun_nsa /= 0) then
        !
        ! Standard atmosphere (air): DO NOT ERASE, this is used for meteo problems!
        !
        !(SM MM)
        cpcoe_nsa = 1004.67_rp
        cvcoe_nsa = 717.5_rp
        rgasc     = 287.17_rp
        adgam_nsa = 1.4_rp
     end if
     
     ! state law setup: compute model constants
     call nsa_stalaw(0,0,densi_nsa,press_nsa,tempe_nsa,dummr,dummr,mowei_nsa,cpcoe_nsa)
     
     if (press_nsa < zensa) then
        messa = &
             '        REFERENCE PRESSURE COMPUTED FROM GIVEN TEMPERATURE AND DENSITY'
        call livinf(0_ip,messa,one)
        call nsa_stalaw(2,0,densi_nsa,press_nsa,tempe_nsa,dummr,dummr,mowei_nsa,cpcoe_nsa)
     else
        if (kfl_refpdt_nsa == 1) then
           messa = &
                '        REFERENCE TEMPERATURE COMPUTED FROM GIVEN PRESSURE AND DENSITY'
           call livinf(0_ip,messa,one)
           call nsa_stalaw(3,0,densi_nsa,press_nsa,tempe_nsa,dummr,dummr,mowei_nsa,cpcoe_nsa)           
        else
           messa = &
                '        REFERENCE PRESSURE COMPUTED FROM GIVEN TEMPERATURE AND DENSITY'
           call livinf(0_ip,messa,one)
           call nsa_stalaw(2,0,densi_nsa,press_nsa,tempe_nsa,dummr,dummr,mowei_nsa,cpcoe_nsa)              
        end if
     end if
     
     if (press_dynamic_nsa < zensa) then
        press_dynamic_nsa= 0.5_rp * densi_nsa * speed_nsa * speed_nsa
     end if
     
     if (kfl_isent_nsa == 1) &
          call nsa_stalaw(3,0,densi_nsa,press_nsa,tempe_nsa,dummr,dummr,mowei_nsa,cpcoe_nsa)
     
     xsoun_nsa    = sqrt(adgam_nsa * rgasc * tinfi)
     rmach_nsa    = speed_nsa / xsoun_nsa
     
     tstag_nsa = tinfi  * &
          (1.0_rp + trefa_nsa*(adgam_nsa-1.0_rp)*rmach_nsa*rmach_nsa*0.5_rp)
     
     rgasc_nsa = rgasc  ! rgasc is a kernel constant, supersedes rgasc_nsa

     
     if (lawvi_nsa == 1_ip) then  !power viscosity law
        vispa_nsa(2) = 0.76_rp 
        vispa_nsa(1) = visco_nsa / tempe_nsa**vispa_nsa(2)
     else  if (lawvi_nsa == 2_ip) then    ! sutherland's viscosity law (parameters for air)
        vispa_nsa(2) = 110.4_rp 
        vispa_nsa(1) = (tempe_nsa + vispa_nsa(2)) * visco_nsa / tempe_nsa**1.5_rp 
     end if
     
     if (speed_nsa < zensa) then
        kfl_zero_initial_velocity_nsa = 1
     end if
     
     ivert_nsa = 1    ! default value to one, always
     do idime= 1,ndime
        if ((gravi_nsa(idime) > 1.0e-10) &
             .or.(gravi_nsa(idime) < -1.0e-10)) then
           ivert_nsa= idime
        end if
     end do
     
     if (visco_nsa < zensa) then
        kfl_visco_nsa = 0                 ! no viscosity -> no viscous terms
     else
!        kfl_visco_nsa = 1
     end if

     !
     ! Update nodal materials: 1=fluid, -1=solid, 0=interface 
     !     
     if(nmate_nsa>1) then
        lmatn_nsa=-1
        do ielem=1,nelem
           if(lmate_nsa(ielem)==1) then
              ielty=ltype(ielem)
              pnode=nnode(ielty)
              do inode=1,pnode
                 ipoin=lnods(inode,ielem)
                 lmatn_nsa(ipoin)=1
              end do
           end if
        end do
        do ielem=1,nelem
           if(lmate_nsa(ielem)==-1) then
              ielty=ltype(ielem)
              pnode=nnode(ielty)
              do inode=1,pnode
                 ipoin=lnods(inode,ielem)
                 if(lmatn_nsa(ipoin)==1) lmatn_nsa(ipoin)=0
              end do
           end if
        end do
     end if

     sofac_nsa = sqrt(adgam_nsa*rgasc)
     rgacv_nsa = rgasc / cvcoe_nsa
     
     rgcou_nsa    = 0.0_rp
     cpcou_nsa    = 0.0_rp
     cvcou_nsa    = 0.0_rp
     if (kfl_coupl_nsa == 1) then
        cvcou_nsa(2) = 1424.0_rp   ! water vapor cv     (COUPLING WITH CHEMIC)
        cpcou_nsa(4) = 4186.0_rp   ! liquid water cp    (COUPLING WITH CHEMIC)
        rgcou_nsa(2) =  461.0_rp   ! water vapor rgasc  (COUPLING WITH CHEMIC)
     else
        nspec = 0                  ! no species (WHEN NOT COUPLED WITH CHEMIC)
     end if


     !
     ! Count the columns for 3D problems: (SM)
     !
     if (ndime > 2 .and. kfl_benme_nsa >= 200) then
        if(nelx_nsa == 0 .or. nely_nsa == 0) then
           print*, ' '
           print*, '      ERROR in *.nsa.dat '
           print*, '      For 3D problems with moisture '
           print*, '      you need to give all three values for:'
           print*, '      NELX, NELY, and NELZ'
           print*, ' '
           print*, '      Correct your *.nsa.dat file and re-execute the code'
           print*, ' '
           stop

        else
           print*,'-| ALYA             nnodex   = ', nx_nsa 
           print*,'-| ALYA             nnodey   = ', ny_nsa 
           print*,'-| ALYA             nnodez   = ', nz_nsa 

           ncol_nsa = nx_nsa*ny_nsa
        end if

     else if (ndime < 3 .and. kfl_benme_nsa >= 200) then
        if(nelz_nsa == 0) then
           print*, ' '
           print*, '     ERROR in *.nsa.dat '
           print*, '     For 2D problems with moisture '
           print*, '     you need to give the value of vertical layers:'
           print*, '     NELZ'
           print*, ' '
           print*, '     Correct your *.nsa.dat file and re-execute the code'
           print*, ' '
           stop

        else
           print*,'-| ALYA             nnodex   = ', nx_nsa 
           print*,'-| ALYA             nnodez   = ', nz_nsa 
        end if
        
     end if
!!     print*,'-| ALYA             ncolumns = ', ncol_nsa


  end if


end subroutine nsa_reaphy
