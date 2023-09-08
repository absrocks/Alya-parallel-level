subroutine tur_reaphy()
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_reaphy
  ! NAME 
  !    tur_reaphy
  ! DESCRIPTION
  !    This routine reads the physical treatment 
  ! USED BY
  !    tur_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_inpout
  use def_master
  use def_turbul
  use def_domain
  use mod_memchk
  use def_kermod, only     :  cmu_st, kfl_kemod_ker, kfl_logva
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip) :: kfl_atmbl, imate
  real(rp)    :: dummr
  logical     :: rwind ! realizable using wind coeffficients
 
  if( INOTSLAVE ) then
     !
     ! Initializations (defaults)
     !
     kfl_model_tur = 0                                   ! No model
     kfl_timei_tur = 0                                   ! Stationary flow
     kfl_cotem_tur = 0                                   ! No temperature coupling
     kfl_advec_tur = 1                                   ! Advection
     kfl_colev_tur = 0                                   ! No coupling with LEVELS
     kfl_ddesm_tur = 0                                   ! No model DDES 
     kfl_sasim_tur = 0                                   ! No SAS model
     kfl_rotat_tur = 0                                   ! No rotating reference frame
     kfl_inifi_tur(1) = 0                                ! Initial fields
     kfl_inifi_tur(2) = 0                                ! Initial fields
     kfl_inifi_tur(3) = 0                                ! Initial fields
     inits_tur     = 0                                   ! Initial time step
     nturb_tur     = 1                                   ! # turbulence variables
     ipara_tur     = 0                                   ! Integer Parameters
     lawde_tur     = 0                                   ! Law for rho
     lawvi_tur     = 0                                   ! Law for mu 
     kfl_kemod     = 0                                   ! k eps modification
     kfl_logva_tur     = 0                                   ! works with logarithm of unknowns
     boube_tur     = 0.0_rp                              ! Beta
     grnor_tur     = 0.0_rp                              ! Gravity norm |g|
     gravi_tur     = 0.0_rp                              ! Gravity vector g
     prtur_tur     = 1.0_rp                              ! Turbulent Prandtl number
     param_tur     = 0.0_rp                              ! Real parameters
     densi_tur     = 0.0_rp                              ! Density
     visco_tur     = 0.0_rp                              ! Viscosity
     densa_tur     = 1.0_rp                              ! Air density
     visca_tur     = 1.0_rp                              ! Air viscosity
     cddes_tur     = 0.65_rp                             ! CDDES - smagorinsky-like constant for DES
     inv_l_max          = 0.0_rp                         ! Mixing length (k-eps, CFDWind2)
     kfl_discd_tur      = 0
     ldiss_material_tur = 0
     !
     ! Local variables
     !
     kfl_atmbl     = 0                                   ! Atmospheric boundary layer model
     rwind         =.false.                              ! k eps realizable model using wind coefficients
     !
     ! Reach the section
     !
     rewind(lisda)
     call ecoute('tur_reaphy')
     do while(words(1)/='PHYSI')
        call ecoute('tur_reaphy')
     end do
     !
     ! Begin to read data
     ! 
     do while(words(1)/='ENDPH')
        call ecoute('tur_reaphy')

        if(words(1)=='PROBL') then
           !
           ! Problem definition data
           !
           call ecoute('tur_reaphy')
           do while(words(1)/='ENDPR')

              if(words(1)=='MODEL') then
                 if(words(2)=='SPALA') then
                    kfl_model_tur =  1
                 else if (words(2)=='KXUCH') then
                    kfl_model_tur =  2
                 else if (words(2)=='STDKE') then
                    kfl_model_tur = 11
                    if( exists('ABL  ') ) kfl_atmbl = 1 ! Neutral ABL
!                    if( exists('NEUTR') ) kfl_atmbl = 1 ! Neutral ABL
!                    if( exists('NONNE') ) kfl_atmbl = 2 ! Non-neutral ABL
                    if (exists('RNG  ') ) kfl_kemod = 1
                    if (exists('REALI') ) then 
                       kfl_kemod = 2
                       if (exists('WIND ') )  rwind=.true.
                    end if
                    if (exists('KEFP ') ) kfl_kemod = 3 !  An improved k - ε model applied to a wind turbine wake in
                    !              atmospheric turbulence, Wind Energ. 2013; 00:1–19 DOI: 10.1002/we
                    if (exists('LOGVA') ) kfl_logva_tur = 1 ! works with logarithmic functions 
                    ! copy variable to kernel 
                    if (kfl_logva_tur==1) kfl_logva=1
                    if (kfl_logva==1)     kfl_logva_tur=1
                    kfl_kemod_ker = kfl_kemod  ! copy variable to kernel 
                 else if (words(2)=='LAUND') then
                    kfl_model_tur = 12
                 else if (words(2)=='CHIEN') then
                    kfl_model_tur = 13
                 else if (words(2)=='LAMBR') then
                    kfl_model_tur = 14
                 else if (words(2)=='RODIT') then
                    kfl_model_tur = 15
                 else if (words(2)=='NATXU') then 
                    kfl_model_tur = 16
                 else if (words(2)=='KEPSV') then
                    kfl_model_tur = 17
                 else if (words(2)=='KEPSJ'.or.words(2)=='JAWHW') then
                    kfl_model_tur = 18
                 else if (words(2)=='KEPSN'.or.words(2)=='NAGAN') then
                    kfl_model_tur = 19
                 else if (words(2)=='KEPSP') then
                    kfl_model_tur = 20
                 else if (words(2)=='KOMEG') then
                    kfl_model_tur = 30
                    if( exists('MODIF') ) ipara_tur(1) = 1
                 else if (words(2)=='BREDB') then
                    kfl_model_tur = 31
                 else if (words(2)=='SSTKO') then
                    kfl_model_tur = 32
                    if(exists('ROTAT'))  kfl_rotat_tur = 1! Rotating Reference frame
                 end if
                 call tur_inivar(2_ip)
                 if(exists('START')) then
                    inits_tur=getint('START',0_ip,'#Start at time step')
                 end if

              else if(words(1)=='TEMPO') then                 ! Temporal evolution
                 if(exists('ON   ')) kfl_timei_tur = 1 
              
              else if(words(1)=='DDES ') then                 ! DDES model flag
                 kfl_ddesm_tur = 1
                 if(exists('IMPRO ')) kfl_ddesm_tur = 2
                 cddes_tur     = getrea('CDDES ',0.65_rp,'#CDDES Smagorinski Constant')

              else if(words(1)=='SAS') then                 ! SAS model flag
                 if(exists('ON   ')) kfl_sasim_tur = 1

              else if(words(1)=='CONVE') then                 ! Advection
                 if(exists('NASTI')) kfl_advec_tur = 1 
                 if(exists('ON   ')) kfl_advec_tur = 1 
                 if(exists('OFF  ')) kfl_advec_tur = 0 

              else if(words(1)=='TEMPE'.and.words(2)=='BOUSS') then
                 kfl_cotem_tur = 1
                 boube_tur    = getrea('BETA ',0.0_rp,'#Beta coefficient')
                 grnor_tur    = getrea('G    ',0.0_rp,'#Gravity norm')
                 gravi_tur(1) = getrea('GX   ',0.0_rp,'#x-component of G')
                 gravi_tur(2) = getrea('GY   ',0.0_rp,'#y-component of G')
                 gravi_tur(3) = getrea('GZ   ',0.0_rp,'#z-component of G')
                 call vecuni(3_ip,gravi_tur,dummr)
                 if (exists('SOGAC')) kfl_cotem_tur = 2 ! Sogachev's model

              else if(words(1)=='LEVEL') then
                 if(words(2)=='ON   ') then
                    kfl_colev_tur=1
                    if(exists('THICK')) then
                       call runend('TUR_REAPHY: now thickness is read by level and stored in a thick that belongs to master')  
                    end if
                    if(exists('STAGG')) kfl_colev_tur=3
                 end if

              else if(words(1)=='GRAVI') then
                 grnor_tur     = getrea('NORM ',0.0_rp,'#Gravity norm')
                 gravi_tur(1)  = getrea('GX   ',0.0_rp,'#x-component of g')
                 gravi_tur(2)  = getrea('GY   ',0.0_rp,'#y-component of g')
                 gravi_tur(3)  = getrea('GZ   ',0.0_rp,'#z-component of g')
                 call vecuni(3_ip,gravi_tur,dummr)
              else if(words(1)=='DISCD') then ! extra dissipation term for actuator disc model
                 kfl_discd_tur= 1
               
                 call ecoute('tur_reaphy')
                 do while( words(1) /= 'ENDDI' )
                    if( words(1) == 'MATER' ) then
                       imate = getint('MATER',1_ip,'#Material force number')
                       if( imate < 1 .or. imate > nmate ) call runend('TUR_REAPHY: WRONG MATERIAL NUMBER')
                       ldiss_material_tur(imate) = 1 ! materials for dissipation
                    else 
                       call runend('TUR_REAPHY: WRONG DISSIPATION FIELD') 
                    end if
                    if (max_mater_tur.lt.nmate) call runend('TUR_REAPHY: too many materials for ldiss_material')
                    call ecoute('tur_reaphy')
                 end do
              end if
              call ecoute('tur_reaphy')
           end do

        else if(words(1)=='PROPE') then
           !
           ! Properties and model constants
           !
           do while(words(1)/='ENDPR')
              call ecoute('tur_reaphy')

              if(words(1)=='DENSI') then                      ! Density (rho)
                 densi_tur(1:10)=param(1:10)

              else if(words(1)=='VISCO') then                 ! Viscosity (mu)
                 visco_tur(1:10)=param(1:10)

              else if( words(1) == 'LMAXI' ) then
                 !
                 ! max. mixing length l_max
                 !
                 inv_l_max = getrea('LMAXI',0.0_rp,'#Max. mixing length (CFDWind model)')
                 ! The inverse of l_max  will be stored
                 if( inv_l_max < 1.0e8_rp.and.abs(inv_l_max)>1.0e-8_rp) then
                    inv_l_max= 1.0_rp/inv_l_max  
                 else if (inv_l_max >1.0e7_rp) then
                    inv_l_max= 0.0_rp
                 else
                    call runend ('tur_reaphy:invalid input data for l_max')
                 end if

              else if(words(1)=='LAWDE') then                 ! Law for rho
                 if(words(2)=='CONST') then
                    lawde_tur=0
                 else if(words(2)=='DENSI') then
                    lawde_tur=1
                 end if

              else if(words(1)=='LAWVI') then                 ! Law for mu
                 if(words(2)=='CONST') then
                    lawvi_tur=0
                 else if(words(2)=='VISCO') then
                    lawvi_tur=1
                 end if

              else if(words(1)=='TURBU') then                 ! turbulent Prandtl number
                 prtur_tur = getrea('TURBU',0.0_rp,'#Turbulent Prandtl number')

              else if(words(1)=='REALP') then                 ! Model real parameters   
                 if(words(2)/='AUTOM') then
                    param_tur(1:11)=param(1:11)
                    if (kfl_kemod==1) then !RNG model
                       param_tur(1)= 0.7179_rp    !sigma_k
                       param_tur(2)= 0.7179_rp    !sigma_e
                       param_tur(3)= 1.42_rp      !c1
                       param_tur(4)= 1.68_rp      !c2
                       param_tur(5)= 0.0_rp       !c3
                       param_tur(6)= 0.085_rp     !cmu
                       cmu_st      = param_tur(6) !copy only in master
                    else if(kfl_kemod==2) then !Realizable 
                       param_tur(1)= 1.0_rp       !sigma_k
                       param_tur(2)= 1.2_rp       !sigma_e
                       param_tur(3)= 0.0_rp       !c1
                       param_tur(4)= 1.90_rp      !c2
                       param_tur(5)= 0.0_rp       !c3
                       param_tur(6)= 0.09_rp      !cmu
                       cmu_st      = param_tur(6) !copy only in master
                       if (rwind) then ! realizable for wind
                          param_tur(6)= 0.0333_rp      !cmu
                          cmu_st      = param_tur(6)   !copy (only in master, later in all subdomains (tur_sendat))
                          ! sigma_e to satisfy wall equilibrium
                          param_tur(2)= 0.41_rp*0.41_rp/( param_tur(4)*sqrt(cmu_st)-0.43_rp &
                               *sqrt(cmu_st/0.09_rp))    
                       end if
                       
                    end if
                 end if

              else if(words(1)=='INTEG') then                 ! Model integer parameters 
  
                 if(words(2)/='AUTOM') then
                    ipara_tur(1:2)=int(param(1:2))
                 end if

              else if(words(1)=='AIRDE') then
                 densa_tur=getrea('AIRDE',1.0_rp,'#Air density')

              else if(words(1)=='AIRVI') then
                 visca_tur=getrea('AIRVI',1.0_rp,'#Air viscosity')

              else if(words(1)=='INITI') then
                 !
                 ! Initial Fields
                 !
                 if(words(2)=='CODES') then ! Initial Fields given by codes
                    call runend('TUR_REAPHY: "CODES" OPTION IN INITIAL_CONDITION IS NOT AVAILABLE')
                 else ! Initial Fields given on nodes
                    call ecoute('tur_reaphy')
                    do while(words(1)/='ENDIN')
                       if (words(1) == 'KEY  ') then
                          kfl_inifi_tur(1) = 1
                          nfiel_tur(1) = -getint('FIELD',1_ip,'#Field Number for kinetic energy')
                       else if (words(1) == 'NUFOR') then
                          kfl_inifi_tur(1) = 1
                          kfl_inifi_tur(2) = 1
                          nfiel_tur(1) = -getint('FIELD',1_ip,'#Field Number for kinetic eddy viscosity')
                       else if (words(1) == 'EPSIL') then
                          kfl_inifi_tur(2) = 1
                          nfiel_tur(2) = -getint('FIELD',1_ip,'#Field Number for second variable (epsilon)')
                       else if (words(1) == 'OMEGA') then
                          kfl_inifi_tur(2) = 2
                          nfiel_tur(2) = -getint('FIELD',1_ip,'#Field Number for second variable (omega)')
                       else if (words(1) == 'THIRD') then ! in case of third variable, not programmed yet
                          kfl_inifi_tur(3) = 1
                          nfiel_tur(3) = -getint('FIELD',1_ip,'#Field Number for third variable')
                       end if
                       call ecoute('tur_reaphy')
                    end do
                    if (kfl_inifi_tur(1) > 0 .and. kfl_inifi_tur(2) == 0) then
                       ! hay que ver aun...
                       call runend('TUR_REAPHY: KEY AND SECOND VARIABLE (EPSILON OR OMEGA) MUST BE GIVEN')
                    else if ((kfl_inifi_tur(1) + kfl_inifi_tur(2) + kfl_inifi_tur(3)) == 0) then
                       call runend('TUR_REAPHY: NO FIELDS WERE GIVEN')
                    end if
                 end if

              end if
           end do
        end if
     end do
     !
     ! Define number of turbulence variables
     !
     if(TUR_ONE_EQUATION_MODEL) then
        nturb_tur=1                                           ! One-equation models
     else if(TUR_K_EPS_V2_F) then
        nturb_tur=4                                           ! k-eps-v2-f
     else if(TUR_K_EPS_PHI_F) then
        nturb_tur=4                                           ! k-eps-phi-f
     else
        nturb_tur=2                                           ! Two-equation models
     end if
     !
     ! ABL model
     !
     if( kfl_atmbl /= 0 ) then
        ipara_tur(1) = kfl_atmbl
     end if

  end if

end subroutine tur_reaphy
