  !-----------------------------------------------------------------------
  !> @addtogroup ExmediInput
  !> @ingroup    Exmedi
  !> @{ 
  !> @file     exm_reaphy
  !> @author   jaz,mv
  !> @brief    This routine reads the physical problem definition.
  !> @details  And sets up the models to use 
  !> @}
  !-----------------------------------------------------------------------
subroutine exm_reaphy

  use      def_parame
  use      def_inpout
  use      def_master
  use      def_domain
  use      mod_memchk
  use      def_exmedi
  use mod_ecoute, only :  ecoute
  use mod_messages, only : messages_live, livinf  
  use mod_opfpos, only: postpr_intto8

  implicit none
  integer(ip) :: imate,iipar,iauxi,nauxi
  integer(ip) :: igrou,ivalu,istim,iall_materials,icoor
  character(300)           :: messa

  thiso_exm = huge(1.0_rp)
  nstrb_exm=0_ip
  kfl_nodif_exm = 0_ip ! compute diffusion terms 
  nvint_exm = 2_ip*mitim


  if( INOTSLAVE ) then

     !
     ! Reach the section
     !
     rewind(lisda)
     call ecoute('exm_reaphy')
     do while(words(1)/='PHYSI')
        call ecoute('exm_reaphy')
     end do

     !
     ! Initializations (defaults)
     !
     kfl_stead_exm = 0_ip                    ! No steady-state
     kfl_timei_exm = 1_ip                     ! No steady-state
     kfl_appli_exm = 0_ip                     ! Applied currents and starting pots. are off
     kfl_gcoup_exm = 0_ip                     ! No geometric coupling with SOLIDZ
     kfl_genal_exm = 1_ip                     ! General algorithm type (EXPLICIT=1, decoupled 
     !    implicit, monolithic implicit, ...)

!!!!! LEAVE THIS FOR RETROCOMPATIBILITY, TO BE ELIMINATED
     kfl_gemod_exm = 0_ip                     ! Default: FHN
     !    (TT,LR,BR, ...) or 0=no-subcell or approximate models (FHN, FENTON...)
!!!!!

     kfl_cemod_exm = 1_ip                     ! Cell model (propagation model): monodomain 1 o bidomain 2
     kfl_ptrig_exm = 0_ip                     ! Stimuli start by time and not by pressure trigger     

     kfl_stree_exm = 0_ip                     ! No streeter fiber field created
     kfl_voini_exm = 0_ip                     ! No initial potential imposed

     kfl_psecg_exm = 0_ip 
     kcopeecg_exm = -1_ip 
     kfl_paced_exm = 1_ip 

     kfl_fract_diffusion_exm = -1_ip ! Fractional diffusion is deactivated by default
     fract_diff_coef_exm = 1.0_rp        ! Default value of Fractional diffusion is 1.
     fract_diff_nintp_exm = 0_ip       ! Default value for integration points 0
     
     modfi_exm     = 0_ip                     ! Field: Fiber model
     modst_exm     = 0_ip                     ! Field: Starting stimuli
     nstim_exm     = 1_ip                     ! One starting stimulus
     nconc_exm = 8_ip                         ! TT MODEL number of ionic concentration
     nauxi_exm = 12_ip                        ! TT MODEL 
     nicel_exm = 17_ip                        ! TT MODEL 
     apval_exm     = 0.0_rp                               ! Applied current intensity
     aplap_exm     = 0.0_rp                               ! Applied current time lapse
     apcen_exm     = 0.0_rp                               ! Applied current center
     aprea_exm     = 0.0_rp                               ! Applied current reach
     aploo_exm     = 0.0_rp                               ! Loops

     ngrou_exm     = 0
     nstis_exm     = 0

     xthri_exm     = 1.0e20
     xthrs_exm     = 1.0e20
     xmccm_exm     = 1.0_rp

     kfl_hetermate_exm = 0_ip                      ! homogeneous cell model
     xmccmmate_exm     = 1.0_rp
     kfl_atbhe_exm = 0_ip 
     iall_materials= 0                                    ! Each material reads its own set of properties

     fiaxe_exm = 1.0_rp
     strbo_exm(1)=1_ip
     strbo_exm(2)=2_ip
     strbo_exm(3)=3_ip
     stran_endo_exm=60.0_rp
     stran_epi_exm=200.0_rp ! This vale is to check if there's an input

     kfl_hfmod_exm = 0_ip 
     kfl_hfmodmate_exm = 0_ip 
     kfl_inaga_exm = 0_ip  
     gdiff_exm = 0.0_rp

     ! FOR SUBCELLULAR IONIC CURRENTS MODELS (TT, LR, BR, ...)
     nauxi_exm = 1    ! Default value of number of variables used for activation/inactivation 
     !modab_exm = 1
     !
     ! Begin to read data
     !
     messa = &
          '        READING PHYSICS...'
     call livinf(0_ip,messa,one)
     messa = &
          '        MODULE EXMEDI USES THESE UNITS:'
     call livinf(0_ip,messa,one)
     messa = &
          '          [D] = cm2 / msec  ;   [C_m] = muF / cm2'
     call livinf(0_ip,messa,one)
     messa = &
          '          [I] = muA         ;   [phi] = mV'
     call livinf(0_ip,messa,one)
     messa = &
          '          [x] = cm          ;   [t]   = s'
     call livinf(0_ip,messa,one)
     messa = &
          '          Stimuli: [Current density] =  muA / mm3      '
     call livinf(0_ip,messa,one)
     messa = &
          '        WARNING: TO OBTAIN t IN SECONDS D IS INTERNALLY MULTIPLIED BY 1000.'
     call livinf(0_ip,messa,one)

     !--><group>
     !-->    <groupName>PHYSICAL_PROBLEM</groupName>
     !
     ! ADOC[0]> $-----------------------------------------------------------------------
     ! ADOC[0]> $ Physical properties definition
     ! ADOC[0]> $-----------------------------------------------------------------------
     ! ADOC[0]> PHYSICAL_PROBLEM
     !
     do while(words(1)/='ENDPH')
        call ecoute('exm_reaphy')
        if(words(1)=='PROBL') then
           !
           ! Problem definition data
           !
           call ecoute('exm_reaphy')
           !
           ! ADOC[1]> PROBLEM_DEFINITION
           !
           do while(words(1)/='ENDPR')
              !
              ! ADOC[2]> ALL_MATERIALS: ON | OFF
              !
              ! ADOC[d]> ALL_MATERIALS: Properties are applied to all the materials in the problem. If the flag is not present, everty 
              if(words(1)=='ALLMA') then   ! When many materials are present, all of them share properties
                 iall_materials = 1
                 if(words(2)=='OFF  ') iall_materials = 0
              else if(words(1)=='APPLI') then               ! Applied currents           

                 call runend("EXM_COMAPP: ALWAYS STARTING POTENTIAL OPTION")

                 call ecoute('exm_reaphy')

                 !
                 ! ADOC[2]> PSEUDOECG: ON
                 !
                 ! ADOC[d]> Compute the pseudo-ecg
                 !
              else if(words(1)=='PSEUD') then
                 if (words(2)=='ON') then
                    nvint_exm = 1
                    kfl_psecg_exm = 1_ip 
                    kcopeecg_exm = kcopeecg_exm + 1
                    messa = &
                         '        CALCULATING PSEUDO-ECG'
                    call livinf(0_ip,messa,one)
                    do while (words(1)/='ENDPS')
                       call ecoute('exm_reaphy')   
                       if (words(1) == 'NUMBE') then
                          nrootecg_exm = int(param(1))     
                          if ( nrootecg_exm > size(coordecg_exm,2,KIND=ip) ) call runend('EXM_REAPHY: TOO MANY ECG ROOTS, MAXVAL IS 256!!')
                       elseif (words(1) == 'FREQU') then
!!!                          frequecg_exm = param(1)
                          nvint_exm = int(param(1))
                       elseif (words(1) == 'COORD') then
                          do icoor = 1,nrootecg_exm
                             call ecoute('exm_reaphy')   
                             coordecg_exm(1:ndime,icoor) = param(1:ndime)
                          end do
                       end if
                    end do
                 end if

                 !
                 ! ADOC[2]> GEO_COUPLING: EULERIAN | LAGRANGIAN 
                 !
                 ! ADOC[d]> GEO_COUPLING: Mesh in wich the excitable media problem is solved.
                 ! ADOC[d]> EULERIAN: The problem is solved in the original fixed mesh.
                 ! ADOC[d]> LAGRANGIAN: The problem is solved in the deformed mesh
                 ! Pseudo-ecg

              else if(words(1)=='GEOCO') then               ! Geometric coupling: 
                 if (words(2)=='LAGRA') then
                    kfl_gcoup_exm = 1_ip    !    is SOLIDZ moving my mesh?
                    messa = &
                         '        RUNNING EXMEDI ON A DEFORMED MESH'
                    call livinf(0_ip,messa,one)
                 else if (words(2)=='NEWLA') then
                    kfl_gcoup_exm = 2_ip    !    is SOLIDZ moving my mesh?
                    messa = &
                         '        RUNNING EXMEDI ON A DEFORMED MESH, WITH A NEW-LAGRANGIAN FORMULATION'
                    call livinf(0_ip,messa,one)
                    call runend('EXM_REAPHY: NEWLAGRANGIAN DEPRECATED. USE LAGRANGIAN.')
                 else if (words(2)=='EULERIAN') then
                    kfl_gcoup_exm = 0_ip 
                    messa = &
                         '        RUNNING EXMEDI ON A NON-DEFORMED MESH'
                    call livinf(0_ip,messa,one)
                 end if
                 !
                 ! ADOC[2]> STARTING POTENTIAL
                 !
              else if(words(1)=='START') then               ! Starting potential
                 kfl_appli_exm = 200_ip 
                 kfl_appty_exm = 1_ip                           ! Decay time considered, good for TT-like models
                 messa = &
                      '        STARTING IMPULSES...'
                 call livinf(0_ip,messa,one)                 
                 !
                 ! ADOC[2]> FLASH: ON | OFF
                 !
                 ! ADOC[d]> FLASH: When ON, no temporal variation of the stimuli is considered. the option LAPSE has no effect.
                 ! ADOC[d]> This option must only be used with Fitzhugh nagumo celular model
                 if(exists('FLASH')) then
                    kfl_appty_exm = 2_ip       ! No decay considered, just a flash, good for FHN
                    if(words(2)=='OFF  ') kfl_appty_exm = 1_ip 
                    messa = &
                         '          TYPE OF INITIAL IMPULSE: FLASH'
                    call livinf(0_ip,messa,one)
                 end if
                 !
                 ! ADOC[2]> RESET: 
                 !
                 ! ADOC[d]> RESET: Vaya usted a saber que significa "reset time for loops"
                 aploo_exm(1) = getrea('RESET',0.0_rp,'#Reset time for loops')  
                 if(exists('RESET')) then
                    messa = &
                         '          RESET AT: '//trim(retost(aploo_exm(1)))//' secs'
                    call livinf(0_ip,messa,one)
                 end if
                 !
                 ! ADOC[2]> LOOPT
                 !
                 ! ADOC[d]> RESET: Vaya usted a saber que significa "loop time"
                 aploo_exm(2) = getrea('LOOPT',0.0_rp,'#Loop time')  
                 if(exists('LOOPT')) then
                    messa = &
                         '          LOOP EVERY: '//trim(retost(aploo_exm(2)))//' secs'
                    call livinf(0_ip,messa,one)
                 end if

                 if(exists('TABLE')) then
                    kfl_appva_exm = 1_ip 
                    messa = &
                         '          READING STARTING POTENTIAL FROM TABLE...'
                    call livinf(0_ip,messa,one)
                    call ecoute('exm_reaphy')
                    if (words(1) == 'NSTIM') then       ! Reach from the center
                       nstim_exm = int(param(1))
                       if (words(2) == 'BOUND') nstim_exm = -1_ip 
                       if (words(3) == 'CURDE') kfl_appva_exm = 2_ip 
                       if (words(3) == 'VOLTA') kfl_appva_exm = 1_ip 
                    else 
                       call runend('EXM_REAPHY: GIVE FIRST THE TOTAL NUMBER OF STARTING STIMULI')
                    end if

                    if(size(apval_exm, KIND=ip) .lt. nstim_exm) then
                       call runend('EXM_REAPHY: apval_exm: MAXIMUM NUMBER OF STARTING STIMULI EXCEEDED ')
                    end if
                    if(size(apval_exm, KIND=ip) .lt. nstim_exm/ndime) then
                       call runend('EXM_REAPHY: apval_exm: MAXIMUM NUMBER OF STARTING STIMULI EXCEEDED ')
                    end if
                    if(size(aprea_exm, KIND=ip) .lt. nstim_exm) then
                       call runend('EXM_REAPHY: aprea_exm: MAXIMUM NUMBER OF STARTING STIMULI EXCEEDED ')
                    end if
                    if(size(aptim, KIND=ip) .lt. nstim_exm) then
                       call runend('EXM_REAPHY: aptim: MAXIMUM NUMBER OF STARTING STIMULI EXCEEDED ')
                    end if

                    nauxi= nstim_exm
                    if (nstim_exm < 0) nauxi= 1                    
                    do istim=1,nauxi
                       ! SIEMPRE DARLAS EN ESTE ORDEN: VALUE, CENTER (DE TAMANYO NDIME), REACH, TSTAR, LAPSE
                       call ecoute('exm_reaphy')

                       !Check if there enough numeric columns in the file
                       if (nnpar < ndime+4) then
                          call runend("EXMEDI: Starting potential table row "//trim(intost(istim))//" has less than "//trim(intost(ndime+4))//" numeric columns")
                       end if

                       apval_exm(istim) = param(1)
                       apcen_exm(1:ndime,istim) = param(2:ndime+1)
                       aprea_exm(istim) = param(ndime+2)
                       aptim(istim) = param(ndime+3)
                       aplap_exm(istim) = param(ndime+4)
                    end do

                 else  if(exists('FIELD')) then
                    call ecoute('exm_reaphy')                       
                    if (exists('CURDE')) then 
                       kfl_appva_exm = 2_ip 
                       modst_exm = -getint('FIELD',1_ip,'#Field Number for the current density starting stimuli') 
                       messa = &
                            '          CURRENT DENSITY INITIAL STIMULI STORED IN FIELD: '//trim(intost(-modst_exm))
                       call livinf(0_ip,messa,one)
                    else if (exists('VOLTA')) then 
                       modst_exm = -getint('FIELD',1_ip,'#Field Number for the voltage starting stimuli') 
                       kfl_appva_exm = 1_ip 
                       messa = &
                            '          VOLTAGE INITIAL STIMULI STORED IN FIELD: '//trim(intost(-modst_exm))
                       call livinf(0_ip,messa,one)
                    end if
                 else       
                    call ecoute('exm_reaphy')
                    messa = &
                         '          READING STARTING POTENTIAL FROM PARAMETERS...'
                    call livinf(0_ip,messa,one)
                    if (words(1) == 'NSTIM') then       ! Reach from the center
                       nstim_exm = getint('NSTIM',-1_ip,'#Number of initial stimuli')
                       if (words(2) == 'BOUND') then
                          nstim_exm = -1_ip 
                          nstis_exm = getint('BOUND', 2_ip,'#Set number for the initial stimuli')
                       end if
                    else 
                       call runend('EXM_REAPHY: GIVE FIRST THE TOTAL NUMBER OF STARTING STIMULI')
                    end if
                    nauxi= nstim_exm
                    if (nstim_exm < 0) nauxi = 1_ip 
                    do while(words(1)/='ENDST')
                       if (words(1) == 'VALUE' .or. words(1) == 'VOLTA') then                 ! Value
                          kfl_appva_exm = 1_ip 
                          apval_exm(1:nauxi) = param(1:nauxi)
                       else if (words(1) == 'CURDE') then            ! Density current            
                          kfl_appva_exm = 2_ip 
                          apval_exm(1:nauxi) = param(1:nauxi)
                       else if (words(1) == 'CENTE') then       ! Center of application
                          do iauxi=1,nauxi
                             iipar=(iauxi-1)*ndime + 1
                             apcen_exm(1:ndime,iauxi) = param(iipar:ndime+iipar-1)
                          end do
                       else if (words(1) == 'LAPSE') then       ! Time lapse of application
                          aplap_exm(1:nauxi) = param(1:nauxi)
                       else if (words(1) == 'TSTAR') then       ! Start time
                          if (words(2) == 'PTRIG') then
                             kfl_ptrig_exm = 1_ip 
                             aptim(1:nauxi) = param(2:nauxi+1)
                          else 
                             aptim(1:nauxi) = param(1:nauxi)
                          end if
                       else if (words(1) == 'REACH') then       ! Reach from the center
                          aprea_exm(1:nauxi) = param(1:nauxi)
                       else if (words(1) == 'NSTIM') then       ! Reach from the center
                          if (nauxi==0) nauxi = int(param(1), kind=ip)
                       end if
                       call ecoute('exm_reaphy')

                    end do

                 end if

              end if
              call ecoute('exm_reaphy')
           end do

        else if(words(1)=='PROPE') then
           !conve= 0_ip
           ! Default: Ca+ concentration is not coming from a cell model, but from time spent after a wave 
           ! has passed, like FHN or Fenton.
           !kfl_voini_exm     = 0_ip           

           imate= 1_ip
           messa = &
                '        READING MATERIALS...'
           call livinf(0_ip,messa,one)
           do while(words(1)/='ENDPR')
              if( words(1) == 'NODIF' ) then
                 kfl_nodif_exm = 1_ip                  ! do not compute diffusion terms
                 messa = &
                      '        WARNING: NO DIFFUSION TERMS COMPUTED, NO MATTER THE DIFFUSION VALUE READ!!!! '
                 call livinf(0_ip,messa,one)
              else if( words(1) == 'MATER' ) then
                 imate = getint('MATER',1_ip,'#Current material')
                 messa = &
                      '        MATERIAL: '//trim(intost(imate))
                 call livinf(0_ip,messa,one)

                 if( imate > nmate_exm ) then
                    call runend('EXM_REAPHY: THIS MODULE HAS MORE MATERIALS THAN THOSE DEFINED IN DOM.DAT')
                 end if

                 if( imate > kfl_exm_max_nmaterials ) then
                    call runend('EXM_REAPHY: THIS MODULE HAS MORE MATERIALS THAN HARDCODED IN EXMEDI. MODIFY kfl_exm_max_nmaterials AND RECOMPILE ALYA.')
                 end if



              else if (words(1)=='CONTI') then         ! continuum model properties 
                 messa = &
                      '        CONTINUUM MODEL '
                 call livinf(0_ip,messa,one)
                 do while(words(1)/='ENDCO')

                    !Fractional  diffusion
                    ! BCAM collaboration (ncusimano / lgerardo [at] bcamath.org )
                    if(words(1)=='FRACT') then
                       messa = '           USING FRACTIONAL DIFFUSION'  
                       call livinf(0_ip,messa,one)

                       do while(words(1)/='ENDFR')
                            call ecoute('exm_reaphy')
                            if(words(1)=='COEFF')then
                               fract_diff_coef_exm(imate) = getrea('COEFF',1.0_rp,'#Fractional diffusion coefficient')  
                            elseif(words(1)=='ITERA')then
                               kfl_fract_diffusion_exm(imate)=getint('ITERA',1_ip,'#Fractional diffusion iterations')
                            elseif(words(1)=='NUMBE')then
                               fract_diff_nintp_exm(imate)=getint('NUMBE',1_ip,'#Number of integration points')
                            endif
                       end do
                        ! Checking parameters for fractional diffusion
                        if (fract_diff_coef_exm(imate).gt.1.0_rp .or. fract_diff_coef_exm(imate).le.0.0_rp) then
                           call runend('FRACTIONAL COEFFICIENT SHOULD BE 0>=S>=1')
                        endif
                        if(kfl_fract_diffusion_exm(imate).eq.-1_ip)then
                           call runend('FRACTIONAL DIFFUSION ITERATIONS CANT BE ZERO')
                        endif
                        if (fract_diff_nintp_exm(imate).le.0_ip) then
                           call runend('NUMBER OF INTEGRATION POINTS SHOULD BE >=0')
                        endif
                    endif

                    if(words(1)=='INCON' .or. words(1)=='CONDU') then  ! we retain INCON for retrocompatibility
                       call runend('EXM_REAPHY: USE INDIF FOR THE INTRACELLULAR DIFFUSIVITY')
                    else if(words(1)=='INDIF') then  
                       gdiff_exm(1,1,imate) = param(1)
                       gdiff_exm(1,2,imate) = param(2)
                       gdiff_exm(1,3,imate) = param(3)                       
                    else if (exists('INITI')) then
                       kfl_voini_exm(imate) = 1_ip 
                       voini_exm(imate) = getrea('INITI',0.0_rp,'#Initial voltage of the material')  
                    end if
                    call ecoute('exm_reaphy')
                 end do
                 messa = &
                      '           DIFFU 1: '//trim(retost(gdiff_exm(1,1,imate)))//' cm2/msec'
                 call livinf(0_ip,messa,one)
                 messa = &
                      '           DIFFU 2: '//trim(retost(gdiff_exm(1,2,imate)))//' cm2/msec'
                 call livinf(0_ip,messa,one)
                 messa = &
                      '           DIFFU 3: '//trim(retost(gdiff_exm(1,3,imate)))//' cm2/msec'
                 call livinf(0_ip,messa,one)
                 messa = &
                      '           INITIAL VOLTAGE: '//trim(retost(voini_exm(imate)))//' mV'
                 call livinf(0_ip,messa,one)

                 ! internally convert diffusion units from ms to 
                 gdiff_exm(1,1,imate) = 1000.0_rp * gdiff_exm(1,1,imate) 
                 gdiff_exm(1,2,imate) = 1000.0_rp * gdiff_exm(1,2,imate) 
                 gdiff_exm(1,3,imate) = 1000.0_rp * gdiff_exm(1,3,imate) 



              else if(words(1)=='CELLM') then  ! cell model properties
                 messa = &
                      '        CELL MODEL '
                 call livinf(0_ip,messa,one)

                 if (words(2)=='FITZH') then
                    messa = words(2)
                    kfl_cellmod(imate) = 1_ip                        
                    !default initial values            
                    !    ALPHA= 0.1_rp   $ old alpha
                    !    BETA = 0.01_rp  $ old epsilon
                    !    DELTA= 0.5_rp   $ old gamma_membrane
                    !    C1FHN= 10.0_rp  $ old c_membrane_const
                    !    C2FHN= 0.5_rp   $ no estaba
                    xmopa_exm(1,imate)=0.13_rp  ! original values of paper: a=0.13;
                    xmopa_exm(2,imate)=0.0035_rp ! original values of paper: b=0.013;
                    xmopa_exm(3,imate)=0.26_rp  ! original values of paper: c1=0.26;
                    xmopa_exm(4,imate)=0.1_rp   ! original values of paper: c2=0.1;
                    xmopa_exm(5,imate)=1.0_rp    ! original values of paper: d=1.0;

                    !
                    ! Reference voltages 
                    ! needed and used only by the FHN model:
                    poref_fhn_exm(1)  = -86.2_rp
                    poref_fhn_exm(2)  =  35.0_rp

                    if (kfl_voini_exm(imate) ==1_ip  ) poref_fhn_exm(1)  = voini_exm(imate)

                    do while(words(1)/='ENDCE')

                       if(words(1)=='CMCON') then               ! Model constant Cm
                          xmccmmate_exm(imate)=param(1)                          
                       end if

                       if(words(1)=='ALPHA') then               ! Model constant 
                          xmopa_exm( 1,imate)=param(1)
                       else if(words(1)=='BETA') then               ! Model constant 
                          xmopa_exm( 2,imate)=param(1)
                       else if(words(1)=='C1FHN') then               ! Model constant 
                          xmopa_exm( 3,imate)=param(1)
                       else if(words(1)=='C2FHN') then               ! Model constant 
                          xmopa_exm( 4,imate)=param(1)
                       else if(words(1)=='DELTA') then               ! Model constant 
                          xmopa_exm( 5,imate)=param(1)                             
                       end if
                       call ecoute('exm_reaphy')
                    end do

                 else if (words(2) =='FENTO') then                      
                    messa = words(2)
                    kfl_cellmod(imate) = 2_ip  
                    kfl_fento_exm = 1_ip 

                    call runend('EXM_REAPHYNEW: FENTON-KARMA MODELS TO BE REPROGRAMMED.')

                    if (words(3) == 'BEELE') kfl_fento_exm = 1_ip 
                    if (words(3) == 'MODBE') kfl_fento_exm = 2_ip 
                    if (words(3) == 'LUORU') kfl_fento_exm = 3_ip 
                    if (words(3) == 'GIROU') kfl_fento_exm = 4_ip 
                    !
                    ! Define some model parameters and dimensions for the fenton model
                    !

                    if (kfl_fento_exm == 1_ip ) then

                       xmopa_exm( 3,imate)  = 0.13_rp    ! phi c
                       xmopa_exm( 4,imate)  = 3.33_rp    ! tau v +
                       xmopa_exm( 5,imate)  = 0.04_rp    ! phi v

                       xmopa_exm( 7,imate)  = 1250.0_rp  ! tau v1 -
                       xmopa_exm( 8,imate)  = 19.6_rp    ! tau v2 -
                       xmopa_exm( 9,imate)  = 33.0_rp    ! tauro             
                       xmopa_exm( 10,imate) = 30.0_rp    ! tausi
                       xmopa_exm( 11,imate) = 12.5_rp    ! tauso
                       xmopa_exm( 12,imate) = 0.85_rp    ! phics
                       xmopa_exm( 13,imate) = 4.0_rp     ! gephi
                       xmopa_exm( 16,imate) = 10.0_rp    ! consk
                       xmopa_exm( 17,imate) = 870.0_rp   ! tau w +
                       xmopa_exm( 18,imate) = 41.0_rp    ! tau w -

                    else if (kfl_fento_exm == 2_ip ) then

                       xmopa_exm( 3,imate)  =  0.13_rp    ! phi c
                       xmopa_exm( 4,imate)  =  3.33_rp    ! tau v +
                       xmopa_exm( 5,imate)  =  0.055_rp   ! phi v

                       xmopa_exm( 7,imate)  =  1000.0_rp  ! tau v1 -
                       xmopa_exm( 8,imate)  =  19.2_rp    ! tau v2 -
                       xmopa_exm( 9,imate)  =  50.0_rp    ! tauro          
                       xmopa_exm( 10,imate) =  45.0_rp    ! tausi
                       xmopa_exm( 11,imate) =  8.3_rp     ! tauso
                       xmopa_exm( 12,imate) =  0.85_rp    ! phics
                       xmopa_exm( 13,imate) =  4.0_rp     ! gephi
                       xmopa_exm( 16,imate) =  10.0_rp    ! consk
                       xmopa_exm( 17,imate) =  667.0_rp   ! tau w +
                       xmopa_exm( 18,imate) =  11.0_rp    ! tau w -

                    else if (kfl_fento_exm == 3_ip ) then

                       xmopa_exm( 3,imate)  =  0.13_rp    ! phi c
                       xmopa_exm( 4,imate)  =  10.0_rp    ! tau v +
                       xmopa_exm( 5,imate)  =  0.0_rp     ! phi v 

                       xmopa_exm( 7,imate)  =  18.2_rp    ! tau v1 - 
                       xmopa_exm( 8,imate)  =  18.2_rp    ! tau v2 -
                       xmopa_exm( 9,imate)  =  130.0_rp   ! tauro                
                       xmopa_exm( 10,imate) =  127.0_rp   ! tausi
                       xmopa_exm( 11,imate) =  12.5_rp    ! tauso 
                       xmopa_exm( 12,imate) =  0.85_rp    ! phics
                       xmopa_exm( 13,imate) =  5.8_rp     ! gephi
                       xmopa_exm( 16,imate) =  10.0_rp    ! consk
                       xmopa_exm( 17,imate) =  1020.0_rp  ! tau w +
                       xmopa_exm( 18,imate) =  80.0_rp    ! tau w -


                    else if (kfl_fento_exm == 4_ip ) then


                       xmopa_exm( 3,imate)  =  0.13_rp    ! phi c
                       xmopa_exm( 4,imate)  =  10.0_rp    ! tau v +
                       xmopa_exm( 5,imate)  =  0.025_rp   ! phi v

                       xmopa_exm( 7,imate)  =  333.0_rp   ! tau v1 -
                       xmopa_exm( 8,imate)  =  40.0_rp    ! tau v2 -
                       xmopa_exm( 9,imate)  =  25.0_rp    ! tauro                
                       xmopa_exm( 10,imate) =  22.0_rp    ! tausi
                       xmopa_exm( 11,imate) =  12.5_rp    ! tauso 
                       xmopa_exm( 12,imate) =  0.85_rp    ! phics
                       xmopa_exm( 13,imate) =  8.7_rp     ! gephi
                       xmopa_exm( 16,imate) =  10.0_rp    ! consk
                       xmopa_exm( 17,imate) =  1000.0_rp  ! tau w +
                       xmopa_exm( 18,imate) =  65.0_rp    ! tau w -
                    end if

                 else if (words(2) =='TENTU') then 
                    messa = words(2)
                    kfl_cellmod(imate) = 3_ip  
                    do while(words(1)/='ENDCE')                       

                       if(words(1)=='CMCON') then               ! Model constant Cm
                          xmccmmate_exm(imate)=param(1)                          
                       end if

                       ! if needed, user-defined model parameters must go here
                       call ecoute('exm_reaphy')
                    end do

                 else if (words(2) =='TTHET') then 
                    messa = words(2)
                    kfl_cellmod(imate) = 4_ip                      
                    nconc_exm = 10_ip                            
                    nauxi_exm = 12_ip                           
                    nicel_exm = 18_ip                                               
                    do while(words(1)/='ENDCE')                       
                       if(words(1)=='CMCON') then               ! Model constant Cm
                          xmccmmate_exm(imate)=param(1)                          
                       end if
                       ! if needed, user-defined model parameters must go here
                       call ecoute('exm_reaphy')
                    end do
                 else if (words(2) =='OHARA') then 
                    messa = words(2)
                    kfl_cellmod(imate) = 5_ip                         
                    nconc_exm = 11_ip                            
                    nauxi_exm = 29_ip                           
                    nicel_exm = 26_ip     
                    if (words(3) == 'INAPA') then
                       messa = words(3)
                       kfl_inaga_exm(imate) = 1_ip 
                    end if     
                    do while(words(1)/='ENDCE')                       
                       if(words(1)=='CMCON') then               ! Model constant Cm
                          xmccmmate_exm(imate)=param(1)                          
                       end if
                       ! if needed, user-defined model parameters must go here
                       call ecoute('exm_reaphy')
                    end do
                 else if (words(2) == 'TOROR') then
                    messa = words(2)
                    kfl_cellmod(imate) = CELL_TORORD_EXMEDI
                    nconc_exm = 14_ip
                    nauxi_exm = 31_ip
                    nicel_exm = 21_ip
                    ! if (words(3) == 'BORDE') then
                    !    messa = words(2)//' '//words(3)
                    !    kfl_borde_exm(imate) = 1_ip
                    !    border_gkrsc_exm(imate) = param(3)
                    ! end if
                 else if (words(2) =='NOMOD') then                      ! no ionic current model
                    messa = words(2)
                    kfl_cellmod(imate) = 0_ip 
                    do while(words(1)/='ENDCE')                       
                       if(words(1)=='CMCON') then               ! Model constant Cm
                          xmccmmate_exm(imate)=param(1)                          
                       end if
                       ! if needed, user-defined model parameters must go here
                       call ecoute('exm_reaphy')
                    end do
                 else if (words(2) =='SCATR') then 
                    messa = words(2)
                    kfl_cellmod(imate) = 6_ip                           
                    nconc_exm = 3_ip                            
                    nauxi_exm = 14_ip                           
                    nicel_exm = 15_ip   
                    do while(words(1)/='ENDCE')   
                       if(words(1)=='CMCON') then               ! Model constant Cm
                          xmccmmate_exm(imate)=param(1)                          
                       end if  
                       if(words(1)=='PACED') then               ! The model can be paced or it can beat on it's own
                          kfl_paced_exm = 1_ip                     
                       end if                    
                       ! if needed, user-defined model parameters must go here
                       call ecoute('exm_reaphy')
                    end do
                 else if (words(2) =='SCVEN') then 
                    messa = words(2)
                    kfl_cellmod(imate) = 7_ip                           
                    nconc_exm = 3_ip                            
                    nauxi_exm = 14_ip                           
                    nicel_exm = 15_ip   
                    do while(words(1)/='ENDCE')    
                       if(words(1)=='CMCON') then               ! Model constant Cm
                          xmccmmate_exm(imate)=param(1)                          
                       end if
                       if(words(1)=='PACED') then               ! The model can be paced or it can beat on it's own
                          kfl_paced_exm = 1_ip   
                          else 
                          kfl_paced_exm = 0_ip                         
                       end if                     
                       ! if needed, user-defined model parameters must go here
                       call ecoute('exm_reaphy')
                    end do
                 end if

                 messa = &
                      '           TYPE: '//adjustl(trim(messa))
                 call livinf(0_ip,messa,one)

              else if(words(1)=='INICE') then                       ! CURRENTS PARAMETERS AND CONSTANTS
                 if (exists('HETER')) then
                    kfl_hetermate_exm(imate) = 1_ip 
                    kfl_heter_exm = 1_ip 
                 end if
                 igrou = 0
                 do while(words(1)/='ENDIN')
                    call ecoute('exm_reaphy')                   
                    if(words(1)=='BEATS') then
                       moneclmate_exm( 1,imate)=int(param(1))             ! Number of Beats to simulate until steady state     
                    else if(words(1)=='CYCLE') then
                       moneclmate_exm( 2,imate)=int(param(1))              ! Heart beat cycle length

                    else if(words(1)=='CELLT') then                        !CELLTYPES = ENDO MID EPI, which celltypes will be used. By default all. This will change for which celltypes the ODE will be solved in exm_oneohr.f90
                       kfl_user_specified_celltypes_exm(imate, 1) = 0_ip !reset
                       kfl_user_specified_celltypes_exm(imate, 2) = 0_ip
                       kfl_user_specified_celltypes_exm(imate, 3) = 0_ip
                       if ( exists('ENDO ') ) then
                           kfl_user_specified_celltypes_exm(imate, 1) = 1_ip          
                       end if
                       if ( exists('MID  ') ) then
                           kfl_user_specified_celltypes_exm(imate, 2) = 1_ip          
                       end if
                       if ( exists('EPI  ') ) then
                           kfl_user_specified_celltypes_exm(imate, 3) = 1_ip          
                       end if

                    !IGNORE_STEADYSTATE = ENDO MID EPI. Write the names of the celltypes (ENDO MID EPI). For the specified celltypes Alya will NOT terminate if steady state  is not reached. Only warning will be displayed
                    else if(words(1)=='IGNOR') then                        
                       if ( exists('ENDO ') ) then
                           kfl_ignore_steadystate_celltypes_exm(imate, 1) = 1_ip          
                       end if
                       if ( exists('MID  ') ) then
                           kfl_ignore_steadystate_celltypes_exm(imate, 2) = 1_ip          
                       end if
                       if ( exists('EPI  ') ) then
                           kfl_ignore_steadystate_celltypes_exm(imate, 3) = 1_ip          
                       end if

                    !STEADY_STATE_VARIABLE = (CALCIUM|VOLTAGE) [TOLERANCE = real]
                    else if(words(1)=='STEAD') then  
                       if (exists('CALCI')) then
                           kfl_steadystate_variable(imate) = EXM_CELL_STEADY_CALCIUM
                       end if
                       
                       if (exists('VOLTA')) then
                           kfl_steadystate_variable(imate) = EXM_CELL_STEADY_VOLTAGE
                       end if

                       if ( exists('TOLER') ) then
                           kfl_steadystate_tolerance(imate) = getrea('TOLER',-1.0_rp,'#TOLERANCE FOR STEADY STATE')
                       end if

                       !call runend("EXM_REAPHY: KEYWORD STEAD IN CELL DEFINITION FOR MATERIAL "//trim(intost(imate))//" IS NOT FOLLOWED BY ANYTHING RECOGNISEABLE.")                           
                                          

   

                    !kfl_hfmodmate_exm 0 or 1 -- flag if hardcoded initial conditions should be used: 0 - uses hardcoded initial conditions for cell, 1 - modified cell type, meaning initial conditions need to be calculated by alya
                    !kfl_hfmod_exm - type of cell: male/female/pig. 
                    else if(words(1)=='MYOCY') then               ! Normal Cell or Heart Failure simulation
                       if(words(2)=='NORMA') then
                          kfl_hfmodmate_exm(imate) = 0_ip
                          ttparmate_exm(:,:,imate) = 1.0_rp
                          !call ecoute('exm_reaphy')
                          if(words(3)=='MALE') then
                             kfl_hfmod_exm(imate) = 2_ip
                             !kfl_hfmodmate_exm(imate) = 1_ip
                          else if(words(3)=='FEMAL') then
                             kfl_hfmod_exm(imate) = 3_ip 
                             !kfl_hfmodmate_exm(imate) = 1_ip 
                          else if(words(3)=='PIG') then
                             kfl_hfmod_exm(imate) = 1_ip                                             
                          end if  
                       else if(words(2)=='MODIF') then
                          kfl_hfmodmate_exm(imate) = 1_ip   
                          kfl_hfmod_exm(imate) = 10_ip  !So that the parameter table is read                     
                          if(words(3)=='MALE') then
                             kfl_hfmod_exm(imate) = 2_ip
                          else if(words(3)=='FEMAL') then
                             kfl_hfmod_exm(imate) = 3_ip 
                          else if(words(3)=='PIG') then
                             kfl_hfmod_exm(imate) = 4_ip                           
                          end if                             
                       end if
                    !else if(words(3)=='MALE') then
                    !         kfl_hfmod_exm(imate) = 2_ip
                    !else if(words(3)=='FEMALE') then
                    !         kfl_hfmod_exm(imate) = 3_ip                          
                        !end if
                       !end if

                       !call ecoute('exm_reaphy')
        
                    !The following case should execute only if NORMAL is specified in celltype, but not PIG,MALE or FEMALE
                    else if(words(1)=='CONDU' .and. kfl_hfmod_exm(imate) >= 1_ip) then
                       call ecoute('exm_reaphy')
                       igrou = igrou+1
                       ivalu = 0
                       if(words(1)=='INAME') then  
                          call ecoute('exm_reaphy')
                          do while(words(1)/='ENDCO') 
                             ivalu=ivalu+1
                             ttparmate_exm(1,ivalu,imate) = param(1)
                             ttparmate_exm(2,ivalu,imate) = param(2)
                             ttparmate_exm(3,ivalu,imate) = param(3)

                             if( (param(2) .eq. 0.0_rp) .or. (param(3) .eq. 0.0_rp)) then
                                !PLEASE, we need to allow blocking some of the currents as not all the experimental models have all the currents
                                call livinf(-17_ip, 'exm_reaphy: SOME CURRENTS ARE BLOCKED IN THE HETEROGENEOUS CELL CONFIGURATION',0_ip)
                             endif
                             call ecoute('exm_reaphy')
                          end do
                       end if


                       if ( ivalu<size(ttparmate_exm,2,KIND=ip) ) then
                           call messages_live("EXM_REAPHY: Read only "//postpr_intto8(ivalu)//" triplets of conductances out of "//postpr_intto8(size(ttparmate_exm,2,KIND=ip))//" from the CONDUCTANCES section for material "//postpr_intto8(imate)//". The unread values will be reset to 1.","WARNING")
                       end if

                       ngrou_exm = igrou
                       !end if
                    else if(words(1)=='ONDRU') then               ! Simulation with or without a drug
                       if(exists('YES  ') .or. exists('ON   ')) then
                          kfl_drugsmate_exm(imate) = 1_ip 
                       else
                          kfl_drugsmate_exm(imate) = 0_ip 
                       end if
                    else if(words(1)=='DOSIS') then
                       call ecoute('exm_reaphy')
                       igrou = igrou+1
                       ivalu = 0
                       if(words(1)=='INAME') then  
                          call ecoute('exm_reaphy')
                          do while(words(1)/='ENDDO') 
                             ivalu=ivalu+1
                             drugdmate_exm(ivalu,imate) = param(1)
                             call ecoute('exm_reaphy')
                          end do
                          !if((ivalu.ne.8_ip) .and. (kfl_drugsmate_exm(imate)==1_ip)) call runend('EXM_REAPHY: DRUGS ARE ON BUT NOT ALL DRUGS ARE DEFINED!')
                       end if
                       ngrou_exm = igrou                            
                    end if
                 end do
              end if

              !call ecoute('exm_reaphy')               				 
              if(words(1)=='ISOCH') then
                 if(words(2)=='HIGHE') then
                    thiso_exm(1)=  1.0_rp         ! isochrones are taken when value becomes higher than threshold
                    thiso_exm(2)= param(2)        ! isochrones threshold
                 else if(words(2)=='LOWER') then
                    thiso_exm(1)= -1.0_rp         ! isochrones are taken when value becomes lower than threshold
                    thiso_exm(2)= param(2)        ! isochrones threshold
                 end if

                 thiso_exm(3)= -1.0_rp

                 if (exists('AUTOS')) then  !automatically determine when to save isochrones
                     thiso_exm(3)= 1.0_rp
                 end if

              else if(words(1)=='IONIZ') then               ! Ionization current function (FHN)
                 continue
              end if

              call ecoute('exm_reaphy')

           end do

        else if (words(1)=='FIBER') then              
           if(words(2)=='XALIG') then
              modfi_exm = 1
           else if(words(2)=='YALIG') then
              modfi_exm = 2
           else if(words(2)=='ZALIG') then
              modfi_exm = 3
           else if(words(2)=='SPHER') then
              modfi_exm = 4
           else
              modfi_exm = -getint('FIELD',1_ip,'#Field Number for fibers')
           end if

          else if (words(1)=='ORTHO') then              
           call ecoute('exm_reaphy')
           do while(words(1).ne.'ENDOR')
              if (words(1)=='SHEET') then             
                 if(words(2)=='XALIG') then
                  modor_exm(1) = 1_ip
                 else if(words(2)=='YALIG') then
                  modor_exm(1) = 2_ip
                 else if(words(2)=='ZALIG') then
                  modor_exm(1) = 3_ip
                 else
                  modor_exm(1) = -getint('FIELD',1_ip,'#Field Number for fibers')
                 end if
              elseif (words(1)=='NORMA') then             
                 if(words(2)=='XALIG') then
                  modor_exm(2) = 1_ip
                 else if(words(2)=='YALIG') then
                  modor_exm(2) = 2_ip
                 else if(words(2)=='ZALIG') then
                  modor_exm(2) = 3_ip
                 else
                  modor_exm(2) = -getint('FIELD',1_ip,'#Field Number for fibers')
                 end if
              else
                call runend('EXM_REAPHY: ORTHO_MODEL FIBER FIELD NOT IDENTIFIED')
              endif
              call ecoute('exm_reaphy')
           enddo

           if((modor_exm(1).eq.0_ip) .or. (modor_exm(2).eq.0_ip) )then
             call runend('EXM_ADARR: ORTHO MODEL DEFINED BUT ONE FIBER FIELD IS MISSING')
           endif


           
        else if (words(1)=='CREAT') then
           call ecoute('exm_reaphy')
           do while(words(1).ne.'ENDCR')
              if(words(1)=='FUNCT') then
                 if (words(2) == 'CUBIC') then
                    kfl_stree_exm=3
                 elseif (words(2) == 'LINEA') then
                    kfl_stree_exm=1
                 else
                    kfl_stree_exm=1 ! Default option is linear
                 endif
              elseif (words(1) == 'VAXIS') then
                 fiaxe_exm(1:3)= param(1:3)

              elseif (words(1) == 'EPICA') then
                 strbo_exm(1) = int(param(1), kind=ip)
                 nstrb_exm=nstrb_exm+1
              elseif(words(1) == 'LEFTE') then
                 strbo_exm(2) = int(param(1), kind=ip)
                 nstrb_exm=nstrb_exm+1
              elseif(words(1) == 'RIGHT') then
                 strbo_exm(3) = int(param(1), kind=ip)
                 nstrb_exm=nstrb_exm+1
              elseif(words(1) == 'ANGLE') then
                 stran_endo_exm = param(1)
                 stran_epi_exm = param(2)
                 if(stran_epi_exm .eq. 0.0_rp) stran_epi_exm = - stran_endo_exm
              endif
              call ecoute('exm_reaphy')
           enddo
        else if (words(1)=='CELLT') then              
           modce_exm = -getint('FIELD',2_ip,'#heterogeneous cell type')
        else if (words(1)=='APEXT') then   
           kfl_atbhe_exm = 1_ip           
           modab_exm = -getint('FIELD',3_ip,'#apex-to-base heterogeneity')
        end if
     end do

     !
     ! Unknowns element-wise dimensions
     !
     ndofn_exm = 1 
     if (kfl_cemod_exm == 2) ndofn_exm = 2_ip 
     ndof2_exm = ndofn_exm*ndofn_exm

     !  FOR NO CELL IONIC CURRENTS MODELS (FHN, ...)
     nevat_exm = mnode    ! <--- explicit (default value)
     if (kfl_genal_exm == 2) nevat_exm = mnode                   ! <--- decoupled implicit
     if (kfl_genal_exm == 3) nevat_exm = ndofn_exm*mnode         ! <--- monolithic implicit

     !
     ! Initialization
     !
     dtinv_exm=0.0_rp

     ! WARNING WARNING: Esta opcion ahora solo tiene sentido para gdiff.
     ! xmopa son variables de FHN, no tiene sentido asignarselas a todos los materiales.
     if (iall_materials == 1) then

        do imate= 1,nmate
           gdiff_exm(1,1,imate) = gdiff_exm(1,1,1) 
           gdiff_exm(1,2,imate) = gdiff_exm(1,2,1) 
           gdiff_exm(1,3,imate) = gdiff_exm(1,3,1) 
           xmopa_exm( 3,imate)  = xmopa_exm(  3,1)  
           xmopa_exm( 4,imate)  = xmopa_exm(  4,1)  
           xmopa_exm( 5,imate)  = xmopa_exm(  5,1)  

           xmopa_exm( 7,imate)  = xmopa_exm(  7,1)  
           xmopa_exm( 8,imate)  = xmopa_exm(  8,1)  
           xmopa_exm( 9,imate)  = xmopa_exm(  9,1)             
           xmopa_exm( 10,imate) = xmopa_exm( 10,1) 
           xmopa_exm( 11,imate) = xmopa_exm( 11,1) 
           xmopa_exm( 12,imate) = xmopa_exm( 12,1) 
           xmopa_exm( 13,imate) = xmopa_exm( 13,1) 
           xmopa_exm( 16,imate) = xmopa_exm( 16,1) 
           xmopa_exm( 17,imate) = xmopa_exm( 17,1) 
           xmopa_exm( 18,imate) = xmopa_exm( 18,1)  
        end do

     end if


     iauxi= -1
     do imate= 1,nmate
        if (kfl_cellmod(imate) .ge. 0_ip ) iauxi= 1
     end do
     if (iauxi == -1) then

        call runend('EXM_REAPHY: FITZHUGH MUST BE EXPLICTILY DEFINED AS CELL=FITZHUGH')

     end if

     !if (kfl_gemod_exm == 1 .or. conve == 1) then
     ! TEN TUSSCHER model  
     ! kfl_cellmod(imate) = 1 --> TenTusscher model 
     !       (Paper 'A model for human ventricular tissue', 2003)
     ! kfl_cellmod(imate) = 2 --> LuoRudyII  model (Paper ..................... , 1994)
     ! kfl_cellmod(imate) = 3 --> BeelerReuter model (Paper ................., 1977)
     !gdiff_exm(1,:,imate) = gdiff_exm(1,:,imate) * 8991997122.56092078050535_rp
     !(1.0_rp / (4.0_rp * 3.1415_rp * 0.00000000000885_rp)) !conversion from SI to CGS


  end if      !! kfl_paral


end subroutine exm_reaphy
