!-----------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_inivar.f90
!> @author  Mariano Vazquez
!> @brief   Initialize data
!> @date    16/11/1966
!> @details Initialize data
!> @} 
!-----------------------------------------------------------------------
subroutine exm_inivar(itask)
  !-----------------------------------------------------------------------
  !****f* Temper/exm_inivar
  ! NAME 
  !    exm_inicar
  ! DESCRIPTION
  !    This routine initializes some variables
  !    ITASK=1 ... When starting the run (from Turnon)
  !    ITASK=2 ... First time step. This is needed as some variables 
  !                are not initialized before
  ! USES
  ! USED BY
  !    exm_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_exmedi 
  use def_solver
  use mod_iofile
  use mod_exm_sld_eccoupling,  only : kfl_exmsld_3Dcou_ecc
  use mod_sld_cou,        only : mod_sld_cou_initvar
  use mod_sld_cou,        only : mod_sld_cou_initexchange
  use mod_arrays, only : arrays_register
  use mod_memory, only : memory_alloca
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: iconc,iauxi,ivari
  integer(ip)             :: imate, icelltype, iionconcentr

  select case(itask)

  case(0_ip)
     !
     ! Register variables
     !
     call arrays_register(  5_ip,(/'TAULO','SCALA','NPOIN','PRIMA'/),taulo         ,ENTITY_POSITION=1_ip)
     call arrays_register(  6_ip,(/'FISOC','SCALA','NPOIN','PRIMA'/),fisoc         ,ENTITY_POSITION=1_ip)
     call arrays_register(  7_ip,(/'KWAVE','SCALA','NPOIN','PRIMA'/),kwave_exm     ,ENTITY_POSITION=1_ip)
     call arrays_register(  8_ip,(/'ISOCC','SCALA','NPOIN','PRIMA'/),isoch_modified,ENTITY_POSITION=1_ip)
     call arrays_register(  9_ip,(/'ELMAG','SCALA','NPOIN','PRIMA'/),elmag         ,ENTITY_POSITION=1_ip,TIME_POSITION=2_ip)
     call arrays_register( 10_ip,(/'VCONC','SCALA','NPOIN','PRIMA'/),vconc         ,ENTITY_POSITION=2_ip,TIME_POSITION=3_ip)
     call arrays_register( 11_ip,(/'QNETT','SCALA','NPOIN','PRIMA'/),qneto_exm     ,ENTITY_POSITION=1_ip)
     call arrays_register( 12_ip,(/'REFHN','SCALA','NPOIN','PRIMA'/),refhn_exm     ,ENTITY_POSITION=1_ip,TIME_POSITION=2_ip)
     call arrays_register( 13_ip,(/'VAUXI','SCALA','NPOIN','PRIMA'/),vauxi_exm     ,ENTITY_POSITION=2_ip,TIME_POSITION=3_ip,COMPONENT_POSITION=1_ip)
     call arrays_register( 14_ip,(/'VICEL','SCALA','NPOIN','PRIMA'/),vicel_exm     ,ENTITY_POSITION=2_ip,TIME_POSITION=3_ip,COMPONENT_POSITION=1_ip)

     !
     ! Postprocess
     !

     ! VARIABLE NAMES:

     postp(1) % wopos ( 1, 1) = 'FIBER' ! Fiber direction
     postp(1) % wopos ( 1, 2) = 'GRAFI'
     postp(1) % wopos ( 1, 4) = 'IOCON'   ! VCONC (ION CONCENTRATIONS), A MULTIDIMENSIONAL NODAL FIELD

     ! UNREGISTERED VARIABLES
     postp(1) % wopos ( 1,15) = 'CALCI'   ! calcium
     postp(1) % wopos ( 1,16) = 'POTAS'   ! calcium

     
     postp(1) % wopos ( 1,20) = 'INTRA'
     postp(1) % wopos ( 1,21) = 'EXTRA'
     postp(1) % wopos ( 1,22) = 'RECOV'
     postp(1) % wopos ( 1,23) = 'SHEET' ! Sheet fiber direction
     postp(1) % wopos ( 1,24) = 'NORMA' ! Normal fiber direction
     postp(1) % wopos ( 1,25) = 'FREE4' ! FREE

     postp(1) % wopos ( 1,26) = 'ISOCH'
     postp(1) % wopos ( 1,27) = 'CECO1'
     postp(1) % wopos ( 1,28) = 'CECO2'
     postp(1) % wopos ( 1,29) = 'CECO3'
     postp(1) % wopos ( 1,30) = 'FREE5'! FREE     
     postp(1) % wopos ( 1,31) = 'CURRE'  ! VICEL (CURRENTS), A MULTIDIMENSIONAL NODAL FIELD
     postp(1) % wopos ( 1,32) = 'QNETO'  ! VICEL (CURRENTS), A MULTIDIMENSIONAL NODAL FIELD

     ! VARIABLE TYPES:

     postp(1) % wopos ( 2, 1) = 'VECTO'
     postp(1) % wopos ( 2, 2) = 'VECTO'
     postp(1) % wopos ( 2, 4) = 'MULTI'

     postp(1) % wopos ( 2,15) = 'SCALA'
     postp(1) % wopos ( 2,16) = 'SCALA'

     postp(1) % wopos ( 2,20) = 'SCALA'
     postp(1) % wopos ( 2,21) = 'SCALA'
     postp(1) % wopos ( 2,22) = 'SCALA'
     postp(1) % wopos ( 2,23) = 'SCALA'
     postp(1) % wopos ( 2,24) = 'SCALA'
     postp(1) % wopos ( 2,25) = 'SCALA'

     postp(1) % wopos ( 2,26) = 'SCALA'
     postp(1) % wopos ( 2,27) = 'VECTO'
     postp(1) % wopos ( 2,28) = 'VECTO'
     postp(1) % wopos ( 2,29) = 'VECTO'
     !+MRV
     postp(1) % wopos ( 2,30) = 'VECTO'     
     !-MRV

     postp(1) % wopos ( 2,31) = 'MULTI'
     postp(1) % wopos ( 2,32) = 'SCALA'


     do iconc= 1,nconc_exm
        ivari= 32 + iconc
        if (iconc.le.9) then
           postp(1) % wopos ( 1,ivari) = 'VCO0'//trim(intost(iconc))
        else
           postp(1) % wopos ( 1,ivari) = 'VCO'//trim(intost(iconc))
        end if
        postp(1) % wopos ( 2,ivari) = 'SCALA'
     end do

     do iauxi= 1,nauxi_exm
        ivari= 32 + nconc_exm + iauxi
        if (iconc.le.9) then
           postp(1) % wopos ( 1,ivari) = 'VAU0'//trim(intost(iconc))
        else
           postp(1) % wopos ( 1,ivari) = 'VAU'//trim(intost(iconc))
        end if
        postp(1) % wopos ( 2,ivari) = 'SCALA'
     end do

     !
     ! Sets variables
     !

     postp(1) % woese (  1) = 'INTRA'
     postp(1) % woese (  2) = 'CALCI'         
     postp(1) % woese (  3) = 'SODIU'         
     postp(1) % woese (  4) = 'POTAS'         

     postp(1) % wonse (  1) = 'INTRA'
     postp(1) % wonse (  2) = 'CALCI'   
     postp(1) % wonse (  3) = 'SODIU'   
     postp(1) % wonse (  4) = 'POTAS'   

     !
     ! Witness variables
     !
     postp(1) % wowit  (1)     = 'INTRA'
     postp(1) % wowit  (2)     = 'COORX'
     postp(1) % wowit  (3)     = 'CAI  '  !VCONC(1:11) ! 2+1
     postp(1) % wowit  (4)     = 'NASS '
     postp(1) % wowit  (5)     = 'KI   '
     postp(1) % wowit  (6)     = 'KISS '
     postp(1) % wowit  (7)     = 'NAI  '
     postp(1) % wowit  (8)     = 'CASS '
     postp(1) % wowit  (9)     = 'CANSR'
     postp(1) % wowit (10)     = 'CAJSR'
     postp(1) % wowit (11)     = 'JRLNP' !Jrelnp
     postp(1) % wowit (12)     = 'JRLP ' !Jrelp
     postp(1) % wowit (13)     = 'CAMKT'
     postp(1) % wowit (14)     = 'INA  '  !VICEL(1:26) ! 13+1
     postp(1) % wowit (15)     = 'INAL '
     postp(1) % wowit (16)     = 'ITO  '
     postp(1) % wowit (17)     = 'ICAL '
     postp(1) % wowit (18)     = 'IKR  '
     postp(1) % wowit (19)     = 'IKS  '
     postp(1) % wowit (20)     = 'IK1  '
     postp(1) % wowit (21)     = 'INCI ' !INaCa_i
     postp(1) % wowit (22)     = 'INCSS' !INaCa_ss
     postp(1) % wowit (23)     = 'INAK '
     postp(1) % wowit (24)     = 'IKB  '
     postp(1) % wowit (25)     = 'INAB '
     postp(1) % wowit (26)     = 'ICAB '
     postp(1) % wowit (27)     = 'IPCA '
     postp(1) % wowit (28)     = 'JDI  ' !Jdiff
     postp(1) % wowit (29)     = 'JDINA' !JdiffNa
     postp(1) % wowit (30)     = 'JDIK '
     postp(1) % wowit (31)     = 'JUP '
     postp(1) % wowit (32)     = 'JLEAK'
     postp(1) % wowit (33)     = 'JTR  '
     postp(1) % wowit (34)     = 'JREL '
     postp(1) % wowit (35)     = 'CAMKA'
     postp(1) % wowit (36)     = 'STIM '
     postp(1) % wowit (37)     = 'ICANA'
     postp(1) % wowit (38)     = 'ICAK '
     postp(1) % wowit (39)     = 'CAMKB'
     

     !
     ! Nullify arrays
     !
     nullify(idima_exm)  
     nullify(kgrfi_exm)        
     nullify(cedif_exm)
     nullify(grafi_exm)
     nullify(fiber_exm)
     nullify(sheet_exm)
     nullify(normal_exm)
     nullify(celty_exm)
     nullify(atbhe_exm)
     nullify(vdiag_exm)     
     nullify(fibe2_exm) 
     nullify(tncod_exm)      
     nullify(tbcod_exm)      
     nullify(amatr_auxi_exm) 
     nullify(appfi_exm)      
     nullify(lmate_exm)    
     nullify(kfl_fixno_exm)
     nullify(kfl_fixbo_exm)  
     nullify(bvess_exm) 
     nullify(bvnat_exm) 
     nullify(vauxi_exm)
     nullify(ticel_exm)    
     nullify(jicel_exm)    
     nullify(vicel_exm)
     nullify(qneto_exm)
     nullify(refhn_exm) 
     nullify(lapno_exm)     
     nullify(kwave_exm)      
     nullify(isoch_modified) 
     !
     ! Pseudo-ecg
     !
     kcopeecg_exm = 0
     !
     ! Solvers
     !     
     call soldef(-1_ip)

     solve(1) % wprob     = 'ACTIVATION_POTENTIAL'
     solve(1) % kfl_solve = 1
     solve(1) % ndofn     = 1

     ! In the case of implicit, these values will be defined in nsa_reanut, when calling reasol
     ! These are default values, corresponding to explicit with global time step
     solve(1) % kfl_algso = SOL_SOLVER_RICHARDSON            
     solve(1) % kfl_preco = SOL_CLOSE_MASS_MATRIX

     ! Time variables
     kfl_timei = 1      ! exmedi is always transient

     !
     ! Materials
     !
     nmate_exm =  nmate
     lmate_exm => lmate
     !lcell_exm => lcell

     epres = 0.0_rp

     sms_conversion_currents_exm = 1000.0_rp

     cpold_exm= 0.0_rp

     do imate=1,size(kfl_user_specified_celltypes_exm, 1, KIND=ip)
        do icelltype=1,size(kfl_user_specified_celltypes_exm, 2, KIND=ip)
           kfl_user_specified_celltypes_exm( imate, icelltype ) = 1_ip   !by deafult the ODE for all celltypes will be executed
           kfl_ignore_steadystate_celltypes_exm( imate, icelltype ) = 0_ip !do not ignore steady state by deafult
        end do
        kfl_steadystate_variable(imate) = EXM_CELL_STEADY_VOLTAGE !use voltage for steady state by default
        kfl_steadystate_tolerance(imate) = -1.0_rp                !let subroutine determine the optimal value
     end do


     !reset ttparmate_exm to 1, afterwards in reaphy it will be filled
     do icelltype=1,size(ttparmate_exm, 1, KIND=ip)
        do iionconcentr=1,size(ttparmate_exm, 2, KIND=ip)
            do imate=1,size(ttparmate_exm, 3, KIND=ip)
                ttparmate_exm(icelltype,iionconcentr,imate) = 1_rp
            end do
        end do
     end do


     !
     ! Flags
     !
     call  mod_sld_cou_initvar()

     kfl_save_convergence_cellmodel = 0_ip
     kfl_save_init_cellmodel = 0_ip

  case ( 1_ip )
     !
     ! More initialisations after reading files
     !
     if(kfl_exmsld_3Dcou_ecc) call mod_sld_cou_initexchange()

     if(kfl_timei_exm==1) ncomp_exm = 2 + kfl_tiacc_exm
     ! 
     !
     ! FOR SUBCELLULAR IONIC CURRENTS MODELS (TT, LR, BR, ...)
     !
     nauxi_exm = 1
     nconc_exm = 1
     nicel_exm = 1
     ngate_exm = 1    ! Default for fhn

     do imate= 1,nmate

        if (kfl_cellmod(imate) == CELL_NOMOD_EXMEDI) then
           nauxi_exm = max(nauxi_exm,1_ip)    ! Default value of number of variables used for activation/inactivation 
           nconc_exm = max(nconc_exm,1_ip)    ! Default value of number of variables used for concentrations
           nicel_exm = max(nicel_exm,1_ip)
        else if (kfl_cellmod(imate) == CELL_FENTON_EXMEDI) then
           call runend("EXM_MEMALL: FENTON KARMA MODELS NOT UPDATED")
           nconc_exm = 2
           ngate_exm = 2
        else if (kfl_cellmod(imate) == CELL_FITZHUGH_EXMEDI) then
           nauxi_exm = max(nauxi_exm,12_ip)  
           nconc_exm = max(nconc_exm, 8_ip)
        else if (kfl_cellmod(imate) == CELL_TT2006_EXMEDI) then
           nauxi_exm = max(nauxi_exm,12_ip)  
           nconc_exm = max(nconc_exm,10_ip)  
           nicel_exm = max(nconc_exm,18_ip) 
        else if (kfl_cellmod(imate) == CELL_OHARA_EXMEDI) then
           nauxi_exm = max(nauxi_exm,29_ip)  
           nconc_exm = max(nconc_exm,11_ip)  
           nicel_exm = max(nconc_exm,26_ip)
        else if (kfl_cellmod(imate) == CELL_TORORD_EXMEDI) then
           nauxi_exm = max(nauxi_exm,31_ip)
           nconc_exm = max(nconc_exm,14_ip)
           nicel_exm = max(nicel_exm,21_ip)
        else if (kfl_cellmod(imate) == CELL_SCATRIA_EXMEDI) then
           nauxi_exm = max(nauxi_exm,14_ip)  
           nconc_exm = max(nconc_exm,3_ip)  
           nicel_exm = max(nconc_exm,15_ip)
        else if (kfl_cellmod(imate) == CELL_SCVENTRI_EXMEDI) then
           nauxi_exm = max(nauxi_exm,14_ip)  
           nconc_exm = max(nconc_exm,3_ip)  
           nicel_exm = max(nconc_exm,15_ip)
           !     else if (kfl_cellmod(imate) == 6) then  !TTINA
           !        nauxi_exm = max(nauxi_exm,14_ip)  
           !        nconc_exm = max(nconc_exm,10_ip)  
           !        nicel_exm = max(nconc_exm,18_ip) 
        end if

     end do
     
  end select

end subroutine exm_inivar
