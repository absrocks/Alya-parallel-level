subroutine tur_inivar(itask)
  !-----------------------------------------------------------------------
  !****f* Temper/tur_inivar
  ! NAME 
  !    tur_inicar
  ! DESCRIPTION
  !    This routine initializes some variables
  !    ITASK=0 ... When reading model (tur_reaphy)
  !    ITASK=1 ... When starting the run (from Turnon)
  !    ITASK=2 ... First time step. This is needed as some variables 
  !                are not initialized before
  ! USES
  ! USED BY
  !    tur_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_solver
  use def_turbul
  use def_domain
  use def_kermod, only : kfl_adj_prob
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: iunkn

  select case(itask)

  case(0_ip)

     !-------------------------------------------------------------------
     ! 
     ! Postprocess Variable names
     !
     !-------------------------------------------------------------------

     postp(1) % wopos (1, 1) = 'KEY  '
     postp(1) % wopos (1, 2) = 'EPSIL'
     postp(1) % wopos (1, 3) = 'OMEGA'
     postp(1) % wopos (1, 4) = 'NUTIL'
     postp(1) % wopos (1, 5) = 'LENGT'
     postp(1) % wopos (1, 6) = 'TURVI'
     postp(1) % wopos (1,10) = ''
     postp(1) % wopos (1,11) = 'YPLUS'
     postp(1) % wopos (1,12) = 'USTAR'
     postp(1) % wopos (1,13) = 'YVSTA'
     postp(1) % wopos (1,14) = 'YSTAR'
     postp(1) % wopos (1,15) = 'RESI1'
     postp(1) % wopos (1,16) = 'RESI2'
     postp(1) % wopos (1,17) = 'RESIT'
     postp(1) % wopos (1,18) = 'PHI  '
     postp(1) % wopos (1,19) = 'F    '
     postp(1) % wopos (1,20) = 'GRSQK'
     postp(1) % wopos (1,21) = 'INTEN'
     postp(1) % wopos (1,22) = 'PROJE'
     postp(1) % wopos (1,23) = 'LIMIT'
     postp(1) % wopos (1,24) = 'RATIO'
     postp(1) % wopos (1,25) = 'UNKN1'
     postp(1) % wopos (1,26) = 'UNKN2'
     postp(1) % wopos (1,27) = 'UNKN3'
     postp(1) % wopos (1,28) = 'UNKN4'
     postp(1) % wopos (1,29) = 'PROJ1'
     postp(1) % wopos (1,30) = 'PROJ2'
     postp(1) % wopos (1,31) = 'PROJ3'
     postp(1) % wopos (1,32) = 'PROJ4'
     postp(1) % wopos (1,33) = 'LINEL'
     postp(1) % wopos (1,34) = 'GROUP'
     postp(1) % wopos (1,35) = 'FDDES'
     postp(1) % wopos (1,36) = 'GDDES'
     postp(1) % wopos (1,37) = 'SSTF1'
     postp(1) % wopos (1,38) = 'SSTF2'
     postp(1) % wopos (1,39) = 'SASSO'
     postp(1) % wopos (1,40) = 'AVKEY'
     postp(1) % wopos (1,41) = 'AVOME'
     postp(1) % wopos (1,42) = 'AVTVI'
     postp(1) % wopos (1,43) = 'VORTI'
     postp(1) % wopos (1,44) = 'CACAS'
     postp(1) % wopos (1,45) = 'ERROR'
     postp(1) % wopos (1,46) = 'FIXT1'    ! I do not call it FIXNO otherwise it will overlap with the one for vel
     postp(1) % wopos (1,47) = 'FIXT2'
     postp(1) % wopos (1,48) = 'SENSM'    ! Mesh sensitivities
     postp(1) % wopos (1,49) = 'MAXLE'    ! Max mixing length     


     postp(1) % wopos (2, 1) = 'SCALA'
     postp(1) % wopos (2, 2) = 'SCALA'
     postp(1) % wopos (2, 3) = 'SCALA'
     postp(1) % wopos (2, 4) = 'SCALA'
     postp(1) % wopos (2, 5) = 'SCALA'
     postp(1) % wopos (2, 6) = 'SCALA'
     postp(1) % wopos (2, 7) = 'SCALA'
     postp(1) % wopos (2, 8) = 'SCALA'
     postp(1) % wopos (2, 9) = 'SCALA'
     postp(1) % wopos (2,10) = ''
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
     postp(1) % wopos (2,41) = 'SCALA'
     postp(1) % wopos (2,42) = 'SCALA'
     postp(1) % wopos (2,43) = 'SCALA'
     postp(1) % wopos (2,44) = 'SCALA'
     postp(1) % wopos (2,45) = 'SCALA'
     postp(1) % wopos (2,46) = 'SCALA'
     postp(1) % wopos (2,47) = 'SCALA'
     postp(1) % wopos (2,48) = 'VECTO'
     postp(1) % wopos (2,49) = 'SCALA'
     !
     ! Witness variables 
     !
     postp(1) % wowit (1)    = 'KEY  ' 
     postp(1) % wowit (2)    = 'EPSIL' 
     !
     ! Nullify pointers
     !
     nullify(kfl_fixno_tur) 
     nullify(bvess_tur)     
     nullify(ustar_tur)        
     nullify(grve2_tur)        
     nullify(grsqk_tur)         
     nullify(dunkn_tur)        
     nullify(unold_tur)       
     nullify(grk12_tur)       
     nullify(grono_tur)       
     nullify(greps_tur)       
     nullify(grphi_tur)      
     nullify(unpro_tur)      
     nullify(unprr_tur)      
     nullify(unpgr_tur)      
     nullify(detur_tur)   
     nullify(vitur_tur)   
     nullify(vorti_tur)   
     nullify(fddes_tur)   
     nullify(gddes_tur)  
     nullify(sstf1_tur)
     nullify(sstf2_tur) 
     nullify(sasso_tur) 
     nullify(avkey_tur)
     nullify(avome_tur)
     nullify(avtvi_tur)
     nullify(produ_tur)   
  

     !-------------------------------------------------------------------
     ! 
     ! Solver
     !
     !-------------------------------------------------------------------

     call soldef(-5_ip)
     solve(1)   % wprob     = '1ST_TURBUL_VARIABLE' ! k
     solve(2)   % wprob     = '2ND_TURBUL_VARIABLE' ! eps, w
     solve(3)   % wprob     = '3RD_TURBUL_VARIABLE' ! v2
     solve(4)   % wprob     = '4TH_TURBUL_VARIABLE' ! f, phi
     solve(5)   % wprob     = 'WALL_DISTANCE'       ! d
     solve(1:5) % kfl_solve = 1

  case(1_ip)
     !
     ! Number of turbulent components 
     !
     if(kfl_timei_tur==1) then
        if(kfl_tisch_tur==1) then
           !
           ! Trapezoidal rule
           !
           ncomp_tur=3
        else if(kfl_tisch_tur==2) then
           !
           ! BDF scheme
           !
           ncomp_tur=2+kfl_tiacc_tur
        end if
     else
        ncomp_tur=2     
     end if
     nprev_tur = min(3_ip,ncomp_tur)  ! Last time step or global iteration
     nbdfp_tur = 2
     
     if(kfl_timei_tur==0) then
        dtinv_tur=1.0_rp              ! Time step
     else
        kfl_timei=1                   ! Time integration
     end if
     ittot_tur = 0                    ! Total iterations
     kfl_tiaor_tur=kfl_tiacc_tur      ! Time accuracy: save original value
     if (kfl_adj_prob == 1) kfl_stead_tur = 0

     dtmax_tur = 0.0_rp ! otherwise full ckeck stops because it is unitialized when doing max in tur_tistep 

  case(2_ip)

     !-------------------------------------------------------------------
     !
     ! Models selections
     !
     !-------------------------------------------------------------------

     TUR_FAMILY_K_EPS         = .false.
     TUR_FAMILY_K_OMEGA       = .false.
     TUR_ONE_EQUATION_MODEL   = .false.  
     TUR_TWO_EQUATION_MODEL   = .false.  

     TUR_SPALART_ALLMARAS     = .false.
     TUR_K_XU_CHIEN           = .false. 
     TUR_K_EPS_STD            = .false. 
     TUR_K_EPS_LAUNDER_SHARMA = .false.
     TUR_K_EPS_CHIEN          = .false.
     TUR_K_EPS_LAM_BREMHORST  = .false.
     TUR_TWO_LAYER_RODI       = .false.
     TUR_TWO_LAYER_XU_CHEN    = .false. 
     TUR_K_EPS_V2_F           = .false. 
     TUR_K_EPS_JAW_HWANG      = .false.
     TUR_K_EPS_NAGANO_TAGAWA  = .false.
     TUR_K_EPS_PHI_F          = .false. 
     TUR_K_OMEGA              = .false.   
     TUR_K_OMEGA_BREDBERG     = .false.  
     TUR_SST_K_OMEGA          = .false. 

     if(kfl_model_tur== 1) TUR_SPALART_ALLMARAS     = .true. 
     if(kfl_model_tur== 2) TUR_K_XU_CHIEN           = .true. 
     if(kfl_model_tur==11) TUR_K_EPS_STD            = .true. 
     if(kfl_model_tur==12) TUR_K_EPS_LAUNDER_SHARMA = .true. 
     if(kfl_model_tur==13) TUR_K_EPS_CHIEN          = .true. 
     if(kfl_model_tur==14) TUR_K_EPS_LAM_BREMHORST  = .true. 
     if(kfl_model_tur==15) TUR_TWO_LAYER_RODI       = .true. 
     if(kfl_model_tur==16) TUR_TWO_LAYER_XU_CHEN    = .true. 
     if(kfl_model_tur==17) TUR_K_EPS_V2_F           = .true. 
     if(kfl_model_tur==18) TUR_K_EPS_JAW_HWANG      = .true. 
     if(kfl_model_tur==19) TUR_K_EPS_NAGANO_TAGAWA  = .true.
     if(kfl_model_tur==20) TUR_K_EPS_PHI_F          = .true.
     if(kfl_model_tur==30) TUR_K_OMEGA              = .true. 
     if(kfl_model_tur==31) TUR_K_OMEGA_BREDBERG     = .true. 
     if(kfl_model_tur==32) TUR_SST_K_OMEGA          = .true. 

     if(kfl_model_tur>=10.and.kfl_model_tur<=29) TUR_FAMILY_K_EPS   = .true.
     if(kfl_model_tur>=30.and.kfl_model_tur<=39) TUR_FAMILY_K_OMEGA = .true.

     if(kfl_model_tur<10) then
        TUR_ONE_EQUATION_MODEL = .true.
     else
        TUR_TWO_EQUATION_MODEL = .true.
     end if

     !-------------------------------------------------------------------
     !
     ! Models constants
     !
     !-------------------------------------------------------------------

     if(TUR_SPALART_ALLMARAS) then
        !
        ! Spalart-Allmars model
        !
        param_tur(1)=0.1355_rp                 ! cb1 
        param_tur(2)=0.622_rp                  ! cb2
        param_tur(3)=0.6670_rp                 ! sig 
        param_tur(4)=7.100_rp                  ! cv1
        param_tur(5)=0.3000_rp                 ! cw2 
        param_tur(6)=2.000_rp                  ! cw3
        param_tur(7)=0.4100_rp                 ! kap
        param_tur(8)=param_tur(1)&             ! cw1=cb1/kap^2+(1+cb2)/sig
             /(param_tur(7)**2)&       
             +(1.0_rp+param_tur(2))&
             /param_tur(3)
        ipara_tur(1)=1                         ! Geuzaine, Delanaye and Liu corrections
        ipara_tur(2)=1                         ! Rotation term

     else if(TUR_K_XU_CHIEN) then
        !
        ! Xu-Chen natural convection
        !
        param_tur(1)=1.00_rp                   ! sig_k 
        param_tur(2)=1.30_rp                   ! sig_e
        param_tur(3)=1.44_rp                   ! Ce1 
        param_tur(4)=1.92_rp                   ! Ce2
        param_tur(5)=1.44_rp                   ! Ce3 
        param_tur(6)=0.09_rp                   ! Cmu 
        param_tur(7)=0.09_rp                   ! Cmu 
        ipara_tur(1)=1                         ! los de chen
        ipara_tur(2)=1                         ! ...


     else if(kfl_model_tur==10) then
        !
        ! Chen k-l
        !

     else if(TUR_K_EPS_STD) then
        !
        ! Standard k-e
        !
        !        if( ipara_tur(1) == 0 ) then
        if (abs(param_tur(1)).lt.1.0e-4_rp)    param_tur(1)=1.00_rp                ! sig_k 
        if (abs(param_tur(2)).lt.1.0e-4_rp)    param_tur(2)=1.30_rp                ! sig_e
        if (abs(param_tur(3)).lt.1.0e-4_rp)    param_tur(3)=1.44_rp                ! Ce1 
        if (abs(param_tur(4)).lt.1.0e-4_rp)    param_tur(4)=1.92_rp                ! Ce2
        if (abs(param_tur(5)).lt.1.0e-4_rp)    param_tur(5)=1.44_rp                ! Ce3 
        if (abs(param_tur(6)).lt.1.0e-4_rp)    param_tur(6)=0.09_rp                ! Cmu 
        if (abs(param_tur(7)).lt.1.0e-4_rp)    param_tur(7)=param_tur(6)           ! beta*=Cmu

 !       else if( ipara_tur(1) == 1 ) then
           !
           ! Neutral ABL model
           !
 !          param_tur(1)=1.00_rp                ! sig_k 
 !          param_tur(2)=1.30_rp                ! sig_e
 !          param_tur(3)=1.176_rp               ! Ce1 
 !          param_tur(4)=1.92_rp                ! Ce2
 !          param_tur(6)=0.0333_rp               ! Cmu 
 !          param_tur(7)=param_tur(6)           ! beta*=Cmu
 !          param_tur(2)=0.41_rp*0.41_rp&
 !               /((param_tur(4)-param_tur(3))*sqrt(param_tur(6)))

 !       else if( ipara_tur(1) == 2 ) then
           !
           ! Non-neutral ABL model
           !
 !          param_tur(1)=0.74_rp                ! sig_k 
 !          param_tur(2)=1.30_rp                ! sig_e
 !          param_tur(3)=1.13_rp                ! Ce1 
 !          param_tur(4)=1.90_rp                ! Ce2
 !          param_tur(6)=0.0256_rp              ! Cmu 
 !          param_tur(7)=param_tur(6)           ! beta*=Cmu
 
 !       end if

     else if(TUR_K_EPS_LAUNDER_SHARMA) then
        !
        ! Launder-Sharma k-e
        !
        param_tur(1)=1.00_rp                   ! sig_k 
        param_tur(2)=1.30_rp                   ! sig_e
        param_tur(3)=1.44_rp                   ! Ce1 
        param_tur(4)=1.92_rp                   ! Ce2
        param_tur(5)=0.80_rp                   ! Ce3 
        param_tur(6)=0.09_rp                   ! Cmu 
        param_tur(7)=param_tur(6)              ! beta*=Cmu

     else if(TUR_K_EPS_CHIEN) then
        !
        ! Chien k-e
        !
        param_tur(1)=1.00_rp                   ! sig_k 
        param_tur(2)=1.30_rp                   ! sig_e
        param_tur(3)=1.35_rp                   ! Ce1 
        param_tur(4)=1.80_rp                   ! Ce2
        param_tur(5)=0.80_rp                   ! Ce3
        param_tur(6)=0.09_rp                   ! Cmu 
        param_tur(7)=param_tur(6)              ! beta*=Cmu

     else if(TUR_K_EPS_LAM_BREMHORST) then
        !
        ! Lam-Bremhorst k-e
        !
        param_tur(1)=1.00_rp                   ! sig_k 
        param_tur(2)=1.30_rp                   ! sig_e
        param_tur(3)=1.44_rp                   ! Ce1 
        param_tur(4)=1.92_rp                   ! Ce2
        param_tur(6)=0.09_rp                   ! Cmu 
        param_tur(7)=param_tur(6)              ! beta*=Cmu

     else if(TUR_K_EPS_JAW_HWANG) then
        !
        ! Jaw-Hwang k-e
        !
        param_tur(1)=1.00_rp                   ! sig_k 
        param_tur(2)=1.2857_rp                 ! sig_e
        param_tur(3)=1.44_rp                   ! Ce1 
        param_tur(4)=1.92_rp                   ! Ce2
        param_tur(5)=0.80_rp                   ! Ce3
        param_tur(6)=0.09_rp                   ! CD
        param_tur(7)=0.09_rp                   ! CD=Cmu

     else if(TUR_K_EPS_NAGANO_TAGAWA) then
        !
        ! Nagano-Tagawa k-e
        !
        param_tur(1)=1.40_rp                   ! sig_k 
        param_tur(2)=1.30_rp                   ! sig_e
        param_tur(3)=1.45_rp                   ! Ce1 
        param_tur(4)=1.90_rp                   ! Ce2
        param_tur(5)=0.80_rp                   ! Ce3
        param_tur(6)=0.09_rp                   ! Cmu
        param_tur(7)=0.09_rp                   ! Cmu

     else if(TUR_K_EPS_V2_F) then
        !
        ! Lien-Durbin k-eps-v2-f
        !
        param_tur(1)=1.00_rp                   ! sig_k 
        param_tur(2)=1.30_rp                   ! sig_e
        param_tur(3)=1.00_rp                   ! sig_phi
        param_tur(6)=0.22_rp                   ! Cmu
        param_tur(7)=0.22_rp                   ! Cmu

     else if(TUR_K_EPS_PHI_F) then
        !
        ! k-eps-phi-f
        !
        param_tur(1)=1.00_rp                   ! sig_k 
        param_tur(2)=1.30_rp                   ! sig_e
        param_tur(3)=1.00_rp                   ! sig_phi
        param_tur(6)=0.22_rp                   ! Cmu
        param_tur(7)=0.22_rp                   ! Cmu

     else if(TUR_K_EPS_PHI_F) then
        !
        ! Laurence-Uribe k-eps-phi-f
        !
        param_tur(1) =1.00_rp                  ! sig_k 
        param_tur(2) =1.30_rp                  ! sig_e
        param_tur(3) =1.00_rp                  ! sig_phi
        param_tur(6) =0.22_rp                  ! Cmu
        param_tur(7) =0.22_rp                  ! Cmu
        param_tur(8) =1.90_rp                  ! Ce2
        param_tur(9) =0.80_rp                  ! Ce3
        param_tur(10)=1.40_rp                  ! c1
        param_tur(11)=0.30_rp                  ! c2
        param_tur(12)=0.25_rp                  ! CL
        param_tur(13)=110.00_rp                ! Cn

     else if(TUR_TWO_LAYER_RODI) then
        !
        ! Rodi two-layer 
        !
        param_tur(1)=1.00_rp                   ! sig_k 
        param_tur(2)=1.30_rp                   ! sig_e
        param_tur(3)=1.44_rp                   ! Ce1 
        param_tur(4)=1.92_rp                   ! Ce2
        param_tur(5)=1.44_rp                   ! Ce3 
        param_tur(6)=0.09_rp                   ! Cmu 
        param_tur(7)=param_tur(6)              ! beta*=Cmu
        param_tur(8)=0.5_rp                    ! mixing slope
        param_tur(9)=80.0_rp                   ! mixing value

     else if(TUR_TWO_LAYER_XU_CHEN) then
        !
        ! Xu-Chen tow-layer natural convection
        !
        param_tur(1)=1.00_rp                   ! sig_k 
        param_tur(2)=1.30_rp                   ! sig_e
        param_tur(3)=1.44_rp                   ! Ce1 
        param_tur(4)=1.92_rp                   ! Ce2
        param_tur(5)=1.44_rp                   ! Ce3 
        param_tur(6)=0.09_rp                   ! Cmu 
        param_tur(7)=0.09_rp                   ! Cmu 
        param_tur(8)=0.5_rp                    ! mixing slope
        param_tur(9)=160.0_rp                  ! mixing value

     else if(TUR_K_OMEGA) then
        !
        ! Standard k-w
        !
        param_tur( 1)=2.0_rp                   ! sig_k
        param_tur( 2)=2.0_rp                   ! sig_w
        param_tur( 3)=5.0_rp/9.0_rp            ! alpha
        param_tur( 4)=3.0_rp/40.0_rp           ! beta
        param_tur( 5)=1.0_rp                   ! alpha*
        param_tur( 6)=9.0_rp/100.0_rp          ! beta*
        param_tur( 7)=6.0_rp                   ! Rk
        param_tur( 8)=27.0_rp/10.0_rp          ! Rw
        param_tur( 9)=8.0_rp                   ! Rbeta
        param_tur(10)=param_tur(4)/3.0_rp      ! alpha0*
        param_tur(11)=1.0_rp/10.0_rp           ! alpha0 

     else if(TUR_K_OMEGA_BREDBERG) then
        !
        ! Bredberg-Peng-Davidson k-w model
        !
        param_tur( 1)=1.00_rp                  ! sig_k
        param_tur( 2)=1.80_rp                  ! sig_w
        param_tur( 3)=0.09_rp                  ! Ck
        param_tur( 4)=1.10_rp                  ! Cw
        param_tur( 5)=0.49_rp                  ! Cw1
        param_tur( 6)=1.00_rp                  ! Cmu
        param_tur( 7)=0.072_rp                 ! Cw2
        param_tur( 9)=8.0_rp                   ! Rbeta

     else if(TUR_SST_K_OMEGA) then
        !
        ! Menter SST k-w
        !
        param_tur( 1)=2.0_rp                   ! sig_k
        param_tur( 2)=2.0_rp                   ! sig_w
        param_tur( 3)=5.0_rp/9.0_rp            ! alpha
        param_tur( 4)=3.0_rp/40.0_rp           ! beta
        param_tur( 5)=1.0_rp                   ! alpha*
        param_tur( 6)=9.0_rp/100.0_rp          ! beta*
        param_tur( 7)=6.0_rp                   ! Rk
        param_tur( 8)=27.0_rp/10.0_rp          ! Rw
        param_tur( 9)=8.0_rp                   ! Rbeta
        param_tur(10)=param_tur(4)/3.0_rp      ! alpha0*
        param_tur(11)=1.0_rp/10.0_rp           ! alpha0    
        ipara_tur(1)=1                         ! Menter SST-2003 version

     end if
     !
     ! Solver fixity
     !
     if( INOTMASTER ) then
        do iunkn = 1,nturb_tur
           solve(iunkn) % bvess     => bvess_tur(:,:,iunkn)
           solve(iunkn) % kfl_fixno => kfl_fixno_tur(:,:,iunkn)
        end do
     end if

  end select

end subroutine tur_inivar
