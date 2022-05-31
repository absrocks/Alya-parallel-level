!------------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    nsi_inivar.f90
!> @author  Guillaume Houzeaux
!> @brief   Initialize some variables
!> @details Initialize some variables\n
!!          ITASK=1 ... When starting the run (from Turnon)\n
!!          ITASK=2 ... First time step. This is needed as some variables
!!                      are not initialized before\n
!!          ITASK=3 ... When starting a time step (from nsi_begste)\n
!> @}
!------------------------------------------------------------------------
subroutine nsi_inivar(itask)
  use def_parame
  use def_master
  use def_nastin
  use def_domain
  use def_solver
  use mod_memchk
  use def_kermod,               only : gasco, kfl_adj_prob
  use mod_nsi_schur_complement, only : nsi_schur_complement_initialization
  use mod_nsi_multi_step_fs,    only : nsi_multi_step_fs_initialization
  use mod_nsi_fractional_step,  only : nsi_fractional_step_initialization
  use mod_arrays,               only : arrays_register
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ivar1,ivar2

  select case(itask)

  case(0_ip)
     !
     ! Modules initialization
     !
     call nsi_schur_complement_initialization()
     call nsi_multi_step_fs_initialization()
     call nsi_fractional_step_initialization()
     !
     ! Postprocess Variable names and types
     !
     call arrays_register(  1_ip,(/'VELOC','VECTO','NPOIN','PRIMA'/),veloc,ENTITY_POSITION=2_ip,TIME_POSITION=3_ip,EXCLUDE_RST=(/2_ip,3_ip/))
     call arrays_register(  2_ip,(/'PRESS','SCALA','NPOIN','PRIMA'/),press,ENTITY_POSITION=1_ip,TIME_POSITION=2_ip,EXCLUDE_RST=(/2_ip,3_ip/))
     call arrays_register(  3_ip,(/'STREA','SCALA','NPOIN','SECON'/))      
     call arrays_register(  4_ip,(/'RESID','SCALA','NPOIN','SECON'/))      
     call arrays_register(  5_ip,(/'VESGS','R3PVE','NELEM','PRIMA'/),vesgs,ENTITY_POSITION=1_ip,TIME_POSITION=3_ip,EXCLUDE_RST=(/2_ip/))
     call arrays_register(  6_ip,(/'BOUND','VECTO','NPOIN','SECON'/))      
     call arrays_register(  7_ip,(/'DENSI','SCALA','NPOIN','PRIMA'/),densi,ENTITY_POSITION=1_ip,TIME_POSITION=2_ip)
     call arrays_register(  8_ip,(/'PMV  ','SCALA','NPOIN','SECON'/))
     call arrays_register(  9_ip,(/'PPD  ','SCALA','NPOIN','SECON'/))
     call arrays_register( 10_ip,(/'MACHN','SCALA','NPOIN','SECON'/))
     call arrays_register( 11_ip,(/'TANGE','VECTO','NPOIN','SECON'/))
     call arrays_register( 12_ip,(/'VORTI','VECTO','NPOIN','SECON'/))
     call arrays_register( 13_ip,(/'MODVO','SCALA','NPOIN','SECON'/))
     call arrays_register( 14_ip,(/'LAMB2','SCALA','NPOIN','SECON'/))
     call arrays_register( 15_ip,(/'GROUP','SCALA','NPOIN','SECON'/))
     call arrays_register( 16_ip,(/'LINEL','SCALA','NPOIN','SECON'/))
     call arrays_register( 17_ip,(/'VISCO','SCALA','NPOIN','SECON'/))
     call arrays_register( 18_ip,(/'FIXPR','SCALA','NPOIN','SECON'/))
     call arrays_register( 19_ip,(/'FIXNO','VECTO','NPOIN','SECON'/))
     call arrays_register( 20_ip,(/'TAU  ','SCALA','NPOIN','SECON'/))                                        ! FREE POSITION
     call arrays_register( 21_ip,(/'AVVEL','VECTO','NPOIN','PRIMA'/),avvel_nsi,ENTITY_POSITION=2_ip)
     call arrays_register( 22_ip,(/'SCHUR','SCALA','NPOIN','SECON'/))         
     call arrays_register( 23_ip,(/'VEPRO','VECTO','NPOIN','PRIMA'/),vepro_nsi,ENTITY_POSITION=2_ip)
     call arrays_register( 24_ip,(/'PRPRO','SCALA','NPOIN','PRIMA'/),prpro_nsi,ENTITY_POSITION=1_ip)
     call arrays_register( 25_ip,(/'YPLUS','SCALA','NPOIN','SECON'/))         
     call arrays_register( 26_ip,(/'LIMIT','SCALA','NPOIN','SECON'/))         
     call arrays_register( 27_ip,(/'GRPRO','VECTO','NPOIN','PRIMA'/),grpro_nsi,ENTITY_POSITION=2_ip)
     call arrays_register( 28_ip,(/'AVPRE','SCALA','NPOIN','PRIMA'/),avpre_nsi,ENTITY_POSITION=1_ip)
     call arrays_register( 29_ip,(/'VEERR','SCALA','NPOIN','SECON'/))
     call arrays_register( 30_ip,(/'PECLE','SCALA','NPOIN','SECON'/))
     call arrays_register( 31_ip,(/'HYDRO','SCALA','NPOIN','SECON'/))
     call arrays_register( 32_ip,(/'DUDT ','VECTO','NPOIN','SECON'/))
     call arrays_register( 33_ip,(/'D2UDT','VECTO','NPOIN','SECON'/))
     call arrays_register( 34_ip,(/'DT1  ','SCALA','NPOIN','SECON'/))
     call arrays_register( 35_ip,(/'DT2  ','SCALA','NPOIN','SECON'/))
     call arrays_register( 36_ip,(/'PRERR','SCALA','NPOIN','SECON'/))
     call arrays_register( 37_ip,(/'LINVE','SCALA','NPOIN','SECON'/))
     call arrays_register( 38_ip,(/'DEFOR','VECTO','NPOIN','SECON'/))
     call arrays_register( 39_ip,(/'NODPR','SCALA','NPOIN','SECON'/))
     call arrays_register( 40_ip,(/'TTRAC','VECTO','NPOIN','SECON'/))                                       ! total traction                              
     call arrays_register( 41_ip,(/'VELO2','VECTO','NPOIN','SECON'/))                                       ! For BDF restart   
     call arrays_register( 42_ip,(/'PRES2','SCALA','NPOIN','SECON'/))                                       ! For BDF restart   
     call arrays_register( 43_ip,(/'VELO3','VECTO','NPOIN','SECON'/))                                       ! For BDF restart   
     call arrays_register( 44_ip,(/'PRES3','SCALA','NPOIN','SECON'/))                                       ! For BDF restart   
     call arrays_register( 45_ip,(/'MODVE','SCALA','NPOIN','SECON'/))                                                                                           
     call arrays_register( 46_ip,(/'AVTAN','VECTO','NPOIN','PRIMA'/),avtan_nsi,ENTITY_POSITION=2_ip)        ! Average tangent component of boundary traction
     call arrays_register( 47_ip,(/'LINVE','SCALA','NPOIN','SECON'/))                                                                                           
     call arrays_register( 48_ip,(/'VELOM','VECTO','NPOIN','SECON'/))                                                                                           
     call arrays_register( 49_ip,(/'FORCF','VECTO','NPOIN','PRIMA'/),forcf,ENTITY_POSITION=2_ip)                                                                      
     call arrays_register( 50_ip,(/'INTFO','VECTO','NPOIN','SECON'/))                                                                                           
     call arrays_register( 51_ip,(/'USTAR','SCALA','NPOIN','SECON'/))                                                                                           
     call arrays_register( 52_ip,(/'TRACT','VECTO','NPOIN','SECON'/))                                       ! traction                                                        
     call arrays_register( 53_ip,(/'FIXRS','SCALA','NPOIN','SECON'/))                                                                                 
     call arrays_register( 54_ip,(/'FVELO','VECTO','NPOIN','SECON'/))                                       ! new velocity for filter output             
     call arrays_register( 55_ip,(/'FPRES','SCALA','NPOIN','SECON'/))                                       ! new pressure for filter output             
     call arrays_register( 56_ip,(/'FTANG','VECTO','NPOIN','SECON'/))                                       ! new wall shear stress for filter output    
     call arrays_register( 57_ip,(/'AVVE2','VECTO','NPOIN','PRIMA'/),avve2_nsi,ENTITY_POSITION=2_ip)        ! Average Vi**2                              
     call arrays_register( 58_ip,(/'AVVXY','VECTO','NPOIN','PRIMA'/),avvxy_nsi,ENTITY_POSITION=2_ip)        ! Average Vx*Vy, Vy*Vz, Vz*Vx                
     call arrays_register( 59_ip,(/'AVPR2','SCALA','NPOIN','PRIMA'/),avpr2_nsi,ENTITY_POSITION=1_ip)        ! Average Pr**2                              
     call arrays_register( 60_ip,(/'MESHR','VECTO','NPOIN','SECON'/))                                       ! Mesh Rotation                              
     call arrays_register( 61_ip,(/'NODEF','SCALA','NPOIN','SECON'/))                                       ! If reaction force are computed on nodes    
     call arrays_register( 62_ip,(/'MOMEN','VECTO','NPOIN','SECON'/))                                       ! Momentum equaiton residual                 
     call arrays_register( 63_ip,(/'FIXPP','SCALA','NPOIN','SECON'/))                                       ! Pressure fixity                            
     call arrays_register( 64_ip,(/'BVNAT','SCALA','NPOIN','SECON'/))                                       ! Algebraic Neumann condition                
     call arrays_register( 65_ip,(/'REACT','VECTO','NPOIN','SECON'/))                                       ! Algebraic reaction (b-Ax)                  
     call arrays_register( 66_ip,(/'CHA01','VECTO','NPOIN','SECON'/))                                                                             
     call arrays_register( 67_ip,(/'CHA02','VECTO','NPOIN','SECON'/))                                                                             
     call arrays_register( 68_ip,(/'CHA03','VECTO','NPOIN','SECON'/))                                                                             
     call arrays_register( 69_ip,(/'CHA04','VECTO','NPOIN','SECON'/))                                                                             
     call arrays_register( 70_ip,(/'CHA05','VECTO','NPOIN','SECON'/))                                                                             
     call arrays_register( 71_ip,(/'WETNO','SCALA','NPOIN','SECON'/))                                                                             
     call arrays_register( 72_ip,(/'AVVRE','VECTO','WHATE','SECON'/))                                       ! For time-averaging restart              
     call arrays_register( 74_ip,(/'TURMU','SCALA','NPOIN','SECON'/))                                       ! LES turbulent viscosity                 
     call arrays_register( 75_ip,(/'AVMUT','SCALA','NPOIN','PRIMA'/),avmut_nsi,ENTITY_POSITION=1_ip)        ! Average turbulent viscosity             
     call arrays_register( 76_ip,(/'AVSTX','VECTO','NPOIN','PRIMA'/),avstx_nsi,ENTITY_POSITION=2_ip)        ! Average stress (mu_t * grad(u))         
     call arrays_register( 77_ip,(/'AVSTY','VECTO','NPOIN','PRIMA'/),avsty_nsi,ENTITY_POSITION=2_ip)        ! Average stress (mu_t * grad(v))         
     call arrays_register( 78_ip,(/'AVSTZ','VECTO','NPOIN','PRIMA'/),avstz_nsi,ENTITY_POSITION=2_ip)        ! Average stress (mu_t * grad(w))         
     call arrays_register( 79_ip,(/'BUBBL','SCALA','NELEM','SECON'/))                                       ! Pressure bubble                         
     call arrays_register( 80_ip,(/'TAUPR','SCALA','NPOIN','SECON'/))                                       ! Projection of tau                       
     call arrays_register( 81_ip,(/'MFCFO','SCALA','NPOIN','SECON'/))                                       ! For mass flow control restart           
     call arrays_register( 82_ip,(/'UBPRE','SCALA','NPOIN','SECON'/))                                       ! For mass flow control restart           
     call arrays_register( 83_ip,(/'NORMA','SCALA','NPOIN','SECON'/))                                       ! Normals                                 
     call arrays_register( 84_ip,(/'SENSM','VECTO','NPOIN','SECON'/))                                       ! Mesh sensitivities                      
     call arrays_register( 85_ip,(/'LAGRA','VECTO','NPOIN','PRIMA'/),lagra_nsi,ENTITY_POSITION=2_ip)        ! Lagrange multiplier                     
     call arrays_register( 86_ip,(/'ENVEL','VECTO','NPOIN','PRIMA'/),envel_nsi,ENTITY_POSITION=2_ip)                                                                   
     call arrays_register( 87_ip,(/'ENPRE','SCALA','NPOIN','PRIMA'/),enpre_nsi,ENTITY_POSITION=1_ip)                                                                   
     call arrays_register( 88_ip,(/'ENVE2','VECTO','NPOIN','PRIMA'/),enve2_nsi,ENTITY_POSITION=2_ip)        ! Ensemble Vi**2                                               
     call arrays_register( 89_ip,(/'ENVXY','VECTO','NPOIN','PRIMA'/),envxy_nsi,ENTITY_POSITION=2_ip)        ! Ensemble Vx*Vy, Vy*Vz, Vz*Vx                                 
     call arrays_register( 90_ip,(/'ENPR2','SCALA','NPOIN','PRIMA'/),enpr2_nsi,ENTITY_POSITION=1_ip)        ! Ensemble Pr**2                                               
     call arrays_register( 91_ip,(/'ENTAN','VECTO','NPOIN','PRIMA'/),entan_nsi,ENTITY_POSITION=2_ip)        ! Ensemble tau                                                 
     call arrays_register( 92_ip,(/'ENMUT','SCALA','NPOIN','PRIMA'/),enmut_nsi,ENTITY_POSITION=1_ip)        ! Ensemble turbulent viscosity                                 
     call arrays_register( 93_ip,(/'ENSTX','VECTO','NPOIN','PRIMA'/),enstx_nsi,ENTITY_POSITION=2_ip)        ! Ensemble stress (mu_t * grad(u))                             
     call arrays_register( 94_ip,(/'ENSTY','VECTO','NPOIN','PRIMA'/),ensty_nsi,ENTITY_POSITION=2_ip)        ! Ensemble stress (mu_t * grad(v))                             
     call arrays_register( 95_ip,(/'ENSTZ','VECTO','NPOIN','PRIMA'/),enstz_nsi,ENTITY_POSITION=2_ip)        ! Ensemble stress (mu_t * grad(w))                             
     call arrays_register( 97_ip,(/'BTRAC','VECTO','NPOIN','PRIMA'/),btrac_nsi,ENTITY_POSITION=2_ip)        ! For two-layer model restart - LES part                       
     call arrays_register( 98_ip,(/'TRACR','VECTO','NPOIN','PRIMA'/),tluav_nsi,ENTITY_POSITION=2_ip)        ! For two-layer model restart - RANS part                      
     call arrays_register( 99_ip,(/'TLUAV','VECTO','NPOIN','PRIMA'/),tracr_nsi,ENTITY_POSITION=2_ip)        ! For two-layer model restart - Average velocity               
     call arrays_register(100_ip,(/'AVVTA','VECTO','WHATE','SECON'/))                                       ! Average viscous part of tangential stres
     call arrays_register(101_ip,(/'VAFOR','VECTO','NPOIN','PRIMA'/),vafor_nsi,ENTITY_POSITION=2_ip)        ! VARIATIONAL FORCE                                             
     call arrays_register(102_ip,(/'AVVAF','VECTO','NPOIN','PRIMA'/),avvaf_nsi,ENTITY_POSITION=2_ip)        ! Average VARIATIONAL FORCE                                     
     call arrays_register(103_ip,(/'NOTRA','VECTO','NPOIN','PRIMA'/),notra_nsi,ENTITY_POSITION=2_ip)        ! VARIATIONAL Tangential traction                               
     call arrays_register(104_ip,(/'AVNTR','VECTO','NPOIN','PRIMA'/),avntr_nsi,ENTITY_POSITION=2_ip)        ! Average VARIATIONAL Tangential traction  - beware running average             
     call arrays_register(105_ip,(/'AVGTR','VECTO','NPOIN','PRIMA'/),avgtr_nsi,ENTITY_POSITION=2_ip)        ! Average gradient based Tangential traction  - beware running average
     call arrays_register(106_ip,(/'FANSW','SCALA','NPOIN','SECON'/))                                       ! Factor for no slip wall law - 4 restart 
     call arrays_register(107_ip,(/'PORFO','VECTO','NPOIN','SECON'/))                                       ! Porous force  
     call arrays_register(108_ip,(/'MASRH','SCALA','NPOIN','PRIMA'/),mass_rho_nsi,ENTITY_POSITION=1_ip,TIME_POSITION=2_ip) ! M * rho for low-Mach                                     
     call arrays_register(109_ip,(/'MOMSK','VECTO','NPOIN','SECON'/))                                       ! Momentum source from partis               
     call arrays_register(110_ip,(/'DRHOD','SCALA','NPOIN','SECON'/))                                       ! unlumped drho/dt for restarting low mach  
     call arrays_register(111_ip,(/'VEOLD','VECTO','NPOIN','PRIMA'/),veold_nsi,     ENTITY_POSITION=2_ip)   ! Old velocity        
     call arrays_register(112_ip,(/'TAUIB','VECTO','NPOIN','PRIMA'/),tauib_nsi,     ENTITY_POSITION=3_ip)   ! tauib_nsi          
     call arrays_register(113_ip,(/'DUNKN','VECTO','NPOIN','PRIMA'/),dunkn_nsi,     ENTITY_POSITION=2_ip)   ! dunkn_nsi    
     call arrays_register(114_ip,(/'DUNKP','SCALA','NPOIN','PRIMA'/),dunkp_nsi,     ENTITY_POSITION=1_ip)   ! dunkp_nsi                                                               
     call arrays_register(115_ip,(/'BUBBL','SCALA','NELEM','PRIMA'/),bubble_nsi    ,ENTITY_POSITION=1_ip)   ! Bubble
     call arrays_register(116_ip,(/'BUAQQ','SCALA','NELEM','PRIMA'/),bubble_aqq_nsi,ENTITY_POSITION=1_ip)   ! Bubble 
     call arrays_register(117_ip,(/'BUAQU','MATRI','NELEM','PRIMA'/),bubble_aqu_nsi,ENTITY_POSITION=2_ip)   ! Bubble 
     call arrays_register(118_ip,(/'BUAQP','MATRI','NELEM','PRIMA'/),bubble_aqp_nsi,ENTITY_POSITION=2_ip)   ! Bubble 
     call arrays_register(119_ip,(/'BUBQ ','SCALA','NELEM','PRIMA'/),bubble_bq_nsi ,ENTITY_POSITION=1_ip)   ! Bubble
     call arrays_register(120_ip,(/'VEPAV','VECTO','NPOIN','SECON'/))                                       ! Velocity at boundaries projected on nodes
     call arrays_register(121_ip,(/'QCRIT','SCALA','NPOIN','SECON'/))                                       ! Q criterion
     call arrays_register(122_ip,(/'AVMOS','VECTO','NPOIN','PRIMA'/),avmos_nsi,ENTITY_POSITION=2_ip)        ! Average momentum source from spray
     !                                                         
     ! Set variables
     !
     postp(1) % woese (1)     = 'VELOC'
     postp(1) % woese (2)     = 'VORTI'
     postp(1) % woese (3)     = 'KINET'
     postp(1) % woese (4)     = 'DIVU2'
     postp(1) % woese (5:7)   = 'MOMEN'

     postp(1) % wobse (1)     = 'MEANP'   ! Mean pressure
     postp(1) % wobse (2)     = 'MASS '   ! Mass rho*u.n
     postp(1) % wobse (3)     = 'FORCE'   ! Viscous force
     postp(1) % wobse (4)     = 'F_v_y'
     postp(1) % wobse (5)     = 'F_v_z'
     postp(1) % wobse (6)     = 'F_p_x'   ! Pressure force
     postp(1) % wobse (7)     = 'F_p_y'
     postp(1) % wobse (8)     = 'F_p_z'
     postp(1) % wobse (9)     = 'TORQU'   ! Viscous torque
     postp(1) % wobse (10)    = 'T_v_y'
     postp(1) % wobse (11)    = 'T_v_z'
     postp(1) % wobse (12)    = 'T_p_x'   ! Pressure torque
     postp(1) % wobse (13)    = 'T_p_y'
     postp(1) % wobse (14)    = 'T_p_z'
     postp(1) % wobse (15)    = 'MEANY'   ! Mean y+
     postp(1) % wobse (16)    = 'MEANV'   ! Mean velocity
     postp(1) % wobse (17)    = 'WETFO'   ! Viscous wet force
     postp(1) % wobse (18)    = 'Fwv_y'
     postp(1) % wobse (19)    = 'Fwv_z'
     postp(1) % wobse (20)    = 'Fwp_x'   ! Pressure wet force
     postp(1) % wobse (21)    = 'Fwp_y'
     postp(1) % wobse (22)    = 'Fwp_z'
     postp(1) % wobse (23)    = 'WETSU'
     postp(1) % wobse (24)    = 'INTFX'   ! Internal force
     postp(1) % wobse (25)    = 'INTFY'
     postp(1) % wobse (26)    = 'INTFZ'
     postp(1) % wobse (27)    = 'REATT'   ! Reattachment
     postp(1) % wobse (28)    = 'REATY'
     postp(1) % wobse (29)    = 'REATZ'
     postp(1) % wobse (30)    = 'REATM'
     postp(1) % wobse (31)    = 'UDOTN'   ! u.n
     postp(1) % wobse (32)    = 'FPORX'   ! Porous force -- ojo no meti como si fuera boundary pero en realidad es volumetrico
     postp(1) % wobse (33)    = 'FPORY'
     postp(1) % wobse (34)    = 'FPORZ'

     postp(1) % wonse (1)     = 'VELOX'   ! x-velocity
     postp(1) % wonse (2)     = 'VELOY'   ! y-velocity
     postp(1) % wonse (3)     = 'VELOZ'   ! z-velocity
     postp(1) % wonse (4)     = 'PRESS'   ! Pressure
     postp(1) % wonse (5)     = 'YPLUS'   ! y+
     postp(1) % wonse (6)     = 'PMV  '   ! PMV - comfort index
     postp(1) % wonse (7)     = 'PPD  '   ! PPD - comfort index
     !
     ! Witness variables
     !
     postp(1) % wowit ( 1)    = 'VELOX'   ! x-velocity
     postp(1) % wowit ( 2)    = 'VELOY'   ! y-velocity
     postp(1) % wowit ( 3)    = 'VELOZ'   ! z-velocity
     postp(1) % wowit ( 4)    = 'PRESS'   ! Pressure
     postp(1) % wowit ( 5)    = 'S11  '   ! Strain rate
     postp(1) % wowit ( 6)    = 'S22  '   ! Strain rate
     postp(1) % wowit ( 7)    = 'S12  '   ! Strain rate
     postp(1) % wowit ( 8)    = 'S33  '   ! Strain rate
     postp(1) % wowit ( 9)    = 'S13  '   ! Strain rate
     postp(1) % wowit (10)    = 'S23  '   ! Strain rate
     postp(1) % wowit (11)    = 'COORX'   ! x-coordinate
     postp(1) % wowit (12)    = 'COORY'   ! y-coordinate
     postp(1) % wowit (13)    = 'COORZ'   ! z-coordinate
     postp(1) % wowit (14)    = 'TURBU'   ! Turbulent viscosity
     !
     ! Solvers
     !
     call soldef(-9_ip)                   ! Allocate memory
     solve(1) % kfl_solve = 1             ! Momentum:   Output flag
     solve(2) % kfl_solve = 1             ! Continuity: Output flag
     !
     ! Pointers: matrices used in Schur complement methods
     !
     Auu_nsi => nulir
     Aup_nsi => nulir
     Apu_nsi => nulir
     App_nsi => nulir
     !
     ! Nullify pointers
     !
     nullify(lforc_material_nsi)
     nullify(xforc_material_nsi)
     nullify(velta_nsi)
     nullify(thrta_nsi)
     nullify(powta_nsi)
     nullify(veave_nsi)
     nullify(radiu_nsi)
     nullify(forcn_nsi)
     nullify(forct_nsi)
     nullify(kfl_fixno_nsi)
     nullify(kfl_fixpr_nsi)
     nullify(kfl_fixpp_nsi)
     nullify(kfl_fixbo_nsi)
     nullify(kfl_fixrs_nsi)
     nullify(kfl_funno_nsi)
     nullify(kfl_funbo_nsi)
     nullify(kfl_funtn_nsi)
     nullify(kfl_funtb_nsi)
     nullify(kfl_wlawf_nsi)
     nullify(bvess_nsi)
     nullify(bvnat_nsi)
     nullify(skcos_nsi)
     nullify(tncod_nsi)
     nullify(tgcod_nsi)
     nullify(tbcod_nsi)
     nullify(veold_nsi)
     nullify(gradv_nsi)
     nullify(unk2n_nsi)
     nullify(bpess_nsi)
     nullify(dunkn_nsi)
     nullify(dunkp_nsi)
     nullify(avvel_nsi)
     nullify(avve2_nsi)
     nullify(avvxy_nsi)
     nullify(avpre_nsi)
     nullify(avpr2_nsi)
     nullify(avtan_nsi)
     nullify(resch_nsi)
     nullify(prope_nsi)
     nullify(norle_nsi)
     nullify(curle_nsi)
     nullify(outflow_mass)
     nullify(intfo_nsi)
     nullify(iboun_nsi)
     nullify(Auu_nsi)
     nullify(Aup_nsi)
     nullify(Apu_nsi)
     nullify(App_nsi)
     nullify(amatr_nsi)
     nullify(lapla_nsi)
     nullify(visco_nsi)
     nullify(cmama_nsi)
     nullify(deltp_nsi)
     nullify(dt_rho_nsi)
     nullify(mass_rho_nsi)
     nullify(tau_nsi)
     nullify(bubble_nsi)
     nullify(bubble_aqq_nsi)
     nullify(bubble_aqu_nsi)
     nullify(bubble_aqp_nsi)
     nullify(bubble_bq_nsi)
     nullify(vepro_nsi)
     nullify(grpro_nsi)
     nullify(prpro_nsi)
     nullify(vepr2_nsi)
     nullify(grpr2_nsi)
     nullify(prpr2_nsi)
     nullify(kfl_fixno_div_nsi)
     nullify(notra_nsi)
     nullify(massb_nsi)
     nullify(avmut_nsi)
     nullify(avstx_nsi)
     nullify(avsty_nsi)
     nullify(avstz_nsi)
     nullify(envel_nsi)
     nullify(enve2_nsi)
     nullify(envxy_nsi)
     nullify(enpre_nsi)
     nullify(enpr2_nsi)
     nullify(entan_nsi)
     nullify(enmut_nsi)
     nullify(enstx_nsi)
     nullify(ensty_nsi)
     nullify(enstz_nsi)
     nullify(btrac_nsi)
     nullify(tracr_nsi)
     nullify(tluav_nsi)
     nullify(vafor_nsi)
     nullify(avvaf_nsi)
     nullify(bupor_nsi)
     nullify(avntr_nsi)
     nullify(avgtr_nsi)
     nullify(drhodt_nsi)
     nullify(avmos_nsi)

     nullify(resdiff_nsi)
     nullify(dcost_dx_nsi)
     nullify(costdiff_nsi)
     
     nullify(lbpse)
     !
     ! Others
     !
     nodpr_nsi           = 0      ! Node where to impose pressure
     kfl_exist_fixi7_nsi = 0      ! exists nodal fixity of type 7
     kfl_exist_fib20_nsi = 0      ! exists boundary fixity of type 20
     pabdf_nsi           = 0.0_rp
     dtinv_nsi           = 0.0_rp

  case(1_ip)
     !
     ! Dimensions
     !
     ivari_nsi = ivari_nsi_mom
     !
     ! Number of velocity components
     !
     if(kfl_timei_nsi==1) then
        if(kfl_tisch_nsi==1) then
           !
           ! Trapezoidal rule
           !
           ncomp_nsi=3
        else if(kfl_tisch_nsi==2) then
           !
           ! BDF scheme
           !
           ncomp_nsi=2+kfl_tiacc_nsi
        else if(kfl_tisch_nsi==3 .or. kfl_tisch_nsi==4) then !=3 adams, =4 RK
           !
           ! Explicit time step
           !
           ncomp_nsi = 4
           ! AB WHAT IS THIS TODO TODO TODO if(kfl_tisch_nsi==4 .and. kfl_fscon_nsi == 0) then !RK extrapolating pressure
           ! AB WHAT IS THIS TODO TODO TODO    ncomp_nsi=4   ! only needed for the pressure unknown 
           ! AB WHAT IS THIS TODO TODO TODO else
           ! AB WHAT IS THIS TODO TODO TODO    ncomp_nsi=3 
           ! AB WHAT IS THIS TODO TODO TODO end if
        end if
        if ( kfl_ini_ts_guess_order_nsi > 1 ) then
           ncomp_nsi = max (ncomp_nsi,kfl_ini_ts_guess_order_nsi + 2 )
        end if
     else
        ncomp_nsi=2
     end if
     nprev_nsi = min(3_ip,ncomp_nsi)  ! Last time step or global iteration
     nbdfp_nsi = 2
     !
     ! First or second order velocity derivatives
     !
     ivar1 = 32
     call posdef(25_ip,ivar1)
     ivar2 = 33
     call posdef(25_ip,ivar2)
     if( kfl_dttyp_nsi == 2 .or. ivar1 /= 0 .or. ivar2 /= 0 ) then
        ncomp_nsi = max(5_ip,ncomp_nsi)
     end if
     !
     ! Time variables
     !
     if( kfl_timei_nsi == 0 ) then
        dtinv_nsi = 1.0_rp
     else
        kfl_timei = 1
     end if
     kfl_stead_nsi = 0
     !
     ! Dimensions
     !
     if( NSI_MONOLITHIC ) then
        solve(1) % ndofn = ndime + 1
        solve(1) % wprob = 'NAVIER_STOKES'
     else
        solve(1) % ndofn = ndime
        solve(2) % ndofn = 1
        solve(1) % wprob = 'MOMENTUM'
        solve(2) % wprob = 'CONTINUITY'
     end if
     if( NSI_SEMI_IMPLICIT ) then
        solve(9) % kfl_solve = 1
        solve(9) % ndofn     = ndime
        solve(9) % wprob     = 'VISCOUS_TERM'        
     end if
     !
     ! Multiple solves using the same matrix
     !
     if( NSI_SCHUR_COMPLEMENT ) then
        solve(1) % num_multiple_solves = 2
        if( kfl_sosch_nsi == 3 ) then
           solve(2) % num_multiple_solves = 2
        end if
     end if
     !
     ! Off-diagonal part of momentum equations
     !
     kfl_rmom2_nsi = 0
     kfl_p1ve2_nsi = 0

     if( fvnoa_nsi     >  0.0_rp )                          kfl_rmom2_nsi = 1
     if( kfl_visco_nsi == 1 .and. fvins_nsi     >= 0.9_rp ) kfl_rmom2_nsi = 1
     if( kfl_advec_nsi == 1 .and. kfl_linea_nsi == 2      ) kfl_rmom2_nsi = 1
     if( kfl_advec_nsi == 1 .and. kfl_adj_prob  == 1      ) kfl_rmom2_nsi = 1
     if( fvnoa_nsi     >  0.0_rp )                          kfl_p1ve2_nsi = 1
     if( kfl_anipo_nsi /= 0      )                          kfl_rmom2_nsi = 1
     
     ndbgs_nsi = ndime * npoin

     if( kfl_initi_nsi == 5 ) then
        solve(3)             = solve(1)                ! Check boundary conditions (same solver as Momentum)
        solve(3) % wprob     = 'BOUNDARY_CONDITIONS'
        solve(3) % ndofn     = ndime
        solve(3) % kfl_solve = 0
        if( kfl_initi_nsi == 5 ) then
           solve(3) % kfl_solve = 1
           solve(3) % kfl_cvgso = 1
        end if
     end if
     !
     ! Mass correction solver
     !
     solve(4) % wprob     = 'MASS_CORRECTION'
     if( kfl_corre_nsi == 3 ) then
        solve(4) % ndofn     = ndime
        solve(4) % kfl_solve = 1                        ! Output flag
     else
        solve(4) % kfl_algso = SOL_SOLVER_RICHARDSON    ! Mass correction: Diagonal solver
        solve(4) % ndofn     = ndime
        solve(4) % xdiag     = 1.0_rp
        if( kfl_corre_nsi == 1 ) then
           solve(4) % kfl_preco = SOL_CLOSE_MASS_MATRIX ! Uses VMASC
        else
           solve(4) % kfl_preco = SOL_MASS_MATRIX       ! Uses VMASS
        end if
        solve(4) % kfl_solve = 0
        solve(4) % kfl_cvgso = 0
     end if

     if( kfl_inipr_nsi == NSI_PDE_HYDROSTATIC_PRESSURE ) then
        solve(5) % wprob     = 'HYDROSTATIC_PRESSURE'
        solve(5) % ndofn     = 1
        solve(5) % kfl_solve = 1         
     end if

     if( kfl_divcorrec_nsi /= 0 ) then
        solve(6) % wprob     = 'ZERO_DIVERGENCE'
        solve(6) % ndofn     = 1
        solve(6) % kfl_solve = 1        
     end if

     solve(7) % kfl_solve = 0
     if( 1 == 0 ) then
        solve(7) % wprob     = 'NORMAL_EXTENSIONS'
        solve(7) % ndofn     = 1
        solve(7) % kfl_solve = 1 
     end if
     
     solve(8) % wprob     = 'CONSISTENT'                  ! Consistent matrix
     solve(8) % kfl_solve = 1
     solve(8) % kfl_iffix = 2
     solve(8) % ndofn     = ndime

     nevat_nsi          = (ndime+1)*mnode
     nschu_nsi          = 0                               ! # Schur complement solves
     nmome_nsi          = 0                               ! # Momentum solves

     kfl_perip_nsi      = 0                               ! Periodicity: if pressure prescribed on periodic nodes
     pcoef_nsi          = 1.0_rp-gasco/sphea_nsi
     gamth_nsi          = sphea_nsi/(sphea_nsi-gasco)     ! Low-Mach: gamma=Cp/(Cp-R)

     prthe(1)           = prthe_nsi                       ! Low-Mach: p0
     prthe(2)           = prthe_nsi                       ! Low-Mach: p0^{n}
     prthe(3)           = prthe_nsi                       ! Low-Mach: p0^{n-1}
     prthe(4)           = prthe_nsi                       ! Low-Mach: p0^{0}

     imass              = tmass_nsi                       ! initial density
     kfl_prthe          = kfl_prthe_nsi                   ! type of thpressure calculation
     tflux              = 0.0_rp                          ! Low-Mach: heat flux
     dpthe              = 0.0_rp                          ! Low-Mach: dp0/dt
     xmass_nsi          = 0.0_rp                          ! Low-Mach: mass
     cputi_nsi          = 0.0_rp
     cputi_assembly_nsi = 0.0_rp
     cpu_ass_sol_nsi    = 0.0_rp
     iteqn_nsi(1)       = 0
     iteqn_nsi(2)       = 0
     ittot_nsi          = 0
     resin_nsi(1)       = 1.0_rp                          ! Algebraic inner residual
     resin_nsi(2)       = 1.0_rp
     resou_nsi(1)       = 1.0_rp                          ! Algebraic outer residual
     resou_nsi(2)       = 1.0_rp
     resss_nsi(1)       = 1.0_rp                          ! Algebraic Steady state residual
     resss_nsi(2)       = 1.0_rp
     reinf_nsi(1)       = 1.0_rp                          ! Linf residual
     reinf_nsi(2)       = 1.0_rp
     relpa_nsi(1)       = 0.0_rp                          ! Relaxation
     relpa_nsi(2)       = 0.0_rp
     corio_nsi          = 0.0_rp                          ! Coriolis force module
     kfl_tiaor_nsi      = kfl_tiacc_nsi                   ! Time accuracy: save original value
     resgs_nsi(1)       = 0.0_rp
     resgs_nsi(2)       = 0.0_rp
     rmsgs_nsi          = 0.0_rp
     dtsgs_nsi          = 0.0_rp
     kfl_wlare_nsi      = 0
     porfo_nsi          = 0.0_rp
     dtmax_nsi          = 0.0_rp                          ! otherwise full ckeck stops because it is unitialized when doing max in nsi_tistep
     actav_nsi          = 0.0_rp                          ! Accumulated time for averaging
     if(kfl_sgsco_nsi==0) misgs_nsi=1                     ! Subgrid scale number of iterations
     if( kfl_sgsac_nsi == -1 ) then                       ! SGS time accuracy default
        if(kfl_tisch_nsi/=2) then
           kfl_sgsac_nsi = kfl_tiacc_nsi
        else
           kfl_sgsac_nsi = 1_ip                           ! For BDF use order 1 for the subscale by default
        end if
     end if
     !
     ! CFD wake
     !
     call nsi_cfdwak(1_ip)
#ifdef matiaslma
     !
     ! for low Mach, never solve pressure
     !
     if (kfl_regim_nsi==3) then
        kfl_confi_nsi = -1   ! do not fix pressure!!!
        fvins_nsi = 2.0_rp     ! Complete form
     end if
     !  pressure stabilization coefficient
     staco_nsi(4) = min(staco_nsi(4), 0.25_rp/staco_nsi(1))
#endif
     
  case(2_ip)
     !
     ! If temperature subgrid scale should be considered
     !
     if( associated(tesgs) .and. (kfl_cotem_nsi == 1 .or. kfl_regim_nsi == 3) ) then !!FER Check if SGS of temper is needed for low mach
        kfl_sgste_nsi = 1
     else
        kfl_sgste_nsi = 0
     end if
     !
     ! Activate flag for Low Mach in combustion (properties on GP or Node)
     !
     kfl_lookg_nsi=0
     if (kfl_regim_nsi==3 .and. (associated(tempe_gp) .or. kfl_coupl(ID_NASTIN,ID_CHEMIC) /= 0)) kfl_lookg_nsi = 1

  case(3_ip)
     !
     ! Before starting a time step
     !
     relax_nsi = 1.0_rp
     relap_nsi = 1.0_rp

  end select

end subroutine nsi_inivar
