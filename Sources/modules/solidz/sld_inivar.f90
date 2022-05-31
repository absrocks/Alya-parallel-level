!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_inivar.f90
!> @author  Mariano Vazquez
!> @date    August, 2006
!>          - Subroutine written
!> @brief   Initialize some variables
!>
!> @details
!>          \verbatim
!>          Initialize some variables:
!>          ITASK = 0 ... When starting the run (from Turnon)
!>                  1 ... When starting the run (from Turnon)
!>                  2 ... Correct fibers when needed
!>          \endverbatim
!>
!> @}
!-----------------------------------------------------------------------

subroutine sld_inivar(itask)

  use def_parame,         only : pi
  use def_kintyp,         only : ip,rp
  use def_master,         only : INOTMASTER,postp
  use def_master,         only : kfl_timei,coupling,solve
  use def_master,         only : displ
  use def_domain,         only : xfiel, lmate, ndime, nmate, nelem
  use def_solver,         only : SOL_LOCAL_MASS_MATRIX,SOL_SOLVER_RICHARDSON
  use mod_communications, only : PAR_MAX
  use mod_arrays,         only : arrays_register
  use def_solidz,         only : lmate_sld, nmate_sld, nvoig_sld, eleng_sld, celen_sld
  use def_solidz,         only : kfl_fixbo_sld, kfl_fixno_sld, kfl_fixrs_sld
  use def_solidz,         only : kfl_funbo_sld, kfl_funno_sld, kfl_funtn_sld, kfl_funtb_sld
  use def_solidz,         only : tbcod_sld, tncod_sld, tgcod_sld
  use def_solidz,         only : jacrot_du_dq_sld, jacrot_dq_du_sld
  use def_solidz,         only : bvess_sld, bvnat_sld
  use def_solidz,         only : veloc_sld, accel_sld, dunkn_sld, ddisp_sld
  use def_solidz,         only : uruk4_sld
  use def_solidz,         only : ndofn_sld, dtcri_sld
  use def_solidz,         only : fint2_sld, fext2_sld, fine2_sld, fnatu_sld
  use def_solidz,         only : finte_sld, fexte_sld, macce_sld, frxid_sld
  use def_solidz,         only : fintt_sld, fextt_sld
  use def_solidz,         only : allie_sld, allwk_sld, allke_sld, etota_sld
  use def_solidz,         only : epsel_sld, lepse_sld
  use def_solidz,         only : caunp_sld, caulo_sld, cause_sld
  use def_solidz,         only : svegm_sld, svaux_sld
  use def_solidz,         only : isoch_sld, iswav_sld
  use def_solidz,         only : nvgij_inv_sld, nvgij_sld, nsvar_sld
  use def_solidz,         only : rorig_sld
  use def_solidz,         only : vmass_sld
  use def_solidz,         only : kfl_xfeme_sld, vxfem_sld, dxfem_sld, axfem_sld
  use def_solidz,         only : cranx_sld, crapx_sld, lnenr_sld, leenr_sld
  use def_solidz,         only : lecoh_sld, treff_sld, cohnx_sld, nocoh_sld
  use def_solidz,         only : dcmax_sld, dceff_sld, dslip_sld
  use def_solidz,         only : nopio_sld
  use def_solidz,         only : iwave_sld, kfl_gdepo, kfl_sdvar_sld
  use def_solidz,         only : lawch_sld, lawst_sld, lawpl_sld
  use def_solidz,         only : relco_sld, relpa_sld
  use def_solidz,         only : kfl_restr_sld, restr_sld
  use def_solidz,         only : kfl_timei_sld, kfl_stead_sld, kfl_tisch_sld
  use def_solidz,         only : ncomp_sld, nprev_sld, nderi_sld, dtinv_sld
  use def_solidz,         only : kfl_moduf_sld, fiemo_sld
  use def_solidz,         only : axis1_sld, axis2_sld, axis3_sld, orien_sld
  use def_solidz,         only : aux_sld
  use def_solidz,         only : dttau_sld, ddism_sld
  use def_solidz,         only : cpu_ass_sol_sld
  use def_solidz,         only : SLD_STATIC_PROBLEM
  use def_solidz,         only : fcont_sld, saved_rhsid
  use mod_sld_fibers,     only : sld_fib_initialise
  use mod_sld_csys,       only : sld_csys_assign_material_axes
  use mod_exm_cou,        only : mod_exm_cou_initvar
  use mod_exm_cou,        only : mod_exm_cou_initexchange
  use mod_exm_sld_eccoupling, only : kfl_exmsld_3Dcou_ecc, has_exmsld_coupling
  use mod_exm_sld_eccoupling, only : exm_sld_ecc_manage_arrays

  implicit none

  integer(ip), intent(in) :: itask       !< when this subroutine is called
  integer(ip)             :: ivari,ielem

  select case(itask)

     ! IMPORTANT: nvarp is the maximum number of wopos. it is defined in def_kintyp.
     ! check it before adding new woposes.

  case(0)
     !
     !Flags
     call mod_exm_cou_initvar()

     !
     ! Postprocess Variable names and types
     !
     call arrays_register(  1_ip, (/'DISPL','VECTO','NPOIN','PRIMA'/),     displ, ENTITY_POSITION= 2_ip, TIME_POSITION= 3_ip )
     call arrays_register(  2_ip, (/'VELOC','VECTO','NPOIN','PRIMA'/), veloc_sld, ENTITY_POSITION= 2_ip, TIME_POSITION= 3_ip )
     call arrays_register(  3_ip, (/'ACCEL','VECTO','NPOIN','PRIMA'/), accel_sld, ENTITY_POSITION= 2_ip, TIME_POSITION= 3_ip )
     call arrays_register(  4_ip, (/'SIGMA','TENSO','NPOIN','SECON'/) )
     call arrays_register(  5_ip, (/'EPSIL','TENSO','NPOIN','SECON'/) )
     call arrays_register(  6_ip, (/'LNEPS','TENSO','NPOIN','SECON'/) )
     call arrays_register(  7_ip, (/'SEQVM','SCALA','NPOIN','SECON'/) )
     call arrays_register(  8_ip, (/'FIBER','VECTO','NPOIN','SECON'/) )
     call arrays_register(  9_ip, (/'BVESS','VECTO','NPOIN','SECON'/) )
     call arrays_register( 10_ip, (/'FIBEL','R3PVE','NPOIN','SECON'/) )
     call arrays_register( 11_ip, (/'DXFEM','VECTO','NPOIN','SECON'/) )
     call arrays_register( 12_ip, (/'VXFEM','VECTO','NPOIN','SECON'/) )
     call arrays_register( 13_ip, (/'AXFEM','VECTO','NPOIN','SECON'/) )
     call arrays_register( 14_ip, (/'LINEL','SCALA','NPOIN','SECON'/) )
     call arrays_register( 15_ip, (/'NSIGN','SCALA','NPOIN','SECON'/) )
     call arrays_register( 16_ip, (/'XYROT','SCALA','NPOIN','SECON'/) )
     call arrays_register( 17_ip, (/'ZXROT','SCALA','NPOIN','SECON'/) )
     call arrays_register( 18_ip, (/'DCOHE','SCALA','NELEM','SECON'/) )
     call arrays_register( 19_ip, (/'LNENR','SCALA','NPOIN','SECON'/) )
     call arrays_register( 20_ip, (/'ERROR','VECTO','NPOIN','SECON'/) )
     call arrays_register( 21_ip, (/'CRACK','SCALA','NPOIN','SECON'/) )
     call arrays_register( 22_ip, (/'GROUP','SCALA','NPOIN','SECON'/) )
     call arrays_register( 23_ip, (/'SIGEI','SCALA','NPOIN','SECON'/) )
     call arrays_register( 24_ip, (/'SDV1E','SCALA','NELEM','SECON'/) )
     call arrays_register( 25_ip, (/'SDV2E','SCALA','NELEM','SECON'/) )
     call arrays_register( 26_ip, (/'SDV3E','SCALA','NELEM','SECON'/) )
     call arrays_register( 27_ip, (/'SDV4E','SCALA','NELEM','SECON'/) )
     call arrays_register( 28_ip, (/'FIXRS','SCALA','NPOIN','SECON'/) )
     call arrays_register( 29_ip, (/'SDV1G','R3P  ','NPOIN','SECON'/) )
     call arrays_register( 30_ip, (/'SDV2G','R3P  ','NPOIN','SECON'/) )
     call arrays_register( 31_ip, (/'SDV3G','R3P  ','NPOIN','SECON'/) )
     call arrays_register( 32_ip, (/'SDV4G','R3P  ','NPOIN','SECON'/) )
     call arrays_register( 33_ip, (/'FIXNO','VECTO','NPOIN','SECON'/) )
     call arrays_register( 34_ip, (/'EXAFO','VECTO','NPOIN','SECON'/) )
     call arrays_register( 35_ip, (/'FORCF','VECTO','NPOIN','SECON'/) )
     call arrays_register( 36_ip, (/'SEGMA','R3PVE','NELEM','SECON'/) )
     call arrays_register( 37_ip, (/'FRXID','VECTO','NPOIN','SECON'/) )
     call arrays_register( 38_ip, (/'SIGNP','R3PVE','NPOIN','SECON'/) )
     call arrays_register( 39_ip, (/'SIGLO','R3PVE','NPOIN','SECON'/) )
     call arrays_register( 40_ip, (/'AXIS1','VECTO','NELEM','SECON'/) )
     call arrays_register( 41_ip, (/'AXIS2','VECTO','NELEM','SECON'/) )
     call arrays_register( 42_ip, (/'AXIS3','VECTO','NELEM','SECON'/) )
     call arrays_register( 43_ip, (/'ORIEN','SCALA','NELEM','SECON'/) )
     call arrays_register( 44_ip, (/'SREAC','VECTO','NPOIN','SECON'/) )
     call arrays_register( 45_ip, (/'NSIGM','VECTO','NPOIN','SECON'/) )
     call arrays_register( 46_ip, (/'PMATE','SCALA','NELEM','SECON'/) )
     call arrays_register( 47_ip, (/'TOUCH','SCALA','NPOIN','SECON'/) )
     call arrays_register( 48_ip, (/'NFIXN','SCALA','NPOIN','SECON'/) )
     call arrays_register( 49_ip, (/'TSIGM','VECTO','NPOIN','SECON'/) )
     call arrays_register( 50_ip, (/'DIISO','SCALA','NPOIN','SECON'/) )
     call arrays_register( 51_ip, (/'VMISO','SCALA','NPOIN','SECON'/) )
     call arrays_register( 52_ip, (/'BOSET','SCALA','NPOIN','SECON'/) )
     call arrays_register( 53_ip, (/'NOMIN','VECTO','NPOIN','SECON'/) )
     call arrays_register( 54_ip, (/'NNSIG','SCALA','NPOIN','SECON'/) )
     call arrays_register( 55_ip, (/'ELSET','SCALA','NELEM','SECON'/) )
     call arrays_register( 56_ip, (/'KDTRE','VECTO','NPOIN','SECON'/) )
     call arrays_register( 57_ip, (/'EPSIS','TENSO','NPOIN','SECON'/) )
     call arrays_register( 58_ip, (/'PK2  ','TENSO','NPOIN','SECON'/) )
     call arrays_register( 59_ip, (/'INTLI','SCALA','NPOIN','SECON'/) )
     call arrays_register( 60_ip, (/'SVEGM','R3PVE','NELEM','PRIMA'/), svegm_sld, ENTITY_POSITION= 1_ip, TIME_POSITION= 3_ip )
     call arrays_register( 61_ip, (/'CELEN','SCALA','NELEM','SECON'/) )
     call arrays_register( 62_ip, (/'FIXBO','SCALA','NPOIN','SECON'/) )
     call arrays_register( 63_ip, (/'BOCOD','SCALA','NPOIN','SECON'/) )
     call arrays_register( 64_ip, (/'ELNOR','VECTO','NELEM','SECON'/) )
     call arrays_register( 65_ip, (/'FEXTE','VECTO','NPOIN','SECON'/) )
     call arrays_register( 66_ip, (/'WETNO','SCALA','NPOIN','SECON'/) )
     call arrays_register( 67_ip, (/'SBVNA','VECTO','NPOIN','SECON'/) )
     call arrays_register( 68_ip, (/'ELSET','SCALA','NELEM','SECON'/) )
     call arrays_register( 69_ip, (/'PARTI','SCALA','NELEM','SECON'/) )
     call arrays_register( 70_ip, (/'ELEPS','TENSO','NPOIN','SECON'/) )
     call arrays_register( 71_ip, (/'ROTM1','VECTO','NPOIN','SECON'/) )
     call arrays_register( 72_ip, (/'ROTM2','VECTO','NPOIN','SECON'/) )
     call arrays_register( 73_ip, (/'ROTM3','VECTO','NPOIN','SECON'/) )
     call arrays_register( 74_ip, (/'FCONT','VECTO','NPOIN','SECON'/) )
     call arrays_register( 75_ip, (/'DAMAG','MULTI','NELEM','SECON'/) )
     call arrays_register( 90_ip, (/'MICNL','SCALA','NELEM','SECON'/) )
     call arrays_register( 91_ip, (/'MICCO','SCALA','NELEM','SECON'/) )
     call arrays_register( 92_ip, (/'MICCV','SCALA','NELEM','SECON'/) )
     call arrays_register( 93_ip, (/'EPRIN','SCALA','NPOIN','SECON'/) )
     !
     ! Register arrays (SLD-EXM-ECC Coupling): these arrays take 5 positions
     ! that is, from narray + 1 count 5 more and DO NOT put any array_register
     !
     if( has_exmsld_coupling() .or. kfl_exmsld_3Dcou_ecc ) then
         call exm_sld_ecc_manage_arrays(0_ip,narray=93_ip)
     endif
     !
     ! Sets variables
     !
     ! Boundary sets
     postp(1) % wobse (1)     = 'FORCE'  ! Force
     postp(1) % wobse (2)     = 'F_y'
     postp(1) % wobse (3)     = 'F_z'
     postp(1) % wobse (4)     = 'NORMA'  ! Normal Force
     postp(1) % wobse (5)     = 'DIBOX'  ! Averaged displacement in x-direction
     postp(1) % wobse (6)     = 'DIBOY'  ! Averaged displacement in y-direction
     postp(1) % wobse (7)     = 'DIBOZ'  ! Averaged displacement in z-direction
     postp(1) % wobse (8)     = 'DIBOU'  ! Averaged displacement (Magnitude)
     postp(1) % wobse (9)     = 'FRBOX'  ! Sum of reaction force in x-direction
     postp(1) % wobse (10)    = 'FRBOY'  ! Sum of reaction force in y-direction
     postp(1) % wobse (11)    = 'FRBOZ'  ! Sum of reaction force in z-direction
     postp(1) % wobse (12)    = 'FRBOU'  ! Sum of reaction force (Magnitude)
     postp(1) % wobse (13)    = 'FCONX'  ! Sum of contact force in x-direction
     postp(1) % wobse (14)    = 'FCONY'  ! Sum of contact force in y-direction
     postp(1) % wobse (15)    = 'FCONZ'  ! Sum of contact force in z-direction
     postp(1) % wobse (16)    = 'FCONT'  ! Sum of contact force (Magnitude)
     postp(1) % wobse (17)    = 'PRESS'  ! Pressure
     postp(1) % wobse (18)    = 'TRNOR'  ! Normal traction
     ! Node sets
     postp(1) % wonse (1)     = 'DISPX'  ! x-displacement
     postp(1) % wonse (2)     = 'DISPY'  ! y-displacement
     postp(1) % wonse (3)     = 'DISPZ'  ! z-displacement
     postp(1) % wonse (4)     = 'VELOX'  ! x-velocity
     postp(1) % wonse (5)     = 'VELOY'  ! y-velocity
     postp(1) % wonse (6)     = 'VELOZ'  ! z-velocity
     postp(1) % wonse (7)     = 'ACCEX'  ! x-acceleration
     postp(1) % wonse (8)     = 'ACCEY'  ! y-acceleration
     postp(1) % wonse (9)     = 'ACCEZ'  ! z-acceleration
     postp(1) % wonse (10)    = 'FRXIX'  ! x-reaction force
     postp(1) % wonse (11)    = 'FRXIY'  ! y-reaction force
     postp(1) % wonse (12)    = 'FRXIZ'  ! z-reaction force
     postp(1) % wonse (13)    = 'SIGXX'  ! Stress-XX
     postp(1) % wonse (14)    = 'SIGYY'  ! Stress-YY
     postp(1) % wonse (15)    = 'SIGZZ'  ! Stress-ZZ
     postp(1) % wonse (16)    = 'SIGYZ'  ! Stress-YZ
     postp(1) % wonse (17)    = 'SIGXZ'  ! Stress-XZ
     postp(1) % wonse (18)    = 'SIGXY'  ! Stress-XY
     postp(1) % wonse (19)    = 'EPSXX'  ! Strain-XX
     postp(1) % wonse (20)    = 'EPSYY'  ! Strain-YY
     postp(1) % wonse (21)    = 'EPSZZ'  ! Strain-ZZ
     postp(1) % wonse (22)    = 'EPSYZ'  ! Strain-YZ
     postp(1) % wonse (23)    = 'EPSXZ'  ! Strain-XZ
     postp(1) % wonse (24)    = 'EPSXY'  ! Strain-XY
     postp(1) % wonse (25)    = 'LEPXX'  ! Logarithmic Strain-XX
     postp(1) % wonse (26)    = 'LEPYY'  ! Logarithmic Strain-YY
     postp(1) % wonse (27)    = 'LEPZZ'  ! Logarithmic Strain-ZZ
     postp(1) % wonse (28)    = 'LEPYZ'  ! Logarithmic Strain-YZ
     postp(1) % wonse (29)    = 'LEPXZ'  ! Logarithmic Strain-XZ
     postp(1) % wonse (30)    = 'LEPXY'  ! Logarithmic Strain-XY
     postp(1) % wonse (31)    = 'COORX'  ! x-coordinate
     postp(1) % wonse (32)    = 'COORY'  ! y-coordinate
     postp(1) % wonse (33)    = 'COORZ'  ! z-coordinate
     postp(1) % wonse (34)    = 'SEQVM'  ! Von Mises Stress
     postp(1) % wonse (35)    = 'FEXTX'  ! External force x-component
     postp(1) % wonse (36)    = 'FEXTY'  ! External force y-component
     postp(1) % wonse (37)    = 'FEXTZ'  ! External force z-component
     postp(1) % wonse (38)    = 'FCONX'  ! Contact force x-component
     postp(1) % wonse (39)    = 'FCONY'  ! Contact force y-component
     postp(1) % wonse (40)    = 'FCONZ'  ! Contact force z-component
     postp(1) % wonse (41)    = 'EPRIN'  ! Principal stretches
     ! Element sets
     postp(1) % woese (1)     = 'VOLUM'  ! Volume  of the set
     postp(1) % woese (2)     = 'INV1E'  ! 1st Invariant of green-lagrange strain tensor
     postp(1) % woese (3)     = 'INV2E'  ! 2nd Invariant of green-lagrange strain tensor
     postp(1) % woese (4)     = 'INV3E'  ! 3rd Invariant of green-lagrange strain tensor
     postp(1) % woese (5)     = 'GRL11'  ! component of green-lagrange strain tensor
     postp(1) % woese (6)     = 'GRL12'  ! component of green-lagrange strain tensor
     postp(1) % woese (7)     = 'GRL13'  ! component of green-lagrange strain tensor
     postp(1) % woese (8)     = 'GRL21'  ! component of green-lagrange strain tensor
     postp(1) % woese (9)     = 'GRL22'  ! component of green-lagrange strain tensor
     postp(1) % woese (10)    = 'GRL23'  ! component of green-lagrange strain tensor
     postp(1) % woese (11)    = 'GRL31'  ! component of green-lagrange strain tensor
     postp(1) % woese (12)    = 'GRL32'  ! component of green-lagrange strain tensor
     postp(1) % woese (13)    = 'GRL33'  ! component of green-lagrange strain tensor
     postp(1) % woese (14)    = 'LEPRT'  ! Averaged log strain in polar cylindrical (eps_rr + eps_tt)
     postp(1) % woese (15)    = 'VELOC'  ! Mean velocity
     postp(1) % woese (16)    = 'TMASS'  ! Total mass
     postp(1) % woese (17)    = 'ALLKE'  ! All kinetic energy
     !
     ! Witness variables
     !
     postp(1) % wowit (1)     = 'DISPX'
     postp(1) % wowit (2)     = 'DISPY'
     postp(1) % wowit (3)     = 'DISPZ'
     postp(1) % wowit (4)     = 'VELOX'
     postp(1) % wowit (5)     = 'VELOY'
     postp(1) % wowit (6)     = 'VELOZ'
     postp(1) % wowit (7)     = 'ACCEX'
     postp(1) % wowit (8)     = 'ACCEY'
     postp(1) % wowit (9)     = 'ACCEZ'
     postp(1) % wowit (10)    = 'FRXIX'
     postp(1) % wowit (11)    = 'FRXIY'
     postp(1) % wowit (12)    = 'FRXIZ'
     postp(1) % wowit (13)    = 'SIGXX'
     postp(1) % wowit (14)    = 'SIGYY'
     postp(1) % wowit (15)    = 'SIGZZ'
     postp(1) % wowit (16)    = 'SIGYZ'
     postp(1) % wowit (17)    = 'SIGXZ'
     postp(1) % wowit (18)    = 'SIGXY'
     postp(1) % wowit (19)    = 'EPSXX'
     postp(1) % wowit (20)    = 'EPSYY'
     postp(1) % wowit (21)    = 'EPSZZ'
     postp(1) % wowit (22)    = 'EPSYZ'
     postp(1) % wowit (23)    = 'EPSXZ'
     postp(1) % wowit (24)    = 'EPSXY'
     postp(1) % wowit (25)    = 'LEPXX'
     postp(1) % wowit (26)    = 'LEPYY'
     postp(1) % wowit (27)    = 'LEPZZ'
     postp(1) % wowit (28)    = 'LEPYZ'
     postp(1) % wowit (29)    = 'LEPXZ'
     postp(1) % wowit (30)    = 'LEPXY'
     postp(1) % wowit (31)    = 'COORX'
     postp(1) % wowit (32)    = 'COORY'
     postp(1) % wowit (33)    = 'COORZ'
     postp(1) % wowit (34)    = 'SEQVM'
     postp(1) % wowit (35)    = 'SIGFI'
     postp(1) % wowit (36)    = 'MICRO'

     !
     ! Solvers
     !
     call soldef(-2_ip)                            ! Allocate memory for solvers
     solve(1) % kfl_solve = 1                      ! Output flag
     solve(1) % wprob     = 'DISPLACEMENT'
     solve(1) % ndofn     = ndime

     solve(1) % kfl_algso = SOL_SOLVER_RICHARDSON
     solve(1) % kfl_preco = SOL_LOCAL_MASS_MATRIX  ! explicit is the default

     solve(1) % kfl_force = 1                      ! Force solution continuity after calling the solver

     solve(2) % kfl_solve = 1                      ! Output flag ! OJO CAMBIAR
     solve(2) % wprob     = 'DISP-ROT'
     solve(2) % ndofn     = 6
     !solve(3) % kfl_solve = 1                     ! Output flag

     !
     ! Materials
     !
     nmate_sld =  nmate
     lmate_sld => lmate

     !
     ! Postprocess rotation origin
     !
     rorig_sld(1) = 0.0_rp
     rorig_sld(2) = 0.0_rp
     rorig_sld(3) = 0.0_rp

     !
     ! Force vector norms
     !
     fext2_sld = 0.0_rp
     fint2_sld = 0.0_rp
     fine2_sld = 0.0_rp
     fnatu_sld = 0.0_rp

     !
     ! Relaxation parameters
     !
     relpa_sld = 0.0_rp ! for solvers
     relco_sld = 1.0_rp ! for coupling

     !
     ! Others
     !
     cpu_ass_sol_sld = 0.0_rp
     dtcri_sld       = 0.0_rp

     !
     !
     ! Nullify pointers
     !
     ! Fixities and bc codes
     nullify(kfl_fixbo_sld)
     nullify(kfl_fixno_sld)
     nullify(kfl_fixrs_sld)
     nullify(kfl_funbo_sld)
     nullify(kfl_funno_sld)
     nullify(kfl_funtn_sld)
     nullify(kfl_funtb_sld)
     nullify(tbcod_sld)
     nullify(tncod_sld)
     nullify(tgcod_sld)
     nullify(bvess_sld)
     nullify(bvnat_sld)
     nullify(jacrot_du_dq_sld)
     nullify(jacrot_dq_du_sld)
     ! Motion variables and mass
     nullify(dunkn_sld)
     nullify(ddisp_sld)
     nullify(vmass_sld)
     nullify(veloc_sld)
     nullify(accel_sld)
     nullify(uruk4_sld)
     ! Force vectors
     nullify(finte_sld)
     nullify(fexte_sld)
     nullify(macce_sld)
     nullify(fintt_sld)
     nullify(fextt_sld)
     nullify(frxid_sld)
     ! Energies
     nullify(allie_sld)
     nullify(allwk_sld)
     nullify(allke_sld)
     nullify(etota_sld)
     ! Stresses and strains
     nullify(epsel_sld)
     nullify(lepse_sld)
     ! Contact
     nullify(fcont_sld)
     nullify(saved_rhsid)
     ! Other
     nullify(celen_sld)
     nullify(ddism_sld)
     nullify(dttau_sld)
     nullify(cause_sld)
     ! X-FEM
     nullify(dxfem_sld)
     nullify(vxfem_sld)
     nullify(axfem_sld)
     nullify(crapx_sld)
     nullify(cranx_sld)
     nullify(eleng_sld)
     nullify(lnenr_sld)
     nullify(leenr_sld)
     nullify(lecoh_sld)

     nullify(dslip_sld)
     nullify(dceff_sld)
     nullify(dcmax_sld)
     nullify(nopio_sld)
     nullify(treff_sld)
     nullify(nocoh_sld)
     nullify(cohnx_sld)

     !
     ! Composite stuff
     !
     nullify(axis1_sld)
     nullify(axis2_sld)
     nullify(axis3_sld)
     nullify(orien_sld)
     nullify(svegm_sld)
     nullify(svaux_sld)
     nullify(caunp_sld)
     nullify(caulo_sld)
     nullify(aux_sld)
     nullify(isoch_sld)
     nullify(iswav_sld)

  case(1)
     !
     ! Voigt dimensions
     !
     nvoig_sld = (ndime + 1) * ndime / 2
     !
     ! Voigt rule conversion
     !
     if ( ndime == 2_ip) then
        !
        ! 2D Voigt rule
        ! \sigma_{ij} \sigma_{a}
        !         ij          a
        !         11          1
        !         22          2
        !         12          3
        !
        nvgij_inv_sld(1,1) = 1_ip
        nvgij_inv_sld(2,2) = 2_ip
        nvgij_inv_sld(1,2) = 3_ip
        ! Symmetric
        nvgij_inv_sld(2,1) = nvgij_inv_sld(1,2)

        nvgij_sld(1,1) = 1
        nvgij_sld(1,2) = 1
        nvgij_sld(2,1) = 2
        nvgij_sld(2,2) = 2
        nvgij_sld(3,1) = 1
        nvgij_sld(3,2) = 2

     else if ( ndime == 3_ip) then
        !
        ! 3D Voigt rule
        ! \sigma_{ij} \sigma_{a}
        !         ij          a
        !         11          1
        !         22          2
        !         33          3
        !         23          4
        !         13          5
        !         12          6
        !
        nvgij_inv_sld(1,1) = 1_ip
        nvgij_inv_sld(2,2) = 2_ip
        nvgij_inv_sld(3,3) = 3_ip
        nvgij_inv_sld(2,3) = 4_ip
        nvgij_inv_sld(1,3) = 5_ip
        nvgij_inv_sld(1,2) = 6_ip
        ! Symmetric
        nvgij_inv_sld(2,1) = nvgij_inv_sld(1,2)
        nvgij_inv_sld(3,1) = nvgij_inv_sld(1,3)
        nvgij_inv_sld(3,2) = nvgij_inv_sld(2,3)

        nvgij_sld(1,1) = 1
        nvgij_sld(1,2) = 1
        nvgij_sld(2,1) = 2
        nvgij_sld(2,2) = 2
        nvgij_sld(3,1) = 3
        nvgij_sld(3,2) = 3
        nvgij_sld(4,1) = 2
        nvgij_sld(4,2) = 3
        nvgij_sld(5,1) = 1
        nvgij_sld(5,2) = 3
        nvgij_sld(6,1) = 1
        nvgij_sld(6,2) = 2

     end if
     !
     ! Dimensions
     !
     ndofn_sld = ndime
     if (kfl_xfeme_sld == 1) ndofn_sld = 2*ndime
     !
     ! Push
     !
     ivari = 7
     call posdef(25_ip,ivari)
     if( ivari /= 0 ) kfl_gdepo = 1
     !
     ! Number of temporal variables
     !
     if( kfl_tisch_sld == 1 ) then
        !
        ! Newmark-Beta Scheme/Central differences
        !
        ncomp_sld = 4_ip

     else if( kfl_tisch_sld == 2 ) then
        !
        ! Dissipative Tchamwa - Wielgosz scheme
        !
        ncomp_sld = 3_ip

     else if( kfl_tisch_sld == 3 ) then
        !
        ! Runge-Kutta Scheme
        !
        ncomp_sld = 3_ip
        nderi_sld = ndime*2
     else
        ncomp_sld = 3_ip
        nderi_sld = ndime*2
     end if
     nprev_sld = min(3_ip,ncomp_sld)  ! Last time step or global iteration
     !
     ! Time variables
     !
     dtinv_sld = 0.0_rp
     kfl_timei = 1_ip
     if ( kfl_timei_sld == SLD_STATIC_PROBLEM ) dtinv_sld = 1.0_rp
     kfl_stead_sld = 0_ip
     !
     ! Dimensions state variables
     !
     if( kfl_sdvar_sld == 1_ip ) then
        nsvar_sld = 0
        do ielem = 1,nelem
           if (     lawch_sld(lmate(ielem)) == 904_ip .or. &
                lawch_sld(lmate(ielem)) == 905_ip ) then
              ! Cohesive elements: 904 and 905
              nsvar_sld = max(nsvar_sld,5_ip)
           else if( lawst_sld(lmate(ielem)) == 152_ip ) then
              ! sm152
              nsvar_sld = max(nsvar_sld,7_ip)
           else if( lawst_sld(lmate(ielem)) == 153_ip ) then
              ! sm153
              nsvar_sld = max(nsvar_sld,3_ip)
           else if( lawst_sld(lmate(ielem)) == 154_ip ) then
              ! sm154
              nsvar_sld = max(nsvar_sld,12_ip)
           else if( lawst_sld(lmate(ielem)) == 0 .and. lawpl_sld(lmate(ielem))==1_ip ) then
              ! smXXX
              nsvar_sld = max(nsvar_sld,20_ip)
           else if( lawst_sld(lmate(ielem)) == 200_ip ) then
              ! sm200
              nsvar_sld = max(nsvar_sld,22_ip)
           end if
        end do
        call PAR_MAX(nsvar_sld)
     end if
     !
     ! Bio-Fiber model
     !
     call sld_fib_initialise()
     !
     ! Stresses field
     !
     if (kfl_restr_sld < 0) then
        restr_sld => xfiel(-kfl_restr_sld) % a(:,:,1)
     end if
     !
     ! Rigidity modulator field
     !
     if( kfl_moduf_sld(1) > 0 ) then
        fiemo_sld => xfiel(kfl_moduf_sld(1)) % a(:,:,1)
     end if
     !
     ! Polarization wave index
     !
     iwave_sld = 0
     !
     ! Interchange data with exmedi if there is coupling
     !
     if(kfl_exmsld_3Dcou_ecc) call  mod_exm_cou_initexchange()

  end select

end subroutine sld_inivar

