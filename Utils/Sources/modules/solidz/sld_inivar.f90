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
  use def_master,         only : fiber,displ
  use def_kermod,         only : kfl_lface
  use def_domain,         only : coord, xfiel, lmate, ndime, nmate, npoin, nelem
  use def_solver,         only : SOL_LOCAL_MASS_MATRIX,SOL_SOLVER_RICHARDSON
  use mod_communications, only : PAR_MAX
  use def_solidz,         only : lmate_sld, nmate_sld, nvoig_sld, eleng_sld, celen_sld
  use def_solidz,         only : kfl_fixbo_sld, kfl_fixno_sld, kfl_fixrs_sld
  use def_solidz,         only : kfl_funbo_sld, kfl_funno_sld
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
  use def_solidz,         only : myofi_sld
  use def_solidz,         only : isoch_sld, iswav_sld
  use def_solidz,         only : nvgij_inv_sld, nvgij_sld, nsvar_sld
  use def_solidz,         only : kfl_indis_sld
  use def_solidz,         only : modfi_sld, rorig_sld
  use def_solidz,         only : vmass_sld
  use def_solidz,         only : kfl_xfeme_sld, vxfem_sld, dxfem_sld, axfem_sld
  use def_solidz,         only : cranx_sld, crapx_sld, lnenr_sld, leenr_sld
  use def_solidz,         only : lecoh_sld, treff_sld, cohnx_sld, nocoh_sld
  use def_solidz,         only : dcmax_sld, dceff_sld, dslip_sld
  use def_solidz,         only : nopio_sld
  use def_solidz,         only : iwave_sld, kacti_sld, kfl_gdepo, kfl_sdvar_sld
  use def_solidz,         only : lawch_sld, lawst_sld, lawpl_sld
  use def_solidz,         only : relco_sld, relpa_sld
  use def_solidz,         only : kfl_restr_sld, restr_sld
  use def_solidz,         only : kfl_timei_sld, kfl_stead_sld, kfl_tisch_sld
  use def_solidz,         only : ncomp_sld, nprev_sld, nderi_sld, dtinv_sld
  use def_solidz,         only : kfl_vofor_sld, vofor_sld
  use def_solidz,         only : kfl_moduf_sld, fiemo_sld
  use def_solidz,         only : kfl_fiber_sld, axis1_sld, axis2_sld, axis3_sld, orien_sld
  use def_solidz,         only : fibts_sld, fibtn_sld, modor_sld
  use def_solidz,         only : aux_sld
  use def_solidz,         only : dttau_sld, ddism_sld
  use def_solidz,         only : cpu_ass_sol_sld
  use def_solidz,         only : SLD_STATIC_PROBLEM
  use def_solidz,         only : fcont_sld
  use mod_sld_csys,       only : sld_csys_assign_material_axes

  implicit none

  integer(ip), intent(in) :: itask      !< when this subroutine is called
  integer(ip)             :: ipoin,idime,imate,ivari,ielem
  real(rp)                :: fimod,dummr(3),vauxi(3),vmodu

  select case(itask)

     ! IMPORTANT: nvarp is the maximum number of wopos. it is defined in def_kintyp.
     ! check it before adding new woposes.

  case(0)
     !
     ! Variable Names
     !
     postp(1) % wopos ( 1, 1) = 'DISPL'  ! Displacement
     postp(1) % wopos ( 1, 2) = 'VELOC'  ! Velocity
     postp(1) % wopos ( 1, 3) = 'ACCEL'  ! Acceleration
     postp(1) % wopos ( 1, 4) = 'SIGMA'  ! Cauchy Stresses
     postp(1) % wopos ( 1, 5) = 'EPSIL'  ! Strains (Green)
     postp(1) % wopos ( 1, 6) = 'LNEPS'  ! Logarithmic strains
     postp(1) % wopos ( 1, 7) = 'SEQVM'  ! Von Mises stress
     postp(1) % wopos ( 1, 8) = 'FIBER'
     postp(1) % wopos ( 1, 9) = 'BVESS'  ! Essential boundary conditions
     postp(1) % wopos ( 1,10) = 'FIBEL'
     postp(1) % wopos ( 1,11) = 'DXFEM'
     postp(1) % wopos ( 1,12) = 'VXFEM'
     postp(1) % wopos ( 1,13) = 'AXFEM'
     postp(1) % wopos ( 1,14) = 'LINEL'
     postp(1) % wopos ( 1,15) = 'NSIGN'
     postp(1) % wopos ( 1,16) = 'XYROT'
     postp(1) % wopos ( 1,17) = 'ZXROT'
     postp(1) % wopos ( 1,18) = 'DCOHE'  ! Cohesive elements: d (Scalar per element)
     postp(1) % wopos ( 1,19) = 'LNENR'
     postp(1) % wopos ( 1,20) = 'ERROR'
     postp(1) % wopos ( 1,21) = 'CRACK'
     postp(1) % wopos ( 1,22) = 'GROUP'
     postp(1) % wopos ( 1,23) = 'SIGEI'
     postp(1) % wopos ( 1,24) = 'SDV1E'  ! sm152: d1 (Scalar per element)
     postp(1) % wopos ( 1,25) = 'SDV2E'  ! sm152: dG (Scalar per element)
     postp(1) % wopos ( 1,26) = 'SDV3E'  ! sm152: dK (Scalar per element)
     postp(1) % wopos ( 1,27) = 'SDV4E'  ! sm152: d6 (Scalar per element)
     postp(1) % wopos ( 1,28) = 'FIXRS'  ! Fixity code for local axes
     postp(1) % wopos ( 1,29) = 'SDV1G'  ! sm152: d1 (Gauss point)
     postp(1) % wopos ( 1,30) = 'SDV2G'  ! sm152: dG (Gauss point)
     postp(1) % wopos ( 1,31) = 'SDV3G'  ! sm152: dK (Gauss point)
     postp(1) % wopos ( 1,32) = 'SDV4G'  ! sm152: d6 (Gauss point)
     postp(1) % wopos ( 1,33) = 'FIXNO'  ! Fixity codes for each DoF
     postp(1) % wopos ( 1,34) = 'EXAFO'
     postp(1) % wopos ( 1,35) = 'FORCF'
     postp(1) % wopos ( 1,36) = 'SEGMA'
     postp(1) % wopos ( 1,37) = 'FRXID'  ! Nodal Reaction Forces
     postp(1) % wopos ( 1,38) = 'SIGNP'
     postp(1) % wopos ( 1,39) = 'SIGLO'  ! Cauchy stress (sigma) according to the local coordinate system
     postp(1) % wopos ( 1,40) = 'AXIS1'  ! sm151 and 152: Fiber direction (1-axis)
     postp(1) % wopos ( 1,41) = 'AXIS2'  ! sm151 and 152: Transverse direction (2-axis)
     postp(1) % wopos ( 1,42) = 'AXIS3'  ! sm151 and 152: Normal deirection (3-axis)
     postp(1) % wopos ( 1,43) = 'ORIEN'  ! sm151 and 152: Ply orientation
     postp(1) % wopos ( 1,44) = 'SREAC'  ! Solver reaction
     postp(1) % wopos ( 1,45) = 'NSIGM'
     postp(1) % wopos ( 1,46) = 'PMATE'  ! Material number
     postp(1) % wopos ( 1,47) = 'TOUCH'
     postp(1) % wopos ( 1,48) = 'NFIXN'
     postp(1) % wopos ( 1,49) = 'TSIGM'
     postp(1) % wopos ( 1,50) = 'DIISO'
     postp(1) % wopos ( 1,51) = 'VMISO'
     postp(1) % wopos ( 1,52) = 'TRACE'
     postp(1) % wopos ( 1,53) = 'NOMIN'
     postp(1) % wopos ( 1,54) = 'NNSIG'
     postp(1) % wopos ( 1,55) = 'NNOMI'
     postp(1) % wopos ( 1,56) = 'KDTRE'
     postp(1) % wopos ( 1,57) = 'EPSIS'  ! Infinitesimal stresses
     postp(1) % wopos ( 1,58) = 'PK2'    ! Second Piola Kirchoff
     postp(1) % wopos ( 1,59) = 'INTLI'  ! PLEPP interior list points
     postp(1) % wopos ( 1,60) = 'SVEGM'  ! State Dependent Variables (SDVs)
     postp(1) % wopos ( 1,61) = 'CELEN'  ! Characteristic element length
     postp(1) % wopos ( 1,62) = 'FIXBO'  ! Fixity code boundaries
     postp(1) % wopos ( 1,63) = 'BOCOD'  ! Boundary codes
     postp(1) % wopos ( 1,64) = 'ELNOR'  ! Element normal for HEX08, PEN06 and QUA04 (Interface element)
     postp(1) % wopos ( 1,65) = 'FEXTE'  ! External nodal force
     postp(1) % wopos ( 1,66) = 'WETNO'  ! wet nodes for solidz
     postp(1) % wopos ( 1,67) = 'SBVNA'  ! Solver natural boundary condition
     postp(1) % wopos ( 1,68) = 'ELSET'  ! Element sets
     postp(1) % wopos ( 1,69) = 'PARTI'  ! Partitions
     postp(1) % wopos ( 1,70) = 'ELEPS'  ! Elastic strain
     postp(1) % wopos ( 1,71) = 'ROTM1'  ! Rotation matrix 1st vector
     postp(1) % wopos ( 1,72) = 'ROTM2'  ! Rotation matrix 2n vector
     postp(1) % wopos ( 1,73) = 'ROTM3'  ! Rotation matrix 3rd vector
     postp(1) % wopos ( 1,74) = 'FCONT'  ! PDN Contact force
     postp(1) % wopos ( 1,75) = 'DAMAG'  ! Damage varialbes for sm152
     postp(1) % wopos ( 1,90) = 'MICNL'  ! MicroPP that enter in non-linear range
     postp(1) % wopos ( 1,91) = 'MICCO'  ! MicroPP cost
     postp(1) % wopos ( 1,92) = 'MICCV'  ! MicroPP that has converged
     !
     ! Variable Types: SCALA,VECTO,R3P,R3PVE,TENSO
     !
     postp(1) % wopos ( 2, 1) = 'VECTO'
     postp(1) % wopos ( 2, 2) = 'VECTO'
     postp(1) % wopos ( 2, 3) = 'VECTO'
     postp(1) % wopos ( 2, 4) = 'TENSO'
     postp(1) % wopos ( 2, 5) = 'TENSO'
     postp(1) % wopos ( 2, 6) = 'TENSO'
     postp(1) % wopos ( 2, 7) = 'SCALA'
     postp(1) % wopos ( 2, 8) = 'VECTO'
     postp(1) % wopos ( 2, 9) = 'VECTO'
     postp(1) % wopos ( 2,10) = 'R3PVE'
     postp(1) % wopos ( 2,11) = 'VECTO'
     postp(1) % wopos ( 2,12) = 'VECTO'
     postp(1) % wopos ( 2,13) = 'VECTO'
     postp(1) % wopos ( 2,14) = 'SCALA'
     postp(1) % wopos ( 2,15) = 'SCALA'
     postp(1) % wopos ( 2,16) = 'SCALA'
     postp(1) % wopos ( 2,17) = 'SCALA'
     postp(1) % wopos ( 2,18) = 'SCALA'
     postp(1) % wopos ( 2,19) = 'SCALA'
     postp(1) % wopos ( 2,20) = 'VECTO'
     postp(1) % wopos ( 2,21) = 'SCALA'
     postp(1) % wopos ( 2,22) = 'SCALA'
     postp(1) % wopos ( 2,23) = 'SCALA'
     postp(1) % wopos ( 2,24) = 'SCALA'
     postp(1) % wopos ( 2,25) = 'SCALA'
     postp(1) % wopos ( 2,26) = 'SCALA'
     postp(1) % wopos ( 2,27) = 'SCALA'
     postp(1) % wopos ( 2,28) = 'SCALA'
     postp(1) % wopos ( 2,29) = 'R3P'
     postp(1) % wopos ( 2,30) = 'R3P'
     postp(1) % wopos ( 2,31) = 'R3P'
     postp(1) % wopos ( 2,32) = 'R3P'
     postp(1) % wopos ( 2,33) = 'VECTO'
     postp(1) % wopos ( 2,34) = 'VECTO'
     postp(1) % wopos ( 2,35) = 'VECTO'
     postp(1) % wopos ( 2,36) = 'R3PVE'
     postp(1) % wopos ( 2,37) = 'VECTO'
     postp(1) % wopos ( 2,38) = 'R3PVE'
     postp(1) % wopos ( 2,39) = 'R3PVE'
     postp(1) % wopos ( 2,40) = 'VECTO'
     postp(1) % wopos ( 2,41) = 'VECTO'
     postp(1) % wopos ( 2,42) = 'VECTO'
     postp(1) % wopos ( 2,43) = 'SCALA'
     postp(1) % wopos ( 2,44) = 'VECTO'
     postp(1) % wopos ( 2,45) = 'VECTO'
     postp(1) % wopos ( 2,46) = 'SCALA'
     postp(1) % wopos ( 2,47) = 'SCALA'
     postp(1) % wopos ( 2,48) = 'SCALA'
     postp(1) % wopos ( 2,49) = 'VECTO'
     postp(1) % wopos ( 2,50) = 'SCALA'
     postp(1) % wopos ( 2,51) = 'SCALA'
     postp(1) % wopos ( 2,52) = 'SCALA'
     postp(1) % wopos ( 2,53) = 'VECTO'
     postp(1) % wopos ( 2,54) = 'SCALA'
     postp(1) % wopos ( 2,55) = 'SCALA'
     postp(1) % wopos ( 2,56) = 'VECTO'
     postp(1) % wopos ( 2,57) = 'TENSO'
     postp(1) % wopos ( 2,58) = 'TENSO'
     postp(1) % wopos ( 2,59) = 'SCALA'
     postp(1) % wopos ( 2,60) = 'R3P'
     postp(1) % wopos ( 2,61) = 'SCALA'
     postp(1) % wopos ( 2,62) = 'SCALA'
     postp(1) % wopos ( 2,63) = 'SCALA'
     postp(1) % wopos ( 2,64) = 'VECTO'
     postp(1) % wopos ( 2,65) = 'VECTO'
     postp(1) % wopos ( 2,66) = 'SCALA'
     postp(1) % wopos ( 2,67) = 'VECTO'
     postp(1) % wopos ( 2,68) = 'SCALA'
     postp(1) % wopos ( 2,69) = 'SCALA'
     postp(1) % wopos ( 2,70) = 'TENSO'
     postp(1) % wopos ( 2,71) = 'VECTO'
     postp(1) % wopos ( 2,72) = 'VECTO'
     postp(1) % wopos ( 2,73) = 'VECTO'
     postp(1) % wopos ( 2,74) = 'VECTO'
     postp(1) % wopos ( 2,75) = 'MULTI'
     postp(1) % wopos ( 2,90) = 'SCALA'
     postp(1) % wopos ( 2,91) = 'SCALA'
     postp(1) % wopos ( 2,92) = 'SCALA'
     !
     ! Variable Locations: NPOIN (default), NELEM
     !
     postp(1) % wopos ( 3,18) = 'NELEM'
     postp(1) % wopos ( 3,36) = 'NELEM'
     postp(1) % wopos ( 3,24) = 'NELEM'
     postp(1) % wopos ( 3,25) = 'NELEM'
     postp(1) % wopos ( 3,26) = 'NELEM'
     postp(1) % wopos ( 3,27) = 'NELEM'
     postp(1) % wopos ( 3,40) = 'NELEM'
     postp(1) % wopos ( 3,41) = 'NELEM'
     postp(1) % wopos ( 3,42) = 'NELEM'
     postp(1) % wopos ( 3,43) = 'NELEM'
     postp(1) % wopos ( 3,46) = 'NELEM'
     postp(1) % wopos ( 3,61) = 'NELEM'
     postp(1) % wopos ( 3,64) = 'NELEM'
     postp(1) % wopos ( 3,68) = 'NELEM'
     postp(1) % wopos ( 3,69) = 'NELEM'
     postp(1) % wopos ( 3,75) = 'NELEM'
     postp(1) % wopos ( 3,90) = 'NELEM'
     postp(1) % wopos ( 3,91) = 'NELEM'
     postp(1) % wopos ( 3,92) = 'NELEM'
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
     postp(1) % woese (15)    = 'ALLKE'  ! Kinetic energy
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
     ! Nullify pointers
     !
     ! Fixities and bc codes
     nullify(kfl_fixbo_sld)
     nullify(kfl_fixno_sld)
     nullify(kfl_fixrs_sld)
     nullify(kfl_funbo_sld)
     nullify(kfl_funno_sld)
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
     ! Rigid body
     !nullify(rbdil_sld)
     !nullify(rbvel_sld)
     !nullify(rbacl_sld)
     !!nullify(rbfor_sld)
     ! Other
     nullify(celen_sld)
     nullify(ddism_sld)
     nullify(dttau_sld)
     nullify(cause_sld)
     nullify(kacti_sld)
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
     ! Material system or axes (fibers)
     !
     if( kfl_fiber_sld > 3 ) call sld_csys_assign_material_axes()
     !
     ! Volume force field
     !
     if( kfl_vofor_sld > 0 ) then
        vofor_sld => xfiel(kfl_vofor_sld) % a(:,:,1)
     end if
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
     ! XFEM ndofn solver
     !
     if (kfl_xfeme_sld == 1) then
        solve(1) % ndofn = 2*ndime
        !solve(1) % kfl_xfeme = 1
        kfl_lface = 1                      ! Requires global list of faces
        !
        ! Second solver: the enrichment field (X-FEM or alike)
        !
        !        solve(2) % ndofn     = ndime
        !        solve(2) % wprob     = 'DISP_XFEM'
     end if
     !
     ! Initial (REFERENCE) displacements field: COORD IS CHANGED when not stressed
     !
     if (kfl_indis_sld(1) < 0) then
        if( associated(xfiel)) then
           if (INOTMASTER) then
              if (kfl_indis_sld(2) == 0) then
                 do ipoin = 1,npoin
                    do idime= 1,ndime
                       coord(idime,ipoin) = coord(idime,ipoin) &
                            + xfiel(-kfl_indis_sld(1)) % a(idime,ipoin,1)
                    end do
!!                    write (1932,100) ipoin,10.0_rp * coord(1:3,ipoin)
                 end do
              end if
           end if
        else
           call runend('SLD_INIVAR: NO INITIAL DISPLACEMENTS FIELD ASSOCIATED.')
        end if
     end if
     !
     ! Fibers
     !
     if ( kfl_fiber_sld > 0 .and. kfl_fiber_sld < 4_ip ) then
        if( associated(xfiel) ) then
           do imate=1,nmate_sld
              if( modfi_sld(imate) < 0 .and. INOTMASTER ) then
                 fiber => xfiel(-modfi_sld(imate)) % a(:,:,1)
                 do ipoin = 1,npoin
                    fimod= 0.0_rp
                    do idime=1,ndime
                       fimod= fimod + fiber(idime,ipoin)*fiber(idime,ipoin)
                    end do
                    ! fibers are always normalized
                    fimod= sqrt(fimod)
                    if (fimod == 0.0_rp) fimod= 1.0_rp
                    do idime=1,ndime
                       fiber(idime,ipoin)= fiber(idime,ipoin) / fimod
                    end do
                 end do
              end if

              if( modor_sld(1,imate) < 0 .and. INOTMASTER ) then
                 fibts_sld => xfiel(-modor_sld(1,imate)) % a(:,:,1)
                 do ipoin=1,npoin
                    fimod= 0.0_rp
                    do idime=1,ndime
                       fimod= fimod + fibts_sld(idime,ipoin)*fibts_sld(idime,ipoin)
                    end do
                    ! fibers are always normalized
                    fimod= sqrt(fimod)
                    if (fimod == 0.0_rp) fimod= 1.0_rp
                    do idime=1,ndime
                       fibts_sld(idime,ipoin)= fibts_sld(idime,ipoin) / fimod
                    end do
                 end do
              end if

              if( modor_sld(2,imate) < 0 .and. INOTMASTER ) then
                 fibtn_sld => xfiel(-modor_sld(2,imate)) % a(:,:,1)
                 do ipoin = 1,npoin
                    fimod= 0.0_rp
                    do idime=1,ndime
                       fimod= fimod + fibtn_sld(idime,ipoin)*fibtn_sld(idime,ipoin)
                    end do
                    ! fibers are always normalized
                    fimod= sqrt(fimod)
                    if (fimod == 0.0_rp) fimod= 1.0_rp
                    do idime=1,ndime
                       fibtn_sld(idime,ipoin)= fibtn_sld(idime,ipoin) / fimod
                    end do
                 end do
              end if
           end do
        end if
     end if

     !
     ! Stresses field
     !
     if (kfl_restr_sld < 0) then
        if( associated(xfiel) ) then
           restr_sld => xfiel(-kfl_restr_sld) % a(:,:,1)
        else
           call runend('SLD_INIVAR: NO STRESS FIELD ASSOCIATED.')
        end if
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

     ! y esto??
     myofi_sld(1,:) = 0.0_rp   !xXBprer
     myofi_sld(2,:) = 0.007_rp * 0.0001_rp    !(cm) xXBpostr
     myofi_sld(3,:) = 1.0_rp   !Nxb
     myofi_sld(4,:) = 2.0579e-7_rp   !XBprer
     myofi_sld(5,:) = 8.7571e-7_rp   !Xbpostr
     myofi_sld(6,:) = (1.0_rp-1.8e-7_rp)  !N_NoXB
     myofi_sld(7,:) = 1.8e-7_rp  !P_NoXB
     myofi_sld(8,:) = 0.0196_rp   !TRPNCaL
     myofi_sld(9,:) = 0.1667_rp   !TRPNCaH
     myofi_sld(10,:) = 1.85_rp* 0.0001_rp  ! SL
     myofi_sld(11,:) = 0.0_rp  !integral of Force
     myofi_sld(12,:) = 0.0_rp  !P initialised within sld_myofil

  case (2)

     !
     ! Correct fibers when required
     !
     if (kfl_fiber_sld == 3 .and. INOTMASTER) then
        dummr(:) = 0.0_rp
        do ipoin = 1,npoin

           vauxi   = 0.0_rp
           vauxi(1)= 1.0_rp
           call vecpro(fiber(1,ipoin),vauxi,fibts_sld(1,ipoin),ndime)
           dummr(1:ndime) = fibts_sld(1:ndime,ipoin)
           vmodu= sqrt(dummr(1)*dummr(1) + dummr(2)*dummr(2) + dummr(3)*dummr(3))
           if (vmodu < 1.0e-10_rp) then
              vauxi   = 0.0_rp
              vauxi(2)= 1.0_rp
              call vecpro(fiber(1,ipoin),vauxi,fibts_sld(1,ipoin),ndime)
              dummr(1:ndime) = fibts_sld(1:ndime,ipoin)
              vmodu= sqrt(dummr(1)*dummr(1) + dummr(2)*dummr(2) + dummr(3)*dummr(3))
           end if

           fibts_sld(1:ndime,ipoin) = fibts_sld(1:ndime,ipoin)/vmodu

           call vecpro(fiber(1,ipoin),fibts_sld(1,ipoin),fibtn_sld(1,ipoin),ndime)
           dummr(1:ndime) = fibtn_sld(1:ndime,ipoin)
           vmodu= sqrt(dummr(1)*dummr(1) + dummr(2)*dummr(2) + dummr(3)*dummr(3))

           fibtn_sld(1:ndime,ipoin) = fibtn_sld(1:ndime,ipoin)/vmodu

        end do

     end if

     !
     ! Initial (REFERENCE) displacements field: DISPL IS SET when stressed
     !
     if (kfl_indis_sld(1) < 0) then
        if( associated(xfiel)) then
           if (INOTMASTER) then
              if (kfl_indis_sld(2) == 1) then
                 do ipoin = 1,npoin
                    do idime= 1,ndime
                       displ(idime,ipoin,nprev_sld) = xfiel(-kfl_indis_sld(1)) % a(idime,ipoin,1)
                    end do
                 end do
              end if
           end if
        else
           call runend('SLD_INIVAR: NO INITIAL DISPLACEMENTS FIELD ASSOCIATED.')
        end if
     end if

  end select

end subroutine sld_inivar

