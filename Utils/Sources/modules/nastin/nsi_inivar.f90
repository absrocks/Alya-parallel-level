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
  use def_kermod, only    :  gasco, kfl_adj_prob, velav_ker
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: ivar1,ivar2

  select case(itask)

  case(0_ip)
     !
     ! Postprocess Variable names and types
     !
     postp(1) % wopos ( 1, 1) = 'VELOC'
     postp(1) % wopos ( 1, 2) = 'PRESS'
     postp(1) % wopos ( 1, 3) = 'STREA'
     postp(1) % wopos ( 1, 4) = 'RESID'
     postp(1) % wopos ( 1, 5) = 'VESGS'
     postp(1) % wopos ( 1, 6) = 'BOUND'
     postp(1) % wopos ( 1, 7) = 'DENSI'
     postp(1) % wopos ( 1, 8) = 'PMV  '
     postp(1) % wopos ( 1, 9) = 'PPD  '
     postp(1) % wopos ( 1,10) = 'MACHN'
     postp(1) % wopos ( 1,11) = 'TANGE'
     postp(1) % wopos ( 1,12) = 'VORTI'
     postp(1) % wopos ( 1,13) = 'MODVO'
     postp(1) % wopos ( 1,14) = 'LAMB2'
     postp(1) % wopos ( 1,15) = 'GROUP'
     postp(1) % wopos ( 1,16) = 'LINEL'
     postp(1) % wopos ( 1,17) = 'VISCO'
     postp(1) % wopos ( 1,18) = 'FIXPR'
     postp(1) % wopos ( 1,19) = 'FIXNO'
     postp(1) % wopos ( 1,20) = 'TAU  '
     postp(1) % wopos ( 1,21) = 'AVVEL'
     postp(1) % wopos ( 1,22) = 'SCHUR'
     postp(1) % wopos ( 1,23) = 'VEPRO'
     postp(1) % wopos ( 1,24) = 'PRPRO'
     postp(1) % wopos ( 1,25) = 'YPLUS'
     postp(1) % wopos ( 1,26) = 'LIMIT'
     postp(1) % wopos ( 1,27) = 'GRPRO'
     postp(1) % wopos ( 1,28) = 'AVPRE'
     postp(1) % wopos ( 1,29) = 'VEERR'
     postp(1) % wopos ( 1,30) = 'PECLE'
     postp(1) % wopos ( 1,31) = 'HYDRO'
     postp(1) % wopos ( 1,32) = 'DUDT '
     postp(1) % wopos ( 1,33) = 'D2UDT'
     postp(1) % wopos ( 1,34) = 'DT1  '
     postp(1) % wopos ( 1,35) = 'DT2  '
     postp(1) % wopos ( 1,36) = 'PRERR'
     postp(1) % wopos ( 1,37) = 'LINVE'
     postp(1) % wopos ( 1,38) = 'DEFOR'
     postp(1) % wopos ( 1,39) = 'NODPR'
     postp(1) % wopos ( 1,40) = 'TTRAC'
     postp(1) % wopos ( 1,41) = 'VELO2'   ! For BDF restart
     postp(1) % wopos ( 1,42) = 'PRES2'   ! For BDF restart
     postp(1) % wopos ( 1,43) = 'VELO3'   ! For BDF restart
     postp(1) % wopos ( 1,44) = 'PRES3'   ! For BDF restart
     postp(1) % wopos ( 1,45) = 'MODVE'
     postp(1) % wopos ( 1,46) = 'AVTAN'   ! Average tau
     postp(1) % wopos ( 1,47) = 'LINVE'
     postp(1) % wopos ( 1,48) = 'VELOM'
     postp(1) % wopos ( 1,49) = 'FORCF'
     postp(1) % wopos ( 1,50) = 'INTFO'
     postp(1) % wopos ( 1,51) = 'USTAR'
     postp(1) % wopos ( 1,52) = 'TRACT'
     postp(1) % wopos ( 1,53) = 'FIXRS'
     postp(1) % wopos ( 1,54) = 'FVELO'   ! new velocity for filter output
     postp(1) % wopos ( 1,55) = 'FPRES'   ! new pressure for filter output
     postp(1) % wopos ( 1,56) = 'FTANG'   ! new wall shear stress for filter output
     postp(1) % wopos ( 1,57) = 'AVVE2'   ! Average Vi**2
     postp(1) % wopos ( 1,58) = 'AVVXY'   ! Average Vx*Vy, Vy*Vz, Vz*Vx
     postp(1) % wopos ( 1,59) = 'AVPR2'   ! Average Pr**2
     postp(1) % wopos ( 1,60) = 'MESHR'   ! Mesh Rotation
     postp(1) % wopos ( 1,61) = 'NODEF'   ! If reaction force are computed on nodes
     postp(1) % wopos ( 1,62) = 'MOMEN'   ! Momentum equaiton residual
     postp(1) % wopos ( 1,63) = 'FIXPP'   ! Pressure fixity
     postp(1) % wopos ( 1,64) = 'BVNAT'   ! Algebraic Neumann condition
     postp(1) % wopos ( 1,65) = 'REACT'   ! Algebraic reaction (b-Ax)
     postp(1) % wopos ( 1,66) = 'CHA01'   !
     postp(1) % wopos ( 1,67) = 'CHA02'   !
     postp(1) % wopos ( 1,68) = 'CHA03'   !
     postp(1) % wopos ( 1,69) = 'CHA04'   !
     postp(1) % wopos ( 1,70) = 'CHA05'   !
     postp(1) % wopos ( 1,71) = 'WETNO'   !
     postp(1) % wopos ( 1,72) = 'AVVRE'   ! For time-averaging restart
     postp(1) % wopos ( 1,73) = 'AVUPO'   ! For time-averaging restart
     postp(1) % wopos ( 1,74) = 'TURMU'   ! LES turbulent viscosity
     postp(1) % wopos ( 1,75) = 'AVMUT'   ! Average turbulent viscosity
     postp(1) % wopos ( 1,76) = 'AVSTX'   ! Average stress (mu_t * grad(u))
     postp(1) % wopos ( 1,77) = 'AVSTY'   ! Average stress (mu_t * grad(v))
     postp(1) % wopos ( 1,78) = 'AVSTZ'   ! Average stress (mu_t * grad(w))
     postp(1) % wopos ( 1,79) = 'BUBBL'   ! Pressure bubble
     postp(1) % wopos ( 1,80) = 'TAUPR'   ! Projection of tau
     postp(1) % wopos ( 1,81) = 'MFCFO'   ! For mass flow control restart
     postp(1) % wopos ( 1,82) = 'UBPRE'   ! For mass flow control restart
     postp(1) % wopos ( 1,83) = 'NORMA'   ! Normals
     postp(1) % wopos ( 1,84) = 'SENSM'   ! Mesh sensitivities
     postp(1) % wopos ( 1,85) = 'LAGRA'   ! Lagrange multiplier
     postp(1) % wopos ( 1,86) = 'ENVEL'
     postp(1) % wopos ( 1,87) = 'ENPRE'
     postp(1) % wopos ( 1,88) = 'ENVE2'   ! Ensemble Vi**2
     postp(1) % wopos ( 1,89) = 'ENVXY'   ! Ensemble Vx*Vy, Vy*Vz, Vz*Vx
     postp(1) % wopos ( 1,90) = 'ENPR2'   ! Ensemble Pr**2
     postp(1) % wopos ( 1,91) = 'ENTAN'   ! Ensemble tau
     postp(1) % wopos ( 1,92) = 'ENMUT'   ! Ensemble turbulent viscosity
     postp(1) % wopos ( 1,93) = 'ENSTX'   ! Ensemble stress (mu_t * grad(u))
     postp(1) % wopos ( 1,94) = 'ENSTY'   ! Ensemble stress (mu_t * grad(v))
     postp(1) % wopos ( 1,95) = 'ENSTZ'   ! Ensemble stress (mu_t * grad(w))
     postp(1) % wopos ( 1,96) = 'VELAV'   ! Velocity at boundaries
     postp(1) % wopos ( 1,97) = 'TITIM'   ! For two-layer model restart - LES part
     postp(1) % wopos ( 1,98) = 'TRACR'   ! For two-layer model restart - RANS part
     postp(1) % wopos ( 1,99) = 'TLUAV'   ! For two-layer model restart - Average velocity
     postp(1) % wopos ( 1,100) = 'AVVTA'  ! Average viscous part of tangential stress at the wall 4 restart - variable avta1_nsw_ker
     postp(1) % wopos ( 1,101) = 'VAFOR'  ! VARIATIONAL FORCE
     postp(1) % wopos ( 1,102) = 'AVVAF'  ! Average VARIATIONAL FORCE
     postp(1) % wopos ( 1,103) = 'NOTRA'  ! VARIATIONAL Tangential traction
     postp(1) % wopos ( 1,104) = 'AVNTR'  ! Average VARIATIONAL Tangential traction  - beware running average
     postp(1) % wopos ( 1,105) = 'AVGTR'  ! Average gradient based Tangential traction  - beware running average
     postp(1) % wopos ( 1,106) = 'FANSW'  ! Factor for no slip wall law - 4 restart 
     postp(1) % wopos ( 1,107) = 'PORFO'  ! Porous force 

     
     postp(1) % wopos ( 2, 1) = 'VECTO'
     postp(1) % wopos ( 2, 2) = 'SCALA'
     postp(1) % wopos ( 2, 3) = 'SCALA'
     postp(1) % wopos ( 2, 4) = 'SCALA'
     postp(1) % wopos ( 2, 5) = 'R3PVE'
     postp(1) % wopos ( 2, 6) = 'VECTO'
     postp(1) % wopos ( 2, 7) = 'SCALA'
     postp(1) % wopos ( 2, 8) = 'SCALA'
     postp(1) % wopos ( 2, 9) = 'SCALA'
     postp(1) % wopos ( 2,10) = 'SCALA'
     postp(1) % wopos ( 2,11) = 'VECTO'
     if(ndime==2) then
        postp(1) % wopos ( 2,12) = 'SCALA'
     else
        postp(1) % wopos ( 2,12) = 'VECTO'
     end if
     postp(1) % wopos ( 2,13) = 'SCALA'
     postp(1) % wopos ( 2,14) = 'SCALA'
     postp(1) % wopos ( 2,15) = 'SCALA'
     postp(1) % wopos ( 2,16) = 'SCALA'
     postp(1) % wopos ( 2,17) = 'SCALA'
     postp(1) % wopos ( 2,18) = 'SCALA'
     postp(1) % wopos ( 2,19) = 'VECTO'
     postp(1) % wopos ( 2,20) = 'SCALA'
     postp(1) % wopos ( 2,21) = 'VECTO'
     postp(1) % wopos ( 2,22) = 'SCALA'
     postp(1) % wopos ( 2,23) = 'VECTO'
     postp(1) % wopos ( 2,24) = 'SCALA'
     postp(1) % wopos ( 2,25) = 'SCALA'
     postp(1) % wopos ( 2,26) = 'SCALA'
     postp(1) % wopos ( 2,27) = 'VECTO'
     postp(1) % wopos ( 2,28) = 'SCALA'
     postp(1) % wopos ( 2,29) = 'SCALA'
     postp(1) % wopos ( 2,30) = 'SCALA'
     postp(1) % wopos ( 2,31) = 'SCALA'
     postp(1) % wopos ( 2,32) = 'VECTO'
     postp(1) % wopos ( 2,33) = 'VECTO'
     postp(1) % wopos ( 2,34) = 'SCALA'
     postp(1) % wopos ( 2,35) = 'SCALA'
     postp(1) % wopos ( 2,36) = 'SCALA'
     postp(1) % wopos ( 2,37) = 'SCALA'
     postp(1) % wopos ( 2,38) = 'VECTO'
     postp(1) % wopos ( 2,39) = 'SCALA'
     postp(1) % wopos ( 2,40) = 'VECTO'
     postp(1) % wopos ( 2,41) = 'VECTO'
     postp(1) % wopos ( 2,42) = 'SCALA'
     postp(1) % wopos ( 2,43) = 'VECTO'
     postp(1) % wopos ( 2,44) = 'SCALA'
     postp(1) % wopos ( 2,45) = 'SCALA'
     postp(1) % wopos ( 2,46) = 'VECTO'
     postp(1) % wopos ( 2,47) = 'SCALA'
     postp(1) % wopos ( 2,48) = 'VECTO'
     postp(1) % wopos ( 2,49) = 'VECTO'
     postp(1) % wopos ( 2,50) = 'VECTO'
     postp(1) % wopos ( 2,51) = 'SCALA'
     postp(1) % wopos ( 2,52) = 'VECTO'
     postp(1) % wopos ( 2,53) = 'SCALA'
     postp(1) % wopos ( 2,54) = 'VECTO'
     postp(1) % wopos ( 2,55) = 'SCALA'
     postp(1) % wopos ( 2,56) = 'VECTO'
     postp(1) % wopos ( 2,57) = 'VECTO'
     postp(1) % wopos ( 2,58) = 'VECTO'
     postp(1) % wopos ( 2,59) = 'SCALA'
     postp(1) % wopos ( 2,60) = 'VECTO'
     postp(1) % wopos ( 2,61) = 'SCALA'
     postp(1) % wopos ( 2,62) = 'VECTO'
     postp(1) % wopos ( 2,63) = 'SCALA'
     postp(1) % wopos ( 2,64) = 'SCALA'
     postp(1) % wopos ( 2,65) = 'VECTO'
     postp(1) % wopos ( 2,66) = 'VECTO'
     postp(1) % wopos ( 2,67) = 'VECTO'
     postp(1) % wopos ( 2,68) = 'VECTO'
     postp(1) % wopos ( 2,69) = 'VECTO'
     postp(1) % wopos ( 2,70) = 'VECTO'
     postp(1) % wopos ( 2,71) = 'SCALA'
     postp(1) % wopos ( 2,72) = 'VECTO'
     postp(1) % wopos ( 2,73) = 'VECTO'
     postp(1) % wopos ( 2,74) = 'SCALA'
     postp(1) % wopos ( 2,75) = 'SCALA'
     postp(1) % wopos ( 2,76) = 'VECTO'
     postp(1) % wopos ( 2,77) = 'VECTO'
     postp(1) % wopos ( 2,78) = 'VECTO'
     postp(1) % wopos ( 2,79) = 'SCALA'
     postp(1) % wopos ( 2,80) = 'SCALA'
     postp(1) % wopos ( 2,81) = 'SCALA'
     postp(1) % wopos ( 2,82) = 'SCALA'
     postp(1) % wopos ( 2,83) = 'SCALA'
     postp(1) % wopos ( 2,84) = 'VECTO'
     postp(1) % wopos ( 2,85) = 'VECTO'
     postp(1) % wopos ( 2,86) = 'VECTO'
     postp(1) % wopos ( 2,87) = 'SCALA'
     postp(1) % wopos ( 2,88) = 'VECTO'   ! Ensemble Vi**2
     postp(1) % wopos ( 2,89) = 'VECTO'   ! Ensemble Vx*Vy, Vy*Vz, Vz*Vx
     postp(1) % wopos ( 2,90) = 'SCALA'   ! Ensemble Pr**2
     postp(1) % wopos ( 2,91) = 'VECTO'   ! Ensemble tau
     postp(1) % wopos ( 2,92) = 'SCALA'   ! Ensemble turbulent viscosity
     postp(1) % wopos ( 2,93) = 'VECTO'   ! Ensemble stress (mu_t * grad(u))
     postp(1) % wopos ( 2,94) = 'VECTO'   ! Ensemble stress (mu_t * grad(v))
     postp(1) % wopos ( 2,95) = 'VECTO'   ! Ensemble stress (mu_t * grad(w))
     postp(1) % wopos ( 2,96) = 'VECTO'
     postp(1) % wopos ( 2,97) = 'SCALA'
     postp(1) % wopos ( 2,98) = 'VECTO'
     postp(1) % wopos ( 2,99) = 'VECTO'
     postp(1) % wopos ( 2,100) = 'VECTO'
     postp(1) % wopos ( 2,101) = 'VECTO'
     postp(1) % wopos ( 2,102) = 'VECTO'
     postp(1) % wopos ( 2,103) = 'VECTO'
     postp(1) % wopos ( 2,104) = 'VECTO'
     postp(1) % wopos ( 2,105) = 'VECTO'
     postp(1) % wopos ( 2,106) = 'SCALA'
     postp(1) % wopos ( 2,107) = 'VECTO'

     postp(1) % wopos ( 3,79) = 'NELEM'
     postp(1) % wopos ( 3,72) = 'WHATE'
     postp(1) % wopos ( 3,96) = 'NBOUN'
     postp(1) % wopos ( 3,100) = 'WHATE'

     !
     ! Set variables
     !
     postp(1) % woese (1)     = 'VELOC'
     postp(1) % woese (2)     = 'VORTI'
     postp(1) % woese (3)     = 'KINET'
     postp(1) % woese (4)     = 'DIVU2'

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
     call soldef(-7_ip)                   ! Allocate memory
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
     nullify(kfl_wlawf_nsi)
     nullify(bvess_nsi)
     nullify(bvnat_nsi)
     nullify(skcos_nsi)
     nullify(tncod_nsi)
     nullify(tgcod_nsi)
     nullify(tbcod_nsi)
     !nullify(itsta_nsi)
     nullify(veold_nsi)
     !nullify(resis_nsi)
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
     nullify(Auu_nsi)
     nullify(Aup_nsi)
     nullify(Apu_nsi)
     nullify(App_nsi)
     nullify(amatr_nsi)
     nullify(lapla_nsi)
     nullify(cmama_nsi)
     nullify(deltp_nsi)
     nullify(dt_rho_nsi)
     nullify(mu_rho_nsi)
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
     nullify(velav_ker)
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

     nullify(resdiff_nsi)
     nullify(dcost_dx_nsi)
     nullify(costdiff_nsi)
     
     nullify(lbpse)
     !
     ! Others
     !
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
        else if(kfl_tisch_nsi==3 .or. kfl_tisch_nsi==4) then
           !
           ! Explicit time step
           !
           if(kfl_tisch_nsi==4 .and. kfl_fscon_nsi == 0) then
              ncomp_nsi=4
           else
              ncomp_nsi=3
           end if
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
     ! Gradient of viscosity needed
     !
     if( kfl_cotur_nsi == 1 ) then
        kfl_grvis_nsi = 1
     else
        kfl_grvis_nsi = 0
     end if
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
     if( kfl_visco_nsi == 1 .and. kfl_grvis_nsi == 1      ) kfl_rmom2_nsi = 1
     if( kfl_advec_nsi == 1 .and. kfl_linea_nsi == 2      ) kfl_rmom2_nsi = 1
     if( kfl_advec_nsi == 1 .and. kfl_adj_prob == 1      )  kfl_rmom2_nsi = 1
     if( fvnoa_nsi     >  0.0_rp )                          kfl_p1ve2_nsi = 1

     if( IMASTER ) then
        ndbgs_nsi = 0
     else
        ndbgs_nsi = ndime * npoin
     end if

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

     nevat_nsi     = (ndime+1)*mnode
     nschu_nsi     = 0                               ! # Schur complement solves
     nmome_nsi     = 0                               ! # Momentum solves

     kfl_perip_nsi = 0                               ! Periodicity: if pressure prescribed on periodic nodes
     pcoef_nsi     = 1.0_rp-gasco/sphea_nsi
     gamth_nsi     = sphea_nsi/(sphea_nsi-gasco)     ! Low-Mach: gamma=Cp/(Cp-R)

     prthe(1)      = prthe_nsi                       ! Low-Mach: p0
     prthe(2)      = prthe_nsi                       ! Low-Mach: p0^{n}
     prthe(3)      = prthe_nsi                       ! Low-Mach: p0^{n-1}
     prthe(4)      = prthe_nsi                       ! Low-Mach: p0^{0}

     imass         = tmass_nsi                       ! initial density
     kfl_prthe     = kfl_prthe_nsi                   ! type of thpressure calculation
     tflux         = 0.0_rp                          ! Low-Mach: heat flux
     dpthe         = 0.0_rp                          ! Low-Mach: dp0/dt
     xmass_nsi     = 0.0_rp                          ! Low-Mach: mass
     cputi_nsi     = 0.0_rp
     cputi_assembly_nsi = 0.0_rp
     cpu_ass_sol_nsi = 0.0_rp
     iteqn_nsi(1)  = 0
     iteqn_nsi(2)  = 0
     ittot_nsi     = 0
     resin_nsi(1)  = 1.0_rp                          ! Algebraic inner residual
     resin_nsi(2)  = 1.0_rp
     resou_nsi(1)  = 1.0_rp                          ! Algebraic outer residual
     resou_nsi(2)  = 1.0_rp
     resss_nsi(1)  = 1.0_rp                          ! Algebraic Steady state residual
     resss_nsi(2)  = 1.0_rp
     reinf_nsi(1)  = 1.0_rp                          ! Linf residual
     reinf_nsi(2)  = 1.0_rp
     relpa_nsi(1)  = 0.0_rp                          ! Relaxation
     relpa_nsi(2)  = 0.0_rp
     corio_nsi     = 0.0_rp                          ! Coriolis force module
     kfl_tiaor_nsi = kfl_tiacc_nsi                   ! Time accuracy: save original value
     resgs_nsi(1) =  0.0_rp
     resgs_nsi(2) =  0.0_rp
     rmsgs_nsi    =  0.0_rp
     dtsgs_nsi    =  0.0_rp
     kfl_wlare_nsi=  0
     if(kfl_sgsco_nsi==0) misgs_nsi=1                ! Subgrid scale number of iterations
     actav_nsi     = 0.0_rp                          ! Accumulated time for averaging
     if( kfl_sgsac_nsi == -1 ) then                  ! SGS time accuracy default
        if(kfl_tisch_nsi/=2) then
           kfl_sgsac_nsi = kfl_tiacc_nsi
        else
           kfl_sgsac_nsi = 1_ip                      ! For BDF use order 1 for the subscale by default
        end if
     end if
     dtmax_nsi = 0.0_rp ! otherwise full ckeck stops because it is unitialized when doing max in nsi_tistep
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

  case(3_ip)
     !
     ! Before starting a time step
     !
     relax_nsi = 1.0_rp
     relap_nsi = 1.0_rp

  end select

end subroutine nsi_inivar
