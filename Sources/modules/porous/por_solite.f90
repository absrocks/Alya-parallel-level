!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_solite.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Solves an iteration of the porous equations - pressure.
!> @details Solves an iteration of the porous equations - pressure.
!> @} 
!------------------------------------------------------------------------
subroutine por_solite()
  use def_parame
  use def_master
  use def_domain
  use def_porous
  implicit none
  !
  ! set solver to presssure solver
  !
  call por_inisol()
  !
  ! Construct the system matrix and right-hand-side
  !
  call por_matrpr()
  !
  ! initialize unkno to pressure for solver
  !
  call por_updunk(9_ip)
  !
  ! Solve the algebraic system
  !
  call solver(rhsid,unkno,amatr,pmatr)

end subroutine por_solite
