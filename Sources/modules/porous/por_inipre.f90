!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_inipre.f90
!> @date    13/08/2013
!> @author  Herbert Owen
!> @brief   Initialices pressure. Very similar to por_solite but only laplacian and gravity term.
!> @details Initialices pressure. Very similar to por_solite but only laplacian and gravity term.
!> @} 
!------------------------------------------------------------------------
subroutine por_inipre()
  use def_parame
  use def_master
  use def_domain
  use def_porous
!  use mod_gradie
  implicit none
  ! 
  ! set solver to presssure solver
  !
  call por_inisol()
  !
  ! Construct the system matrix and right-hand-side
  !
  call por_matipr()
  !
  ! initialize unkno to pressure for solver
  !
  call por_updunk(9_ip)
  !
  ! Solve the algebraic system
  !
  call solver(rhsid,unkno,amatr,pmatr)
  !
  ! press(:,1) = unkno
  !
  call por_updunk(3_ip)

end subroutine por_inipre
