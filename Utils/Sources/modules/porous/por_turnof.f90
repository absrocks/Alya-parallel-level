!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_turnof.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Closes the run for the porous equations.
!> @details Closes the run for the porous equations.
!> @} 
!------------------------------------------------------------------------
subroutine por_turnof()
  use def_parame
  use def_master
  use def_solver
  use def_porous
  implicit none
  !
  ! Output latex file
  !
!  call por_outlat(2_ip)
 
end subroutine por_turnof

