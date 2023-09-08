!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_inisol.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Loads the solver data for  the porous equation
!> @details Loads the solver data for  the porous equation
!!            In general, it may change from time step to time step or even
!!            from iteration to iteration.
!> @} 
!------------------------------------------------------------------------
subroutine por_inisol()
  use def_domain
  use def_porous
  use def_master
  use def_solver
  implicit none

  solve_sol => solve(kprsa_por:)                              ! Solver type
  if( dtinv /= 0.0_rp ) solve_sol(1)%xdiag = 1.0_rp / dtinv   ! Not used with LOC_DIAGONAL precond.

  call inisol()

end subroutine por_inisol
