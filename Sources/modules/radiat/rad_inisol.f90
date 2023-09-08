subroutine rad_inisol()
  !-----------------------------------------------------------------------
  !   
  ! This routine loads the solver data for the radiation equation.
  ! In general, it may change from time step to time step or even
  ! from iteration to iteration.
  !
  !-----------------------------------------------------------------------
  use def_domain
  use def_radiat
  use def_master
  use def_solver
  implicit none

  solve_sol => solve                       ! Solver type
  solve_sol(1)%xdiag = 1.0_rp/dtinv

end subroutine rad_inisol
 
