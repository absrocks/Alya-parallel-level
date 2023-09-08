subroutine tem_inisol()
  !-----------------------------------------------------------------------
  !   
  ! This routine loads the solver data for the temperature equation.
  ! In general, it may change from time step to time step or even
  ! from iteration to iteration.
  !
  !-----------------------------------------------------------------------
  use def_domain
  use def_temper
  use def_master
  use def_solver
  implicit none

  solve_sol => solve                       ! Solver type
  if( dtinv /= 0.0_rp ) solve_sol(1)%xdiag = 1.0_rp/dtinv

end subroutine tem_inisol
 
