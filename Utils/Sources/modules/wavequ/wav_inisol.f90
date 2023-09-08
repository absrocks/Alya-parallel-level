subroutine wav_inisol()
  !-----------------------------------------------------------------------
  !   
  ! This routine loads the solver data for the temperature equation.
  ! In general, it may change from time step to time step or even
  ! from iteration to iteration.
  !
  !-----------------------------------------------------------------------
  use def_master
  use def_solver
  implicit none

  solve_sol => solve

end subroutine wav_inisol
