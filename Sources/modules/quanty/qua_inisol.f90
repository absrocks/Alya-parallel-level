subroutine qua_inisol()
  !-----------------------------------------------------------------------
  !   
  ! This routine loads the solver data for the Shrodinger equation.
  !
  !-----------------------------------------------------------------------
  use def_domain
  use def_quanty
  use def_solver
  implicit none


  
   eigen_sol => eigen_qua                       ! Eigen Solver type
   solve_sol => solve_qua                       ! Eigen Solver type

end subroutine qua_inisol
 
