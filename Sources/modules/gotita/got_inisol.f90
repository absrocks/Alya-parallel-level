subroutine got_inisol(itask)
  !-----------------------------------------------------------------------
  !****f* Gotita/got_inisol
  ! NAME 
  !    got_inisol
  ! DESCRIPTION
  !    This routine loads the solver data for Gotita equations.
  !    In general, it may change from time step to time step or even
  !    from iteration to iteration.
  ! USED BY
  !    got_begite
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_gotita
  use def_solver
  implicit none
  integer(ip), intent(in) :: itask
  !
  ! Which equation(s) is(are) to be solved
  !
  if(itask==1.or.itask==2) then
     ivari_got  = 1                                           ! Momentum (+Continuity) equation  
  else if(itask==3) then
     ivari_got  = 2                                           ! Continuity
  end if

  solve_sol => solve(ivari_got:)                      ! Solver type

end subroutine got_inisol
