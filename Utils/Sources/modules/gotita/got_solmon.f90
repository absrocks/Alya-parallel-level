subroutine got_solmon()
  !-----------------------------------------------------------------------
  !****f* Gotita/got_solmon
  ! NAME 
  !    got_solmon
  ! DESCRIPTION
  !    This routine solves an iteration for the incompressible NS equations
  !    using a monolitic scheme.
  ! USES
  !    got_matrix
  !    solver
  ! USED BY
  !    got_solite
  !***
  !-----------------------------------------------------------------------
  use def_master
  use def_gotita
  implicit none
  !
  ! Construct the system matrix and right-hand-side
  !
  call got_matrix(1_ip)
  !
  ! Solve the algebraic system
  !
  call solver(rhsid,unkno,amatr,pmatr)
  
end subroutine got_solmon
