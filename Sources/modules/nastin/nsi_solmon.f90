subroutine nsi_solmon()
  !-----------------------------------------------------------------------
  !****f* Nastin/nsi_solmon
  ! NAME 
  !    nsi_solmon
  ! DESCRIPTION
  !    This routine solves an iteration for the incompressible NS equations
  !    using a monolitic scheme.
  ! USES
  !    nsi_ifconf
  !    nsi_matrix
  !    Soldir
  !    Solite
  !    nsi_rotunk
  ! USED BY
  !    nsi_solite
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_nastin
  implicit none
  !
  ! Construct the system matrix and right-hand-side
  !
  call nsi_inisol(2_ip)                                 ! Initialize solver
  call nsi_matrix()
  !call nsi_rescon()                                     ! Algebraic residual
  !
  ! Transform from global to local systems 
  !
  call nsi_rotunk(one,unkno)
  !
  ! Solve the algebraic system
  !
  call solver(rhsid,unkno,amatr,pmatr)
  !
  ! Transform from local to global systems 
  !
  call nsi_rotunk(two,unkno)

end subroutine nsi_solmon
