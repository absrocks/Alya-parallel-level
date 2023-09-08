subroutine rad_solite()
  !-----------------------------------------------------------------------
  !****f* Radiat/rad_solite
  ! NAME 
  !    rad_solite
  ! DESCRIPTION
  !    This routine solves an iteration of the radiation heat transfer equations.
  ! USES
  !    rad_matrix
  !    Soldir
  !    Solite
  ! USED BY
  !    rad_doiter
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_radiat
  use mod_gradie

  implicit none

  integer(ip) :: kfl_advec_old,kfl_timei_old
  integer(ip) :: i,j,k
  !
  ! Update inner iteration counter and write headings in the solver
  ! file.
  !
  itinn(modul) = itinn(modul) + 1
  ittot_rad    = ittot_rad + 1
  !
  ! Update boundary conditions
  !
  call rad_updbcs(three)
!!$  !
!!$  ! Update density
!!$  !
!!$  call rad_updunk(7_ip)
  !
  ! Construct the system matrix and right-hand-side
  !
  call rad_matrix()
  !
  ! Solve the algebraic system
  !

  call solver(rhsid,unkno,amatr,pmatr)


end subroutine rad_solite
