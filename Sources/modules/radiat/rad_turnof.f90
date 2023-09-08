subroutine rad_turnof()
  !-----------------------------------------------------------------------
  !****f* Radiat/rad_turnof
  ! NAME 
  !    rad_turnof
  ! DESCRIPTION
  !    This routine closes the run for the radiation heat transfer
  ! USES
  !    rad_outcpu
  !    rad_output
  ! USED BY
  !    Radiat
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_solver
  use def_radiat
  implicit none
  !
  ! Output latex file
  !
  !call rad_outlat(2_ip)
 
end subroutine rad_turnof

