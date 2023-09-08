subroutine got_restar(itask)
  !------------------------------------------------------------------------
  !****f* Gotita/got_restar
  ! NAME 
  !    got_restar
  ! DESCRIPTION
  !    This routine reads the initial values from the restart file
  ! USES
  ! USED BY
  !    got_turnon
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_gotita
  use mod_iofile
  implicit none
  integer(ip), intent(in) :: itask

end subroutine got_restar
