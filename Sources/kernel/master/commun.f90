subroutine commun(messa)
  !------------------------------------------------------------------------
  !****f* master/commun
  ! NAME 
  !    commun
  ! DESCRIPTION
  !    This subroutine communicates with Alya
  !
  ! OUTPUT
  ! USES
  ! USED BY
  !    Turnon
  !    Turnof
  !***
  !------------------------------------------------------------------------
  use def_master
  implicit none
  character(100), intent(out) :: messa

  messa= ' '

  !rewind(lun_commu)
  !read(lun_commu,'(a)',end=10,err=10) messa

end subroutine commun
