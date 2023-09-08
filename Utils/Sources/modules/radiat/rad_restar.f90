subroutine rad_restar(itask)
  !------------------------------------------------------------------------
  !****f* Radiat/rad_restar
  ! NAME 
  !    rad_restar
  ! DESCRIPTION
  !    This routine:
  !    ITASK = 1 ... Reads the initial values from the restart file
  !            2 ... Writes restart file
  ! USES
  ! USED BY
  !    rad_turnon
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_radiat
  use mod_postpr
  use mod_memchk
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: icomp,iwopo,kfl_gores
  !
  ! Check if restart file should be read or written
  !
  call respre(itask,kfl_gores)
  if( kfl_gores == 0 ) return

  if( itask == READ_RESTART_FILE ) then
     icomp = 3
  else
     icomp = 1
  end if

  !----------------------------------------------------------------------
  !
  ! Average radiation intensity
  !
  !----------------------------------------------------------------------

  iwopo = 1
  if( INOTMASTER ) gesca => radav_rad(:,icomp)
  call postpr(gesca,postp(1)%wopos(1:2,iwopo),ittim,cutim)

  !----------------------------------------------------------------------
  !
  ! Finish
  !
  !----------------------------------------------------------------------

  call respre(3_ip,kfl_gores)

end subroutine rad_restar
 
