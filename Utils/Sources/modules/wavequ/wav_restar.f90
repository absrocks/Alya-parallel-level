subroutine wav_restar(itask)
!------------------------------------------------------------------------
!****f* Wavequ/wav_restar
! NAME 
!    wav_restar
! DESCRIPTION
!    This routine writes in and reads from the restart file
! USES
! USED BY
!    wav_turnon
!***
!------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_wavequ
  use mod_postpr
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: icomp,iwopo,kfl_gores
  !
  ! Check if restrt file should be read or written
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
  ! Wave amplitude
  !
  !----------------------------------------------------------------------

  iwopo = 1
  if( INOTMASTER ) gesca => wavam(:,icomp)
  call postpr(gesca,postp(1)%wopos(1:2,iwopo),ittim,cutim)

  !----------------------------------------------------------------------
  !
  ! Finish
  !
  !----------------------------------------------------------------------

  call respre(3_ip,kfl_gores)

end subroutine wav_restar
