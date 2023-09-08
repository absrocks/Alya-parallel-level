subroutine rad_outerr()
!------------------------------------------------------------------------
!****f* Radiat/rad_outerr
! NAME 
!    rad_outerr
! DESCRIPTION
!    This routine checks if there are errros and warnings
! USES
! USED BY
!    rad_turnon
!***
!------------------------------------------------------------------------
  use def_master
  use def_domain
  use def_radiat
  use mod_outfor, only : outfor
  implicit none
  integer(ip)   :: ierro=0,iwarn=0
  integer(ip)   :: ivara
  character(20) :: messa

!!$  !
!!$  ! Exact solution
!!$  !
!!$  if(kfl_sourc_rad==1.and.kfl_exacs_rad/=0) then
!!$     iwarn=iwarn+1
!!$     kfl_sourc_rad=0
!!$     call outfor(2_ip,momod(modul)%lun_outpu,&
!!$          'SOURCE TERM WAS AUTOMATICALLY SET TO ZERO TO SOLVE AN EXACT SOLUTION')     
!!$  end if

  !----------------------------------------------------------------------
  !
  ! Postprocess
  !
  !----------------------------------------------------------------------
  !
  ! Orthogonal projection
  !
  ivara = 10
  call posdef(25_ip,ivara)
  if( ivara /= 0 .and. kfl_ortho_rad == 0 ) then
     call posdef(26_ip,ivara)
     iwarn = iwarn + 1
     call outfor(2_ip,momod(modul)%lun_outpu,'CANNOT POSTPROCESS ORTHOGONAL PROJECTION')
  end if
  !----------------------------------------------------------------------
  !
  ! ERROR MESSAGE
  !
  !----------------------------------------------------------------------

  call errors(3_ip,ierro,iwarn,' ')

end subroutine rad_outerr
