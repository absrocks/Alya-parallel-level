subroutine Doiter
!-----------------------------------------------------------------------
!****f* master/Doiter
! NAME
!    Doiter
! DESCRIPTION
!    This routine calls the different problems to be solved within
!    one iteration      
! USES
!    Nastin
!    Temper
!    Codire
!    Alefor
! USED BY
!    Alya
!***
!-----------------------------------------------------------------------
  use def_kintyp,        only : ip
  use def_master,        only : itinn
  use def_master,        only : ittim
  use def_master,        only : ITASK_DOITER
  use mod_ker_detection, only : ker_detection_doiter
  use def_kermod,        only : kfl_detection
  use mod_messages, only : livinf
  implicit none

  call livinf(5_ip,' ',0_ip)
  call livinf(6_ip,' ',0_ip)

  itinn(0) = ittim

  call moduls(ITASK_DOITER)

  !
  ! Detect non-converged modules
  !
  !if( kfl_detection /= 0 ) call ker_detection_doiter()

end subroutine Doiter
