subroutine ada_turnon
!-----------------------------------------------------------------------
!****f* adapti/ada_turnon
! NAME 
!    ada_turnon
! DESCRIPTION
!    This routine starts adapti
! USES
!
! USED BY
!
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_adapti
  implicit none

  !
  ! Open files
  !
  call ada_openfi(one)
  !
  ! Read the strategy parameters
  !
  call ada_reastr(one)
  !
  ! Compute element graph when needed
  !
  call ada_elmgra
  !
  ! Allocate work space
  !
  call ada_memall(one)
  !
  ! Read the strategy vectors
  !
  call ada_reastr(two)
  !
  ! Close input data file
  !
  call ada_openfi(ten)


end subroutine ada_turnon
