subroutine cgnser
!-----------------------------------------------------------------------
!****f* Services/cgnser
! NAME
!    cgnser
! DESCRIPTION
!    This routine manages CGNS errors
! USED BY
!    cgns24
!***
!-----------------------------------------------------------------------
  use      def_parame
  implicit none
  include  'cgnslib_f.h'
  character(150) :: message

  call cg_get_error_f(message) 
  call runend(trim(message))

end subroutine cgnser
