subroutine Cgns24(order)
!-----------------------------------------------------------------------
!****f* Services/cgns24
! NAME
!    cgns24
! DESCRIPTION
!    This routine
! USED BY
!    outdom
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      mod_postpr
  use      mod_memchk
  implicit none
  integer(ip) :: order

  call runend('CGNS24: CANNOT USE THIS FORMAT: COMPILE THIS SERVICE')

end subroutine Cgns24
