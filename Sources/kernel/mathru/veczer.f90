subroutine veczer(n,v,rmod)

!-----------------------------------------------------------------------
!
! This routine initializes a vector
!
!-----------------------------------------------------------------------
  use      def_kintyp, only : ip,rp
  implicit none
  integer(ip), intent(in) :: n
  real(rp), intent(in)    :: rmod
  real(rp), intent(out)   :: v(n)
  integer(ip)             :: i

  do i=1,n
     v(i) = rmod
  end do
 
end subroutine veczer
      
