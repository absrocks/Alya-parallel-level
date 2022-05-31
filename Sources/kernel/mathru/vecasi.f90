subroutine vecasi(n,v1,v2)

!-----------------------------------------------------------------------
!
! Vector assign:    v2(i) = v1(i)   i=1..n
!
!-----------------------------------------------------------------------
  use      def_kintyp
  implicit none
  integer(ip), intent(in)  :: n
  real(rp),    intent(in)  :: v1(n)
  real(rp),    intent(out) :: v2(n)
  integer(ip)              :: i

  do i=1,n
     v2(i)=v1(i)
  end do
  
end subroutine vecasi


