subroutine shaga0(ngaus,shaga,ierro)

!-----------------------------------------------------------------------
!
! This routine evaluates shape functions associated to gauss points
!
!-----------------------------------------------------------------------
  use      def_kintyp
  implicit none
  integer(ip), intent(in)  :: ngaus
  integer(ip), intent(out) :: ierro
  real(rp)                 :: s
  real(rp),    intent(out) :: shaga(ngaus)

  ierro    = 0
  shaga(1) = 1.0_rp
  
end subroutine shaga0

