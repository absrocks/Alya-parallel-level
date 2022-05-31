subroutine shaga1(    s,ngaus,shaga,ierro)

!-----------------------------------------------------------------------
!
! This routine evaluates shape functions associated to gauss points
! for 1D : NGAUS = 1, 2, 3
!
!-----------------------------------------------------------------------
  use      def_kintyp
  implicit none
  integer(ip), intent(in)  :: ngaus
  integer(ip), intent(out) :: ierro
  real(rp),    intent(in)  :: s
  real(rp),    intent(out) :: shaga(ngaus)

  ierro=0
  if(ngaus==1) then
     shaga(1)=1.0_rp
  else if(ngaus==2) then
     shaga(1)= 0.5_rp*sqrt(3.0_rp)*(1.0_rp/sqrt(3.0_rp)-s)
     shaga(2)= 0.5_rp*sqrt(3.0_rp)*(1.0_rp/sqrt(3.0_rp)+s)
  else if(ngaus==3) then
     shaga(1)= 5.0_rp/6.0_rp*(s-sqrt(0.6_rp))*s
     shaga(2)=-5.0_rp/3.0_rp*(s*s-0.60_rp)
     shaga(3)= 5.0_rp/6.0_rp*(s+sqrt(0.6_rp))*s
  else
     ierro=1
  end if
  
end subroutine shaga1
