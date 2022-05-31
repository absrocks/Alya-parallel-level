subroutine ruyclo(ndime,ngaus,posgp,weigp,ierro)

  !-----------------------------------------------------------------------
  ! 
  ! There is no closed quadrature rule for the pyramid as the Jacobian
  ! is zero at the apex!    
  ! 
  !----------------------------------------------------------------------- 
  use      def_kintyp
  implicit none
  integer(ip), intent(in)  :: ndime,ngaus
  integer(ip), intent(out) :: ierro
  real(rp),    intent(out) :: posgp(ndime,ngaus),weigp(ngaus)
  integer(ip)              :: ii,m
  real(rp)                 :: g1,g2,g3,g4,g5,gg2(2)
  real(rp)                 :: j,k,wg9(3)
  real(rp)                 :: w1,w2,w3,w4,ww2(2)
  real(rp)                 :: jk(2,20),jk4(2,20),jk9(2,9)

  jk(1,1)  = -1.0_rp 
  jk(2,1)  = -1.0_rp
  jk(1,2)  =  1.0_rp
  jk(2,2)  = -1.0_rp
  jk(1,3)  =  1.0_rp
  jk(2,3)  =  1.0_rp
  jk(1,4)  = -1.0_rp
  jk(2,4)  =  1.0_rp

  ierro=1
  posgp=0.0_rp
  weigp=0.0_rp

  if( ngaus == 5 ) then
     g1 = 8.0_rp*sqrt(2.0_rp/15.0_rp)/5.0_rp
     do ii = 1,4
        j = jk(1,ii)
        k = jk(2,ii)
        posgp(1,ii) = j
        posgp(2,ii) = k
        posgp(3,ii) =-2.0_rp/3.0_rp
        weigp(  ii) = 81.0_rp/100.0_rp           
     end do
     posgp(1,5) = 0.0_rp
     posgp(2,5) = 0.0_rp
     posgp(3,5) = 2.0_rp/5.0_rp
     weigp(  5) = 125.0_rp/27.0_rp  
  end if


end subroutine ruyclo
