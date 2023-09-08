subroutine chm_velfun(npoin,coord,vefun)
  !-----------------------------------------------------------------------
  !****f* Temper/chm_velfun
  ! NAME
  !   chm_velfun
  ! DESCRIPTION
  !   Compute velocity according to the function number 
  ! INPUT
  !   KFL_ADVEC_CHM ... Function
  !   COORD ........... Coordinates
  !   NPOIN ........... Number of nodes or Gauss points
  ! OUTPUT 
  !   VEFUN ........... Velocity
  ! USES
  ! USED BY
  !    chm_elmope
  !***
  !-----------------------------------------------------------------------
  use def_parame, only       :  pi
  use def_kintyp, only       :  ip,rp 
  use def_chemic, only       :  kfl_advec_chm
  use def_master, only       :  veloc,cutim
  use def_domain, only       :  ndime
  implicit none
  integer(ip), intent(in)    :: npoin
  real(rp),    intent(in)    :: coord(ndime,npoin)
  real(rp),    intent(out)   :: vefun(ndime,npoin)
  integer(ip)                :: idime,ipoin
  real(rp)                   :: x,y,theta,w

  if( kfl_advec_chm == 0 ) then
     do ipoin=1,npoin
        do idime=1,ndime
           vefun(idime,ipoin)= 0.0_rp
        end do
     end do
     
  else if( kfl_advec_chm == 1 ) then
     do ipoin=1,npoin
        do idime=1,ndime
           vefun(idime,ipoin)= veloc(idime,ipoin,1)
        end do
     end do
     
  else if( kfl_advec_chm == 2 ) then
     do ipoin=1,npoin
        x=coord(1,ipoin)
        y=coord(2,ipoin)
        vefun(1,ipoin)= 0.5_rp*(1.0_rp-x*x)*(1.0_rp+y)
        vefun(2,ipoin)=-0.5_rp*x*(4.0_rp-(1.0_rp+y)**2)
     end do

  else if( kfl_advec_chm == 3 ) then
     do ipoin=1,npoin
        vefun(    1,ipoin) =  1.0_rp
        vefun(    2,ipoin) = -1.0_rp   
        vefun(ndime,ipoin) = -1.0_rp   
     end do

  else if( kfl_advec_chm == 4 ) then
     do ipoin=1,npoin
        x=coord(1,ipoin)
        y=coord(2,ipoin)
        vefun(1,ipoin)=  pi/2.0_rp*(y-0.2_rp)
        vefun(2,ipoin)= -pi/2.0_rp*(x-0.8_rp)
     end do

  else if( kfl_advec_chm == 5 ) then
     do ipoin=1,npoin
        x=coord(1,ipoin)
        y=coord(2,ipoin)
        vefun(1,ipoin)=  sin(pi*x)**2*sin(2.0_rp*pi*y)*cos(pi*cutim /5.0_rp)
        vefun(2,ipoin)= -sin(pi*y)**2*sin(2.0_rp*pi*x)*cos(pi*cutim /5.0_rp) 
     end do

  else if( kfl_advec_chm == 6 ) then
     w = 2.0_rp*pi/1000.0_rp
     do ipoin=1,npoin
        x=coord(1,ipoin)
        y=coord(2,ipoin)
        vefun(1,ipoin)=  -w*(y-0.5_rp*(2.0_rp))
        vefun(2,ipoin)=   w*(x-0.5_rp*(2.0_rp))
     end do
  else if( kfl_advec_chm == 7 ) then
     do ipoin=1,npoin
        vefun(1,ipoin)=  0.1_rp
        vefun(2,ipoin)=   0.0_rp
     end do
  end if

end subroutine chm_velfun
