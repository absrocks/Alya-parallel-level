program unitt_maths_linear_regressions 
  
   use def_kintyp_basic,      only : ip, rp
   use mod_maths_basic,       only : maths_weighted_linear_regression
   use mod_maths_basic,       only : maths_linear_regression
  
  implicit none
 
   real(rp),    pointer     :: X(:)
   real(rp),    pointer     :: Y(:)
   real(rp),    pointer     :: W(:)
   real(rp)                 :: a, b
   real(rp), parameter :: a0=-0.621636457587501_rp
   real(rp), parameter :: b0=7.99188885777447_rp
   real(rp), parameter :: a1=-0.544011033702302_rp 
   real(rp), parameter :: b1=12.9875870849100_rp 


   nullify(X,Y,W)
   allocate(X(8),Y(8),W(8))

   X = (/ 1.5_rp,2.0_rp,5.7_rp,3.2_rp,11.3_rp,12.4_rp,100.12339890_rp,12.3234234_rp /) 
   Y = (/ 7.5_rp,25.0_rp,4.7_rp,1.2_rp,13.3_rp,1.4_rp,10.12339890_rp,111.3234234_rp /) 
   W = (/ 3.1_rp,5.3232_rp,11.9_rp,24.0_rp,1.111_rp,0.04343_rp,4.1232_rp,11.0434_rp /)

   call maths_weighted_linear_regression(X,Y,8_ip,a,b,W,6_ip)
   if( abs(a-a0) > 1.e-12 .or. abs(b-b0) > 1.e-12) stop 3;
   
   call maths_linear_regression(X,Y,8_ip,a,b,6_ip)
   if( abs(a-a1) > 1.e-12 .or. abs(b-b1) > 1.e-12) stop 4;

   deallocate(X,Y,W)
   stop
  
end program unitt_maths_linear_regressions
