subroutine shaga2(s,t,ltopo,ngaus,shaga,ierro)

  !-----------------------------------------------------------------------
  !
  ! This routine evaluates shape functions associated to Gauss points
  ! for 2D. The cases available so far are:
  !
  ! BRICKS    -->   NGAUS =   1   4   9  
  ! SIMPLICES -->   NGAUS =   1   3   4  6  7
  !
  !-----------------------------------------------------------------------
  use      def_kintyp
  implicit none
  integer(ip), intent(in)  :: ngaus,ltopo
  integer(ip), intent(out) :: ierro
  real(rp),    intent(in)  :: s,t
  real(rp),    intent(out) :: shaga(ngaus)
  !
  ! Quadrilateral and hexahedral elements
  !
  if(ngaus==1) then
     shaga(1)=1.0_rp
  else if(ltopo==0 .and. ngaus==4) then
     shaga(1)= .75_rp*(s-1.0_rp/sqrt(3.0_rp))*(t-1.0_rp/sqrt(3.0_rp))
     shaga(2)=-.75_rp*(s-1.0_rp/sqrt(3.0_rp))*(t+1.0_rp/sqrt(3.0_rp))
     shaga(3)=-.75_rp*(s+1.0_rp/sqrt(3.0_rp))*(t-1.0_rp/sqrt(3.0_rp))
     shaga(4)= .75_rp*(s+1.0_rp/sqrt(3.0_rp))*(t+1.0_rp/sqrt(3.0_rp))
  else if(ltopo==0 .and. ngaus==9) then
     shaga(1)= 25.0_rp/36.0_rp*(s-sqrt(.6))*(t-sqrt(.6))*s*t
     shaga(2)=-25.0_rp/18.0_rp*(s-sqrt(.6))*(t*t-.6)*s
     shaga(3)= 25.0_rp/36.0_rp*(s-sqrt(.6))*(t+sqrt(.6))*s*t
     shaga(4)=-25.0_rp/18.0_rp*(s*s-.6)*(t-sqrt(.6))*t
     shaga(5)= 25.0_rp/9.0_rp*(s*s-.6)*(t*t-.6)
     shaga(6)=-25.0_rp/18.0_rp*(s*s-.6)*(t+sqrt(.6))*t
     shaga(7)= 25.0_rp/36.0_rp*(s+sqrt(.6))*(t-sqrt(.6))*s*t
     shaga(8)=-25.0_rp/18.0_rp*(s+sqrt(.6))*(t*t-.6)*s
     shaga(9)= 25.0_rp/36.0_rp*(s+sqrt(.6))*(t+sqrt(.6))*s*t
     !
     ! Triangular and tetrahedral elements
     !        
  else if (ltopo==1 .and. ngaus==3) then
     shaga(1)=2.0_rp*s-1.0_rp/3.
     shaga(2)=2.0_rp*t-1.0_rp/3.
     shaga(3)=2.0_rp*(1.0_rp-s-t)-1.0_rp/3.0_rp
  else if (ltopo==1 .and. ngaus==4) then
     shaga(1)=(-45.0_rp*(s+t)+225.0_rp*s*t+9.0_rp)/4.0_rp
     shaga(2)=(  5.0_rp*(s+t)- 75.0_rp*s*t+5.0_rp)/4.0_rp
     shaga(3)=(25.0_rp*s+15.0_rp*t-75.0_rp*s*t-5.0_rp)/4.0_rp
     shaga(4)=(15.0_rp*s+25.0_rp*t-75.0_rp*s*t-5.0_rp)/4.0_rp
  else if (ltopo==1 .and. ngaus==6) then
     shaga( 1)= 0.13855958741e+00_rp-0.19897337353e+01_rp*s&
          + 0.16597397329e+00_rp*t-0.16597397329e+00_rp*s*t&
          + 0.37248334214e+01_rp*s*s-0.16597397329e+00_rp*t*t
     shaga( 2)= 0.18736592735e+01_rp-0.54599331075e+01_rp*s&
          -0.54599331075e+01_rp*t+ 0.76156408161e+01_rp*s*t&
          + 0.37248334214e+01_rp*s*s+ 0.37248334214e+01_rp*t*t
     shaga( 3)= 0.13855958741e+00_rp+ 0.16597397329e+00_rp*s&
          -0.19897337353e+01_rp*t-0.16597397329e+00_rp*s*t&
          -0.16597397329e+00_rp*s*s+ 0.37248334214e+01_rp*t*t
     shaga( 4)=-0.63855958741e+00_rp+ 0.40859490520e+00_rp*s&
          + 0.79963036869e+01_rp*t-0.79963036869e+01_rp*s*t&
          + 0.35630540870e+00_rp*s*s-0.79963036869e+01_rp*t*t
     shaga( 5)= 0.12634072649e+00_rp-0.11212057226e+01_rp*s&
          -0.11212057226e+01_rp*t+ 0.87089145043e+01_rp*s*t&
          + 0.35630540870e+00_rp*s*s+ 0.35630540870e+00_rp*t*t
     shaga( 6)=-0.63855958741e+00_rp+ 0.79963036869e+01_rp*s&
          + 0.40859490520e+00_rp*t-0.79963036869e+01_rp*s*t&
          -0.79963036869e+01_rp*s*s+ 0.35630540870e+00_rp*t*t
  else if (ltopo==1 .and. ngaus==7 ) then
     shaga( 1)= 0.46660223518e-14_rp-0.22500000000e+01_rp*s&
          -0.22500000000e+01_rp*t+ 0.49500000000e+02_rp*s*t&
          + 0.22500000000e+01_rp*s*s+ 0.22500000000e+01_rp*t*t&
          -0.47250000000e+02_rp*(s*s*t+s*t*t)
     shaga( 2)=-0.76639777949e+00_rp+ 0.27797375097e+01_rp*s&
          + 0.87162291828e+01_rp*t-0.32857759237e+02_rp*s*t&
          -0.21106850116e+01_rp*s*s-0.87162291828e+01_rp*t*t&
          + 0.24141530054e+02_rp*(s*s*t+s*t*t)
     shaga( 3)=-0.76639777949e+00_rp+ 0.87162291828e+01_rp*s&
          + 0.27797375097e+01_rp*t-0.32857759237e+02_rp*s*t&
          -0.87162291828e+01_rp*s*s-0.21106850116e+01_rp*t*t&
          + 0.24141530054e+02_rp*(s*s*t+s*t*t)
     shaga( 4)=-0.97345281425e-01_rp+ 0.14416325135e+01_rp*s&
          + 0.14416325135e+01_rp*t-0.19646670894e+02_rp*s*t&
          -0.21106850116e+01_rp*s*s-0.21106850116e+01_rp*t*t&
          + 0.24141530054e+02_rp*(s*s*t+s*t*t)
     shaga( 5)= 0.26639777949e+00_rp-0.30297375097e+01_rp*s&
          -0.96622918276e+00_rp*t+ 0.93577592368e+01_rp*s*t&
          + 0.48606850116e+01_rp*s*s+ 0.96622918276e+00_rp*t*t&
          -0.83915300541e+01_rp*(s*s*t+s*t*t)
     shaga( 6)= 0.26639777949e+00_rp-0.96622918276e+00_rp*s&
          -0.30297375097e+01_rp*t+ 0.93577592368e+01_rp*s*t&
          + 0.96622918276e+00_rp*s*s+ 0.48606850116e+01_rp*t*t&
          -0.83915300541e+01_rp*(s*s*t+s*t*t)
     shaga( 7)= 0.20973452814e+01_rp-0.66916325135e+01_rp*s&
          -0.66916325135e+01_rp*t+ 0.17146670894e+02_rp*s*t&
          + 0.48606850116e+01_rp*s*s+ 0.48606850116e+01_rp*t*t&
          -0.83915300541e+01_rp*(s*s*t+s*t*t)
  else
     ierro=1
     !call runend('SHAGA2: INTERPOLATION NOT AVAILABLE')
  end if

end subroutine shaga2
      
