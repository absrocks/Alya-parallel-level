subroutine rad_exacso(&
     itask,gpcod,gpden,gpdif,gprea,gpgrd,&
     gpvel,gpsou,extem,exteg,gprhs)
  !-----------------------------------------------------------------------
  !****f* Radiat/rad_exacso
  ! NAME 
  !    rad_exacso
  ! DESCRIPTION
  !    This routine computes the exact solution at the integration points
  !    ITASK = 1 ... Computes temp in and gradients 
  !    ITASK = 2 ... Computes force vector in to be applied as a source
  ! USES
  ! USED BY
  !    rad_elmope
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master, only       :  cutim
  use def_domain, only       :  ndime,kfl_naxis
  use def_radiat, only       :  kfl_exacs_rad,expar_rad
  implicit none
  integer(ip), intent(in)    :: itask
  real(rp),    intent(in)    :: gpden,gpdif,gprea,gpgrd(ndime)
  real(rp),    intent(in)    :: gpsou,gpvel(ndime),gpcod(ndime)
  real(rp),    intent(out)   :: extem,exteg(3)
  real(rp),    intent(inout) :: gprhs
  real(rp)                   :: x,y,z,t,u,v,w,r,k,a,dtdx,dtdy,dtdz,zetem
  real(rp)                   :: d2tdy2,d2tdx2,d2tdz2,dtdt,Q,dkdx,dkdy,dkdz
  real(rp)                   :: rmayo,rmeno
!!$  !
!!$  ! Initializations
!!$  ! 
!!$  zetem = 1.0e-12_rp
!!$  x     = gpcod(1)
!!$  dtdz  = 0.0_rp
!!$  if(ndime>=2) y = gpcod(2)
!!$  if(ndime==3) z = gpcod(3)
!!$
!!$  if(itask==2) then
!!$     dtdt    = 0.0_rp
!!$     dtdx    = 0.0_rp
!!$     dtdy    = 0.0_rp
!!$     dtdz    = 0.0_rp
!!$     d2tdx2  = 0.0_rp
!!$     d2tdy2  = 0.0_rp
!!$     d2tdz2  = 0.0_rp
!!$     dkdx    = gpgrd(1)
!!$     dkdy    = gpgrd(2)
!!$     dkdz    = 0.0_rp
!!$     if(ndime==3) dkdy = gpgrd(3)
!!$     u       = gpvel(1)
!!$     v       = gpvel(2)
!!$     w       = 0.0_rp
!!$     if(ndime==3) w = gpvel(3)
!!$     k       = gpdif
!!$     r       = gprea
!!$     Q       = gpsou
!!$  end if
!!$
!!$  if (kfl_exacs_rad==1)then
!!$
!!$     !
!!$     ! T=sin(pi*x)*sin(pi*y)*exp(x*y) in [0,1]x[0,1]
!!$     !
!!$     t      = sin(pi*x)*sin(pi*y)*exp(x*y)
!!$     dtdx   = sin(pi*y)*exp(x*y)*(pi*cos(pi*x)+y*sin(pi*x))
!!$     dtdy   = sin(pi*x)*exp(x*y)*(pi*cos(pi*y)+x*sin(pi*y))
!!$     d2tdx2 = sin(pi*y)*exp(x*y)*(2.0_rp*y*pi*cos(pi*x)&
!!$          &   -pi*pi*sin(pi*x)+y*y*sin(pi*x))
!!$     d2tdy2 = sin(pi*x)*exp(x*y)*(2.0_rp*x*pi*cos(pi*y)&
!!$          &   -pi*pi*sin(pi*y)+x*x*sin(pi*y))
!!$
!!$  else if (kfl_exacs_rad==2)then
!!$
!!$     a=expar_rad(1)
!!$     k=expar_rad(2)
!!$     r=expar_rad(3)
!!$     Q=expar_rad(4)
!!$     !
!!$     ! Solution of Advection-Diffusion-Reaction equation in [0,1]
!!$     !
!!$     if(r<=zetem .and. k>zetem .and. a>zetem) then
!!$        !
!!$        ! AD   :  T = 1/a*(exp(a/k)-1)*[(1-exp(a*x/k))+x/a]*Q
!!$        !
!!$        t=((1.0_rp-cosh(0.5_rp*a*x/k)-sinh(0.5_rp*x/k))/               &
!!$        &   (a*(cosh(0.5_rp*a/k)+sinh(0.5_rp*a/k)-1.0_rp))+x/a)*Q
!!$        dtdx=-(0.5_rp*a/k)*(cosh(0.5_rp*a*x/k)+sinh(0.5_rp*a*x/k))  &
!!$        &    /(a*(cosh(0.5_rp*a/k)+sinh(0.5_rp*a/k)-1.0_rp))+1.0_rp/a
!!$        d2tdx2=-0.25_rp*a*(cosh(0.5_rp*a*x/k)+sinh(0.5_rp*a*x/k))&
!!$        &    /(k*k*(cosh(0.5_rp*a/k)+sinh(0.5_rp*a/k)-1.0_rp))
!!$
!!$     else if(r<=zetem .and. k>zetem .and. a<=zetem) then
!!$        !
!!$        ! D    : T = [ x/(2*k)*(1-x) ]*Q
!!$        !
!!$        t      = (x/(2.0_rp*k)*(1.0_rp-x))*Q
!!$        dtdx   = 1.0_rp/(2.0_rp*k)-x/k 
!!$        d2tdx2 = -1.0_rp/k
!!$
!!$     else if(r<=zetem .and. k<=zetem .and. a>zetem) then
!!$        !
!!$        ! A   : T = [x/a]*Q
!!$        !
!!$        t      = (x/a)*Q
!!$        dtdx   = 1.0_rp/a
!!$        d2tdx2 = 0.0_rp
!!$
!!$     else if(k<=zetem .and. r>zetem .and. a>zetem ) then
!!$        !
!!$        ! AR  : T = [1/r*(1-exp(-r*x/a))]*Q
!!$        !
!!$        t=((1.0_rp-cosh(r*x/a)+sinh(r*x/a))/r)*Q
!!$        dtdx=(cosh(r*x/a)-sinh(r*x/a))/a
!!$        d2tdx2=(sinh(r*x/a)-cosh(r*x/a))*r/a
!!$
!!$     else if(k<=zetem .and. r>zetem .and. a<=zetem) then
!!$        !
!!$        ! R   : T = Q/s
!!$        !
!!$        t    = (1.0_rp/r)*Q
!!$        dtdx = 0.0_rp
!!$        d2tdx2=0.0_rp
!!$
!!$     else if(a<=zetem .and. k>zetem .and. r>zetem) then
!!$        !
!!$        ! DR  : T = [1/(r*(exp(r/k)^(1/2)-exp(r/k)^(1/2))*(exp((r/k)^(1/2)*x)-exp(-(r/k)^(1/2)*x))-1/r]*Q
!!$        !
!!$        t=((sinh(sqrt(r/k)*x)/(r*sinh(sqrt(r/k)))-1.0_rp/r))*Q
!!$        dtdx=(sqrt(r/k)*cosh(sqrt(r/k)*x)/(r*sinh(sqrt(r/k))))
!!$        d2tdx2=sinh(sqrt(r/k)*x)/(k*sinh(sqrt(r/k)))
!!$
!!$     else
!!$        !
!!$        ! ADR : T= 1/(s*a)*[(exp(r+)-1)/(exp(r-)-exp(r+)) * exp((r-)*x) -(exp(r-)-1)/(exp(r-)-exp(r+))*exp((r+)*x)+1]*Q
!!$        !
!!$        rmayo = a/(k*2.0_rp)+sqrt( (a/(k*2.0_rp))**2.0_rp + r/k )
!!$        rmeno = a/(k*2.0_rp)-sqrt( (a/(k*2.0_rp))**2.0_rp + r/k )
!!$         t=((((cosh(rmayo)+sinh(rmayo)-1.0_rp)*(cosh(rmeno*x)+sinh(rmeno*x))-       &
!!$         &  (cosh(rmeno)+sinh(rmeno)-1.0_rp)*(cosh(rmayo*x)+sinh(rmayo*x)))/       &
!!$         &  (cosh(rmeno)+sinh(rmeno)-cosh(rmayo)-sinh(rmayo))+1.0_rp)*(1.0_rp/(r*a)))*Q
!!$        dtdx=(((cosh(rmayo)+sinh(rmayo)-1.0_rp)*rmeno*(cosh(rmeno*x)+sinh(rmeno*x))-  &
!!$        &   (cosh(rmeno)+sinh(rmeno)-1.0_rp)*rmayo*(cosh(rmayo*x)+sinh(rmayo*x)))/    &
!!$        &  ( cosh(rmeno)+sinh(rmeno)-cosh(rmayo)-sinh(rmayo)))*(1.0_rp/(r*a))
!!$        d2tdx2=(1.0_rp/(r*a))*(((cosh(rmayo)+sinh(rmayo)-1.0_rp)*((cosh(rmeno*x)+sinh(rmeno*x))*rmeno*rmeno) &
!!$        &   -(cosh(rmeno)+sinh(rmeno)-1.0_rp)*((cosh(rmayo*x)+sinh(rmayo*x))*rmayo*rmayo))/      &
!!$        &   (cosh(rmeno)+sinh(rmeno)-cosh(rmayo)-sinh(rmayo)))
!!$
!!$     end if
!!$
!!$  else if(kfl_exacs_rad==3)then
!!$     !
!!$     ! T=2*x+3*y in [0,1]x[0,1]
!!$     !
!!$     t      = 2.0_rp*x+3.0_rp*y
!!$     dtdx   = 2.0_rp
!!$     dtdy   = 3.0_rp
!!$
!!$  else if(kfl_exacs_rad==4)then
!!$     !
!!$     ! T=2.0*exp(-(20.0*(x-0.5))**2)*sin(x)*sin(pi*x)*sin(y)*sin(pi*y)   
!!$     !
!!$     t      =  2.0_rp*exp(-(20.0_rp*(x-0.5_rp))**2.0_rp)*sin(x)*sin(pi*x)*sin(y)*sin(pi*y) 
!!$
!!$     dtdx   =  2.0_rp*(-800.0_rp*x+400.0_rp)*exp((-400.0_rp*((x-0.5_rp))**2))*sin(x)*sin(pi*x)*sin(y)*sin(pi*y)&
!!$          &   +2.0_rp*exp((-400.0_rp*((x-0.5))**2))*cos(x)*sin(pi*x)*sin(y)*sin(pi*y)    &
!!$          &   +2.0_rp*exp((-400.0_rp*((x-0.5_rp))**2))*sin(x)*cos(pi*x)*pi*sin(y)*sin(pi*y)
!!$
!!$     dtdy   =  2.0_rp*exp(-(20.0_rp*(x-0.5_rp))**2.0_rp)*sin(x)*sin(pi*x)&
!!$          &    *(cos(y)*sin(pi*y)+pi*sin(y)*cos(pi*y))
!!$
!!$     d2tdx2 = -1602.0_rp*exp((-400.0_rp*((x-0.5))**2))*sin(x)*sin(pi*x)*sin(y)*sin(pi*y)                            &
!!$          &   +2.0_rp*((-800.0_rp*x+400.0_rp))**2 *exp((-400.0_rp*((x-0.5))**2))*sin(x)*sin(pi*x)*sin(y)*sin(pi*y)  &
!!$          &   +4.0_rp*(-800.0_rp*x+400.0_rp)*exp((-400.0_rp*((x-0.5))**2))*cos(x)*sin(pi*x)*sin(y)*sin(pi*y)        &
!!$          &   +4.0_rp*(-800.0_rp*x+400.0_rp)*exp((-400.0_rp*((x-0.5))**2))*sin(x)*cos(pi*x)*pi*sin(y)*sin(pi*y)     &
!!$          &   +4.0_rp*exp((-400.0_rp*((x-0.5))**2))*cos(x)*cos(pi*x)*pi*sin(y)*sin(pi*y)                            &
!!$          &   -2.0_rp*exp((-400.0_rp*((x-0.5))**2))*sin(x)*sin(pi*x)*pi*pi*sin(y)*sin(pi*y)
!!$
!!$     d2tdy2 = -2.0_rp*exp((-400.0_rp*((x-0.5_rp))**2))*sin(x)*sin(pi*x)*sin(y)*sin(pi*y)          &
!!$          &   +4.0_rp*exp((-400.0_rp*((x-0.5))**2))* sin(x)*sin(pi*x)*cos(y)*cos(pi*y)            &
!!$          &   *pi-2.0_rp*exp((-400.0_rp*((x-0.5))**2))*sin(x)*sin(pi*x)*pi**2*sin(y)* sin(pi*y)
!!$
!!$  else if (kfl_exacs_rad==5)then
!!$     !
!!$     ! T=2x+3
!!$     !
!!$     t      = 2.0_rp*x*x+1.0_rp
!!$     dtdx   = 4.0_rp*x
!!$
!!$  else if (kfl_exacs_rad==6)then
!!$     !
!!$     ! T=t*(2x+3y+4)
!!$     !
!!$     t      = cutim*(2.0_rp*x+3.0_rp*y+4.0_rp)
!!$     dtdx   = cutim*2.0_rp
!!$     dtdy   = cutim*3.0_rp
!!$     dtdt   = 2.0_rp*x+3.0_rp*y+4.0_rp
!!$
!!$ else if (kfl_exacs_rad==5)then
!!$     !
!!$     ! T=2x^2+1
!!$     !
!!$     t      = 2.0_rp*x*x+1.0_rp
!!$     dtdx   = 4.0_rp*x
!!$     d2tdx2 = 4.0_rp
!!$
!!$  else if (kfl_exacs_rad==7)then
!!$     !
!!$     ! T=2x^2+3x^2
!!$     !
!!$     t      = 2.0_rp*x*x+3.0_rp*y*y
!!$     dtdx   = 4.0_rp*x
!!$     d2tdx2 = 4.0_rp
!!$     dtdy   = 6.0_rp*y
!!$     d2tdy2 = 6.0_rp
!!$
!!$  else if (kfl_exacs_rad==8)then
!!$     !
!!$     ! T=2x^3+3y^3
!!$     !
!!$     t      = 2.0_rp*x*x*x+3.0_rp*y*y*y
!!$     dtdx   = 6.0_rp*x*x
!!$     d2tdx2 = 12.0_rp*x
!!$     dtdy   = 9.0_rp*y*y
!!$     d2tdy2 = 18.0_rp*y
!!$
!!$  else if(kfl_exacs_rad==9)then
!!$
!!$     !
!!$     ! T=x^3+y^3
!!$     !
!!$     t      = x*x*x+y*y*y
!!$     dtdx   = 3.0_rp*x*x
!!$     d2tdx2 = 6.0_rp*x
!!$     dtdy   = 3.0_rp*y*y
!!$     d2tdy2 = 6.0_rp*y
!!$
!!$  else if(kfl_exacs_rad==10)then
!!$     !
!!$     ! T=x+2*y+3*z
!!$     !
!!$     t      = x+2.0_rp*y+3.0_rp*z
!!$     dtdx   = 1.0_rp
!!$     dtdy   = 2.0_rp
!!$     dtdz   = 3.0_rp
!!$
!!$  end if
!!$
!!$  if(itask==1) then
!!$     !
!!$     ! Exact temp and gradients
!!$     !
!!$     extem    = t
!!$     exteg(1) = dtdx
!!$     exteg(2) = dtdy
!!$     exteg(3) = dtdz
!!$
!!$  else if(itask==2) then
!!$     !
!!$     ! Force term GPRHS= rho*cp*[dT/dt+u.grad(T)]-div[k*grad(T)]+r*T
!!$     !
!!$     gprhs = gprhs + gpden*( dtdt + u*dtdx + v*dtdy + w*dtdz ) &
!!$          &  - dkdx*dtdx - dkdy*dtdy - dkdz*dtdz &
!!$          &  - k*(d2tdx2+d2tdy2+d2tdz2) + r*t 
!!$     !
!!$     ! Axisymmetric problem: -k/x*dT/dx
!!$     !
!!$     if(kfl_naxis==1.and.abs(x)/=0.0_rp) gprhs = gprhs - k/x*dtdx
!!$
!!$  end if

end subroutine rad_exacso
