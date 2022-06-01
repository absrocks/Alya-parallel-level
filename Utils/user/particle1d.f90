!-----------------------------------------------------------------------
!> @{
!> @file    particle1d.f90
!> @author  Guillaume Houzeaux
!> @date    26/10/2017
!> @brief   Particle trajectory
!> @details Compute a particle trajectory submitted to a drag law.
!>
!>          d up
!>          ----
!>           dt
!>          Main assumptions:
!>          - uf is constant
!>
!>          Input
!>          -----
!>          rhof ... Constant fluid density
!>          muf .... Constant fluid viscosity
!>          uf ..... Constant fluid velocity
!>          dp ..... Particle diameter
!>          rhop ... Particle density
!>          ap ..... Particle acceleration
!>          up ..... Particle velocity
!>          xp ..... Particle position
!>          dt ..... Time step
!>
!>          Output (each column):
!>          ---------------------
!>          1. Time
!>          2. Particle velocity
!>          3. tau
!>          6. Stokes number
!>          A particle with a low Stokes number follows fluid streamlines
!>          (perfect advection), while a particle with a large Stokes number
!>          is dominated by its inertia and continues along its initial trajectory.
!>          where is theparticle density, is the particle diameter and is the
!>          gas dynamic viscosity.
!> @}
!-----------------------------------------------------------------------

program particle1d
  implicit none
  integer, parameter  :: ip = 4             ! 4-byte integer
  integer, parameter  :: rp = 8             ! Double precision

  integer(ip) :: ii
  real(rp)    :: alphad,Cd,Re,CdRe,dt,psip
  real(rp)    :: dCdRedRe,Du,dp,t,ap,eps,tauinv
  real(rp)    :: tau,up(2),uf,rhop,muf,rhof,tf
  real(rp)    :: rper1,rper2,percent,Stk,uu,exact,up0
  real(rp)    :: tauStk

  open(unit=10,file='u.txt',status='unknown')
  write(10,1)

  t      = 0.0_rp  
  eps    = epsilon(1.0_rp)
  rhof   = 1.2047_rp
  muf    = 1.82e-05_rp
  uf     = 1.0_rp

  up     = 0.0_rp
  psip   = 1.0_rp
  !
  ! 1 Micro
  ! 
  dp     = 1.0e-6_rp
  rhop   = 1.0e3_rp
  dt     = 1.0e-11_rp
  tf     = 1.0e-4_rp
  !
  ! 100 micros
  !
  dp     = 1.0e-4_rp
  rhop   = 1.0e3_rp
  dt     = 1.0e-8_rp
  tf     = 1.0e0_rp 
  !
  ! 1 Nano
  !
  !dp     = 1.0e-9_rp
  !rhop   = 1.2047_rp
  !dt     = 1.0e-17_rp
  !tf     = 1.0e-12_rp 
  !
  ! up^{n+1} = up^n + (uf-up^n} * dt / tau
  !
  alphad  = 0.75_rp * muf / (rhop*dp**2)
  uu      = 1.0_rp ! max(uf,up(1))
  rper2   = 0.0_rp
  percent = 1.0e-3_rp 
  !percent = 0.0_rp
  tauStk  = rhop*dp**2/(18.0_rp*muf)
  up0     = uf-up(1)

  do while( t < tf )
     Du     = abs(up(1)-uf)     
     call dragfo(2_ip,Du,muf,rhof,dp,psip,CdRe,Cd,Re,dCdRedRe)     
     tauinv = alphad * CdRe
     tau    = 1.0_rp / max(eps,tauinv)
     ap     = ( uf - up(1) ) * tauinv
     Stk    = tau * uf / dp
     up(2)  = up(1) + ap * dt
     t      = t + dt
     exact  = uf - up0 * exp(-t/tauStk)

     rper1 = 100.0_rp * t / tf
     if( rper1 >= rper2 + percent ) then
        rper2 = rper2 + percent
        write(10,'(10(1x,e13.6))') t,abs(up(2)-uf)/uf,tau,Re,CdRe,Stk,exact/uu,alphad
     end if

     up(1)  = up(2)

  end do

  close(unit=10)

1 format('#            t',&
       & '            up',&
       & '           tau',&
       & '            Re',&
       & '         Cd*Re',&
       & '        Stokes',&
       & '         Exact')
end program particle1d


subroutine dragfo(imode,u,mu,rho,d,psi,CdRe,Cd,Re,dCdRedRe)
  !-----------------------------------------------------------------------
  !****f* dragfo/dragfo
  ! NAME
  !    dragfo
  ! DESCRIPTION
  !    By definition:
  !    Fd = 3/4 * rho_a * u^2 * V/d * Cd
  !    By using V = (4/3)*pi*(d/2)^3, we end up with:
  !
  !    Fd = rho_a * u^2 * pi * d^2/8 * Cd
  !       = ( 1/2 * rho_a * u^2 ) * ( pi * (d/2)^2 ) * Cd
  !
  !    To be used for a drag model of the kind:
  !
  !    Fd(i) = - pi * mu * d / 8.0_rp * Cd * Re * ( u(i)-uf(i) )
  !    ad(i) =   Fd(i) / (rho_p*pi*d^3*)/6)
  !           =  visfl / rho_p * 0.75 * Cd * Re / d^2 * ( u(i)-uf(i) )
  !
  !    or, equivalently:
  !
  !    Fd = pi * mu * d/8 * (Cd*Re) * u
  !
  !    The associated terminal velocity for a sphere is:
  !    ut = sqrt( 4*g*d/(3*Cd) * (rho_p-rho_a)/rho_a )
  !
  !    CD FORMULA:
  !    -----------
  !
  !    1. CHENG model.
  !       Nian-Sheng Cheng, 'Comparison of formulas for drag coefficient
  !       and settling velocity of spherical particles',
  !       Powder Technology, 2008.
  !
  !       Cd = 24/Re * ( 1 + 0.27 * Re ) ** 0.43 + 0.4 * ( 1 - exp(-0.04*Re**0.38) )
  !
  !    2. GANSER model.
  !       Ganser, G., 1993. A rational approach to drag prediction of spherical and nonspherical
  !       particles. Powder Technology 77, 143-152.
  !
  !       Cd = 24/(Re*k1) * ( 1 + 0.1118*(Re*k1*k2)^0.6567 ) + 0.4305*k2/ [ 1 + 3305/(Re*k1*k2) ]
  !       k1 = 3/(1+2*p^-0.5)
  !       k2 = 10^(1.84148*(-log10(p))^0.5743)
  !       p  = sphericity (=1 for a sphere)
  !
  !    3. ARASTOOPOUR model.
  !       Arastoopour, H., Wang, C., Weil, S., 1982. Particle-particle interaction
  !       force in a dilute gas-solid system. Chemical Engineering Science 37 (9), 1379-1386.
  !
  !       if( Re <= 1000 ) then
  !           Cd = 24.0/Re * ( 1 + 0.15 * Re ** 0.687 )
  !       else
  !           Cd = 0.44
  !       end if
  !
  !    4. WILSON model.
  !       Wilson, L., Huang, T., 1979. The influence of shape on the atmospheric settling velocity
  !       of volcanic ash particles. Earth and Planetary Science Letters 44, 311-324.
  !
  !       if( Re <= 100 ) then
  !           Cd = 24/Re*p^-0.828 + 2 sqrt( 1.0-p )
  !       else if( 100 <= Re <= 1000 ) then
  !           Cd = 1 - ( 1-Cd|Re=10^2 )/900 * (10^3-Re)
  !       else
  !           Cd = 1
  !       end if
  !       p = (b+c)/2a = particle aspect ratio (a>b>c are the particle semi-axes); p=1 for sphere
  !
  !    5. TURTON and LEVENSPIEL model.
  !       Turton, R., and O. Levenspiel. 1986. A short note on drag correlation for spheres. Powder Technolo-
  !       gy Journal, 47, 83.
  !       Very similar to Cheng's. Valid for Re < 2.6x10^5
  !
  !
  !
  !                                 Cd for different formula
  !    Cd
  !      10000 ++-+---+-+--+---+-++--+--+-++--+---+-++--+--+-++--+---+-++-+---+-++
  !            +        +         +        +         +        +    GANSER ****** +
  !            **                                                  WILSON ###### +
  !       1000 +***                                                             ++
  !            +   ****                                                          +
  !            |      ***                                                        |
  !            +         ***                                                     +
  !        100 ++           ***                                                 ++
  !            +              ****                                               +
  !            +                 #***                                            +
  !         10 ++                   ****                                        ++
  !            +                      ##***                                      +
  !            |                          #****                                  |
  !            +                            ##*****                              +
  !          1 ++                              ### ******  #######################
  !            +                                  ### ##**************************
  !            +        +         +        +        ###       +         +        +
  !        0.1 ++-+---+-+--+---+-++--+--+-++--+---+-++--+--+-++--+---+-++-+---+-++
  !           0.01     0.1        1        10       100      1000     10000    100000
  !                                            Re
  !
  !
  ! USED BY
  !    domain
  !***
  !-----------------------------------------------------------------------
  implicit none
  integer, parameter  :: ip = 4             ! 4-byte integer
  integer, parameter  :: rp = 8             ! Double precision
  integer(ip), intent(in)  :: imode
  real(rp),    intent(in)  :: u          !< Relative velocity |u_fluid-u_particle|
  real(rp),    intent(in)  :: mu         !< Fluid viscosity
  real(rp),    intent(in)  :: rho        !< Fluid density
  real(rp),    intent(in)  :: d          !< Particle diameter
  real(rp),    intent(in)  :: psi        !< Sphericity
  real(rp),    intent(out) :: CdRe       !< Drag * Reynolds
  real(rp),    intent(out) :: Cd         !< Drag coefficient
  real(rp),    intent(out) :: Re         !< Reynolds number
  real(rp),    intent(out) :: dCdRedRe   !< d(Cd*Re) / dRe
  real(rp)                 :: k1,k2,Re2
  real(rp)                 :: zeror
  !
  ! Particle Reynolds number
  !
  zeror    = epsilon(1.0_rp)
  Re       = rho * u * d / mu + zeror
  dCdRedRe = 0.0_rp

  select case ( imode )

  case ( 1_ip )
     !
     ! Cheng model
     !
     CdRe =   24.0_rp * ( 1.0_rp + 0.27_rp * Re )**0.43_rp &
          &  + 0.4_rp * Re * ( 1.0_rp - exp(-0.04_rp*Re**0.38_rp) )

  case ( 2_ip )
     !
     ! Ganser model
     !
     Re2  = min(1.0e5_rp,Re)
     k1   = 3.0_rp / ( 1.0_rp + 2.0_rp/sqrt(psi) )
     k2   = 10.0_rp**( 1.84148_rp * (-log10(psi))**0.5743_rp )
     CdRe =   24.0_rp/k1 * ( 1.0_rp + 0.1118_rp*(Re2*k1*k2)**0.6567_rp ) &
          & + 0.4305_rp*Re2*Re2*k2/ ( Re2 + 3305.0_rp/(k1*k2) )
     if( Re > zeror ) then
        dCdRedRe =   24.0_rp/k1 * ( 0.1118_rp*0.6567_rp*(k1*k2)**0.6567_rp*Re2**(1.0_rp-0.6567_rp) ) &
             &    + 0.4305_rp * Re2 * 2.0_rp * k2 / ( Re2 + 3305.0_rp/(k1*k2) ) &
             &    - 0.4305_rp * Re2 * Re2    * k2 / ( Re2 + 3305.0_rp/(k1*k2) ) ** 2
     end if

  case ( 3_ip )
     !
     ! Arastoopour model
     !
     if( Re <= 1000.0_rp ) then
        CdRe   = 24.0_rp * ( 1.0_rp + 0.15_rp * Re ** 0.687_rp )
        if( Re > zeror ) then
           dCdRedRe = -24.0_rp / Re * ( 1.0_rp + 0.15_rp * Re ** 0.687_rp ) &
                + 24.0_rp * ( 0.15_rp * 0.687_rp * Re ** (0.687_rp-1.0_rp) )
        end if
     else
        CdRe     = 0.44_rp * Re
        dCdRedRe = 0.0_rp
     end if

  case ( 4_ip )
     !
     ! Wilson model
     !
     if( Re <= 100_rp ) then
        CdRe = 24.0_rp*psi**(-0.828_rp) + 2.0_rp * Re * sqrt( 1.0_rp-psi )
     else if( Re > 100.0_rp .and. Re <= 1000.0_rp ) then
        k1   = 24.0_rp/(100.0_rp)*psi**(-0.828_rp) + 2.0_rp * sqrt( 1.0_rp-psi )
        CdRe = Re * ( 1.0_rp - ( 1.0_rp-k1 )/900 * (1000.0_rp-Re) )
     else
        CdRe = Re
     end if

  case ( 5_ip )
     !
     ! Turton and Levenspiel model
     !
     CdRe = 24.0_rp * (1.0_rp+0.173_rp*Re**0.657_rp) + 0.413_rp*Re/(1.0_rp+11630.0_rp*Re**(-1.09_rp))

  case ( 6_ip )
     !
     ! Stokes model, dCdRedRe = dCd / dRe * Re
     !
     CdRe     =  24.0_rp
     dCdRedRe =   0.0_rp

  case ( 7_ip )
     !
     ! Dallavalle
     !
     CdRe     = ( sqrt(Re)*0.63_rp + 4.8_rp )**2
     dCdRedRe = 0.63_rp * Re + 4.8_rp * sqrt(Re)

  case ( 8_ip )
     !
     ! Stokes
     !
     Re2  = min(1.0e5_rp,Re)
     CdRe = 24.0_rp

  case default
     !
     ! Others
     !
     stop

  end select

  if( Re /= 0.0_rp ) then
     Cd = CdRe/Re
  else
     Cd = 0.0_rp
  end if

end subroutine dragfo

!-----------------------------------------------------------------------
!
!   #----------------------
!   # PLOT CD USING GNUPLOT
!   #----------------------
!   #
!   reset
!   set xrange[0.01:100000]
!   set title 'Cd for different formula'
!   set log x
!   set log y
!   set xlabel 'Re'
!   set ylabel 'Cd'
!   #
!   # CHENG
!   #
!   c1(x)  = 24.0/x * ( 1.0 + 0.27 * x )**0.43 +  0.4 * ( 1 - exp(-0.04*x**0.38) )
!   #
!   # ARASTOOPOUR
!   #
!   c21(x) = 24.0/x * ( 1.0 + 0.15 * x ** 0.687 )
!   c22(x) = 0.44
!   c2(x) = ( x < 1000 ? c21(x):c22(x))
!   #
!   # GANSER
!   #
!   p  = 1.0
!   k1 = 3.0/(1.0+2.0*p**(-0.5))
!   k2 = 10.0**(1.84148*(-log10(p))**0.5743)
!   c3(x) = 24.0/(x*k1) * ( 1.0 + 0.1118*(x*k1*k2)**0.6567 ) + 0.4305*k2/( 1.0 + 3305.0/(x*k1*k2) )
!   #
!   # WILSON
!   #
!   p      = 1.0
!   c41(x) = 24.0/x*p**(-0.828)   + 2.0 *( 1.0-p )**0.5
!   cd100  = 24.0/100.0*p**(-0.828) + 2.0 *( 1.0-p )**0.5
!   c42(x) = 1.0 - (( 1-cd100 )/900.0 * (1000.0-x))
!   c43(x) = 1.0
!   c44(x) = ( x < 100  ? c41(x):c42(x))
!   c4(x)  = ( x < 1000 ? c44(x):c43(x))
!   #
!   # TURTON and LEVENSPIEL
!   #
!   c5(x) = 24.0/x * (1.0+0.173*x**0.657) + 0.413/(1.0+11630.0*x**(-1.09))
!
!   plot c1(x) t 'CHENG',c2(x) t 'ARASTOOPOUR',c3(x) t 'GANSER',c4(x) t 'WILSON',c5(x) t 'TURTON'
!
!
!   #------------------------
!   # PLOT CDRe USING GNUPLOT
!   #------------------------
!   #
!   reset
!   set xrange[0.01:100000]
!   set title 'CdRe for different formula'
!   set log x
!   set log y
!   set xlabel 'Re'
!   set ylabel 'CdRe'
!   #
!   # CHENG
!   #
!   c1(x)  = 24.0 * ( 1.0 + 0.27 * x )**0.43 + 0.4 * x * ( 1.0 - exp(-0.04*x**0.38) )
!   #
!   # ARASTOOPOUR
!   #
!   c21(x) = 24.0 * ( 1.0 + 0.15 * x ** 0.687 )
!   c22(x) = 0.44 * x
!   c2(x) = ( x < 1000 ? c21(x):c22(x))
!   #
!   # GANSER
!   #
!   p  = 1.0
!   k1 = 3.0/(1.0+2.0*p**(-0.5))
!   k2 = 10.0**(1.84148*(-log10(p))**0.5743)
!   c3(x) =  24.0/k1 * ( 1.0 + 0.1118*(x*k1*k2)**0.6567 ) + 0.4305*x*x*k2/ ( x + 3305.0/(k1*k2) )
!   #
!   # WILSON
!   #
!   p      = 1.0
!   c41(x) = 24.0*p**(-0.828) + 2.0 * x * ( 1.0-p )**0.5
!   cd100  = 24.0/(100.0)*p**(-0.828) + 2.0 * ( 1.0-p )**0.5
!   c42(x) = x * (1.0 - (( 1-cd100 )/900.0 * (1000.0-x)))
!   c43(x) = x
!   c44(x) = ( x < 100  ? c41(x):c42(x))
!   c4(x)  = ( x < 1000 ? c44(x):c43(x))
!   #
!   # TURTON and LEVENSPIEL
!   #
!   c5(x) = 24.0 * (1.0+0.173*x**0.657) + 0.413*x/(1.0+11630.0*x**(-1.09))
!
!
!   plot c1(x) t 'CHENG',c2(x) t 'ARASTOOPOUR',c3(x) t 'GANSER',c4(x) t 'WILSON',c5(x) t 'TURTON'
!
!-----------------------------------------------------------------------
