subroutine dragfo(imode,u,mu,rho,d,dista,psi,CdRe,Cd,Re,dCdRedRe)
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
  !    8. Wall correction: Zeng et al. The effects of near wall corrections to hydrodynamic forces on particle deposition and transport in vertical...    
  !          Cd ={ 1+ 0.15[1-exp(-sqrt(delta)]Re_p^[0.687+0.313*exp(-2*sqrt(delta)]}
  !          C_DO = [ 1.028 - 0.07/(1+4*delta^2)-8/15*ln(270*delta/(135+256*delta))]*(24/Re_p)
  !          delta = L/d_p - 0.5  L:distance to the wall, d_p: particle diameter
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
  use def_kintyp, only     :  ip,rp
  implicit none
  integer(ip), intent(in)  :: imode
  real(rp),    intent(in)  :: u          !< Relative velocity |u_fluid-u_particle|
  real(rp),    intent(in)  :: mu         !< Fluid viscosity
  real(rp),    intent(in)  :: rho        !< Fluid density
  real(rp),    intent(in)  :: d          !< Particle diameter
  real(rp),    intent(in)  :: dista      !< Distance to the wall     
  real(rp),    intent(in)  :: psi        !< Sphericity
  real(rp),    intent(out) :: CdRe       !< Drag * Reynolds
  real(rp),    intent(out) :: Cd         !< Drag coefficient
  real(rp),    intent(out) :: Re         !< Reynolds number
  real(rp),    intent(out) :: dCdRedRe   !< d(Cd*Re) / dRe
  real(rp)                 :: k1,k2,Re2
  real(rp)                 :: zeror,delta,CDORe
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
     CdRe = 24.0_rp * (1.0_rp+0.173_rp*Re**0.657_rp) + 0.413_rp*Re/(1.0_rp+11630.0_rp*Re**(-1.09_rp)) !!!!OJO..es la derivada cero!??

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
     ! Wall Correction: Zeng
     !
     delta    = dista/d-0.5_rp
     CDORe    = ( 1.028_rp - 0.07_rp/(1.0_rp+4.0_rp*delta**2)- 8.0_rp/15.0_rp * log( 270.0_rp*delta/(135.0_rp+256.0_rp*delta) ) )*24.0_rp
     CdRe     = ( 1.0_rp + 0.15_rp*(1.0-exp(-sqrt(delta)))*Re**(0.687_rp+0.313_rp*exp(-2.0_rp*sqrt(delta))))*CDORe
     dCdRedRe = ( 1.0_rp + 0.15_rp*(1.0-exp(-sqrt(delta)))*(0.687_rp+0.313_rp*exp(-2.0_rp*sqrt(delta))*&
                 Re**(0.687_rp+0.313_rp*exp(-2.0_rp*sqrt(delta))-1.0_rp)) )*CDORe
  case default
     !
     ! Others
     !
     call runend('DRAGFO: NON-EXISTING MODEL')

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
!
!   Visulization of the function phi of the Jacobian (only drag):
!   -------------------------------------------------------------
!
!                                  Dui Duj 
!   Jij = -(1+k) * dij - k * phi * -------
!                                  ||Du||^2
!
!   with k      = gamma * alphaD * Cd * Rep * dt
!        alphaD = 3/4 muf/(rhop*dp^2*Cslip)
!        phi    = 1/Cd * d(RepCd)/dRe
!        det(J) = -(-1-k)^{nd-1} * (1+k+k*phi), nd=problem dimension
!
!   reset
!   psi         = 1.0
!   k1          = 3.0 / ( 1.0 + 2.0/sqrt(psi) )
!   k2          = 10.0**( 1.84148 * (-log10(psi))**0.5743 )
!   CdRe(x)     = 24.0/k1 * ( 1.0 + 0.1118*(x*k1*k2)**0.6567 ) + 0.4305*x*x*k2/ ( x + 3305.0/(k1*k2) )
!   Cd(x)       = 1.0/x*(24.0/k1 * ( 1.0 + 0.1118*(x*k1*k2)**0.6567 ) + 0.4305*x*x*k2/ ( x + 3305.0/(k1*k2) ))
!   dCdRedRe(x) = 24.0/k1 * ( 0.1118*0.6567*(k1*k2)**0.6567*x**(1.0-0.6567) ) \
!                 + 0.4305 * x * 2.0 * k2 / ( x + 3305.0/(k1*k2) ) \
!                 - 0.4305 * x * x * k2 / ( x + 3305.0/(k1*k2) ) ** 2
!   phi(x)      = dCdRedRe(x)/Cd(x)
!   set xrange[1.0e-3:1e6]
!   set log x
!   set log y
!   plot phi(x)
!
!-----------------------------------------------------------------------
