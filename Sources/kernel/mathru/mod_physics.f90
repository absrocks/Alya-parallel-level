!-----------------------------------------------------------------------
!> @addtogroup Physics
!> @{
!> @file    mod_physics.f90
!> @author  houzeaux
!> @date    2018-09-22
!> @brief   Physics
!> @details Subroutines for physical laws
!>
!-----------------------------------------------------------------------

module mod_physics

  use def_kintyp,   only : ip,rp
  use def_parame,   only : pi,oneo3
  use mod_messages, only : messages_live 
  implicit none
  
  !---------------------------------
  ! Structure for liquid parameters 
  !---------------------------------
  type liquid_params
      !
      ! New method for Daubert & Danner (1985) Data compilation tables of properties of pure compounds
      !
      real(rp)      :: rho(6)
      real(rp)      :: cp(6)
      real(rp)      :: Lv(6)
      real(rp)      :: psat(6)
      integer(ip)   :: irho
      integer(ip)   :: icp
      integer(ip)   :: iLv
      integer(ip)   :: ipsat
      
  end type liquid_params

  !----------------------------
  ! Structure for liquid state 
  !----------------------------
  type liquid_state
      character(5)         :: name             ! Name of liquid 
      integer(ip)          :: ID               ! ID of liquid
      real(rp)             :: W                ! Molar mass 
      real(rp)             :: T                ! Temperature 
      real(rp)             :: Pc               ! Critical pressure
      real(rp)             :: Tc               ! Critical temperature 
      real(rp)             :: P                ! Pressure
      real(rp)             :: rho              ! Density
      real(rp)             :: cp               ! Specific heat
      real(rp)             :: Lv               ! Heat of vaporization
      real(rp)             :: psat             ! Saturation pressure
      real(rp)             :: Tsat             ! Saturation temperature

      type(liquid_params)  :: param            ! Parameters 
  end type liquid_state

  real(rp), parameter :: zeror                  = epsilon(1.0_rp)
  real(rp), parameter :: universal_gas_constant = 8.13144598_rp
  real(rp), parameter :: k_boltzmann            = 1.3806503e-23_rp

  private 

  public :: universal_gas_constant
  
  public :: physics_drag_force
  public :: physics_lift_force
  public :: physics_sphere_diameter
  public :: physics_sphere_mass
  public :: physics_sphere_Nusselt
  public :: physics_sphere_Sherwood
  public :: physics_air_water_vapor_diffusion_coefficient
  public :: physics_Cunningham
  public :: physics_humidity
  public :: liquid_state
  public :: physics_initialize_liquid
  public :: physics_set_liquid_temperature
  public :: physics_set_liquid_pressure
  public :: physics_mass_2_w
  public :: physics_mole_2_mass
  public :: physics_T_2_HCp
  public :: physics_H_2_TCp
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-09-22
  !> @brief   Drag force
  !> @details Drag force correlations
  !>          By definition:
  !>          Fd = 3/4 * rho_a * u^2 * V/d * Cd
  !>          By using V = (4/3)*pi*(d/2)^3, we end up with:
  !>      
  !>          Fd = rho_a * u^2 * pi * d^2/8 * Cd
  !>             = ( 1/2 * rho_a * u^2 ) * ( pi * (d/2)^2 ) * Cd
  !>      
  !>          To be used for a drag model of the kind:
  !>      
  !>          Fd(i) = - pi * mu * d / 8.0_rp * Cd * Re * ( u(i)-uf(i) )
  !>          ad(i) =   Fd(i) / (rho_p*pi*d^3*)/6)
  !>                 =  visfl / rho_p * 0.75 * Cd * Re / d^2 * ( u(i)-uf(i) )
  !>      
  !>          or, equivalently:
  !>      
  !>          Fd = pi * mu * d/8 * (Cd*Re) * u
  !>      
  !>          The associated terminal velocity for a sphere is:
  !>          ut = sqrt( 4*g*d/(3*Cd) * (rho_p-rho_a)/rho_a )
  !>      
  !>          CD FORMULA:
  !>          -----------
  !>      
  !>          1. CHENG model.
  !>             Nian-Sheng Cheng, 'Comparison of formulas for drag coefficient
  !>             and settling velocity of spherical particles',
  !>             Powder Technology, 2008.
  !>      
  !>             Cd = 24/Re * ( 1 + 0.27 * Re ) ** 0.43 + 0.4 * ( 1 - exp(-0.04*Re**0.38) )
  !>      
  !>          2. GANSER model.
  !>             Ganser, G., 1993. A rational approach to drag prediction of spherical and nonspherical
  !>             particles. Powder Technology 77, 143-152.
  !>      
  !>             Cd = 24/(Re*k1) * ( 1 + 0.1118*(Re*k1*k2)^0.6567 ) + 0.4305*k2/ [ 1 + 3305/(Re*k1*k2) ]
  !>             k1 = 3/(1+2*p^-0.5)
  !>             k2 = 10^(1.84148*(-log10(p))^0.5743)
  !>             p  = sphericity (=1 for a sphere)
  !>      
  !>          3. ARASTOOPOUR model.
  !>             Arastoopour, H., Wang, C., Weil, S., 1982. Particle-particle interaction
  !>             force in a dilute gas-solid system. Chemical Engineering Science 37 (9), 1379-1386.
  !>             AB: I don't know why not call it Schiller and Neumann model (1933)
  !>      
  !>             if( Re <= 1000 ) then
  !>                 Cd = 24.0/Re * ( 1 + 0.15 * Re ** 0.687 )
  !>             else
  !>                 Cd = 0.44
  !>             end if
  !>      
  !>          4. WILSON model.
  !>             Wilson, L., Huang, T., 1979. The influence of shape on the atmospheric settling velocity
  !>             of volcanic ash particles. Earth and Planetary Science Letters 44, 311-324.
  !>      
  !>             if( Re <= 100 ) then
  !>                 Cd = 24/Re*p^-0.828 + 2 sqrt( 1.0-p )
  !>             else if( 100 <= Re <= 1000 ) then
  !>                 Cd = 1 - ( 1-Cd|Re=10^2 )/900 * (10^3-Re)
  !>             else
  !>                 Cd = 1
  !>             end if
  !>             p = (b+c)/2a = particle aspect ratio (a>b>c are the particle semi-axes); p=1 for sphere
  !>      
  !>          5. TURTON and LEVENSPIEL model.
  !>             Turton, R., and O. Levenspiel. 1986. A short note on drag correlation for spheres. Powder Technolo-
  !>             gy Journal, 47, 83.
  !>             Very similar to Cheng's. Valid for Re < 2.6x10^5
  !>      
  !>          8. Wall correction: Zeng et al. The effects of near wall corrections to hydrodynamic
  !>             forces on particle deposition and transport in vertical...    
  !>             Cd ={ 1+ 0.15[1-exp(-sqrt(delta)]Re_p^[0.687+0.313*exp(-2*sqrt(delta)]}
  !>             C_DO = [ 1.028 - 0.07/(1+4*delta^2)-8/15*ln(270*delta/(135+256*delta))]*(24/Re_p)
  !>             delta = L/d_p - 0.5  L:distance to the wall, d_p: particle diameter
  !>      
  !>      
  !>                                       Cd for different formula
  !>          Cd
  !>            10000 ++-+---+-+--+---+-++--+--+-++--+---+-++--+--+-++--+---+-++-+---+-++
  !>                  +        +         +        +         +        +    GANSER ****** +
  !>                  **                                                  WILSON ###### +
  !>             1000 +***                                                             ++
  !>                  +   ****                                                          +
  !>                  |      ***                                                        |
  !>                  +         ***                                                     +
  !>              100 ++           ***                                                 ++
  !>                  +              ****                                               +
  !>                  +                 #***                                            +
  !>               10 ++                   ****                                        ++
  !>                  +                      ##***                                      +
  !>                  |                          #****                                  |
  !>                  +                            ##*****                              +
  !>                1 ++                              ### ******  #######################
  !>                  +                                  ### ##**************************
  !>                  +        +         +        +        ###       +         +        +
  !>              0.1 ++-+---+-+--+---+-++--+--+-++--+---+-++--+--+-++--+---+-++-+---+-++
  !>                 0.01     0.1        1        10       100      1000     10000    100000
  !>                                                  Re
  !-----------------------------------------------------------------------

  subroutine physics_drag_force(imode,u,mu,rho,d,Cd,Re,CdRe,dCdRedRe,SPHERICITY,DISTANCE)

    integer(ip), intent(in)            :: imode               !< Drag model
    real(rp),    intent(in)            :: u                   !< Relative velocity |u_fluid-u_particle|
    real(rp),    intent(in)            :: mu                  !< Fluid viscosity
    real(rp),    intent(in)            :: rho                 !< Fluid density
    real(rp),    intent(in)            :: d                   !< Particle diameter
    real(rp),    intent(out)           :: Cd                  !< Drag coefficient
    real(rp),    intent(out)           :: Re                  !< Reynolds number
    real(rp),    intent(out)           :: CdRe                !< Drag * Reynolds
    real(rp),    intent(out)           :: dCdRedRe            !< d(Cd*Re) / dRe
    real(rp),    intent(in),  optional :: DISTANCE            !< Distance to the wall     
    real(rp),    intent(in),  optional :: SPHERICITY          !< Sphericity
    real(rp)                           :: k1,k2,Re2,psi,dista
    real(rp)                           :: zeror,delta,CDORe
    !
    ! Optional
    !
    if( present(SPHERICITY) ) then
       psi = SPHERICITY
    else
       psi = 1.0_rp
    end if
    if( present(DISTANCE) ) then
       dista = DISTANCE
    else
       dista = 1.0e6_rp
    end if
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
       ! Arastoopour model  A.K.A.  Schiller & Neumann model
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
       if( .not. present(DISTANCE) ) call runend('PHYSICS_DRAG_FORCE: DISTANCE REQUIRED')
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

  end subroutine physics_drag_force

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

  function physics_funReZeng(L,Re,d) result (x)
    
    real(rp),  intent(in)  :: L
    real(rp),  intent(in)  :: Re
    real(rp),  intent(in)  :: d
    real(rp)               :: x
    real(rp)               :: f0, f1, Clt0, L2

    f0 = 1.0_rp + 0.329_rp*Re + 0.00485_rp*Re*Re
    f1 = -0.9_rp * tanh( 0.022_rp*Re )

    L2 = (L*Re)/d 

    if( L2 < 10_rp ) then

       Clt0 = ( 9.0_rp/8.0_rp + 5.78e-6_rp * L2 ) * exp( -0.292_rp*L2 ) 

    else

       Clt0 = 8.94_rp*L2**(-2.09_rp)

    end if

    x = f0 * Clt0 * L**f1

  end function physics_funReZeng

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-09-22
  !> @brief   Lift force
  !> @details Lift coefficent for a particle immersed  into a viscous fluid
  !>          See:
  !>          Zeng, L., Najjar, F., Balachandar, S., Fischer, P., 2009. Forces
  !>          on a finite-sized particle located close to a wall in a linear
  !>          shear flow. Phys. Fluids 21, 033302.
  !> 
  !-----------------------------------------------------------------------

  subroutine physics_lift_force(imode,ur,wf,mu,rho,d,dista,Cls,Re,Res)

    integer(ip), intent(in)  :: imode       !< Model
    real(rp),    intent(in)  :: ur          !< Relative velocity |u_fluid-u_particle|
    real(rp),    intent(in)  :: wf          !< Vorticity of velocity 
    real(rp),    intent(in)  :: mu          !< Fluid viscosity
    real(rp),    intent(in)  :: rho         !< Fluid density
    real(rp),    intent(in)  :: d           !< Particle diameter
    real(rp),    intent(in)  :: dista       !< Distance to the wall     
    real(rp),    intent(out) :: Cls         !< Lift coefficient
    real(rp),    intent(out) :: Re          !< Reynolds number
    real(rp),    intent(out) :: Res         !< Shear Reynolds number
    real(rp)                 :: alpha,beta
    real(rp)                 :: lambda,delta
    real(rp)                 :: Clsw,funcRe
    real(rp)                 :: zeror,g
    !
    ! Particle Reynolds number
    !
    Re       =  rho * ur * d / mu + zeror
    Res      =  rho * d*d * wf / mu

    select case ( imode )
       
    case ( 1_ip )
       !
       ! Saffman Mei
       !
       if (Re == 0.0_rp) then
          beta   =  0.0_rp
       else
          beta   =  0.5_rp * Res/Re
       end if

       if(Re > 40.0_rp) then
          funcRe = 0.0524_rp * sqrt(beta*Re)
       else
          funcRe = (1.0_rp-0.3314_rp * sqrt(beta)) * exp(-Re/10.0_rp) + 0.3314_rp * sqrt(beta)
       end if

       if( Res /= 0.0_rp ) then
          Cls = 4.1126_rp * funcRe/sqrt(Res)
       else
          Cls = 0.0_rp
       end if
       
    case ( 2_ip )
       !
       ! Zeng, L., Najjar, F., Balachandar, S., Fischer, P., 2009. Forces on a finite-sized particle 
       ! located close to a wall in a linear shear flow. Phys. Fluids 21, 033302.
       !
       ! Is applicable in circumstances when a stationary particle is positioned in a wall-bounded linear shear flow for 1 < Res < 200 and
       ! even when the particle touches the wall. 
       !
       delta  = dista/d-0.5_rp      
       Clsw   = 3.663_rp / (Res*Res + 0.1173_rp)**0.22_rp 
       alpha  = - exp ( -0.3_rp + 0.025_rp*Res) 
       beta   = 0.8_rp + 0.01_rp * Res
       lambda = (1.0_rp - exp(-delta))**(2.5_rp)
       Cls    = Clsw * exp( -0.5_rp*delta  * (Res/250_rp)**(4.0_rp/3.0_rp) ) * (exp( alpha * delta**beta) - lambda) 

    case ( 3_ip )
       !
       ! Zeng, L., Najjar, F., Balachandar, S., Fischer, P., 2009. Forces on a finite-sized particle 
       ! located close to a wall in a linear shear flow. Phys. Fluids 21, 033302.
       !
       ! Is applicable with translation parallel to the nearby plane wall induced lift force, for 0 < Re < 100 and 0 < L* < 300. 
       !
       ! L* = LRe/d
       !
       ! where L is s the distance from the center of the particle to the nearby wall.
       !
       delta  = dista/d-0.5_rp      

       Clsw   = 0.313_rp + 0.812_rp*exp( -0.125_rp*Re**0.77_rp ) 
       g      = 3.0_rp * exp( -0.17_rp*Re**0.7_rp )
       Cls    = physics_funReZeng(dista,Re,d) + (Clsw - physics_funReZeng(0.5_rp,Re,d))*exp( -11.0_rp*(delta/g)**1.2_rp )

    case default
       !
       ! Others
       !
       call runend('LIFTFO: NON-EXISTING MODEL')

    end select


  end subroutine physics_lift_force

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-03-18
  !> @brief   Give the diameter of a sphere 
  !> @details Give the diameter of a sphere given its mass and density
  !> 
  !-----------------------------------------------------------------------

  real(rp) function physics_sphere_diameter(mass,density)
    
    real(rp), intent(in) :: mass
    real(rp), intent(in) :: density
    real(rp)             :: term1, term2

    !
    ! To avoid precision problems:
    !
    term1= ( 6.0_rp * mass ) ** oneo3
    term2= ( density * pi  ) ** oneo3
    
    physics_sphere_diameter = term1 / term2
    
  end function physics_sphere_diameter

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-03-18
  !> @brief   Give the mass of a sphere 
  !> @details Give the mass of a sphere given its diameter and density
  !> 
  !-----------------------------------------------------------------------

  real(rp) function physics_sphere_mass(diameter,density)
    
    real(rp), intent(in) :: diameter !< Mass
    real(rp), intent(in) :: density  !< Density
    
    physics_sphere_mass = density * diameter ** 3 * pi / 6.0_rp 
    
  end function physics_sphere_mass
  
  
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-03-18
  !> @brief   Sphere Nusselt
  !> @details Correlations for Nusselt number of a sphere
  !>          - Model 1:
  !>            F. P. Incropera and D. P. De Witt. Fundamentals of Heat and Mass
  !>            Transfer. John Wiley and Sons, New York, 4th edition, 1996. 69, 82
  !> 
  !-----------------------------------------------------------------------
  
  real(rp) function physics_sphere_Nusselt(model,Re,k,Cp,mu,Prandtl)

    integer(ip),           intent(in) :: model    !< Specific model
    real(rp),              intent(in) :: Re       !< Particle Reynolds number 
    real(rp),              intent(in) :: k        !< Thermal conductivity of medium [ W/(m.K) ]
    real(rp),    optional, intent(in) :: Cp       !< Specific heat                  [ J/(kg.K) ]
    real(rp),    optional, intent(in) :: mu       !< Viscosity                      [ kg/(m.s) ]
    real(rp),    optional, intent(in) :: Prandtl  !< Prandtl number
    real(rp)                          :: Nu,Pr
    !
    ! Prandtl number
    !
    if( present(Prandtl) ) then
       Pr = Prandtl
    else if( present(cp) .and. present(mu) ) then 
       Pr = cp * mu / k
    else
       call runend('PARAMETERS ARE MISSING TO COMPUTE HEAT TRANSFER COEFFICIENT')
    end if

    select case ( model )

    case ( 1_ip )
       !
       ! Incropera and Witt  = Ranz & Mashsall model
       !
       physics_sphere_Nusselt = 2.0_rp + 0.6_rp * sqrt(Re) * Pr ** oneo3 ! Nusselt number
       
    case ( 2_ip )
       !
       ! Bear and Pruppacher model
       !
       if( sqrt(Re) * Pr ** oneo3 < 1.4_rp ) then
          physics_sphere_Nusselt = 2.0_rp + 0.2160_rp * ( sqrt(Re) * Pr ** oneo3 )
       else
          physics_sphere_Nusselt = 1.560_rp + 0.6160_rp *( sqrt(Re) * Pr ** oneo3 )
       end if
       
    case default
       !
       ! Others
       !
       call runend('NON-EXISTING HEAT TRANSFER COEFFICIENT MODEL')

    end select
    
  end function physics_sphere_Nusselt
  

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-03-18
  !> @brief   Mass transfer coefficient
  !> @details Correlations for mass transfer coefficient [m/s]
  !>          - Model 1:
  !>            F. P. Incropera and D. P. De Witt. Fundamentals of Heat and Mass
  !>            Transfer. John Wiley and Sons, New York, 4th edition, 1996. 69, 82
  !> 
  !-----------------------------------------------------------------------
  
  real(rp) function physics_sphere_Sherwood(model,Re,D,mu,rho,Schmidt)

    integer(ip),           intent(in) :: model    !< Specific model
    real(rp),              intent(in) :: Re       !< Particle Reynolds number 
    real(rp),              intent(in) :: D        !< Binary diffusion coef between particle and medium [ m^2/s ]
    real(rp),    optional, intent(in) :: mu       !< Viscosity of medium                               [ kg/(m^2.s) ]
    real(rp),    optional, intent(in) :: rho      !< Density of medium                                 [ kg/m^3 ]
    real(rp),    optional, intent(in) :: Schmidt  !< Schmidt number
    real(rp)                          :: Sc,Sh
    !
    ! Schmidt number
    !
    if( present(Schmidt) ) then
       Sc = Schmidt
    else if( present(mu) .and. present(rho) ) then
       Sc = mu / ( rho * D )
    else
       call runend('PARAMETERS ARE MISSING TO COMPUTE MASS TRANSFER COEFFICIENT')
    end if
    
    select case ( model )

    case ( 1_ip )
       !
       ! Incropera and Witt = Ranz & Mashsall model
       !
       physics_sphere_Sherwood = 2.0_rp + 0.6_rp * sqrt(Re) * Sc ** oneo3  
       
    case ( 2_ip )
       !
       ! Bear and Pruppacher model
       !
       if( sqrt(Re) * Sc ** oneo3 < 1.4_rp ) then
          physics_sphere_Sherwood = 2.0_rp + 0.2160_rp * ( sqrt(Re) * Sc ** oneo3 )
       else
          physics_sphere_Sherwood = 1.560_rp + 0.6160_rp * ( sqrt(Re) * Sc ** oneo3 )
       end if

     case default
       !
       ! Others
       !
       call runend('NON-EXISTING MASS TRANSFER COEFFICIENT MODEL')

    end select
    
  end function physics_sphere_Sherwood

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-03-18
  !> @brief   Air-water vapor diffusion coefficient
  !> @details Air-water vapor diffusion coefficient in [m^2/s]
  !> 
  !-----------------------------------------------------------------------

  real(rp) function physics_air_water_vapor_diffusion_coefficient(model,T)

    integer(ip), intent(in) :: model    !< Specific model
    real(rp),    intent(in) :: T        !< Temperature    [K]

    select case ( model )

    case ( 1_ip )
       !
       ! Bolz and Tuve (1976)
       !
       physics_air_water_vapor_diffusion_coefficient = &
            - 2.775e-6_rp                  &
            + 4.479e-8_rp  * T             &
            + 1.656e-10_rp * T * T

    case ( 2_ip )
       !
       ! Approximate fit (valid up to a temperature of 40C)
       !
       physics_air_water_vapor_diffusion_coefficient = &
            22.5e-6_rp * (T/273.15_rp)**1.8_rp
       
    case default
       !
       ! Others
       !
       call runend('NON-EXISTING BINARY DIFFUSION COEFFICIENT MODEL')

    end select

  end function physics_air_water_vapor_diffusion_coefficient

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-03-18
  !> @brief   liquid equilibrium vapour mass fraction
  !> @details liquid equilibrium vapour mass fraction, the liquid is the
  !>          droplet
  !>                                            xl                  
  !>           Vapor fraction yl = -----------------------------    
  !>                                  +-    w_fluid -+   w_fluid    
  !>                               xl | 1 - -------  | + -------
  !>                                  +-       w    -+      w
  !>
  !>                    +-    H_vapor w    +- 1    1 -+ -+
  !>           xl = exp | ---------------- |  -  - -  |  |
  !>                    +-      R          +- T_b  T -+ -+
  !>
  !-----------------------------------------------------------------------

!! AB not used   !MMG! !Change name to: liquid equilibrium vapor mass fraction
!! AB not used   real(rp) function physics_vapor_mass_fraction(H_vapor,w,R,T_boiling,T,w_fluid)
!! AB not used 
!! AB not used     real(rp), intent(in) :: H_vapor   !< Latent heat of vaporisation of the liquid [ J/kg ]
!! AB not used     real(rp), intent(in) :: w         !< Molecular weight of the liquid            [ kg/mol ]
!! AB not used     real(rp), intent(in) :: R         !< Universal gas constant                    [ J / (mol.K) ]
!! AB not used     real(rp), intent(in) :: T_boiling !< Boiling temperature of the liquid         [ K ]
!! AB not used     real(rp), intent(in) :: T         !< Temperature of the droplet                [ K ]
!! AB not used     real(rp), intent(in) :: w_fluid   !< Molecular weight of the fluid             [ kg/mol ]
!! AB not used     real(rp)             :: xl
!! AB not used     
!! AB not used     xl = exp( H_vapor * w / R * ( 1.0_rp / T_boiling - 1.0_rp / T ) )
!! AB not used     physics_vapor_mass_fraction = xl / ( xl*(1.0_rp-w_fluid/w) + w_fluid/w )
!! AB not used 
!! AB not used   end function physics_vapor_mass_fraction


  !
  ! Mass fraction from mole fraction and molar masses
  !
  real(rp) function physics_mole_2_mass(x,wx,wo)
    real(rp), intent(in) :: x
    real(rp), intent(in) :: wx
    real(rp), intent(in) :: wo

    !             Wk * Xk
    ! Yk = --------------------------
    !       Wk * Xk + Wother * (1-Xk)
    !       
    physics_mole_2_mass = x / ( x + (1.0_rp -x) * wo / wx )

  end function physics_mole_2_mass

  !
  ! Mean molar mass from mass fractions
  !
  real(rp) function physics_mass_2_w(y,wy,wo)
    real(rp), intent(in) :: y
    real(rp), intent(in) :: wy
    real(rp), intent(in) :: wo

    !               1
    !  W = -----------------------   
    !         Yk/Wk + (1-Yk)/Wo
    !   
    physics_mass_2_w = 1.0_rp / ( y/wy + (1.0_rp -y) / wo )

  end function physics_mass_2_w


  !
  ! Calculate enthalpy and specific heat at a given temperature
  ! Using the NASA polynomials
  !
  subroutine physics_T_2_HCp(T,cpcoef,h,cp)
    real(rp), intent(in) :: T
    real(rp), intent(in) :: cpcoef(6,2)
    real(rp), intent(out):: h
    real(rp), intent(out):: cp

    real(rp)    :: Tclip
    real(rp)    :: cpc(6)
    integer(ip) :: ii

    Tclip =  min(3000.0_rp, max(200.0_rp, T))
    !
    ! First is low temperature, second is high temperature in table
    !
    if (Tclip < 1000.0_rp) then
        cpc = cpcoef(:,1)
    else
        cpc = cpcoef(:,2)
    endif

    cp = cpc(5)
    h  = cpc(5)/5.0_rp

    do ii = 4,1,-1
        cp = cp * Tclip
        h  = h  * Tclip
        cp = cp + cpc(ii)
        h  = h  + cpc(ii)/real(ii, KIND=rp)
    enddo
    h = h * Tclip
    h = h + cpc(6)

  end subroutine physics_T_2_HCp


  !
  ! Calculate temperature and specific heat at a given enthalpy
  ! Using the NASA polynomials
  !
  subroutine physics_H_2_TCp(h,cpcoef,T,cp)
    real(rp), intent(in)    :: h
    real(rp), intent(in)    :: cpcoef(6,2)
    real(rp), intent(inout) :: T
    real(rp), intent(out)   :: cp

    real(rp)    :: error, h_guess, cp_guess
    integer(ip) :: ii

    T       =  min(3000.0_rp, max(200.0_rp, T))
    error   = 1.0e10_rp
    ii      = 0_ip

    do while( (abs(error) > 1.0e-6_rp) .and. (ii<100_ip) ) 
        call physics_T_2_HCp(T,cpcoef,h_guess,cp_guess)
        error   = h_guess - h
        T       = T - error / cp_guess
        ii      = ii + 1_ip
    enddo

    !
    ! Do a last clipping here
    !
    T       =  min(3000.0_rp, max(200.0_rp, T))
    call physics_T_2_HCp(T,cpcoef,h_guess,cp_guess)
    cp = cp_guess
  
  end subroutine physics_H_2_TCp



  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-03-18
  !> @brief   Cunningham slip correction factor
  !> @details Compute the Cunningham slip correction factor to be applied
  !>          to the drag:
  !>
  !>          Cunningham slip correction factor
  !>          Some values of airborne particles at standrad conditions (T=293K,P=101kPa)
  !>
  !>          Particle diameter   Slip Correction
  !>          d(m)                factor Cc
  !>          -----------------------------------
  !>          10^-9               224.332
  !>          10^-8                22.976
  !>          10^-7                 2.928
  !>          10^-6                 1.155
  !>          10^-5                 1.015
  !>          10^-4                 1.002
  !>
  !-----------------------------------------------------------------------

  real(rp) function physics_Cunningham(mean_free_path,diameter)

    real(rp), intent(in) :: mean_free_path !< Mean free path 
    real(rp), intent(in) :: diameter       !< Particle diameter
    
    if( mean_free_path > 0.0_rp ) then
       physics_Cunningham = 1.0_rp+2.0_rp*mean_free_path/diameter &
            *(1.257_rp+0.4_rp*exp(-0.55_rp*diameter/mean_free_path))
    else
       physics_Cunningham = 1.0_rp       
    end if
    
  end function physics_Cunningham

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-03-19
  !> @brief   Real humidity
  !> @details Compute real humidity as a function of temperature and
  !>          relative humidity
  !> 
  !-----------------------------------------------------------------------

  !MMG! !Change name to: vapor mass fraction.
  real(rp) function physics_humidity(T,RelHum)

    real(rp), intent(in) :: T
    real(rp), intent(in) :: Relhum
    !
    !Approximation of absolute humidity (of water) of air at that temperature
    !
    physics_humidity = -37.562_rp + 0.7491_rp*T + (-.0059263)*T**2 &
         + 0.00002331_rp*T**3 + (-4.5683e-8_rp)*T**4 + 3.5756e-11_rp*T**5
    !
    ! Computation of real humidity based on absolute humidity and relative humidity
    !
    physics_humidity = RelHum * physics_humidity

  end function physics_humidity


  
    !
    ! AB: function to initialize liquid state 
    ! 
    subroutine physics_initialize_liquid(liq,T,P, Lv, T_boiling, cp, rho, W)
        type(liquid_state),   intent(inout) :: liq
        real(rp),                intent(in) :: T
        real(rp),                intent(in) :: P

        real(rp), optional,      intent(in) :: Lv 
        real(rp), optional,      intent(in) :: T_boiling
        real(rp), optional,      intent(in) :: cp
        real(rp), optional,      intent(in) :: rho
        real(rp), optional,      intent(in) :: W
         
        liq % T            = T
        liq % P            = P
        liq % param % rho  = 0.0_rp
        liq % param % cp   = 0.0_rp
        liq % param % Lv   = 0.0_rp
        liq % param % psat = 0.0_rp


        select case (liq % name)
            
            case ('USER ') 
               if (present( W )) then
                   liq % W     = W
               else
                   call runend('mod_physics: W not given.')
               endif

               if (present( rho )) then
                  liq % param % irho    = 100_ip
                  liq % param % rho(1)  = rho
               else
                   call runend('mod_physics: rho not given.')
               endif

               if (present( cp )) then
                  liq % param % icp     = 100_ip
                  liq % param % cp(1)   = cp 
               else
                   call runend('mod_physics: cp not given.')
               endif

               liq % Tc    = 99999.0_rp ! K, critical 
               
               liq % param % iLv     = 100_ip
               liq % param % Lv(1)   = Lv 

               !
               ! Clausius-Clapeyron
               !
               liq % param % ipsat   = 1_ip
               liq % param % psat(1) = P
               liq % param % psat(2) = T_boiling
               liq % param % psat(3) = Lv
               liq % param % psat(4) = W


            
            case ('NHEPT')
               liq % W  = 0.10020592_rp 
               liq % Tc = 540.2_rp ! K, critical
               
               liq % param % irho    = 105_ip
               liq % param % rho(1)  = 61.38396836_rp
               liq % param % rho(2)  = 0.26211_rp
               liq % param % rho(3)  = liq % Tc
               liq % param % rho(4)  = 0.28141_rp
               
               liq % param % icp     = 114_ip
               liq % param % cp(1)   = 6.11976102401216_rp
               liq % param % cp(2)   = 3137.69909384855_rp
               liq % param % cp(3)   = 182.274175063868_rp
               liq % param % cp(4)   = -254.530511150515_rp
               
               liq % param % iLv     = 106_ip
               liq % param % Lv(1)   = 499121.791545248_rp
               liq % param % Lv(2)   = 0.38795_rp

               liq % param % ipsat   = 101_ip
               liq % param % psat(1) = 87.829_rp
               liq % param % psat(2) = -6996.4_rp
               liq % param % psat(3) = -9.8802_rp
               liq % param % psat(4) = 7.2099e-6_rp
               liq % param % psat(5) = 2.0_rp 

            case ('NDODE')
               liq % W  = 0.170338_rp 
               liq % Tc = 658.0_rp ! K, critical

               liq % param % irho    = 105_ip
               liq % param % rho(1)  = 60.53982858_rp 
               liq % param % rho(2)  = 0.25511_rp
               liq % param % rho(3)  = liq % Tc
               liq % param % rho(4)  = 0.29368_rp
               
               liq % param % icp     = 100_ip
               liq % param % cp(1)   = 2983.53861146661_rp
               liq % param % cp(2)   = -8.0352006011577_rp
               liq % param % cp(3)   = 0.018207916025784_rp
               
               liq % param % iLv     = 106_ip
               liq % param % Lv(1)   = 454020.829174935_rp 
               liq % param % Lv(2)   = 0.40681_rp
              
               liq % param % ipsat   = 101_ip
               liq % param % psat(1) = 137.47_rp
               liq % param % psat(2) = -11976.0_rp
               liq % param % psat(3) = -16.698_rp
               liq % param % psat(4) = 8.0906e-06_rp
               liq % param % psat(5) = 2.0_rp
                
            case ('H2O  ')
               liq % W  = 0.018015_rp 
               liq % Tc = 647.13_rp ! K, critical

               liq % param % irho    = 105_ip
               liq % param % rho(1)  = 98.343885_rp 
               liq % param % rho(2)  = 0.30542_rp
               liq % param % rho(3)  = liq % Tc
               liq % param % rho(4)  = 0.081_rp
               
               liq % param % icp     = 100_ip
               liq % param % cp(1)   = 15341.1046350264_rp
               liq % param % cp(2)   = -116.019983347211_rp
               liq % param % cp(3)   = 0.451013044684985_rp
               liq % param % cp(4)   = -0.000783569247849015_rp
               liq % param % cp(5)   = 5.20127671384957e-07_rp
               
               liq % param % iLv     = 106_ip
               liq % param % Lv(1)   = 2889425.47876769_rp 
               liq % param % Lv(2)   = 0.3199_rp
               liq % param % Lv(3)   = -0.212_rp
               liq % param % Lv(4)   = 0.25795_rp
              
               liq % param % ipsat   = 101_ip
               liq % param % psat(1) = 73.649_rp
               liq % param % psat(2) = -7258.2_rp
               liq % param % psat(3) = -7.3037_rp
               liq % param % psat(4) = 4.1653e-06_rp
               liq % param % psat(5) = 2.0_rp
                
            case default
                call runend('mod_physics: Liquid name is not recognized.')
            
        end select

        !
        ! Get critical pressure
        !
        call physics_set_liquid_temperature(liq, liq % Tc)
        liq % Pc = liq % psat
        call physics_set_liquid_temperature(liq,T)


        call physics_set_liquid_pressure(liq,P)
        call physics_set_liquid_temperature(liq,T)

    end subroutine physics_initialize_liquid





    !
    ! AB: function to calculate liquid properties 
    ! 
    subroutine physics_set_liquid_temperature(liq,T)
        type(liquid_state),      intent(inout) :: liq
        real(rp),                intent(in) :: T
        !
        ! Heptane coeffs
        !
        real(rp) :: Tnorm, eta, dummr
        real(rp) :: rho_a, rho_b, rho_c
        real(rp) :: cp_a, cp_b, cp_c, cp_d
        real(rp) :: Lv_a, Lv_b
        real(rp) :: psat_a, psat_b, psat_c, psat_d

        liq % T = max(200.0_rp,  T)

        ! 
        ! Density
        !
        liq % rho = physics_indexed_property_functions(liq % param % irho, &
                                                       liq % T,            &
                                                       liq % param % rho,  &
                                                       liq % Tc)
        !
        ! Specific heat
        !
        liq % cp  = physics_indexed_property_functions(liq % param % icp,  &
                                                       liq % T,            &
                                                       liq % param % cp,   &
                                                       liq % Tc)
        !
        ! Heat of vaporization
        !
        liq % Lv  = physics_indexed_property_functions(liq % param % iLv,  &
                                                       liq % T,            &
                                                       liq % param % Lv,   &
                                                       liq % Tc)
        !
        ! Saturation pressure
        !
        liq % psat= physics_indexed_property_functions(liq % param % ipsat,&
                                                       liq % T,            &
                                                       liq % param % psat, &
                                                       liq % Tc)


    end subroutine physics_set_liquid_temperature


    function physics_indexed_property_functions(itask,x,param,xcrit) result (y)
        integer(ip),          intent(in)    :: itask
        real(rp),             intent(in)    :: x
        real(rp),             intent(in)    :: param(:)
        real(rp),   intent(in), optional    :: xcrit
        real(rp)                            :: y

        !real(rp)                            :: xnorm, eta, xbound

        if (itask >= 100_ip) then
            !
            ! Daubert & Danner (1985) Data compilation tables of properties of pure compounds
            !
            y = physics_NSRDS_functions(itask,x,param,xcrit)
        else
            select case(itask)
                case(1)
                    !
                    ! Clausius-Clapeyron equation for saturation pressure of fluids 
                    ! param(1): P_boiling
                    ! param(2): T_boiling
                    ! param(3): Lv
                    ! param(4): W
                    ! psat =  P_boiling * exp( Lv * W / universal_gas_constant  *  (min(0.0_rp, T - T_boiling) ) / (T_boiling * T ) )
                    !
                    y =  param(1) * exp( param(3) * param(4) / universal_gas_constant  *  (min(0.0_rp,x - param(2)) ) / ( param(2) * x ) )
            end select
        endif
    end function physics_indexed_property_functions


    function physics_NSRDS_functions(itask,x,param,xcrit) result (y)
        integer(ip),          intent(in)    :: itask
        real(rp),             intent(in)    :: x
        real(rp),             intent(in)    :: param(:)
        real(rp),   intent(in), optional    :: xcrit
        real(rp)                            :: y

        real(rp)                            :: xnorm, eta, xbound


        xbound = max(1.0e-8_rp,x)
        y = 0.0_rp
        select case(itask)
           case(100)
              y = ((((param(6)*x + param(5))*x + param(4))*x + param(3))*x + param(2))*x + param(1)
           case(101)
              y      = exp(param(1) + param(2)/xbound + param(3)*log(xbound) + param(4)*xbound**param(5))
           case(102)
              y = param(1)* (xbound**param(2))/(1.0_rp + param(3)/xbound + param(5)/(xbound**2.0_rp))
           case(103)
              y = param(1) + param(2)*exp(-1.0_rp * param(3)/(xbound**param(4))) 
           case(104)
              y = param(1) + param(2)/xbound + param(3)/(xbound**3.0_rp) &
                 + param(4)/(xbound**8.0_rp) + param(5)/(xbound**9.0_rp)
           case(105)
              y = param(1) / (param(2)**(1.0_rp + (max(1.0e-8_rp,(1.0_rp - x/param(3)))**param(4))))
           case(106) 
              xnorm = min(xbound,xcrit-1.0e-8_rp)/xcrit 
              eta   = 1.0_rp - xnorm
              y     = param(1) * eta**(((param(5)*xnorm + param(4))*xnorm + param(3))*xnorm + param(2)) 
           case(107)
              y = param(1) + param(2)*((param(3)/xbound)/sinh(param(3)/xbound))**2.0_rp & 
                 +param(4)*((param(5)/xbound)/cosh(param(5)/xbound))**2.0_rp 
           case(114)
              xnorm = min(xbound,xcrit-1.0e-8_rp)/xcrit 
              eta   = 1.0_rp - xnorm
              y     =param(1)*param(1)/(eta + 1.0e-8_rp) + param(2) - eta &
                     *(                                                   &
                       2.0_rp*param(1)*param(3)                           &
                       + eta*(param(1)*param(4)                           &
                         + eta*(param(3)*param(3)/3.0_rp                  &
                           + eta*(0.5_rp*param(3)*param(4)                &
                             + 0.2_rp*param(4)*param(4)*eta)))            &
                      )   


        end select

    end function physics_NSRDS_functions


    !
    ! AB: function to calculate liquid properties 
    ! 
    subroutine physics_set_liquid_pressure(liq,P)
        type(liquid_state),   intent(inout) :: liq
        real(rp),             intent(in)    :: P

        real(rp)                :: T0, T1, T2, er, psat1, psat2, grad
        integer(ip)             :: iiter
        integer(ip), parameter  :: niter = 100_ip 
        
        liq % P = P
       
        if (liq % P >= liq % Pc) then
           liq % Tsat = liq % Tc
           call messages_live('    SUPERCRITICAL AMBIENT PRESSURE, P_CRIT= ', REAL_NUMBER = liq%Pc )
        else
           !
           ! Initialize
           !
           liq % Tsat = min(liq % Tc, 3000.0_rp)

           !
           ! Save temperature
           !
           T0 = liq % T
           er = 1.0e10_rp
           iiter = 0_ip
           call physics_set_liquid_temperature(liq, liq % Tc/2.0_rp)

           !
           ! Newton-Raphson's to calculate saturation temperature
           ! 
           do while ( (iiter < niter) .and. (abs(er) > 1.0e-6_rp) )
               
               er = (liq % psat - liq % P)
               T1 = liq % T
               psat1 = liq % psat
               
               call physics_set_liquid_temperature(liq, T1+1.0_rp)
               T2 = liq % T
               psat2 = liq % psat

               grad = (T2-T1)/max(psat2-psat1,1e-8_rp)
               
               call physics_set_liquid_temperature(liq, T1 - 0.5_rp * er * grad)

               iiter = iiter + 1_ip
           end do

           liq % Tsat = liq % T

           if (liq % Tsat > liq % Tc) then
              call messages_live('    SUPERCRITICAL AMBIENT PRESSURE, P_CRIT= ', REAL_NUMBER = liq%Pc )
           endif

           !
           ! Restore original state
           !
           call physics_set_liquid_temperature(liq, T0)

        endif

    end subroutine physics_set_liquid_pressure

  
end module mod_physics
!> @}
