program alya_terminal_velocity
  
  !-----------------------------------------------------------------------
  !> @{
  !> @file    alya_terminal_velocity.f90
  !> @author  Guillaume Houzeaux
  !> @date    22/02/2017
  !> @brief   Terminal velocity
  !> @details Compute the terminal velocity u of a particle in a fluid
  !>          by equalizing the gravity, buoyancy and drag forces,
  !>          by solving the following non-linear equation for u:
  !>
  !>          g * (rhop-rho) = 3/4 * rho * u^2 * Cd / d
  !>
  !>          See dragfo.f90 in Alya for the different available drag models.
  !>          To compile, execute ./alya-terminal-velocity.sh
  !>
  !>          Special case of Stokes flow:
  !>          Cd   = 24/Re
  !>          u    = g * (rhop-rho) * d^2 / (18 mu)
  !>          u(t) = (rhop-rho)*V*g/b * [1-exp(-b*t/m)], b = 3*pi*mu*d
  !>
  !> @} 
  !-----------------------------------------------------------------------
  
  use def_kintyp, only :  ip,rp
  implicit none
  real(rp), parameter  :: pi = 3.141592653589793238462643383279502884197_rp
  real(rp)             :: rhop,rho,mu,d,Cd,Re,g,V,b,m,u
  real(rp)             :: unew,eps,dCdRedRe,CdRe,psi,Ar,td,Ri
  integer(ip)          :: ii,imodel
  character(32)        :: cval

  imodel = 2
  
  rho    = 1.1839_rp                      ! Fluid density
  mu     = 1.86e-5_rp                     ! Fluid viscosity
  rhop   = 2500.0_rp                      ! Particle density
  d      = 1.0e-4_rp                      ! Particle diameter
  g      = 9.81_rp                        ! Gravity
  psi    = 1.0_rp                         ! Particle sphericity
  
  write(*,'("Drag model        (default=",i13,") = ")',ADVANCE='no') imodel
  read(*,'(a)') cval 
  if( trim(cval) /= '' ) read(cval,*) imodel

  write(*,'("Fluid density     (default=",e13.6,") = ")',ADVANCE='no') rho
  read(*,'(a)') cval 
  if( trim(cval) /= '' ) read(cval,*) rho

  write(*,'("Fluid viscosity   (default=",e13.6,") = ")',ADVANCE='no') mu
  read(*,'(a)') cval 
  if( trim(cval) /= '' ) read(cval,*) mu

  write(*,'("Particle density  (default=",e13.6,") = ")',ADVANCE='no') rhop
  read(*,'(a)') cval 
  if( trim(cval) /= '' ) read(cval,*) rhop

  write(*,'("Particle diameter (default=",e13.6,") = ")',ADVANCE='no') d
  read(*,'(a)') cval 
  if( trim(cval) /= '' ) read(cval,*) d

  write(*,'("Gravity           (default=",e13.6,") = ")',ADVANCE='no') g
  read(*,'(a)') cval 
  if( trim(cval) /= '' ) read(cval,*) g
  
  Re    = 1.0_rp
  u     = Re / ( rho * d ) * mu
  eps   = 1.0_rp
  ii    = 0

  do while( ii < 1000 .and. eps > 1.0e-6_rp )
     call dragfo(imodel,u,mu,rho,d,psi,CdRe,Cd,Re,dCdRedRe)     
     unew = sqrt( (rhop-rho)*g*4.0_rp/3.0_rp/rho*d/Cd)
     eps  = abs(unew-u)/unew
     u    = unew     
     ii   = ii + 1    
  end do
  !
  ! Derived parameters
  !
  V      = 4.0_rp/3.0_rp*pi*(d/2.0_rp)**3   ! Particle volume
  m      = rhop*V                           ! Particle mass
  b      = 3.0_rp*pi*mu*d
  Ar     = d**3*rho*(rhop-rho)*g/(mu**2)    ! Archimede's number
  Ri     = Ar/Re**2                         ! Richardson's number
  td     = (rhop*d**2) / (0.75_rp*mu*CdRe)  ! Relaxation time, td=m/b for Stokes flow
  
  if( ii >= 1000 ) then
     print*,'Did not converge'
  else
     write(6,'(a)')         ''
     write(6,'(a)')         '-------------------------------------------'
     write(6,'(a)')         ''
     write(6,'(a,i13)')     'Number of iterations= ',ii
     write(6,'(a,e13.6)')   'Mass=                 ',m
     write(6,'(a,e13.6)')   'Volume=               ',V
     write(6,'(a,e13.6)')   'Cd=                   ',Cd
     write(6,'(a,e13.6,a)') 'Re=                   ',Re,' (inertia/viscous)'
     write(6,'(a,e13.6,a)') 'Ar=                   ',Ar,' (gravity/viscous)'
     write(6,'(a,e13.6,a)') 'Ri=Ar/Re^2            ',Ri,' (gravity/interia)'
     write(6,'(a,e13.6,a)') 'Relaxation time=      ',td,' (time to reach terminal velocity)'
     write(6,'(a,e13.6)')   'Terminal velocity=    ',u
     write(6,'(a)')         ''
  end if

end program alya_terminal_velocity
