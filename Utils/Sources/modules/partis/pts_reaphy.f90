!------------------------------------------------------------------------
!> @addtogroup Partis
!! @{
!> @name    Partis read physical problem
!! @file    pts_reaphy.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!! @brief   This routine reads the physical problem
!! @details This routine reads the physical problem
!> @}
!------------------------------------------------------------------------

subroutine pts_reaphy()
  use def_parame,    only : xmaxint4
  use def_kintyp
  use def_master
  use def_kermod
  use def_inpout
  use def_domain
  use def_partis
  use mod_ecoute, only :  ecoute
  implicit none
  integer(ip) :: itype,iprop

  if( INOTSLAVE ) then

     mlagr     =  20                                        ! Lagrangian particles
     dtmin_pts = 1.0e-12_rp                                 ! Minimum time step
     dimin_pts = 0.0_rp                                     ! Minimum distance
     mean_free_path_pts = -69e-9_rp                         ! Molecular mean free path
     do itype = 1,mtyla
        parttyp(itype) % kfl_exist = 0                      ! Type does not exist
        parttyp(itype) % kfl_modla = 1                      ! Transport model
        parttyp(itype) % kfl_grafo = 0                      ! Gravity
        parttyp(itype) % kfl_buofo = 0                      ! Buoyancy
        parttyp(itype) % kfl_drafo = 0                      ! Drag force model
        parttyp(itype) % kfl_extfo = 0                      ! External force
        parttyp(itype) % kfl_brown = 0                      ! Brownian motion
        parttyp(itype) % kfl_saffm = 0                      ! Saffman force
        parttyp(itype) % kfl_schem = 0                      ! Time integration scheme: Analytical or numerical integration
        parttyp(itype) % kfl_tstep = 0                      ! Time step
        parttyp(itype) % denpa     = 0.0_rp                 ! Density particle
        parttyp(itype) % spher     = 1.0_rp                 ! Particle sphericity
        parttyp(itype) % diame     = 0.0_rp                 ! Particle diameter
        parttyp(itype) % diffu     = 0.0_rp                 ! Particle diffusion coefficient [m^2/s]
        parttyp(itype) % dtime     = 1.0_rp                 ! Default time step
        parttyp(itype) % safet     = 1.0_rp                 ! Safety factor
        parttyp(itype) % chale     = -1.0_rp                ! Characteristic length
        do iprop = 1,mlapr
           parttyp(itype) % prope(iprop) = 0.0_rp           ! Particle type properties
        end do
     end do

     !-------------------------------------------------------------------
     !
     ! Read/write unit
     !
     !-------------------------------------------------------------------
     lispa = 0
     lisda = momod(modul) % lun_pdata ! Reading file
     lisre = momod(modul) % lun_outpu ! Writing file

     !
     ! ADOC[0]> $-----------------------------------------------------------------------
     ! ADOC[0]> $ Physical properties definition
     ! ADOC[0]> $-----------------------------------------------------------------------
     ! ADOC[0]> PHYSICAL_PROBLEM
     !
     call ecoute('pts_reaphy')
     do while( words(1) /= 'PHYSI' )
        call ecoute('pts_reaphy')
     end do
     call ecoute('pts_reaphy')

     do while(words(1) /= 'ENDPH' )

        if( words(1) == 'MAXIM' ) then

           mlagr = getint('MAXIM',1_ip,'#MAXIMUM NB OF LAGRANGIAN PARTICLE')
           if( ip == 4 .and. (mlagr > xmaxint4 .or. mlagr < 0.0_rp) ) call runend('PTS_REAPHY: GOT TO 8 BYTES INTEGERS')

        else if( words(1) == 'DEPOS' ) then
           !
           ! ADOC[1]> DEPOSITION_DISTANCE: real                                                       $ Deposition distance
           ! ADOC[d]> DEPOSITION:
           ! ADOC[d]> Depopsition distance, under which the particle is considered as deposited
           !
           dimin_pts = getrea('DEPOS',0.0_rp,'#DEPOSITION DISTANCE TO WALL')

        else if( words(1) == 'TYPE ' ) then
           !
           ! ADOC[1]> TYPE: int                                                                       $ Definition of particle type int.
           ! ADOC[d]> TYPE:
           ! ADOC[d]> Definition of particle type int.
           !
           if( .not. exists('OFF  ') ) then

              itype = getint('TYPE ',1_ip,'#TYPE NUMBER OF LAGRANGIAN PARTICLE')
              if( itype < 1 .or. itype > mtyla ) call runend('PTS_REAPHY: WRONG PARTICLE TYPE')
              parttyp(itype) % kfl_exist = 1
              call ecoute('pts_reaphy')
              do while(words(1) /= 'ENDTY' )
                 !
                 ! ADOC[2]> MODEL : FORCE | VELOCITY ?                                               $ Definition of model particle   
                 ! ADOC[d]> MODEL:
                 ! ADOC[d]> Definition of model particle 
                 !
                 if( words(1) == 'MODEL' ) then
                    if( words(2) == 'FORCE' ) then
                       parttyp(itype) % kfl_modla = 2
                    else
                       parttyp(itype) % kfl_modla = 1
                    end if
                    !
                    ! ADOC[2]> TIME_STRATEGY : PRESCRIBED[VALUE=real] | EXTERIOR | ADAPTATIVE [ERROR|POSITION|FROM_CRITIQUAL|VELOCITY]
                    ! ADOC[d]> TIME_STRATEGY :
                    ! ADOC[d]> Definition of time strategy to solve the particle motion, prescribed: you fix the value, exterior: take the time step of the fluid 
                    ! ADOC[d]> adaptative is ....? 
                    !   
                 else if( words(1) == 'TIMES' ) then
                    if(      words(2) == 'PRESC' ) then
                       parttyp(itype) % kfl_tstep = -1
                       parttyp(itype) % dtime     = getrea('VALUE',1.0e-6_rp,'#TIME STEP OF PARTICLES')                       
                    else if( words(2) == 'EXTER' ) then
                       parttyp(itype) % kfl_tstep =  0
                    else if( words(2) == 'ADAPT' ) then
                       if(      exists('ERROR') ) then
                          parttyp(itype) % kfl_tstep =  1
                       else if( exists('POSIT') ) then
                          parttyp(itype) % kfl_tstep =  2
                       else if( exists('FROMC') ) then
                          parttyp(itype) % kfl_tstep =  3
                       else if( exists('VELOC') ) then
                          parttyp(itype) % kfl_tstep =  4
                       end if
                    end if
                    !
                    ! ADOC[3]> CHARACTERISTIC LENGTH : DIAMETER | ELEMENT | VISCOSITY | TAU | real                 $ Definition of characteristic length    
                    ! ADOC[d]> CHARACTERISTIC LENGTH :
                    ! ADOC[d]> Definition of characteristic length of the particle 
                    !
                    if( exists('CHARA') ) then
                       if( exists('DIAME') ) then
                          parttyp(itype) % chale =  0.0_rp
                       else if( exists('ELEME') ) then
                          parttyp(itype) % chale = -1.0_rp
                       else if( exists('VISCO') ) then
                          parttyp(itype) % chale = -2.0_rp
                       else if( exists('TAU  ') ) then
                          parttyp(itype) % chale = -3.0_rp
                       else
                          parttyp(itype) % chale = getrea('CHARA',1.0_rp,'#CHARACTERISTIC LENGTH')
                       end if
                    end if
                    !
                    ! ADOC[3]> SAFETY FACTOR : real                                                                 $ Definition of safety factor    
                    ! ADOC[d]> SAFETY FACTOR :
                    ! ADOC[d]> Definition of safety factor of the particle 
                    !
                    if( exists('SAFET') ) then
                       parttyp(itype) % safet = getrea('SAFET',1.0_rp,'#SAFETY FACTOR')
                    end if

                 else if( words(1) == 'MINIM' ) then
                    !
                    ! ADOC[2]> MINIMUM_TIME_STEP: real                                                         $ Minimum time step
                    ! ADOC[d]> MINIMUM_TIME_STEP:
                    ! ADOC[d]> Minimum time step for particles under which they disappear from the simulation
                    !
                    dtmin_pts = getrea('MINIM',1.0e-12_rp,'#MINIMUM TIME STEP')
                    !
                    ! ADOC[2]> GRAVITY : ON | OFF                                                         $ gravity force
                    ! ADOC[d]> GRAVITY :
                    ! ADOC[d]> add gravity force to the particle
                    !
                 else if( words(1) == 'GRAVI' ) then
                    if (words(2) == 'ON   ') parttyp(itype) % kfl_grafo = 1
                    !
                    ! ADOC[2]> BUOYANCY : ON | OFF                                                         $ buoyancy force
                    ! ADOC[d]> BUOYANCY :
                    ! ADOC[d]> add buoyancy force to the particle
                    !
                 else if( words(1) == 'BUOYA' ) then
                    if( words(2) == 'ON   ') parttyp(itype) % kfl_buofo = 1
                    !
                    ! ADOC[2]> DRAG : ON or GANSER , real(SPHERICITY OF THE PARTICLE TYPE)                  $ drag force
                    ! ADOC[d]> DRAG :
                    ! ADOC[d]> add drag force to the particle
                    !
                 else if( words(1) == 'DRAG ' .or. words(1) == 'DRAGF' ) then
                    if( words(2) == 'ON   ' .or. words(2) == 'GANSE' ) then
                       parttyp(itype) % kfl_drafo = 2
                       parttyp(itype) % spher = getrea('SPHER',1.0_rp,'#SPHERICITY OF THE PARTICLE TYPE')
                    end if

                 else if( words(1) == 'MEANF' ) then
                    !
                    ! ADOC[2]> MEAN_FREE_PATH: real                                                            $ Mean free path of medium
                    ! ADOC[d]> MEAN_FREE_PATH:
                    ! ADOC[d]> Mean free path lambda required for Cunningham correction factor. Put a positive value
                    ! ADOC[d]> to activate Cunningham correction.
                    ! ADOC[d]> <p>
                    ! ADOC[d]> http://en.wikipedia.org/wiki/Cunningham_correction_factor
                    ! ADOC[d]> <p>
                    ! ADOC[d]> In fluid dynamics, the Cunningham correction factor or Cunningham slip correction
                    ! ADOC[d]> factor is used to account for noncontinuum effects when calculating the drag on small
                    ! ADOC[d]> particles. The derivation of Stokes Law, which is used to calculate the drag force on
                    ! ADOC[d]> small particles, assumes a No-slip condition which is no longer correct at high
                    ! ADOC[d]> Knudsen number. The Cunningham slip correction factor allows predicting the drag force
                    ! ADOC[d]> on a particle moving a fluid with Knudsen number between the continuum regime and
                    ! ADOC[d]> free molecular flow. C=1+2*lambda/d*(1.257+0.4*exp(-0.55*d/lambda).
                    ! ADOC[d]> <p>
                    ! ADOC[d]> <pre>
                    ! ADOC[d]> Vacuum range           Pressure in hPa (mbar)      Molecules / cm3   Molecules / m3   Mean free path
                    ! ADOC[d]> ----------------------------------------------------------------------------------------------------
                    ! ADOC[d]> Ambient pressure       1013                        2.7 x 1019        2.7 x 1025       68 nm
                    ! ADOC[d]> Low vacuum             300 - 1[dubious - discuss]  1019 - 1016       1025 - 1022      0.1 - 100 micro-m
                    ! ADOC[d]> Medium vacuum          1 - 10-3                    1016 - 1013       1022 - 1019      0.1 - 100 mm
                    ! ADOC[d]> High vacuum            10-3 - 10-7                 1013 - 109        1019 - 1015      10 cm - 1 km
                    ! ADOC[d]> Ultra high vacuum      10-7 - 10-12                109 - 104         1015 - 1010      1 km - 105 km
                    ! ADOC[d]> Extremely high vacuum  <10-12                      <104              <1010            >105 km
                    ! ADOC[d]> </pre>
                    !
                    mean_free_path_pts = getrea('MEANF',69e-9_rp,'#MEAN FREE PATH')

                 else if( words(1) == 'EXTER' ) then
                    !
                    ! ADOC[2]> EXTERIOR: FUNCTION: int                                                     $ exterior force acting on the particle 
                    ! ADOC[d]> EXTERIOR:
                    ! ADOC[d]> exterior force acting on the particle 
                    ! 
                    if( words(2) == 'FUNCT') &
                         parttyp(itype) % kfl_extfo = getint('FUNCT',1_ip,'#EXTERNAL FORCE FUNCTION NUMBER')

                 else if( words(1) == 'TIMEI' ) then
                    !
                    ! ADOC[2]> TIME_INTEGRATION: ANALY | NUMER                                                   $ Analytical or numerical integration
                    ! ADOC[d]> TIME_INTEGRATION:
                    ! ADOC[d]> Analytical or numerical (Newmark + Newton-Raphson) integration
                    !
                    if( words(2) == 'ANALY' ) then
                       parttyp(itype) % kfl_schem = 1
                    else if( words(2) == 'NUMER' ) then
                       parttyp(itype) % kfl_schem = 0
                    else if( words(2) == 'RK4  ' ) then
                       parttyp(itype) % kfl_schem = 2
                    else if( words(2) == 'RK45 ' ) then
                       parttyp(itype) % kfl_schem = 3
                    end if

                 else if( words(1) == 'FORCE' ) then
                    !
                    ! ADOC[2]> FORCE: ALL , FUNCTION [int] | SPHERICITY [real] | CHENG | ARAST | STOKE | DALLA      $ force acting on the particle 
                    ! ADOC[2]> FORCE: GRAVI , BUOYA, SAFFM , DRAG | SPHERICITY [real] | CHENG | ARAST | STOKE | DALLA, EXTER [int] 
                    ! ADOC[d]> FORCE:
                    ! ADOC[d]> force acting on the particle 
                    ! 
                    if( words(2) == 'ALL  ' ) then

                       parttyp(itype) % kfl_grafo = 1
                       parttyp(itype) % kfl_buofo = 1
                       parttyp(itype) % kfl_drafo = 2
                       parttyp(itype) % kfl_extfo = 1
                       parttyp(itype) % kfl_saffm = 1
                       if( exists('FUNCT') ) &
                            parttyp(itype) % kfl_extfo = getint('FUNCT',1_ip,'#EXTERNAL FORCE FUNCTION NUMBER')
                       if( exists('SPHER') ) &
                            parttyp(itype) % spher = getrea('SPHER',1.0_rp,'#SPHERICITY OF THE PARTICLE TYPE')
                       if( exists('CHENG') ) &
                            parttyp(itype) % kfl_drafo = 1
                       if( exists('ARAST') ) &
                            parttyp(itype) % kfl_drafo = 3
                       if( exists('STOKE') ) &
                            parttyp(itype) % kfl_drafo = 6 ! Stokes
                       if( exists('DALLA') ) &
                            parttyp(itype) % kfl_drafo = 7 ! Dalavalle
                       if( exists('WALLC') ) &
                            parttyp(itype) % kfl_drafo = 8 ! Wall Correction Zeng et al:The effects of near wall corrections to
                                                           ! hydrodynamic forces on particle deposition and transport in vertical...
                    else
                       if( exists('GRAVI') ) parttyp(itype) % kfl_grafo = 1
                       if( exists('BUOYA') ) parttyp(itype) % kfl_buofo = 1
                       if( exists('SAFFM') ) then
                          parttyp(itype) % kfl_saffm = 1 ! Standard Saffman lift force
                          if( exists('WALLS') ) &
                               parttyp(itype) % kfl_saffm = 2  ! wall effects for stationary particle, Zeng et al:The effects of near wall...
                          if( exists('WALLC') ) &
                               parttyp(itype) % kfl_saffm = 3  ! wall effects for a moving particle, Zeng et al:The effects of near wall...
                       end if
                       if( exists('DRAG ') .or. exists('DRAGF') ) then
                          parttyp(itype) % kfl_drafo = 2
                          parttyp(itype) % spher = getrea('SPHER',1.0_rp,'#SPHERICITY OF THE PARTICLE TYPE')
                          if( exists('CHENG') ) &
                               parttyp(itype) % kfl_drafo = 1
                          if( exists('ARAST') ) &
                               parttyp(itype) % kfl_drafo = 3
                          if( exists('STOKE') ) &
                               parttyp(itype) % kfl_drafo = 6 ! Stokes
                          if( exists('DALLA') ) &
                               parttyp(itype) % kfl_drafo = 7 ! Dalavalle
                          if( exists('WALLC') ) &
                               parttyp(itype) % kfl_drafo = 8 ! Wall Correction Zeng et al:The effects of near wall corrections to
                                                           ! hydrodynamic forces on particle deposition and transport in vertical...
                       end if
                       if( exists('EXTER') ) then
                          parttyp(itype) % kfl_extfo = getint('EXTER',1_ip,'#EXTERNAL FORCE FUNCTION NUMBER')
                       end if

                    end if
                    !
                    ! ADOC[2]> DENSITY: real                                                $ density of the particle 
                    ! ADOC[d]> DENSITY:
                    ! ADOC[d]> density of the particle 
                    ! 
                 else if( words(1) == 'DENSI' ) then
                    parttyp(itype) % denpa = getrea('DENSI',1.0_rp,'#PARTICLE DENSITY')
                    !
                    ! ADOC[2]> DIAMETER: real                                                $ diameter of the particle 
                    ! ADOC[d]> DIAMETER:
                    ! ADOC[d]> diameter of the particle 
                    ! 
                 else if( words(1) == 'DIAME' ) then
                    parttyp(itype) % diame = getrea('DIAME',1.0_rp,'#PARTICLE DIAMETER')
                    !
                    ! ADOC[2]> SPHERICITY : real                                              $ sphericity of the particle
                    ! ADOC[d]> SPHERICITY :
                    ! ADOC[d]> sphericity of the particle 
                    ! 
                 else if( words(1) == 'SPHER' ) then
                    parttyp(itype) % spher = getrea('SPHER',1.0_rp,'#SPHERICITY OF THE PARTICLE TYPE')
                    !
                    ! ADOC[2]> BROWNIAN_MOTION:   ON [,DIFFUSION = real] | FORCE [,DIFFUSION = real]                      $ Brownian motion
                    ! ADOC[d]> BROWNIAN_MOTION:
                    ! ADOC[d]> Add a Brownian motion to the particles after advecting them, using as DIFFUSION
                    ! ADOC[d]> coefficient real. If the option is not present, then the diffusion coefficient
                    ! ADOC[d]> is computed using Einstein's law: D = kb * T / (6.0_rp*pi*mu*r).
                    ! ADOC[d]> The viscosity should therefore be defined in Kermod and the temperature solved.
                    ! ADOC[d]> If the temperature is not available, then a temperature of 20C is used.
                    ! ADOC[d]> 2 options are available ....
                    !                   
                 else if( words(1) == 'BROWN' ) then                  
                    if( words(2) == 'ON   ' .or. words(2) == 'DISPL' ) then
                       parttyp(itype) % kfl_brown = 1
                       parttyp(itype) % diffu = getrea('DIFFU',0.0_rp,'#PARTICLE DIFFUSION COEFFICIENT M**2/S')
                    else if (words(2) == 'FORCE ') then
                       parttyp(itype) % kfl_brown = 2
                       parttyp(itype) % diffu = getrea('DIFFU',0.0_rp,'#PARTICLE DIFFUSION COEFFICIENT M**2/S')
                    end if
                 end if
                 call ecoute('pts_reaphy')
              end do
           end if
           
        end if

        call ecoute('pts_reaphy')
     end do

  end if
end subroutine pts_reaphy
