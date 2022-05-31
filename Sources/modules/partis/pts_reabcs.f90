!-----------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @file    pts_reabcs.f90
!> @author  Guillaume Houzeaux
!> @date    16/11/1966   
!> @brief   Read boundary conditions 
!> @details The boundary codes are:
!>          KFL_FIXBOU(IBOUN) = 0 ... Outflow
!>                              1 ... Wall
!>                              3 ... Slip
!>                              4 ... bounce 
!>
!>  To read prescribed velocity fields, these need to be specified in dom.dat and here in pts.dat. Periodic fields are supported. 
!>  Linear interpolation in time is also supported
!>  boundary_conditions
!>     VELOC = field_id
!> 
!>
!>  For bouncing boundaries, WALLD must not be used. Kernel needs only solver for normal propagation
!>
!>  For slip boundairies, WALL_NORMAL nees to be added to ker.dat, potentially for partis
!>  also WALL_DISTANCE needs to be added. Sample definitions can be like this:
!>     NUMERICAL_TREATMENT
!>       WALL_DISTANCE                                
!>         ALGEBRAIC_SOLVER    DEFLATED_CG, ITERA= 5000,TOLER= 1e-9, ADAPTIVE, RATIO=1e-9
!>         PRECONDITIONING     DIAGONAL               
!>         CODES, BOUNDARIES
!>           1 3
!>         END_CODES
!>       END_WALL_DISTANCE 
!>     
!>       WALL_NORMAL
!>         ALGEBRAIC_SOLVER    DEFLATED_CG, ITERA= 5000,TOLER= 1e-9, ADAPTIVE, RATIO=1e-9
!>         PRECONDITIONING     DIAGONAL               
!>         CODES, BOUNDARIES
!>           1 3
!>         END_CODES
!>       END_WALL_NORMAL
!>     END_NUMERICAL_TREATMENT
!>
!>   For everything else:
!>
!>       PHYSICAL_PROBLEM
!>         MAXIMUM int    $maximum number of particles to preallocate memory, if it is smaller than the number of injected particles, the memory will be reallocated
!>         TYPE= int  $1,2,...
!>           MODEL= FORCE
!>           FORCES= DRAG, GRAVITY, BUOYANCY
!>           DENSITY= real
!>           DIAMETER= real
!>         END_TYPE
!>       END_PHYSICAL_PROBLEM
!>
!>
!>       BOUNDARY_CONDITIONS
!>         INJECTOR= int
!>           GEOMETRY: injector_name, PARAM= injector dependent numbers
!>           TYPE:      int | ALL   $ type of the particles to inject (provided in the PHYSICAL_PROBLEM section
!>           DISTRIBUTION:         UNICA | UNIPO   $ spatial particle disrtibution, optional
!>           NPARTICLES: ASIS | [ MASS_FLOW:  CONSTANT= mdot | TIME_FUNCTION=NAMEFUN ] $ mass flow overwrites the number of particles given in the geometry  
!>           STOCA: ON  $random particles, otherwise nonrandom (no clue what is the nonrandom)
!>           INITIAL_TIME: real  $initial time of injection
!>           FINAL_TIME: real  $final time of injection
!>           PERIOD_TIME: real  $time step of injection
!>           VELOCITY: SPRAY, PARAM=velocity, angle_degrees
!>           SIZE_DISTRIBUTION: CONSTANT | UNIFORM PARAM= dmin, dmax | ROSIN_RAMMLER PARAM= dmin, dmax, dmean, nparam $ size distribution overwrites the diameter from the type 
!>         END_INJECTOR
!>         CODES, BOUNDARIES
!>           ...
!>           int1 1                                    $   int1=id of the boundary, 1 - wall
!>           ...
!>         END_CODES
!>       END_BOUNDARY_CONDITIONS
!>
!>       Injector types:
!>           GEOMETRY: CIRCLE, PARAMS=x_center, y_cenetr, z_center, radius, x_normal, y_normal, z_normal, number_of_particles
!>
!>           The normal has to point in the direction of injection.
!>
!>           The real number of particles injected by default is number_of_particles^2, 
!>           unless NPARTICLES: ASIS is specified, then number_of_particles is injected.
!>
!>           DISTRIBUTION parameter (works only for CIRCLE injector for 3D meshes) specifies 
!>           the spatial distrbution of the particles, by default it's uniform in polar
!>           coordinates UNIPO (as it was implmented). Specifying UNICA will produce uniform
!>           distribution in cartesian coordinates.
!>        
!>        Other injector types to be added.
!> @} 
!-----------------------------------------------------------------------

subroutine pts_reabcs()

  use def_kintyp
  use def_master
  use def_kermod
  use def_inpout
  use def_domain
  use def_partis
  use mod_opebcs
  use mod_memory
  use mod_ecoute,   only : ecoute
  use mod_opebcs,   only : boundary_conditions_read_boundary_codes
  use mod_opebcs,   only : opebcs_initialization_structure
  use mod_opebcs,   only : opebcs_initialization_variable
  use mod_messages, only : livinf
  use mod_ker_space_time_function
  implicit none
  character(len=8) :: fmt,num1,str2
  integer(ip)      :: ipara,iinj,ii
  real(rp)         :: normlen
  character(5) :: wfname
  !
  ! Allocate memory for boundary codes
  !
  call opebcs_initialization_structure(1_ip,tbcod_pts)     ! Velocity 
  call opebcs_initialization_variable (1_ip,tbcod_pts)

  if( INOTSLAVE ) then

     kfl_injve_pts                  = 0
     kfl_injte_pts                  = 0
     kfl_injma_pts                  = 0
     kfl_injsd_pts                  = 0
     kfl_injVelocFunction_pts       = 0
     kfl_injMassFunction_pts        = 0
     ninj_pts                       = 0
     tinla_pts                      = 0.0_rp                   
     tfila_pts                      = 1.6_rp                   
     tpela_pts                      = 1.0e6_rp
     kfl_injla_pts                  = 0
     kfl_injty_pts                  = 0
     kfl_random_pts                 = 0
     kfl_boundary_injection         = 0
     codbo_pts                      = 0_ip 
     parla_pts                      = 0.0_rp
     injector_particle_distribution = 0
     injector_npts_asis             = 0
     !
     ! Reach section
     !
     call ecoute('pts_reabcs')
     do while( words(1) /= 'BOUND' )
        call ecoute('pts_reabcs')
     end do
     call ecoute('pts_reabcs')
     !
     ! Read boundary conditions field
     !
     do while( words(1) /= 'ENDBO' )

        if( words(1) == 'CODES'.and.exists('BOUND') ) then
           ! 
           ! User-defined codes on boundaries
           !          
           kfl_fixbo => kfl_fixbo_pts
           bvnat     => bvnat_pts
           tbcod     => tbcod_pts(1:)
           call boundary_conditions_read_boundary_codes('PARTICLES')
           !call reacod(2_ip)

        else if( words(1) == 'INJEC') then
           !
           ! Injection field
           !
           iinj = getint('INJEC',1_ip,'#INJECTOR')
           if( iinj > pts_minj ) call runend('PTS_REABCS: WRONG INJECTOR NUMBER, INCREASE PTS_MINJ IN DEF_PARTIS... SORRY FOR THAT!')
           call ecoute('pts_reabcs')
           do while( words(1) /= 'ENDIN' )
              if( words(1) == 'GEOME' .or. words(1) == 'BOUND' ) then
                 !
                 ! Get boundary number
                 !
                 if( words(1) == 'BOUND' ) then
                    kfl_boundary_injection = 1_ip
                    codbo_pts(iinj) = getint('BOUND ',0_ip,'BOUNDARY CODE')
                    if( words(3) == 'PARTI') then
                       parla_pts(iinj,8) = getrea('PARTI',0.0_rp,'#PARTICLES')                      
                    end if
                 end if
                 !
                 ! Geometry of the injector
                 !
                 do ipara = 1,min(mpala,nnpar)
                    parla_pts(iinj,ipara) = param(ipara+2)
                 end do
                 if( words(2) == 'SQUAR' ) then
                    kfl_injla_pts(iinj) = 1
                 else if( words(2) == 'SPHER' ) then
                    kfl_injla_pts(iinj) = 2
                 else if( words(2) == 'SEMIS' ) then
                    kfl_injla_pts(iinj) = 3
                 else if( words(2) == 'CIRCL' ) then
                    kfl_injla_pts(iinj) = 4
                    !
                    ! Normalize normal vector
                    !
                    normlen = sqrt(parla_pts(iinj, 5)**2 + parla_pts(iinj, 6)**2 + parla_pts(iinj, 7)**2)
                    parla_pts(iinj, 5) = parla_pts(iinj, 5)/normlen 
                    parla_pts(iinj, 6) = parla_pts(iinj, 6)/normlen
                    parla_pts(iinj, 7) = parla_pts(iinj, 7)/normlen
                 else if( words(2) == 'RECTA' ) then
                    kfl_injla_pts(iinj) = 5
                 else if( words(2) == 'POINT' ) then
                    kfl_injla_pts(iinj) = 6
                 else if( words(2) == 'CONE'  ) then
                    kfl_injla_pts(iinj) = 7
                 else if( words(2) == 'RANDO' ) then
                    kfl_injla_pts(iinj) = 8
                 else if( words(2) == 'SEGME' ) then
                    kfl_injla_pts(iinj) = 9
                 else if( words(2) == 'ANNUL' ) then
                    kfl_injla_pts(iinj) = 11
                    !
                    ! Normalize normal vector
                    !
                    normlen = sqrt(parla_pts(iinj, 6)**2 + parla_pts(iinj, 7)**2 + parla_pts(iinj, 8)**2)
                    parla_pts(iinj, 6) = parla_pts(iinj, 6)/normlen 
                    parla_pts(iinj, 7) = parla_pts(iinj, 7)/normlen
                    parla_pts(iinj, 8) = parla_pts(iinj, 8)/normlen
                 end if
              

              else if ( words(1) == 'COORD' ) then
                 injection_pts(iinj) % number_particles=getint('NUMBE',1_ip,'*PARTICLES')
                 call livinf(-9_ip,'INJECTED PARTICLES FROM FILE   ',injection_pts(iinj) % number_particles )
                 kfl_injla_pts(iinj) = 10
                 call memory_alloca(mem_modul(1:2,modul),'INJECTION_PTS(iinj) % COORD_PARTICLES','pts_reabcs',injection_pts(iinj) % coord_particles,ndime*injection_pts(iinj) % number_particles) 
                 call ecoute('pts_reabcs')
                 ii = 1
                 do while( words(1) /= 'ENDCO' )
                    if( ii > injection_pts(iinj) % number_particles*ndime ) call runend('PTS_REABCS: WRONG NUMBER OF PARTICLES')
                    injection_pts(iinj) % coord_particles(ii:(ii-1)+ndime ) = param(1:ndime)
                    call ecoute('pts_reabcs')
                    ii = ii + ndime
                 end do

              else if( words(1) == 'DISTR' ) then
                 if( words(2) == 'UNICA' ) then ! uniform distribution in cartesian coordinates
                     injector_particle_distribution(iinj) = 1
                 else if( words(2) == 'UNIPO') then ! uniform distribution in polar coordinates
                     injector_particle_distribution(iinj) = 0_ip
                 else
                    call runend("pts_reabcs: Particle distribution " // words(2) // " is not implemented")
                 end if
                  
              else if( words(1) == 'NPART' ) then
                 if ( words(2) == 'MASSF' ) then
                    !
                    ! Read type of mass flow 
                    !
                    injector_npts_asis(iinj) = 2

                    if ( exists('CONST') ) then
                       param_mass_pts(1, iinj) = getrea('CONST',0.0_rp,'#MASS FLOW OF INJECTOR')
                    elseif ( exists('SPACE') ) then
                       wfname                 = getcha('SPACE','NULL ','#Space/time Function name')
                       kfl_injMassFunction_pts(iinj) = space_time_function_number(wfname)        
                    else
                       call runend('In pts_reabcs: Mass flow rate should be specified')
                    endif

                 elseif ( words(2) == 'ASIS ' ) then
                     injector_npts_asis(iinj) = 1
                 else
                     injector_npts_asis(iinj) = 0
                 end if

              else if( (words(1) == 'STOCA') .or. (words(1) == 'STOCH') ) then
                 kfl_random_pts(iinj) = 1
              else if( words(1) == 'PIPE' ) then
                 kfl_random_pts(iinj) = 2

              else if( words(1) == 'TYPE ' ) then
                 !
                 ! On which type to apply
                 !
                 if( words(2) == 'ALL  ' ) then
                    kfl_injty_pts(iinj) = 0
                 else
                    kfl_injty_pts(iinj) =  getint('TYPE ',0_ip,'#TYPE ON WHICH TO APPLY INJECTOR')
                 end if
              
              else if( words(1) == 'INITI' ) then
                 !
                 ! Initial injection time
                 !
                 tinla_pts = getrea('INITI',1.0_rp,'#INITIAL INJECTION TIME LAGRANGIAN PARTICLE')
              else if( words(1) == 'PERIO' ) then
                 !
                 ! Adding a velocity field (unifrom or gaussian) in injected particles
                 !
                 tpela_pts = getrea('PERIO',1.0_rp,'#PERIOD INJECTION TIME LAGRANGIAN PARTICLE')
              else if( words(1) == 'FINAL' ) then
                 !
                 ! Final injection time
                 !
                 tfila_pts = getrea('FINAL',1.0_rp,'#FINAL INJECTION TIME LAGRANGIAN PARTICLE')
                 
              else if( words(1) == 'VELOC' ) then
                 !
                 ! Injection velocity
                 !
                 if(     words(2) == 'ZERO ' ) then
                    kfl_injve_pts(iinj) = -1
                 else if( words(2) == 'FLUID' ) then
                    kfl_injve_pts(iinj) =  0
                 else if( words(2) == 'UNIFO' ) then
                    kfl_injve_pts(iinj) =  1
                 else if( words(2) == 'GAUSS' ) then
                    kfl_injve_pts(iinj) =  2
                 else if( words(2) == 'CONIC' ) then
                    kfl_injve_pts(iinj) =  3
                 else if( words(2) == 'SPRAY' ) then
                    kfl_injve_pts(iinj) =  4
                 else if( words(2) == 'CONST' ) then
                    kfl_injve_pts(iinj) =  5
                 end if

                 do ipara = 1,min(size(param_veloc_pts,DIM=1,KIND=ip),nnpar)
                    param_veloc_pts(ipara,iinj) = param(ipara+2)
                 end do

                 if ( exists('SPACE') ) then
                     wfname                 = getcha('SPACE','NULL ','#Space/time Function name')
                     kfl_injVelocFunction_pts(iinj) = space_time_function_number(wfname)        
                 endif

              else if( words(1) == 'TEMPE' ) then
                 !
                 ! Injection tempreature
                 !
                 if(      words(2) == 'FLUID' ) then
                    kfl_injte_pts(iinj) =  0
                 else if( words(2) == 'CONST' ) then
                    kfl_injte_pts(iinj) =  5
                 end if

                 do ipara = 1,min(size(param_tempe_pts,DIM=1,KIND=ip),nnpar)
                    param_tempe_pts(ipara,iinj) = param(ipara+2)
                 end do

              else if( words(1) == 'SIZED' ) then
                 !
                 ! Size distribution
                 !
                 if(      words(2) == 'CONST' ) then
                    kfl_injsd_pts(iinj) =  0
                 else if( words(2) == 'UNIFO' ) then
                    kfl_injsd_pts(iinj) =  1
                 else if( words(2) == 'ROSIN' ) then
                    kfl_injsd_pts(iinj) =  2
                 end if

                 do ipara = 1,min(size(param_sized_pts,DIM=1,KIND=ip),nnpar)
                    param_sized_pts(ipara,iinj) = param(ipara+2)
                 end do
                
                
              end if
              call ecoute('pts_reabcs')
           end do

        else if( words(1)(1:4) == 'INJE') then
           !
           ! Inject particles with a saqure, sphere, semi-sphere or a circle
           !
           fmt = '(I1.1)'
           num1 = words(1)(5:5)
           read (num1, fmt) iinj
           if( iinj > 0 .and. iinj <= pts_minj) then
              if( words(2) == 'SQUAR' ) then
                 kfl_injla_pts(iinj) = 1
              else if( words(2) == 'SPHER' ) then
                 kfl_injla_pts(iinj) = 2
              else if( words(2) == 'SEMIS' ) then
                 kfl_injla_pts(iinj) = 3
              else if( words(2) == 'CIRCL' ) then
                 kfl_injla_pts(iinj) = 4
              else if( words(2) == 'RECTA' ) then
                 kfl_injla_pts(iinj) = 5
              else if( words(2) == 'POINT' ) then
                 kfl_injla_pts(iinj) = 6
              else if( words(2) == 'CONE'  ) then
                 kfl_injla_pts(iinj) = 7
              else if( words(2) == 'RANDO' ) then
                 kfl_injla_pts(iinj) = 8
              else if( words(2) == 'SEGME' ) then
                 kfl_injla_pts(iinj) = 9
              else if( words(2) == 'ANNUL' ) then
                 kfl_injla_pts(iinj) = 11
              end if
              do ipara = 1,min(mpala,nnpar)
                 parla_pts(iinj,ipara) = param(ipara+2)
              end do
           else
              write(str2,fmt) pts_minj
              call runend("NUMBER OF INJECTIONS MUST BE BETWEEN 1 AND "//str2)
           end if


        else if( words(1) == 'INITI' ) then
           !
           ! Initial injection time
           !
           tinla_pts = getrea('INITI',1.0_rp,'#INITIAL INJECTION TIME LAGRANGIAN PARTICLE')
        else if( words(1) == 'FINAL' ) then
           !
           ! Final injection time
           !
           tfila_pts = getrea('FINAL',1.0_rp,'#FINAL INJECTION TIME LAGRANGIAN PARTICLE')
        else if( words(1) == 'PERIO' ) then
           !
           ! Adding a velocity field (unifrom or gaussian) in injected particles
           !
           tpela_pts = getrea('PERIO',1.0_rp,'#PERIOD INJECTION TIME LAGRANGIAN PARTICLE')
        else if( words(1) == 'VELOC' ) then
           !
           ! Injection velocity
           !
           call runend('PTS_REABCS: OBSOLETE WAY OF DEFINING INJECTION VELOCITY')

        end if

        call ecoute('pts_reabcs')
     end do

  end if

end subroutine pts_reabcs
