!------------------------------------------------------------------------
!> @addtogroup Partis
!> @{
!> @name   Partis inner iteration
!> @file    pts_soltie.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!> @brief   This routine integrates particle paths
!> @details
!>
!>    ALGORITHM
!>    ---------
!>
!>       We are at time n+1:
!>
!>       do while( no han llegado todas a su tiempo )
!>
!>          Loop over particles
!>            do while( t < tf )
!>              Start from last element: ielem
!>              x^n+1 = x^n + dt * u
!>              Look for neighbors of ielem
!>              - if found
!>                  Validate new position
!>                  x^n = x^n+1
!>                  Save element
!>                  if( fringe element ) save particle to be sent
!>              - else
!>                  dt <= dt / 2
!>              - endif
!>              t = t + dt
!>            end do
!>          end do
!>          Send what I have to send
!>          Receive what I have to receive
!>       end do
!>
!>    NEWMARK
!>    -------
!>
!>    The Newmark scheme has a small inconsistency. The mag force is computed at k
!>    and not k+1
!>
!>    v^{k+1}_NS => F_drag                             x^{k+1} = f(x^k, a^{k+1}, a^k)
!>                           => F^{k+1} => a^{k+1} =>
!>    x^k        => F_mag                              u^{k+1} = f(a^{k+1}, a^k)
!>
!>    x^{n+1} = x^n + u^n dt + 1/2 a^n dt^2 + 1/6 a^n dt^3
!>    a^{n+1} = a^n + a'^n dt + 1/2 a''^n dt^2 =>
!>    a'^n dt = ( a^{n+1} - a^n ) / dt =>
!>
!>    u^{n+1} = u^n + a^n dt + gamma dt ( a^{n+1} - a^n )
!>    x^{n+1} = x^n + u^n dt + 1/2 a^n dt^2 + beta dt^2 ( a^{n+1} - a^n )
!>
!>    x^{n+1} = x^n + dt [ u^{n+1} - a^n dt - gamma * dt ( a^{n+1} - a^n ) ] + 1/2 a^n dt^2 + beta dt^2 ( a^{n+1} - a^n )
!>    x^{n+1} = x^n + dt u^{n+1} - a^n dt^2 - gamma * dt^2 ( a^{n+1} - a^n ) ] + 1/2 a^n dt^2 + beta dt^2 ( a^{n+1} - a^n )
!>    x^{n+1} = x^n + dt u^{n+1} + dt^2 ( - a^n  - gamma * ( a^{n+1} - a^n ) + 1/2 a^n + beta ( a^{n+1} - a^n ) ]
!>    x^{n+1} = x^n + dt u^{n+1} + dt^2 [ (-1/2+gamma-beta) a^n  + (-gamma+beta)  a^{n+1} ]
!>
!>                            gamma   beta
!>    Fox Goodwin             1 / 2   1 / 12    conditionnaly stable
!>    Linear acceleration     1 / 2   1 / 6     conditionnaly stable
!>    Average acceleration    1 / 2   1 / 4     unconditionnaly stable
!>    Nosotros                0.75    0.390625  diffusive
!>    Nosotros                1.0     0.5625    super diffusive
!>    External                1 / 2   0         explicit
!>
!>    Relation between beta and gamma: beta = 0.25*(gamma+1/2)^2
!>
!>    RANDOM WALK
!>    -----------
!>
!>    It is also possible to simulate the diffusion process by using a random walk
!>    algorithm. Once the new position of the particle has been found, the
!>    position is corrected as:
!>
!>    x^{n+1} <= x^{n+1} + eps * (2D*dt)^1/2,
!>
!>    where D is the diffusion coefficient [m^2/s] and eps follows a normal distribution.
!>
!>    The diffusion is given by
!>           k*T
!>    D = ---------
!>        6*pi*mu*r
!>
!>    where T is temperature, k is Boltzmann's constant, mu viscosity and r ths solute radius.
!>
!>    The process is the following:
!>
!>    1. Generate two random numbers U1 and U2 in [0,1]
!>    2. Perform a box-muller transformation to generate from U1 and U2 a normal distribution
!>       eps1 = sqrt( -2 *log(U1) ) * cos(2*pi*U2)
!>       eps2 = sqrt( -2 *log(U1) ) * sin(2*pi*U2)
!>    3. Update the position
!>       x <= x + eps1 * (2D*dt)^1/2
!>       y <= y + eps1 * (2D*dt)^1/2
!>
!>    A simple test can be peformed to check the normal distribution:
!>
!>       The probablity P(r,t) is the normal defined as
!>
!>       P(r,t) = 1/(4*D*pi*t)^(d/2) * exp( -r^2/(4*D*t) )
!>
!>       for a d-dimensional diffusion problem.
!>
!>       k is related to the diffusion as k^2 = 2D and P satisfies the
!>       2-dimensional Poisson equation:
!>       dP/dt = D ( d^2P/dx^2 + d^2P/dy^2 )
!>
!>    1. Inject particles at 0.05,0.05 in a domain [0,0] x [0.01,0.01].
!>    2. Choose zero convection velocity; D = 5.35*10^-5 [m^2/s];
!>       dt = 10^-5 [s]. With dt = 10^-4 [s], statistics is bad.
!>    3. At a given time, say t=0.01 [s], count the number of particles
!>       with radius in ranges from [0:R] where R goes from 0 to 0.01.
!>    4. Divide the results by the total number of particles.
!>    5. Compare the results with the following 2D normal distribution:
!>
!>       f(R) = \int_0^{2pi} dtheta \int_0^R r P(r,t) dr
!>            = 1 - exp( -R^2/(4*D*t) )
!>
!>       using \int x*exp(-c*x^2) dx = -1/2c * exp(-c*x^2)
!>
!>       Note that in 2D we have \int_0^2pi dtheta \int_0^infty r P(r,t) dr = 1.
!>
!>       In 3D, ther result is:
!>
!>       f(R) = \int_0^{2pi} dtheta \int_0^pi sin(phi) dphi \int_0^R r^2 P(r,t) dr
!>            = 4*pi/(4*D*pi*t)^3/2 \int_0^R r^2 * exp( -r^2/(4*D*t) ) dr
!>
!>    INJECTION
!>    ---------
!>
!>       KFL_MODLA ......... Model for Lagrangian transport (0=no particle,1=velocity,1=drag)
!>       KFL_INJLA ......... Injection model
!>       TINLA ............. Initial time of injection
!>       TPELA ............. Time period of injection
!>       MPALA ............. Maximum number of parameters
!>       PARLA(MPALA) ...... Parameters for the injection
!>
!>    OTHERS
!>    ------
!>
!>       MLAGR ............. Max. number of particle in each subdomain
!>       NLAGR ............. Number of particles in each subdomain (just needed for info)
!>       NLACC_PTS ......... Total number of existing and disappeared particles in all subdomains
!>       NLAGR_EXISTING_PTS ........... Total number of existing particles in all subdomains
!>       NLAGR_NON_MIGRATING_PTS ........... Particles going from one subdomain to another
!>       NLAGR_GOING_OUT_PTS ........... Particles deposited: going out of the computational domain
!>       NLAGR_ZERO_TIME_PTS ........... Particles that vanish because of zero time step
!>       NLAGR_DEPOSITED_PTS ........... Particles deposited, bue boundary not found
!>       PELEL_2,LELEL_2 ... Element connectivity linked list including
!>                           the neighbors element in parallel
!>
!>    LAGRTYP definition
!>    ------------------
!>
!>       Transport of Lagrangian particles. For particle ilagr
!>       lagrtyp(ilagr) % ilagr     =  i ......... Particle number
!>       lagrtyp(ilagr) % kfl_exist = -1 ......... I am the owner
!>                                  = -2 ......... Particle is deposited: go out of domain
!>                                  = -3 ......... Particle vanishes: time step too small
!>                                  = -4 ......... Particle is out of the flow
!>                                  = -5 ......... Particle has just been injected
!>                                  =  i ......... Send particle to neighbor i
!>                                  =  0 ......... No particle at that position ILAGR
!>       lagrtyp(ilagr) % coord(3)  = x,y,z ...... Coordinates
!>       lagrtyp(ilagr) % veloc(3)  = vx,vy,vz ... Old velocity
!>       lagrtyp(ilagr) % t         = t .......... Current time
!>
!> @}
!------------------------------------------------------------------------

subroutine pts_solite()
  use def_master,         only : vorti
  use def_master,         only : momod,mem_modul,modul,igene
  use def_master,         only : cutim,dtime,dtinv,ittim
  use def_master,         only : nparr,tempe,zeror,pard1
  use def_master,         only : advec,gisca,CPU_ASSEMBLY
  use def_master,         only : npasr
  use def_master,         only : ISLAVE,IMASTER,IPARALL,kfl_paral
  use def_master,         only : INOTMASTER,INOTSLAVE
  use def_master,         only : momen,kfl_coupl
  use def_master,         only : lelbf,leinv_loc,cpu_modul
  use def_master,         only : ID_NASTIN,ID_PARTIS
  use def_parame,         only : pi
  use def_kermod,         only : gravi,grnor,kfl_detection,relse
  use def_domain,         only : lnods,coord,ltype,nnode,npoin
  use def_domain,         only : lelel_2,pelel_2,nelem,leldo
  use def_domain,         only : ndime,mnode,walln
  use def_domain,         only : element_bin
  use def_domain,         only : element_bin_boxes
  use def_domain,         only : walld
  use def_elmtyp,         only : TET04,TRI03
  use mod_ker_proper,     only : ker_proper
  use mod_random,         only : random_grnd
  use mod_elmgeo,         only : elmgeo_natural_coordinates
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_memory,         only : memory_resize
  use mod_maths,          only : maths_mapping_coord_to_3d
  use mod_parall,         only : PAR_CODE_SIZE
  use mod_communications, only : PAR_ALLGATHER
  use mod_communications, only : PAR_SEND_RECEIVE
  use mod_communications, only : PAR_DEFINE_COMMUNICATOR
  use mod_communications, only : PAR_COMM_RANK_AND_SIZE
  use mod_communications, only : PAR_SUM
  use mod_communications, only : PAR_ALLGATHER
  use mod_communications, only : PAR_ALLGATHERV
  use mod_communications, only : PAR_MAX
  use mod_communications, only : PAR_SEND_RECEIVE_TO_ALL
  use mod_ker_detection,  only : ker_events_particle_not_converged
  use mod_maths,          only : maths_local_orthonormal_basis
  use mod_maths,          only : maths_vector_to_new_basis
  use mod_maths,          only : maths_vector_from_new_basis
  use mod_ker_timeline,   only : ker_timeline

  use def_partis
  use mod_pts_injection,  only : pts_injection_injectors
  
  use mod_messages,       only : livinf
#ifdef _OPENMP
  use omp_lib
#endif
  implicit none
  integer(ip)               :: ielem,ielel,jelem,pelty,pnode,inode,ithr
  integer(ip)               :: ipoin,idime,jdime,ifoun,ineig,dummi(3),ipars
  integer(ip)               :: mcros,mcror,nneig,icror,ilagr,iparr,new_size
  integer(ip)               :: itype,iprop,nvar1,comcont,ilagr_local
  integer(ip)               :: nlagr_free,nlagr_new,nlagr_last,iboun
  integer(ip)               :: melel,dumm0,nlagr_local,ii,jj,kk
  integer(ip),  allocatable :: ncros(:)
  integer(ip),  allocatable :: ncror(:)
  integer(ip),  pointer     :: permu_nlagr(:)
  real(rp),     pointer     :: parrs_pts(:)
  real(rp),     pointer     :: parre_pts(:)
  integer(ip),  pointer     :: pari1_pts(:)
  integer(ip),  pointer     :: paris_pts(:)
  integer(ip),  pointer     :: parig_pts(:)
  real(rp)                  :: elcod(ndime,mnode),toler,U1,U2
  real(rp)                  :: coloc(3),deriv(ndime,mnode),shapf(mnode)
  real(rp)                  :: dummr(2),coord_kp1(3),xxd(3),eps_br(3)
  real(rp)                  :: t,tf,dtc,hleng,venor,xfact,dt_elm
  real(rp)                  :: grafo,buofo,Re,dtg,mass
  real(rp)                  :: beta2,strex,ovstr,ti,itmaxr
  real(rp)                  :: accel_kp1(3),veloc_kp1(3)                ! a^{k+1},u^{k+1},x^{k+1}
  real(rp)                  :: vefl1(3),vefl2(3)                        ! Fluid velocity at n,n^+1
  real(rp)                  :: v_fluid_k(3),dt_k,alpha_k
  real(rp)                  :: v_fluid_km1(3),dt_km1,alpha_km1
  real(rp)                  :: v_fluid_km2(3),dt_km2,alpha_km2
  real(rp)                  :: dt012,dt12,dt01

  integer(ip)               :: iwall,itint,itmax0,itmax
  real(rp)                  :: xinte(3),dista
  
  real(rp)                  :: visfl,denfl,diame,Du,CdRe,Cd             ! Drag
  real(rp)                  :: spher,time1,time2,time_total,time0
  real(rp)                  :: time_max,time_ave,load_balance,tau_p
  real(rp)                  :: g(3),denpa,Cc,lambda

  integer(ip)               :: iiter,niter                              ! Newton Raphson for drag
  real(rp)                  :: xerro,xdeno
  real(rp)                  :: xnume
  real(rp)                  :: ff(3),df(3),dRedu
  real(rp)                  :: dCddRe,deltu(3),veloi(3)

  logical(lg)               :: inscont
  integer(ip),  save        :: ipass = 0

  real(rp)                  :: eps1,eps2,eps3,D,mu,r                   ! Random walk

  real(rp)                  :: alpha_str,h,uu,alpha
  logical(lg)               :: accept_time_step

  real(rp)                  :: tau,nu,tauinv,Stk(2),uf			! Stk(1): instantaneous Stk, Stk(2): effective Stk

  real(rp)                  :: eps(ndime,ndime)                         ! Saffman
  real(rp)                  :: K                                        ! Saffman constant coefficient
  real(rp)                  :: saff_deno                                ! Denominator in Saffman force
  real(rp)                  :: wf,ur,beta,funcRe,Cls,Res                ! Saffman Mei 
  real(rp)                  :: urel(ndime),vorti_fl(ndime  )            ! relative velocity, voriticity

  logical(lg)               :: newmark_converged


  logical(lg)               :: local_axes
  real(rp)                  :: basis(ndime,ndime),bnorm

  logical(lg)               :: debugmode
  integer(4)                :: PAR_COMM_TO_USE4
  type(comm_data_par), pointer :: commu
  !
  ! Detection
  !
  integer(4),   pointer     :: par_nlagr_4(:)
  integer(4)                :: my_nlagr_4,ipart4,ndime4,my_rank4
  integer(ip)               :: klagr
  real(rp),     pointer     :: coord_nlagr_4(:,:)
  real(rp),     pointer     :: my_coord_nlagr_4(:,:)
  real(rp),     pointer     :: my_veloc_nlagr_4(:,:)
  !
  ! Automatic tests
  !
  logical(lg)               :: test_brownian
  
#ifdef EVENT
  call mpitrace_user_function(1)
#endif

  test_brownian = .false.
  nullify(commu)
  !
  ! Define communication array
  !
  if( IPARALL ) call PAR_DEFINE_COMMUNICATOR('IN MY ZONE',PAR_COMM_TO_USE4,commu)

  !----------------------------------------------------------------------
  !
  ! USEFUL FOR DEBUGGING IN PARALLEL:
  !
  ! Reorder elements linked list to go through the neighboring elements
  ! in the same order in sequential and parallel
  !
  !----------------------------------------------------------------------
  debugmode = .false.
  if( debugmode .and. INOTMASTER ) then
     do ielem = 1,nelem
        dumm0 = pelel_2(ielem+1) - pelel_2(ielem)
        call memgen(1_ip,dumm0,0_ip)
        inode = 0
        do ielel = pelel_2(ielem),pelel_2(ielem+1)-1
           jelem = lelel_2(ielel)
           inode = inode + 1
           gisca(inode) = leinv_loc(jelem)
        end do
        call heapsorti2(1_ip,inode,gisca,lelel_2(pelel_2(ielem)))
        call memgen(3_ip,dumm0,0_ip)
     end do
  end if
  
  !----------------------------------------------------------------------
  !
  ! Numerical constants and definitions
  !
  !----------------------------------------------------------------------

  ipass     =  ipass + 1
  toler     =  relse(1)                 ! Tolerance for element search
  nvar1     =  33                       ! Number of variables required when sending one particle to another subd.
  if( kfl_posla_pts == 4 ) then         ! Add number of property variables managed by modules
     nvar1 = nvar1 + nlapr
  end if

  !----------------------------------------------------------------------
  !
  ! Inject particles
  !
  !----------------------------------------------------------------------

  nlagr = nlacc_pts

  call pts_injection_injectors()
  
  !----------------------------------------------------------------------
  !
  ! Global number for particles: NLAGR, starting from NLACC_PTS
  !
  !----------------------------------------------------------------------
  !
  ! Allocate memory to communicate with my neighbors
  !
  call Parall(704_ip)
  nneig = pard1
  nullify(parrs_pts)
  nullify(parre_pts)
  nullify(pari1_pts)
  nullify(paris_pts)
  nullify(parig_pts)
  nullify(permu_nlagr)
  nullify(par_nlagr_4)
  nullify(coord_nlagr_4)
  nullify(my_coord_nlagr_4)
  nullify(my_veloc_nlagr_4)
  !
  ! Temporary array to send/receive to my neighbors
  !
  !nthr    =  omp_get_num_threads()

  if( ISLAVE ) then
     allocate(ncror(nneig))
     allocate(ncros(nneig))
  else
     allocate(ncror(1))
     allocate(ncros(1))
     ncror = 0
     ncros = 0
  end if
  !----------------------------------------------------------------------
  !
  ! Initialization
  !
  !----------------------------------------------------------------------

  time_total     = 0.0_rp             ! Total CPU time
  particles_sent = 0                  ! # of particles received
  particles_recv = 0                  ! # of particles sent
  comm_loops_pts = 0                  ! # communication loops
  mcros          = 0                  ! Size of receiving array (Parall)
  mcror          = 0                  ! Size of sending array (Parall)
  comcont        = 0                  ! Communication loop active
  nlagr_existing_pts        = 0                  ! Number of existing particles
  nlagr_local    = 0
  ti             = cutim - dtime      ! Initial time
  tf             = cutim              ! Final time: t^{n+1}
  ovstr          = 1.2_rp             ! Inverse stretching factor
  strex          = 1.0_rp / ovstr     ! Stretching factor
  g(1)           = grnor * gravi(1)   ! Gravity gx
  g(2)           = grnor * gravi(2)   ! Gravity gy
  g(3)           = grnor * gravi(3)   ! Gravity gz
  K              = 2.594_rp           ! Constant Saffman coefficient
  newmark_converged = .true.          ! Flag to set if newmark did not converge  

  !----------------------------------------------------------------------
  !
  ! NLAGR_EXISTING_PTS= Number of existing particles (should be equal to nlacc_pts)
  !
  !----------------------------------------------------------------------

  do ilagr = 1,mlagr
     !do ilagr_local = 1,nlagr_local
     !ilagr = permu_nlagr(ilagr_local)
     if( lagrtyp(ilagr) % kfl_exist == -1 ) then
        nlagr_existing_pts = nlagr_existing_pts + 1
        lagrtyp(ilagr) % t  = ti
     end if
  end do
  if( INOTMASTER ) then
     call memory_alloca(mem_modul(1:2,modul),'PERMU_NLAGR','pts_solite',permu_nlagr,mlagr)
  end if

  call PAR_SUM(nlagr_existing_pts,'IN MY ZONE')

  call livinf(-9_ip,'TRANSPORT LAGRANGIAN PARTICLES= ',nlagr_existing_pts)

  call ker_timeline('INI_SOLVER',nlagr_existing_pts)
  
  !----------------------------------------------------------------------
  !
  ! Loop over communication iterations
  !
  !----------------------------------------------------------------------
  do while( comcont /= -1 )
     !
     ! Communcation variables
     !
     call cputim(time1)
     comcont = comcont + 1
     comm_loops_pts = comm_loops_pts + 1
     nlagr_non_migrating_pts = 0
     if( ISLAVE ) then
        mcros = 0
        mcror = 0
        do ineig = 1,nneig
           ncros(ineig) = 0
           ncror(ineig) = 0
        end do
     end if
     !
     ! List of existing particles
     !
     if( INOTMASTER ) then
        nlagr_local = 0
        nlagr_free  = 0
        do ilagr = 1,mlagr
           if(      lagrtyp(ilagr) % kfl_exist ==  0 ) then
              nlagr_free = nlagr_free + 1
           else if( lagrtyp(ilagr) % kfl_exist == -1 ) then
              nlagr_local = nlagr_local + 1
              permu_nlagr(nlagr_local) = ilagr
           else if( lagrtyp(ilagr) % kfl_exist <= -2 ) then
              nlagr_non_migrating_pts = nlagr_non_migrating_pts + 1
           end if
        end do
     end if
     !
     ! Loop over existing particles
     !


     !----------------------------------------------------------------------
     !
     ! Loop over particles
     !
     !----------------------------------------------------------------------

     !$OMP PARALLEL DO SCHEDULE (DYNAMIC,1000)                                                                &
     !$OMP DEFAULT      (NONE)                                                                                &
     !$OMP PRIVATE      (itype,dt_k,alpha_str,itint,Cc,dtg,ielem,iboun,inode,ipoin,pelty,pnode,t,             &
     !$OMP              xfact,dt_km1,dt_km2,veloc_kp1,coord_kp1,accel_kp1,elcod,beta2,alpha_k,                &
     !$OMP              alpha_km1,dt012,dt01,dt12,alpha_km2,visfl,denfl,denpa,spher,grafo,buofo,nu,eps,eps1,  &
     !$OMP              eps2,eps3,veloi,xerro,niter,iiter,Du,df,alpha,tauinv,tau,dRedu,saff_deno,             &
     !$OMP              xdeno,xnume,ff,deltu,accept_time_step,uu,vv,h,D,mu,venor,dtc,iwall,                   &
     !$OMP              vefl1,vefl2,v_fluid_km1,v_fluid_km2,v_fluid_k,diame,ifoun,Stk,uf,                     &
     !$OMP              dista,melel,inscont,jelem,ineig,ielel,ilagr,hleng,r,tau_p,itmaxr,itmax0,              &
     !$OMP              ilagr_local,dt_elm,idime,jdime,lambda,CdRe,Re,Cd,dCddRe,coloc,deriv,shapf,            &
     !$OMP              dummr,ii,jj,kk,xxd,xinte,local_axes,bnorm,basis,mass,time1,                           &
     !$OMP              Res,wf,ur,beta,funcRe,Cls,urel,vorti,vorti_fl,newmark_converged)                                        &
     !$OMP SHARED       (nlagr_local,parttyp,mean_free_path_pts,dtmin_pts,dtime,ltype,nnode,                  &
     !$OMP              lnods,coord,toler,lagrtyp,lelbf,lelel_2,element_bin,itmax,                            &
     !$OMP              ti,tf,defor_pts,gamma_pts,K,g,strex,beta_pts,dtime_pts,ovstr,tempe,                   &
     !$OMP              lboue_pts,dimin_pts,walln,pelel_2,kfl_usbin_pts,kfl_resid_pts,                        &
     !$OMP              nelem,kfl_depos_pts,leinv_loc,advec,kfl_exacs_pts,leldo,dtinv,kfl_walld_pts,          &
     !$OMP              permu_nlagr,hleng_pts,element_bin_boxes,kfl_fixbo_pts,dumm0,kfl_order_pts,            &
#ifndef NDIMEPAR
     !$OMP              ndime,                                                                                &
#endif
     !$OMP              bvnat_pts,leleboun_pts,bouno_pts,kfl_slip_wall_pts,walld_slip_pts,resid_pts,momen,    &
     !$OMP              kfl_coupl,walld,kfl_paral)                                                                            &
     !$OMP REDUCTION    (+:nlagr_non_migrating_pts,ncros)


     do ilagr_local = 1,nlagr_local
        ilagr = permu_nlagr(ilagr_local)
        iwall = 0 !reset here because for some particles the other iwall=0 never gets executed
        !
        ! Update solution at k+1
        !
        ! k-2          k-1            k            k+1
        !  o-------------o-------------o-------------o------>
        !      dt^k-2         dt^k-1         dt^k
        !
        inscont   = .true.
        itype     = lagrtyp(ilagr) % itype      ! Particle type
        t         = lagrtyp(ilagr) % t          ! Particle time
        dt_k      = lagrtyp(ilagr) % dt_k       ! Particle time step guess: t^k+1 - t^k
        alpha_str = lagrtyp(ilagr) % stret      ! Stretching
        itint     = 0                           ! Total number time step number (including non-accepted ones)
        diame     = parttyp(itype) % diame      ! Particle diameter
        r         = 0.5_rp * diame              ! Particle radius
        denpa     = parttyp(itype) % denpa      ! Particle density
        mass      = 4.0_rp / 3.0_rp * pi * r**3 ! Particle mass
        Cd        = 0.0_rp                      ! Drag Coefficient inicalization
        Re        = 0.0_rp                      ! Reynold's Particle initialization
        v_fluid_k = 0.0_rp
        itmax     = 200 !!!!! OJO 100000
        !
        ! Cunningham slip correction factor
        ! Some values of airborne particles at standrad conditions (T=293K,P=101kPa)
        !
        ! Particle diameter   Slip Correction
        ! d(m)                factor Cc
        ! -----------------------------------
        ! 10^-9               224.332
        ! 10^-8                22.976
        ! 10^-7                 2.928
        ! 10^-6                 1.155
        ! 10^-5                 1.015
        ! 10^-4                 1.002
        !
        if( mean_free_path_pts > 0.0_rp ) then
           lambda = mean_free_path_pts
           Cc     = 1.0_rp+2.0_rp*lambda/diame*(1.257_rp+0.4_rp*exp(-0.55_rp*diame/lambda))
        else
           Cc     = 1.0_rp
        end if
        !-------------------------------------------------------------------
        !
        ! Loop over time
        !
        ! ti               t <--dtk-->    tf
        ! o----------------|----------|---o
        ! n                              n+1
        ! <------------------------------->
        !              dtime
        !
        !-------------------------------------------------------------------
        do while( t < tf-zeror .and. inscont )
           ielem = lagrtyp(ilagr) % ielem
           hleng = hleng_pts(ielem)
           !
           ! Modify time step
           !
           if( parttyp(itype) % kfl_tstep < 0 ) then
              dt_k = max(dtmin_pts,min(dt_k,dtime_pts))
           else
              dt_k = max(dtmin_pts,dt_k*alpha_str)
           end if
           !dt_elm = (hleng / (sqrt(v_fluid_k(1)**2 + v_fluid_k(2)**2 + v_fluid_k(3)**2)+zeror))*0.5_rp
           !dt_k   = min(dt_k,dt_elm)
           dt_k   = min(dtime,dt_k)
!print*,'-------------dt_k---------',dt_k
           dtg    = dt_k
           !
           ! Time and predicted time step
           !
           t      = lagrtyp(ilagr) % t
           !
           ! Synchronize with tf when dt is too large
           !
           if( t + dt_k >= tf-zeror ) then
              dt_k = tf - t

              !if(t>8.0e-2_rp) 
!print*,'mierda=',t,dt_k,tf,tf-zeror
              !if(t>8.0e-2_rp)print*,'mierda'
              !call flush(200)

           end if
           !
           ! Time step is too small, consider we have arrived
           !
           if( dt_k < zeror ) then
              inscont = .false.
              t       = tf
!if(t>8.0e-2_rp)print*,'time step too small'
              goto 20
           end if
           !
           ! Advance in time
           !
           t     = t + dt_k
           itint = itint + 1
! if(t>8.0e-2_rp)print*,'b=',t,dt_k,tf,tf-zeror

           !
           ! Particle is in element IELEM. Get value of shape function SHAPF in it
           ! Compute element length HLENG
           !
           pelty = ltype(ielem)
           pnode = nnode(pelty)
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              elcod(1:ndime,inode) = coord(1:ndime,ipoin)
           end do
!if(t>8.0e-2_rp)then
!print*,'cacaaaaaaaa',lagrtyp(ilagr) % coord,ielem,toler
!print*,'elcod1',elcod(1:3,1)
!print*,'elcod2',elcod(1:3,2)
!print*,'elcod3',elcod(1:3,3)
!print*,'elcod4',elcod(1:3,4)
!print*,'elcod5',elcod(1:3,5)
!print*,'elcod6',elcod(1:3,6)
!end if
           call elmgeo_natural_coordinates(             &
                ndime,pelty,pnode,elcod,shapf,deriv,    &
                lagrtyp(ilagr) % coord,coloc,ifoun,toler)
!if(t>8.0e-2_rp)print*,'salgo con ifoun again',ifoun
           !
           ! A particle did not find its element! This is very strange
           !
           if( ifoun == 0 ) then
              lagrtyp(ilagr) % kfl_exist = -3
              inscont = .false.
              t       = tf
              print *, 'particle did not find its element a=',lagrtyp(ilagr) % ilagr
              print *, 'particle did not find its element b=',lagrtyp(ilagr) % coord(1:ndime)
              print *, 'particle did not find its element c=',leinv_loc(ielem),pelty,pnode
              !do inode = 1,pnode
              !   write(6,'(A,i2,10(1x,e12.6))') 'particle did not find its element d=',inode,elcod(1:ndime,inode)
              !end do
              !call runend('A PARTICLE DID NOT FIND ITS ELEMENT')
              goto 20
           end if
           !
           ! Interpolate fluid velocity uf(t^k,x^k)
           !
           ! Fluid velocity at new particle position
           !
           ! VEFL1 = u^n+1 at tf
           ! VEFL2 = u^n   at ti
           !
           ! ti              t^k            tf
           ! o----------------o-------------o
           ! VEFL2           x^k           VEFL1
           !
           !
           vefl1(1:ndime) = 0.0_rp
           vefl2(1:ndime) = 0.0_rp
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              vefl1(1:ndime) = vefl1(1:ndime) + shapf(inode) * advec(1:ndime,ipoin,1)
              vefl2(1:ndime) = vefl2(1:ndime) + shapf(inode) * advec(1:ndime,ipoin,3)
           end do
           xfact = dtinv * ( lagrtyp(ilagr) % t - ti )
           lagrtyp(ilagr) % v_fluid_k(1:ndime) = (1.0_rp-xfact) * vefl2(1:ndime) + xfact * vefl1(1:ndime)
           !
           !  Fluid velocity at previous times, a priori unknown at k+1
           !
           !   n     k-2     k-1     k     k+1     n+1
           !   o------|-------|------|------|-------o
           !   ti      dt_km2  dt_km1  dt_k         tf
           !
           !   uf^k   = uf(x^k,t^k)     = v_fluid_k(1:3)   = lagrtyp(ilagr) % v_fluid_k(1:3)
           !   uf^k-1 = uf(x^k-1,t^k-1) = v_fluid_km1(1:3) = lagrtyp(ilagr) % v_fluid_km1(1:3)
           !   uf^k-2 = uf(x^k-2,t^k-2) = v_fluid_km2(1:3) = lagrtyp(ilagr) % v_fluid_km2(1:3)
           !
           !   dt^k   = t^k+1 - t^k
           !   dt^k-1 = t^k   - t^k-1
           !   dt^k-2 = t^k-1 - t^k-2
           !
           v_fluid_k(1:3)   = lagrtyp(ilagr) % v_fluid_k(1:3)
           v_fluid_km1(1:3) = lagrtyp(ilagr) % v_fluid_km1(1:3)
           v_fluid_km2(1:3) = lagrtyp(ilagr) % v_fluid_km2(1:3)
           dt_km1           = lagrtyp(ilagr) % dt_km1
           dt_km2           = lagrtyp(ilagr) % dt_km2
           !
           ! Fluid properties
           !
           call ker_proper('VISCO','IGAUS',dumm0,ielem,dummr,pnode,1_ip,shapf)   ! Fluid viscosity
           visfl = dummr(1)
           call ker_proper('DENSI','IGAUS',dumm0,ielem,dummr,pnode,1_ip,shapf)   ! Fluid density
           denfl = dummr(1)
           !
           ! Difussion Coefficeint values
           !
           if( parttyp(itype) % kfl_brown /= 0 ) then
              D = parttyp(itype) % diffu
              call ker_proper('VISCO','IGAUS',dumm0,ielem,dummr,pnode,1_ip,shapf)
              mu = dummr(1)
              call pts_brown(itint,ielem,itype,pnode,shapf,            &
                   &              dt_k,r,mu,Cc,D,eps1,eps2,eps3)
           else
              D = 0.0_rp
           end if

           !----------------------------------------------------------
           !
           !     ----------------------------------------------
           !     UPDATE ACCELERATION, VELOCITY, POSITION AT K+1
           !     ----------------------------------------------
           !
           !----------------------------------------------------------

           if( kfl_exacs_pts /= 0 ) then
              !if( kfl_exacs_pts == 0 ) then
              ! Exact acceleration
              call pts_exacso(3_ip,t,accel_kp1,dummr,dummr)
              do idime = 1,ndime
                 veloc_kp1(idime) = lagrtyp(ilagr) % veloc(idime) &
                      + dt_k * ( (1.0_rp-gamma_pts) * accel_kp1(idime) + gamma_pts * lagrtyp(ilagr) % accel(idime))
                 coord_kp1(idime) = lagrtyp(ilagr) % coord(idime) + dt_k * lagrtyp(ilagr) % veloc(idime) &
                      + dt_k * dt_k * ( 0.5_rp*lagrtyp(ilagr) % accel(idime) &
                      + beta_pts *(accel_kp1(idime)-lagrtyp(ilagr) % accel(idime)) )
              end do

           else if( parttyp(itype) % kfl_modla == 1 ) then
              ! Analytic solution
              if ( 0 == 3 ) then
                 call pts_analytics(ielem, dt_k, lagrtyp(ilagr) % coord, coord_kp1)
              end if
              ! Pure transport: second order
              ! this method is known as a two-step method. More precisely, is is known as the second-order
              ! Adams-Bashforth method (or AB method) dating back to 1883
              ! If dt is constant, it gives:
              ! x^n+1 = x^n + 3/2*dt*u^n -1/2*dt*u^n-1
              !
              !     k-1     k      k+1
              !  o---|------|-------|--------o
              !       <----> <------>
              !        dto      dt
              !
              ! From Taylor series:
              ! (1)  x^{k+1} = x^k + dt * uk + 1/2 a^k * dt^2 + O(dt^3)
              ! (2)  u^{k-1} = u^k - a^k dto + O(dto^2)
              ! where u^k = u(t^k,x^k)
              !
              ! From (2) we have:
              ! a^k = ( u^k - u ^{k-1} ) / dto + O(dto)
              ! Substitute this in (1):
              ! x^{k+1} = x^k + dt * uk + 1/2 ( u^k - u^{k-1} ) / dto * dt^2 + O(dto*dt^2) + O(dt^3)
              ! x^{k+1} = x^k + ( dt + 1/2 dt^2 / dto ) u^k - 1/2 dt^2 / dto u^{k-1}
              !
              ! veloc and accel are always one iteration before, at k, for implementation purpose
              !
              if( kfl_order_pts < 0 .and. ( pelty == TET04 .or. pelty == TRI03 ) ) then
                 !
                 ! Analytical integration
                 !
                 call pts_analytics(ielem, dt_k, lagrtyp(ilagr) % coord, coord_kp1)

              else if( abs(kfl_order_pts) == 2 ) then
                 !
                 ! Second order AB
                 !
                 beta2              =  dt_k / dt_km1
                 alpha_k            =  dt_k * ( 1.0_rp + 0.5_rp * beta2 )
                 alpha_km1          = -dt_k * 0.5_rp * beta2
                 veloc_kp1(1:ndime) =  v_fluid_k(1:ndime)
                 accel_kp1(1:ndime) =  ( v_fluid_k(1:ndime) - v_fluid_km1(1:ndime) ) / dt_km1
                 coord_kp1(1:ndime) =  lagrtyp(ilagr) % coord(1:ndime)   &
                      &                + alpha_k   * v_fluid_k(1:ndime)  &
                      &                + alpha_km1 * v_fluid_km1(1:ndime)

              else if( abs(kfl_order_pts) == 3 ) then
                 !
                 ! Third order AB
                 !
                 dt012              = dt_km2 + dt_km1 + dt_k
                 dt01               = dt_km1 + dt_k
                 dt12               = dt_km2 + dt_km1
                 alpha_k            = ( 0.25_rp*(dt012*dt01*(dt012+dt01)-dt12*dt_km1*(dt12+dt_km1))&
                      &               + 1.0_rp/12.0_rp*(-dt012**3-dt01**3+dt12**3+dt_km1**3))/( dt_km1*dt12)
                 alpha_km1          = ( 0.25_rp*(dt012*dt_k *(dt012+dt_k))&
                      &               + 1.0_rp/12.0_rp*(-dt012**3-dt_k**3 +dt12**3))/(-dt_km2*dt_km1)
                 alpha_km2          = ( 0.25_rp*(dt01 *dt_k *(dt01 +dt_k))&
                      &               + 1.0_rp/12.0_rp*(-dt01**3 -dt_k**3 +dt_km1**3))/( dt_km2*dt12)
                 veloc_kp1(1:ndime) = v_fluid_k(1:ndime)
                 accel_kp1(1:ndime) = ( v_fluid_k(1:ndime) - v_fluid_km1(1:ndime) ) / dt_km1
                 coord_kp1(1:ndime) = lagrtyp(ilagr) % coord(1:ndime) + alpha_k   * v_fluid_k(1:ndime)   &
                      &                                               + alpha_km1 * v_fluid_km1(1:ndime) &
                      &                                               + alpha_km2 * v_fluid_km2(1:ndime)

              else
                 call runend('PTS_SOLITE: WRONG ORDER')
              end if

           else if( parttyp(itype) % kfl_modla == 2) then
              !
              ! Force model
              !
              ! Properties taken from previous time step at old position
              !
              spher      = parttyp(itype) % spher                                     ! Particle sphericity
              grafo      = real( parttyp(itype) % kfl_grafo, rp )                     ! Gravity  force = 1.0
              buofo      = real( parttyp(itype) % kfl_buofo, rp )                     ! Buoyancy force = 1.0
              nu         = visfl / denfl                                              ! Fluid kinematic viscosity
              local_axes = .false.                                                    ! If particle equation is in a local basis
              !
              ! Check if equations should be solved in local axes ! NEW
              !
              if( kfl_slip_wall_pts > 0 ) then
                 dista      = dot_product(shapf(1:pnode),walld_slip_pts(lnods(1:pnode,ielem)))
                 if( dista <= 1.0_rp*diame .or. lboue_pts(ielem) == PTS_SLIP_CONDITION ) then
                    local_axes = .true.
                    bnorm = 0.0_rp
                    do idime = 1,ndime
                       basis(idime,1) = dot_product(shapf(1:pnode),walln(idime,lnods(1:pnode,ielem)))
                       bnorm          = bnorm + basis(idime,1) * basis(idime,1)
                    end do
                    basis(1:ndime,1) = basis(1:ndime,1) / sqrt(bnorm)
                    call maths_local_orthonormal_basis(ndime,basis)
                 else
                    basis(1:ndime,1) = 0.0_rp
                 end if
              end if
              !
              ! Drag force, initial guess
              !
              veloc_kp1(1:ndime) = lagrtyp(ilagr) % veloc(1:ndime)
!write(94,*) t-dt_k, sqrt(veloc_kp1(1)**2+ veloc_kp1(2)**2+ veloc_kp1(3)**2)
              !
              ! Velocity deformation useful for Saffman's force
              !
              if( parttyp(itype) % kfl_saffm /= 0 ) then
                 !eps = 0.0_rp
                 !if (ndime == 2) then
                 !   do inode = 1,pnode
                 !      ipoin = lnods(inode,ielem)
                 !      eps(1,1) = eps(1,1) + 0.5_rp * shapf(inode) * defor_pts(1,ipoin)
                 !      eps(2,2) = eps(2,2) + 0.5_rp * shapf(inode) * defor_pts(2,ipoin)
                 !      eps(1,2) = eps(1,2) + 0.5_rp * shapf(inode) * defor_pts(3,ipoin)
                 !   end do
                 !   eps(2,1) = eps(1,2)
                 !else if(ndime == 3) then
                 !   do inode = 1,pnode
                 !      ipoin = lnods(inode,ielem)
                 !      eps(1,1) = eps(1,1) + 0.5_rp * shapf(inode) * defor_pts(1,ipoin)
                 !      eps(2,2) = eps(2,2) + 0.5_rp * shapf(inode) * defor_pts(2,ipoin)
                 !      eps(1,2) = eps(1,2) + 0.5_rp * shapf(inode) * defor_pts(3,ipoin)
                 !      eps(3,3) = eps(3,3) + 0.5_rp * shapf(inode) * defor_pts(4,ipoin)
                 !      eps(3,1) = eps(3,2) + 0.5_rp * shapf(inode) * defor_pts(5,ipoin)
                 !      eps(3,2) = eps(3,2) + 0.5_rp * shapf(inode) * defor_pts(6,ipoin)
                 !   end do
                 !   eps(2,1) = eps(1,2)
                 !   eps(2,3) = eps(3,2)
                 !   eps(1,3) = eps(3,1)

                    !
                    ! Saffman Mei
                    !
                    vorti_fl(1:ndime) = 0.0_rp
                    do inode = 1,pnode
                       ipoin = lnods(inode,ielem)
                       do idime = 1,ndime
                          vorti_fl(idime) = vorti_fl(idime) + shapf(inode) * vorti(idime,ipoin)
                       end do
                    end do
                 !end if
              end if

              if( 0 == 1 ) then
                 !
                 ! Runge-Kutta 4th-order
                 !
                 veloi = veloc_kp1
                 call pts_rungk4(                                                          &
                      ndime,parttyp(itype) % kfl_drafo,parttyp(itype) % kfl_extfo,         &
                      grafo,buofo,g,v_fluid_k,visfl,denfl,veloi,lagrtyp(ilagr) % coord,    &
                      denpa,diame,spher,t,dt_k,veloc_kp1,coord_kp1)
                 
              else if( 1 == 2 ) then
                 !
                 ! Runge-Kutta 4th/5th-order with time step control
                 !
                 veloi = veloc_kp1
                 call pts_rungk45(                                                          &
                      ndime,parttyp(itype) % kfl_drafo,parttyp(itype) % kfl_extfo,          &
                      grafo,buofo,g,v_fluid_k,visfl,denfl,veloi,lagrtyp(ilagr) % coord,     &
                      denpa,diame,spher,t,dt_k,veloc_kp1,xxd)

              else

                 if( parttyp(itype) % kfl_schem == 1 .and. parttyp(itype) % kfl_drafo /= 0 .and. 1 == 0) then

                    ! Analytical time integration
                    ! Only exact for Stokes approximation, Re <= 1, for which tau = cst
                    ! du/dt  = 1/tau*(uf-u) + a_ext; let u* = u - uf
                    ! du*/dt = -u*/tau - duf/dt + a_ext
                    !        = -u*/tau + a_total; a_total is rewritten so that we recover the pure transport
                    !           equation when there is no drag
                    !
!!$                       call dragfo(parttyp(itype) % kfl_drafo,Du,visfl,denfl,diame,spher,CdRe,Cd,Re,dCddRe)
!!$                       tau      =  denpa * diame * diame * Cc / ( 0.75_rp * visfl * CdRe )
!!$                       if( dt_k/tau > 100.0_rp ) then
!!$                          expdttau =  exp(-100.0_rp)
!!$                       else
!!$                          expdttau =  exp(-dt_k/tau)
!!$                       end if
!!$                       beta2    =  dt / dt_km1
!!$                       alph3    =  dt * ( 1.0_rp + 0.5_rp * beta2 )
!!$                       alph4    = -dt * 0.5_rp * beta2
!!$                       do idime = 1,ndime
!!$                          lagrtyp(ilagr) % acceg(idime) = g(idime) * ( grafo - buofo * denfl / denpa )
!!$                          accel_kp1(idime)                  = lagrtyp(ilagr) % acceg(idime)
!!$                          veloc_kp1(idime) =  lagrtyp(ilagr) % veloc(idime) * expdttau &
!!$                               +          ( lagrtyp(ilagr) % v_fluid_k(idime) + accel_kp1(idime) * tau ) * ( 1.0_rp-expdttau )
!!$                          accel_kp1(idime) =  1.0_rp / tau * ( lagrtyp(ilagr) % v_fluid_k(idime) - veloc_kp1(idime) ) + lagrtyp(ilagr) % acceg(idime)
!!$                          coord_kp1(idime)    =  coord_k(idime) + tau * ( lagrtyp(ilagr) % veloc(idime) &
!!$                               &          - lagrtyp(ilagr) % v_fluid_k(idime) - accel_kp1(idime) * tau ) * ( 1.0_rp - expdttau ) &
!!$                               &          + accel_kp1(idime) * tau * dt &
!!$                               &          + alph3 * lagrtyp(ilagr) % v_fluid_k(idime) + alph4 * lagrtyp(ilagr) % v_fluid_k(idime)
!!$
!!$                       end do

                 else
                    !
                    ! Numerical time integration: Newmark + Newton-Raphson
                    !
                    xerro = 1.0_rp
                    niter = 100
                    iiter = 0
                    if( parttyp(itype) % kfl_drafo + parttyp(itype) % kfl_extfo == 0 ) niter = 1
                    do while( iiter < niter .and. xerro > 1.0e-12_rp )
                       iiter  = iiter + 1
                       Du     = 0.0_rp
                       do idime = 1,ndime
                          Du               =  Du + ( veloc_kp1(idime) - v_fluid_k(idime) ) ** 2
                          df(idime)        = -1.0_rp
                          accel_kp1(idime) =  0.0_rp
                       end do
                       Du = sqrt( Du ) + zeror                                             ! Relative velocity
                       !
                       ! Drag force:
                       ! a^{n+1} = 1/tau * ( u_fluid - u )
                       ! tau is referred to as relaxation time and tau = tau(u)
                       ! tau =  ( rho_p * d^2 ) / ( 3/4 * mu * Cd * Re )
                       !
                       if( parttyp(itype) % kfl_drafo /= 0 ) then
                          if(parttyp(itype) % kfl_drafo == 8 ) then
                              dista = dot_product(shapf(1:pnode),walld(lnods(1:pnode,ielem)))
                           else
                              dista =0.0_rp
                           end if
                          call dragfo(parttyp(itype) % kfl_drafo,Du,visfl,denfl,diame,dista,spher,CdRe,Cd,Re,dCddRe)
                          alpha  =  0.75_rp * visfl / ( denpa * diame * diame * Cc )
                          tauinv =  alpha * CdRe
                          tau    =  1.0_rp / ( tauinv + zeror )
                          uf     =  sqrt(dot_product(v_fluid_k(1:ndime),v_fluid_k(1:ndime)))
                          tau_p  = (denpa * diame * diame)/(18.0_rp*visfl) 
			  Stk    =  tau * uf / diame
			 do idime = 1,ndime
                             lagrtyp(ilagr) % acced(idime) = alpha * CdRe * ( v_fluid_k(idime) - veloc_kp1(idime) )
                             accel_kp1(idime)              = accel_kp1(idime) + lagrtyp(ilagr) % acced(idime)
                             dRedu                         = diame / nu * ( veloc_kp1(idime) - v_fluid_k(idime) ) / Du
                             df(idime)                     = df(idime) - dt_k * gamma_pts * alpha * ( CdRe + ( veloc_kp1(idime) - v_fluid_k(idime) ) * dRedu * dCddRe )
                          end do
                       else
                          Stk = 0.0_rp
                       end if
                       !print*,'caca=',du,cd,re,cdre,dot_product(( v_fluid_k(1:ndime) - veloc_kp1(1:3) ),( v_fluid_k(1:3) - veloc_kp1(1:3) )),dredu
                       !
                       ! Saffman
                       !
                       if( parttyp(itype) % kfl_saffm /= 0 ) then
                         !Saffman old
                         
!!$                          saff_deno = 0.0_rp
!!$                          do iime = 1,ndime
!!$                             do jdime = 1,ndime
!!$                                saff_deno = saff_deno + eps(idime,jdime)*eps(jdime,idime)
!!$                             end do
!!$                          end do
!!$                          saff_deno = (saff_deno) ** (0.25_rp)
!!$                          if(saff_deno /= 0.0_rp )then
!!$                             do idime = 1,ndime
!!$                                do jdime = 1,ndime
!!$                                   !if (abs( (veloi(jdime) - vv(jdime))) > zeror ) then
!!$                                   accel_kp1(idime) = accel_kp1(idime) + ((2.0_rp * K * sqrt(visfl/denfl) *  eps(idime,jdime) * &
!!$                                        denfl ) * (v_fluid_k(jdime) - veloi(jdime))) / (saff_deno * diame * denpa)
!!$                                   !end if
!!$                                end do
!!$                             end do
!!$                             do idime = 1,ndime
!!$                                df(idime) = df(idime) - (dt_k * gamma_pts * 2.0_rp * K *  sqrt(visfl/denfl) *  eps(idime,idime) * &
!!$                                     denfl ) / (saff_deno * diame *denpa)
!!$                             end do
!!$                          end if
                          
                          !
                          !Saffman Mei
                          !                          

                          wf            =  sqrt(dot_product(vorti_fl(1:ndime),vorti_fl(1:ndime)))
                          urel(1:ndime) = v_fluid_k(1:ndime)-veloi(1:ndime)
                          ur            =  sqrt(dot_product(urel(1:ndime),urel(1:ndime)))


                          call liftfo(parttyp(itype) % kfl_saffm,ur,wf,visfl,denfl,diame,dista,Cls,Re,Res) 

                          accel_kp1(1)  =  accel_kp1(1)  + denfl * pi/8.0_rp * diame * diame * diame * Cls * ( urel(2) * vorti_fl(3) - urel(3) * vorti_fl(2) )
                          accel_kp1(2)  =  accel_kp1(2)  + denfl * pi/8.0_rp * diame * diame * diame * Cls * ( urel(3) * vorti_fl(1) - urel(1) * vorti_fl(3) )
                          accel_kp1(3)  =  accel_kp1(3)  + denfl * pi/8.0_rp * diame * diame * diame * Cls * ( urel(1) * vorti_fl(2) - urel(2) * vorti_fl(1) )

                       end if
                       !
                       ! External force
                       !
                       if( parttyp(itype) % kfl_extfo /= 0 ) then
                          call extefo(&
                               parttyp(itype) % kfl_extfo,lagrtyp(ilagr) % coord,&
                               denpa,spher,denfl,visfl,t,lagrtyp(ilagr) % accee)
                          do idime = 1,ndime
                             accel_kp1(idime) = accel_kp1(idime) + lagrtyp(ilagr) % accee(idime)
                          end do

                       end if
                       !
                       ! Gravity and buoyancy
                       !
                       do idime = 1,ndime
                          lagrtyp(ilagr) % acceg(idime) = - g(idime) * buofo * denfl / denpa
                          accel_kp1(idime)              = accel_kp1(idime) + lagrtyp(ilagr) % acceg(idime)
                          accel_kp1(idime)              = accel_kp1(idime) + g(idime) * grafo
                       end do
                       !
                       ! Brownian force
                       !
                       if( parttyp(itype) % kfl_brown == 2 ) then
                          if( ndime == 2 ) then
                             accel_kp1(1) = accel_kp1(1) + eps1
                             accel_kp1(2) = accel_kp1(2) + eps2
                          else
                             accel_kp1(1) = accel_kp1(1) + eps1
                             accel_kp1(2) = accel_kp1(2) + eps2
                             accel_kp1(3) = accel_kp1(3) + eps3
                          end if
                       end if
                       !
                       ! Update velocity and position: given an+1
                       ! u^{n+1} = u^n + dt*( (1-gamma)*a^n + gamma*a^{n+1} )
                       ! x^{n+1} = x^n + u^n dt + dt^2 * [ 1/2 a^n + beta * ( a^{n+1} - a^n ) ]
                       !
                       ! u^{n+1} is computed using a Newton-Raphson method
                       !
                       ! f(u) = 0 => f^{k+1} = f^k + df/du Du.
                       ! Let f^{k+1} = 0, then
                       ! u^{n+1} = u^{n+1}(k) - f(u) / df(u)/du
                       ! where     f(u) = -u + u^n + dt*( (1-gamma)*a^n + gamma*a^{n+1} )
                       !       df(u)/du = -1 + dt*gamma*da^{n+1}/du
                       !
                       !print*,'b=',veloc_kp1(1),accel_kp1(1),veloc_kp1(1)+dt_k*accel_kp1(1),CdRe,Re
                       !stop
!print*,'avant deltu: ',deltu(3),df(3),ff(3)
                       do idime = 1,ndime
                          ff(idime)       =  ( -veloc_kp1(idime) + lagrtyp(ilagr) % veloc(idime) &
                               &             + dt_k * ( (1.0_rp-gamma_pts) * lagrtyp(ilagr) % accel(idime) + gamma_pts * accel_kp1(idime) ) )
                          deltu(idime)     = - ff(idime) / ( df(idime) + zeror )
                          veloc_kp1(idime) =   veloc_kp1(idime) + deltu(idime)
                       end do
!print*,'apres deltu: ',deltu(3),df(3),ff(3)
                       !
                       ! Local basis: rotate and cancel outflow normal velocity
                       !
                       if( local_axes ) then
                          call maths_vector_to_new_basis(ndime,basis,veloc_kp1)
                          deltu(1)     = min(veloc_kp1(1),0.0_rp) - veloc_kp1(1)
                          veloc_kp1(1) = min(veloc_kp1(1),0.0_rp)
                          call maths_vector_from_new_basis(ndime,basis,veloc_kp1)
                       end if
                       !
                       ! Residual
                       !
                       xnume = dot_product(deltu(1:ndime),deltu(1:ndime))
                       xdeno = dot_product(veloc_kp1(1:ndime),veloc_kp1(1:ndime))
                       xnume = sqrt(xnume)
                       xdeno = sqrt(xdeno) + zeror
                       xerro = xnume / xdeno
                       !print*,'NR: ',iiter,xnume,xdeno
                       !print*,'NR: ',iiter,ff(3),df(3),deltu(3),lagrtyp(ilagr) % coord(3)
     
                    end do
                    !
                    ! Update new position
                    !
                    if( iiter == niter .and. niter /= 1 ) then
                       alpha_str        =  strex
                       accept_time_step = .false.
                       t = lagrtyp(ilagr) % t
                       newmark_converged = .false.  !temporarily set to throw an error in 20
                       !print *, 'Newmark didnt converge',lagrtyp(ilagr) % ilagr,iiter,niter,' dt=',dt_k, ', v=',veloc_kp1(1:ndime),', |v|=', xdeno, kfl_paral
                       goto 20
                    else
                       do idime = 1,ndime
                          coord_kp1(idime) = lagrtyp(ilagr) % coord(idime) + dt_k * lagrtyp(ilagr) % veloc(idime)  &
                               &         + dt_k * dt_k * ( 0.5_rp*lagrtyp(ilagr) % accel(idime) &
                               &         + beta_pts *(accel_kp1(idime)-lagrtyp(ilagr) % accel(idime)) )
                       end do
                       ! Compute instantaneous & effective Stokes numbers (Nicolaou & Zaki, JAS, 2016)
                       call pts_stk_local(ielem, tau_p, lagrtyp(ilagr) % t_inject, Stk(1), Stk(2))
                    end if
                 end if

              end if

           end if

           !----------------------------------------------------------
           !
           !                  -----------
           !                  RANDOM WALK
           !                  -----------
           !
           !----------------------------------------------------------

           if( parttyp(itype) % kfl_brown == 1 ) then
              if( ndime == 2 ) then
                 coord_kp1(1) = coord_kp1(1) + eps1
                 coord_kp1(2) = coord_kp1(2) + eps2
              else
                 coord_kp1(1) = coord_kp1(1) + eps1
                 coord_kp1(2) = coord_kp1(2) + eps2
                 coord_kp1(3) = coord_kp1(3) + eps3
              end if
           end if

           !----------------------------------------------------------
           !
           !                    ---------------
           !                    WALL DEPOSITION
           !                    ---------------
           !
           ! Check if trajectory crosses the wall
           ! IWALL = 0 ... stay in domain, slip, bouncing
           !       = 1 ... outflow, wall
           !
           !----------------------------------------------------------

           iwall = 0
           if( lboue_pts(ielem) > 0 .or. kfl_walld_pts > 0 ) then
              call pts_cross_wall(&
                   ielem,pnode,ilagr,diame,hleng,toler,shapf,deriv,iwall,&
                   iboun,coord_kp1,veloc_kp1,accel_kp1,xinte,t,dt_k)
  
              !print *,'^^^ part=',lagrtyp(ilagr) % ilagr, 'iwall after cross=',iwall,'; kfl_paral=',kfl_paral
              if( iwall == 1 ) then
                 !
                 ! Hit a wall: go out of time loop
                 !
                 inscont = .false.
                 !print *,'part=',lagrtyp(ilagr) % ilagr, 'stopping the loop'
                 lagrtyp(ilagr) % iboun          = iboun
                 lagrtyp(ilagr) % coord(1:ndime) = xinte(1:ndime)
                 lagrtyp(ilagr) % t              = t
                 goto 20

              end if
           end if

           !----------------------------------------------------------
           !
           ! Adaptive time step
           !
           !----------------------------------------------------------

           accept_time_step = .true.
           if( parttyp(itype) % kfl_tstep > 0 .and. parttyp(itype) % kfl_modla == 2 ) then
              if( parttyp(itype) % chale == -3.0_rp ) then
                 !
                 ! Tau, h = tau * |uf|
                 !
                 uu = 0.0_rp
                 do idime = 1,ndime
                    uu = uu + veloc_kp1(idime)*veloc_kp1(idime)
                 end do
                 h = denpa * diame * diame / ( 18.0_rp * visfl ) * (sqrt(uu)+zeror)

              else if( parttyp(itype) % chale == -2.0_rp ) then
                 !
                 ! Re = 1
                 !
                 h = visfl / ( denfl * ( Du + zeror ) )

              else if( parttyp(itype) % chale == -1.0_rp ) then
                 !
                 ! Mesh size
                 !
                 h = hleng

              else if( parttyp(itype) % chale == 0.0_rp ) then
                 !
                 ! Particle diameter
                 !
                 h = sqrt(diame)

              else
                 !
                 ! User-defined
                 !
                 h = parttyp(itype) % chale

              end if
              call pts_adapti(&
                   parttyp(itype) % kfl_tstep,h,coord_kp1,lagrtyp(ilagr)%coord,veloc_kp1,lagrtyp(ilagr)%veloc,&
                   v_fluid_k,accel_kp1,lagrtyp(ilagr)%accel,tau,dt_k,lagrtyp(ilagr)%dt_km1,parttyp(itype) % safet,&
                   alpha_str,accept_time_step,parttyp(itype) % kfl_modla)
              if( .not. accept_time_step ) then
                 t = lagrtyp(ilagr) % t
                 goto 20
              end if
           else
              alpha_str = ovstr
           end if
           !
           ! Critical time step based on element minimum length HLENG
           !
           venor = sqrt(veloc_kp1(1)*veloc_kp1(1) + veloc_kp1(2)*veloc_kp1(2) + veloc_kp1(3)*veloc_kp1(3))
           if( venor /= 0.0_rp .or. D /= 0.0_rp ) then
              dtc = 1.0_rp / ( venor / hleng + 2.0_rp * D / hleng**2 )
           else
              dtc = 1.0e12_rp
           end if
           !
           ! Maximum number of iterations inside loop
           !
           if( itint == 1 .and. 1==2) then
              !uu = 0.0_rp
              !vv = 0.0_rp
              !do idime = 1,ndime
              !   uu = uu + veloc_kp1(idime)*veloc_kp1(idime)
              !end do
              !do idime = 1,ndime
              !   vv = vv + lagrtyp(ilagr) % v_fluid_k(idime)*lagrtyp(ilagr) % v_fluid_k(idime)
              !end do
              !if( parttyp(itype) kfl_tstep <= 0 .and. parttyp(itype) % kfl_modla /= 2 ) h = 1.0_rp
              !tau_p  = denpa * diame * diame / ( 18.0_rp * visfl )  ! Particle response time
              !itmaxr = ( abs(uu-vv) /( uu*tau_p*h + zeror))*dtime
              !itmax0 =  int(itmaxr)
              !if( itmax0 < 2000 ) itmax0 = 2000
              !itmax = min(itmax,itmax0)
           end if

           !----------------------------------------------------------
           !
           !                      -------------------
           !                      SEARCH HOST ELEMENT
           !                      -------------------
           !
           !----------------------------------------------------------
           !
           ! Start by checking in IELEM... Is it really likely?
           !
           ifoun = 0
           jelem = ielem
           pelty = ltype(jelem)
           pnode = nnode(pelty)
           do inode = 1,pnode
              ipoin = lnods(inode,jelem)
              elcod(1:ndime,inode) = coord(1:ndime,ipoin)
           end do
!print*,'otro elmgeo',toler
           call elmgeo_natural_coordinates(          &
                ndime,pelty,pnode,elcod,shapf,deriv, &
                coord_kp1,coloc,ifoun,toler)
!print*,'Holaaaaaaaaaa-busco coordenada del elemento',ifoun,coord_kp1,jelem
           if( ifoun == 0 ) then
!print*,'holaaaaaaaaaa-bin?',kfl_usbin_pts
              if( kfl_usbin_pts == 1 ) then
                 !
                 ! Use element neighboring bin to reduce the search
                 !
                 call maths_mapping_coord_to_3d(&
                      ndime,element_bin_boxes,element_bin(ielem) % comin,&
                      element_bin(ielem) % comax,coord_kp1,ii,jj,kk)
                 if( ii*jj*kk /= 0 ) then
                    ielel = 1
                    melel = element_bin(ielem) % bin_size(ii,jj,kk)
                    do while( ifoun == 0 .and. ielel <= melel )
                       jelem = element_bin(ielem) % list_elements(ii,jj,kk) % l(ielel)
                       if( jelem /= ielem ) then
                          pelty = ltype(jelem)
                          pnode = nnode(pelty)
                          do inode = 1,pnode
                             ipoin = lnods(inode,jelem)
                             elcod(1:ndime,inode) = coord(1:ndime,ipoin)
                          end do
!print*,'mas elmgeo',toler
                          call elmgeo_natural_coordinates(         &
                               ndime,pelty,pnode,elcod,shapf,deriv,&
                               coord_kp1,coloc,ifoun,toler)
                          !print*,jelem,ifoun
                       end if
                       ielel = ielel + 1
                    end do
                 end if
              else
                 !
                 ! Look for host element JELEM in neigboring list of IELEM
                 !
!if(t>8.0e-2_rp)
!print*,'estoy en ifoun=0',coord_kp1
             ! do ielel = pelel_2(ielem),pelel_2(ielem+1)
             !    if( ielel == pelel_2(ielem+1) ) then
             !       jelem = ielem
             !    else
             !       jelem = lelel_2(ielel)
             !    end if
             !   print*,'-----',jelem,lnods(1:nnode(ltype(jelem)),jelem)
             !    !do inode = 1,nnode(ltype(jelem))
             !    !   write(6,*) inode,lnods(inode,jelem),coord(1:ndime,lnods(inode,jelem))
                 !end do
             ! end do
                 ielel = pelel_2(ielem)
                 melel = pelel_2(ielem+1)-1
                 do while( ifoun == 0 .and. ielel <= melel )
                    jelem = lelel_2(ielel)
                    pelty = ltype(jelem)
                    pnode = nnode(pelty)
                    do inode = 1,pnode
                       ipoin = lnods(inode,jelem)
                       elcod(1:ndime,inode) = coord(1:ndime,ipoin)
                    end do
!if(t>8.0e-2_rp)
!print*,'elmgeo-en look host',jelem
                    call elmgeo_natural_coordinates(          &
                         ndime,pelty,pnode,elcod,shapf,deriv, &
                         coord_kp1,coloc,ifoun,toler)
!print*,'salgo con',ifoun,ielel,melel
                    ielel = ielel + 1
                 end do
!if(t>8.0e-2_rp)
!print*,'salgo de aqui con ifoun',ifoun,'y jeleme',jelem
              end if
           end if

           if( ifoun /= 0 ) then
!!$              print*,'c=',leinv_loc(ielem),lagrtyp(ilagr) % coord(1:ndime)
!!$              print*,'d=',leinv_loc(ielem),coord_kp1(1:ndime)
!!$              do ielel = pelel_2(ielem),pelel_2(ielem+1)
!!$                 if( ielel == pelel_2(ielem+1) ) then
!!$                    jelem = ielem
!!$                 else
!!$                    jelem = lelel_2(ielel)
!!$                 end if
!!$                print*,'-----',jelem,ielem,lelel_2(ielel)
!!$                 do inode = 1,nnode(ltype(jelem))
!!$                    write(6,*) inode,lnods(inode,jelem),coord(1:ndime,lnods(inode,jelem))
!!$                 end do
!!$              end do
              !
              ! Residence time
              !
              if( kfl_resid_pts > 0 ) then
                 if( jelem == ielem ) then
                    !$OMP ATOMIC
                    resid_pts(itype,ielem) = resid_pts(itype,ielem) + dt_k
                 else
                    !
                    ! Intersection between segment and face should be computed in mod_elmgeo
                    ! If JELEM is a ghost element, info will be further sent in pts_endite()
                    !
                    !$OMP ATOMIC
                    resid_pts(itype,ielem) = resid_pts(itype,ielem) + 0.5_rp * dt_k
                    !$OMP ATOMIC
                    resid_pts(itype,jelem) = resid_pts(itype,jelem) + 0.5_rp * dt_k
                 end if
              end if
              !
              ! Coupling with Nastin
              !
              if( kfl_coupl(ID_NASTIN,ID_PARTIS) /= 0 ) then
                 do inode = 1,pnode
                    ipoin = lnods(inode,jelem)
                    do idime = 1,ndime
                       !$OMP ATOMIC
                       momen(idime,ipoin) = momen(idime,ipoin)        &
                            &  - mass * shapf(inode) * dt_k / dtime * &
                            &  ( lagrtyp(ilagr) % acced(idime)      + &
                            &    lagrtyp(ilagr) % acceg(idime)        )
                    end do
                 end do
              end if

              if( jelem > nelem ) then
                 !
                 ! JELEM > NELEM: particle goes to neighboring subdomain INEIG
                 !
                 ineig                      = leldo(1,jelem-nelem)         ! Which neighbor holds this element
                 lagrtyp(ilagr) % kfl_exist = ineig                        ! Neighbor
!if(lagrtyp(ilagr) % kfl_exist==-3)print*,'parti con exist -3',lagrtyp(ilagr) % ilagr,lagrtyp(ilagr) % kfl_exist
                 inscont                    = .false.                      ! Go out of internal time loop
                 ncros(ineig)               = ncros(ineig) + 1             ! Number of particles to send to neighbor INEIG
                 lagrtyp(ilagr) % ielem     = leldo(2,jelem-nelem)         ! Local numbering of JELEM in my neighbor
!if(t>8.0e-2_rp)print*,'pillo el elem-1',leinv_loc(lagrtyp(ilagr) % ielem)
              else
                 !
                 ! Element found: accept time step update new position
                 !
                 lagrtyp(ilagr) % ielem = jelem                            ! JELEM remains in my subdomain
!if(t>8.0e-2_rp)print*,'pillo el elem-2',leinv_loc(lagrtyp(ilagr) % ielem)
              end if
              lagrtyp(ilagr) % ittim                = lagrtyp(ilagr) % ittim + 1 ! One new time step
              lagrtyp(ilagr) % Cd                   = Cd                         ! Accept drag
              lagrtyp(ilagr) % Re                   = Re                         ! Accept Reynolds
              lagrtyp(ilagr) % Stk                  = Stk                        ! Stokes number ! (Stk(1): instantaneous, Stk(2): effective)
              lagrtyp(ilagr) % stret                = alpha_str                  ! Stretching
              lagrtyp(ilagr) % t                    = t                          ! Accept new time
              lagrtyp(ilagr) % dt_k                 = dtg                        ! Guess for next time step
              lagrtyp(ilagr) % dt_km2               = lagrtyp(ilagr) % dt_km1    ! Save time step
              lagrtyp(ilagr) % dt_km1               = dt_k                       ! Save time step
              lagrtyp(ilagr) % coord(1:ndime)       = coord_kp1(1:ndime)         ! Accept position
              lagrtyp(ilagr) % veloc(1:ndime)       = veloc_kp1(1:ndime)         ! Accept velocity
              lagrtyp(ilagr) % accel(1:ndime)       = accel_kp1(1:ndime)         ! Accept acceleration
              lagrtyp(ilagr) % v_fluid_km2(1:ndime) = lagrtyp(ilagr) % v_fluid_km1(1:ndime)
              lagrtyp(ilagr) % v_fluid_km1(1:ndime) = lagrtyp(ilagr) % v_fluid_k(1:ndime)

           else
              !
              ! Element not found: reduce time step
              !
              t         = lagrtyp(ilagr) % t     ! Go back to previous time
!print*,'ELEMENT NOT FOUND',t
              alpha_str = strex                  ! Decrease time step
              accept_time_step = .false.         ! Do not accept time step
           end if
           !
           ! Very small dt: remove particle
           !
20         continue
           !print *,'*** particle=',lagrtyp(ilagr) % ilagr, '; iwall=',iwall, '; kfl_parall=', kfl_paral
           if( ( dt_k < 0.1_rp * dtmin_pts .or. itint > itmax ) .and. iwall /= 1 ) then
              inscont = .false.                    ! Go out of time loop

              if(.not. newmark_converged) then                 
                 print *, 'Newmark did not converge after all iteraitons. Particle id=', lagrtyp(ilagr) % ilagr,'; coords=',lagrtyp(ilagr) % coord(1:ndime),', xerro=',xerro  
              end if

              if( lelbf(ielem) % n /= 0 ) then
                 !print*,'paso por aqui-A'
                 !lagrtyp(ilagr) % kfl_exist = -2   ! Assume it hits the wall
                 lagrtyp(ilagr) % kfl_exist = -3   ! It is lost
                 !print *, 'Particle hits wall 2', lagrtyp(ilagr) % ilagr,lagrtyp(ilagr) % coord(1:ndime),ielem
              else
                 lagrtyp(ilagr) % kfl_exist = -3   ! Assume it is lost for numerical reasons... :o(
                 !print *, 'Particle lost for numerical reasons', lagrtyp(ilagr) % ilagr,lagrtyp(ilagr) % coord(1:ndime),', xerro=',xerro  
              end if
              !t = tf
           end if

           newmark_converged = .true. !reset for the next iteration
           !print *,'ilagr=',lagrtyp(ilagr) % ilagr,'t=',t,' ;dt_k=',dt_k,' ; inscont=',inscont
        end do
        !end if
        if( lagrtyp(ilagr) % kfl_exist < 0 ) nlagr_non_migrating_pts = nlagr_non_migrating_pts + 1

     end do
     !$OMP END PARALLEL DO
     !
     ! All reduce total number of non-crossing particles
     ! if NLAGR_EXISTING_PTS /= NLAGR_NON_MIGRATING_PTS, some particles migrate to other subdomains
     !


     call cputim(time2)
     time_total = time_total + time2 - time1

     call PAR_SUM(nlagr_non_migrating_pts,'IN MY ZONE')
     !print *,'after par_sum: nlagr_existing_pts=',nlagr_existing_pts,', nlagr_non_migrating_pts=',nlagr_non_migrating_pts
     if( nlagr_non_migrating_pts > nlagr_existing_pts ) nlagr_non_migrating_pts = nlagr_existing_pts    ! Careful, this is a temporal patch.

     if( nlagr_existing_pts == nlagr_non_migrating_pts ) then
        !
        ! Communication is over
        !
        comcont = -1

     else if( ISLAVE ) then
        !
        ! Send/receive number of particles
        !
        !call PAR_SEND_RECEIVE_TO_ALL(ncros,ncror,'IN MY ZONE')
        do ineig = 1,nneig
           call PAR_SEND_RECEIVE(1_ip,1_ip,ncros(ineig:ineig),ncror(ineig:ineig),'IN MY ZONE',commu % neights(ineig) )
        end do
        !
        ! Count total [number of particles to receive] - [number of particles to send]
        ! UPS!!! I cannot take into account the send particles because the place is freed
        ! as we go through the neigbors
        !
        do ineig = 1,nneig
           particles_sent = particles_sent + ncros(ineig)
           particles_recv = particles_recv + ncror(ineig)
        end do
        nlagr_new = particles_recv !- particles_sent

        if( nlagr_new > nlagr_free ) then
           new_size = int(1.2_rp*real(mlagr+nlagr_new-nlagr_free,rp),ip)
           call pts_reallocate(new_size)
           call memory_resize(mem_modul(1:2,modul),'PERMU_NLAGR','pts_solite',permu_nlagr,mlagr)
        end if
        !
        ! Send/Recv crossing particles to neighbors: coord, tinit, ilagr, etc.
        ! This is the minimum information for the particle to be correctly tracked by
        ! the neighboring subdomain
        !
        do ineig = 1,nneig

           nullify(parrs_pts)
           nullify(parre_pts)
           if( ncros(ineig) > 0 ) then
              npasr = ncros(ineig) * nvar1
              call memory_alloca(mem_modul(1:2,modul),'PARRS_PTS','pts_solite',parrs_pts,npasr,'DO_NOT_INITIALIZE')
              ipars = 0
              !do ilagr = 1,mlagr
              do ilagr_local = 1,nlagr_local
                 ilagr = permu_nlagr(ilagr_local)
                 if( lagrtyp(ilagr) % kfl_exist == ineig ) then
                    ipars            = ipars + 1                      ! 1
                    parrs_pts(ipars) = lagrtyp(ilagr) % coord(1)
                    ipars            = ipars + 1                      ! 2
                    parrs_pts(ipars) = lagrtyp(ilagr) % coord(2)
                    ipars            = ipars + 1                      ! 3
                    parrs_pts(ipars) = lagrtyp(ilagr) % coord(3)
                    ipars            = ipars + 1                      ! 4
                    parrs_pts(ipars) = lagrtyp(ilagr) % veloc(1)
                    ipars            = ipars + 1                      ! 5
                    parrs_pts(ipars) = lagrtyp(ilagr) % veloc(2)
                    ipars            = ipars + 1                      ! 6
                    parrs_pts(ipars) = lagrtyp(ilagr) % veloc(3)
                    ipars            = ipars + 1                      ! 7
                    parrs_pts(ipars) = lagrtyp(ilagr) % accel(1)
                    ipars            = ipars + 1                      ! 8
                    parrs_pts(ipars) = lagrtyp(ilagr) % accel(2)
                    ipars            = ipars + 1                      ! 9
                    parrs_pts(ipars) = lagrtyp(ilagr) % accel(3)
                    ipars            = ipars + 1                      ! 10: Useful for nastin coupling
                    parrs_pts(ipars) = lagrtyp(ilagr) % acced(1)
                    ipars            = ipars + 1                      ! 11
                    parrs_pts(ipars) = lagrtyp(ilagr) % acced(2)
                    ipars            = ipars + 1                      ! 12
                    parrs_pts(ipars) = lagrtyp(ilagr) % acced(3)
                    ipars            = ipars + 1                      ! 13
                    parrs_pts(ipars) = lagrtyp(ilagr) % stret
                    ipars            = ipars + 1                      ! 14
                    parrs_pts(ipars) = lagrtyp(ilagr) % t
                    ipars            = ipars + 1                      ! 15
                    parrs_pts(ipars) = real(lagrtyp(ilagr) % ilagr,rp)
                    ipars            = ipars + 1                      ! 16
                    parrs_pts(ipars) = real(lagrtyp(ilagr) % itype,rp)
                    ipars            = ipars + 1                      ! 17
                    parrs_pts(ipars) = real(lagrtyp(ilagr) % ielem,rp)
                    ipars            = ipars + 1                      ! 18
                    parrs_pts(ipars) = real(lagrtyp(ilagr) % ittim,rp)
                    ipars            = ipars + 1                      ! 19
                    parrs_pts(ipars) = lagrtyp(ilagr) % dt_k
                    ipars            = ipars + 1                      ! 20
                    parrs_pts(ipars) = lagrtyp(ilagr) % dt_km1
                    ipars            = ipars + 1                      ! 21
                    parrs_pts(ipars) = lagrtyp(ilagr) % dt_km2
                    ipars            = ipars + 1                      ! 22: Prediction time step
                    parrs_pts(ipars) = lagrtyp(ilagr) % dtg
                    ipars            = ipars + 1                      ! 23
                    parrs_pts(ipars) = lagrtyp(ilagr) % v_fluid_k(1)
                    ipars            = ipars + 1                      ! 24
                    parrs_pts(ipars) = lagrtyp(ilagr) % v_fluid_k(2)
                    ipars            = ipars + 1                      ! 25
                    parrs_pts(ipars) = lagrtyp(ilagr) % v_fluid_k(3)
                    ipars            = ipars + 1                      ! 26
                    parrs_pts(ipars) = lagrtyp(ilagr) % v_fluid_km1(1)
                    ipars            = ipars + 1                      ! 27
                    parrs_pts(ipars) = lagrtyp(ilagr) % v_fluid_km1(2)
                    ipars            = ipars + 1                      ! 28
                    parrs_pts(ipars) = lagrtyp(ilagr) % v_fluid_km1(3)
                    ipars            = ipars + 1                      ! 29
                    parrs_pts(ipars) = lagrtyp(ilagr) % v_fluid_km2(1)
                    ipars            = ipars + 1                      ! 30
                    parrs_pts(ipars) = lagrtyp(ilagr) % v_fluid_km2(2)
                    ipars            = ipars + 1                      ! 31
                    parrs_pts(ipars) = lagrtyp(ilagr) % v_fluid_km2(3)
		    ipars            = ipars + 1                      ! 32
		    parrs_pts(ipars) = lagrtyp(ilagr) % t_inject
		    ipars            = ipars + 1                      ! 33
		    parrs_pts(ipars) = lagrtyp(ilagr) % Stk(2)		
                    do iprop = 1,nlapr
                       ipars            = ipars + 1
                       parrs_pts(ipars) = lagrtyp(ilagr) % prope(iprop)
                    enddo
                    lagrtyp(ilagr) % kfl_exist = 0                ! Take off particle from my list
                 end if
              end do
           end if

           if( ncror(ineig) > 0 ) then
              nparr = ncror(ineig) * nvar1
              call memory_alloca(mem_modul(1:2,modul),'PARRE_PTS','pts_solite',parre_pts,nparr)
           end if

           call PAR_SEND_RECEIVE(parrs_pts,parre_pts,'IN MY ZONE',commu % neights(ineig) )
           !
           ! Allocate new particles
           !
           ilagr      = 1
           iparr      = 0
           nlagr_last = 0

           do icror = 1,ncror(ineig)
              !
              ! Look for a free space in LAGRTYP to save particle
              ! If there is no space, allocate more memory
              !
              ifoun = 0
              ilagr = nlagr_last
              loop_find_position: do while( ilagr < mlagr )
                 ilagr = ilagr + 1
                 if( lagrtyp(ilagr) % kfl_exist == 0 ) then
                    ifoun = 1
                    exit loop_find_position
                 end if
              end do loop_find_position
              nlagr_last = ilagr
              if( ifoun == 0 ) then
                 call runend('PTS_SOLITE: WE ARE IN TROUBLE!')
                 !   new_size =  int(1.2_rp*real(nlagr_new,rp),ip)
                 !   ilagr = mlagr + 1
                 !   call pts_lagdef(2_ip,0_ip)
                 !   call pts_reallocate(
                 !   call memory_resize(mem_modul(1:2,modul),'PERMU_NLAGR','pts_solite',permu_nlagr,mlagr)
              end if
              lagrtyp(ilagr) = lagrtyp_init
              iparr                           = iparr + 1             ! 1
              lagrtyp(ilagr) % coord(1)       = parre_pts(iparr)
              iparr                           = iparr + 1             ! 2
              lagrtyp(ilagr) % coord(2)       = parre_pts(iparr)
              iparr                           = iparr + 1             ! 3
              lagrtyp(ilagr) % coord(3)       = parre_pts(iparr)
              iparr                           = iparr + 1             ! 4
              lagrtyp(ilagr) % veloc(1)       = parre_pts(iparr)
              iparr                           = iparr + 1             ! 5
              lagrtyp(ilagr) % veloc(2)       = parre_pts(iparr)
              iparr                           = iparr + 1             ! 6
              lagrtyp(ilagr) % veloc(3)       = parre_pts(iparr)
              iparr                           = iparr + 1             ! 7
              lagrtyp(ilagr) % accel(1)       = parre_pts(iparr)
              iparr                           = iparr + 1             ! 8
              lagrtyp(ilagr) % accel(2)       = parre_pts(iparr)
              iparr                           = iparr + 1             ! 9
              lagrtyp(ilagr) % accel(3)       = parre_pts(iparr)
              iparr                           = iparr + 1             ! 10: Useful for nastin co
              lagrtyp(ilagr) % acced(1)       = parre_pts(iparr)
              iparr                           = iparr + 1             ! 11
              lagrtyp(ilagr) % acced(2)       = parre_pts(iparr)
              iparr                           = iparr + 1             ! 12
              lagrtyp(ilagr) % acced(3)       = parre_pts(iparr)
              iparr                           = iparr + 1             ! 13
              lagrtyp(ilagr) % stret          = parre_pts(iparr)
              iparr                           = iparr + 1             ! 14
              lagrtyp(ilagr) % t              = parre_pts(iparr)
              iparr                           = iparr + 1             ! 15
              lagrtyp(ilagr) % ilagr          = int(parre_pts(iparr),ip)
              iparr                           = iparr + 1             ! 16
              lagrtyp(ilagr) % itype          = int(parre_pts(iparr),ip)
              iparr                           = iparr + 1             ! 17
              lagrtyp(ilagr) % ielem          = int(parre_pts(iparr),ip)
              iparr                           = iparr + 1             ! 18
              lagrtyp(ilagr) % ittim          = int(parre_pts(iparr),ip)
              iparr                           = iparr + 1             ! 19
              lagrtyp(ilagr) % dt_k           = parre_pts(iparr)
              iparr                           = iparr + 1             ! 20
              lagrtyp(ilagr) % dt_km1         = parre_pts(iparr)
              iparr                           = iparr + 1             ! 21
              lagrtyp(ilagr) % dt_km2         = parre_pts(iparr)
              iparr                           = iparr + 1             ! 22
              lagrtyp(ilagr) % dtg            = parre_pts(iparr)
              iparr                           = iparr + 1             ! 23
              lagrtyp(ilagr) % v_fluid_k(1)   = parre_pts(iparr)
              iparr                           = iparr + 1             ! 24
              lagrtyp(ilagr) % v_fluid_k(2)   = parre_pts(iparr)
              iparr                           = iparr + 1             ! 25
              lagrtyp(ilagr) % v_fluid_k(3)   = parre_pts(iparr)
              iparr                           = iparr + 1             ! 26
              lagrtyp(ilagr) % v_fluid_km1(1) = parre_pts(iparr)
              iparr                           = iparr + 1             ! 27
              lagrtyp(ilagr) % v_fluid_km1(2) = parre_pts(iparr)
              iparr                           = iparr + 1             ! 28
              lagrtyp(ilagr) % v_fluid_km1(3) = parre_pts(iparr)
              iparr                           = iparr + 1             ! 29
              lagrtyp(ilagr) % v_fluid_km2(1) = parre_pts(iparr)
              iparr                           = iparr + 1             ! 30
              lagrtyp(ilagr) % v_fluid_km2(2) = parre_pts(iparr)
              iparr                           = iparr + 1             ! 31
              lagrtyp(ilagr) % v_fluid_km2(3) = parre_pts(iparr)
	      iparr                           = iparr + 1             ! 32
	      lagrtyp(ilagr) % t_inject       = parre_pts(iparr)
	      iparr                           = iparr + 1             ! 33
	      lagrtyp(ilagr) % Stk(2)         = parre_pts(iparr)
              do iprop = 1,nlapr
                 iparr = iparr + 1
                 lagrtyp(ilagr) % prope(iprop) = parre_pts(iparr)
              enddo
              lagrtyp(ilagr) % kfl_exist = -1
           end do
           !
           ! Deallocate
           !
           call memory_deallo(mem_modul(1:2,modul),'PARRE_PTS','pts_solite',parre_pts)
           call memory_deallo(mem_modul(1:2,modul),'PARRS_PTS','pts_solite',parrs_pts)
           ncros(ineig) = 0
           ncror(ineig) = 0

        end do

     end if
  end do

  !
  ! Waht I send and receive from neighbors
  !
  deallocate(ncror)
  deallocate(ncros)
  !
  ! Deallocate list of existing particles
  !
  if( INOTMASTER ) then
     call memory_deallo(mem_modul(1:2,modul),'PERMU_NLAGR','pts_solite',permu_nlagr)
     !
     ! Save particles on deposition map
     !
     if( kfl_depos_pts == 1 .OR. kfl_depos_surface_pts /= 0) then
        do ilagr = 1,mlagr
           !do ilagr_local = 1,nlagr_local
           !ilagr = permu_nlagr(ilagr_local)
           if( lagrtyp(ilagr) % kfl_exist == -2 ) then
              ielem = lagrtyp(ilagr) % ielem
              iboun = lagrtyp(ilagr) % iboun
              if( iboun > 0 ) then
                 !depoe_pts(lagrtyp(ilagr) % itype,ielem) = depoe_pts(lagrtyp(ilagr) % itype,ielem) + 1.0_rp
                 depob_pts(lagrtyp(ilagr) % itype,iboun) = depob_pts(lagrtyp(ilagr) % itype,iboun) + 1.0_rp
              end if
           end if
        end do
     end if
  end if


  
  !----------------------------------------------------------------------
  !
  ! Particles out of the domain: put a negative particle number
  ! Particles with zero time step: take them off
  !
  !----------------------------------------------------------------------
  
  nlagr_going_out_pts = 0
  nlagr_zero_time_pts = 0
  nlagr_deposited_pts = 0
  
  if( INOTMASTER ) then
     do ilagr = 1,mlagr
        !do ilagr_local = 1,nlagr_local
        !ilagr = permu_nlagr(ilagr_local)
        if( lagrtyp(ilagr) % kfl_exist == -2 .or. lagrtyp(ilagr) % kfl_exist == -4 ) then
           !
           ! Deposited or outflow particle
           !
           nlagr_going_out_pts = nlagr_going_out_pts + 1
           lagrtyp(ilagr) % ilagr = -abs(lagrtyp(ilagr) % ilagr) ! Mark as deposited
           if( lagrtyp(ilagr) % boundary_set == -1 ) nlagr_deposited_pts = nlagr_deposited_pts + 1
           
        else if( lagrtyp(ilagr) % kfl_exist == -3 ) then
           !
           ! Zero time step
           !
           nlagr_zero_time_pts = nlagr_zero_time_pts + 1
        end if
        
     end do

  end if
  my_nlagr_4 = int(nlagr_zero_time_pts,4)
  dummi(1)   = nlagr_going_out_pts
  dummi(2)   = nlagr_zero_time_pts
  dummi(3)   = nlagr_deposited_pts
 
  call PAR_SUM(3_ip,dummi,'IN MY ZONE')

  nlagr_going_out_pts  = dummi(1)
  nlagr_zero_time_pts  = dummi(2)
  nlagr_deposited_pts  = dummi(3)
  
  if( nlagr_going_out_pts > 0 ) call livinf(-9_ip,'PARTICLES GOING OUT OF COMPUTATIONAL DOMAIN= ',nlagr_going_out_pts)
  if( nlagr_zero_time_pts > 0 ) call livinf(-9_ip,'PARTICLES WITH ZERO TIME STEP= ',nlagr_zero_time_pts)
  if( nlagr_deposited_pts > 0 ) call livinf(-9_ip,'PARTICLES DEPOSITED OUT OF BOUNDARY= ',nlagr_deposited_pts)

  !----------------------------------------------------------------------
  !
  ! Events: particles are lost
  !
  !----------------------------------------------------------------------
  if( nlagr_zero_time_pts > 0 .and. kfl_detection /= 0 ) then
     !if( kfl_detection /= 0 ) then

     !my_nlagr_4 = 0
     !do ilagr = 1,mlagr
     !   if( lagrtyp(ilagr) % kfl_exist == -1 ) then
     !      my_nlagr_4 = my_nlagr_4 + 1
     !   end if
     !end do
     !nlagr_4 = my_nlagr_4
     !call PAR_SUM(nlagr_4,'IN MY ZONE')

     if( IPARALL ) then
        !
        ! All gather particle coordinates
        !
        call memory_alloca(mem_modul(1:2,modul),'PAR_NLAGR_4','pts_solite',par_nlagr_4,int(PAR_CODE_SIZE,4),'INITIALIZE',0_4)
        call PAR_ALLGATHER(my_nlagr_4,par_nlagr_4,1_4,'IN MY CODE')
        if( INOTMASTER ) then
           call memory_alloca(mem_modul(1:2,modul),'MY_COORD_NLAGR_4','pts_solite',my_coord_nlagr_4,ndime,int(my_nlagr_4,ip),'INITIALIZE')
           call memory_alloca(mem_modul(1:2,modul),'MY_VELOC_NLAGR_4','pts_solite',my_veloc_nlagr_4,ndime,int(my_nlagr_4,ip),'INITIALIZE')
           my_nlagr_4 = 0
           do ilagr = 1,mlagr
              !do ilagr_local = 1,nlagr_local
              !ilagr = permu_nlagr(ilagr_local)
              if( lagrtyp(ilagr) % kfl_exist == -3 ) then
                 !if( lagrtyp(ilagr) % kfl_exist == -1 ) then
                 my_nlagr_4 = my_nlagr_4 + 1
                 my_coord_nlagr_4(1:ndime,my_nlagr_4) = lagrtyp(ilagr) % coord(1:ndime)
                 my_veloc_nlagr_4(1:ndime,my_nlagr_4) = lagrtyp(ilagr) % veloc(1:ndime)
              end if
           end do
        end if
        ndime4 = int(ndime,4)
        do ipart4 = 0,int(PAR_CODE_SIZE,4)-1_4
           par_nlagr_4(ipart4) = par_nlagr_4(ipart4) * ndime4
        end do
        call memory_alloca(mem_modul(1:2,modul),'COORD_NLAGR_4','pts_solite',coord_nlagr_4,ndime,nlagr_zero_time_pts,'INITIALIZE')
        call PAR_ALLGATHERV(my_coord_nlagr_4,coord_nlagr_4,par_nlagr_4,'IN MY CODE')
        !
        ! Output event
        !
        call PAR_COMM_RANK_AND_SIZE(my_rank4,'IN MY ZONE')
        klagr = 0
        do ipart4 = 0,int(PAR_CODE_SIZE,4)-1_4
           par_nlagr_4(ipart4) = par_nlagr_4(ipart4) / ndime4
           do ilagr = 1,par_nlagr_4(ipart4)
              klagr = klagr + 1
              if( ipart4 == my_rank4 ) then
                 call ker_events_particle_not_converged(ipart4,coord_nlagr_4(1:ndime,klagr),ndime,advec,my_veloc_nlagr_4(1:ndime,klagr))
              else
                 call ker_events_particle_not_converged(ipart4,coord_nlagr_4(1:ndime,klagr),ndime,advec)
              end if
           end do
        end do

        call memory_deallo(mem_modul(1:2,modul),'COORD_NLAGR_4',   'pts_solite',coord_nlagr_4)
        call memory_deallo(mem_modul(1:2,modul),'MY_VELOC_NLAGR_4','pts_solite',my_veloc_nlagr_4)
        call memory_deallo(mem_modul(1:2,modul),'MY_COORD_NLAGR_4','pts_solite',my_coord_nlagr_4)
        call memory_deallo(mem_modul(1:2,modul),'PAR_NLAGR_4',     'pts_solite',par_nlagr_4)

     end if
  end if
  !print*,'salgo de solite, rank=',my_rank4

  !----------------------------------------------------------------------
  !
  ! Statistics law: test for random walk/Brownian motion
  !
  !----------------------------------------------------------------------
  if( test_brownian ) then
     if( ittim == 1000 ) then
        do jj = 1,1000
           toler = real(jj-1,rp)/999.0_rp*0.005_rp
           ii = 0
           do ilagr = 1,mlagr
              if( lagrtyp(ilagr) % kfl_exist == -1 ) then
                 xerro = 0.0_rp
                 do idime = 1,ndime
                    xerro = xerro + (lagrtyp(ilagr) % coord(idime)-0.005_rp)**2
                 end do
                 xerro = sqrt(xerro)
                 if( xerro < toler ) ii = ii + 1
              end if
           end do
           write(90,*) toler,ii ; flush(90)
        end do
     end if
  end if

  !----------------------------------------------------------------------
  !
  ! Timing
  !
  !----------------------------------------------------------------------

  cpu_modul(CPU_ASSEMBLY,modul) = cpu_modul(CPU_ASSEMBLY,modul) + time_total
  time_max = time_total
  time_ave = time_total

  call PAR_MAX(time_max)
  call PAR_SUM(time_ave)
  call PAR_MAX(particles_sent)
  call PAR_MAX(particles_recv)

  time_ave     = time_ave / max(1.0_rp,real(PAR_CODE_SIZE-1,rp))
  load_balance = time_ave / max(zeror,time_max)

  if( INOTSLAVE ) then
     write(momod(modul) % lun_conve,111) cutim,nlacc_pts,nlagr_existing_pts,nlagr_going_out_pts,nlagr_zero_time_pts,time_max,time_ave,load_balance,particles_sent,particles_recv,comm_loops_pts
     flush(momod(modul) % lun_conve)
  end if

  call ker_timeline('END_SOLVER',nlagr_existing_pts)
  !
  ! Formats
  !
111 format((e12.6,2x,4(2x,i7),3(2x,e12.6),3(2x,i7)))

#ifdef EVENT
  call mpitrace_user_function(0)
#endif

end subroutine pts_solite
!
!-----------------------------------------------------------------------
!
!   BROWNIAN MOTION
!   Activate option: test_brownian
!
!   #--------------
!   # USING GNUPLOT FOR 2D CASE
!   #--------------
!   #
!   D=5.35E-06
!   t=0.01
!   set xrange[0:0.0015]
!   plot 'fort.90' u 1:($2/2500.0) w p,  1.0 - exp( -x**2/(4.0*D*t) )
!
!
!
!
!
! Ratios
!
!dtacc  = 1.0e6_rp
!dt1    = dto       ! dt^n
!dt2    = dt2       ! dt^n-1
!dt1dt2 = dt1/dt2   ! dt^n/dt^n-1
!dt2dt1 = dt2/dt1   ! dt^n-1/dt^n
!if( xx(1) /= 0.0_rp ) then
!   e3e1  = x2(1) / xx(1)
!   e2e1  = x1(1) / xx(1)
!   xfac1 = 1.0_rp - (1.0_rp+dt1dt2)*e2e1 + dt1dt2*e3e1
!else
!   e3e1  = 1.0_rp
!   e2e1  = 1.0_rp
!   xfac1 = 0.0_rp
!end if
!!
!! Previous time step values
!!
!u1      = xx(1)
!u2      = x1(1)
!u3      = x2(1)
!dudt    = (   u1*(dt2dt1+2.0_rp) -   u2*(dt1dt2+dt2dt1+2.0_rp) +   u3*dt1dt2) / (dt1+dt2)
!dudtb   = u1 * ( (dt2dt1+2.0_rp) - e2e1*(dt1dt2+dt2dt1+2.0_rp) + e3e1*dt1dt2) / (dt1+dt2)
!d2udt2  = (   u1*dt2dt1 -   u2*(1.0_rp+dt2dt1) +   u3) / (0.5_rp*dt2*(dt1+dt2))
!d2udt2b = u1 * ( 1.0_rp - (1.0_rp+dt1dt2)*e2e1 + dt1dt2*e3e1 ) / ( 0.5_rp * dt1 * ( dt1 + dt2 ) )
!if( abs(xfac1) > zeror .and. abs(d2udt2b) > zeror ) then
!   err1 = 0.5_rp * dt1 * ( dt1 + dt2 ) / xfac1        ! |e| / |d^2e/dt^2| (local truncation)
!   err2 =  ( dt2dt1 + (dt1dt2-dt2dt1)*e2e1 &          ! |de/dt| / |d^2e/dt^2| (high order terms)
!        &  - dt1dt2*e3e1 ) / xfac1 * 0.5_rp * dt1
!   dtht = sqrt(2.0_rp * 0.03_rp * abs(err1))
!   dttr = 2.0_rp * 0.03_rp * abs(err2)
!else
!   dtht = 1.0_rp / strex * dto
!   dttr = 1.0_rp / strex * dto
!end if
!!
!! Criteria
!!
!dtall = min(dtht,dttr)
!dtacc = dtall
!!
!! Calculation of accuracy reached with precious timestep dtold:
!! facc = max(dt_trunc,dt_HT) / dtold
!!
!facc = dtacc / dto
!!
!! Prediction of next time step dtn+1 to be used in the solver
!! High damping: more stiff
!! Low damping:  more smooth
!!
!a = facc
!!d = dampi_chm
!d = 2.0_rp
!if( a >= 1.0_rp ) then
!   s = 1.0_rp/strex
!else
!   s = strex
!end if
!xfact = (s-1.0_rp)*tanh( (a-1.0_rp)/(d*(s-1.0_rp)) ) + 1.0_rp
!!stop
!dt    = dto * xfact
!!
!! NEW
!!
!!$
!!$             if( parttyp(itype) % kfl_modla == 2 .and. itrex == 1000 ) then
!!$                 Du = 0.0_rp
!!$                 Da = 0.0_rp
!!$                 do idime = 1,ndime
!!$                    Du = Du + ( lagrtyp(ilagr) % veloc(idime) ) ** 2
!!$                    Da = Da + ( lagrtyp(ilagr) % accel(idime)-lagrtyp(ilagr) % acceo(idime) ) ** 2
!!$                 end do
!!$                 Du = sqrt( Du )
!!$                 Da = sqrt( Da )
!!$                 if( Du > zeror ) then
!!$                    facc = abs(dt*dt*(1.0_rp/6.0_rp-beta_pts)*Da)/Du
!!$                 else
!!$                    facc = 1.0e-04_rp
!!$                 end if
!!$                 write(200,'(4(1x,e12.6))') t,facc,Du,Da
!!$                 if( facc > 1.0e-12_rp ) then
!!$                    facc = 1.0e-04_rp / facc
!!$                 else
!!$                    facc = 1.0e-04_rp / 1.0e-12
!!$                 end if
!!$                 strec = 1.5_rp
!!$                 dampi = 1.2_rp
!!$                 !dt    = dt * min(strec,max(1.0_rp/strec,log(1.0_rp + (dampi-1.0_rp)*facc) / log(dampi)))
!!$                 dt    = dt * min(strec,log(1.0_rp + (dampi-1.0_rp)*facc) / log(dampi))
!!$                 dt    = max(dt,dtmin_pts)
!!$                 dt    = min(dt,dtime)
!!$                 !dt_trn = dtmin_pts
!!$                 !dt_err = dtmin_pts
!!$                 !d2udt2 = 0.0_rp
!!$                 !dudt   = 0.0_rp
!!$                 !u      = 0.0_rp
!!$                 !do idime = 1,ndime
!!$                 !   d2udt2 = d2udt2 + lagrtyp(ilagr) % accel(idime) ** 2
!!$                 !   dudt   = dudt   + lagrtyp(ilagr) % veloc(idime) ** 2
!!$                 !   u      = u      + lagrtyp(ilagr) % veloc(idime) ** 2
!!$                 !end do
!!$                 !d2udt2 = sqrt(d2udt2)
!!$                 !dudt   = sqrt(dudt)
!!$                 !u      = sqrt(u)
!!$                 !if( d2udt2 > 0.0_rp ) then
!!$                 !   dt_trn = max(dt_trn,2.0_rp*1.0e-3_rp*dudt/d2udt2)
!!$                 !   dt_err = max(dt_err,sqrt(2.0_rp*1.0e-3_rp*u/d2udt2))
!!$                 !end if
!!$                 !facc  = max(dt_err,dt_trn) / dt
!!$                 !strec = 1.5_rp
!!$                 !dampi = 1.2_rp
!!$                 !dt    = dt * min(strec,log(1.0_rp + (dampi-1.0_rp)*facc) / log(dampi))
!!$                 !dt    = max(dt,dtmin_pts)
!!$                 !dt    = min(dt,dtime)
!!$              end if
