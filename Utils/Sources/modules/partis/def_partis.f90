module def_partis
  !------------------------------------------------------------------------
  !****f* Partis/def_partis
  ! NAME 
  !    def_partis
  ! DESCRIPTION
  !    Heading for the Partis routines
  ! USES
  ! USED BY
  !    Almost all
  !***
  !------------------------------------------------------------------------
  use def_kintyp
  use def_inpout, only : maxwp
  USE, INTRINSIC :: ISO_C_BINDING

  !------------------------------------------------------------------------
  !
  ! Lagrangian particle, needed by Master for coupling
  !
  !------------------------------------------------------------------------

  integer(ip)               :: mlagr             ! Lagrangian particles: Max. number of particle in each subdomain
  integer(ip), parameter    :: mlapr = 10        ! Maximum number of properties in Lagrangian particles
  integer(ip)               :: nlapr             ! Current number of properties requested by modules
  integer(ip), parameter    :: mtyla = 100       ! Max # of lagrangian particle types
  integer(ip), parameter    :: pts_advec_narrays=3 !number of arrays of advec(:,:,:) to use, partis uses only advec(:,:,1) and advec(:,:,3). Needed when nastin is off and timsteps are loaded
  
!  character(LEN=255)        :: veloc_field_prefix ='' !template for velocity field files for each timestep, i.e. for files like "V/VELOC00000001.mpio.bin" this wuld be 'V/VELOC'
!  integer(4)                :: veloc_field_Ndigits !number of digits in the template of velocity field filename
!  integer(ip)               :: veloc_field_Nsteps  !number of timesteps in the velocity field, or number of files. 
  integer(ip)                :: veloc_field_id = -1 !field id of the velocity field

  type typlatyp
     integer(ip)            :: kfl_exist         ! If type exists
     integer(ip)            :: kfl_modla         ! Transport model
     integer(ip)            :: kfl_grafo         ! Gravity
     integer(ip)            :: kfl_buofo         ! Buoyancy
     integer(ip)            :: kfl_drafo         ! Drag force model
     integer(ip)            :: kfl_brown         ! Brownian force
     integer(ip)            :: kfl_extfo         ! External force
     integer(ip)            :: kfl_saffm         ! Saffman force
     integer(ip)            :: kfl_schem         ! Time integration scheme
     integer(ip)            :: kfl_tstep         ! Time step strategy
     real(rp)               :: denpa             ! Density particle
     real(rp)               :: spher             ! Particle sphericity
     real(rp)               :: diame             ! Particle diameter
     real(rp)               :: calor             ! Calorific capacity
     real(rp)               :: emisi             ! Particle radiation emissivity
     real(rp)               :: scatt             ! Particle radiation scattering factor
     real(rp)               :: diffu             ! Particle diffusion coefficient [m^2/s]
     real(rp)               :: dtime             ! Time step
     real(rp)               :: safet             ! Safety factor
     real(rp)               :: chale             ! Characteristic length
     real(rp)               :: prope(mlapr)      ! Default value of properties for each type
     integer(ip)            :: prova(mlapr)      ! ID of variable in properties vector
  end type typlatyp
  type latyp
     integer(ip)            :: ilagr             ! Absolute ID
     integer(ip)            :: itype             ! Type of particle
     integer(ip)            :: kfl_exist         ! If I have it
     integer(ip)            :: ielem             ! Last interpolation element
     integer(ip)            :: iboun             ! Boundary element where particle is deposited
     integer(ip)            :: ittim             ! Time iteration
     integer(ip)            :: boundary_set      ! Boundary set where wall intersection
     real(rp)               :: coord(3)          ! Coordinates
     real(rp)               :: veloc(3)          ! Velocity
     real(rp)               :: accel(3)          ! Acceleration
     real(rp)               :: coord_k(3)        ! Coordinates at time k
     real(rp)               :: coord_km1(3)      ! Coordinates at time k-1
     real(rp)               :: v_fluid_k(3)      ! Fluid velocity at time k
     real(rp)               :: v_fluid_km1(3)    ! Fluid velocity at time k-1
     real(rp)               :: v_fluid_km2(3)    ! Fluid velocity at time k-2
     real(rp)               :: acced(3)          ! Drag acceleration
     real(rp)               :: accee(3)          ! External acceleration
     real(rp)               :: acceg(3)          ! Gravity/buoyancy acceleration
     real(rp)               :: stret             ! Stretching factor
     real(rp)               :: t_inject	       ! Time particle was injected
     real(rp)               :: t                 ! Time
     real(rp)               :: dt_k              ! Time step: t^k+1-t^k
     real(rp)               :: dt_km1            ! Time step: t^k  -t^k-1
     real(rp)               :: dt_km2            ! Time step: t^k-1-t^k-2
     real(rp)               :: dtg               ! Guessed time step
     real(rp)               :: Cd                ! Drag coefficient
     real(rp)               :: Re                ! Reynolds number
     real(rp)               :: Stk(2)            ! Stokes number (instantaneous & effective) 
     real(rp)               :: dista             ! Distance
     real(rp)               :: coord1d           ! 1D coordinate
     real(rp)               :: sign              ! sign used for 1D coordinate
     real(rp)               :: prope(mlapr)      ! Various properties (managed by modules)
  end type latyp
  integer(ip)               :: nlagr
  type(latyp),    pointer   :: lagrtyp(:)
  type(typlatyp)            :: parttyp(mtyla)    ! Particle types

  !------------------------------------------------------------------------
  !
  ! Parameters
  !
  !------------------------------------------------------------------------

  integer(ip), parameter :: kfl_version_pts = 1 ! 0=old, 1=new
  integer(ip), parameter   ::      &
       PTS_OUTFLOW_CONDITION  = 0, &
       PTS_WALL_CONDITION     = 1, &
       PTS_SLIP_CONDITION     = 3, &
       PTS_BOUNCING_CONDITION = 4
  
  integer(ip), parameter   ::      &
       PTS_PARTICLE_EXISTS    = -1,&
       PTS_PARTICLE_HITS_WALL = -2,&
       PTS_PARTICLE_ZERO_TIME = -3,&
       PTS_PARTICLE_OUTFLOW   = -4

  integer(ip), parameter   :: &
       mpala     = 20,        & ! Max # parameters for lagrangian particles
       pts_minj  = 100,       & ! Max # injections
       mvarp_pts = 70           ! Max # of postprocess variables
  
  integer(ip), parameter   :: &
       lun_resul_pts = 1725,  & ! Result file unit
       lun_oudep_pts = 1727,  & ! Deposition map file unit
       lun_depsu_pts = 1729     ! Deposition surface file unit

  character(150)           :: &
       fil_resul_pts,         & ! Result file 
       fil_oudep_pts,         & ! Deposition map file
       fil_depsu_pts            ! Deposition surface file

  real(rp),    parameter   :: &
       zepts = epsilon(1.0_rp)   ! Zero

  !------------------------------------------------------------------------
  !
  ! Physical problem: read in pts_reaphy
  !
  !------------------------------------------------------------------------
  !--BEGIN REA GROUP

  real(rp)                 :: &
       dtmin_pts,             &      ! Minimum time step 
       dimin_pts,             &      ! Minimum distance to wall under which particle is deposited
       mean_free_path_pts            ! Mean free path of medium

  !------------------------------------------------------------------------
  !
  ! Numerical problem: read in pts_reanut
  !
  !------------------------------------------------------------------------

  integer(ip)              :: &
       kfl_adapt_pts,         &      ! Adaptive time step
       kfl_usbin_pts,         &      ! Use element bin for element search
       kfl_order_pts,         &      ! Order of integration of AB
       kfl_walld_pts                 ! Activate walld calculation for corss wall test
  real(rp)                 :: & 
       gamma_pts,             &      ! Lagrangian particles: Initial time of injection
       beta_pts,              &      ! Lagrangian particles: Time period of injection
       chale_pts,             &      ! Characteristic length for adaptive time step
       safet_pts,             &      ! Safety factor
       dtime_pts                     ! Particle time step
  
  !------------------------------------------------------------------------
  !
  ! Output and Postprocess: read in pts_reaous
  !
  !------------------------------------------------------------------------

  type, BIND(C) :: dbparams
     integer(C_INT) :: kfl_db_deep            !number of cube levels in the daba base
     character  :: kfl_db_url_conn(150)       !url database connection data
     integer(C_INT) :: port                   !port connection
     integer(C_INT) :: maxInserts             !number maximum of parall inserts
     integer(C_INT) :: numParticlesBatch      !number maximum of parall inserts    
  end type dbparams

  type (dbparams) dbSettings

  integer(ip)              :: &     
       kfl_posla_pts,         &              ! Lagrangian particles: Postprocess 
       kfl_oudep_pts,         &              ! Output of deposition map
       kfl_depos_surface_pts, &              ! Postprocess deposition surface
       kfl_oufre_pts,         &              ! Output frequency
       kfl_dbfre_pts,         &              ! Output DB frequency
       kfl_exacs_pts                         ! Exact solution
  integer(ip)               :: kfl_ptsres_binary     ! save pts.res in binary
  integer(ip)               :: kfl_ptsres_split      ! save separate pts.res per timestep

  character(100)           :: wordsdb(maxwp) !database load word

  !------------------------------------------------------------------------
  !
  ! Physical problem: read in pts_reabcs
  !
  !------------------------------------------------------------------------

  integer(ip)                   :: &  
       kfl_injla_pts(pts_minj) ,   &        ! Lagrangian particles: Injection model
       kfl_injty_pts(pts_minj) ,   &        ! Lagrangian particles: type on which injector is applied
       kfl_random_pts(pts_minj),   &        ! Lagrangian particles: random mode
       kfl_injve_pts,              &        ! Lagrangian particles: Velocity injection model
       kfl_boundary_injection ,    &        ! Inject particles from boundary
       codbo_pts(pts_minj),        &        ! Code of the boundary to inject particles
       injector_particle_distribution(pts_minj), & ! Code for the particle distribution, 1-uniform in cartesian coords, 0-uniform in polar coords
       injector_npts_asis(pts_minj)                ! 1 - Inject exactly the requested number of particles, not square, cube or whatever else, 0-default behaviour
  real(rp)                      :: & 
       tinla_pts,                  &        ! Lagrangian particles: Initial time of injection
       tfila_pts,                  &        ! Lagrangian particles: Final time of injection
       tpela_pts,                  &        ! Lagrangian particles: Time period of injection
       parla_pts(pts_minj,mpala),  &        ! Parameters for the injection
       parla2_pts(mpala)
  integer(ip),   pointer  ::  &
       kfl_fixbo_pts(:)              ! Element boundary fixity     
  real(rp),      pointer  ::  &
       bvnat_pts(:,:)                ! Natural bc values
  type(bc_bound),pointer  :: &     
       tbcod_pts(:)                  ! Boundary code type  

  type typ_injector_coord
     integer(ip)              :: number_particles
     !!!!!!!real(rp),        pointer :: coord_particles(:,:)  !para hacer el broadcast de esto cambio la estructura
     real(rp),        pointer :: coord_particles(:)
     real(rp),        pointer :: parameters(:)
  end type typ_injector_coord
  type(typ_injector_coord) :: injection_pts(mpala)

  !--END REA GROUP

  !------------------------------------------------------------------------
  !
  ! Others
  !
  !------------------------------------------------------------------------

  integer(ip)              :: &
       nlacc_pts,             &        ! Number of accumulated particles
       kfl_depos_pts,         &        ! Postprocess deposition map
       kfl_resid_pts,         &        ! Postprocess residence time
       ntyla_pts,             &        ! # Number of effective particle types (max type used)     
       number_types_pts,      &        ! # Number of particle types
       ninj_pts,              &        !
       kfl_slip_wall_pts,     &        ! If there are slip walls
       kfl_bouncing_wall_pts, &        ! If there are bouncing walls
       kfl_injec,             &        ! If particles have been injected
       nlagr_local_pts,       &        ! Number of existing local particles
       nlagr_free_pts,        &        ! Number of free positions in type 
       nlagr_existing_pts,    &        ! # total number of existing particles in all subdomains _1
       nlagr_non_migrating_pts,&       ! # particles going from one subdomain to another        _2
       nlagr_going_out_pts,   &        ! # particles going out of the computational domain      _3
       nlagr_zero_time_pts,   &        ! # particles that disappear because of zero time step   _4
       nlagr_deposited_pts,   &        ! # particles that deposited without nboun               _5
       nlagr_hits_wall_pts,   &        ! # particles that deposited on wall
       nvarp_pts,             &        ! # Number of postprocess variables
       nvard_pts                       ! # Number of deposition variables
  real(rp)                 :: &   
       cutla_pts,             &        ! Current injection time
       xc,yc,zc,nx,ny,nz,     &        ! Center of the injected particles
       vv, sigma,             &        ! Velocity and sigma of injector
       time_transport_pts              ! CPU time for transport
  integer(ip), pointer     :: &
       lboue_pts(:),          &        ! List of element touching boundaries
       kfl_fixno_walld_slip_pts(:,:),& ! Fixity wall distance to slip boundaries
       kfl_fixno_walld_bouncing_pts(:,:),& ! Fixity wall distance to bouncing boundaries
       permu_nlagr_pts(:)              ! Permutation
  real(rp),    pointer     :: &
       resid_pts(:,:),        &        ! Residence time
       depoe_pts(:,:),        &        ! Deposition map over elements
       depob_pts(:,:),        &        ! Deposition map over boundaries
       defor_pts(:,:),        &        ! Velocity deformation tensor
       hleng_pts(:),          &        ! Element characteristic length
       bouno_pts(:,:),        &        ! Boundary outwards normal
       walld_slip_pts(:),     &        ! Wall distance to slip boundaries
       walld_bouncing_pts(:), &        ! Wall distance to bouncing boundaries
       friction_pts(:)                 ! Friction coefficient
  type(i1p),   pointer     :: &
       leleboun_pts(:)                 ! List of element boundaries
  type(latyp)              :: &
       lagrtyp_init                    ! Initialization particle
  type(latyp), pointer     :: &
       lagrtyp_tmp(:)                  ! Copy of particle type
  integer(ip)              :: &
       particles_sent,        &        ! Number of sent particles
       particles_recv,        &        ! Number of received particles
       comm_loops_pts                  ! Number of communicaiton loops
  integer(ip)              :: &
       kfl_rstar_pts                   ! Independent restart of particles
  !
  ! Postprocess
  !
  logical(lg)              :: &
       postprocess_var_pts(mvarp_pts)  ! Postprocess variables
  character(5)             :: &
       postprocess_name_pts(mvarp_pts) ! Postprocess variable names
  logical(lg)              :: &
       deposition_var_pts(mvarp_pts)   ! Deposition variables
  logical(lg)              :: &
       migrated_variables_pts(mvarp_pts) ! Variables to migrate
  integer(ip)              :: &
       number_migrated_variables_pts   ! Variables to migrate
  integer(ip),  pointer    :: &
       postprocess_list_pts(:)         ! List of postprocess variables
  integer(ip),  pointer    :: &
       deposition_list_pts(:)          ! List of deposition variables

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-10-03
  !> @brief   Find a variable ID
  !> @details According to a string find the variable ID
  !> 
  !-----------------------------------------------------------------------

  integer(ip) function pts_name_to_variable_number(wname)

    implicit none
    character(*), intent(in) :: wname
    integer(ip)              :: ii

    pts_name_to_variable_number = 0
    do ii = 1,mvarp_pts
       if( trim(wname) == trim(postprocess_name_pts(ii)) ) then
          pts_name_to_variable_number = ii
          return
       end if
    end do
    
  end function pts_name_to_variable_number
  
end module def_partis
