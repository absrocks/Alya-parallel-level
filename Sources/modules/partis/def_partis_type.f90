module def_partis_type

  use def_kintyp_basic, only : ip,rp
  use def_master,       only : dtime
  use mod_physics,      only : liquid_state
  use mod_interp_tab,   only : typ_lookup_framework
  use mod_interp_tab,   only : typ_tab_coord, typ_lookup_table
  implicit none
  private

  integer(ip), parameter    :: mlapr = 10        ! Maximum number of properties in Lagrangian particles
  integer(ip), parameter    :: mtyla = 100       ! Max # of lagrangian particle types
  !
  ! Particle type
  !
  type typlatyp
     integer(ip)            :: kfl_exist         ! If type exists
     integer(ip)            :: kfl_modla         ! Transport model
     integer(ip)            :: kfl_therm         ! Thermodynamic model
     integer(ip)            :: kfl_heattr_corr   ! Correction of heat transfer due to evaporation
     integer(ip)            :: kfl_mass_pot      ! Mass transfer potential formulation
     integer(ip)            :: kfl_grafo         ! Gravity
     integer(ip)            :: kfl_buofo         ! Buoyancy
     integer(ip)            :: kfl_drafo         ! Drag force model
     integer(ip)            :: kfl_brown         ! Brownian force
     integer(ip)            :: kfl_turbu         ! Turbulent diffusion
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
     real(rp)               :: tursc             ! Turbulent schmidt number
     real(rp)               :: prope(mlapr)      ! Default value of properties for each type
     integer(ip)            :: prova(mlapr)      ! ID of variable in properties vector
     !
     ! Thermodynamic
     !
     real(rp)               :: d_min             ! Thermodynamic model: minimum diameter [m] 
     real(rp)               :: n_drop            ! Number of droplets represented by a parcel
     real(rp)               :: L_vapor           ! Thermodynamic model: enthalpy of vaporization          [ J/kg ] 
     real(rp)               :: cp                ! Thermodynamic model: specific heat capacity            [ J / (kg K) ]        
     real(rp)               :: w                 ! Thermodynamic model: molecular weight                  [ kg / mol ]
     real(rp)               :: T_boiling         ! Thermodynamic model: Boiling temperature               [ K ]
     real(rp)               :: cpcoef_v_chm(6,2) ! Coefficients of NASA polynomial of fuel in the lookup table 
     real(rp)               :: weight_seen       ! Thermodynamic model: weighting factor of seen gas properties for evaluating mean properties, default: 1/3
     type(liquid_state)     :: liq               ! Liquid state
     integer(ip)                        :: kfl_tab_fw              ! index of lookup framework
     type(typ_lookup_framework),pointer :: table_fw                ! lookup framework
     type(typ_tab_coord),   pointer     :: h_scaling_coords(:)     ! control variable discretization for enthalpy scaling
     type(typ_tab_coord),   pointer     :: yc_scaling_coords(:)    ! control variable discretization for progress variable scaling
     type(typ_tab_coord),   pointer     :: table_coords(:)         ! control variable discretization for lookup table
     type(typ_lookup_table),pointer     :: h_scaling_tab        ! enthalpy scaling table
     type(typ_lookup_table),pointer     :: yc_scaling_tab       ! progress variable scaling table
     type(typ_lookup_table),pointer     :: table_tab            ! lookup table
  end type typlatyp
  !
  ! Particle
  !
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
     real(rp)               :: t_inject          ! Time particle was injected
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
     !
     ! Thermodynamic model
     !
     real(rp)               :: tempe_k           ! Temperature at time k
     real(rp)               :: tempe_km1         ! Temperature at time k-1
     real(rp)               :: tempe_km2         ! Temperature at time k-2
     real(rp)               :: mass_k            ! Mass at time k
     real(rp)               :: mass_km1          ! Mass at time k-1
     real(rp)               :: mass_km2          ! Mass at time k-2
     real(rp)               :: diam_k            ! Diameter for postprocessing
     real(rp)               :: Temp_fluid_k      ! Seen temperature for postprocessing
     real(rp)               :: Yvap_fluid_k      ! Seen vapour mass fraction for postprocessing
     !
     ! Redistribution
     !
     integer(ip)            :: mpi_rank          ! MPI rank

   contains

     procedure, pass        :: init              ! Initialize

  end type latyp

  type(latyp),    pointer   :: lagrtyp(:)
  type(typlatyp)            :: parttyp(mtyla)    ! Particle types
  
  public :: mlapr
  public :: latyp
  public :: mtyla
  public :: typlatyp
  public :: lagrtyp
  public :: parttyp
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-07-10
  !> @brief   Initialization
  !> @details Particle initialization
  !> 
  !-----------------------------------------------------------------------

  subroutine init(particle)

    class(latyp) :: particle

    particle % ilagr          = 0
    particle % itype          = 1
    particle % kfl_exist      = 0
    particle % ielem          = 0
    particle % ittim          = 0
    particle % boundary_set   = 0
    particle % iboun          = 0
    particle % coord          = 0.0_rp
    particle % veloc          = 0.0_rp
    particle % accel          = 0.0_rp
    particle % coord_k        = 0.0_rp
    particle % coord_km1      = 0.0_rp
    particle % v_fluid_k      = 0.0_rp
    particle % v_fluid_km1    = 0.0_rp
    particle % v_fluid_km2    = 0.0_rp
    particle % acced          = 0.0_rp
    particle % accee          = 0.0_rp
    particle % acceg          = 0.0_rp 
    particle % stret          = 1.0_rp
    particle % t_inject       = 0.0_rp
    particle % t              = 0.0_rp
    particle % dt_k           = dtime
    particle % dt_km1         = dtime
    particle % dt_km2         = dtime
    particle % dtg            = dtime
    particle % Cd             = 0.0_rp
    particle % Re             = 0.0_rp
    particle % Stk            = 0.0_rp
    particle % dista          = 0.0_rp
    particle % coord1d        = 0.0_rp
    particle % sign           = 1.0_rp
    particle % prope(1:mlapr) = parttyp(1) % prope(1:mlapr)
    particle % tempe_k        = 0.0_rp
    particle % tempe_km1      = 0.0_rp
    particle % tempe_km2      = 0.0_rp
    particle % mass_k         = 0.0_rp
    particle % mass_km1       = 0.0_rp
    particle % mass_km2       = 0.0_rp
    particle % diam_k         = 0.0_rp
    particle % Temp_fluid_k   = 0.0_rp
    particle % Yvap_fluid_k   = 0.0_rp

  end subroutine init

end module def_partis_type
