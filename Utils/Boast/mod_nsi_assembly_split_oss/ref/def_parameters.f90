module def_parameters

  use def_kintyp
  !
  ! Parameters that do not change
  !
  integer(ip), parameter  :: ndime = 3               ! Spatial dimension
  real(rp),    parameter  :: zeror = epsilon(1.0_rp) ! Almost zero
  integer(ip), parameter  :: TET04 = 30              ! 3D
  integer(ip), parameter  :: TET10 = 31              ! 3D
  integer(ip), parameter  :: PYR05 = 32              ! 3D
  integer(ip), parameter  :: PYR14 = 33              ! 3D
  integer(ip), parameter  :: PEN06 = 34              ! 3D
  integer(ip), parameter  :: PEN15 = 35              ! 3D
  integer(ip), parameter  :: PEN18 = 36              ! 3D
  integer(ip), parameter  :: HEX08 = 37              ! 3D
  integer(ip), parameter  :: HEX20 = 38              ! 3D
  integer(ip), parameter  :: HEX27 = 39              ! 3D
  integer(ip), parameter  :: HEX64 = 40              ! 3D
  integer(ip), parameter  :: SHELL = 51              ! 3D shell element
  integer(ip), parameter  :: BAR3D = 52              ! 3D bar element
  !
  ! Size of vector: "VECTOR_SIZE" elements are assembled at the same time
  !
  integer(ip), parameter  :: VECTOR_SIZE = 4         ! Size for vectorization
  !
  ! Nastin parameters
  !
  integer(ip),   parameter :: &
       NSI_INCREMENTAL_PROJECTION          =  2, &
       NSI_PREDICTOR_CORRECTOR             =  3, &
       NSI_BLOCK_GAUSS_SEIDEL              =  4, &
       NSI_MOMENTUM                        =  2, &
       NSI_CONTINUITY                      =  3, &
       NSI_MOMENTUM_AND_CONTINUITY         =  1, &
       NSI_INCOMPRESSIBLE                  =  0, &
       NSI_COMPRESSIBLE                    =  1, &
       NSI_LOW_MACH                        =  3, &
       NSI_ANALYTICAL_HYDROSTATIC_PRESSURE =  1, &
       NSI_PDE_HYDROSTATIC_PRESSURE        =  2, &
       NSI_GALERKIN                        = -1, &
       NSI_ASGS                            =  0, &
       NSI_OSS                             =  1, &
       NSI_SPLIT_OSS                       =  2, &
       NSI_ALGEBRAIC_SPLIT_OSS             =  3, &
       NSI_CONVECTION_NON_CONSERVATIVE     =  0, &
       NSI_CONVECTION_CONSERVATIVE         =  1, &
       NSI_CONVECTION_SKEW                 =  2, &
       NSI_CONVECTION_EMAC                 =  3
  !
  ! Internal variables
  !
  logical(lg)              :: &
       NSI_MONOLITHIC,        &      ! Monolithic algorithm
       NSI_SCHUR_COMPLEMENT,  &      ! Schur complement algorithm
       NSI_FRACTIONAL_STEP

end module def_parameters
