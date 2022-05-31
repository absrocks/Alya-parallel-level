!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    mod_AMR.f90
!> @author  houzeaux
!> @date    2020-03-07
!> @brief   Adaptive mesh refinement
!> @details All tools for adaptive mesh refinement
!-----------------------------------------------------------------------

module def_AMR

  use def_kintyp_basic, only : ip,rp
  use def_coupli,       only : typ_color_coupling
  implicit none

  private

  !----------------------------------------------------------------------
  !
  ! Background mesh
  !
  !----------------------------------------------------------------------

  integer(ip), parameter   :: AMR_BACKGROUND_BIN           = 0
  integer(ip), parameter   :: AMR_BACKGROUND_OCTREE        = 1
  integer(ip), parameter   :: AMR_BACKGROUND_ORIGINAL_MESH = 2
  integer(ip), parameter   :: AMR_BACKGROUND_OCTBIN        = 3
  
  !----------------------------------------------------------------------
  !
  ! Interpolation strategies
  !
  !----------------------------------------------------------------------

  integer(ip), parameter   :: AMR_ON_NODE_ELEMENT_INTERPOLAITON = 0
  integer(ip), parameter   :: AMR_ON_NODE_NEAREST_NODE          = 1
  integer(ip), parameter   :: AMR_ON_ELEMENT                    = 2
  integer(ip), parameter   :: AMR_BOUNDARY                      = 3
  
  !----------------------------------------------------------------------
  !
  ! Input variables
  !
  !----------------------------------------------------------------------

  integer(ip)              :: kfl_amr                 ! AMR
  integer(ip)              :: kfl_amr_post            ! Postprocess technique for AMR
  integer(ip)              :: kfl_amr_freq            ! Frequency
  integer(ip)              :: nelem_amr               ! Target number of elements
  integer(ip)              :: maxit_amr               ! Maximum number of iterations
  integer(ip)              :: kfl_amr_varia           ! Variable
  integer(ip)              :: kfl_amr_background      ! Background mesh
  integer(ip)              :: limit_amr_background    ! Size background mesh
  integer(ip)              :: boxes_amr_background    ! Size background mesh
  integer(ip)              :: kfl_size_amr            ! Size strategy
  real(rp)                 :: min_size_amr            ! Min size fot eh mesh
  real(rp)                 :: max_size_amr            ! Max size for the mesh
  
  !----------------------------------------------------------------------
  !
  ! Local variables
  !
  !----------------------------------------------------------------------

  type(typ_color_coupling), target :: coupling_AMR_npoin         ! Node volume interpolation
  type(typ_color_coupling), target :: coupling_AMR_npoin_nearest ! Nereast node vale
  type(typ_color_coupling), target :: coupling_AMR_nelem         ! Element volume interpolation
  type(typ_color_coupling), target :: coupling_AMR_nboun         ! Boundary volume interpolation
 
  integer(ip)              :: num_amr                    ! Number of mesh refinement steps
  integer(ip)              :: npoin_new
  integer(ip)              :: nelem_new
  integer(ip)              :: nboun_new

  public :: AMR_BACKGROUND_BIN          
  public :: AMR_BACKGROUND_OCTREE       
  public :: AMR_BACKGROUND_ORIGINAL_MESH
  public :: AMR_BACKGROUND_OCTBIN      

  public :: num_amr
  public :: kfl_amr
  public :: kfl_amr_post
  public :: kfl_amr_freq
  public :: nelem_amr
  public :: maxit_amr
  public :: kfl_amr_varia
  public :: kfl_amr_background
  public :: limit_amr_background
  public :: boxes_amr_background
  public :: kfl_size_amr 
  public :: min_size_amr 
  public :: max_size_amr 
  
  public :: coupling_AMR_npoin
  public :: coupling_AMR_npoin_nearest
  public :: coupling_AMR_nelem
  public :: coupling_AMR_nboun

  public :: npoin_new
  public :: nelem_new
  public :: nboun_new

end module def_AMR
!> @}

