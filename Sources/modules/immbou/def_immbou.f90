module def_immbou
  !------------------------------------------------------------------------
  !****f* Immbou/def_immbou
  ! NAME 
  !    def_immbou
  ! DESCRIPTION
  !    Heading for the Immbou routines
  ! USES
  ! USED BY
  !    Almost all
  !***
  !------------------------------------------------------------------------
  use def_kintyp
  use def_master, only : netyp

  !------------------------------------------------------------------------
  !
  ! Parameters
  !
  !------------------------------------------------------------------------

  integer(ip), parameter :: &
       lun_mshib_ibm = 2110, lun_resib_ibm = 2111, lun_mshi2_ibm = 2112, &
       lun_resi2_ibm = 2113, lun_outpu_ibm = 2114
 
  real(rp),      parameter :: &
       zeibm = epsilon(1.0_rp)
 
  integer(ip), parameter      :: &
       IBM_EMBEDDED      = 0
!--BEGIN REA GROUP
  !------------------------------------------------------------------------
  !
  ! Physical problem: read in ibm_readim
  !
  !------------------------------------------------------------------------

  
  !------------------------------------------------------------------------
  !
  ! Physical problem: read in ibm_reaphy
  !
  !------------------------------------------------------------------------

  integer(ip)                           :: &
       kfl_timei_ibm,                      & ! Existence of time behavior
       kfl_rotib_ibm,                      & ! IB enable rotation
       kfl_linib_ibm,                      & ! IB enable linear motion
       nstro_ibm(3),                       & ! IB start rotation motion
       nstli_ibm(3),                       & ! IB start linear motion
       kfl_mvext_ibm,                      & ! IB move exterior idem COG in x & y
       kfl_ralei_ibm,                      & ! IB enable Raleigh damping
       nstra_ibm,                          & ! IB start eliminating Raleigh damping
       nenra_ibm,                          & ! IB end eliminating Raleigh damping
       kfl_staib_ibm,                      & ! IB start
       kfl_colli_ibm,                      & ! IB collisions
       kfl_grafo_ibm,                      & ! Gravity force
       kfl_catfo_ibm,                      & ! Catamaran force
       kfl_buofo_ibm,                      & ! Buoyancy force
       kfl_drafo_ibm,                      & ! Drag force
       kfl_extfo_ibm                         ! External force
  real(rp)                              :: &
       staib_ibm,                          & ! Starting time
       spher_ibm,                          & ! Sphericity 
       xline_ibm(3),                       & ! Linear motion factor
       xrota_ibm(3),                       & ! Angular motion factor
       ralei_ibm                             ! Raleigh damping factor

  !------------------------------------------------------------------------
  !
  ! Numerical problem: read in ibm_reanut
  !
  !------------------------------------------------------------------------

  integer(ip)                           :: &
       kfl_inter_ibm                         ! Method of interpolation/ALE
  integer(ip)                           :: &
       kfl_nforc_ibm                         ! Numerical force calculation
  real(rp)                              :: &
       safet_ibm,                          & ! Safety factor
       beta_ibm,                           & ! Beta for Newmark
       gamma_ibm                             ! Gamma for Newmark

  !------------------------------------------------------------------------
  !
  ! Output and Postprocess: read in ibm_reaous
  !
  !------------------------------------------------------------------------

!--END REA GROUP
  !------------------------------------------------------------------------
  !
  ! Others
  !
  !------------------------------------------------------------------------

  type bxtyp
     integer(ip)          :: nnode             ! Number of nodes in the box 
     integer(ip), pointer :: nodes(:)          ! Nodes of the box
  end type bxtyp 

  integer(ip),allocatable :: r_dom_2(:),c_dom_2(:)

  type(bxtyp), pointer    :: boxes_ibm(:,:,:) 
  integer(ip)             :: nubox_ibm         ! Number of boxes

  real(rp),    pointer    :: xmima_ibm(:,:)    ! Boundary box of subdomains including their neibors  

  real(rp),    pointer    :: lndib_ibm(:)      ! Signed distances to the nearest faces
  integer(ip), pointer    :: lnfib_ibm(:)      ! List of valid fringe nodes (for remove deformed fringe nodes)
  integer(ip)             :: ncoll_ibm         ! number of collitions
  integer(ip), pointer    :: colli_ibm(:,:)    ! Particle List of ollitions


  !
  ! LELCH
  !
  integer(ip), pointer    :: lelch_ibm(:)      ! Immbou element characteristics
  integer(ip), pointer    :: lnoch_ibm(:)      ! Immbou node characteristics


  !
  ! Internal variables
  !
  integer(ip)             :: kfl_embed_ibm     ! If embedded exist
  integer(ip)             :: kfl_diric_ibm     ! Existence of Dirichlet IB
  integer(ip)             :: kfl_force_ibm     ! Existence of Force IB
  integer(ip)             :: kfl_stead_ibm     ! Steady-state has been reached 
  real(rp)                :: dtcri_ibm         ! Critical time step
  real(rp)                :: bacdt_ibm         ! Critical time step backup
  real(rp)                :: dtinv_ibm         ! 1/dt
  real(rp)                :: dtime_ibm         ! dt used for collisions
  real(rp)                :: cutim_ibm         ! Current time
  !
  ! Walls
  !
  integer(ip)             :: nwaib
  type waltypib
     integer(ip)          :: nboun             ! Number of boundaries
     integer(ip)          :: npoin             ! Number of nodes
     integer(ip)          :: numbe             ! Global numbering
     integer(ip), pointer :: lnodb(:,:)        ! Boundary connectivity
     integer(ip), pointer :: ltypb(:)          ! Type of boundaries
     integer(ip), pointer :: lninv(:)          ! Permutation
     real(rp),    pointer :: coord(:,:)        ! Coordinates
     real(rp),    pointer :: bouno(:,:)        ! Normal
     real(rp),    pointer :: fabox(:,:,:)      ! Tree structure
     real(rp),    pointer :: sabox(:,:,:)      ! ->
     integer(ip), pointer :: blink(:)          !
     integer(ip), pointer :: stru2(:)          !
     real(rp),    pointer :: ldist(:)          ! 
     real(rp)             :: bobox(3,2)        ! <- Bounding box
     type(netyp), pointer :: lnele(:)          ! Node to element graph

     real(rp)             :: cotim             ! Estimated time of collision
  end type waltypib 
  type(waltypib), pointer :: twall_ibm(:)      ! Wall structure

end module def_immbou
