module def_dodeme
  !-----------------------------------------------------------------------
  !****f* Dodeme/def_dodeme
  ! NAME
  !    def_dodeme
  ! DESCRIPTION
  !    Dodeme module
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_master, only       :  netyp
  !
  ! File units
  !
  integer(ip),   parameter   :: lun_pdata_dod = 5701
  integer(ip),   parameter   :: lun_outpu_dod = 5702
!  real(rp),      parameter   :: increm = epsilon(1.0_rp) 
   !
  ! Parameters
  !
  integer(ip),     parameter :: minte_dod      = 4 !!OJO antes 3
  integer(ip),     parameter :: DOD_CHIMERA    = 1
  integer(ip),     parameter :: DOD_PATCH      = 2               
  integer(ip),     parameter :: DOD_PRESCRIBED = 3
  integer(ip),     parameter :: DOD_HOLED_PATCH = 4               
  integer(ip),     parameter :: DOD_FREE       = 0
  integer(ip),     parameter :: DOD_SOLID      = 2
  !
  ! Interface arrays
  !
  integer(ip)                :: ictop_dod(minte_dod)     ! Interface topology
  integer(ip)                :: ismoo_dod(minte_dod)     ! Smoothing interface is needed
  integer(ip), pointer       :: intyp_dod(:,:)           ! Interface type

  !
  ! Subdomain topologies
  !
  integer(ip),     pointer   :: ihole_dod(:)             ! Hole cutting is necessary
  integer(ip),     pointer   :: ipatc_dod(:)             ! This is a patch
  integer(ip),     pointer   :: ipres_dod(:)             ! This is a prescribed interface
  !
  ! Extension nodes
  !
  integer(ip)                :: number_fringe_nodes
  type typ_extension
     integer(ip)             :: number_candidates
     integer(ip),  pointer   :: lnods(:,:)
     integer(ip),  pointer   :: ltype(:)
     integer(ip),  pointer   :: candidates(:)
     integer(ip),  pointer   :: lelch(:)
     integer(ip),  pointer   :: lmate(:)
     integer(ip),  pointer   :: lelez(:)
     integer(ip),  pointer   :: lesub(:)
  end type typ_extension
  type(typ_extension), pointer :: lpext(:) ! Remove this
  !
  ! Various
  !  
  real(rp)                   :: shape_dod(64)            ! For Elsest
  real(rp)                   :: deriv_dod(192)           ! For Elsest
  integer(ip)                :: nelem_dod                ! Number of extension elements
 ! integer(ip)                :: nsubd                    ! Number of subdomains
  integer(ip)                :: kfl_smoot_dod           ! Type of smooth strategy
  integer(ip),     pointer   :: lsubd_npoin(:)           ! List of subdomain   for nodes
  integer(ip),     pointer   :: lsubd_nelem(:)           ! List of subdomain   for elements
  integer(ip),     pointer   :: lsubd_nboun(:)           ! List of subdomain   for boundaries
  integer(ip),     pointer   :: linvp_npoin(:)           ! Inverse permutation for nodes
  integer(ip),     pointer   :: linvp_nelem(:)           ! Inverse permutation for elements
  integer(ip),     pointer   :: linvp_nboun(:)           ! Inverse permutation for boundaries
  integer(ip),     pointer   :: lmatn_dod(:)             ! List of materials per node   
  integer(ip),     pointer   :: lpoiz_dod(:)             ! List of zones per node   
  integer(ip),     pointer   :: prescribed_boundaries(:) ! Prescribed interface boundaries
  !
  ! Subdomain arrays
  ! 
  type typ_subdomain2
     !
     ! Dimensions
     !
     integer(ip)             :: npoin                    ! local number nodes 
     integer(ip)             :: nelem                    ! local number elements
     integer(ip)             :: nboun                    ! local number boundaries 
     !
     ! Transmission conditions
     ! 
     integer(ip),  pointer   :: lsubd_npoin(:)           ! Node subdomain
     integer(ip),  pointer   :: lsubd_nelem(:)           ! Element subdomain
     integer(ip),  pointer   :: lsubd_nboun(:)           ! Boundary subdomain
     integer(ip),  pointer   :: host_element(:)          ! Host elements for hole nodes
     !
     ! Permutation arrays
     !
     integer(ip),  pointer   :: lnper(:)                 ! Global => local numbering for nodes
     integer(ip),  pointer   :: leper(:)                 ! Global => local numbering for elements
     integer(ip),  pointer   :: lbper(:)                 ! Global => local numbering for boundaries
     !
     ! Geometrical arrays
     !
     real(rp),     pointer   :: coord(:,:)
     integer(ip),  pointer   :: lnods(:,:)
     integer(ip),  pointer   :: lelch(:)
     integer(ip),  pointer   :: lnnod(:)
     integer(ip),  pointer   :: ltype(:)
     integer(ip),  pointer   :: lboch(:)
     integer(ip),  pointer   :: lnodb(:,:)
     integer(ip),  pointer   :: ltypb(:)
     integer(ip),  pointer   :: lboel(:,:)
     integer(ip),  pointer   :: lelbo(:)
     real(rp),     pointer   :: bouno(:,:)
     real(rp)                :: embox(3,2)
     !    
     ! Element-element graph: for recursive dod_gallet
     !
     integer(ip),  pointer   :: pelel(:)
     integer(ip),  pointer   :: lelel(:)
     !    
     ! Element-node graph: for holyfying free nodes surrounded by hole nodes
     !
     integer(ip),  pointer   :: pelpo(:)
     integer(ip),  pointer   :: lelpo(:)
     !    
     ! Boundary-boundary graph: to detect outer boundaries of patches automatically
     !
     integer(ip),  pointer   :: pbobo(:)
     integer(ip),  pointer   :: lbobo(:)
     !    
     ! Node-node graph: for smoothing
     !
     integer(ip),  pointer   :: ppopo(:)
     integer(ip),  pointer   :: lpopo(:)
     !
     ! Node-boundary graph: for dod_boucex, loop over boundaries connected to fringe nodes
     !
     integer(ip),  pointer   :: pbopo(:)
     integer(ip),  pointer   :: lbopo(:)
     !
     ! SKD-Tree structure
     ! 
     real(rp),     pointer   :: fabox(:,:,:)
     real(rp),     pointer   :: sabox(:,:,:)
     integer(ip),  pointer   :: blink(:)
     integer(ip),  pointer   :: stru2(:)
     real(rp),     pointer   :: ldist(:)
     real(rp)                :: bobox(3,2)
     type(netyp),  pointer   :: lnele(:)
     integer(ip)             :: npoin_bou
     integer(ip),  pointer   :: lnodb_bou(:,:)           ! Boundary connectivity using boundary nodes
     real(rp),     pointer   :: coord_bou(:,:)           ! Boundary node coordinates
     !
     ! Extension type
     !
     type(typ_extension), pointer :: extension(:)

  end type typ_subdomain2

  type(typ_subdomain2), pointer :: subdomain(:)
  type(typ_subdomain2), pointer :: current_subdomain
  type(typ_subdomain2), pointer :: neighbor_subdomain
  !
  ! CPU times
  !
  real(rp)                     ::     &
       cpu_dod_copy_mesh,             &
       cpu_dod_graphs,                &
       cpu_dod_kdtree,                &
       cpu_dod_elsest,                &
       cpu_dod_holcut,                & ! =>
       cpu_dod_holcut_marknodes,      & !
       cpu_dod_holcut_inversehole,    & !
       cpu_dod_holcut_markelements,   & ! HOLE CUTTING
       cpu_dod_holcut_holeboundary,   & !
       cpu_dod_holcut_fringenodes,    & !
       cpu_dod_graphs_kdtree,         & ! <=
       cpu_dod_patch_boundaries,      &
       cpu_dod_prescribed_boundaries, &
       cpu_dod_candidate_nodes,       &
       cpu_dod_extension,             &
       cpu_dod_merge_new_elements
  !
  ! Global boundary needed to creat extensions
  !
  integer(ip)                  :: nboun_global
  integer(ip)                  :: npoin_global
  integer(ip)                  :: nbopo_global
  integer(ip),         pointer :: lnodb_global(:,:)
  integer(ip),         pointer :: lboel_global(:)
  integer(ip),         pointer :: lelbo_global(:)
  integer(ip),         pointer :: lboch_global(:)
  integer(ip),         pointer :: ltypb_global(:)
  integer(ip),         pointer :: pbopo_global(:)
  integer(ip),         pointer :: lbopo_global(:)
  integer(ip),         pointer :: nbpoi_global(:)
  integer(ip),         pointer :: ledbo_global(:,:)
  type(i1p),           pointer :: ledge_npoin_global(:)
  
end module def_dodeme
