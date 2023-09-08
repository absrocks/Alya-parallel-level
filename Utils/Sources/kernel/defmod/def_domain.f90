module def_domain

  !-----------------------------------------------------------------------
  !****f* defmod/def_domain
  ! NAME
  !   def_domain
  ! DESCRIPTION
  !   This module is the header of the domain
  !***
  !-----------------------------------------------------------------------
  
  use def_kintyp
  use mod_htable,    only : hash_t    

  !------------------------------------------------------------------------
  ! Parameters and units
  !------------------------------------------------------------------------

  integer(ip)              :: &
       lun_pdata_dom ,        &      ! Domain data file unit
       lun_outpu_dom,         &      ! Output domain file unit
       lun_elsta_dom,         &      ! Elsest statistics
       lun_elmsh_dom,         &      ! Elsest mesh
       lun_elres_dom                 ! Elsest results
  integer(ip),   parameter :: &
       nelty=60,              &      ! # of element types
       mfree=200,             &      ! Max # free surfaces
       mnode_max=64,          &      ! Max # nodes per element
       mfiel=100,             &      ! Maximum number of fields
       interval_funno=1000           ! kfl_funno now has diferent meaning if it is >1000 or not.

  !------------------------------------------------------------------------
  ! Dimensions: read in readim
  !------------------------------------------------------------------------

  integer(ip)              :: &
       kfl_autbo,             &      ! Automatic boundaries
       npoin,                 &      ! # of nodal points
       nelem,                 &      ! # of elements
       necnt,                 &      ! # of contact elements
       nncnt,                 &      ! # of contact nodes
       nboun,                 &      ! # of boundary elements
       nperi,                 &      ! # periodic nodes
       lexis(nelty),          &      ! List of existing elements
       nexis,                 &      ! Number of different element types in the geometry file
       utype,                 &      ! Value of the unique type, -1 if several types
       npoib,                 &      ! # immersed nodes
       nboib,                 &      ! # immersed boundaries
       lexib(nelty),          &      ! List of existing IB elements
       nhang,                 &      ! # hanging nodes
       nimbo,                 &      ! # IB
       nrbod,                 &      ! # RB
       nzone,                 &      ! Number of zones
       nsubd,                 &      ! Number of subdomains
       nmate,                 &      ! # of materials
       nfiel,                 &      ! Number of fields
       kfl_field(7,mfiel),    &      ! Fields dimensions
       mcodb=99                         ! Max # codes

#ifdef NDIMEPAR
  ! much more comfortable this way; you can have a forder unix2d with -DTWODIM in the config.in
#ifdef TWODIM
  integer(ip), parameter   :: ndime = 2  ! # of space dimensions
#else
  integer(ip), parameter   :: ndime = 3  ! # of space dimensions
#endif
#else
  integer(ip)              :: ndime  ! # of space dimensions
#endif

  !------------------------------------------------------------------------
  ! Strategy: read in reastr
  !------------------------------------------------------------------------

  integer(ip), target      :: &
       ngaus(nelty),          &      ! # of Gauss points per element
       ngaib(nelty)                  ! IB: # of Gauss points per element
  integer(ip)              :: &
       lquad(nelty),          &      ! List of quadrature
       kfl_ngrou,             &      ! Strategy to construct groups
       ngrou_dom,             &      ! Groups (for deflated CG)
       ngrou_boxes_coarse,    &      ! Number of coarse boxes for SFC
       ngrou_boxes_fine,      &      ! Number of fines boxes for SFC
       kfl_savda,             &      ! Save element data base
       lquib(nelty),          &      ! IB: List of quadrature
       kfl_geome,             &      ! Geometrical normals should be computed
       kfl_convx,             &      ! What to do with convex nodes
       kfl_frees,             &      ! Freestream criterion
       kfl_extra,             &      ! Extrapolate from boundary to nodes
       npbcs(8),              &      ! # parameters geometrical bcs
       lsbcs(100,8),          &      ! # list of geometrical bcs
       kfl_chege,             &      ! Check geometry
       kfl_naxis,             &      ! Axi-symmetry
       kfl_spher,             &      ! Spherical
       kfl_bouel,             &      ! Boundary-element connectivity
       kfl_divid,             &      ! Divide element into TET04
       curvatureDataField,    &      ! The field that occupies the curved data for mesh division
       curvatureField,        &      ! The field that occupies the curved geometry for mesh division
       materials_nlaye(5),    &      ! Automatic generation of materials
       materials_icode(5),    &      ! Automatic generation of materials
       materials_imate(5)            ! Automatic generation of materials

  real(rp)                 :: &
       xscal(3),              &      ! Geometric scale factors
       trans(3)                      ! Geometric translation factors
  real(rp)                 :: &
       awind,                 &      ! Wind angle (for freestream condition)
       tolan,                 &      ! Tolerance used to define inflow from freestream
       geoan                         ! Geometrical angle

  !------------------------------------------------------------------------
  ! Geometry: read in reageo
  !------------------------------------------------------------------------

  integer(ip), pointer     :: &
       lnods(:,:),            &      ! Interior element connectivity
       ltype(:),              &      ! List of element types
       lesub(:),              &      ! List of element subdomains
       lgaus(:),              &      ! List of element Gauss points
       lnodb(:,:),            &      ! Boundary element connectivity
       ltypb(:),              &      ! List of boundary types
       lboch(:),              &      ! List of boundary characteristics
       lelbo(:),              &      ! List of boudnary elements (old LBOEL)
       lnnod(:),              &      ! Element number of nodes
       lelch(:),              &      ! Element characteristic
       lmate(:),              &      ! Materials (elements)
       lnoch(:),              &      ! List of node characteristics
       lmast(:),              &      ! List of masters for periodicity
       lgrou_dom(:),          &      ! List of groups (deflated CG)
       lperi(:,:)                    ! List of Master/Slave
  real(rp),    pointer     :: &
       coord(:,:),            &      ! Coordinates
       skcos(:,:,:)                  ! Cosine matrices of skew systems
  type(i1p),   pointer     :: &
       lhang(:)                      ! List of hanging nodes
  type(r3p), pointer       :: &
       xfiel(:)                      ! Fields
  type(r1p)                :: &
       time_field(mfiel)             ! Time for Fields for more than 1 step

  !------------------------------------------------------------------------
  ! Sets: read in reaset
  !------------------------------------------------------------------------

  integer(ip)              :: &
       neset,                 &      ! # of element sets
       nbset,                 &      ! # of boundary sets
       nnset                         ! # of node sets
  integer(ip), pointer     :: &
       leset(:),&                    ! List of element sets
       lbset(:),&                    ! List of boundary sets
       lnset(:)                      ! List of node sets

  !------------------------------------------------------------------------
  ! Sets: read in reabcs
  !------------------------------------------------------------------------

  integer(ip)              :: &
       kfl_icodn,             &      ! If codes on nodes
       kfl_icodb,             &      ! If codes on boundaries
       mcono                         ! Max # codes per nodes
  integer(ip),   pointer   :: &
       kfl_codno(:,:),        &      ! Node codes
       kfl_codbo(:)                  ! Boundary codes

  !------------------------------------------------------------------------
  ! Global variables
  !------------------------------------------------------------------------
  !
  ! Dimensions
  !
  integer(ip)              :: &
       nedge,                 &      ! Number of edges
       medge,                 &      ! Maximum number of edges
       mecnt,                 &      ! Total number of contact elements
       mnode,                 &      ! Max # of nodes per element
       mnoga,                 &      ! Max mnode and mgaus
       mnoib,                 &      ! Max # of nodes per IB (boundary)
       mnodi,                 &      ! Max # of nodes per IB (volume)
       mnodb,                 &      ! Max # of nodes per boundary
       mgaus,                 &      ! Max # of Gauss points per element
       mgaib,                 &      ! Max # of Gauss points per IB (surface)
       mgaui,                 &      ! Max # of Gauss points per IB (volume)
       mgaub,                 &      ! Max # of Gauss points per boundary
       mlapl,                 &      ! Max # of llapl(ielty), ielty=1,nelty
       nfacs,                 &      ! Number of faces
       ntens,                 &      ! # of components of symmetric tensors
       niner,                 &      ! # of components of inertial tensor
       nrule,                 &      ! # of integration rule
       nbopo,                 &      ! # of boundary points
       ndimb,                 &      ! # of space dimensions-1
       nzdom,                 &      ! # of nonzero elements in the mesh graph
       nzdom_own,             &      ! # of nonzero elements in the own mesh graph
       nzbou,                 &      ! # of nonzero elements in the boundary mesh graph
       nzsky,                 &      ! # of comp. in the skyline matrix of the graph
       nzsol,                 &      ! = nzdom + 2*nslav: Components due to slaves
       nzmat,                 &      ! = max(nzsol,nzsky,nzexp): Components of A
       nzmbt,                 &      ! Components of B
       nzrhs,                 &      ! = max(nzsol,nzsky,nzexp): Components of RHS
       nzpre,                 &      ! Components of Preconditioner
       neige,                 &      ! Size of eigen value vector
       neiva,                 &      ! Number of eigenvalues
       nzerr,                 &      ! Components of the Error Estimator
       elmin,                 &      ! Element with minimum volume
       elmax,                 &      ! Element with maximum volume
       nzsym,                 &      ! # of nonzero elements in the symmetric mesh graph
       nxexa,                 &      ! # x-element for example
       nyexa,                 &      ! # y-element for example
       nzexa,                 &      ! # z-element for example
       nxexp,                 &      ! # x-partition for example
       nyexp,                 &      ! # y-partition for example
       nzexp,                 &      ! # z-partition for example
       kfl_parex,             &      ! Automatic parallelization of domain
       nbono                         ! Number of boundary nodes
  integer(ip)              :: &
       nzmax,                 &      ! = max(nzsol,nzsky,nzexp): Components of A
       nzrhx,                 &      ! = max(nzsol,nzsky,nzexp): Components of RHS
       nzprx                         ! Components of Preconditioner

  integer(8)               :: &
       memor_dom(2)                  ! Memory counter
  integer(ip),     pointer :: &
       lpoib(:),              &      ! IB List of point IB
       nmatn(:),              &      ! Number of nodes per material
       nodemat(:),            &      ! Nodal material (computed from elements material)
       lbono(:),              &      ! List of nodes attached to boundaries
       lnnob(:),              &      ! Boundary number of nodes
       lesec(:),              &      ! Element set connectivity
       lbsec(:),              &      ! Boundary set connectivity
       lnsec(:)                      ! Node set connectivity
  type(i1p),   pointer     :: &
       lmatn(:)                      ! Materials (nodes)
  !
  ! Domain properties
  !
  integer(ip)              :: &
       bandw_dom,             &      ! Bandwidth
       naved_dom,             &      ! Average number of edges
       nmied_dom,             &      ! Min number of edges
       nmaed_dom                     ! Max number of edges
  real(rp)                 :: &
       vodom,                 &      ! Measure of the domain
       vomin,                 &      ! Minimum element volume
       vomax,                 &      ! Maximum element volume
       voave,                 &      ! Averaged element volume
       profi_dom                     ! Profile
  real(rp),        target  :: &
       xmima(2,3)                    ! Bounding box
  !
  ! Hanging nodes
  !
  integer(ip),    pointer  :: &
       nehan(:)                      ! Number of hanging nodes per elements
  type(i1pi1p),   pointer  :: &
       lehan(:)                      ! List of hanging elements
  !
  ! Reals
  !
  real(rp)                 :: &
       permx(9)                      ! Rotation matrix relating periodic faces
  !
  ! Graph
  !
  integer(ip)              :: &
       mepoi,                 &      ! Max. # of elements by node
       mpopo,                 &      ! Max. # of point-point connectivity
       kfl_crbou,             &      ! If boundary graph has been computed
       kfl_pelel,             &      ! If element graph has been computed
       kfl_domar,             &      ! If geometrical arrays should be recomputed (then goto domarr)
       npoin_ii,              &      ! Interior nodes
       npoin_bb                      ! Boundary nodes
  integer(ip), pointer :: &
       nepoi(:),              &      ! # of neighbor elements
       pelpo(:),              &      ! Pointer node/element connectivity
       lelpo(:),              &      ! List node/element connectivities
       pelpo_2(:),            &      ! Pointer node/element extended connectivity
       lelpo_2(:),            &      ! List node/element extended connectivities
       pelel(:),              &      ! Pointer element/element connectivity
       lelel(:),              &      ! List element/element connectivities
       pelel_2(:),            &      ! Pointer element/element extended connectivity
       lelel_2(:),            &      ! List element/element extended connectivities
       lezdo(:,:,:),          &      ! Lnods to graph array
       lbzdo(:,:,:)                  ! Lnodb to graph array

  integer(ip), pointer     :: &
       r_sol(:),              &      ! Row array for the matrix CSR storage
       c_sol(:)                      ! Column array for the matrix CSR storage
  integer(ip), pointer     :: &
       r_dom(:),              &      ! Row array for the domain CSR storage
       c_dom(:),              &      ! Column array for the domain CSR storage
       r_bou(:),              &      ! Row array for the boundary domain CSR storage
       c_bou(:),              &      ! Column array for the boundary domain CSR storage
       r_dom_aii(:),          &      ! Row array for the domain CSR storage
       permr_aii(:),          &      ! Permutation for Aii
       invpr_aii(:),          &      ! Inverse permutation for Aii
       c_dom_aii(:),          &      ! Column array for the domain CSR storage
       r_dom_aib(:),          &      ! Row array for the domain CSR storage
       c_dom_aib(:),          &      ! Column array for the domain CSR storage
       r_dom_abi(:),          &      ! Row array for the domain CSR storage
       c_dom_abi(:),          &      ! Column array for the domain CSR storage
       r_dom_abb(:),          &      ! Row array for the domain CSR storage
       c_dom_abb(:),          &      ! Column array for the domain CSR storage
       permr_abb(:),          &      ! Permutation for Abb
       invpr_abb(:),          &      ! Inverse permutation for Abb
       r_dom_prec(:),         &      ! Graph for Schur preconditioner
       c_dom_prec(:),         &      ! Graph for Schur preconditioner
       r_dom_own(:),          &      ! Full row graph
       c_dom_own(:),          &      ! Full row graph
       permr_prec(:),         &      ! Permutation for Schur preconditioner
       invpr_prec(:),         &      ! Inverse permutation for Schur preconditioner
       r_sym(:),              &      ! Row array for the matrix CSR symmetric storage
       c_sym(:)                      ! Column array for the matrix CSR symmetric storage

  !------OMPS DOMAIN--------
  type(ompss_domain),   pointer :: ompss_domains(:)
  type(ompss_domain),   pointer :: ompss_boundaries(:)
  !------END OMPSS DOMAIN------

  !
  ! Integer Arrays
  !
  integer(ip), pointer     :: &
       lpoty(:),              &      ! List of point types
       lfacs(:,:),            &      ! List faces
       lboel(:,:),            &      ! List of boundary elements
       lgaib(:),              &      ! IB: number of Gauss points on IB
       lrenn(:),              &      ! Node renumbering
       lfcnt(:,:),            &      ! List of contact element faces
       lncnt(:,:),            &      ! List of contact element nodes
       lessl(:,:)                    ! Groups in parallel
  integer(ip)              :: &
       lnuty(nelty),          &      ! Number of elements for each ielty
       lnuib(nelty)                  ! Number of IB types
  !
  ! Node and boundary codes
  !
  integer(ip), pointer     :: &
       kfl_fixno(:,:),        &      ! General node fixity
       kfl_fixbo(:),          &      ! General boundary fixity
       kfl_funno(:),          &      ! Function number for nodes
       kfl_funbo(:),          &      ! Function number for boundaries
       kfl_fixrs(:),          &      ! Axes
       kfl_geobo(:),          &      ! Geometrical boundary b.c
       kfl_geono(:),          &      ! Geometrical node b.c.
       wallo(:),              &      ! Order of the Index of the nearest wall node
       lpoin(:)                      ! Type of point (for geometrical arrays)
  integer(ip)              :: &
       iffun,                 &      ! If function should be read
       ifloc,                 &      ! If axis should be read
       ifbop,                 &      ! If bc are imposed on boundary nodes
       ifbes                         ! If value should be assigned
  real(rp),    pointer     :: &
       bvess(:,:),            &      ! General node fixity
       bvnat(:,:)                    ! General boundary fixity
  type(bc_nodes), pointer  :: &
       tncod(:)                      ! Node code type
  type(bc_nodes), pointer  :: &
       tgcod(:)                      ! Geometrical Node code type
  type(bc_bound), pointer  :: &
       tbcod(:)                      ! Boundary code type
  !
  ! Real Arrays
  !
  real(rp),    pointer     :: &
       exnor(:,:,:),          &      ! Exterior normal
       vmass(:),              &      ! Lumped mass matrix
       vmasc(:),              &      ! Mass matrix with close rule
       cmass(:),              &      ! Consistent mass matrix
       cmass_weighted(:),     &      ! Consistent weighted mass matrix
       walld(:),              &      ! Distance to the wall
       wallcoor(:,:),         &      ! Coordinates of the to the nearest wall points
       walln(:,:),            &      ! Normal to the wall
       rough(:),              &      ! Roughness
       canhe(:),              &      ! Canopy height
       heiov(:),              &      ! Height over terrain
       canla(:),              &      ! Canopy Leaf Area Density
       ywalb(:),              &      ! Distance to the wall at each boundary for variable wall distance
       ywalp(:),              &      ! Projection of boundary wall distance onto boundary nodes
       ywale(:)                      ! Minimum of ywalb for each element
  !
  ! Element shape functions and derivatives
  !
  type(elm),   pointer     :: &
       elmar(:)                      ! Element data base
  type(elmgp), pointer     :: &
       elmda(:)                      ! Element Gauss point data base
  real(rp), pointer        :: &
       elmda_gpvol(:,:),      &
       elmda_gpcar(:,:,:,:)

  real(rp)                 :: &
       hnatu(nelty)                  ! Natural element length
  integer(ip)              :: &
       lenex(mnode_max+1,nelty),&    ! List of next element node
       ldime(nelty),          &      ! List of element dimensions
       ltopo(nelty),          &      ! List of element topology
       llapl(nelty),          &      ! List of element Laplacian
       lrule(nelty),          &      ! List of element integration rules
       lruib(nelty),          &      ! List of IB integration rules
       lorde(nelty),          &      ! List of element order
       nnode(-nelty:nelty),   &      ! List of element # of nodes
       nface(nelty),          &      ! List of element # of faces
       needg(nelty),          &      ! List of element # of edges
       leedg(2,20,nelty),     &      ! List of element of edges
       iesta_dom,             &      ! Where element starts
       iesto_dom,             &      ! Where element stops
       ibsta_dom,             &      ! Where boundary starts
       ibsto_dom,             &      ! Where boundary stops
       kfl_horde,             &      ! IF high order element exist
       kfl_elcoh,             &      ! Cohesive elements
       mface                         ! Maximum number of faces
  type(i1p),   target      :: &
       ltypf(nelty),          &      ! List of faces type
       nnodf(nelty)                  ! Number of node for each face
  type(i2p),   target      :: &
       lface(nelty)                  ! List of face nodes
  character(7)             :: &
       cenal(nelty)                  ! List of element names (lower case)
  character(13)            :: &
       cetop(nelty)                  ! List of element topology name
  character(5)             :: &
       cenam(nelty)                  ! List of element names
  !
  ! File name
  !
  character(150)           :: &
       fil_outpu_dom                 ! Output domain mesh file
  !
  ! Old mesh data
  !
  integer(ip)              :: &
       npoin_old,             &      ! # of nodal points
       nelem_old,             &      ! # of elements
       nboun_old                     ! # of boundary elements
  !
  ! Special arrays for Level Set reinitialization
  !
  integer(ip)              :: &
       nelwh                         ! # of total elements in the whole mesh
  integer(ip), pointer     :: &
       pefpo(:),              &      ! Pointer node/element (where interface must be sought) connectivity
       lefpo(:),              &      ! List node/element (where interface must be sought) connectivities
       lnuew(:)                      ! Numeration of an element in the whole mesh
  !
  ! Parall service
  !
  integer(ip)              :: &
       nelem_2,               &      ! nelem + fringe elements
       nboun_2,               &      ! nboun + fringe boundaries
       npoin_2,               &      ! npoin + fringe nodes
       npoin_own,             &      ! Own nodes
       npoin_halo,            &      ! Number of nodes up to halo nodes
       npoin_total,           &      ! Number total of points (boundary replicated)
       nelem_total,           &      ! Number total of elements
       nboun_total                   ! Number total of boundaries
  integer(ip), pointer     :: &
       leldo(:,:)                    ! Fringe elements: subdomains/local numbering
  !
  ! Mesh multiplication
  !
  integer(ip),  pointer    :: &
       facel(:,:,:)                  ! List of faces

  !------------------------------------------------------------------------
  !
  ! Mesh structures
  !
  !------------------------------------------------------------------------

  integer(ip),   pointer  :: lnlev(:)        ! List of node level
  integer(ip),   pointer  :: lelev(:)        ! List of element level
  integer(ip),   pointer  :: lblev(:)        ! List of boundary level

  type linno_type
     integer(ip), pointer :: l(:)
     integer(ip)          :: n
  end type linno_type
  type mesh_type
     !
     ! GEOMETRY
     !
     integer(ip)               :: ndime                     ! Dimension number
     integer(ip)               :: ntens                     ! Number Hessian components
     integer(ip)               :: npoin                     ! Number of nodes
     integer(ip)               :: nelem                     ! Number of elements
     integer(ip)               :: nboun                     ! Number of boundaries
     integer(ip)               :: nedge                     ! Number of edges
     integer(ip)               :: mnode                     ! Max number of node per element
     integer(ip)               :: mnodb                     ! Max number of node per boundary
     integer(ip)               :: mgaus                     ! Max number of Gauss points
     integer(ip)               :: medge                     ! Max number of edge per element
     integer(ip)               :: nbopo                     ! Number of boundary nodes
     integer(ip)               :: npoi1                     ! Number of interior node
     integer(ip)               :: npoi2                     ! First own boundary node
     integer(ip)               :: npoi3                     ! Last own boundary node
     integer(ip)               :: npoin_2                   ! Number of nodes including halo nodes
     integer(ip)               :: nelem_2                   ! Number of elements including halo nodes
     integer(ip)               :: nboun_2                   ! Number of boundaries including halo nodes
     integer(ip)               :: npoin_own                 ! Own nodes = npoi3
     integer(ip)               :: npoin_halo                ! Number of nodes up to halo nodes
     integer(ip),      pointer :: npoin_par(:)
     integer(ip),      pointer :: nelem_par(:)
     integer(ip),      pointer :: nboun_par(:)
     integer(ip)               :: npoin_total
     integer(ip)               :: nelem_total
     integer(ip)               :: nboun_total
     integer(ip),      pointer :: leinv_loc(:)              ! NELEM
     integer(ip),      pointer :: lnods(:,:)                ! NELEM
     integer(ip),      pointer :: ltype(:)                  ! NELEM
     integer(ip),      pointer :: lelch(:)                  ! NELEM
     integer(ip),      pointer :: lnnod(:)                  ! NELEM
     integer(ip),      pointer :: lesub(:)                  ! NELEM
     integer(ip),      pointer :: lmate(:)                  ! NELEM
     integer(ip),      pointer :: lgaus(:)                  ! NELEM
     integer(ip),      pointer :: lbinv_loc(:)              ! NBOUN
     integer(ip),      pointer :: lnodb(:,:)                ! NBOUN
     integer(ip),      pointer :: lboel(:,:)                ! NBOUN
     integer(ip),      pointer :: lelbo(:)                  ! NBOUN
     integer(ip),      pointer :: ltypb(:)                  ! NBOUN
     integer(ip),      pointer :: lboch(:)                  ! NBOUN
     integer(ip),      pointer :: lnnob(:)                  ! NBOUN
     real(rp),         pointer :: coord(:,:)                ! NDIME,NPOIN: node coordinates
     integer(ip),      pointer :: lnoch(:)                  ! NPOIN: List of node characteristics
     integer(ip),      pointer :: lmast(:)                  ! NPOIN: List of master nodes
     integer(ip),      pointer :: lpoty(:)                  ! NPOIN: list of boudnary nodes
     !
     ! Derived arrays
     !
     integer(ip),      pointer :: lninv_loc(:)              ! NPOIN
     integer(ip),      pointer :: edge_to_node(:,:)         ! NEDGE
     integer(ip),      pointer :: local_to_global_edge(:,:) ! NEDGE
     integer(ip),      pointer :: ledgs(:,:)                ! NELEM: element edge connectivity  (see LNODS)
     integer(ip),      pointer :: ledgb(:,:)                ! NBOUN: boundary edge connectivity (see LNODB)
     integer(ip),      pointer :: lnned(:)                  ! NELEM: Number of edge by element  (see LNNOD)
     integer(ip),      pointer :: lnneb(:)                  ! NBOUN: Number of edge by element  (see LNNOB)
     !
     ! SETS
     !
     integer(ip),      pointer :: leset(:)                  ! NELEM
     integer(ip),      pointer :: lbset(:)                  ! NBOUN
     integer(ip),      pointer :: lnset(:)                  ! NPOIN
     !
     ! BOUNDARY CONDITIONS
     !
     integer(ip),      pointer :: kfl_codno(:,:)            ! NPOIN
     integer(ip),      pointer :: kfl_codbo(:)              ! NBOUN
     !
     ! GROUPS
     !
     integer(ip),      pointer :: lgrou_dom(:)              ! NPOIN
     !
     ! INTERPOLATION
     !
     type(linno_type), pointer :: linno(:)                  ! NPOIN
     !
     ! Graphs
     !
     integer(ip)               :: nzdom                     ! Size of node graph
     integer(ip)               :: nzsym                     ! Size of node symmetric graph
     integer(ip)               :: nzedg                     ! Size of edge graph
     integer(ip)               :: nzelm_2                   ! Size of element graph (including halos)
     integer(ip)               :: nzdom_own                 ! Size of node graph of own nodes
     integer(ip)               :: nzdom_ell                 ! Size of ell graph
     integer(ip),      pointer :: c_dom(:)                  ! Node CSR graph
     integer(ip),      pointer :: r_dom(:)                  ! Node CSR graph
     integer(ip),      pointer :: coo_rows(:)               ! Node COO graph
     integer(ip),      pointer :: coo_cols(:)               ! Node COO graph
     integer(ip),      pointer :: ell_cols(:,:)             ! Node ELL graph
     integer(ip),      pointer :: c_sym(:)                  ! Symmetric node graph
     integer(ip),      pointer :: r_sym(:)                  ! Symmetric node graph
     integer(ip),      pointer :: c_edg(:)                  ! Edge graph
     integer(ip),      pointer :: r_edg(:)                  ! Edge graph
     integer(ip),      pointer :: c_elm_2(:)                ! Element graph including halos
     integer(ip),      pointer :: r_elm_2(:)                ! Element graph including halos
     integer(ip),      pointer :: c_dom_2(:)                ! Node CSR graph including halos
     integer(ip),      pointer :: r_dom_2(:)                ! Node CSR graph including halos
     integer(ip),      pointer :: r_dom_own(:)              ! Full row graph
     integer(ip),      pointer :: r_dom_end(:)              ! Full row graph
     integer(ip),      pointer :: r_dom_ini(:)              ! Full row graph
     integer(ip),      pointer :: c_dom_own(:)              ! Full row graph
     !
     ! Geometrical arrays
     !
     real(rp),         pointer :: exnor(:,:,:)              ! NDIME,NDIME,NBOPO: Local basis
     real(rp),         pointer :: vmass(:)                  ! NPOIN: Lumped mass matrix
     real(rp),         pointer :: vmasc(:)                  ! NPOIN: Mass matrix with closed rule
  end type mesh_type

  type(mesh_type), pointer :: meshe(:)
  integer(ip),     pointer :: lpmsh(:)
  integer(ip),     pointer :: lemsh(:)
  integer(ip),     pointer :: lbmsh(:)
  !
  ! Cut elements structure
  !
  type subel_type
     integer(ip)               :: inout      ! Type of elemen: -1 if is inside, 1 1 if is outside.
     real(rp),pointer          :: elcod(:,:) ! Coordinates of a subelement
  end type subel_type
  type subbo_type
     real(rp),pointer          :: bocod(:,:) ! Coordinates of the boundary
  end type subbo_type

  type cutel_type
     integer(ip)               ::  nelem     ! List of cut elements
     integer(ip)               ::  iimbo     ! Id of the cut particle
     type(subel_type), pointer ::  l(:)      ! List of subelements inside each cut element
     integer(ip)               ::  nboun     ! List of boundaries
     type(subbo_type), pointer ::  lb(:)     ! List of boundary elements that cut an element
     integer(ip),      pointer ::  linou(:)  ! Determina if a gauss point is inside or outside the particle
  end type cutel_type

  type(cutel_type),    pointer :: cutel(:)
  !
  ! Element bin
  !
  integer(ip)                    :: element_bin_boxes(3)
  type(typ_element_bin), pointer :: element_bin(:)
  !
  ! Additional vectors recalculated every time step for transient fields
  !
  real(rp)                  :: x_tran_fiel(mfiel)      ! interpolation for transient field
  integer(ip)               :: k_tran_fiel(mfiel)      ! begining of interval for transient field. Always = 1 for fields loaded on demand
  integer(ip)               :: k_tran_fiel_real(mfiel) ! real value of the k_tran_fiel. For checking if the timestep changed and files need to be loaded
  integer(ip), parameter    :: nsteps_fiel_ondemand=2  ! how many timesteps allocate for the fields loaded on demand (kfl_field(6,ifiel) == 1)
  integer(ip)               :: kexist_tran_fiel        ! do transient fields exist?
  !
  ! Edge codes
  !
  integer(ip),      pointer :: &
       kfl_coded(:,:)                                  ! Edge codes
  !
  ! Finite volume arrays
  !
  real(rp),         pointer ::   &
       fv_center_coord(:,:) ,    &                     ! Coordinates of centroides
       fv_cell_volume(:)    ,    &                     ! Volume of cells
       fv_face_area(:)      ,    &                     ! Face area
       fv_face_normal(:,:)  ,    &                     ! Face normals
       fv_face_orientation(:),   &                     ! Face orientation
       fv_center_distance(:),    &                     ! Distance bewteen centroids
       fv_center_vector(:,:),    &                     ! Vector between centroids
       fv_center_face(:,:,:)                           ! Vector from centroide to face centroide
  integer(ip),      pointer ::   &
       fv_face_boundary(:),      &                     ! Correspondance face booundary
       fv_face_graph(:,:),       &                     ! Correspondance face to graph position
       fv_graph_diag(:)                                ! Diagonal position of element in graph
  !
  ! Hash table to performa PAR_GLOBAL_TO_LOCAL_NODE
  !
  type(hash_t)              :: htable_lninv_loc        ! Hash table for global->local

end module def_domain
