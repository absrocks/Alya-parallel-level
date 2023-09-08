!------------------------------------------------------------------------
!> @defgroup Parall_Toolbox
!> Toolbox for parallelization tools
!> @{
!> @name    Parallelization toolbox
!> @file    mod_parall.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!> @brief   "The easiest way of making software scalable is to make it sequentially inefficient."
!>          - Bill Gropp, 1999
!>
!>          Michael E. Henderson, Christopher Radcliff Anderson, Stephen L. Lyons
!>          Object Oriented Methods for Interoperable Scientific and Engineering Computing: Proceedings of the 1998 SIAM Workshop
!>          SIAM, Jan 1, 1999 - Technology & Engineering - 321 pages
!>
!> @details ToolBox and definitions for parall.
!>
!>          The variables
!>          -------------
!>
!>          \verbatim
!>
!>          PAR_COMM_UNIVERSE........................................ Communicator of the universe  (Always MPI_COMM_WORLD)
!>          PAR_COMM_MY_CODE ........................................ Communicator of current code
!>          PAR_COMM_MY_CODE_WM4 .................................... Communicator of current code without master
!>          PAR_COMM_WORLD .......................................... Alya world communicator (MPI_COMM_WORLD or application split of MPI_COMM_WORLD)
!>          PAR_WORLD_SIZE .......................................... World size
!>          PAR_MY_UNIVERSE_RANK..................................... My rank in the universe
!>          PAR_MY_WORLD_RANK_WM .................................... My rank in the world without master
!>          PAR_MY_WORLD_RANK ....................................... My rank in the world
!>          PAR_UNIVERSE_SIZE........................................ Size of the universe
!>          PAR_CODE_SIZE ........................................... Size of my code communicator
!>          PAR_MY_CODE_RANK ........................................ My rank in my code
!>          PAR_MY_PARMETIS_RANK .................................... My rank in the parallel partitioners communicator
!>          PAR_MY_PARMETI2_RANK .................................... My rank in the parallel partitioners communicator + master
!>          PAR_INTEGER ............................................. Length of integer for MPI (ip)
!>
!>          I_AM_IN_COLOR(ICOLO) .................................... if I have color ICOLO (TRUE/FALSE)
!>          PAR_COMM_COLOR(0:MCOLO,0:MCOLO) ......................... Intercolor communicator
!>          PAR_COMM_COLOR_ARRAY(0:MCOLO) ........................... Color communication arrays
!>          PAR_CPU_TO_COLOR(0:PAR_WORLD_SIZE) % L(:) ............... List of colors for each world partition
!>          PAR_COLOR_TO_CPU(0:MCOLO) % L(:) ........................ List of world partitions for each color
!>          PAR_COMM_COLOR_PERM(0:MCOLO,0:MCOLO,0:PAR_WORLD_SIZE) ... Ranks for each communicator
!>          PAR_COMM_WORLD_TO_CODE_PERM(2,0:PAR_WORLD_SIZE) ......... Rank permutation from world to code
!>
!>          MCOLO ................................................... Nb of colors (codes,zones,subds) in the world
!>          MCODE ................................................... Nb of codes in the world
!>          MZONE ................................................... Nb of zones in the world
!>          MSUBD ................................................... Nb of subds in the world
!>          ncolo ................................................... Nb of colors in my code
!>          MAPPS ................................................... Nb of applications in the universe (MPI_COMM_WORLD)
!>
!>          \endverbatim
!>
!>          The functions
!>          -------------
!>
!>          \verbatim
!>
!>          PAR_WORLD_RANK_OF_A_CODE_NEIGHBOR ....................... World rank of a code neighbor
!>          PAR_WORLD_MASTER_RANK_OF_A_CODE ......................... Rank of the master of a code
!>          PAR_COLOR_COUPLING_RANK_TO_WORLD ........................ My rank in the world given my rank in a coupling
!>          PAR_CODE_ZONE_SUBD_TO_COLOR ............................. Mapping (code,zone,subd) => color
!>          PAR_COLOR_TO_CODE ....................................... Mapping color => code
!>          PAR_COLOR_TO_ZONE ....................................... Mapping color => zone
!>          PAR_COLOR_TO_SUBD ....................................... Mapping color => subd
!>          PAR_COLOR_TO_CODE_ZONE_SUBD(3) .......................... Mapping color => (code,zone,subd)
!>          PAR_PART_IN_COLOR ....................................... Check if a partition (rank in world) is in color
!>          PAR_THIS_NODE_IS_MINE ................................... If a node is interior or own boundary
!>          PAR_INITIALIZE_COMMUNICATION_ARRAY ...................... Initialize/Nullify a communication array
!>          PAR_COPY_COMMUNICATION_ARRAY ............................ Copy a communication array to another
!>          PAR_NODE_NUMBER_OF_PARTITIONS ........................... Number of partitions sharing this node
!>          PAR_GLOBAL_TO_LOCAL_NODE ................................ Local to global numbering of a ndoe
!>          par_omp_coloring ........................................ Color the elements for openmp
!>
!>          \endverbatim
!>
!>          Communication structure
!>          -----------------------
!>
!>          COMMD % NNEIG .................. Number of neigbors
!>          COMMD % NEIGHTS(:) ............. List of neighbor's ranks
!>          COMMD % BOUND_DIM ............. Total size of interface (repeated nodes)
!>          COMMD % BOUND_SIZE(:) ......... Linked lisT to exchange nodes with my neighbors
!>                                          Number of nodes to exchange with ISUBD: COMMD % BOUND_SIZE(ISUBD+1)-COMMD % BOUND_SIZE(ISUBD)
!>          COMMD % BOUND_PERM(:) ......... Linked list to exchange nodes with my neighbors
!>                                          commd % BOUND_PERM(COMMD % BOUND_SIZE(ISUBD):COMMD % BOUND_SIZE(ISUBD+1)-1) is the
!>                                          list of nodes to exchange with my neighbor ISUBD
!>          COMMD % BOUND_OWNER_RANK(:) ... Rank of the subdomain that owns the node
!>
!>
!------------------------------------------------------------------------

module mod_parall

  use def_kintyp,         only : ip,rp,lg,i1p,comm_data_par,ompss_domain
  use def_domain,         only : npoin,nelem,nboun,npoin_2
  use def_domain,         only : mesh_type
  use def_domain,         only : htable_lninv_loc
  use mod_graphs,         only : graphs_number_to_linked_list
  use def_master,         only : ISEQUEN,ISLAVE,IMASTER,INOTMASTER
  use def_master,         only : npoi1,npoi2,npoi3,lninv_loc
  use def_master,         only : intost,gisca,igene,lun_outpu
  use def_master,         only : ioutp
  use mod_memory,         only : memory_copy
  use mod_memory,         only : memory_alloca
  use mod_memory,         only : memory_deallo
  use mod_memory,         only : memory_size
  use mod_graphs,         only : graphs_eleele
  use mod_graphs,         only : graphs_dealep
  use mod_graphs,         only : graphs_coloring
  use mod_graphs,         only : graphs_coloring_greedy
  use mod_graphs,         only : graphs_deallocate
  use mod_htable,         only : htalid
  
  implicit none

  private

#ifndef MPI_OFF
  include  'mpif.h'
  integer(4) :: status(MPI_STATUS_SIZE)
#endif

  integer(ip)                   :: PAR_COMM_CURRENT                  ! Current communicator
  integer(ip)                   :: PAR_COMM_MY_CODE                  ! Communicator of current code (defined at beginning)
  integer(4)                    :: PAR_COMM_MY_CODE4                 ! Communicator of current code (defined at beginning)
  integer(4)                    :: PAR_COMM_MY_CODE_WM4              ! Communicator of current code without master
  integer(4)                    :: PAR_COMM_VIZ                      ! Communicator for insitu vizulization
  integer(4)                    :: PAR_COMM_SIM                      ! Communicator for simulation (same as alya world)
  integer(4)                    :: PAR_COMM_SIMVIZ                   ! Communicator combined simulation and viz
  integer(ip)                   :: PAR_COMM_UNIVERSE                 ! Universe communicator
  integer(ip)                   :: PAR_UNIVERSE_SIZE                 ! Size of the universe communicator
  integer(ip)                   :: PAR_MY_UNIVERSE_RANK              ! My rank in the universe
  integer(ip)                   :: PAR_COMM_WORLD                    ! Alya World communicator
  integer(ip)                   :: PAR_WORLD_SIZE                    ! Size of the world communicator
  integer(ip)                   :: PAR_MY_WORLD_RANK                 ! My communicator in the world
  integer(ip)                   :: PAR_MY_WORLD_RANK_WM              ! My communicator in the world without master
  integer(ip)                   :: PAR_CODE_SIZE                     ! Size of the code communicator
  integer(ip)                   :: PAR_MY_CODE_RANK                  ! My communicator in the world
  integer(ip)                   :: PAR_MY_PARMETIS_RANK              ! My communicator in the partitioners world
  integer(ip)                   :: PAR_MY_PARMETI2_RANK              ! My communicator in the partitioners world + master
  integer(4)                    :: PAR_COMM_SFC                      ! SFC partition communicator
  integer(4)                    :: PAR_COMM_SFC_WM                   ! SFC partition communicator without master
  integer(4)                    :: PAR_MY_SFC_RANK_WM                ! SFC partition communicator without master
  integer(ip)                   :: PAR_COMM_MPIO                     ! // READING communicator
  integer(ip)                   :: PAR_COMM_MPIO_WM                  ! // READING communicator without master
  integer(ip)                   :: PAR_COMM_MPIO_WM_SIZE             ! // READING communicator size without master
  integer(ip)                   :: PAR_COMM_MPIO_RANK_WM             ! // READING rank without master
  integer(4)                    :: PAR_INTEGER                       ! Integer types for METIS
  integer(ip),         pointer  :: PAR_COMM_COLOR(:,:)               ! Intercolor MPI communicator
  integer(ip),         pointer  :: PAR_COMM_COLOR_PERM(:,:,:)        ! Communicator permutation array
  integer(ip),         pointer  :: PAR_COMM_WORLD_TO_CODE_PERM(:,:)  ! Permutation world to code and rank
  logical(lg),         pointer  :: I_AM_IN_COLOR(:)                  ! If I am in color

  type(comm_data_par), pointer  :: PAR_COMM_COLOR_ARRAY(:)           ! Intercolor communicator
  type(comm_data_par), pointer  :: PAR_COMM_MY_CODE_ARRAY(:)         ! Intercolor communicator
  type(comm_data_par), pointer  :: PAR_COMM_ZONE_ARRAY(:)            ! Interzone communicator
  type(comm_data_par), pointer  :: commd                             ! Generic communication array
  !
  ! Color structure
  !
  type(i1p),           pointer  :: PAR_CPU_TO_COLOR(:)               ! Given a CPU, gives the colors and associated rank
  type(i1p),           pointer  :: PAR_COLOR_TO_CPU(:)               ! Given a color, gives the CPUs
  type(i1p),           pointer  :: PAR_COLOR_BIN(:)                  ! Bin structure for colors
  !
  ! Current CODE/ZONE/SUBDOMAIN
  !
  integer(ip)                   :: mcolo                             ! Maximum number of colors
  integer(ip)                   :: mcode                             ! Maximum number of codes (over all processes)
  integer(ip)                   :: mzone                             ! Maximum number of zones (over all processes)
  integer(ip)                   :: msubd                             ! Maximum number of subds (over all processes)
  integer(ip)                   :: ncolo                             ! Number of colors of my code
  integer(ip)                   :: mapps                             ! Number of applications in the universe (MPI_COMM_WORLD)
  !
  ! Bin structure to perform subdomain geometrical queries
  ! and subdomain bounding box
  !
  integer(ip)                   :: par_bin_boxes(3)                  ! # boxes in each direction
  real(rp)                      :: par_bin_comin(3)                  ! Minimum box coordinates
  real(rp)                      :: par_bin_comax(3)                  ! Maximum box coordinates
  integer(ip),         pointer  :: par_bin_size(:,:,:)               ! # partitions per box
  type(i1p),           pointer  :: par_bin_part(:,:,:)               ! Bin structure of world partition
  real(rp),            pointer  :: par_part_comin(:,:)               ! Subdomain minimum coordinates
  real(rp),            pointer  :: par_part_comax(:,:)               ! Subdomain maximum coordinates
  type typ_bin_structure
     integer(ip)                :: boxes(3)                          ! # boxes in each direction
     real(rp)                   :: comin(3)                          ! Minimum box coordinates
     real(rp)                   :: comax(3)                          ! Maximum box coordinates
     integer(ip),      pointer  :: size(:,:,:)                       ! # partitions per box
     type(i1p),        pointer  :: part(:,:,:)                       ! Bin structure of world partition
     real(rp),         pointer  :: part_comin(:,:)                   ! Subdomain minimum coordinates
     real(rp),         pointer  :: part_comax(:,:)                   ! Subdomain maximum coordinates
  end type typ_bin_structure
  !
  ! Others
  !
  integer(8)                    :: par_memor(2)                      ! Memory counter
  integer(ip)                   :: color_target                      ! Target color when using coupling
  integer(ip)                   :: color_source                      ! Source color when using coupling
  integer(ip),        parameter :: lun_outpu_par = 5502              ! Output
  integer(ip),        parameter :: lun_parti_msh = 5507              ! Partition mesh file
  integer(ip),        parameter :: lun_parti_res = 5508              ! Partition result file
  integer(ip),        parameter :: lun_matri_msh = 5509              ! Matrix mesh file
  integer(ip),        parameter :: lun_matri_res = 5510              ! Matrix result file
  integer(4)                    :: PARMETIS_COMM, PARMETIS_INTERCOMM ! parmetis communicators
  integer(4)                    :: PARMETI2_COMM                     ! parmetis communicators + master
  !
  ! OpenMP
  !
  integer(ip)                   :: par_omp_num_blocks                ! Size of blocks to compute chunks n/par_omp_num_blocks
  integer(ip)                   :: par_omp_granularity               ! Granularity to compute par_omp_num_blocks=par_omp_num_threads*par_omp_granularity
  integer(ip)                   :: par_omp_nelem_chunk   = 0         ! Element loop chunk size
  integer(ip)                   :: par_omp_npoin_chunk               ! Node loop chunk size
  integer(ip)                   :: par_omp_nboun_chunk               ! Boundary loop chunk size
  integer(ip)                   :: par_omp_num_threads               ! Number of openmp threads
  integer(ip)                   :: par_omp_coloring_alg              ! Coloring algorithm
  integer(ip)                   :: par_omp_partition_alg             ! Partitioning algorithm for OmpSs (use same nomenclature for mesh partitioning)

  integer(ip)                   :: par_omp_num_colors                ! Element: Number of colors
  integer(ip),        pointer   :: par_omp_list_colors(:)            ! Element: Element colors
  integer(ip),        pointer   :: par_omp_ia_colors(:)              ! Element: Linked list IA for colors
  integer(ip),        pointer   :: par_omp_ja_colors(:)              ! Element: Linked list JA for colors

  integer(ip)                   :: par_omp_nboun_num_colors          ! Boundary: Number of colors
  integer(ip),        pointer   :: par_omp_nboun_list_colors(:)      ! Boundary: Element colors
  integer(ip),        pointer   :: par_omp_nboun_ia_colors(:)        ! Boundary: Linked list IA for colors
  integer(ip),        pointer   :: par_omp_nboun_ja_colors(:)        ! Boundary: Linked list JA for colors
  !
  ! Hybrid parallelization (OpenMP/OmpSs) and vectorization
  !
  integer(ip),        parameter        :: PAR_HYBRID_OFF         = 0
  integer(ip),        parameter        :: PAR_OPENMP_COLORING    = 1
  integer(ip),        parameter        :: PAR_OPENMP_NO_COLORING = 2
  integer(ip),        parameter        :: PAR_OMPSS              = 3
  integer(ip)                          :: par_hybrid
  type typ_list_elements_par
     type(i1p),                pointer :: packs(:)
  end type typ_list_elements_par

  integer(ip)                          :: num_subd_par               ! Element loop: For loops with race condition
  integer(ip),                 pointer :: num_pack_par(:)
  type(typ_list_elements_par), pointer :: list_elements_par(:)

  integer(ip)                          :: num_subd_norace_par        ! Element loop: For loops without race condition
  integer(ip),                 pointer :: num_pack_norace_par(:)
  type(typ_list_elements_par), pointer :: list_elements_norace_par(:)

  integer(ip)                          :: num_subd_nboun_par         ! Boundary loop: For loops with race condition
  integer(ip),                 pointer :: num_pack_nboun_par(:)
  type(typ_list_elements_par), pointer :: list_boundaries_par(:)

  integer(ip)                          :: num_subd_norace_nboun_par  ! Boundary loop: For loops without race condition
  integer(ip),                 pointer :: num_pack_norace_nboun_par(:)
  type(typ_list_elements_par), pointer :: list_boundaries_norace_par(:)
  !
  ! Topology
  !
  integer(ip)                          :: par_topo_num_nodes                ! Number of computing nodes
  integer(ip)                          :: par_topo_num_cores_per_node       ! Number of cores per node
  !
  ! OPMSs
  !
!  type ompss_domain
!     integer(ip), allocatable          :: neighbours(:)
!     integer(ip), allocatable          :: elements(:)
!     integer(ip)                       :: neighIdx
!     integer(ip)                       :: elemIdx
!  end type ompss_domain
  !
  ! Partitioning method
  !
  integer(ip),        parameter        :: PAR_METIS4               = 0
  integer(ip),        parameter        :: PAR_SFC                  = 1
  integer(ip),        parameter        :: PAR_ORIENTED_BIN         = 2
  integer(ip),        parameter        :: PAR_FRONTAL              = 3 
  integer(ip),        parameter        :: PAR_NUMBERING            = 4
  integer(ip),        parameter        :: PAR_USING_RANK           = 5
  
  integer(ip),        parameter        :: PAR_SEQUENTIAL_PARTITION = 0
  integer(ip),        parameter        :: PAR_PARALLEL_PARTITION   = 1

  integer(ip),        parameter        :: PAR_WEIGHT_GAUSS         =  0
  integer(ip),        parameter        :: PAR_WEIGHT_OFF           = -1
  integer(ip),        parameter        :: PAR_WEIGHT_ELEMENT       = -2
  integer(ip),        parameter        :: PAR_WEIGHT_SQUARE        = -3
  integer(ip),        parameter        :: PAR_WEIGHT_MATERIAL      = -4
  !
  ! Files
  !
  character(150)                       :: fil_parti_msh                     ! Partition mesh
  character(150)                       :: fil_parti_res                     ! Partition result
  !
  ! Interfaces
  !
  interface PAR_INITIALIZE_COMMUNICATION_ARRAY
     module procedure PAR_INITIALIZE_COMMUNICATION_ARRAY_s, &
          &           PAR_INITIALIZE_COMMUNICATION_ARRAY_1
  end interface PAR_INITIALIZE_COMMUNICATION_ARRAY

  interface NODE_IN_NEIGHBOR
     module procedure NODE_IN_NEIGHBOR_P,NODE_IN_NEIGHBOR_s
  end interface NODE_IN_NEIGHBOR

  public :: PAR_INITIALIZE_COMMUNICATION_ARRAY
  public :: PAR_DEALLOCATE_COMMUNICATION_ARRAY
  public :: PAR_TRANSPOSE_COMMUNICATION_ARRAY
  public :: NODE_IN_NEIGHBOR
  !
  ! Subroutines
  !
  public :: par_code_zone_subd_to_color
  public :: par_color_to_code
  public :: par_color_to_zone
  public :: par_color_to_subd
  public :: par_part_in_color
  public :: par_color_coupling_rank_to_world
  public :: par_world_rank_of_a_code_neighbor
  public :: PAR_THIS_NODE_IS_MINE
  public :: PAR_GLOBAL_TO_LOCAL_NODE
  public :: PAR_SUB_COMMUNICATION_ARRAY
  public :: PAR_NODE_NUMBER_OF_PARTITIONS
  public :: PAR_COPY_COMMUNICATION_ARRAY
  public :: PAR_POINT_COMMUNICATION_ARRAY
  !
  ! Types
  !
  public :: typ_bin_structure
  !
  ! Variables
  !
  public :: PAR_COMM_CURRENT
  public :: PAR_COMM_MY_CODE
  public :: PAR_COMM_MY_CODE4
  public :: PAR_COMM_MY_CODE_WM4
  public :: PAR_COMM_VIZ
  public :: PAR_COMM_SIM
  public :: PAR_COMM_SIMVIZ
  public :: PAR_COMM_UNIVERSE
  public :: PAR_UNIVERSE_SIZE
  public :: PAR_MY_UNIVERSE_RANK
  public :: PAR_COMM_WORLD
  public :: PAR_WORLD_SIZE
  public :: PAR_MY_WORLD_RANK
  public :: PAR_MY_WORLD_RANK_WM
  public :: PAR_CODE_SIZE
  public :: PAR_MY_CODE_RANK
  public :: PAR_MY_PARMETIS_RANK
  public :: PAR_MY_PARMETI2_RANK
  public :: PAR_INTEGER
  public :: PAR_COMM_SFC
  public :: PAR_COMM_SFC_WM
  public :: PAR_MY_SFC_RANK_WM
  public :: PAR_COMM_COLOR
  public :: PAR_COMM_COLOR_PERM
  public :: PAR_COMM_WORLD_TO_CODE_PERM
  public :: I_AM_IN_COLOR
  public :: PAR_COMM_COLOR_ARRAY
  public :: PAR_COMM_MY_CODE_ARRAY
  public :: PAR_COMM_ZONE_ARRAY
  public :: commd
  public :: PAR_CPU_TO_COLOR
  public :: PAR_COLOR_TO_CPU
  public :: PAR_COLOR_BIN
  public :: mcolo
  public :: mcode
  public :: mzone
  public :: msubd
  public :: ncolo
  public :: mapps
  public :: par_bin_boxes
  public :: par_bin_comin
  public :: par_bin_comax
  public :: par_bin_size
  public :: par_bin_part
  public :: par_part_comin
  public :: par_part_comax
  public :: par_memor
  public :: color_target
  public :: color_source
  public :: lun_outpu_par
  public :: lun_parti_msh
  public :: lun_matri_res
  public :: lun_matri_msh
  public :: lun_parti_res
  public :: PARMETIS_COMM
  public :: PARMETI2_COMM
  public :: PARMETIS_INTERCOMM
  public :: par_omp_num_blocks
  public :: par_omp_granularity
  public :: par_omp_nelem_chunk
  public :: par_omp_npoin_chunk
  public :: par_omp_nboun_chunk
  public :: par_omp_num_threads

  public :: par_omp_num_colors
  public :: par_omp_list_colors
  public :: par_omp_ia_colors
  public :: par_omp_ja_colors

  public :: par_omp_nboun_num_colors
  public :: par_omp_nboun_list_colors
  public :: par_omp_nboun_ia_colors
  public :: par_omp_nboun_ja_colors

  public :: par_omp_coloring_alg
  public :: par_omp_partition_alg
  public :: ompss_domain
  public :: par_topo_num_nodes                ! Number of computing nodes
  public :: par_topo_num_cores_per_node       ! Number of cores per node

  public :: PAR_METIS4
  public :: PAR_SFC
  public :: PAR_ORIENTED_BIN
  public :: PAR_NUMBERING
  public :: PAR_USING_RANK
  public :: PAR_SEQUENTIAL_PARTITION
  public :: PAR_PARALLEL_PARTITION
  public :: PAR_WEIGHT_GAUSS     
  public :: PAR_WEIGHT_OFF       
  public :: PAR_WEIGHT_ELEMENT   
  public :: PAR_WEIGHT_SQUARE    
  public :: PAR_WEIGHT_MATERIAL

  public :: PAR_COMM_MPIO                       ! // READING communicator
  public :: PAR_COMM_MPIO_WM                    ! // READING communicator without master
  public :: PAR_COMM_MPIO_RANK_WM               ! // READING rank without master
  public :: PAR_COMM_MPIO_WM_SIZE

  public :: fil_parti_msh                     ! Partition mesh
  public :: fil_parti_res                     ! Partition result
  !
  ! Hybrid parallelization
  !
  public :: PAR_OPENMP_COLORING
  public :: PAR_OPENMP_NO_COLORING
  public :: PAR_OMPSS
  public :: PAR_HYBRID_OFF
  public :: par_hybrid

  public :: typ_list_elements_par
  public :: num_subd_par
  public :: num_pack_par
  public :: list_elements_par
  public :: num_subd_norace_par
  public :: num_pack_norace_par
  public :: list_elements_norace_par

  public :: num_subd_nboun_par
  public :: num_pack_nboun_par
  public :: list_boundaries_par
  public :: num_subd_norace_nboun_par
  public :: num_pack_norace_nboun_par
  public :: list_boundaries_norace_par

contains

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/03/2014
  !> @brief   Give the world rank of a code neighbor
  !> @details Given a neighbor INEIG, returns its rank in the MPI_COMM_WORLD
  !
  !----------------------------------------------------------------------

  function par_world_rank_of_a_code_neighbor(ipart,icode)
    integer(ip), intent(in) :: ipart
    integer(ip), intent(in) :: icode
    integer(ip)             :: par_world_rank_of_a_code_neighbor
    integer(ip)             :: ipart_world,ipart_code,jcode

    ipart_world = 0
    par_world_rank_of_a_code_neighbor = -1
    do while( ipart_world <= PAR_WORLD_SIZE-1 )
       jcode      = PAR_COMM_WORLD_TO_CODE_PERM(1,ipart_world)
       ipart_code = PAR_COMM_WORLD_TO_CODE_PERM(2,ipart_world)
       if( ipart_code == ipart .and. icode == jcode ) then
          par_world_rank_of_a_code_neighbor = ipart_world
          ipart_world = PAR_WORLD_SIZE-1
       end if
       ipart_world = ipart_world + 1
    end do
    if( par_world_rank_of_a_code_neighbor == -1 ) &
         call runend('par_world_rank_of_a_code_neighbor: WE ARE IN TROUBLE')

  end function par_world_rank_of_a_code_neighbor

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/03/2014
  !> @brief   Give the rank of the master of a code
  !> @details Given a code ICODE, returns the rank of the master of the
  !>          code in the MPI_COMM_WORLD
  !
  !----------------------------------------------------------------------

  function par_world_master_rank_of_a_code(icode)
    integer(ip), intent(in) :: icode
    integer(ip)             :: par_world_master_rank_of_a_code
    integer(ip)             :: icolo,isize

    icolo = par_code_zone_subd_to_color(icode,0_ip,0_ip)
    isize = 0
    do while( isize <= PAR_WORLD_SIZE-1 )
       par_world_master_rank_of_a_code = PAR_COMM_COLOR_PERM(icolo,icolo,isize)
       if( par_world_master_rank_of_a_code == 0 ) then
          isize = PAR_WORLD_SIZE
       end if
       isize = isize + 1
    end do

  end function par_world_master_rank_of_a_code

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    13/03/2014
  !> @brief   Give my rank in the world given my rank in a coupling
  !> @details Give my rank in the world given my rank in a coupling
  !
  !----------------------------------------------------------------------

  function par_color_coupling_rank_to_world(my_rank,icolo,jcolo)
    integer(ip), intent(in) :: my_rank
    integer(ip), intent(in) :: icolo
    integer(ip), intent(in) :: jcolo
    integer(ip)             :: par_color_coupling_rank_to_world
    integer(ip)             :: isize

    isize = 0
    do while( isize <= PAR_WORLD_SIZE-1 )
       par_color_coupling_rank_to_world = PAR_COMM_COLOR_PERM(icolo,jcolo,isize)
       if( par_color_coupling_rank_to_world == my_rank ) then
          isize = PAR_WORLD_SIZE
       end if
       isize = isize + 1
    end do

  end function par_color_coupling_rank_to_world

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Mapping (code,zone,subd) => color
  !> @details Given a code, zone and subd, compute the corresponding
  !>          color, where:
  !>          code: 0 to mcode
  !>          zone: 0 to mzone
  !>          subd: 0 to msubd
  !>          3D <=> 1D mapping:
  !>          i = x + Lx * (y+Ly*z)
  !>          x = modulo(i,Lx)
  !>          y = modulo(i/Lx,Ly)
  !>          z = i/(Lx*Ly)
  !
  !----------------------------------------------------------------------

  function par_code_zone_subd_to_color(icode,izone,isubd)
    integer(ip), intent(in) :: icode,izone,isubd
    integer(ip)             :: par_code_zone_subd_to_color

    par_code_zone_subd_to_color = isubd+(msubd+1)*(izone+icode*(mzone+1))

  end function par_code_zone_subd_to_color

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Mapping color => code
  !> @details Given a color, returns the code number
  !
  !----------------------------------------------------------------------

  function par_color_to_code(icolo)
    integer(ip), intent(in) :: icolo
    integer(ip)             :: par_color_to_code

    par_color_to_code = icolo / ((msubd+1)*(mzone+1))

  end function par_color_to_code

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Mapping color => code
  !> @details Given a color, returns the zone number
  !
  !----------------------------------------------------------------------

  function par_color_to_zone(icolo)
    integer(ip), intent(in) :: icolo
    integer(ip)             :: par_color_to_zone

    par_color_to_zone = modulo(icolo/(msubd+1),mzone+1)

  end function par_color_to_zone

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Mapping color => code
  !> @details Given a color, returns the subdomain number
  !
  !----------------------------------------------------------------------

  function par_color_to_subd(icolo)
    integer(ip), intent(in) :: icolo
    integer(ip)             :: par_color_to_subd

    par_color_to_subd = modulo(icolo,msubd+1)

  end function par_color_to_subd

  function par_color_to_code_zone_subd(icolo)
    implicit none
    integer(ip), intent(in) :: icolo
    integer(ip)             :: par_color_to_code_zone_subd(3)

    par_color_to_code_zone_subd(1) = icolo / ((msubd+1)*(mzone+1))
    par_color_to_code_zone_subd(2) = modulo(icolo/(msubd+1),mzone+1)
    par_color_to_code_zone_subd(3) = modulo(icolo,msubd+1)

  end function par_color_to_code_zone_subd

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    28/06/2012
  !> @brief   Check if a partition IPART is in color ICOLO
  !> @details IPART should be given in the world
  !>                                                       ICOLO
  !>                                           +---+---+---+---+---+---+
  !>          PAR_CPU_TO_COLOR(IPART) % L(:) = | 1 | 4 | 7 | 8 | 9 |10 |
  !>                                           +---+---+---+---+---+---+
  !>                                           KCOLO =>
  !>
  !----------------------------------------------------------------------

  function par_part_in_color(ipart,icolo)
    integer(ip), intent(in) :: ipart
    integer(ip), intent(in) :: icolo
    logical(lg)             :: par_part_in_color
    integer(ip)             :: kcolo,jcolo,ksize

    par_part_in_color = .false.
    if( associated(PAR_CPU_TO_COLOR(ipart) % l) ) then
       jcolo = 1
       kcolo = PAR_CPU_TO_COLOR(ipart) % l(jcolo)
       ksize = size(PAR_CPU_TO_COLOR(ipart) % l)
       do while( jcolo <= ksize .and. kcolo < icolo )
          kcolo = PAR_CPU_TO_COLOR(ipart) % l(jcolo)
          jcolo = jcolo + 1
       end do
       if( kcolo == icolo ) par_part_in_color = .true.
    end if

  end function par_part_in_color

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    06/03/2014
  !> @brief   Returns if a node is interior or inside my own boundary
  !> @details Returns if a node is interior or inside my own boundary
  !>
  !----------------------------------------------------------------------

  function PAR_THIS_NODE_IS_MINE(ipoin,where)
    integer(ip),  intent(in)           :: ipoin
    character(*), intent(in), optional :: where
    logical(lg)                        :: PAR_THIS_NODE_IS_MINE

    PAR_THIS_NODE_IS_MINE = .false.

    if( ISEQUEN ) then
       if( ipoin <= npoin ) PAR_THIS_NODE_IS_MINE = .true.
    else if( INOTMASTER ) then
       if( .not. present(where) ) then
          if( ipoin <= npoi3 ) then
             PAR_THIS_NODE_IS_MINE = .true.
          else
             PAR_THIS_NODE_IS_MINE = .false.
          end if
       else
          call runend('PAR_THIS_NODE_IS_MINE: NOT CODED')
       end if
    end if

  end function PAR_THIS_NODE_IS_MINE

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    06/03/2014
  !> @brief   Returns the local node number of a global number
  !> @details Returns the local node number of a global number
  !>          and zero if I do not own the node
  !
  !----------------------------------------------------------------------

  function PAR_GLOBAL_TO_LOCAL_NODE(ipoin_global,wherein,mask)
    
    integer(ip),  intent(in)                    :: ipoin_global
    character(*), intent(in), optional          :: wherein
    logical(lg),  intent(in), optional, pointer :: mask(:)
    integer(ip)                                 :: ipoin
    integer(ip)                                 :: PAR_GLOBAL_TO_LOCAL_NODE
    integer(ip)                                 :: npoin_end

    npoin_end = npoin 
    if( present(wherein) ) then
       if( trim(wherein) == 'INCLUDING HALOS' ) then
          npoin_end = npoin_2
       end if
    end if

    PAR_GLOBAL_TO_LOCAL_NODE = 0

    if( ISEQUEN ) then
       PAR_GLOBAL_TO_LOCAL_NODE = ipoin_global
    else if( ISLAVE ) then
       if( htable_lninv_loc % sizet /= 0 ) then
          ipoin = htalid(htable_lninv_loc,ipoin_global)
          if( ipoin < 1 .or. ipoin > npoin_end ) ipoin = 0
          PAR_GLOBAL_TO_LOCAL_NODE = ipoin
          if( present(mask) ) then
             if( .not. mask(ipoin) ) PAR_GLOBAL_TO_LOCAL_NODE = 0
          end if
          return
       else
          if( present(mask) ) then
             do ipoin = 1,npoin_end
                if( mask(ipoin) ) then
                   if( lninv_loc(ipoin) == ipoin_global ) then
                      PAR_GLOBAL_TO_LOCAL_NODE = ipoin
                      return
                   end if
                end if
             end do
          else
             do ipoin = 1,npoin_end
                if( lninv_loc(ipoin) == ipoin_global ) then
                   PAR_GLOBAL_TO_LOCAL_NODE = ipoin
                   return
                end if
             end do
          end if
       end if
    end if

  end function PAR_GLOBAL_TO_LOCAL_NODE

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    06/03/2014
  !> @brief   Number of partitions for nodes
  !> @details Gives the number of neighbors for a given node or all
  !>          the nodes. This includes my partition. In next example,
  !>          the middle node (o) has 4 partitions.
  !>
  !>          x----------x----------x
  !>          |          |          |
  !>          |   CPU1   |   CPU2   |
  !>          |          |          |
  !>          x----------o----------x
  !>          |          |          |
  !>          |   CPU3   |   CPU4   |
  !>          |          |          |
  !>          x----------x----------x
  !
  !----------------------------------------------------------------------

  subroutine PAR_NODE_NUMBER_OF_PARTITIONS(npoin,lneig,ipoin)
    integer(ip), intent(in)           :: npoin
    integer(ip), intent(out)          :: lneig(*)
    integer(ip), intent(in), optional :: ipoin
    integer(ip)                       :: kpoin
    integer(ip)                       :: ineig,jj

    if( present(ipoin) ) then
       lneig(1) = 1
       do jj = 1,commd % bound_dim
          kpoin = commd % bound_perm(jj)
          if( kpoin == ipoin ) lneig(ipoin) = lneig(ipoin) + 1
       end do
    else
       do kpoin = 1,npoin
          lneig(kpoin) = 1
       end do
       do jj = 1,commd % bound_dim
          kpoin = commd % bound_perm(jj)
          lneig(kpoin) = lneig(kpoin) + 1
       end do
    end if

  end subroutine PAR_NODE_NUMBER_OF_PARTITIONS

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    27/03/2018
  !> @brief   Deallocate communication array
  !> @details Deallocate communication array
  !
  !----------------------------------------------------------------------

  subroutine PAR_DEALLOCATE_COMMUNICATION_ARRAY(COMM,memor_opt,PAR_COMM_OPT,OPTION)

    type(comm_data_par), intent(inout)           :: COMM          !< Input communicator
    integer(8),          intent(inout), optional :: memor_opt(2)  !< Memory counter
    logical(lg),         intent(in),    optional :: PAR_COMM_OPT
    character(*),        intent(in),    optional :: OPTION
    integer(8)                                   :: memor(2)
    logical(lg)                                  :: PAR_COMM
    integer(ip)                                  :: ii
    logical(lg)                                  :: only_nodes

    only_nodes = .false.
    PAR_COMM   = .true.

    if( present(PAR_COMM_OPT) ) PAR_COMM = PAR_COMM_OPT
    if( present(option) ) THEN
       if( trim(OPTION) == 'ONLY NODES' ) only_nodes = .true.
    end if
    if( present(memor_opt) ) then
       memor = memor_opt
    else
       memor = 0_8
    end if

    COMM % nneig               = 0
    COMM % npoi1               = 0
    COMM % npoi2               = 0
    COMM % npoi3               = 0
    COMM % bound_dim           = 0
    call memory_deallo(memor,'NEIGHTS'             ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % neights      )
    call memory_deallo(memor,'BOUND_SIZE'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bound_size   )
    call memory_deallo(memor,'BOUND_PERM'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bound_perm   )

    if( .not. only_nodes ) then
       COMM % nedg1               = 0
       COMM % nedg2               = 0
       COMM % nedg3               = 0

       COMM % bedge_dim           = 0
       COMM % lsend_dim           = 0
       COMM % lrecv_dim           = 0
       COMM % lscat_dim           = 0
       COMM % matrix_nzdom        = 0
       COMM % bface_dim           = 0
       COMM % full_row_send_dim   = 0
       COMM % full_row_recv_nneig = 0
       COMM % full_row_recv_dim   = 0
       COMM % ghost_send_node_dim = 0
       COMM % ghost_recv_node_dim = 0
       COMM % ghost_send_elem_dim = 0
       COMM % ghost_recv_elem_dim = 0
       COMM % ghost_send_boun_dim = 0
       COMM % ghost_recv_boun_dim = 0
       if(  PAR_COMM ) COMM % PAR_COMM_WORLD      = 0
       COMM % RANK4               = -1_4

       call memory_deallo(memor,'BOUND_SCAL'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bound_scal   )
       call memory_deallo(memor,'BOUND_SCAL'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bound_multiplicity )
       call memory_deallo(memor,'BOUND_SCAL'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bound_owner_rank )

       call memory_deallo(memor,'BFACE_SIZE'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bface_size   )
       call memory_deallo(memor,'BFACE_PERM'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bface_perm   )

       call memory_deallo(memor,'BFACE_PERM'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % node_number_in_owner)

       if( associated(COMM % bound_matrix) ) then
          do ii = 1,size(COMM % bound_matrix)
             call memory_deallo(memor,'BFACE_PERM'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bound_matrix(ii) % ja)
             call memory_deallo(memor,'BFACE_PERM'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bound_matrix(ii) % nzdom_ii)
             COMM % bound_matrix(ii) % nzdom_ii = 0
          end do
          nullify(COMM % bound_matrix)
       end if
       if( associated(COMM % bound_mat_halo_send) ) then
          do ii = 1,size(COMM % bound_mat_halo_send)
             call memory_deallo(memor,'BFACE_PERM'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bound_mat_halo_send(ii) % ja)
             call memory_deallo(memor,'BFACE_PERM'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bound_mat_halo_send(ii) % nzdom_ii)
             COMM % bound_mat_halo_send(ii) % nzdom_ii = 0
          end do
          nullify(COMM % bound_mat_halo_send)
       end if
       if( associated(COMM % bound_mat_halo_recv) ) then
          do ii = 1,size(COMM % bound_mat_halo_recv)
             call memory_deallo(memor,'BFACE_PERM'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bound_mat_halo_recv(ii) % ja)
             call memory_deallo(memor,'BFACE_PERM'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bound_mat_halo_recv(ii) % nzdom_ii)
             COMM % bound_mat_halo_recv(ii) % nzdom_ii = 0
          end do
          nullify(COMM % bound_mat_halo_recv)
       end if
       call memory_deallo(memor,'BFACE_PERM'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bedge_size)
       call memory_deallo(memor,'BFACE_PERM'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bedge_perm)
       call memory_deallo(memor,'BFACE_PERM'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bedge_adja)
       call memory_deallo(memor,'BFACE_PERM'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bedge_scal)
       call memory_deallo(memor,'BFACE_PERM'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bedge_multiplicity)
       call memory_deallo(memor,'BFACE_PERM'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % bedge_owner_rank)
       call memory_deallo(memor,'BFACE_PERM'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % lsend_size)
       call memory_deallo(memor,'BFACE_PERM'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % lsend_perm)
       call memory_deallo(memor,'BFACE_PERM'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % lrecv_size)
       call memory_deallo(memor,'BFACE_PERM'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % lrecv_perm)
       call memory_deallo(memor,'BFACE_PERM'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % lscat_perm)
       call memory_deallo(memor,'BFACE_PERM'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % matrix_ia)
       call memory_deallo(memor,'BFACE_PERM'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % matrix_ja)
       call memory_deallo(memor,'BFACE_PERM'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % matrix_aa)
       call memory_deallo(memor,'BFACE_PERM'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % full_row_send_neights)
       call memory_deallo(memor,'BFACE_PERM'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % full_row_send_size)
       call memory_deallo(memor,'BFACE_PERM'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % full_row_send_perm)
       call memory_deallo(memor,'BFACE_PERM'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % full_row_recv_neights)
       call memory_deallo(memor,'BFACE_PERM'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % full_row_recv_size)
       call memory_deallo(memor,'BFACE_PERM'          ,'PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % full_row_recv_perm)

       call memory_deallo(memor,'GHOST_SEND_NODE_SIZE','PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % ghost_send_node_size   )
       call memory_deallo(memor,'GHOST_SEND_NODE_PERM','PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % ghost_send_node_perm   )
       call memory_deallo(memor,'GHOST_RECV_NODE_SIZE','PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % ghost_recv_node_size   )
       call memory_deallo(memor,'GHOST_RECV_NODE_PERM','PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % ghost_recv_node_perm   )
       call memory_deallo(memor,'GHOST_SEND_ELEM_SIZE','PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % ghost_send_elem_size   )
       call memory_deallo(memor,'GHOST_SEND_ELEM_PERM','PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % ghost_send_elem_perm   )
       call memory_deallo(memor,'GHOST_RECV_ELEM_SIZE','PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % ghost_recv_elem_size   )
       call memory_deallo(memor,'GHOST_RECV_ELEM_PERM','PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % ghost_recv_elem_perm   )
       call memory_deallo(memor,'GHOST_SEND_BOUN_SIZE','PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % ghost_send_boun_size   )
       call memory_deallo(memor,'GHOST_SEND_BOUN_PERM','PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % ghost_send_boun_perm   )
       call memory_deallo(memor,'GHOST_RECV_BOUN_SIZE','PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % ghost_recv_boun_size   )
       call memory_deallo(memor,'GHOST_RECV_BOUN_PERM','PAR_DEALLOCATE_COMMUNICATION_ARRAY',COMM % ghost_recv_boun_perm   )

       if( present(memor_opt) ) then
          memor_opt = memor
       end if
    end if

  end subroutine PAR_DEALLOCATE_COMMUNICATION_ARRAY

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    06/03/2014
  !> @brief   Initialize communication array
  !> @details Initialize communication array
  !
  !----------------------------------------------------------------------

  subroutine PAR_INITIALIZE_COMMUNICATION_ARRAY_s(COMM,PAR_COMM_OPT)
    type(comm_data_par), intent(inout)        :: COMM     !< Input communicator
    logical(lg),         intent(in), optional :: PAR_COMM_OPT
    logical(lg)                               :: PAR_COMM

    PAR_COMM = .true.
    if( present(PAR_COMM_OPT) ) PAR_COMM = PAR_COMM_OPT

    COMM % nneig               =  0
    COMM % npoi1               =  0
    COMM % npoi2               =  0
    COMM % npoi3               = -1
    COMM % nedg1               =  0
    COMM % nedg2               =  0
    COMM % nedg3               = -1
    COMM % bound_dim           =  0
    COMM % bedge_dim           =  0
    COMM % lsend_dim           =  0
    COMM % lrecv_dim           =  0
    COMM % lscat_dim           =  0
    COMM % matrix_nzdom        =  0
    COMM % bface_dim           =  0

    COMM % full_row_send_nneig =  0
    COMM % full_row_send_dim   = -1
    COMM % full_row_recv_nneig =  0
    COMM % full_row_recv_dim   = -1

    COMM % ghost_send_elem_dim = -1
    COMM % ghost_recv_elem_dim = -1
    COMM % ghost_send_node_dim = -1
    COMM % ghost_recv_node_dim = -1
    COMM % ghost_send_boun_dim = -1
    COMM % ghost_recv_boun_dim = -1

    COMM % RANK4               = -1_4

    if( PAR_COMM ) COMM % PAR_COMM_WORLD = -1
    nullify(COMM % neights)

    nullify(COMM % bound_size)
    nullify(COMM % bound_perm)
    nullify(COMM % bound_scal)
    nullify(COMM % bound_multiplicity)
    nullify(COMM % bound_owner_rank)
    nullify(COMM % node_number_in_owner)
    nullify(COMM % bound_matrix)
    nullify(COMM % bound_mat_halo_send)
    nullify(COMM % bound_mat_halo_recv)

    nullify(COMM % bedge_size)
    nullify(COMM % bedge_perm)
    nullify(COMM % bedge_adja)
    nullify(COMM % bedge_scal)
    nullify(COMM % bedge_multiplicity)
    nullify(COMM % bedge_owner_rank)

    nullify(COMM % lsend_size)
    nullify(COMM % lrecv_size)
    nullify(COMM % lsend_perm)
    nullify(COMM % lrecv_perm)
    nullify(COMM % lscat_perm)
    nullify(COMM % matrix_ia)
    nullify(COMM % matrix_ja)
    nullify(COMM % matrix_aa)
    nullify(COMM % bface_size)
    nullify(COMM % bface_perm)

    nullify(COMM % full_row_send_neights)
    nullify(COMM % full_row_send_size)
    nullify(COMM % full_row_send_perm)
    nullify(COMM % full_row_recv_neights)
    nullify(COMM % full_row_recv_size)
    nullify(COMM % full_row_recv_perm)

    nullify(COMM % ghost_send_elem_size)
    nullify(COMM % ghost_send_elem_perm)
    nullify(COMM % ghost_recv_elem_size)
    nullify(COMM % ghost_recv_elem_perm)
    nullify(COMM % ghost_send_node_size)
    nullify(COMM % ghost_send_node_perm)
    nullify(COMM % ghost_recv_node_size)
    nullify(COMM % ghost_recv_node_perm)
    nullify(COMM % ghost_send_boun_size)
    nullify(COMM % ghost_send_boun_perm)
    nullify(COMM % ghost_recv_boun_size)
    nullify(COMM % ghost_recv_boun_perm)

  end subroutine PAR_INITIALIZE_COMMUNICATION_ARRAY_s

  subroutine PAR_INITIALIZE_COMMUNICATION_ARRAY_1(COMM,PAR_COMM_OPT)
    type(comm_data_par), pointer, intent(inout)          :: COMM(:)     !< Input communicator
    logical(lg),                  intent(in),   optional :: PAR_COMM_OPT
    integer(ip)                                          :: icomm

    if( associated(COMM) ) then
       do icomm = lbound(COMM,1),ubound(COMM,1)
          call PAR_INITIALIZE_COMMUNICATION_ARRAY_s(COMM(icomm),PAR_COMM_OPT)
       end do
    end if

  end subroutine PAR_INITIALIZE_COMMUNICATION_ARRAY_1

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    06/03/2014
  !> @brief   Copy a communication array
  !> @details Copy a communication array to another
  !
  !----------------------------------------------------------------------

  subroutine PAR_POINT_COMMUNICATION_ARRAY(COMM_IN,COMM_OUT)
    type(comm_data_par), intent(in)    :: COMM_IN(:)     !< Input communicator
    type(comm_data_par), intent(inout) :: COMM_OUT(:)    !< Output communicator

    COMM_OUT(1) % nneig                 =  COMM_IN(1) % nneig
    COMM_OUT(1) % npoi1                 =  COMM_IN(1) % npoi1
    COMM_OUT(1) % npoi2                 =  COMM_IN(1) % npoi2
    COMM_OUT(1) % npoi3                 =  COMM_IN(1) % npoi3
    COMM_OUT(1) % bound_dim             =  COMM_IN(1) % bound_dim
    COMM_OUT(1) % lsend_dim             =  COMM_IN(1) % lsend_dim
    COMM_OUT(1) % lrecv_dim             =  COMM_IN(1) % lrecv_dim
    COMM_OUT(1) % lscat_dim             =  COMM_IN(1) % lscat_dim
    COMM_OUT(1) % matrix_nzdom          =  COMM_IN(1) % matrix_nzdom
    COMM_OUT(1) % bface_dim             =  COMM_IN(1) % bface_dim
    COMM_OUT(1) % ghost_send_elem_dim   =  COMM_IN(1) % ghost_send_elem_dim
    COMM_OUT(1) % ghost_recv_elem_dim   =  COMM_IN(1) % ghost_recv_elem_dim
    COMM_OUT(1) % ghost_send_node_dim   =  COMM_IN(1) % ghost_send_node_dim
    COMM_OUT(1) % ghost_recv_node_dim   =  COMM_IN(1) % ghost_recv_node_dim
    COMM_OUT(1) % ghost_send_boun_dim   =  COMM_IN(1) % ghost_send_boun_dim
    COMM_OUT(1) % ghost_recv_boun_dim   =  COMM_IN(1) % ghost_recv_boun_dim

    COMM_OUT(1) % PAR_COMM_WORLD        =  COMM_IN(1) % PAR_COMM_WORLD
    COMM_OUT(1) % RANK4                 =  COMM_IN(1) % RANK4
    COMM_OUT(1) % neights               => COMM_IN(1) % neights
    COMM_OUT(1) % bound_size            => COMM_IN(1) % bound_size
    COMM_OUT(1) % bound_perm            => COMM_IN(1) % bound_perm
    COMM_OUT(1) % bound_scal            => COMM_IN(1) % bound_scal
    COMM_OUT(1) % bound_multiplicity    => COMM_IN(1) % bound_multiplicity
    COMM_OUT(1) % bound_owner_rank      => COMM_IN(1) % bound_owner_rank
    COMM_OUT(1) % node_number_in_owner  => COMM_IN(1) % node_number_in_owner
    COMM_OUT(1) % bound_matrix          => COMM_IN(1) % bound_matrix
    COMM_OUT(1) % bound_mat_halo_send   => COMM_IN(1) % bound_mat_halo_send
    COMM_OUT(1) % bound_mat_halo_recv   => COMM_IN(1) % bound_mat_halo_recv

    COMM_OUT(1) % bedge_size            => COMM_IN(1) % bedge_size
    COMM_OUT(1) % bedge_perm            => COMM_IN(1) % bedge_perm
    COMM_OUT(1) % bedge_adja            => COMM_IN(1) % bedge_adja
    COMM_OUT(1) % bedge_scal            => COMM_IN(1) % bedge_scal
    COMM_OUT(1) % bedge_dim             =  COMM_IN(1) % bedge_dim
    COMM_OUT(1) % bedge_multiplicity    => COMM_IN(1) % bedge_multiplicity
    COMM_OUT(1) % bedge_owner_rank      => COMM_IN(1) % bedge_owner_rank

    COMM_OUT(1) % lsend_size            => COMM_IN(1) % lsend_size
    COMM_OUT(1) % lrecv_size            => COMM_IN(1) % lrecv_size
    COMM_OUT(1) % lsend_perm            => COMM_IN(1) % lsend_perm
    COMM_OUT(1) % lrecv_perm            => COMM_IN(1) % lrecv_perm
    COMM_OUT(1) % lscat_perm            => COMM_IN(1) % lrecv_perm

    COMM_OUT(1) % matrix_ia             => COMM_IN(1) % matrix_ia
    COMM_OUT(1) % matrix_ja             => COMM_IN(1) % matrix_ja
    COMM_OUT(1) % matrix_aa             => COMM_IN(1) % matrix_aa

    COMM_OUT(1) % bface_size            => COMM_IN(1) % bface_size
    COMM_OUT(1) % bface_perm            => COMM_IN(1) % bface_perm

    COMM_OUT(1) % ghost_send_elem_size  => COMM_IN(1) % ghost_send_elem_size
    COMM_OUT(1) % ghost_send_elem_perm  => COMM_IN(1) % ghost_send_elem_perm
    COMM_OUT(1) % ghost_recv_elem_size  => COMM_IN(1) % ghost_recv_elem_size
    COMM_OUT(1) % ghost_recv_elem_perm  => COMM_IN(1) % ghost_recv_elem_perm
    COMM_OUT(1) % ghost_send_node_size  => COMM_IN(1) % ghost_send_node_size
    COMM_OUT(1) % ghost_send_node_perm  => COMM_IN(1) % ghost_send_node_perm
    COMM_OUT(1) % ghost_recv_node_size  => COMM_IN(1) % ghost_recv_node_size
    COMM_OUT(1) % ghost_recv_node_perm  => COMM_IN(1) % ghost_recv_node_perm

    COMM_OUT(1) % ghost_send_boun_size  => COMM_IN(1) % ghost_send_boun_size
    COMM_OUT(1) % ghost_send_boun_perm  => COMM_IN(1) % ghost_send_boun_perm
    COMM_OUT(1) % ghost_recv_boun_size  => COMM_IN(1) % ghost_recv_boun_size
    COMM_OUT(1) % ghost_recv_boun_perm  => COMM_IN(1) % ghost_recv_boun_perm

  end subroutine PAR_POINT_COMMUNICATION_ARRAY

  subroutine PAR_COPY_COMMUNICATION_ARRAY(COMM_IN,COMM_OUT,memor_opt)

    use def_master, only : kfl_paral
    type(comm_data_par), intent(inout)           :: COMM_IN       !< Input communicator
    type(comm_data_par), intent(inout)           :: COMM_OUT      !< Output communicator
    integer(8),          intent(inout), optional :: memor_opt(2)  !< Memory counter
    integer(8)                                   :: memor(2)

    if( present(memor_opt) ) then
       memor = memor_opt
    else
       memor = 0_8
    end if

    COMM_OUT % nneig                 =  COMM_IN % nneig
    COMM_OUT % npoi1                 =  COMM_IN % npoi1
    COMM_OUT % npoi2                 =  COMM_IN % npoi2
    COMM_OUT % npoi3                 =  COMM_IN % npoi3
    COMM_OUT % bound_dim             =  COMM_IN % bound_dim
    COMM_OUT % lsend_dim             =  COMM_IN % lsend_dim
    COMM_OUT % lrecv_dim             =  COMM_IN % lrecv_dim
    COMM_OUT % lscat_dim             =  COMM_IN % lscat_dim
    COMM_OUT % matrix_nzdom          =  COMM_IN % matrix_nzdom
    COMM_OUT % bface_dim             =  COMM_IN % bface_dim
    COMM_OUT % ghost_send_elem_dim   =  COMM_IN % ghost_send_elem_dim
    COMM_OUT % ghost_recv_elem_dim   =  COMM_IN % ghost_recv_elem_dim
    COMM_OUT % ghost_send_node_dim   =  COMM_IN % ghost_send_node_dim
    COMM_OUT % ghost_recv_node_dim   =  COMM_IN % ghost_recv_node_dim
    COMM_OUT % ghost_send_boun_dim   =  COMM_IN % ghost_send_boun_dim
    COMM_OUT % ghost_recv_boun_dim   =  COMM_IN % ghost_recv_boun_dim
                                               
    COMM_OUT % PAR_COMM_WORLD        =  COMM_IN % PAR_COMM_WORLD
    COMM_OUT % RANK4                 =  COMM_IN % RANK4
                                               
    COMM_OUT % bedge_dim             =  COMM_IN % bedge_dim

    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % neights               , COMM_OUT % neights,'DO_NOT_DEALLOCATE')             
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bound_size            , COMM_OUT % bound_size,'DO_NOT_DEALLOCATE')          
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bound_perm            , COMM_OUT % bound_perm,'DO_NOT_DEALLOCATE')          
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bound_scal            , COMM_OUT % bound_scal,'DO_NOT_DEALLOCATE')          
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bound_multiplicity    , COMM_OUT % bound_multiplicity,'DO_NOT_DEALLOCATE')  
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bound_owner_rank      , COMM_OUT % bound_owner_rank,'DO_NOT_DEALLOCATE')    
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % node_number_in_owner  , COMM_OUT % node_number_in_owner,'DO_NOT_DEALLOCATE')

    if( associated(COMM_IN % bound_matrix) )        call runend('YOU SHOULD CODE THESE LINES...')
    if( associated(COMM_IN % bound_mat_halo_send) ) call runend('YOU SHOULD CODE THESE LINES...')
    if( associated(COMM_IN % bound_mat_halo_recv) ) call runend('YOU SHOULD CODE THESE LINES...')
    !call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bound_matrix          , COMM_OUT % bound_matrix,'DO_NOT_DEALLOCATE')        
    !call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bound_mat_halo_send   , COMM_OUT % bound_mat_halo_send,'DO_NOT_DEALLOCATE') 
    !call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bound_mat_halo_recv   , COMM_OUT % bound_mat_halo_recv,'DO_NOT_DEALLOCATE') 

    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bedge_size            , COMM_OUT % bedge_size,'DO_NOT_DEALLOCATE')          
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bedge_perm            , COMM_OUT % bedge_perm,'DO_NOT_DEALLOCATE')          
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bedge_adja            , COMM_OUT % bedge_adja,'DO_NOT_DEALLOCATE')          
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bedge_scal            , COMM_OUT % bedge_scal,'DO_NOT_DEALLOCATE')          
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bedge_multiplicity    , COMM_OUT % bedge_multiplicity,'DO_NOT_DEALLOCATE')  
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bedge_owner_rank      , COMM_OUT % bedge_owner_rank,'DO_NOT_DEALLOCATE')  

    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % lsend_size            , COMM_OUT % lsend_size,'DO_NOT_DEALLOCATE')          
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % lrecv_size            , COMM_OUT % lrecv_size,'DO_NOT_DEALLOCATE')          
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % lsend_perm            , COMM_OUT % lsend_perm,'DO_NOT_DEALLOCATE')          
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % lrecv_perm            , COMM_OUT % lrecv_perm,'DO_NOT_DEALLOCATE')          
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % lscat_perm            , COMM_OUT % lrecv_perm,'DO_NOT_DEALLOCATE')          

    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % matrix_ia             , COMM_OUT % matrix_ia,'DO_NOT_DEALLOCATE')           
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % matrix_ja             , COMM_OUT % matrix_ja,'DO_NOT_DEALLOCATE')           
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % matrix_aa             , COMM_OUT % matrix_aa,'DO_NOT_DEALLOCATE')           

    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bface_size            , COMM_OUT % bface_size,'DO_NOT_DEALLOCATE')          
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % bface_perm            , COMM_OUT % bface_perm,'DO_NOT_DEALLOCATE')          

    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % ghost_send_elem_size  , COMM_OUT % ghost_send_elem_size,'DO_NOT_DEALLOCATE')
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % ghost_send_elem_perm  , COMM_OUT % ghost_send_elem_perm,'DO_NOT_DEALLOCATE')
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % ghost_recv_elem_size  , COMM_OUT % ghost_recv_elem_size,'DO_NOT_DEALLOCATE')
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % ghost_recv_elem_perm  , COMM_OUT % ghost_recv_elem_perm,'DO_NOT_DEALLOCATE')
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % ghost_send_node_size  , COMM_OUT % ghost_send_node_size,'DO_NOT_DEALLOCATE')
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % ghost_send_node_perm  , COMM_OUT % ghost_send_node_perm,'DO_NOT_DEALLOCATE')
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % ghost_recv_node_size  , COMM_OUT % ghost_recv_node_size,'DO_NOT_DEALLOCATE')
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % ghost_recv_node_perm  , COMM_OUT % ghost_recv_node_perm,'DO_NOT_DEALLOCATE')

    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % ghost_send_boun_size  , COMM_OUT % ghost_send_boun_size,'DO_NOT_DEALLOCATE')
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % ghost_send_boun_perm  , COMM_OUT % ghost_send_boun_perm,'DO_NOT_DEALLOCATE')
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % ghost_recv_boun_size  , COMM_OUT % ghost_recv_boun_size,'DO_NOT_DEALLOCATE')
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_COPY_COMMUNICATION_ARRAY',COMM_IN % ghost_recv_boun_perm  , COMM_OUT % ghost_recv_boun_perm,'DO_NOT_DEALLOCATE')
    
    if( present(memor_opt) ) memor_opt = memor

  end subroutine PAR_COPY_COMMUNICATION_ARRAY

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    06/03/2014
  !> @brief   Create a sub communication array
  !> @details Given a communicaiton array structure, create a sub
  !           communication structure for a subset of nodes
  !
  !----------------------------------------------------------------------

  subroutine PAR_SUB_COMMUNICATION_ARRAY(COMM_IN,COMM_OUT,mask)
    type(comm_data_par),          intent(in)    :: COMM_IN        !< Input communicator
    type(comm_data_par),          intent(inout) :: COMM_OUT       !< Output communicator
    integer(ip),         pointer, intent(in)    :: mask(:)        !< mask (=1 to consider node)
    integer(ip),         pointer                :: bound_perm(:)
    integer(ip),         pointer                :: bound_size(:)
    integer(ip),         pointer                :: neights(:)
    integer(ip)                                 :: nneig,bound_dim
    integer(ip)                                 :: ineig,jj,jneig
    integer(ip)                                 :: ipoin

    COMM_OUT % PAR_COMM_WORLD = COMM_IN % PAR_COMM_WORLD
    COMM_OUT % RANK4          = COMM_IN % RANK4

    if( COMM_IN % bound_dim > 0 ) then

       nullify(bound_perm)
       nullify(bound_size)
       nullify(neights)

       allocate( bound_perm(COMM_IN % bound_dim) )
       allocate( bound_size(COMM_IN % nneig)     )
       allocate( neights   (COMM_IN % nneig)     )

       nneig     = 0
       bound_dim = 0

       do jj = 1,COMM_IN % bound_dim
          bound_perm(jj) = COMM_IN % bound_perm(jj)
       end do
       do ineig = 1,COMM_IN % nneig
          bound_size(ineig) = 0
       end do

       do ineig = 1,COMM_IN % nneig
          do jj = COMM_IN % bound_size(ineig),COMM_IN % bound_size(ineig+1)-1
             ipoin = COMM_IN % bound_perm(jj)
             if( mask(ipoin) > 0 ) then
                bound_size(ineig) = bound_size(ineig) + 1
             else
                bound_perm(jj) = 0
             end if
          end do
          bound_dim = bound_dim + bound_size(ineig)

          if( bound_size(ineig) > 0 ) nneig = nneig + 1
       end do
       !
       ! Allocate interzone communication array
       !
       COMM_OUT % nneig     = nneig
       COMM_OUT % bound_dim = bound_dim
       COMM_OUT % npoi1     = COMM_IN % npoi1
       COMM_OUT % npoi2     = COMM_IN % npoi2
       COMM_OUT % npoi3     = COMM_IN % npoi3
       allocate( COMM_OUT % neights(nneig) )
       allocate( COMM_OUT % bound_perm(bound_dim) )
       allocate( COMM_OUT % bound_size(nneig+1) )
       !
       ! Permutation array BOUND_PERM(1:BOUND_DIM)
       !
       jneig = 0
       bound_dim = 0
       do ineig = 1,COMM_IN % nneig
          if( bound_size(ineig) > 0 ) then

             jneig = jneig + 1
             COMM_OUT % neights(jneig) = COMM_IN % neights(ineig)

             do jj = COMM_IN % bound_size(ineig),COMM_IN % bound_size(ineig+1)-1
                ipoin = bound_perm(jj)
                if( ipoin > 0 ) then
                   bound_dim = bound_dim + 1
                   COMM_OUT % bound_perm(bound_dim) = ipoin
                end if
             end do

          end if
       end do
       !
       ! Construct linked list BOUND_SIZE(1:NNEIG+1)
       !
       COMM_OUT % bound_size(1) = 1
       jneig = 0
       do ineig = 1,COMM_IN % nneig
          if( bound_size(ineig) > 0 ) then
             jneig = jneig + 1
             COMM_OUT % bound_size(jneig+1) = &
                  COMM_OUT % bound_size(jneig) + bound_size(ineig)
          end if
       end do

       deallocate( bound_perm )
       deallocate( bound_size )
       deallocate( neights    )

    end if

  end subroutine PAR_SUB_COMMUNICATION_ARRAY

  function NODE_IN_NEIGHBOR_P(ipoin,ineig,commu)
    implicit none
    integer(ip),                  intent(in) :: ipoin     !< Node
    integer(ip),                  intent(in) :: ineig     !< Neighbor
    type(comm_data_par), pointer, intent(in) :: commu(:)  !< Generic communication array
    logical(lg)                              :: NODE_IN_NEIGHBOR_P
    integer(ip)                              :: ii


    NODE_IN_NEIGHBOR_P = .false.
    loop_nodes: do ii = commu(1) % bound_size(ineig),commu(1) % bound_size(ineig+1)-1
       if( ipoin == commu(1) % bound_perm(ii) )then
          NODE_IN_NEIGHBOR_P = .true.
          exit loop_nodes
       end if
    end do loop_nodes

  end function NODE_IN_NEIGHBOR_P

  function NODE_IN_NEIGHBOR_s(ipoin,ineig,commu)
    implicit none
    integer(ip),         intent(in) :: ipoin  !< Node
    integer(ip),         intent(in) :: ineig  !< Neighbor
    type(comm_data_par), intent(in) :: commu  !< Generic communication array
    logical(lg)                     :: NODE_IN_NEIGHBOR_s
    integer(ip)                     :: ii

    NODE_IN_NEIGHBOR_s = .false.
    loop_nodes: do ii = commu % bound_size(ineig),commu % bound_size(ineig+1)-1
       if( ipoin == commu % bound_perm(ii) )then
          NODE_IN_NEIGHBOR_s = .true.
          exit loop_nodes
       end if
    end do loop_nodes

  end function NODE_IN_NEIGHBOR_s

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    08/02/2019
  !> @brief   Transpose communication array
  !> @details Transpose communication array
  !
  !----------------------------------------------------------------------

  subroutine PAR_TRANSPOSE_COMMUNICATION_ARRAY(COMM,memor_opt)

    type(comm_data_par), intent(inout)           :: COMM          !< Input communicator
    integer(8),          intent(inout), optional :: memor_opt(2)  !< Memory counter
    integer(8)                                   :: memor(2)
    integer(ip)                                  :: lsend_dim_copy
    integer(ip),  pointer                        :: lsend_size_copy(:)
    integer(ip),  pointer                        :: lsend_perm_copy(:)

    nullify(lsend_size_copy)
    nullify(lsend_perm_copy)
    
    if( present(memor_opt) ) then
       memor = memor_opt
    else
       memor = 0_8
    end if
    !
    ! SEND => SEND_COPY 
    !
    lsend_dim_copy = COMM % lsend_dim
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_TRANSPOSE_COMMUNICATION_ARRAY',COMM % lsend_size,lsend_size_copy)
    call memory_copy(memor,'COMM % LSEND_PERM','PAR_TRANSPOSE_COMMUNICATION_ARRAY',COMM % lsend_perm,lsend_perm_copy)
    !
    ! RECV => SEND
    !
    COMM % lsend_dim = COMM % lrecv_dim
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_TRANSPOSE_COMMUNICATION_ARRAY',COMM % lrecv_size,COMM % lsend_size)
    call memory_copy(memor,'COMM % LSEND_PERM','PAR_TRANSPOSE_COMMUNICATION_ARRAY',COMM % lrecv_perm,COMM % lsend_perm)
    !
    ! SEND_COPY => RECV
    !    
    COMM % lrecv_dim = lsend_dim_copy
    call memory_copy(memor,'COMM % LSEND_SIZE','PAR_TRANSPOSE_COMMUNICATION_ARRAY',lsend_size_copy,COMM % lrecv_size)
    call memory_copy(memor,'COMM % LSEND_PERM','PAR_TRANSPOSE_COMMUNICATION_ARRAY',lsend_perm_copy,COMM % lrecv_perm)
       
    if( present(memor_opt) ) memor_opt = memor

  end subroutine PAR_TRANSPOSE_COMMUNICATION_ARRAY

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    06/03/2014
  !> @brief   Ghost geometry communicator
  !> @details Ghost geometry communicator
  !>
  !----------------------------------------------------------------------

  !subroutine par_ghost_communication_array(COMM)
  !  type(comm_data_par), intent(inout) :: COMM             !< Input communicator
  !  integer(ip)                        :: PAR_CURRENT_RANK
  !  integer(ip)                        :: PAR_CURRENT_SIZE
  !  integer(ip)                        :: dom_i,dom_j
  !  integer(ip)                        :: ghost_elem_dim_ineig,ii
  !  integer(ip),         pointer       :: permu(:)
  !  PAR_COMM_CURRENT = COMM % PAR_COMM_WORLD
  !  call PAR_COMM_RANK_AND_SIZE(PAR_COMM_CURRENT,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)
  !end subroutine par_ghost_communication_array

end module mod_parall
!
!> @}
!-----------------------------------------------------------------------
