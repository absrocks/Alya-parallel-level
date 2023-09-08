!------------------------------------------------------------------------
!> @defgroup Communication_Toolbox
!> Toolbox for MPI communication, bridge to MPI
!> @{
!> @name    Parallelization toolbox
!> @file    mod_communications.f90
!> @author  Guillaume Houzeaux
!> @date    28/06/2012
!> @brief   ToolBox for parallel communications
!> @details ToolBox for parallel communications
!------------------------------------------------------------------------

module mod_communications

  use def_kintyp, only : ip,rp,lg,r1p
  use def_domain, only : npoin
  use def_domain, only : npoin_2
  use def_domain, only : nelem_2
  use def_domain, only : nboun_2
  use def_master, only : npoi1,kfl_paral
  use def_master, only : IPARALL,IMASTER,ISLAVE,INOTSLAVE,ISEQUEN
  use def_master, only : INOTMASTER
  use def_master, only : comm_data_par
  use def_master, only : current_code,current_zone,current_subd
  use mod_maths,  only : maths_heap_sort
  use mod_memory, only : memory_alloca
  use mod_memory, only : memory_deallo
  use mod_memory, only : memory_size
  use mod_parall, only : par_code_zone_subd_to_color
  use mod_parall, only : PAR_COMM_COLOR
  use mod_parall, only : PAR_COMM_MY_CODE
  use mod_parall, only : PAR_COMM_MY_CODE_WM4
  use mod_parall, only : PAR_COMM_COLOR_ARRAY
  use mod_parall, only : PAR_COMM_MY_CODE_ARRAY
  use mod_parall, only : PAR_MY_CODE_RANK, PAR_MY_WORLD_RANK
  use mod_parall, only : PAR_COMM_WORLD
  use mod_parall, only : PAR_COMM_CURRENT
  use mod_parall, only : PAR_INTEGER
  use mod_parall, only : PAR_COMM_SFC
  use mod_parall, only : PAR_COMM_SFC_WM
  use mod_parall, only : PAR_COMM_MPIO
  use mod_parall, only : PAR_COMM_MPIO_WM
  use mod_parall, only : commd
  use mod_parall, only : color_target
  use mod_parall, only : color_source
  use mod_parall, only : par_memor
  use mod_parall, only : PAR_CODE_SIZE
  use mod_std

#ifdef EXTRAE
  use extrae_module
#endif

  implicit none

  private

#ifndef MPI_OFF
  include 'mpif.h'
  integer(4) :: status(MPI_STATUS_SIZE)
#endif

  type non_blocking_typ
     integer(4),  pointer  :: request4(:)
     integer(4)            :: count4
  end type non_blocking_typ
  logical(lg),            allocatable :: tmp_lsend(:)
  logical(lg),            allocatable :: tmp_lrecv(:)
  integer(ip),            allocatable :: tmp_isend(:)
  integer(ip),            allocatable :: tmp_irecv(:)
  real(rp),               allocatable :: tmp_rsend(:)
  real(rp),               allocatable :: tmp_rrecv(:)
  integer(4),             allocatable :: ireq4(:)
  integer(4)                          :: ireq41(1)
  type(non_blocking_typ), pointer     :: non_blocking(:)
  integer(ip)                         :: inonblocking
  real(rp),               pointer     :: yy_non_blocking(:)
  integer(ip),            parameter   :: ON_NODES = 1
  integer(ip),            parameter   :: ON_EDGES = 2
  !
  ! All reduce operations:
  !
  ! PAR_SUM
  ! PAR_MAX
  ! PAR_MIN
  !
  interface PAR_SUM
     module procedure PAR_SUM_IP_0,PAR_SUM_IP_s,PAR_SUM_IP_1,PAR_SUM_IP_2,PAR_SUM_IP_3,PAR_SUM_IP_03,&
          &           PAR_SUM_RP_0,PAR_SUM_RP_0b,PAR_SUM_RP_s,PAR_SUM_RP_1,PAR_SUM_RP_2,PAR_SUM_RP_3,&
          &           PAR_SUM_RP_02,&
          &           PAR_SUM_CX_0,PAR_SUM_CX_s,PAR_SUM_CX_1,PAR_SUM_CX_2
  end interface PAR_SUM
  interface PAR_MAX
     module procedure PAR_MAX_4_s,PAR_MAX_8_s,PAR_MAX_IP_0,PAR_MAX_IP_1,PAR_MAX_IP_2,&
          &           PAR_MAX_RP_s,PAR_MAX_RP_0,PAR_MAX_RP_1,PAR_MAX_RP_2,PAR_MAX_RP_3,&
          &           PAR_MAX_CX_s,PAR_MAX_CX_0,PAR_MAX_CX_1,PAR_MAX_CX_2
  end interface PAR_MAX
  interface PAR_MIN
     module procedure PAR_MIN_4_s ,PAR_MIN_8_s ,PAR_MIN_IP_0,PAR_MIN_IP_1,PAR_MIN_IP_2,&
          &           PAR_MIN_RP_s,PAR_MIN_RP_0,PAR_MIN_RP_1,PAR_MIN_RP_2,&
          &           PAR_MIN_CX_s,PAR_MIN_CX_0,PAR_MIN_CX_1,PAR_MIN_CX_2
  end interface PAR_MIN
  interface PAR_OR
     module procedure PAR_OR_LG_s
  end interface PAR_OR
  interface PAR_AVERAGE
     module procedure PAR_AVERAGE_IP_s,PAR_AVERAGE_IP_0,&
          &           PAR_AVERAGE_RP_s,PAR_AVERAGE_RP_0
  end interface PAR_AVERAGE
  interface PAR_LOAD_BALANCE
     module procedure PAR_LOAD_BALANCE_RP_0,PAR_LOAD_BALANCE_RP_s
  end interface PAR_LOAD_BALANCE
  !
  ! Reduction operations in the Alya world (PAR_COMM_WORLD) including masters
  !
  interface PAR_MAX_ALL
     module procedure PAR_MAX_ALL_IP_s, PAR_MAX_ALL_IP_1
  end interface PAR_MAX_ALL
  interface PAR_SUM_ALL
     module procedure PAR_SUM_ALL_IP_1, PAR_SUM_ALL_IP_3
  end interface PAR_SUM_ALL
  !
  ! All to all
  !
  interface PAR_ALLTOALL
     module procedure PAR_ALLTOALL_IP_0,&
          &           PAR_ALLTOALL_IP_1
  end interface PAR_ALLTOALL
  !
  ! Exchange array between neighbors
  !
  interface PAR_ARRAY_EXCHANGE
     module procedure PAR_ARRAY_EXCHANGE_IP
  end interface PAR_ARRAY_EXCHANGE
  !
  ! Operation on arrays between all partitions
  !
  interface PAR_ALL_TO_ALL_ARRAY_OPERATION
     module procedure PAR_ALL_TO_ALL_ARRAY_OPERATION_IP
  end interface PAR_ALL_TO_ALL_ARRAY_OPERATION
  !
  ! Exchange between slaves:
  ! type:  1. nodes  2. faces    3. fringe nodes
  ! kind:  1. real   2. integer  3. Complex
  ! dim:   1. x(:)   2. x(:,:)   3. x(:,:,:)   4. n,x(n)
  ! what:  1. SUM    2. MAX      3. MIN
  ! where: in my code, in my zone
  !
  interface PAR_INTERFACE_NODE_EXCHANGE
     module procedure PAR_INTERFACE_NODE_EXCHANGE_IP_0,  &
          &           PAR_INTERFACE_NODE_EXCHANGE_IP_1,  &
          &           PAR_INTERFACE_NODE_EXCHANGE_IP_2,  &
          &           PAR_INTERFACE_NODE_EXCHANGE_IP_3,  &
          &           PAR_INTERFACE_NODE_EXCHANGE_IP_2b, &
          &           PAR_INTERFACE_NODE_EXCHANGE_RP_00, &
          &           PAR_INTERFACE_NODE_EXCHANGE_RP_0,  &
          &           PAR_INTERFACE_NODE_EXCHANGE_RP_1,  &
          &           PAR_INTERFACE_NODE_EXCHANGE_RP_2,  &
          &           PAR_INTERFACE_NODE_EXCHANGE_RP_3,  &
          &           PAR_INTERFACE_NODE_EXCHANGE_RP_2b, &
          &           PAR_INTERFACE_NODE_EXCHANGE_LG_1
  end interface PAR_INTERFACE_NODE_EXCHANGE
  interface PAR_INTERFACE_EDGE_EXCHANGE
     module procedure PAR_INTERFACE_EDGE_EXCHANGE_IP_0,  &
          &           PAR_INTERFACE_EDGE_EXCHANGE_IP_1,  &
          &           PAR_INTERFACE_EDGE_EXCHANGE_IP_2,  &
          &           PAR_INTERFACE_EDGE_EXCHANGE_IP_3,  &
          &           PAR_INTERFACE_EDGE_EXCHANGE_IP_2b, &
          &           PAR_INTERFACE_EDGE_EXCHANGE_RP_00, &
          &           PAR_INTERFACE_EDGE_EXCHANGE_RP_0,  &
          &           PAR_INTERFACE_EDGE_EXCHANGE_RP_1,  &
          &           PAR_INTERFACE_EDGE_EXCHANGE_RP_2,  &
          &           PAR_INTERFACE_EDGE_EXCHANGE_RP_3,  &
          &           PAR_INTERFACE_EDGE_EXCHANGE_RP_2b, &
          &           PAR_INTERFACE_EDGE_EXCHANGE_LG_1
  end interface PAR_INTERFACE_EDGE_EXCHANGE
  !
  ! Exchange on own nodes
  !
  interface PAR_INTERFACE_OWN_NODE_EXCHANGE
     module procedure PAR_INTERFACE_OWN_NODE_EXCHANGE_IP_0,  &
          &           PAR_INTERFACE_OWN_NODE_EXCHANGE_IP_1,  &
          &           PAR_INTERFACE_OWN_NODE_EXCHANGE_IP_2,  &
          &           PAR_INTERFACE_OWN_NODE_EXCHANGE_RP_00, &
          &           PAR_INTERFACE_OWN_NODE_EXCHANGE_RP_0,  &
          &           PAR_INTERFACE_OWN_NODE_EXCHANGE_RP_1,  &
          &           PAR_INTERFACE_OWN_NODE_EXCHANGE_RP_2
  end interface PAR_INTERFACE_OWN_NODE_EXCHANGE
  !
  ! Interface node matrix exchange
  !
  interface PAR_INTERFACE_MATRIX_EXCHANGE
     module procedure PAR_INTERFACE_MATRIX_EXCHANGE_WHERE,&
          &           PAR_INTERFACE_MATRIX_EXCHANGE_COMM, &
          &           PAR_INTERFACE_MATRIX_EXCHANGE_COMM_SHAPE
  end interface PAR_INTERFACE_MATRIX_EXCHANGE
  !
  ! Interface node matrix with halos exchange
  !
  interface PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE
     module procedure PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE_RP_1,&
          &           PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE_RP_2
  end interface PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE
  !
  ! Exchange between assymetric
  !
  interface PAR_COUPLING_NODE_EXCHANGE
     module procedure PAR_COUPLING_NODE_EXCHANGE_RP_0,  &
          &           PAR_COUPLING_NODE_EXCHANGE_RP_1,  &
          &           PAR_COUPLING_NODE_EXCHANGE_RP_1b, &
          &           PAR_COUPLING_NODE_EXCHANGE_RP_2,  &
          &           PAR_COUPLING_NODE_EXCHANGE_RP_3,  &
          &           PAR_COUPLING_NODE_EXCHANGE_RP_2b
  end interface PAR_COUPLING_NODE_EXCHANGE
  !
  ! Exchange between slaves on ghost elements/boundaries/nodes
  !
  ! +----+----+----+....o
  !             ||   /|
  !             \/   ||
  !           o....+----+----+----+
  !
  !
  interface PAR_GHOST_ELEMENT_EXCHANGE
     module procedure PAR_GHOST_ELEMENT_EXCHANGE_IP_0,  &
          &           PAR_GHOST_ELEMENT_EXCHANGE_IP_1,  &
          &           PAR_GHOST_ELEMENT_EXCHANGE_IP_2,  &
          &           PAR_GHOST_ELEMENT_EXCHANGE_IP_3,  &
          &           PAR_GHOST_ELEMENT_EXCHANGE_IP_2b, &
          &           PAR_GHOST_ELEMENT_EXCHANGE_RP_00, &
          &           PAR_GHOST_ELEMENT_EXCHANGE_RP_0,  &
          &           PAR_GHOST_ELEMENT_EXCHANGE_RP_1,  &
          &           PAR_GHOST_ELEMENT_EXCHANGE_RP_2,  &
          &           PAR_GHOST_ELEMENT_EXCHANGE_RP_3,  &
          &           PAR_GHOST_ELEMENT_EXCHANGE_RP_2b
  end interface PAR_GHOST_ELEMENT_EXCHANGE
  interface PAR_GHOST_BOUNDARY_EXCHANGE
     module procedure PAR_GHOST_BOUNDARY_EXCHANGE_IP_0,  &
          &           PAR_GHOST_BOUNDARY_EXCHANGE_IP_1,  &
          &           PAR_GHOST_BOUNDARY_EXCHANGE_IP_2,  &
          &           PAR_GHOST_BOUNDARY_EXCHANGE_IP_3,  &
          &           PAR_GHOST_BOUNDARY_EXCHANGE_IP_2b, &
          &           PAR_GHOST_BOUNDARY_EXCHANGE_RP_00, &
          &           PAR_GHOST_BOUNDARY_EXCHANGE_RP_0,  &
          &           PAR_GHOST_BOUNDARY_EXCHANGE_RP_1,  &
          &           PAR_GHOST_BOUNDARY_EXCHANGE_RP_2,  &
          &           PAR_GHOST_BOUNDARY_EXCHANGE_RP_3,  &
          &           PAR_GHOST_BOUNDARY_EXCHANGE_RP_2b
  end interface PAR_GHOST_BOUNDARY_EXCHANGE
  interface PAR_GHOST_NODE_EXCHANGE
     module procedure PAR_GHOST_NODE_EXCHANGE_IP_0,  &
          &           PAR_GHOST_NODE_EXCHANGE_IP_1,  &
          &           PAR_GHOST_NODE_EXCHANGE_IP_2,  &
          &           PAR_GHOST_NODE_EXCHANGE_IP_3,  &
          &           PAR_GHOST_NODE_EXCHANGE_IP_2b, &
          &           PAR_GHOST_NODE_EXCHANGE_RP_00, &
          &           PAR_GHOST_NODE_EXCHANGE_RP_0,  &
          &           PAR_GHOST_NODE_EXCHANGE_RP_1,  &
          &           PAR_GHOST_NODE_EXCHANGE_RP_2,  &
          &           PAR_GHOST_NODE_EXCHANGE_RP_3,  &
          &           PAR_GHOST_NODE_EXCHANGE_RP_2b
  end interface PAR_GHOST_NODE_EXCHANGE
  !
  ! Exchange between slaves on ghost elements/boundaries/nodes
  !
  ! +----+----+----+....o
  !             /|   ||
  !             ||   \/
  !           o....+----+----+----+
  !
  !
  interface PAR_FROM_GHOST_ELEMENT_EXCHANGE
     module procedure PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_00, &
          &           PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_0,  &
          &           PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_1,  &
          &           PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_2,  &
          &           PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_3,  &
          &           PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_2b, &
          &           PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_00, &
          &           PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_0,  &
          &           PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_1,  &
          &           PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_2,  &
          &           PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_3,  &
          &           PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_2b
  end interface PAR_FROM_GHOST_ELEMENT_EXCHANGE
  interface PAR_FROM_GHOST_BOUNDARY_EXCHANGE
     module procedure PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_00, &
          &           PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_0,  &
          &           PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_1,  &
          &           PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_2,  &
          &           PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_3,  &
          &           PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_2b, &
          &           PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_00, &
          &           PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_0,  &
          &           PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_1,  &
          &           PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_2,  &
          &           PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_3,  &
          &           PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_2b
  end interface PAR_FROM_GHOST_BOUNDARY_EXCHANGE
  interface PAR_FROM_GHOST_NODE_EXCHANGE
     module procedure PAR_FROM_GHOST_NODE_EXCHANGE_IP_00, &
          &           PAR_FROM_GHOST_NODE_EXCHANGE_IP_0,  &
          &           PAR_FROM_GHOST_NODE_EXCHANGE_IP_1,  &
          &           PAR_FROM_GHOST_NODE_EXCHANGE_IP_2,  &
          &           PAR_FROM_GHOST_NODE_EXCHANGE_IP_3,  &
          &           PAR_FROM_GHOST_NODE_EXCHANGE_IP_2b
  end interface PAR_FROM_GHOST_NODE_EXCHANGE
  !
  ! Send
  !
  interface PAR_SEND
     module procedure PAR_SEND_IP_s,&
          &           PAR_SEND_RP_s
  end interface PAR_SEND
  !
  ! Receive
  !
  interface PAR_RECEIVE
     module procedure PAR_RECEIVE_IP_s,&
          &           PAR_RECEIVE_RP_s
  end interface PAR_RECEIVE
  !
  ! Send/receive
  !
  interface PAR_SEND_RECEIVE
     module procedure PAR_SEND_RECEIVE_IP_s,&
          &           PAR_SEND_RECEIVE_IP_0,&
          &           PAR_SEND_RECEIVE_IP_1,&
          &           PAR_SEND_RECEIVE_IP_2,&
          &           PAR_SEND_RECEIVE_IP_3,&
          &           PAR_SEND_RECEIVE_RP_s,&
          &           PAR_SEND_RECEIVE_RP_0,&
          &           PAR_SEND_RECEIVE_RP_1,&
          &           PAR_SEND_RECEIVE_RP_2,&
          &           PAR_SEND_RECEIVE_RP_3
  end interface PAR_SEND_RECEIVE
  interface PAR_SEND_RECEIVE_TO_ALL
     module procedure PAR_SEND_RECEIVE_TO_ALL_RP_s,&
          &           PAR_SEND_RECEIVE_TO_ALL_RP_0,&
          &           PAR_SEND_RECEIVE_TO_ALL_RP_1,&
          &           PAR_SEND_RECEIVE_TO_ALL_RP_1c,&
          &           PAR_SEND_RECEIVE_TO_ALL_RP_2,&
          &           PAR_SEND_RECEIVE_TO_ALL_RP_3,&
          &           PAR_SEND_RECEIVE_TO_ALL_IP_1,&
          &           PAR_SEND_RECEIVE_TO_ALL_IP_1c,&
          &           PAR_SEND_RECEIVE_TO_ALL_IP_2,&
          &           PAR_SEND_RECEIVE_TO_ALL_IP_0
  end interface PAR_SEND_RECEIVE_TO_ALL
  !
  ! Operations with neighbors
  !
  interface PAR_POINT_TO_POINT_ARRAY_OPERATION
     module procedure PAR_POINT_TO_POINT_ARRAY_OPERATION_IP_1
  end interface PAR_POINT_TO_POINT_ARRAY_OPERATION
  !
  ! Split
  !
  interface PAR_COMM_SPLIT
     module procedure PAR_COMM_SPLIT4,PAR_COMM_SPLIT8
  end interface PAR_COMM_SPLIT
  !
  ! Bridge to broadcast data from master to slaves
  !
  interface PAR_EXCHANGE
     module procedure PAR_EXCHANGE_IP_s,PAR_EXCHANGE_IP_0,PAR_EXCHANGE_IP_1,PAR_EXCHANGE_IP_2,&
          &           PAR_EXCHANGE_RP_s,PAR_EXCHANGE_RP_0,PAR_EXCHANGE_RP_02,PAR_EXCHANGE_RP_1,PAR_EXCHANGE_RP_2,&
          &           PAR_EXCHANGE_CH,PAR_EXCHANGE_CH_1,&
          &           PAR_EXCHANGE_LG_s,PAR_EXCHANGE_LG_0,PAR_EXCHANGE_LG_1
  end interface PAR_EXCHANGE
  !
  ! Communicator operations
  !
  interface PAR_COMM_RANK_AND_SIZE
     module procedure PAR_COMM_RANK_AND_SIZE_4,PAR_COMM_RANK_AND_SIZE_4W,&
          &           PAR_COMM_RANK_AND_SIZE_41W,&
          &           PAR_COMM_RANK_AND_SIZE_8,PAR_COMM_RANK_AND_SIZE_8W
  end interface PAR_COMM_RANK_AND_SIZE
  !
  ! GATHER
  !
  interface PAR_GATHER
     module procedure PAR_GATHER_CHARACTER,&
          &           PAR_GATHER_IP_s4,PAR_GATHER_IP_s8,PAR_GATHER_IP_s48,PAR_GATHER_IP_14,PAR_GATHER_IP_18,&
          &           PAR_GATHER_IP_12,PAR_GATHER_IP_23,&
          &           PAR_GATHER_RP_s,PAR_GATHER_RP_1,PAR_GATHER_RP_12
  end interface PAR_GATHER
  !
  ! GATHERV
  !
  interface PAR_GATHERV
     module procedure PAR_GATHERV_RP_1,PAR_GATHERV_RP_21,PAR_GATHERV_RP_22,PAR_GATHERV_RP_0,&
          &           PAR_GATHERV_IP_1,PAR_GATHERV_IP_21, PAR_GATHERV_IP_22,&
                      PAR_GATHERV_RP_21_SEND,PAR_GATHERV_RP_22_SEND,&
          &           PAR_GATHERV_IP_1_SEND,PAR_GATHERV_IP_21_SEND,PAR_GATHERV_IP_22_SEND
  end interface PAR_GATHERV
  !
  ! ALLGATHERV
  !
  interface PAR_ALLGATHERV
     module procedure PAR_ALLGATHERV_IP4  ,PAR_ALLGATHERV_IP8  ,&
          &           PAR_ALLGATHERV_IP4_2,&
          &           PAR_ALLGATHERV_RP_1,&
          &           PAR_ALLGATHERV_RP_18,&
          &           PAR_ALLGATHERV_RP_2,&
          &           PAR_ALLGATHERV_RP_3
  end interface PAR_ALLGATHERV
  !
  ! ALLGATHER
  !
  interface PAR_ALLGATHER
     module procedure PAR_ALLGATHER_s4,PAR_ALLGATHER_s8,&
          &           PAR_ALLGATHER_IP_14,PAR_ALLGATHER_IP_18,&
          &           PAR_ALLGATHER_RP_0,PAR_ALLGATHER_RP_2,&
          &           PAR_ALLGATHER_RP_02,PAR_ALLGATHER_LG
  end interface PAR_ALLGATHER
  !
  ! SCATTER
  !
  interface PAR_SCATTER
     module procedure PAR_SCATTER_IP_s
  end interface PAR_SCATTER
  !
  ! SCATTERV
  !
  interface PAR_SCATTERV
     module procedure PAR_SCATTERV_IP_1,PAR_SCATTERV_IP_2,&
          &           PAR_SCATTERV_RP_1,PAR_SCATTERV_RP_2,PAR_SCATTERV_RP_0,&
          &           PAR_SCATTERV_IP_1_RCV,PAR_SCATTERV_IP_2_RCV,&
          &           PAR_SCATTERV_R4_2_RCV, PAR_SCATTERV_R8_2_RCV, PAR_SCATTERV_R16_2_RCV
  end interface PAR_SCATTERV
  !
  ! BORADCAST
  !
  interface PAR_BROADCAST_IP
     module procedure PAR_BROADCAST_I4,PAR_BROADCAST_I8
  end interface PAR_BROADCAST_IP

  interface PAR_BROADCAST
     module procedure PAR_BROADCAST_IP_04,PAR_BROADCAST_IP_08,&
          &           PAR_BROADCAST_I4_s,PAR_BROADCAST_I8_s,&
          &           PAR_BROADCAST_I4_1,PAR_BROADCAST_I8_1,&
          &           PAR_BROADCAST_R8_s,PAR_BROADCAST_R8_0,PAR_BROADCAST_R8_1,&
          &           PAR_BROADCAST_R4_s,PAR_BROADCAST_R4_0,PAR_BROADCAST_R4_1,&
          &           PAR_BROADCAST_R16_s,PAR_BROADCAST_R16_0,PAR_BROADCAST_R16_1,&
          &           PAR_BROADCAST_LG_s,PAR_BROADCAST_LG_1,&
          &           PAR_BROADCAST_CH
  end interface PAR_BROADCAST
  !
  ! ALL_TO_ALLV
  !
  interface PAR_ALLTOALLV
     module procedure PAR_ALLTOALLV_IP_1,&
          &           PAR_ALLTOALLV_IP_2,&
          &           PAR_ALLTOALLV_RP_1,&
          &           PAR_ALLTOALLV_RP_2
  end interface PAR_ALLTOALLV
  !
  ! Start non-blocking communications
  !
  interface PAR_START_NON_BLOCKING_COMM
     module procedure PAR_START_NON_BLOCKING_COMM_4,&
          &           PAR_START_NON_BLOCKING_COMM_4_OPT,&
          &           PAR_START_NON_BLOCKING_COMM_8
  end interface PAR_START_NON_BLOCKING_COMM
  !
  ! Public functions
  !
  public :: PAR_SUM                            ! AllReduce SUM
  public :: PAR_MAX                            ! AllReduce MAX
  public :: PAR_MIN                            ! AllReduce MIN
  public :: PAR_AVERAGE                        ! Average value
  public :: PAR_LOAD_BALANCE                   ! Load balance: average / max
  public :: PAR_SUM_ALL                        ! AllReduce SUM in the Alya world (PAR_COMM_WORLD)
  public :: PAR_MAX_ALL                        ! AllReduce MAX in the Alya world (PAR_COMM_WORLD)
  public :: PAR_ALLTOALL                       ! Alltoall
  public :: PAR_INTERFACE_NODE_EXCHANGE        ! Interface nodal exchange (send/receive)
  public :: PAR_INTERFACE_OWN_NODE_EXCHANGE    ! Interface OWN node exchange (send/receive)
  public :: PAR_INTERFACE_EDGE_EXCHANGE        ! Interface edge exchange (send/receive)
  public :: PAR_INTERFACE_MATRIX_EXCHANGE      ! Matrix exchange on interface nodes
  public :: PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE ! Matrix with halos exchange on interface nodes
  public :: PAR_COUPLING_NODE_EXCHANGE         ! Coupling nodal exchange (send/receive)
  public :: PAR_GHOST_ELEMENT_EXCHANGE         ! Ghost element exchange (send/receive): from element to ghost
  public :: PAR_GHOST_NODE_EXCHANGE            ! Ghost nodal exchange (send/receive)
  public :: PAR_FROM_GHOST_ELEMENT_EXCHANGE    ! Ghost element exchange (send/receive): from ghost to element
  public :: PAR_FROM_GHOST_BOUNDARY_EXCHANGE   ! Ghost element exchange (send/receive): from ghost to boundary
  public :: PAR_FROM_GHOST_NODE_EXCHANGE       ! Ghost node exchange (send/receive): from ghost to node
  public :: PAR_GHOST_BOUNDARY_EXCHANGE        ! Ghost boundary exchange (send/receive)
  public :: PAR_DEFINE_COMMUNICATOR            ! Define the communicator according to some keywords
  public :: PAR_COMM_SPLIT                     ! Split a communicator
  public :: PAR_SEND                           ! Send arrays to a specific partition
  public :: PAR_RECEIVE                        ! Receive arrays to a specific partition
  public :: PAR_SEND_RECEIVE                   ! Send and receive arrays to a specific partition
  public :: PAR_SEND_RECEIVE_IP                ! Send and receive arrays to a specific partition
  public :: PAR_SEND_RECEIVE_RP                ! Send and receive arrays to a specific partition
  public :: PAR_SEND_RECEIVE_TO_ALL            ! Send and receive arrays to all
  public :: PAR_EXCHANGE                       ! Bridge to broadcast data from master to slave
  public :: PAR_COMM_RANK_AND_SIZE             ! Give rank (and size) of a communicator
  public :: PAR_GATHER                         ! Gather of integers, reals, characters
  public :: PAR_GATHERV                        ! Gather of integers, reals, characters
  public :: PAR_SCATTER                        ! Scatter of integers, reals, characters
  public :: PAR_SCATTERV                       ! Scatter of integers, reals, characters
  public :: PAR_ALLGATHERV                     ! All Gatherv
  public :: PAR_ALLGATHER                      ! All Gather
  public :: PAR_INIT                           ! Initialize MPI
  public :: PAR_BROADCAST                      ! Broadcast
  public :: PAR_LENGTH_INTEGER                 ! Length of integers
  public :: PAR_IMASTER_IN_COMMUNICATOR        ! If I am a master of a given communicator
  public :: PAR_BARRIER                        ! Barrier
  public :: PAR_INITIALIZE_NON_BLOCKING_COMM   ! Initialize non-blocking communications variables
  public :: PAR_START_NON_BLOCKING_COMM        ! Start non-blocking communications
  public :: PAR_END_NON_BLOCKING_COMM          ! End non-blocking communications
  public :: PAR_SET_NON_BLOCKING_COMM_NUMBER   ! Set the non-blocking communicator number
  public :: PAR_WAITALL                        ! Waitall
  public :: PAR_ALL_TO_ALL_ARRAY_OPERATION     ! Array operations between all partitions of communicators
  public :: PAR_POINT_TO_POINT_ARRAY_OPERATION ! Array operations between all partitions of communicators
  public :: PAR_INTERFACE_NODE_EXCHANGE_VALUE  ! Interface nodal exchange only send receive, used in GPU directs.
  public :: COMMUNICATOR_VALUE                 ! getting communicator pointers and size
  public :: PAR_COMM_SPLIT_SIMVIZ              ! communicator splitting for viz ans im processes
  public :: PAR_ALLTOALLV                      ! All to all V
  public :: PAR_MPI_ERROR_TO_MESSAGE           ! MPI error number to message

contains

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-12-29
  !> @brief   Define communcator and communication arrays
  !> @details Define the communicator according to a keyword
  !>
  !-----------------------------------------------------------------------

  subroutine PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    implicit none

    character(*),        optional,          intent(in)  :: wherein
    integer(4),                             intent(out) :: PAR_COMM_TO_USE
    type(comm_data_par), optional, pointer, intent(inout) :: commu
    integer(ip)                                         :: icolo,jcolo

    !if( IPARALL ) then
    if( present(wherein) ) then
       if( trim(wherein) == 'IN THE UNIVERSE' ) then
          !
          ! In the universe
          !
#ifndef MPI_OFF
          PAR_COMM_TO_USE = MPI_COMM_WORLD
#endif
       else if( trim(wherein) == 'IN THE WORLD' ) then
          !
          ! In the world
          !
          PAR_COMM_TO_USE = int(PAR_COMM_WORLD,4_ip)   ! Alya world
          if( present(commu) .and. associated(PAR_COMM_MY_CODE_ARRAY) )   commu => commd

       else if( trim(wherein) == 'IN MY CODE' ) then
          !
          ! In my code
          !
          !icolo = par_code_zone_subd_to_color(current_code,0_ip,0_ip)
          PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,4_ip)
          !PAR_COMM_COLOR(icolo,icolo)
          if( present(commu) .and. associated(PAR_COMM_MY_CODE_ARRAY) ) commu => PAR_COMM_MY_CODE_ARRAY(1)

       else if( trim(wherein) == 'IN MY CODE WITHOUT MASTER' ) then
          !
          ! In my code without master
          !
          PAR_COMM_TO_USE = PAR_COMM_MY_CODE_WM4
          if( present(commu) .and. associated(PAR_COMM_MY_CODE_ARRAY) ) commu => PAR_COMM_MY_CODE_ARRAY(1)

       else if( trim(wherein) == 'IN MY ZONE' .or. trim(wherein) == 'IN CURRENT ZONE' ) then
          !
          ! In my current zone
          !
          icolo = par_code_zone_subd_to_color(current_code,current_zone,0_ip)
          PAR_COMM_TO_USE = int(PAR_COMM_COLOR(icolo,icolo),4)
          if( present(commu) .and. associated(PAR_COMM_MY_CODE_ARRAY) ) commu => PAR_COMM_COLOR_ARRAY(icolo)

       else if( trim(wherein) == 'IN MY SUBD' ) then
          !
          ! In my current subd
          !
          icolo = par_code_zone_subd_to_color(current_code,0_ip,current_subd)
          PAR_COMM_TO_USE = int(PAR_COMM_COLOR(icolo,icolo),4)
          if( present(commu) .and. associated(PAR_COMM_MY_CODE_ARRAY) ) commu => PAR_COMM_COLOR_ARRAY(icolo)

       else if( trim(wherein) == 'IN CURRENT COUPLING' ) then
          !
          ! In my current coupling
          !
          icolo = color_target
          jcolo = color_source
          ! REVISAR IN CURRENT COUPLING
          PAR_COMM_TO_USE = int(PAR_COMM_COLOR(icolo,jcolo),4)
          ! print*,"DEBUG: color_target, color_source ", color_target, color_source, PAR_COMM_TO_USE
          if( present(commu) ) then
             call runend('PAR_DEFINE_COMMUNICATOR: WRONG OPTION 1')
          end if

       else if( trim(wherein) == 'IN CURRENT COLOR' .or. trim(wherein) == 'IN CURRENT TARGET COLOR' ) then
          !
          ! In my current target color
          !
          icolo = color_target
          PAR_COMM_TO_USE = int(PAR_COMM_COLOR(icolo,icolo),4)
          if( present(commu) .and. associated(PAR_COMM_MY_CODE_ARRAY) ) commu => PAR_COMM_COLOR_ARRAY(icolo)

       else if( trim(wherein) == 'IN CURRENT SOURCE COLOR' ) then
          !
          ! In my current source color
          !
          icolo = color_source
          PAR_COMM_TO_USE = int(PAR_COMM_COLOR(icolo,icolo),4)
          if( present(commu) .and. associated(PAR_COMM_MY_CODE_ARRAY) ) commu => PAR_COMM_COLOR_ARRAY(icolo)

       else if( trim(wherein) == 'IN CURRENT' ) then
          !
          ! Uses current communicator
          !
          PAR_COMM_TO_USE = int(PAR_COMM_CURRENT,4)

          if( present(commu) ) then
             call runend('PAR_DEFINE_COMMUNICATOR: WRONG OPTION 2')
          end if
       else if( trim(wherein) == 'IN SFC PARTITION' ) then

          PAR_COMM_TO_USE = int(PAR_COMM_SFC_WM,4)
       else if( trim(wherein) == 'IN SFC PARTITION WITH MASTER' ) then

          PAR_COMM_TO_USE = int(PAR_COMM_SFC,4)
       else if( trim(wherein) == 'IN MPIO' ) then

          PAR_COMM_TO_USE = int(PAR_COMM_MPIO_WM,4)
       else if( trim(wherein) == 'IN MPIO WITH MASTER' ) then

          PAR_COMM_TO_USE = int(PAR_COMM_MPIO,4)

       else

          call runend('PAR DEFINE COMMUNICATOR: INVALID COMMUNICATOR OPTION')

       end if
    else
       PAR_COMM_TO_USE =  int(PAR_COMM_WORLD,4)
       commu            => commd
    end if
    !else
    !   PAR_COMM_TO_USE =  0
    !   commu            => commd
    !end if

  end subroutine PAR_DEFINE_COMMUNICATOR

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_INTERFACE_NODE_EXCHANGE_IP
  !
  !----------------------------------------------------------------------

  subroutine PAR_INTERFACE_NODE_EXCHANGE_IP_0(n,xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: n
    integer(ip),                    intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),         optional, intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE .or. npoin == 0 ) return
    ndofn = n
    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
    end if
    call PAR_INTERFACE_NODE_EXCHANGE_IP(ON_NODES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_NODE_EXCHANGE_IP_0

  subroutine PAR_INTERFACE_NODE_EXCHANGE_IP_1(xx,what,wherein,wsynch,dom_k,COMM)
    implicit none
    integer(ip),          pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),         optional, intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    type(comm_data_par),  optional, intent(in)    :: COMM
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    integer(4)                                    :: PAR_COMM_TO_USE

   if( INOTSLAVE .or. npoin == 0 ) return
    ndofn = 1
    if( size(xx,1) /= npoin .and. size(xx,1) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE IP_1')
    if( present(COMM) ) then
       PAR_COMM_TO_USE = int(COMM % PAR_COMM_WORLD,4)
       call PAR_INTERFACE_NODE_EXCHANGE_IP(ON_NODES,ndofn,xx,what,COMM,PAR_COMM_TO_USE,wsynch,dom_k)
    else
       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
       end if
       call PAR_INTERFACE_NODE_EXCHANGE_IP(ON_NODES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_INTERFACE_NODE_EXCHANGE_IP_1

  subroutine PAR_INTERFACE_NODE_EXCHANGE_IP_2(xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),         optional, intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE .or. npoin == 0 ) return
    ndofn = size(xx,1)
    if( size(xx,2) /= npoin .and. size(xx,2) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE IP_2')
    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
    end if
    call PAR_INTERFACE_NODE_EXCHANGE_IP(ON_NODES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_NODE_EXCHANGE_IP_2

  subroutine PAR_INTERFACE_NODE_EXCHANGE_IP_3(xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),          pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),         optional, intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE .or. npoin == 0 ) return
    if( size(xx,3) <= 2 ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= npoin .and. size(xx,2) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE IP_3')
    else
       ndofn = size(xx,1)*size(xx,2)
       if( size(xx,3) /= npoin .and. size(xx,3) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE IP_3')
    end if
    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
    end if
    call PAR_INTERFACE_NODE_EXCHANGE_IP(ON_NODES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_NODE_EXCHANGE_IP_3

  subroutine PAR_INTERFACE_NODE_EXCHANGE_IP_2b(xx,what,commu,wsynch,dom_k)
    implicit none
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),  pointer,  intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE .or. npoin == 0 ) return
    ndofn = size(xx,1)
    if( size(xx,2) /= npoin .and. size(xx,2) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
    PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
    call PAR_INTERFACE_NODE_EXCHANGE_IP(ON_NODES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_NODE_EXCHANGE_IP_2b

  !----------------------------------------------------------------------
  !
  ! PAR_INTERFACE_NODE_EXCHANGE_IP: NODE ASSEMBLY FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_INTERFACE_NODE_EXCHANGE_IP(onwhat,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: onwhat
    integer(ip),                    intent(in)    :: ndofn
    integer(ip),                    intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    integer(4),                     intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,nsize,jj,dom_i,kdofn
    integer(ip)                                   :: ipoin,ini,kk,idofn,jdofn
    integer(4)                                    :: istat4,nsize4,count4
    integer(4)                                    :: dom_i4,my_rank4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    integer(4),                        pointer    :: status4(:,:)
    integer(ip)                                   :: dom_j
    integer(ip)                                   :: bound_dim
    integer(ip),                       pointer    :: bound_perm(:)
    integer(ip),                       pointer    :: bound_size(:)

#ifndef MPI_OFF
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! On nodes or edges
       !
       if( onwhat == ON_NODES ) then
          bound_dim  =  commu % bound_dim
          bound_perm => commu % bound_perm
          bound_size => commu % bound_size
       else
          bound_dim  =  commu % bedge_dim
          bound_perm => commu % bedge_perm
          bound_size => commu % bedge_size
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_isend(bound_dim * ndofn))
             allocate(tmp_irecv(bound_dim * ndofn))
             !
             ! Save in temp_send
             !
             if( trim(what) == 'DISTRIBUTE' ) then
                do jj = 1,bound_dim
                   ipoin = bound_perm(jj)
                   if( ipoin < commu % npoi2 .or. ipoin > commu % npoi3 ) xx(1:ndofn,ipoin) = 0_ip
                end do
             end if

             kk = 0
             do jj = 1,bound_dim
                ipoin = bound_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_isend(kk) = xx(idofn,ipoin)
                   tmp_irecv(kk) = 0
                end do
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             istat4 = 0_4
             kk = 0
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   ini   = ndofn * ( bound_size(ii)   - 1 ) + 1
                   nsize = ndofn * ( bound_size(ii+1) - 1 ) + 1 - ini

                   nsize4 = int(nsize,4)
                   if( asynch ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           tmp_isend(ini:ini+nsize-1), nsize4, &
                           PAR_INTEGER,  dom_i4, 0_4,          &
                           PAR_COMM_TO_USE, ireq4(kk), istat4 )
                      kk = kk + 1
                      call MPI_Irecv(&
                           tmp_irecv(ini:ini+nsize-1), nsize4, &
                           PAR_INTEGER,  dom_i4, 0_4,          &
                           PAR_COMM_TO_USE, ireq4(kk), istat4 )
                   else
                      call MPI_Sendrecv(                       &
                           tmp_isend(ini:), nsize4,            &
                           PAR_INTEGER, dom_i4, 0_4,           &
                           tmp_irecv(ini:), nsize4,            &
                           PAR_INTEGER, dom_i4, 0_4,           &
                           PAR_COMM_TO_USE, status, istat4    )
                   end if
                   if( istat4 /= 0_4 ) call runend('PAR_INTERFACE_NODE_EXCHANGE_IP: MPI ERROR')
                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE,2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if( trim(what) == 'SUM' .or. trim(what) == 'ASSEMBLY' .or. trim(what) == 'DISTRIBUTE' ) then
                !
                ! SUM
                !
                kk = 0
                do jj = 1,bound_dim
                   ipoin = bound_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = xx(idofn,ipoin) + tmp_irecv(kk)
                   end do
                end do

             else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' ) then
                !
                ! MAX
                !
                kk = 0
                do jj = 1,bound_dim
                   ipoin = bound_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = max(xx(idofn,ipoin),tmp_irecv(kk))
                   end do
                end do

             else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' ) then
                !
                ! MIN
                !
                kk = 0
                do jj = 1,bound_dim
                   ipoin = bound_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = min(xx(idofn,ipoin),tmp_irecv(kk))
                   end do
                end do

             else if( trim(what) == 'TAKE MIN' ) then
                !
                ! TAKE MIN
                !
                call runend('TAKE MIN NOT CODED')
                call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank4)
                kk = 0
                do ii = 1,commu % nneig
                   dom_i = commu % neights(ii)
                   if( my_rank4 > dom_i ) then
                      ini   = ndofn * ( bound_size(ii)   - 1 ) + 1
                      nsize = ndofn * ( bound_size(ii+1) - 1 ) + 1 - ini
                      do jj = bound_size(ii),bound_size(ii+1)-1
                         ipoin = bound_perm(jj)
                         do idofn = 1,ndofn
                            kk = kk + 1
                            xx(idofn,ipoin) = tmp_irecv(kk)
                         end do
                      end do
                   end if
                end do

             else if( trim(what) == 'MIN RANK OR NEGATIVE' ) then
                !
                ! MIN RANK OR NEGATIVE
                !
                call runend('MIN RANK OR NEGATIVE NOT CODED')
                call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank4)
                kk = 0
                do ii = 1,commu % nneig
                   dom_i = commu % neights(ii)
                   ini   = ndofn * ( bound_size(ii)   - 1 ) + 1
                   nsize = ndofn * ( bound_size(ii+1) - 1 ) + 1 - ini
                   do jj = bound_size(ii),bound_size(ii+1)-1
                      ipoin = bound_perm(jj)
                      do idofn = 1,ndofn
                         kk = kk + 1
                         if( xx(idofn,ipoin) > 0 .and. tmp_irecv(kk) > 0 ) then
                            if( my_rank4 > dom_i ) then
                               !pard1 = 1
                               xx(idofn,ipoin) = -abs(xx(idofn,ipoin))
                            end if
                         end if
                      end do
                   end do
                end do

             else if( trim(what) == 'DIF' .or. trim(what) == 'DIFFERENCE' .or. trim(what) == 'DIFF' ) then
                !
                ! DIF
                !
                kk = 0
                do jj = 1,bound_dim
                   ipoin = bound_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      if( abs(xx(idofn,ipoin)-tmp_irecv(kk)) > 0 ) then
                         if( ipoin <= npoi1 .or. ( ipoin >= commu % npoi2 .and. ipoin <= commu % npoi3 ) ) write(6,*) 'DIF=',xx(idofn,ipoin)-tmp_irecv(kk)
                         xx(idofn,ipoin) = huge(1_ip)
                      end if
                   end do
                end do

             else if( trim(what) == 'MERGE' ) then
                !
                ! MERGE
                !
                kk = 0

                do jj = 1,bound_dim
                   ipoin = bound_perm(jj)

                   loop_kdofn1: do kdofn = 1,ndofn
                      if( xx(kdofn,ipoin) == 0 ) exit loop_kdofn1
                   end do loop_kdofn1
                   kdofn = min(kdofn,ndofn)

                   do idofn = 1,ndofn
                      kk = kk + 1
                      jdofn = 0
                      do while( jdofn < ndofn )
                         jdofn = jdofn + 1
                         if( tmp_irecv(kk) == xx(jdofn,ipoin) ) jdofn = 2*ndofn
                      end do
                      if( jdofn /= 2*ndofn ) then
                         xx(kdofn,ipoin) = tmp_irecv(kk)
                      end if
                   end do
                end do

                do jj = 1,bound_dim
                   ipoin = bound_perm(jj)
                   loop_kdofn2: do kdofn = 1,ndofn
                      if( xx(kdofn,ipoin) == 0 ) exit loop_kdofn2
                   end do loop_kdofn2
                   kdofn = min(kdofn,ndofn)
                   if( kdofn > 0 ) call maths_heap_sort(1_ip,kdofn,xx(1:kdofn,ipoin))
                end do

             else
                call runend('UNKNOWN ORDER')
             end if

             ipass = 0
             deallocate(tmp_irecv)
             deallocate(tmp_isend)

          end if

       end if

    end if
#endif

  end subroutine PAR_INTERFACE_NODE_EXCHANGE_IP

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_INTERFACE_NODE_EXCHANGE_LG
  !
  !----------------------------------------------------------------------

  subroutine PAR_INTERFACE_NODE_EXCHANGE_LG_1(xx,what,wherein,wsynch,dom_k)
    implicit none
    logical(lg),          pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),         optional, intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    integer(4)                                    :: PAR_COMM_TO_USE

   if( INOTSLAVE ) return
    ndofn = 1
    if( size(xx,1) /= npoin .and. size(xx,1) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
    end if
    call PAR_INTERFACE_NODE_EXCHANGE_LG(ON_NODES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_NODE_EXCHANGE_LG_1

  !----------------------------------------------------------------------
  !
  ! PAR_INTERFACE_NODE_EXCHANGE_IP: NODE ASSEMBLY FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_INTERFACE_NODE_EXCHANGE_LG(onwhat,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: onwhat
    integer(ip),                    intent(in)    :: ndofn
    logical(lg),                    intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    integer(4),                     intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,nsize,jj,dom_i
    integer(ip)                                   :: ipoin,ini,kk,idofn
    integer(4)                                    :: istat4,nsize4,count4
    integer(4)                                    :: dom_i4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    integer(4),                        pointer    :: status4(:,:)
    integer(ip)                                   :: dom_j
    integer(ip)                                   :: bound_dim
    integer(ip),                       pointer    :: bound_perm(:)
    integer(ip),                       pointer    :: bound_size(:)

#ifndef MPI_OFF
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! On nodes or edges
       !
       if( onwhat == ON_NODES ) then
          bound_dim  =  commu % bound_dim
          bound_perm => commu % bound_perm
          bound_size => commu % bound_size
       else
          bound_dim  =  commu % bedge_dim
          bound_perm => commu % bedge_perm
          bound_size => commu % bedge_size
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_lsend(bound_dim * ndofn))
             allocate(tmp_lrecv(bound_dim * ndofn))
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,bound_dim
                ipoin = bound_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_lsend(kk) = xx(idofn,ipoin)
                   tmp_lrecv(kk) = .false.
                end do
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             istat4 = 0_4
             kk = 0
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   ini   = ndofn * ( bound_size(ii)   - 1 ) + 1
                   nsize = ndofn * ( bound_size(ii+1) - 1 ) + 1 - ini

                   nsize4 = int(nsize,4)
                   if( asynch ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           tmp_lsend(ini:ini+nsize-1), nsize4, &
                           MPI_LOGICAL,  dom_i4, 0_4,          &
                           PAR_COMM_TO_USE, ireq4(kk), istat4 )
                      kk = kk + 1
                      call MPI_Irecv(&
                           tmp_lrecv(ini:ini+nsize-1), nsize4, &
                           MPI_LOGICAL,  dom_i4, 0_4,          &
                           PAR_COMM_TO_USE, ireq4(kk), istat4 )
                   else
                      call MPI_Sendrecv(                       &
                           tmp_lsend(ini:), nsize4,            &
                           MPI_LOGICAL, dom_i4, 0_4,           &
                           tmp_lrecv(ini:), nsize4,            &
                           PAR_INTEGER, dom_i4, 0_4,           &
                           PAR_COMM_TO_USE, status, istat4    )
                   end if
                   if( istat4 /= 0_4 ) call runend('PAR_INTERFACE_NODE_EXCHANGE_IP: MPI ERROR')
                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE,2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if( trim(what) == 'OR' ) then
                !
                ! OR
                !
                kk = 0
                do jj = 1,bound_dim
                   ipoin = bound_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = xx(idofn,ipoin) .or. tmp_lrecv(kk)
                   end do
                end do

             else if( trim(what) == 'AND' ) then
                !
                ! AND
                !
                kk = 0
                do jj = 1,bound_dim
                   ipoin = bound_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = xx(idofn,ipoin) .and. tmp_lrecv(kk)
                   end do
                end do

             else
                call runend('UNKNOWN ORDER')
             end if

             ipass = 0
             deallocate(tmp_lrecv)
             deallocate(tmp_lsend)

          end if

       end if

    end if
#endif

  end subroutine PAR_INTERFACE_NODE_EXCHANGE_LG

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_INTERFACE_NODE_EXCHANGE_RP
  !
  !----------------------------------------------------------------------

  subroutine PAR_INTERFACE_NODE_EXCHANGE_RP_00(n,xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(*)
    character(*),                   intent(in)    :: what
    character(*),         optional, intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
    end if
    call PAR_INTERFACE_NODE_EXCHANGE_RP(ON_NODES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_NODE_EXCHANGE_RP_00

  subroutine PAR_INTERFACE_NODE_EXCHANGE_RP_0(n,xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),         optional, intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
    end if
    call PAR_INTERFACE_NODE_EXCHANGE_RP(ON_NODES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_NODE_EXCHANGE_RP_0

  subroutine PAR_INTERFACE_NODE_EXCHANGE_RP_1(xx,what,wherein,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),         optional, intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE .or. npoin == 0 ) return
    ndofn = 1
    if( what(len(what):len(what)) /= 'I' ) then
       if( size(xx,1) /= npoin .and. size(xx,1) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE RP_1')
    end if
    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
    end if
    call PAR_INTERFACE_NODE_EXCHANGE_RP(ON_NODES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_NODE_EXCHANGE_RP_1

  subroutine PAR_INTERFACE_NODE_EXCHANGE_RP_2(xx,what,wherein,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),         optional, intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE .or. npoin == 0 ) return
    ndofn = size(xx,1)
    if( what(len(what):len(what)) /= 'I' ) then
       if( size(xx,2) /= npoin .and. size(xx,2) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE RP_2')
    end if
    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
    end if
    call PAR_INTERFACE_NODE_EXCHANGE_RP(ON_NODES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_NODE_EXCHANGE_RP_2

  subroutine PAR_INTERFACE_NODE_EXCHANGE_RP_3(xx,what,wherein,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),         optional, intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE .or. npoin == 0 ) return
    if( what(len(what):len(what)) /= 'I' ) then
       if( size(xx,3) <= 2 ) then
          ndofn = size(xx,1)
          if( size(xx,2) /= npoin .and. size(xx,2) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE IP_3')
       else
          ndofn = size(xx,1)*size(xx,2)
          if( size(xx,3) /= npoin .and. size(xx,3) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE RP_3')
       end if
    end if
    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE,commu)
    end if
    call PAR_INTERFACE_NODE_EXCHANGE_RP(ON_NODES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_NODE_EXCHANGE_RP_3

  subroutine PAR_INTERFACE_NODE_EXCHANGE_RP_2b(xx,what,commu,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE .or. npoin == 0 ) return
    ndofn = size(xx,1)
    if( what(len(what):len(what)) /= 'I' ) then
       if( size(xx,2) /= npoin .and. size(xx,2) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
    end if
    PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
    call PAR_INTERFACE_NODE_EXCHANGE_RP(ON_NODES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_NODE_EXCHANGE_RP_2b

  !----------------------------------------------------------------------
  !
  ! PAR_INTERFACE_NODE_EXCHANGE_RP: NODE ASSEMBLY FOR INTEGERS
  ! If last letter of what is I (meaning I interface) then it is assumed
  ! that only the interface part of the vector is passed, from npoi1+1
  ! to npoin.
  !
  !              Interior nodes                    Interface nodes
  ! <----------------------------------------><----------------------->
  ! +-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
  ! |  1  |     |     |     |     |     |NPOI1|     |     |     |NPOIN|
  ! +-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
  !
  !----------------------------------------------------------------------

  subroutine PAR_INTERFACE_NODE_EXCHANGE_RP(onwhat,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: onwhat
    integer(ip),                    intent(in)    :: ndofn
    real(rp),                       intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    integer(4),                     intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,nsize,jj,dom_i,dom_min
    integer(ip)                                   :: ipoin,ini,kk,idofn
    integer(ip)                                   :: offset,ipoin_offset
    integer(4)                                    :: istat4,nsize4,count4
    integer(4)                                    :: dom_i4,my_rank4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    integer(4),                        pointer    :: status4(:,:)
    integer(ip)                                   :: dom_j
    integer(ip)                                   :: bound_dim
    integer(ip),                       pointer    :: bound_perm(:)
    integer(ip),                       pointer    :: bound_size(:)
    real(rp),                          pointer    :: bound_diff(:,:)

#ifndef MPI_OFF
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! On nodes or edges
       !
       if( onwhat == ON_NODES ) then
          bound_dim  =  commu % bound_dim
          bound_perm => commu % bound_perm
          bound_size => commu % bound_size
       else
          bound_dim  =  commu % bedge_dim
          bound_perm => commu % bedge_perm
          bound_size => commu % bedge_size
       end if
       !
       ! Offset
       !
       if( what(len(what):len(what)) /= 'I' ) then
          offset = 0
       else
          offset = -npoi1
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_rsend(bound_dim * ndofn))
             allocate(tmp_rrecv(bound_dim * ndofn))
             !
             ! Distribute values among neighbors
             !
             if( trim(what) == 'DISTRIBUTE' ) then
                do jj = 1,bound_dim
                   ipoin = bound_perm(jj)
                   if( ipoin < commu % npoi2 .or. ipoin > commu % npoi3 ) xx(1:ndofn,ipoin) = 0.0_rp
                end do
             end if
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,bound_dim
                ipoin = bound_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_rsend(kk) = xx(idofn,ipoin+offset)
                   tmp_rrecv(kk) = 0.0_rp
                end do
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             istat4 = 0_4
             kk = 0
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   ini   = ndofn * ( bound_size(ii)   - 1 ) + 1
                   nsize = ndofn * ( bound_size(ii+1) - 1 ) + 1 - ini

                   nsize4 = int(nsize,4)
                   if( asynch ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           tmp_rsend(ini:ini+nsize-1), nsize4, &
                           MPI_DOUBLE_PRECISION,  dom_i4, 0_4, &
                           PAR_COMM_TO_USE, ireq4(kk), istat4 )
                      kk = kk + 1
                      call MPI_Irecv(&
                           tmp_rrecv(ini:ini+nsize-1), nsize4, &
                           MPI_DOUBLE_PRECISION,  dom_i4, 0_4, &
                           PAR_COMM_TO_USE, ireq4(kk), istat4 )
                   else
                      call MPI_Sendrecv(                       &
                           tmp_rsend(ini:), nsize4,            &
                           MPI_DOUBLE_PRECISION, dom_i4, 0_4,  &
                           tmp_rrecv(ini:), nsize4, MPI_DOUBLE_PRECISION, dom_i4, 0_4,&
                           & PAR_COMM_TO_USE, status, istat4    )
                   end if
                   if( istat4 /= 0_4 ) call runend('PAR_INTERFACE_NODE_EXCHANGE_RP: MPI ERROR')

                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE,2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if(    trim(what) == 'SUM'  .or. trim(what) == 'ASSEMBLY' .or. trim(what) == 'SUMI'&
                  & .or. trim(what) == 'DISTRIBUTE' ) then
                !
                ! SUM
                !
                kk = 0
                do jj = 1,bound_dim
                   ipoin = bound_perm(jj)
                   ipoin_offset = ipoin + offset
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin_offset) = xx(idofn,ipoin_offset) + tmp_rrecv(kk)
                   end do
                end do

             else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' .or. trim(what) == 'MAXI' )&
                  & then
                !
                ! MAX
                !
                kk = 0
                do jj = 1,bound_dim
                   ipoin = bound_perm(jj)
                   ipoin_offset = ipoin + offset
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin_offset) = max(xx(idofn,ipoin_offset),tmp_rrecv(kk))
                   end do
                end do

             else if( trim(what) == 'TAKE MIN' ) then
                !
                ! TAKE MIN
                !
                call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank4)
                kk = 0
                dom_min = min(int(my_rank4,ip),minval(commu % neights))
                do ii = 1,commu % nneig
                   dom_i = commu % neights(ii)
                   if( dom_i == dom_min ) then
                      ini   = ndofn * ( bound_size(ii)   - 1 ) + 1
                      nsize = ndofn * ( bound_size(ii+1) - 1 ) + 1 - ini
                      do jj = bound_size(ii),bound_size(ii+1)-1
                         ipoin = bound_perm(jj)
                         do idofn = 1,ndofn
                            kk = kk + 1
                            xx(idofn,ipoin) = tmp_rrecv(kk)
                         end do
                      end do
                   end if
                end do

             else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' .or. trim(what) == 'MINI' ) then
                !
                ! MIN
                !
                kk = 0
                do jj = 1,bound_dim
                   ipoin = bound_perm(jj)
                   ipoin_offset = ipoin + offset
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin_offset) = min(xx(idofn,ipoin_offset),tmp_rrecv(kk))
                   end do
                end do

             else if( trim(what) == 'DIF' .or. trim(what) == 'DIFFERENCE' .or. trim(what) == 'DIFF' ) then
                !
                ! DIF
                !
                kk = 0
                do jj = 1,bound_dim
                   ipoin = bound_perm(jj)
                   ipoin_offset = ipoin + offset
                   do idofn = 1,ndofn
                      kk = kk + 1
                      if( abs(xx(idofn,ipoin_offset)-tmp_rrecv(kk)) / (abs(xx(idofn,ipoin_offset))+epsilon(1.0_rp)) > 1.0e-12_rp ) then
                         if( ipoin <= npoi1 .or. ( ipoin >= commu % npoi2 .and. ipoin <= commu % npoi3 ) ) write(6,*) 'DIF=',xx(idofn,ipoin_offset)-tmp_rrecv(kk)
                         xx(idofn,ipoin_offset) = huge(1.0_rp)
                      end if
                   end do
                end do

             else if( trim(what) == 'MAX_DIF' .or. trim(what) == 'MAX_DIFFERENCE' .or. trim(what) == 'MAX_DIFF' ) then
                !
                ! MAX_DIF
                !
                nullify(bound_diff)
                call memory_alloca(par_memor,'bound_diff','par_send_receive_to_all_rp_0',bound_diff,ndofn,npoin)
                bound_diff(1:ndofn,1:npoin) = xx(1:ndofn,1:npoin)
                xx(1:ndofn,1:npoin) = 0.0_rp

                kk = 0
                do jj = 1,bound_dim
                   ipoin = bound_perm(jj)
                   ipoin_offset = ipoin + offset
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin_offset) = max( abs(bound_diff(idofn,ipoin)-tmp_rrecv(kk)),xx(idofn,ipoin_offset))
                   end do
                end do
                call memory_deallo(par_memor,'bound_diff','par_send_receive_to_all_rp_0',bound_diff)

             else
                call runend('UNKNOWN ORDER')
             end if

             ipass = 0
             deallocate(tmp_rrecv)
             deallocate(tmp_rsend)

          end if

       end if

    end if
#endif

  end subroutine PAR_INTERFACE_NODE_EXCHANGE_RP

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_INTERFACE_OWN_NODE_EXCHANGE
  !
  ! KPOIN enables to limit the positions in the receiving arrays.
  ! It can be used for example not to write on halos, for example if XX
  ! is only defined until own+boundary nodes
  !
  !----------------------------------------------------------------------

  subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_IP_0(ndofn,xx,message,kpoin)
    integer(ip),                    intent(in)    :: ndofn
    integer(ip),                    intent(inout) :: xx(ndofn,*)
    character(*),         optional, intent(in)    :: message
    integer(ip),          optional, intent(in)    :: kpoin
    call PAR_INTERFACE_OWN_NODE_EXCHANGE_IP(ndofn,xx,message,kpoin)
  end subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_IP_0

  subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_IP_1(xx,message,kpoin)
    integer(ip),          pointer,  intent(inout) :: xx(:)
    character(*),         optional, intent(in)    :: message
    integer(ip),          optional, intent(in)    :: kpoin
    integer(ip)                                   :: ndofn
    if( .not. associated(xx) ) return
    ndofn = 1
    call PAR_INTERFACE_OWN_NODE_EXCHANGE_IP(ndofn,xx,message,kpoin)
  end subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_IP_1

  subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_IP_2(xx,message,kpoin)
    implicit none
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),         optional, intent(in)    :: message
    integer(ip),          optional, intent(in)    :: kpoin
    integer(ip)                                   :: ndofn
    if( .not. associated(xx) ) return
    ndofn = size(xx,1)
    call PAR_INTERFACE_OWN_NODE_EXCHANGE_IP(ndofn,xx,message,kpoin)
  end subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_IP_2

  subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_RP_00(ndofn,xx,message,kpoin)
    integer(ip),                    intent(in)    :: ndofn
    real(rp),                       intent(inout) :: xx(*)
    character(*),         optional, intent(in)    :: message
    integer(ip),          optional, intent(in)    :: kpoin
    call PAR_INTERFACE_OWN_NODE_EXCHANGE_RP(ndofn,xx,message,kpoin)
  end subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_RP_00

  subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_RP_0(ndofn,xx,message,kpoin)
    integer(ip),                    intent(in)    :: ndofn
    real(rp),                       intent(inout) :: xx(ndofn,*)
    character(*),         optional, intent(in)    :: message
    integer(ip),          optional, intent(in)    :: kpoin
    call PAR_INTERFACE_OWN_NODE_EXCHANGE_RP(ndofn,xx,message,kpoin)
  end subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_RP_0

  subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_RP_1(xx,message,kpoin)
    real(rp),             pointer,  intent(inout) :: xx(:)
    character(*),         optional, intent(in)    :: message
    integer(ip),          optional, intent(in)    :: kpoin
    integer(ip)                                   :: ndofn
    if( .not. associated(xx) ) return
    ndofn = 1
    call PAR_INTERFACE_OWN_NODE_EXCHANGE_RP(ndofn,xx,message,kpoin)
  end subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_RP_1

  subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_RP_2(xx,message,kpoin)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),         optional, intent(in)    :: message
    integer(ip),          optional, intent(in)    :: kpoin
    integer(ip)                                   :: ndofn
    if( .not. associated(xx) ) return
    ndofn = size(xx,1)
    call PAR_INTERFACE_OWN_NODE_EXCHANGE_RP(ndofn,xx,message,kpoin)
  end subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_RP_2

  subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_RP(ndofn,xx,message,kpoin)
    implicit none
    integer(ip),                    intent(in)    :: ndofn
    real(rp),                       intent(inout) :: xx(ndofn,*)
    character(*), optional,         intent(in)    :: message
    integer(ip),  optional,         intent(in)    :: kpoin
    integer(ip)                                   :: ii,nsize,jj,dom_i
    integer(ip)                                   :: ipoin,kk,idofn
    integer(ip)                                   :: ii_recv_size,ii_send_size
    integer(ip)                                   :: ii_recv,ii_send,ineig
    integer(ip)                                   :: size_non_blocking
    real(rp),                          target     :: tmp_rsend_ok(2)
    real(rp),                          target     :: tmp_rrecv_ok(2)
    logical(lg)                                   :: do_send_receive
    logical(lg)                                   :: do_wait_and_assemble

    do_send_receive      = .false.
    do_wait_and_assemble = .false.

    if( present(message) ) then
       if( trim(message) == 'SEND RECEIVE'              ) do_send_receive      = .true.
       if( trim(message) == 'WAIT AND ASSEMBLE'         ) do_wait_and_assemble = .true.
       if( trim(message) == 'SEND RECEIVE AND ASSEMBLE' ) then
          do_send_receive      = .true.
          do_wait_and_assemble = .true.
       end if
    else
       do_send_receive      = .true.
       do_wait_and_assemble = .true.
    end if

    if( do_send_receive ) then
       !
       ! Fill-in sending array
       !
       allocate(tmp_rsend(max(commd % full_row_send_dim * ndofn,1_ip)))
       allocate(tmp_rrecv(max(commd % full_row_recv_dim * ndofn,1_ip)))

       do ii = 1,commd % full_row_send_dim
          ipoin = commd % full_row_send_perm(ii)
          do idofn = 1,ndofn
             !if( (ii-1)*ndofn+idofn > size(tmp_rsend) ) call runend('TROUBLE 0')
             tmp_rsend((ii-1)*ndofn+idofn) = xx(idofn,ipoin)
          end do
       end do
       !
       ! Maximum number of non-blocking operations
       !
       size_non_blocking = commd % full_row_send_nneig + commd % full_row_recv_nneig
       !
       ! Send and receive
       !
       call PAR_START_NON_BLOCKING_COMM(1_ip,size_non_blocking)     ! Set number of requests
       call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)                  ! Set non-blocking communication number to 1

       do ineig = 1,commd % full_row_send_nneig
          dom_i        =   commd % full_row_send_neights(ineig)
          ii_send_size = ( commd % full_row_send_size(ineig+1)-commd % full_row_send_size(ineig) ) * ndofn
          ii_send      = ( commd % full_row_send_size(ineig)-1)*ndofn+1
          ii_recv_size = 0
          !if(ii_send < 1 .or. ii_send > size(tmp_rsend) )   call runend('TROUBLE 1')
          !if(ii_send + ii_send_size - 1 > size(tmp_rsend) ) call runend('TROUBLE 2')
          call PAR_SEND_RECEIVE(ii_send_size,ii_recv_size,tmp_rsend(ii_send:),tmp_rrecv_ok,'IN MY CODE',dom_i,'ASYNCHRONOUS')
       end do
       do ineig = 1,commd % full_row_recv_nneig
          dom_i        =   commd % full_row_recv_neights(ineig)
          ii_recv_size = ( commd % full_row_recv_size(ineig+1)-commd % full_row_recv_size(ineig) ) * ndofn
          ii_recv      = ( commd % full_row_recv_size(ineig)-1)*ndofn+1
          ii_send_size = 0
          !if(ii_recv < 1 .or. ii_recv > size(tmp_rrecv) )   call runend('TROUBLE 3')
          !if(ii_recv + ii_recv_size - 1 > size(tmp_rrecv) ) call runend('TROUBLE 4')
          call PAR_SEND_RECEIVE(ii_send_size,ii_recv_size,tmp_rsend_ok,tmp_rrecv(ii_recv:),'IN MY CODE',dom_i,'ASYNCHRONOUS')
       end do

    end if

    if( do_wait_and_assemble ) then
       !
       ! Wait and assemble
       !
       call PAR_END_NON_BLOCKING_COMM(1_ip)

       if( present(kpoin) ) then
          do ii = 1,commd % full_row_recv_dim
             ipoin = commd % full_row_recv_perm(ii)
             if( ipoin <= kpoin ) then
                do idofn = 1,ndofn
                   xx(idofn,ipoin) = tmp_rrecv((ii-1)*ndofn+idofn)
                end do
             end if
          end do
       else
          do ii = 1,commd % full_row_recv_dim
             ipoin = commd % full_row_recv_perm(ii)
             do idofn = 1,ndofn
                !if( (ii-1)*ndofn+idofn > size(tmp_rrecv) ) call runend('TROUBLE 5')
                xx(idofn,ipoin) = tmp_rrecv((ii-1)*ndofn+idofn)
             end do
          end do
       end if

       deallocate(tmp_rsend)
       deallocate(tmp_rrecv)

    end if

  end subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_RP

  subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_IP(ndofn,xx,message,kpoin)
    implicit none
    integer(ip),                    intent(in)    :: ndofn
    integer(ip),                    intent(inout) :: xx(ndofn,*)
    character(*), optional,         intent(in)    :: message
    integer(ip),  optional,         intent(in)    :: kpoin
    integer(ip)                                   :: ii,nsize,jj,dom_i
    integer(ip)                                   :: ipoin,kk,idofn
    integer(ip)                                   :: ii_recv_size,ii_send_size
    integer(ip)                                   :: ii_recv,ii_send,ineig
    integer(ip)                                   :: size_non_blocking
    integer(ip),                       target     :: tmp_isend_ok(2)
    integer(ip),                       target     :: tmp_irecv_ok(2)
    logical(lg)                                   :: do_send_receive
    logical(lg)                                   :: do_wait_and_assemble

    do_send_receive      = .false.
    do_wait_and_assemble = .false.

    if( present(message) ) then
       if( trim(message) == 'SEND RECEIVE'      ) do_send_receive      = .true.
       if( trim(message) == 'WAIT AND ASSEMBLE' ) do_wait_and_assemble = .true.
    else
       do_send_receive      = .true.
       do_wait_and_assemble = .true.
    end if

    if( do_send_receive ) then
       !
       ! Send receive
       !
       allocate(tmp_isend(max(commd % full_row_send_dim * ndofn,1_ip)))
       allocate(tmp_irecv(max(commd % full_row_recv_dim * ndofn,1_ip)))

       do ii = 1,commd % full_row_send_dim
          ipoin = commd % full_row_send_perm(ii)
          do idofn = 1,ndofn
             tmp_isend((ii-1)*ndofn+idofn) = xx(idofn,ipoin)
          end do
       end do
       !
       ! Maximum number of non-blocking operations
       !
       size_non_blocking = commd % full_row_send_nneig + commd % full_row_recv_nneig
       !
       ! Send and receive
       !
       call PAR_START_NON_BLOCKING_COMM(1_ip,size_non_blocking)     ! Set number of requests
       call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)                  ! Set non-blocking communication number to 1

       do ineig = 1,commd % full_row_send_nneig
          dom_i        =   commd % full_row_send_neights(ineig)
          ii_send_size = ( commd % full_row_send_size(ineig+1)-commd % full_row_send_size(ineig) ) * ndofn
          ii_send      = ( commd % full_row_send_size(ineig)-1)*ndofn+1
          ii_recv_size = 0
          call PAR_SEND_RECEIVE(ii_send_size,ii_recv_size,tmp_isend(ii_send:),tmp_irecv_ok,'IN MY CODE',dom_i,'ASYNCHRONOUS')
       end do
       do ineig = 1,commd % full_row_recv_nneig
          dom_i        =   commd % full_row_recv_neights(ineig)
          ii_recv_size = ( commd % full_row_recv_size(ineig+1)-commd % full_row_recv_size(ineig) ) * ndofn
          ii_recv      = ( commd % full_row_recv_size(ineig)-1)*ndofn+1
          ii_send_size = 0
          call PAR_SEND_RECEIVE(ii_send_size,ii_recv_size,tmp_isend_ok,tmp_irecv(ii_recv:),'IN MY CODE',dom_i,'ASYNCHRONOUS')
       end do

    end if

    if( do_wait_and_assemble ) then
       !
       ! Wait and assemble
       !
       call PAR_END_NON_BLOCKING_COMM(1_ip)

       if( present(kpoin) ) then
          do ii = 1,commd % full_row_recv_dim
             ipoin = commd % full_row_recv_perm(ii)
             if( ipoin <= kpoin ) then
                do idofn = 1,ndofn
                   xx(idofn,ipoin) = tmp_irecv((ii-1)*ndofn+idofn)
                end do
             end if
          end do
       else
          do ii = 1,commd % full_row_recv_dim
             ipoin = commd % full_row_recv_perm(ii)
             do idofn = 1,ndofn
                xx(idofn,ipoin) = tmp_irecv((ii-1)*ndofn+idofn)
             end do
          end do
       end if

       deallocate(tmp_isend)
       deallocate(tmp_irecv)

    end if

  end subroutine PAR_INTERFACE_OWN_NODE_EXCHANGE_IP

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_INTERFACE_NODE_EXCHANGE_RP
  !
  !----------------------------------------------------------------------

  subroutine PAR_COUPLING_NODE_EXCHANGE_RP_00(n,xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_COUPLING_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_COUPLING_NODE_EXCHANGE_RP_00

  subroutine PAR_COUPLING_NODE_EXCHANGE_RP_0(n,xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_COUPLING_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_COUPLING_NODE_EXCHANGE_RP_0

  subroutine PAR_COUPLING_NODE_EXCHANGE_RP_1(xx,what,wherein,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = 1
    if( size(xx,1) /= npoin .and. size(xx,1) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE RP_1')
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_COUPLING_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_COUPLING_NODE_EXCHANGE_RP_1

  subroutine PAR_COUPLING_NODE_EXCHANGE_RP_1b(xx,what,commu,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    ndofn = 1
    if( size(xx,1) /= npoin .and. size(xx,1) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE RP_1b')
    PAR_COMM_TO_USE = int(commu % PAR_COMM_WORLD,4)
    call PAR_COUPLING_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_COUPLING_NODE_EXCHANGE_RP_1b

  subroutine PAR_COUPLING_NODE_EXCHANGE_RP_2(xx,what,wherein,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = size(xx,1)
    if( size(xx,2) /= npoin .and. size(xx,2) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE RP_2')
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_COUPLING_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_COUPLING_NODE_EXCHANGE_RP_2

  subroutine PAR_COUPLING_NODE_EXCHANGE_RP_3(xx,what,wherein,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    if( size(xx,3) <= 2 ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= npoin .and. size(xx,2) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE IP_3')
    else
       ndofn = size(xx,1)*size(xx,2)
       if( size(xx,3) /= npoin .and. size(xx,3) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE RP_3')
    end if
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_COUPLING_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_COUPLING_NODE_EXCHANGE_RP_3

  subroutine PAR_COUPLING_NODE_EXCHANGE_RP_2b(xx,what,commu,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    ndofn = size(xx,1)
    if( size(xx,2) /= npoin .and. size(xx,2) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
    PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
    call PAR_COUPLING_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_COUPLING_NODE_EXCHANGE_RP_2b

  !----------------------------------------------------------------------
  !
  ! PAR_COUPLING_NODE_EXCHANGE_RP: NODE ASSEMBLY FOR REALS
  !
  ! From the receiver point of view
  ! +---+
  ! | 1 |
  ! +---+                                                    +---+  /|
  ! | 2 |                                                    |   |  ||
  ! +---+   lscat_dim                                        +---+  ||
  ! | 3 |                                                    |   |  ||
  ! +---+      /\  +---+    +---+---+---+---+---+---+---+    +---+  || lrecv_size(2)-lrecv_size(1)
  ! | 4 |      ||  |   |    |   |   |   |   |   |   |   |    |   |  ||
  ! +---+      ||  +---+    +---+---+---+---+---+---+---+    +---+  ||
  ! | 5 |      ||  |   |  = |   |   |   |   |   |   |   |    |   |  ||
  ! +---+      ||  +---+    +---+---+---+---+---+---+---+    +---+  \/
  ! | 6 |      ||  |   |    |   |   |   |   |   |   |   |    +---+  /|
  ! +---+      \/  +---+    +---+---+---+---+---+---+---+    |   |  ||
  ! | 7 |            |             matrix_ia(:)              +---+  ||
  ! +---+  <---------+             matrix_ja(:)              |   |  || lrecv_size(3)-lrecv_size(2)
  ! | 8 |    lscat_perm(:)         matrix_aa(:)              +---+  ||
  ! +---+                          matrix_nzdom              |   |  ||
  ! | 9 |                                                    +---+  \/
  ! +---+
  !
  ! Nodal array
  !
  !----------------------------------------------------------------------

  subroutine PAR_COUPLING_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

    implicit none
    integer(ip),                    intent(in)    :: ndofn
    real(rp),                       intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    integer(4),                     intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,nsize_recv,jj,dom_i,ll
    integer(ip)                                   :: ipoin,ini_recv,kk,idofn
    integer(4)                                    :: istat4,nsize4_send,nsize4_recv,count4
    integer(4)                                    :: ini_send,nsize_send
    integer(4)                                    :: dom_i4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    integer(4),                        pointer    :: status4(:,:)
    integer(ip)                                   :: dom_j

#ifndef MPI_OFF
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_rsend(commu % lsend_dim * ndofn))
             allocate(tmp_rrecv(commu % lrecv_dim * ndofn))
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,commu % lsend_dim
                ipoin = commu % lsend_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_rsend(kk) = xx(idofn,ipoin)
                end do
             end do
             do kk = 1,commu % lrecv_dim * ndofn
                tmp_rrecv(kk) = 0.0_rp
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             istat4 = 0_4
             kk = 0
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   ini_send   = ndofn * ( commu % lsend_size(ii)   - 1 ) + 1
                   nsize_send = ndofn * ( commu % lsend_size(ii+1) - 1 ) + 1 - ini_send
                   ini_recv   = ndofn * ( commu % lrecv_size(ii)   - 1 ) + 1
                   nsize_recv = ndofn * ( commu % lrecv_size(ii+1) - 1 ) + 1 - ini_recv

                   nsize4_send = int(nsize_send,4)
                   nsize4_recv = int(nsize_recv,4)

                   if( asynch ) then
                      !if( nsize_send > 0 ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           tmp_rsend(ini_send:ini_send+nsize_send-1), nsize4_send, &
                           MPI_DOUBLE_PRECISION,  dom_i4, 0_4, &
                           PAR_COMM_TO_USE, ireq4(kk), istat4 )
                      !end if
                      !if( nsize_recv > 0 ) then
                      kk = kk + 1
                      call MPI_Irecv(&
                           tmp_rrecv(ini_recv:ini_recv+nsize_recv-1), nsize4_recv, &
                           MPI_DOUBLE_PRECISION,  dom_i4, 0_4, &
                           PAR_COMM_TO_USE, ireq4(kk), istat4 )
                      !end if
                   else
                      call MPI_Sendrecv(                       &
                           tmp_rsend(ini_send:), nsize4_send,  &
                           MPI_DOUBLE_PRECISION, dom_i4, 0_4,  &
                           tmp_rrecv(ini_recv:), nsize4_recv,  &
                           MPI_DOUBLE_PRECISION, dom_i4, 0_4,  &
                           PAR_COMM_TO_USE, status, istat4    )
                   end if
                   if( istat4 /= 0_4 ) call runend('PAR_COUPLING_NODE_EXCHANGE_RP: MPI ERROR')

                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE,2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if( trim(what) == 'MATRIX' ) then
                !
                ! MATRIX
                !
                kk = 0
                do ii = 1,commu % lscat_dim
                   ipoin = commu % lscat_perm(ii)
                   do kk = commu % matrix_ia(ii),commu % matrix_ia(ii+1)-1
                      jj = commu % matrix_ja(kk)
                      xx(1:ndofn,ipoin) = xx(1:ndofn,ipoin) + commu % matrix_aa(kk) * tmp_rrecv( jj:jj )
                   end do
                end do

             else
                call runend('UNKNOWN ORDER')
             end if

             ipass = 0
             deallocate(tmp_rrecv)
             deallocate(tmp_rsend)

          end if

       end if

    end if
#endif

  end subroutine PAR_COUPLING_NODE_EXCHANGE_RP

  !----------------------------------------------------------------------
  !
  ! Norm
  !
  !----------------------------------------------------------------------

  !subroutine PAR_NORM_s(xx,xnorm,wherein)
  !  implicit none
  !  real(rp),     intent(in)  :: xx
  !  real(rp),     intent(out) :: xnorm
  !  character(*),             :: wherein
  !end subroutine PAR_NORM_s

  !----------------------------------------------------------------------
  !
  ! SUM FOR REALS
  !
  !----------------------------------------------------------------------

  subroutine PAR_SUM_RP_s(xx,wherein,commu)
    implicit none
    real(rp),                      intent(inout) :: xx
    character(*),        optional                :: wherein
    type(comm_data_par), optional, intent(in)    :: commu
    integer(ip)                                  :: nsize
    integer(4)                                   :: istat4,nsize4
    integer(4)                                   :: PAR_COMM_TO_USE
    real(rp)                                     :: yy
    character(30)                                :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       my_wherein = 'IN MY CODE'
       if( present(commu) ) then
          PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
       else if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       if( PAR_COMM_TO_USE /= -1 ) then
          nsize  = 1
          nsize4 = int(nsize,4)
          if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
             yy = 0.0_rp
          else
             yy = xx
          end if
          call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_PRECISION,&
               MPI_SUM,PAR_COMM_TO_USE,istat4)
          if( istat4 /= 0_4 ) call runend('PAR_SUM_RP_s: MPI ERROR')
       end if
    end if
#endif

  end subroutine PAR_SUM_RP_s

  subroutine PAR_SUM_RP_0(n,xx,wherein,wsynch,PAR_COMM_IN4)
    implicit none
    integer(ip),  intent(in)    :: n
    real(rp),     intent(inout) :: xx(n)
    character(*), optional      :: wherein
    character(*), optional      :: wsynch
    integer(4),   optional      :: PAR_COMM_IN4
    integer(ip)                 :: ii,nsize
    integer(4)                  :: istat4,nsize4
    integer(4)                  :: PAR_COMM_TO_USE
    logical(lg)                 :: asynch
    character(30)               :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL .and. n > 0 ) then
       my_wherein = 'IN MY CODE'
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE=PAR_COMM_IN4
       else
          if( present(wherein) ) then
             my_wherein = trim(wherein)
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
          end if
       endif
       asynch = .false.
       if( present(wsynch) ) then
          if( trim(wsynch) == 'NON BLOCKING' ) asynch = .true.
       end if
       nsize  = n
       nsize4 = int(nsize,4)
       allocate(yy_non_blocking(n))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do ii = 1,nsize
             yy_non_blocking(ii) = 0.0_rp
          end do
       else
          do ii = 1,nsize
             yy_non_blocking(ii) = xx(ii)
          end do
       end if
       if( asynch ) then
#ifdef MPI3
          call MPI_IAllReduce(yy_non_blocking,xx,nsize4,MPI_DOUBLE_PRECISION,&
               MPI_SUM,PAR_COMM_TO_USE,ireq41(1),istat4)
#else
          call runend('THUS FUNCTIONALITY REQUIRES MPI3: COMPILE ALYA WITH -DMPI3')
#endif
       else
          call MPI_AllReduce(yy_non_blocking,xx,nsize4,MPI_DOUBLE_PRECISION,&
               MPI_SUM,PAR_COMM_TO_USE,istat4)
          deallocate(yy_non_blocking)
       end if
       if( istat4 /= 0_4 ) call runend('PAR_SUM_RP_0: MPI ERROR')
    end if
#endif

  end subroutine PAR_SUM_RP_0

  subroutine PAR_SUM_RP_0b(xx,PAR_COMM_IN4)
    implicit none
    real(rp),     intent(inout) :: xx
    integer(4),   intent(in)    :: PAR_COMM_IN4
    real(rp)                    :: aa(1)

    aa(1)=xx
    call PAR_SUM_RP_0(1_ip,aa,PAR_COMM_IN4=PAR_COMM_IN4)
    xx=aa(1)

  end subroutine PAR_SUM_RP_0b

  subroutine PAR_SUM_RP_02(n1,n2,xx,wherein)
    implicit none
    integer(ip),  intent(in)    :: n1
    integer(ip),  intent(in)    :: n2
    real(rp),     intent(inout) :: xx(n1,n2)
    character(*), optional      :: wherein
    integer(ip)                 :: ii,jj,nsize
    integer(4)                  :: istat4,nsize4
    integer(4)                  :: PAR_COMM_TO_USE
    real(rp),     allocatable   :: yy(:,:)
    character(30)               :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    nsize  = n1*n2
    if( IPARALL .and. nsize > 0 ) then
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize4 = int(nsize,4)
       allocate(yy(n1,n2))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do jj = 1,n2
             do ii = 1,n1
                yy(ii,jj) = 0.0_rp
             end do
          end do
       else
          do jj = 1,n2
             do ii = 1,n1
                yy(ii,jj) = xx(ii,jj)
             end do
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_PRECISION,&
            MPI_SUM,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_SUM_RP_0: MPI ERROR')
       deallocate(yy)
    end if
#endif

  end subroutine PAR_SUM_RP_02

  subroutine PAR_SUM_RP_1(xx,wherein,wsynch,PAR_COMM_IN4)
    implicit none
    real(rp),     pointer, intent(inout) :: xx(:)
    character(*), optional,intent(in)    :: wherein
    character(*), optional,intent(in)    :: wsynch
    integer(4),   optional, intent(in)   :: PAR_COMM_IN4
    integer(ip)                          :: ii,nsize
    integer(4)                           :: istat4,nsize4
    integer(4)                           :: PAR_COMM_TO_USE
    logical(lg)                          :: asynch
    character(30)                        :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       my_wherein = 'IN MY CODE'
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             my_wherein = trim(wherein)
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
          end if
       endif
       asynch = .false.
       if( present(wsynch) ) then
          if( trim(wsynch) == 'NON BLOCKING' ) asynch = .true.
       end if
       nsize  = size(xx,1)
       nsize4 = int(nsize,4)
       allocate( yy_non_blocking(lbound(xx,1):ubound(xx,1)) )
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do ii = lbound(yy_non_blocking,1),ubound(yy_non_blocking,1)
             yy_non_blocking(ii) = 0.0_rp
          end do
       else
          do ii = lbound(yy_non_blocking,1),ubound(yy_non_blocking,1)
             yy_non_blocking(ii) = xx(ii)
          end do
       end if
      if( asynch ) then
#ifdef MPI3
          call MPI_IAllReduce(yy_non_blocking,xx,nsize4,MPI_DOUBLE_PRECISION,&
               MPI_SUM,PAR_COMM_TO_USE,ireq41(1),istat4)
#else
          call runend('THUS FUNCTIONALITY REQUIRES MPI3: COMPILE ALYA WITH -DMPI3')
#endif
       else
          call MPI_AllReduce(yy_non_blocking,xx,nsize4,MPI_DOUBLE_PRECISION,&
               MPI_SUM,PAR_COMM_TO_USE,istat4)
          deallocate( yy_non_blocking )
       end if
       if( istat4 /= 0_4 ) call runend('PAR_SUM_RP_1: MPI ERROR')
    end if
#endif

  end subroutine PAR_SUM_RP_1

  subroutine PAR_SUM_RP_2(xx,wherein)
    implicit none
    real(rp),     pointer, intent(inout)         :: xx(:,:)
    character(*),          optional              :: wherein
    integer(ip)                                  :: ii,jj,nsize1,nsize2,nsize
    integer(4)                                   :: istat4,nsize4
    integer(4)                                   :: PAR_COMM_TO_USE
    real(rp),     pointer                        :: yy(:,:)
    character(30)                                :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize1 = size(xx,1)
       nsize2 = size(xx,2)
       nsize  = nsize1*nsize2
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1),lbound(xx,2):ubound(xx,2)) )
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = 0.0_rp
             end do
          end do
       else
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = xx(ii,jj)
             end do
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_PRECISION,&
            MPI_SUM,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_SUM_RP_2: MPI ERROR')
       deallocate( yy )
    end if
#endif

  end subroutine PAR_SUM_RP_2

  subroutine PAR_SUM_RP_3(xx,wherein,who,PAR_COMM_IN4)
    implicit none
    real(rp),     pointer, intent(inout)         :: xx(:,:,:)
    character(*),          optional              :: wherein
    character(*), optional                       :: who
    integer(4),   optional, intent(in)           :: PAR_COMM_IN4
    integer(ip)                                  :: ii,jj
    integer(4)                                   :: istat4,nsize4
    integer(4)                                   :: PAR_COMM_TO_USE
    real(rp),     pointer                        :: yy(:,:,:)
    logical(lg)                                  :: not_include_master
    character(30)                                :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       my_wherein = 'IN MY CODE'
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             my_wherein = trim(wherein)
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
          end if
       end if
       not_include_master = .true.
       if( present(who) ) then
          if( trim(who) == 'INCLUDE MASTER' ) not_include_master = .false.
       end if
       nsize4 = int(size(xx),4)
       allocate( yy(lbound(xx,1):ubound(xx,1),lbound(xx,2):ubound(xx,2),lbound(xx,3):ubound(xx,3)) )
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. not_include_master ) then
          yy = 0.0_rp
       else
          yy = xx
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_PRECISION,&
            MPI_SUM,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_SUM_RP_2: MPI ERROR')
       deallocate( yy )
    end if
#endif

  end subroutine PAR_SUM_RP_3

  !----------------------------------------------------------------------
  !
  ! SUM FOR COMPLEX
  !
  !----------------------------------------------------------------------

  subroutine PAR_SUM_CX_s(xx,wherein)
    implicit none
    complex(rp),  intent(inout) :: xx
    character(*), optional      :: wherein
    integer(ip)                 :: nsize
    integer(4)                  :: istat4,nsize4
    integer(4)                  :: PAR_COMM_TO_USE
    complex(rp)                 :: yy
    character(30)               :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize  = 1
       nsize4 = int(nsize,4)
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          yy = 0.0_rp
       else
          yy = xx
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_COMPLEX,&
            MPI_SUM,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_SUM_CX_s: MPI ERROR')
    end if
#endif

  end subroutine PAR_SUM_CX_s

  subroutine PAR_SUM_CX_0(n,xx,wherein)
    implicit none
    integer(ip),  intent(in)    :: n
    complex(rp),  intent(inout) :: xx(n)
    character(*), optional      :: wherein
    integer(ip)                 :: ii,nsize
    integer(4)                  :: istat4,nsize4
    integer(4)                  :: PAR_COMM_TO_USE
    complex(rp),  allocatable   :: yy(:)
    character(30)               :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL .and. n > 0 ) then
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize  = n
       nsize4 = int(nsize,4)
       allocate(yy(n))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do ii = 1,nsize
             yy(ii) = 0.0_rp
          end do
       else
          do ii = 1,nsize
             yy(ii) = xx(ii)
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_COMPLEX,&
            MPI_SUM,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_SUM_CX_0: MPI ERROR')
       deallocate(yy)
    end if
#endif

  end subroutine PAR_SUM_CX_0

  subroutine PAR_SUM_CX_1(xx,wherein)
    implicit none
    complex(rp),  pointer, intent(inout) :: xx(:)
    character(*),          optional      :: wherein
    integer(ip)                          :: ii,nsize
    integer(4)                           :: istat4,nsize4
    integer(4)                           :: PAR_COMM_TO_USE
    complex(rp),  pointer                :: yy(:)
    character(30)                        :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize  = size(xx,1)
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1)) )
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) = 0.0_rp
          end do
       else
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) = xx(ii)
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_COMPLEX,&
            MPI_SUM,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_SUM_CX_1: MPI ERROR')
       deallocate( yy )
    end if
#endif

  end subroutine PAR_SUM_CX_1

  subroutine PAR_SUM_CX_2(xx,wherein)
    implicit none
    complex(rp),  pointer, intent(inout) :: xx(:,:)
    character(*), optional               :: wherein
    integer(ip)                          :: ii,jj,nsize1,nsize2,nsize
    integer(4)                           :: istat4,nsize4
    integer(4)                           :: PAR_COMM_TO_USE
    complex(rp),  pointer                :: yy(:,:)
    character(30)                        :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize1 = size(xx,1)
       nsize2 = size(xx,2)
       nsize  = nsize1*nsize2
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1),lbound(xx,2):ubound(xx,2)))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = (0.0_rp,0.0_rp)
             end do
          end do
       else
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = xx(ii,jj)
             end do
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_COMPLEX,&
            MPI_SUM,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_SUM_CX_2: MPI ERROR')
       deallocate( yy )
    end if
#endif

  end subroutine PAR_SUM_CX_2

  !----------------------------------------------------------------------
  !
  ! AVERAGE FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_AVERAGE_IP_s(xx,wherein)
    implicit none
    integer(ip),  intent(inout) :: xx
    character(*), optional      :: wherein
    integer(4)                  :: PAR_COMM_TO_USE
    integer(4)                  :: my_rank4,comm_size4

    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
    end if
    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank4,comm_size4)
    if( present(wherein) ) then
       call PAR_SUM(xx,wherein)
    else
       call PAR_SUM(xx)
    end if
    if( IPARALL ) xx = xx / (max(1_ip,comm_size4-1_ip))

  end subroutine PAR_AVERAGE_IP_s

  subroutine PAR_AVERAGE_IP_0(n,xx,wherein)
    implicit none
    integer(ip),  intent(in)    :: n
    integer(ip),  intent(inout) :: xx(n)
    character(*), optional      :: wherein
    integer(4)                  :: PAR_COMM_TO_USE
    integer(4)                  :: my_rank4,comm_size4

    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
    end if
    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank4,comm_size4)
    if( present(wherein) ) then
       call PAR_SUM(n,xx,wherein)
    else
       call PAR_SUM(n,xx)
    end if
    if( IPARALL ) xx(1:n) = xx(1:n) / (max(1_ip,comm_size4-1_ip))

  end subroutine PAR_AVERAGE_IP_0

  !----------------------------------------------------------------------
  !
  ! AVERAGE FOR REAL
  !
  !----------------------------------------------------------------------

  subroutine PAR_AVERAGE_RP_s(xx,wherein)
    implicit none
    real(rp),     intent(inout) :: xx
    character(*), optional      :: wherein
    integer(4)                  :: PAR_COMM_TO_USE
    integer(4)                  :: my_rank4,comm_size4

    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
    end if
    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank4,comm_size4)
    if( present(wherein) ) then
       call PAR_SUM(xx,wherein)
    else
       call PAR_SUM(xx)
    end if
    if( IPARALL ) xx = xx / real(max(1_ip,comm_size4-1_ip),rp)

  end subroutine PAR_AVERAGE_RP_s

  subroutine PAR_AVERAGE_RP_0(n,xx,wherein)
    implicit none
    integer(ip),  intent(in)    :: n
    real(rp),     intent(inout) :: xx(n)
    character(*), optional      :: wherein
    integer(4)                  :: PAR_COMM_TO_USE
    integer(4)                  :: my_rank4,comm_size4

    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
    end if
    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank4,comm_size4)
    if( present(wherein) ) then
       call PAR_SUM(n,xx,wherein)
    else
       call PAR_SUM(n,xx)
    end if
    if( IPARALL ) xx(1:n) = xx(1:n) / real(max(1_ip,comm_size4-1_ip),rp)

  end subroutine PAR_AVERAGE_RP_0

 !----------------------------------------------------------------------
  !
  ! AVERAGE FOR REAL
  !
  !----------------------------------------------------------------------

  subroutine PAR_LOAD_BALANCE_RP_s(xx,wherein)
    implicit none
    real(rp),     intent(inout) :: xx
    character(*), optional      :: wherein
    integer(4)                  :: PAR_COMM_TO_USE
    integer(4)                  :: my_rank4,comm_size4
    real(rp)                    :: xx_max

    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
    end if
    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank4,comm_size4)
    call PAR_SUM(xx,wherein)
    xx_max = xx
    if( present(wherein) ) then
       call PAR_MAX(xx_max,wherein)
    else
       call PAR_MAX(xx_max)
    end if
    if( IPARALL ) xx = xx / real(max(1_ip,comm_size4-1_ip),rp)
    xx = xx / (xx_max + epsilon(1.0_rp))

  end subroutine PAR_LOAD_BALANCE_RP_s

  subroutine PAR_LOAD_BALANCE_RP_0(n,xx,wherein)
    implicit none
    integer(ip),  intent(in)    :: n
    real(rp),     intent(inout) :: xx(n)
    character(*), optional      :: wherein
    integer(4)                  :: PAR_COMM_TO_USE
    integer(4)                  :: my_rank4,comm_size4
    real(rp)                    :: xx_max(n)

    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
    end if
    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank4,comm_size4)
    if( present(wherein) ) then
       call PAR_SUM(n,xx,wherein)
    else
       call PAR_SUM(n,xx)
    end if
    xx_max = xx
    if( IPARALL ) xx(1:n) = xx(1:n) / real(max(1_ip,comm_size4-1_ip),rp)
     xx(1:n) = xx(1:n) / (xx_max(1:n) + epsilon(1.0_rp))

  end subroutine PAR_LOAD_BALANCE_RP_0

  !----------------------------------------------------------------------
  !
  ! SUM FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_SUM_IP_s(xx,wherein,who,PAR_COMM_IN4)
    implicit none
    integer(ip),            intent(inout) :: xx
    character(*), optional                :: wherein
    character(*), optional                :: who
    integer(4),   optional, intent(in)    :: PAR_COMM_IN4
    integer(ip)                           :: nsize
    integer(4)                            :: istat4,nsize4
    integer(4)                            :: PAR_COMM_TO_USE
    integer(ip)                           :: yy
    logical(lg)                           :: not_include_master
    character(30)                         :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       my_wherein = 'IN MY CODE'
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             my_wherein = trim(wherein)
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
          end if
       end if
       not_include_master = .true.
       if( present(who) ) then
          if( trim(who) == 'INCLUDE MASTER' ) not_include_master = .false.
       end if
       nsize  = 1
       nsize4 = int(nsize,4)
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. not_include_master ) then
          yy = 0_ip
       else
          yy = xx
       end if
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_SUM,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_SUM_IP_s: MPI ERROR')
    end if
#endif

  end subroutine PAR_SUM_IP_s

  subroutine PAR_SUM_IP_0(n,xx,wherein,who,PAR_COMM_IN4)
    implicit none
    integer(ip),  intent(in)            :: n
    integer(ip),  intent(inout)         :: xx(n)
    character(*), intent(in),  optional :: wherein
    character(*), intent(in),  optional :: who
    integer(4),   intent(in),  optional :: PAR_COMM_IN4
    integer(ip)                         :: ii,nsize
    integer(4)                          :: istat4,nsize4
    integer(4)                          :: PAR_COMM_TO_USE
    integer(ip),  allocatable           :: yy(:)
    logical(lg)                         :: not_include_master
    character(30)                       :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL .and. n > 0 ) then
       my_wherein = 'IN MY CODE'
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             my_wherein = trim(wherein)
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
          end if
       end if
       not_include_master = .true.
       if( present(who) ) then
          if( trim(who) == 'INCLUDE MASTER' ) not_include_master = .false.
       end if
       nsize  = n
       nsize4 = int(nsize,4)
       allocate(yy(n))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. not_include_master ) then
          yy = 0_ip
       else
          yy = xx
       end if
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_SUM,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_SUM_IP_0: MPI ERROR')
       deallocate(yy)
    end if
#endif

  end subroutine PAR_SUM_IP_0

  subroutine PAR_SUM_IP_03(n1,n2,n3,xx,wherein,who,PAR_COMM_IN4)
    implicit none
    integer(ip),  intent(in)           :: n1,n2,n3
    integer(ip),  intent(inout)        :: xx(n1,n2,n3)
    character(*), optional             :: wherein
    character(*), optional             :: who
    integer(4),   optional, intent(in) :: PAR_COMM_IN4
    integer(ip)                        :: ii
    integer(4)                         :: istat4,nsize4
    integer(4)                         :: PAR_COMM_TO_USE
    integer(ip),  allocatable          :: yy(:,:,:)
    logical(lg)                        :: not_include_master
    character(30)                      :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL .and. n1*n2*n3 > 0 ) then
       my_wherein = 'IN MY CODE'
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             my_wherein = trim(wherein)
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
          end if
       end if
       not_include_master = .true.
       if( present(who) ) then
          if( trim(who) == 'INCLUDE MASTER' ) not_include_master = .false.
       end if
       nsize4 = int(n1*n2*n3,4)
       allocate(yy(n1,n2,n3))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. not_include_master ) then
          yy = 0_ip
       else
          yy = xx
       end if
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_SUM,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_SUM_IP_0: MPI ERROR')
       deallocate(yy)
    end if
#endif

  end subroutine PAR_SUM_IP_03

  subroutine PAR_SUM_IP_1(xx,wherein,who,PAR_COMM_IN4)
    implicit none
    integer(ip),     pointer,  intent(inout) :: xx(:)
    character(*),    optional                :: wherein
    character(*),    optional                :: who
    integer(4),      optional, intent(in)    :: PAR_COMM_IN4
    integer(ip)                              :: ii,nsize
    integer(4)                               :: istat4,nsize4
    integer(4)                               :: PAR_COMM_TO_USE
    integer(ip),     pointer                 :: yy(:)
    logical(lg)                              :: not_include_master
    character(30)                            :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       my_wherein = 'IN MY CODE'
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             my_wherein = trim(wherein)
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
          end if
       end if
       not_include_master = .true.
       if( present(who) ) then
          if( trim(who) == 'INCLUDE MASTER' ) not_include_master = .false.
       end if
       nsize  = size(xx,1)
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1)) )
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. not_include_master ) then
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) = 0_ip
          end do
       else
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) = xx(ii)
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_SUM,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_SUM_IP_1: MPI ERROR')
       deallocate( yy )
    end if
#endif

  end subroutine PAR_SUM_IP_1

  subroutine PAR_SUM_IP_2(xx,wherein,who)
    implicit none
    integer(ip),     pointer, intent(inout) :: xx(:,:)
    character(*),    optional               :: wherein
    character(*),    optional               :: who
    integer(ip)                             :: ii,jj,nsize1,nsize2,nsize
    integer(4)                              :: istat4,nsize4
    integer(4)                              :: PAR_COMM_TO_USE
    integer(ip),     pointer                :: yy(:,:)
    logical(lg)                             :: not_include_master
    character(30)                           :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if

       not_include_master = .true.
       if( present(who) ) then
          if( trim(who) == 'INCLUDE MASTER' ) not_include_master = .false.
       end if

       nsize1 = size(xx,1)
       nsize2 = size(xx,2)
       nsize  = nsize1*nsize2
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1),lbound(xx,2):ubound(xx,2)))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. not_include_master ) then
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = 0_ip
             end do
          end do
       else
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = xx(ii,jj)
             end do
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_SUM,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_SUM_IP_2: MPI ERROR')
       deallocate( yy )
    end if
#endif

  end subroutine PAR_SUM_IP_2

  subroutine PAR_SUM_IP_3(xx,wherein,who,PAR_COMM_IN4)
    implicit none
    integer(ip),    pointer, intent(inout) :: xx(:,:,:)
    character(*),   optional               :: wherein
    character(*),   optional               :: who
    integer(4),     optional, intent(in)   :: PAR_COMM_IN4
    integer(4)                             :: istat4,nsize4
    integer(4)                             :: PAR_COMM_TO_USE
    integer(ip),    pointer                :: yy(:,:,:)
    logical(lg)                            :: not_include_master
    character(30)                          :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       my_wherein = 'IN MY CODE'
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             my_wherein = trim(wherein)
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
          end if
       end if
       not_include_master = .true.
       if( present(who) ) then
          if( trim(who) == 'INCLUDE MASTER' ) not_include_master = .false.
       end if

       nsize4 = int(size(xx),4)
       allocate( yy(lbound(xx,1):ubound(xx,1),lbound(xx,2):ubound(xx,2),lbound(xx,3):ubound(xx,3)) )
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. not_include_master ) then
          yy = 0
       else
          yy = xx
       end if
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_SUM,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_SUM_IP_3: MPI ERROR')
       deallocate( yy )
    end if
#endif

  end subroutine PAR_SUM_IP_3

  !----------------------------------------------------------------------
  !
  ! MAX FOR REALS
  !
  !----------------------------------------------------------------------

  subroutine PAR_MAX_RP_s(xx,wherein,who,rank_max_owner,PAR_COMM_IN4)
    implicit none
    real(rp),     intent(inout) :: xx
    character(*), optional      :: wherein
    character(*), optional      :: who
    integer(ip),  optional      :: rank_max_owner
    integer(4),   optional, intent(in) :: PAR_COMM_IN4
    integer(ip)                 :: nsize
    integer(4)                  :: istat4,nsize4
    integer(4)                  :: PAR_COMM_TO_USE
    integer(4)                  :: rank_min
    integer(4)                  :: rank_max_owner4
    real(rp)                    :: yy
    logical(lg)                 :: not_include_master
    character(30)               :: my_wherein

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
        if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             my_wherein = trim(wherein)
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             my_wherein = 'IN MY CODE'
             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
          end if
       end if
       not_include_master = .true.
       if( present(who) ) then
          if( trim(who) == 'INCLUDE MASTER' ) not_include_master = .false.
       end if

       nsize  = 1
       nsize4 = int(nsize,4)
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. not_include_master ) then
          yy = -huge(1.0_rp)
       else
          yy = xx
       end if

       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_PRECISION,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)

       if( istat4 /= 0_4 ) call runend('PAR_MAX_RP_s: MPI ERROR')
       if( present(rank_max_owner) ) then
          rank_min = huge(1_4)
          if( abs(yy-xx) < epsilon(1.0_rp) ) call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,rank_min)
          call MPI_AllReduce(rank_min,rank_max_owner4,nsize4,MPI_INTEGER,&
               MPI_MIN,PAR_COMM_TO_USE,istat4)
          rank_max_owner = int(rank_max_owner4,ip)
       end if
#endif
    else
       if( present(rank_max_owner) ) then
          rank_max_owner = -1
       end if
    end if

  end subroutine PAR_MAX_RP_s

  subroutine PAR_MAX_RP_0(n,xx,wherein,who,PAR_COMM_IN4)
    implicit none
    integer(ip),  intent(in)           :: n
    real(rp),     intent(inout)        :: xx(n)
    character(*), optional, intent(in) :: wherein
    character(*), optional, intent(in) :: who
    integer(4),   optional, intent(in) :: PAR_COMM_IN4
    integer(ip)                        :: ii,nsize
    integer(4)                         :: istat4,nsize4
    integer(4)                         :: PAR_COMM_TO_USE
    real(rp),     allocatable          :: yy(:)
    logical(lg)                        :: not_include_master
    character(30)                      :: my_wherein

    if( IPARALL .and. n > 0 ) then
#ifndef MPI_OFF
       istat4 = 0_4
       my_wherein = 'IN MY CODE'
        if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             my_wherein = trim(wherein)
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
          end if
       end if
       not_include_master = .true.
       if( present(who) ) then
          if( trim(who) == 'INCLUDE MASTER' ) not_include_master = .false.
       end if

       nsize  = n
       nsize4 = int(nsize,4)
       allocate(yy(n))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. not_include_master ) then
          do ii = 1,nsize
             yy(ii) = -huge(1.0_rp)
          end do
       else
          do ii = 1,nsize
             yy(ii) = xx(ii)
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_PRECISION,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_MAX_RP_0: MPI ERROR')
       deallocate(yy)
#endif
    end if

  end subroutine PAR_MAX_RP_0

  subroutine PAR_MAX_RP_1(xx,wherein,who,PAR_COMM_IN4)
    implicit none
    real(rp),     pointer, intent(inout) :: xx(:)
    character(*), optional, intent(in)   :: wherein
    character(*), optional, intent(in)   :: who
    integer(4),   optional, intent(in)   :: PAR_COMM_IN4
    integer(ip)                          :: ii,nsize
    integer(4)                           :: istat4,nsize4
    integer(4)                           :: PAR_COMM_TO_USE
    real(rp),     pointer                :: yy(:)
    logical(lg)                          :: not_include_master
    character(30)                        :: my_wherein

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       my_wherein = 'IN MY CODE'
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
          end if
       end if
       not_include_master = .true.
       if( present(who) ) then
          if( trim(who) == 'INCLUDE MASTER' ) not_include_master = .false.
       end if

       nsize  = size(xx,1)
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1)) )
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) = -huge(1.0_rp)
          end do
       else
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) = xx(ii)
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_PRECISION,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_MAX_RP_1: MPI ERROR')
       deallocate( yy )
#endif
    end if

  end subroutine PAR_MAX_RP_1

  subroutine PAR_MAX_RP_2(xx,wherein)
    implicit none
    real(rp),     pointer, intent(inout) :: xx(:,:)
    character(*),          optional      :: wherein
    integer(ip)                          :: ii,jj,nsize1,nsize2,nsize
    integer(4)                           :: istat4,nsize4
    integer(4)                           :: PAR_COMM_TO_USE
    real(rp),     pointer                :: yy(:,:)
    character(30)                        :: my_wherein

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          my_wherein = 'IN MY CODE'
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize1 = size(xx,1)
       nsize2 = size(xx,2)
       nsize  = nsize1*nsize2
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1),lbound(xx,2):ubound(xx,2)))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = -huge(1.0_rp)
             end do
          end do
       else
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = xx(ii,jj)
             end do
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_PRECISION,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_MAX_RP_2: MPI ERROR')
       deallocate( yy )
#endif
    end if

  end subroutine PAR_MAX_RP_2

  subroutine PAR_MAX_RP_3(xx,wherein)
    implicit none
    real(rp),     pointer, intent(inout) :: xx(:,:,:)
    character(*),          optional      :: wherein
    integer(ip)                          :: ii,jj,kk,nsize1,nsize2,nsize3,nsize
    integer(4)                           :: istat4,nsize4
    integer(4)                           :: PAR_COMM_TO_USE
    real(rp),     pointer                :: yy(:,:,:)
    character(30)                        :: my_wherein

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          my_wherein = 'IN MY CODE'
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize1 = size(xx,1)
       nsize2 = size(xx,2)
       nsize3 = size(xx,3)
       nsize  = nsize1*nsize2*nsize3
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1),lbound(xx,2):ubound(xx,2),lbound(xx,3):ubound(xx,3)))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do kk = lbound(yy,3),ubound(yy,3)
             do jj = lbound(yy,2),ubound(yy,2)
                do ii = lbound(yy,1),ubound(yy,1)
                   yy(ii,jj,kk) = -huge(1.0_rp)
                end do
             end do
          end do
       else
          do kk = lbound(yy,3),ubound(yy,3)
             do jj = lbound(yy,2),ubound(yy,2)
                do ii = lbound(yy,1),ubound(yy,1)
                   yy(ii,jj,kk) = xx(ii,jj,kk)
                end do
             end do
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_PRECISION,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_MAX_RP_2: MPI ERROR')
       deallocate( yy )
#endif
    end if

  end subroutine PAR_MAX_RP_3

  !----------------------------------------------------------------------
  !
  ! MAX FOR COMPLEX
  !
  !----------------------------------------------------------------------

  subroutine PAR_MAX_CX_s(xx,wherein)
    implicit none
    complex(rp),  intent(inout) :: xx
    character(*), optional      :: wherein
    integer(ip)                 :: nsize
    integer(4)                  :: istat4,nsize4
    integer(4)                  :: PAR_COMM_TO_USE
    complex(rp)                 :: yy
    character(30)               :: my_wherein

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          my_wherein = 'IN MY CODE'
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize  = 1
       nsize4 = int(nsize,4)
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          yy = -huge(1.0_rp)
       else
          yy = xx
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_COMPLEX,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_MAX_CX_s: MPI ERROR')
#endif
    end if

  end subroutine PAR_MAX_CX_s

  subroutine PAR_MAX_CX_0(n,xx,wherein)
    implicit none
    integer(ip),  intent(in)    :: n
    complex(rp),  intent(inout) :: xx(n)
    character(*), optional      :: wherein
    integer(ip)                 :: ii,nsize
    integer(4)                  :: istat4,nsize4
    integer(4)                  :: PAR_COMM_TO_USE
    complex(rp),  allocatable   :: yy(:)
    character(30)               :: my_wherein

    if( IPARALL .and. n > 0 ) then
#ifndef MPI_OFF
       istat4 = 0_4
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          my_wherein = 'IN MY CODE'
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize  = n
       nsize4 = int(nsize,4)
       allocate(yy(n))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do ii = 1,nsize
             yy(ii) = -huge(1.0_rp)
          end do
       else
          do ii = 1,nsize
             yy(ii) = xx(ii)
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_COMPLEX,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_MAX_CX_0: MPI ERROR')
       deallocate(yy)
#endif
    end if

  end subroutine PAR_MAX_CX_0

  subroutine PAR_MAX_CX_1(xx,wherein)
    implicit none
    complex(rp),  pointer, intent(inout) :: xx(:)
    character(*),          optional      :: wherein
    integer(ip)                          :: ii,nsize
    integer(4)                           :: istat4,nsize4
    integer(4)                           :: PAR_COMM_TO_USE
    complex(rp),  pointer                :: yy(:)
    character(30)                        :: my_wherein

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          my_wherein = 'IN MY CODE'
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize  = size(xx,1)
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1)))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) = -huge(1.0_rp)
          end do
       else
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) = xx(ii)
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_COMPLEX,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_MAX_CX_1: MPI ERROR')
       deallocate( yy )
#endif
    end if

  end subroutine PAR_MAX_CX_1

  subroutine PAR_MAX_CX_2(xx,wherein)
    implicit none
    complex(rp),  pointer, intent(inout) :: xx(:,:)
    character(*), optional               :: wherein
    integer(ip)                          :: ii,jj,nsize1,nsize2,nsize
    integer(4)                           :: istat4,nsize4
    integer(4)                           :: PAR_COMM_TO_USE
    complex(rp),  pointer                :: yy(:,:)
    character(30)                        :: my_wherein

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          my_wherein = 'IN MY CODE'
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize1 = size(xx,1)
       nsize2 = size(xx,2)
       nsize  = nsize1*nsize2
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1),lbound(xx,2):ubound(xx,2)))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = -huge(1.0_rp)
             end do
          end do
       else
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = xx(ii,jj)
             end do
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_COMPLEX,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_MAX_CX_2: MPI ERROR')
       deallocate( yy )
#endif
    end if

  end subroutine PAR_MAX_CX_2

  !----------------------------------------------------------------------
  !
  ! MAX FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_MAX_4_s(xx,wherein,rank_max_owner,PAR_COMM_IN4,INCLUDE_ROOT)
    implicit none
    integer(4),   intent(inout)           :: xx
    character(*), intent(in),    optional :: wherein
    integer(4),   intent(inout), optional :: rank_max_owner
    integer(4),   intent(in),    optional :: PAR_COMM_IN4
    logical(lg),  intent(in),    optional :: INCLUDE_ROOT
    integer(ip)                           :: nsize
    integer(4)                            :: istat4,nsize4
    integer(4)                            :: PAR_COMM_TO_USE
    integer(4)                            :: rank_min
    integer(4)                            :: rank_max_owner4
    integer(4)                            :: yy
    logical(lg)                           :: if_include_root
    character(30)                         :: my_wherein

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       my_wherein = 'IN MY CODE'
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       if( present(INCLUDE_ROOT) ) then
          if_include_root = INCLUDE_ROOT
       else
          if_include_root = .false.
       end if
       nsize  = 1
       nsize4 = int(nsize,4)
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. .not. if_include_root ) then
          yy = -huge(1_4)
       else
          yy = xx
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_INTEGER,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_MAX_4_s: MPI ERROR')

       if( present(rank_max_owner) ) then
          rank_min = huge(1_4)
          if( xx == yy ) call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,rank_min)
          call MPI_AllReduce(rank_min,rank_max_owner4,nsize4,MPI_INTEGER,&
               MPI_MIN,PAR_COMM_TO_USE,istat4)
          rank_max_owner = int(rank_max_owner4,ip)
       end if

#endif
    else
       if( present(rank_max_owner) ) then
          rank_max_owner = -1
       end if
    end if

  end subroutine PAR_MAX_4_s

  subroutine PAR_MAX_8_s(xx,wherein,rank_max_owner,PAR_COMM_IN4,INCLUDE_ROOT)
    implicit none
    integer(8),   intent(inout)           :: xx
    character(*), intent(in),    optional :: wherein
    integer(8),   intent(inout), optional :: rank_max_owner
    integer(4),   intent(in),    optional :: PAR_COMM_IN4
    logical(lg),  intent(in),    optional :: INCLUDE_ROOT
    integer(ip)                           :: nsize
    integer(4)                            :: istat4,nsize4
    integer(4)                            :: PAR_COMM_TO_USE
    integer(4)                            :: rank_min
    integer(4)                            :: rank_max_owner4
    integer(8)                            :: yy
    logical(lg)                           :: if_include_root
    character(30)                         :: my_wherein

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       my_wherein = 'IN MY CODE'
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       if( present(INCLUDE_ROOT) ) then
          if_include_root = INCLUDE_ROOT
       else
          if_include_root = .false.
       end if
       nsize  = 1
       nsize4 = int(nsize,4)
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. .not. if_include_root ) then
          yy = -huge(1_8)
       else
          yy = xx
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_INTEGER8,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_MAX_IP_s: MPI ERROR')
       if( present(rank_max_owner) ) then
          rank_min = huge(1_4)
          if( yy == xx ) call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,rank_min)
          call MPI_AllReduce(rank_min,rank_max_owner4,nsize4,MPI_INTEGER,&
               MPI_MIN,PAR_COMM_TO_USE,istat4)
          rank_max_owner = int(rank_max_owner4,8)
       end if
#endif
    else
       if( present(rank_max_owner) ) then
          rank_max_owner = -1
       end if
    end if

  end subroutine PAR_MAX_8_s

  subroutine PAR_MAX_IP_0(n,xx,wherein,PAR_COMM_IN4)
    implicit none
    integer(ip),  intent(in)    :: n
    integer(ip),  intent(inout) :: xx(n)
    character(*), optional      :: wherein
    integer(4),   optional      :: PAR_COMM_IN4
    integer(ip)                 :: ii,nsize
    integer(4)                  :: istat4,nsize4
    integer(4)                  :: PAR_COMM_TO_USE
    integer(ip),  allocatable   :: yy(:)
    character(30)               :: my_wherein

    if( IPARALL .and. n > 0 ) then
#ifndef MPI_OFF
       istat4 = 0_4
       my_wherein = 'IN MY CODE'
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             my_wherein = trim(wherein)
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
          end if
       end if
       nsize  = n
       nsize4 = int(nsize,4)
       allocate(yy(n))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do ii = 1,nsize
             yy(ii) = -huge(1_ip)
          end do
       else
          do ii = 1,nsize
             yy(ii) = xx(ii)
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_MAX_IP_0: MPI ERROR')
       deallocate(yy)
#endif
    end if

  end subroutine PAR_MAX_IP_0

  subroutine PAR_MAX_IP_1(xx,wherein)
    implicit none
    integer(ip),     pointer, intent(inout) :: xx(:)
    character(*),    optional               :: wherein
    integer(ip)                             :: ii,nsize
    integer(4)                              :: istat4,nsize4
    integer(4)                              :: PAR_COMM_TO_USE
    integer(ip),     pointer                :: yy(:)
    character(30)                           :: my_wherein

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize  = size(xx,1)
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1)) )
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) = -huge(1_ip)
          end do
       else
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) = xx(ii)
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_MAX_IP_1: MPI ERROR')
       deallocate( yy )
#endif
    end if

  end subroutine PAR_MAX_IP_1

  subroutine PAR_MAX_IP_2(xx,wherein)
    implicit none
    integer(ip),     pointer, intent(inout) :: xx(:,:)
    character(*),    optional               :: wherein
    integer(ip)                             :: ii,jj,nsize1,nsize2,nsize
    integer(4)                              :: istat4,nsize4
    integer(4)                              :: PAR_COMM_TO_USE
    integer(ip),     pointer                :: yy(:,:)
    character(30)                           :: my_wherein

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize1 = size(xx,1)
       nsize2 = size(xx,2)
       nsize  = nsize1*nsize2
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1),lbound(xx,2):ubound(xx,2)))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = -huge(1_ip)
             end do
          end do
       else
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = xx(ii,jj)
             end do
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_MAX_IP_2: MPI ERROR')
       deallocate( yy )
#endif
    end if

  end subroutine PAR_MAX_IP_2

  !----------------------------------------------------------------------
  !
  ! OR FOR REALS
  !
  !----------------------------------------------------------------------

  subroutine PAR_OR_LG_s(xx_lg,wherein)
    implicit none
    logical(lg),  intent(inout) :: xx_lg
    character(*), optional      :: wherein
    integer(4)                  :: istat4,nsize4
    integer(4)                  :: PAR_COMM_TO_USE
    integer(ip)                 :: yy,xx
    character(30)               :: my_wherein

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize4 = 1_4
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          yy =  0
       else
          if( xx_lg ) then
             yy = 1
          else
             yy = 0
          end if
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_INTEGER,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_MIN_RP_s: MPI ERROR')
       if( xx == 1 ) then
          xx_lg = .true.
       else
          xx_lg = .false.
       end if
#endif
    end if

  end subroutine PAR_OR_LG_s

  !----------------------------------------------------------------------
  !
  ! MIN FOR REALS
  !
  !----------------------------------------------------------------------

  subroutine PAR_MIN_RP_s(xx,wherein,who,rank_max_owner)
    implicit none
    real(rp),     intent(inout) :: xx
    character(*), optional      :: wherein
    character(*), optional      :: who
    integer(ip),  optional      :: rank_max_owner
    integer(ip)                 :: nsize
    integer(4)                  :: istat4,nsize4
    integer(4)                  :: PAR_COMM_TO_USE
    integer(4)                  :: rank_min
    integer(4)                  :: rank_max_owner4
    real(rp)                    :: yy
    logical(lg)                 :: not_include_master
    character(30)               :: my_wherein

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       not_include_master = .true.
       if( present(who) ) then
          if( trim(who) == 'INCLUDE MASTER' ) not_include_master = .false.
       end if

       nsize  = 1
       nsize4 = int(nsize,4)
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. not_include_master ) then
          yy =  huge(1.0_rp)
       else
          yy = xx
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_PRECISION,&
            MPI_MIN,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_MIN_RP_s: MPI ERROR')
       if( present(rank_max_owner) ) then
          rank_min = huge(1_4)
          if( abs(yy-xx) < epsilon(1.0_rp) ) call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,rank_min)
          call MPI_AllReduce(rank_min,rank_max_owner4,nsize4,MPI_INTEGER,&
               MPI_MIN,PAR_COMM_TO_USE,istat4)
          rank_max_owner = int(rank_max_owner4,ip)
       end if
#endif
    else
       if( present(rank_max_owner) ) then
          rank_max_owner = -1
       end if
    end if

  end subroutine PAR_MIN_RP_s

  subroutine PAR_MIN_RP_0(n,xx,wherein,who,PAR_COMM_IN4)
    implicit none
    integer(ip),  intent(in)           :: n
    real(rp),     intent(inout)        :: xx(n)
    character(*), optional, intent(in) :: wherein
    character(*), optional, intent(in) :: who
    integer(4),   optional, intent(in) :: PAR_COMM_IN4
    integer(ip)                        :: ii,nsize
    integer(4)                         :: istat4,nsize4
    integer(4)                         :: PAR_COMM_TO_USE
    real(rp),     allocatable          :: yy(:)
    logical(lg)                        :: not_include_master
    character(30)                      :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    my_wherein = 'IN MY CODE'
    if( IPARALL .and. n > 0 ) then
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             my_wherein = trim(wherein)
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
          end if
       end if
       not_include_master = .true.
       if( present(who) ) then
          if( trim(who) == 'INCLUDE MASTER' ) not_include_master = .false.
       end if

       nsize  = n
       nsize4 = int(nsize,4)
       allocate(yy(n))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. not_include_master ) then
          do ii = 1,nsize
             yy(ii) =  huge(1.0_rp)
          end do
       else
          do ii = 1,nsize
             yy(ii) = xx(ii)
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_PRECISION,&
            MPI_MIN,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_MIN_RP_0: MPI ERROR')
       deallocate(yy)
    end if
#endif

  end subroutine PAR_MIN_RP_0

  subroutine PAR_MIN_RP_1(xx,wherein,who,PAR_COMM_IN4)
    implicit none
    real(rp),     pointer,  intent(inout) :: xx(:)
    character(*), optional, intent(in)    :: wherein
    character(*), optional, intent(in)    :: who
    integer(4),   optional, intent(in)    :: PAR_COMM_IN4
    integer(ip)                           :: ii,nsize
    integer(4)                            :: istat4,nsize4
    integer(4)                            :: PAR_COMM_TO_USE
    real(rp),     pointer                 :: yy(:)
    logical(lg)                           :: not_include_master
    character(30)                         :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    my_wherein = 'IN MY CODE'
    if( IPARALL ) then
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             my_wherein = trim(wherein)
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
          end if
       end if
       not_include_master = .true.
       if( present(who) ) then
          if( trim(who) == 'INCLUDE MASTER' ) not_include_master = .false.
       end if
       nsize  = size(xx,1)
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1)) )
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. not_include_master ) then
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) =  huge(1.0_rp)
          end do
       else
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) = xx(ii)
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_PRECISION,&
            MPI_MIN,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_MIN_RP_1: MPI ERROR')
       deallocate( yy )
    end if
#endif

  end subroutine PAR_MIN_RP_1

  subroutine PAR_MIN_RP_2(xx,wherein)
    implicit none
    real(rp),     pointer, intent(inout) :: xx(:,:)
    character(*),          optional      :: wherein
    integer(ip)                          :: ii,jj,nsize1,nsize2,nsize
    integer(4)                           :: istat4,nsize4
    integer(4)                           :: PAR_COMM_TO_USE
    real(rp),     pointer                :: yy(:,:)
    character(30)                        :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    my_wherein = 'IN MY CODE'
    if( IPARALL ) then
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize1 = size(xx,1)
       nsize2 = size(xx,2)
       nsize  = nsize1*nsize2
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1),lbound(xx,2):ubound(xx,2)))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) =  huge(1.0_rp)
             end do
          end do
       else
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = xx(ii,jj)
             end do
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_PRECISION,&
            MPI_MIN,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_MIN_RP_2: MPI ERROR')
       deallocate( yy )
    end if
#endif

  end subroutine PAR_MIN_RP_2

  !----------------------------------------------------------------------
  !
  ! MIN FOR COMPLEX
  !
  !----------------------------------------------------------------------

  subroutine PAR_MIN_CX_s(xx,wherein)
    implicit none
    complex(rp),  intent(inout) :: xx
    character(*), optional      :: wherein
    integer(ip)                 :: nsize
    integer(4)                  :: istat4,nsize4
    integer(4)                  :: PAR_COMM_TO_USE
    complex(rp)                 :: yy
    character(30)               :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    my_wherein = 'IN MY CODE'
    if( IPARALL ) then
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize  = 1
       nsize4 = int(nsize,4)
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          yy = huge(1.0_rp)
       else
          yy = xx
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_COMPLEX,&
            MPI_MIN,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_MIN_CX_s: MPI ERROR')
    end if
#endif

  end subroutine PAR_MIN_CX_s

  subroutine PAR_MIN_CX_0(n,xx,wherein)
    implicit none
    integer(ip),  intent(in)    :: n
    complex(rp),  intent(inout) :: xx(n)
    character(*), optional      :: wherein
    integer(ip)                 :: ii,nsize
    integer(4)                  :: istat4,nsize4
    integer(4)                  :: PAR_COMM_TO_USE
    complex(rp),  allocatable   :: yy(:)
    character(30)               :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    my_wherein = 'IN MY CODE'
    if( IPARALL .and. n > 0 ) then
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize  = n
       nsize4 = int(nsize,4)
       allocate(yy(n))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do ii = 1,nsize
             yy(ii) = huge(1.0_rp)
          end do
       else
          do ii = 1,nsize
             yy(ii) = xx(ii)
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_COMPLEX,&
            MPI_MIN,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_MIN_CX_0: MPI ERROR')
       deallocate(yy)
    end if
#endif

  end subroutine PAR_MIN_CX_0

  subroutine PAR_MIN_CX_1(xx,wherein)
    implicit none
    complex(rp),  pointer, intent(inout) :: xx(:)
    character(*),          optional      :: wherein
    integer(ip)                          :: ii,nsize
    integer(4)                           :: istat4,nsize4
    integer(4)                           :: PAR_COMM_TO_USE
    complex(rp),  pointer                :: yy(:)
    character(30)                        :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize  = size(xx,1)
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1)) )
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) = huge(1.0_rp)
          end do
       else
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) = xx(ii)
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_COMPLEX,&
            MPI_MIN,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_MIN_CX_1: MPI ERROR')
       deallocate( yy )
    end if
#endif

  end subroutine PAR_MIN_CX_1

  subroutine PAR_MIN_CX_2(xx,wherein)
    implicit none
    complex(rp),  pointer, intent(inout) :: xx(:,:)
    character(*), optional               :: wherein
    integer(ip)                          :: ii,jj,nsize1,nsize2,nsize
    integer(4)                           :: istat4,nsize4
    integer(4)                           :: PAR_COMM_TO_USE
    complex(rp),  pointer                :: yy(:,:)
    character(30)                        :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize1 = size(xx,1)
       nsize2 = size(xx,2)
       nsize  = nsize1*nsize2
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1),lbound(xx,2):ubound(xx,2)))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = huge(1.0_rp)
             end do
          end do
       else
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = xx(ii,jj)
             end do
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_DOUBLE_COMPLEX,&
            MPI_MIN,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_MIN_CX_2: MPI ERROR')
       deallocate( yy )
    end if
#endif

  end subroutine PAR_MIN_CX_2

  !----------------------------------------------------------------------
  !
  ! MIN FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_MIN_4_s(xx,wherein,rank_max_owner,PAR_COMM_IN4,INCLUDE_ROOT)
    implicit none
    integer(4),   intent(inout)        :: xx
    character(*),             optional :: wherein
    integer(4),               optional :: rank_max_owner
    integer(4),   intent(in), optional :: PAR_COMM_IN4
    logical(lg),  intent(in), optional :: INCLUDE_ROOT
    integer(ip)                        :: nsize
    integer(4)                         :: istat4,nsize4
    integer(4)                         :: PAR_COMM_TO_USE
    integer(4)                         :: yy
    integer(4)                         :: rank_min
    integer(4)                         :: rank_max_owner4
    logical(lg)                        :: if_include_root
    character(30)                      :: my_wherein

    istat4 = 0_4
    if( IPARALL ) then
#ifndef MPI_OFF
       my_wherein = 'IN MY CODE'
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       if( present(INCLUDE_ROOT) ) then
          if_include_root = INCLUDE_ROOT
       else
          if_include_root = .false.
       end if

       nsize  = 1
       nsize4 = int(nsize,4)
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. .not. if_include_root ) then
          yy =  huge(4_4)
       else
          yy = xx
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_INTEGER,&
            MPI_MIN,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_MIN_4_S: MPI ERROR')
       if( present(rank_max_owner) ) then
          rank_min = huge(1_4)
          if( abs(yy-xx) == 0_4 ) call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,rank_min)
          call MPI_AllReduce(rank_min,rank_max_owner4,nsize4,MPI_INTEGER,&
               MPI_MIN,PAR_COMM_TO_USE,istat4)
          rank_max_owner = rank_max_owner4
       end if
#endif
    else
       if( present(rank_max_owner) ) then
          rank_max_owner = -1
       end if
    end if

  end subroutine PAR_MIN_4_s

  subroutine PAR_MIN_8_s(xx,wherein,rank_max_owner,PAR_COMM_IN4,INCLUDE_ROOT)
    implicit none
    integer(8),   intent(inout)        :: xx
    character(*), optional             :: wherein
    integer(8),   optional             :: rank_max_owner
    integer(4),   intent(in), optional :: PAR_COMM_IN4
    logical(lg),  intent(in), optional :: INCLUDE_ROOT
    integer(ip)                        :: nsize
    integer(4)                         :: istat4,nsize4
    integer(4)                         :: PAR_COMM_TO_USE
    integer(8)                         :: yy
    integer(4)                         :: rank_min
    integer(4)                         :: rank_max_owner4
    logical(lg)                        :: if_include_root
    character(30)                      :: my_wherein

    istat4 = 0_4
    if( IPARALL ) then
#ifndef MPI_OFF
       my_wherein = 'IN MY CODE'
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       if( present(INCLUDE_ROOT) ) then
          if_include_root = INCLUDE_ROOT
       else
          if_include_root = .false.
       end if
       nsize  = 1
       nsize4 = int(nsize,4)
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' .and. .not. if_include_root ) then
          yy =  huge(1_8)
       else
          yy = xx
       end if
       call MPI_AllReduce(yy,xx,nsize4,MPI_INTEGER8,&
            MPI_MIN,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_MIN_8_S: MPI ERROR')
       if( present(rank_max_owner) ) then
          rank_min = huge(1_4)
          if( abs(yy-xx) == 0_8 ) call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,rank_min)
          call MPI_AllReduce(rank_min,rank_max_owner4,nsize4,MPI_INTEGER,&
               MPI_MIN,PAR_COMM_TO_USE,istat4)
          rank_max_owner = int(rank_max_owner4,8)
       end if
#endif
    else
       if( present(rank_max_owner) ) then
          rank_max_owner = -1
       end if
    end if

  end subroutine PAR_MIN_8_s


  subroutine PAR_MIN_IP_0(n,xx,wherein)
    implicit none
    integer(ip),  intent(in)    :: n
    integer(ip),  intent(inout) :: xx(n)
    character(*), optional      :: wherein
    integer(ip)                 :: ii,nsize
    integer(4)                  :: istat4,nsize4
    integer(4)                  :: PAR_COMM_TO_USE
    integer(ip),  allocatable   :: yy(:)
    character(30)               :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL .and. n > 0 ) then
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize  = n
       nsize4 = int(nsize,4)
       allocate(yy(n))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do ii = 1,nsize
             yy(ii) =  huge(1_ip)
          end do
       else
          do ii = 1,nsize
             yy(ii) = xx(ii)
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_MIN,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_MIN_IP_0: MPI ERROR')
       deallocate(yy)
    end if
#endif

  end subroutine PAR_MIN_IP_0

  subroutine PAR_MIN_IP_1(xx,wherein)
    implicit none
    integer(ip),     pointer, intent(inout) :: xx(:)
    character(*),    optional               :: wherein
    integer(ip)                             :: ii,nsize
    integer(4)                              :: istat4,nsize4
    integer(4)                              :: PAR_COMM_TO_USE
    integer(ip),     pointer                :: yy(:)
    character(30)                           :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize  = size(xx,1)
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1)) )
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) =  huge(1_ip)
          end do
       else
          do ii = lbound(yy,1),ubound(yy,1)
             yy(ii) = xx(ii)
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_MIN,PAR_COMM_TO_USE,istat4)
       deallocate( yy )
       if( istat4 /= 0_4 ) call runend('PAR_MIN_IP_1: MPI ERROR')
    end if
#endif

  end subroutine PAR_MIN_IP_1

  subroutine PAR_MIN_IP_2(xx,wherein)
    implicit none
    integer(ip),     pointer, intent(inout) :: xx(:,:)
    character(*),    optional               :: wherein
    integer(ip)                             :: ii,jj,nsize1,nsize2,nsize
    integer(4)                              :: istat4,nsize4
    integer(4)                              :: PAR_COMM_TO_USE
    integer(ip),     pointer                :: yy(:,:)
    character(30)                           :: my_wherein

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       my_wherein = 'IN MY CODE'
       if( present(wherein) ) then
          my_wherein = trim(wherein)
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       nsize1 = size(xx,1)
       nsize2 = size(xx,2)
       nsize  = nsize1*nsize2
       nsize4 = int(nsize,4)
       allocate( yy(lbound(xx,1):ubound(xx,1),lbound(xx,2):ubound(xx,2)))
       if( PAR_MY_CODE_RANK == 0 .and. trim(my_wherein) /= 'IN THE UNIVERSE' ) then
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) =  huge(1_ip)
             end do
          end do
       else
          do jj = lbound(yy,2),ubound(yy,2)
             do ii = lbound(yy,1),ubound(yy,1)
                yy(ii,jj) = xx(ii,jj)
             end do
          end do
       end if
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_MIN,PAR_COMM_TO_USE,istat4)
       deallocate( yy )
       if( istat4 /= 0_4 ) call runend('PAR_MIN_IP_2: MPI ERROR')
    end if
#endif

  end subroutine PAR_MIN_IP_2
  !
  ! SUM for integers in the Alya world (PAR_COMM_WORLD) including masters
  !
  subroutine PAR_SUM_ALL_IP_1(xx)

    implicit none
    integer(ip),  pointer, intent(inout) :: xx(:)
    integer(ip)                          :: ii,nsize
    integer(4)                           :: istat4,nsize4
    integer(4)                           :: PAR_COMM_TO_USE
    integer(ip),  pointer                :: yy(:)

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       call PAR_DEFINE_COMMUNICATOR('IN THE WORLD',PAR_COMM_TO_USE)
       nsize = size(xx,1)
       nsize4 = int(nsize,4)
       allocate(yy(lbound(xx,1):ubound(xx,1)))
       do ii =  lbound(yy,1), ubound(yy,1)
          yy(ii) = xx(ii)
       end do
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_SUM,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_SUM_ALL_IP_1: MPI ERROR')
       deallocate(yy)
    end if
#endif

  end subroutine PAR_SUM_ALL_IP_1

  subroutine PAR_SUM_ALL_IP_3(xx)

    implicit none
    integer(ip),  pointer, intent(inout) :: xx(:,:,:)
    integer(ip)                          :: ii,jj,kk,nsize1,nsize2,nsize3,nsize
    integer(4)                           :: istat4,nsize4
    integer(4)                           :: PAR_COMM_TO_USE
    integer(ip),  pointer                :: yy(:,:,:)

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       call PAR_DEFINE_COMMUNICATOR('IN THE WORLD',PAR_COMM_TO_USE)
       nsize1  = size(xx,1)
       nsize2  = size(xx,2)
       nsize3  = size(xx,3)
       nsize   = nsize1 * nsize2 * nsize3
       nsize4 = int(nsize,4)
       allocate(yy(lbound(xx,1):ubound(xx,1),lbound(xx,2):ubound(xx,2),lbound(xx,3):ubound(xx,3)))
       do kk = lbound(yy,3), ubound(yy,3)
          do jj = lbound(yy,2), ubound(yy,2)
             do ii =  lbound(yy,1), ubound(yy,1)
                yy(ii,jj,kk) = xx(ii,jj,kk)
             end do
          end do
       end do
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_SUM,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_SUM_ALL_IP_3: MPI ERROR')
       deallocate(yy)
    end if
#endif

  end subroutine PAR_SUM_ALL_IP_3
  !
  ! MAX for integers in the Alya world (PAR_COMM_WORLD) including masters
  !
  subroutine PAR_MAX_ALL_IP_s(xx)

    implicit none
    integer(ip),  intent(inout) :: xx
    integer(4)                  :: istat4
    integer(4)                  :: PAR_COMM_TO_USE
    integer(ip)                 :: yy

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       call PAR_DEFINE_COMMUNICATOR('IN THE WORLD',PAR_COMM_TO_USE)
       yy = xx
       call MPI_AllReduce(yy,xx,1_ip,PAR_INTEGER,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_MAX_ALL_IP_s: MPI ERROR')
    end if
#endif
  end subroutine PAR_MAX_ALL_IP_s

  subroutine PAR_MAX_ALL_IP_1(xx)

    implicit none
    integer(ip),  pointer, intent(inout) :: xx(:)
    integer(ip)                          :: ii,nsize
    integer(4)                           :: istat4,nsize4
    integer(4)                           :: PAR_COMM_TO_USE
    integer(ip),  pointer                :: yy(:)

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       call PAR_DEFINE_COMMUNICATOR('IN THE WORLD',PAR_COMM_TO_USE)
       nsize = size(xx,1)
       nsize4 = int(nsize,4)
       allocate(yy(lbound(xx,1):ubound(xx,1)))
       do ii =  lbound(yy,1), ubound(yy,1)
          yy(ii) = xx(ii)
       end do
       call MPI_AllReduce(yy,xx,nsize4,PAR_INTEGER,&
            MPI_MAX,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_SUM_ALL_IP_1: MPI ERROR')
       deallocate(yy)
    end if
#endif

  end subroutine PAR_MAX_ALL_IP_1

  !----------------------------------------------------------------------
  !
  ! Gather to rank = 0
  !
  !----------------------------------------------------------------------

  subroutine PAR_GATHER_CHARACTER(character_in,character_out,wherein)
    implicit none
    character(*),           intent(in)    :: character_in
    character(*), pointer,  intent(inout) :: character_out(:)
    character(*), optional, intent(in)    :: wherein
    character(1), pointer                 :: dummi(:)
    integer(4)                            :: PAR_COMM_TO_USE
    integer(4)                            :: l_character_in
    integer(4)                            :: l_character_out
    integer(4)                            :: istat4

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if

       l_character_in  = len(character_in)

       if( .not. associated(character_out) ) then
          allocate(dummi(l_character_in))
          call MPI_Gather(character_in, l_character_in, MPI_CHARACTER,&
               &          dummi,        l_character_in, MPI_CHARACTER,&
               &          0_4,PAR_COMM_TO_USE,istat4)
          deallocate(dummi)
       else
          l_character_out = len(character_out(1))
          if( l_character_in /= l_character_out ) call runend('PAR_GATHER: WRONG CHARACTER LENGTH')

          call MPI_Gather(character_in, l_character_in, MPI_CHARACTER,&
               &          character_out,l_character_out,MPI_CHARACTER,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       end if
       if( istat4 /= 0_4 ) call runend('PAR_GATHER_CHARACTER: MPI ERROR')
#endif
    else
       character_out = character_in
    end if

  end subroutine PAR_GATHER_CHARACTER

  subroutine PAR_GATHER_IP_s4(sendbuf,recvbuf,wherein,PAR_COMM_IN4)
    implicit none
    integer(4),             intent(in)    :: sendbuf
    integer(4),   pointer,  intent(inout) :: recvbuf(:)
    character(*), optional, intent(in)    :: wherein
    integer(4),   optional, intent(in)    :: PAR_COMM_IN4
    integer(4)                            :: dummi(2)
    integer(4)                            :: PAR_COMM_TO_USE
    integer(4)                            :: istat4,my_rank
    integer(ip)                           :: rl

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       call MPI_Comm_rank(PAR_COMM_TO_USE,my_rank,istat4)
       if( .not. associated(recvbuf) ) then
          if( my_rank == 0 ) call runend('PAR_GATHER_IP_s: NOT ASSOCIATED')
          call MPI_Gather(sendbuf,      1_4, MPI_INTEGER,&
               &          dummi,        1_4, MPI_INTEGER,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       else
          rl = lbound(recvbuf,1)
          call MPI_Gather(sendbuf,      1_4, MPI_INTEGER,&
               &          recvbuf(rl:), 1_4, MPI_INTEGER,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       end if
       if( istat4 /= 0_4 ) call runend('PAR_GATHER_IP_s4: MPI ERROR')
#endif
    else
       recvbuf = sendbuf
    end if

  end subroutine PAR_GATHER_IP_s4

  subroutine PAR_GATHER_IP_s48(sendbuf,recvbuf,wherein,PAR_COMM_IN4)
    implicit none
    integer(8),             intent(in)    :: sendbuf
    integer(4),   pointer,  intent(inout) :: recvbuf(:)
    character(*), optional, intent(in)    :: wherein
    integer(4),   optional, intent(in)    :: PAR_COMM_IN4
    integer(4)                            :: dummi(2)
    integer(4)                            :: PAR_COMM_TO_USE
    integer(4)                            :: istat4,my_rank
    integer(ip)                           :: rl
    integer(4)                            :: sendbuf4

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       sendbuf4 = int(sendbuf,4)
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       call MPI_Comm_rank(PAR_COMM_TO_USE,my_rank,istat4)
       if( .not. associated(recvbuf) ) then
          if( my_rank == 0 ) call runend('PAR_GATHER_IP_s: NOT ASSOCIATED')
          call MPI_Gather(sendbuf4,     1_4, MPI_INTEGER,&
               &          dummi,        1_4, MPI_INTEGER,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       else
          rl = lbound(recvbuf,1)
          call MPI_Gather(sendbuf4,     1_4, MPI_INTEGER,&
               &          recvbuf(rl:), 1_4, MPI_INTEGER,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       end if
       if( istat4 /= 0_4 ) call runend('PAR_GATHER_IP_s4: MPI ERROR')
#endif
    else
       recvbuf = int(sendbuf,4)
    end if

  end subroutine PAR_GATHER_IP_s48

  subroutine PAR_GATHER_IP_s8(sendbuf,recvbuf,wherein,PAR_COMM_IN4)
    implicit none
    integer(8),             intent(in)    :: sendbuf
    integer(8),   pointer,  intent(inout) :: recvbuf(:)
    character(*), optional, intent(in)    :: wherein
    integer(4),   optional, intent(in)    :: PAR_COMM_IN4
    integer(8)                            :: dummi(2)
    integer(4)                            :: PAR_COMM_TO_USE
    integer(4)                            :: istat4,my_rank
    integer(ip)                           :: rl

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       call MPI_Comm_rank(PAR_COMM_TO_USE,my_rank,istat4)
       if( .not. associated(recvbuf) ) then
          if( my_rank == 0 ) call runend('PAR_GATHER_IP_s: NOT ASSOCIATED')
          call MPI_Gather(sendbuf,      1_4, MPI_INTEGER8,&
               &          dummi,        1_4, MPI_INTEGER8,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       else
          rl = lbound(recvbuf,1)
          call MPI_Gather(sendbuf,      1_4, MPI_INTEGER8,&
               &          recvbuf(rl:), 1_4, MPI_INTEGER8,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       end if
       if( istat4 /= 0_4 ) call runend('PAR_GATHER_IP_s8: MPI ERROR')
#endif
    else
       recvbuf = sendbuf
    end if

  end subroutine PAR_GATHER_IP_s8

  subroutine PAR_GATHER_IP_14(sendbuf,recvbuf,wherein)
    implicit none
    integer(4),  pointer,  intent(in)    :: sendbuf(:)
    integer(4),  pointer,  intent(inout) :: recvbuf(:)
    character(*),          intent(in)    :: wherein
    integer(4)                           :: dummi(2)
    integer(4)                           :: PAR_COMM_TO_USE
    integer(4)                           :: istat4,my_rank
    integer(4)                           :: sendcount4
    integer(ip)                          :: sl,rl,lbouns,lbounr

    if( IPARALL ) then
#ifndef MPI_OFF
    istat4 = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call MPI_Comm_rank(PAR_COMM_TO_USE,my_rank,istat4)

       sendcount4 = size(sendbuf)
       sl = lbound(sendbuf,1)
       if( .not. associated(recvbuf) ) then
          if( my_rank == 0 ) call runend('PAR_GATHER_IP_s: NOT ASSOCIATED')
          call MPI_Gather(sendbuf(sl:), sendcount4,MPI_INTEGER,&
               &          dummi,        sendcount4,MPI_INTEGER,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       else
          rl = lbound(recvbuf,1)
          call MPI_Gather(sendbuf(sl:), sendcount4,MPI_INTEGER,&
               &          recvbuf(rl:), sendcount4,MPI_INTEGER,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       end if
       if( istat4 /= 0_4 ) call runend('PAR_GATHER_IP_14: MPI ERROR')
#endif
    else
       lbouns  = lbound(sendbuf,1)
       lbounr  = lbound(recvbuf,1)
       sl      = size(sendbuf)
       recvbuf(lbounr:lbounr+sl-1) = sendbuf(lbouns:lbouns+sl-1)
    end if

  end subroutine PAR_GATHER_IP_14

  subroutine PAR_GATHER_IP_18(sendbuf,recvbuf,wherein)
    implicit none
    integer(8),  pointer,  intent(in)    :: sendbuf(:)
    integer(8),  pointer,  intent(inout) :: recvbuf(:)
    character(*),          intent(in)    :: wherein
    integer(8)                           :: dummi(2)
    integer(4)                           :: PAR_COMM_TO_USE
    integer(4)                           :: istat4,my_rank
    integer(4)                           :: sendcount4
    integer(ip)                          :: sl,rl,lbouns,lbounr

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call MPI_Comm_rank(PAR_COMM_TO_USE,my_rank,istat4)

       sendcount4 = size(sendbuf)
       sl = lbound(sendbuf,1)
       if( .not. associated(recvbuf) ) then
          if( my_rank == 0 ) call runend('PAR_GATHER_IP_s: NOT ASSOCIATED')
          call MPI_Gather(sendbuf(sl:), sendcount4,MPI_INTEGER8,&
               &          dummi,        sendcount4,MPI_INTEGER8,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       else
          rl = lbound(recvbuf,1)
          call MPI_Gather(sendbuf(sl:), sendcount4,MPI_INTEGER8,&
               &          recvbuf(rl:), sendcount4,MPI_INTEGER8,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       end if
       if( istat4 /= 0_4 ) call runend('PAR_GATHER_IP_18: MPI ERROR')
#endif
    else
       lbouns  = lbound(sendbuf,1)
       lbounr  = lbound(recvbuf,1)
       sl      = size(sendbuf)
       recvbuf(lbounr:lbounr+sl-1) = sendbuf(lbouns:lbouns+sl-1)
    end if

  end subroutine PAR_GATHER_IP_18

  subroutine PAR_GATHER_IP_12(sendbuf,recvbuf,wherein)
    implicit none
    integer(ip),  pointer, intent(in)    :: sendbuf(:)
    integer(ip),  pointer, intent(inout) :: recvbuf(:,:)
    character(*),          intent(in)    :: wherein
    integer(4)                           :: dummr(2)
    integer(4)                           :: PAR_COMM_TO_USE
    integer(4)                           :: istat4,my_rank
    integer(4)                           :: sendcount4
    integer(ip)                          :: rl1,rl2,lbouns,sl,ii

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank)

       sendcount4 = size(sendbuf)
       if( .not. associated(recvbuf) ) then
          if( my_rank == 0 ) call runend('PAR_GATHER_IP_s: NOT ASSOCIATED')
          call MPI_Gather(sendbuf,  sendcount4,PAR_INTEGER,&
               &          dummr,    sendcount4,PAR_INTEGER,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       else
          rl1 = lbound(recvbuf,1)
          rl2 = lbound(recvbuf,2)
          call MPI_Gather(sendbuf,           sendcount4,PAR_INTEGER,&
               &          recvbuf(rl1:,rl2:),sendcount4,PAR_INTEGER,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       end if
       if( istat4 /= 0_4 ) call runend('PAR_GATHER_IP_12: MPI ERROR')
#endif
    else
       lbouns = lbound(sendbuf,1)
       rl1    = lbound(recvbuf,1)
       rl2    = lbound(recvbuf,2)
       sl     = size(sendbuf)
       do ii = rl2,ubound(recvbuf,2)
          recvbuf(rl1:rl1+sl-1,ii) = sendbuf(lbouns:lbouns+sl-1)
       end do
    end if

  end subroutine PAR_GATHER_IP_12

  subroutine PAR_GATHER_IP_23(sendbuf,recvbuf,wherein)
    implicit none
    integer(ip),  pointer, intent(in)    :: sendbuf(:,:)
    integer(ip),  pointer, intent(inout) :: recvbuf(:,:,:)
    character(*),          intent(in)    :: wherein
    integer(4)                           :: dummr(2)
    integer(4)                           :: PAR_COMM_TO_USE
    integer(4)                           :: istat4,my_rank
    integer(4)                           :: sendcount4
    integer(ip)                          :: rl1,rl2,rl3

#ifndef MPI_OFF
    istat4 = 0_4
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank)

    sendcount4 = size(sendbuf)
    if( .not. associated(recvbuf) ) then
       if( my_rank == 0 ) call runend('PAR_GATHER_IP_s: NOT ASSOCIATED')
       call MPI_Gather(sendbuf,  sendcount4,PAR_INTEGER,&
            &          dummr,    sendcount4,PAR_INTEGER,&
            &          0_4,PAR_COMM_TO_USE,istat4)
    else
       rl1 = lbound(recvbuf,1)
       rl2 = lbound(recvbuf,2)
       rl3 = lbound(recvbuf,3)
       call MPI_Gather(sendbuf,                sendcount4,PAR_INTEGER,&
            &          recvbuf(rl1:,rl2:,rl3:),sendcount4,PAR_INTEGER,&
            &          0_4,PAR_COMM_TO_USE,istat4)
    end if
    if( istat4 /= 0_4 ) call runend('PAR_GATHER_IP_12: MPI ERROR')
#endif

  end subroutine PAR_GATHER_IP_23

  subroutine PAR_GATHER_RP_s(sendbuf,recvbuf,wherein)
    implicit none
    real(rp),              intent(in)    :: sendbuf
    real(rp),     pointer, intent(inout) :: recvbuf(:)
    character(*),          intent(in)    :: wherein
    real(rp)                             :: dummr(2)
    integer(4)                           :: PAR_COMM_TO_USE
    integer(4)                           :: istat4,my_rank
    integer(ip)                          :: rl

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call MPI_Comm_rank(PAR_COMM_TO_USE,my_rank,istat4)
       if( .not. associated(recvbuf) ) then
          if( my_rank == 0 ) call runend('PAR_GATHER_IP_s: NOT ASSOCIATED')
          call MPI_Gather(sendbuf,      1_4, MPI_DOUBLE_PRECISION,&
               &          dummr,        1_4, MPI_DOUBLE_PRECISION,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       else
          rl = lbound(recvbuf,1)
          call MPI_Gather(sendbuf,      1_4, MPI_DOUBLE_PRECISION,&
               &          recvbuf(rl:), 1_4, MPI_DOUBLE_PRECISION,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       end if
       if( istat4 /= 0_4 ) call runend('PAR_GATHER_IP_s: MPI ERROR')
#endif
    else
       recvbuf = sendbuf
    end if

  end subroutine PAR_GATHER_RP_s

  subroutine PAR_GATHER_RP_1(sendbuf,recvbuf,wherein)
    implicit none
    real(rp),     pointer, intent(in)    :: sendbuf(:)
    real(rp),     pointer, intent(inout) :: recvbuf(:)
    character(*),          intent(in)    :: wherein
    real(rp)                             :: dummr(2)
    integer(4)                           :: PAR_COMM_TO_USE
    integer(4)                           :: istat4,my_rank
    integer(4)                           :: sendcount4
    integer(ip)                          :: sl,rl,lbouns,lbounr

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call MPI_Comm_rank(PAR_COMM_TO_USE,my_rank,istat4)

       sendcount4 = size(sendbuf)
       if( .not. associated(recvbuf) ) then
          if( my_rank == 0 ) call runend('PAR_GATHER_IP_s: NOT ASSOCIATED')
          call MPI_Gather(sendbuf, sendcount4,MPI_DOUBLE_PRECISION,&
               &          dummr,   sendcount4,MPI_DOUBLE_PRECISION,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       else
          rl = lbound(recvbuf,1)
          call MPI_Gather(sendbuf,      sendcount4,MPI_DOUBLE_PRECISION,&
               &          recvbuf(rl:), sendcount4,MPI_DOUBLE_PRECISION,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       end if
       if( istat4 /= 0_4 ) call runend('PAR_GATHER_RP_1: MPI ERROR')
#endif
    else
       lbouns  = lbound(sendbuf,1)
       lbounr  = lbound(recvbuf,1)
       sl      = size(sendbuf)
       recvbuf(lbounr:lbounr+sl-1) = sendbuf(lbouns:lbouns+sl-1)
    end if

  end subroutine PAR_GATHER_RP_1

  subroutine PAR_GATHER_RP_12(sendbuf,recvbuf,wherein)
    implicit none
    real(rp),     pointer, intent(in)    :: sendbuf(:)
    real(rp),     pointer, intent(inout) :: recvbuf(:,:)
    character(*),          intent(in)    :: wherein
    real(rp)                             :: dummr(2)
    integer(4)                           :: PAR_COMM_TO_USE
    integer(4)                           :: istat4,my_rank
    integer(4)                           :: sendcount4
    integer(ip)                          :: rl1,rl2,lbouns,sl

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank)

       sendcount4 = size(sendbuf)
       if( .not. associated(recvbuf) ) then
          if( my_rank == 0 ) call runend('PAR_GATHER_IP_s: NOT ASSOCIATED')
          call MPI_Gather(sendbuf,  sendcount4,MPI_DOUBLE_PRECISION,&
               &          dummr,    sendcount4,MPI_DOUBLE_PRECISION,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       else
          rl1 = lbound(recvbuf,1)
          rl2 = lbound(recvbuf,2)
          call MPI_Gather(sendbuf,           sendcount4,MPI_DOUBLE_PRECISION,&
               &          recvbuf(rl1:,rl2:),sendcount4,MPI_DOUBLE_PRECISION,&
               &          0_4,PAR_COMM_TO_USE,istat4)
       end if
       if( istat4 /= 0_4 ) call runend('PAR_GATHER_RP_12: MPI ERROR')
#endif
    else
       lbouns = lbound(sendbuf,1)
       rl1    = lbound(recvbuf,1)
       rl2    = lbound(recvbuf,2)
       sl     = size(sendbuf)
       recvbuf(rl1:rl1+sl-1,rl2) = sendbuf(lbouns:lbouns+sl-1)
    end if

  end subroutine PAR_GATHER_RP_12

  !----------------------------------------------------------------------
  !
  ! Scatter from rank = 0
  !
  !----------------------------------------------------------------------

  subroutine PAR_SCATTER_IP_s(xx_in,xx_out,wherein)
    implicit none
    integer(ip),  pointer, intent(in)  :: xx_in(:)
    integer(ip),           intent(out) :: xx_out
    character(*),          intent(in)  :: wherein
    integer(ip)                        :: dummi(2)
    integer(4)                         :: PAR_COMM_TO_USE
    integer(4)                         :: istat4

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       if( .not. associated(xx_in) ) then
          call MPI_Scatter(dummi,1_4,PAR_INTEGER,xx_out,1_4,PAR_INTEGER,0_4,PAR_COMM_TO_USE,istat4)
       else
          call MPI_Scatter(xx_in,1_4,PAR_INTEGER,xx_out,1_4,PAR_INTEGER,0_4,PAR_COMM_TO_USE,istat4)
       end if
       if( istat4 /= 0_4 ) call runend('PAR_SCATTER_IP_s: MPI ERROR 1')
    end if
#endif

  end subroutine PAR_SCATTER_IP_s

  subroutine PAR_SCATTERV_IP_1(sendbuf,recvbuf,sendcount4,wherein,displs4)
    implicit none
    integer(ip),  pointer, intent(in)           :: sendbuf(:)           !< Send buffer
    integer(ip),  pointer, intent(inout)        :: recvbuf(:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)        !< Recv counts
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: recvcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rl,sl,scl,dl,jpart
    integer(4),   pointer                       :: sendcount4_tmp(:)
    integer(4),   target                        :: sendcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       root_rank4      = 0_4
       sendcount4_null = 0_4
       displs4_null    = 0_4
       recvcount4      = 0_4

       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)

       if( associated(recvbuf) ) then

          recvcount4 = size(recvbuf)
          rl         = lbound(recvbuf,1_ip)

       end if

       if( associated(sendbuf)    ) sl = lbound(sendbuf,1_ip)

!       if( associated(sendcount4) ) then
       if( my_rank == root_rank4 ) then

          scl = lbound(sendcount4,1_4)
          sendcount4_tmp => sendcount4(scl:)

       else

          sendcount4_tmp => sendcount4_null

       end if

       if( present(displs4) ) then

          if( associated(displs4) ) then

             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)

          else

             displs4_tmp => displs4_null

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,              PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,                  PAR_INTEGER,&
                     &            root_rank4,  PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                  PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                      PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

       else

          if( my_rank == root_rank4 ) then

             allocate( my_displs4(0_4:comm_size-1_4) )
             jpart = lbound(sendcount4,1_4)
             my_displs4(0_4) = 0_4
             do ipart = 1,comm_size-1

                my_displs4(ipart) = my_displs4(ipart-1) + sendcount4(jpart)
                jpart = jpart + 1

             end do

          else

             allocate( my_displs4(1_4) )
             my_displs4(1_4) = 0_4

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,             PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,                 PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                 PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                     PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

          deallocate( my_displs4 )

       end if

       if( istat4 /= 0 ) call runend('PAR_SCATTERV_IP_1: MPI ERROR')

    end if

#endif

  end subroutine PAR_SCATTERV_IP_1

  subroutine PAR_SCATTERV_IP_2(sendbuf,recvbuf,sendcount4,wherein,displs4)
    implicit none
    integer(ip),  pointer, intent(in)           :: sendbuf(:,:)           !< Send buffer
    integer(ip),  pointer, intent(inout)          :: recvbuf(:,:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)        !< Recv counts
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: recvcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rl,sl,scl,dl,jpart
    integer(4),   pointer                       :: sendcount4_tmp(:)
    integer(4),   target                        :: sendcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       root_rank4      = 0_4
       sendcount4_null = 0_4
       displs4_null    = 0_4
       recvcount4      = 0_4

       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)

       if( associated(recvbuf) ) then

          recvcount4 = size(recvbuf)
          rl         = lbound(recvbuf,1_ip)

       end if

       if( associated(sendbuf)    ) sl = lbound(sendbuf,1_ip)

!       if( associated(sendcount4) ) then
       if( my_rank == root_rank4 ) then

          scl = lbound(sendcount4,1_4)
          sendcount4_tmp => sendcount4(scl:)

       else

          sendcount4_tmp => sendcount4_null

       end if

       if( present(displs4) ) then

          if( associated(displs4) ) then

             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)

          else

             displs4_tmp => displs4_null

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,              PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,                  PAR_INTEGER,&
                     &            root_rank4,  PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                  PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                      PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

       else

          if( my_rank == root_rank4 ) then

             allocate( my_displs4(0_4:comm_size-1_4) )
             jpart = lbound(sendcount4,1_4)
             my_displs4(0_4) = 0_4
             do ipart = 1,comm_size-1

                my_displs4(ipart) = my_displs4(ipart-1) + sendcount4(jpart)
                jpart = jpart + 1

             end do

          else

             allocate( my_displs4(1_4) )
             my_displs4(1_4) = 0_4

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,             PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,                 PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                 PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                     PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

          deallocate( my_displs4 )

       end if

       if( istat4 /= 0 ) call runend('PAR_SCATTERV_IP_2: MPI ERROR')

    end if

#endif

  end subroutine PAR_SCATTERV_IP_2

  subroutine PAR_SCATTERV_RP_0(sendbuf,recvbuf,sendcount4,recvcount4,wherein,displs4)
    implicit none
    real(rp),     pointer, intent(in)           :: sendbuf(:)           !< Send buffer
    real(rp),              intent(inout)        :: recvbuf(*)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)        !< Send counts
    integer(4),            intent(in)           :: recvcount4           !< Recv count
    character(*),          intent(in)           :: wherein              !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4),   pointer                       :: my_displs4(:)
    real(rp)                                    :: sendbuf_tmp(2)
    real(rp)                                    :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rl,sl,scl,dl,jpart
    integer(4),   pointer                       :: sendcount4_tmp(:)
    integer(4),   target                        :: sendcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       root_rank4      = 0_4
       sendcount4_null = 0_4
       displs4_null    = 0_4

       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)

       rl = 1
       if( associated(sendbuf)    ) sl = lbound(sendbuf,1_ip)

!       if( associated(sendcount4) ) then
       if( my_rank == root_rank4 ) then

          scl = lbound(sendcount4,1_4)
          sendcount4_tmp => sendcount4(scl:)

       else

          sendcount4_tmp => sendcount4_null

       end if

       if( present(displs4) ) then

          if( associated(displs4) ) then

             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)

          else

             displs4_tmp => displs4_null

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,              MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,                  MPI_DOUBLE_PRECISION,&
                     &            root_rank4,  PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                  MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                      MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

       else

          if( my_rank == root_rank4 ) then

             allocate( my_displs4(0_4:comm_size-1_4) )
             jpart = lbound(sendcount4,1_4)
             my_displs4(0_4) = 0_4
             do ipart = 1,comm_size-1

                my_displs4(ipart) = my_displs4(ipart-1) + sendcount4(jpart)
                jpart = jpart + 1

             end do

          else

             allocate( my_displs4(1_4) )
             my_displs4(1_4) = 0_4

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,             MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,                 MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                 MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                     MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

          deallocate( my_displs4 )

       end if

       if( istat4 /= 0 ) call runend('PAR_SCATTERV_RP_1: MPI ERROR')

    end if

#endif

  end subroutine PAR_SCATTERV_RP_0

  subroutine PAR_SCATTERV_RP_1(sendbuf,recvbuf,sendcount4,wherein,displs4)
    implicit none
    real(rp),     pointer, intent(in)           :: sendbuf(:)           !< Send buffer
    real(rp),     pointer, intent(inout)        :: recvbuf(:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)        !< Send counts
    character(*),          intent(in)           :: wherein              !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: recvcount4
    integer(4),   pointer                       :: my_displs4(:)
    real(rp)                                    :: sendbuf_tmp(2)
    real(rp)                                    :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rl,sl,scl,dl,jpart
    integer(4),   pointer                       :: sendcount4_tmp(:)
    integer(4),   target                        :: sendcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       root_rank4      = 0_4
       sendcount4_null = 0_4
       displs4_null    = 0_4
       recvcount4      = 0_4

       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)

       if( associated(recvbuf) ) then

          recvcount4 = size(recvbuf)
          rl         = lbound(recvbuf,1_ip)

       end if

       if( associated(sendbuf)    ) sl = lbound(sendbuf,1_ip)

!       if( associated(sendcount4) ) then
       if( my_rank == root_rank4 ) then

          scl = lbound(sendcount4,1_4)
          sendcount4_tmp => sendcount4(scl:)

       else

          sendcount4_tmp => sendcount4_null

       end if

       if( present(displs4) ) then

          if( associated(displs4) ) then

             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)

          else

             displs4_tmp => displs4_null

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,              MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,                  MPI_DOUBLE_PRECISION,&
                     &            root_rank4,  PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                  MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                      MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

       else

          if( my_rank == root_rank4 ) then

             allocate( my_displs4(0_4:comm_size-1_4) )
             jpart = lbound(sendcount4,1_4)
             my_displs4(0_4) = 0_4
             do ipart = 1,comm_size-1

                my_displs4(ipart) = my_displs4(ipart-1) + sendcount4(jpart)
                jpart = jpart + 1

             end do

          else

             allocate( my_displs4(1_4) )
             my_displs4(1_4) = 0_4

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,             MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,                 MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                 MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                     MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

          deallocate( my_displs4 )

       end if

       if( istat4 /= 0 ) call runend('PAR_SCATTERV_RP_1: MPI ERROR')

    end if

#endif

  end subroutine PAR_SCATTERV_RP_1

  subroutine PAR_SCATTERV_RP_2(sendbuf,recvbuf,sendcount4,wherein,displs4)
    implicit none
    real(rp),  pointer, intent(in)              :: sendbuf(:,:)           !< Send buffer
    real(rp),  pointer, intent(inout)             :: recvbuf(:,:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)          !< Recv counts
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)             !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: recvcount4
    integer(4),   pointer                       :: my_displs4(:)
    real(rp)                                    :: sendbuf_tmp(2)
    real(rp)                                    :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rl,sl,scl,dl,jpart
    integer(4),    pointer                      :: sendcount4_tmp(:)
    integer(4),    target                       :: sendcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       root_rank4      = 0_4
       sendcount4_null = 0_4
       displs4_null    = 0_4
       recvcount4      = 0_4

       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)

       if( associated(recvbuf) ) then

          recvcount4 = size(recvbuf)
          rl         = lbound(recvbuf,1_ip)

       end if

       if( associated(sendbuf)    ) sl = lbound(sendbuf,1_ip)

       if( my_rank == root_rank4 ) then

          scl = lbound(sendcount4,1_4)
          sendcount4_tmp => sendcount4(scl:)

       else

          sendcount4_tmp => sendcount4_null

       end if

       if( present(displs4) ) then

          if( associated(displs4) ) then

             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)

          else

             displs4_tmp => displs4_null

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,              MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,                  MPI_DOUBLE_PRECISION,&
                     &            root_rank4,  PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                  MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                      MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

       else

          if( my_rank == root_rank4 ) then

             allocate( my_displs4(0_4:comm_size-1_4) )
             jpart = lbound(sendcount4,1_4)
             my_displs4(0_4) = 0_4
             do ipart = 1,comm_size-1

                my_displs4(ipart) = my_displs4(ipart-1) + sendcount4(jpart)
                jpart = jpart + 1

             end do

          else

             allocate( my_displs4(1_4) )
             my_displs4(1_4) = 0_4

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,             MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,                 MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                 MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                     MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

          deallocate( my_displs4 )

       end if

       if( istat4 /= 0 ) call runend('PAR_SCATTERV_RP_1: MPI ERROR')

    end if

#endif

  end subroutine PAR_SCATTERV_RP_2

  subroutine PAR_SCATTERV_IP_1_RCV(sendbuf,recvbuf,sendcount4,recvcount4,wherein,displs4)
    implicit none
    integer(ip),  pointer, intent(in)           :: sendbuf(:)           !< Send buffer
    integer(ip),  pointer, intent(inout)        :: recvbuf(:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)        !< Recv counts
    integer(4),            intent(in)           :: recvcount4
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rl,sl,scl,dl,jpart
    integer(4),   pointer                       :: sendcount4_tmp(:)
    integer(4),   target                        :: sendcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       root_rank4      = 0_4
       sendcount4_null = 0_4
       displs4_null    = 0_4

       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)

       if( associated(recvbuf) ) then

          rl         = lbound(recvbuf,1_ip)

       end if

       if( associated(sendbuf)    ) sl = lbound(sendbuf,1_ip)

       if( my_rank == root_rank4 ) then

          scl = lbound(sendcount4,1_4)
          sendcount4_tmp => sendcount4(scl:)

       else

          sendcount4_tmp => sendcount4_null

       end if

       if( present(displs4) ) then

          if( associated(displs4) ) then

             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)

          else

             displs4_tmp => displs4_null

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,              PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,                  PAR_INTEGER,&
                     &            root_rank4,  PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                  PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                      PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

       else

          if( my_rank == root_rank4 ) then

             allocate( my_displs4(0_4:comm_size-1_4) )
             jpart = lbound(sendcount4,1_4)
             my_displs4(0_4) = 0_4
             do ipart = 1,comm_size-1

                my_displs4(ipart) = my_displs4(ipart-1) + sendcount4(jpart)
                jpart = jpart + 1

             end do

          else

             allocate( my_displs4(1_4) )
             my_displs4(1_4) = 0_4

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,             PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,                 PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                 PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                     PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

          deallocate( my_displs4 )

       end if

       if( istat4 /= 0 ) call runend('PAR_SCATTERV_IP_1: MPI ERROR')

    end if

#endif

  end subroutine PAR_SCATTERV_IP_1_RCV

  subroutine PAR_SCATTERV_IP_2_RCV(sendbuf,recvbuf,sendcount4,recvcount4,wherein,displs4)
    implicit none
    integer(ip),  pointer, intent(in)           :: sendbuf(:,:)         !< Send buffer
    integer(ip),  pointer, intent(inout)        :: recvbuf(:,:)         !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)        !< Recv counts
    character(*),          intent(in)           :: wherein              !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4),            intent(in)           :: recvcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rl,sl,scl,dl,jpart
    integer(4),   pointer                       :: sendcount4_tmp(:)
    integer(4),   target                        :: sendcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       root_rank4      = 0_4
       sendcount4_null = 0_4
       displs4_null    = 0_4

       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)

       if( associated(recvbuf) ) then

          rl         = lbound(recvbuf,1_ip)

       end if

       if( associated(sendbuf)    ) sl = lbound(sendbuf,1_ip)

       if( my_rank == root_rank4 ) then

          scl = lbound(sendcount4,1_4)
          sendcount4_tmp => sendcount4(scl:)

       else

          sendcount4_tmp => sendcount4_null

       end if

       if( present(displs4) ) then

          if( associated(displs4) ) then

             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)

          else

             displs4_tmp => displs4_null

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,              PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,                  PAR_INTEGER,&
                     &            root_rank4,  PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                  PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                      PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

       else

          if( my_rank == root_rank4 ) then

             allocate( my_displs4(0_4:comm_size-1_4) )
             jpart = lbound(sendcount4,1_4)
             my_displs4(0_4) = 0_4
             do ipart = 1,comm_size-1

                my_displs4(ipart) = my_displs4(ipart-1) + sendcount4(jpart)
                jpart = jpart + 1

             end do

          else

             allocate( my_displs4(1_4) )
             my_displs4(1_4) = 0_4

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,             PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf_tmp, recvcount4,                 PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                 PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, PAR_INTEGER,&
                     &            recvbuf, recvcount4,                     PAR_INTEGER,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

          deallocate( my_displs4 )

       end if

       if( istat4 /= 0 ) call runend('PAR_SCATTERV_IP_2: MPI ERROR')

    end if

#endif

  end subroutine PAR_SCATTERV_IP_2_RCV

  subroutine PAR_SCATTERV_RP_1_RCV(sendbuf,recvbuf,sendcount4,recvcount4,wherein,displs4)
    implicit none
    real(rp),  pointer, intent(in)              :: sendbuf(:)           !< Send buffer
    real(rp),  pointer, intent(inout)           :: recvbuf(:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)        !< Recv counts
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4),    intent(in)                   :: recvcount4
    integer(4),   pointer                       :: my_displs4(:)
    real(rp)                                    :: sendbuf_tmp(2)
    real(rp)                                    :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rl,sl,scl,dl,jpart
    integer(4),    pointer                      :: sendcount4_tmp(:)
    integer(4),    target                       :: sendcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       root_rank4      = 0_4
       sendcount4_null = 0_4
       displs4_null    = 0_4

       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)

       if( associated(recvbuf) ) then

          rl         = lbound(recvbuf,1_ip)

       end if

       if( associated(sendbuf)    ) sl = lbound(sendbuf,1_ip)

       if( my_rank == root_rank4 ) then

          scl = lbound(sendcount4,1_4)
          sendcount4_tmp => sendcount4(scl:)

       else

          sendcount4_tmp => sendcount4_null

       end if

       if( present(displs4) ) then

          if( associated(displs4) ) then

             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)

          else

             displs4_tmp => displs4_null

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,              MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,                  MPI_DOUBLE_PRECISION,&
                     &            root_rank4,  PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                  MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                      MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

       else

          if( my_rank == root_rank4 ) then

             allocate( my_displs4(0_4:comm_size-1_4) )
             jpart = lbound(sendcount4,1_4)
             my_displs4(0_4) = 0_4
             do ipart = 1,comm_size-1

                my_displs4(ipart) = my_displs4(ipart-1) + sendcount4(jpart)
                jpart = jpart + 1

             end do

          else

             allocate( my_displs4(1_4) )
             my_displs4(1_4) = 0_4

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,             MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,                 MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                 MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                     MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

          deallocate( my_displs4 )

       end if

       if( istat4 /= 0 ) call runend('PAR_SCATTERV_RP_1: MPI ERROR')

    end if

#endif

  end subroutine PAR_SCATTERV_RP_1_RCV

  subroutine PAR_SCATTERV_R8_2_RCV(sendbuf,recvbuf,sendcount4,recvcount4,wherein,displs4)
    implicit none
    real(8),  pointer, intent(in)               :: sendbuf(:,:)           !< Send buffer
    real(8),  pointer, intent(inout)            :: recvbuf(:,:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)          !< Recv counts
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)             !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4),     intent(in)                  :: recvcount4
    integer(4),   pointer                       :: my_displs4(:)
    real(8)                                    :: sendbuf_tmp(2)
    real(8)                                    :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rl,sl,scl,dl,jpart
    integer(4),    pointer                      :: sendcount4_tmp(:)
    integer(4),    target                       :: sendcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       root_rank4      = 0_4
       sendcount4_null = 0_4
       displs4_null    = 0_4

       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)

       if( associated(recvbuf) ) then

          rl         = lbound(recvbuf,1_ip)

       end if

       if( associated(sendbuf)    ) sl = lbound(sendbuf,1_ip)

       if( my_rank == root_rank4 ) then

          scl = lbound(sendcount4,1_4)
          sendcount4_tmp => sendcount4(scl:)

       else

          sendcount4_tmp => sendcount4_null

       end if

       if( present(displs4) ) then

          if( associated(displs4) ) then

             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)

          else

             displs4_tmp => displs4_null

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,              MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,                  MPI_DOUBLE_PRECISION,&
                     &            root_rank4,  PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                  MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                      MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

       else

          if( my_rank == root_rank4 ) then

             allocate( my_displs4(0_4:comm_size-1_4) )
             jpart = lbound(sendcount4,1_4)
             my_displs4(0_4) = 0_4
             do ipart = 1,comm_size-1

                my_displs4(ipart) = my_displs4(ipart-1) + sendcount4(jpart)
                jpart = jpart + 1

             end do

          else

             allocate( my_displs4(1_4) )
             my_displs4(1_4) = 0_4

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,             MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf_tmp, recvcount4,                 MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                 MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, MPI_DOUBLE_PRECISION,&
                     &            recvbuf, recvcount4,                     MPI_DOUBLE_PRECISION,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

          deallocate( my_displs4 )

       end if

       if( istat4 /= 0 ) call runend('PAR_SCATTERV_RP_1: MPI ERROR')

    end if

#endif

  end subroutine PAR_SCATTERV_R8_2_RCV

  subroutine PAR_SCATTERV_R4_2_RCV(sendbuf,recvbuf,sendcount4,recvcount4,wherein,displs4)
    implicit none
    real(4),  pointer, intent(in)               :: sendbuf(:,:)           !< Send buffer
    real(4),  pointer, intent(inout)            :: recvbuf(:,:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)          !< Recv counts
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)             !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4),     intent(in)                  :: recvcount4
    integer(4),   pointer                       :: my_displs4(:)
    real(4)                                    :: sendbuf_tmp(2)
    real(4)                                    :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rl,sl,scl,dl,jpart
    integer(4),    pointer                      :: sendcount4_tmp(:)
    integer(4),    target                       :: sendcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       root_rank4      = 0_4
       sendcount4_null = 0_4
       displs4_null    = 0_4

       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)

       if( associated(recvbuf) ) then

          rl         = lbound(recvbuf,1_ip)

       end if

       if( associated(sendbuf)    ) sl = lbound(sendbuf,1_ip)

       if( my_rank == root_rank4 ) then

          scl = lbound(sendcount4,1_4)
          sendcount4_tmp => sendcount4(scl:)

       else

          sendcount4_tmp => sendcount4_null

       end if

       if( present(displs4) ) then

          if( associated(displs4) ) then

             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)

          else

             displs4_tmp => displs4_null

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, MPI_REAL4,&
                     &            recvbuf_tmp, recvcount4,              MPI_REAL4,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, MPI_REAL4,&
                     &            recvbuf_tmp, recvcount4,                  MPI_REAL4,&
                     &            root_rank4,  PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, displs4_tmp, MPI_REAL4,&
                     &            recvbuf, recvcount4,                  MPI_REAL4,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, displs4_tmp, MPI_REAL4,&
                     &            recvbuf, recvcount4,                      MPI_REAL4,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

       else

          if( my_rank == root_rank4 ) then

             allocate( my_displs4(0_4:comm_size-1_4) )
             jpart = lbound(sendcount4,1_4)
             my_displs4(0_4) = 0_4
             do ipart = 1,comm_size-1

                my_displs4(ipart) = my_displs4(ipart-1) + sendcount4(jpart)
                jpart = jpart + 1

             end do

          else

             allocate( my_displs4(1_4) )
             my_displs4(1_4) = 0_4

          end if

          if( recvcount4 == 0_4 ) then

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, MPI_REAL4,&
                     &            recvbuf_tmp, recvcount4,             MPI_REAL4,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, MPI_REAL4,&
                     &            recvbuf_tmp, recvcount4,                 MPI_REAL4,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          else

             if( associated(sendbuf) ) then

                call MPI_SCATTERV(sendbuf, sendcount4_tmp, my_displs4, MPI_REAL4,&
                     &            recvbuf, recvcount4,                 MPI_REAL4,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             else

                call MPI_SCATTERV(sendbuf_tmp, sendcount4_tmp, my_displs4, MPI_REAL4,&
                     &            recvbuf, recvcount4,                     MPI_REAL4,&
                     &            root_rank4, PAR_COMM_TO_USE, istat4)

             end if

          end if

          deallocate( my_displs4 )

       end if

       if( istat4 /= 0 ) call runend('PAR_SCATTERV_RP_1: MPI ERROR')

    end if

#endif

  end subroutine PAR_SCATTERV_R4_2_RCV

#ifdef __PGI
  subroutine PAR_SCATTERV_R16_2_RCV(dummy)
    implicit none
    type(pgi_dummy1), intent(in)                :: dummy
  end subroutine PAR_SCATTERV_R16_2_RCV
#else

  subroutine PAR_SCATTERV_R16_2_RCV(sendbuf,recvbuf,sendcount4,recvcount4,wherein,displs4)
    implicit none
    real(16),  pointer, intent(in)              :: sendbuf(:,:)           !< Send buffer
    real(16),  pointer, intent(inout)           :: recvbuf(:,:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)          !< Recv counts
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)             !< Displacement
    integer(4),     intent(in)                  :: recvcount4

#ifndef MPI_OFF
    if( IPARALL ) then
        call runend('PAR_BROADCAST: RP 16 IS NOT COMPATIBLE WITH MPI')
    end if

#endif

  end subroutine PAR_SCATTERV_R16_2_RCV

#endif

  !----------------------------------------------------------------------
  !
  ! RANK and SIZE of a communicator
  !
  !----------------------------------------------------------------------

  subroutine PAR_COMM_RANK_AND_SIZE_4(PAR_COMM_ORIGINAL,my_rank,comm_size)
    implicit none
    integer(4),            intent(in)  :: PAR_COMM_ORIGINAL  !< Communicator
    integer(4),            intent(out) :: my_rank            !< Rank in this communicator
    integer(4), optional,  intent(out) :: comm_size          !< Size of this communicator
    integer(4)                         :: istat4

#ifndef MPI_OFF
    istat4 = 0_4
    if( present(comm_size) ) then
       call MPI_Comm_size(PAR_COMM_ORIGINAL,comm_size,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_COMM_RANK_AND_SIZE_4: MPI ERROR 1')
    end if
    call MPI_Comm_rank(PAR_COMM_ORIGINAL,my_rank,istat4)
    if( istat4 /= 0_4 ) call runend('PAR_COMM_RANK_AND_SIZE_4: MPI ERROR 2')
#else
    my_rank = 0
    if( present(comm_size) ) comm_size = 0
#endif

  end subroutine PAR_COMM_RANK_AND_SIZE_4

  subroutine PAR_COMM_RANK_AND_SIZE_8(PAR_COMM_ORIGINAL,my_rank,comm_size)
    implicit none
    integer(8),            intent(in)  :: PAR_COMM_ORIGINAL !< Communicator
    integer(8),            intent(out) :: my_rank           !< Rank in this communicator
    integer(8), optional,  intent(out) :: comm_size         !< Size of this communicator
    integer(4)                         :: istat4
    integer(4)                         :: my_rank4
    integer(4)                         :: comm_size4
    integer(4)                         :: PAR_COMM_ORIGINAL4

#ifndef MPI_OFF
    istat4 = 0_4
    PAR_COMM_ORIGINAL4 = int(PAR_COMM_ORIGINAL,4)
    if( present(comm_size) ) then
       call MPI_Comm_size(PAR_COMM_ORIGINAL4,comm_size4,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_COMM_RANK_AND_SIZE_8: MPI ERROR 1')
       comm_size = int(comm_size4,8)
    end if
    call MPI_Comm_rank(PAR_COMM_ORIGINAL4,my_rank4,istat4)
    if( istat4 /= 0_4 ) call runend('PAR_COMM_RANK_AND_SIZE_8: MPI ERROR 2')
    my_rank = int(my_rank4,8)
#else
    my_rank = 0
    if( present(comm_size) ) comm_size = 0
#endif

  end subroutine PAR_COMM_RANK_AND_SIZE_8

  subroutine PAR_COMM_RANK_AND_SIZE_4W(my_rank,comm_size,wherein)
    implicit none
    integer(4),             intent(out) :: my_rank           !< Rank in this communicator
    integer(4),             intent(out) :: comm_size         !< Size of this communicator
    character(*),           intent(in)  :: wherein             !< Wherein
    integer(4)                          :: istat4
    integer(4)                          :: PAR_COMM_ORIGINAL

#ifndef MPI_OFF
    istat4 = 0_4
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_ORIGINAL)
    call MPI_Comm_size(PAR_COMM_ORIGINAL,comm_size,istat4)
    if( istat4 /= 0_4 ) call runend('PAR_COMM_RANK_AND_SIZE_4W: MPI ERROR 1')
    call MPI_Comm_rank(PAR_COMM_ORIGINAL,my_rank,istat4)
    if( istat4 /= 0_4 ) call runend('PAR_COMM_RANK_AND_SIZE_4W: MPI ERROR 2')
#else
    my_rank   = 0
    comm_size = 0
#endif

  end subroutine PAR_COMM_RANK_AND_SIZE_4W

  subroutine PAR_COMM_RANK_AND_SIZE_41W(my_rank,wherein)
    implicit none
    integer(4),             intent(out) :: my_rank           !< Rank in this communicator
    character(*),           intent(in)  :: wherein             !< Wherein
    integer(4)                          :: istat4
    integer(4)                          :: PAR_COMM_ORIGINAL

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_ORIGINAL)
       call MPI_Comm_rank(PAR_COMM_ORIGINAL,my_rank,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_COMM_RANK_AND_SIZE_4W: MPI ERROR 2')
#endif
    else
       my_rank = -1
    end if

  end subroutine PAR_COMM_RANK_AND_SIZE_41W

  subroutine PAR_COMM_RANK_AND_SIZE_8W(my_rank,comm_size,wherein)
    implicit none
    integer(8),   intent(out) :: my_rank           !< Rank in this communicator
    integer(8),   intent(out) :: comm_size         !< Size of this communicator
    character(*), intent(in)  :: wherein           !< Wherein
    integer(4)                :: istat4
    integer(4)                :: PAR_COMM_ORIGINAL
    integer(4)                :: my_rank4
    integer(4)                :: comm_size4

#ifndef MPI_OFF
    istat4 = 0_4
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_ORIGINAL)
    call MPI_Comm_size(PAR_COMM_ORIGINAL,comm_size4,istat4)
    if( istat4 /= 0_4 ) call runend('PAR_COMM_RANK_AND_SIZE_8W: MPI ERROR 1')
    call MPI_Comm_rank(PAR_COMM_ORIGINAL,my_rank4,istat4)
    if( istat4 /= 0_4 ) call runend('PAR_COMM_RANK_AND_SIZE_8W: MPI ERROR 2')
    my_rank = int(my_rank4,8)
    comm_size = int(comm_size4,8)
#else
    my_rank   = 0
    comm_size = 0
#endif

  end subroutine PAR_COMM_RANK_AND_SIZE_8W

  !----------------------------------------------------------------------
  !
  ! SPLIT COMMUNICATOR
  ! IKEY should have the rank
  !
  !----------------------------------------------------------------------

  subroutine PAR_COMM_SPLIT4(icolor,PAR_COMM_FINAL,my_new_rank,wherein)
    implicit none
    integer(4),   intent(in)           :: icolor
    integer(4),   intent(out)          :: PAR_COMM_FINAL
    integer(4),   intent(out)          :: my_new_rank
    character(*), intent(in), optional :: wherein
    integer(4)                         :: ikey4,istat4,jcolor4
    integer(4)                         :: PAR_COMM_INITIAL

#ifndef MPI_OFF
    istat4 = 0_4
    if( icolor == 0 ) then
       jcolor4 = MPI_UNDEFINED
    else
       jcolor4 = int(icolor,4)
    end if
    if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_INITIAL)
    else
       call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_INITIAL)
    end if
    call MPI_COMM_RANK(PAR_COMM_INITIAL,ikey4,istat4)

    if( istat4 /= 0_4 ) call runend('PAR_COMM_SPLIT4: MPI ERROR 1')
    call MPI_COMM_SPLIT(PAR_COMM_INITIAL,jcolor4,ikey4,PAR_COMM_FINAL,istat4)
    if( istat4 /= 0_4 ) call runend('PAR_COMM_SPLIT4: MPI ERROR 2')
    if( PAR_COMM_FINAL /= MPI_COMM_NULL ) then
       call MPI_COMM_RANK (PAR_COMM_FINAL,my_new_rank,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_COMM_SPLIT4: MPI ERROR 3')
    else
       my_new_rank = -1
    end if
#else
    PAR_COMM_FINAL = 0
#endif

  end subroutine PAR_COMM_SPLIT4

  subroutine PAR_COMM_SPLIT8(icolor,PAR_COMM_FINAL,my_new_rank,wherein)
    implicit none
    integer(8),   intent(in)           :: icolor
    integer(8),   intent(out)          :: PAR_COMM_FINAL
    integer(8),   intent(out)          :: my_new_rank
    character(*), intent(in), optional :: wherein
    integer(4)                         :: my_new_rank4
    integer(4)                         :: ikey4,istat4,jcolor4
    integer(4)                         :: PAR_COMM_INITIAL4
    integer(4)                         :: PAR_COMM_FINAL4

#ifndef MPI_OFF
    istat4 = 0_4
    if( icolor == 0 ) then
       jcolor4 = MPI_UNDEFINED
    else
       jcolor4 = int(icolor,4)
    end if

    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_INITIAL4)
    call MPI_COMM_RANK (PAR_COMM_INITIAL4,ikey4,istat4)

    if( istat4 /= 0_4 ) call runend('PAR_COMM_SPLIT8: MPI ERROR 1')
    call MPI_COMM_SPLIT(PAR_COMM_INITIAL4,jcolor4,ikey4,PAR_COMM_FINAL4,istat4)
    if( istat4 /= 0_4 ) call runend('PAR_COMM_SPLIT8: MPI ERROR 2')
    if( PAR_COMM_FINAL4 /= MPI_COMM_NULL ) then
       call MPI_COMM_RANK (PAR_COMM_FINAL4,my_new_rank4,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_COMM_SPLIT8: MPI ERROR 3')
       PAR_COMM_FINAL = int(PAR_COMM_FINAL4,8)
       my_new_rank    = int(my_new_rank4,8)
    else
       my_new_rank    = -1
    end if
#else
    PAR_COMM_FINAL = 0
#endif

  end subroutine PAR_COMM_SPLIT8

  !----------------------------------------------------------------------
  !
  ! BROADCAST FOR INTEGERS
  !
  !----------------------------------------------------------------------
  subroutine PAR_BROADCAST_I8_s(xx,wherein,PAR_COMM_IN4)
    implicit none
    integer(8),   intent(inout)        :: xx
    character(*),             optional :: wherein
    integer(4),   intent(in), optional :: PAR_COMM_IN4
    integer(8)                         :: n
    integer(8)                         :: yy(2)
    if( ISEQUEN ) return
    yy(1) = xx
    n     = 1
    if( present(PAR_COMM_IN4) ) then
       call PAR_BROADCAST_IP(n,yy,PAR_COMM_IN4=PAR_COMM_IN4)
    else
       call PAR_BROADCAST_IP(n,yy,wherein)
    end if
    xx    = yy(1)
  end subroutine PAR_BROADCAST_I8_s

  subroutine PAR_BROADCAST_I4_s(xx,wherein,PAR_COMM_IN4)
    implicit none
    integer(4),   intent(inout)        :: xx
    character(*),             optional :: wherein
    integer(4),   intent(in), optional :: PAR_COMM_IN4
    integer(4)                         :: n
    integer(4)                         :: yy(2)
    if( ISEQUEN ) return
    yy(1) = xx
    n     = 1
    if( present(PAR_COMM_IN4) ) then
       call PAR_BROADCAST_IP(n,yy,PAR_COMM_IN4=PAR_COMM_IN4)
    else
       call PAR_BROADCAST_IP(n,yy,wherein)
    end if
    xx    = yy(1)
  end subroutine PAR_BROADCAST_I4_s

  subroutine PAR_BROADCAST_I4_1(xx,wherein,PAR_COMM_IN4)

    integer(4),  pointer, intent(inout)         :: xx(:)
    character(*),                      optional :: wherein
    integer(4),            intent(in), optional :: PAR_COMM_IN4
    integer(4)                                  :: n

    if( ISEQUEN ) return

    if( associated(xx) ) then
       n = size(xx)
       if( n > 0 ) then
          if(     present(PAR_COMM_IN4) ) then
             call PAR_BROADCAST_IP(n,xx,PAR_COMM_IN4=PAR_COMM_IN4)
          else if( present(wherein) ) then
             call PAR_BROADCAST_IP(n,xx,wherein)
          else
             call PAR_BROADCAST_I4(n,xx,'IN MY CODE')
          end if
       end if
    end if
  end subroutine PAR_BROADCAST_I4_1

  subroutine PAR_BROADCAST_I8_1(xx,wherein,PAR_COMM_IN4)

    integer(8),   pointer, intent(inout)        :: xx(:)
    character(*),                      optional :: wherein
    integer(4),            intent(in), optional :: PAR_COMM_IN4
    integer(8)                                  :: n
    if( ISEQUEN ) return

    if( associated(xx) ) then
       n = size(xx,KIND=8)
       if( n > 0 ) then
          if(     present(PAR_COMM_IN4) ) then
             call PAR_BROADCAST_IP(n,xx,PAR_COMM_IN4=PAR_COMM_IN4)
          else if( present(wherein) ) then
             call PAR_BROADCAST_IP(n,xx,wherein)
          else
             call PAR_BROADCAST_IP(n,xx,'IN MY CODE')
          end if
       end if
    end if
  end subroutine PAR_BROADCAST_I8_1

  subroutine PAR_BROADCAST_IP_04(n,xx,wherein,PAR_COMM_IN4)

    integer(4),          intent(in)           :: n
    integer(4),          intent(inout)        :: xx(*)
    character(*),                    optional :: wherein
    integer(4),          intent(in), optional :: PAR_COMM_IN4

    if( ISEQUEN ) return

    if( n > 0 ) then
       if(      present(PAR_COMM_IN4) ) then
          call PAR_BROADCAST_IP(n,xx,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( present(wherein) ) then
          call PAR_BROADCAST_IP(n,xx,wherein)
       else
          call PAR_BROADCAST_IP(n,xx,'IN MY CODE')
       end if
    end if

  end subroutine PAR_BROADCAST_IP_04

  subroutine PAR_BROADCAST_IP_08(n,xx,wherein,PAR_COMM_IN4)
    implicit none
    integer(8),          intent(in)           :: n
    integer(8),          intent(inout)        :: xx(*)
    character(*),                    optional :: wherein
    integer(4),          intent(in), optional :: PAR_COMM_IN4

    if( ISEQUEN .or. n <= 0 ) then
       return
    else
       if(      present(PAR_COMM_IN4) ) then
          call PAR_BROADCAST_IP(n,xx,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( present(wherein) ) then
          call PAR_BROADCAST_IP(n,xx,wherein)
       else
          call PAR_BROADCAST_IP(n,xx,'IN MY CODE')
       end if
    end if

  end subroutine PAR_BROADCAST_IP_08

  subroutine PAR_BROADCAST_I4(n,xx,wherein,root_rank4,PAR_COMM_IN4)
    implicit none
    integer(4),   intent(in)           :: n
    integer(4),   intent(inout)        :: xx(n)
    character(*),             optional :: wherein
    integer(4),   intent(in), optional :: PAR_COMM_IN4
    integer(4),   intent(in), optional :: root_rank4
    integer(4)                         :: istat4,n4,PAR_COMM_TO_USE
    if( ISEQUEN ) return

    if( n > 0 ) then
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       n4 = int(n,4)
#ifndef MPI_OFF
       istat4 = 0_4
       if( present(root_rank4) ) then
          call MPI_Bcast(xx(1:n),n4,MPI_INTEGER,root_rank4,PAR_COMM_TO_USE,istat4)
       else
          call MPI_Bcast(xx(1:n),n4,MPI_INTEGER,0_4,PAR_COMM_TO_USE,istat4)
       end if
       if( istat4 /= 0_4 ) call runend('PAR_BROADCAST_IP: MPI ERROR')
#endif
    end if
  end subroutine PAR_BROADCAST_I4

  subroutine PAR_BROADCAST_I8(n,xx,wherein,root_rank4,PAR_COMM_IN4)
    implicit none
    integer(8),  intent(in)            :: n
    integer(8),  intent(inout)         :: xx(n)
    character(*),             optional :: wherein
    integer(4),   intent(in), optional :: root_rank4
    integer(4),   intent(in), optional :: PAR_COMM_IN4
    integer(4)                         :: istat4,n4,PAR_COMM_TO_USE
    if( ISEQUEN ) return

    if( n > 0 ) then
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       n4 = int(n,4)
#ifndef MPI_OFF
       istat4 = 0_4
       if( present(root_rank4) ) then
          call MPI_Bcast(xx(1:n),n4,MPI_INTEGER8,root_rank4,PAR_COMM_TO_USE,istat4)
       else
          call MPI_Bcast(xx(1:n),n4,MPI_INTEGER8,0_4,PAR_COMM_TO_USE,istat4)
       end if
       if( istat4 /= 0_4 ) call runend('PAR_BROADCAST_I8: MPI ERROR')
#endif
    end if
  end subroutine PAR_BROADCAST_I8


  subroutine PAR_BROADCAST_R8_s(xx,wherein)
    implicit none
    real(8),     intent(inout) :: xx
    character(*), optional      :: wherein
    integer(ip)                 :: n
    real(8)                    :: yy(2)
    yy(1) = xx
    n     = 1
    call PAR_BROADCAST_R8(n,yy,wherein)
    xx    = yy(1)
  end subroutine PAR_BROADCAST_R8_s

  subroutine PAR_BROADCAST_R8_1(xx,wherein)
    implicit none
    real(8),     pointer, intent(inout) :: xx(:)
    character(*),          optional      :: wherein
    integer(ip)                          :: n

    if( ISEQUEN ) return
    if( associated(xx) ) then
       n = size(xx)
       if( n > 0 ) then
          if( present(wherein) ) then
             call PAR_BROADCAST_R8(n,xx,wherein)
          else
             call PAR_BROADCAST_R8(n,xx,'IN MY CODE')
          end if
       end if
    end if

  end subroutine PAR_BROADCAST_R8_1

  subroutine PAR_BROADCAST_R8_0(n,xx,wherein,root_rank)
    implicit none
    integer(ip),  intent(in)           :: n
    real(8),     intent(inout)        :: xx(n)
    character(*), intent(in), optional :: wherein
    integer(ip),  intent(in), optional :: root_rank

    if( present(wherein) ) then
       if( present(root_rank) ) then
          call PAR_BROADCAST_R8(n,xx,wherein,root_rank)
       else
          call PAR_BROADCAST_R8(n,xx,wherein)
       end if
    else
       if( present(root_rank) ) then
          call PAR_BROADCAST_R8(n,xx,'IN MY CODE',root_rank)
       else
          call PAR_BROADCAST_R8(n,xx,'IN MY CODE')
       end if
    end if

  end subroutine PAR_BROADCAST_R8_0

  subroutine PAR_BROADCAST_R8(n,xx,wherein,root_rank)
    implicit none
    integer(ip),  intent(in)           :: n
    real(8),     intent(inout)        :: xx(n)
    character(*), intent(in), optional :: wherein
    integer(ip),  intent(in), optional :: root_rank
    integer(4)                         :: istat4,n4,PAR_COMM_TO_USE,root_rank4

    if( ISEQUEN ) return
    if( n > 0 ) then
       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       n4 = int(n,4)
       if( present(root_rank) ) then
          root_rank4 = int(root_rank,4)
       else
          root_rank4 = 0_4
       end if
#ifndef MPI_OFF
       istat4 = 0_4
       call MPI_Bcast(xx(1:n),n4,MPI_DOUBLE_PRECISION,root_rank4,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_BROADCAST_RP: MPI ERROR')
#endif
    end if
  end subroutine PAR_BROADCAST_R8

  !--4

    subroutine PAR_BROADCAST_R4_s(xx,wherein)
    implicit none
    real(4),     intent(inout) :: xx
    character(*), optional      :: wherein
    integer(ip)                 :: n
    real(4)                    :: yy(2)
    yy(1) = xx
    n     = 1
    call PAR_BROADCAST_R4(n,yy,wherein)
    xx    = yy(1)
  end subroutine PAR_BROADCAST_R4_s

  subroutine PAR_BROADCAST_R4_1(xx,wherein)
    implicit none
    real(4),     pointer, intent(inout) :: xx(:)
    character(*),          optional      :: wherein
    integer(ip)                          :: n

    if( ISEQUEN ) return
    if( associated(xx) ) then
       n = size(xx)
       if( n > 0 ) then
          if( present(wherein) ) then
             call PAR_BROADCAST_R4(n,xx,wherein)
          else
             call PAR_BROADCAST_R4(n,xx,'IN MY CODE')
          end if
       end if
    end if

  end subroutine PAR_BROADCAST_R4_1

  subroutine PAR_BROADCAST_R4_0(n,xx,wherein,root_rank)
    implicit none
    integer(ip),  intent(in)           :: n
    real(4),     intent(inout)        :: xx(n)
    character(*), intent(in), optional :: wherein
    integer(ip),  intent(in), optional :: root_rank

    if( present(wherein) ) then
       if( present(root_rank) ) then
          call PAR_BROADCAST_R4(n,xx,wherein,root_rank)
       else
          call PAR_BROADCAST_R4(n,xx,wherein)
       end if
    else
       if( present(root_rank) ) then
          call PAR_BROADCAST_R4(n,xx,'IN MY CODE',root_rank)
       else
          call PAR_BROADCAST_R4(n,xx,'IN MY CODE')
       end if
    end if

  end subroutine PAR_BROADCAST_R4_0

  subroutine PAR_BROADCAST_R4(n,xx,wherein,root_rank)
    implicit none
    integer(ip),  intent(in)           :: n
    real(4),     intent(inout)        :: xx(n)
    character(*), intent(in), optional :: wherein
    integer(ip),  intent(in), optional :: root_rank
    integer(4)                         :: istat4,n4,PAR_COMM_TO_USE,root_rank4

    if( ISEQUEN ) return
    if( n > 0 ) then
       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       n4 = int(n,4)
       if( present(root_rank) ) then
          root_rank4 = int(root_rank,4)
       else
          root_rank4 = 0_4
       end if
#ifndef MPI_OFF
       istat4 = 0_4
       call MPI_Bcast(xx(1:n),n4,MPI_REAL4,root_rank4,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_BROADCAST_RP: MPI ERROR')
#endif
    end if
  end subroutine PAR_BROADCAST_R4

    !--16

#ifdef __PGI
  subroutine PAR_BROADCAST_R16_s(dummy)
    implicit none
    type(pgi_dummy1), intent(in)                :: dummy
  end subroutine

  subroutine PAR_BROADCAST_R16_1(dummy)
    implicit none
    type(pgi_dummy2), intent(in)                :: dummy
  end subroutine

  subroutine PAR_BROADCAST_R16_0(dummy)
    implicit none
    type(pgi_dummy3), intent(in)                :: dummy
  end subroutine

  subroutine PAR_BROADCAST_R16(dummy)
    implicit none
    type(pgi_dummy5), intent(in)                :: dummy
  end subroutine

#else

  subroutine PAR_BROADCAST_R16_s(xx,wherein)
    implicit none
    real(16),     intent(inout) :: xx
    character(*), optional      :: wherein
    integer(ip)                 :: n
    real(16)                    :: yy(2)
    yy(1) = xx
    n     = 1
    call PAR_BROADCAST_R16(n,yy,wherein)
    xx    = yy(1)
  end subroutine PAR_BROADCAST_R16_s

  subroutine PAR_BROADCAST_R16_1(xx,wherein)
    implicit none
    real(16),     pointer, intent(inout) :: xx(:)
    character(*),          optional      :: wherein
    integer(ip)                          :: n

    if( ISEQUEN ) return
    if( associated(xx) ) then
       n = size(xx)
       if( n > 0 ) then
          if( present(wherein) ) then
             call PAR_BROADCAST_R16(n,xx,wherein)
          else
             call PAR_BROADCAST_R16(n,xx,'IN MY CODE')
          end if
       end if
    end if

  end subroutine PAR_BROADCAST_R16_1

  subroutine PAR_BROADCAST_R16_0(n,xx,wherein,root_rank)
    implicit none
    integer(ip),  intent(in)           :: n
    real(16),     intent(inout)        :: xx(n)
    character(*), intent(in), optional :: wherein
    integer(ip),  intent(in), optional :: root_rank

    if( present(wherein) ) then
       if( present(root_rank) ) then
          call PAR_BROADCAST_R16(n,xx,wherein,root_rank)
       else
          call PAR_BROADCAST_R16(n,xx,wherein)
       end if
    else
       if( present(root_rank) ) then
          call PAR_BROADCAST_R16(n,xx,'IN MY CODE',root_rank)
       else
          call PAR_BROADCAST_R16(n,xx,'IN MY CODE')
       end if
    end if

  end subroutine PAR_BROADCAST_R16_0

  subroutine PAR_BROADCAST_R16(n,xx,wherein,root_rank)
    implicit none
    integer(ip),  intent(in)           :: n
    real(16),     intent(inout)        :: xx(n)
    character(*), intent(in), optional :: wherein
    integer(ip),  intent(in), optional :: root_rank
    integer(4)                         :: istat4,n4,PAR_COMM_TO_USE,root_rank4

    if( ISEQUEN ) return
    if( n > 0 ) then
       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       end if
       n4 = int(n,4)
       if( present(root_rank) ) then
          root_rank4 = int(root_rank,4)
       else
          root_rank4 = 0_4
       end if
#ifndef MPI_OFF
       istat4 = 0_4
       call runend('PAR_BROADCAST: RP 16 IS NOT COMPATIBLE WITH MPI')
       !call MPI_Bcast(xx(1:n),n4,MPI_REAL4,root_rank4,PAR_COMM_TO_USE,istat4)
       !if( istat4 /= 0_4 ) call runend('PAR_BROADCAST_RP: MPI ERROR')
#endif
    end if
  end subroutine PAR_BROADCAST_R16

#endif

  subroutine PAR_BROADCAST_CH(n,xx,wherein)
    implicit none
    integer(ip),  intent(in)           :: n
    character(*), intent(inout)        :: xx
    character(*), intent(in)           :: wherein
    integer(4)                         :: istat4,n4,PAR_COMM_TO_USE

    if( IPARALL .and. n > 0 ) then
#ifndef MPI_OFF
       !if( present(wherein) ) then
       !   call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       !else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE)
       !end if
       istat4 = 0_4
       n4 = int(n,4)
       call MPI_Bcast(xx(1:n),n4,MPI_CHARACTER,0_4,PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_BROADCAST_CP: MPI ERROR')
#endif
    end if

  end subroutine PAR_BROADCAST_CH

  subroutine PAR_BROADCAST_LG_s(xx,wherein)
    implicit none
    logical(lg),  intent(inout) :: xx
    character(*), optional      :: wherein
    integer(4)                  :: istat4,PAR_COMM_TO_USE4

    if( IPARALL ) then
#ifndef MPI_OFF
       istat4 = 0_4
       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE4)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE4)
       end if
       call MPI_Bcast(xx,1_4,MPI_LOGICAL,0_4,PAR_COMM_TO_USE4,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_BROADCAST_CP: MPI ERROR')
#endif
    end if
  end subroutine PAR_BROADCAST_LG_s

  subroutine PAR_BROADCAST_LG_1(xx,wherein)
    implicit none
    logical(lg),  pointer, intent(inout) :: xx(:)
    character(*), optional               :: wherein
    integer(4)                           :: istat4,PAR_COMM_TO_USE4,n4

    n4 = int(memory_size(xx),4_4)
    if( IPARALL .and. n4 > 0 ) then
#ifndef MPI_OFF
       istat4 = 0_4
       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE4)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE4)
       end if
       call MPI_Bcast(xx,n4,MPI_LOGICAL,0_4,PAR_COMM_TO_USE4,istat4)
       if( istat4 /= 0_4 ) call runend('PAR_BROADCAST_CP: MPI ERROR')
#endif
    end if
  end subroutine PAR_BROADCAST_LG_1

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_SEND and PAR_RECEIVE
  !
  !----------------------------------------------------------------------

  subroutine PAR_RECEIVE_IP_s(xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN4)
    implicit none
    integer(ip),           intent(out)          :: xx_recv
    character(*),          intent(in), optional :: wherein
    integer(ip),           intent(in), optional :: dom_i
    character(*),          intent(in), optional :: wsynch
    integer(4),            intent(in), optional :: PAR_COMM_IN4
    integer(ip)                                 :: yy_send(2)
    integer(ip)                                 :: yy_recv(2)

    if( ISEQUEN ) return
    if( present(wsynch) ) then
       call PAR_SEND_RECEIVE_IP(0_ip,1_ip,yy_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN4=PAR_COMM_IN4)
    else
       call PAR_SEND_RECEIVE_IP(0_ip,1_ip,yy_send,yy_recv,wherein,dom_i,PAR_COMM_IN4=PAR_COMM_IN4)
    end if
    xx_recv    = yy_recv(1)

  end subroutine PAR_RECEIVE_IP_s

  subroutine PAR_RECEIVE_RP_s(xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN4)
    implicit none
    real(rp),              intent(out)          :: xx_recv
    character(*),          intent(in), optional :: wherein
    integer(ip),           intent(in), optional :: dom_i
    character(*),          intent(in), optional :: wsynch
    integer(4),            intent(in), optional :: PAR_COMM_IN4
    real(rp)                                    :: yy_send(2)
    real(rp)                                    :: yy_recv(2)

    if( ISEQUEN ) return
    if( present(wsynch) ) then
       call PAR_SEND_RECEIVE_RP(0_ip,1_ip,yy_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN4=PAR_COMM_IN4)
    else
       call PAR_SEND_RECEIVE_RP(0_ip,1_ip,yy_send,yy_recv,wherein,dom_i,PAR_COMM_IN4=PAR_COMM_IN4)
    end if
    xx_recv    = yy_recv(1)

  end subroutine PAR_RECEIVE_RP_s

  subroutine PAR_SEND_IP_s(xx_send,wherein,dom_i,wsynch,PAR_COMM_IN4)
    implicit none
    integer(ip),           intent(in)           :: xx_send
    character(*),          intent(in), optional :: wherein
    integer(ip),           intent(in), optional :: dom_i
    character(*),          intent(in), optional :: wsynch
    integer(4),            intent(in), optional :: PAR_COMM_IN4
    integer(ip)                                 :: yy_send(2)
    integer(ip)                                 :: yy_recv(2)

    if( ISEQUEN ) return
    yy_send(1) = xx_send
    if( present(wsynch) ) then
       call PAR_SEND_RECEIVE_IP(1_ip,0_ip,yy_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN4=PAR_COMM_IN4)
    else
       call PAR_SEND_RECEIVE_IP(1_ip,0_ip,yy_send,yy_recv,wherein,dom_i,PAR_COMM_IN4=PAR_COMM_IN4)
    end if

  end subroutine PAR_SEND_IP_s

  subroutine PAR_SEND_RP_s(xx_send,wherein,dom_i,wsynch,PAR_COMM_IN4)
    implicit none
    real(rp),              intent(in)           :: xx_send
    character(*),          intent(in), optional :: wherein
    integer(ip),           intent(in), optional :: dom_i
    character(*),          intent(in), optional :: wsynch
    integer(4),            intent(in), optional :: PAR_COMM_IN4
    real(rp)                                    :: yy_send(2)
    real(rp)                                    :: yy_recv(2)

    if( ISEQUEN ) return
    yy_send(1) = xx_send
    if( present(wsynch) ) then
       call PAR_SEND_RECEIVE_RP(1_ip,0_ip,yy_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN4=PAR_COMM_IN4)
    else
       call PAR_SEND_RECEIVE_RP(1_ip,0_ip,yy_send,yy_recv,wherein,dom_i,PAR_COMM_IN4=PAR_COMM_IN4)
    end if

  end subroutine PAR_SEND_RP_s
  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_SEND_RECEIVE_IP
  !
  !----------------------------------------------------------------------

  subroutine PAR_SEND_RECEIVE_IP_s(xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN4)
    implicit none
    integer(ip),           intent(in)           :: xx_send
    integer(ip),           intent(out)          :: xx_recv
    character(*),          intent(in), optional :: wherein
    integer(ip),           intent(in), optional :: dom_i
    character(*),          intent(in), optional :: wsynch
    integer(4),            intent(in), optional :: PAR_COMM_IN4
    integer(ip)                                 :: yy_send(2)
    integer(ip)                                 :: yy_recv(2)
    if( ISEQUEN ) return
    yy_send(1) = xx_send
    if( present(wsynch) ) then
       call PAR_SEND_RECEIVE_IP(1_ip,1_ip,yy_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN4=PAR_COMM_IN4)
    else
       call PAR_SEND_RECEIVE_IP(1_ip,1_ip,yy_send,yy_recv,wherein,dom_i,PAR_COMM_IN4=PAR_COMM_IN4)
    end if
    xx_recv    = yy_recv(1)
  end subroutine PAR_SEND_RECEIVE_IP_s
  subroutine PAR_SEND_RECEIVE_IP_0(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN4)
    implicit none
    integer(ip),           intent(in)           :: nsend
    integer(ip),           intent(in)           :: nrecv
    integer(ip),           intent(in)           :: xx_send(*)
    integer(ip),           intent(out)          :: xx_recv(*)
    character(*),          intent(in), optional :: wherein
    integer(ip),           intent(in)           :: dom_i
    character(*),          intent(in), optional :: wsynch
    integer(4),            intent(in), optional :: PAR_COMM_IN4
    integer(ip)                                 :: yy_send(2)
    integer(ip)                                 :: yy_recv(2)
    if( ISEQUEN ) return
    if( present(wsynch) ) then
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN4=PAR_COMM_IN4)
       else
          return
       end if
    else
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,PAR_COMM_IN4=PAR_COMM_IN4)
       else
          return
       end if
    end if
  end subroutine PAR_SEND_RECEIVE_IP_0
  subroutine PAR_SEND_RECEIVE_IP_1(xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN4)
    implicit none
    integer(ip),  pointer, intent(in)           :: xx_send(:)
    integer(ip),  pointer, intent(inout)          :: xx_recv(:)
    character(*),          intent(in), optional :: wherein
    integer(ip),           intent(in), optional :: dom_i
    character(*),          intent(in), optional :: wsynch
    integer(4),            intent(in), optional :: PAR_COMM_IN4
    integer(ip)                                 :: nsend,nrecv
    integer(ip)                                 :: yy_send(2)
    integer(ip)                                 :: yy_recv(2)
    if( ISEQUEN ) return
    if( .not. associated(xx_send) ) then
       nsend = 0
    else
       nsend = size(xx_send)
    end if
    if( .not. associated(xx_recv) ) then
       nrecv = 0
    else
       nrecv = size(xx_recv)
    end if
    if( present(wsynch) ) then
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN4=PAR_COMM_IN4)
       else
          return
       end if
    else
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,PAR_COMM_IN4=PAR_COMM_IN4)
       else
          return
       end if
    end if
  end subroutine PAR_SEND_RECEIVE_IP_1
  subroutine PAR_SEND_RECEIVE_IP_2(xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN4)
    implicit none
    integer(ip),  pointer, intent(in)           :: xx_send(:,:)
    integer(ip),  pointer, intent(inout)          :: xx_recv(:,:)
    character(*),          intent(in), optional :: wherein
    integer(ip),           intent(in), optional :: dom_i
    character(*),          intent(in), optional :: wsynch
    integer(4),            intent(in), optional :: PAR_COMM_IN4
    integer(ip)                                 :: nsend,nrecv
    integer(ip)                                 :: yy_send(2)
    integer(ip)                                 :: yy_recv(2)

    if( ISEQUEN ) return

    if( .not. associated(xx_send) ) then
       nsend = 0
    else
       nsend = size(xx_send,1)*size(xx_send,2)
    end if
    if( .not. associated(xx_recv) ) then
       nrecv = 0
    else
       nrecv = size(xx_recv,1)*size(xx_recv,2)
    end if
    if( present(wsynch) ) then
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN4=PAR_COMM_IN4)
       else
          return
       end if
    else
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,PAR_COMM_IN4=PAR_COMM_IN4)
       else
          return
       end if
    end if
  end subroutine PAR_SEND_RECEIVE_IP_2
  subroutine PAR_SEND_RECEIVE_IP_3(xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN4)
    implicit none
    integer(ip),  pointer, intent(in)           :: xx_send(:,:,:)
    integer(ip),  pointer, intent(inout)          :: xx_recv(:,:,:)
    character(*),          intent(in), optional :: wherein
    integer(ip),           intent(in), optional :: dom_i
    character(*),          intent(in), optional :: wsynch
    integer(4),            intent(in), optional :: PAR_COMM_IN4
    integer(ip)                                 :: nsend,nrecv
    integer(ip)                                 :: yy_send(2)
    integer(ip)                                 :: yy_recv(2)
    if( ISEQUEN ) return
    if( .not. associated(xx_send) ) then
       nsend = 0
    else
       nsend = size(xx_send,1)*size(xx_send,2)*size(xx_send,3)
    end if
    if( .not. associated(xx_recv) ) then
       nrecv = 0
    else
       nrecv = size(xx_recv,1)*size(xx_recv,2)*size(xx_recv,3)
    end if
    if( present(wsynch) ) then
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN4=PAR_COMM_IN4)
       else
          return
       end if
    else
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,PAR_COMM_IN4=PAR_COMM_IN4)
       else
          return
       end if
    end if
  end subroutine PAR_SEND_RECEIVE_IP_3

  !----------------------------------------------------------------------
  !
  ! PAR_SEND_RECEIVE_IP
  !
  !----------------------------------------------------------------------

  subroutine PAR_SEND_RECEIVE_IP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN4)
    implicit none
    integer(ip),            intent(in)  :: nsend
    integer(ip),            intent(in)  :: nrecv
    integer(ip),            intent(in)  :: xx_send(*)
    integer(ip),            intent(out) :: xx_recv(*)
    character(*), optional, intent(in)  :: wherein
    integer(ip),            intent(in)  :: dom_i
    character(*), optional, intent(in)  :: wsynch
    integer(4),   optional, intent(in)  :: PAR_COMM_IN4
    integer(ip)                         :: kk
    integer(4)                          :: istat4,nsend4,nrecv4,dom_i4
    integer(4)                          :: PAR_COMM_TO_USE
    logical(lg)                         :: asynch
    !
    ! Define communicator
    !
    if( present(PAR_COMM_IN4) ) then
       PAR_COMM_TO_USE = PAR_COMM_IN4
    else if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
    else
       PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
    end if
    !
    ! BlockinG/non blocking
    !
    if( present(wsynch) ) then
       if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
          asynch = .false.
       else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
          asynch = .true.
          if( .not. associated(non_blocking(inonblocking) % request4) ) then
             call runend('NON-BLOCKING SEND/RECEIVE SHOULD BE STARTED')
          end if
       else
          call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
       end if
    else
       asynch = .false.
    end if
    !
    ! Synchronous Send/receive
    !
    nsend4 = int(nsend,4)
    nrecv4 = int(nrecv,4)
    dom_i4 = int(dom_i,4)
    istat4 = 0

#ifndef MPI_OFF
    if( asynch ) then
       if( nrecv /= 0 ) then
          non_blocking(inonblocking) % count4 = non_blocking(inonblocking) % count4 + 1
          kk = non_blocking(inonblocking) % count4
          call MPI_IRecv(                                 &
               xx_recv(1:nrecv), nrecv4,                  &
               PAR_INTEGER, dom_i4, 0_4,                  &
               PAR_COMM_TO_USE,                           &
               non_blocking(inonblocking) % request4(kk), &
               istat4                                     )
       end if
       if( nsend /= 0 ) then
          non_blocking(inonblocking) % count4 = non_blocking(inonblocking) % count4 + 1
          kk = non_blocking(inonblocking) % count4
          call MPI_ISend(                                 &
               xx_send(1:nsend), nsend4,                  &
               PAR_INTEGER, dom_i4, 0_4,                  &
               PAR_COMM_TO_USE,                           &
               non_blocking(inonblocking) % request4(kk), &
               istat4                                     )
       end if
    else
       if( nrecv /= 0 .and. nsend == 0 ) then
          call MPI_Recv(                          &
               xx_recv(1:nrecv), nrecv4,          &
               PAR_INTEGER, dom_i4, 0_4,          &
               PAR_COMM_TO_USE, status, istat4    )

       else if( nrecv == 0 .and. nsend /= 0 ) then
          call MPI_Send(                          &
               xx_send(1:nsend), nsend4,          &
               PAR_INTEGER, dom_i4, 0_4,          &
               PAR_COMM_TO_USE, istat4            )

       else if( nrecv /= 0 .and. nsend /= 0 ) then
          call MPI_Sendrecv(                      &
               xx_send(1:nsend), nsend4,          &
               PAR_INTEGER, dom_i4, 0_4,          &
               xx_recv(1:nrecv), nrecv4,          &
               PAR_INTEGER, dom_i4, 0_4,          &
               PAR_COMM_TO_USE, status, istat4    )

       end if
    end if

#endif
    if( istat4 /= 0_4 ) call runend('PAR_SEND_RECEIVE_IP: MPI ERROR')

  end subroutine PAR_SEND_RECEIVE_IP

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_SEND_RECEIVE_RP
  !
  !----------------------------------------------------------------------

  subroutine PAR_SEND_RECEIVE_RP_s(xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN4)
    implicit none
    real(rp),            intent(in)           :: xx_send
    real(rp),            intent(out)          :: xx_recv
    character(*),        intent(in)           :: wherein
    integer(ip),         intent(in)           :: dom_i
    character(*),        intent(in), optional :: wsynch
    integer(4),          intent(in), optional :: PAR_COMM_IN4
    real(rp)                                  :: yy_send(2)
    real(rp)                                  :: yy_recv(2)
    if( ISEQUEN ) return
    yy_send(1) = xx_send
    if( present(wsynch) ) then
       call PAR_SEND_RECEIVE_RP(1_ip,1_ip,yy_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN4=PAR_COMM_IN4)
    else
       call PAR_SEND_RECEIVE_RP(1_ip,1_ip,yy_send,yy_recv,wherein,dom_i,PAR_COMM_IN4=PAR_COMM_IN4)
    end if
    xx_recv    = yy_recv(1)
  end subroutine PAR_SEND_RECEIVE_RP_s
  subroutine PAR_SEND_RECEIVE_RP_0(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN4)
    implicit none
    integer(ip),         intent(in)           :: nsend
    integer(ip),         intent(in)           :: nrecv
    real(rp),            intent(in)           :: xx_send(*)
    real(rp),            intent(out)          :: xx_recv(*)
    character(*),        intent(in)           :: wherein
    integer(ip),         intent(in)           :: dom_i
    character(*),        intent(in), optional :: wsynch
    integer(4),          intent(in), optional :: PAR_COMM_IN4
    real(rp)                                  :: yy_send(2)
    real(rp)                                  :: yy_recv(2)
    if( ISEQUEN ) return
    if( present(wsynch) ) then
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN4=PAR_COMM_IN4)
       else
          return
       end if
    else
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,PAR_COMM_IN4=PAR_COMM_IN4)
       else
          return
       end if
    end if
  end subroutine PAR_SEND_RECEIVE_RP_0
  subroutine PAR_SEND_RECEIVE_RP_1(xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN4)
    implicit none
    real(rp),     pointer, intent(in)           :: xx_send(:)
    real(rp),     pointer, intent(inout)        :: xx_recv(:)
    character(*),          intent(in)           :: wherein
    integer(ip),           intent(in)           :: dom_i
    character(*),          intent(in), optional :: wsynch
    integer(4),            intent(in), optional :: PAR_COMM_IN4
    integer(ip)                                 :: nsend,nrecv
    real(rp)                                    :: yy_send(2)
    real(rp)                                    :: yy_recv(2)
    if( ISEQUEN ) return
    if( .not. associated(xx_send) ) then
       nsend = 0
    else
       nsend = size(xx_send)
    end if
    if( .not. associated(xx_recv) ) then
       nrecv = 0
    else
       nrecv = size(xx_recv)
    end if
    if( present(wsynch) ) then
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN4=PAR_COMM_IN4)
       else
          return
       end if
    else
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,PAR_COMM_IN4=PAR_COMM_IN4)
       else
          return
       end if
    end if
  end subroutine PAR_SEND_RECEIVE_RP_1
  subroutine PAR_SEND_RECEIVE_RP_2(xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN4)
    implicit none
    real(rp),     pointer, intent(in)           :: xx_send(:,:)
    real(rp),     pointer, intent(inout)        :: xx_recv(:,:)
    character(*),          intent(in)           :: wherein
    integer(ip),           intent(in)           :: dom_i
    character(*),          intent(in), optional :: wsynch
    integer(4),            intent(in), optional :: PAR_COMM_IN4
    integer(ip)                                 :: nsend,nrecv
    real(rp)                                    :: yy_send(2)
    real(rp)                                    :: yy_recv(2)
    if( ISEQUEN ) return
    if( .not. associated(xx_send) ) then
       nsend = 0
    else
       nsend = size(xx_send,1)*size(xx_send,2)
    end if
    if( .not. associated(xx_recv) ) then
       nrecv = 0
    else
       nrecv = size(xx_recv,1)*size(xx_recv,2)
    end if
    if( present(wsynch) ) then
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN4=PAR_COMM_IN4)
       else
          return
       end if
    else
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,PAR_COMM_IN4=PAR_COMM_IN4)
       else
          return
       end if
    end if
  end subroutine PAR_SEND_RECEIVE_RP_2
  subroutine PAR_SEND_RECEIVE_RP_3(xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN4)
    implicit none
    real(rp),     pointer, intent(in)           :: xx_send(:,:,:)
    real(rp),     pointer, intent(inout)        :: xx_recv(:,:,:)
    character(*),          intent(in)           :: wherein
    integer(ip),           intent(in)           :: dom_i
    character(*),          intent(in), optional :: wsynch
    integer(4),            intent(in), optional :: PAR_COMM_IN4
    integer(ip)                                 :: nsend,nrecv
    real(rp)                                    :: yy_send(2)
    real(rp)                                    :: yy_recv(2)
    if( ISEQUEN ) return
    if( .not. associated(xx_send) ) then
       nsend = 0
    else
       nsend = size(xx_send,1)*size(xx_send,2)*size(xx_send,3)
    end if
    if( .not. associated(xx_recv) ) then
       nrecv = 0
    else
       nrecv = size(xx_recv,1)*size(xx_recv,2)*size(xx_recv,3)
    end if
    if( present(wsynch) ) then
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,wsynch,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN4=PAR_COMM_IN4)
       else
          return
       end if
    else
       if( nsend == 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,yy_send,xx_recv,wherein,dom_i,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv == 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,yy_recv,wherein,dom_i,PAR_COMM_IN4=PAR_COMM_IN4)
       else if( nsend /= 0 .and. nrecv /= 0 ) then
          call PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,PAR_COMM_IN4=PAR_COMM_IN4)
       else
          return
       end if
    end if
  end subroutine PAR_SEND_RECEIVE_RP_3

  !----------------------------------------------------------------------
  !
  ! PAR_SEND_RECEIVE_RP
  !
  !----------------------------------------------------------------------

  subroutine PAR_SEND_RECEIVE_RP(nsend,nrecv,xx_send,xx_recv,wherein,dom_i,wsynch,PAR_COMM_IN4)
    implicit none
    integer(ip),            intent(in)  :: nsend
    integer(ip),            intent(in)  :: nrecv
    real(rp),               intent(in)  :: xx_send(*)
    real(rp),               intent(out) :: xx_recv(*)
    character(*),           intent(in)  :: wherein
    integer(ip),            intent(in)  :: dom_i
    character(*), optional, intent(in)  :: wsynch
    integer(4),   optional, intent(in)  :: PAR_COMM_IN4
    integer(ip)                         :: kk
    integer(4)                          :: istat4,nsend4,nrecv4,dom_i4
    integer(4)                          :: PAR_COMM_TO_USE
    logical(lg)                         :: asynch

    if( IPARALL ) then
       !
       ! Define communicator
       !
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       end if
       !
       ! BlockinG/non blocking
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
             if( .not. associated(non_blocking(inonblocking) % request4) ) then
                call runend('NON-BLOCKING SEND/RECEIVE SHOULD BE STARTED')
             end if
          else
             call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if
       !
       ! Synchronous Send/receive
       !
       nsend4 = int(nsend,4)
       nrecv4 = int(nrecv,4)
       dom_i4 = int(dom_i,4)
       istat4 = 0_4

#ifndef MPI_OFF
       if( asynch ) then
          if( nrecv /= 0 ) then
             non_blocking(inonblocking) % count4 = non_blocking(inonblocking) % count4 + 1
             kk = non_blocking(inonblocking) % count4
             call MPI_IRecv(                                 &
                  xx_recv(1:nrecv), nrecv4,                  &
                  MPI_DOUBLE_PRECISION, dom_i4, 0_4,         &
                  PAR_COMM_TO_USE,                           &
                  non_blocking(inonblocking) % request4(kk), &
                  istat4                                     )
          end if
          if( nsend /= 0 ) then
             non_blocking(inonblocking) % count4 = non_blocking(inonblocking) % count4 + 1
             kk = non_blocking(inonblocking) % count4
             call MPI_ISend(                                 &
                  xx_send(1:nsend), nsend4,                  &
                  MPI_DOUBLE_PRECISION, dom_i4, 0_4,         &
                  PAR_COMM_TO_USE,                           &
                  non_blocking(inonblocking) % request4(kk), &
                  istat4                                     )
          end if
       else
          if( nrecv /= 0 .and. nsend == 0 ) then
             call MPI_Recv(                          &
                  xx_recv(1:nrecv), nrecv4,          &
                  MPI_DOUBLE_PRECISION, dom_i4, 0_4, &
                  PAR_COMM_TO_USE, status, istat4   )
          else if( nrecv == 0 .and. nsend /= 0 ) then
             call MPI_Send(                          &
                  xx_send(1:nsend), nsend4,          &
                  MPI_DOUBLE_PRECISION, dom_i4, 0_4, &
                  PAR_COMM_TO_USE, istat4           )
          else if( nrecv /= 0 .and. nsend /= 0 ) then
             call MPI_Sendrecv(                      &
                  xx_send(1:nsend), nsend4,          &
                  MPI_DOUBLE_PRECISION, dom_i4, 0_4, &
                  xx_recv(1:nrecv), nrecv4,          &
                  MPI_DOUBLE_PRECISION, dom_i4, 0_4, &
                  PAR_COMM_TO_USE, status, istat4   )
          end if
       end if
#endif
       if( istat4 /= 0_4 ) call runend('PAR_SEND_RECEIVE_RP: MPI ERROR')

    end if

  end subroutine PAR_SEND_RECEIVE_RP

  !-----------------------------------------------------------------------
  !
  !> @brief   Bridge to MPI_ALLGATHERV
  !> @details Bridge to MPI_ALLGATHERV. If the displacement is not
  !>          prescribed the recvbuf are put one after the other
  !>          automatically
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

    subroutine PAR_GATHERV_RP_1(sendbuf,recvbuf,recvcount4,wherein,displs4)
    implicit none
    real(rp),     pointer, intent(in)           :: sendbuf(:)           !< Send buffer
    real(rp),     pointer, intent(inout)        :: recvbuf(:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in), optional :: wherein              !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rcl,dl,jpart
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

    if( ISEQUEN ) then
       sendcount4 = size(sendbuf)
       recvbuf    = sendbuf
    else if( IPARALL ) then
#ifndef MPI_OFF
       root_rank4      = 0_4
       recvcount4_null = 0_4
       displs4_null    = 0_4
       sendcount4      = 0_4

       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          PAR_COMM_TO_USE = PAR_COMM_MY_CODE
       end if
       if( associated(sendbuf) ) then
          sendcount4 = size(sendbuf)
       end if
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                    MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                    MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
       else
          call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)
          if( my_rank == root_rank4 ) then
             allocate( my_displs4(0:comm_size-1) )
             jpart = lbound(recvcount4,1)
             my_displs4(0) = 0
             do ipart = 1,comm_size-1
                my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
                jpart = jpart + 1
             end do
          else
             allocate( my_displs4(1) )
             my_displs4 = 0
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                   MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                   MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= 0 ) call runend('PAR_GATHERV_RP_1: MPI ERROR')
#endif
    end if

  end subroutine PAR_GATHERV_RP_1

  subroutine PAR_GATHERV_RP_0(sendbuf,recvbuf,sendcount4,recvcount4,wherein,displs4)
    implicit none
    real(rp),              intent(in)           :: sendbuf(*)           !< Send buffer
    real(rp),     pointer, intent(inout)        :: recvbuf(:)           !< Recv buffer
    integer(4),            intent(in)           :: sendcount4           !< Send counts
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in)           :: wherein              !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rcl,dl,jpart
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

    if( ISEQUEN ) then
       recvbuf    = sendbuf(1:sendcount4)
    else if( IPARALL ) then
#ifndef MPI_OFF
       root_rank4      = 0_4
       recvcount4_null = 0_4
       displs4_null    = 0_4

       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                    MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                    MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
       else
          call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)
          if( my_rank == root_rank4 ) then
             allocate( my_displs4(0:comm_size-1) )
             jpart = lbound(recvcount4,1)
             my_displs4(0) = 0
             do ipart = 1,comm_size-1
                my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
                jpart = jpart + 1
             end do
          else
             allocate( my_displs4(1) )
             my_displs4 = 0
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                   MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                   MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= 0 ) call runend('PAR_GATHERV_RP_1: MPI ERROR')
#endif
    end if

  end subroutine PAR_GATHERV_RP_0

  subroutine PAR_GATHERV_RP_21(sendbuf,recvbuf,recvcount4,wherein,displs4)
    implicit none
    real(rp),     pointer, intent(in)           :: sendbuf(:,:)         !< Send buffer
    real(rp),     pointer, intent(inout)        :: recvbuf(:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in)           :: wherein              !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rcl,dl,jpart,ii,jj,kk
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

    if( ISEQUEN ) then
       sendcount4 = size(sendbuf)
       kk = lbound(recvbuf,1)
       do jj = lbound(sendbuf,2),ubound(sendbuf,2)
          do ii = lbound(sendbuf,1),ubound(sendbuf,1)
             recvbuf(kk) = sendbuf(ii,jj)
             kk = kk + 1
          end do
       end do
    else if( IPARALL ) then
#ifndef MPI_OFF
       root_rank4      = 0_4
       recvcount4_null = 0_4
       displs4_null    = 0_4
       sendcount4      = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       if( associated(sendbuf) ) then
          sendcount4 = size(sendbuf)
       end if
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                    MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                    MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
       else
          call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)
          if( my_rank == root_rank4 ) then
             allocate( my_displs4(0:comm_size-1) )
             jpart = lbound(recvcount4,1)
             my_displs4(0) = 0
             do ipart = 1,comm_size-1
                my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
                jpart = jpart + 1
             end do
          else
             allocate( my_displs4(1) )
             my_displs4(1) = 0
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                   MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                   MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= 0 ) call runend('PAR_GATHERV_RP_2: MPI ERROR')
#endif
    end if

  end subroutine PAR_GATHERV_RP_21

  subroutine PAR_GATHERV_RP_22(sendbuf,recvbuf,recvcount4,wherein,displs4)
    implicit none
    real(rp),     pointer, intent(in)           :: sendbuf(:,:)         !< Send buffer
    real(rp),     pointer, intent(inout)          :: recvbuf(:,:)         !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rcl,dl,jpart,ii,jj,kk
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

    if( ISEQUEN ) then
       sendcount4 = size(sendbuf)
       kk = lbound(recvbuf,1)
       do jj = lbound(sendbuf,2),ubound(sendbuf,2)
          do ii = lbound(sendbuf,1),ubound(sendbuf,1)
             recvbuf(ii,jj) = sendbuf(ii,jj)
             kk = kk + 1
          end do
       end do
    else if( IPARALL ) then
#ifndef MPI_OFF
       root_rank4      = 0_4
       recvcount4_null = 0_4
       displs4_null    = 0_4
       sendcount4      = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       if( associated(sendbuf) ) then
          sendcount4 = size(sendbuf)
       end if
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                    MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                    MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
       else
          call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)
          if( my_rank == root_rank4 ) then
             allocate( my_displs4(0:comm_size-1) )
             jpart = lbound(recvcount4,1)
             my_displs4(0) = 0
             do ipart = 1,comm_size-1
                my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
                jpart = jpart + 1
             end do
          else
             allocate( my_displs4(1) )
             my_displs4(1) = 0
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                   MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                   MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= 0 ) call runend('PAR_GATHERV_RP_2: MPI ERROR')
#endif
    end if

  end subroutine PAR_GATHERV_RP_22

  subroutine PAR_GATHERV_IP_1(sendbuf,recvbuf,recvcount4,wherein,displs4,PAR_COMM_IN4)
    implicit none
    integer(ip),  pointer, intent(in)           :: sendbuf(:)           !< Send buffer
    integer(ip),  pointer, intent(inout)        :: recvbuf(:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in), optional :: wherein              !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4),   pointer, intent(in), optional :: PAR_COMM_IN4         !< Communicator
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rl,sl,rcl,dl,jpart
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    if( IPARALL ) then
       root_rank4      = 0_4
       recvcount4_null = 0_4
       displs4_null    = 0_4
       sendcount4      = 0_4
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
       end if

       if( associated(sendbuf) ) then
          sendcount4 = size(sendbuf)
          sl = lbound(sendbuf,1)
       end if
       if( associated(recvbuf) )    rl  = lbound(recvbuf,1)
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                    PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                    PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
       else
          call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)
          if( my_rank == root_rank4 ) then
             allocate( my_displs4(0:comm_size-1) )
             jpart = lbound(recvcount4,1)
             my_displs4(0) = 0
             do ipart = 1,comm_size-1
                my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
                jpart = jpart + 1
             end do
          else
             allocate( my_displs4(1) )
             my_displs4(1) = 0
           end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                   PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                   PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= 0 ) call runend('PAR_GATHERV_IP_1: MPI ERROR')
    end if
#endif

  end subroutine PAR_GATHERV_IP_1

  subroutine PAR_GATHERV_IP_21(sendbuf,recvbuf,recvcount4,wherein,displs4)
    implicit none
    integer(ip),  pointer, intent(in)           :: sendbuf(:,:)         !< Send buffer
    integer(ip),  pointer, intent(inout)          :: recvbuf(:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in), optional :: wherein              !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rl,rcl,dl
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)
    integer(ip)                                 :: jpart

#ifndef MPI_OFF
    if( IPARALL ) then
       root_rank4      = 0_4
       recvcount4_null = 0_4
       displs4_null    = 0_4
       sendcount4      = 0_4
       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
       end if
       if( associated(sendbuf) ) then
          sendcount4 = size(sendbuf,1)*size(sendbuf,2)
       end if
       if( associated(recvbuf) )    rl  = lbound(recvbuf,1)
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if
       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                    PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                    PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
       else
          call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)
          if( my_rank == root_rank4 ) then
             allocate( my_displs4(0:comm_size-1) )
             jpart = lbound(recvcount4,1)
             my_displs4(0) = 0
             do ipart = 1,comm_size-1
                my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
                jpart = jpart + 1
             end do
          else
             allocate( my_displs4(1) )
             my_displs4(1) = 0
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                   PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                   PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= 0 ) call runend('PAR_GATHERV_IP_2: MPI ERROR')
    end if
#endif

  end subroutine PAR_GATHERV_IP_21

  subroutine PAR_GATHERV_IP_22(sendbuf,recvbuf,recvcount4,wherein,displs4)
    implicit none
    integer(ip),  pointer, intent(in)           :: sendbuf(:,:)         !< Send buffer
    integer(ip),  pointer, intent(inout)        :: recvbuf(:,:)         !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in), optional :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rcl,dl,jpart,ii,jj,kk
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

    if( ISEQUEN ) then
       sendcount4 = size(sendbuf)
       kk = lbound(recvbuf,1)
       do jj = lbound(sendbuf,2),ubound(sendbuf,2)
          do ii = lbound(sendbuf,1),ubound(sendbuf,1)
             recvbuf(ii,jj) = sendbuf(ii,jj)
             kk = kk + 1
          end do
       end do
    else if( IPARALL ) then
#ifndef MPI_OFF
       root_rank4      = 0_4
       recvcount4_null = 0_4
       displs4_null    = 0_4
       sendcount4      = 0_4
       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
       end if
       if( associated(sendbuf) ) then
          sendcount4 = size(sendbuf)
       end if
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                    PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                    PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
       else
          call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)
          if( my_rank == root_rank4 ) then
             allocate( my_displs4(0:comm_size-1) )
             jpart = lbound(recvcount4,1)
             my_displs4(0) = 0
             do ipart = 1,comm_size-1
                my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
                jpart = jpart + 1
             end do
          else
             allocate( my_displs4(1) )
             my_displs4(1) = 0
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                   PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                   PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= 0 ) call runend('PAR_GATHERV_IP_22: MPI ERROR')
#endif
    end if

  end subroutine PAR_GATHERV_IP_22

  subroutine PAR_GATHERV_RP_21_SEND(sendbuf,recvbuf,sendcount4,recvcount4,wherein,displs4)
    implicit none
    real(rp),     pointer, intent(in)           :: sendbuf(:,:)         !< Send buffer
    real(rp),     pointer, intent(inout)          :: recvbuf(:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4), intent(in)                      :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rcl,dl,jpart,ii,jj,kk
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

    if( ISEQUEN ) then
       kk = lbound(recvbuf,1)
       do jj = lbound(sendbuf,2),ubound(sendbuf,2)
          do ii = lbound(sendbuf,1),ubound(sendbuf,1)
             recvbuf(kk) = sendbuf(ii,jj)
             kk = kk + 1
          end do
       end do
    else if( IPARALL ) then
#ifndef MPI_OFF
       root_rank4      = 0_4
       recvcount4_null = 0_4
       displs4_null    = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                    MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                    MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
       else
          call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)
          if( my_rank == root_rank4 ) then
             allocate( my_displs4(0:comm_size-1) )
             jpart = lbound(recvcount4,1)
             my_displs4(0) = 0
             do ipart = 1,comm_size-1
                my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
                jpart = jpart + 1
             end do
          else
             allocate( my_displs4(1) )
             my_displs4(1) = 0
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                   MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                   MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= 0 ) call runend('PAR_GATHERV_RP_2: MPI ERROR')
#endif
    end if

  end subroutine PAR_GATHERV_RP_21_SEND

  subroutine PAR_GATHERV_RP_22_SEND(sendbuf,recvbuf,sendcount4,recvcount4,wherein,displs4)
    implicit none
    real(rp),     pointer, intent(in)           :: sendbuf(:,:)         !< Send buffer
    real(rp),     pointer, intent(inout)          :: recvbuf(:,:)         !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4), intent(in)                      :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rcl,dl,jpart,ii,jj,kk
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

    if( ISEQUEN ) then
       kk = lbound(recvbuf,1)
       do jj = lbound(sendbuf,2),ubound(sendbuf,2)
          do ii = lbound(sendbuf,1),ubound(sendbuf,1)
             recvbuf(ii,jj) = sendbuf(ii,jj)
             kk = kk + 1
          end do
       end do
    else if( IPARALL ) then
#ifndef MPI_OFF
       root_rank4      = 0_4
       recvcount4_null = 0_4
       displs4_null    = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                    MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                    MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
       else
          call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)
          if( my_rank == root_rank4 ) then
             allocate( my_displs4(0:comm_size-1) )
             jpart = lbound(recvcount4,1)
             my_displs4(0) = 0
             do ipart = 1,comm_size-1
                my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
                jpart = jpart + 1
             end do
          else
             allocate( my_displs4(1) )
             my_displs4(1) = 0
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                   MPI_DOUBLE_PRECISION,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                   MPI_DOUBLE_PRECISION,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= 0 ) call runend('PAR_GATHERV_RP_2: MPI ERROR')
#endif
    end if

  end subroutine PAR_GATHERV_RP_22_SEND

  subroutine PAR_GATHERV_IP_1_SEND(sendbuf,recvbuf,sendcount4,recvcount4,wherein,displs4)
    implicit none
    integer(ip),  pointer, intent(in)           :: sendbuf(:)           !< Send buffer
    integer(ip),  pointer, intent(inout)          :: recvbuf(:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in), optional :: wherein              !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4), intent(in)                      :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rl,sl,rcl,dl,jpart
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    if( IPARALL ) then
       root_rank4      = 0_4
       recvcount4_null = 0_4
       displs4_null    = 0_4
       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
       end if
       if( associated(sendbuf) ) then
          sl = lbound(sendbuf,1)
       end if
       if( associated(recvbuf) )    rl  = lbound(recvbuf,1)
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                    PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                    PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
       else
          call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)
          if( my_rank == root_rank4 ) then
             allocate( my_displs4(0:comm_size-1) )
             jpart = lbound(recvcount4,1)
             my_displs4(0) = 0
             do ipart = 1,comm_size-1
                my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
                jpart = jpart + 1
             end do
          else
             allocate( my_displs4(1) )
             my_displs4(1) = 0
           end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                   PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                   PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= 0 ) call runend('PAR_GATHERV_IP_1: MPI ERROR')
    end if
#endif

  end subroutine PAR_GATHERV_IP_1_SEND

  subroutine PAR_GATHERV_IP_21_SEND(sendbuf,recvbuf,sendcount4,recvcount4,wherein,displs4)
    implicit none
    integer(ip),  pointer, intent(in)           :: sendbuf(:,:)         !< Send buffer
    integer(ip),  pointer, intent(inout)          :: recvbuf(:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in), optional :: wherein              !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4), intent(in)                      :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rl,rcl,dl
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)
    integer(ip)                                 :: jpart

#ifndef MPI_OFF
    if( IPARALL ) then
       root_rank4      = 0_4
       recvcount4_null = 0_4
       displs4_null    = 0_4
       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
       end if
       if( associated(recvbuf) )    rl  = lbound(recvbuf,1)
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if
       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                    PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                    PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
       else
          call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)
          if( my_rank == root_rank4 ) then
             allocate( my_displs4(0:comm_size-1) )
             jpart = lbound(recvcount4,1)
             my_displs4(0) = 0
             do ipart = 1,comm_size-1
                my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
                jpart = jpart + 1
             end do
          else
             allocate( my_displs4(1) )
             my_displs4(1) = 0
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                   PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                   PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= 0 ) call runend('PAR_GATHERV_IP_2: MPI ERROR')
    end if
#endif

  end subroutine PAR_GATHERV_IP_21_SEND


  subroutine PAR_GATHERV_IP_22_SEND(sendbuf,recvbuf,sendcount4,recvcount4,wherein,displs4)
    implicit none
    integer(ip),     pointer, intent(in)        :: sendbuf(:,:)         !< Send buffer
    integer(ip),     pointer, intent(inout)       :: recvbuf(:,:)         !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4), intent(in)                      :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(ip)                                 :: recvbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4)                                  :: root_rank4
    integer(4)                                  :: my_rank
    integer(ip)                                 :: rcl,dl,jpart,ii,jj,kk
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

    if( ISEQUEN ) then
       kk = lbound(recvbuf,1)
       do jj = lbound(sendbuf,2),ubound(sendbuf,2)
          do ii = lbound(sendbuf,1),ubound(sendbuf,1)
             recvbuf(ii,jj) = sendbuf(ii,jj)
             kk = kk + 1
          end do
       end do
    else if( IPARALL ) then
#ifndef MPI_OFF
       root_rank4      = 0_4
       recvcount4_null = 0_4
       displs4_null    = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,                PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                    PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,displs4_tmp,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                    PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
       else
          call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank,comm_size)
          if( my_rank == root_rank4 ) then
             allocate( my_displs4(0:comm_size-1) )
             jpart = lbound(recvcount4,1)
             my_displs4(0) = 0
             do ipart = 1,comm_size-1
                my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
                jpart = jpart + 1
             end do
          else
             allocate( my_displs4(1) )
             my_displs4(1) = 0
          end if
          if( sendcount4 == 0 ) then
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf_tmp,sendcount4,               PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          else
             if( associated(recvbuf) ) then
                call MPI_GATHERV(sendbuf,sendcount4,                   PAR_INTEGER,&
                     &           recvbuf,recvcount4_tmp,my_displs4,    PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             else
                call MPI_GATHERV(sendbuf,sendcount4,                   PAR_INTEGER,&
                     &           recvbuf_tmp,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                     &           root_rank4,PAR_COMM_TO_USE,istat4)
             end if
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= 0 ) call runend('PAR_GATHERV_IP_22: MPI ERROR')
#endif
    end if

  end subroutine PAR_GATHERV_IP_22_SEND

  !-----------------------------------------------------------------------
  !
  !> @brief   Bridge to MPI_ALLGATHERV
  !> @details Bridge to MPI_ALLGATHERV. If the displacement is not
  !>          prescribed the recvbuf are put one after the other
  !>          automatically
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine PAR_ALLGATHERV_IP4(sendbuf,recvbuf,recvcount4,wherein,displs4,PAR_COMM_IN4)
    implicit none
    integer(ip),  pointer, intent(in)           :: sendbuf(:)           !< Send buffer
    integer(ip),  pointer, intent(inout)        :: recvbuf(:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in), optional :: wherein              !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4),            intent(in), optional :: PAR_COMM_IN4         !< Communicator
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(4)                                  :: ipart,lbouns,lbounr,lbounc
    integer(ip)                                 :: rl,rcl,dl,jpart
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

    if( IPARALL ) then
#ifndef MPI_OFF
       recvcount4_null = 0_4
       displs4_null    = 0_4
       sendcount4      = 0_4
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
       end if
       if( associated(sendbuf) ) then
          sendcount4 = size(sendbuf)
       end if
       rl = lbound(recvcount4,1)
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             CALL MPI_ALLGATHERV(sendbuf_tmp,sendcount4,            PAR_INTEGER,&
                  &              recvbuf,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                  &              PAR_COMM_TO_USE,istat4)
          else
             CALL MPI_ALLGATHERV(sendbuf,sendcount4,                PAR_INTEGER,&
                  &              recvbuf,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                  &              PAR_COMM_TO_USE,istat4)
          end if
       else
          call MPI_Comm_size(PAR_COMM_TO_USE,comm_size,istat4)
          allocate( my_displs4(0:comm_size-1) )
          jpart = lbound(recvcount4,1)
          my_displs4(0) = 0
          do ipart = 1,comm_size-1
             my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
             jpart = jpart + 1
          end do
          if( sendcount4 == 0 ) then
             CALL MPI_ALLGATHERV(sendbuf_tmp,sendcount4,           PAR_INTEGER,&
                  &              recvbuf,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                  &              PAR_COMM_TO_USE,istat4)
          else
             CALL MPI_ALLGATHERV(sendbuf,sendcount4,               PAR_INTEGER,&
                  &              recvbuf,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                  &              PAR_COMM_TO_USE,istat4)
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= 0_4 ) call runend('PAR_ALLGATHERV_IP4: MPI ERROR')
#endif
    else
       lbounr = lbound(recvbuf,1)
       lbouns = lbound(sendbuf,1)
       lbounc = lbound(recvcount4,1)
       recvbuf(lbounr:lbounr+recvcount4(lbounc)-1) = sendbuf(lbouns:lbouns+recvcount4(lbounc)-1)
    end if

  end subroutine PAR_ALLGATHERV_IP4

  subroutine PAR_ALLGATHERV_IP8(sendbuf,recvbuf,recvcount8,wherein,displs8,PAR_COMM_IN4)
    implicit none
    integer(ip),  pointer, intent(in)           :: sendbuf(:)           !< Send buffer
    integer(ip),  pointer, intent(inout)        :: recvbuf(:)           !< Recv buffer
    integer(8),   pointer, intent(in)           :: recvcount8(:)        !< Recv counts
    character(*),          intent(in), optional :: wherein              !< Wherein
    integer(8),   pointer, intent(in), optional :: displs8(:)           !< Displacement
    integer(4),            intent(in), optional :: PAR_COMM_IN4         !< Communicator
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(4),   pointer                       :: displs4(:)
    integer(4),   pointer                       :: recvcount4(:)
    integer(ip)                                 :: jpart
    integer(4)                                  :: lbounr,lbouns,lbounc
    integer(4)                                  :: dumm4

    if( IPARALL ) then
#ifndef MPI_OFF
       sendcount4 = 0_4
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,4)
       end if
       if( associated(sendbuf) ) then
          sendcount4 = size(sendbuf)
       end if
       allocate( recvcount4(max(1_ip,size(recvcount8,KIND=ip))) )
       recvcount4 = int(recvcount8,4)
       if( present(displs8) ) then
          if( .not. associated(displs8) ) then
             allocate( displs4(1) )
             displs4 = 0
          else
             allocate( displs4(size(displs8)) )
             displs4 = int(displs8,4)
          end if
          if( sendcount4 == 0 ) then
             CALL MPI_ALLGATHERV(sendbuf_tmp,sendcount4,    PAR_INTEGER,&
                  &              recvbuf,recvcount4,displs4,PAR_INTEGER,&
                  &              PAR_COMM_TO_USE,istat4)
          else
             CALL MPI_ALLGATHERV(sendbuf,sendcount4,        PAR_INTEGER,&
                  &              recvbuf,recvcount4,displs4,PAR_INTEGER,&
                  &              PAR_COMM_TO_USE,istat4)
          end if
          deallocate(displs4)
       else
          call MPI_Comm_size(PAR_COMM_TO_USE,comm_size,istat4)
          allocate( my_displs4(0:comm_size-1) )
          jpart = lbound(recvcount8,1)
          my_displs4(0) = 0
          do ipart = 1,comm_size-1
             my_displs4(ipart) = my_displs4(ipart-1) + int(recvcount8(jpart),4)
             jpart = jpart + 1
          end do
          if( sendcount4 == 0 ) then
             CALL MPI_ALLGATHERV(sendbuf_tmp,sendcount4,       PAR_INTEGER,&
                  &              recvbuf,recvcount4,my_displs4,PAR_INTEGER,&
                  &              PAR_COMM_TO_USE,istat4)
          else
             CALL MPI_ALLGATHERV(sendbuf,sendcount4,           PAR_INTEGER,&
                  &              recvbuf,recvcount4,my_displs4,PAR_INTEGER,&
                  &              PAR_COMM_TO_USE,istat4)
          end if
          deallocate( my_displs4 )
       end if
       deallocate( recvcount4 )
       if( istat4 /= 0_4 ) call runend('PAR_ALLGATHERV_IP8: MPI ERROR')
#endif
    else
       lbounr = lbound(recvbuf,1)
       lbouns = lbound(sendbuf,1)
       lbounc = lbound(recvcount8,1)
       dumm4  = int(recvcount8(lbounc),4)
       recvbuf(lbounr:lbounr+dumm4-1) = sendbuf(lbouns:lbouns+dumm4-1)
    end if

  end subroutine PAR_ALLGATHERV_IP8

  subroutine PAR_ALLGATHERV_IP4_2(sendbuf,recvbuf,recvcount4,wherein,displs4)
    implicit none
    integer(ip),  pointer, intent(in)           :: sendbuf(:,:)         !< Send buffer
    integer(ip),  pointer, intent(inout)          :: recvbuf(:,:)         !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(ip)                                 :: rl,rcl,dl,jpart
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    if( IPARALL ) then
       recvcount4_null = 0_4
       displs4_null    = 0_4
       sendcount4      = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       if( associated(sendbuf) ) then
          sendcount4 = size(sendbuf)
       end if
       rl = lbound(recvcount4,1)
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             CALL MPI_ALLGATHERV(sendbuf_tmp,sendcount4,            PAR_INTEGER,&
                  &              recvbuf,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                  &              PAR_COMM_TO_USE,istat4)
          else
             CALL MPI_ALLGATHERV(sendbuf,sendcount4,                PAR_INTEGER,&
                  &              recvbuf,recvcount4_tmp,displs4_tmp,PAR_INTEGER,&
                  &              PAR_COMM_TO_USE,istat4)
          end if
       else
          call MPI_Comm_size(PAR_COMM_TO_USE,comm_size,istat4)
          allocate( my_displs4(0:comm_size-1) )
          jpart = lbound(recvcount4,1)
          my_displs4(0) = 0
          do ipart = 1,comm_size-1
             my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
             jpart = jpart + 1
          end do
          if( sendcount4 == 0 ) then
             CALL MPI_ALLGATHERV(sendbuf_tmp,sendcount4,           PAR_INTEGER,&
                  &              recvbuf,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                  &              PAR_COMM_TO_USE,istat4)
          else
             CALL MPI_ALLGATHERV(sendbuf,sendcount4,               PAR_INTEGER,&
                  &              recvbuf,recvcount4_tmp,my_displs4,PAR_INTEGER,&
                  &              PAR_COMM_TO_USE,istat4)
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= 0_4 ) call runend('PAR_ALLGATHERV_IP4_2: MPI ERROR')
    end if
#endif

  end subroutine PAR_ALLGATHERV_IP4_2


 subroutine PAR_ALLGATHERV_RP_1(sendbuf,recvbuf,recvcount4,wherein,displs4)
    implicit none
    real(rp),     pointer, intent(in)           :: sendbuf(:)           !< Send buffer
    real(rp),     pointer, intent(inout)          :: recvbuf(:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(ip)                                 :: rl,rcl,dl,jpart
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    if( IPARALL ) then
       recvcount4_null = 0_4
       displs4_null    = 0_4
       sendcount4      = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       if( associated(sendbuf) ) then
          sendcount4 = size(sendbuf)
       end if
       rl = lbound(recvcount4,1)
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             CALL MPI_ALLGATHERV(sendbuf_tmp,sendcount4,            MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          else
             CALL MPI_ALLGATHERV(sendbuf,sendcount4,                MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          end if
       else
          call MPI_Comm_size(PAR_COMM_TO_USE,comm_size,istat4)
          allocate( my_displs4(0:comm_size-1) )
          jpart = lbound(recvcount4,1)
          my_displs4(0) = 0
          do ipart = 1,comm_size-1
             my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
             jpart = jpart + 1
          end do
          if( sendcount4 == 0 ) then
             CALL MPI_ALLGATHERV(sendbuf_tmp,sendcount4,           MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          else
             CALL MPI_ALLGATHERV(sendbuf,sendcount4,               MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= 0_4 ) call runend('PAR_ALLGATHERV_RP1: MPI ERROR')
    end if
#endif

  end subroutine PAR_ALLGATHERV_RP_1

 subroutine PAR_ALLGATHERV_RP_18(sendbuf,recvbuf,recvcount8,wherein,displs8)
    implicit none
    real(rp),     pointer, intent(in)           :: sendbuf(:)           !< Send buffer
    real(rp),     pointer, intent(inout)          :: recvbuf(:)           !< Recv buffer
    integer(8),   pointer, intent(in)           :: recvcount8(:)        !< Recv counts
    integer(4),   pointer                       :: recvcount4(:)        !< Recv counts integer 4
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(8),   pointer, intent(in), optional :: displs8(:)           !< Displacement
    integer(4),   pointer                       :: displs4(:)           !< Displacement integer 4
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(ip)                                 :: rl,rcl,dl,jpart
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    if( IPARALL ) then
       recvcount4_null = 0_4
       displs4_null    = 0_4
       sendcount4      = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       if( associated(sendbuf) ) then
          sendcount4 = size(sendbuf)
       end if
       allocate( recvcount4(max(1_ip,size(recvcount8,KIND=ip))) )
       recvcount4 = int(recvcount8,4)
       rl = lbound(recvcount4,1)
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if
      if( present(displs8) ) then
          if( .not. associated(displs8) ) then
             allocate( displs4(1) )
             displs4 = 0
          else
             allocate( displs4(size(displs8)) )
             displs4 = int(displs8,4)
          end if
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             CALL MPI_ALLGATHERV(sendbuf_tmp,sendcount4,            MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          else
             CALL MPI_ALLGATHERV(sendbuf,sendcount4,                MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          end if
       else
          call MPI_Comm_size(PAR_COMM_TO_USE,comm_size,istat4)
          allocate( my_displs4(0:comm_size-1) )
          jpart = lbound(recvcount4,1)
          my_displs4(0) = 0
          do ipart = 1,comm_size-1
             my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
             jpart = jpart + 1
          end do
          if( sendcount4 == 0 ) then
             CALL MPI_ALLGATHERV(sendbuf_tmp,sendcount4,           MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          else
             CALL MPI_ALLGATHERV(sendbuf,sendcount4,               MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= 0_4 ) call runend('PAR_ALLGATHERV_RP1: MPI ERROR')
    end if
#endif
  end subroutine PAR_ALLGATHERV_RP_18

  subroutine PAR_ALLGATHERV_RP_2(sendbuf,recvbuf,recvcount4,wherein,displs4)
    implicit none
    real(rp),     pointer, intent(in)           :: sendbuf(:,:)         !< Send buffer
    real(rp),     pointer, intent(inout)          :: recvbuf(:,:)         !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(ip)                                 :: rl,rcl,dl,jpart
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    if( IPARALL ) then
       recvcount4_null = 0_4
       displs4_null    = 0_4
       sendcount4      = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       if( associated(sendbuf) ) then
          sendcount4 = size(sendbuf)
       end if
       rl = lbound(recvcount4,1)
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             CALL MPI_ALLGATHERV(sendbuf_tmp,sendcount4,            MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          else
             CALL MPI_ALLGATHERV(sendbuf,sendcount4,                MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          end if
       else
          call MPI_Comm_size(PAR_COMM_TO_USE,comm_size,istat4)
          allocate( my_displs4(0:comm_size-1) )
          jpart = lbound(recvcount4,1)
          my_displs4(0) = 0
          do ipart = 1,comm_size-1
             my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
             jpart = jpart + 1
          end do
          if( sendcount4 == 0 ) then
             CALL MPI_ALLGATHERV(sendbuf_tmp,sendcount4,           MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          else
             CALL MPI_ALLGATHERV(sendbuf,sendcount4,               MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= 0_4 ) call runend('PAR_ALLGATHERV_IP4: MPI ERROR')
    end if
#endif

  end subroutine PAR_ALLGATHERV_RP_2

  subroutine PAR_ALLGATHERV_RP_3(sendbuf,recvbuf,recvcount4,wherein,displs4)
    implicit none
    real(rp),     pointer, intent(in)           :: sendbuf(:,:,:)       !< Send buffer
    real(rp),     pointer, intent(inout)          :: recvbuf(:,:,:)       !< Recv buffer
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in)           :: wherein                !< Wherein
    integer(4),   pointer, intent(in), optional :: displs4(:)           !< Displacement
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4),   pointer                       :: my_displs4(:)
    integer(ip)                                 :: sendbuf_tmp(2)
    integer(4)                                  :: ipart
    integer(ip)                                 :: rl,rcl,dl,jpart
    integer(4),   pointer                       :: recvcount4_tmp(:)
    integer(4),   target                        :: recvcount4_null(2)
    integer(4),   pointer                       :: displs4_tmp(:)
    integer(4),   target                        :: displs4_null(2)

#ifndef MPI_OFF
    if( IPARALL ) then
       recvcount4_null = 0_4
       displs4_null    = 0_4
       sendcount4      = 0_4
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       if( associated(sendbuf) ) then
          sendcount4 = size(sendbuf)
       end if
       rl = lbound(recvcount4,1)
       if( associated(recvcount4) ) then
          rcl = lbound(recvcount4,1)
          recvcount4_tmp => recvcount4(rcl:)
       else
          recvcount4_tmp => recvcount4_null
       end if

       if( present(displs4) ) then
          if( associated(displs4) ) then
             dl = lbound(displs4,1)
             displs4_tmp => displs4(dl:)
          else
             displs4_tmp => displs4_null
          end if
          if( sendcount4 == 0 ) then
             CALL MPI_ALLGATHERV(sendbuf_tmp,sendcount4,            MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          else
             CALL MPI_ALLGATHERV(sendbuf,sendcount4,                MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,displs4_tmp,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          end if
       else
          call MPI_Comm_size(PAR_COMM_TO_USE,comm_size,istat4)
          allocate( my_displs4(0:comm_size-1) )
          jpart = lbound(recvcount4,1)
          my_displs4(0) = 0
          do ipart = 1,comm_size-1
             my_displs4(ipart) = my_displs4(ipart-1) + recvcount4(jpart)
             jpart = jpart + 1
          end do
          if( sendcount4 == 0 ) then
             CALL MPI_ALLGATHERV(sendbuf_tmp,sendcount4,           MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          else
             CALL MPI_ALLGATHERV(sendbuf,sendcount4,               MPI_DOUBLE_PRECISION,&
                  &              recvbuf,recvcount4_tmp,my_displs4,MPI_DOUBLE_PRECISION,&
                  &              PAR_COMM_TO_USE,istat4)
          end if
          deallocate( my_displs4 )
       end if
       if( istat4 /= 0_4 ) call runend('PAR_ALLGATHERV_IP4: MPI ERROR')
    end if
#endif

  end subroutine PAR_ALLGATHERV_RP_3

  !-----------------------------------------------------------------------
  !
  !> @brief   Bridge to MPI_ALLGATHER
  !> @details Bridge to MPI_ALLGATHER
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine PAR_ALLGATHER_s4(sendbuf,recvbuf,recvcount4_opt,wherein,PAR_COMM_IN4)
    implicit none
    integer(4),            intent(in)           :: sendbuf              !< Send buffer
    integer(4),  pointer,  intent(inout)        :: recvbuf(:)           !< Recv buffer
    integer(4),            intent(in), optional :: recvcount4_opt       !< Recv counts
    character(*),          intent(in), optional :: wherein              !< Wherein
    integer(4),            intent(in), optional :: PAR_COMM_IN4         !< Wherein
    integer(4)                                  :: istat4,lbounr
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4)                                  :: recvcount4

    if( IPARALL ) then
#ifndef MPI_OFF
       if( present(recvcount4_opt) ) then
          recvcount4 = recvcount4_opt
       else
          recvcount4 = 1_4
       end if
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
          end if
       end if
       sendcount4 = 1_4
       call MPI_ALLGATHER(sendbuf,sendcount4,MPI_INTEGER,&
            &             recvbuf,recvcount4,MPI_INTEGER,&
            &             PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0 ) call runend('PAR_ALLGATHER_s: MPI ERROR')
#endif
    else
       if( associated(recvbuf) ) then
          lbounr = lbound(recvbuf,1)
          recvbuf(lbounr) = sendbuf
       end if
    end if

  end subroutine PAR_ALLGATHER_s4

  subroutine PAR_ALLGATHER_s8(sendbuf,recvbuf,recvcount4_opt,wherein,PAR_COMM_IN4)
    implicit none
    integer(8),            intent(in)           :: sendbuf              !< Send buffer
    integer(8),  pointer,  intent(inout)        :: recvbuf(:)           !< Recv buffer
    integer(4),            intent(in), optional :: recvcount4_opt       !< Recv counts
    character(*),          intent(in), optional :: wherein              !< Wherein
    integer(4),            intent(in), optional :: PAR_COMM_IN4         !< Wherein
    integer(4)                                  :: istat4,lbounr
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: sendcount4
    integer(4)                                  :: recvcount4

    if( IPARALL ) then
#ifndef MPI_OFF
       if( present(recvcount4_opt) ) then
          recvcount4 = recvcount4_opt
       else
          recvcount4 = 1_4
       end if
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
          end if
       end if
       sendcount4 = 1_4
       call MPI_ALLGATHER(sendbuf,sendcount4,MPI_INTEGER8,&
            &             recvbuf,recvcount4,MPI_INTEGER8,&
            &             PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0 ) call runend('PAR_ALLGATHER_s: MPI ERROR')
#endif
    else
       if( associated(recvbuf) ) then
          lbounr = lbound(recvbuf,1)
          recvbuf(lbounr) = sendbuf
       end if
    end if

  end subroutine PAR_ALLGATHER_s8

  subroutine PAR_ALLGATHER_IP_14(sendbuf,recvbuf,recvcount4,wherein,PAR_COMM_IN4)
    implicit none
    integer(4),  pointer,  intent(in)             :: sendbuf(:)           !< Send buffer
    integer(4),  pointer,  intent(inout)          :: recvbuf(:)           !< Recv buffer
    integer(4),            intent(in)             :: recvcount4           !< Recv counts
    character(*),          intent(in),   optional :: wherein              !< Wherein
    integer(4),            intent(in),   optional :: PAR_COMM_IN4         !< Communicator
    integer(4)                                    :: istat4
    integer(4)                                    :: lbouns,lbounr
    integer(4)                                    :: PAR_COMM_TO_USE
    integer(4)                                    :: sendcount4
    integer(ip)                                   :: sendnul(2)

    if( IPARALL ) then
#ifndef MPI_OFF
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
          end if
       end if
       if( associated(sendbuf) ) then
          sendcount4 = size(sendbuf)
          call MPI_ALLGATHER(sendbuf,sendcount4,MPI_INTEGER,&
               &             recvbuf,recvcount4,MPI_INTEGER,&
               &             PAR_COMM_TO_USE,istat4)
       else
          sendcount4 = 0_4
          sendnul    = 0_ip
          call MPI_ALLGATHER(sendnul,sendcount4,MPI_INTEGER,&
               &             recvbuf,recvcount4,MPI_INTEGER,&
               &             PAR_COMM_TO_USE,istat4)
       end if
       if( istat4 /= 0 ) call runend('PAR_ALLGATHER_IP: MPI ERROR')
#endif
    else
       lbounr = lbound(recvbuf,1)
       lbouns = lbound(sendbuf,1)
       recvbuf(lbounr:lbounr+recvcount4-1) = sendbuf(lbouns:lbouns+recvcount4-1)
    end if

  end subroutine PAR_ALLGATHER_IP_14

  subroutine PAR_ALLGATHER_IP_18(sendbuf,recvbuf,recvcount4,wherein,PAR_COMM_IN4)
    implicit none
    integer(8),  pointer,  intent(in)             :: sendbuf(:)           !< Send buffer
    integer(8),  pointer,  intent(inout)          :: recvbuf(:)           !< Recv buffer
    integer(4),            intent(in)             :: recvcount4           !< Recv counts
    character(*),          intent(in),   optional :: wherein              !< Wherein
    integer(4),            intent(in),   optional :: PAR_COMM_IN4         !< Communicator
    integer(4)                                    :: istat4
    integer(4)                                    :: lbouns,lbounr
    integer(4)                                    :: PAR_COMM_TO_USE
    integer(4)                                    :: sendcount4
    integer(ip)                                   :: sendnul(2)

    if( IPARALL ) then
#ifndef MPI_OFF
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else
          if( present(wherein) ) then
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
          else
             PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
          end if
       end if
       if( associated(sendbuf) ) then
          sendcount4 = size(sendbuf)
          call MPI_ALLGATHER(sendbuf,sendcount4,MPI_INTEGER8,&
               &             recvbuf,recvcount4,MPI_INTEGER8,&
               &             PAR_COMM_TO_USE,istat4)
       else
          sendcount4 = 0_4
          sendnul    = 0_ip
          call MPI_ALLGATHER(sendnul,sendcount4,MPI_INTEGER8,&
               &             recvbuf,recvcount4,MPI_INTEGER8,&
               &             PAR_COMM_TO_USE,istat4)
       end if
       if( istat4 /= 0 ) call runend('PAR_ALLGATHER_IP: MPI ERROR')
#endif
    else
       lbounr = lbound(recvbuf,1)
       lbouns = lbound(sendbuf,1)
       recvbuf(lbounr:lbounr+recvcount4-1) = sendbuf(lbouns:lbouns+recvcount4-1)
    end if

  end subroutine PAR_ALLGATHER_IP_18

  subroutine PAR_ALLGATHER_RP_0(sendbuf,recvbuf,recvcount4,wherein)
    implicit none
    real(rp),              intent(in)  :: sendbuf(*)           !< Send buffer
    real(rp),              intent(out) :: recvbuf(*)           !< Recv buffer
    integer(4),            intent(in)  :: recvcount4           !< Recv counts
    character(*),          intent(in)  :: wherein                !< Wherein
    integer(4)                         :: istat4
    integer(4)                         :: PAR_COMM_TO_USE

#ifndef MPI_OFF
    if( IPARALL ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call MPI_ALLGATHER(sendbuf,recvcount4,MPI_DOUBLE_PRECISION,&
            &             recvbuf,recvcount4,MPI_DOUBLE_PRECISION,&
            &             PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0 ) call runend('PAR_ALLGATHER_IP: MPI ERROR')
    end if
#endif

  end subroutine PAR_ALLGATHER_RP_0

  subroutine PAR_ALLGATHER_RP_2(sendbuf,recvbuf,recvcount4,wherein)
    implicit none
    real(rp),     pointer, intent(in)    :: sendbuf(:,:)         !< Send buffer
    real(rp),     pointer, intent(inout) :: recvbuf(:,:)         !< Recv buffer
    integer(4),            intent(in)    :: recvcount4           !< Recv counts
    character(*),          intent(in)    :: wherein                !< Wherein
    integer(4)                           :: istat4
    integer(4)                           :: PAR_COMM_TO_USE

#ifndef MPI_OFF
    if( IPARALL ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       call MPI_ALLGATHER(sendbuf,recvcount4,MPI_DOUBLE_PRECISION,&
            &             recvbuf,recvcount4,MPI_DOUBLE_PRECISION,&
            &             PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0 ) call runend('PAR_ALLGATHER_IP: MPI ERROR')
    end if
#endif

  end subroutine PAR_ALLGATHER_RP_2

  subroutine PAR_ALLGATHER_RP_02(sendbuf,recvbuf,recvcount4,wherein,PAR_COMM_IN4)
    implicit none
    real(rp),              intent(in)    :: sendbuf(*)           !< Send buffer
    real(rp),     pointer, intent(inout) :: recvbuf(:,:)         !< Recv buffer
    integer(4),            intent(in)    :: recvcount4           !< Recv counts
    character(*), optional,intent(in)    :: wherein              !< Wherein
    integer(4),   optional,intent(in)    :: PAR_COMM_IN4         !< Communicator
    integer(4)                           :: istat4,lbounr,lbouns
    integer(4)                           :: PAR_COMM_TO_USE

    if( IPARALL ) then
#ifndef MPI_OFF
       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       end if
       call MPI_ALLGATHER(sendbuf,recvcount4,MPI_DOUBLE_PRECISION,&
            &             recvbuf,recvcount4,MPI_DOUBLE_PRECISION,&
            &             PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0 ) call runend('PAR_ALLGATHER_IP: MPI ERROR')
#endif
    else
       lbounr = lbound(recvbuf,2)
       recvbuf(1:recvcount4,lbounr) = sendbuf(1:recvcount4)
    end if

  end subroutine PAR_ALLGATHER_RP_02

  subroutine PAR_ALLGATHER_LG(sendbuf,recvbuf,recvcount4,wherein)
    implicit none
    logical(lg),           intent(in)  :: sendbuf              !< Send buffer
    logical(lg),  pointer, intent(inout) :: recvbuf(:)           !< Recv buffer
    integer(4),            intent(in)  :: recvcount4           !< Recv counts
    character(*),          intent(in)  :: wherein                !< Wherein
    integer(4)                         :: istat4
    integer(4)                         :: PAR_COMM_TO_USE
    integer(4)                         :: sendcount4

    if( IPARALL ) then
#ifndef MPI_OFF
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       sendcount4 = 1_4
       call MPI_ALLGATHER(sendbuf,sendcount4,MPI_LOGICAL,&
            &             recvbuf,recvcount4,MPI_LOGICAL,&
            &             PAR_COMM_TO_USE,istat4)
       if( istat4 /= 0 ) call runend('PAR_ALLGATHER_s: MPI ERROR')
#endif
    end if

  end subroutine PAR_ALLGATHER_LG

  !-----------------------------------------------------------------------
  !
  !> @brief   Initialize MPI
  !> @details Initialize MPI
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine PAR_INIT()
    implicit none
    integer(4) :: istat4

#ifndef MPI_OFF
    istat4 = 0
    call MPI_Init(istat4)
    if( istat4 /= 0_4 ) call runend('COULD NOT INITIALIZE MPI')
#endif
#ifdef EXTRAE
    call extrae_shutdown()
#endif

  end subroutine PAR_INIT

  !-----------------------------------------------------------------------
  !
  !> @brief   Define length of integers
  !> @details Define length of integers
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine PAR_LENGTH_INTEGER()
    implicit none

#ifndef MPI_OFF
    if( ip == 4 ) then
       PAR_INTEGER = MPI_INTEGER
    else
       PAR_INTEGER = MPI_INTEGER8
    end if
#endif

  end subroutine PAR_LENGTH_INTEGER

  !-----------------------------------------------------------------------
  !
  !> @brief   Send receive
  !> @details Send and receive to all my neghbors within communicator COMMU
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine PAR_SEND_RECEIVE_TO_ALL_IP_1(xx_send,xx_recv,wherein,wsynch,lsend_perm,PERMUTE_SEND,PERMUTE_RECV)
    implicit none
    integer(ip), pointer,           intent(in)    :: xx_send(:)
    integer(ip), pointer,           intent(inout) :: xx_recv(:)
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip), pointer, optional, intent(in)    :: lsend_perm(:)
    logical(lg),          optional, intent(in)    :: PERMUTE_SEND
    logical(lg),          optional, intent(in)    :: PERMUTE_RECV

    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    integer(4)                                    :: PAR_COMM_TO_USE
    integer(ip),                       pointer    :: xx_send_perm(:)
    integer(ip),                       pointer    :: xx_recv_perm(:)
    integer(ip)                                   :: icont
    logical(lg)                                   :: if_permute_send
    logical(lg)                                   :: if_permute_recv
    integer(ip)                                   :: xx_send_tmp(2)
    integer(ip)                                   :: xx_recv_tmp(2)

    if( ISEQUEN ) return
    ndofn = 1_ip
    if_permute_send = .false.
    if_permute_recv = .false.
    if( present(PERMUTE_SEND) .or. present(lsend_perm) ) if_permute_send = PERMUTE_SEND
    if( present(PERMUTE_RECV) ) if_permute_recv = PERMUTE_RECV
    nullify(xx_send_perm)
    nullify(xx_recv_perm)

    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)

    if( if_permute_send .and. commu % lsend_dim > 0 ) then
       call memory_alloca(par_memor,'xx_send_perm','par_send_receive_to_all_rp_0',xx_send_perm,max(1_ip,commu % lsend_dim))
       do icont = 1,commu % lsend_dim
          xx_send_perm(icont) = xx_send(commu % lsend_perm(icont))
       end do
       if( .not. associated(xx_recv) ) then
          call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send_perm,xx_recv_tmp,commu,wsynch)
       else
          call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send_perm,xx_recv,commu,wsynch)
       end if
    else
      if( .not. associated(xx_send) .and. .not. associated(xx_recv) ) then
          return
       else if( .not. associated(xx_send) ) then
          call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send_tmp,xx_recv,commu,wsynch)
       else if( .not. associated(xx_recv) ) then
          call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send,xx_recv_tmp,commu,wsynch)
       else
          call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send,xx_recv,commu,wsynch)
       end if
    endif

    if( if_permute_send .and. commu % lsend_dim > 0 ) then
       call memory_deallo(par_memor,'xx_send_perm','par_send_receive_to_all_rp_0',xx_send_perm)
    end if

    if( if_permute_recv .and. commu % lrecv_dim > 0  ) then
       call runend('PAR_SEND_RECEIVE_TO_ALL_RP_0: NOT CODED')
    end if

  end subroutine PAR_SEND_RECEIVE_TO_ALL_IP_1

  subroutine PAR_SEND_RECEIVE_TO_ALL_IP_1c(xx_send,xx_recv,commu,wsynch,lsend_perm,PERMUTE_SEND,PERMUTE_RECV)

    integer(ip), pointer,           intent(in)    :: xx_send(:)
    integer(ip), pointer,           intent(inout) :: xx_recv(:)
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip), pointer, optional, intent(in)    :: lsend_perm(:)
    logical(lg),          optional, intent(in)    :: PERMUTE_SEND
    logical(lg),          optional, intent(in)    :: PERMUTE_RECV

    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    integer(ip),                       pointer    :: xx_send_perm(:)
    integer(ip),                       pointer    :: xx_recv_perm(:)
    integer(ip)                                   :: icont,dom_i,ineig
    logical(lg)                                   :: if_permute_send
    logical(lg)                                   :: if_permute_recv
    integer(ip)                                   :: xx_send_tmp(2)
    integer(ip)                                   :: xx_recv_tmp(2)
    logical(lg)                                   :: alltoallv
    integer(ip)                                   :: PAR_CURRENT_RANK,PAR_CURRENT_SIZE
    integer(4),                        pointer    :: sendcount4(:)
    integer(4),                        pointer    :: recvcount4(:)

    if( ISEQUEN ) return
    ndofn = 1_ip
    if_permute_send = .false.
    if_permute_recv = .false.
    alltoallv       = .false.

    if( present(PERMUTE_SEND) .or. present(lsend_perm) ) if_permute_send = PERMUTE_SEND
    if( present(PERMUTE_RECV) ) if_permute_recv = PERMUTE_RECV
    nullify(xx_send_perm)
    nullify(xx_recv_perm)

    if( present(wsynch) ) then
       if( trim(wsynch) == 'ALLTOALLV' ) then
          alltoallv = .true.
       end if
    end if

    if( if_permute_send ) then
       call memory_alloca(par_memor,'xx_send_perm','par_send_receive_to_all_rp_0',xx_send_perm,max(1_ip,commu % lsend_dim))
       do icont = 1,commu % lsend_dim
          xx_send_perm(icont) = xx_send(commu % lsend_perm(icont))
       end do

       if( alltoallv ) then

          call PAR_COMM_RANK_AND_SIZE(COMMU % PAR_COMM_WORLD,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)
          allocate(sendcount4(0:PAR_CURRENT_SIZE-1))
          allocate(recvcount4(0:PAR_CURRENT_SIZE-1))
          sendcount4 = 0_4
          recvcount4 = 0_4
          do ineig = 1,commu % nneig
             dom_i             = commu % neights(ineig)
             sendcount4(dom_i) = commu % lsend_size(ineig+1) - commu % lsend_size(ineig)
             recvcount4(dom_i) = commu % lrecv_size(ineig+1) - commu % lrecv_size(ineig)
          end do
          call PAR_ALLTOALLV_IP_1(xx_send_perm,xx_recv,sendcount4,recvcount4,PAR_COMM_IN4=int(commu % PAR_COMM_WORLD,4))
          deallocate(sendcount4)
          deallocate(recvcount4)

       else

          if( commu % lsend_dim > 0 ) then
             if( .not. associated(xx_recv) ) then
                call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send_perm,xx_recv_tmp,commu,wsynch)
             else
                call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send_perm,xx_recv,commu,wsynch)
             end if
          else
             if( .not. associated(xx_recv) ) then
                call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send_tmp,xx_recv_tmp,commu,wsynch)
             else
                call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send_tmp,xx_recv,commu,wsynch)
             end if
          end if
       end if

    else
       if( .not. associated(xx_send) .and. .not. associated(xx_recv) ) then
          return
       else if( .not. associated(xx_send) ) then
          call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send_tmp,xx_recv,commu,wsynch)
       else if( .not. associated(xx_recv) ) then
          call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send,xx_recv_tmp,commu,wsynch)
       else
          call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send,xx_recv,commu,wsynch)
       end if
    endif

    if( if_permute_send .and. commu % lsend_dim > 0 ) then
       call memory_deallo(par_memor,'xx_send_perm','par_send_receive_to_all_rp_0',xx_send_perm)
    end if

    if( if_permute_recv .and. commu % lrecv_dim > 0  ) then
       call runend('PAR_SEND_RECEIVE_TO_ALL_RP_0: NOT CODED')
    end if

  end subroutine PAR_SEND_RECEIVE_TO_ALL_IP_1c

  subroutine PAR_SEND_RECEIVE_TO_ALL_IP_2(xx_send,xx_recv,commu,wsynch,lsend_perm,PERMUTE_SEND,PERMUTE_RECV)
    implicit none
    integer(ip),          pointer,  intent(in)    :: xx_send(:,:)
    integer(ip),          pointer,  intent(inout) :: xx_recv(:,:)
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip), pointer, optional, intent(in)    :: lsend_perm(:)
    logical(lg),          optional                :: PERMUTE_SEND
    logical(lg),          optional                :: PERMUTE_RECV
    integer(ip)                                   :: yy_recv(2)
    integer(ip)                                   :: yy_send(2)
    integer(ip),             pointer              :: xx_send_perm(:,:)
    integer(ip),             pointer              :: xx_recv_perm(:,:)
    integer(ip)                                   :: ndofn,icont,dom_i,ineig
    logical(lg)                                   :: if_permute_send
    logical(lg)                                   :: if_permute_recv
    logical(lg)                                   :: alltoallv
    integer(ip)                                   :: PAR_CURRENT_RANK,PAR_CURRENT_SIZE
    integer(4),                        pointer    :: sendcount4(:)
    integer(4),                        pointer    :: recvcount4(:)

    if( ISEQUEN ) return

    if(      associated(xx_send) ) then
       ndofn = size(xx_send,1)
    else if( associated(xx_recv) ) then
       ndofn = size(xx_recv,1)
    else
       return
    end if

    if_permute_send = .false.
    if_permute_recv = .false.
    alltoallv       = .false.

    if( present(PERMUTE_SEND) .or. present(lsend_perm) ) if_permute_send = PERMUTE_SEND
    if( present(PERMUTE_RECV) ) if_permute_recv = PERMUTE_RECV
    nullify(xx_send_perm)
    nullify(xx_recv_perm)

    if( present(wsynch) ) then
       if( trim(wsynch) == 'ALLTOALLV' ) then
          alltoallv = .true.
       end if
    end if

    if( alltoallv ) then

       call PAR_COMM_RANK_AND_SIZE(COMMU % PAR_COMM_WORLD,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)
       allocate(sendcount4(0:PAR_CURRENT_SIZE-1))
       allocate(recvcount4(0:PAR_CURRENT_SIZE-1))
       sendcount4 = 0_4
       recvcount4 = 0_4
       do ineig = 1,commu % nneig
          dom_i             = commu % neights(ineig)
          sendcount4(dom_i) = commu % lsend_size(ineig+1) - commu % lsend_size(ineig)
          recvcount4(dom_i) = commu % lrecv_size(ineig+1) - commu % lrecv_size(ineig)
       end do
       sendcount4 = sendcount4 * ndofn
       recvcount4 = recvcount4 * ndofn

    end if

    if( if_permute_recv .and. associated(xx_recv) ) then
       call memory_alloca(par_memor,'xx_recv_perm','par_send_receive_to_all_ip_1',xx_recv_perm,ndofn,max(1_ip,commu % lrecv_dim))
    else
       xx_recv_perm => xx_recv
    endif

    if( if_permute_send .and. associated(xx_send) ) then
       nullify(xx_send_perm)
       call memory_alloca(par_memor,'xx_send_perm','par_send_receive_to_all_ip_1',xx_send_perm,ndofn,max(1_ip,commu % lsend_dim))
       do icont= 1,commu % lsend_dim
          xx_send_perm(1:ndofn,icont) = xx_send(1:ndofn,commu % lsend_perm(icont))
       end do
       if( alltoallv ) then
          call PAR_ALLTOALLV_IP_2(xx_send_perm,xx_recv_perm,sendcount4,recvcount4,PAR_COMM_IN4=int(commu % PAR_COMM_WORLD,4))
       else
          if( .not. associated(xx_recv) ) then
             call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send_perm,yy_recv,commu,wsynch)
          else
             call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send_perm,xx_recv_perm,commu,wsynch)
          end if
       end if
       call memory_deallo(par_memor,'xx_send_perm','par_send_receive_to_all_ip_01',xx_send_perm)
    else
       if( alltoallv ) then
          call PAR_ALLTOALLV_IP_2(xx_send,xx_recv_perm,sendcount4,recvcount4,PAR_COMM_IN4=int(commu % PAR_COMM_WORLD,4))
       else
          if( .not. associated(xx_recv) ) then
             call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send,yy_recv,commu,wsynch)
          else if( .not. associated(xx_send) ) then
             call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,yy_send,xx_recv_perm,commu,wsynch)
          else
             call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send,xx_recv_perm,commu,wsynch)
          end if
       end if
    endif

    if( if_permute_recv .and. associated(xx_recv) ) then
       do icont= 1,commu % lrecv_dim
          xx_recv(1:ndofn,commu % lrecv_perm(icont)) = xx_recv_perm(1:ndofn,icont)
       end do
    end if

    if( alltoallv ) then
       deallocate(sendcount4)
       deallocate(recvcount4)
    end if

  end subroutine PAR_SEND_RECEIVE_TO_ALL_IP_2

  subroutine PAR_SEND_RECEIVE_TO_ALL_IP_0(ndofn,xx_send,xx_recv,commu,wsynch,lsend_perm,PERMUTE_SEND,PERMUTE_RECV)
    implicit none
    integer(ip),                    intent(in)  :: ndofn
    integer(ip),                    intent(in)  :: xx_send(ndofn,*)
    integer(ip),                    intent(out) :: xx_recv(ndofn,*)
    type(comm_data_par),            intent(in)  :: commu
    character(*),         optional, intent(in)  :: wsynch
    integer(ip), pointer, optional, intent(in)  :: lsend_perm(:)
    logical(lg),          optional              :: PERMUTE_SEND
    logical(lg),          optional              :: PERMUTE_RECV
    integer(ip),                       pointer  :: xx_send_perm(:,:)
    integer(ip),                       pointer  :: xx_recv_perm(:,:)
    integer(ip)                                 :: icont
    logical(lg)                                 :: if_permute_send
    logical(lg)                                 :: if_permute_recv

    if( ISEQUEN ) return

    if_permute_send = .false.
    if_permute_recv = .false.
    if( present(PERMUTE_SEND) .or. present(lsend_perm) ) if_permute_send = PERMUTE_SEND
    if( present(PERMUTE_RECV) ) if_permute_recv = PERMUTE_RECV
    nullify(xx_send_perm)
    nullify(xx_recv_perm)

    if( if_permute_send ) then
       nullify(xx_send_perm)
       call memory_alloca(par_memor,'xx_send_perm','par_send_receive_to_all_ip_0',xx_send_perm,ndofn,max(1_ip,commu % lsend_dim))
       do icont=1, commu % lsend_dim
          xx_send_perm(1:ndofn,icont) = xx_send(1:ndofn,commu % lsend_perm(icont))
       end do
       call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send_perm,xx_recv,commu,wsynch)
       call memory_deallo(par_memor,'xx_send_perm','par_send_receive_to_all_ip_0',xx_send_perm)
    else
       call PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send,xx_recv,commu,wsynch)
    end if

  end subroutine PAR_SEND_RECEIVE_TO_ALL_IP_0

  subroutine PAR_SEND_RECEIVE_TO_ALL_IP(ndofn,xx_send,xx_recv,commu,wsynch)
    implicit none
    integer(ip),                    intent(in)  :: ndofn
    integer(ip),                    intent(in)  :: xx_send(*)
    integer(ip),                    intent(out) :: xx_recv(*)
    type(comm_data_par),            intent(in)  :: commu
    character(*),         optional, intent(in)  :: wsynch
    integer(ip)                                 :: inise,finse,inire,finre,counti
    integer(ip)                                 :: dom_i,ineig,nsend,nrecv,kk
    integer(4)                                  :: istat4,nsend4,nrecv4,dom_i4
    integer(4)                                  :: PAR_COMM_TO_USE,count4
    logical(lg)                                 :: asynch
    integer(4),           allocatable           :: ireqq4(:)
    integer(4),           pointer               :: status4(:,:)
    !
    ! Define communicator
    !
    PAR_COMM_TO_USE = int(commu % PAR_COMM_WORLD,4)
    !
    ! Asynchronous
    !
    if( present(wsynch) ) then
       if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
          asynch = .false.
       else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
          asynch = .true.
       else
          call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
       end if
    else
       asynch = .false.
    end if
    if( asynch ) allocate(ireqq4(commu % nneig*2))
    kk = 0
    !
    ! Synchronous Send/receive
    !
    do ineig = 1,commu % nneig

       dom_i  = commu % neights(ineig)

       inise  = ndofn * ( commu % lsend_size(ineig)   - 1 ) + 1
       nsend  = ndofn * ( commu % lsend_size(ineig+1) - 1 ) + 1 - inise
       finse  = inise + nsend - 1
       inire  = ndofn * ( commu % lrecv_size(ineig)   - 1 ) + 1
       nrecv  = ndofn * ( commu % lrecv_size(ineig+1) - 1 ) + 1 - inire
       finre  = inire + nrecv - 1

       nsend4 = int(nsend,4)
       nrecv4 = int(nrecv,4)
       dom_i4 = int(dom_i,4)

#ifndef MPI_OFF
       if( asynch .and. commu % nneig /= 0_ip ) then
          if( nsend > 0_ip ) then
             kk = kk + 1
             call MPI_Isend(                          &
                  xx_send(inise:finse), nsend4,       &
                  PAR_INTEGER,  dom_i4, 0_4,          &
                  PAR_COMM_TO_USE, ireqq4(kk), istat4 )
          end if
          if( nrecv > 0_ip ) then
             kk = kk + 1
             call MPI_Irecv(                          &
                  xx_recv(inire:finre), nrecv4,       &
                  PAR_INTEGER,  dom_i4, 0_4,          &
                  PAR_COMM_TO_USE, ireqq4(kk), istat4 )
          end if
       else
          if( nrecv /= 0 .and. nsend == 0 ) then
             call MPI_Recv(                       &
                  xx_recv(inire:finre), nrecv4,   &
                  PAR_INTEGER, dom_i4, 0_4,       &
                  PAR_COMM_TO_USE, status, istat4 )

          else if( nrecv == 0 .and. nsend /= 0 ) then
             call MPI_Send(                       &
                  xx_send(inise:finse), nsend4,   &
                  PAR_INTEGER, dom_i4, 0_4,       &
                  PAR_COMM_TO_USE, istat4         )

          else if( nrecv /= 0 .and. nsend /= 0 ) then
             call MPI_Sendrecv(                   &
                  xx_send(inise:finse), nsend4,   &
                  PAR_INTEGER, dom_i4, 0_4,       &
                  xx_recv(inire:finre), nrecv4,   &
                  PAR_INTEGER, dom_i4, 0_4,       &
                  PAR_COMM_TO_USE, status, istat4 )
          end if
       end if
#endif

    end do
    !
    ! Wait in case of asynchronous, count is the total of communications actually performe
    ! and is smaller or equal than the asynchronous communications requested (ireqq)
    !
#ifndef MPI_OFF
    if( asynch .and. commu % nneig /= 0 ) then
       counti = kk
       count4 = int(counti,4)
       allocate( status4(MPI_STATUS_SIZE,counti) )
       CALL MPI_WAITALL(count4,ireqq4,status4,istat4)
       deallocate( status4 )
       deallocate( ireqq4  )
    end if
#endif

  end subroutine PAR_SEND_RECEIVE_TO_ALL_IP

  subroutine PAR_SEND_RECEIVE_TO_ALL_RP_s(xx_send,xx_recv,commu,wsynch)
    implicit none
    real(rp),                       intent(in)  :: xx_send
    real(rp),                       intent(out) :: xx_recv
    type(comm_data_par),            intent(in)  :: commu
    character(*),         optional, intent(in)  :: wsynch
    real(rp)                                    :: yy_send(2)
    real(rp)                                    :: yy_recv(2)
    if( ISEQUEN ) return
    yy_send(1) = xx_send
    call PAR_SEND_RECEIVE_TO_ALL_RP(1_ip,yy_send,yy_recv,commu,wsynch)
    xx_recv    = yy_recv(1)
  end subroutine PAR_SEND_RECEIVE_TO_ALL_RP_s

 subroutine PAR_SEND_RECEIVE_TO_ALL_RP_0(ndofn,xx_send,xx_recv,commu,wsynch,lsend_perm,PERMUTE_SEND,PERMUTE_RECV)
    implicit none
    integer(ip),                    intent(in)  :: ndofn
    real(rp),                       intent(in)  :: xx_send(ndofn,*)
    real(rp),                       intent(out) :: xx_recv(ndofn,*)
    type(comm_data_par),            intent(in)  :: commu
    integer(ip), pointer, optional, intent(in)  :: lsend_perm(:)
    logical(lg),          optional, intent(in)  :: PERMUTE_SEND
    logical(lg),          optional, intent(in)  :: PERMUTE_RECV
    character(*),         optional, intent(in)  :: wsynch
    real(rp),                          pointer  :: xx_send_perm(:,:)
    real(rp),                          pointer  :: xx_recv_perm(:,:)
    integer(ip)                                 :: icont,idofn,itotn,jcont,jtotn
    logical(lg)                                 :: if_permute_send
    logical(lg)                                 :: if_permute_recv

    if( ISEQUEN ) return

    if_permute_send = .false.
    if_permute_recv = .false.
    if( present(PERMUTE_SEND) .or. present(lsend_perm) ) if_permute_send = PERMUTE_SEND
    if( present(PERMUTE_RECV) ) if_permute_recv = PERMUTE_RECV
    nullify(xx_send_perm)
    nullify(xx_recv_perm)

    if( if_permute_send .and. commu % lsend_dim > 0 ) then
       call memory_alloca(par_memor,'xx_send_perm','par_send_receive_to_all_rp_0',xx_send_perm,ndofn,max(1_ip,commu % lsend_dim))
       do icont = 1,commu % lsend_dim
          xx_send_perm(1:ndofn,icont) = xx_send(1:ndofn,commu % lsend_perm(icont))
       end do
       call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send_perm,xx_recv,commu,wsynch)
    else
       call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send,xx_recv,commu,wsynch)
    endif

    if( if_permute_send .and. commu % lsend_dim > 0 ) then
       call memory_deallo(par_memor,'xx_send_perm','par_send_receive_to_all_rp_0',xx_send_perm)
    end if

    if( if_permute_recv .and. commu % lrecv_dim > 0  ) then
       call runend('PAR_SEND_RECEIVE_TO_ALL_RP_0: NOT CODED')
    end if

 end subroutine PAR_SEND_RECEIVE_TO_ALL_RP_0

  subroutine PAR_SEND_RECEIVE_TO_ALL_RP_1(xx_send,xx_recv,wherein,wsynch,lsend_perm,PERMUTE_SEND,PERMUTE_RECV)
    implicit none
    real(rp),    pointer,           intent(in)    :: xx_send(:)
    real(rp),    pointer,           intent(inout) :: xx_recv(:)
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip), pointer, optional, intent(in)    :: lsend_perm(:)
    logical(lg),          optional, intent(in)    :: PERMUTE_SEND
    logical(lg),          optional, intent(in)    :: PERMUTE_RECV

    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    integer(4)                                    :: PAR_COMM_TO_USE
    real(rp),                          pointer    :: xx_send_perm(:)
    real(rp),                          pointer    :: xx_recv_perm(:)
    integer(ip)                                   :: icont
    logical(lg)                                   :: if_permute_send
    logical(lg)                                   :: if_permute_recv
    real(rp)                                      :: xx_send_tmp(2)
    real(rp)                                      :: xx_recv_tmp(2)

    if( ISEQUEN ) return
    ndofn = 1_ip
    if_permute_send = .false.
    if_permute_recv = .false.
    if( present(PERMUTE_SEND) .or. present(lsend_perm) ) if_permute_send = PERMUTE_SEND
    if( present(PERMUTE_RECV) ) if_permute_recv = PERMUTE_RECV
    nullify(xx_send_perm)
    nullify(xx_recv_perm)

    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)

    if( if_permute_send .and. commu % lsend_dim > 0 ) then
       call memory_alloca(par_memor,'xx_send_perm','par_send_receive_to_all_rp_0',xx_send_perm,max(1_ip,commu % lsend_dim))
       do icont = 1,commu % lsend_dim
          xx_send_perm(icont) = xx_send(commu % lsend_perm(icont))
       end do
       if( .not. associated(xx_recv) ) then
          call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send_perm,xx_recv_tmp,commu,wsynch)
       else
          call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send_perm,xx_recv,commu,wsynch)
       end if
    else
      if( .not. associated(xx_send) .and. .not. associated(xx_recv) ) then
          return
       else if( .not. associated(xx_send) ) then
          call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send_tmp,xx_recv,commu,wsynch)
       else if( .not. associated(xx_recv) ) then
          call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send,xx_recv_tmp,commu,wsynch)
       else
          call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send,xx_recv,commu,wsynch)
       end if
    endif

    if( if_permute_send .and. commu % lsend_dim > 0 ) then
       call memory_deallo(par_memor,'xx_send_perm','par_send_receive_to_all_rp_0',xx_send_perm)
    end if

    if( if_permute_recv .and. commu % lrecv_dim > 0  ) then
       call runend('PAR_SEND_RECEIVE_TO_ALL_RP_0: NOT CODED')
    end if

  end subroutine PAR_SEND_RECEIVE_TO_ALL_RP_1

  subroutine PAR_SEND_RECEIVE_TO_ALL_RP_1c(xx_send,xx_recv,commu,wsynch,lsend_perm,PERMUTE_SEND,PERMUTE_RECV)

    real(rp),    pointer,           intent(in)    :: xx_send(:)
    real(rp),    pointer,           intent(inout) :: xx_recv(:)
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip), pointer, optional, intent(in)    :: lsend_perm(:)
    logical(lg),          optional, intent(in)    :: PERMUTE_SEND
    logical(lg),          optional, intent(in)    :: PERMUTE_RECV

    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    real(rp),                          pointer    :: xx_send_perm(:)
    real(rp),                          pointer    :: xx_recv_perm(:)
    integer(ip)                                   :: icont,dom_i,ineig
    logical(lg)                                   :: if_permute_send
    logical(lg)                                   :: if_permute_recv
    real(rp)                                      :: xx_send_tmp(2)
    real(rp)                                      :: xx_recv_tmp(2)
    logical(lg)                                   :: alltoallv
    integer(ip)                                   :: PAR_CURRENT_RANK,PAR_CURRENT_SIZE
    integer(4),                        pointer    :: sendcount4(:)
    integer(4),                        pointer    :: recvcount4(:)

    if( ISEQUEN ) return
    ndofn = 1_ip
    if_permute_send = .false.
    if_permute_recv = .false.
    alltoallv       = .false.

    if( present(PERMUTE_SEND) .or. present(lsend_perm) ) if_permute_send = PERMUTE_SEND
    if( present(PERMUTE_RECV) ) if_permute_recv = PERMUTE_RECV
    nullify(xx_send_perm)
    nullify(xx_recv_perm)

    if( present(wsynch) ) then
       if( trim(wsynch) == 'ALLTOALLV' ) then
          alltoallv = .true.
       end if
    end if

    if( alltoallv ) then

       call PAR_COMM_RANK_AND_SIZE(COMMU % PAR_COMM_WORLD,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)
       allocate(sendcount4(0:PAR_CURRENT_SIZE-1))
       allocate(recvcount4(0:PAR_CURRENT_SIZE-1))
       sendcount4 = 0_4
       recvcount4 = 0_4
       do ineig = 1,commu % nneig
          dom_i             = commu % neights(ineig)
          sendcount4(dom_i) = commu % lsend_size(ineig+1) - commu % lsend_size(ineig)
          recvcount4(dom_i) = commu % lrecv_size(ineig+1) - commu % lrecv_size(ineig)
       end do

    end if

    if( if_permute_send .and. commu % lsend_dim > 0 ) then
       call memory_alloca(par_memor,'xx_send_perm','par_send_receive_to_all_rp_0',xx_send_perm,max(1_ip,commu % lsend_dim))
       do icont = 1,commu % lsend_dim
          xx_send_perm(icont) = xx_send(commu % lsend_perm(icont))
       end do
       if( alltoallv ) then
          call PAR_ALLTOALLV_RP_1(xx_send_perm,xx_recv,sendcount4,recvcount4,PAR_COMM_IN4=int(commu % PAR_COMM_WORLD,4))
       else
          if( .not. associated(xx_recv) ) then
             call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send_perm,xx_recv_tmp,commu,wsynch)
          else
             call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send_perm,xx_recv,commu,wsynch)
          end if
       end if
    else
       if( .not. associated(xx_send) .and. .not. associated(xx_recv) ) then
          return
       else if( .not. associated(xx_send) ) then
          call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send_tmp,xx_recv,commu,wsynch)
       else if( .not. associated(xx_recv) ) then
          call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send,xx_recv_tmp,commu,wsynch)
       else
          call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send,xx_recv,commu,wsynch)
       end if
    endif

    if( if_permute_send .and. commu % lsend_dim > 0 ) then
       call memory_deallo(par_memor,'xx_send_perm','par_send_receive_to_all_rp_0',xx_send_perm)
    end if

    if( if_permute_recv .and. commu % lrecv_dim > 0  ) then
       call runend('PAR_SEND_RECEIVE_TO_ALL_RP_0: NOT CODED')
    end if

    if( alltoallv ) then
       deallocate(sendcount4)
       deallocate(recvcount4)
    end if

  end subroutine PAR_SEND_RECEIVE_TO_ALL_RP_1c

  subroutine PAR_SEND_RECEIVE_TO_ALL_RP_2(xx_send,xx_recv,commu,wsynch,lsend_perm,PERMUTE_SEND,PERMUTE_RECV)
    implicit none
    real(rp),             pointer,  intent(in)    :: xx_send(:,:)
    real(rp),             pointer,  intent(inout) :: xx_recv(:,:)
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip), pointer, optional, intent(in)    :: lsend_perm(:)
    logical(lg),          optional                :: PERMUTE_SEND
    logical(lg),          optional                :: PERMUTE_RECV
    integer(ip)                                   :: ndofn,icont,dom_i,ineig
    real(rp)                                      :: yy_send(2)
    real(rp)                                      :: yy_recv(2)
    real(rp),             pointer                 :: xx_send_perm(:,:)
    real(rp),             pointer                 :: xx_recv_perm(:,:)
    logical(lg)                                   :: if_permute_send
    logical(lg)                                   :: if_permute_recv
    real(rp)                                      :: xx_send_tmp(1,1)
    real(rp)                                      :: xx_recv_tmp(1,1)
    logical(lg)                                   :: alltoallv
    integer(ip)                                   :: PAR_CURRENT_RANK,PAR_CURRENT_SIZE
    integer(4),                        pointer    :: sendcount4(:)
    integer(4),                        pointer    :: recvcount4(:)

    if_permute_send = .false.
    if_permute_recv = .false.
    alltoallv       = .false.

    if( present(PERMUTE_SEND) .or. present(lsend_perm) ) if_permute_send = PERMUTE_SEND
    if( present(PERMUTE_RECV) ) if_permute_recv = PERMUTE_RECV
    nullify(xx_send_perm)
    nullify(xx_recv_perm)

    if( present(wsynch) ) then
       if( trim(wsynch) == 'ALLTOALLV' ) then
          alltoallv = .true.
       end if
    end if

    if( if_permute_send .and. associated(xx_send) ) then
       ndofn = size(xx_send,1)
       call memory_alloca(par_memor,'xx_send_perm','par_send_receive_to_all_rp_0',xx_send_perm,ndofn,max(1_ip,commu % lsend_dim))
       do icont = 1,commu % lsend_dim
          xx_send_perm(1:ndofn,icont) = xx_send(1:ndofn,commu % lsend_perm(icont))
       end do
    else
       xx_send_perm => xx_send
    end if

    if( if_permute_recv ) then
       call memory_alloca(par_memor,'xx_recv_perm','par_send_receive_to_all_rp_0',xx_recv_perm,ndofn,max(1_ip,commu % lrecv_dim))
   else
       xx_recv_perm => xx_recv
    end if

    if( .not. if_permute_send ) then
       if( associated(xx_send_perm) ) then
          if( size(xx_send_perm,2) < npoin ) call runend('WRONG SIZE: CANNOT INTERPOLATE')
       end if
    end if

    if( alltoallv ) then

       call PAR_COMM_RANK_AND_SIZE(COMMU % PAR_COMM_WORLD,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)
       allocate(sendcount4(0:PAR_CURRENT_SIZE-1))
       allocate(recvcount4(0:PAR_CURRENT_SIZE-1))
       sendcount4 = 0_4
       recvcount4 = 0_4
       do ineig = 1,commu % nneig
          dom_i             = commu % neights(ineig)
          sendcount4(dom_i) = commu % lsend_size(ineig+1) - commu % lsend_size(ineig)
          recvcount4(dom_i) = commu % lrecv_size(ineig+1) - commu % lrecv_size(ineig)
       end do

    end if

    if( .not. associated(xx_send_perm) .and. associated(xx_recv_perm) ) then
       ndofn = size(xx_recv_perm,1)
       if( alltoallv ) then
          sendcount4 = sendcount4 * ndofn
          recvcount4 = recvcount4 * ndofn
          call PAR_ALLTOALLV_RP_2(xx_send_perm,xx_recv_perm,sendcount4,recvcount4,PAR_COMM_IN4=int(commu % PAR_COMM_WORLD,4))
       else
          call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,yy_send,xx_recv_perm,commu,wsynch)
       end if

    else if( .not. associated(xx_send_perm) .and. .not. associated(xx_recv_perm) ) then
       continue

    else if( associated(xx_send_perm) .and. .not. associated(xx_recv_perm) ) then
       ndofn = size(xx_send_perm,1)
       if( alltoallv ) then
          sendcount4 = sendcount4 * ndofn
          recvcount4 = recvcount4 * ndofn
          call PAR_ALLTOALLV_RP_2(xx_send_perm,xx_recv_perm,sendcount4,recvcount4,PAR_COMM_IN4=int(commu % PAR_COMM_WORLD,4))
       else
          call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send_perm,yy_recv,commu,wsynch)
       end if

    else
       ndofn = size(xx_send_perm,1)
       if( ndofn /= size(xx_recv_perm,1) ) call runend('WRONG SIZE: CANNOT INTERPOLATE')
       if( alltoallv ) then
          sendcount4 = sendcount4 * ndofn
          recvcount4 = recvcount4 * ndofn
          call PAR_ALLTOALLV_RP_2(xx_send_perm,xx_recv_perm,sendcount4,recvcount4,PAR_COMM_IN4=int(commu % PAR_COMM_WORLD,4))
       else
          call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send_perm,xx_recv_perm,commu,wsynch)
       end if

    end if

   if( if_permute_send ) then
      call memory_deallo(par_memor,'xx_send_perm','par_send_receive_to_all_rp_0',xx_send_perm)
   end if

   if( if_permute_recv .and. associated(xx_recv) ) then
      do icont = 1,commu % lrecv_dim
         xx_recv(1:ndofn,commu % lrecv_perm(icont)) = xx_recv_perm(1:ndofn,icont)
      end do
      call memory_deallo(par_memor,'xx_send_perm','par_send_receive_to_all_rp_0',xx_recv_perm)
   end if

   if( alltoallv ) then
      deallocate(sendcount4)
      deallocate(recvcount4)
   end if

  end subroutine PAR_SEND_RECEIVE_TO_ALL_RP_2

  subroutine PAR_SEND_RECEIVE_TO_ALL_RP_3(xx_send,xx_recv,commu,wsynch,PERMUTE_SEND)
    implicit none
    real(rp),             pointer,  intent(in)    :: xx_send(:,:,:)
    real(rp),             pointer,  intent(inout) :: xx_recv(:,:,:)
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    logical(lg),          optional                :: PERMUTE_SEND
    integer(ip)                                   :: ndofn,ndof1,ndof2,icont,dom_i,ineig
    real(rp)                                      :: yy_send(2)
    real(rp)                                      :: yy_recv(2)
    real(rp),             pointer                 :: xx_send_perm(:,:,:)
    logical(lg)                                   :: if_permute_send
    logical(lg)                                   :: alltoallv
    integer(ip)                                   :: PAR_CURRENT_RANK,PAR_CURRENT_SIZE
    integer(4),                        pointer    :: sendcount4(:)
    integer(4),                        pointer    :: recvcount4(:)

    if_permute_send = .false.
    if( present(PERMUTE_SEND) ) if_permute_send = PERMUTE_SEND
    alltoallv = .false.

    if( present(wsynch) ) then
       if( trim(wsynch) == 'ALLTOALLV' ) then
          alltoallv = .true.
       end if
    end if

    if( alltoallv ) then

       call PAR_COMM_RANK_AND_SIZE(COMMU % PAR_COMM_WORLD,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)
       allocate(sendcount4(0:PAR_CURRENT_SIZE-1))
       allocate(recvcount4(0:PAR_CURRENT_SIZE-1))
       sendcount4 = 0_4
       recvcount4 = 0_4
       do ineig = 1,commu % nneig
          dom_i             = commu % neights(ineig)
          sendcount4(dom_i) = commu % lsend_size(ineig+1) - commu % lsend_size(ineig)
          recvcount4(dom_i) = commu % lrecv_size(ineig+1) - commu % lrecv_size(ineig)
       end do

    end if


    if( if_permute_send ) then

       nullify(xx_send_perm)
       if(      associated(xx_send) ) then
          ndof1 = size(xx_send,1)
          ndof2 = size(xx_send,2)
       else if( associated(xx_recv) ) then
          ndof1 = size(xx_recv,1)
          ndof2 = size(xx_recv,2)
       else
          return
       end if
       ndofn = ndof1*ndof2

       call memory_alloca(par_memor,'xx_send_perm','par_send_receive_to_all_rp_0',xx_send_perm,ndof1,ndof2,max(1_ip,commu % lsend_dim))
        do icont = 1, commu % lsend_dim
          xx_send_perm(1:ndof1,1:ndof2,icont) = xx_send(1:ndof1,1:ndof2,commu % lsend_perm(icont))
       end do

       if( alltoallv ) then
          sendcount4 = sendcount4 * ndofn
          recvcount4 = recvcount4 * ndofn
          call PAR_ALLTOALLV_RP_3(xx_send_perm,xx_recv,sendcount4,recvcount4,PAR_COMM_IN4=int(commu % PAR_COMM_WORLD,4))
       else
          if( .not. associated(xx_recv) ) then
             call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send_perm,yy_recv,commu,wsynch)
          else
             call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send_perm,xx_recv,commu,wsynch)
          end if
       end if
       call memory_deallo(par_memor,'xx_send_perm','par_send_receive_to_all_rp_0',xx_send_perm)

    else

       if( .not. associated(xx_send) .and. associated(xx_recv) ) then
          ndofn = size(xx_recv,1)*size(xx_recv,2)
          if( alltoallv ) then
             sendcount4 = sendcount4 * ndofn
             recvcount4 = recvcount4 * ndofn
             call PAR_ALLTOALLV_RP_3(xx_send,xx_recv,sendcount4,recvcount4,PAR_COMM_IN4=int(commu % PAR_COMM_WORLD,4))
          else
             call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,yy_send,xx_recv,commu,wsynch)
          end if
       else if( .not. associated(xx_send) .and. .not. associated(xx_recv) ) then
          ndofn = 1
          if( alltoallv ) then
             sendcount4 = sendcount4
             recvcount4 = recvcount4
             call PAR_ALLTOALLV_RP_3(xx_send,xx_recv,sendcount4,recvcount4,PAR_COMM_IN4=int(commu % PAR_COMM_WORLD,4))
          else
             call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,yy_send,yy_recv,commu,wsynch)
          end if
       else if( associated(xx_send) .and. .not. associated(xx_recv) ) then
          ndofn = size(xx_send,1)* size(xx_send,2)
          if( alltoallv ) then
             sendcount4 = sendcount4 * ndofn
             recvcount4 = recvcount4 * ndofn
             call PAR_ALLTOALLV_RP_3(xx_send,xx_recv,sendcount4,recvcount4,PAR_COMM_IN4=int(commu % PAR_COMM_WORLD,4))
          else
             call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send,yy_recv,commu,wsynch)
          end if
       else
          ndofn = size(xx_send,1)*size(xx_send,2)
          if( ndofn /= size(xx_recv,1)*size(xx_recv,2) ) call runend('WRONG SIZE: CANNOT INTERPOLATE')
          if( alltoallv ) then
             sendcount4 = sendcount4 * ndofn
             recvcount4 = recvcount4 * ndofn
             call PAR_ALLTOALLV_RP_3(xx_send,xx_recv,sendcount4,recvcount4,PAR_COMM_IN4=int(commu % PAR_COMM_WORLD,4))
          else
             call PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send,xx_recv,commu,wsynch)
          end if
       end if
    end if

  end subroutine PAR_SEND_RECEIVE_TO_ALL_RP_3

  subroutine PAR_SEND_RECEIVE_TO_ALL_RP(ndofn,xx_send,xx_recv,commu,wsynch)
    implicit none
    integer(ip),                    intent(in)  :: ndofn
    real(rp),                       intent(in)  :: xx_send(*)
    real(rp),                       intent(out) :: xx_recv(*)
    type(comm_data_par),            intent(in)  :: commu
    character(*),         optional, intent(in)  :: wsynch
    integer(ip)                                 :: inise,finse,inire,finre,counti
    integer(ip)                                 :: dom_i,ineig,nsend,nrecv,kk
    integer(4)                                  :: istat4,nsend4,nrecv4,dom_i4
    integer(4)                                  :: PAR_COMM_TO_USE,count4
    logical(lg)                                 :: asynch
    integer(4),           allocatable           :: ireqq4(:)
    integer(4),           pointer               :: status4(:,:)
    !
    ! Define communicator
    !
    PAR_COMM_TO_USE = int(commu % PAR_COMM_WORLD,4)
    !
    ! Asynchronous
    !
    if( present(wsynch) ) then
       if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
          asynch = .false.
       else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
          asynch = .true.
       else
          call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
       end if
    else
       asynch = .false.
    end if
    if( asynch .and. commu % nneig /= 0_ip ) allocate(ireqq4(commu % nneig*2_ip))

    kk = 0
    istat4 = 0_4
    !
    ! Synchronous Send/receive
    !
    do ineig = 1,commu % nneig

       dom_i  = commu % neights(ineig)
       inise  = ndofn * ( commu % lsend_size(ineig)   - 1 ) + 1
       nsend  = ndofn * ( commu % lsend_size(ineig+1) - 1 ) + 1 - inise
       finse  = inise + nsend - 1
       inire  = ndofn * ( commu % lrecv_size(ineig)   - 1 ) + 1
       nrecv  = ndofn * ( commu % lrecv_size(ineig+1) - 1 ) + 1 - inire
       finre  = inire + nrecv - 1
       nsend4 = int(nsend,4)
       nrecv4 = int(nrecv,4)
       dom_i4 = int(dom_i,4)

#ifndef MPI_OFF
       if( asynch .and. commu % nneig /= 0_ip ) then

          if( nsend > 0 ) then
             kk = kk + 1
             call MPI_Isend(                            &
                  xx_send(inise:finse), nsend4,         &
                  MPI_DOUBLE_PRECISION,  dom_i4, 0_4,   &
                  PAR_COMM_TO_USE, ireqq4(kk), istat4   )
          end if
          if( nrecv > 0 ) then
             kk = kk + 1
             call MPI_Irecv(                            &
                  xx_recv(inire:finre), nrecv4,         &
                  MPI_DOUBLE_PRECISION,  dom_i4, 0_4,   &
                  PAR_COMM_TO_USE, ireqq4(kk), istat4   )
          end if
       else
          if( nrecv /= 0 .and. nsend == 0 ) then
             call MPI_Recv(                          &
                  xx_recv(inire:finre), nrecv4,      &
                  MPI_DOUBLE_PRECISION, dom_i4, 0_4, &
                  PAR_COMM_TO_USE, status, istat4    )
          else if( nrecv == 0 .and. nsend /= 0 ) then
             call MPI_Send(                          &
                  xx_send(inise:finse), nsend4,      &
                  MPI_DOUBLE_PRECISION, dom_i4, 0_4, &
                  PAR_COMM_TO_USE, istat4            )
          else if( nrecv /= 0 .and. nsend /= 0 ) then
             call MPI_Sendrecv(                      &
                  xx_send(inise:finse), nsend4,      &
                  MPI_DOUBLE_PRECISION, dom_i4, 0_4, &
                  xx_recv(inire:finre), nrecv4,      &
                  MPI_DOUBLE_PRECISION, dom_i4, 0_4, &
                  PAR_COMM_TO_USE, status, istat4    )
          end if
       end if
#endif
    end do
    !
    ! Wait in case of asynchronous, count is the total of communications actually performe
    ! and is smaller or equal than the asynchronous communications requested (ireqq)
    !
#ifndef MPI_OFF
    if( asynch .and. commu % nneig /= 0_ip ) then
       counti = kk
       count4 = int(counti,4)
       allocate( status4(MPI_STATUS_SIZE,counti) )
       CALL MPI_WAITALL(count4,ireqq4,status4,istat4)
       deallocate( status4 )
       deallocate( ireqq4  )
    end if
#endif

  end subroutine PAR_SEND_RECEIVE_TO_ALL_RP

  !-----------------------------------------------------------------------
  !
  !> @brief   Bridge to broadcast data
  !> @details Bridge to broadcast data
  !> @author  Guillaume Houzeaux
  !
  !-----------------------------------------------------------------------

  subroutine PAR_EXCHANGE_IP_s(xx,xarra,icoun,ipass)
    implicit none
    integer(ip), intent(inout)          :: xx
    integer(ip), intent(inout), pointer :: xarra(:)
    integer(ip), intent(out)            :: icoun
    integer(ip), intent(in)             :: ipass
    icoun = icoun + 1
    if( ipass == 2 ) then
       if( IMASTER ) xarra(icoun) = xx
       if( ISLAVE  ) xx           = xarra(icoun)
    end if
  end subroutine PAR_EXCHANGE_IP_s
  subroutine PAR_EXCHANGE_IP_0(n,xx,xarra,icoun,ipass)
    integer(ip), intent(in)             :: n
    integer(ip), intent(inout)          :: xx(:)
    integer(ip), intent(inout), pointer :: xarra(:)
    integer(ip), intent(out)            :: icoun
    integer(ip), intent(in)             :: ipass
    integer(ip)                         :: ii
    do ii = 1,n
       icoun = icoun + 1
       if( ipass == 2 ) then
          if( IMASTER ) xarra(icoun) = xx(ii)
          if( ISLAVE  ) xx(ii)       = xarra(icoun)
       end if
    end do
  end subroutine PAR_EXCHANGE_IP_0
  subroutine PAR_EXCHANGE_IP_1(xx,xarra,icoun,ipass)
    implicit none
    integer(ip), intent(inout), pointer :: xx(:)
    integer(ip), intent(inout), pointer :: xarra(:)
    integer(ip), intent(out)            :: icoun
    integer(ip), intent(in)             :: ipass
    integer(ip)                         :: ii
    do ii = 1,size(xx)
       icoun = icoun + 1
       if( ipass == 2 ) then
          if( IMASTER ) xarra(icoun) = xx(ii)
          if( ISLAVE  ) xx(ii)       = xarra(icoun)
       end if
    end do
  end subroutine PAR_EXCHANGE_IP_1
  subroutine PAR_EXCHANGE_IP_2(xx,xarra,icoun,ipass)
    implicit none
    integer(ip), intent(inout), pointer :: xx(:,:)
    integer(ip), intent(inout), pointer :: xarra(:)
    integer(ip), intent(out)            :: icoun
    integer(ip), intent(in)             :: ipass
    integer(ip)                         :: ii,jj
    do jj = 1,size(xx,2)
       do ii = 1,size(xx,1)
          icoun = icoun + 1
          if( ipass == 2 ) then
             if( IMASTER ) xarra(icoun) = xx(ii,jj)
             if( ISLAVE  ) xx(ii,jj)    = xarra(icoun)
          end if
       end do
    end do
  end subroutine PAR_EXCHANGE_IP_2

  subroutine PAR_EXCHANGE_RP_s(xx,xarra,icoun,ipass)
    implicit none
    real(rp),    intent(inout)          :: xx
    real(rp),    intent(inout), pointer :: xarra(:)
    integer(ip), intent(out)            :: icoun
    integer(ip), intent(in)             :: ipass
    icoun = icoun + 1
    if( ipass == 2 ) then
       if( IMASTER ) xarra(icoun) = xx
       if( ISLAVE  ) xx           = xarra(icoun)
    end if
  end subroutine PAR_EXCHANGE_RP_s
  subroutine PAR_EXCHANGE_RP_0(n,xx,xarra,icoun,ipass)
    integer(ip), intent(in)             :: n
    real(rp),    intent(inout)          :: xx(:)
    real(rp),    intent(inout), pointer :: xarra(:)
    integer(ip), intent(out)            :: icoun
    integer(ip), intent(in)             :: ipass
    integer(ip)                         :: ii
    do ii = 1,n
       icoun = icoun + 1
       if( ipass == 2 ) then
          if( IMASTER ) xarra(icoun) = xx(ii)
          if( ISLAVE  ) xx(ii)       = xarra(icoun)
       end if
    end do
  end subroutine PAR_EXCHANGE_RP_0
  subroutine PAR_EXCHANGE_RP_02(n1,n2,xx,xarra,icoun,ipass)
    integer(ip), intent(in)             :: n1
    integer(ip), intent(in)             :: n2
    real(rp),    intent(inout)          :: xx(:,:)
    real(rp),    intent(inout), pointer :: xarra(:)
    integer(ip), intent(out)            :: icoun
    integer(ip), intent(in)             :: ipass
    integer(ip)                         :: i1,i2
    do i2 = 1,n2
       do i1 = 1,n1
          icoun = icoun + 1
          if( ipass == 2 ) then
             if( IMASTER ) xarra(icoun) = xx(i1,i2)
             if( ISLAVE  ) xx(i1,i2)    = xarra(icoun)
          end if
       end do
    end do
  end subroutine PAR_EXCHANGE_RP_02
  subroutine PAR_EXCHANGE_RP_1(xx,xarra,icoun,ipass)
    implicit none
    real(rp),    intent(inout), pointer :: xx(:)
    real(rp),    intent(inout), pointer :: xarra(:)
    integer(ip), intent(out)            :: icoun
    integer(ip), intent(in)             :: ipass
    integer(ip)                         :: ii
    do ii = 1,size(xx)
       icoun = icoun + 1
       if( ipass == 2 ) then
          if( IMASTER ) xarra(icoun) = xx(ii)
          if( ISLAVE  ) xx(ii)       = xarra(icoun)
       end if
    end do
  end subroutine PAR_EXCHANGE_RP_1
  subroutine PAR_EXCHANGE_RP_2(xx,xarra,icoun,ipass)
    real(rp),    intent(inout), pointer :: xx(:,:)
    real(rp),    intent(inout), pointer :: xarra(:)
    integer(ip), intent(out)            :: icoun
    integer(ip), intent(in)             :: ipass
    integer(ip)                         :: ii,jj
    do jj = 1,size(xx,2)
       do ii = 1,size(xx,1)
          icoun = icoun + 1
          if( ipass == 2 ) then
             if( IMASTER ) xarra(icoun) = xx(ii,jj)
             if( ISLAVE  ) xx(ii,jj)    = xarra(icoun)
          end if
       end do
    end do
  end subroutine PAR_EXCHANGE_RP_2

  subroutine PAR_EXCHANGE_CH(xx,xarra,icoun,ipass)
    character(*), intent(inout) :: xx
    character(*), intent(inout) :: xarra
    integer(ip),  intent(out)   :: icoun
    integer(ip),  intent(in)    :: ipass
    integer(ip)                 :: iposi,len_xarra,len_xx
    len_xx    = len(xx)
    len_xarra = len(xarra)
    iposi     = icoun
    icoun     = icoun + len_xx
    if( ipass == 2 ) then
       if( iposi+len_xx > len_xarra ) call runend('CEXCHA: MAXIMUM PARCH LENGTH EXCEEDED')
       if( IMASTER ) xarra(iposi+1:iposi+len_xx) = xx(1:len_xx)
       if( ISLAVE  ) xx(1:len_xx) = xarra(iposi+1:iposi+len_xx)
    end if
  end subroutine PAR_EXCHANGE_CH

  subroutine PAR_EXCHANGE_CH_1(ncoun,xx,xarra,icoun,ipass)

    integer(ip),      intent(in)             :: ncoun
    character(len=*), intent(inout)          :: xx
    character(len=:), intent(inout), pointer :: xarra
    integer(ip),      intent(out)            :: icoun
    integer(ip),      intent(in)             :: ipass
    integer(ip)                              :: iposi,len_xarra,len_xx
    len_xx    = ncoun !len(xx)
    len_xarra = len(xarra)
    iposi     = icoun
    icoun     = icoun + len_xx
    if( ipass == 2 ) then
       if( iposi+len_xx > len_xarra ) call runend('CEXCHA: MAXIMUM PARCH LENGTH EXCEEDED')
       if( IMASTER ) xarra(iposi+1:iposi+len_xx) = xx(1:len_xx)
       if( ISLAVE  ) xx(1:len_xx) = xarra(iposi+1:iposi+len_xx)
    end if
    
  end subroutine PAR_EXCHANGE_CH_1

  subroutine PAR_EXCHANGE_LG_s(xx,xarra,icoun,ipass)
    logical(lg), intent(inout)          :: xx
    logical(lg), intent(inout), pointer :: xarra(:)
    integer(ip), intent(out)            :: icoun
    integer(ip), intent(in)             :: ipass
    icoun = icoun + 1
    if( ipass == 2 ) then
       if( IMASTER ) xarra(icoun) = xx
       if( ISLAVE  ) xx           = xarra(icoun)
    end if
  end subroutine PAR_EXCHANGE_LG_s
  subroutine PAR_EXCHANGE_LG_0(n,xx,xarra,icoun,ipass)
    integer(ip), intent(in)             :: n
    logical(lg), intent(inout)          :: xx(:)
    logical(lg), intent(inout), pointer :: xarra(:)
    integer(ip), intent(out)            :: icoun
    integer(ip), intent(in)             :: ipass
    integer(ip)                         :: ii
    do ii = 1,n
       icoun = icoun + 1
       if( ipass == 2 ) then
          if( IMASTER ) xarra(icoun) = xx(ii)
          if( ISLAVE  ) xx(ii)       = xarra(icoun)
       end if
    end do
  end subroutine PAR_EXCHANGE_LG_0
  subroutine PAR_EXCHANGE_LG_1(xx,xarra,icoun,ipass)
    logical(lg), intent(inout), pointer :: xx(:)
    logical(lg), intent(inout), pointer :: xarra(:)
    integer(ip), intent(out)            :: icoun
    integer(ip), intent(in)             :: ipass
    integer(ip)                         :: ii
    do ii = 1,memory_size(xx)
       icoun = icoun + 1
       if( ipass == 2 ) then
          if( IMASTER ) xarra(icoun) = xx(ii)
          if( ISLAVE  ) xx(ii)       = xarra(icoun)
       end if
    end do
  end subroutine PAR_EXCHANGE_LG_1

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    29/04/2014
  !> @brief   Define if I am the master in this comunicatior
  !> @details Answer is TRUE if my rank is zero in communicator
  !
  !----------------------------------------------------------------------

  function PAR_IMASTER_IN_COMMUNICATOR(PAR_COMM_TO_USE)
    integer(4),  intent(in) :: PAR_COMM_TO_USE              !< Communicator
    integer(4)              :: my_rank
    logical(lg)             :: PAR_IMASTER_IN_COMMUNICATOR

    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank)
    if( my_rank == 0_4 ) then
       PAR_IMASTER_IN_COMMUNICATOR = .true.
    else
       PAR_IMASTER_IN_COMMUNICATOR = .false.
    end if

  end function PAR_IMASTER_IN_COMMUNICATOR

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    16/05/2014
  !> @brief   MPI Barrier
  !> @details MPI Barrier
  !
  !----------------------------------------------------------------------

  subroutine PAR_BARRIER(wherein,PAR_COMM_IN4)
    character(*), optional, intent(in) :: wherein
    integer(4),   optional, intent(in) :: PAR_COMM_IN4
    integer(4)                         :: PAR_COMM_TO_USE4
    integer(4)                         :: istat4

#ifndef MPI_OFF
    if( IPARALL ) then
       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE4=PAR_COMM_IN4
       else
          if( present(wherein) ) then
             call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE4)
          else
             call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE4)
          end if
       endif
       call MPI_Barrier( PAR_COMM_TO_USE4, istat4 )
    end if
#endif

  end subroutine PAR_BARRIER

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    16/05/2014
  !> @brief   Initialize non-blocking communications
  !> @details Initialize non-blocking communications
  !
  !----------------------------------------------------------------------

  subroutine PAR_INITIALIZE_NON_BLOCKING_COMM()
    integer(ip) :: ii
    nullify(non_blocking)
    allocate(non_blocking(10))
    do ii = 1,size(non_blocking)
       non_blocking(ii) % count4 = 0
       nullify(non_blocking(ii) % request4)
    end do
  end subroutine PAR_INITIALIZE_NON_BLOCKING_COMM

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    16/05/2014
  !> @brief   Start non-blocking communications
  !> @details Start non-blocking communications
  !>
  !>          KNONBLOCKING ............. communication number
  !>          MAX_NUMBER_NEIGHBORS*2 ... maximum number of requests
  !
  !----------------------------------------------------------------------

  subroutine PAR_START_NON_BLOCKING_COMM_8(knonblocking,max_number_neighbors)
    integer(ip), intent(in)           :: knonblocking
    integer(8),  intent(in)           :: max_number_neighbors

    if( knonblocking <= 0 .or. knonblocking > size(non_blocking) ) then
       call runend('PAR_START_NON_BLOCKING_COMM: WRONG NON-BLOCKING STRUCTURE NUMBER')
    else
       if( non_blocking(knonblocking) % count4 /= 0 ) then
          call runend('PAR_START_NON_BLOCKING_COMM: AN ASYNCHRONOUS SEND/RECEIVE IS ALREADY STARTED')
       else
          if (associated(non_blocking(knonblocking) % request4) ) then
             call runend('PAR_START_NON_BLOCKING_COMM: A NON-BLOCKING COMMUNICATOR NUMBER IS ALREADY ASSOCIATED')
          else
             if( max_number_neighbors > 0 ) &
                  allocate(non_blocking(knonblocking) % request4(max_number_neighbors*2))
             non_blocking(knonblocking) % count4 = 0
          end if
       end if
    end if
  end subroutine PAR_START_NON_BLOCKING_COMM_8

  subroutine PAR_START_NON_BLOCKING_COMM_4(knonblocking,max_number_neighbors)
    integer(ip), intent(in)           :: knonblocking
    integer(4),  intent(in)           :: max_number_neighbors

    if( knonblocking <= 0 .or. knonblocking > size(non_blocking) ) then
       write(*,*)'knonblocking,size(non_blocking)',knonblocking,size(non_blocking)
       call runend('PAR_START_NON_BLOCKING_COMM: WRONG NON-BLOCKING STRUCTURE NUMBER')
    else
       if( non_blocking(knonblocking) % count4 /= 0 ) then
          call runend('PAR_START_NON_BLOCKING_COMM: AN ASYNCHRONOUS SEND/RECEIVE IS ALREADY STARTED')
       else
          if (associated(non_blocking(knonblocking) % request4) ) then
             call runend('PAR_START_NON_BLOCKING_COMM: A NON-BLOCKING COMMUNICATOR NUMBER IS ALREADY ASSOCIATED')
          else
             allocate(non_blocking(knonblocking) % request4(max_number_neighbors*2))
             non_blocking(knonblocking) % count4 = 0
          end if
       end if
    end if
  end subroutine PAR_START_NON_BLOCKING_COMM_4

  subroutine PAR_START_NON_BLOCKING_COMM_4_OPT(knonblocking)
    integer(ip), intent(in)           :: knonblocking
    integer(ip)                       :: max_number_neighbors
    integer(ip)                       :: PAR_CURRENT_RANK,PAR_CURRENT_SIZE

    call PAR_COMM_RANK_AND_SIZE(commd % PAR_COMM_WORLD,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)
    max_number_neighbors = PAR_CURRENT_SIZE

    if( knonblocking <= 0 .or. knonblocking > size(non_blocking) ) then
       write(*,*)'knonblocking,size(non_blocking)',knonblocking,size(non_blocking)
       call runend('PAR_START_NON_BLOCKING_COMM: WRONG NON-BLOCKING STRUCTURE NUMBER')
    else
       if( non_blocking(knonblocking) % count4 /= 0 ) then
          call runend('PAR_START_NON_BLOCKING_COMM: AN ASYNCHRONOUS SEND/RECEIVE IS ALREADY STARTED')
       else
          if (associated(non_blocking(knonblocking) % request4) ) then
             call runend('PAR_START_NON_BLOCKING_COMM: A NON-BLOCKING COMMUNICATOR NUMBER IS ALREADY ASSOCIATED')
          else
             allocate(non_blocking(knonblocking) % request4(max_number_neighbors*2))
             non_blocking(knonblocking) % count4 = 0
          end if
       end if
    end if
  end subroutine PAR_START_NON_BLOCKING_COMM_4_OPT

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    16/05/2014
  !> @brief   End non-blocking communications
  !> @details End non-blocking communications waiting all messages have
  !>          arrived
  !>
  !>          KNONBLOCKING ............. communication number
  !
  !----------------------------------------------------------------------

  subroutine PAR_END_NON_BLOCKING_COMM(knonblocking)
    integer(ip), intent(in) :: knonblocking
    integer(4),  pointer    :: status4(:,:)
    integer(4)              :: istat4

#ifndef MPI_OFF
    if( non_blocking(knonblocking) % count4 > 0 ) then
       allocate( status4(MPI_STATUS_SIZE,non_blocking(knonblocking) % count4) )
       CALL MPI_WAITALL(non_blocking(knonblocking) % count4,non_blocking(knonblocking) % request4,status4,istat4)
       if( istat4 /= 0 ) call runend('NON BLOCKING SEND/RECEIVE COULD NOT BE COMPLETED')
       deallocate( status4 )
       deallocate( non_blocking(knonblocking) % request4 )
       non_blocking(knonblocking) % count4 = 0
    end if
    if( associated(non_blocking(knonblocking) % request4) ) then
       deallocate( non_blocking(knonblocking) % request4  )
    end if
#endif

  end subroutine PAR_END_NON_BLOCKING_COMM

  !----------------------------------------------------------------------
  !
  !> @author  Guillaume Houzeaux
  !> @date    16/05/2014
  !> @brief   Set non-blocking communication number
  !> @details Set non-blocking communication number
  !>
  !>          KNONBLOCKING ............. current communication number
  !
  !----------------------------------------------------------------------

  subroutine PAR_SET_NON_BLOCKING_COMM_NUMBER(knonblocking)
    integer(ip), intent(in) :: knonblocking

    if( knonblocking < 0 .or. knonblocking > size(non_blocking) ) then
       call runend('PAR_START_NON_BLOCKING_COMM: WRONG NON-BLOCKING STRUCTURE NUMBER')
    else
       inonblocking = knonblocking
    end if

  end subroutine PAR_SET_NON_BLOCKING_COMM_NUMBER

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_GHOST_NODE_EXCHANGE_RP
  !
  !----------------------------------------------------------------------

  subroutine PAR_GHOST_NODE_EXCHANGE_RP_00(n,xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_GHOST_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_GHOST_NODE_EXCHANGE_RP_00

  subroutine PAR_GHOST_NODE_EXCHANGE_RP_0(n,xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_GHOST_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_GHOST_NODE_EXCHANGE_RP_0

  subroutine PAR_GHOST_NODE_EXCHANGE_RP_1(xx,what,wherein,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = 1
    if( associated(xx) ) then
       if( size(xx,1) /= npoin_2 ) call runend('PAR_GHOST_NODE_EXCHANGE_RP_1: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_NODE_EXCHANGE_RP_1

  subroutine PAR_GHOST_NODE_EXCHANGE_RP_2(xx,what,wherein,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = size(xx,1)
    if( associated(xx) ) then
       if( size(xx,2) /= npoin_2 ) call runend('PAR_GHOST_NODE_EXCHANGE_RP_2: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_NODE_EXCHANGE_RP_2

  subroutine PAR_GHOST_NODE_EXCHANGE_RP_3(xx,what,wherein,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = size(xx,1)*size(xx,2)
       if( size(xx,3) /= npoin_2 ) call runend('PAR_GHOST_NODE_EXCHANGE_RP_3: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_NODE_EXCHANGE_RP_3

  subroutine PAR_GHOST_NODE_EXCHANGE_RP_2b(xx,what,commu,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= npoin_2 ) call runend('PAR_GHOST_NODE_EXCHANGE_RP_2b: WRONG SIZE')
       PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
       call PAR_GHOST_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_NODE_EXCHANGE_RP_2b

  !----------------------------------------------------------------------
  !
  ! PAR_GHOST_NODE_EXCHANGE_RP: NODE ASSEMBLY FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_GHOST_NODE_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: ndofn
    real(rp),                       intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    integer(4),                     intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,jj,dom_i
    integer(ip)                                   :: nsize_send,nsize_recv
    integer(ip)                                   :: ipoin,ini_send,ini_recv,kk,idofn
    integer(4)                                    :: nsize_send4,nsize_recv4
    integer(4)                                    :: istat4,count4
    integer(4)                                    :: dom_i4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    integer(4),                        pointer    :: status4(:,:)
    integer(ip)                                   :: dom_j

#ifndef MPI_OFF
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( commu % ghost_send_node_dim == -1 .or. commu % ghost_recv_node_dim == -1 ) &
                  call runend('FRINGE NODE EXCHANGE NOT COMPUTED')
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_rsend(commu % ghost_send_node_dim * ndofn))
             allocate(tmp_rrecv(commu % ghost_recv_node_dim * ndofn))
             istat4 = 0_4
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,commu % ghost_send_node_dim
                ipoin = commu % ghost_send_node_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_rsend(kk) = xx(idofn,ipoin)
                end do
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             kk = 0
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   dom_i      = commd % neights(ii)
                   ini_send   = ndofn * ( commd % ghost_send_node_size(ii)   -1 ) + 1
                   nsize_send = ndofn * ( commd % ghost_send_node_size(ii+1) -1 ) + 1 - ini_send
                   ini_recv   = ndofn * ( commd % ghost_recv_node_size(ii)   -1 ) + 1
                   nsize_recv = ndofn * ( commd % ghost_recv_node_size(ii+1) -1 ) + 1 - ini_recv

                   nsize_send4 = int(nsize_send,4)
                   nsize_recv4 = int(nsize_recv,4)

                   if( asynch ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           tmp_rsend(ini_send:ini_send+nsize_send-1), nsize_send4, &
                           MPI_DOUBLE_PRECISION,  dom_i4, 0_4,                     &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                      kk = kk + 1
                      call MPI_Irecv(&
                           tmp_rrecv(ini_recv:ini_recv+nsize_recv-1), nsize_recv4, &
                           MPI_DOUBLE_PRECISION,  dom_i4, 0_4,                     &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                   else
                      if( nsize_recv /= 0 .and. nsize_send == 0 ) then
                         call MPI_Recv(                          &
                              tmp_rrecv(ini_recv:), nsize_recv4, &
                              MPI_DOUBLE_PRECISION, dom_i4, 0_4, &
                              PAR_COMM_TO_USE, status, istat4    )
                      else if( nsize_recv == 0 .and. nsize_send /= 0 ) then
                         call MPI_Send(                          &
                              tmp_rsend(ini_send:), nsize_send4, &
                              MPI_DOUBLE_PRECISION, dom_i4, 0_4, &
                              PAR_COMM_TO_USE, istat4            )
                      else if( nsize_recv /= 0 .and. nsize_send /= 0 ) then
                         call MPI_Sendrecv(                         &
                              tmp_rsend(ini_send:), nsize_send4,    &
                              MPI_DOUBLE_PRECISION, dom_i4, 0_4,    &
                              tmp_rrecv(ini_recv:), nsize_recv4,    &
                              MPI_DOUBLE_PRECISION, dom_i4, 0_4,    &
                              PAR_COMM_TO_USE, status, istat4       )
                      end if
                   end if
                   if( istat4 /= 0_4 ) call runend('PAR_GHOST_NODE_EXCHANGE_RP: MPI ERROR')

                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE,2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             if( istat4 /= 0 ) call runend('WRONG SEND/RECEIVE')
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if( trim(what) == 'SUM' .or. trim(what) == 'ASSEMBLY' ) then
                !
                ! SUM
                !
                kk = 0
                do jj = 1,commu % ghost_recv_node_dim
                   ipoin = commu % ghost_recv_node_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = xx(idofn,ipoin) + tmp_rrecv(kk)
                   end do
                end do

             else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' ) then
                !
                ! MAX
                !
                kk = 0
                do jj = 1,commu % ghost_recv_node_dim
                   ipoin = commu % ghost_recv_node_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = max(xx(idofn,ipoin),tmp_rrecv(kk))
                   end do
                end do

             else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' ) then
                !
                ! MIN
                !
                kk = 0
                do jj = 1,commu % ghost_recv_node_dim
                   ipoin = commu % ghost_recv_node_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = min(xx(idofn,ipoin),tmp_rrecv(kk))
                   end do
                end do

             else if( trim(what) == 'REPLACE' .or. trim(what) == 'SUBSTITUTE' ) then
                !
                ! Replace value on my fringe nodes according to what I have received
                !
                kk = 0
                do jj = 1,commu % ghost_recv_node_dim
                   ipoin = commu % ghost_recv_node_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = tmp_rrecv(kk)
                   end do
                end do

             else
                call runend('UNKNOWN ORDER')
             end if

             ipass = 0
             deallocate(tmp_rrecv)
             deallocate(tmp_rsend)

          end if

       end if

    end if
#endif

  end subroutine PAR_GHOST_NODE_EXCHANGE_RP

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_GHOST_NODE_EXCHANGE_IP
  !
  !----------------------------------------------------------------------

  subroutine PAR_GHOST_NODE_EXCHANGE_IP_0(n,xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: n
    integer(ip),                    intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_GHOST_NODE_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_GHOST_NODE_EXCHANGE_IP_0

  subroutine PAR_GHOST_NODE_EXCHANGE_IP_1(xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),          pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    integer(4)                                    :: PAR_COMM_TO_USE

   if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = 1
       if( size(xx,1) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_NODE_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_NODE_EXCHANGE_IP_1

  subroutine PAR_GHOST_NODE_EXCHANGE_IP_2(xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_NODE_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_NODE_EXCHANGE_IP_2

  subroutine PAR_GHOST_NODE_EXCHANGE_IP_3(xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),          pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = size(xx,1)*size(xx,2)
       if( size(xx,3) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_NODE_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_NODE_EXCHANGE_IP_3

  subroutine PAR_GHOST_NODE_EXCHANGE_IP_2b(xx,what,commu,wsynch,dom_k)
    implicit none
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),  pointer,  intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= npoin .and. size(xx,2) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
       call PAR_GHOST_NODE_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_NODE_EXCHANGE_IP_2b

 !----------------------------------------------------------------------
  !
  ! PAR_GHOST_NODE_EXCHANGE_RP: NODE ASSEMBLY FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_GHOST_NODE_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: ndofn
    integer(ip),                    intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    integer(4),                     intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,jj,dom_i
    integer(ip)                                   :: nsize_send,nsize_recv
    integer(ip)                                   :: ipoin,ini_send,ini_recv,kk,idofn
    integer(4)                                    :: nsize_send4,nsize_recv4
    integer(4)                                    :: istat4,count4
    integer(4)                                    :: dom_i4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    integer(4),                        pointer    :: status4(:,:)
    integer(ip)                                   :: dom_j

#ifndef MPI_OFF
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( commu % ghost_send_node_dim == -1 .or. commu % ghost_recv_node_dim == -1 ) &
                  call runend('FRINGE NODE EXCHANGE NOT COMPUTED')
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_isend(commu % ghost_send_node_dim * ndofn))
             allocate(tmp_irecv(commu % ghost_recv_node_dim * ndofn))
             istat4 = 0_4
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,commu % ghost_send_node_dim
                ipoin = commu % ghost_send_node_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_isend(kk) = xx(idofn,ipoin)
                end do
             end do
             istat4 = 0_4
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             kk = 0
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   dom_i      = commd % neights(ii)
                   ini_send   = ndofn * ( commd % ghost_send_node_size(ii)   -1 ) + 1
                   nsize_send = ndofn * ( commd % ghost_send_node_size(ii+1) -1 ) + 1 - ini_send
                   ini_recv   = ndofn * ( commd % ghost_recv_node_size(ii)   -1 ) + 1
                   nsize_recv = ndofn * ( commd % ghost_recv_node_size(ii+1) -1 ) + 1 - ini_recv

                   nsize_send4 = int(nsize_send,4)
                   nsize_recv4 = int(nsize_recv,4)

                   if( asynch ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           tmp_isend(ini_send:ini_send+nsize_send-1), nsize_send4, &
                           PAR_INTEGER,  dom_i4, 0_4,                              &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                      kk = kk + 1
                      call MPI_Irecv(&
                           tmp_irecv(ini_recv:ini_recv+nsize_recv-1), nsize_recv4, &
                           PAR_INTEGER,  dom_i4, 0_4,                              &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                   else
                      if( nsize_recv /= 0 .and. nsize_send == 0 ) then
                         call MPI_Recv(                          &
                              tmp_irecv(ini_recv:), nsize_recv4, &
                              PAR_INTEGER, dom_i4, 0_4,          &
                              PAR_COMM_TO_USE, status, istat4    )
                      else if( nsize_recv == 0 .and. nsize_send /= 0 ) then
                         call MPI_Send(                          &
                              tmp_isend(ini_send:), nsize_send4, &
                              PAR_INTEGER, dom_i4, 0_4,          &
                              PAR_COMM_TO_USE, istat4            )
                      else if( nsize_recv /= 0 .and. nsize_send /= 0 ) then
                         call MPI_Sendrecv(                      &
                              tmp_isend(ini_send:), nsize_send4, &
                              PAR_INTEGER, dom_i4, 0_4,          &
                              tmp_irecv(ini_recv:), nsize_recv4, &
                              PAR_INTEGER, dom_i4, 0_4,          &
                              PAR_COMM_TO_USE, status, istat4    )
                      end if
                   end if
                   if( istat4 /= 0_4 ) call runend('PAR_GHOST_NODE_EXCHANGE_RP: MPI ERROR')

                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE,2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             if( istat4 /= 0 ) call runend('WRONG SEND/RECEIVE')
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if( trim(what) == 'SUM' .or. trim(what) == 'ASSEMBLY' ) then
                !
                ! SUM
                !
                kk = 0
                do jj = 1,commu % ghost_recv_node_dim
                   ipoin = commu % ghost_recv_node_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = xx(idofn,ipoin) + tmp_irecv(kk)
                   end do
                end do

             else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' ) then
                !
                ! MAX
                !
                kk = 0
                do jj = 1,commu % ghost_recv_node_dim
                   ipoin = commu % ghost_recv_node_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = max(xx(idofn,ipoin),tmp_irecv(kk))
                   end do
                end do

             else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' ) then
                !
                ! MIN
                !
                kk = 0
                do jj = 1,commu % ghost_recv_node_dim
                   ipoin = commu % ghost_recv_node_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = min(xx(idofn,ipoin),tmp_irecv(kk))
                   end do
                end do

             else if( trim(what) == 'REPLACE' .or. trim(what) == 'SUBSTITUTE' ) then
                !
                ! Replace value on my fringe nodes according to what I have received
                !
                kk = 0
                do jj = 1,commu % ghost_recv_node_dim
                   ipoin = commu % ghost_recv_node_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = tmp_irecv(kk)
                   end do
                end do

             else
                call runend('UNKNOWN ORDER')
             end if

             ipass = 0
             deallocate(tmp_irecv)
             deallocate(tmp_isend)

          end if

       end if

    end if
#endif

  end subroutine PAR_GHOST_NODE_EXCHANGE_IP

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_GHOST_ELEMENT_EXCHANGE_IP
  !
  !----------------------------------------------------------------------

  subroutine PAR_GHOST_ELEMENT_EXCHANGE_IP_0(n,xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: n
    integer(ip),                    intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_GHOST_ELEMENT_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_GHOST_ELEMENT_EXCHANGE_IP_0

  subroutine PAR_GHOST_ELEMENT_EXCHANGE_IP_1(xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),          pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = 1
       if( size(xx,1) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_ELEMENT_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_ELEMENT_EXCHANGE_IP_1

  subroutine PAR_GHOST_ELEMENT_EXCHANGE_IP_2(xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_ELEMENT_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_ELEMENT_EXCHANGE_IP_2

  subroutine PAR_GHOST_ELEMENT_EXCHANGE_IP_3(xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),          pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = size(xx,1)*size(xx,2)
       if( size(xx,3) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_ELEMENT_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_ELEMENT_EXCHANGE_IP_3

  subroutine PAR_GHOST_ELEMENT_EXCHANGE_IP_2b(xx,what,commu,wsynch,dom_k)
    implicit none
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),  pointer,  intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
       call PAR_GHOST_ELEMENT_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_ELEMENT_EXCHANGE_IP_2b

 !----------------------------------------------------------------------
  !
  ! PAR_GHOST_ELEMENT_EXCHANGE_RP: NODE ASSEMBLY FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_GHOST_ELEMENT_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: ndofn
    integer(ip),                    intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    integer(4),                     intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,jj,dom_i
    integer(ip)                                   :: nsize_send,nsize_recv
    integer(ip)                                   :: ipoin,ini_send,ini_recv,kk,idofn
    integer(4)                                    :: nsize_send4,nsize_recv4
    integer(4)                                    :: istat4,count4
    integer(4)                                    :: dom_i4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    integer(4),                        pointer    :: status4(:,:)
    integer(ip)                                   :: dom_j

#ifndef MPI_OFF
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_ELEM_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( commu % ghost_send_elem_dim == -1 .or. commu % ghost_recv_elem_dim == -1 ) &
                  call runend('FRINGE NODE EXCHANGE NOT COMPUTED')
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_isend(commu % ghost_send_elem_dim * ndofn))
             allocate(tmp_irecv(commu % ghost_recv_elem_dim * ndofn))
             istat4 = 0_4
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,commu % ghost_send_elem_dim
                ipoin = commu % ghost_send_elem_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_isend(kk) = xx(idofn,ipoin)
                end do
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             kk = 0
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   dom_i      = commd % neights(ii)
                   ini_send   = ndofn * ( commd % ghost_send_elem_size(ii)   -1 ) + 1
                   nsize_send = ndofn * ( commd % ghost_send_elem_size(ii+1) -1 ) + 1 - ini_send
                   ini_recv   = ndofn * ( commd % ghost_recv_elem_size(ii)   -1 ) + 1
                   nsize_recv = ndofn * ( commd % ghost_recv_elem_size(ii+1) -1 ) + 1 - ini_recv

                   nsize_send4 = int(nsize_send,4)
                   nsize_recv4 = int(nsize_recv,4)

                   if( asynch ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           tmp_isend(ini_send:ini_send+nsize_send-1), nsize_send4, &
                           PAR_INTEGER,  dom_i4, 0_4,                              &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                      kk = kk + 1
                      call MPI_Irecv(&
                           tmp_irecv(ini_recv:ini_recv+nsize_recv-1), nsize_recv4, &
                           PAR_INTEGER,  dom_i4, 0_4,                              &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                   else
                      if( nsize_recv /= 0 .and. nsize_send == 0 ) then
                         call MPI_Recv(                          &
                              tmp_irecv(ini_recv:), nsize_recv4, &
                              PAR_INTEGER, dom_i4, 0_4,          &
                              PAR_COMM_TO_USE, status, istat4    )
                      else if( nsize_recv == 0 .and. nsize_send /= 0 ) then
                         call MPI_Send(                          &
                              tmp_isend(ini_send:), nsize_send4, &
                              PAR_INTEGER, dom_i4, 0_4,          &
                              PAR_COMM_TO_USE, istat4            )
                      else if( nsize_recv /= 0 .and. nsize_send /= 0 ) then
                         call MPI_Sendrecv(                      &
                              tmp_isend(ini_send:), nsize_send4, &
                              PAR_INTEGER, dom_i4, 0_4,          &
                              tmp_irecv(ini_recv:), nsize_recv4, &
                              PAR_INTEGER, dom_i4, 0_4,          &
                              PAR_COMM_TO_USE, status, istat4    )
                      end if
                   end if
                   if( istat4 /= 0_4 ) call runend('PAR_GHOST_ELEMENT_EXCHANGE_RP: MPI ERROR')

                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE,2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             if( istat4 /= 0 ) call runend('WRONG SEND/RECEIVE')
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if( trim(what) == 'SUM' .or. trim(what) == 'ASSEMBLY' ) then
                !
                ! SUM
                !
                kk = 0
                do jj = 1,commu % ghost_recv_elem_dim
                   ipoin = commu % ghost_recv_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = xx(idofn,ipoin) + tmp_irecv(kk)
                   end do
                end do

             else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' ) then
                !
                ! MAX
                !
                kk = 0
                do jj = 1,commu % ghost_recv_elem_dim
                   ipoin = commu % ghost_recv_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = max(xx(idofn,ipoin),tmp_irecv(kk))
                   end do
                end do

             else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' ) then
                !
                ! MIN
                !
                kk = 0
                do jj = 1,commu % ghost_recv_elem_dim
                   ipoin = commu % ghost_recv_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = min(xx(idofn,ipoin),tmp_irecv(kk))
                   end do
                end do

             else if( trim(what) == 'REPLACE' .or. trim(what) == 'SUBSTITUTE' ) then
                !
                ! Replace value on my fringe nodes according to what I have received
                !
                kk = 0
                do jj = 1,commu % ghost_recv_elem_dim
                   ipoin = commu % ghost_recv_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = tmp_irecv(kk)
                   end do
                end do

             else
                call runend('UNKNOWN ORDER')
             end if

             ipass = 0
             deallocate(tmp_irecv)
             deallocate(tmp_isend)

          end if

       end if

    end if
#endif

  end subroutine PAR_GHOST_ELEMENT_EXCHANGE_IP

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_GHOST_ELEMENT_EXCHANGE_RP
  !
  !----------------------------------------------------------------------

  subroutine PAR_GHOST_ELEMENT_EXCHANGE_RP_00(n,xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_GHOST_ELEMENT_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_GHOST_ELEMENT_EXCHANGE_RP_00

  subroutine PAR_GHOST_ELEMENT_EXCHANGE_RP_0(n,xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_GHOST_ELEMENT_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_GHOST_ELEMENT_EXCHANGE_RP_0

  subroutine PAR_GHOST_ELEMENT_EXCHANGE_RP_1(xx,what,wherein,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = 1
       if( size(xx,1) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_ELEMENT_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_ELEMENT_EXCHANGE_RP_1

  subroutine PAR_GHOST_ELEMENT_EXCHANGE_RP_2(xx,what,wherein,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = size(xx,1)
    if( size(xx,2) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_GHOST_ELEMENT_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_GHOST_ELEMENT_EXCHANGE_RP_2

  subroutine PAR_GHOST_ELEMENT_EXCHANGE_RP_3(xx,what,wherein,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = size(xx,1)*size(xx,2)
       if( size(xx,3) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_ELEMENT_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_ELEMENT_EXCHANGE_RP_3

  subroutine PAR_GHOST_ELEMENT_EXCHANGE_RP_2b(xx,what,commu,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
       call PAR_GHOST_ELEMENT_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_ELEMENT_EXCHANGE_RP_2b

  !----------------------------------------------------------------------
  !
  ! PAR_GHOST_ELEMENT_EXCHANGE_RP: NODE ASSEMBLY FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_GHOST_ELEMENT_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: ndofn
    real(rp),                       intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    integer(4),                     intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,jj,dom_i
    integer(ip)                                   :: nsize_send,nsize_recv
    integer(ip)                                   :: ipoin,ini_send,ini_recv,kk,idofn
    integer(4)                                    :: nsize_send4,nsize_recv4
    integer(4)                                    :: istat4,count4
    integer(4)                                    :: dom_i4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    integer(4),                        pointer    :: status4(:,:)
    integer(ip)                                   :: dom_j

#ifndef MPI_OFF
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_ELEMENT_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( commu % ghost_send_elem_dim == -1 .or. commu % ghost_recv_elem_dim == -1 ) &
                  call runend('FRINGE NODE EXCHANGE NOT COMPUTED')
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_rsend(commu % ghost_send_elem_dim * ndofn))
             allocate(tmp_rrecv(commu % ghost_recv_elem_dim * ndofn))
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,commu % ghost_send_elem_dim
                ipoin = commu % ghost_send_elem_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_rsend(kk) = xx(idofn,ipoin)
                end do
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             kk = 0
             istat4 = 0_4
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   dom_i      = commd % neights(ii)
                   ini_send   = ndofn * ( commd % ghost_send_elem_size(ii)   -1 ) + 1
                   nsize_send = ndofn * ( commd % ghost_send_elem_size(ii+1) -1 ) + 1 - ini_send
                   ini_recv   = ndofn * ( commd % ghost_recv_elem_size(ii)   -1 ) + 1
                   nsize_recv = ndofn * ( commd % ghost_recv_elem_size(ii+1) -1 ) + 1 - ini_recv

                   nsize_send4 = int(nsize_send,4)
                   nsize_recv4 = int(nsize_recv,4)

                   if( asynch ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           tmp_rsend(ini_send:ini_send+nsize_send-1), nsize_send4, &
                           MPI_DOUBLE_PRECISION,  dom_i4, 0_4,                     &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                      kk = kk + 1
                      call MPI_Irecv(&
                           tmp_rrecv(ini_recv:ini_recv+nsize_recv-1), nsize_recv4, &
                           MPI_DOUBLE_PRECISION,  dom_i4, 0_4,                     &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                   else
                      if( nsize_recv /= 0 .and. nsize_send == 0 ) then
                         call MPI_Recv(                          &
                              tmp_rrecv(ini_recv:), nsize_recv4, &
                              MPI_DOUBLE_PRECISION, dom_i4, 0_4, &
                              PAR_COMM_TO_USE, status, istat4    )
                      else if( nsize_recv == 0 .and. nsize_send /= 0 ) then
                         call MPI_Send(                          &
                              tmp_rsend(ini_send:), nsize_send4, &
                              MPI_DOUBLE_PRECISION, dom_i4, 0_4, &
                              PAR_COMM_TO_USE, istat4            )
                      else if( nsize_recv /= 0 .and. nsize_send /= 0 ) then
                         call MPI_Sendrecv(                         &
                              tmp_rsend(ini_send:), nsize_send4,    &
                              MPI_DOUBLE_PRECISION, dom_i4, 0_4,    &
                              tmp_rrecv(ini_recv:), nsize_recv4,    &
                              MPI_DOUBLE_PRECISION, dom_i4, 0_4,    &
                              PAR_COMM_TO_USE, status, istat4       )
                      end if
                   end if
                   if( istat4 /= 0_4 ) call runend('PAR_GHOST_ELEMENT_EXCHANGE_RP: MPI ERROR')

                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE,2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             if( istat4 /= 0 ) call runend('WRONG SEND/RECEIVE')
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if( trim(what) == 'SUM' .or. trim(what) == 'ASSEMBLY' ) then
                !
                ! SUM
                !
                kk = 0
                do jj = 1,commu % ghost_recv_elem_dim
                   ipoin = commu % ghost_recv_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = xx(idofn,ipoin) + tmp_rrecv(kk)
                   end do
                end do

             else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' ) then
                !
                ! MAX
                !
                kk = 0
                do jj = 1,commu % ghost_recv_elem_dim
                   ipoin = commu % ghost_recv_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = max(xx(idofn,ipoin),tmp_rrecv(kk))
                   end do
                end do

             else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' ) then
                !
                ! MIN
                !
                kk = 0
                do jj = 1,commu % ghost_recv_elem_dim
                   ipoin = commu % ghost_recv_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = min(xx(idofn,ipoin),tmp_rrecv(kk))
                   end do
                end do

             else if( trim(what) == 'REPLACE' .or. trim(what) == 'SUBSTITUTE' ) then
                !
                ! Replace value on my fringe nodes according to what I have received
                !
                kk = 0
                do jj = 1,commu % ghost_recv_elem_dim
                   ipoin = commu % ghost_recv_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = tmp_rrecv(kk)
                   end do
                end do

             else
                call runend('UNKNOWN ORDER')
             end if

             ipass = 0
             deallocate(tmp_rrecv)
             deallocate(tmp_rsend)

          end if

       end if

    end if
#endif

  end subroutine PAR_GHOST_ELEMENT_EXCHANGE_RP

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_GHOST_BOUNDARY_EXCHANGE_IP
  !
  !----------------------------------------------------------------------

  subroutine PAR_GHOST_BOUNDARY_EXCHANGE_IP_0(n,xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: n
    integer(ip),                    intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_GHOST_BOUNDARY_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_GHOST_BOUNDARY_EXCHANGE_IP_0

  subroutine PAR_GHOST_BOUNDARY_EXCHANGE_IP_1(xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),          pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    ndofn = 1
    if( associated(xx) ) then
       !
       ! Boundary arrays are generally allocated to a minimum of 1
       ! for robustness reasons, therefore we have not SIZE(XX) = NBOUN_2
       ! when NBOUN_2 = 0
       !
       if( size(xx,1) /= nboun_2 .and. nboun_2 /= 0 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_BOUNDARY_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_BOUNDARY_EXCHANGE_IP_1

  subroutine PAR_GHOST_BOUNDARY_EXCHANGE_IP_2(xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    ndofn = size(xx,1)
    if( associated(xx) ) then
       !
       ! Boundary arrays are generally allocated to a minimum of 1
       ! for robustness reasons, therefore we have not SIZE(XX) = NBOUN_2
       ! when NBOUN_2 = 0
       !
       if( size(xx,2) /= nboun_2 .and. nboun_2 /= 0 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_BOUNDARY_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_BOUNDARY_EXCHANGE_IP_2

  subroutine PAR_GHOST_BOUNDARY_EXCHANGE_IP_3(xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),          pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    ndofn = size(xx,1)*size(xx,2)
    if( associated(xx) ) then
       !
       ! Boundary arrays are generally allocated to a minimum of 1
       ! for robustness reasons, therefore we have not SIZE(XX) = NBOUN_2
       ! when NBOUN_2 = 0
       !
       if( size(xx,3) /= nboun_2 .and. nboun_2 /= 0 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_BOUNDARY_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_BOUNDARY_EXCHANGE_IP_3

  subroutine PAR_GHOST_BOUNDARY_EXCHANGE_IP_2b(xx,what,commu,wsynch,dom_k)
    implicit none
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),  pointer,  intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    ndofn = size(xx,1)
    if( associated(xx) ) then
       !
       ! Boundary arrays are generally allocated to a minimum of 1
       ! for robustness reasons, therefore we have not SIZE(XX) = NBOUN_2
       ! when NBOUN_2 = 0
       !
       if( size(xx,2) /= nboun_2 .and. nboun_2 /= 0 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
       call PAR_GHOST_BOUNDARY_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_BOUNDARY_EXCHANGE_IP_2b

 !----------------------------------------------------------------------
  !
  ! PAR_GHOST_BOUNDARY_EXCHANGE_RP: NODE ASSEMBLY FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_GHOST_BOUNDARY_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: ndofn
    integer(ip),                    intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    integer(4),                     intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,jj,dom_i,dummi
    integer(ip)                                   :: nsize_send,nsize_recv
    integer(ip)                                   :: ipoin,ini_send,ini_recv,kk,idofn
    integer(4)                                    :: nsize_send4,nsize_recv4
    integer(4)                                    :: istat4,count4
    integer(4)                                    :: dom_i4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    integer(4),                        pointer    :: status4(:,:)
    integer(ip)                                   :: dom_j

#ifndef MPI_OFF
    istat4 = 0_4
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_BOUN_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( commu % ghost_send_boun_dim == -1 .or. commu % ghost_recv_boun_dim == -1 ) &
                  call runend('FRINGE NODE EXCHANGE NOT COMPUTED')
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_isend(commu % ghost_send_boun_dim * ndofn))
             allocate(tmp_irecv(commu % ghost_recv_boun_dim * ndofn))
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,commu % ghost_send_boun_dim
                ipoin = commu % ghost_send_boun_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_isend(kk) = xx(idofn,ipoin)
                end do
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             kk = 0
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   dom_i      = commd % neights(ii)
                   ini_send   = ndofn * ( commd % ghost_send_boun_size(ii)   -1 ) + 1
                   nsize_send = ndofn * ( commd % ghost_send_boun_size(ii+1) -1 ) + 1 - ini_send
                   ini_recv   = ndofn * ( commd % ghost_recv_boun_size(ii)   -1 ) + 1
                   nsize_recv = ndofn * ( commd % ghost_recv_boun_size(ii+1) -1 ) + 1 - ini_recv

                   nsize_send4 = int(nsize_send,4)
                   nsize_recv4 = int(nsize_recv,4)

                   if( asynch ) then
                      kk = kk + 1

                      if( nsize_send > 0 ) then
                         call MPI_Isend(&
                              dummi, nsize_send4,                                     &
                              PAR_INTEGER,  dom_i4, 0_4,                              &
                              PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                      else
                         call MPI_Isend(&
                              tmp_isend(ini_send:ini_send+nsize_send-1), nsize_send4, &
                              PAR_INTEGER,  dom_i4, 0_4,                              &
                              PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                      end if

                      kk = kk + 1

                      if( nsize_recv > 0 ) then
                         call MPI_Irecv(&
                              tmp_irecv(ini_recv:ini_recv+nsize_recv-1), nsize_recv4, &
                              PAR_INTEGER,  dom_i4, 0_4,                              &
                              PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                      else
                         call MPI_Irecv(&
                              dummi, nsize_recv4,                                     &
                              PAR_INTEGER,  dom_i4, 0_4,                              &
                              PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                      end if

                   else
                      if( nsize_recv > 0 .and. nsize_send <= 0 ) then
                         call MPI_Recv(                          &
                              tmp_irecv(ini_recv:), nsize_recv4, &
                              PAR_INTEGER, dom_i4, 0_4,          &
                              PAR_COMM_TO_USE, status, istat4    )
                      else if( nsize_recv <= 0 .and. nsize_send > 0 ) then
                         call MPI_Send(                          &
                              tmp_isend(ini_send:), nsize_send4, &
                              PAR_INTEGER, dom_i4, 0_4,          &
                              PAR_COMM_TO_USE, istat4            )
                      else if( nsize_recv > 0 .and. nsize_send > 0 ) then
                         call MPI_Sendrecv(                      &
                              tmp_isend(ini_send:), nsize_send4, &
                              PAR_INTEGER, dom_i4, 0_4,          &
                              tmp_irecv(ini_recv:), nsize_recv4, &
                              PAR_INTEGER, dom_i4, 0_4,          &
                              PAR_COMM_TO_USE, status, istat4    )
                      end if
                   end if
                   if( istat4 /= 0_4 ) call runend('PAR_GHOST_BOUNDARY_EXCHANGE_IP: MPI ERROR')

                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE,2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             if( istat4 /= 0 ) call runend('WRONG SEND/RECEIVE')
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if( trim(what) == 'SUM' .or. trim(what) == 'ASSEMBLY' ) then
                !
                ! SUM
                !
                kk = 0
                do jj = 1,commu % ghost_recv_boun_dim
                   ipoin = commu % ghost_recv_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = xx(idofn,ipoin) + tmp_irecv(kk)
                   end do
                end do

             else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' ) then
                !
                ! MAX
                !
                kk = 0
                do jj = 1,commu % ghost_recv_boun_dim
                   ipoin = commu % ghost_recv_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = max(xx(idofn,ipoin),tmp_irecv(kk))
                   end do
                end do

             else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' ) then
                !
                ! MIN
                !
                kk = 0
                do jj = 1,commu % ghost_recv_boun_dim
                   ipoin = commu % ghost_recv_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = min(xx(idofn,ipoin),tmp_irecv(kk))
                   end do
                end do

             else if( trim(what) == 'REPLACE' .or. trim(what) == 'SUBSTITUTE' ) then
                !
                ! Replace value on my fringe nodes according to what I have received
                !
                kk = 0
                do jj = 1,commu % ghost_recv_boun_dim
                   ipoin = commu % ghost_recv_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = tmp_irecv(kk)
                   end do
                end do

             else
                call runend('UNKNOWN ORDER')
             end if

             ipass = 0
             deallocate(tmp_irecv)
             deallocate(tmp_isend)

          end if

       end if

    end if
#endif

  end subroutine PAR_GHOST_BOUNDARY_EXCHANGE_IP

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_GHOST_BOUNDARY_EXCHANGE_RP
  !
  !----------------------------------------------------------------------

  subroutine PAR_GHOST_BOUNDARY_EXCHANGE_RP_00(n,xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_GHOST_BOUNDARY_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_GHOST_BOUNDARY_EXCHANGE_RP_00

  subroutine PAR_GHOST_BOUNDARY_EXCHANGE_RP_0(n,xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_GHOST_BOUNDARY_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_GHOST_BOUNDARY_EXCHANGE_RP_0

  subroutine PAR_GHOST_BOUNDARY_EXCHANGE_RP_1(xx,what,wherein,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = 1
    if( associated(xx) ) then
       if( size(xx,1) /= nboun_2 .and. nboun_2 /= 0 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_BOUNDARY_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_BOUNDARY_EXCHANGE_RP_1

  subroutine PAR_GHOST_BOUNDARY_EXCHANGE_RP_2(xx,what,wherein,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = size(xx,1)
    if( associated(xx) ) then
       if( size(xx,2) /= nboun_2 .and. nboun_2 /= 0 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_BOUNDARY_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_BOUNDARY_EXCHANGE_RP_2

  subroutine PAR_GHOST_BOUNDARY_EXCHANGE_RP_3(xx,what,wherein,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = size(xx,1)*size(xx,2)
    if( associated(xx) ) then
       if( size(xx,3) /= nboun_2 .and. nboun_2 /= 0 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_GHOST_BOUNDARY_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_BOUNDARY_EXCHANGE_RP_3

  subroutine PAR_GHOST_BOUNDARY_EXCHANGE_RP_2b(xx,what,commu,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    ndofn = size(xx,1)
    if( associated(xx) ) then
       if( size(xx,2) /= nboun_2 .and. nboun_2 /= 0 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
       call PAR_GHOST_BOUNDARY_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_GHOST_BOUNDARY_EXCHANGE_RP_2b

  !----------------------------------------------------------------------
  !
  ! PAR_GHOST_BOUNDARY_EXCHANGE_RP: NODE ASSEMBLY FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_GHOST_BOUNDARY_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: ndofn
    real(rp),                       intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    integer(4),                     intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,jj,dom_i
    integer(ip)                                   :: nsize_send,nsize_recv
    integer(ip)                                   :: ipoin,ini_send,ini_recv,kk,idofn
    integer(4)                                    :: nsize_send4,nsize_recv4
    integer(4)                                    :: istat4,count4
    integer(4)                                    :: dom_i4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    integer(4),                        pointer    :: status4(:,:)
    integer(ip)                                   :: dom_j

#ifndef MPI_OFF
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_BOUNDARY_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( commu % ghost_send_boun_dim == -1 .or. commu % ghost_recv_boun_dim == -1 ) &
                  call runend('FRINGE NODE EXCHANGE NOT COMPUTED')
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_rsend(commu % ghost_send_boun_dim * ndofn))
             allocate(tmp_rrecv(commu % ghost_recv_boun_dim * ndofn))
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,commu % ghost_send_boun_dim
                ipoin = commu % ghost_send_boun_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_rsend(kk) = xx(idofn,ipoin)
                end do
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             istat4 = 0_4
             kk = 0
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   dom_i      = commd % neights(ii)
                   ini_send   = ndofn * ( commd % ghost_send_boun_size(ii)   -1 ) + 1
                   nsize_send = ndofn * ( commd % ghost_send_boun_size(ii+1) -1 ) + 1 - ini_send
                   ini_recv   = ndofn * ( commd % ghost_recv_boun_size(ii)   -1 ) + 1
                   nsize_recv = ndofn * ( commd % ghost_recv_boun_size(ii+1) -1 ) + 1 - ini_recv

                   nsize_send4 = int(nsize_send,4)
                   nsize_recv4 = int(nsize_recv,4)

                   if( asynch ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           tmp_rsend(ini_send:ini_send+nsize_send-1), nsize_send4, &
                           MPI_DOUBLE_PRECISION,  dom_i4, 0_4,                     &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                      kk = kk + 1
                      call MPI_Irecv(&
                           tmp_rrecv(ini_recv:ini_recv+nsize_recv-1), nsize_recv4, &
                           MPI_DOUBLE_PRECISION,  dom_i4, 0_4,                     &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                   else
                      if( nsize_recv /= 0 .and. nsize_send == 0 ) then
                         call MPI_Recv(                          &
                              tmp_rrecv(ini_recv:), nsize_recv4, &
                              MPI_DOUBLE_PRECISION, dom_i4, 0_4, &
                              PAR_COMM_TO_USE, status, istat4    )
                      else if( nsize_recv == 0 .and. nsize_send /= 0 ) then
                         call MPI_Send(                          &
                              tmp_rsend(ini_send:), nsize_send4, &
                              MPI_DOUBLE_PRECISION, dom_i4, 0_4, &
                              PAR_COMM_TO_USE, istat4            )
                      else if( nsize_recv /= 0 .and. nsize_send /= 0 ) then
                         call MPI_Sendrecv(                         &
                              tmp_rsend(ini_send:), nsize_send4,    &
                              MPI_DOUBLE_PRECISION, dom_i4, 0_4,    &
                              tmp_rrecv(ini_recv:), nsize_recv4,    &
                              MPI_DOUBLE_PRECISION, dom_i4, 0_4,    &
                              PAR_COMM_TO_USE, status, istat4       )
                      end if
                   end if
                   if( istat4 /= 0_4 ) call runend('PAR_GHOST_BOUNDARY_EXCHANGE_RP: MPI ERROR')

                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE,2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             if( istat4 /= 0 ) call runend('WRONG SEND/RECEIVE')
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if( trim(what) == 'SUM' .or. trim(what) == 'ASSEMBLY' ) then
                !
                ! SUM
                !
                kk = 0
                do jj = 1,commu % ghost_recv_boun_dim
                   ipoin = commu % ghost_recv_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = xx(idofn,ipoin) + tmp_rrecv(kk)
                   end do
                end do

             else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' ) then
                !
                ! MAX
                !
                kk = 0
                do jj = 1,commu % ghost_recv_boun_dim
                   ipoin = commu % ghost_recv_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = max(xx(idofn,ipoin),tmp_rrecv(kk))
                   end do
                end do

             else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' ) then
                !
                ! MIN
                !
                kk = 0
                do jj = 1,commu % ghost_recv_boun_dim
                   ipoin = commu % ghost_recv_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = min(xx(idofn,ipoin),tmp_rrecv(kk))
                   end do
                end do

             else if( trim(what) == 'REPLACE' .or. trim(what) == 'SUBSTITUTE' ) then
                !
                ! Replace value on my fringe nodes according to what I have received
                !
                kk = 0
                do jj = 1,commu % ghost_recv_boun_dim
                   ipoin = commu % ghost_recv_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = tmp_rrecv(kk)
                   end do
                end do

             else
                call runend('UNKNOWN ORDER')
             end if

             ipass = 0
             deallocate(tmp_rrecv)
             deallocate(tmp_rsend)

          end if

       end if

    end if
#endif

  end subroutine PAR_GHOST_BOUNDARY_EXCHANGE_RP

  !----------------------------------------------------------------------
  !
  ! WAITALL
  !
  !----------------------------------------------------------------------

  subroutine PAR_WAITALL()
    integer(4) :: istat4
    integer(4) :: count4
#ifndef MPI_OFF
    integer(4) :: status4(MPI_STATUS_SIZE,1_ip)
#endif

#ifndef MPI_OFF
    if( IPARALL ) then
       count4 = 1_4
       call MPI_WAITALL(count4,ireq41,status4,istat4)
       deallocate(yy_non_blocking)
    end if
#endif
  end subroutine PAR_WAITALL

  !----------------------------------------------------------------------
  !
  ! PAR_ARRAY_EXCHANGE_IP: FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_ARRAY_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: ndofn
    integer(ip),                    intent(inout) :: xx(*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    integer(4),                     intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,nsize,jj,dom_i
    integer(ip)                                   :: ipoin,ini,kk,idofn
    integer(4)                                    :: istat4,nsize4,count4
    integer(4)                                    :: dom_i4,my_rank4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    integer(4),                        pointer    :: status4(:,:)
    integer(ip)                                   :: dom_j

#ifndef MPI_OFF
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_isend(commu % nneig * ndofn))
             allocate(tmp_irecv(commu % nneig * ndofn))
             !
             ! Save in temp_send
             !
             kk = 0
             do ii = 1,commu % nneig
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_isend(kk) = xx(idofn)
                   tmp_irecv(kk) = 0
                end do
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             istat4 = 0_4
             kk     = 0
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   ini   = ( ii - 1 ) * ndofn + 1
                   nsize = ndofn

                   nsize4 = int(nsize,4)
                   if( asynch ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           tmp_isend(ini:ini+nsize-1), nsize4, &
                           PAR_INTEGER,  dom_i4, 0_4,          &
                           PAR_COMM_TO_USE, ireq4(kk), istat4 )
                      kk = kk + 1
                      call MPI_Irecv(&
                           tmp_irecv(ini:ini+nsize-1), nsize4, &
                           PAR_INTEGER,  dom_i4, 0_4,          &
                           PAR_COMM_TO_USE, ireq4(kk), istat4 )
                   else
                      call MPI_Sendrecv(                       &
                           tmp_isend(ini:), nsize4,            &
                           PAR_INTEGER, dom_i4, 0_4,           &
                           tmp_irecv(ini:), nsize4,            &
                           PAR_INTEGER, dom_i4, 0_4,           &
                           PAR_COMM_TO_USE, status, istat4    )
                   end if
                   if( istat4 /= 0_4 ) call runend('PAR_INTERFACE_NODE_EXCHANGE_IP: MPI ERROR')
                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE,2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if( trim(what) == 'SUM' .or. trim(what) == 'ASSEMBLY' ) then
                !
                ! SUM
                !
                kk = 0
                do ii = 1,commu % nneig
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn) = xx(idofn) + tmp_irecv(kk)
                   end do
                end do

             else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' ) then
                !
                ! MAX
                !
                kk = 0
                do ii = 1,commu % nneig
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn) = max(xx(idofn),tmp_irecv(kk))
                   end do
                end do

             else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' ) then
                !
                ! MIN
                !
                kk = 0
                do ii = 1,commu % nneig
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn) = min(xx(idofn),tmp_irecv(kk))
                   end do
                end do

             else if( trim(what) == 'TAKE MIN' ) then
                !
                ! TAKE MIN
                !
                call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank4)
                kk = 0
                do ii = 1,commu % nneig
                   dom_i = commu % neights(ii)
                   if( my_rank4 < dom_i ) then
                      do idofn = 1,ndofn
                         kk = kk + 1
                         xx(idofn) = tmp_irecv(kk)
                      end do
                   else
                      kk = kk + ndofn
                   end if
                end do

             else
                call runend('UNKNOWN ORDER')
             end if

             ipass = 0
             deallocate(tmp_irecv)
             deallocate(tmp_isend)

          end if

       end if

    end if
#endif

  end subroutine PAR_ARRAY_EXCHANGE_IP

  !----------------------------------------------------------------------
  !
  ! PAR_ALL_TO_ALL_ARRAY_OPERATION_IP: FOR INTEGERS
  ! INPUT:  XX(NDOFN) for all slaves
  ! OUTPUT: XX(IDOFN) = XX(IDOFN) if my ranks is the max who have XX(IDOFN) /= 0
  !                   = 0 otherwise
  !         For the master XX(IDOFN) = rank of partition who has XX(IDOFN) /= 0
  !
  !
  !----------------------------------------------------------------------

  subroutine PAR_ALL_TO_ALL_ARRAY_OPERATION_IP(ndofn,xx,what,wherein)
    implicit none
    integer(ip),         intent(in)    :: ndofn
    integer(ip),         intent(inout) :: xx(*)
    character(*),        intent(in)    :: what
    character(*),        optional      :: wherein
    integer(ip)                        :: idofn
    integer(4)                         :: my_rank4
    integer(4)                         :: PAR_COMM_TO_USE4
    integer(ip),         pointer       :: lranks(:)

#ifndef MPI_OFF
    if( IPARALL ) then

       if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE4)
       else
          call PAR_DEFINE_COMMUNICATOR('IN MY CODE',PAR_COMM_TO_USE4)
       end if
       call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE4,my_rank4)

       allocate( lranks(ndofn) )

       if( trim(what) == 'CHOOSE ONLY ONE' ) then
          if( INOTMASTER ) then
             do idofn = 1,ndofn
                if( xx(idofn) > 0 ) then
                   lranks(idofn) = int(my_rank4,ip)
                else
                   lranks(idofn) = 0_ip
                end if
             end do
          end if

          if( present(wherein) ) then
             call PAR_MAX(ndofn,lranks,wherein)
          else
             call PAR_MAX(ndofn,lranks)
          end if

          if( INOTMASTER ) then
             do idofn = 1,ndofn
                if( my_rank4 /= lranks(idofn) ) xx(idofn) = 0
             end do
          else
             do idofn = 1,ndofn
                xx(idofn) = lranks(idofn)
             end do
          end if

       end if

       deallocate( lranks )

    end if
#endif

  end subroutine PAR_ALL_TO_ALL_ARRAY_OPERATION_IP

  !----------------------------------------------------------------------
  !
  ! PAR_POINT_TO_POINT_ARRAY_OPERATION_IP: FOR INTEGERS
  ! INPUT/OUTPUT:  XX_SEND_RECV(NDOFN) for all slaves
  ! Operate on arrays between neighbors in the communicator.
  !
  !----------------------------------------------------------------------

  subroutine PAR_POINT_TO_POINT_ARRAY_OPERATION_IP_1(xx_send_recv,wherein,what)
    implicit none
    integer(ip),          pointer,  intent(inout) :: xx_send_recv(:)
    character(*),                   intent(in)    :: wherein
    character(*),                   intent(in)    :: what
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    integer(4)                                    :: PAR_COMM_TO_USE

    if( ISEQUEN ) return

    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)

    if( .not. associated(xx_send_recv) ) then
       return
    else
       ndofn = size(xx_send_recv)
       if( ndofn > 0 ) call PAR_POINT_TO_POINT_ARRAY_OPERATION_IP(ndofn,xx_send_recv,PAR_COMM_TO_USE,commu,what)
    end if

  end subroutine PAR_POINT_TO_POINT_ARRAY_OPERATION_IP_1

  subroutine PAR_POINT_TO_POINT_ARRAY_OPERATION_IP(ndofn,xx_send_recv,PAR_COMM_TO_USE,commu,what)
    implicit none
    integer(ip),                    intent(in)    :: ndofn
    integer(ip),                    intent(inout) :: xx_send_recv(*)
    integer(4),                     intent(in)    :: PAR_COMM_TO_USE
    type(comm_data_par),            intent(in)    :: commu
    character(*),                   intent(in)    :: what
    integer(ip)                                   :: dom_i,ineig,nsend,nrecv,ii
    integer(ip),          pointer                 :: xx_recv(:)
    integer(4)                                    :: istat4,dom_i4,my_rank4
    integer(4)                                    :: ndofn4
    !
    ! Define communicator
    !
    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank4)
    !
    ! Allocate
    !
    nullify(xx_recv)
    allocate(xx_recv(ndofn))
    ndofn4 = int(ndofn,4)
    !
    ! Send/receive and operate
    !
    do ineig = 1,commu % nneig

       dom_i  = commu % neights(ineig)
       dom_i4 = int(dom_i,4)

#ifndef MPI_OFF
       call MPI_Sendrecv(                      &
            xx_send_recv(1:ndofn), ndofn4,     &
            PAR_INTEGER, dom_i4, 0_4,          &
            xx_recv(1:ndofn), ndofn4,          &
            PAR_INTEGER, dom_i4, 0_4,          &
            PAR_COMM_TO_USE, status, istat4    )
#endif
       !
       ! Operate on array
       !
       if(      trim(what) == 'MAX' ) then
          !
          ! MAX
          !
          do ii = 1,ndofn
             xx_send_recv(ii) = max(xx_send_recv(ii),xx_recv(ii))
          end do
       else if( trim(what) == 'MIN' ) then
          !
          ! MIN
          !
          do ii = 1,ndofn
             xx_send_recv(ii) = min(xx_send_recv(ii),xx_recv(ii))
          end do
       else if( trim(what) == 'SUM' ) then
          !
          ! SUM
          !
          do ii = 1,ndofn
             xx_send_recv(ii) = xx_send_recv(ii) + xx_recv(ii)
          end do
       else if( trim(what) == 'MIN RANK OR NEGATIVE' ) then
          !
          ! At the end, only one partition will have the positive sign
          ! if the value is positive
          !
          if( my_rank4 > dom_i ) then
             do ii = 1,ndofn
                if( xx_send_recv(ii) > 0 .and. xx_recv(ii) > 0 ) then
                   xx_send_recv(ii) = -abs(xx_send_recv(ii))
                end if
             end do
          end if
       end if

    end do
    !
    ! Deallocate
    !
    deallocate(xx_recv)

  end subroutine PAR_POINT_TO_POINT_ARRAY_OPERATION_IP

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP
  !
  !----------------------------------------------------------------------

  subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_00(n,xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_00

  subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_0(n,xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_0

  subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_1(xx,what,wherein,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = 1
       if( size(xx,1) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_1

  subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_2(xx,what,wherein,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = size(xx,1)
    if( size(xx,2) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_2

  subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_3(xx,what,wherein,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = size(xx,1)*size(xx,2)
       if( size(xx,3) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_3

  subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_2b(xx,what,commu,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    if( associated(xx) ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
       call PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP_2b

  !----------------------------------------------------------------------
  !
  ! PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP: NODE ASSEMBLY FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: ndofn
    real(rp),                       intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    integer(4),                     intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,jj,dom_i
    integer(ip)                                   :: nsize_send,nsize_recv
    integer(ip)                                   :: ipoin,ini_send,ini_recv,kk,idofn
    integer(4)                                    :: nsize_send4,nsize_recv4
    integer(4)                                    :: istat4,count4
    integer(4)                                    :: dom_i4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    integer(4),                        pointer    :: status4(:,:)
    integer(ip)                                   :: dom_j

#ifndef MPI_OFF
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_ELEMENT_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( commu % ghost_send_elem_dim == -1 .or. commu % ghost_recv_elem_dim == -1 ) &
                  call runend('FRINGE NODE EXCHANGE NOT COMPUTED')
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_rsend(commu % ghost_recv_elem_dim * ndofn))
             allocate(tmp_rrecv(commu % ghost_send_elem_dim * ndofn))
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,commu % ghost_recv_elem_dim
                ipoin = commu % ghost_recv_elem_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_rsend(kk) = xx(idofn,ipoin)
                end do
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             kk = 0
             istat4 = 0_4
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   dom_i      = commd % neights(ii)
                   ini_send   = ndofn * ( commd % ghost_recv_elem_size(ii)   -1 ) + 1
                   nsize_send = ndofn * ( commd % ghost_recv_elem_size(ii+1) -1 ) + 1 - ini_send
                   ini_recv   = ndofn * ( commd % ghost_send_elem_size(ii)   -1 ) + 1
                   nsize_recv = ndofn * ( commd % ghost_send_elem_size(ii+1) -1 ) + 1 - ini_recv

                   nsize_send4 = int(nsize_send,4)
                   nsize_recv4 = int(nsize_recv,4)

                   if( asynch ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           tmp_rsend(ini_send:ini_send+nsize_send-1), nsize_send4, &
                           MPI_DOUBLE_PRECISION,  dom_i4, 0_4,                     &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                      kk = kk + 1
                      call MPI_Irecv(&
                           tmp_rrecv(ini_recv:ini_recv+nsize_recv-1), nsize_recv4, &
                           MPI_DOUBLE_PRECISION,  dom_i4, 0_4,                     &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                   else
                      if( nsize_recv /= 0 .and. nsize_send == 0 ) then
                         call MPI_Recv(                          &
                              tmp_rrecv(ini_recv:), nsize_recv4, &
                              MPI_DOUBLE_PRECISION, dom_i4, 0_4, &
                              PAR_COMM_TO_USE, status, istat4    )
                      else if( nsize_recv == 0 .and. nsize_send /= 0 ) then
                         call MPI_Send(                          &
                              tmp_rsend(ini_send:), nsize_send4, &
                              MPI_DOUBLE_PRECISION, dom_i4, 0_4, &
                              PAR_COMM_TO_USE, istat4            )
                      else if( nsize_recv /= 0 .and. nsize_send /= 0 ) then
                         call MPI_Sendrecv(                         &
                              tmp_rsend(ini_send:), nsize_send4,    &
                              MPI_DOUBLE_PRECISION, dom_i4, 0_4,    &
                              tmp_rrecv(ini_recv:), nsize_recv4,    &
                              MPI_DOUBLE_PRECISION, dom_i4, 0_4,    &
                              PAR_COMM_TO_USE, status, istat4       )
                      end if
                   end if
                   if( istat4 /= 0_4 ) call runend('PAR_GHOST_ELEMENT_EXCHANGE_RP: MPI ERROR')

                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE,2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             if( istat4 /= 0 ) call runend('WRONG SEND/RECEIVE')
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if( trim(what) == 'SUM' .or. trim(what) == 'ASSEMBLY' ) then
                !
                ! SUM
                !
                kk = 0
                do jj = 1,commu % ghost_send_elem_dim
                   ipoin = commu % ghost_send_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = xx(idofn,ipoin) + tmp_rrecv(kk)
                   end do
                end do

             else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' ) then
                !
                ! MAX
                !
                kk = 0
                do jj = 1,commu % ghost_send_elem_dim
                   ipoin = commu % ghost_send_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = max(xx(idofn,ipoin),tmp_rrecv(kk))
                   end do
                end do

             else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' ) then
                !
                ! MIN
                !
                kk = 0
                do jj = 1,commu % ghost_send_elem_dim
                   ipoin = commu % ghost_send_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = min(xx(idofn,ipoin),tmp_rrecv(kk))
                   end do
                end do

             else if( trim(what) == 'REPLACE' .or. trim(what) == 'SUBSTITUTE' ) then
                !
                ! Replace value on my fringe nodes according to what I have received
                !
                kk = 0
                do jj = 1,commu % ghost_send_elem_dim
                   ipoin = commu % ghost_send_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = tmp_rrecv(kk)
                   end do
                end do

             else
                call runend('UNKNOWN ORDER')
             end if

             ipass = 0
             deallocate(tmp_rrecv)
             deallocate(tmp_rsend)

          end if

       end if

    end if
#endif

  end subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_RP

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_INTERFACE_EDGE_EXCHANGE_IP
  !
  !----------------------------------------------------------------------

  subroutine PAR_INTERFACE_EDGE_EXCHANGE_IP_0(n,xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: n
    integer(ip),                    intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_INTERFACE_NODE_EXCHANGE_IP(ON_EDGES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_EDGE_EXCHANGE_IP_0

  subroutine PAR_INTERFACE_EDGE_EXCHANGE_IP_1(xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),          pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    integer(4)                                    :: PAR_COMM_TO_USE

   if( INOTSLAVE ) return
    ndofn = 1
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_INTERFACE_NODE_EXCHANGE_IP(ON_EDGES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_EDGE_EXCHANGE_IP_1

  subroutine PAR_INTERFACE_EDGE_EXCHANGE_IP_2(xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    ndofn = size(xx,1)
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_INTERFACE_NODE_EXCHANGE_IP(ON_EDGES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_EDGE_EXCHANGE_IP_2

  subroutine PAR_INTERFACE_EDGE_EXCHANGE_IP_3(xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),          pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    if( size(xx,3) <= 2 ) then
       ndofn = size(xx,1)
    else
       ndofn = size(xx,1)*size(xx,2)
    end if
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_INTERFACE_NODE_EXCHANGE_IP(ON_EDGES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_EDGE_EXCHANGE_IP_3

  subroutine PAR_INTERFACE_EDGE_EXCHANGE_IP_2b(xx,what,commu,wsynch,dom_k)
    implicit none
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),  pointer,  intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    ndofn = size(xx,1)
    PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
    call PAR_INTERFACE_NODE_EXCHANGE_IP(ON_EDGES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_EDGE_EXCHANGE_IP_2b

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_INTERFACE_EDGE_EXCHANGE_RP
  !
  !----------------------------------------------------------------------

  subroutine PAR_INTERFACE_EDGE_EXCHANGE_RP_00(n,xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_INTERFACE_NODE_EXCHANGE_RP(ON_EDGES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_EDGE_EXCHANGE_RP_00

  subroutine PAR_INTERFACE_EDGE_EXCHANGE_RP_0(n,xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_INTERFACE_NODE_EXCHANGE_RP(ON_EDGES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_EDGE_EXCHANGE_RP_0

  subroutine PAR_INTERFACE_EDGE_EXCHANGE_RP_1(xx,what,wherein,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = 1
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_INTERFACE_NODE_EXCHANGE_RP(ON_EDGES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_EDGE_EXCHANGE_RP_1

  subroutine PAR_INTERFACE_EDGE_EXCHANGE_RP_2(xx,what,wherein,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = size(xx,1)
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_INTERFACE_NODE_EXCHANGE_RP(ON_EDGES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_EDGE_EXCHANGE_RP_2

  subroutine PAR_INTERFACE_EDGE_EXCHANGE_RP_3(xx,what,wherein,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    if( what(len(what):len(what)) /= 'I' ) then
       if( size(xx,3) <= 2 ) then
          ndofn = size(xx,1)
       else
          ndofn = size(xx,1)*size(xx,2)
       end if
    end if
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_INTERFACE_NODE_EXCHANGE_RP(ON_EDGES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_EDGE_EXCHANGE_RP_3

  subroutine PAR_INTERFACE_EDGE_EXCHANGE_RP_2b(xx,what,commu,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    ndofn = size(xx,1)
    PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
    call PAR_INTERFACE_NODE_EXCHANGE_RP(ON_EDGES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_EDGE_EXCHANGE_RP_2b

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_INTERFACE_EDGE_EXCHANGE_LG
  !
  !----------------------------------------------------------------------

  subroutine PAR_INTERFACE_EDGE_EXCHANGE_LG_1(xx,what,wherein,wsynch,dom_k)
    implicit none
    logical(lg),          pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    type(comm_data_par),  pointer                 :: commu
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE ) return
    ndofn = 1
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_INTERFACE_NODE_EXCHANGE_LG(ON_EDGES,ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_INTERFACE_EDGE_EXCHANGE_LG_1

  !----------------------------------------------------------------------
  !
  ! PAR_INTERFACE_MATRIX_EXCHANGE: EXCHANGE MATRIX ON INTERFACE NODES
  !
  !----------------------------------------------------------------------

  subroutine PAR_INTERFACE_MATRIX_EXCHANGE_WHERE(ndofn,aa,wherein)
    implicit none
    integer(ip),          intent(in)    :: ndofn
    real(rp),             intent(inout) :: aa(ndofn,ndofn,*)
    character(*),         intent(in)    :: wherein
    integer(4)                          :: PAR_COMM_TO_USE
    type(comm_data_par),  pointer       :: commu
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_INTERFACE_MATRIX_EXCHANGE_RP(ndofn,aa,commu)
  end subroutine PAR_INTERFACE_MATRIX_EXCHANGE_WHERE

  subroutine PAR_INTERFACE_MATRIX_EXCHANGE_COMM_SHAPE(ndofn,aa,commu)
    implicit none
    integer(ip),          intent(in)    :: ndofn
    real(rp),             intent(inout) :: aa(*)
    type(comm_data_par),  intent(in)    :: commu
    call PAR_INTERFACE_MATRIX_EXCHANGE_RP(ndofn,aa,commu)
  end subroutine PAR_INTERFACE_MATRIX_EXCHANGE_COMM_SHAPE

  subroutine PAR_INTERFACE_MATRIX_EXCHANGE_COMM(ndofn,aa,commu)
    implicit none
    integer(ip),          intent(in)    :: ndofn
    real(rp),             intent(inout) :: aa(ndofn,ndofn,*)
    type(comm_data_par),  intent(in)    :: commu
    call PAR_INTERFACE_MATRIX_EXCHANGE_RP(ndofn,aa,commu)
  end subroutine PAR_INTERFACE_MATRIX_EXCHANGE_COMM

  subroutine PAR_INTERFACE_MATRIX_EXCHANGE_RP(ndofn,aa,commu)
    integer(ip),          intent(in)    :: ndofn
    real(rp),             intent(inout) :: aa(ndofn,ndofn,*)
    type(comm_data_par),  intent(in)    :: commu
    integer(ip)                         :: ineig,dom_i,bsize,jj,kk,ll
    integer(ip)                         :: idofn,jdofn,ii,iz
    integer(ip)                         :: PAR_CURRENT_RANK
    integer(ip)                         :: PAR_CURRENT_SIZE
    type(r1p),            pointer       :: loc_sparr1(:),loc_rparr1(:)
    !
    ! Allocate memory
    !
    nullify(loc_sparr1)
    nullify(loc_rparr1)
    call memory_alloca(par_memor,'LOC_SPARR1','PAR_INTERFACE_MATRIX_EXCHANGE',loc_sparr1,commu % nneig)
    call memory_alloca(par_memor,'LOC_RPARR1','PAR_INTERFACE_MATRIX_EXCHANGE',loc_rparr1,commu % nneig)
    do ineig = 1,commu % nneig
       bsize = ndofn * ndofn *  commu % bound_matrix(ineig) % nzdom
       call memory_alloca(par_memor,'LOC_SPARR1','PAR_INTERFACE_MATRIX_EXCHANGE',loc_sparr1(ineig) % a,bsize)
       call memory_alloca(par_memor,'LOC_RPARR1','PAR_INTERFACE_MATRIX_EXCHANGE',loc_rparr1(ineig) % a,bsize)
    end do
    !
    ! Exchange matrix coefficients
    !
    call PAR_COMM_RANK_AND_SIZE(COMMU % PAR_COMM_WORLD,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)
    call PAR_START_NON_BLOCKING_COMM(1_ip,PAR_CURRENT_SIZE)
    call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
    do ineig = 1,commu % nneig
       dom_i = commu % neights(ineig)
       bsize = ndofn * ndofn *  commu % bound_matrix(ineig) % nzdom
       jj = 0
       kk = 0
       do ii = commu % bound_size(ineig),commu % bound_size(ineig+1)-1
          jj = jj + 1
          do ll = 1,commu % bound_matrix(ineig) % nzdom_ii(jj)
             iz = commu % bound_matrix(ineig) % ja(jj) % l(ll)
             do idofn = 1,ndofn
                do jdofn = 1,ndofn
                   kk = kk + 1
                   loc_sparr1(ineig) % a(kk) = aa(idofn,jdofn,iz)
                end do
             end do
          end do
       end do
       if( bsize > 0 ) call PAR_SEND_RECEIVE(bsize,bsize,loc_sparr1(ineig)%a,loc_rparr1(ineig)%a,'IN MY CODE',dom_i,'NON BLOCKING')
       !if( bsize > 0 ) call PAR_SEND_RECEIVE(bsize,bsize,loc_sparr1(ineig)%a,loc_rparr1(ineig)%a,'IN MY CODE',dom_i)
    end do
    call PAR_END_NON_BLOCKING_COMM(1_ip)
    !
    ! Assemble the contribution of my neighbor
    !
    do ineig = 1,commu % nneig
       jj = 0
       kk = 0
       do ii = commu % bound_size(ineig),commu % bound_size(ineig+1)-1
          jj = jj + 1
          do ll = 1,commu % bound_matrix(ineig) % nzdom_ii(jj)
             iz = commu % bound_matrix(ineig) % ja(jj) % l(ll)
             do idofn = 1,ndofn
                do jdofn = 1,ndofn
                   kk = kk + 1
                   aa(idofn,jdofn,iz) = aa(idofn,jdofn,iz) + loc_rparr1(ineig) % a(kk)
                end do
             end do
          end do
       end do

    end do

    call memory_deallo( par_memor,'LOC_SPARR1','PAR_INTERFACE_MATRIX_EXCHANGE',loc_sparr1)
    call memory_deallo( par_memor,'LOC_RPARR1','PAR_INTERFACE_MATRIX_EXCHANGE',loc_rparr1)

  end subroutine PAR_INTERFACE_MATRIX_EXCHANGE_RP

  !----------------------------------------------------------------------
  !
  ! PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE: EXCHANGE MATRIX WITH HALOS ON INTERFACE NODES
  !
  !----------------------------------------------------------------------
  ! ojo aa lo tengo que partir en 2 un aasend y un aarecv  el primero es sin halos y el segundo con !!

  subroutine PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE_RP_2(ndofn,aa_send,aa_recv,commu)
    implicit none
    integer(ip),          intent(in)             :: ndofn
    real(rp),             intent(in)             :: aa_send(ndofn,ndofn,*)
    real(rp),             intent(inout)          :: aa_recv(ndofn,ndofn,*)
    type(comm_data_par),  intent(in),   optional :: commu

    if( present(commu) ) then
       call PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE_RP(ndofn,aa_send,aa_recv,commu)
    else
       call PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE_RP(ndofn,aa_send,aa_recv,PAR_COMM_MY_CODE_ARRAY(1))
    end if

  end subroutine PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE_RP_2

  subroutine PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE_RP_1(ndofn,aa_send,aa_recv,commu)
    implicit none
    integer(ip),          intent(in)             :: ndofn
    real(rp),             intent(in)             :: aa_send(*)
    real(rp),             intent(inout)          :: aa_recv(*)
    type(comm_data_par),  intent(in),   optional :: commu

    if( present(commu) ) then
       call PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE_RP(ndofn,aa_send,aa_recv,commu)
    else
       call PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE_RP(ndofn,aa_send,aa_recv,PAR_COMM_MY_CODE_ARRAY(1))
    end if

  end subroutine PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE_RP_1


  subroutine PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE_RP(ndofn,aa_send,aa_recv,commu)

    integer(ip),          intent(in)    :: ndofn
    real(rp),             intent(in)    :: aa_send(ndofn,ndofn,*)
    real(rp),             intent(inout) :: aa_recv(ndofn,ndofn,*)
    type(comm_data_par),  intent(in)    :: commu
    integer(ip)                         :: ineig,dom_i,bsize,jj,kk,ll
    integer(ip)                         :: idofn,jdofn,ii,iz,ipoin
    integer(ip)                         :: PAR_CURRENT_RANK
    integer(ip)                         :: PAR_CURRENT_SIZE
    type(r1p),            pointer       :: loc_sparr1(:),loc_rparr1(:)
    integer(ip)                         :: ssize,rsize

    !
    ! Allocate memory
    !
    nullify(loc_sparr1)
    nullify(loc_rparr1)
    call memory_alloca(par_memor,'LOC_SPARR1','PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE',loc_sparr1,commu % nneig)
    call memory_alloca(par_memor,'LOC_RPARR1','PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE',loc_rparr1,commu % nneig)
    do ineig = 1,commu % nneig
       ssize = ndofn * ndofn * commu % bound_mat_halo_send(ineig) % nzdom
       rsize = ndofn * ndofn * commu % bound_mat_halo_recv(ineig) % nzdom
       call memory_alloca(par_memor,'LOC_SPARR1','PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE',loc_sparr1(ineig) % a,max(1_ip,ssize))
       call memory_alloca(par_memor,'LOC_RPARR1','PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE',loc_rparr1(ineig) % a,max(1_ip,rsize))
    end do
    !
    ! Exchange matrix coefficients
    !
    call PAR_COMM_RANK_AND_SIZE(COMMU % PAR_COMM_WORLD,PAR_CURRENT_RANK,PAR_CURRENT_SIZE)
    call PAR_START_NON_BLOCKING_COMM(1_ip,PAR_CURRENT_SIZE)
    call PAR_SET_NON_BLOCKING_COMM_NUMBER(1_ip)
    do ineig = 1,commu % nneig
       dom_i = commu % neights(ineig)
       ssize = ndofn * ndofn *  commu % bound_mat_halo_send(ineig) % nzdom
       rsize = ndofn * ndofn *  commu % bound_mat_halo_recv(ineig) % nzdom
       jj = 0
       kk = 0
       do ii = commu % bound_size(ineig),commu % bound_size(ineig+1)-1     ! These are all the nodes that are in the boundary with my neighbour ineig
          jj = jj + 1
          do ll = 1,commu % bound_mat_halo_send(ineig) % nzdom_ii(jj)      ! These are the number of conections to a certain node
             iz = commu % bound_mat_halo_send(ineig) % ja(jj) % l(ll)
             do idofn = 1,ndofn
                do jdofn = 1,ndofn
                   kk = kk + 1
                   loc_sparr1(ineig) % a(kk) = aa_send(idofn,jdofn,iz)
                end do
             end do
          end do
       end do
       if( ssize > 0 .or. rsize > 0 ) call PAR_SEND_RECEIVE(ssize,rsize,loc_sparr1(ineig)%a,loc_rparr1(ineig)%a,'IN MY CODE',dom_i,'NON BLOCKING')
    end do
    call PAR_END_NON_BLOCKING_COMM(1_ip)
    !
    ! Assemble the contribution of my neighbor
    !
    do ineig = 1,commu % nneig
       jj = 0
       kk = 0
       do ii = commu % bound_size(ineig),commu % bound_size(ineig+1)-1
          jj = jj + 1
          do ll = 1,commu % bound_mat_halo_recv(ineig) % nzdom_ii(jj)
             iz = commu % bound_mat_halo_recv(ineig) % ja(jj) % l(ll)
             do idofn = 1,ndofn
                do jdofn = 1,ndofn
                   kk = kk + 1
                   aa_recv(idofn,jdofn,iz) = aa_recv(idofn,jdofn,iz) + loc_rparr1(ineig) % a(kk)
                end do
             end do
          end do
       end do

    end do

    call memory_deallo( par_memor,'LOC_SPARR1','PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE',loc_sparr1)
    call memory_deallo( par_memor,'LOC_RPARR1','PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE',loc_rparr1)

  end subroutine PAR_INTERFACE_MATRIX_W_HALOS_EXCHANGE_RP



  subroutine PAR_INTERFACE_NODE_EXCHANGE_VALUE(ndofn,exlen,sendb,recvb,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: ndofn
    integer(4)                                    :: exlen
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    real(rp), pointer                             :: sendb(:),recvb(:)
    logical(lg)                                   :: asynch
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu
    integer(4),                        pointer    :: status4(:,:)
    integer(4)                                    :: dom_i,istat4,dom_i4,dom_j,ii,kk,count4,nsize4
    integer(ip)                                   :: bound_dim,nsize,ini
    integer(ip),                       pointer    :: bound_perm(:)
    integer(ip),                       pointer    :: bound_size(:)


#ifndef MPI_OFF

    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)

    if( present(dom_k) ) then
       dom_j = dom_k
    else
       dom_j = 0
    end if

    ! On nodes or edges (this one only on nodes )
    !if( onwhat == ON_NODES ) then
    bound_dim  =  commu % bound_dim
    bound_perm => commu % bound_perm
    bound_size => commu % bound_size
    !else
    !   bound_dim  =  commu % bedge_dim
    !   bound_perm => commu % bedge_perm
    !   bound_size => commu % bedge_size
    !end if

    if( present(wsynch) ) then
       if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
          asynch = .false.
       else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
          asynch = .true.
       else
          call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
       end if
    else
       asynch = .false.
    end if
    kk = 0
    if( ISLAVE ) then

       if( asynch ) allocate(ireq4(commu % nneig*2))
       ! Send    temp_send
       ! Receive temp_recv
       !
       istat4 = 0_4

       do ii = 1,commu % nneig

          dom_i  = commu % neights(ii)
          dom_i4 = int(dom_i,4)

          if( dom_j == 0 .or. dom_j == dom_i ) then

             ini   = ndofn * ( bound_size(ii)   - 1 ) + 1
             nsize = ndofn * ( bound_size(ii+1) - 1 ) + 1 - ini

             nsize4 = int(nsize,4)

             !print *,kfl_paral,'--->',ini,nsize

             if( asynch ) then
                kk = kk + 1
                call MPI_Isend(&
                     sendb(ini), nsize4, &
                     MPI_DOUBLE_PRECISION,  dom_i4, 0_4, &
                     PAR_COMM_TO_USE, ireq4(kk), istat4 )
                kk = kk + 1
                call MPI_Irecv(&
                     recvb(ini), nsize4, &
                     MPI_DOUBLE_PRECISION,  dom_i4, 0_4, &
                     PAR_COMM_TO_USE, ireq4(kk), istat4 )
             else
                call MPI_Sendrecv(                       &
                     sendb(ini), nsize4, &
                     MPI_DOUBLE_PRECISION, dom_i4, 0_4,  &
                     recvb(ini), nsize4, &
                     MPI_DOUBLE_PRECISION, dom_i4, 0_4,  &
                     PAR_COMM_TO_USE, status, istat4    )
             end if
             if( istat4 /= 0_4 ) call runend('PAR_INTERFACE_NODE_EXCHANGE_RP: MPI ERROR')

          end if

       end do

    end if

    if( asynch ) then
       count4 = 2*int(commu % nneig,4)
       allocate( status4(MPI_STATUS_SIZE,2*commu % nneig) )
       CALL MPI_WAITALL(count4,ireq4,status4,istat4)
       deallocate( status4 )
       deallocate(ireq4)
    end if

#endif

  end subroutine PAR_INTERFACE_NODE_EXCHANGE_VALUE

  subroutine COMMUNICATOR_VALUE(bound_dim,bound_perm,bound_size,wherein)
    implicit none
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu
    integer(ip)                                   :: bound_dim
    integer(ip),                       pointer    :: bound_perm(:)
    integer(ip),                       pointer    :: bound_size(:)
    character(*),                   intent(in)    :: wherein

#ifndef MPI_OFF

    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)


    ! On nodes or edges (this one only on nodes )
    !if( onwhat == ON_NODES ) then
    bound_dim  =  commu % bound_dim
    bound_perm => commu % bound_perm
    bound_size => commu % bound_size
    !else
    !   bound_dim  =  commu % bedge_dim
    !   bound_perm => commu % bedge_perm
    !   bound_size => commu % bedge_size
    !end if

#endif

  end subroutine COMMUNICATOR_VALUE

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_FROM_GHOST_NODE_EXCHANGE_RP
  !
  !----------------------------------------------------------------------

  subroutine PAR_FROM_GHOST_NODE_EXCHANGE_IP_00(n,xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: n
    integer(ip),                    intent(inout) :: xx(*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) then
       return
    else if( n > 0 ) then
       ndofn = n
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_NODE_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_NODE_EXCHANGE_IP_00

  subroutine PAR_FROM_GHOST_NODE_EXCHANGE_IP_0(n,xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: n
    integer(ip),                    intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_FROM_GHOST_NODE_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_FROM_GHOST_NODE_EXCHANGE_IP_0

  subroutine PAR_FROM_GHOST_NODE_EXCHANGE_IP_1(xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),          pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) then
       return
    else if( associated(xx) ) then
       ndofn = 1
       if( size(xx,1) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_NODE_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_NODE_EXCHANGE_IP_1

  subroutine PAR_FROM_GHOST_NODE_EXCHANGE_IP_2(xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) then
       return
    else if( associated(xx) ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_NODE_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_NODE_EXCHANGE_IP_2

  subroutine PAR_FROM_GHOST_NODE_EXCHANGE_IP_3(xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),          pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) then
       return
    else if( associated(xx) ) then
       ndofn = size(xx,1)*size(xx,2)
       if( size(xx,3) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_NODE_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_NODE_EXCHANGE_IP_3

  subroutine PAR_FROM_GHOST_NODE_EXCHANGE_IP_2b(xx,what,commu,wsynch,dom_k)
    implicit none
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE ) then
       return
    else if( associated(xx) ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= npoin_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
       call PAR_FROM_GHOST_NODE_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_NODE_EXCHANGE_IP_2b

  !----------------------------------------------------------------------
  !
  ! PAR_FROM_GHOST_NODE_EXCHANGE_IP: NODE ASSEMBLY FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_FROM_GHOST_NODE_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

    implicit none
    integer(ip),                    intent(in)    :: ndofn
    integer(ip),                    intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    integer(4),                     intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,jj,dom_i
    integer(ip)                                   :: nsize_send,nsize_recv
    integer(ip)                                   :: ipoin,ini_send,ini_recv,kk,idofn
    integer(4)                                    :: nsize_send4,nsize_recv4
    integer(4)                                    :: istat4,count4
    integer(4)                                    :: dom_i4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    integer(4),                        pointer    :: status4(:,:)
    integer(ip)                                   :: dom_j

#ifndef MPI_OFF
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_ELEMENT_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( commu % ghost_send_node_dim == -1 .or. commu % ghost_recv_node_dim == -1 ) &
                  call runend('FRINGE NODE EXCHANGE NOT COMPUTED')
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_isend(commu % ghost_send_node_dim * ndofn))
             allocate(tmp_irecv(commu % ghost_recv_node_dim * ndofn))
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,commu % ghost_send_node_dim
                ipoin = commu % ghost_send_node_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_isend(kk) = xx(idofn,ipoin)
                end do
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             kk = 0
             istat4 = 0_4
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   dom_i      = commd % neights(ii)
                   ini_send   = ndofn * ( commd % ghost_send_node_size(ii)   -1 ) + 1
                   nsize_send = ndofn * ( commd % ghost_send_node_size(ii+1) -1 ) + 1 - ini_send
                   ini_recv   = ndofn * ( commd % ghost_recv_node_size(ii)   -1 ) + 1
                   nsize_recv = ndofn * ( commd % ghost_recv_node_size(ii+1) -1 ) + 1 - ini_recv

                   nsize_send4 = int(nsize_send,4)
                   nsize_recv4 = int(nsize_recv,4)

                   if( asynch ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           tmp_isend(ini_send:ini_send+nsize_send-1), nsize_send4, &
                           PAR_INTEGER, dom_i4, 0_4,                               &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                      kk = kk + 1
                      call MPI_Irecv(&
                           tmp_irecv(ini_recv:ini_recv+nsize_recv-1), nsize_recv4, &
                           PAR_INTEGER,  dom_i4, 0_4,                              &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                   else
                      if( nsize_recv /= 0 .and. nsize_send == 0 ) then
                         call MPI_Recv(                          &
                              tmp_irecv(ini_recv:), nsize_recv4, &
                              PAR_INTEGER, dom_i4, 0_4,          &
                              PAR_COMM_TO_USE, status, istat4    )
                      else if( nsize_recv == 0 .and. nsize_send /= 0 ) then
                         call MPI_Send(                          &
                              tmp_isend(ini_send:), nsize_send4, &
                              PAR_INTEGER, dom_i4, 0_4,          &
                              PAR_COMM_TO_USE, istat4            )
                      else if( nsize_recv /= 0 .and. nsize_send /= 0 ) then
                         call MPI_Sendrecv(                         &
                              tmp_isend(ini_send:), nsize_send4,    &
                              PAR_INTEGER, dom_i4, 0_4,             &
                              tmp_irecv(ini_recv:), nsize_recv4,    &
                              PAR_INTEGER, dom_i4, 0_4,             &
                              PAR_COMM_TO_USE, status, istat4       )
                      end if
                   end if
                   if( istat4 /= 0_4 ) call runend('PAR_FROM_GHOST_NODE_EXCHANGE_RP: MPI ERROR')

                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE,2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             if( istat4 /= 0 ) call runend('WRONG SEND/RECEIVE')
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if( trim(what) == 'SUM' .or. trim(what) == 'ASSEMBLY' ) then
                !
                ! SUM
                !
                kk = 0
                do jj = 1,commu % ghost_recv_node_dim
                   ipoin = commu % ghost_recv_node_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = xx(idofn,ipoin) + tmp_irecv(kk)
                   end do
                end do

             else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' ) then
                !
                ! MAX
                !
                kk = 0
                do jj = 1,commu % ghost_recv_node_dim
                   ipoin = commu % ghost_recv_node_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = max(xx(idofn,ipoin),tmp_irecv(kk))
                   end do
                end do

             else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' ) then
                !
                ! MIN
                !
                kk = 0
                do jj = 1,commu % ghost_recv_node_dim
                   ipoin = commu % ghost_recv_node_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = min(xx(idofn,ipoin),tmp_irecv(kk))
                   end do
                end do

             else if( trim(what) == 'REPLACE' .or. trim(what) == 'SUBSTITUTE' ) then
                !
                ! Replace value on my fringe nodes according to what I have received
                !
                kk = 0
                do jj = 1,commu % ghost_recv_node_dim
                   ipoin = commu % ghost_recv_node_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = tmp_irecv(kk)
                   end do
                end do

             else
                call runend('UNKNOWN ORDER')
             end if

             ipass = 0
             deallocate(tmp_irecv)
             deallocate(tmp_isend)

          end if

       end if

    end if
#endif

  end subroutine PAR_FROM_GHOST_NODE_EXCHANGE_IP

  subroutine PAR_COMM_SPLIT_SIMVIZ(PAR_COMM_SIM,PAR_COMM_VIZ,PAR_COMM_SIMVIZ,rank,mpisplit)

    implicit none
    integer(4) :: PAR_COMM_SIM, PAR_COMM_VIZ, PAR_COMM_SIMVIZ, rank, mpisplit
    integer(4) :: vizcolor,simcolor,istat4

#ifndef MPI_OFF

    if ( mod(rank+1,mpisplit) .ne. 0_ip ) then
        simcolor = 1
        vizcolor = MPI_UNDEFINED
     else
        simcolor = MPI_UNDEFINED
        vizcolor = 1
     end if

     call MPI_COMM_SPLIT(PAR_COMM_SIMVIZ,simcolor,rank,PAR_COMM_SIM,istat4)
     if( istat4 /= MPI_SUCCESS ) call runend('PAR_COMM_SPLIT_simviz: MPI ERROR 2')
     call MPI_COMM_SPLIT(PAR_COMM_SIMVIZ,vizcolor,rank,PAR_COMM_VIZ,istat4)
     if( istat4 /= MPI_SUCCESS ) call runend('PAR_COMM_SPLIT_simviz: MPI ERROR 2')

#endif

   end subroutine PAR_COMM_SPLIT_SIMVIZ


  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP
  !
  !----------------------------------------------------------------------

  subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_00(n,xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( ISLAVE ) then
       ndofn = n
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_00

  subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_0(n,xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: n
    real(rp),                       intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( ISLAVE ) then
       ndofn = n
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_0

  subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_1(xx,what,wherein,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( ISLAVE .and. associated(xx) ) then
       ndofn = 1
       if( size(xx,1) /= nboun_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_1

  subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_2(xx,what,wherein,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( ISLAVE .and. associated(xx) ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= nboun_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_2

  subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_3(xx,what,wherein,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( ISLAVE .and. associated(xx) ) then
       ndofn = size(xx,1)*size(xx,2)
       if( size(xx,3) /= nboun_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_3

  subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_2b(xx,what,commu,wsynch,dom_k)
    implicit none
    real(rp),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE

    if( ISLAVE .and. associated(xx) ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= nboun_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
       call PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP_2b

  !----------------------------------------------------------------------
  !
  ! PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP: NODE ASSEMBLY FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: ndofn
    real(rp),                       intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    integer(4),                     intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,jj,dom_i
    integer(ip)                                   :: nsize_send,nsize_recv
    integer(ip)                                   :: ipoin,ini_send,ini_recv,kk,idofn
    integer(4)                                    :: nsize_send4,nsize_recv4
    integer(4)                                    :: istat4,count4
    integer(4)                                    :: dom_i4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    integer(4),                        pointer    :: status4(:,:)
    integer(ip)                                   :: dom_j

#ifndef MPI_OFF
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_BOUNDARY_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( commu % ghost_send_boun_dim == -1 .or. commu % ghost_recv_boun_dim == -1 ) &
                  call runend('FRINGE NODE EXCHANGE NOT COMPUTED')
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_rsend(commu % ghost_recv_boun_dim * ndofn))
             allocate(tmp_rrecv(commu % ghost_send_boun_dim * ndofn))
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,commu % ghost_recv_boun_dim
                ipoin = commu % ghost_recv_boun_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_rsend(kk) = xx(idofn,ipoin)
                end do
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             kk = 0
             istat4 = 0_4
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   dom_i      = commd % neights(ii)
                   ini_send   = ndofn * ( commd % ghost_recv_boun_size(ii)   -1 ) + 1
                   nsize_send = ndofn * ( commd % ghost_recv_boun_size(ii+1) -1 ) + 1 - ini_send
                   ini_recv   = ndofn * ( commd % ghost_send_boun_size(ii)   -1 ) + 1
                   nsize_recv = ndofn * ( commd % ghost_send_boun_size(ii+1) -1 ) + 1 - ini_recv

                   nsize_send4 = int(nsize_send,4)
                   nsize_recv4 = int(nsize_recv,4)

                   if( asynch ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           tmp_rsend(ini_send:ini_send+nsize_send-1), nsize_send4, &
                           MPI_DOUBLE_PRECISION,  dom_i4, 0_4,                     &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                      kk = kk + 1
                      call MPI_Irecv(&
                           tmp_rrecv(ini_recv:ini_recv+nsize_recv-1), nsize_recv4, &
                           MPI_DOUBLE_PRECISION,  dom_i4, 0_4,                     &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                   else
                      if( nsize_recv /= 0 .and. nsize_send == 0 ) then
                         call MPI_Recv(                          &
                              tmp_rrecv(ini_recv:), nsize_recv4, &
                              MPI_DOUBLE_PRECISION, dom_i4, 0_4, &
                              PAR_COMM_TO_USE, status, istat4    )
                      else if( nsize_recv == 0 .and. nsize_send /= 0 ) then
                         call MPI_Send(                          &
                              tmp_rsend(ini_send:), nsize_send4, &
                              MPI_DOUBLE_PRECISION, dom_i4, 0_4, &
                              PAR_COMM_TO_USE, istat4            )
                      else if( nsize_recv /= 0 .and. nsize_send /= 0 ) then
                         call MPI_Sendrecv(                         &
                              tmp_rsend(ini_send:), nsize_send4,    &
                              MPI_DOUBLE_PRECISION, dom_i4, 0_4,    &
                              tmp_rrecv(ini_recv:), nsize_recv4,    &
                              MPI_DOUBLE_PRECISION, dom_i4, 0_4,    &
                              PAR_COMM_TO_USE, status, istat4       )
                      end if
                   end if
                   if( istat4 /= 0_4 ) call runend('PAR_GHOST_BOUNDARY_EXCHANGE_RP: MPI ERROR')

                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE,2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             if( istat4 /= 0 ) call runend('WRONG SEND/RECEIVE')
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if( trim(what) == 'SUM' .or. trim(what) == 'ASSEMBLY' ) then
                !
                ! SUM
                !
                kk = 0
                do jj = 1,commu % ghost_send_boun_dim
                   ipoin = commu % ghost_send_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = xx(idofn,ipoin) + tmp_rrecv(kk)
                   end do
                end do

             else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' ) then
                !
                ! MAX
                !
                kk = 0
                do jj = 1,commu % ghost_send_boun_dim
                   ipoin = commu % ghost_send_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = max(xx(idofn,ipoin),tmp_rrecv(kk))
                   end do
                end do

             else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' ) then
                !
                ! MIN
                !
                kk = 0
                do jj = 1,commu % ghost_send_boun_dim
                   ipoin = commu % ghost_send_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = min(xx(idofn,ipoin),tmp_rrecv(kk))
                   end do
                end do

             else if( trim(what) == 'REPLACE' .or. trim(what) == 'SUBSTITUTE' ) then
                !
                ! Replace value on my fringe nodes according to what I have received
                !
                kk = 0
                do jj = 1,commu % ghost_send_boun_dim
                   ipoin = commu % ghost_send_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = tmp_rrecv(kk)
                   end do
                end do

             else
                call runend('UNKNOWN ORDER')
             end if

             ipass = 0
             deallocate(tmp_rrecv)
             deallocate(tmp_rsend)

          end if

       end if

    end if
#endif

  end subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_RP

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP
  !
  !----------------------------------------------------------------------

  subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_00(n,xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: n
    integer(ip),                    intent(inout) :: xx(*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( ISLAVE  ) then
       ndofn = n
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_00

  subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_0(n,xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: n
    integer(ip),                    intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( ISLAVE ) then
       ndofn = n
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_0

  subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_1(xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),          pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( ISLAVE .and. associated(xx) ) then
       ndofn = 1
       if( size(xx,1) /= nboun_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_1

  subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_2(xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( ISLAVE .and. associated(xx) ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= nboun_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_2

  subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_3(xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),          pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( ISLAVE .and. associated(xx) ) then
       ndofn = size(xx,1)*size(xx,2)
       if( size(xx,3) /= nboun_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_3

  subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_2b(xx,what,commu,wsynch,dom_k)
    implicit none
    integer(ip),          pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE

    if( ISLAVE .and. associated(xx) ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= nboun_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
       call PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP_2b

  !----------------------------------------------------------------------
  !
  ! PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP: NODE ASSEMBLY FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: ndofn
    integer(ip),                    intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    integer(4),                     intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,jj,dom_i
    integer(ip)                                   :: nsize_send,nsize_recv
    integer(ip)                                   :: ipoin,ini_send,ini_recv,kk,idofn
    integer(4)                                    :: nsize_send4,nsize_recv4
    integer(4)                                    :: istat4,count4
    integer(4)                                    :: dom_i4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    integer(4),                        pointer    :: status4(:,:)
    integer(ip)                                   :: dom_j

#ifndef MPI_OFF
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_BOUNDARY_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( commu % ghost_send_boun_dim == -1 .or. commu % ghost_recv_boun_dim == -1 ) &
                  call runend('FRINGE NODE EXCHANGE NOT COMPUTED')
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_isend(commu % ghost_recv_boun_dim * ndofn))
             allocate(tmp_irecv(commu % ghost_send_boun_dim * ndofn))
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,commu % ghost_recv_boun_dim
                ipoin = commu % ghost_recv_boun_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_isend(kk) = xx(idofn,ipoin)
                end do
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             kk = 0
             istat4 = 0_4
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   dom_i      = commd % neights(ii)
                   ini_send   = ndofn * ( commd % ghost_recv_boun_size(ii)   -1 ) + 1
                   nsize_send = ndofn * ( commd % ghost_recv_boun_size(ii+1) -1 ) + 1 - ini_send
                   ini_recv   = ndofn * ( commd % ghost_send_boun_size(ii)   -1 ) + 1
                   nsize_recv = ndofn * ( commd % ghost_send_boun_size(ii+1) -1 ) + 1 - ini_recv

                   nsize_send4 = int(nsize_send,4)
                   nsize_recv4 = int(nsize_recv,4)

                   if( asynch ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           tmp_isend(ini_send:ini_send+nsize_send-1), nsize_send4, &
                           PAR_INTEGER,  dom_i4, 0_4,                     &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                      kk = kk + 1
                      call MPI_Irecv(&
                           tmp_irecv(ini_recv:ini_recv+nsize_recv-1), nsize_recv4, &
                           PAR_INTEGER,  dom_i4, 0_4,                     &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                   else
                      if( nsize_recv /= 0 .and. nsize_send == 0 ) then
                         call MPI_Recv(                          &
                              tmp_irecv(ini_recv:), nsize_recv4, &
                              PAR_INTEGER, dom_i4, 0_4, &
                              PAR_COMM_TO_USE, status, istat4    )
                      else if( nsize_recv == 0 .and. nsize_send /= 0 ) then
                         call MPI_Send(                          &
                              tmp_isend(ini_send:), nsize_send4, &
                              PAR_INTEGER, dom_i4, 0_4, &
                              PAR_COMM_TO_USE, istat4            )
                      else if( nsize_recv /= 0 .and. nsize_send /= 0 ) then
                         call MPI_Sendrecv(                         &
                              tmp_isend(ini_send:), nsize_send4,    &
                              PAR_INTEGER, dom_i4, 0_4,    &
                              tmp_irecv(ini_recv:), nsize_recv4,    &
                              PAR_INTEGER, dom_i4, 0_4,    &
                              PAR_COMM_TO_USE, status, istat4       )
                      end if
                   end if
                   if( istat4 /= 0_4 ) call runend('PAR_GHOST_BOUNDARY_EXCHANGE_IP: MPI ERROR')

                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE,2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             if( istat4 /= 0 ) call runend('WRONG SEND/RECEIVE')
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if( trim(what) == 'SUM' .or. trim(what) == 'ASSEMBLY' ) then
                !
                ! SUM
                !
                kk = 0
                do jj = 1,commu % ghost_send_boun_dim
                   ipoin = commu % ghost_send_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = xx(idofn,ipoin) + tmp_irecv(kk)
                   end do
                end do

             else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' ) then
                !
                ! MAX
                !
                kk = 0
                do jj = 1,commu % ghost_send_boun_dim
                   ipoin = commu % ghost_send_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = max(xx(idofn,ipoin),tmp_irecv(kk))
                   end do
                end do

             else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' ) then
                !
                ! MIN
                !
                kk = 0
                do jj = 1,commu % ghost_send_boun_dim
                   ipoin = commu % ghost_send_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = min(xx(idofn,ipoin),tmp_irecv(kk))
                   end do
                end do

             else if( trim(what) == 'REPLACE' .or. trim(what) == 'SUBSTITUTE' ) then
                !
                ! Replace value on my fringe nodes according to what I have received
                !
                kk = 0
                do jj = 1,commu % ghost_send_boun_dim
                   ipoin = commu % ghost_send_boun_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = tmp_irecv(kk)
                   end do
                end do

             else
                call runend('UNKNOWN ORDER')
             end if

             ipass = 0
             deallocate(tmp_irecv)
             deallocate(tmp_isend)

          end if

       end if

    end if
#endif

  end subroutine PAR_FROM_GHOST_BOUNDARY_EXCHANGE_IP

  !----------------------------------------------------------------------
  !
  ! Bridges to PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP
  !
  !----------------------------------------------------------------------

  subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_00(n,xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: n
    integer(ip),                       intent(inout) :: xx(*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_00

  subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_0(n,xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: n
    integer(ip),                       intent(inout) :: xx(n,*)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE ) return
    ndofn = n
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_0

  subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_1(xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),             pointer,  intent(inout) :: xx(:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE .or. .not. associated(xx) ) return
    if( associated(xx) ) then
       ndofn = 1
       if( size(xx,1) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_1

  subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_2(xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE .or. .not. associated(xx) ) return
    ndofn = size(xx,1)
    if( size(xx,2) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
    call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
    call PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)

  end subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_2

  subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_3(xx,what,wherein,wsynch,dom_k)
    implicit none
    integer(ip),             pointer,  intent(inout) :: xx(:,:,:)
    character(*),                   intent(in)    :: what
    character(*),                   intent(in)    :: wherein
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE
    type(comm_data_par),               pointer    :: commu

    if( INOTSLAVE .or. .not. associated(xx) ) return
    if( associated(xx) ) then
       ndofn = size(xx,1)*size(xx,2)
       if( size(xx,3) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE,commu)
       call PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_3

  subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_2b(xx,what,commu,wsynch,dom_k)
    implicit none
    integer(ip),             pointer,  intent(inout) :: xx(:,:)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ndofn
    integer(4)                                    :: PAR_COMM_TO_USE

    if( INOTSLAVE .or. .not. associated(xx) ) return
    if( associated(xx) ) then
       ndofn = size(xx,1)
       if( size(xx,2) /= nelem_2 ) call runend('PAR_COMMUNICATIONS: WRONG SIZE')
       PAR_COMM_TO_USE = commu % PAR_COMM_WORLD
       call PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    end if

  end subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP_2b

  !----------------------------------------------------------------------
  !
  ! PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP: NODE ASSEMBLY FOR INTEGERS
  !
  !----------------------------------------------------------------------

  subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP(ndofn,xx,what,commu,PAR_COMM_TO_USE,wsynch,dom_k)
    implicit none
    integer(ip),                    intent(in)    :: ndofn
    integer(ip),                       intent(inout) :: xx(ndofn,*)
    character(*),                   intent(in)    :: what
    type(comm_data_par),            intent(in)    :: commu
    integer(4),                     intent(in)    :: PAR_COMM_TO_USE
    character(*),         optional, intent(in)    :: wsynch
    integer(ip),          optional, intent(in)    :: dom_k
    integer(ip)                                   :: ii,jj,dom_i
    integer(ip)                                   :: nsize_send,nsize_recv
    integer(ip)                                   :: ipoin,ini_send,ini_recv,kk,idofn
    integer(4)                                    :: nsize_send4,nsize_recv4
    integer(4)                                    :: istat4,count4
    integer(4)                                    :: dom_i4
    logical(lg)                                   :: asynch
    integer(ip),                       save       :: ipass = 0
    integer(4),                        pointer    :: status4(:,:)
    integer(ip)                                   :: dom_j

#ifndef MPI_OFF
    if( IPARALL ) then
       !
       ! Passes
       !
       ipass = ipass + 1
       if( present(dom_k) ) then
          dom_j = dom_k
       else
          dom_j = 0
       end if
       !
       ! Synchronous or asynchronous
       !
       if( present(wsynch) ) then
          if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
             asynch = .false.
          else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
             asynch = .true.
          else
             call runend('PAR_ELEMENT_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
          end if
       else
          asynch = .false.
       end if

       if( ISLAVE ) then

          if( ipass == 1 ) then
             !
             ! Allocate memory
             !
             if( commu % ghost_send_elem_dim == -1 .or. commu % ghost_recv_elem_dim == -1 ) &
                  call runend('FRINGE NODE EXCHANGE NOT COMPUTED')
             if( asynch ) allocate(ireq4(commu % nneig*2))
             allocate(tmp_isend(commu % ghost_recv_elem_dim * ndofn))
             allocate(tmp_irecv(commu % ghost_send_elem_dim * ndofn))
             !
             ! Save in temp_send
             !
             kk = 0
             do jj = 1,commu % ghost_recv_elem_dim
                ipoin = commu % ghost_recv_elem_perm(jj)
                do idofn = 1,ndofn
                   kk = kk + 1
                   tmp_isend(kk) = xx(idofn,ipoin)
                end do
             end do
             !
             ! Send    temp_send
             ! Receive temp_recv
             !
             kk = 0
             istat4 = 0_4
             do ii = 1,commu % nneig

                dom_i  = commu % neights(ii)
                dom_i4 = int(dom_i,4)

                if( dom_j == 0 .or. dom_j == dom_i ) then

                   dom_i      = commd % neights(ii)
                   ini_send   = ndofn * ( commd % ghost_recv_elem_size(ii)   -1 ) + 1
                   nsize_send = ndofn * ( commd % ghost_recv_elem_size(ii+1) -1 ) + 1 - ini_send
                   ini_recv   = ndofn * ( commd % ghost_send_elem_size(ii)   -1 ) + 1
                   nsize_recv = ndofn * ( commd % ghost_send_elem_size(ii+1) -1 ) + 1 - ini_recv

                   nsize_send4 = int(nsize_send,4)
                   nsize_recv4 = int(nsize_recv,4)

                   if( asynch ) then
                      kk = kk + 1
                      call MPI_Isend(&
                           tmp_isend(ini_send:ini_send+nsize_send-1), nsize_send4, &
                           PAR_INTEGER,  dom_i4, 0_4,                     &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                      kk = kk + 1
                      call MPI_Irecv(&
                           tmp_irecv(ini_recv:ini_recv+nsize_recv-1), nsize_recv4, &
                           PAR_INTEGER,  dom_i4, 0_4,                     &
                           PAR_COMM_TO_USE, ireq4(kk), istat4                      )
                   else
                      if( nsize_recv /= 0 .and. nsize_send == 0 ) then
                         call MPI_Recv(                          &
                              tmp_irecv(ini_recv:), nsize_recv4, &
                              PAR_INTEGER, dom_i4, 0_4, &
                              PAR_COMM_TO_USE, status, istat4    )
                      else if( nsize_recv == 0 .and. nsize_send /= 0 ) then
                         call MPI_Send(                          &
                              tmp_isend(ini_send:), nsize_send4, &
                              PAR_INTEGER, dom_i4, 0_4, &
                              PAR_COMM_TO_USE, istat4            )
                      else if( nsize_recv /= 0 .and. nsize_send /= 0 ) then
                         call MPI_Sendrecv(                         &
                              tmp_isend(ini_send:), nsize_send4,    &
                              PAR_INTEGER, dom_i4, 0_4,    &
                              tmp_irecv(ini_recv:), nsize_recv4,    &
                              PAR_INTEGER, dom_i4, 0_4,    &
                              PAR_COMM_TO_USE, status, istat4       )
                      end if
                   end if
                   if( istat4 /= 0_4 ) call runend('PAR_GHOST_ELEMENT_EXCHANGE_IP: MPI ERROR')

                end if

             end do

          end if
          !
          ! sum,max,min on temp_recv
          !
          if( asynch .and. ipass == 2 ) then
             count4 = 2*int(commu % nneig,4)
             allocate( status4(MPI_STATUS_SIZE,2*commu % nneig) )
             CALL MPI_WAITALL(count4,ireq4,status4,istat4)
             if( istat4 /= 0 ) call runend('WRONG SEND/RECEIVE')
             deallocate( status4 )
             deallocate(ireq4)
          end if

          if( ( asynch .and. ipass == 2 ) .or. ( .not. asynch .and. ipass == 1 ) ) then

             if( trim(what) == 'SUM' .or. trim(what) == 'ASSEMBLY' ) then
                !
                ! SUM
                !
                kk = 0
                do jj = 1,commu % ghost_send_elem_dim
                   ipoin = commu % ghost_send_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = xx(idofn,ipoin) + tmp_irecv(kk)
                   end do
                end do

             else if( trim(what) == 'MAX' .or. trim(what) == 'MAXIMUM' ) then
                !
                ! MAX
                !
                kk = 0
                do jj = 1,commu % ghost_send_elem_dim
                   ipoin = commu % ghost_send_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = max(xx(idofn,ipoin),tmp_irecv(kk))
                   end do
                end do

             else if( trim(what) == 'MIN' .or. trim(what) == 'MINIMUM' ) then
                !
                ! MIN
                !
                kk = 0
                do jj = 1,commu % ghost_send_elem_dim
                   ipoin = commu % ghost_send_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = min(xx(idofn,ipoin),tmp_irecv(kk))
                   end do
                end do

             else if( trim(what) == 'REPLACE' .or. trim(what) == 'SUBSTITUTE' ) then
                !
                ! Replace value on my fringe nodes according to what I have received
                !
                kk = 0
                do jj = 1,commu % ghost_send_elem_dim
                   ipoin = commu % ghost_send_elem_perm(jj)
                   do idofn = 1,ndofn
                      kk = kk + 1
                      xx(idofn,ipoin) = tmp_irecv(kk)
                   end do
                end do

             else
                call runend('UNKNOWN ORDER')
             end if

             ipass = 0
             deallocate(tmp_irecv)
             deallocate(tmp_isend)

          end if

       end if

    end if
#endif

  end subroutine PAR_FROM_GHOST_ELEMENT_EXCHANGE_IP

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-06-04
  !> @brief   All to all
  !> @details MPI_ALLTOALL
  !>
  !-----------------------------------------------------------------------

!!$  subroutine PAR_ALLATOALL_IP_1(xx_send,xx_recv,wherein,wsynch,PAR_COMM_IN4)
!!$
!!$    integer(ip),  pointer,  intent(in)  :: xx_send(:)
!!$    integer(ip),  pointer,  intent(inout) :: xx_recv(:)
!!$    character(*), optional, intent(in)  :: wherein
!!$    character(*), optional, intent(in)  :: wsynch
!!$    integer(4),   optional, intent(in)  :: PAR_COMM_IN4
!!$    integer(ip)                         :: nsend,nrecv
!!$    integer(ip)                         :: yy_send(2)
!!$    integer(ip)                         :: yy_recv(2)
!!$
!!$    nsend = memor_size(xx_send)
!!$    nrecv = memor_size(xx_recv)
!!$
!!$    if(      nsend == 0 .and. nrecv /= 0 ) then
!!$       call PAR_ALLATOALL_IP(nsend,nrecv,yy_send,xx_recv,wherein,wsynch,PAR_COMM_IN4)
!!$    else if( nsend /= 0 .and. nrecv == 0 ) then
!!$       call PAR_ALLATOALL_IP(nsend,nrecv,xx_send,yy_recv,wherein,wsynch,PAR_COMM_IN4)
!!$    else if( nsend == 0 .and. nrecv == 0 ) then
!!$       return
!!$    else
!!$       call PAR_ALLATOALL_IP(nsend,nrecv,xx_send,xx_recv,wherein,wsynch,PAR_COMM_IN4)
!!$    end if
!!$
!!$  end subroutine PAR_ALLATOALL_IP_1

  subroutine PAR_ALLTOALL_IP_0(nsend,nrecv,xx_send,xx_recv,wherein,wsynch,PAR_COMM_IN4)

    integer(ip),            intent(in)  :: nsend
    integer(ip),            intent(in)  :: nrecv
    integer(ip),            intent(in)  :: xx_send(*)
    integer(ip),            intent(out) :: xx_recv(*)
    character(*), optional, intent(in)  :: wherein
    character(*), optional, intent(in)  :: wsynch
    integer(4),   optional, intent(in)  :: PAR_COMM_IN4
    integer(4)                          :: istat4,nsend4,nrecv4
    integer(4)                          :: PAR_COMM_TO_USE
    logical(lg)                         :: asynch
    !
    ! Define communicator
    !
    if( present(PAR_COMM_IN4) ) then
       PAR_COMM_TO_USE = PAR_COMM_IN4
    else if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
    else
       PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
    end if
    !
    ! BlockinG/non blocking
    !
    if( present(wsynch) ) then
       if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
          asynch = .false.
       else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
          asynch = .true.
          if( .not. associated(non_blocking(inonblocking) % request4) ) then
             call runend('NON-BLOCKING SEND/RECEIVE SHOULD BE STARTED')
          end if
       else
          call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
       end if
    else
       asynch = .false.
    end if
    if( asynch ) call runend('ALL TO ALL NOT CODED')
    !
    ! Synchronous Send/receive
    !
    nsend4 = int(nsend,4)
    nrecv4 = int(nrecv,4)
    istat4 = 0

#ifndef MPI_OFF
    call MPI_ALLTOALL(xx_send,nsend4,PAR_INTEGER,xx_recv,nrecv4,PAR_INTEGER,PAR_COMM_TO_USE,istat4)
#endif
    if( istat4 /= 0_4 ) call runend('PAR_SEND_RECEIVE_IP: MPI ERROR')

  end subroutine PAR_ALLTOALL_IP_0

  subroutine PAR_ALLTOALL_IP_1(xx_send,xx_recv,wherein,wsynch,PAR_COMM_IN4)

    integer(ip),            intent(in),    pointer :: xx_send(:)
    integer(ip),            intent(inout), pointer :: xx_recv(:)
    character(*), optional, intent(in)             :: wherein
    character(*), optional, intent(in)             :: wsynch
    integer(4),   optional, intent(in)             :: PAR_COMM_IN4
    integer(ip)                                    :: nsend,nrecv
    integer(ip)                                    :: lboun_send,lboun_recv
    integer(4)                                     :: istat4,nsend4,nrecv4
    integer(4)                                     :: PAR_COMM_TO_USE
    integer(4)                                     :: my_rank4,my_size4
    logical(lg)                                    :: asynch
    !
    ! Define communicator
    !
    if( present(PAR_COMM_IN4) ) then
       PAR_COMM_TO_USE = PAR_COMM_IN4
    else if( present(wherein) ) then
       call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
    else
       PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
    end if
    call PAR_COMM_RANK_AND_SIZE(PAR_COMM_TO_USE,my_rank4,my_size4)
    !
    ! BlockinG/non blocking
    !
    if( present(wsynch) ) then
       if( trim(wsynch) == 'SYNCHRONOUS' .or. trim(wsynch) == 'BLOCKING' ) then
          asynch = .false.
       else if( trim(wsynch) == 'ASYNCHRONOUS' .or. trim(wsynch) == 'NON BLOCKING' ) then
          asynch = .true.
          if( .not. associated(non_blocking(inonblocking) % request4) ) then
             call runend('NON-BLOCKING SEND/RECEIVE SHOULD BE STARTED')
          end if
       else
          call runend('PAR_NODE_ASSMEMBLY: UNKNOWN COMMUNICATION TYPE')
       end if
    else
       asynch = .false.
    end if
    if( asynch ) call runend('ALL TO ALL NOT CODED')
    !
    ! Synchronous Send/receive
    !
    nsend      = memory_size(xx_send)/int(my_size4,ip)
    nrecv      = memory_size(xx_recv)/int(my_size4,ip)
    lboun_send = lbound(xx_send,1_ip)
    lboun_recv = lbound(xx_recv,1_ip)
    nsend4     = int(nsend,4)
    nrecv4     = int(nrecv,4)
    istat4     = 0

#ifndef MPI_OFF
    call MPI_ALLTOALL(xx_send(lboun_send:),nsend4,PAR_INTEGER,xx_recv(lboun_recv:),nrecv4,PAR_INTEGER,PAR_COMM_TO_USE,istat4)
#endif
    if( istat4 /= 0_4 ) call runend('PAR_SEND_RECEIVE_IP: MPI ERROR')

  end subroutine PAR_ALLTOALL_IP_1

  !-----------------------------------------------------------------------
  !>
  !> @author  Ricard Borell
  !> @date    2018-12-04
  !> @brief   All to all v
  !> @details All to all v
  !>
  !-----------------------------------------------------------------------

  subroutine PAR_ALLTOALLV_IP_1(sendbuf,recvbuf,sendcount4,recvcount4,wherein,PAR_COMM_IN4)

    integer(ip),  pointer, intent(in)           :: sendbuf(:)           !< Send buffer
    integer(ip),  pointer, intent(inout)        :: recvbuf(:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)        !< Send counts
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in), optional :: wherein              !< Wherein
    integer(4),            intent(in), optional :: PAR_COMM_IN4         !< Communicator
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: ipart
    integer(4)                                  :: mpi_sumsend
    integer(4)                                  :: mpi_sumrecv
    integer(4),   pointer                       :: mpi_sdispls(:)
    integer(4),   pointer                       :: mpi_rdispls(:)
    integer(ip),  target                        :: send_null(2)
    integer(ip),  target                        :: recv_null(2)
    integer(ip),  pointer                       :: send_tmp(:)
    integer(ip),  pointer                       :: recv_tmp(:)

#ifndef MPI_OFF

    if( IPARALL ) then

       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
       end if
       call MPI_Comm_size(PAR_COMM_TO_USE,comm_size,istat4)

       allocate(mpi_sdispls(0:comm_size-1))
       allocate(mpi_rdispls(0:comm_size-1))

       mpi_sumsend = 0_4
       mpi_sumrecv = 0_4
       do ipart = 0,comm_size-1
          mpi_sdispls(ipart) = mpi_sumsend
          mpi_rdispls(ipart) = mpi_sumrecv
          mpi_sumsend        = mpi_sumsend + sendcount4(ipart)
          mpi_sumrecv        = mpi_sumrecv + recvcount4(ipart)
       end do

       if( mpi_sumsend == 0 ) then
          send_tmp => send_null
       else
          send_tmp => sendbuf
       end if
       if( mpi_sumrecv == 0 ) then
          recv_tmp => recv_null
       else
          recv_tmp => recvbuf
       end if

       call MPI_Alltoallv(&
            send_tmp,sendcount4,mpi_sdispls, PAR_INTEGER,  &
            recv_tmp,recvcount4,mpi_rdispls, PAR_INTEGER,  &
            PAR_COMM_TO_USE,istat4)

       deallocate(mpi_sdispls)
       deallocate(mpi_rdispls)

    end if

#endif

  end subroutine PAR_ALLTOALLV_IP_1

  subroutine PAR_ALLTOALLV_IP_2(sendbuf,recvbuf,sendcount4,recvcount4,wherein,PAR_COMM_IN4)

    integer(ip),  pointer, intent(in)           :: sendbuf(:,:)         !< Send buffer
    integer(ip),  pointer, intent(inout)        :: recvbuf(:,:)         !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)        !< Send counts
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in), optional :: wherein              !< Wherein
    integer(4),            intent(in), optional :: PAR_COMM_IN4         !< Communicator
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: ipart
    integer(4)                                  :: mpi_sumsend
    integer(4)                                  :: mpi_sumrecv
    integer(4),   pointer                       :: mpi_sdispls(:)
    integer(4),   pointer                       :: mpi_rdispls(:)
    integer(ip),  target                        :: send_null(2,2)
    integer(ip),  target                        :: recv_null(2,2)
    integer(ip),  pointer                       :: send_tmp(:,:)
    integer(ip),  pointer                       :: recv_tmp(:,:)

#ifndef MPI_OFF

    if( IPARALL ) then

       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
       end if
       call MPI_Comm_size(PAR_COMM_TO_USE,comm_size,istat4)

       allocate(mpi_sdispls(0:comm_size-1))
       allocate(mpi_rdispls(0:comm_size-1))

       mpi_sumsend = 0_4
       mpi_sumrecv = 0_4
       do ipart = 0,comm_size-1
          mpi_sdispls(ipart) = mpi_sumsend
          mpi_rdispls(ipart) = mpi_sumrecv
          mpi_sumsend        = mpi_sumsend + sendcount4(ipart)
          mpi_sumrecv        = mpi_sumrecv + recvcount4(ipart)
       end do

       if( mpi_sumsend == 0 ) then
          send_tmp => send_null
       else
          send_tmp => sendbuf
       end if
       if( mpi_sumrecv == 0 ) then
          recv_tmp => recv_null
       else
          recv_tmp => recvbuf
       end if

       call MPI_Alltoallv(&
            send_tmp,sendcount4,mpi_sdispls, PAR_INTEGER,  &
            recv_tmp,recvcount4,mpi_rdispls, PAR_INTEGER,  &
            PAR_COMM_TO_USE,istat4)

       deallocate(mpi_sdispls)
       deallocate(mpi_rdispls)

    end if

#endif

  end subroutine PAR_ALLTOALLV_IP_2

  subroutine PAR_ALLTOALLV_RP_1(sendbuf,recvbuf,sendcount4,recvcount4,wherein,PAR_COMM_IN4)

    real(rp),     pointer, intent(in)           :: sendbuf(:)           !< Send buffer
    real(rp),     pointer, intent(inout)        :: recvbuf(:)           !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)        !< Send counts
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in), optional :: wherein              !< Wherein
    integer(4),            intent(in), optional :: PAR_COMM_IN4         !< Communicator
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: ipart
    integer(4)                                  :: mpi_sumsend
    integer(4)                                  :: mpi_sumrecv
    integer(4),   pointer                       :: mpi_sdispls(:)
    integer(4),   pointer                       :: mpi_rdispls(:)
    real(rp),     target                        :: send_null(2)
    real(rp),     target                        :: recv_null(2)
    real(rp),     pointer                       :: send_tmp(:)
    real(rp),     pointer                       :: recv_tmp(:)

#ifndef MPI_OFF

    if( IPARALL ) then

       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
       end if
       call MPI_Comm_size(PAR_COMM_TO_USE,comm_size,istat4)

       allocate(mpi_sdispls(0:comm_size-1))
       allocate(mpi_rdispls(0:comm_size-1))

       mpi_sumsend = 0_4
       mpi_sumrecv = 0_4
       do ipart = 0,comm_size-1
          mpi_sdispls(ipart) = mpi_sumsend
          mpi_rdispls(ipart) = mpi_sumrecv
          mpi_sumsend        = mpi_sumsend + sendcount4(ipart)
          mpi_sumrecv        = mpi_sumrecv + recvcount4(ipart)
       end do

       if( mpi_sumsend == 0 ) then
          send_tmp => send_null
       else
          send_tmp => sendbuf
       end if
       if( mpi_sumrecv == 0 ) then
          recv_tmp => recv_null
       else
          recv_tmp => recvbuf
       end if

       call MPI_Alltoallv(&
            send_tmp,sendcount4,mpi_sdispls, MPI_DOUBLE_PRECISION,  &
            recv_tmp,recvcount4,mpi_rdispls, MPI_DOUBLE_PRECISION,  &
            PAR_COMM_TO_USE,istat4)

       deallocate(mpi_sdispls)
       deallocate(mpi_rdispls)

    end if

#endif

  end subroutine PAR_ALLTOALLV_RP_1

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2019-01-02
  !> @brief   Transform an MPI error code into a string
  !> @details Transform an MPI error code into a string
  !>
  !-----------------------------------------------------------------------

  subroutine PAR_MPI_ERROR_TO_MESSAGE(ierror,length,message)

    integer(ip),      intent(out)       :: length    !< Length of the message
    integer(4),       intent(out)       :: ierror    !< Error code
    character(len=*), intent(out)       :: message   !< Message
    integer(4)                          :: length4,temp4,ierror4

#ifndef MPI_OFF
    character(len=MPI_MAX_ERROR_STRING) :: message_mpi

    if( ierror4 /= MPI_SUCCESS ) then
       call MPI_Error_string(ierror4,message_mpi,length4,temp4)
       length4 = min(len(message),length4)
       message = message_mpi(1:length4)
    end if
#endif

  end subroutine PAR_MPI_ERROR_TO_MESSAGE

  subroutine PAR_ALLTOALLV_RP_2(sendbuf,recvbuf,sendcount4,recvcount4,wherein,PAR_COMM_IN4)

    real(rp),     pointer, intent(in)           :: sendbuf(:,:)         !< Send buffer
    real(rp),     pointer, intent(inout)        :: recvbuf(:,:)         !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)        !< Send counts
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in), optional :: wherein              !< Wherein
    integer(4),            intent(in), optional :: PAR_COMM_IN4         !< Communicator
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: ipart
    integer(4)                                  :: mpi_sumsend
    integer(4)                                  :: mpi_sumrecv
    integer(4),   pointer                       :: mpi_sdispls(:)
    integer(4),   pointer                       :: mpi_rdispls(:)
    real(rp),     target                        :: send_null(2,2)
    real(rp),     target                        :: recv_null(2,2)
    real(rp),     pointer                       :: send_tmp(:,:)
    real(rp),     pointer                       :: recv_tmp(:,:)

#ifndef MPI_OFF

    if( IPARALL ) then

       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
       end if
       call MPI_Comm_size(PAR_COMM_TO_USE,comm_size,istat4)

       allocate(mpi_sdispls(0:comm_size-1))
       allocate(mpi_rdispls(0:comm_size-1))

       mpi_sumsend = 0_4
       mpi_sumrecv = 0_4
       do ipart = 0,comm_size-1
          mpi_sdispls(ipart) = mpi_sumsend
          mpi_rdispls(ipart) = mpi_sumrecv
          mpi_sumsend        = mpi_sumsend + sendcount4(ipart)
          mpi_sumrecv        = mpi_sumrecv + recvcount4(ipart)
       end do

       if( mpi_sumsend == 0 ) then
          send_tmp => send_null
       else
          send_tmp => sendbuf
       end if
       if( mpi_sumrecv == 0 ) then
          recv_tmp => recv_null
       else
          recv_tmp => recvbuf
       end if

       call MPI_Alltoallv(&
            send_tmp,sendcount4,mpi_sdispls, MPI_DOUBLE_PRECISION,  &
            recv_tmp,recvcount4,mpi_rdispls, MPI_DOUBLE_PRECISION,  &
            PAR_COMM_TO_USE,istat4)

       deallocate(mpi_sdispls)
       deallocate(mpi_rdispls)

    end if

#endif

  end subroutine PAR_ALLTOALLV_RP_2

  subroutine PAR_ALLTOALLV_RP_3(sendbuf,recvbuf,sendcount4,recvcount4,wherein,PAR_COMM_IN4)

    real(rp),     pointer, intent(in)           :: sendbuf(:,:,:)       !< Send buffer
    real(rp),     pointer, intent(inout)        :: recvbuf(:,:,:)       !< Recv buffer
    integer(4),   pointer, intent(in)           :: sendcount4(:)        !< Send counts
    integer(4),   pointer, intent(in)           :: recvcount4(:)        !< Recv counts
    character(*),          intent(in), optional :: wherein              !< Wherein
    integer(4),            intent(in), optional :: PAR_COMM_IN4         !< Communicator
    integer(4)                                  :: istat4
    integer(4)                                  :: comm_size
    integer(4)                                  :: PAR_COMM_TO_USE
    integer(4)                                  :: ipart
    integer(4)                                  :: mpi_sumsend
    integer(4)                                  :: mpi_sumrecv
    integer(4),   pointer                       :: mpi_sdispls(:)
    integer(4),   pointer                       :: mpi_rdispls(:)
    real(rp),     target                        :: send_null(2,2,2)
    real(rp),     target                        :: recv_null(2,2,2)
    real(rp),     pointer                       :: send_tmp(:,:,:)
    real(rp),     pointer                       :: recv_tmp(:,:,:)

#ifndef MPI_OFF

    if( IPARALL ) then

       if( present(PAR_COMM_IN4) ) then
          PAR_COMM_TO_USE = PAR_COMM_IN4
       else if( present(wherein) ) then
          call PAR_DEFINE_COMMUNICATOR(wherein,PAR_COMM_TO_USE)
       else
          PAR_COMM_TO_USE = int(PAR_COMM_MY_CODE,ip)
       end if
       call MPI_Comm_size(PAR_COMM_TO_USE,comm_size,istat4)

       allocate(mpi_sdispls(0:comm_size-1))
       allocate(mpi_rdispls(0:comm_size-1))

       mpi_sumsend = 0_4
       mpi_sumrecv = 0_4
       do ipart = 0,comm_size-1
          mpi_sdispls(ipart) = mpi_sumsend
          mpi_rdispls(ipart) = mpi_sumrecv
          mpi_sumsend        = mpi_sumsend + sendcount4(ipart)
          mpi_sumrecv        = mpi_sumrecv + recvcount4(ipart)
       end do

       if( mpi_sumsend == 0 ) then
          send_tmp => send_null
       else
          send_tmp => sendbuf
       end if
       if( mpi_sumrecv == 0 ) then
          recv_tmp => recv_null
       else
          recv_tmp => recvbuf
       end if

       call MPI_Alltoallv(&
            send_tmp,sendcount4,mpi_sdispls, MPI_DOUBLE_PRECISION,  &
            recv_tmp,recvcount4,mpi_rdispls, MPI_DOUBLE_PRECISION,  &
            PAR_COMM_TO_USE,istat4)

       deallocate(mpi_sdispls)
       deallocate(mpi_rdispls)

    end if

#endif

  end subroutine PAR_ALLTOALLV_RP_3

end module mod_communications
!> @}


