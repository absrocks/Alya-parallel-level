!-----------------------------------------------------------------------
!> @addtogroup Kinds_and_types
!> @{
!> @file    def_kintyp_comm.g90
!> @author  houzeaux
!> @date    2020-04-04
!> @brief   Functions
!> @details Communications
!-----------------------------------------------------------------------

module def_kintyp_comm

  use def_kintyp_basic

  type comm_bound_matrix
     type(i1p),   pointer :: ja(:)
     integer(ip), pointer :: nzdom_ii(:)
     integer(ip)          :: nzdom
  end type comm_bound_matrix
  !
  ! Automatic partitioning based on parallel efficiency
  !
  type typ_optimum_npart
     integer(ip) :: kfl_method
     integer(ip) :: frequency
     integer(ip) :: modul
     integer(ip) :: min_cores
     integer(ip) :: max_cores
     real(rp)    :: min_criterion     
     real(rp)    :: max_criterion
     real(rp)    :: change_rate
  end type typ_optimum_npart
  
  type comm_data_par_basic
     ! Neighbors
     integer(ip)                      :: nneig           ! Number of neigbors
     integer(ip), pointer             :: neights(:)      ! List of neighbors
     ! Interface node communication
     integer(ip)                      :: bound_dim       ! Size of interface
     integer(ip), pointer             :: bound_size(:)   ! Size of interface node esxchange
     integer(ip), pointer             :: bound_perm(:)   ! List of subdomain interface nodes
     ! Offsets
     integer(ip)                      :: offset_npoin    ! Offset of nodes
     integer(ip)                      :: offset_nelem    ! Offset of elements
     ! Communicator
     integer(ip)                      :: PAR_COMM_WORLD  ! Communicator
     integer(4)                       :: SIZE4           ! Size 
     integer(4)                       :: RANK4           ! Rank
  end type comm_data_par_basic
  
  type, extends(comm_data_par_basic)  :: comm_data_par
     ! Ordered neighbors
     integer(ip)                      :: nneig_1
     integer(ip)                      :: nneig_2
     integer(ip), pointer             :: neights_ordered(:)
     ! Interior, own and other boundary nodes
     integer(ip)                      :: npoi1
     integer(ip)                      :: npoi2
     integer(ip)                      :: npoi3
     integer(ip)                      :: npoin
     ! Interior, own and other boundary edges
     integer(ip)                      :: nedg1
     integer(ip)                      :: nedg2
     integer(ip)                      :: nedg3
     ! Offsets
     integer(ip)                      :: offset_nboun    ! Offset of boundaries     
     ! Interface node communication
     integer(ip), pointer             :: bound_invp(:)
     integer(ip), pointer             :: bound_scal(:)
     integer(ip), pointer             :: bound_multiplicity(:)
     integer(ip), pointer             :: bound_owner_rank(:)
     integer(ip), pointer             :: node_number_in_owner(:)
     integer(ip), pointer             :: perm_ordered(:)
     type(comm_bound_matrix), pointer :: bound_matrix(:)
     type(comm_bound_matrix), pointer :: bound_mat_halo_send(:)
     type(comm_bound_matrix), pointer :: bound_mat_halo_recv(:)
     ! Interface edge communication
     integer(ip), pointer             :: bedge_size(:)
     integer(ip), pointer             :: bedge_perm(:)
     integer(ip), pointer             :: bedge_adja(:)
     integer(ip), pointer             :: bedge_scal(:)
     integer(ip)                      :: bedge_dim
     integer(ip), pointer             :: bedge_multiplicity(:)
     integer(ip), pointer             :: bedge_owner_rank(:)
     ! Non-symmetric send receive
     integer(ip), pointer             :: lsend_size(:)
     integer(ip), pointer             :: lsend_perm(:)
     integer(ip)                      :: lsend_dim
     integer(ip), pointer             :: lrecv_size(:)
     integer(ip), pointer             :: lrecv_perm(:)
     integer(ip)                      :: lrecv_dim
     integer(ip), pointer             :: lscat_perm(:)
     integer(ip)                      :: lscat_dim
     integer(ip)                      :: matrix_nzdom
     integer(ip), pointer             :: matrix_ia(:)
     integer(ip), pointer             :: matrix_ja(:)
     real(rp),    pointer             :: matrix_aa(:)
     ! Face communication
     integer(ip), pointer             :: bface_size(:)
     integer(ip), pointer             :: bface_perm(:)
     integer(ip)                      :: bface_dim
     ! Full row matrix communicator
     integer(ip)                      :: full_row_send_nneig
     integer(ip), pointer             :: full_row_send_neights(:)
     integer(ip), pointer             :: full_row_send_size(:)
     integer(ip), pointer             :: full_row_send_perm(:)
     integer(ip)                      :: full_row_send_dim
     integer(ip)                      :: full_row_recv_nneig
     integer(ip), pointer             :: full_row_recv_neights(:)
     integer(ip), pointer             :: full_row_recv_size(:)
     integer(ip), pointer             :: full_row_recv_perm(:)
     integer(ip)                      :: full_row_recv_dim
     ! Ghost node communication
     integer(ip), pointer             :: ghost_send_node_size(:)
     integer(ip), pointer             :: ghost_send_node_perm(:)
     integer(ip)                      :: ghost_send_node_dim
     integer(ip), pointer             :: ghost_recv_node_size(:)
     integer(ip), pointer             :: ghost_recv_node_perm(:)
     integer(ip)                      :: ghost_recv_node_dim
     ! Ghost element communication
     integer(ip), pointer             :: ghost_send_elem_size(:)
     integer(ip), pointer             :: ghost_send_elem_perm(:)
     integer(ip)                      :: ghost_send_elem_dim
     integer(ip), pointer             :: ghost_recv_elem_size(:)
     integer(ip), pointer             :: ghost_recv_elem_perm(:)
     integer(ip)                      :: ghost_recv_elem_dim
     ! Ghost boundary communication
     integer(ip), pointer             :: ghost_send_boun_size(:)
     integer(ip), pointer             :: ghost_send_boun_perm(:)
     integer(ip)                      :: ghost_send_boun_dim
     integer(ip), pointer             :: ghost_recv_boun_size(:)
     integer(ip), pointer             :: ghost_recv_boun_perm(:)
     integer(ip)                      :: ghost_recv_boun_dim
  end type comm_data_par

   type :: tAdj_par
     integer(ip)            :: node1
     integer(ip)            :: node2
  end type tAdj_par
  type comm_data_level_par
     integer(ip), pointer   :: &
          neighDom(:),           &
          xadjDom(:),            &
          adjDom(:),             &
          translDual(:),         &
          iaDual(:),             &
          jaDual(:),             &
          colours(:),            &
          lnpar_par(:),          &
          lneig_par(:),          &
          lcomm_par(:,:),        &
          ngrou_par(:),          &
          lgrou(:),              &
          domli(:,:),            &
          ndomi(:),              &
          badj(:),               &
          bdom(:),               &
          bpoin(:),              &
                                ! Gather
          lbig(:),               &
          displ(:),              &
          lcoun(:)
     integer(4), pointer    ::   &
          disp4(:),              &
          lcou4(:)
     real(rp),   pointer    ::   &
          xsmall(:),             &
          xbig(:)
     integer(ip)            ::   &
          nbcol,                 &
          ngrou_total,           &
          ngrou,                 &
          gni,                   &
          gnb,                   &
          nneig,                 &
          nbig,                  &
          nsmall

     type(comm_data_par), pointer :: commd
  end type comm_data_level_par
  
end module def_kintyp_comm
!> @}
