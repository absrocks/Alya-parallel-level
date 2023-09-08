module def_parall

  !-----------------------------------------------------------------------
  !****f* Parall/def_parall
  ! NAME
  !    def_parall
  ! DESCRIPTION
  !    Heading for the PARALL service
  ! USED BY
  !    Almost all parall subroutines
  !***
  !-----------------------------------------------------------------------
  
  use def_kintyp
  use def_domain, only : nelty
  
  !------------------------------------------------------------------------
  ! Types
  !------------------------------------------------------------------------

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
  integer(4),  pointer      :: PAR_COMM_ZONES(:)
  type(comm_data_par), pointer :: commc      ! Coarse grid communication arrays
  type(comm_data_par), pointer :: commz      ! Zone-wise communication arrays
  type(comm_data_par), pointer :: commd_glo  ! Global communication arrays

  !------------------------------------------------------------------------
  ! Parameters
  !------------------------------------------------------------------------

  integer(ip), parameter   :: &
       lun_domai_par = 5503,  &    ! Graph partition
       lun_trace_par = 5504,  &    ! Trace
       lun_conve_par = 5506,  &    ! Time statistics
       lun_rstar_par = 5516,  &    ! Restart
       lun_aonlp_par = 10000       ! additional number for ONLY_PREPROCESS
  integer(ip),   parameter :: &
       nvarp_par=20                ! # postprocess variables

  !------------------------------------------------------------------------
  ! File names
  !------------------------------------------------------------------------

  character(150)           :: &
       fil_rstar_par               ! Restart file

  !------------------------------------------------------------------------
  ! Reapro
  !------------------------------------------------------------------------

  integer(ip)                            :: &
       npart_par,                           &    ! Number of subdomains
       kfl_ascii_par,                       &    ! Restart SCII format(=1)
       kfl_bytes_par,                       &    ! Integer bytes for files
       kfl_parti_par,                       &    ! Partition type
       kfl_fileh_par,                       &    ! File hierarchy file
       kfl_filio_par,                       &    ! Open and close files in preprocess
       kfl_weigh_par,                       &    ! Weighting of Metis graph
       kfl_virfi_par,                       &    ! Virtual files
       kfl_global_numbering_par,            &    ! Global numbering strategy
       nsire_par,                           &    ! Number of simultaneous reading (restart)
       lzone_par(10),                       &    ! Number of METIS zones
       kfl_matri_par,                       &    ! Global matrix postprocess
       kfl_partition_par,                   &    ! Partition method
       kfl_interface_parti_par,             &    ! Interface partitioning strategy
       kfl_parseq_par,                      &    ! Parallel or sequential partitioning
       kfl_cores_per_gpu,                   &    ! Cores per GPU
       kfl_streams_per_gpu,                 &    ! Streams per GPU
       boxes_coarse_par(3),                 &    ! Number of boxes coarse bin
       boxes_fine_par(3),                   &    ! Number of boxes fine bin
       sfc_criteria,                        &
       sfc_check,                           &
       sfc_dim_bin_core                      
 
  real(rp)                               :: &
       rmbyt_par,                           &    ! Max number of Gb for virtual files
       vect_partition_par(3),               &    ! Direction of partition for oriented bin
       mpio_val_hybrid_threshold,           &    ! Minimum size (in MB) per process to enable parallel IO (disabled if 0)
       weights_elements_par(nelty),         &    ! Weights for element types
       weights_materials_par(10_ip)              ! Weights for materials 

  character(20)                          :: &
       method_redistribution_par                 ! Redistribution strategy
  
  !------------------------------------------------------------------------
  ! Others
  !------------------------------------------------------------------------

  !
  ! Partition scalar data
  !
  integer(ip)              :: &
       nproc_par,             &    ! Number of processes
       iproc_par,             &    ! My process
       nbcol,                 &    ! Number of colours
       nneig,                 &    ! Local number of neighbour domains
       slfbo,                 &    ! Init of self boundary
       iproc_part_par              ! My order into the partition slaves

  real(rp) :: &
       cpu_paral(50)

  character(5)             :: &
       wopos_par(2,nvarp_par)      ! Name and character of the postprocess variables
  !
  ! Mesh graph
  !
  integer(ip), pointer     :: &
       padja_par(:),          &
       ladja_par(:)

  integer(ip), pointer     :: &
       lepar_par(:),          & ! Domain of every element
       lnpar_par(:),          & ! Domain of every node
       lbpar_par(:),          & ! Domain of every boundary element
       leper_par(:),          &
       lbper_par(:),          & ! Boundary permutation
       lneig_par(:),          & ! number of neightbours of every domain
       ginde_par(:,:),        &
       lcomm_par(:,:),        &
       nskew_par(:),          &
       ngive_par(:),          &
       slfbo_par(:),          &
       leind_par(:),          &
       lbind_par(:),          &
       lnods_par(:,:),        &
       nhang_par(:),          &
       lsubz_par(:),          & ! Subdomain zones
       nelew_par(:)             ! Number of weights

  integer(ip)              :: &
       gnbop_loc

  integer(ip), pointer     :: &
       xlnin_loc(:),          &
       exnpe_loc(:)

  integer(ip), pointer     :: &
       neighDom(:),           &
       xadjDom(:),            &
       adjDom(:),             &
       translDual(:),         &
       iaDual(:),             &
       jaDual(:),             &
       colours(:)

  integer(ip), pointer     :: &
       badj(:),               &
       bdom(:),               &
       bpoin(:)

  integer(ip), pointer     :: &
       leinv_par(:),          &  ! inverse perm of elements
       lbinv_par(:),          &  ! inverse perm of boundaries
       lnper_par(:),          &  ! perm of nodes
       lninv_par(:)              ! inverse perm of nodes (1..npoin)
  !
  ! Asynchronous communications
  !
  integer(4),  pointer     :: &
       ireq4(:)                  ! Request for asynchronous communication
  integer(ip)              :: &
       ipass_par                 ! Request for asynchronous communication
  real(rp),    pointer     :: &
       parws(:)                  ! Send array
  real(rp),    pointer     :: &
       parwr(:)                  ! Receive array
  !
  ! mesh multiplication
  !
  integer(ip), pointer     :: &
       lowns_par(:),          &  ! List of my own nodes
       lownr_par(:)              ! List of my neighbor's own nodes
  !
  ! Sets
  !
  integer(ip), pointer     :: &
       lnsec_par(:,:)            ! Node set
  integer(ip), pointer     :: &
       nnset_par(:),          &  ! Node set per subdomain
       nwitn_par(:)              ! Witness points per subdomain
  !
  ! Fringe geometry communication arrays
  !
  type fringe_comm
     integer(ip)               :: nbcos
     integer(ip)               :: nbcor
     integer(ip), pointer      :: lbcos(:)
     integer(ip), pointer      :: lbcor(:)
  end type fringe_comm
  type(fringe_comm), pointer :: frcom(:)

  type(comm_data_par)       :: gro_commd
  type(comm_data_level_par) :: comle(10)

end module def_parall
