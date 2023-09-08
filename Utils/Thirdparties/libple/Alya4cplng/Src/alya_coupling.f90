module alya_coupling
  use mpi
  use iso_c_binding 
  implicit none

  integer, parameter :: lchar = 16 

!
!typedef struct 
!{
!  int          status;    /* Status flag for synchronization info */
!  int          root_rank; /* Application root rank in MPI_COMM_WORLD */
!  int          n_ranks;   /* Number of ranks associated with application */
!  const char  *app_type;  /* Application type name (may be empty) */
!  const char  *app_name;  /* Application instance name (may be empty) */
!
!} ple_coupling_mpi_set_info_t;
!
  type ple_coupling_mpi_set_info_t
    sequence 
    integer                        :: status
    integer                        :: root_rank 
    integer                        :: n_ranks
    character(len= lchar), pointer :: app_type
    character(len= lchar), pointer :: app_name
  end type

!
!struct _ple_coupling_mpi_set_t
!{
!  int       n_apps;       /* Number of distinct applications 
!  int       app_id;       /* Id of the local application in the application info 
!  int       app_names_l;  /* Length of application names array 
!  int      *app_info;     /* For each application, 4 integers: 
!                             associated root in base_comm, n_ranks, and indexes in app_names 
!  char     *app_names;    /* Array of application type names and instance names 
!  int      *app_status;   /* Synchronization status for each application 
!  double   *app_timestep; /* Current time step for each application 
!  MPI_Comm  base_comm;    /* Handle to base communicator 
!  MPI_Comm  app_comm;     /* Handle to local application communicator 
!};
!
  type alya_coupling_mpi_set_t
    sequence 
    integer                          :: n_apps
    integer                          :: app_id
    integer                          :: app_names_l
    integer,   pointer               :: app_info
   !>integer,   pointer               :: app_info(:)
   !>character(len=16*lchar)          :: app_names
    character(len=16*lchar), pointer :: app_names
   !>character(len=16*lchar, kind=c_char), pointer        :: app_names
    integer,   pointer               :: app_status
    real,      pointer               :: app_timestep
   !>integer,   pointer               :: app_status(:);   
   !>real,      pointer               :: app_timestep(:)
    integer                          :: base_comm;    
    integer                          :: app_comm;    
  end type


contains 
  subroutine alya_coupling_mpi_set_create2(app_type, & 
                                           app_name, &
                                           base_comm, &
                                           app_comm, & 
                                           s)
   !integer                       ::  sync_flag 
    character(len= 2*lchar)       ::  app_type
    character(len= 2*lchar)       ::  app_name
    character(len=16*lchar)       ::  app_names, sendbuf 
    integer                       ::  base_comm
    integer                       ::  app_comm
    integer                       ::  app_rank
    integer                       ::  n_app_ranks
    integer                       ::  set_id
    integer                       ::  error, i, j, k, start, app_ok 
    integer                       ::  root_marker
    integer                       ::  info_tag 
    integer                       ::  l_rank_info(5) = (/-7, -11, -13, -17, -19/)
    integer                       ::  counter(2)     = (/0, 0/) 
    integer, pointer              ::  rank_info(:)
    integer                       ::  app_type_size, app_name_size, msg_len, name_tag, sendbuf_l
    integer                       ::  recvbuf_tag, recvbuf_l
    type(alya_coupling_mpi_set_t) :: s

    info_tag    =  1
    name_tag    =  2
    recvbuf_tag =  3
    app_rank    = -1
    set_id      = -1
    n_app_ranks =  0 
    root_marker =  0
    sendbuf_l   =  0 
    sendbuf     = ""
    app_names   = ""

    !> Initialization
    s%n_apps      =  0
    s%app_id      = -1
    s%app_names_l =  0
    s%base_comm   =  base_comm
    s%app_comm    =  app_comm

    call MPI_Comm_rank(base_comm, set_id, error)

    if(app_comm /= 0) then
      call MPI_Comm_rank(app_comm,   app_rank, error)
      call MPI_Comm_size(app_comm, n_app_ranks, error)
    else
      app_rank    = 0
      n_app_ranks = 1
    endif 

    l_rank_info(2) = set_id
    l_rank_info(3) = n_app_ranks
    if(len_trim(app_type)>0) l_rank_info(4) = len_trim(app_type)
    
    l_rank_info(5) = 0 
    if(len_trim(app_name)>1) l_rank_info(5) = len_trim(app_name) 

    !> Root rank of base_comm counts applications and receives info 
    if(app_rank == 0) root_marker = 1
    call MPI_Reduce(root_marker, counter(1), 1, MPI_INTEGER, MPI_SUM, 0, base_comm, error)
    call MPI_Bcast(counter, 2, MPI_INTEGER, 0, base_comm, error)

    allocate( rank_info(counter(1)*5) ) 
    rank_info = -1 
    
    !> Root of base_comm collects all info
    if(set_id == 0) then 
      app_ok = 0 
      if(app_rank==0) app_ok = 1

      start = 1
      if(app_ok == 1) start = 2
      !> Root of base_comm collects all info 
      !> Use of different tags for info and strings is important
      !> here as we use MPI_ANY_SOURCE and messages could be mixed
      do i = start,counter(1)
        call MPI_Recv(rank_info((i-1)*5+1:), 5, MPI_INTEGER, MPI_ANY_SOURCE, info_tag, base_comm, error) 
      enddo 
      rank_info(1:5) = l_rank_info(1:5)

      !> Convert rank_info count to index values
      do i = 1,counter(1)
        counter(2) = counter(2) + rank_info((i-1)*5 + 4) + rank_info((i-1)*5 + 5)
      enddo

    !> Other root ranks send info to root
    else if(app_rank == 0) then 
      call MPI_Send(l_rank_info, 5, MPI_INTEGER, 0, info_tag, base_comm, error)
    endif 

    call MPI_Bcast(counter, 2, MPI_INTEGER, 0, base_comm, error)

    !> Now root broadcasts application info
    call MPI_Bcast(rank_info, 5*counter(1),   MPI_INTEGER, 0, base_comm, error) 
    call MPI_Bcast(app_names,   counter(2), MPI_CHARACTER, 0, base_comm, error)

    !> Set global values
    s%n_apps      = counter(1)
    s%app_names_l = counter(2)
    !s%app_names   = app_names
!    allocate( s%app_info(     counter(1)*4) )
!    allocate( s%app_status(   counter(1)  ) )
!    allocate( s%app_timestep( counter(1)  ) )

    do i = 1,s%n_apps 
      if(s%app_id < 0) then 
        do j = 1,4
        !  s%app_info((i-1)*4 + j) = rank_info((i-1)*5 + j+1)
        enddo
        !s%app_status(i) = rank_info((i-1)*5 +1)
        !s%app_timestep(i) = 0.0
      endif 
    enddo 

    !> Set rank set to that of the application root for matching
    call MPI_Bcast(set_id, 1, MPI_INTEGER, 0, app_comm, error)

    do i = 1,s%n_apps
      if(s%app_id < 0) then 
!        if(s%app_info((i-1)*4+1) == set_id) s%app_id = i-1
      endif 
    enddo 
  end subroutine

end module
