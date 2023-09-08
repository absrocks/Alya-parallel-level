program main
  use mpi
  use alya_coupling
  implicit none 

  integer local_comm
  integer world_comm
  integer commij
  integer size_comm 

  integer app_id
  integer error
  integer world_rank, world_size
  integer local_rank, local_size
  integer rankij, sizeij

  integer local_range(2)
  integer rangei(2), rangej(2)
  integer appi, appj
  integer napps, ncplng, app_root 
  integer i, j

  character(len=lchar)               app_type, app_name
  type(alya_coupling_mpi_set_t)      cplng
  type(ple_coupling_mpi_set_info_t), pointer :: info(:)

  integer, pointer :: rangeij(:,:)
  real,    pointer :: timestep(:)  
  
  !> World
  world_comm = MPI_COMM_WORLD

  call MPI_Init(error)
  call MPI_Comm_size(world_comm, world_size, error)
  call MPI_Comm_rank(world_comm, world_rank, error)

  !> Split 
  app_id = -1
  call getarg(1, app_type)
  call alya_coupling_mpi_name_to_id(world_comm, app_type, app_id)

  local_comm = MPI_COMM_NULL
  if(app_id > -1) then 
    call MPI_Comm_split(world_comm, app_id, world_rank, local_comm, error)
  else
    local_comm = MPI_COMM_WORLD
  endif 

  call MPI_Comm_size(local_comm, local_size, error)
  call MPI_Comm_rank(local_comm, local_rank, error)

  !> Count apps amount 
  app_root =  0 
  napps    = -66
  if(local_rank == 0) app_root = 1
  call MPI_Reduce(app_root, napps, 1, MPI_INTEGER, MPI_SUM, 0, world_comm, error)
  call MPI_Bcast(napps, 1, MPI_INTEGER, 0, world_comm, error)

  allocate( timestep(napps) ) 
  call MPI_Barrier(world_comm, error) 

  cplng%n_apps      = -1 
  cplng%app_id      = -1
  cplng%app_names_l = -1
  cplng%app_comm    = -2
  cplng%base_comm   = -3
  allocate( cplng%app_names ) 
  cplng%app_names   =  "-" 
  allocate( cplng%app_timestep ) 
  !allocate( cplng%app_timestep(napps) ) 
  !cplng%app_timestep = -1 
  allocate( cplng%app_info ) 
  allocate( cplng%app_status ) 
  !cplng%app_info     = -1 
  !print *, cplng%app_info

  app_name = "xyz"
  call alya_coupling_mpi_set_create3(0, app_type, app_name, world_comm, local_comm, cplng)
  call MPI_Barrier(world_comm, error) 
  !print *, world_rank, cplng%app_id, cplng%n_apps, cplng%app_names, cplng%app_status!, cplng%app_info 
  !print *, world_rank, cplng%app_id, cplng%n_apps, cplng%app_status!, cplng%app_info 
  !print *, world_rank, cplng%app_id, cplng%n_apps, cplng%app_info 

  !napps = -1
  !call alya_coupling_mpi_set_get_app_id(cplng, napps) 
  !call alya_coupling_mpi_set_n_apps(cplng, napps) 
  
  !timestep = -3.14159
  !call alya_coupling_mpi_set_get_timestep(cplng, timestep) 
  !print *, "f2: ", timestep 
  !print *, cplng%app_timestep

  !napps = -3
  !call alya_coupling_mpi_get_app_comm(cplng, napps)
  !print *, napps 

  !call alya_coupling_mpi_set_get_info(cplng, app_id, infoij) 
  !print *, world_rank, cplng%app_id, app_id, local_rank, infoij%root_rank, infoij%app_type
  !call MPI_Barrier(world_comm, error) 

  !> get all information 
  call alya_coupling_mpi_set_n_apps(cplng, napps)
  allocate( info(napps) ) 
  allocate( rangeij(napps,2) ) 

  do appj = 1,napps
    call alya_coupling_mpi_set_get_info(cplng, appj-1, info(appj) ) 
    !print *, info(appj)%root_rank, info(appj)%n_ranks, info(appj)%app_type
    rangeij(appj,1) = info(appj)%root_rank
    rangeij(appj,2) = info(appj)%root_rank + info(appj)%n_ranks
  enddo 

  !> interaction: 1<->2 
  if(app_id==1 .or. app_id==0) then 
    appi = app_id 
    appj = mod(app_id+1,2)
    rangei(:) = rangeij(appi+1,:)
    rangej(:) = rangeij(appj+1,:)
    call alya_coupling_mpi_intracomm_create(local_comm, rangei, rangej, commij) 
  endif

  !> reduction commij -> world_comm
  call MPI_Comm_size(commij, sizeij, error)
  call MPI_Comm_rank(commij, rankij, error)

  app_root = -1 
  call MPI_Reduce(app_root, napps, 1, MPI_INTEGER, MPI_SUM, 0, commij, error)
  call MPI_Bcast(napps, 1, MPI_INTEGER, 0, world_comm, error)  
  print *, world_rank, local_rank, rankij, napps 

  !> End 
  call MPI_Barrier(world_comm, error)
  !call alya_coupling_mpi_set_destroy(cplng)
  call MPI_Finalize(error)
end program main 
