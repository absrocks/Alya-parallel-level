program main_wrap
  use iso_c_binding
  !use mpi 
  implicit none
  include 'mpif.h'

  integer      error
  character*8  app_name, my_name, namei, namej
  integer      world_comm, world_rank, world_size
  integer      local_comm, local_rank, local_size
  integer      commij, rankij
  character*256  filename
  character*32  my_surname 
  integer size_commij, who_recv, who_send
  integer  idx

  integer stag, rtag
  integer send, recv, n_send, n_recv
  integer lstring

  world_comm = MPI_COMM_WORLD
  local_comm = MPI_COMM_WORLD
  
  call MPI_Init(error)
  call MPI_Comm_rank(world_comm, world_rank, error); 
  call MPI_Comm_size(world_comm, world_size, error); 

  my_surname = "XXXX_XXX"
  call getarg(1, my_name)

  local_comm  = -1
  size_commij = -1
  call commdom_create()
  call commdom_set_names(trim(my_surname), len_trim(my_surname), trim(my_name), len_trim(my_name))
  call commdom_create_commij(world_comm, local_comm)
  call MPI_Comm_rank(local_comm, local_rank, error)
  call MPI_Comm_size(local_comm, local_size, error)
   
  call commdom_get_commij_size(size_commij)

  !================================================================| COMMS |===!
  call MPI_Barrier(local_comm, error); 
  
  recv = -1
  send = -1
  stag =  123
  rtag =  456
  who_recv = -1
  who_send = -1


  call MPI_Barrier(local_comm, error); 
  print *, "wrank/val f:", world_rank+1, send


  call commdom_delete()
  call MPI_Finalize(error) 
  if(local_rank==0) print *, "OK!!"

end program 
