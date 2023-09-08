!>
!>  @page pagename03 Fortran
!>  @brief  Communication domains [safe fluids interchange]  
!>  @author J. M. Zavala Aké
!>  @see     main_loader.f90
!>  @include main_loader.f90
!>
program main_wrap
  use iso_c_binding
  !use mpi 
  implicit none
  include 'mpif.h'

  interface
    function f_stripperdLength(STRNG, LNGTH) result(STRIPPERD)  bind(c, name='stripperdLength') 
      use, intrinsic :: iso_c_binding
      implicit none
      character(kind=c_char), dimension(*) :: STRNG
      integer(kind=c_int), value           :: LNGTH
      integer(kind=c_int)                  :: STRIPPERD
    end function
  end interface


  integer      error
  character*8  app_name, my_name
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
  integer app_ok, i
  integer min_recv 

  character*32  namei, namej
  character*33  char_send, char_recv 
  character*32  token


  world_comm = MPI_COMM_WORLD
  local_comm = MPI_COMM_WORLD

  call commdom_create()

 
  call MPI_Init(error)
  call MPI_Comm_rank(world_comm, world_rank, error); 
  call MPI_Comm_size(world_comm, world_size, error); 


  my_surname = "ALYA_CFD"
  call getarg(1, my_name)
  !lstring = f_stripperdLength(my_name, 8) 
  !print *, my_name, lstring


  !filename="/home/jmake/z2014/Cplng/MeshCplng/MeshCoarse/disk"
  !idx = 0
  !
  !call c_loader_alya_create_array()
  !call c_loader_alya_set_ple_data(filename, idx) 
  !call c_loader_alya_delete()  
  !
  !filename="/home/jmake/z2014/Cplng/MeshCplng/MeshCoarse/square"
  !call c_loader_alya_create_array()
  !call c_loader_alya_set_ple_data(filename, idx) 
  !call c_loader_alya_delete()  

  local_comm  = -1
  size_commij = -1
  call commdom_create()

  do i = 1, command_argument_count() !iargc()
    token = ""
    call getarg(i, token)
    call commdom_set_argvs(trim(token), len_trim(token)) 
  enddo  

  token = "--name"
  call commdom_analyse_argvs(trim(token), len_trim(token))
  call commdom_get_argvs(my_name); 


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

  min_recv = -1

  !ok = -1 
  !surnamei = my_surname
  !call commdom_my_surname(trim(surnamei), len_trim(surnamei), ok)
!  if(trim(my_name)=="CCCC") then 
!    namej = "AAAA" !WhoIsSending?
  namei = "CCCC"
  call commdom_who_iam(namei, len_trim(namei), app_ok)
  if(app_ok==1) then
    namej = "AAAA"
    call commdom_who_areyou(namej, len_trim(namej), app_ok)
    if(app_ok==1) then
      call commdom_get_commij(trim(namej), len_trim(namej), commij); print *, ""

      send   = 1
      n_send = 0
      n_recv = 1

      call commdom_reduce_min_int(send, min_recv, local_comm, commij)

      call commdom_sendrecv_int(recv, n_send, send, n_recv, local_comm, commij)
      call MPI_Bcast(send,  1, MPI_INTEGER, 0, local_comm, error) 

      char_send = namei
      char_recv = ""
      n_send = 32
      n_recv = 0
      call commdom_sendrecv_char(char_send, n_send, char_recv, n_recv, local_comm, commij)
    endif 
  endif 


!  if(trim(my_name)=="AAAA") then 
!    namej = "CCCC" !WhoIsRecv?
  namei = "AAAA"
  call commdom_who_iam(namei, len_trim(namei), app_ok)
  if(app_ok==1) then
    namej = "CCCC"
    call commdom_who_areyou(namej, len_trim(namej), app_ok)
    if(app_ok==1) then

      call commdom_get_commij(trim(namej), len_trim(namej), commij);

      if(local_rank==0) send = 3

      call commdom_reduce_min_int(send, min_recv, local_comm, commij)

      call MPI_Bcast(send, 1, MPI_INTEGER, 0, local_comm, error) 

      n_send = 1
      n_recv = 0
      call commdom_sendrecv_int(send, n_send, recv, n_recv, local_comm, commij)

      char_send = ""
      char_recv = ""
      n_send = 0
      n_recv = 32
      call commdom_sendrecv_char(char_send, n_send, char_recv, n_recv, local_comm, commij)
      !print *, char_recv, "<---"
      

      if(size_commij==2) then 
        namej = "BBBB" !WhoIsRecv?
        call commdom_who_areyou(namej, len_trim(namej), app_ok)
        if(app_ok==1) then
          call commdom_get_commij(trim(namej), len_trim(namej), commij); print *, ""

          if(local_rank==0) send = 3
          call MPI_Bcast(send,  1, MPI_INTEGER, 0, local_comm, error) 

          who_recv = 0
          if(local_rank==0) call MPI_Send(send, 1, MPI_INTEGER, who_recv, rtag, commij, error) 
        endif 
      endif
      
    endif 
  endif 


!  if(trim(my_name)=="BBBB") then 
!    namej = "AAAA" !WhoIsSending?
  namei = "BBBB"
  call commdom_who_iam(namei, len_trim(namei), app_ok)
  if(app_ok==1) then
    namej = "AAAA"
    call commdom_who_areyou(namej, len_trim(namej), app_ok)
    if(app_ok==1) then
      call commdom_get_commij(trim(namej), len_trim(namej), commij); print *, ""

      send  = 2
      who_recv = 0
      if(local_rank==0) call MPI_Recv(send, 1, MPI_INTEGER, who_send, rtag, commij, MPI_STATUS_IGNORE, error) 

      call MPI_Bcast(send, 1, MPI_INTEGER, 0, local_comm, error) 
    endif 
  endif 


!  if(trim(my_name)=="4TRN") then 
!    namej = "PYTH" !WhoIsSending?
  namei = "4TRN"
  call commdom_who_iam(namei, len_trim(namei), app_ok)
  if(app_ok==1) then
    namej = "PYTH"
    call commdom_who_areyou(namej, len_trim(namej), app_ok)
    if(app_ok==1) then
      call commdom_get_commij(trim(namej), len_trim(namej), commij); print *, ""

      send   = 69
      n_send = 1
      n_recv = 0

      call commdom_sendrecv_int(send, n_send, recv, n_recv, local_comm, commij)
      call MPI_Bcast(send,  1, MPI_INTEGER, 0, local_comm, error) 
    endif 
  endif 


!  if(trim(my_name)=="COUPLING") then
!    namej = "Gases" !WhoIsSending?
  namei = "COUPLING"
  call commdom_who_iam(namei, len_trim(namei), app_ok)
  if(app_ok==1) then
    namej = "Temper01"
    call commdom_who_areyou(namej, len_trim(namej), app_ok)
    if(app_ok==1) then
      call commdom_get_commij(trim(namej), len_trim(namej), commij); print *, ""

      print *, ""
      print *, "Coupling OK!!. kill!!..."
      print *, ""

      char_send = "coupling:type:b  "
      char_recv = ""
      n_send = 33
      n_recv = 33
      call commdom_sendrecv_char(char_send, n_send, char_recv, n_recv, local_comm, commij)
      print *, "from_syrthes:", char_recv, "<-----"
 
    endif 
  endif


  call MPI_Barrier(local_comm, error); 
  print *, "wrank/val/min:", world_rank+1, send, min_recv


  call commdom_delete()
  call MPI_Finalize(error) 
  if(local_rank==0) print *, "OK!!"

end program 
