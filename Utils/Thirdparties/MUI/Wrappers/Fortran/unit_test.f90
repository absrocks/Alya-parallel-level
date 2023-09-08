program main
  use iso_c_binding, only : c_ptr,c_null_char, c_double
  implicit none
  include  'mpif.h'
!  integer, pointer :: uniface
  character(len=1024) arg
  external mui_create_uniface3d_f, mui_push_f, mui_destroy_uniface3d_f

  integer(4)       :: iarg
  character(256)   :: token
  real(8) :: zero(1) 
!  real(c_double) :: zero=0.0_c_double
  type(c_ptr),target :: uniface
  integer(4)    :: world_comm

  call mui_create() 

  do iarg = 1, command_argument_count() !iargc()
    token = ""
    call getarg(iarg, token)
    call mui_set_argvs(trim(token), len_trim(token))
  enddo

  token = "--name"
  call mui_analyse_argvs(trim(token), len_trim(token))

  arg = '' 
  call mui_get_argvs(arg)
!!print*, "-->'", trim(arg), "'" 

  world_comm = MPI_COMM_NULL  
  call mui_create_uniface3d_f(arg, world_comm)

!!call mui_push_f("position", 0., 0., 0., 0.)
  call mui_push_f(uniface, "position", zero(1), zero(1), zero(1), zero(1))


  call mui_destroy_uniface3d_f()

  call mui_delete() 

end program main
