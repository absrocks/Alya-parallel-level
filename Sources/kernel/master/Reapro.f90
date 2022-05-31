subroutine reapro()
  !-----------------------------------------------------------------------
  !****f* master/reapro
  ! NAME
  !    Reapro
  ! DESCRIPTION
  !    This routine starts reading data.
  ! USES
  !    inirun   to perform some initializations.
  !    openfi   to get file names and open them.
  !    rrudat   to read run data.
  !    rproda   to read general problem data.
  !    cputim
  !    Nastin
  !    Temper
  !    Codire
  !    Alefor
  ! USED BY
  !    Alya
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use mod_ker_timeline,       only : ker_timeline_synchronization
  use mod_messages,           only : livinf
  use mod_messages,           only : messages_live
  use mod_messages,           only : messages_header
  use mod_output_postprocess, only : output_postprocess_allocate
  implicit none
  !
  ! Check if process was initiated by MPI
  !
  call par_initialize_mpi()
  !
  ! Message header
  !
  call messages_header()
  !
  ! Get data file name and open it
  !
  call openfi(1_ip)
  !
  ! Read run data
  !
  call rrudat()
  !
  ! Get result file names and open them
  !
  call openfi(2_ip)
  !
  ! Live information
  !
  call livinf(1_ip,' ',zero)
  call messages_live('READ PROBLEM DATA')
  !
  ! Checkpoint for Parall communication
  !
  call par_checkpoint()
  !
  ! Checkpoint for OpenMP
  !
  call par_initialize_omp()

  !----------------------------------------------------------------------
  !
  ! Read main data file *.dat
  !
  !----------------------------------------------------------------------
  !
  ! Read general problem data
  !
  call readat()
  !
  ! Read MPI IO data
  !
  call reampio()
  !
  ! Modules: read data
  !
  do modul = 1,mmodu
     call read_module_options()
  end do
  modul = 0
  !
  ! Parallelization
  !
  call par_turnon()
  
  !----------------------------------------------------------------------
  !
  ! Some initializations
  !
  !----------------------------------------------------------------------
  !
  ! Define some memory output options and open files
  !
  call openfi(9_ip)
  !
  ! Block ordering
  !
  call modser()
  !
  ! Allocate and initialize postp type
  !
  call output_postprocess_allocate()
  !
  ! Initial variable
  !
  call inivar(1_ip)
  !
  ! Check errors
  !
  call outerr(0_ip)
  !
  ! Synchronization of timeline
  !
  call ker_timeline_synchronization()

end subroutine reapro
