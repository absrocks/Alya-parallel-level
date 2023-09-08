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
  use mod_ker_timeline, only : ker_timeline_synchronization
  use mod_messages,     only : livinf
  use mod_messages,     only : messages_header
  implicit none
  !
  ! Initializations
  !
  call inirun()
  !
  ! Check if process was initiated by MPI
  !
  call Parall(-1_ip)
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
  ! Live information
  !
  call livinf(1_ip,' ',zero)
  !
  ! Checkpoint for OpenMP
  !
  call par_initialize_omp()
  !
  ! Checkpoint for Parall communication
  !
  call Parall(-2_ip)
  !
  ! Get result file names and open them
  !
  call livinf(2_ip,' ',zero)
  call openfi(2_ip)
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
     call reamod()
  end do
  modul = 0
  !
  ! Services: read data
  !
  call Dodeme(0_ip)
  call Hdfpos(0_ip)
  call Optsol(ITASK_REAPRO)
  call Parall(0_ip) 
  call Adapti(ITASK_REAPRO)
  !
  ! Define some memory output options and open files
  !
  call openfi(9_ip)
  !
  ! Block ordering
  !
  call modser()
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
