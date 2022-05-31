module mod_clock 
  use def_kintyp,           only: ip,rp
  use def_master,           only: mmodu
  use def_master,           only: iblok, ittim, itcou
  use def_master,           only: modul, current_code
  use def_master,           only: title, INOTSLAVE
  use def_master,           only: namod, mmodu
  use def_master,           only: ITASK_INIUNK, ITASK_TURNOF
  use def_master,           only: ITASK_TIMSTE, ITASK_ENDSTE
  use def_master,           only: ITASK_BEGZON, ITASK_ENDZON
  use def_master,           only: ITASK_AFTER,  ITASK_BEFORE
  use def_master,           only: ITASK_BEGSTE, ITASK_CONBLK
  use def_master,           only: ITASK_TURNON
  use def_master,           only: ID_KERMOD, ITASK_DOITER, ITASK_CONCOU
  use def_master,           only: ID_NASTIN, ID_ALEFOR, ID_SOLIDZ, ID_TEMPER
  use mod_outfor,           only : outfor

  use mod_parall, only : PAR_COMM_MY_CODE

  use mod_parall,         only : PAR_CODE_SIZE

  implicit none
#ifndef MPI_OFF
  include 'mpif.h'
#endif
  !-----------------------------------------------------------------------||---!
  !
  logical(ip), parameter :: debug    = .false.  
  ! 
  character(6) :: name_task(20) = (/ 'REAPRO', 'TURNON', 'INIUNK', &
                                     'TIMSTE', 'BEGSTE', 'DOITER', &
                                     'CONCOU', 'CONBLK', 'NEWMSH', &
                                     'ENDSTE', 'FILTER', 'OUTPUT', &
                                     'TURNOF', 'BEGITE', 'ENDITE', &
                                     'MATRIX', 'DOOPTI', 'ENDOPT', &
                                     'BEGZON', 'ENDZON' /)
  !
  real(rp)    :: click, clack   
  real(rp)    :: aux_reducer(MMODU)
  real(rp)    :: ALYA_END(1) = huge(0_rp) 
  !
  integer(4)  :: MPI_RANK, MPI_SIZE  
  !
  logical(ip) :: inited = .false. 
  logical(ip) :: ended  = .false.
  !
  character(100)     :: filename
  integer, parameter :: filehandle = 4321  
  !
  type CLOCK_REGISTER
    real(rp)    :: dtime = 0.0 !huge(0_rp)
    integer(ip) :: calls = 0   !huge(0_ip)
  end type 

  integer(ip),     parameter :: N_DT_TASKS_BY_RANK = 3  
  type(CLOCK_REGISTER), save ::   DT_TASKS_BY_RANK(MMODU,N_DT_TASKS_BY_RANK)
  integer(ip),     parameter :: N_CLOCK_COMMDOM    = 2 
  type(CLOCK_REGISTER), save ::      CLOCK_COMMDOM(N_CLOCK_COMMDOM)

  interface
    subroutine func_template00(       )
      implicit none
      ! 
    end subroutine
  end interface
  !
  interface
    subroutine func_template01( itask )
      use def_kintyp,         only : ip
      implicit none
      integer(ip), intent(in) :: itask
    end subroutine  
  end interface  
  !
  interface
    subroutine func_template02( CPLNG )
      use mod_commdom_alya,     only: COMMDOM_COUPLING
      implicit none
      type(COMMDOM_COUPLING), intent(inout) :: CPLNG
    end subroutine
  end interface


  !-----------------------------------------------------------------------||---!
  !
  interface set_dtime  
     module procedure &
                      set_time_params_0,   & 
                      set_time_params_1, & 
                      set_time_params_2   
  end interface  
  !
  interface get_dtime
     module procedure &
                      get_dtime_00 !,   &
  end interface 


  public :: set_dtime 
  public :: get_dtime
  public :: clock_driver 

  !=============================================================| contains |===!
  contains

!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  subroutine set_time_params_0( )
  implicit none

  click = huge(1_ip)  
  call cputim( click )

  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  subroutine get_dtime_00( idx, dt )
  implicit none
  integer(ip),         intent(in   ) :: idx
  real(rp), optional,  intent(inout) :: dt
  real(rp) :: clack, dtime

  clack = huge(1_ip)
  call cputim( clack )

  dtime = clack - click

  if(present(dt)) dt = dtime

  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  subroutine set_time_params_1( f_param_0 )
  use mod_commdom_alya,     only: COMMDOM_COUPLING
  implicit none
  procedure(func_template00)            :: f_param_0 
  !
  real(rp) :: time01=huge(0_rp), time02=huge(0_rp), dt
  !-----------------------------------------------------------------------||---!
  ! 
  call cputim( time01 )
  call f_param_0()
  call cputim( time02 )
  dt = time02 - time01 
  ! 
  CLOCK_COMMDOM(1)%dtime = CLOCK_COMMDOM(1)%dtime + dt     
  CLOCK_COMMDOM(1)%calls = CLOCK_COMMDOM(1)%calls +  1 
  !
  if(INOTSLAVE) print *, "[set_time_params_0]", CLOCK_COMMDOM(1) 
  !
  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  subroutine set_time_params_2( f_param_1, i_param_2 )
  use mod_commdom_alya,     only: COMMDOM_COUPLING
  implicit none
  procedure(func_template02)            :: f_param_1
  type(COMMDOM_COUPLING), intent(inout) :: i_param_2
  !
  real(rp) :: time01=huge(0_rp), time02=huge(0_rp), dt
  !-----------------------------------------------------------------------||---!
  ! 
  call cputim( time01 )
  call f_param_1( i_param_2 )
  call cputim( time02 )
  dt = time02 - time01
  ! 
  CLOCK_COMMDOM(2)%dtime = CLOCK_COMMDOM(2)%dtime + dt
  CLOCK_COMMDOM(2)%calls = CLOCK_COMMDOM(2)%calls +  1
  !
  if(INOTSLAVE) print *, "[set_time_params_2]", CLOCK_COMMDOM(2)
  !
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  subroutine clock_driver( current_when, current_task ) 
  implicit none
  integer(ip),  intent(in)  :: current_when
  integer(ip),  intent(in)  :: current_task
  !
  character(16)            :: saux(3) = ' '
  character(64)            :: sms     = ' '
  character( 4), parameter :: frmt = '(I2)'
  ! 
  real(rp)                 :: dt 
  integer(ip)              :: imodul
  !-----------------------------------------------------------------------||---!
  if(current_when==ITASK_BEFORE) then
    saux(1) = "+"
    saux(2) = "-"
  else&
  if(current_when==ITASK_AFTER ) then
    saux(1) = "-"
    saux(2) = "+"
  endif
  !
  write(saux(3), frmt) MPI_RANK  

  sms = "'"//trim(title) &
           //"."//trim(saux(3)) &
           //"."//trim(namod(modul)) &
           //"."//trim(saux(1))//trim(name_task(current_task))//saux(2) &
           //"'"
  ! 
  !-----------------------------------------------------------------------||---!
  call clock_turnon( current_when, current_task, sms ) 
  call clock_turnof( current_when, current_task, sms )  

  call clock_set_data( current_when, current_task, ITASK_DOITER, DT_TASKS_BY_RANK(1:MMODU,1), sms )
  call clock_set_data( current_when, current_task, ITASK_BEGZON, DT_TASKS_BY_RANK(1:MMODU,2), sms )
  call clock_set_data( current_when, current_task, ITASK_ENDZON, DT_TASKS_BY_RANK(1:MMODU,3), sms )
 !call clock_set_data( current_when, current_task, ITASK_BEGSTE, DT_TASKS_BY_RANK(1:MMODU,4), sms )

  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  end subroutine


  !-----------------------------------------------------------------------||---!
  subroutine clock_turnon( current_when, current_task, sms )
  implicit none
  integer(ip),   intent(in)  :: current_when
  integer(ip),   intent(in)  :: current_task
  character(64), intent(in)  :: sms  
  integer(4) :: istat4
  !-----------------------------------------------------------------------||---!
  if(current_task==ITASK_TURNON) then
    if(.not.inited) then 
      if(debug.and.INOTSLAVE) print *, sms
#ifndef MPI_OFF
      call MPI_Comm_size(PAR_COMM_MY_CODE, MPI_SIZE, istat4)
      call MPI_Comm_rank(PAR_COMM_MY_CODE, MPI_RANK, istat4)
#endif
      MPI_RANK = MPI_RANK+1
      ! 
      if(INOTSLAVE) then 
        write(filename,'(a,"_",i6.6,".dat")') trim(title), MPI_RANK 
        open(filehandle, file=trim(filename), STATUS='REPLACE')
        if(INOTSLAVE) write(filehandle,*) "# DT_TASKSj = \sum_i (DT_TASK/MODULEi); i=NASTIN, TEMPER, ... "
        if(INOTSLAVE) write(filehandle,*) "# TOTAL_TD  = \sum_j (DT_TASKj);        j=DOITER, BEGZON, ENDZON, BEGSTE"
        if(INOTSLAVE) write(filehandle,*) "#                                         SETMESH(=31), SETDT(=32) "
        if(INOTSLAVE) write(filehandle,*) "#    1        2      3                   4            5         6           7 " 
        if(INOTSLAVE) write(filehandle,*) "# CODE, MODULEi, TASKj, <DT_TASKj/MODULEi>, <DT_TASKSj>, TOTAL_TD, TOTAL_ALYA "
      endif 
      ! 
      inited = .true. 
     !print *, "[clock_turnon]", MPI_RANK, PAR_CODE_SIZE
    endif 
  endif
  !-----------------------------------------------------------------------||---!
  end subroutine


  !-----------------------------------------------------------------------||---!
  subroutine clock_turnof( current_when, current_task, sms )
  use mod_communications, only : PAR_SUM
  use mod_parall,         only : PAR_CODE_SIZE
  use def_master,         only : cpu_initi  
  implicit none
  integer(ip),   intent(in)  :: current_when
  integer(ip),   intent(in)  :: current_task
  character(64), intent(in)  :: sms
  !
  real(rp)    :: AVRG_DT_TASKS(N_DT_TASKS_BY_RANK), TOTAL_DT  
  real(rp)    :: aux02(MMODU) 
  integer(ip) :: iCLOCKS, iCODE 
  integer(4)  :: istat4
  CHARACTER(LEN=*), PARAMETER  :: FMT2 = '( 3I5, 4E14.7 )' 
  !-----------------------------------------------------------------------||---!
  if(current_task==ITASK_TURNOF) then
    ! 
    if((current_when==ITASK_BEFORE).and.(.not.modul==ID_KERMOD)) then
   !if(current_when==ITASK_AFTER) then
      if(INOTSLAVE) write(filehandle,*) "# "//namod(modul)//"_"//name_task(current_task)
      !
      iCODE = current_code
#ifdef COMMDOM 
      if(title(1:5)=='DIRIC') iCODE=1 
      if(title(1:5)=='NEUMA') iCODE=1 
#endif 
      ! 
      TOTAL_DT = 0.0
      ! 
      if(.not.ended) then 
        ! All is necessary only ONCE ...  
        call cputim( ALYA_END(1) )
        ALYA_END(1) = ALYA_END(1) - cpu_initi
        call clock_reduce_data_n( ALYA_END, 1_ip )    
        call clock_reduce_data_n( CLOCK_COMMDOM(1:N_CLOCK_COMMDOM)%dtime, N_CLOCK_COMMDOM )
        !
        ended = .true.
      endif 
      !
      do iCLOCKS = 1,N_CLOCK_COMMDOM
        TOTAL_DT = TOTAL_DT + CLOCK_COMMDOM(iCLOCKS)%dtime
      enddo
      ! 
      !! [DOITER, BEGSTE, ENDSTE] = [ DOITER, SENDRECV ]  
      !! AVRG_DT_TASKS(i) = \sum_j DT_TASKS_BY_RANK(j,i). j=1,...,MMODU 
      AVRG_DT_TASKS = 0.0
      do iCLOCKS = 1,N_DT_TASKS_BY_RANK
        call clock_reduce_data_mmodu( DT_TASKS_BY_RANK(1:MMODU,iCLOCKS), AVRG_DT_TASKS(iCLOCKS) ) 
      enddo
      ! 
      TOTAL_DT = TOTAL_DT + sum(AVRG_DT_TASKS) 
      !
      do iCLOCKS = 1,N_DT_TASKS_BY_RANK
        if(INOTSLAVE) write(         *,FMT2) iCODE,   modul, iCLOCKS, DT_TASKS_BY_RANK(modul,iCLOCKS)%dtime, AVRG_DT_TASKS(iCLOCKS), TOTAL_DT, ALYA_END(1)
        if(INOTSLAVE) write(filehandle,FMT2) iCODE,   modul, iCLOCKS, DT_TASKS_BY_RANK(modul,iCLOCKS)%dtime, AVRG_DT_TASKS(iCLOCKS), TOTAL_DT, ALYA_END(1) 
      enddo 
      !
      !! [SET_MESH]   
      do iCLOCKS = 1,N_CLOCK_COMMDOM 
        if(INOTSLAVE)  write(         *,FMT2) iCODE, MMODU+iCLOCKS, iCLOCKS+N_DT_TASKS_BY_RANK, CLOCK_COMMDOM(iCLOCKS)%dtime, CLOCK_COMMDOM(iCLOCKS)%dtime, TOTAL_DT, ALYA_END(1)
        if(INOTSLAVE)  write(filehandle,FMT2) iCODE, MMODU+iCLOCKS, iCLOCKS+N_DT_TASKS_BY_RANK, CLOCK_COMMDOM(iCLOCKS)%dtime, CLOCK_COMMDOM(iCLOCKS)%dtime, TOTAL_DT, ALYA_END(1)
      enddo 
      if(INOTSLAVE)  write(filehandle,"(/)")  
      ! 
    else if(current_when==ITASK_AFTER) then
     !if(INOTSLAVE) close(filehandle)
    endif 
    !
  endif
  !-----------------------------------------------------------------------||---!
  end subroutine


  !-----------------------------------------------------------------------||---!
  subroutine clock_outfor( dummy )
  use def_master, only : lun_outpu 
  implicit none
  real(rp) :: dummy(2)
  character(128), parameter :: frmt = '("cosa:", f10.2, " ", f10.2)'
  !-----------------------------------------------------------------------||---!

  write(lun_outpu,frmt) dummy  

  !-----------------------------------------------------------------------||---!
  end subroutine


  !-----------------------------------------------------------------------||---!
  subroutine clock_set_data( current_when, current_task, itask, CLOCK, sms )
  implicit none
  integer(ip),   intent(in)  :: current_when
  integer(ip),   intent(in)  :: current_task
  integer(ip),   intent(in)  :: itask   
  character(64), intent(in)  :: sms
  type(CLOCK_REGISTER), intent(inout) :: CLOCK(:) 
  !
  real(rp)                 :: dt
  !-----------------------------------------------------------------------||---!
  if(current_task==ITASK) then
    !
    if(current_when==ITASK_BEFORE) call set_dtime()
    if(current_when==ITASK_AFTER ) then
      dt = huge(1_rp)
      call get_dtime(-1_ip, dt)
      CLOCK(modul)%dtime = CLOCK(modul)%dtime + dt
      CLOCK(modul)%calls = CLOCK(modul)%calls + 1

      if(debug.and.INOTSLAVE) print *, trim(sms), CLOCK(modul)%calls, CLOCK(modul)%dtime
    endif
    ! 
  endif
  !-----------------------------------------------------------------------||---!
  end subroutine


  !-----------------------------------------------------------------------||---!
  subroutine clock_reduce_data_mmodu( CLOCK, CLOCK_SUM )
  implicit none 
  type(CLOCK_REGISTER), intent(inout) :: CLOCK(:)
  real(rp),             intent(  out) :: CLOCK_SUM  
  integer(4) :: istat4
  !-----------------------------------------------------------------------||---!
    aux_reducer(1:MMODU) = huge(0_rp)
#ifndef MPI_OFF
    call MPI_Allreduce( CLOCK(1:MMODU)%dtime, aux_reducer(1:MMODU), MMODU, MPI_DOUBLE_PRECISION, MPI_SUM, PAR_COMM_MY_CODE, istat4);
#endif
    CLOCK(1:MMODU)%dtime = aux_reducer(1:MMODU) / MPI_SIZE

    CLOCK_SUM = sum( CLOCK(1:MMODU)%dtime )
  !-----------------------------------------------------------------------||---!
  end subroutine

  !-----------------------------------------------------------------------||---!
  subroutine clock_reduce_data_n( ARRAY, ARRAY_SIZE, ARRAY_SUM )
  implicit none
  integer(ip),          intent(in   ) :: ARRAY_SIZE
  real(rp),             intent(inout) :: ARRAY(:)
  real(rp), optional,   intent(inout) :: ARRAY_SUM
  !
  real(rp), pointer :: aux(:) => null()
  integer(4)        :: istat4
  !-----------------------------------------------------------------------||---!
    allocate( aux(ARRAY_SIZE) ) 
    aux(1:ARRAY_SIZE)   = huge(0_rp)
#ifndef MPI_OFF
    call MPI_Allreduce( ARRAY(1:ARRAY_SIZE), aux(1:ARRAY_SIZE), ARRAY_SIZE, MPI_DOUBLE_PRECISION, MPI_SUM, PAR_COMM_MY_CODE, istat4);
#endif
    ARRAY(1:ARRAY_SIZE) = aux(1:ARRAY_SIZE) / MPI_SIZE
    deallocate( aux )

    if(present(ARRAY_SUM)) ARRAY_SUM = sum( ARRAY(1:ARRAY_SIZE) )
  !-----------------------------------------------------------------------||---!
  end subroutine

  !-----------------------------------------------------------------------||---!
  subroutine clock_reduce_data_0( SCALAR )
  implicit none
  real(rp),             intent(inout) :: SCALAR  
  !
  integer(4)        :: istat4
  real(rp)  :: aux = huge(0_rp) 
  !-----------------------------------------------------------------------||---!
#ifndef MPI_OFF
    call MPI_Allreduce( SCALAR, aux, 1_4, MPI_DOUBLE_PRECISION, MPI_SUM, PAR_COMM_MY_CODE, istat4 );
#endif
    SCALAR = aux / MPI_SIZE  
  !-----------------------------------------------------------------------||---!
  end subroutine

!-------------------------------------------------------------------------||---!
end module mod_clock  

!
! defmod/def_master.f90 
! 
!  real(rp)                 :: &
!       cpu_initi,             &      ! Initial CPU time
!       cpu_start(4),          &      ! CPU for starting operations
!       cpu_outpu,             &      ! CPU for output operations
!       cpu_other(30),         &      ! CPU time
!       ini_tim,               &      ! Initial time to be used in time measurment with cputim
!       fin_tim,               &      ! Final time to be used in time measurment with cputim
!       rate_time                     ! Time rate computed with system_clock


!kernel/defmod/def_master.f90
!  cpu_modul(40,mmodu),   &      ! Module CPU time


!   call cputim(time1)
! |_ Turnon 
  ! Compute time rate
  !
!  call system_clock(count_rate=count_rate8)
!  rate_time = 1.0_rp / max(real(count_rate8,rp),zeror)

! |_ Turnon 
!   |_ Reapro 
!     |_ inirun
!       |_ call cputim(cpu_initi)         ! Initial CPU time

! |_ cputim 
!   |_ call system_clock(itim8)
!      rtime = real(itim8,rp) * rate_time  


! |_ Turnof
!   |_ outrut/outcpu 
!     |_ 
!       call cputim(cpu_refer)
!       ALYA_END(1) =   cpu_refer - cpu_initi
!       ...
!       call outfor(18_ip,lun_outpu,' ')         ! lun_outpu=12 [Output (log) file unit] 
!       ...
!       call outfor(-19_ip,0_ip,' ') -19?? 
!  


