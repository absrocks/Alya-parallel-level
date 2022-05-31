!
! 2016DIC24. 
!          FROM /home/bsc21/bsc21704/z2016/REPOSITORY/ALYAs/ALYA_2016/Sources/kernel
! 2017JAN07 
!          FROM ~/z2016/REPOSITORY/ALYAs/ALYA_2016DEC23/Sources/kernel/coupli/mod_measurements.f90 
! 
!
module mod_measurements 
  use def_kintyp,           only: ip,rp
  use def_master,           only: mmodu
  use def_master,           only: iblok, ittim, itcou
  use def_master,           only: modul, current_code
  use def_master,           only: title, INOTSLAVE, ISEQUEN
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
  use mod_parall,           only: PAR_COMM_MY_CODE, PAR_CODE_SIZE
  use mod_std

  implicit none
#ifndef MPI_OFF
  include 'mpif.h'
#endif
  !-----------------------------------------------------------------------||---!
  !
  logical(ip), parameter :: STATISTICS = .False.                                       !<- ON|OFF. 2017JAN07  
  logical(ip), parameter :: debug      = .False.  
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
  real(rp)    :: ALYA_END(1) = huge(0.0_rp) 
  !
  integer(4)  :: MPI_RANK, MPI_SIZE  
  !
  logical(ip) :: launched = .false. 
  logical(ip) :: ended    = .false.
  !
  character(100)     :: fileName, fileFormat 
  integer, parameter :: fileHandle = 991  
  !
  type CLOCK_REGISTER
    real(rp)    :: dtime =  0.0_rp !huge(0_rp)
    integer(ip) :: calls =  0_ip   !huge(0_ip) 
    integer(ip) ::  when = -1_ip
    integer(ip) :: modul = -1_ip
    integer(ip) ::  task = -1_ip 
  end type 
  ! 
  integer(ip),     parameter :: N_DT_TASKS_BY_RANK    = 3  
  type(CLOCK_REGISTER), save ::   DT_TASKS_BY_RANK(MMODU,N_DT_TASKS_BY_RANK)
  integer(ip),     parameter :: N_DT_FUNCTION_BY_RANK = 2 
  type(CLOCK_REGISTER), save ::   DT_FUNCTION_BY_RANK(N_DT_FUNCTION_BY_RANK)
  !
  real(rp)          :: toSend
  real(rp), pointer :: toRecv(:) => null()
  !
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
  !
  interface measurements_set_function  
     module procedure &
                      measurements_set_function_01  ,   &
                      measurements_set_function_02 !,   &
  end interface
  !
  !-----------------------------------------------------------------------||---!
  !
  public :: measurements_driver 
  public :: measurements_set_function  

  !=============================================================| contains |===!
  contains

!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  subroutine measurements_set_function0( )
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
  subroutine measurements_driver( current_when, current_task ) 
  implicit none
  integer(ip),  intent(in)  :: current_when
  integer(ip),  intent(in)  :: current_task
  !
  character(16)            :: saux(3) = ''
  character(64)            :: sms     = '?'
  character( 4), parameter :: frmt = '(I2)'
  ! 
  real(rp)                 :: dt 
  integer(ip)              :: imodul
  ! 
if(STATISTICS) then
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
           //"."//trim(saux(1))//trim(name_task(current_task))//trim(saux(2)) &
           //"'"
  sms = trim(sms)
  ! 
  !-----------------------------------------------------------------------||---!
  call measurements_turnon( current_when, current_task, sms ) 
  call measurements_turnof( current_when, current_task, sms )  

  call measurements_set_task( current_when, current_task, ITASK_DOITER, DT_TASKS_BY_RANK(1:MMODU,1), sms )
  call measurements_set_task( current_when, current_task, ITASK_BEGZON, DT_TASKS_BY_RANK(1:MMODU,2), sms )
  call measurements_set_task( current_when, current_task, ITASK_ENDZON, DT_TASKS_BY_RANK(1:MMODU,3), sms )
 !call measurements_set_task( current_when, current_task, ITASK_BEGSTE, DT_TASKS_BY_RANK(1:MMODU,4), sms )
  !-----------------------------------------------------------------------||---!
endif
  ! 
  end subroutine


  !-----------------------------------------------------------------------||---!
  subroutine measurements_turnon( current_when, current_task, sms )
  use def_domain,           only: npoin
  implicit none
  integer(ip),   intent(in)  :: current_when
  integer(ip),   intent(in)  :: current_task
  character(64), intent(in)  :: sms
  character(64) :: aux  
  integer(4)    :: istat4
  real(rp)      :: toSend 
  !-----------------------------------------------------------------------||---!
  if(current_task==ITASK_TURNON) then
    if(.not.launched) then 
      if(debug.and.INOTSLAVE) print *, sms
      ! 
#ifndef MPI_OFF
      call MPI_Comm_size(PAR_COMM_MY_CODE, MPI_SIZE, istat4) !< 2016JAN09 
      call MPI_Comm_rank(PAR_COMM_MY_CODE, MPI_RANK, istat4)
#endif
      MPI_RANK = MPI_RANK+1
      ! 
      if(INOTSLAVE.or.ISEQUEN) then
        write(aux,'(I6)') MPI_SIZE  
        fileFormat = trim( "(" //"5I8,"//trim(aux)//"E14.7"//")" )  

        write(fileName,'(a,"_",i6.6,".tms")') trim(title), MPI_SIZE
        open(filehandle, file=trim(filename), STATUS='REPLACE')
        if(INOTSLAVE) write(filehandle,*) "# iblok, ittim, modul, current_task, current_when, MPI_SIZEx(F) "
        ! 
        if(.not.associated(toRecv)) allocate( toRecv(MPI_SIZE) )
        toRecv   = -1
      else
        if(.not.associated(toRecv)) allocate( toRecv(0) )
      endif
      ! 
      launched = .true. 
      if(debug.and.INOTSLAVE) print *, "[measurements_turnon]", MPI_SIZE  !PAR_CODE_SIZE
      ! 
      toSend   = npoin*1.0_rp 
      aux(1:8) = "npart"  
      call measurements_gatherToFile02( toSend, aux(1:8) )
      !
    endif
  endif
  !-----------------------------------------------------------------------||---!
  end subroutine


  !-----------------------------------------------------------------------||---!
  subroutine measurements_turnof( current_when, current_task, sms )
  use mod_communications, only : PAR_SUM
  use mod_parall,         only : PAR_CODE_SIZE
  use def_master,         only : cpu_initi  
! mod_commdom_driver, only : CNT_CPLNG
  implicit none
  integer(ip),   intent(in)  :: current_when
  integer(ip),   intent(in)  :: current_task
  character(64), intent(in)  :: sms
  !
  real(rp)    :: AVRG_DT_TASKS(N_DT_TASKS_BY_RANK), TOTAL_DT  
  real(rp)    :: aux02(MMODU) 
  integer(ip) :: iCLOCKS !, iCODE 
  integer(4)  :: istat4
  CHARACTER(LEN=*), PARAMETER  :: FMT2 = '( 3I5, 4E14.7 )' 
  !-----------------------------------------------------------------------||---!
  if(current_task==ITASK_TURNOF) then
    ! 
    if((current_when==ITASK_BEFORE).and.(.not.modul==ID_KERMOD)) then
      !
!      iCODE = current_code
#ifdef COMMDOM 
!      iCODE = CNT_CPLNG%current_code 
#endif 
      ! 
      TOTAL_DT = 0.0_rp
      ! 
      if(.not.ended) then 
        ! All is necessary only ONCE ...  
        call cputim( ALYA_END(1) )
        ALYA_END(1) = ALYA_END(1) - cpu_initi
        !
       !call measurements_gatherToFile( DT_FUNCTION_BY_RANK(1)%dtime, current_modul=mmodu+1_ip, current_task=21_ip )
        call measurements_gatherToFile( DT_FUNCTION_BY_RANK(2)%dtime, current_modul=mmodu+2_ip, current_task=22_ip )
        call measurements_gatherToFile(                  ALYA_END(1), current_modul=mmodu+3_ip, current_task=23_ip )
        !
        ended = .true.
      endif 
      !
      if(INOTSLAVE) write(filehandle,*) "# "//namod(modul)//"_"//name_task(current_task), modul 
      !
    else if(current_when==ITASK_AFTER) then
     !if(INOTSLAVE) close(filehandle)
    endif 
    !
  endif
  !-----------------------------------------------------------------------||---!
  end subroutine


  !-----------------------------------------------------------------------||---!
  subroutine measurements_set_function_01( f_param_0, current_when, current_task, current_modul )
  use mod_commdom_alya,     only: COMMDOM_COUPLING
  use def_master,           only: mmodu
  implicit none
  procedure(func_template00)            :: f_param_0
  integer(ip),   intent(in)  :: current_when
  integer(ip),   intent(in)  :: current_task
  integer(ip),   intent(in)  :: current_modul 
  !
  real(rp) :: time01=huge(0_rp), time02=huge(0_rp), dt
  !
if(STATISTICS) then
  !-----------------------------------------------------------------------||---!
  call cputim( time01 )
  call f_param_0()
  call cputim( time02 )
  dt = time02 - time01
  ! 
  DT_FUNCTION_BY_RANK(1)%dtime = dt
  DT_FUNCTION_BY_RANK(1)%calls = DT_FUNCTION_BY_RANK(1)%calls +  1 
  DT_FUNCTION_BY_RANK(1)%when  = current_when  
  DT_FUNCTION_BY_RANK(1)%task  = current_task
  DT_FUNCTION_BY_RANK(1)%modul = current_modul 
  !
  call measurements_gatherToFile( DT_FUNCTION_BY_RANK(1)%dtime, current_when, current_task, current_modul )
  !
  if(INOTSLAVE) write(*,'("[measurements_set_func_01] ",A,E14.7,I10)') "'"//trim(title)//"'", DT_FUNCTION_BY_RANK(1)
  !-----------------------------------------------------------------------||---!
else 
  !
  call f_param_0()
  ! 
endif
  ! 
  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  subroutine measurements_set_function_02( f_param_1, i_param_2, current_when, current_task )
  use mod_commdom_alya,     only: COMMDOM_COUPLING
  implicit none
  procedure(func_template02)            :: f_param_1
  type(COMMDOM_COUPLING), intent(inout) :: i_param_2
  integer(ip), optional, intent(in)  :: current_when
  integer(ip), optional, intent(in)  :: current_task
  !
  real(rp) :: time01=huge(0_rp), time02=huge(0_rp), dt
  ! 
if(STATISTICS) then
  !-----------------------------------------------------------------------||---!
  call cputim( time01 )
  call f_param_1( i_param_2 )
  call cputim( time02 )
  dt = time02 - time01
  ! 
  DT_FUNCTION_BY_RANK(2)%dtime = dt 
  DT_FUNCTION_BY_RANK(2)%calls = DT_FUNCTION_BY_RANK(2)%calls +  1
  DT_FUNCTION_BY_RANK(2)%when  = current_when
  DT_FUNCTION_BY_RANK(2)%task  = current_task
 !DT_FUNCTION_BY_RANK(2)%modul = current_modul
  !
  if(INOTSLAVE) write(*,'("[measurements_set_func_02] ",A,E14.7,I10)') "'"//trim(title)//"'", DT_FUNCTION_BY_RANK(2)
  !-----------------------------------------------------------------------||---!
else
  !
  call f_param_1( i_param_2 )
  ! 
endif
  !
  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  subroutine measurements_set_task( current_when, current_task, itask, CLOCK, sms )
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
    if(current_when==ITASK_BEFORE) call measurements_set_function0  
    if(current_when==ITASK_AFTER ) then
      dt = huge(1_rp)
      call get_dtime_00(-1_ip, dt)
      CLOCK(modul)%dtime = dt
!      CLOCK(modul)%dtime = CLOCK(modul)%dtime + dt <- ERROR!! 
      CLOCK(modul)%calls = CLOCK(modul)%calls + 1
      !
      call measurements_gatherToFile( CLOCK(modul)%dtime, current_when, current_task )
      !
      if(debug.and.INOTSLAVE) write(*,'("[measurements_set_task]",A,I10,E14.7)') trim(sms), CLOCK(modul)%calls, CLOCK(modul)%dtime
      !
    endif
    ! 
  endif
  !-----------------------------------------------------------------------||---!
  end subroutine

  !-----------------------------------------------------------------------||---!
  subroutine measurements_reduce_data_n( ARRAY, ARRAY_SIZE, ARRAY_SUM )
  implicit none
  integer(ip),          intent(in   ) :: ARRAY_SIZE
  real(rp),             intent(inout) :: ARRAY(:)
  real(rp), optional,   intent(inout) :: ARRAY_SUM
  !
  real(rp), pointer :: aux(:) => null()
  integer(4)        :: istat4
#ifndef MPI_OFF
  !-----------------------------------------------------------------------||---!
    allocate( aux(ARRAY_SIZE) ) 
    aux(1:ARRAY_SIZE)   = huge(0_rp)
    call MPI_Allreduce( ARRAY(1:ARRAY_SIZE), aux(1:ARRAY_SIZE), ARRAY_SIZE, MPI_DOUBLE_PRECISION, MPI_SUM, PAR_COMM_MY_CODE, istat4);
    ARRAY(1:ARRAY_SIZE) = aux(1:ARRAY_SIZE) / MPI_SIZE
    deallocate( aux )

    if(present(ARRAY_SUM)) ARRAY_SUM = sum( ARRAY(1:ARRAY_SIZE) )
  !-----------------------------------------------------------------------||---!
#endif 
  end subroutine

  !-----------------------------------------------------------------------||---!
  subroutine measurements_gatherToFile( toSend, current_when, current_task, current_modul )
  use mod_communications,   only: PAR_GATHER 
  use def_master,           only: iblok, modul, ittim
  implicit none
  real(rp),              intent(in)  :: toSend
  integer(ip), optional, intent(in)  :: current_when
  integer(ip), optional, intent(in)  :: current_task
  integer(ip), optional, intent(in)  :: current_modul
  !
  integer(ip) :: ii, currentModul, currentTask, currentWhen  
  !-----------------------------------------------------------------------||---!
  !
  currentWhen  = -1 
  if(present(current_when))  currentWhen  =  current_when 
  currentTask  = -1  
  if(present(current_task))  currentTask  =  current_task 
  currentModul = modul 
  if(present(current_modul)) currentModul =  current_modul  
  !   
  !-----------------------------------------------------------------------||---!
  call PAR_GATHER(toSend, toRecv,'IN MY CODE')
  !
  if(INOTSLAVE.or.ISEQUEN) then
   !write(*,'(A, 3F12.2)') " [measurements_gatherToFile] '"//trim(title)//"' ", minval(toRecv,mask=toRecv>0), sum(toRecv,mask=toRecv>0)/count(toRecv>0), maxval(toRecv)  
    if(launched) then
      write(*,'("[measurements_gatherToFile]",A,3E14.7)') "'"//trim(title)//"' ", minval(toRecv,mask=toRecv>0), sum(toRecv,mask=toRecv>0)/count(toRecv>0,KIND=ip), maxval(toRecv) 
      !
      write(filehandle,trim(fileFormat)) iblok, ittim, currentModul, currentTask, currentWhen, ( toRecv(ii), ii=1,MPI_SIZE ) 
    else 
      call runend("ERROR: [measurements_gatherToFile] NOT launched")  
    endif 
    ! 
  endif
  !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  subroutine measurements_gatherToFile02( toSend, extension )
  use mod_communications,   only: PAR_GATHER
  use def_master,           only: iblok, modul, ittim
  implicit none
  real(rp),              intent(in)  :: toSend
  character(8),          intent(in)  :: extension  
  !
  integer(ip) :: ii
  character(100)              :: filename
  integer,          parameter :: fHandle = 997  
  !-----------------------------------------------------------------------||---!
  !
  if(INOTSLAVE.or.ISEQUEN)  toRecv = -1  
  call PAR_GATHER(toSend, toRecv,'IN MY CODE')
  !
  if(INOTSLAVE.or.ISEQUEN) then
    if(launched) then
      write(filename,'(a,"_",i6.6,".",a)') trim(title), MPI_SIZE, trim(extension)
       open(fHandle, file=trim(filename), STATUS='REPLACE')
      write(fHandle,'("# ",3E14.7)') minval(toRecv,mask=toRecv>0), sum(toRecv,mask=toRecv>0)/count(toRecv>0,KIND=ip), maxval(toRecv)
      do ii=1,MPI_SIZE
        write(fHandle,*) ii, int( toRecv(ii) )
      enddo
      close(fHandle)
    else
      call runend("ERROR: [measurements_gatherToFile02] NOT launched")
    endif
    ! 
  endif
  !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


!-------------------------------------------------------------------------||---!
end module mod_measurements  

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
!       rate_time                     ! Time rate computed with system_measurements


!kernel/defmod/def_master.f90
!  cpu_modul(40,mmodu),   &      ! Module CPU time


!   call cputim(time1)
! |_ Turnon 
  ! Compute time rate
  !
!  call system_measurements(count_rate=count_rate8)
!  rate_time = 1.0_rp / max(real(count_rate8,rp),zeror)

! |_ Turnon 
!   |_ Reapro 
!     |_ inirun
!       |_ call cputim(cpu_initi)         ! Initial CPU time

! |_ cputim 
!   |_ call system_measurements(itim8)
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


