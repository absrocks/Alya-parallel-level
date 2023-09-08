!<
!< 2017ABR05  
!< 2017SEP22  
!<
module mod_coupling_timer 
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
  use def_master,           only: cpu_initi

  use mod_parall,           only: PAR_COMM_MY_CODE, PAR_CODE_SIZE
  use mod_std

  implicit none
#ifndef MPI_OFF
  include 'mpif.h'
#endif
  !-----------------------------------------------------------------------||---!
  !
  integer(ip)            :: COU_STATISTICS = 0_ip                              !  ON OFF  
  logical(ip), parameter :: debug          = .False.  
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
  real(rp)    :: ttime(3) = 0_rp
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
    integer(ip) ::   idx = -1_ip
    character(64) :: sms = '' 
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
    subroutine func_template02( xx, icolo, jcolo, coupling )
      use def_kintyp,         only : ip, rp 
      use def_coupli,         only : typ_color_coupling
      implicit none
      integer(ip),              intent(in)    :: icolo
      integer(ip),              intent(in)    :: jcolo
      real(rp),    pointer,     intent(in)    :: xx(:,:)
      type(typ_color_coupling), intent(inout) :: coupling
    end subroutine
  end interface
  !
!  interface coupling_timer_set_function  
!     module procedure &
!                      coupling_timer_set_function_01  ,   &
!                      coupling_timer_set_function_02 !,   &
!  end interface
  !
  !-----------------------------------------------------------------------||---!
  !
  public :: coupling_timer_driver 
!  public :: coupling_timer_set_function  

  public :: coupling_readat   
  public :: coupling_timer_dynamic_statistics
  public :: coupling_timer_set_clik
  public :: coupling_timer_set_clak

  !=============================================================| contains |===!
  contains

!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  subroutine coupling_readat( )
  implicit none

  COU_STATISTICS = 1_ip   

  if(INOTSLAVE.and.COU_STATISTICS>0) write(*,'("[coupling_readat] ",A,E14.7)') "'"//trim(title)//"'"//" STATISTICS ON!!"

  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  subroutine coupling_timer_set_function0( )
  implicit none

  click = huge(1.0_rp)  
  call cputim( click )

  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  subroutine get_dtime_00( idx, dt )
  implicit none
  integer(ip),         intent(in   ) :: idx
  real(rp), optional,  intent(inout) :: dt
  real(rp) :: clack, dtime

  clack = huge(1.0_rp)
  call cputim( clack )

  dtime = clack - click

  if(present(dt)) dt = dtime

  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  subroutine coupling_timer_driver( current_when, current_task ) 
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
if(COU_STATISTICS>0) then
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
  !
  sms = "'"//trim(title) &
           //"."//trim(saux(3)) &
           //"."//trim(namod(modul)) &
           //"."//trim(saux(1))//trim(name_task(current_task))//trim(saux(2)) &
           //"'"
  sms = trim(sms)
  ! 
  !-----------------------------------------------------------------------||---!
  call coupling_timer_turnon( current_when, current_task, sms ) 
  call coupling_timer_turnof( current_when, current_task, sms )  
  !
  call coupling_timer_set_task( current_when, current_task, ITASK_DOITER, DT_TASKS_BY_RANK(1:MMODU,1), sms )
  call coupling_timer_set_task( current_when, current_task, ITASK_BEGZON, DT_TASKS_BY_RANK(1:MMODU,2), sms )
  call coupling_timer_set_task( current_when, current_task, ITASK_ENDZON, DT_TASKS_BY_RANK(1:MMODU,3), sms )
 !call coupling_timer_set_task( current_when, current_task, ITASK_BEGSTE, DT_TASKS_BY_RANK(1:MMODU,4), sms )
  ! 
  !-----------------------------------------------------------------------||---!
endif
  ! 
  end subroutine


  !-----------------------------------------------------------------------||---!
  subroutine coupling_timer_turnon( current_when, current_task, sms )
  use def_domain,           only: npoin
  use def_coupli,           only:  coupling_type
  implicit none
  integer(ip),   intent(in)  :: current_when
  integer(ip),   intent(in)  :: current_task
  character(64), intent(in)  :: sms
  character(64) :: aux  
  integer(4)    :: istat4
  real(rp)      :: toSend 

  logical(ip)             :: code_j = .False. 
  integer(ip)             :: n_wets = -1 
!  integer(ip), intent(in) :: icoup
  integer(ip)             :: icoup = -1 

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
      if(debug.and.INOTSLAVE) print *, "[coupling_timer_turnon]", MPI_SIZE  !PAR_CODE_SIZE
      ! 
      toSend   = npoin*1.0_rp 
      aux(1:8) = "npart"  
      call coupling_timer_gatherToFile02( toSend, aux(1:8) )
      !
!      code_j = current_code == coupling_type(icoup) % code_target
!      if(code_j) n_wets = coupling_type(icoup)%geome%npoin_wet
!      call coupling_timer_dynamic_statistics( -1_ip, -1_ip )
      !
    endif
  endif
  !-----------------------------------------------------------------------||---!
  end subroutine


  !-----------------------------------------------------------------------||---!
  subroutine coupling_timer_turnof( current_when, current_task, sms )
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
      TOTAL_DT = 0.0
      ! 
      if(.not.ended) then 
        ! All is necessary only ONCE ...  
        call cputim( ALYA_END(1) )
        ALYA_END(1) = ALYA_END(1) - cpu_initi
        !
       !call coupling_timer_gatherToFile( DT_FUNCTION_BY_RANK(1)%dtime, current_modul=mmodu+1_ip, current_task=21_ip )
        call coupling_timer_gatherToFile( DT_FUNCTION_BY_RANK(2)%dtime, current_modul=mmodu+2_ip, current_task=22_ip )
        call coupling_timer_gatherToFile(                  ALYA_END(1), current_modul=mmodu+3_ip, current_task=23_ip )
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
  subroutine coupling_timer_set_clik(  )
  implicit none
  !
  call coupling_timer_set_function0
  !
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  subroutine coupling_timer_set_clak( idx, sms)
  implicit none
  integer(ip),   intent(in)  :: idx   
  character(*), intent(in)  :: sms
  real(rp)  :: dt
  !
  if(COU_STATISTICS>0) then  !< 2017SEP22    
    !
    if (idx>N_DT_FUNCTION_BY_RANK) call runend("ERROR:[coupling_timer_set_clak] idx>N_DT_FUNCTION_BY_RANK")
    !
    dt = huge(1.0_rp)
    call get_dtime_00(-1_ip, dt)
    !
    DT_FUNCTION_BY_RANK(idx)%idx   = idx
    DT_FUNCTION_BY_RANK(idx)%dtime = dt  
    !
    if(INOTSLAVE) write(*,'("[coupling_timer_set_clak] ",A,I10,E14.7)') trim(sms), idx, dt  
    !  
  endif 
  !
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  subroutine coupling_timer_set_task( current_when, current_task, itask, CLOCK, sms )
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
  if(.True.) then !< 2017SEP22
    ttime(1) = ttime(2)
    ttime(2) = huge(1.0_rp)
    call cputim( ttime(2) )
    ttime(2) = ttime(2) - cpu_initi
    ttime(3) = ttime(2) - ttime(1)

    !if( INOTSLAVE.and.(current_task>ITASK_TIMSTE.and.current_task<ITASK_ENDSTE).or.(current_task==ITASK_BEGZON.or.current_task==ITASK_ENDZON) ) &
    write(*,'("|_",A, 3F12.2)') trim(sms), ttime(:)

    call coupling_timer_gatherToFile( ttime(3), current_when, current_task )
  endif
    !
  if(.False.)then   !< 2017SEP22
    if(current_when==ITASK_BEFORE) call coupling_timer_set_function0  
    if(current_when==ITASK_AFTER ) then
      dt = huge(1.0_rp)
      call get_dtime_00(-1_ip, dt)
      CLOCK(modul)%dtime = dt
      CLOCK(modul)%calls = CLOCK(modul)%calls + 1
      !
      call coupling_timer_gatherToFile( CLOCK(modul)%dtime, current_when, current_task )
      !
      if(debug.and.INOTSLAVE) write(*,'("[coupling_timer_set_task]",A,I10,E14.7)') trim(sms), CLOCK(modul)%calls, CLOCK(modul)%dtime
      !
    endif
  endif  
    ! 
  endif
  !-----------------------------------------------------------------------||---!
  end subroutine

  !-----------------------------------------------------------------------||---!
  subroutine coupling_timer_reduce_data_n( ARRAY, ARRAY_SIZE, ARRAY_SUM )
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
    aux(1:ARRAY_SIZE)   = huge(0.0_rp)
    call MPI_Allreduce( ARRAY(1:ARRAY_SIZE), aux(1:ARRAY_SIZE), ARRAY_SIZE, MPI_DOUBLE_PRECISION, MPI_SUM, PAR_COMM_MY_CODE, istat4);
    ARRAY(1:ARRAY_SIZE) = int(aux(1:ARRAY_SIZE),rp) / MPI_SIZE
    deallocate( aux )

    if(present(ARRAY_SUM)) ARRAY_SUM = sum( ARRAY(1:ARRAY_SIZE) )
  !-----------------------------------------------------------------------||---!
#endif 
  end subroutine

  !-----------------------------------------------------------------------||---!
  subroutine coupling_timer_gatherToFile( toSend, current_when, current_task, current_modul )
  use def_coupli,  only: coupling_driver_iteration, mcoup 

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
   !write(*,'(A, 3F12.2)') " [coupling_timer_gatherToFile] '"//trim(title)//"' ", minval(toRecv,mask=toRecv>0), sum(toRecv,mask=toRecv>0)/count(toRecv>0,KIND=ip), maxval(toRecv)  
    if(launched) then  !< 2017SEP22
      !!write(*,'("[coupling_timer_gatherToFile]",A,3E14.7)') "'"//trim(title)//"' ", minval(toRecv,mask=toRecv>0), sum(toRecv,mask=toRecv>0)/count(toRecv>0,KIND=ip), maxval(toRecv) 
      !write(*,"('+',A,3I8, ES14.3)") "'"//trim(title)//"' ", ittim, currentModul, currentTask, sum(toRecv,mask=toRecv>0)/count(toRecv>0,KIND=ip)
      !
      if (mcoup > 0) then  
        write(filehandle,trim(fileFormat)) coupling_driver_iteration(iblok), ittim, currentModul, currentTask, currentWhen, ( toRecv(ii), ii=1,MPI_SIZE ) 
      else 
        write(filehandle,trim(fileFormat))                            iblok, ittim, currentModul, currentTask, currentWhen, ( toRecv(ii), ii=1,MPI_SIZE )
      endif
      !
    else 
      call runend("ERROR: [coupling_timer_gatherToFile] NOT launched")  
    endif 
    ! 
  endif
  !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  subroutine coupling_timer_gatherToFile02( toSend, extension )
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
      call runend("ERROR: [coupling_timer_gatherToFile02] NOT launched")
    endif
    ! 
  endif
  !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!



  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!  
  subroutine coupling_timer_dynamic_statistics( n_ij, n_ji )
  use mod_communications, only : PAR_SUM, PAR_MAX, PAR_MIN, PAR_GATHER
  use mod_parall,         only : PAR_COMM_MY_CODE, PAR_CODE_SIZE
  use def_master,         only : title, INOTSLAVE
  use def_domain,         only : npoin

  implicit none
  integer(ip), intent(in   )   :: n_ij ! PLEPP%n_ij   = n_send 
  integer(ip), intent(in   )   :: n_ji ! PLEPP%n_ji   = n_recv 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  integer(4)  :: MPI_RANK, MPI_SIZE, istat4
  integer(ip) :: Nsend, Nrecv, Ncoupled, Npts, Ntrscts
  integer(ip) :: Msend, Mrecv, Mtrscts
  integer(ip) :: ii
  !
  real(rp)    :: Asend, Arecv, Apts, Atrscts
  !
  real(rp)    :: toSend
  real(rp), pointer :: toRecv(:) => null()
  ! 
  character(len=*), parameter :: FMT2 = '(A, 2I5, 4F12.2)'
  character(100)     :: filename
  integer, parameter :: filehandle = 4321
  ! 
  !-----------------------------------------------------------------------||---!
  ! 
if(COU_STATISTICS>0) then
  MPI_RANK = 1
  MPI_SIZE = 1
#ifndef MPI_OFF
  call MPI_Comm_size(PAR_COMM_MY_CODE, MPI_SIZE, istat4)
  call MPI_Comm_rank(PAR_COMM_MY_CODE, MPI_RANK, istat4)
#endif

  Npts = npoin
  call PAR_SUM(Npts,'IN MY CODE')

  Ncoupled = 0
  if( (INOTSLAVE.or.ISEQUEN).and.n_ij > 0) Ncoupled=1 ! coupled subdomain  
  call PAR_SUM(Ncoupled,'IN MY CODE')

  Msend = 0
  if(INOTSLAVE.or.ISEQUEN) Msend = n_ij
  Mrecv = n_ji
  call PAR_MAX(Msend,'IN MY CODE')
  call PAR_MAX(Mrecv,'IN MY CODE')

  Nsend = 0
  if(INOTSLAVE.or.ISEQUEN) Nsend = n_ij
  Nrecv = n_ji
  call PAR_SUM(Nsend,'IN MY CODE')
  call PAR_SUM(Nrecv,'IN MY CODE')

  Ntrscts = huge( 0_ip )
!
!#ifdef COMMDOM 
!  call commdom_locator_get_get_n_intersects( Ntrscts )
!#endif 
!
  call PAR_SUM(Ntrscts,'IN MY CODE')
  call PAR_MAX(Mtrscts,'IN MY CODE')

  if(ISEQUEN) then
    Apts    = Npts
    Asend   = Nsend
    Arecv   = Nrecv
    Atrscts = Ntrscts
  else
    Apts    = Npts    / (MPI_SIZE-1)
!    Asend   = Nsend   / Ncoupled
!    Arecv   = Nrecv   / Ncoupled
!    Atrscts = Ntrscts / Ncoupled
  endif

 !if(INOTSLAVE) then
 !  print *, "[commdom_dynamic_statistics] '", trim(title),"' ",  MPI_SIZE,
 !  Ncoupled, Asend, Arecv, Atrscts, Apts
 !endif 

  if(INOTSLAVE.or.ISEQUEN) then
    if(.not.associated(toRecv)) allocate( toRecv(MPI_SIZE) )
    toRecv   = -1
  else
    if(.not.associated(toRecv)) allocate( toRecv(0) )
  endif

  toSend = n_ij  
  if(INOTSLAVE.or.ISEQUEN)   toSend   = 0  
  call PAR_GATHER(toSend, toRecv,'IN MY CODE')

  if(INOTSLAVE.or.ISEQUEN) then
!    write(*,FMT2) " [commdom_dynamic_statistics] '"//trim(title)//"' ", MPI_SIZE, Ncoupled, Atrscts, Asend, Arecv, Apts
    write(filename,'(a,"_",i6.6,".ncou")') trim(title), MPI_SIZE
    open(filehandle, file=trim(filename), STATUS='REPLACE')
    write(filehandle,*) "##MPI_SIZE, Ncoupled, Atrscts, Asend, Arecv, Apts"
   !write(filehandle,*) MPI_SIZE, Ncoupled, Atrscts, Asend, Arecv, Apts 
    write(filehandle,'(A, 2I6.6, 3F12.2)') " ##", MPI_SIZE, Ncoupled, Atrscts,minval(toRecv,mask=toRecv>0), Asend, maxval(toRecv) !< 2017JAN09 
    do ii=1,MPI_SIZE
      write(filehandle,*) ii, int( toRecv(ii) )
    enddo
    close(filehandle)
  endif

!  deallocate( toRecv )
endif
  !-----------------------------------------------------------------------||---!  
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!  


!-------------------------------------------------------------------------||---!
end module mod_coupling_timer 

!-------------------------------------------------------------------------||---!
!
! |_domain
!   |_cou_initialize_coupling  
!     |_cou_init_interpolate_points_values
!   |_cou_define_wet_geometry                                             | cou_define_wet_geometry.f90
!     |_ coupling_type(icoup) % geome % npoin_wet | Number of wet nodes
!  
! 
! |_moduls
!   |_mod_coupling_driver                                 1.0 
!     |_XXXXXX 
!       |_xxx_plugin                                      1.0
!         |_ COU_INTERPOLATE_NODAL_VALUES                 1.0             |  mod_couplings.f90 +1342
!           |_ COU_GET_INTERPOLATE_POINTS_VALUES          0.1 13.7/15.7   |  mod_couplings.f90 +1385 
!              |_ cou_get_interpolate_points_values_rp_0  0.1 13.7/15.7   |  mod_interpolation.f90  
!           |_ COU_UPDATE_POINTS_VALUES                   0.1 01.6/15.7   |  mod_couplings.f90 +1390 
!               |_ if IQNLS_SCHEME                                        |  mod_couplings.f90 +771
!                 |_ compute_alpha                        0.1 01.6/15.7
!                   |_par_max                             0.1 01.6/15.7   
!
! |_mod_coupling_driver
!   |_ coupling_timer_driver
!     |_ call coupling_timer_turnon( current_when, current_task, 
!     |_ call coupling_timer_turnof( current_when, current_task, 
!     |_ call coupling_timer_set_task( current_when, current_task, ITASK_DOITER,
!     |_ call coupling_timer_set_task( current_when, current_task, ITASK_BEGZON,
!     |_ call coupling_timer_set_task( current_when, current_task, ITASK_ENDZON,
!
! 
! |_ domain 
!   |_ cou_initialize_coupling 
!     |_ COU_INIT_INTERPOLATE_POINTS_VALUES 
!       |_ 
!
!
!-------------------------------------------------------------------------||---!
