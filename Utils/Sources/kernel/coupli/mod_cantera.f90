!!
!! 2018NOV12. 
!!
module mod_cantera 
  use def_kintyp,           only: ip,rp
  use def_master,           only: mmodu
  use def_master,           only: iblok, ittim, itcou
  use def_master,           only: modul, current_code
  use def_master,           only: namod, mmodu
  use def_master,           only: TITLE, INOTSLAVE, ISEQUEN
  use def_master,           only: ITASK_INIUNK, ITASK_TURNOF
  use def_master,           only: ITASK_TIMSTE, ITASK_ENDSTE
  use def_master,           only: ITASK_BEGZON, ITASK_ENDZON
  use def_master,           only: ITASK_AFTER,  ITASK_BEFORE
  use def_master,           only: ITASK_BEGSTE, ITASK_CONBLK
  use def_master,           only: ITASK_TURNON
  use def_master,           only: ID_KERMOD, ITASK_DOITER, ITASK_CONCOU
  use def_master,           only: ID_NASTIN, ID_ALEFOR, ID_SOLIDZ, ID_TEMPER
  use mod_parall,           only: PAR_COMM_MY_CODE, PAR_CODE_SIZE
  ! 
#ifdef CANTERA  
  use cantera 
#endif 
  !
  implicit none
  !  
#ifndef MPI_OFF
  include 'mpif.h'
#endif
  ! 
  !-----------------------------------------------------------------------||---!
  !
  logical(ip), parameter :: DEBUG      = .True. !.False.  
  ! 
  character(6) :: name_task(20) = (/ 'REAPRO', 'TURNON', 'INIUNK', &
                                     'TIMSTE', 'BEGSTE', 'DOITER', &
                                     'CONCOU', 'CONBLK', 'NEWMSH', &
                                     'ENDSTE', 'FILTER', 'OUTPUT', &
                                     'TURNOF', 'BEGITE', 'ENDITE', &
                                     'MATRIX', 'DOOPTI', 'ENDOPT', &
                                     'BEGZON', 'ENDZON' /)
  !
  integer(4)  :: MPI_RANK, MPI_SIZE  
  !
  logical(ip) :: launched = .false. 
  logical(ip) :: ended    = .false.
  !
  real(rp), pointer :: toRecv(:) => null()
  !
  !-----------------------------------------------------------------------||---!
  !
  public :: cantera_driver 

  !=============================================================| contains |===!
  contains

  !-----------------------------------------------------------------------||---!
  subroutine cantera_driver( current_when, current_task ) 
  implicit none  
  integer(ip),  intent(in)  :: current_when  
  integer(ip),  intent(in)  :: current_task  
  character(16)             :: saux(3) = ''  
  character(64)             :: sms     = '?'  
  character( 4), parameter  :: frmt    = '(I2)'  
  !-----------------------------------------------------------------------||---!
  !  
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
!  if(INOTSLAVE) print *, "[cantera_turnon]", sms
  ! 
  !-----------------------------------------------------------------------||---!
  !
  call cantera_turnon( current_when, current_task, sms ) 
  call cantera_turnof( current_when, current_task, sms )  
  !  
  !-----------------------------------------------------------------------||---!
  end subroutine


  !-----------------------------------------------------------------------||---!
  subroutine cantera_turnon( current_when, current_task, sms )
  use def_domain,           only: npoin
  implicit none
  integer(ip),   intent(in)  :: current_when
  integer(ip),   intent(in)  :: current_task
  character(64), intent(in)  :: sms
#ifdef CANTERA 
  type(phase_t) gas
  double precision :: t, p
  !-----------------------------------------------------------------------||---!
  if(current_task==ITASK_TURNON) then
  
    gas = importPhase('h2o2.cti','ohmech')
    t   = 1200.0     ! K
    p   = 101325.0   ! Pa
    call setState_TPX(gas, t, p, 'H2:1, O2:1, AR:2')

    call cantera_alya_create()    

    if(DEBUG.and.INOTSLAVE) print *, "[cantera_turnon]", sms

  endif
#endif
  !-----------------------------------------------------------------------||---!
  end subroutine


  !-----------------------------------------------------------------------||---!
  subroutine cantera_turnof( current_when, current_task, sms )
  use mod_communications, only : PAR_SUM
  use mod_parall,         only : PAR_CODE_SIZE
  use def_master,         only : cpu_initi  
  implicit none
  integer(ip),   intent(in)  :: current_when
  integer(ip),   intent(in)  :: current_task
  character(64), intent(in)  :: sms
  !-----------------------------------------------------------------------||---!
  if(current_task==ITASK_TURNOF) then
    ! 
    !   
  endif
  !-----------------------------------------------------------------------||---!
  end subroutine

  !=============================================================| contains |===!
end module mod_cantera  
