!==============================================================================!
!  I am your father...
!
!< 2014Sep05
!< 2015Abr14: add  'n_max_its' 'dynamic'
!< 2015Abr28: add  'commdom_driver_nodes_to_reaction' and 'commdom_bvnat'
!             change inivar.f90 + line 305 'endif'
!
!< 2015May10
!< 2015May11 -> unstick
!< 2015Jul02 -> segmentation fixed!!
!< 2015Jul03 -> exchanged, avoid overlapping
!< 2015JUL17 ->
!< 2016Feb15 ->
!< 2016Feb24 ->
!< 2016MAR29 ->
!< 2016MAR30 ->
!< 2016ABR04 ->
!< 2016ABR12 ->
!< 2016MAR26 ->
!< 2016MAY27 ->
!< 2016JUN02 -> (COMMDOM==2)||(COMMDOM==-2)
!< 2016JUN06 ->
!< 2017JAN07 -> (COMMDOM==3)||(COMMDOM==-3)
!< 2017JAN09 ->
!< 2017JAN10 ->
!< 2017JAN11 ->
!< 2017JUN22 -> commented line to avoid errors due to the new compilation flags for the testsuite
!
!==============================================================================!
  !-----------------------------------------------------------------------||---!
  !   + current_code                                      ___________current_task
  !   |_Alya                                       ______|_____
  !     |_call Turnon()                            ITASK_TURNON  02
  !     |_call Iniunk()                            ITASK_INIUNK  03
  !     |_time: do while
  !       |_call Timste()                          ITASK_TIMSTE  04
  !       |_reset: do
  !       | |_call Begste()                        ITASK_BEGSTE  05
  !       |    |_block: do while
  !       |       |_coupling: do while
  !       |         |_call Begzon()                ITASK_BEGZON  19   _
  !       |         |_modules: do while                              / TASK_BEGITE  14
  !       |           |_call Doiter()              ITASK_DOITER  06-|
  !       |           |_call Concou()              ITASK_CONCOU  07  \_ITASK_ENDITE 15
  !       |         |_call Endzon()                ITASK_ENDZON  20
  !       |       |_call Conblk()                  ITASK_CONBLK  08
  !       |_call Endste()                          ITASK_ENDSTE  10
  !   |_call Turnof()                              ITASK_TURNOF  13
  !          __
  ! BLOCK 3_   |
  !   1 X   |  |--current_block  -> CPLNG%blocks_list
  !   2 Y Z |  |
  !   3 W  _|-----current_module -> CPLNG%moduls_list
  ! END_BLOCK__|
  !
  !-----------------------------------------------------------------------||---!
  !
  ! <code, block, modul, task, when, send|recv>
  !
  !-----------------------------------------------------------------------||---!
module mod_commdom_driver
  use mod_commdom_alya,     only: INONE
  use def_parame,           only: ip, rp
  use def_master,           only: inotmaster, imaster, isequen, title, inotslave
  use def_domain,           only: coord, mnode, ndime, npoin
  use def_domain,           only: ltype, lnods
  use mod_commdom_alya,     only: COMMDOM_COUPLING
  use def_kintyp,           only: soltyp
  use def_master,           only: momod, modul
  use def_kermod,           only: kfl_conta
  use mod_messages,         only: livinf
  use mod_std
#ifdef COMMDOM
  use mod_commdom_plepp,    only: PLEPP_CPLNG
  use mod_commdom_plepp,    only: commdom_plepp_compare_dtinv
  use mod_commdom_plepp,    only: commdom_plepp_set_source_nodes
  use mod_commdom_dynamic,  only: commdom_dynamic_check_fixno
  use mod_commdom_dynamic,  only: commdom_dynamic_set_mesh
  use mod_commdom_dynamic,  only: commdom_dynamic_deallocate
  use mod_commdom_dynamic,  only: commdom_dynamic_exchange02
  use mod_commdom_dynamic,  only: CPLNG_PROPS
#endif
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  implicit none
  !
  logical(ip)  :: debug = .false.
  !
  character(6) :: name_task(20) = (/ 'REAPRO', 'TURNON', 'INIUNK', &
                                     'TIMSTE', 'BEGSTE', 'DOITER', &
                                     'CONCOU', 'CONBLK', 'NEWMSH', &
                                     'ENDSTE', 'FILTER', 'OUTPUT', &
                                     'TURNOF', 'BEGITE', 'ENDITE', &
                                     'MATRIX', 'DOOPTI', 'ENDOPT', &
                                     'BEGZON', 'ENDZON' /)

  character(6) :: name_when(2) = (/ 'BEFORE', 'AFTERE'/)
  !
 !integer(ip), parameter   :: n_times = 9                                      !< 2016MAR26. 7 -> 8. 2016JUN06 8 -> 9
 !integer(ip), parameter   :: n_times = 11                                     !< 2017JAN07. 9 -> 11
  integer(ip), parameter   :: n_times = 12_ip                                  !< 2017Nov19. 11 -> 12
  !
  logical(ip)              :: CNT_SENDRECV(n_times) = .false.
  character(64)            :: CNT_SMS               = ' '
  type(COMMDOM_COUPLING),save  :: CNT_CPLNG
  !
  integer(ip)              :: n_max_its
  logical(ip)              :: dynamic
  !
  logical(ip) :: residual   = .false.
  logical(ip) :: exchanged  = .false.
  !
  type(soltyp), pointer    :: solve(:)
  !
  integer(ip) :: modules_i(3) = -1  !< 2016JUN06, 2016MAR29
  integer(ip) :: modules_j(3) = -1
  !
  integer(ip)              :: N_DOITER  = 0       !< 2017JAN07
  logical(ip)              :: INITIATED = .False. !< 2017JAN07
  !
  real(rp)                 :: dt_inv = -1_ip
  !
  private
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  !
  public :: CNT_SENDRECV
  public :: CNT_SMS
  public :: CNT_CPLNG
  public :: commdom_driver_init_cht
  public :: commdom_driver_init_contact
  public :: commdom_driver_init_fsi        !< 2016MAR29
  public :: commdom_driver_sendrecv
  public :: commdom_driver_exchange02
  public :: commdom_driver_n_fixno
  public :: commdom_driver_get_total_flux
  public :: commdom_driver_get_residual
  public :: commdom_driver_coupling_driver !< 2016Feb15
  public :: commdom_driver_domain_create   !< 2016Feb15
  public :: commdom_driver_domain_init     !< 2016Feb15
  public :: commdom_driver_init            !< 2016MAR30
  public :: commdom_driver_set_mesh        !< 2016MAR30
  public :: commdom_driver_init_cht_parallelMPMD !< 2017JAN07
  public :: N_DOITER, INITIATED                  !< 2017JAN07
  !
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!

  contains
!-------------------------------------------------------------------------||---!
!-----------------------------------------------------------------| PUBLIC |---!
!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  subroutine commdom_driver_init_cht_parallelMPMD( CPLNG )                                  !< 2017JAN07.
  use def_master,       only: ITASK_AFTER,  ITASK_BEFORE
  use def_master,       only: ITASK_BEGZON, ITASK_ENDZON
  use def_coupli,       only: UNKNOWN, RESIDUAL
  use def_coupli,       only: mcoup
  use def_master,       only: ID_TEMPER
  use def_master,       only: current_code
  implicit none
  type(COMMDOM_COUPLING), intent(inout) :: CPLNG
  !-----------------------------------------------------------------------||---!
  !-----------------------------------------------------| DIRIC | HEATF--> |---!
  CPLNG%code_i       =  1_ip        !< CODEi
  CPLNG%what_i       = -RESIDUAL    !<---- never
  !
  CPLNG%send_task_i  =  ITASK_ENDZON
  CPLNG%recv_task_i  =  ITASK_BEGZON
  CPLNG%send_when_i  =  ITASK_BEFORE
  CPLNG%recv_when_i  =  ITASK_BEFORE
  !
  CPLNG%send_modul_i =  ID_TEMPER
  CPLNG%recv_modul_i =  ID_TEMPER
  CPLNG%module_i     =  ID_TEMPER
  !-----------------------------------------------------| NEUMA | TEMPE--> |---!
  CPLNG%code_j       =  2_ip        !< CODEj
  CPLNG%what_j       =  RESIDUAL    !< 'physical' (<0)  or 'numerical'(>0) coupling
  !
  CPLNG%send_task_j  =  ITASK_ENDZON
  CPLNG%recv_task_j  =  ITASK_BEGZON
  CPLNG%send_when_j  =  ITASK_BEFORE
  CPLNG%recv_when_j  =  ITASK_BEFORE
  !
  CPLNG%send_modul_j =  ID_TEMPER
  CPLNG%recv_modul_j =  ID_TEMPER
  CPLNG%module_j     =  ID_TEMPER
  !
  CPLNG%tolerance    =  1e-3_rp
  n_max_its          =  1_ip
  !-----------------------------------------------------------------------||---!
  debug              = .false.
  dynamic            = .false.
  CPLNG%n_dof        =  1_ip        !< D.O.F.
  !-----------------------------------------------------------------------||---!
#ifdef COMMDOM
  if( (PLEPP_CPLNG%app_name(1:5) == 'NEUMA').and.(PLEPP_CPLNG%app_name(1:5) == title(1:5)).and.(current_code==CPLNG%code_j) ) then
  else&
  if( (PLEPP_CPLNG%app_name(1:5) == 'DIRIC').and.(PLEPP_CPLNG%app_name(1:5) == title(1:5)).and.(current_code==CPLNG%code_i) ) then
  else
    print *, "[", trim(title),"] ERROR: 'DIRIC' or 'NEUMA' not found!!", " app_name:'", trim(PLEPP_CPLNG%app_name),"' title:'", trim(title) , "' "
    print *, "or DIRIC==1, NEUMA==2 CODE=", current_code
    stop
  endif
#endif
  !
  if( .not.( (CPLNG%what_i/=RESIDUAL).and.(CPLNG%what_j==RESIDUAL) ) ) then
    print *, "[", trim(title),"] ERROR: ", "[commdom_driver_init_fsi] SET 'what_i==RESIDUAL and what_j/=RESIDUAL' "
    print *, " 'what_i==RESIDUAL and what_j/=RESIDUAL' <-> DIRIC: kfl_bvnat=0, kfl_react=1; NEUMA: kfl_bvnat=1, kfl_react=0"
    stop
  endif
  !-----------------------------------------------------------------------||---!
  if(CPLNG%code_i==current_code) CPLNG%current_fixbo = CPLNG%fixbo_i
  if(CPLNG%code_j==current_code) CPLNG%current_fixbo = CPLNG%fixbo_j
  !
  CPLNG%current_code = current_code !< *.dat CODE: ID_CODE
  current_code       = 1_ip         !< trick!! -> PAR_COLOR_COMMUNICATORS: CODE NUMBER EXCEED NUMBER OF CODES
  mcoup              = 0_ip         !< evoid cou_turnon
  !
  if(IMASTER.or.ISEQUEN) print*, "[commdom_driver_init_cht_parallelMPMD]"
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_driver_init_cht( CPLNG )
  use def_coupli,       only: UNKNOWN, RESIDUAL
  use def_coupli,       only: mcoup
  use def_master,       only: ID_TEMPER
  use def_master,       only: current_code
  implicit none
  type(COMMDOM_COUPLING), intent(inout) :: CPLNG
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  !
  CPLNG%code_i       =  1_ip        !< CODEi
  CPLNG%module_i     =  ID_TEMPER   !< MODULEi
  CPLNG%fixbo_i      = -1_ip
  CPLNG%what_i       = -RESIDUAL    !<---- never
  !
  CPLNG%code_j       =  2_ip        !< CODEj
  CPLNG%module_j     =  ID_TEMPER   !< MODULEj
  CPLNG%fixbo_j      = -1_ip
  CPLNG%what_j       =  RESIDUAL    !< 'physical' (<0)  or 'numerical'(>0) coupling
  !
  CPLNG%tolerance    =  1e-3_rp
  n_max_its          =  1_ip
  !
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  dynamic            = .false.
  CPLNG%n_dof        =  1_ip        !< D.O.F.
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if(CPLNG%code_i==current_code) CPLNG%current_fixbo = CPLNG%fixbo_i
  if(CPLNG%code_j==current_code) CPLNG%current_fixbo = CPLNG%fixbo_j
  !
  CPLNG%current_code = current_code !< *.dat CODE: ID_CODE
  current_code       = 1_ip         !< trick!! -> PAR_COLOR_COMMUNICATORS: CODE NUMBER EXCEED NUMBER OF CODES
  mcoup              = 0_ip         !< evoid cou_turnon
  !
  !-----------------------------------------------------------------------||---!
  !print *, trim(PLEPP_CPLNG%namej), trim(PLEPP_CPLNG%namei), trim(title)
  !-----------------------------------------------------------------------||---!
  if(IMASTER.or.ISEQUEN) print*, "[commdom_driver_init_cht]"
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------
  !> @author
  !> @date
  !> @brief   Initializacion PDN contact
  !> @details Initializacion PDN contact
  !>
  !-----------------------------------------------------------------------

  subroutine commdom_driver_init_contact(CPLNG)

    use def_coupli,       only : UNKNOWN, RESIDUAL
    use def_coupli,       only : mcoup
    use def_master,       only : ID_NASTIN, ID_TEMPER, ID_SOLIDZ
    use def_master,       only : current_code

    implicit none

    type(COMMDOM_COUPLING), intent(inout) :: CPLNG

    !-----------------------------------------------------------------------||---!
    !                                                                            !
    !-----------------------------------------------------------------------||---!
    CPLNG%code_i       =  1_ip        !< CODEi
    CPLNG%module_i     =  ID_SOLIDZ   !< MODULEi
    CPLNG%fixbo_i      = -1_ip
    CPLNG%what_i       =  RESIDUAL    !< 'physical'  or 'numerical' coupling
    !
    CPLNG%code_j       =  2_ip        !< CODEj
    CPLNG%module_j     =  ID_SOLIDZ   !< MODULEj
    CPLNG%fixbo_j      = -1_ip
    CPLNG%what_j       = -RESIDUAL    !<---- never
    !
    CPLNG%tolerance    = -1e-4_rp
    !
   !n_max_its          =  1_ip        !< 2015JUL17
    n_max_its          =  20_ip       !< 2017Nov19
    !
#ifdef COMMDOM
    PLEPP_CPLNG%tol    = 1.0e-4_rp    !< Contact tolerance
#endif
    !
    !-----------------------------------------------------------------------||---!
    !                                                                            !
    !-----------------------------------------------------------------------||---!
    dynamic            = .true.
    CPLNG%n_dof        =  ndime       !< D.O.F.
#ifdef COMMDOM
    CPLNG_PROPS%ndime_props =  1_ip   !< props @ vertex_coord_j
#endif
    !-----------------------------------------------------------------------||---!
    !                                                                            !
    !-----------------------------------------------------------------------||---!
    if(CPLNG%code_i==current_code) CPLNG%current_fixbo = CPLNG%fixbo_i
    if(CPLNG%code_j==current_code) CPLNG%current_fixbo = CPLNG%fixbo_j
    !
    CPLNG%current_code = current_code !< *.dat CODE: ID_CODE
    current_code       = 1_ip         !< trick!! -> PAR_COLOR_COMMUNICATORS: CODE NUMBER EXCEED NUMBER OF CODES
    mcoup              = 0_ip         !< evoid cou_turnon
    !
    if(IMASTER) print*, "[commdom_driver_init_contact]"
    !-----------------------------------------------------------------------||---!
    !                                                                            !
    !-----------------------------------------------------------------------||---!
  end subroutine commdom_driver_init_contact

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_driver_init_fsi(CPLNG)                                    !< 2016MAR29
  use def_coupli,       only: UNKNOWN, RESIDUAL
  use def_coupli,       only: mcoup
  use def_master,       only: ID_NASTIN, ID_TEMPER, ID_SOLIDZ, ID_ALEFOR
  use def_master,       only: current_code
  implicit none
  type(COMMDOM_COUPLING), intent(inout) :: CPLNG
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  !
  ! DIRIC
  CPLNG%code_i       =  2_ip        !< CODEi                                   !< 2016MAR26 1 -> 2
  CPLNG%module_i     =  ID_SOLIDZ   !< MODULEi -->
  modules_i(2)       = -ID_SOLIDZ   !< MODULEi <--
  modules_i(3)       = -ID_TEMPER                                              !< 2016JUN06
  CPLNG%fixbo_i      = -1_ip
  !
  ! NEUMA
  CPLNG%code_j       =  1_ip        !< CODEj                                   !< 2016MAR26 2 -> 1
  CPLNG%module_j     =  ID_NASTIN   !< MODULEj <-- !< 2016MAR30
  modules_j(2)       =  ID_ALEFOR   !< MODULEj -->
  modules_j(3)       = -ID_TEMPER                                              !< 2016JUN06
  CPLNG%fixbo_j      = -1_ip
  !
  ! 'ALGEBRAIC' DIRIC-NEUMA system activation...   !< 2016MAR30
  CPLNG%what_i       =  RESIDUAL    !<
  CPLNG%what_j       = -RESIDUAL    !< 'physical' (<0) or 'numerical' (>0) coupling
  !
  CPLNG%tolerance    = -1e-4_rp
  n_max_its          =  1_ip       !< 2015JUL17
  !
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  debug              = .false.
  dynamic            = .false.
  CPLNG%n_dof        =  ndime      !< D.O.F.                                   !< 2016ABR04
 !CPLNG%n_dof        =  ndime + 1  !< D.O.F.                                   !< 2016JUN06
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
#ifdef COMMDOM
  if( (PLEPP_CPLNG%app_name(1:5) == 'NEUMA').and.(PLEPP_CPLNG%app_name(1:5) == title(1:5)).and.(current_code==CPLNG%code_j) ) then !< 2016MAR26 i -> j
  else&
  if( (PLEPP_CPLNG%app_name(1:5) == 'DIRIC').and.(PLEPP_CPLNG%app_name(1:5) == title(1:5)).and.(current_code==CPLNG%code_i) ) then !< 2016MAR26 j -> i
  else
    print *, "[", trim(title),"] ERROR: 'DIRIC' or 'NEUMA' not found!!", " app_name:'", trim(PLEPP_CPLNG%app_name),"' title:'", trim(title) , "' "
    print *, "or DIRIC==2, NEUMA==1 CODE=", current_code
    stop
  endif
#endif
  !
  !
  if( .not.( (CPLNG%what_i==RESIDUAL).and.(CPLNG%what_j/=RESIDUAL) ) ) then
    print *, "[", trim(title),"] ERROR: ", "[commdom_driver_init_fsi] SET 'what_i==RESIDUAL and what_j/=RESIDUAL' "
    print *, " 'what_i==RESIDUAL and what_j/=RESIDUAL' <-> DIRIC: kfl_bvnat=0, kfl_react=1; NEUMA: kfl_bvnat=1, kfl_react=0"
    stop
  endif
  !
  !-----------------------------------------------------------------------||---!
  if(CPLNG%code_i==current_code) CPLNG%current_fixbo = CPLNG%fixbo_i
  if(CPLNG%code_j==current_code) CPLNG%current_fixbo = CPLNG%fixbo_j
  !
  CPLNG%current_code = current_code !< *.dat CODE: ID_CODE
  current_code       = 1_ip         !< trick!! -> PAR_COLOR_COMMUNICATORS: CODE NUMBER EXCEED NUMBER OF CODES
  mcoup              = 0_ip         !< evoid cou_turnon
  !
  if(IMASTER) print*, "[commdom_driver_init_contact]"
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_driver_sendrecv(CPLNG, current_when, current_task)
  use def_master,           only: iblok, ittim
  use def_master,           only: modul, current_code
  use def_master,           only: nblok
  use def_master,           only: namod, mmodu
  use def_master,           only: ITASK_INIUNK, ITASK_TURNOF
  use def_master,           only: ITASK_TIMSTE, ITASK_ENDSTE
  use def_master,           only: ITASK_BEGZON, ITASK_ENDZON
  use def_master,           only: ITASK_AFTER,  ITASK_BEFORE
  use def_master,           only: ITASK_BEGSTE, ITASK_CONBLK
  use def_master,           only: ID_KERMOD, ITASK_DOITER, ITASK_CONCOU
  use def_master,           only: ID_NASTIN, ID_ALEFOR, ID_SOLIDZ, ID_TEMPER
  !
  use def_master,           only: dtinv
  !
  !< ITERATIONs
  use def_coupli,           only: coupling_driver_iteration
  use mod_measurements,     only : measurements_set_function !< 2017JAN07
  !
  implicit none
  type(COMMDOM_COUPLING), intent(inout) :: CPLNG
  !
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  integer(ip),  intent(in)  :: current_when
  integer(ip),  intent(in)  :: current_task
  !
 !integer(ip), parameter   :: n_times = 6
  integer(ip)              :: itime
  integer(ip)              ::       now(8_ip)
  integer(ip)              ::  the_time(8_ip,n_times) = -1_ip
  logical(ip)              :: sendrecv(n_times) = .false.
  !
  character(16)            :: saux(8_ip) = ' '
  character(64)            :: sms        = ' '
  character( 4), parameter :: frmt = '(I2)'
  !
  CNT_SMS = '+-+-+-+-+-+-+-+-'
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  !
  !  ittim: Current time step
  !  itcou: Current global iteration
  !  iblok: Current block
  !  is_when = itcou==micou(iblok)
  !
  now = (/ CPLNG%current_code, iblok, modul, current_task, current_when, ittim, -1_ip, -1_ip /)
  !
  !-----------------------------------------------------------------------||---!
  write(saux(1), frmt) CPLNG%current_code
  write(saux(2), frmt) iblok
  write(saux(3), frmt) ittim
  !write(saux(4), frmt) coupling_driver_iteration(iblok) !< 2017JUN22
  !
  if(current_when==ITASK_BEFORE) then
    saux(5) = "+"
    saux(6) = "-"
  else&
  if(current_when==ITASK_AFTER ) then
    saux(5) = "-"
    saux(6) = "+"
  endif
  !
  sms = "'"//trim(saux(1))//&
        "."//namod(modul)//&
        "."//trim(saux(5))//name_task(current_task)//trim(saux(6))// &
       !"."//name_when(current_when)//&
        ".B"//trim(saux(2))// &
        ".T"//trim(saux(3))// &
        ".I"//trim(saux(4))//"'"
  !-----------------------------------------------------------------------||---!
  if(nblok>1) then
    print *, "[commdom_driver_sendrecv] ", sms
    print *, "[commdom_driver_sendrecv] ", "nblok==1 !!"
    call runend("EXIT!!")
  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  !-----------------------------------------------------------| ITERATIONS |---!
  !   +
  !   |_Alya
  !     |_call Turnon()
  !     |_call Iniunk()
  !     |_time: do while
  !       |_call Timste()                (-1+)
  !       |_reset: do
  !         |_call Begste()              (+2-)
  !           |_block: do while
  !             |_coupling: do while
  !               |_call Begzon()        (+4-)  AFTER<-never INTO the module: into 'coupli/mod_coupling_driver.f90'
  !                                                                           add 'if(current_when==ITASK_AFTER) call Temper(-1_ip)'
  !               |_modules: do while
  !                 |_call Doiter()
  !                 |_call Concou()      (+7-)
  !               |_call Endzon()        (-5+)
  !             |_call Conblk()          (-6+)
  !       |_call Endste()
  !   |_call Turnof()                    (-3+)
  !
  !-----------------------------------------------------------------------||---!
  if(CPLNG%current_code==CPLNG%code_i) then
    the_time(:,1_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_i, ITASK_TIMSTE, ITASK_AFTER,  ittim, -1_ip, -1_ip /)
    the_time(:,2_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_i, ITASK_BEGSTE, ITASK_BEFORE, ittim, -1_ip, -1_ip /)
    the_time(:,3_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_i, ITASK_TURNOF, ITASK_AFTER,  ittim, -1_ip, -1_ip /)
    the_time(:,4_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_i, ITASK_BEGZON, ITASK_AFTER,  ittim, -1_ip, -1_ip /)
    the_time(:,5_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_i, ITASK_ENDZON, ITASK_AFTER,  ittim, -1_ip, -1_ip /)
    the_time(:,6_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_i, ITASK_CONBLK, ITASK_BEFORE, ittim, -1_ip, -1_ip /)
    !
   !the_time(:,7_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_i, ITASK_ENDZON, ITASK_AFTER,  ittim, -1_ip, -1_ip /)  !< 2017JAN09
    the_time(:,7_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_i, ITASK_ENDZON, ITASK_BEFORE, ittim, -1_ip, -1_ip /)  !< 2017Nov19
    the_time(:,8_ip)  = (/ CPLNG%current_code, iblok,   modules_i(2), ITASK_ENDZON, ITASK_AFTER,  ittim, -1_ip, -1_ip /)  !< 2016MAR26
    the_time(:,9_ip)  = (/ CPLNG%current_code, iblok,   modules_i(3), ITASK_ENDZON, ITASK_AFTER,  ittim, -1_ip, -1_ip /)  !< 2016JUN06
    !
    the_time(:,10  )  = (/ CPLNG%current_code, iblok,  CPLNG%send_modul_i, CPLNG%send_task_i, CPLNG%send_when_i, ittim, -1_ip, -1_ip /)  !< 2017JAN07
    the_time(:,11  )  = (/ CPLNG%current_code, iblok,  CPLNG%recv_modul_i, CPLNG%recv_task_i, CPLNG%recv_when_i, ittim, -1_ip, -1_ip /)  !< 2017JAN07
    !
    the_time(:,12_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_i, ITASK_ENDZON, ITASK_AFTER,  ittim, -1_ip, -1_ip /) !< 2017Nov19
  endif
  if(CPLNG%current_code==CPLNG%code_j) then
    the_time(:,1_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_j, ITASK_TIMSTE, ITASK_AFTER,  ittim, -1_ip, -1_ip /)
    the_time(:,2_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_j, ITASK_BEGSTE, ITASK_BEFORE, ittim, -1_ip, -1_ip /)
    the_time(:,3_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_j, ITASK_TURNOF, ITASK_AFTER,  ittim, -1_ip, -1_ip /)
    the_time(:,4_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_j, ITASK_BEGZON, ITASK_AFTER,  ittim, -1_ip, -1_ip /)
    the_time(:,5_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_j, ITASK_ENDZON, ITASK_AFTER,  ittim, -1_ip, -1_ip /)
    the_time(:,6_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_j, ITASK_CONBLK, ITASK_BEFORE, ittim, -1_ip, -1_ip /)
    !
    the_time(:,7_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_j, ITASK_BEGZON, ITASK_BEFORE, ittim, -1_ip, -1_ip /)  !< 2017JAN09
    the_time(:,8_ip)  = (/ CPLNG%current_code, iblok,   modules_j(2), ITASK_BEGZON, ITASK_BEFORE, ittim, -1_ip, -1_ip /)  !< 2016MAR26
    the_time(:,9_ip)  = (/ CPLNG%current_code, iblok,   modules_j(3), ITASK_BEGZON, ITASK_BEFORE, ittim, -1_ip, -1_ip /)  !< 2016JUN06
    !
    the_time(:,10  )  = (/ CPLNG%current_code, iblok,  CPLNG%send_modul_j, CPLNG%send_task_j, CPLNG%send_when_j, ittim, -1_ip, -1_ip /)  !< 2017JAN07
    the_time(:,11  )  = (/ CPLNG%current_code, iblok,  CPLNG%recv_modul_j, CPLNG%recv_task_j, CPLNG%recv_when_j, ittim, -1_ip, -1_ip /)  !< 2017JAN07
    !
    the_time(:,12_ip)  = (/ CPLNG%current_code, iblok, CPLNG%module_j, ITASK_ENDZON, ITASK_AFTER, ittim, -1_ip, -1_ip /)  !< 2017Nov19
  endif
  !
  sendrecv = (/ (all(the_time(:,itime)==now), itime=1,n_times) /)
  !
  do itime = 1,n_times
    if( sendrecv(itime) ) then
      if(debug) print *, " [", trim(title),"] ",  trim(sms)
    endif
  enddo
  !
  CNT_SENDRECV = sendrecv                                                      !< 2016MAR29
  if(any(sendrecv)) CNT_SMS = trim(sms)
  !
  call commdom_driver_init_exchange( current_when, current_task )              !< 2017JAN07
  !call commdom_driver_time_exchange( current_when, current_task )             !< 2017JAN11
  !
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-------------------------------------------------------------| -TIMSTE+ |---!
#ifdef COMMDOM
  if( sendrecv(1_ip) ) then
#if COMMDOM!=2
     call commdom_localize(CPLNG)     !< 2017Nov19 !< Bilateral
#endif
  endif
  !-------------------------------------------------------------| +BEGSTE- |---!
  if( sendrecv(2_ip) ) then
    call commdom_driver_begste()
  !!call measurements_set_function( commdom_driver_begste, current_when, current_task, modul ) !< 2017JAN07
  endif
  !-------------------------------------------------------------| -TURNOF+ |---!
  if( sendrecv(3_ip) ) then
#if COMMDOM!=-3
  call commdom_plepp_compare_dtinv(dtinv)                                      !< 2017JAN07
#endif
  endif
  !-------------------------------------------------------------| -BEGZON+ |---!
  if( sendrecv(4_ip) ) then
    call commdom_driver_begzon()
  endif
  !-------------------------------------------------------------| -ENDZON+ |---!
  if( sendrecv(5_ip) ) then !.and.exchanged ) then
    call commdom_driver_endzon()
!    exchanged = .false.                                                        !< 2015Jul03
  endif
  !-------------------------------------------------------------| +CONBLK- |---!
  if( sendrecv(6_ip) ) then
    call commdom_driver_conblk()
    !-----------------------------------------------------------| or here? |---!
    !---------------------------------------------------------------------||---!
  endif
  !-------------------------------------------------------------| EXCHANGE |---!
  if( any(sendrecv(7_ip:n_times)) ) then !.and.(.not.exchanged) ) then            !< 2017Nov19
 !if( sendrecv(7_ip) ) then !.and.(.not.exchanged) ) then
    !
#if COMMDOM==2
     if ( kfl_conta == 1 ) then          ! Unilateral
         if ( sendrecv(7_ip) ) call commdom_localize_unilateral(CPLNG)
     else if ( kfl_conta == 2 ) then     ! Bilateral
        if ( sendrecv(7_ip) .or. sendrecv(12_ip) ) call commdom_localize(CPLNG)  !< 2017Nov19 !< Bilateral
     end if
#endif
    !< AFTER<-never INTO the module !!
    if(current_when==ITASK_AFTER) then                                         !< 2016MAR31
     !print *, "["//trim(title)//"] ==", trim(CNT_SMS)
      if( modul==ID_SOLIDZ ) call Solidz( -1_ip )
      if( modul==ID_NASTIN ) call Nastin( -1_ip )
      if( modul==ID_ALEFOR ) call Alefor( -1_ip )
      if( modul==ID_TEMPER ) call Temper( -1_ip )
    endif
    !
  endif
  !-----------------------------------------------------------------------||---!
#endif
  !-----------------------------------------------------------------------||---!
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------
  !>
  !> @author  M. Rivero
  !> @date
  !> @brief   Localization nodes contact and deform mesh
  !> @details
  !>          It worked for unilateral only implicit time integration scheme
  !>          not for explicit cases. It is particularly designed for bilateral
  !>          cases.
  !-----------------------------------------------------------------------

  subroutine commdom_localize(CPLNG)

    use def_master,           only : displ
    use def_master,           only : ITER_K

    implicit none

    type(COMMDOM_COUPLING), intent(inout) :: CPLNG

#ifdef COMMDOM
    if ( DYNAMIC ) then

       if ( INOTMASTER ) then
          !
          ! Unilateral contact part
          !
          if ( CNT_SENDRECV(7_ip) ) then   ! solamente para el unilateral

             ! Localizacion i body (A1) with deformed shape
             if ( CNT_CPLNG%current_code==CNT_CPLNG%code_i ) then ! Neumann detection (identer)

                coord(1:ndime,1:npoin) = coord(1:ndime,1:npoin) + displ(1:ndime,1:npoin,ITER_K)

             end if

             ! Localization j body (A2) with undeformed shape
             if ( CNT_CPLNG%current_code==CNT_CPLNG%code_j ) then ! Dirichlet detection (block)
                coord(1:ndime,1:npoin) = coord(1:ndime,1:npoin)
             end if

          end if
          !
          ! Neumann contact part
          !
          if ( CNT_SENDRECV(12_ip) ) then

             ! Localization for both bodies
             coord(1:ndime,1:npoin) = coord(1:ndime,1:npoin) + displ(1:ndime,1:npoin,ITER_K)

          end if

       end if

       !
       ! Deallocate and deform mesh
       ! <------ Relax
       if ( CNT_SENDRECV(7_ip ) ) PLEPP_CPLNG%tol   =  1.0e-10_rp      !< localization tolerance update (defined in mod_commdom_plepp.f90)
       if ( CNT_SENDRECV(12_ip) ) PLEPP_CPLNG%tol   =  1.0e-3_rp       !< localization tolerance update (defined in mod_commdom_plepp.f90)
       ! ------> Relax
       call commdom_dynamic_deallocate( PLEPP_CPLNG )
       call commdom_dynamic_set_mesh( CPLNG%current_fixbo, CPLNG%n_dof)

       if ( INOTMASTER ) then
          !
          ! Unilateral contact part
          !
          if ( CNT_SENDRECV(7_ip) ) then

             ! Localizacion i body (A1) with deformed shape
             if ( CNT_CPLNG%current_code==CNT_CPLNG%code_i ) then ! Neumann detection (identer)

                coord(1:ndime,1:npoin) = coord(1:ndime,1:npoin) - displ(1:ndime,1:npoin,ITER_K)

             end if

             ! Localization j body (A2) with undeformed shape
             if ( CNT_CPLNG%current_code==CNT_CPLNG%code_j ) then ! Dirichlet detection (block)

                coord(1:ndime,1:npoin) = coord(1:ndime,1:npoin)

             end if

          end if
          !
          ! Neumann contact part
          !
          if ( CNT_SENDRECV(12_ip) ) then

             ! Localization for both bodies
             coord(1:ndime,1:npoin) = coord(1:ndime,1:npoin) - displ(1:ndime,1:npoin,ITER_K)

          end if

       end if

       if ( RESIDUAL ) then                                                       !< 2015Abr28
          !      call commdom_bvnat(n_dirichlet, n_neumann, init=.false., debug=.true.)  !< 2015May11
       else
          !
          call commdom_driver_nodes_to_reaction( CPLNG )                          !< 2015May10
          !      call commdom_bvnat(n_dirichlet, n_neumann, init=.true., debug=.true.)   !< 2015May11
          residual = .true.
       endif                                                           !< 2015Abr28

    else

       if ( RESIDUAL ) then                                                       !< 2016Feb24
       else
          call commdom_driver_nodes_to_reaction( CPLNG )                          !< moved here from kernel/master/inivar.f90
          residual = .true.
       endif

    endif
#endif

  end subroutine commdom_localize

  !-----------------------------------------------------------------------
  !>
  !> @author  Gerard Guillamet
  !> @date    January, 2019
  !> @brief   Localization nodes contact for Unilateral cases
  !> @details Localization nodes contact for Unilateral cases
  !>
  !-----------------------------------------------------------------------

  subroutine commdom_localize_unilateral(CPLNG)

    use def_master, only : displ
    use def_master, only : ITER_K

    implicit none

    type(COMMDOM_COUPLING), intent(inout) :: CPLNG

#ifdef COMMDOM
    if ( DYNAMIC ) then

       if ( INOTMASTER ) then
          !
          ! Unilateral/Neumann contact partd
          !
          coord(1:ndime,1:npoin) = coord(1:ndime,1:npoin) + displ(1:ndime,1:npoin,ITER_K)

       end if

       ! Dynamic deallocate
       call commdom_dynamic_deallocate( PLEPP_CPLNG )

       ! Deform mesh
       call commdom_dynamic_set_mesh( CPLNG%current_fixbo, CPLNG%n_dof)

       if ( INOTMASTER ) then
          !
          ! Unilateral/Neumann contact parts
          !
          coord(1:ndime,1:npoin) = coord(1:ndime,1:npoin) - displ(1:ndime,1:npoin,ITER_K)

       end if

       !
       ! Residual
       !
       if ( RESIDUAL ) then
          !      call commdom_bvnat(n_dirichlet, n_neumann, init=.false., debug=.true.)  !< 2015May11
       else
          !call commdom_driver_nodes_to_reaction( CPLNG )
          !      call commdom_bvnat(n_dirichlet, n_neumann, init=.true., debug=.true.)   !< 2015May11
          !residual = .true.
       endif

    else

       if ( RESIDUAL ) then
       else
          call commdom_driver_nodes_to_reaction( CPLNG )
          residual = .true.
       endif

    endif
#endif

  end subroutine commdom_localize_unilateral

  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_driver_exchange02( CPLNG, debug)
  implicit none
  type(COMMDOM_COUPLING), intent(inout) :: CPLNG
  logical(ip), optional,  intent(in   ) :: debug
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
#ifdef COMMDOM
  if( present(debug) ) then
    call commdom_dynamic_exchange02( CPLNG%var_ij(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), &
                                     CPLNG%var_ji(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), &
                                     CPLNG%n_dof,                                     &
                                     debug                                            &
                                   )
  else
    call commdom_dynamic_exchange02( CPLNG%var_ij(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), &
                                     CPLNG%var_ji(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), &
                                     CPLNG%n_dof,                                     &
                                     .True.                                           &
                                   )
  endif
#endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_driver_domain_init()
  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_driver_domain_create()
  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_driver_coupling_driver() !current_when,current_task)
  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  subroutine commdom_driver_init( CPLNG )
  use mod_precice,        only: precice_create
  use mod_precice_driver, only: PRCC                                           !< 2016ABR12
  implicit none
  type(COMMDOM_COUPLING), intent(inout) :: CPLNG
  !-----------------------------------------------------------------------||---!
#ifdef COMMDOM
   !!
#if   COMMDOM==1
!!call commdom_alya_init_cht(       CPLNG )
#elif (COMMDOM==2)||(COMMDOM==-2)
  call commdom_driver_init_contact( CPLNG )
#elif COMMDOM==3
  call commdom_driver_init_cht(     CPLNG )
#elif COMMDOM==-3
  call commdom_driver_init_cht_parallelMPMD( CPLNG )
#elif COMMDOM==4
  call commdom_driver_init_fsi(     CPLNG )
#endif
#endif
  !
#ifdef PRECICE
  call precice_create( PRCC )
#endif
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  subroutine commdom_driver_set_mesh( CPLNG )
  use mod_precice,          only: precice_initialize
  use mod_precice_driver,   only: PRCC                                         !< 2016ABR12
  use mod_commdom_alya,     only: commdom_alya_memall
  implicit none
  type(COMMDOM_COUPLING), intent(inout) :: CPLNG
  !-----------------------------------------------------------------------||---!
#ifdef COMMDOM
  !
#if   COMMDOM==1
!!call commdom_plepp_set_mesh(   CPLNG%current_fixbo, CPLNG%n_dof )
#elif (COMMDOM==2)||(COMMDOM==-2)
!!  call commdom_dynamic_set_mesh( CPLNG%current_fixbo, CPLNG%n_dof )   !< 2016JUN01
#elif (COMMDOM==3)||(COMMDOM==-3)
  call commdom_dynamic_set_mesh( CPLNG%current_fixbo, CPLNG%n_dof )
#elif COMMDOM==4
  call commdom_dynamic_set_mesh( CPLNG%current_fixbo, CPLNG%n_dof )
#endif
#if   COMMDOM==1
!!call commdom_cht_memall(  CPLNG )
#elif (COMMDOM==2)||(COMMDOM==-2)
  call commdom_alya_memall( CPLNG )
#elif (COMMDOM==3)||(COMMDOM==-3)
  call commdom_alya_memall( CPLNG )
#elif COMMDOM==4
  call commdom_alya_memall( CPLNG )
#endif
  !
#elif PRECICE
  call precice_initialize( PRCC )
#endif
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!

!-------------------------------------------------------------------------||---!
!----------------------------------------------------------------| PRIVATE |---!
!-------------------------------------------------------------------------||---!

subroutine  commdom_driver_nodes_to_reaction( CPLNG )
  use def_parame,         only: ip, rp
  use def_master,         only: inotmaster
  use def_domain,         only: npoin, mcono, mcodb,  kfl_codno
  use def_kintyp,         only: soltyp
  use mod_commdom_alya,   only: COMMDOM_COUPLING
  use def_master,         only: momod, modul
  use def_coupli,         only: RESIDUAL
  implicit none
  type(COMMDOM_COUPLING), intent(inout) :: CPLNG
  logical(ip) :: sendrecv(5_ip)
  !
  integer(ip) :: ipoin, ndofn, ncono
  logical(ip), pointer :: touched(:) => null()
  logical(ip) :: codno(mcono)
  logical(ip), pointer ::   fixno(:)
  !
  !
  type(soltyp), pointer :: solve_sol(:)
  solve_sol => momod(modul) % solve(1_ip:)
  ndofn = solve_sol(1)%ndofn
#ifdef COMMDOM
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if(inotmaster) then
    allocate( touched(npoin) )
    allocate(   fixno(ndofn) )
  else
    allocate( touched(1) )
    allocate(   fixno(1) )
  endif
  !
  touched     = .false.
  !
  if(inotmaster) then
    do ipoin = 1,npoin
      codno(1:mcono) = kfl_codno(1:mcono,ipoin) /= mcodb+1                     !< Is it at 'boundary'?
      ncono          = count( codno(1:mcono) , KIND=ip)                        !< How many 'codes' are?
      !
      if( ncono>0 ) then                                                       !< Is a                   'boundary'?
        !fixno(1:ndofn) = abs( solve_sol(1)%kfl_fixno(1:ndofn,ipoin) ) == 0    !< Is it fixed? !<2019June17 <GGU>
        if( all(fixno(1:ndofn)) ) touched(ipoin) = .true.
      endif
      !
    enddo
  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  CPLNG%current_what = .false.
  CPLNG%current_what(1_ip) = (CPLNG%what_i==RESIDUAL)
  CPLNG%current_what(2_ip) = (CPLNG%what_j==RESIDUAL)
  CPLNG%current_what(3_ip) = CPLNG%current_what(1_ip).and.CPLNG%current_what(2_ip)
  CPLNG%current_what(4_ip) = CPLNG%current_what(1_ip).or. CPLNG%current_what(2_ip)

  sendrecv = .false.
  sendrecv(1_ip) = (CPLNG%code_i==CPLNG%current_code).and.(CPLNG%module_i==modul).and.(CPLNG%current_what(4_ip)).and.INOTMASTER
  sendrecv(2_ip) = (CPLNG%code_j==CPLNG%current_code).and.(CPLNG%module_j==modul).and.(CPLNG%current_what(4_ip)).and.INOTMASTER
  sendrecv(3_ip) = sendrecv(1_ip).and.sendrecv(2_ip)
  sendrecv(4_ip) = sendrecv(1_ip).or. sendrecv(2_ip)
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( sendrecv(4_ip) ) then                   ! OR
    !
    if( sendrecv(2_ip) ) then                 !  CODE_J
      if( CPLNG%current_what(2_ip) ) then     ! BVNAT_J
        call allocate_react_bvnat(  momod(modul), modul, .false., .true.)
      else                                    ! REACT_J
        call allocate_react_bvnat(  momod(modul), modul, .true., .false.)
        !
        call commdom_plepp_set_source_nodes( solve_sol(1)%lpoin_reaction(1:npoin)  )  ! Mark the nodes where reaction is required   . 2016Feb24 reactived!!
       !solve_sol(1)%lpoin_reaction(1:npoin) = touched(1:npoin)                       ! I do not know why, but this doesnt work!! :(. 2016Feb24
        !
        call commdom_set_lpoin_reaction() !sets to true all the lpoin_reaction nodes which doesn't have fixno > 0

        call allocate_block_system( momod(modul), solve_sol(1)%lpoin_reaction(1:npoin) )
      endif
      if(INOTSLAVE) print *, "[commdom_plepp_inivar]", " 'RESIDUALj'", count( solve_sol(1)%lpoin_reaction(1:npoin) , KIND=ip )
    else&
    if( sendrecv(1_ip) ) then                 !  CODE_I
      if( CPLNG%current_what(1_ip) ) then     ! BVNAT_I
        call allocate_react_bvnat(  momod(modul), modul, .false., .true.)
      else                                    ! REACT_I
        call allocate_react_bvnat(  momod(modul), modul, .true., .false.)
        !
        call commdom_plepp_set_source_nodes( solve_sol(1)%lpoin_reaction(1:npoin)  )  ! Mark the nodes where reaction is required   . 2016Feb24 reactived!!
       !solve_sol(1)%lpoin_reaction(1:npoin) = touched(1:npoin)                       ! I do not know why, but this doesnt work!! :(. 2016Feb24
        !
        call commdom_set_lpoin_reaction() !sets to true all the lpoin_reaction nodes which doesn't have fixno > 0

        call allocate_block_system( momod(modul), solve_sol(1)%lpoin_reaction(1:npoin) )
      endif
      if(INOTSLAVE) print *, "[commdom_plepp_inivar]", " 'RESIDUALi'"
    endif
    !
  else
  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  deallocate( touched )
  deallocate( fixno   )
#endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
end subroutine
  !-------------------------------------------------------------------------||---!
  !                                                                              !
  !------------------------------------------------------------------------------!

  !-----------------------------------------------------------------------||---!
  subroutine commdom_set_lpoin_reaction()
#ifdef COMMDOM
#if   COMMDOM==2
  use def_domain,          only: npoin
  implicit none
  integer(ip)                :: ipoin
  solve => momod(modul) % solve(1:)

  do ipoin = 1,npoin
    if ( .not. maxval(solve(1) % kfl_fixno(:,ipoin)) > 0_ip ) solve(1) % lpoin_reaction(1:npoin) = .true.
  end do
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
#endif
#endif
  end subroutine commdom_set_lpoin_reaction
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_bvnat( n_dirichlet, n_neumann, init, debug )
  use def_domain
  use def_master
  implicit none
  integer(ip), optional, intent(out) :: n_dirichlet(:)
  integer(ip), optional, intent(out) :: n_neumann(:)
  logical(ip), optional, intent(in ) :: init
  logical(ip), optional, intent(in ) :: debug
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  integer(ip) :: idofn, ipoin, itouch
  integer(ip) :: ncono   = -1
  integer(ip) :: ndofn   = -1
  integer(ip) :: n_touch = -1
  integer(ip), pointer :: n_fixno(:) => null()
  integer(ip), pointer :: n_react(:) => null()
  integer(ip), pointer :: n_bvnat(:) => null()
  !
  logical(ip) ::  init_aux = .false.
  logical(ip) :: kfl_react = .false.
  logical(ip) :: kfl_bvnat = .false.
  logical(ip) :: codno(mcono)
  logical(ip), pointer ::   fixno(:)
  logical(ip), pointer :: touched(:) => null()
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  ndofn = solve(1)%ndofn
  idofn = 1
  !
  if(inotmaster) then
    allocate( touched(npoin) )
    allocate(   fixno(ndofn) )
    allocate( n_fixno(ndofn) )
    allocate( n_react(ndofn) )
    allocate( n_bvnat(ndofn) )
  else
    allocate( touched(1) )
    allocate(   fixno(1) )
    allocate( n_fixno(1) )
    allocate( n_react(1) )
    allocate( n_bvnat(1) )
  endif
  !
  touched     = .false.
  !
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------| cach error |---!
  if( present(n_dirichlet).and.( size(n_dirichlet,1)<ndofn) ) then
    print *, "[commdom_bvnat] ", " shape(n_dirichlet)/=ndofn ", shape(n_dirichlet), "/=", ndofn
    call runend("[commdom_bvnat] ERROR!!")
  endif
  !-----------------------------------------------------------| cach error |---!
  if( present(n_neumann).and.( size(n_neumann,1)<ndofn) ) then
    print *, "[commdom_bvnat] ", " shape(n_neumann)/=ndofn ", shape(n_neumann),  "/=", ndofn
    call runend("[commdom_bvnat] ERROR!!")
  endif
  !-----------------------------------------------------------------------||---!
  if( present(n_dirichlet) ) n_dirichlet = 0
  if( present(n_neumann  ) ) n_neumann   = 0
  if( present(init       ) ) init_aux    = init
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
#ifdef COMMDOM
  inotmaster02: &
  if(inotmaster) then
    !---------------------------------------------------------------------||---!
    kfl_react  = solve(1)%kfl_react==1
    kfl_bvnat  = solve(1)%kfl_bvnat==1
    !
    if(kfl_react) then
      !-------------------------------------------------------------------||---!
      call commdom_dynamic_check_fixno(solve(1)%kfl_fixno, 1_ip, 1_ip, .True.) ! fixno, idofn, fixval, ToDo
      call commdom_dynamic_check_fixno(solve(1)%kfl_fixno, 2_ip, 1_ip, .True.) ! fixno, idofn, fixval, ToDo
      call commdom_dynamic_check_fixno(solve(1)%kfl_fixno, 3_ip, 1_ip, .True.) ! fixno, idofn, fixval, ToDo
      !-------------------------------------------------------------------||---!
    endif
    !
    call commdom_plepp_set_source_nodes( touched(1:npoin)  )                   !< Is it 'touched'?
    n_touch = count( touched(1:npoin) , KIND=ip)
    !
    !---------------------------------------------------------------------||---!
    n_fixno(1:ndofn) = 0
    n_react(1:ndofn) = 0
    n_bvnat(1:ndofn) = 0
    do ipoin = 1,npoin
      codno(1:mcono) = kfl_codno(1:mcono,ipoin) /= mcodb+1                     !< Is it at 'boundary'?
      ncono          = count( codno(1:mcono) , KIND=ip)                                 !< How many 'codes' are?
      !
      fixno(1:ndofn) = abs( solve(1)%kfl_fixno(1:ndofn,ipoin) ) /= 0           !< Is it fixed?
      !
      if( ncono>0 ) then                                                       !< Is a                   'boundary'?
        if( touched(ipoin).OR.init_aux ) then                                  !< Is a           'touched boundary'?
          n_fixno(1:ndofn) = n_fixno(1:ndofn) + 1                              !< Is a             'fixed boundary'?
          where(     fixno(1:ndofn) ) n_react(1:ndofn) = n_react(1:ndofn) + 1  !< Is a     'fixed touched boundary'?
          where(.not.fixno(1:ndofn) ) n_bvnat(1:ndofn) = n_bvnat(1:ndofn) + 1  !< Is a not 'fixed touched boundary'?
        endif
      endif
      !
    enddo
    !---------------------------------------------------------| cach error |---!
    !
    ! mca_oob_tcp_msg_recv: readv failed: Connection reset by peer (104) ??
    ! ORTE_ERROR_LOG:  A message is attempting to be sent to... ??
    !
    if( .not.all(n_bvnat(1:ndofn) + n_react(1:ndofn) == n_fixno(1:ndofn)) ) then
      print *, "[commdom_bvnat] ", " 'n_bvnat+n_react /= n_fixno' ", n_bvnat(1:ndofn) + n_react(1:ndofn),"/=", n_fixno(1:ndofn)
      call runend("[commdom_bvnat] ERROR!!")
    endif
    !---------------------------------------------------------------------||---!
    !                                                                          !
    !---------------------------------------------------------------------||---!
    if(kfl_react) then
      !-------------------------------------------------------| cach error |---!
      if( (any(n_bvnat(1:ndofn)>0).and.(.not.init_aux) ) ) then
        itouch = 0
        do ipoin = 1,npoin
          fixno(1:ndofn) = abs( solve(1)%kfl_fixno(1:ndofn,ipoin) ) /= 0
          if( any( .not.fixno(1:ndofn) ).and.( touched(ipoin) ) ) then
            itouch = itouch + 1
            print *, "[commdom_react] ", "codno->", itouch, "|", kfl_codno(1:mcono,ipoin), "|", solve(1)%kfl_fixno(1:ndofn,ipoin)
          endif
        enddo
        !
        print *, "[commdom_react] ", "'"//trim(title)//"'", " 'n_bvnat>0' (", n_bvnat(1:ndofn), ") ", count( solve(1)%lpoin_reaction(1:npoin) , KIND=ip), &
                                    ", 'CHANGE fixno -> 1'"
        call runend("[commdom_react] ERROR!!")
      endif
      !------------------------------------------------------| n_dirichlet |---!
      !
      if(debug) print *, "[kfl_react] ", "'"//trim(title)//"'", n_react(1:ndofn), count( solve(1)%lpoin_reaction(1:npoin) , KIND=ip)
      !
      if( present(n_dirichlet) ) n_dirichlet(1:ndofn) = n_react(1:ndofn)
      !
      !-------------------------------------------------------------------||---!
    endif
    !---------------------------------------------------------------------||---!
    if(kfl_bvnat) then
      !-------------------------------------------------------| cach error |---!
      if( (any(n_react(1:ndofn)>0).and.(.not.init_aux)) ) then
        do ipoin = 1,npoin
          fixno(1:ndofn) = abs( solve(1)%kfl_fixno(1:ndofn,ipoin) ) /= 0
          if( any( fixno(1:ndofn) ).and.( touched(ipoin) ) ) &
            print *, "[commdom_bvnat] ", "codno->", kfl_codno(1:mcono,ipoin), "|", solve(1)%kfl_fixno(1:ndofn,ipoin)
        enddo
        !
        print *, "[commdom_bvnat] ", "'"//trim(title)//"'", " 'n_react>0' ", n_react(1:ndofn)
        call runend("[commdom_bvnat] ERROR!!")
      endif
      !------------------------------------------------------------| error |---!
      if( count(kfl_fixbo(1:nboun)>0,KIND=ip) > 0 ) then
        print *, "[commdom_bvnat] ", "'"//trim(title)//"'", " 'kfl_fixbo>0', unset 'CODES, BOUNDARIES' "
        call runend("[commdom_bvnat] ERROR!!")
      endif
      !--------------------------------------------------------| n_neumann |---!
      !
      if(debug) print *, "[kfl_bvnat] ", "'"//trim(title)//"'", n_bvnat(1:ndofn), count( solve(1)%lpoin_reaction(1:npoin) , KIND=ip)
      !
      if( present(n_neumann) ) n_neumann(1:ndofn) = n_bvnat(1:ndofn)
      !
      !-------------------------------------------------------------------||---!
    endif
    !---------------------------------------------------------------------||---!
  endif inotmaster02
#endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  deallocate( touched )
  deallocate( fixno   )
  deallocate( n_fixno )
  deallocate( n_react )
  deallocate( n_bvnat )
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_driver_begste()
  use def_master,   only: iblok
  use def_coupli,   only: coupling_driver_iteration
  use def_coupli,   only: coupling_driver_number_couplings
  use def_coupli,   only: coupling_driver_max_iteration
  use def_coupli,   only: max_block_cou
  use mod_messages, only: livinf
  !
  use def_master,        only: dtinv, cutim, dtime
  !
  implicit none
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  coupling_driver_iteration(1:max_block_cou) = 0
  coupling_driver_number_couplings(iblok)    = 1
  coupling_driver_max_iteration(iblok)       = n_max_its
  !-----------------------------------------------------------------------||---!
#ifdef COMMDOM
#if COMMDOM==-3
  !
#else
  call commdom_plepp_compare_dtinv(dtinv)                                      !< 2017JAN07
  !                                                                            !< 2017JAN10
  cutim  = cutim - dtime
  call setgts(2_ip)
  call livinf(201_ip, ' ',1_ip)
#endif
#endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_driver_begzon()
  use def_master,    only : iblok
  use def_master,    only : mmodu
  use def_master,    only : lmord
  use def_master,    only : itinn
  use def_coupli,    only : coupling_driver_iteration
  use def_coupli,    only : coupling_driver_number_couplings
  use mod_messages, only : livinf
  implicit none
  integer(ip) :: iorde,imodu
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
     if( coupling_driver_number_couplings(iblok) /= 0 .and. coupling_driver_iteration(iblok) == 0 ) then
        call livinf(-6_ip,'ZONAL COUPLING FOR BLOCK ', iblok)
     end if
     !
     ! Put inner iterations to zero
     !
     do iorde = 1,mmodu
        imodu = lmord(iorde,iblok)
        itinn(imodu) = 0
     end do
     !
     ! Iteration counter
     !
     coupling_driver_iteration(iblok) = coupling_driver_iteration(iblok) + 1
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_driver_endzon()
  use def_master,  only: iblok
  use def_master,  only: kfl_gocou
  use def_coupli,  only: coupling_driver_iteration
  use def_coupli,  only: coupling_driver_number_couplings
  use def_coupli,  only: coupling_driver_max_iteration
  use def_coupli,  only: kfl_gozon
  implicit none
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
     if( coupling_driver_number_couplings(iblok) /= 0 ) then
        !
        kfl_gozon = 0
        if( coupling_driver_iteration(iblok) >= coupling_driver_max_iteration(iblok) ) then
 kfl_gozon = 0 !< kernel/coupli/mod_couplings.f90
!          return
        else
! if( resid_cou(1,icoup) > coupling_driver_tolerance(iblok) ) kfl_gozon = 1 !< kernel/coupli/mod_couplings.f90
          kfl_gozon = 1
        endif
        !
        if( kfl_gozon == 1 ) kfl_gocou = 1
     else
        kfl_gozon = 0
     end if
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_driver_conblk()
  use def_master,   only: iblok
  use def_coupli,   only: coupling_driver_number_couplings
  use def_coupli,   only : coupling_driver_iteration
  use mod_messages, only : livinf
  implicit none
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( coupling_driver_number_couplings(iblok) /= 0 ) then
    call livinf(-13_ip,'END ZONAL COUPLING: ', coupling_driver_iteration(iblok) )
  end if
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  subroutine commdom_driver_init_exchange( current_when, current_task )        !< 2017JAN07
  use def_master,       only: ittim, mitim
  use def_master,       only: ITASK_AFTER,  ITASK_BEFORE
  use def_master,       only: ITASK_DOITER, ITASK_ENDZON, ITASK_TURNOF
  implicit none
  integer(ip),  intent(in)    :: current_when
  integer(ip),  intent(in)    :: current_task
  logical(ip) :: initSend = .false.
  logical(ip) :: initRecv = .false.
  !-----------------------------------------------------------------------||---!
#ifdef COMMDOM
#if COMMDOM==-3
  !
  if(current_task==ITASK_DOITER.and.current_when==ITASK_AFTER) then
    N_DOITER = N_DOITER + 1
  endif
  !
  !-----------------------------------------------------------------------||---!
  if( any(CNT_SENDRECV) ) then
    !
    if( CNT_SENDRECV(5) ) N_DOITER = 0
    !
    if( CNT_CPLNG%current_code==CNT_CPLNG%code_i ) then
      !
      if( CNT_SENDRECV(10) ) then ! SENDi
        if( (ittim==1    ).and.(N_DOITER >0) ) initSend  = .True.
        call commdom_check_parallelMPMD( current_when, current_task, initSend, initRecv )
      endif
      if( CNT_SENDRECV(11) ) then ! RECVi
        !
      endif
    endif
    !
    if( CNT_CPLNG%current_code==CNT_CPLNG%code_j ) then
      !
      if( CNT_SENDRECV(10) ) then ! SENDj
        if( (ittim==mitim).and.(N_DOITER >0) ) INITIATED = .FALSE.
      endif
      if( CNT_SENDRECV(11) ) then ! RECVj
        if( (ittim==1    ).and.(N_DOITER==0) ) initRecv  = .True.
        call commdom_check_parallelMPMD( current_when, current_task, initSend, initRecv )
      endif
      !
      if(current_task==ITASK_TURNOF.and.current_when==ITASK_AFTER) then        !< 2017JAN09
        call par_finali(1_ip)
        if(INOTSLAVE) print *, " '"//trim(title)//"' "//trim(CNT_SMS)
        stop
      endif
      !
    endif
    !---------------------------------------------------------------------||---!
  endif
  !
#endif
#endif
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  subroutine commdom_check_parallelMPMD( current_when, current_task, initSend, initRecv ) !< 2017JAN07
  use def_master,       only: ITASK_AFTER,  ITASK_BEFORE
  use def_master,       only: ITASK_DOITER, ITASK_ENDZON
  use def_master,       only: ittim
  use def_master,       only: title, inotslave
  use mod_communications, only: PAR_BROADCAST
  implicit none
  integer(ip),  intent(in)    :: current_when
  integer(ip),  intent(in)    :: current_task
  logical(ip),  intent(inout) :: initSend
  logical(ip),  intent(inout) :: initRecv
  !-----------------------------------------------------------------------||---!
  !
  real(rp)      :: send, recv
  integer(ip)   :: n_send, n_recv
  character(16) :: msg
  !
  !-----------------------------------------------------------------------||---!
#ifdef COMMDOM
  !
  send   =  real(ittim,rp)
  recv   = -1
  n_send =  1
  n_recv =  1
  !
  !if(INOTSLAVE) print *, trim(CNT_SMS)//"'"//trim(title)//"'", N_DOITER, ittim, initSend, initRecv,  INITIATED, recv
  !
  if( (PLEPP_CPLNG%commij /= -1).or.(PLEPP_CPLNG%commji /= -1) ) then
    if(.not.(PLEPP_CPLNG%commij==-1)) then
      if(initSend) call commdom_sendrecv_real(send, n_send, recv, n_recv, PLEPP_CPLNG%local_comm, PLEPP_CPLNG%commij)
    endif
    if(.not.(PLEPP_CPLNG%commji==-1)) then
      if(initRecv) call commdom_sendrecv_real(send, n_send, recv, n_recv, PLEPP_CPLNG%local_comm, PLEPP_CPLNG%commji)
    endif
  else
    if(INOTSLAVE) print *, trim(CNT_SMS)//"'"//trim(title)//"'", N_DOITER, ittim, initSend, initRecv, INITIATED
    call runend("[commdom_check_parallelMPMD] (PLEPP_CPLNG%commij /= -1).or.(PLEPP_CPLNG%commji /= -1)  EXIT!!")
  endif
  call PAR_BROADCAST(recv,'IN MY CODE')
  !
  if( recv>=0 ) INITIATED = .True.
  !
  if(initSend) msg = "-->"
  if(initRecv) msg = "<--"
  if(INOTSLAVE.and.(initSend.or.initRecv)) print *, "[commdom_check_parallelMPMD] "//"'"//trim(title)//"'", "  ON ", trim(msg), INITIATED
  !
  initSend = .false.
  initRecv = .false.
  !
#endif
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  subroutine commdom_driver_time_exchange( current_when, current_task )        !< 2017JAN07
  use def_master,       only: ITASK_AFTER,  ITASK_BEFORE
  use def_master,       only: ITASK_DOITER, ITASK_ENDZON, ITASK_TURNOF
  implicit none
  integer(ip),  intent(in)    :: current_when
  integer(ip),  intent(in)    :: current_task
  !-----------------------------------------------------------------------||---!
#ifdef COMMDOM
#if COMMDOM==-3
  !-----------------------------------------------------------------------||---!
  if( any(CNT_SENDRECV) ) then
    !
    if( CNT_CPLNG%current_code==CNT_CPLNG%code_i ) then
      !
      if( CNT_SENDRECV(10) ) then ! SENDi
        if(INITIATED) call commdom_change_dt( EXCHANGE=.True., SET=.False.)
      endif
      if( CNT_SENDRECV( 2) ) then ! BEGSTEi
!        if(INITIATED) call commdom_change_dt( EXCHANGE=.True., SET=.False.)
!        if(INITIATED) call commdom_change_dt( EXCHANGE=.False., SET=.True.)
      endif
    endif
    !
    if( CNT_CPLNG%current_code==CNT_CPLNG%code_j ) then
      !
      if( CNT_SENDRECV(10) ) then ! SENDj
!        if(INITIATED) call commdom_change_dt( EXCHANGE=.True., SET=.False.)
     endif
      if( CNT_SENDRECV( 2) ) then ! BEGSTEj
        if(N_DOITER==0) &
call commdom_change_dt( EXCHANGE=.True., SET=.False.)
        if(N_DOITER==0) &
call commdom_change_dt( EXCHANGE=.False., SET=.True.)
      endif
      !
    endif
    !---------------------------------------------------------------------||---!
  endif
  !
#endif
#endif
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  subroutine commdom_change_dt( EXCHANGE, SET )                                                 !< 2017JAN10
  use def_master,         only : dtinv, dtime, cutim, oltim, ittim , dtold
  use mod_couplings,      only : THERE_EXISTS_A_ZONE_COUPLING
  use mod_communications, only : par_min, par_max
  implicit none
  logical(ip), intent(in)  :: EXCHANGE
  logical(ip), intent(in)  :: SET
  !-----------------------------------------------------------------------||---!
  !
#ifdef COMMDOM
  !
  if( EXCHANGE ) then
    !                                                                           !< dt^{ij}_{n} = min( dt^{i}, dt^{j} )
    dt_inv = dtinv

!print *, dtinv
!if(INOTSLAVE) print *, " '"//trim(title)//"' "//trim(CNT_SMS), dtinv, dtime, cutim, oltim, dt_inv, "->"

    call commdom_plepp_compare_dtinv( dt_inv )
! HAY un error muy, pero muy raro en el PRIMER paso de tiempo!!!
if(INOTSLAVE) print *, " '"//trim(title)//"' "//trim(CNT_SMS), dtinv, dtime, cutim, oltim, dt_inv, "<-"

  else &
  if( SET ) then
    if(ittim==1) then
      !                                                                         !< t_{0} = t_{-1}  + dt^{ij}_{n}
!      if( CNT_CPLNG%current_code==CNT_CPLNG%code_i ) &
!        cutim  = 1/dt_inv
!      if( CNT_CPLNG%current_code==CNT_CPLNG%code_j ) &
!        cutim  = 1/dt_inv
!!if(INITIATED) then
!if( CNT_CPLNG%current_code==CNT_CPLNG%code_j ) then

      cutim   = 0.0_rp !oltim - dtime + 1/dt_inv
      dtinv   = dt_inv
      dtime   = 0.0_rp
      oltim   = 0.0_rp
     dtold(:) = 0.0_rp
!endif
    else
      !                                                                         !< t_{n} = t_{n-1} + dt^{ij}_{n}
      dtinv  = dt_inv
      cutim  = cutim - dtime        !
    endif
    call setgts(2_ip)               ! |__ Begste, if(kfl_reset == 1)
    call livinf(201_ip, ' ',1_ip)   !
    !
    if(INOTSLAVE) print *, " '"//trim(title)//"' "//trim(CNT_SMS), dtinv, dtime, cutim, oltim
  else
    call runend("[commdom_change_dt] ERROR dt !!")
  endif
  !
#endif
  !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!

!-------------------------------------------------------------------------||---!
!------------------------------------------------------------------| TOOLs |---!
!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_driver_get_total_flux( dU, U, Ut )
  use def_domain
  use def_master,   only: ittim
  use mod_communications, only: PAR_SUM
  implicit none
  !-----------------------------------------------------------------------||---!
  !
  !     /
  ! U = | dU dA  ->  U = dU_i \delta A_i \Psi_i  ??
  !     /
  !
  !  [U] = units       <- Flux        : the rate of 'dU' transfer through a given surface A, per unit 'surface'
  ! [dU] = units/m^2   <- Flux density: the rate per unit area
  !
  !-----------------------------------------------------------------------||---!
  real(rp), intent(inout)             :: dU(npoin)
  real(rp), intent(inout)             ::  U(npoin)
  real(rp), optional, intent(inout)   ::  Ut(1)
  real(rp)                            ::  Area, dA
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  real(rp)                 :: baloc(ndime,ndime)
  real(rp)                 :: bocod(ndime,mnodb)
  real(rp)                 :: elcod(ndime,mnode)
  real(rp)                 :: xjaci(ndime,ndime)
  real(rp)                 :: xjacm(ndime,ndime)
  real(rp)                 :: gpcar(ndime,mnode,mgaus)
  real(rp)                 :: eucta, detjm
  integer(ip)              :: ielem, inode, ipoin
  integer(ip)              :: igaus, igaub, iboun, inodb, pblty
  integer(ip)              :: pnodb, pnode, pelty, pgaus, pgaub
  type(soltyp), pointer    :: solve_sol(:)
  solve_sol => momod(modul) % solve(1_ip:)
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  !
  gpcar = 0.0_rp
  Area  = 0.0_rp
  !
  if( INOTMASTER ) then
  !-----------------------------------------------------------------------||---!
  !
  U(1:npoin) = 0.0_rp
  !
  boundaries: &
  do iboun = 1,nboun
    !---------------------------------------------------------------------||---!
    pblty = ltypb(iboun)
    pnodb = lnnob(iboun)
    pgaub = ngaus(pblty)
    ielem = lelbo(iboun)
    pelty = ltype(ielem)
    pnode = nnode(pelty)
    pgaus = ngaus(pelty)
    !---------------------------------------------------------------------||---!
    do inodb = 1,pnodb
      ipoin = lnodb(inodb,iboun)
      bocod(1:ndime,inodb) = coord(1:ndime,ipoin)
    end do
    !
    do inode = 1,pnode
      ipoin = lnods(inode,ielem)
      elcod(1:ndime,inode) = coord(1:ndime,ipoin)
    end do
    !
    do igaus = 1,pgaus
      call elmder(pnode, ndime, elmar(pelty)%deriv(1,1,igaus), elcod, gpcar(1,1,igaus), detjm, xjacm, xjaci)
    end do
    !---------------------------------------------------------------------||---!
    gauss_points: &
    do igaub = 1,pgaub
      !
      call bouder(pnodb, ndime, ndimb, elmar(pblty)%deriv(1,1,igaub), bocod, baloc, eucta)
      call chenor(pnode, baloc, bocod, elcod)
      dA   = elmar(pblty)%weigp(igaub) * eucta
      !
      do inodb = 1,pnodb
        ipoin    = lnodb(inodb,iboun)
        U(ipoin) = U(ipoin) + dU(ipoin) * elmar(pblty)%shape(inodb,igaub) * dA
      end do
      !
    end do gauss_points
    !---------------------------------------------------------------------||---!
  end do boundaries
  !-----------------------------------------------------------------------||---!
  if( present(Ut) ) then                                                       !< 2016ABR04
    if( .not.associated(solve_sol(1)%lpoin_reaction) .and. (ittim>0) ) &
      call runend("ERROR: [commdom_driver_get_total_flux] .not.associated(lpoin_reaction) ")
    !
    Ut(1) = 0.0_rp
    if( size(solve_sol(1)%lpoin_reaction(:))==npoin ) Ut(1) = sum( U, mask=solve_sol(1)%lpoin_reaction(:) )
    !
  endif
  !-----------------------------------------------------------------------||---!
  endif
  !call PAR_SUM(1_ip, Ut, 'IN MY CODE')
  !print *, "Ut", Ut
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  end subroutine commdom_driver_get_total_flux
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_driver_get_residual()
  use def_master, only: gesca, gevec
  implicit none
  integer(ip) :: ndofn, nblocks
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if(inotmaster) then
    !
    solve => momod(modul) % solve(1:)
    ndofn   = solve(1) % ndofn
    nblocks = solve(1) % num_blocks
    !
    if( (nblocks>1).and.(solve(1)%kfl_bvnat /= 0) ) then
      print *," [commdom_driver_get_residual] nblocks>1:", nblocks, solve(1)%kfl_bvnat
      print *, " gesca(1:npoin) = solve_sol(1) % block_array(2) % bvnat(1,1:npoin) "
      stop
    endif
    !
!print *, solve(1)%kfl_bvnat, solve(1)%kfl_react, size(solve(1)%reaction(ndofn,:)), size(solve(2)%reaction(ndofn,:))
    !
    if(ndofn==ndime) then
      call memgen(0_ip,ndime,npoin)
      gevec(1:ndime,1:npoin) = 0.0_rp !                                     _____________ !< 2016ABR04
      !                                                                    /
      if(solve(1)%kfl_bvnat == 1 .and.(size(solve(1)%bvnat(ndofn,:))==npoin) ) then
        gevec(1:ndofn,1:npoin) =  solve(1)%bvnat(   1:ndofn,1:npoin)
      endif
      if(solve(1)%kfl_react == 1 .and.(size(solve(1)%reaction(ndofn,:))==npoin) ) then
        gevec(1:ndofn,1:npoin) =  solve(1)%reaction(1:ndofn,1:npoin)
      endif
    else&
    if(ndofn==1_ip) then
      call memgen(0_ip,npoin,0_ip)
      gesca(1:npoin) = 0.0_rp
      !
      if(solve(1)%kfl_bvnat == 1 .and.(size(solve(1)%bvnat(1,:))==npoin) ) then
        gesca(1:npoin) =  solve(1)%bvnat(   1,1:npoin)
      endif
      if(solve(1)%kfl_react == 1 .and.(size(solve(1)%reaction(1,:))==npoin) ) then
        gesca(1:npoin) =  solve(1)%reaction(1,1:npoin)
      endif
    else
      print *, "[commdom_driver_get_residual] ", " 'residue d.o.f='", solve(1) % ndofn
      call runend("[commdom_driver_get_residual] ERROR!!")
    endif
    !
  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine commdom_driver_get_residual
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine  commdom_driver_n_fixno( n_fixno )
    use def_domain,         only: mcono, mcodb,  kfl_codno
    implicit none
    !
    real(rp), intent(out) :: n_fixno(npoin)
    !
    integer(ip)          :: ipoin, ndofn, ncono
    logical(ip)          :: codno(mcono)
    logical(ip)          :: dirich, neumann
    !
    type(soltyp), pointer :: solve_sol(:)
    solve_sol => momod(modul) % solve(1_ip:)
    ndofn = solve_sol(1)%ndofn
    !-----------------------------------------------------------------------||---!
    !                                                                            !
    !-----------------------------------------------------------------------||---!
    if(inotmaster) then
       n_fixno(1:npoin) = 0.0_rp
       do ipoin = 1,npoin
          codno(1:mcono) = kfl_codno(1:mcono,ipoin) /= mcodb+1                     !< Is it at 'boundary'?
          ncono          = count( codno(1:mcono) , KIND=ip)                                 !< How many 'codes' are?
          !
          if( ncono>0 ) then                                                       !< Is a  'boundary'?
             !n_fixno(ipoin) = count( abs(solve_sol(1)%kfl_fixno(1:ndofn,ipoin)) >= 0 , KIND=ip)   !< How many 'boundary' are fixed?
             !
             !< How many 'boundary' are fixed?
             dirich  = all( solve_sol(1)%kfl_fixno(1:ndofn,ipoin) /= 0 )
             neumann = all( solve_sol(1)%kfl_fixno(1:ndofn,ipoin) == 0 )
             n_fixno(ipoin) = real(count( abs(solve_sol(1)%kfl_fixno(1:ndofn,ipoin))>= 0 , KIND=ip),rp)
             if(neumann)  n_fixno(ipoin) = -n_fixno(ipoin)
             !
          endif
          !
       enddo
    endif
    !-----------------------------------------------------------------------||---!
    !                                                                            !
    !-----------------------------------------------------------------------||---!
  end subroutine commdom_driver_n_fixno
  !-----------------------------------------------------------------------||---!

!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!
end module mod_commdom_driver
!==============================================================================!
!==============================================================================!
