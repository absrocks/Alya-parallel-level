!==============================================================================!
!  I am your father...
! 
!< 2014Sep05
!< 2014Sep20, gfortran test
!< 2014Sep23, add 'calc_block'
!< 2014Dic10, from 'ITASK_TIMSTE+ITASK_AFTER' to 'ITASK_BEGSTE+ITASK_BEFORE'
!< 2015Jan29, current_what
!< 2015Feb09, debug: mem_modul<- Subscript #2 of the array MEM_MODUL has value 0 which is less than the lower bound of 1 
!< 2015Jul03, add 'tolerance' and 'tolerance' 
! 
!==============================================================================!
module mod_commdom_alya
  !=================================================================| init |===!
  !-----------------------------------------------------------------------||---!
  !< Sources/kernel/coupli/def_coupli.f90
  !<          o----o----o----o Source
  !<          x------x-------x Target
  !-----------------------------------------------------------------------||---!
  use def_parame,    only: ip, rp
  use def_master,    only: inotmaster, imaster, isequen, islave, inotslave, iparall
  use def_master,    only: kfl_gocou !, mem_modul
  use def_master,    only: ITASK_TURNON, MODUL 
  use def_domain,    only: coord, mnode
  use def_domain,    only: ltype, lnods
  use mod_couplings, only: COU_INTERPOLATE_NODAL_VALUES
  use def_coupli,    only: coupling_type, mcoup, typ_color_coupling
  use mod_ecoute,    only :  ecoute
  !use def_domain,   only: nelem, ndime, npoin, nnode, ngaus
  !use mod_memchk,   only: memchk
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  implicit none 
!#ifndef MPI_OFF
!  include  'mpif.h'
!#endif
  !
  interface
    subroutine func_template(itask)  
      use def_parame, only: ip
      implicit none 
      integer(ip), intent(in) :: itask 
    end subroutine func_template
  end interface 
  !
  type COMMDOM_COUPLING 
    integer(ip) :: n_couplings      =  2_ip
    integer(ip) :: coupling_i       = -1_ip !< COUPLING
    integer(ip) :: coupling_j       = -1_ip !< COUPLING
    !
    integer(ip) :: code_i           = -1_ip !< CODE.   coupling%code_source
    integer(ip) :: code_j           = -1_ip 
    !
    integer(ip) :: module_i         = -1_ip !< MODULE. coupling%module_source
    integer(ip) :: module_j         = -1_ip
    !
    !integer(ip) :: block_i          = -1_ip
    !integer(ip) :: block_j          = -1_ip
    !
    integer(ip) :: send_task_i      = -1_ip !< ITASK_BEGITE|ITASK_TURNON|ITASK_DOITER
    integer(ip) :: send_task_j      = -1_ip
    integer(ip) :: recv_task_i      = -1_ip
    integer(ip) :: recv_task_j      = -1_ip
    !
    integer(ip) :: send_when_i      = -1_ip !< ITASK_BEFORE|ITASK_AFTER
    integer(ip) :: send_when_j      = -1_ip
    integer(ip) :: recv_when_i      = -1_ip
    integer(ip) :: recv_when_j      = -1_ip
    !
    integer(ip) :: send_block_i     = -1_ip !< master/Conblk -> iblok += 1; if(iblok > nblok) kfl_goblk = 0
    integer(ip) :: send_block_j     = -1_ip
    integer(ip) :: recv_block_i     = -1_ip
    integer(ip) :: recv_block_j     = -1_ip
    integer(ip) :: calc_block_i     = -1_ip
    integer(ip) :: calc_block_j     = -1_ip
    !
    integer(ip) :: send_modul_i     = -1_ip
    integer(ip) :: send_modul_j     = -1_ip
    integer(ip) :: recv_modul_i     = -1_ip
    integer(ip) :: recv_modul_j     = -1_ip
    !
    integer(ip) :: fixbo_i          = -1_ip !< WHERE_TYPE: CODE,NUMBER
    integer(ip) :: fixbo_j          = -1_ip
    integer(ip) :: current_fixbo    = -1_ip
    !
    integer(ip) :: what_i             = -1_ip
    integer(ip) :: what_j             = -1_ip
    logical(ip) :: current_what(5_ip) = .false.
    !
    ! <code, block, modul, task, when, send|recv>
    integer(ip) ::      now(     6_ip)  = -1_ip 
    integer(ip) :: the_time(2_ip,6_ip)  = -1_ip
    !
    integer(ip) :: current_code     = -1_ip
    integer(ip) :: current_block    = -1_ip
    integer(ip) :: current_module   = -1_ip
    integer(ip) :: current_task     = -1_ip
    integer(ip) :: current_when     = -1_ip
    integer(ip) :: current_sendrecv = -1_ip
    !
    integer(ip) :: n_dof            = -1_ip
    integer(ip) :: n_pts            = -1_ip
    !
    integer(ip) :: n_wet_pts        = -1_ip
    !
    integer(ip) :: iters_cou        = -1_ip
    integer(ip) :: counter          = -1_ip
    !
    logical(ip) :: sendrecv( 2,7)   = .false. ! 7<-6  
    logical(ip) :: setgetvar(2,3)   = .false.
    !
    integer(ip) :: sendrecv_code    = -1_ip
    integer(ip) :: sendrecv_idx(99) = -1_ip
    integer(ip) :: sendrecv_order   = -1_ip
    !
    integer(ip) :: n_blocks         = -1_ip
    !
    integer(ip) :: n_moduls         = -1_ip
    integer(ip), pointer :: moduls_list(:)
    integer(ip), pointer :: blocks_list(:)
    !
    real(rp), pointer    :: var_ji(:,:)
    real(rp), pointer    :: var_ij(:,:)
    real(rp), pointer    ::  dummy(:,:)
    !
    integer(ip) ::  converged       = 0          !< 2015Jul03
    real(rp)    ::  tolerance(8)    = tiny(1.0)  !< 2015Jul03 
    !
  end type COMMDOM_COUPLING
  ! 
  integer(ip), parameter :: N_CPLNG = 5_ip
  type(COMMDOM_COUPLING), save :: CPLNG(N_CPLNG)
  !
  !
  interface
    subroutine func_template01(itask)  
      use def_parame,       only: ip
      implicit none 
      integer(ip), intent(in) :: itask 
      !type(COMMDOM_COUPLING) :: CPLNG
    end subroutine func_template01
  end interface 
  !procedure(func_template), pointer :: FUNCi01 => null()
  !
  !
  interface commdom_alya_turnon
    module procedure commdom_alya_init
  end interface 
  !
  integer(ip), parameter :: INONE = 0_ip, ISEND = 1_ip, IRECV = 2_ip, ISENDRECV = 3_ip 
  integer(ip), parameter :: SEND_I = 11_ip, RECV_I = 12_ip, SENDRECV_I = 13_ip
  integer(ip), parameter :: SEND_J = 21_ip, RECV_J = 22_ip, SENDRECV_J = 23_ip
  !
  integer(ip), parameter :: CHT = 1
  integer(ip), parameter ::  CC = 4
  !
  private
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
    !
    ! +
    ! |_kernel/defmod/def_master
    !   ID_KERNEL=0, ID_NASTIN=1, ID_TEMPER=2, ID_ALEFOR=7, ID_SOLIDZ=10
    !   ITASK_BEGSTE=5, ITASK_DOITER=6, ITASK_ENDSTE=10
    !   ITASK_BEFORE=1, ITASK_AFTER=2
    !
    ! + 
    ! |_/kernel/master/Concou
    public :: commdom_alya_concou
    !
    ! + 
    ! |_/kernel/coupli/cou_readat
    public :: commdom_alya_turnon
    !
    public :: commdom_alya_memall
    !
public :: commdom_alya_coupling_driver_i
public :: commdom_alya_coupling_driver_j

    public :: COMMDOM_COUPLING, ISEND, IRECV, ISENDRECV

    public :: commdom_alya_exchange

    public :: commdom_alya_set_sendrescv_block_type_extreme
    public :: commdom_alya_sendrecv_driver
    public :: commdom_alya_calculate_driver
    !
    public :: SENDRECV_I, SEND_I, RECV_I
    public :: SENDRECV_J, SEND_J, RECV_J
    public :: INONE
    public :: CPLNG, N_CPLNG, CHT, CC 
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!

  !=======================================================================||===!



  !=============================================================| contains |===!
  contains
!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------| INIT |---!
!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !> @author  JM Zavala-Ake
  !> @date    
  !> @brief   
  !> @details 
  !-----------------------------------------------------------------------||---!
  subroutine commdom_alya_memall(CPLNG)
  use def_domain, only: nelem, ndime, npoin, nnode, ngaus
  use def_master, only: mem_modul
  use def_master, only: ID_NASTIN, ID_SOLIDZ
  use def_master, only: ITASK_BEFORE, ITASK_AFTER
  use def_master, only: current_code, modul
  use mod_memory, only: memory_alloca
  implicit none
  type(COMMDOM_COUPLING), intent(inout) :: CPLNG
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  !
  if(CPLNG%n_dof<=0) then 
    print *, "[commdom_alya_memall] ERROR: n_dof<=0 !!", CPLNG%n_dof
    call runend('EXIT!!')
  endif 
  !
  CPLNG%n_pts          = 0_ip
  if(INOTMASTER) then
!    CPLNG%n_dof        = 1_ip !ndime 
    CPLNG%n_pts        = npoin
  endif   
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!  
  ! se atora!!
  !call memory_alloca( mem_modul(1:2,modul), 'var_ji', 'commdom_memall', CPLNG%var_ji, CPLNG%n_dof, CPLNG%n_pts)
  !call memory_alloca( mem_modul(1:2,modul), 'var_ij', 'commdom_memall', CPLNG%var_ij, CPLNG%n_dof, CPLNG%n_pts)
  !call memory_alloca( mem_modul(1:2,modul),  'dummy', 'commdom_memall',  CPLNG%dummy, CPLNG%n_dof, CPLNG%n_pts)
  !
  if( .not.associated(CPLNG%var_ij) ) then
    allocate( CPLNG%var_ij(CPLNG%n_dof,CPLNG%n_pts) )
  else
     call runend('[commdom_alya_memall] ERROR: code_i==code_j!!')
  endif
  if( .not.associated(CPLNG%var_ji) ) allocate( CPLNG%var_ji(CPLNG%n_dof,CPLNG%n_pts) )
  if( .not.associated(CPLNG%dummy)  ) allocate( CPLNG%dummy( CPLNG%n_dof,CPLNG%n_pts) )
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!  
  CPLNG%var_ji(1:CPLNG%n_dof,1:CPLNG%n_pts) =    0.0_rp
  CPLNG%var_ij(1:CPLNG%n_dof,1:CPLNG%n_pts) =    0.0_rp
  CPLNG%dummy( 1:CPLNG%n_dof,1:CPLNG%n_pts) = -666.666_rp
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if(IMASTER) print*, "[commdom_alya_memall]"
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !> @author  JM Zavala-Ake
  !> @date    
  !> @brief   
  !> @details 
  !-----------------------------------------------------------------------||---!
  subroutine commdom_alya_set_sendrescv_block_type_extreme(CPLNG)
  use def_coupli,           only: mcoup
  use def_master,           only: lmord, iblok, mmodu, mblok, nblok, ITASK_ENDSTE, ITASK_BEGSTE, ITASK_INIUNK, ITASK_AFTER
  use mod_memory, only: memory_alloca
  use def_master, only: mem_modul
  implicit none
  type(COMMDOM_COUPLING), intent(inout) :: CPLNG
  !
  integer(ip)   :: i_block, i_order, i_modul
  integer(ip)   :: first_modul=-1, last_modul=-1
  integer(ip)   :: first_block=-1, last_block=-1
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  !
  !                  __________ITASK_BEGSTE|ITASK_BEFORE <- recv 
  !                 /         _ITASK_ENDSTE|ITASK_AFTER  <- send 
  ! BLOCK 3   \    /         / 
  !   1 X     |   | [1 2 2 3]   |--current_block  -> CPLNG%blocks_list 
  !   2 Y Z   |-->|                        
  !   3 W     |   | [X Y Z W]   |--current_module -> CPLNG%moduls_list 
  ! END_BLOCK /    | |  |  |  |
  !                | |  |  |  |-END TIME STEP
  !                | |  |  |-END BLOCK
  !                | |  |  |-START BLOCK 3 
  !                | |  |-END BLOCK
  !                | |  |-START BLOCK 2 
  !                | |-END BLOCK
  !                | |-START BLOCK 1 
  !                |- START TIME STEP
  !
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!  
  call memory_alloca( mem_modul(1:2,modul), 'moduls_list', 'commdom_memall', CPLNG%moduls_list, mmodu*mblok) !< 2015Feb09
  call memory_alloca( mem_modul(1:2,modul), 'blocks_list', 'commdom_memall', CPLNG%blocks_list, mmodu*mblok) !< 2015Feb09
 !allocate( CPLNG%moduls_list(mmodu*mblok)  )
 !allocate( CPLNG%blocks_list(mmodu*mblok)  )

  CPLNG%moduls_list = -1
  CPLNG%blocks_list = -1
  !
  CPLNG%n_moduls = 0
  do i_block = 1,mblok
    do i_order = 1,mmodu
       i_modul = lmord(i_order,i_block)
       if(i_modul>0) then
          CPLNG%n_moduls = CPLNG%n_moduls + 1
          CPLNG%moduls_list(CPLNG%n_moduls) = i_modul
          CPLNG%blocks_list(CPLNG%n_moduls) = i_block
          if((CPLNG%code_i==CPLNG%current_code).and.(CPLNG%module_i==i_modul)) CPLNG%calc_block_i=i_block
          if((CPLNG%code_j==CPLNG%current_code).and.(CPLNG%module_j==i_modul)) CPLNG%calc_block_j=i_block
       endif
    enddo
  enddo
  !
  if(nblok==1) then
    if(CPLNG%code_i==CPLNG%current_code) CPLNG%calc_block_i = 1
    if(CPLNG%code_j==CPLNG%current_code) CPLNG%calc_block_j = 1
  endif
  !
  first_modul = CPLNG%moduls_list(1)
   last_modul = CPLNG%moduls_list(CPLNG%n_moduls)
  first_block = CPLNG%blocks_list(1)
   last_block = CPLNG%blocks_list(CPLNG%n_moduls)
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if(CPLNG%code_i==CPLNG%code_j) call runend('ERROR: code_i==code_j!!')
  !
  if((CPLNG%code_i==CPLNG%current_code).and.(CPLNG%code_j/=CPLNG%current_code)) then
    CPLNG%recv_modul_i   = first_modul
    CPLNG%send_modul_i   = last_modul
    !
    CPLNG%recv_block_i   = first_block
    CPLNG%send_block_i   = last_block
    ! 
    if(CPLNG%send_task_i==ITASK_ENDSTE) CPLNG%send_block_i = nblok+1 !< Conblk-> iblok = iblok+1   
    !
    if(CPLNG%recv_task_i==CPLNG%send_task_i) then
      if(CPLNG%recv_task_i==ITASK_BEGSTE) then
        CPLNG%send_modul_i = CPLNG%recv_modul_i
        CPLNG%send_block_i = CPLNG%recv_block_i
      endif
      if(CPLNG%recv_task_i==ITASK_ENDSTE) then
        CPLNG%recv_modul_i = CPLNG%send_modul_i
        CPLNG%recv_block_i = CPLNG%send_block_i
      endif
    endif
  endif
  if((CPLNG%code_j==CPLNG%current_code).and.(CPLNG%code_i/=CPLNG%current_code)) then
    CPLNG%recv_modul_j   = first_modul
    CPLNG%send_modul_j   = last_modul
    !
    CPLNG%recv_block_j   = first_block
    CPLNG%send_block_j   = last_block
    !
    if(CPLNG%send_task_j==ITASK_ENDSTE) CPLNG%send_block_j = nblok+1 !< Conblk-> iblok = iblok+1  
    !
    if(CPLNG%recv_task_j==CPLNG%send_task_j) then
      if(CPLNG%recv_task_j==ITASK_BEGSTE) then
        CPLNG%send_modul_j = CPLNG%recv_modul_j
        CPLNG%send_block_j = CPLNG%recv_block_j
        CPLNG%send_when_j  = CPLNG%recv_when_j
      endif
      if(CPLNG%recv_task_j==ITASK_ENDSTE) then
        CPLNG%recv_modul_j = CPLNG%send_modul_j
        CPLNG%recv_block_j = CPLNG%send_block_j
        CPLNG%recv_when_j  = CPLNG%send_when_j
      endif
    endif
  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !> @author  JM Zavala-Ake
  !> @date    
  !> @brief   
  !> @details 
  !-----------------------------------------------------------------------||---!
  subroutine commdom_alya_sendrecv_driver(CPLNG, current_when, current_task)
  use def_master,           only: iblok, nblok, ITASK_INIUNK, ITASK_AFTER, ITASK_TIMSTE, ITASK_BEFORE 
  use def_master,           only: ID_KERMOD, ITASK_BEGSTE
  implicit none
  type(COMMDOM_COUPLING), intent(inout) :: CPLNG
  integer(ip),  intent(in)  :: current_when
  integer(ip),  intent(in)  :: current_task
  integer(ip) :: n_sendrecv
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  !   + current_code                                      ___________current_task 
  !   |_Alya                                       ______|_____
  !     |_call Turnon()                            ITASK_TURNON 02  
  !     |_call Iniunk()                            ITASK_INIUNK 03 
  !     |_time: do while
  !       |_call Timste()                          ITASK_TIMSTE 04 
  !       |_do 
  !       | |_call Begste()                        ITASK_BEGSTE 05 
  !       |    |_block: do while                                     _
  !       |       |_coupling_modules: do while                      / TASK_BEGITE  14 
  !       |       | |_call Doiter()                ITASK_DOITER 06-|  
  !       |       | |_call Concou()                ITASK_CONCOU 07  \_ITASK_ENDITE 15 
  !       |       |_call Conblk()                  ITASK_CONBLK 08 
  !       |       |_call Newmsh()                  ITASK_NEWMSH 09 
  !       |_call Endste()                          ITASK_ENDSTE 10 
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
  ! CPLNG%sendrecv([1,2,3], 1) <- send
  ! CPLNG%sendrecv([1,2,3], 2) <- recv
  ! CPLNG%sendrecv([1,2,3], 3) <- send.and.recv
  ! CPLNG%sendrecv([1,2,3], 4) <- send.or.recv
  ! CPLNG%sendrecv([1,2,3], 5) <- ITASK_AFTER|ITASK_INIUNK|nblok|modul 
  ! CPLNG%sendrecv([1,2,3], 6) <- ITASK_AFTER|ITASK_TIMSTE|nblok+1|modul Begste.iniste.iblok=1 
  !
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  !CPLNG%current_code    <-- commdom_alya_init_XXXXX
  CPLNG%current_block    =  iblok 
  CPLNG%current_module   =  modul
  CPLNG%current_task     =  current_task
  CPLNG%current_when     =  current_when
  CPLNG%current_sendrecv = -1_ip 
  !
  CPLNG%sendrecv = .false. 
  !
  CPLNG%now = (/ CPLNG%current_code, iblok, modul, current_task, current_when, -1_ip /) 
  !
  !--------------------------------------------------| [send|recv]_block_i |---!
  CPLNG%sendrecv(1,1) = (CPLNG%send_task_i ==current_task).and.(CPLNG%send_when_i ==current_when).and.&
                        (CPLNG%send_block_i==iblok       ).and.(CPLNG%send_modul_i==       modul)
  CPLNG%sendrecv(1,2) = (CPLNG%recv_task_i ==current_task).and.(CPLNG%recv_when_i ==current_when).and.&
                        (CPLNG%recv_block_i==iblok       ).and.(CPLNG%recv_modul_i==       modul)
  CPLNG%sendrecv(1,3) =  CPLNG%sendrecv(1,1).and.CPLNG%sendrecv(1,2)
  CPLNG%sendrecv(1,4) =  CPLNG%sendrecv(1,1).or. CPLNG%sendrecv(1,2)
  !--------------------------------------------------| [send|recv]_block_j |---!
  CPLNG%sendrecv(2,1) = (CPLNG%send_task_j ==current_task).and.(CPLNG%send_when_j ==current_when).and.&
                        (CPLNG%send_block_j==       iblok).and.(CPLNG%send_modul_j==       modul)
  CPLNG%sendrecv(2,2) = (CPLNG%recv_task_j ==current_task).and.(CPLNG%recv_when_j ==current_when).and.&
                        (CPLNG%recv_block_j==       iblok).and.(CPLNG%recv_modul_j==       modul)
  CPLNG%sendrecv(2,3) =  CPLNG%sendrecv(2,1).and.CPLNG%sendrecv(2,2)
  CPLNG%sendrecv(2,4) =  CPLNG%sendrecv(2,1).or. CPLNG%sendrecv(2,2)
  !---------------------------------| ITASK_AFTER|ITASK_INIUNK|nblok|modul |---!
  CPLNG%sendrecv(1,5) = (current_task==ITASK_INIUNK).and.(      current_when==ITASK_AFTER).and.&
                        (       iblok==nblok       ).and.(CPLNG%send_modul_i==modul      )
  CPLNG%sendrecv(2,5) = (current_task==ITASK_INIUNK).and.(      current_when==ITASK_AFTER).and.&
                        (       iblok==nblok       ).and.(CPLNG%send_modul_j==modul)
  !---------------------------------| ITASK_AFTER|ITASK_TIMSTE|nblok|modul |---!
  ! |_/master/Begste
  !   |_iniste
  !     |_iblok=1 
!
!  CPLNG%sendrecv(1,6) = (current_task==ITASK_TIMSTE).and.(      current_when==ITASK_AFTER).and.&
!                        (       iblok==nblok+1     ).and.(CPLNG%send_modul_i==modul      )
!  CPLNG%sendrecv(2,6) = (current_task==ITASK_TIMSTE).and.(      current_when==ITASK_AFTER).and.&
!                        (       iblok==nblok+1     ).and.(CPLNG%send_modul_j==modul) 
!
  ! |_Alya
  !   |_Timste 
  !   | |_ moduls(ITASK_TIMSTE)
  !   | |_ setgts(2_ip)        <-- minimum of critical time steps. NUNCA PASA POR MODULS!!
  !   |_Begste                 <-- aqui ya se tiene el minimo general de todos los modulos!!
  CPLNG%sendrecv(1,6) = (current_task==ITASK_BEGSTE).and.(      current_when==ITASK_BEFORE).and.&
                        (       iblok==1           ).and.(CPLNG%send_modul_i==modul       )
  CPLNG%sendrecv(2,6) = (current_task==ITASK_BEGSTE).and.(      current_when==ITASK_BEFORE).and.&
                        (       iblok==1           ).and.(CPLNG%send_modul_j==modul)

  !----------------------------------------------------------|TIMSTE.BEFORE|---!
  !
  CPLNG%the_time(1,:) = (/ CPLNG%current_code, nblok+1, CPLNG%send_modul_i, ITASK_TIMSTE, ITASK_BEFORE, -1_ip /)
  CPLNG%the_time(2,:) = (/ CPLNG%current_code, nblok+1, CPLNG%send_modul_j, ITASK_TIMSTE, ITASK_BEFORE, -1_ip /)
  CPLNG%sendrecv(1,7) = all( CPLNG%now(:) == CPLNG%the_time(1,:) ) 
  CPLNG%sendrecv(2,7) = all( CPLNG%now(:) == CPLNG%the_time(2,:) ) 
  !
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  n_sendrecv = 0
  CPLNG%sendrecv_code = 0
  !
  if( CPLNG%sendrecv(1,4) ) then     ! i 
      if( CPLNG%sendrecv(1,3) ) then ! sendrecv 
        !if(IMASTER.or.ISEQUEN) print *, "sendrecv01_i",  CPLNG%current_code, iblok, modul, current_task, current_when
        CPLNG%sendrecv_code = SENDRECV_I
        n_sendrecv = n_sendrecv + 1
      else&
      if( CPLNG%sendrecv(1,2) ) then ! recv 
        !if(IMASTER.or.ISEQUEN) print *, "recv01_i", CPLNG%current_code, iblok, modul, current_task, current_when
        CPLNG%sendrecv_code = RECV_I
        n_sendrecv = n_sendrecv + 1
      else&
      if( CPLNG%sendrecv(1,1) ) then ! send 
        !if(IMASTER.or.ISEQUEN) print *, "send01_i", CPLNG%current_code, iblok, modul, current_task, current_when
        CPLNG%sendrecv_code = SEND_I
        n_sendrecv = n_sendrecv + 1
      endif
  else &
  if( CPLNG%sendrecv(2,4) ) then     ! j 
      if( CPLNG%sendrecv(2,3) ) then ! sendrecv 
        !if(IMASTER.or.ISEQUEN) print *, "sendrecv01_j",  CPLNG%current_code, iblok, modul, current_task, current_when
        CPLNG%sendrecv_code = SENDRECV_J
        n_sendrecv = n_sendrecv + 1
      else&
      if( CPLNG%sendrecv(2,2) ) then ! recv 
        !if(IMASTER.or.ISEQUEN) print *, "recv01_j",  CPLNG%current_code, iblok, modul, current_task, current_when
        CPLNG%sendrecv_code = RECV_J
        n_sendrecv = n_sendrecv + 1
      else&
      if( CPLNG%sendrecv(2,1) ) then ! send 
        !if(IMASTER.or.ISEQUEN) print *, "send01_j",  CPLNG%current_code, iblok, modul, current_task, current_when
        CPLNG%sendrecv_code = SEND_J
        n_sendrecv = n_sendrecv + 1
      endif
  else &
  if( CPLNG%sendrecv(1,5) ) then     ! init_sendrecv_i
        !if(IMASTER.or.ISEQUEN) print *, "sendrecv00_i",  CPLNG%current_code, iblok, modul, current_task, current_when
        CPLNG%sendrecv_code = -SENDRECV_I
        n_sendrecv = n_sendrecv + 1
  else &
  if( CPLNG%sendrecv(2,5) ) then     ! init_sendrecv_j 
        !if(IMASTER.or.ISEQUEN) print *, "sendrecv00_j",  CPLNG%current_code, iblok, modul, current_task, current_when
        CPLNG%sendrecv_code = -SENDRECV_J
        n_sendrecv = n_sendrecv + 1
  else &
  if( CPLNG%sendrecv(1,6) ) then 
        !if(IMASTER.or.ISEQUEN) print *, "sendrecv00_i",  CPLNG%current_code, iblok, modul, current_task, current_when
        !CPLNG%sendrecv_code = -SENDRECV_I
        n_sendrecv = n_sendrecv + 1
  else &
  if( CPLNG%sendrecv(2,6) ) then 
        !if(IMASTER.or.ISEQUEN) print *, "sendrecv00_j",  CPLNG%current_code, iblok, modul, current_task, current_when
        !CPLNG%sendrecv_code = -SENDRECV_J
        n_sendrecv = n_sendrecv + 1
  else
      !call runend('driver: not coded')
        n_sendrecv = -1
  endif
  !
  if( n_sendrecv>1 ) then
    print *, CPLNG%sendrecv(1,:)
    print *, CPLNG%sendrecv(2,:)
    print *, "n_sendrecv:", n_sendrecv
    call runend('ERROR: send|recv options must be only 1 !!')
  else
  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !> @author  JM Zavala-Ake
  !> @date    
  !> @brief   
  !> @details 
  !-----------------------------------------------------------------------||---!
  subroutine commdom_alya_calculate_driver(CPLNG, restriction) 
  use def_master,           only: ITASK_CONCOU, ITASK_BEGSTE, ITASK_DOITER, ITASK_INIUNK
  use def_master,           only: itcou, micou, iblok
  implicit none
  type(COMMDOM_COUPLING), intent(inout) :: CPLNG
  logical(ip),  optional, intent(in)    :: restriction(2_ip)
  !integer(ip) :: n_sendrecv
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  logical(ip) :: is_when, is_block, is_modul, is_task, is_code
  logical(ip) :: set_var_ij, get_var_ji, ini_var_ij
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !----------------------------------------------------------| CALC_MODULE |---!
  CPLNG%setgetvar(:,:) = .false.
  set_var_ij = .false.
  get_var_ji = .false.
  ini_var_ij = .false. 
  !
  is_code    = CPLNG%current_code   == CPLNG%code_i
  is_modul   = CPLNG%current_module == CPLNG%module_i 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !---------------------------------------------------------| ITASK_CONCOU |---!
  ! <code, block, modul, task, when, send|recv>
  is_block   = CPLNG%current_block == CPLNG%calc_block_i !iblok 
  is_task    = CPLNG%current_task  == ITASK_CONCOU
  is_when    = itcou==micou(iblok)
  set_var_ij = is_code.and.is_block.and.is_modul.and.is_when.and.is_task
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !---------------------------------------------------------| ITASK_BEGSTE |---!
  ! <code, block, modul, task, when, send|recv>
  is_block   = CPLNG%current_block == 1_ip
  is_task    = CPLNG%current_task  == ITASK_BEGSTE
  is_when    = .true.
  get_var_ji = is_code.and.is_block.and.is_modul.and.is_when.and.is_task
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !---------------------------------------------------------| ITASK_BEGSTE |---!
  ! <code, block, modul, task, when, send|recv>
  is_block   = CPLNG%current_block == CPLNG%calc_block_i !iblok 
  is_task    = CPLNG%current_task  == ITASK_INIUNK
  is_when    = .true.
  ini_var_ij = is_code.and.is_block.and.is_modul.and.is_when.and.is_task
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( present(restriction) ) then
    ini_var_ij = ini_var_ij.and.restriction(CPLNG%code_i)
    set_var_ij = set_var_ij.and.restriction(CPLNG%code_i)
    get_var_ji = get_var_ji.and.restriction(CPLNG%code_i)
  endif
  if(IMASTER.or.ISEQUEN) then
    if(ini_var_ij) print *, "[commdom_plepp_code_i]", " ini_var_ij" !, "-->", CPLNG%setgetvar
    if(set_var_ij) print *, "[commdom_plepp_code_i]", " set_var_ij" !, "-->", CPLNG%setgetvar
    if(get_var_ji) print *, "[commdom_plepp_code_i]", " get_var_ji" !, "-->", CPLNG%setgetvar
  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  CPLNG%setgetvar(CPLNG%code_i,1) = ini_var_ij
  CPLNG%setgetvar(CPLNG%code_i,2) = set_var_ij
  CPLNG%setgetvar(CPLNG%code_i,3) = get_var_ji
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  set_var_ij = .false.
  get_var_ji = .false.
  ini_var_ij = .false.
  !
  is_code    = CPLNG%current_code   == CPLNG%code_j
  is_modul   = CPLNG%current_module == CPLNG%module_j
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !---------------------------------------------------------| ITASK_CONCOU |---!
  ! <code, block, modul, task, when, send|recv>
  is_block   = CPLNG%current_block == CPLNG%calc_block_j !iblok 
  is_task    = CPLNG%current_task  == ITASK_CONCOU
  is_when    = itcou==micou(iblok)
  set_var_ij = is_code.and.is_block.and.is_modul.and.is_when.and.is_task
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !---------------------------------------------------------| ITASK_BEGSTE |---!
  ! <code, block, modul, task, when, send|recv>
  is_block   = CPLNG%current_block == 1_ip
  is_task    = CPLNG%current_task  == ITASK_BEGSTE
  is_when    = .true.
  get_var_ji = is_code.and.is_block.and.is_modul.and.is_when.and.is_task
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !---------------------------------------------------------| ITASK_BEGSTE |---!
  ! <code, block, modul, task, when, send|recv>
  is_block   = CPLNG%current_block == CPLNG%calc_block_j !iblok 
  is_task    = CPLNG%current_task  == ITASK_INIUNK
  is_when    = .true.
  ini_var_ij = is_code.and.is_block.and.is_modul.and.is_when.and.is_task
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( present(restriction) ) then
    ini_var_ij = ini_var_ij.and.restriction(CPLNG%code_j)
    set_var_ij = set_var_ij.and.restriction(CPLNG%code_j)
    get_var_ji = get_var_ji.and.restriction(CPLNG%code_j)
  endif
  if(IMASTER.or.ISEQUEN) then
    if(ini_var_ij) print *, "[commdom_plepp_code_j]", " ini_var_ij"
    if(set_var_ij) print *, "[commdom_plepp_code_j]", " set_var_ij"
    if(get_var_ji) print *, "[commdom_plepp_code_j]", " get_var_ji"
  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  CPLNG%setgetvar(CPLNG%code_j,1) = ini_var_ij
  CPLNG%setgetvar(CPLNG%code_j,2) = set_var_ij
  CPLNG%setgetvar(CPLNG%code_j,3) = get_var_ji
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


!-------------------------------------------------------------------------||---!
!-----------------------------------------------------------------| DRIVER |---!
!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !> @author  JM Zavala-Ake
  !> @date    
  !> @brief   
  !> @details 
  !-----------------------------------------------------------------------||---!
  subroutine commdom_alya_coupling_driver_i(               &
                                            CPLNG,         &
                                            sendrecv_code, &
                                            current_when,  &
                                            f_before_ini,  &
                                            f_before_end,  &
                                            f_after_ini,   &
                                            f_after_end,   &
                                            dummy          &
                                            )
  use def_master, only: ITASK_AFTER, ITASK_BEFORE
  implicit none
  type(COMMDOM_COUPLING), intent(inout) :: CPLNG
  integer(ip), intent(in)               :: sendrecv_code
  integer(ip), intent(in),  optional    :: current_when
  procedure(func_template), optional    :: f_before_ini
  procedure(func_template), optional    :: f_before_end
  procedure(func_template), optional    :: f_after_ini
  procedure(func_template), optional    :: f_after_end
  integer(ip)             , optional    :: dummy
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  if(CPLNG%module_i /= CPLNG%current_module) then 
    print *, "module_i /= current_module", CPLNG%module_i, CPLNG%current_module
    call runend('[commdom_alya_coupling_driver_i] module_i /= current_module!!')
  endif 
  !
  if( present(current_when) ) then
    !
    select case(current_when)
      case(ITASK_BEFORE)
        !
        if( present(f_before_ini) ) call f_before_ini(dummy) 
        select case(sendrecv_code) 
          case(ISEND)
            call commdom_alya_exchange(CPLNG, CPLNG%module_j, ISEND) 
          case(IRECV)
            call commdom_alya_exchange(CPLNG, CPLNG%module_i, IRECV)   
        end select
        if( present(f_before_end) ) call f_before_end(dummy) 
        ! 
      case(ITASK_AFTER) 
        !
        if( present(f_after_ini) ) call f_after_ini(dummy) 
        select case(sendrecv_code) 
          case(ISEND)
            call commdom_alya_exchange(CPLNG, CPLNG%module_j, ISEND) 
          case(IRECV)
            call commdom_alya_exchange(CPLNG, CPLNG%module_i, IRECV)   
if(INOTMASTER) print *, "var_ij02:", sum( CPLNG%var_ji(1_ip,1:CPLNG%n_pts) ) /CPLNG%n_pts

        end select
        if( present(f_after_end) ) call f_after_end(dummy) 
        ! 
      case default
      ! 
    end select
    !
  else  
    !
    select case(sendrecv_code) 
      case(ISEND)
        call commdom_alya_exchange(CPLNG, CPLNG%module_j, ISEND) 
      case(IRECV)
        call commdom_alya_exchange(CPLNG, CPLNG%module_i, IRECV)   
    end select 
    !
  endif 
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !> @author  JM Zavala-Ake
  !> @date    
  !> @brief   
  !> @details 
  !-----------------------------------------------------------------------||---!
  subroutine commdom_alya_coupling_driver_j(               &
                                            CPLNG,         &
                                            sendrecv_code, &
                                            current_when,  &
                                            f_before_ini,  &
                                            f_before_end,  &
                                            f_after_ini,   &
                                            f_after_end,   &
                                            dummy          &
                                            )
  use def_master, only: ITASK_AFTER, ITASK_BEFORE
  implicit none
  type(COMMDOM_COUPLING), intent(inout) :: CPLNG
  integer(ip), intent(in)               :: sendrecv_code
  integer(ip), intent(in),  optional    :: current_when
  procedure(func_template), optional    :: f_before_ini
  procedure(func_template), optional    :: f_before_end
  procedure(func_template), optional    :: f_after_ini
  procedure(func_template), optional    :: f_after_end
  integer(ip)             , optional    :: dummy
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  if(CPLNG%module_j /= CPLNG%current_module) then 
    print *, "module_j /= current_module", CPLNG%module_i, CPLNG%current_module
    call runend('[commdom_alya_coupling_driver_j] module_j /= current_module!!')
  endif 
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  if( present(current_when) ) then
    !
    select case(current_when)
      case(ITASK_BEFORE)
        !
        if( present(f_before_ini) ) call f_before_ini(dummy) 
        select case(sendrecv_code) 
          case(ISEND)
            call commdom_alya_exchange(CPLNG, CPLNG%module_i, ISEND) 
          case(IRECV)
            call commdom_alya_exchange(CPLNG, CPLNG%module_j, IRECV)   
        end select
        if( present(f_before_end) ) call f_before_end(dummy) 
        ! 
      case(ITASK_AFTER) 
        !
        if( present(f_after_ini) ) call f_after_ini(dummy) 
        select case(sendrecv_code) 
          case(ISEND)
!if(INOTMASTER) print *, "var_ij03:", sum( CPLNG%var_ij(1_ip,1:CPLNG%n_pts) ) /CPLNG%n_pts, CPLNG%n_pts
!print *, "var_ji01:", sum( CPLNG%var_ji(1_ip,1:CPLNG%n_pts) ) /CPLNG%n_pts
            call commdom_alya_exchange(CPLNG, CPLNG%module_i, ISEND) 
          case(IRECV)
            call commdom_alya_exchange(CPLNG, CPLNG%module_j, IRECV)   
        end select
        if( present(f_after_end) ) call f_after_end(dummy) 
        ! 
      case default
      ! 
    end select
    !
  else  
    !
    select case(sendrecv_code) 
      case(ISEND)
        call commdom_alya_exchange(CPLNG, CPLNG%module_i, ISEND) 
      case(IRECV)
        call commdom_alya_exchange(CPLNG, CPLNG%module_j, IRECV)   
    end select 
    !
  endif 
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


!-------------------------------------------------------------------------||---!
!---------------------------------------------------------------| SENDRECV |---!
!-------------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !> @author  JM Zavala-Ake
  !> @date    
  !> @brief   
  !> @details 
  !-----------------------------------------------------------------------||---!
  subroutine commdom_alya_exchange(CPLNG, WhoIam, sendrecv_code)
  use def_master, only: current_code
  implicit none
  type(COMMDOM_COUPLING), intent(inout) :: CPLNG
  integer(ip),  intent(in)  :: WhoIam, sendrecv_code
  integer(ip)   :: WhoID = -1
  character(16) :: done='', who, str(3)
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if(WhoIam == CPLNG%module_i) then 
     WhoID = CPLNG%code_i
     who = '_i'
  else& 
  if(WhoIam == CPLNG%module_j) then 
     WhoID = CPLNG%code_j
     who = '_j'
  else
     print * , "[commdom_alya_exchange] ERROR: module name", WhoIam  
     call runend('EXIT!!')
  endif 
  !-----------------------------------------------------------------------||---!  
  !                                                                            !
  !-----------------------------------------------------------------------||---!  
  select case(sendrecv_code) 
    case(ISEND)
      done = '_send'
      CPLNG%counter = CPLNG%counter + 1
    case(IRECV)
      done = '_recv'
      CPLNG%counter = CPLNG%counter + 1
    case(ISENDRECV)
      print * , "[commdom_alya_exchange] ERROR: USE 'ISEND'|'IRECV' "
      call runend('EXIT!!')
  end select 
  !-----------------------------------------------------------------------||---!  
  !                                                                            !
  !-----------------------------------------------------------------------||---!  
  select case(sendrecv_code) 
    case(ISEND)
      call COU_INTERPOLATE_NODAL_VALUES(WhoID, &
                                        CPLNG%n_dof,    &
                                        CPLNG%dummy( 1:CPLNG%n_dof,1:CPLNG%n_pts), & !< [-------] 
                                        CPLNG%var_ij(1:CPLNG%n_dof,1:CPLNG%n_pts))   !< [to_send]
if(INOTMASTER) print *, "  var_IJ:", sum( CPLNG%var_ij(1_ip,1:CPLNG%n_pts) )/CPLNG%n_pts !, sum(CPLNG%dummy(1_ip,1:CPLNG%n_pts))/CPLNG%n_pts
    case(IRECV)
      call COU_INTERPOLATE_NODAL_VALUES(WhoID, &
                                        CPLNG%n_dof,    &
                                        CPLNG%var_ji(1:CPLNG%n_dof,1:CPLNG%n_pts), & !< [to_recv] 
                                        CPLNG%dummy( 1:CPLNG%n_dof,1:CPLNG%n_pts))   !< [-------]
if(INOTMASTER) print *, "  var_JI:", sum( CPLNG%var_ji(1_ip,1:CPLNG%n_pts) )/CPLNG%n_pts !, sum(CPLNG%dummy(1_ip,1:CPLNG%n_pts))/CPLNG%n_pts
    case(ISENDRECV)
    !  call COU_INTERPOLATE_NODAL_VALUES(WhoSend, &
    !                                    CPLNG%n_dof,    &
    !                                    CPLNG%dummy(1:CPLNG%n_dof,1:CPLNG%n_pts), & !< [to_recv] 
    !                                    CPLNG%dummy(1:CPLNG%n_dof,1:CPLNG%n_pts))   !< [to_send]
  end select 
  !-----------------------------------------------------------------------||---!  
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if(.not.INOTMASTER) then 
    !
    write(str(1),'(I5.5)') WhoIam 
    write(str(2),'(I5.5)') CPLNG%counter 
    write(str(3),'(I5.5)') current_code 
    !
    print *, '[commdom_alya_exchange'//trim(done)//trim(who)//'] '& 
              //' Current:'//adjustl(str(3))&
              //'Module:'  //adjustl(str(1))&
              //' Order:'  //adjustl(str(2))
    !
  endif 
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!

!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !> @author  JM Zavala-Ake
  !> @date    
  !> @brief   
  !> @details 
  !-----------------------------------------------------------------------||---!
  subroutine commdom_alya_Concou( CPLNG, gocou ) !< switch kfl_gocou
  use mod_parall, only: I_AM_IN_COLOR
  implicit none
  type(COMMDOM_COUPLING), intent(inout) :: CPLNG
  integer(ip), intent(out) :: gocou
  integer(ip)   :: icoup
  integer(ip)   :: color_ij
  character(16) :: str(1)
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  gocou     = 1                           !< set coupling_modules ?? 
  CPLNG%iters_cou = CPLNG%iters_cou + 1   !< number of invocations of subroutine
  !
  do icoup = 1,CPLNG%n_couplings 
      call commdom_alya_combination(CPLNG, icoup, color_ij)
      if(color_ij == 1) then 
!        if(.not.(CPLNG%iters_cou < CPLNG%n_iters_cou)) then
          gocou           = 0
          CPLNG%iters_cou = 0
!        end if
      endif
      !
  enddo 
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  if(.not.INOTMASTER) then 
    write(str(1),'(I6.5)') CPLNG%iters_cou+1
    !if(gocou==1) print *, "[commdom_alya_concou] ", adjustl(str(1)) 
  endif 
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  subroutine commdom_alya_combination( CPLNG, icoup, I_AM_IN_COUPLING ) 
  use mod_parall, only: I_AM_IN_COLOR
  implicit none
  type(COMMDOM_COUPLING), intent(in) :: CPLNG
  integer(ip), intent(in)  :: icoup
  integer(ip), intent(out) :: I_AM_IN_COUPLING
  integer(ip) :: module_source
  integer(ip) :: color_i=-1, color_j=-1
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  module_source = coupling_type(icoup) % module_source
  color_i       = coupling_type(icoup) % color_source
  color_j       = coupling_type(icoup) % color_target
  !
  if((module_source == CPLNG%module_i) .or. (module_source == CPLNG%module_j)) then
  !if((icoup == CPLNG%code_i) .or. (icoup == CPLNG%code_j)) then
    I_AM_IN_COUPLING = -1
    if(I_AM_IN_COLOR(color_i).or.I_AM_IN_COLOR(color_j)) I_AM_IN_COUPLING = 1
  endif 
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!

!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !> @author  JM Zavala-Ake
  !> @date    
  !> @brief   
  !> @details 
  !-----------------------------------------------------------------------||---!
  subroutine commdom_alya_init() 
  use def_kintyp
  use def_master, only: ITASK_DOITER, ITASK_AFTER, ITASK_BEFORE
  use def_inpout, only: words, getint
  implicit none
  integer :: icoup, where_number
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
 !character(*), intent(in),  optional :: wzero
!    coupling_type(icoup) % module_target = idmod(words(3))
!    coupling_type(icoup) % code_target   = int(param(2),ip)

!    coupling_type(icoup) % module_source = idmod(words(3))
!    coupling_type(icoup) % code_source   = int(param(2),ip)

!    coupling_type(icoup) % where_number = ! #Number of the set or the field
!    coupling_type(icoup) % where_type   = ! ON_SET|ON_FIELD|ON_CODE|ON_WHOLE
!    coupling_type(icoup) % what         = ! UNKNOWN|DIRICHLET_IMPLICIT|DIRICHLET_EXPLICIT|RESIDUAL
!    coupling_type(icoup) % task_compute_and_send  ! ITASK_DOITER|ITASK_BEGITE|ITASK_ENDITE|ITASK_TURNON|ITASK_TURNOF
!    coupling_type(icoup) % task_recv_and_assemble ! 
!    coupling_type(icoup) % when_compute_and_send  ! ITASK_AFTER|ITASK_BEFORE
!    coupling_type(icoup) % when_recv_and_assemble   
!
!    coupling_type(icoup) % itera
!  
!  +
!  |_Sources/kernel/coupli/cou_parall
!    xxxxx = target|source 
!    coupling_type(icoup) % color_xxxxx     = par_code_zone_subd_to_color(code_xxxxx, zone_xxxxx, subdomain_xxxxx)
!    coupling_type(icoup) % zone_xxxxx      
!    coupling_type(icoup) % subdomain_xxxxx 
!
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  call ecoute('cou_readat')
  do while( words(1) /= 'ENDPH' )
      if( words(1) == 'NUMBE' ) mcoup = getint('NUMBE',1_ip,'#NUMBER OF COUPLINGS')
     
      if( words(1) == 'COUPL' ) then
           icoup = getint('COUPL',1_ip,'#Number of the set or the field')
           !if( icoup < 1 .or. icoup > mcoup )    call runend('COU_READAT: WRONG COUPLING NUMBER')
           !if( .not. associated(coupling_type) ) call runend('COU_READAT: NUMBER_COUPLING TYPE IS MISSING')

           call ecoute('cou_readat')
           do while( words(1) /= 'ENDCO' )
              if(      words(1) == 'TARGE' .and. words(2) == 'MODUL' ) then
              endif 

              if( words(1) == 'WHERE' ) then
                 where_number = getint('NUMBE',-1_ip,'#Number of the set or the field')
              endif 

              if( words(1) == 'SENDA' )then
              !   if( words(2) == 'DOITE ') coupling_type(icoup) % task_compute_and_send = ITASK_DOITER
              endif 

              if( words(1) == 'RECEI' )then
              !   if( words(2) == 'DOITE ') coupling_type(icoup) % task_recv_and_assemble = ITASK_DOITER
              !   if( words(3) == 'BEFOR' ) coupling_type(icoup) % when_compute_and_send = ITASK_BEFORE
              !   if( words(3) == 'AFTER' ) coupling_type(icoup) % when_compute_and_send = ITASK_AFTER
              endif
              call ecoute('cou_readat') 
           enddo 
      endif 
      call ecoute('cou_readat')
  enddo
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  if(.not.INOTMASTER) then 
    print *, "[commdom_alya_init] ", "where", where_number
  endif 
  !-----------------------------------------------------------------------||---!
  !
  mcoup = 0
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!



  !-----------------------------------------------------------------------||---!
  !> @author  JM Zavala-Ake
  !> @date    
  !> @brief   
  !> @details 
  !-----------------------------------------------------------------------||---!
  !subroutine commdom_alya_XXX()
  !implicit none
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  !end subroutine
  !-----------------------------------------------------------------------||---!

  !=============================================================| contains |===!
end module mod_commdom_alya
!==============================================================================!
!==============================================================================!
