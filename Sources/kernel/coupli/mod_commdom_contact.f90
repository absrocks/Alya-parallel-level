!==============================================================================!
!
!< 2015Mar20 -> created (from 'mod_commdom_alya_contact')
!
!==============================================================================!
module mod_commdom_contact
  use def_parame,        only: ip, rp
  use def_master,        only: inotmaster, imaster,ITASK_TIMSTE
  use def_coupli,        only: mcoup  
  use mod_commdom_alya,  only: COMMDOM_COUPLING
#ifdef COMMDOM 
  use mod_commdom_plepp,   only: commdom_plepp_coupling_send_msg
  use mod_commdom_plepp,   only: commdom_plepp_exchange02
  use mod_commdom_dynamic, only: commdom_dynamic_exchange02
#endif   
  implicit none 
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  type(COMMDOM_COUPLING),save :: CNT_CPLNG   
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  private
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!

  !=============================================================| contains |===!
contains
!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------| INIT |---!
!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_contact_init(CPLNG)
  use def_coupli,       only: UNKNOWN, RESIDUAL  
  use def_master,       only: ID_NASTIN, ID_TEMPER
  use def_master,       only: ITASK_BEFORE, ITASK_AFTER
  use def_master,       only: ITASK_DOITER, ITASK_ENDSTE, ITASK_BEGSTE, ITASK_BEGZON, ITASK_ENDZON
  use def_master,       only: current_code, modul
  use mod_commdom_alya, only: commdom_alya_memall
  use mod_commdom_alya, only: commdom_alya_set_sendrescv_block_type_extreme
  implicit none
  type(COMMDOM_COUPLING), intent(inout) :: CPLNG
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  CPLNG%code_i         =  1_ip        !< CODEi
  CPLNG%module_i       =  ID_TEMPER   !< MODULEi
  CPLNG%fixbo_i        = -1_ip
  CPLNG%what_i         = -RESIDUAL
  !
  CPLNG%code_j         =  2_ip        !< CODEj
  CPLNG%module_j       =  ID_TEMPER   !< MODULEj 
  CPLNG%fixbo_j        = -1_ip
  CPLNG%what_j         = -RESIDUAL  
  !
  CPLNG%n_dof          =  1_ip        !< D.O.F.
  !
  CPLNG%send_task_i    =  ITASK_BEGZON
  CPLNG%send_when_i    =  ITASK_AFTER
  CPLNG%recv_task_i    =  ITASK_BEGZON 
  CPLNG%recv_when_i    =  ITASK_AFTER  
  !
  CPLNG%send_task_j    =  ITASK_ENDZON !ITASK_BEGZON ITASK_BEGSTE !<-- ITASK_DOITER
  CPLNG%send_when_j    =  ITASK_BEFORE
  CPLNG%recv_task_j    =  ITASK_ENDZON !ITASK_BEGZON ITASK_BEGSTE !<-- ITASK_DOITER
  CPLNG%recv_when_j    =  ITASK_BEFORE
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  CPLNG%current_code   = current_code !< *.dat CODE: ID_CODE 
  call commdom_alya_set_sendrescv_block_type_extreme(CPLNG)
  !
  if(CPLNG%code_i==CPLNG%current_code) CPLNG%current_fixbo = CPLNG%fixbo_i
  if(CPLNG%code_j==CPLNG%current_code) CPLNG%current_fixbo = CPLNG%fixbo_j
  !
#ifdef COMMDOM 
  current_code = 1_ip         !< trick!! -> PAR_COLOR_COMMUNICATORS: CODE NUMBER EXCEED NUMBER OF CODES 
  mcoup        = 0_ip         !< evoid cou_turnon 
#endif 
  !
  if(IMASTER) print*, "[commdom_contact_init]"
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
  subroutine commdom_contact_memall(CPLNG)
  use mod_commdom_alya,     only:  commdom_alya_memall
  implicit none
  type(COMMDOM_COUPLING), intent(inout) :: CPLNG
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
#ifdef COMMDOM
  call commdom_alya_memall(CPLNG)
#else
  call commdom_alya_memall(CPLNG)
#endif 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---! 
  CPLNG%sendrecv_idx(2) = 1
  CPLNG%sendrecv_order  = CPLNG%sendrecv_order + 1 !< 2
  if(IMASTER) print*, "[commdom_contact_memall]"
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_contact_driver_plepp(CPLNG, current_when, current_task)
  use mod_commdom_alya,     only: SENDRECV_I, SENDRECV_J
  use mod_commdom_alya,     only: commdom_alya_sendrecv_driver
#ifdef COMMDOM
  use mod_commdom_plepp,    only: commdom_plepp_compare_dtinv
  use mod_commdom_dynamic,  only: commdom_dynamic_set_mesh, commdom_dynamic_deallocate
  use mod_commdom_plepp,    only: PLEPP_CPLNG
!  
use mod_commdom_aitken, only: RELAXATION
use mod_commdom_aitken, only: commdom_aitken_set_prop
!

#endif 
  use mod_commdom_alya,     only: IRECV, ISEND
  use def_master,           only: dtinv, cutim, dtime  
!
use def_master,           only: displ
use def_master,           only: ittim
use def_domain,           only: coord
  use mod_messages, only : livinf
!  
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  implicit none
  integer(ip),            intent(in)    :: current_when
  integer(ip),            intent(in)    :: current_task
  type(COMMDOM_COUPLING), intent(inout) :: CPLNG
  !
  character(32) :: cdummy
  character(8)  :: cdummy_send, cdummy_recv
  integer(ip)   :: cdummy_len
  cdummy_send=''
  cdummy_recv=''
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  call commdom_alya_sendrecv_driver(CPLNG, current_when, current_task)
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
#ifdef COMMDOM 
  !-----------------------------------------------------------------------||---!
  if(CPLNG%sendrecv(1,7)) then 
    !
    coord = coord + displ(:,:,1)
    !
    call commdom_dynamic_deallocate( PLEPP_CPLNG ) 
    call commdom_dynamic_set_mesh( CPLNG%current_fixbo, CPLNG%n_dof) 
    !
    coord = coord - displ(:,:,1)
    !
  else&
  if(CPLNG%sendrecv(2,7)) then
    coord = coord + displ(:,:,1)
    !
    call commdom_dynamic_deallocate( PLEPP_CPLNG ) 
    call commdom_dynamic_set_mesh( CPLNG%current_fixbo, CPLNG%n_dof) 
    !
    coord = coord - displ(:,:,1)
    !
  endif
  !-----------------------------------------------------------------------||---!
  if(CPLNG%sendrecv(1,6)) then 
    call commdom_plepp_compare_dtinv(dtinv) 
    !
    cutim  = cutim - dtime       
    call setgts(ITASK_TIMSTE)             
    call livinf(201_ip, ' ',1_ip)
    !
  else&
  if(CPLNG%sendrecv(2,6)) then
    call commdom_plepp_compare_dtinv(dtinv)
    !
    cutim  = cutim - dtime        
    call setgts(ITASK_TIMSTE)             
    call livinf(201_ip, ' ',1_ip) 
    !
  endif
  !-----------------------------------------------------------------------||---!
  select case(CPLNG%sendrecv_code)
    case(-SENDRECV_I)
!      call commdom_dynamic_exchange02( CPLNG%var_ij(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), &
!                                     CPLNG%var_ji(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), CPLNG%n_dof )
    case(-SENDRECV_J)
!      call commdom_dynamic_exchange02( CPLNG%var_ij(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), &
!                                     CPLNG%var_ji(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), CPLNG%n_dof )
    case( SENDRECV_I )
      CPLNG%current_sendrecv = SENDRECV_I 
!      call commdom_dynamic_exchange02( CPLNG%var_ij(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), &
!                                     CPLNG%var_ji(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), CPLNG%n_dof )
      !
!      call commdom_aitken_set_prop(RELAXATION, &
!                                   CPLNG%var_ji(1:CPLNG%n_dof,PLEPP_CPLNG%interior_list_j), &
!                                   CPLNG%current_iter ) 
      !
    case( SENDRECV_J )
      CPLNG%current_sendrecv = SENDRECV_J  
!      call commdom_dynamic_exchange02( CPLNG%var_ij(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), &
!                                     CPLNG%var_ji(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), CPLNG%n_dof ) 
      !
!      call commdom_aitken_set_prop(RELAXATION, &
!                                   CPLNG%var_ji(1:CPLNG%n_dof,PLEPP_CPLNG%interior_list_j), &
!                                   CPLNG%current_iter ) 
      !
    case default
      
  end select
#endif 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!

!-------------------------------------------------------------------------||---!
!----------------------------------------------------------------| ??????? |---!
!-------------------------------------------------------------------------||---!
end module mod_commdom_contact
!==============================================================================!
!==============================================================================! 
