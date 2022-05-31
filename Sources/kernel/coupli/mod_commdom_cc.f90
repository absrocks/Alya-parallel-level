!==============================================================================!
!
!< 2015Feb02 
! +modules/nastal/nsa_upcons.f90 
!  solver_solve 
!
! +kernel/master/soldef.f90 
!  if( modul == ID_NASTAL )  solve_sol(ivari) % kfl_version = 1 
!
!
!==============================================================================!
!
!==============================================================================!
module mod_commdom_cc
  !=================================================================| init |===!
  use def_parame,    only: ip, rp
  use def_master,    only: inotmaster, imaster, isequen, islave, inotslave, iparall,ITASK_TIMSTE
  use def_master,    only: kfl_gocou !, mem_modul
  use def_master,    only: ITASK_TURNON, MODUL 
  use def_domain,    only: coord, mnode
  use def_domain,    only: ltype, lnods
  use mod_couplings, only: COU_INTERPOLATE_NODAL_VALUES
  use def_coupli,    only: coupling_type, mcoup  
  use mod_commdom_alya, only: COMMDOM_COUPLING
  use mod_messages, only : livinf
#ifdef COMMDOM 
  use mod_commdom_plepp, only: commdom_plepp_coupling_send_msg
  use mod_commdom_plepp, only: commdom_plepp_exchange01
  use mod_commdom_plepp, only: commdom_plepp_exchange02
#endif   
  implicit none 
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  private
    ! + 
    ! |_Alya
    !   |_kernel/domain/domain
          public:: commdom_cc_init
          public:: commdom_cc_memall
    !   |
    !   |_kernel/coupli/mod_coupling_driver
          public:: commdom_cc_driver_plepp
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
  subroutine commdom_cc_init(CPLNG)
  use def_coupli,       only: UNKNOWN, RESIDUAL  
  use def_master,       only: ID_NASTIN, ID_TEMPER, ID_NASTAL
  use def_master,       only: ITASK_BEFORE, ITASK_AFTER
  use def_master,       only: ITASK_DOITER, ITASK_ENDSTE, ITASK_BEGSTE
  use def_master,       only: current_code, modul
use def_domain, only: ndime 
  use mod_commdom_alya, only: commdom_alya_memall
  use mod_commdom_alya, only: commdom_alya_set_sendrescv_block_type_extreme
!
  use def_kintyp,           only: soltyp
  use def_master,           only: momod, modul
!
  implicit none
  type(COMMDOM_COUPLING), intent(inout) :: CPLNG 
  type(soltyp), pointer :: solve(:)
  !-----------------------------------------------------------------------||---!
  !
  CPLNG%coupling_i     = 1_ip       
  CPLNG%coupling_j     = 2_ip
  !
  !< SOURCEi, TARGETj
  CPLNG%code_i         =  1_ip       !< CODE.   coupling%code_source
  CPLNG%module_i       =  ID_NASTAL  !< MODULE. coupling%module_source
  CPLNG%fixbo_i        =  1_ip
  CPLNG%what_i         =  RESIDUAL   !< kfl_bvnat==1,  <-  react 
  !
  !< SOURCEj, TARGETi
  CPLNG%code_j         =  2_ip       !< CODE.   coupling% code_target 
  CPLNG%module_j       =  ID_NASTAL  !< MODULE. coupling% module_target 
  CPLNG%fixbo_j        =  1_ip
  CPLNG%what_j         = -RESIDUAL   !< kfl_react==1,   -> react 
  !
  CPLNG%n_dof          =  ndime+2 
!
!print *, "modul", modul
!  solve => momod(ID_NASTAL) % solve(1:)
!print *, "ndofn", solve(1) % ndofn
!
  !
  !< SEND_AND_COMPUTE|RECEIVE_ASSEMBLE
  !< tasks: ITASK_BEGITE|ITASK_TURNON|ITASK_DOITER 
  !<  when: ITASK_BEFORE|ITASK_AFTER
  CPLNG%send_task_i    = ITASK_ENDSTE
  CPLNG%send_when_i    = ITASK_AFTER
  CPLNG%recv_task_i    = ITASK_ENDSTE !ITASK_BEGSTE 
  CPLNG%recv_when_i    = ITASK_AFTER  !ITASK_BEFORE
  !
  CPLNG%send_task_j    = ITASK_DOITER
  CPLNG%send_when_j    = ITASK_BEFORE
  CPLNG%recv_task_j    = ITASK_DOITER
  CPLNG%recv_when_j    = ITASK_BEFORE
  !
  CPLNG%current_code   = current_code !< *.dat CODE: ID_CODE 

!print *, "modul", modul 

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
  CPLNG%sendrecv_idx(1) = 1
  CPLNG%sendrecv_order  = 1
  !
  if(IMASTER) print*, "[commdom_cc_init]"
  !
  !< nblok                  ! Number of blocks
  !< mmodu                  ! Max. # of modules
  !< mblok                  ! Max. # of blocks
  !< lmord(mmodu, mblok)    ! Order of modules iterations
  !< micou(mblok)           ! Max. # of global iter. for each block
  !< kfl_coupl(imodu,jmodu) !'jmodu->imodu 2' means that jmodu modifies imodu using procedure number 2
  !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_cc_memall(CPLNG)
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
  if(IMASTER) print*, "[commdom_cc_memall]"
  end subroutine
  !-----------------------------------------------------------------------||---!



  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_cc_driver_plepp(CPLNG, current_when, current_task)
  use mod_commdom_alya,  only: SENDRECV_I, SENDRECV_J
  use mod_commdom_alya,  only: commdom_alya_sendrecv_driver
#ifdef COMMDOM
  use mod_commdom_plepp,  only: commdom_plepp_compare_dtinv
#endif 
  use mod_commdom_alya,  only: IRECV, ISEND
  use def_master,        only: dtinv, cutim, dtime  
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
  if(CPLNG%sendrecv(1,6)) then 
   !cdummy_send = "dt_i"
   !cdummy_len  = len_trim( cdummy_send )
   !call commdom_plepp_coupling_send_msg(cdummy_len, trim(cdummy_send), cdummy_recv) 
    call commdom_plepp_compare_dtinv(dtinv) 
    !
    cutim  = cutim - dtime          ! \
   !call iniste(2_ip)               ! |__ 2014Dic10, <-- Begste, if(kfl_reset == 1) 
    call setgts(ITASK_TIMSTE)               ! |
    call livinf(201_ip, ' ',1_ip)   ! /
    !
  else&
  if(CPLNG%sendrecv(2,6)) then
   !cdummy_send = "dt_j"
   !cdummy_len  = len_trim( cdummy_send )
   !call commdom_plepp_coupling_send_msg(cdummy_len, trim(cdummy_send), cdummy_recv) 
    call commdom_plepp_compare_dtinv(dtinv)
    !
    cutim  = cutim - dtime          ! \
   !call iniste(2_ip)               ! |__ 2014Dic10, <-- Begste, if(kfl_reset == 1)  
    call setgts(ITASK_TIMSTE)               ! |
    call livinf(201_ip, ' ',1_ip)   ! /
    !
  endif
  !-----------------------------------------------------------------------||---!
  select case(CPLNG%sendrecv_code)
    case(-SENDRECV_I)
      cdummy_send = "XXxxXXxx"
      cdummy_len  = len_trim( cdummy_send )
      call commdom_plepp_coupling_send_msg(cdummy_len, trim(cdummy_send), cdummy_recv)
      call commdom_plepp_exchange02( CPLNG%var_ij(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), &
                                     CPLNG%var_ji(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), CPLNG%n_dof )
    case(-SENDRECV_J)
      cdummy_send = "YYyyYYyy"
      cdummy_len  = len_trim( cdummy_send )
      call commdom_plepp_coupling_send_msg(cdummy_len, trim(cdummy_send), cdummy_recv)
      call commdom_plepp_exchange02( CPLNG%var_ij(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), &
                                     CPLNG%var_ji(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), CPLNG%n_dof )
    case( SENDRECV_I)
      cdummy_send = "XXxxXXxx"
      cdummy_len  = len_trim( cdummy_send )
      call commdom_plepp_coupling_send_msg(cdummy_len, trim(cdummy_send), cdummy_recv)
      CPLNG%current_sendrecv = SENDRECV_I 
      call commdom_plepp_exchange02( CPLNG%var_ij(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), &
                                     CPLNG%var_ji(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), CPLNG%n_dof )
    case( SENDRECV_J)
      cdummy_send = "YYyyYYyy"
      cdummy_len  = len_trim( cdummy_send )
      call commdom_plepp_coupling_send_msg(cdummy_len, trim(cdummy_send), cdummy_recv)
      CPLNG%current_sendrecv = SENDRECV_J  
      call commdom_plepp_exchange02( CPLNG%var_ij(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), &
                                     CPLNG%var_ji(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), CPLNG%n_dof )
    case default
     !call runend('EXIT!!')
  end select
#endif 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  CPLNG%sendrecv_idx(3) = 1
  !if(IMASTER.and.CPLNG%sendrecv_code/=0) print*, "[commdom_cc_driver_plepp]"
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!

!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!

  !=============================================================| contains |===!
end module mod_commdom_cc
!==============================================================================!
!==============================================================================! 

  !-----------------------------------------------------------------------||---!
  !> @author  JM Zavala-Ake
  !> @date    
  !> @brief   
  !> @details 
  !-----------------------------------------------------------------------||---!
  !subroutine commdom_alya_XXX()
  !implicit none
  !end subroutine
  !-----------------------------------------------------------------------||---!

!==============================================================================!
!==============================================================================! 
