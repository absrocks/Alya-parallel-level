!==============================================================================!
!
!< 2014Sep26. creation  
!
!==============================================================================!
module mod_commdom_lmc 
  use def_parame,    only: ip, rp
  use def_master,    only: inotmaster, imaster, isequen, islave, inotslave, iparall
  use def_master,    only: kfl_gocou !, mem_modul
  use def_master,    only: ITASK_TURNON, MODUL 
  use def_domain,    only: coord, mnode
  use def_domain,    only: ltype, lnods
  use mod_couplings, only: COU_INTERPOLATE_NODAL_VALUES
  use def_coupli,    only: coupling_type, mcoup  
  use mod_commdom_alya, only: COMMDOM_COUPLING
#ifdef COMMDOM 
  use mod_commdom_plepp, only: commdom_plepp_coupling_send_msg
  use mod_commdom_plepp, only: commdom_plepp_exchange01
  use mod_commdom_plepp, only: commdom_plepp_exchange02
#endif   
  !use def_domain,   only: nelem, ndime, npoin, nnode, ngaus
  !use mod_memchk,   only: memchk
  implicit none 
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  type(COMMDOM_COUPLING), save :: LMC_CPLNG  !< Low Mach-Compressible 
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  private
    public:: commdom_lmc_init
    public:: commdom_lmc_memall
    public:: commdom_lmc_driver_plepp
    public:: LMC_CPLNG
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
  subroutine commdom_lmc_init(CPLNG)
  use def_master,       only: ID_NASTIN, ID_NASTAL
  use def_master,       only: ITASK_BEFORE, ITASK_AFTER
  use def_master,       only: ITASK_DOITER, ITASK_ENDSTE, ITASK_BEGSTE
  use def_master,       only: current_code, modul
  use mod_commdom_alya, only: commdom_alya_memall
  use mod_commdom_alya, only: commdom_alya_set_sendrescv_block_type_extreme
  implicit none
  type(COMMDOM_COUPLING), intent(inout) :: CPLNG
  !-----------------------------------------------------------------------||---!
  !
  CPLNG%coupling_i     = 1_ip       
  CPLNG%coupling_j     = 2_ip
  !
  !< SOURCEi, TARGETj
  CPLNG%code_i         =  1_ip       !< CODE.   coupling%code_source
  CPLNG%module_i       =  ID_NASTAL  !< MODULE. coupling%module_source
  CPLNG%fixbo_i        =  5_ip
  !
  !< SOURCEj, TARGETi
  CPLNG%code_j         =  2_ip       !< CODE.   coupling% code_target 
  CPLNG%module_j       =  ID_NASTIN  !< MODULE. coupling% module_target 
  CPLNG%fixbo_j        =  1_ip
  !
  CPLNG%n_dof          =  3_ip
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
  CPLNG%sendrecv_order  = 1
  if(IMASTER.or.ISEQUEN) print*, "[commdom_lmc_init]"
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
  subroutine commdom_lmc_memall(CPLNG)
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
  CPLNG%sendrecv_order  = CPLNG%sendrecv_order + 1 !< 2
  if(IMASTER) print*, "[commdom_lmc_memall]"
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_lmc_driver_alya(CPLNG, current_when, current_task)
  use def_master,           only: modul, kfl_goblk, imaster, kfl_gocou
  use mod_commdom_alya,     only: commdom_alya_sendrecv_driver
  use def_coupli,           only: mcoup
  use def_master,           only: lmord, iblok, mmodu, mblok, nblok, ITASK_ENDSTE, ITASK_BEGSTE, ITASK_INIUNK, ITASK_AFTER
  use mod_commdom_alya,     only: SENDRECV_I, SENDRECV_J
  use mod_commdom_alya,     only: IRECV, ISEND
  use mod_commdom_alya,     only: commdom_alya_exchange
  implicit none
  integer(ip),  intent(in)  :: current_when
  integer(ip),  intent(in)  :: current_task
  type(COMMDOM_COUPLING), intent(inout) :: CPLNG
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  call commdom_alya_sendrecv_driver(CPLNG, current_when, current_task)
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if(IMASTER.and.CPLNG%sendrecv_code/=0) print*, "[commdom_lmc_driver_alya]"
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_lmc_driver_plepp(CPLNG, current_when, current_task)
  use mod_commdom_alya,  only: SENDRECV_I, SENDRECV_J
  use mod_commdom_alya,  only: commdom_alya_sendrecv_driver
#ifdef COMMDOM
  use mod_commdom_plepp,  only: commdom_plepp_compare_dtinv
#endif 
  use mod_commdom_alya,  only: IRECV, ISEND
  use def_master,        only: dtinv 
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
  integer(ip)   :: cdummy_len, idx_exchange
  cdummy_send=''
  cdummy_recv=''
  idx_exchange=1_ip 
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
    call commdom_plepp_compare_dtinv(dtinv) 
  else&
  if(CPLNG%sendrecv(2,6)) then
    call commdom_plepp_compare_dtinv(dtinv)
  endif
  !-----------------------------------------------------------------------||---!
  select case(CPLNG%sendrecv_code)
    case(-SENDRECV_I)
      cdummy_send = "XXxxXXxx"
      cdummy_len  = len_trim( cdummy_send )
      call commdom_plepp_coupling_send_msg(cdummy_len, trim(cdummy_send), cdummy_recv)
      call commdom_plepp_exchange02( CPLNG%var_ij(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), CPLNG%var_ji(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), CPLNG%n_dof )
    case(-SENDRECV_J)
      cdummy_send = "YYyyYYyy"
      cdummy_len  = len_trim( cdummy_send )
      call commdom_plepp_coupling_send_msg(cdummy_len, trim(cdummy_send), cdummy_recv)
      call commdom_plepp_exchange02( CPLNG%var_ij(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), CPLNG%var_ji(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), CPLNG%n_dof )
    case( SENDRECV_I)
      cdummy_send = "XXxxXXxx"
      cdummy_len  = len_trim( cdummy_send )
      call commdom_plepp_coupling_send_msg(cdummy_len, trim(cdummy_send), cdummy_recv)
      call commdom_plepp_exchange02( CPLNG%var_ij(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), CPLNG%var_ji(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), CPLNG%n_dof )
    case( SENDRECV_J)
      cdummy_send = "YYyyYYyy"
      cdummy_len  = len_trim( cdummy_send )
      call commdom_plepp_coupling_send_msg(cdummy_len, trim(cdummy_send), cdummy_recv)
      call commdom_plepp_exchange02( CPLNG%var_ij(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), CPLNG%var_ji(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), CPLNG%n_dof )
    case default
     !call runend('EXIT!!')
  end select
#endif 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  !if(IMASTER.and.CPLNG%sendrecv_code/=0) print*, "[commdom_lmc_driver_plepp]"
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


!-------------------------------------------------------------------------||---!
!----------------------------------------------------------------| ??????? |---!
!-------------------------------------------------------------------------||---!



  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
end module mod_commdom_lmc 
!==============================================================================!
!==============================================================================! 
