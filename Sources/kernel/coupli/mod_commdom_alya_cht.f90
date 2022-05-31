!==============================================================================!
!
!< 2014Sep24 -> CPLNG%current_fixbo  
!< 2014Nov04 -> commdom_cht_get_vals
!<           -> commdom_plepp_exchange02
!<           -> n_dof 
!< 2014Dic10 -> reseting the time step as it is done by 'Begste' 
!< 2014Dic14 -> relax
!
!==============================================================================!
!
!  +./Sources/kernel/master/Timste.f90
!< |_ subroutine Timste()
!  |_ ...
!< |_ call moduls(ITASK_TIMSTE)
!< |_ call par_commdom_compare_dt(dtinv)
!  |_ ...
!
!  +./Sources/modules/temper/tem_commdom.f90
!  |_ ...
!< |_ k*grad(T).n = qr + ar*(T-Tr)
!< |_ ar: 'arobi = bvnat_tem(1,iboun,1)' <- heat exchange 
!< |_ Tr: 'trobi = bvnat_tem(2,iboun,1)' <- reference temperature 
!< |_ qr: 'qrobi = bvnat_tem(3,iboun,1)' <- heat flux 
!  |_ ...
!
!  + *.tem.dat 
!< |_ CODES, BOUNDARIES
!< |_   icode  2  ar Tr  qr ->  kfl_fixbo_tem(nboun) = 2
!< |_ END_CODES_BOUNDARIES
!
!  +./Sources/modules/temper/Temper.f90 
!  |_ ...
!  |_ call commdom_alya_cht_nodes2bvnat( heat_flux(1:npoin), bvnat_tem(3,1:nboun,1) )
!  |_ ...
!  |_ call tem_doiter()
!  |_ ...
!
!==============================================================================!
module mod_commdom_alya_cht
  !=================================================================| init |===!
  !-----------------------------------------------------------------------||---!
  !< Sources/kernel/coupli/def_coupli.f90
  !<          o----o----o----o Source
  !<          x------x-------x Target
  !-----------------------------------------------------------------------||---!
  use def_parame,    only: ip, rp
  use def_master,    only: inotmaster, imaster, isequen, islave, inotslave, iparall,ITASK_TIMSTE
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
  !include 'mpif.h'  
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!

  type(COMMDOM_COUPLING),save :: CHT_CPLNG  !< Conjugate Heat Transfer 
 !type(COMMDOM_COUPLING) :: LMC_CPLNG  !< Low Mach-Compressible 

  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  private
    ! + 
    ! |_Nastin
    !   |_nsi_turnon
    !     |_nsi_memall
    public:: commdom_alya_cht_driver
    public:: commdom_alya_init_cht 
    public:: commdom_cht_memall
    public:: commdom_cht_driver_alya  
    public:: commdom_cht_driver_plepp
    public:: commdom_cht_get_vals
    public:: CHT_CPLNG
!    public:: LMC_CPLNG
!    public:: commdom_alya_init_lowmach_compressible
!    public:: commdom_alya_cht_lowmach_after_endste
!    public:: commdom_alya_cht_tempe_before_doiter
!    public:: commdom_alya_temper_driver
!    public:: commdom_alya_lowmach_driver
  !
  ! public:: 
  !
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
  !
  ! +
  ! |_NASTIN+TEMPER 
  !   |_ITASK_ENDSTE
  !                 \_ Flux,   T  <-- send, recv 
  !
  ! +
  ! |_TEMPER         _    T, q_r  <-- send, recv 
  !   |             / 
  !   |_ITASK_DOITER
  !
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_alya_init_cht(CPLNG)
  use def_coupli,       only: UNKNOWN, RESIDUAL  
  use def_master,       only: ID_NASTIN, ID_TEMPER
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
  CPLNG%module_i       =  ID_TEMPER  !< MODULE. coupling%module_source
  CPLNG%fixbo_i        =  1_ip
  CPLNG%what_i         = -RESIDUAL
  !
  !< SOURCEj, TARGETi
  CPLNG%code_j         =  2_ip       !< CODE.   coupling% code_target 
  CPLNG%module_j       =  ID_TEMPER  !< MODULE. coupling% module_target 
  CPLNG%fixbo_j        =  1_ip
  CPLNG%what_j         = -RESIDUAL  
  !
  CPLNG%n_dof          =  1_ip
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
  CPLNG%sendrecv_idx(1) = 1
  CPLNG%sendrecv_order  = 1
  !
  if(IMASTER) print*, "[commdom_cht_init]"
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
  subroutine commdom_cht_memall(CPLNG)
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
  if(IMASTER) print*, "[commdom_cht_memall]"
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_cht_driver_alya(CPLNG, current_when, current_task)
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
  if(CPLNG%sendrecv(1,6)) then 
  else&
  if(CPLNG%sendrecv(2,6)) then  
  endif 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  select case(CPLNG%sendrecv_code)
    case(-SENDRECV_I)
      CPLNG%var_ij = 1.0
      call commdom_alya_exchange(CPLNG, CPLNG%module_i, ISEND)
    case(-SENDRECV_J)
      CPLNG%var_ji = 0.0
      call commdom_alya_exchange(CPLNG, CPLNG%module_i, IRECV)
    case( SENDRECV_I)
      CPLNG%var_ij = 1.0
      call commdom_alya_exchange(CPLNG, CPLNG%module_i, ISEND)
    case( SENDRECV_J)
      CPLNG%var_ji = 0.0
      call commdom_alya_exchange(CPLNG, CPLNG%module_i, IRECV)
    case default
     !call runend('EXIT!!')
  end select
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  CPLNG%sendrecv_idx(3) = 1
  if(IMASTER.and.CPLNG%sendrecv_code/=0) print*, "[commdom_cht_driver_alya]"
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_cht_driver_plepp(CPLNG, current_when, current_task)
  use mod_commdom_alya,  only: SENDRECV_I, SENDRECV_J
  use mod_commdom_alya,  only: commdom_alya_sendrecv_driver
#ifdef COMMDOM
  use mod_commdom_plepp,  only: commdom_plepp_compare_dtinv
#endif 
  use mod_commdom_alya,  only: IRECV, ISEND
  use def_master,        only: dtinv, cutim, dtime  
  use mod_messages, only : livinf
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
!!      call commdom_plepp_exchange01(CPLNG%var_ij(1_ip,1_ip:CPLNG%n_pts), CPLNG%var_ji(1_ip,1_ip:CPLNG%n_pts))
    case( SENDRECV_J)
      cdummy_send = "YYyyYYyy"
      cdummy_len  = len_trim( cdummy_send )
      call commdom_plepp_coupling_send_msg(cdummy_len, trim(cdummy_send), cdummy_recv)
      CPLNG%current_sendrecv = SENDRECV_J  
      call commdom_plepp_exchange02( CPLNG%var_ij(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), &
                                     CPLNG%var_ji(1_ip:CPLNG%n_dof,1_ip:CPLNG%n_pts), CPLNG%n_dof )
!      call commdom_plepp_exchange01(CPLNG%var_ij(1_ip,1_ip:CPLNG%n_pts), CPLNG%var_ji(1_ip,1_ip:CPLNG%n_pts))
    case default
     !call runend('EXIT!!')
  end select
#endif 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  CPLNG%sendrecv_idx(3) = 1
  !if(IMASTER.and.CPLNG%sendrecv_code/=0) print*, "[commdom_cht_driver_plepp]"
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_cht_get_vals(input, output, relax_in)
  use def_domain,         only: npoin
#ifdef COMMDOM
  use mod_commdom_plepp,  only: commdom_plepp_set_vals
#endif 
  implicit none 
  real(rp),  intent( in)    ::  input(npoin) 
  real(rp),  intent(out)    :: output(npoin) 
  !type(COMMDOM_COUPLING), intent(inout) :: CPLNG
  !
  real(rp) :: relax = 1.0                    !  == 2014Dic14, relax
  real(rp), optional, intent(in) :: relax_in ! \
  if( present(relax_in) ) relax = relax_in   ! /
  !
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
#ifdef COMMDOM
  call commdom_plepp_set_vals(input, output, relax) 
#else
  output = 0.0 
  print *, "[commdom_cht_set_vals] ERROR: no alya set!!"
  call runend("EXIT!!")
#endif
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
  subroutine commdom_cht_lowmach_driver(CPLNG, current_task, current_when)
  use mod_commdom_alya,  only: SENDRECV_I, SENDRECV_J
  use mod_commdom_alya,  only: IRECV, ISEND
  use def_master,        only: ITASK_BEGSTE, ITASK_CONCOU 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  implicit none
  type(COMMDOM_COUPLING), intent(in)           :: CPLNG
  integer(ip),            intent(in)           :: current_task
  integer(ip),            intent(in), optional :: current_when
  !
  if( present(current_when) ) then
  endif
  ! 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!




!-------------------------------------------------------------------------||---!
!----------------------------------------------------------------| LOWMACH |---!
!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !> @author  JM Zavala-Ake
  !> @date    
  !> @brief   
  !> @details 
  !
  !  case(ITASK_ENDSTE)
  !    call tem_endste()
  !    call commdom_alya_cht_lowmach_after_endste() 
  ! 
  !-----------------------------------------------------------------------||---!
  subroutine commdom_alya_lowmach_driver_XXX(current_when) 

!  use def_master,           only: ITASK_AFTER, ITASK_BEFORE 
!  use def_domain,           only: kfl_codno
!  use def_temper,           only: kfl_regim_tem
!  use def_nastin,           only: kfl_regim_nsi
!  use mod_commdom_alya,     only: ISEND, IRECV, ISENDRECV
!  use mod_commdom_alya,     only: commdom_alya_coupling_driver_i
  implicit none
  integer(ip), intent(in), optional  :: current_when
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
!  if((CHT_CPLNG%module_i==CHT_CPLNG%current_module).and.(kfl_regim_tem==3).and.(kfl_regim_nsi==3)) then 
!    call commdom_alya_coupling_driver_i(CHT_CPLNG, ISEND, ITASK_AFTER, f_after_ini=commdom_alya_cht_lowmach_after_endste, dummy=ISEND)
!    call commdom_alya_coupling_driver_i(CHT_CPLNG, IRECV, ITASK_AFTER, f_after_end=commdom_alya_cht_lowmach_after_endste, dummy=IRECV)
!  endif 
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
  !
  !  case(ITASK_ENDSTE)
  !    call tem_endste()
  !    call commdom_alya_cht_lowmach_after_endste() 
  ! 
  !-----------------------------------------------------------------------||---!
  subroutine commdom_alya_cht_lowmach_after_endste(sendrecv_code) 
!  use def_temper,           only: kfl_regim_tem, bvess_tem, bvnat_tem   
!  use def_master,           only: tempe, current_code
!  use def_domain,           only: npoin, nboun 
!  use def_coupli,           only: typ_color_coupling, coupling_type
!  use mod_commdom_alya,     only: ISEND, IRECV, ISENDRECV
  implicit none
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  integer(ip), intent(in) :: sendrecv_code 

!  type(typ_color_coupling), pointer :: alya_coupling => null() 

!  integer(ip), pointer :: wets_status(:)   => null() 
!  integer(ip), pointer :: wets_lpoin(:)    => null() 
!  real(rp),    pointer :: wets_coords(:,:) => null() 
!  integer(ip)   :: n_wets=-1, n_bnds=-1, i_wet, i_pts
!  integer(ip)   :: ipoin, kpoin, icoup, npoin_wet
!  integer(ip), pointer  :: wets_list(:)
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
!  select case(sendrecv_code) 
  !-----------------------------------------------------------------------||---!
!  case(ISEND)
!  if(kfl_regim_tem==3) then
!    if(INOTMASTER) then
!      CHT_CPLNG%var_ij(1_ip,1_ip:CHT_CPLNG%n_pts) = 0.0
!      call commdom_alya_cht_node_flux( CHT_CPLNG%var_ij(1_ip,1_ip:CHT_CPLNG%n_pts) ) 
!      CHT_CPLNG%var_ij(1_ip,1_ip:CHT_CPLNG%n_pts) = -CHT_CPLNG%var_ij(1_ip,1_ip:CHT_CPLNG%n_pts)
!    end if
!  endif 
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
!  case(IRECV)
!  if(kfl_regim_tem==3) then
    !
!    icoup = CHT_CPLNG%coupling_j
!    if(CHT_CPLNG%current_module == coupling_type(icoup) % module_target .and. INOTMASTER ) then
!      wets_lpoin => coupling_type(icoup) % geome % lpoin_wet
!      print *, minval(wets_lpoin), maxval(wets_lpoin), npoin 
!    endif
    ! 
!    if(INOTMASTER) then
!      print *, "CHT_CPLNG%var_ji:", sum( CHT_CPLNG%var_ji(1_ip,1_ip:CHT_CPLNG%n_pts) )
!      bvess_tem(1,1:npoin,1) = CHT_CPLNG%var_ji(1_ip,1_ip:CHT_CPLNG%n_pts) 
     !bvess_tem(1,1:npoin,2) = CHT_CPLNG%var_ji(1_ip,1_ip:CHT_CPLNG%n_pts) ¿segmentation fault?
!      tempe(1:npoin,1)       = CHT_CPLNG%var_ji(1_ip,1_ip:CHT_CPLNG%n_pts) 
!      tempe(1:npoin,2)       = CHT_CPLNG%var_ji(1_ip,1_ip:CHT_CPLNG%n_pts) 
!      tempe(1:npoin,3)       = CHT_CPLNG%var_ji(1_ip,1_ip:CHT_CPLNG%n_pts) 
!    endif 
!  endif
  !-----------------------------------------------------------------------||---!
!    case default 
  !-----------------------------------------------------------------------||---!
!    call runend('[cht_lowmach_after_endste] sendrecv_code==ISEND|IRECV')
  !-----------------------------------------------------------------------||---!
!  end select
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  !< to sent:    CHT_CPLNG%var_ij
  !< to recieve: CHT_CPLNG%var_ji
  !call par_commdom_locator_exchange01( CHT_CPLNG%var_ij(1:npoin,1_ip), CHT_CPLNG%var_ji(1:npoin) ) !> send, recv 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!

!-------------------------------------------------------------------------||---!
!-----------------------------------------------------------------| TEMPER |---!
!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  subroutine commdom_alya_temper_driver(current_when) 
!  use def_master,           only: ITASK_AFTER, ITASK_BEFORE 
!  use def_domain,           only: kfl_codno
!  use def_temper,           only: kfl_regim_tem
!  use mod_commdom_alya,     only: ISEND, IRECV
!  use mod_commdom_alya,     only: commdom_alya_coupling_driver_j
  implicit none
  integer(ip), intent(in), optional  :: current_when
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
!  if((CHT_CPLNG%module_j==modul).and.(kfl_regim_tem==0)) then 
!    call commdom_alya_coupling_driver_j(CHT_CPLNG, IRECV, ITASK_AFTER, f_after_end=commdom_alya_cht_tempe_before_doiter, dummy=IRECV)
!    call commdom_alya_coupling_driver_j(CHT_CPLNG, ISEND, ITASK_AFTER, f_after_ini=commdom_alya_cht_tempe_before_doiter, dummy=ISEND)
!  endif 
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
  !
  !  case(ITASK_DOITER)
  !    call commdom_alya_cht_temper_before_doite() 
  !    call tem_doiter()
  ! 
  !<  k*grad(T).n = qr+ar*(T-Tr)
  !< 'qrobi = bvnat_tem(3,iboun,1)'
  !
  !< CODES, BOUNDARIES
  !<   icode  2  ar Tr  qr
  !< END_CODES_BOUNDARIES
  !CHT_CPLNG%var_ji(1:npoin) = -CHT_CPLNG%var_ji(1:npoin)
  !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_alya_cht_tempe_before_doiter(sendrecv_code) 
!  use def_temper,           only: kfl_regim_tem, bvess_tem, bvnat_tem   
!  use def_master,           only: tempe
!  use def_domain,           only: npoin, nboun  
!  use mod_commdom_alya,     only: ISEND, IRECV, ISENDRECV
  implicit none
  integer(ip), intent(in) :: sendrecv_code 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  !< to sent:    CHT_CPLNG%var_ij
  !< to recieve: CHT_CPLNG%var_ji
  !call par_commdom_locator_exchange01( CHT_CPLNG%var_ij(1:npoin,1_ip), CHT_CPLNG%var_ji(1:npoin) ) !> send, recv 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
!  select case(sendrecv_code) 
  !-----------------------------------------------------------------------||---!
!  case(ISEND)
!  if(kfl_regim_tem==0) then
!    if(INOTMASTER) then 
!      CHT_CPLNG%var_ij(1_ip,1:CHT_CPLNG%n_pts) = tempe(1:npoin,1)  
!      print *, "CHT_CPLNG%var_ij:", sum( CHT_CPLNG%var_ij(1_ip,1:CHT_CPLNG%n_pts) )/CHT_CPLNG%n_pts, CHT_CPLNG%n_pts
!    endif 
!  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
!  case(IRECV)
!  if(kfl_regim_tem==0) then
!    if(INOTMASTER) then 
!      bvnat_tem(:,1:nboun,:) = 0.0 
!      call commdom_alya_cht_nodes2bvnat( CHT_CPLNG%var_ji(1_ip,1_ip:CHT_CPLNG%n_pts), bvnat_tem(3,1:nboun,1) )
!    endif 
!  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
!  case default 
!  call runend('[cht_lowmach_after_endste] sendrecv_code==ISEND|IRECV')
  !-----------------------------------------------------------------------||---!
!  end select
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
  !
  !  case(ITASK_ENDSTE)
  !    call tem_endste()
  !    call commdom_alya_cht_lowmach_after_endste() 
  ! 
  !-----------------------------------------------------------------------||---!
  subroutine commdom_alya_cht_nodes2bvnat02(n_prop, b_prop, i_fixbo) 
  use def_domain
!  use def_master, only: INOTMASTER
!  use def_temper, only: kfl_fixbo_tem
!  use def_elmtyp, only: TRI03, TRI06, TET04, HEX08 
  implicit none
  real(rp),   intent( in) :: n_prop(npoin), i_fixbo 
  real(rp),   intent(in) :: b_prop(nboun) 

!  real(rp)    :: bprop(mnodb)
!  real(rp)    :: xbprop(mgaus) 
!  real(rp)    :: gbsur
!  real(rp)    :: bocod(ndime,mnodb)
!  real(rp)    :: elrhs(mnode)
!  real(rp)    :: xbocod(ndime,mgaus)

!  integer(ip) :: boidx(mnodb)
!  integer(ip) :: igaub, iboun, pgaub, pnodb, pblty
  !
  !
!  if(INOTMASTER) then
  !
!  boundaries: &
!  do iboun = 1,nboun
    !
!    b_prop(iboun) = 0.0
    !
!    if(kfl_fixbo_tem(iboun) == i_fixbo) then
      !
!      pblty = ltypb(iboun)
      !
!      tria03: & 
!      if(pblty == TRI03) then 
!        pnodb          = nnode(pblty)
!        pgaub          = ngaus(pblty) !< pgaub == 1
!        boidx(1:pnodb) = lnodb(1:pnodb,iboun)
        !
!        bprop(1:pnodb)  = n_prop( boidx(1:pnodb) ) 
!        do igaub = 1,pgaub
!          xbprop(igaub) = dot_product( bprop(1:pnodb), elmar(pblty) % shape(1:pnodb,igaub) )
!        enddo 
        !
!        b_prop(iboun) = sum( xbprop(1:pgaub) )/pgaub
        !
!      endif tria03
      !
!    endif 
    !
!  end do boundaries
!  endif
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !> @author  JM Zavala-Ake
  !> @date    
  !> @brief   
  !> @details 
  !
  !  case(ITASK_ENDSTE)
  !    call tem_endste()
  !    call commdom_alya_cht_lowmach_after_endste() 
  ! 
  !-----------------------------------------------------------------------||---!
  subroutine commdom_alya_cht_nodes2bvnat(prop, h_flux) !< ok  
!  use def_parame
!  use def_master
  use def_domain
!  use def_temper
!  use def_elmtyp, only: TRI03, TRI06, TET04, HEX08 

  implicit none
  real(rp),   intent( in) :: prop(npoin) 
  real(rp),   intent( in) :: h_flux(nboun) 

!  real(rp)    :: bprop(mnodb)
!  real(rp)    :: xbprop(mgaus) 
!  real(rp)    :: gbsur
!  real(rp)    :: bocod(ndime,mnodb)
!  real(rp)    :: elrhs(mnode)
!  real(rp)    :: xbocod(ndime,mgaus)

!  integer(ip) :: elidx(mnode), boidx(mnodb)
!  integer(ip) :: igaub, iboun, pgaub, pnodb, pblty
  !
  !
  !print *, sum( kfl_fixbo_tem(1:nboun), kfl_fixbo_tem(1:nboun)==2 )/2, & 
  !         sum( leset(1:nelem),  leset(1:nelem) == -1 )  
  !
  !
!  if(INOTMASTER) then
  !
!  boundaries: &
!  do iboun = 1,nboun
    !
!    h_flux(iboun) = 0.0
    !
!    if(kfl_fixbo_tem(iboun) == 2) then
      !
!      pblty = ltypb(iboun)
      !
!      tria03: & 
!      if(pblty == TRI03) then 
!        pnodb = nnode(pblty)
!        pgaub = ngaus(pblty) !< pgaub == 1
        !
!        boidx(1:pnodb)         = lnodb(1:pnodb,iboun)
        !
!        bocod(1:ndime,1:pnodb) = coord(1:ndime, boidx(1:pnodb) )
!        do igaub = 1,pgaub
!          xbocod(1:ndime,igaub) = matmul( bocod(1:ndime,1:pnodb), elmar(pblty) % shape(1:pnodb,igaub) ) 
!        enddo
        !
!        bprop(1:pnodb)  = prop( boidx(1:pnodb) ) 
!        do igaub = 1,pgaub
!          xbprop(igaub) = dot_product( bprop(1:pnodb), elmar(pblty) % shape(1:pnodb,igaub) )
!        enddo 
        !
!        h_flux(iboun) = sum( xbprop(1:pgaub) )/pgaub
        !
!      endif tria03
      !
      !ielem = lelbo(iboun)
      !pelty = ltype(ielem)
      !pnode = nnode(pelty)
      !pgaus = ngaus(pelty)
      !print *, leset(ielem), ielem
      !elidx(1:pnode)         = lnods(1:pnode,ielem)
      !elcod(1:ndime,1:pnode) = coord(1:ndime, elidx(1:pnode) )
      !
      !bprop(1:mnodb) = 0.0_rp 
      !if(leset(ielem) == -1) then
      !
!    endif 
    !
!  end do boundaries
!  endif
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !> @author  JM Zavala-Ake
  !> @date    
  !> @brief   
  !> @details 
  ! 
  !-----------------------------------------------------------------------||---!
  subroutine commdom_alya_cht_node_flux( prop ) !< ok 
!  use def_parame
!  use def_master
  use def_domain
!  use def_temper
!  use mod_memchk
!  use mod_postpr
!  use mod_gradie
  implicit none
  real(rp),   intent(in) :: prop(npoin) 
!  integer(ip)             :: ipoin,ibopo,idime
!  integer(4)              :: istat
!  real(rp), allocatable   :: gradt(:,:)
  !
  ! Allocate memory
  ! 
!  allocate( gradt(ndime,npoin), stat=istat)
!  call memchk(zero,istat,mem_modul(1:2,modul),'GRADT','tem_bouflux',gradt)
  !
  ! Compute temperature gradients
  ! 
!  call tem_heatfl(gradt)
 
  !
  ! Compute heat flux
  !
!  do ipoin = 1,npoin
!    prop(ipoin) = 0.0 
!    ibopo = lpoty(ipoin)
!    if(ibopo >= 1) then
!        do idime = 1,ndime
!          prop(ipoin) = prop(ipoin) + gradt(idime,ipoin) * exnor(idime,1,ibopo)
!        end do
!    endif 
!  end do
  !
  ! Deallocate memory
  !
!  call memchk(two,istat,mem_modul(1:2,modul),'GRADT','tem_outhfl',gradt)
!  deallocate(gradt,stat=istat)
!  if(istat/=0) call memerr(two,'GRADT','tem_outhfl',0_ip)
  !-----------------------------------------------------------------------||---!
  end subroutine commdom_alya_cht_node_flux
  !-----------------------------------------------------------------------||---!


!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !> @author  JM Zavala-Ake
  !> @date    
  !> @brief   
  !> @details 
  !
  ! call commdom_alya_cht_driver( icoup, 
  !                                coupling_type(current_code) % module_source, & 
  !                                coupling_type(icoup) % module_source, & 
  !                                IRECV, 
  !                                coupling_type(icoup))
  !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_alya_cht_driver(i_code,         &  
                                     active_module,  &   
                                     passive_module, &  
                                     sendrecv_type,  &
                                     passive_coupling, CPLNG)
!  use mod_commdom_alya, only: ISEND, IRECV, ISENDRECV
!  use mod_commdom_alya, only: commdom_alya_coupling_driver_i, commdom_alya_coupling_driver_j, commdom_alya_exchange
  use def_coupli,       only: typ_color_coupling
!  use def_temper, only: kfl_regim_tem
!  use def_nastin, only: kfl_regim_nsi
  implicit none
  integer(ip),  intent(in)  :: i_code
 integer(ip),  intent(in)  :: active_module, passive_module
  integer(ip),  intent(in)  :: sendrecv_type
  type(typ_color_coupling), optional, intent(in) :: passive_coupling
  type(COMMDOM_COUPLING), intent(in) :: CPLNG
  !
!  character(5) :: sendrecv_chosed = ''
!  integer(ip)   :: n_wets=-1, n_bnds=-1, i_wet, i_pts
  !
!  integer(ip), pointer :: wets_status(:)   => null() 
!  real(rp),    pointer :: wets_coords(:,:) => null() 
  !integer(ip) :: calculed=-1, sent=-1, received=-1, i_wets
  !
  !CPLNG%var_ij = -666.66
  !CPLNG%var_ji = -777.77
  !
  !-----------------------------------------------------------------------||---!
!  if( present(passive_coupling) ) then !.and. (active_module /= passive_module)) then 
!    n_bnds      =  passive_coupling % geome % nboun_wet
!    n_wets      =  passive_coupling % geome % npoin_wet
!    wets_status => passive_coupling % geome % status 
!    wets_coords => passive_coupling % geome % coord_wet
!  endif 
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!  
!  select case(sendrecv_type)   
!    case(ISEND)
!      sendrecv_chosed = 'ISEND'
!    case(IRECV)
!      sendrecv_chosed = 'IRECV'
!    case(ISENDRECV)
!      sendrecv_chosed = 'ISENDRECV'
!  end select 
  !
  !if(IMASTER) print*, "[active, passive]", active_module, passive_module, sendrecv_chosed
  !-----------------------------------------------------------------------||---!
  !
!  if(active_module == passive_module) then
    !
    !CPLNG%var_ij =  1.0
    !CPLNG%var_ji =  0.0 !-1.0/CPLNG%n_pts
    !
    !if(IMASTER) print*, "[module_i, module_i]", active_module, passive_module, "-->"
    !call commdom_alya_exchange(active_module, ISEND)
    !
!    if(active_module == CPLNG%module_i) then
!      call commdom_alya_exchange(CPLNG, CPLNG%module_i, sendrecv_type) 
!      if(INOTMASTER) print*, "[module_i, module_i]", sum(CPLNG%var_ji(1,:), dim=1), n_wets, CPLNG%n_pts
!    endif 
    !
!    if(active_module == CPLNG%module_j) then 
!      call commdom_alya_exchange(CPLNG, CPLNG%module_j, sendrecv_type) 
!      if(INOTMASTER) print*, "[module_j, module_j]", sum(CPLNG%var_ji(1,:), dim=1), n_wets, CPLNG%n_pts
!    endif 
    ! 
!  else&
!  if((active_module == CPLNG%module_i).and.(passive_module == CPLNG%module_j)) then
    !
    !CPLNG%var_ij = -2.0
    !CPLNG%var_ji =  0.0
    !
!    call commdom_alya_exchange(CPLNG, CPLNG%module_j, sendrecv_type) !
    !
    !if(IMASTER) print*, "[module_i, module_j]", active_module, passive_module, kfl_regim_tem, kfl_regim_nsi
    !if(INOTMASTER) print*, "[module_i, module_j]", sum(CPLNG%var_ij(1,:), dim=1), sum(CPLNG%var_ji(1,:), dim=1) 
!    if(INOTMASTER) print*, "[module_i, module_j]", sum(CPLNG%var_ji(1,:), dim=1), n_wets, CPLNG%n_pts
    !
    !if(INOTMASTER) then 
    !  !
    !  do i_wet = 1,n_wets
    !    if(wets_status(i_wet)>0) then 
    !  !    print *, wets_status(i_wet), wets_coords(:,i_wet), CPLNG%var_ji(1,wets_status(i_wet)) 
    !    endif 
    !  enddo 
    !  !
    !endif 
    !
!  else&
!  if((active_module == CPLNG%module_j).and.(passive_module == CPLNG%module_i)) then
    !
    !CPLNG%var_ij = -3.0
    !CPLNG%var_ji =  0.0
    !
!    call commdom_alya_exchange(CPLNG, CPLNG%module_i, sendrecv_type) ! 
    !
    !if(IMASTER) print*, "[module_j, module_i]", active_module, passive_module, kfl_regim_tem, kfl_regim_nsi
    !if(INOTMASTER) print*, "[module_j, module_i]", sum(CPLNG%var_ij(1,:), dim=1), sum(CPLNG%var_ji(1,:), dim=1) 
!    if(INOTMASTER) print*, "[module_j, module_i]", sum(CPLNG%var_ji(1,:), dim=1), n_wets, CPLNG%n_pts
    !
    !if(INOTMASTER) then 
    !  !
    !  do i_wet = 1,n_wets
    !    if(wets_status(i_wet)>0) then 
    !      print *, wets_status(i_wet), wets_coords(:,i_wet), CPLNG%var_ji(1,wets_status(i_wet)) 
    !    endif 
    !  enddo 
    !  !
    !endif 
    !
!  endif
  !
  !CPLNG%var_ij = -666.66
  !CPLNG%var_ji = -777.77
  !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  !
  ! +
  ! |_NASTIN+TEMPER
  !   |_ITASK_DOITER
  !                 \_ P, T, V  <- send 
  !
  ! +                _ P, T, V  <- recv
  ! |_NASTAL        /
  !   |_ITASK_DOITER
  !
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  subroutine commdom_alya_init_lowmach_compressible(CPLNG)
  implicit none
  type(COMMDOM_COUPLING), intent(in) :: CPLNG
  end subroutine
  !-----------------------------------------------------------------------||---!

  !=============================================================| contains |===!
end module mod_commdom_alya_cht
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
