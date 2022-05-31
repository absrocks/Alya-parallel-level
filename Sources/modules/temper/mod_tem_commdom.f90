!-------------------------------------------------------------------------||---!
!
!< 2014Dic10, choose a type of regime 
!< 2014Dic14, relax
!< 2015Jan29, current_what 
!< 2015Feb18, changes
!< 2015Abr14, conflic with 'OPTION: ZERO_FIXITY' 
!< 2017JAN07  
!
!-------------------------------------------------------------------------||---!
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  !
  ! +
  ! |_LOWMACH (NASTIN+TEMPER) 
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
!-------------------------------------------------------------------| INIT |---!
!-------------------------------------------------------------------------||---!
module mod_tem_commdom
  use def_parame,           only: ip, rp
  use def_master,           only: IMASTER, INOTMASTER, ISEQUEN
  use def_domain,           only: npoin, nboun
  use mod_commdom_alya_cht, only: CHT_CPLNG, commdom_cht_get_vals
  use mod_commdom_alya,     only: ISEND, IRECV
  use mod_commdom_alya,     only: commdom_alya_calculate_driver 
  use def_kintyp,           only: soltyp
  use def_master,           only: momod, modul
#ifdef COMMDOM
  use mod_commdom_plepp,    only: commdom_plepp_set_vals, commdom_plepp_check_fixno, commdom_plepp_reduce_sum
#endif
  implicit none
  type(soltyp), pointer :: solve(:)
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  private
    ! |_tem_concou                 _
    !   |_tem_coupli(ITASK_CONCOU)  |
    ! |_tem_begste                 _|  
    !   |_tem_coupli(ITASK_BEGSTE)  |
    ! |_tem_iniunk                 _|   
    !   |_tem_coupli(ITASK_INIUNK)  |
    !  _____________________________| 
    ! |_tem_coupling 
    !   |_case(ITASK_CONCOU|ITASK_BEGSTE|ITASK_INIUNK)  
    !public:: 
    !public:: tem_commdom_lm2_code_j
    !
    public:: commdom_alya_cht_node_flux, commdom_alya_cht_nodes2bvnat, temp2enth !< 2016Mar12    
    public:: tem_commdom_plugin 
    !
    interface tem_commdom_lm2_code_i
      module procedure tem_commdom_cht_lowmach 
    end interface
    public :: tem_commdom_lm2_code_i
    !
    interface tem_commdom_lm2_code_j
      module procedure tem_commdom_cht_temper 
    end interface
    public :: tem_commdom_lm2_code_j
    !
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!

  !=============================================================| contains |===!
contains
!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!

#include "tem_commdom.incl" 

  !-----------------------------------------------------------------------||---!
  subroutine tem_commdom_plugin()
#ifdef COMMDOM
#if   COMMDOM==-3 
  call tem_commdom_plugin3_parallelMPMD()
#endif 
#if   COMMDOM==3 
  call tem_commdom_plugin3()
#endif 
#if   COMMDOM==4 
  call tem_commdom_plugin4()
#endif 
#endif 
  end subroutine
  !-----------------------------------------------------------------------||---!

#ifdef COMMDOM
#if   COMMDOM==4  
  !-----------------------------------------------------------------------||---!
  subroutine tem_commdom_plugin4()
  use mod_commdom_alya,    only: COMMDOM_COUPLING
  use mod_commdom_driver,  only: CNT_SMS, CNT_CPLNG, CNT_SENDRECV, commdom_driver_exchange02
  use def_master,          only: title
  use def_domain,          only: ndime
  implicit none
  integer(ip) :: idime, ipoin, idof, icomp
  real(rp) :: d_relax = 1.0
  real(rp) :: n_relax = 1.0
  !
  real(rp) :: residual2(3,2) = 0.0
  !
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( any(CNT_SENDRECV) ) then
    !-----------------------------------------------------------| code_i==1 |---!
    code_i: &
    if(CNT_CPLNG%current_code==CNT_CPLNG%code_j) then                           !< 2016MAY09. j -> i  
       if( CNT_SENDRECV(9) ) then
        !-----------------------------------------------------------| U--> |---!
        to_send01: &
        if(inotmaster) then
          call commdom_alya_cht_node_flux( CNT_CPLNG%var_ij(ndime+1,1:npoin) )
          CNT_CPLNG%var_ij(ndime+1,1:npoin) = -CNT_CPLNG%var_ij(ndime+1,1:npoin)
        endif to_send01
        !-----------------------------------------------------------------||---!
        !
       !print *, "["//trim(title)//"] ", "TEMPER --> ", CNT_SMS
       !call commdom_driver_exchange02( CNT_CPLNG, debug=.false.)               !< 2016MAY09  
        !
        !--------------------------------------------------------| dUdn<-- |---!
        to_recv01: &
        if(inotmaster) then
          CNT_CPLNG%var_ji(ndime+1,1:npoin) = 0.0
        endif to_recv01
        !-----------------------------------------------------------------||---!
      endif
    endif code_i
    !---------------------------------------------------------------------||---!
    !----------------------------------------------------------| code_j==2 |---!
    code_j: &
    if(CNT_CPLNG%current_code==CNT_CPLNG%code_i) then                          !< 2016MAY09. i -> j
      if( CNT_SENDRECV(9) ) then                                        
        print *, "ERROR: ", "["//trim(title)//"]", CNT_SMS, "!!!"
        stop
        !--------------------------------------------------------| dUdn--> |---!
        to_send02: &
        if(inotmaster) then
        endif to_send02
        !-----------------------------------------------------------------||---!
        print *, "["//trim(title)//"] ", "TEMPER <-- ", CNT_SMS
      endif
    endif code_j
    !
  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!
#endif
#endif


  !-----------------------------------------------------------------------||---!
  !   +
  !   |_Alya                                       
  !     |_call Turnon()                            
  !     |_call Iniunk()                             
  !     |_time: do while
  !       |_call Timste()                         
  !       |_reset: do 
  !         |_call Begste()              TEMPE-->, HEATF<--0  
  !           |_block: do while                          
  !             |_coupling: do while                    
  !               |_call Begzon()        TEMPE-->, HEATF<--0  [sendrecv]
  !               |_modules: do while                           
  !                 |_call Doiter()                
  !                 |_call Concou()                
  !               |_call Endzon()        TEMPE<--, HEATF-->                  
  !                                                                           
  !             |_call Conblk()                             
  !       |_call Endste()                                     
  !   |_call Turnof()                    
  !
  !-----------------------------------------------------------------------||---!
  !
  ! COMMDOM==2
  !   
  !   uncoment 'commdom_nmodes_to_reaction'  kernel/master/inivar.f90 +184 
  !            'Temper'                      kernel/coupli/mod_coupling_driver.f90 +141
  !            'commdom_driver_init_cht'     kernel/domain/domain.f90 +68
  !
  !-----------------------------------------------------------------------||---!
  subroutine tem_commdom_plugin3()                                             !<  2016MAY11 
  use mod_commdom_alya,    only: COMMDOM_COUPLING
  use mod_commdom_driver,  only: CNT_SENDRECV, CNT_SMS
  use mod_commdom_driver,  only: CNT_CPLNG, commdom_driver_exchange02
  !
  use def_temper,          only: kfl_regim_tem, bvess_tem, bvnat_tem, kfl_fixno_tem
  use def_temper,          only: kfl_plepp_tem
  use def_master,          only: therm, title
  !
#ifdef COMMDOM
  use mod_commdom_dynamic, only: commdom_dynamic_check_fixno
  use mod_commdom_dynamic, only: commdom_dynamic_set_vals
#endif 
  !
  implicit none
  real(rp) :: relax_temp = 1.0 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( any(CNT_SENDRECV) ) then
    !
    !print *, "[tem_commdom_plugin]", CNT_SMS
    !
#ifdef COMMDOM
#if   COMMDOM==3  
    !----------------------------------------------------------| 2015Abr14 |---!
    if( CNT_CPLNG%current_what(2_ip) ) then 
      if(solve(1) % kfl_iffix /= 1) then 
        print *, "[mod_tem_commdom] '", trim(title) ,"' Set 'OPTION: FIXITY' before 'END_ALGEBRAIC_SOLVER', kfl_iffix:", solve(1) % kfl_iffix
        call runend("[mod_tem_commdom] ERROR!!") 
      endif 
    else
      if(solve(1) % kfl_iffix /= 0) then 
        print *, "[mod_tem_commdom] '", trim(title) ,"' Unset 'OPTION: ZERO_FIXITY' or 'FIXITY' before 'END_ALGEBRAIC_SOLVER', kfl_iffix:", solve(1) % kfl_iffix
        call runend("[mod_tem_commdom] ERROR!!") 
      endif 
    endif 
    !---------------------------------------------------------------------||---!
    !---------------------------------------------------------------------||---!
    if(CNT_CPLNG%current_code==CNT_CPLNG%code_i) then
      if( CNT_SENDRECV(7) ) then ! -ENDZON+
        !-------------------------------------------------------| HEATF--> |---!
        if(inotmaster) then 
          !
          CNT_CPLNG%var_ij(1,1:npoin) = 0.0
          if( CNT_CPLNG%current_what(2_ip) ) then 
            CNT_CPLNG%var_ij(1,1:npoin) = solve(1)%reaction(1,1:npoin) 
          else 
            call commdom_alya_cht_node_flux( CNT_CPLNG%var_ij(1,1:npoin) )
            CNT_CPLNG%var_ij(1,1:npoin) = -CNT_CPLNG%var_ij(1,1:npoin)
          endif 
          !
        endif 
        !-----------------------------------------------------------------||---!
        !
        call commdom_driver_exchange02( CNT_CPLNG )
        !
        !-------------------------------------------------------| TEMPE<-- |---!
        if(inotmaster) then 
          call temp2enth(  CNT_CPLNG%var_ji(1,1:npoin) )
          !
          !if( CNT_CPLNG%current_what(2_ip) ) relax_temp = 0.25 !< 2016Feb24  
          !
          call commdom_dynamic_set_vals( CNT_CPLNG%var_ji(1,1:npoin),     therm(  1:npoin,1), relax_op=relax_temp ) ! T(n)
          call commdom_dynamic_set_vals( CNT_CPLNG%var_ji(1,1:npoin),     therm(  1:npoin,2), relax_op=relax_temp ) ! x?
          call commdom_dynamic_set_vals( CNT_CPLNG%var_ji(1,1:npoin),     therm(  1:npoin,3), relax_op=relax_temp ) ! x?
          call commdom_dynamic_set_vals( CNT_CPLNG%var_ji(1,1:npoin), bvess_tem(1,1:npoin,1), relax_op=relax_temp ) ! T(n)
        endif
        !-----------------------------------------------------------------||---!
      endif
    endif
    !---------------------------------------------------------------------||---!
    if(CNT_CPLNG%current_code==CNT_CPLNG%code_j) then
      if( CNT_SENDRECV(7) ) then ! +BEGZON-
        !-------------------------------------------------------| TEMPE--> |---!
        if(inotmaster) then 
          if (kfl_regim_tem == 4) call tem_clippi()
          CNT_CPLNG%var_ij(1,1:npoin) = therm(1:npoin,1)
        endif
        !-----------------------------------------------------------------||---!
        !
        call commdom_driver_exchange02( CNT_CPLNG )
        !
        !-------------------------------------------------------| HEATF<-- |---!
        if(inotmaster) then 
          if( CNT_CPLNG%current_what(2_ip) ) then 
            call commdom_dynamic_set_vals( CNT_CPLNG%var_ji(1,1:npoin), solve(1)%bvnat(1,1:npoin), relax_op=1.0_rp ) 
          else 
            call commdom_alya_cht_nodes2bvnat( CNT_CPLNG%var_ji(1,1:npoin), bvnat_tem(3,1:nboun,1) )
          endif 
          CNT_CPLNG%var_ji(1,1:npoin) = 0.0_rp                              ! <--- reset values for the next step...
        endif 
        !-----------------------------------------------------------------||---!
      endif
    endif    
    !---------------------------------------------------------------------||---!
#endif
#endif 
    !
  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !   +
  !   |_Alya                                       
  !     |_call Turnon()                            
  !     |_call Iniunk()                             
  !     |_time: do while
  !       |_call Timste()                         
  !       |_reset: do 
  !         |_call Begste()              TEMPE-->, HEATF<--0  
  !           |_block: do while                          
  !             |_coupling: do while                    
  !               |_call Begzon()        TEMPE-->, HEATF<--0  [sendrecv]
  !               |_modules: do while                           
  !                 |_call Doiter()                
  !                 |_call Concou()                
  !               |_call Endzon()        TEMPE<--, HEATF-->                  
  !                                                                           
  !             |_call Conblk()                             
  !       |_call Endste()                                     
  !   |_call Turnof()                    
  !
  !-----------------------------------------------------------------------||---!
  subroutine tem_commdom_plugin_xxx()
  use mod_commdom_alya,    only: COMMDOM_COUPLING
  use mod_commdom_driver,  only: CNT_SENDRECV, CNT_SMS
  use mod_commdom_driver,  only: CNT_CPLNG
  !
  use def_temper,          only: kfl_regim_tem, bvess_tem, bvnat_tem, kfl_fixno_tem
  use def_temper,          only: kfl_plepp_tem
  use def_master,          only: therm
  !
#ifdef COMMDOM
  use mod_commdom_dynamic, only: commdom_dynamic_check_fixno
  use mod_commdom_dynamic, only: commdom_dynamic_set_vals
#endif 
  !
  implicit none
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( any(CNT_SENDRECV) ) then
    !
    print *, "[tem_commdom_plugin]", CNT_SMS
    !
#ifdef COMMDOM
#if   COMMDOM==2 
    !
    if(CNT_CPLNG%current_code==CNT_CPLNG%code_i) then
      if( inotmaster.and.( CNT_SENDRECV(4).or.CNT_SENDRECV(2) ) ) then
        !---------------------------------------------------------| TEMPE--> |---!
        if (kfl_regim_tem == 4) call tem_clippi()
        CNT_CPLNG%var_ij(1,1:npoin) = therm(1:npoin,1)
        !---------------------------------------------------------| HEATF<-- |---!
        call commdom_alya_cht_nodes2bvnat( CNT_CPLNG%var_ji(1_ip,1_ip:npoin), bvnat_tem(3,1:nboun,1) )
        CNT_CPLNG%var_ji(1_ip,1_ip:npoin) = 0.0_rp                              ! <--- reset values for the next step...
        !-------------------------------------------------------------------||---!
      endif
    endif
    if(CNT_CPLNG%current_code==CNT_CPLNG%code_j) then
      if( inotmaster.and.CNT_SENDRECV(5) ) then
        !---------------------------------------------------------| TEMPE<-- |---!
        call temp2enth(CNT_CPLNG%var_ji(1,1:npoin))
        call commdom_dynamic_set_vals( CNT_CPLNG%var_ji(1,1:npoin),       therm(1:npoin,1) ) !<--- T(n)
        call commdom_dynamic_set_vals( CNT_CPLNG%var_ji(1,1:npoin),       therm(1:npoin,2) ) ! x?
        call commdom_dynamic_set_vals( CNT_CPLNG%var_ji(1,1:npoin),       therm(1:npoin,3) ) ! x?
        call commdom_dynamic_set_vals( CNT_CPLNG%var_ji(1,1:npoin), bvess_tem(1,1:npoin,1) ) !<--- T(n)
        !---------------------------------------------------------| HEATF--> |---!
        CNT_CPLNG%var_ij(1,1:npoin) = 0.0
        call commdom_alya_cht_node_flux( CNT_CPLNG%var_ij(1,1:npoin) )
        CNT_CPLNG%var_ij(1,1:npoin) = -CNT_CPLNG%var_ij(1,1:npoin)
        !-------------------------------------------------------------------||---!
      endif
    endif
    ! 
#endif
#endif 
    !
  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!

!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine tem_commdom_cht_lowmach(sendrecv_code, when)
  use def_temper,           only: kfl_regim_tem, bvess_tem, bvnat_tem, kfl_fixno_tem
  use def_temper,           only: kfl_plepp_tem
  use def_master,           only: therm
  use def_domain,           only: npoin, nboun
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  implicit none
  integer(ip), intent(in) :: sendrecv_code, when
  logical(ip) :: set_var_ij, get_var_ji, ini_var_ij 
  logical(ip) :: is_lowmach
  real(rp)    :: relax_fac = 0.5, cploc(6,2),dummr,tenew
  integer(ip) :: ipoin,ivalu
  real(rp)    :: prop_in=1.0, prop_out=-1.0
  solve => momod(modul) % solve(1:)
#ifdef COMMDOM 
#if COMMDOM==1
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  !
  if(kfl_plepp_tem < 0) then
    print *, "[mod_tem_commdom] set 'PLEPP: CONDU' before 'END_NUMERICAL_TREATMENT' in order to choose one coupling regime"
    call runend("[mod_tem_commdom] ERROR!!") 
  endif
  !
  if(solve(1) % kfl_iffix /= 0) then 
    print *, "[mod_tem_commdom] unset 'OPTION: ZERO_FIXITY' before 'END_ALGEBRAIC_SOLVER' ", "kfl_iffix:", solve(1) % kfl_iffix /= 0
    call runend("[mod_tem_commdom] ERROR!!") 
  endif 
  !
  is_lowmach = kfl_regim_tem==kfl_plepp_tem                                    !<-- TEMPER-LOWMACH 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  call commdom_alya_calculate_driver(CHT_CPLNG) 
  ini_var_ij = CHT_CPLNG%setgetvar(CHT_CPLNG%code_i,1).and.is_lowmach
  set_var_ij = CHT_CPLNG%setgetvar(CHT_CPLNG%code_i,2).and.is_lowmach
  get_var_ji = CHT_CPLNG%setgetvar(CHT_CPLNG%code_i,3).and.is_lowmach
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if(INOTMASTER) then
    if(ini_var_ij) then  !< ITASK_INIUNK -->
      !-----------------------------------------------------------------||---!
      !                     Dirichlet nodes : fixval==0 ---V
      call commdom_plepp_check_fixno(kfl_fixno_tem, 1_ip, 1_ip, .True.) ! fixno, idofn, fixval, ToDo
      !-----------------------------------------------------------------||---!
      CHT_CPLNG%var_ij(1_ip,1_ip:npoin) = 0.0
      if( CHT_CPLNG%current_what(2_ip) ) then ! what_j == -RESIDUAL 
        CHT_CPLNG%var_ij(1_ip,1_ip:npoin) = -solve(1)%reaction(1_ip,1_ip:npoin) 
      else
        call commdom_alya_cht_node_flux( CHT_CPLNG%var_ij(1_ip,1_ip:npoin) ) 
      endif 
      CHT_CPLNG%var_ij(1_ip,1_ip:npoin) = -CHT_CPLNG%var_ij(1_ip,1_ip:npoin)
      !-----------------------------------------------------------------||---!
      else &
    if(get_var_ji) then !< ITASK_BEGSTE <--
      !-----------------------------------------------------------------||---!
!      call commdom_plepp_set_vals( CHT_CPLNG%var_ji(1,1:npoin),       tempe(1:npoin,1)  )  \
!      call commdom_plepp_set_vals( CHT_CPLNG%var_ji(1,1:npoin),       tempe(1:npoin,2)  )  |___2015Feb18
!      call commdom_plepp_set_vals( CHT_CPLNG%var_ji(1,1:npoin),       tempe(1:npoin,3)  )  |
!      call commdom_plepp_set_vals( CHT_CPLNG%var_ji(1,1:npoin), bvess_tem(1,1:npoin,1)  )  /
      !-----------------------------------------------------------------||---!
    else &  
    if(set_var_ij) then !< ITASK_CONCOU -->
      !-----------------------------------------------------------------||---!
      CHT_CPLNG%var_ij(1_ip,1_ip:npoin) = 0.0
      if( CHT_CPLNG%current_what(2_ip) ) then ! what_j == -RESIDUAL 
        CHT_CPLNG%var_ij(1_ip,1_ip:npoin) = -solve(1)%reaction(1_ip,1_ip:npoin) 
      else
        call commdom_alya_cht_node_flux( CHT_CPLNG%var_ij(1_ip,1_ip:npoin) ) 
      endif  
      CHT_CPLNG%var_ij(1_ip,1_ip:npoin) = -CHT_CPLNG%var_ij(1_ip,1_ip:npoin)
      !-----------------------------------------------------------------||---!
      call temp2enth(CHT_CPLNG%var_ji(1,1:npoin))

      call commdom_plepp_set_vals( CHT_CPLNG%var_ji(1,1:npoin),       therm (1:npoin,1) ) !<--- T(n)
      call commdom_plepp_set_vals( CHT_CPLNG%var_ji(1,1:npoin),       therm(1:npoin,2) ) ! x
      call commdom_plepp_set_vals( CHT_CPLNG%var_ji(1,1:npoin),       therm(1:npoin,3) ) ! x
      call commdom_plepp_set_vals( CHT_CPLNG%var_ji(1,1:npoin), bvess_tem(1,1:npoin,1) ) !<--- T(n)
      !!! call commdom_plepp_set_vals( CHT_CPLNG%var_ji(1,1:npoin), bvess_tem(1,1:npoin,2) ) !<--- T(n)
      !-----------------------------------------------------------------||---!
    endif
    !-------------------------------------------------------------------||---!
  endif ! INOTMASTER
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  !if(set_var_ij) call commdom_plepp_reduce_sum(prop_in, prop_out)
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
#endif
#endif
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  subroutine tem_commdom_cht_temper(sendrecv_code, when)
  use def_temper,           only: kfl_regim_tem, bvess_tem, bvnat_tem, kfl_fixno_tem
  use def_master,           only: therm  
  use def_domain,           only: npoin, nboun
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  implicit none
  integer(ip), intent(in) :: sendrecv_code, when
  logical(ip) :: set_var_ij, get_var_ji, ini_var_ij
  logical(ip) :: is_temper
  real(rp)    :: prop_in=1.0, prop_out=-1.0
  solve => momod(modul) % solve(1:)
#ifdef COMMDOM 
#if COMMDOM==1
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  call commdom_alya_calculate_driver(CHT_CPLNG)
  is_temper  = kfl_regim_tem==0
  ini_var_ij = CHT_CPLNG%setgetvar(CHT_CPLNG%code_j,1).and.is_temper 
  set_var_ij = CHT_CPLNG%setgetvar(CHT_CPLNG%code_j,2).and.is_temper
  get_var_ji = CHT_CPLNG%setgetvar(CHT_CPLNG%code_j,3).and.is_temper 
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if(INOTMASTER) then
    if(ini_var_ij) then !< ITASK_INIUNK --> 
      !-------------------------------------------------------------------||---!
      !                 Free nodes (Neumann): fixval==0 ---V
      call commdom_plepp_check_fixno(kfl_fixno_tem, 1_ip, 0_ip, .True.) 
      !-------------------------------------------------------------------||---!
      if (kfl_regim_tem == 4) call tem_clippi()
      CHT_CPLNG%var_ij(1,1:npoin) = therm(1:npoin,1)                           ! T(0) --> 
      !-------------------------------------------------------------------||---!
    else &
    if(get_var_ji) then !< ITASK_BEGSTE <--
      !-------------------------------------------------------------------||---!
      if (kfl_regim_tem == 4) call tem_clippi()
      CHT_CPLNG%var_ij(1,1:npoin) = therm(1:npoin,1)                           ! T(n) -->
      !-------------------------------------------------------------------||---!s
      if( CHT_CPLNG%current_what(2_ip) ) then ! what_j == -RESIDUAL 
        call commdom_plepp_set_vals( CHT_CPLNG%var_ji(1,1:npoin), solve(1)%bvnat(1,1:npoin), relax_in=1.0_rp) 
      else
        call commdom_alya_cht_nodes2bvnat( CHT_CPLNG%var_ji(1_ip,1_ip:npoin), bvnat_tem(3,1:nboun,1) )
      endif 
      !-------------------------------------------------------------------||---!
!      if( CHT_CPLNG%current_what(2_ip) ) then ! what_j == -RESIDUAL                                      \
!        solve(1)%bvnat(1_ip,1_ip:npoin) = 0.0_rp                                                         |
!        solve(1)%bvnat(1_ip,1_ip:npoin) = CHT_CPLNG%var_ji(1_ip,1_ip:npoin)                              |
!      else                                                                                               |___2015Feb18
!        bvnat_tem(:,1:nboun,:) = 0.0                                                                     |
!        call commdom_alya_cht_nodes2bvnat( CHT_CPLNG%var_ji(1_ip,1_ip:npoin), bvnat_tem(3,1:nboun,1) )   |
!      endif                                                                                              /
      !-------------------------------------------------------------------||---!
    else & 
    if(set_var_ij) then !< ITASK_CONCOU --> 
      !-------------------------------------------------------------------||---! \
!      CHT_CPLNG%var_ij(1,1:npoin) = tempe(1:npoin,1)                            |___2015Feb18
      !-------------------------------------------------------------------||---! /
    endif
  endif ! INOTMASTER
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  !if(get_var_ji) call commdom_plepp_reduce_sum(prop_in, prop_out)
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
#endif
#endif
  end subroutine
  !-----------------------------------------------------------------------||---!

!-------------------------------------------------------------------------||---!
!----------------------------------------------------------------| PRIVATE |---!
!-------------------------------------------------------------------------||---!



  !-----------------------------------------------------------------------||---!
  !> @author Dani mira 
  !> @date    
  !> @brief   
  !> @details 
  !
  !-----------------------------------------------------------------------||---!
  subroutine temp2enth(temperature) !< ok  
  use def_temper,           only: kfl_regim_tem, kfl_fixno_tem
  use def_master,           only: sphec  
  use def_domain,           only: npoin, nboun 
  use mod_physics,          only: physics_T_2_HCp
  implicit none 
  real(rp), intent(out) :: temperature(npoin) 
  real(rp) :: dummr = 0.0_rp, cploc(6,2) = 0.0_rp, hnew=0.0_rp  
  integer(ip) :: ipoin  
 
  if(kfl_regim_tem==4) then
     do ipoin=1,npoin
        if(kfl_fixno_tem(1,ipoin)==1) then
           cploc(1:6,1:2) = sphec(ipoin,1:6,1:2)
           call physics_T_2_HCp(temperature(ipoin),cploc,hnew,dummr)
           temperature(ipoin) = hnew 
        endif 
     enddo 
  end if
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
  subroutine commdom_alya_cht_nodes2bvnat(prop, h_flux) !< ok  
  use def_parame
  use def_master
  use def_domain
  use def_temper
  use def_elmtyp, only: TRI03, TRI06, TET04, HEX08 
  implicit none
  real(rp),   intent( in   ) :: prop(npoin) 
  real(rp),   intent( inout) :: h_flux(nboun) 
  !
  real(rp)    :: bprop(mnodb)
  real(rp)    :: xbprop(mgaus) 
  real(rp)    :: gbsur
  real(rp)    :: bocod(ndime,mnodb)
  real(rp)    :: elrhs(mnode)
  real(rp)    :: xbocod(ndime,mgaus)
  integer(ip) :: elidx(mnode), boidx(mnodb)
  integer(ip) :: igaub, iboun, pgaub, pnodb, pblty
  integer(ip) :: n_fixbo 
  !
  !
  !print *, sum( kfl_fixbo_tem(1:nboun), kfl_fixbo_tem(1:nboun)==2 )/2, & 
  !         sum( leset(1:nelem),  leset(1:nelem) == -1 )  
  !
  !
  if(INOTMASTER) then
  !
  boundaries: &
  do iboun = 1,nboun
    !
    h_flux(iboun) = 0.0
    !
    if(kfl_fixbo_tem(iboun) == 2) then
      !
      pblty = ltypb(iboun)
      !
!      tria03: & 
!      if(pblty == TRI03) then 
        pnodb = nnode(pblty)
        pgaub = ngaus(pblty) !< pgaub == 1
        !
        boidx(1:pnodb)         = lnodb(1:pnodb,iboun)
        !
        bocod(1:ndime,1:pnodb) = coord(1:ndime, boidx(1:pnodb) )
        do igaub = 1,pgaub
          xbocod(1:ndime,igaub) = matmul( bocod(1:ndime,1:pnodb), elmar(pblty) % shape(1:pnodb,igaub) ) 
        enddo
        !
        bprop(1:pnodb)  = prop( boidx(1:pnodb) ) 
        do igaub = 1,pgaub
          xbprop(igaub) = dot_product( bprop(1:pnodb), elmar(pblty) % shape(1:pnodb,igaub) )
        enddo 
        !
        h_flux(iboun) = sum( xbprop(1:pgaub) )/pgaub
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
    endif 
    !
  end do boundaries
  endif
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
  use def_parame
  use def_master
  use def_domain
  use def_temper
  use mod_memchk
  use mod_postpr
  use mod_gradie
  implicit none
  real(rp), intent(inout) :: prop(npoin) 
  integer(ip)             :: ipoin,ibopo,idime
  integer(4)              :: istat
  real(rp), allocatable   :: gradt(:,:)
  !
  ! Allocate memory
  allocate( gradt(ndime,npoin), stat=istat)
  call memchk(zero,istat,mem_modul(1:2,modul),'GRADT','tem_bouflux',gradt)
  !
  ! Compute temperature gradients
  call tem_heatfl( gradt )

  ! Compute heat flux
  do ipoin = 1,npoin
    prop(ipoin) = 0.0 
    ibopo = lpoty(ipoin)
    if(ibopo >= 1) then
        do idime = 1,ndime
          prop(ipoin) = prop(ipoin) + gradt(idime,ipoin) * exnor(idime,1,ibopo)
        end do
    endif 
  end do
  !
  ! Deallocate memory
  !
  call memchk(two,istat,mem_modul(1:2,modul),'GRADT','tem_outhfl',gradt)
  deallocate(gradt,stat=istat)
  if(istat/=0) call memerr(two,'GRADT','tem_outhfl',0_ip)
  !-----------------------------------------------------------------------||---!
  end subroutine commdom_alya_cht_node_flux
  !-----------------------------------------------------------------------||---!

  !=============================================================| contains |===!
end module mod_tem_commdom 
!==============================================================================!
!==============================================================================!
