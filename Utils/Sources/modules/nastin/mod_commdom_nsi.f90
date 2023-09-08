!
!< 2016MAR29 -> created 
!< 2016MAR30 
!
!==============================================================================!
  !-----------------------------------------------------------------------||---!
  !-----------------------------------------------------------| ITERATIONS |---!
  !   +
  !   |_Alya                                       
  !     |_call Turnon()                            
  !     |_call Iniunk()                             
  !     |_time: do while
  !       |_call Timste()                         
  !       |_reset: do                            [i]             [j]
  !         |_call Begste()                                      
  !           |_block: do while          
  !             |_coupling: do while     
  !               |_call Begzon()         (7)   [ +r-]          [ +r-]   (3) 
  !               |_modules: do while     
  !                 |_call Doiter()       (8)   [+dU-]          [ +U-]   (4) 
  !                 |_call Concou()       (1)   [ +U-]          [+dU-]   (5) 
  !               |_call Endzon()         (2)   [ -s+]          [ -s+]   (6) 
  !             |_call Conblk()                             
  !       |_call Endste()                                     
  !   |_call Turnof()                    
  !
  !-----------------------------------------------------------------------||---!
!==============================================================================!
module mod_commdom_nsi 
  use def_parame,           only: ip, rp
  use def_master,           only: IMASTER, INOTMASTER, ISEQUEN, title 
  use def_domain,           only: npoin, nboun, ndime, coord
  !
  use def_kintyp,           only: soltyp
  use def_master,           only: momod, modul
  !
  use mod_commdom_alya,     only: INONE
#ifdef COMMDOM
  use mod_commdom_dynamic,  only: commdom_dynamic_check_fixno
  use mod_commdom_dynamic,  only: commdom_dynamic_set_vals
  use mod_commdom_dynamic,  only: commdom_dynamic_reduce_sum
#endif 
  !
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  implicit none
  !
  !type(soltyp), pointer :: solve(:)
  private
  public :: commdom_nsi_plugin

  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
contains
!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  subroutine commdom_nsi_plugin()
#ifdef COMMDOM
#if    COMMDOM==4 
  call commdom_nsi_plugin4() 
#endif 
#endif 
  end subroutine
  !-----------------------------------------------------------------------||---!


#ifdef COMMDOM
#if   COMMDOM==4  
  !-----------------------------------------------------------------------||---!
  subroutine commdom_nsi_plugin4()
  use mod_commdom_alya,    only: COMMDOM_COUPLING
  use mod_commdom_driver,  only: CNT_SMS, CNT_CPLNG, CNT_SENDRECV, commdom_driver_exchange02
  use def_master,          only: displ, unkno
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
    if(CNT_CPLNG%current_code==CNT_CPLNG%code_i) then
      if( CNT_SENDRECV(7) ) then
        print *, "ERROR: ", "["//trim(title)//"]", CNT_SMS, "!"
        stop
        !-----------------------------------------------------------| U--> |---!
        to_send01: &
        if(inotmaster) then
          CNT_CPLNG%var_ij(1:ndime,1:npoin) = 0.0
        endif to_send01 
        !-----------------------------------------------------------------||---!
        !
        call commdom_driver_exchange02( CNT_CPLNG, debug=.false.)               !< 
        !
        !--------------------------------------------------------| dUdn<-- |---!
        to_recv01: &
        if(inotmaster) then
          CNT_CPLNG%var_ji(1:ndime,1:npoin) = 0.0
        endif to_recv01 
        !-----------------------------------------------------------------||---!
      endif
    endif code_i
    !---------------------------------------------------------------------||---!
    !----------------------------------------------------------| code_j==2 |---!
    code_j: & 
    if(CNT_CPLNG%current_code==CNT_CPLNG%code_j) then
      if( CNT_SENDRECV(7) ) then                                               !< 2016MAR30. 8 -> 7  
        !--------------------------------------------------------| dUdn--> |---! ! (1) NASTIN --> REACT --> SOLID  
        to_send02: &
        if(inotmaster) then
         !CNT_CPLNG%var_ij(1:ndime,1:npoin) = 1.0
          call commdom_nsi_get_reaction( CNT_CPLNG%var_ij(1:ndime,1:npoin), algebraic=.true. ) !< 2016MAR30  
        endif to_send02 
        !-----------------------------------------------------------------||---!
        !
        print *, "["//trim(title)//"] ", "NASTIN --> REACT", CNT_SMS
        call commdom_driver_exchange02( CNT_CPLNG , debug=.true.)              !< 2016MAR30
        !
        !-----------------------------------------------------------| U<-- |---!
        to_recv02: &
        if(inotmaster) then
        endif to_recv02   
        !-----------------------------------------------------------------||---!
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
  subroutine commdom_nsi_get_reaction( var_ij, algebraic )
  implicit none
  real(rp),               intent(inout) :: var_ij(ndime,npoin)
  logical(ip), optional,  intent(in   ) :: algebraic
  !
  type(soltyp), pointer :: solve_sol(:)
  solve_sol => momod(modul) % solve(1_ip:)
  !--------------------------------------------------------------| dUdn--> |---!
  if(inotmaster) then
    if( (solve_sol(1)%kfl_bvnat/=1).and.(solve_sol(1)%kfl_react/=1) ) then 
      print *, " ERROR: [commdom_nsi_get_reaction]"
      print *, "kfl_bvnat, kfl_react: ", solve_sol(1)%kfl_bvnat, solve_sol(1)%kfl_react 
      call runend(" kfl_bvnat/=1 .and. kfl_react/=1  !!")
    endif 
    !
    if( (solve_sol(1)%kfl_bvnat/=0).and.(solve_sol(1)%kfl_react/=1) ) then
      print *, " ERROR: [commdom_nsi_get_reaction]"
      print *, " DIRIC <-> kfl_react==1 ", solve_sol(1)%kfl_react 
      call runend(" kfl_bvnat/=0 .and. kfl_react/=1  !!")
    endif 
    !
    var_ij(1:ndime,1:npoin) = 0.0_rp
    if( present(algebraic).and.algebraic ) then
     ! 
     !if(solve(1)%kfl_bvnat == 1) then
     !solve_sol(1) % block_array(2) % bvnat
     !
      var_ij(1:ndime,1:npoin) = solve_sol(1)%reaction(1:ndime,1:npoin)
    else
      !
    endif
    !
  endif
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!



!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!
end module mod_commdom_nsi 
!==============================================================================!
!==============================================================================!
