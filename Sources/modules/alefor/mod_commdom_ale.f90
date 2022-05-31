!
! < !< 2016MAR29 -> created 
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
module mod_commdom_ale 
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
  type(soltyp), pointer :: solve(:)
  private
  public :: commdom_ale_plugin

  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
contains
!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  subroutine commdom_ale_plugin()
#ifdef COMMDOM
#if   COMMDOM==4 
  call commdom_ale_plugin4()
#endif 
#endif 
  end subroutine
  !-----------------------------------------------------------------------||---!


#ifdef COMMDOM
#if   COMMDOM==4  
  !-----------------------------------------------------------------------||---!
  subroutine commdom_ale_plugin4()
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
!print *, "["//trim(title)//"] ->", CNT_SMS

    !-----------------------------------------------------------| code_i==1 |---!
    code_i: &
    if(CNT_CPLNG%current_code==CNT_CPLNG%code_i) then
      if( CNT_SENDRECV(7) ) then
        print *, "ERROR: ", "["//trim(title)//"]", CNT_SMS, "!!!"
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
         !call commdom_dynamic_set_vals( CNT_CPLNG%var_ji(idime,1:npoin),     veloc(idime,1:npoin,1), relax_op=n_relax, res2=residual2(3,1:2), debug=.false. )
          CNT_CPLNG%var_ji(1:ndime,1:npoin) = 0.0
        endif to_recv01 
        !-----------------------------------------------------------------||---!
      endif
    endif code_i
    !---------------------------------------------------------------------||---!
    !----------------------------------------------------------| code_j==2 |---!
    code_j: & 
    if(CNT_CPLNG%current_code==CNT_CPLNG%code_j) then
      if( CNT_SENDRECV(8) ) then                                               !< 2016MAR30. 7 -> 8  
        !--------------------------------------------------------| dUdn--> |---!
        to_send02: &
        if(inotmaster) then
        !  CNT_CPLNG%var_ij(1:ndime,1:npoin) = 0.0
        endif to_send02 
        !-----------------------------------------------------------------||---!
        !
        print *, "["//trim(title)//"] ", "ALEFOR <-- UNKNO", CNT_SMS
       !call commdom_driver_exchange02( CNT_CPLNG , debug=.true.)              !< 2016MAR30  
        !
        !-----------------------------------------------------------| U<-- |---! ! (3) ALEFOR <-- UNKNO <-- SOLID  
        call commdom_ale_get_unkno( CNT_CPLNG%var_ji(1:ndime,1:npoin), relax=1.0_rp ) 
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

!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!

  !-----------------------------------------------------------------------||---!
  subroutine commdom_ale_get_unkno( var_ji, relax )
  use def_master,        only :  bvess_ale
  implicit none
  real(rp),               intent(inout) :: var_ji(ndime,npoin)
  real(rp),               intent(in   ) :: relax
  integer(ip) :: idime 
#ifdef COMMDOM
  !--------------------------------------------------------------| dUdn--> |---!
  if(inotmaster) then
    !
    do idime = 1,ndime 
      call commdom_dynamic_set_vals( var_ji(idime,1:npoin), bvess_ale(idime,1:npoin,1), relax_op=relax, debug=.false. )
    enddo  
    var_ji(1:ndime,1:npoin) = 0.0  
    !
  endif
#endif 
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!



!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!
end module mod_commdom_ale 
!==============================================================================!
!==============================================================================!
