subroutine ale_doiter
  !-----------------------------------------------------------------------
  !****f* Alefor/ale_doiter
  ! NAME 
  !    ale_doiter
  ! DESCRIPTION
  !    Check if ALEFOR should be solved.
  !    It is not solved if:
  !    - All d.o.f. are prescribed
  !    - All prescribed nodes are prescribed to zero
  ! USES
  !    ale_begite
  !    ale_solite
  !    ale_endite
  ! USED BY
  !    Aelfor
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_solver
  use def_alefor
  implicit none
  integer(ip) :: ipoin,idime,kfl_fmale_ale,kfl_bvess_ale ! ,kfl_solve_ale
    
  kfl_fmale_ale = 0_ip ! kfl_fmale_ale initialization
  kfl_solve_ale = 1_ip ! kfl_solve_ale initialization, alefor will always solve when called
  kfl_bvess_ale = 1_ip ! kfl_bvess_ale_initialization, always use the boundary conditions in bvess_ale


  !------------------------------------------------------
  !CALL RIGID BODY
  !
  !
  if (kfl_rigid_ale /= 0_ip) call ale_solrbo()
  !_________________________________________________________


  !-------------------------------------------------------------------
  !
  ! Solve mesh deformation and smoothing
  !
  !-------------------------------------------------------------------
  
  call ale_begite()
  call ale_smodef()
  
  !-------------------------------------------------------------------
  !
  ! Recompute some domain variables and output domain mesh
  !
  !-------------------------------------------------------------------
  
  if( kfl_fmale_ale == 0_ip ) kfl_domar = 1_ip


  !-------------------------------------------------------------------
  !
  ! End of iterations
  !
  !-------------------------------------------------------------------
  
  call ale_endite()
  
!   print *, "ppp111", dispm(:,1,1)
!   print *, "ppp222", dispm(:,2,1)
!   print *, "ppp333", dispm(:,3,1)
!   print *, "ppp444", dispm(:,4,1)

end subroutine ale_doiter

 
