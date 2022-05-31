subroutine chm_matrix()
  !-----------------------------------------------------------------------
  !****f* Chemic/chm_matrix
  ! NAME 
  !    chm_matrix
  ! DESCRIPTION
  !    This routine assembles matrix and RHS
  ! USES
  !    chm_matrix
  !    Soldir
  !    Solite
  ! USED BY
  !    chm_solite
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_chemic
  use def_domain
  use def_solver
  use mod_postpr
  use mod_messages, only : livinf
  use mod_timings,  only : timings_assembly
  implicit none
  real(rp)    :: time1,time2,time3,time4,dtinv_tmp,time_elem,time_boun
  !
  ! Solver initializations and preconditioners
  !
  call inisol()

  solve_sol(1) % xdiag = 1.0_rp/dtinv_chm
  if ( kfl_model_chm /= 4 ) then
     call chm_massma(pmatr) !!!!!!!!!!!!! IS THIS ONLY FOR IMPLICIT?
  end if
  
  call cputim(time1)
  !
  ! Element assembly
  !
  if( INOTMASTER ) then

     call chm_elmope_all(1_ip)

  else
     dtinv_chm = 1.0e6_rp
  end if
  call cputim(time2)

  !
  ! Boundary assembly
  !
  call cputim(time3)
  if( INOTMASTER ) call chm_bouope() 
  call cputim(time4)
  
  !
  ! Split assembly
  !
  if( solve_sol(1)%kfl_algso == 9 .or. &
      solve_sol(1)%kfl_algso == 10 ) then
     dtinv_tmp = dtinv_chm
     dtinv_chm = 0.0_rp     
  end if
  if( solve_sol(1)%kfl_algso == 9 .or. &
      solve_sol(1)%kfl_algso == 10 ) then
     dtinv_chm = dtinv_tmp
  end if

  time_elem = time2-time1
  time_boun = time4-time3
  call timings_assembly(time_elem,time_boun,TYPE_OF_ASSEMBLY='ELEMENT, BOUNDARY')
  cputi_chm(1) = cputi_chm(1) + time_elem
  cputi_chm(2) = cputi_chm(2) + time_boun

end subroutine chm_matrix

