subroutine soldir(rhsid,unkno,amatr)

  !-----------------------------------------------------------------------
  !
  ! This routine solves linear systems using a direct method. It is
  ! assumed that the mesh graph has been properly renumbered.      
  !
  !-----------------------------------------------------------------------
  use def_solver
  use def_master, only    :  IPARALL
  implicit none

  real(rp), intent(in)    :: amatr(*)
  real(rp), intent(out)   :: unkno(*)
  real(rp), intent(inout) :: rhsid(*)
  integer(ip)             :: ntotv
  real(rp)                :: cpu_refe1,cpu_refe2
  integer(ip), save       :: ipass=0

  if( IPARALL ) then
     call runend('SOLDIR: DIRECT SOLVER NOT VALID IN PARALLEL')
  end if
 !if(solve_sol(1) % kfl_full_rows == 1 ) then
 ! ia => solve_sol(1) % ia_full
 ! ja => solve_sol(1) % ja_full
 !else
 ! ia => solve_sol % ia
 ! ja => solve_sol % ja
 !end if
  !
  ! Initializations
  !
  call cputim(cpu_refe1)
  if(ipass==0) then
     smemo(1) = 0
     smemo(2) = 0
     ipass    = 1
  end if
  ntotv=solve_sol(1)%ndofn*solve_sol(1)%nequa
  !
  ! Call solver routine
  !      
  call skylin(ntotv,amatr,rhsid,unkno)
  !
  ! Solver statistics
  !
  solve_sol(1)%nsolv = solve_sol(1)%nsolv + 1
  !
  ! Compute CPU time 
  !
  call cputim(cpu_refe2)

end subroutine soldir
