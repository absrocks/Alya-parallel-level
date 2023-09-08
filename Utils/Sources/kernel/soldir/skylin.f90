subroutine skylin(ntotv,amatr,rhsid,unkno)

  !-----------------------------------------------------------------------
  !
  ! This routine solves the set of linear equations 
  !
  !-----------------------------------------------------------------------
  use def_parame
  use def_solver
  implicit none
  integer(ip), intent(in)    :: ntotv
  real(rp),    intent(in)    :: amatr(*)
  real(rp),    intent(inout) :: rhsid(ntotv,solve_sol(1)%nrhss)
  real(rp),    intent(out)   :: unkno(ntotv,solve_sol(1)%nrhss)
  integer(ip)                :: in_up,in_lo,isist,ii
  !
  ! Change the RHS according to the renumbering strategy
  ! (original -> optimal). UNKNO is used as bridge.
  !
  do ii=1,ntotv
     unkno(ii,1)=rhsid(ii,1)
  end do

  call skyren(solve_sol(1)%nrhss,solve_sol(1)%ndofn,&
       solve_sol(1)%nequa,ntotv,unkno,rhsid,lpntn,one)
  !
  ! Compute IN_UP, IN_LO.
  !
  !call skyplu(ntotv,solve_sol(1)%lpdof,in_up,in_lo)
  !
  ! Factorize matrix, if neccesary
  !
  !call skytri(amatr(1),amatr(in_lo),amatr(in_up),solve_sol(1)%lpdof,&
  !     ntotv,solve_sol(1)%lun_solve,solve_sol(1)%kfl_solve)
  !
  ! Solve the equations 
  !
  !do isist=1,solve_sol(1)%nrhss
  !   call skybak(amatr(1),amatr(in_lo),amatr(in_up),solve_sol(1)%lpdof,&
  !        ntotv,rhsid(1,isist),solve_sol(1)%lun_solve,solve_sol(1)%kfl_solve)
  !end do
  !
  ! Change the unknown (stored in RHSID) according to the renumbering
  ! strategy (optimal -> original). 
  !     
  !call skyren(solve_sol(1)%nrhss,solve_sol(1)%ndofn,solve_sol(1)%nequa,&
  !     ntotv,unkno,rhsid,lpntn,two)

end subroutine skylin
