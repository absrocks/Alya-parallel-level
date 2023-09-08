subroutine coafin(nbvar,alpha,xx,yy)
  !
  ! Solve y = y + alpha * (W.A-1.W^t) x
  !
  use def_kintyp,        only :  ip,rp
  use def_master,        only :  INOTMASTER,kfl_paral
  use def_solver,        only :  solve_sol
  use def_domain,        only :  npoin
  use mod_csrdir
  use mod_direct_solver, only : direct_solver_solution
  implicit none
  integer(ip), intent(in)    :: nbvar
  real(rp),    intent(in)    :: alpha
  real(rp),    intent(inout) :: xx(*)
  real(rp),    intent(out)   :: yy(*)
  integer(ip)                :: nrows,ngrou,ii
  real(rp),    pointer       :: mu(:),murhs(:)
  real(rp),    pointer       :: mubig(:)

!return

  ngrou = solve_sol(1) % ngrou
  nrows = npoin * nbvar

  allocate( mu(ngrou*nbvar),murhs(ngrou*nbvar) ) 
  !
  ! murhs = W^T.p
  !
  call wtvect(npoin,ngrou,nbvar,murhs,xx)
  !
  ! mu = A'-1.murhs
  !         
  if( INOTMASTER ) then
     call direct_solver_solution(solve_sol(1) % direct_solver_coarse,murhs,mu) 
     !
     ! x = W.mu
     !
     allocate( mubig(nrows) )
     call wvect(npoin,nbvar,mu,mubig) 
     !
     ! y = y + x: coarse correction
     !
     do ii = 1,nrows
        yy(ii) = yy(ii) + alpha * mubig(ii)
     end do
     deallocate( mubig )
  end if
  deallocate( mu,murhs ) 

end subroutine coafin

subroutine fincoa(nbvar,alpha,xx,yy)
  !
  ! Solve y = y + alpha * (W.A-1.W^t) x
  !
  use def_kintyp, only       :  ip,rp
  use def_master, only       :  INOTMASTER,kfl_paral
  use def_solver, only       :  solve_sol
  use def_domain, only       :  npoin
  use mod_solver, only       :  solver_smvp
  use mod_csrdir
  implicit none
  integer(ip), intent(in)    :: nbvar
  real(rp),    intent(in)    :: alpha
  real(rp),    intent(inout) :: xx(*)
  real(rp),    intent(out)   :: yy(*)
  integer(ip)                :: nrows,ngrou,ii
  real(rp),    pointer       :: mu(:),murhs(:)
  real(rp),    pointer       :: mubig(:)

!return
  call runend('FINCOA: NOT USED')

  ngrou = solve_sol(1) % ngrou
  nrows = npoin * nbvar

  allocate( mu(ngrou*nbvar),murhs(ngrou*nbvar) ) 
  !
  ! murhs = W^T.p
  !
  call wtvect(npoin,ngrou,nbvar,murhs,xx)
  !
  ! A' murhs
  !
  if( INOTMASTER ) then
     !call solver_smvp(ngrou,nbvar,nbvar,solve_sol(1) % askylpredef,&
     !     solve_sol(1) % JLpredef,solve_sol(1) % ILpredef,murhs,mu)
     !
     ! x = W.mu
     !
     allocate( mubig(nrows) )
     call wvect(npoin,nbvar,mu,mubig) 
     !
     ! y = y + x
     !
     do ii = 1,nrows
        yy(ii) = yy(ii) + alpha * mubig(ii)                                   ! Coarse correction
     end do
     deallocate( mubig )
  end if

  deallocate( mu,murhs ) 

end subroutine fincoa
