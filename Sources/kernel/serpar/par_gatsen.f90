subroutine par_gatsen()
  !-------------------------------------------------------------------------------
  !****f* parall/par_gatsen
  ! NAME
  !    par_gatsen
  ! DESCRIPTION
  !    All gather
  ! INPUT
  ! OUTPUT
  ! USED BY
  !***
  !-------------------------------------------------------------------------------
  use def_master
  use def_domain
  use def_solver
  use def_parall
  use mod_parall, only : PAR_COMM_MY_CODE4
  implicit none
#ifdef MPI_OFF
#else
  include  'mpif.h' 
#endif
  integer(4)  :: istat,nsmal4
  integer(ip) :: ii,kk,nsmall,nbig

  nsmall = solve_sol(1) % lcou4(kfl_paral+1)
  nbig   = solve_sol(1) % nbig
  nsmal4 = int(nsmall,4)

  if( ISLAVE ) then 
     !
     ! Slave scatter
     !
     kk = solve_sol(1) % disp4(kfl_paral+1)
     do ii = 1,nsmall
        kk                        = kk + 1
        solve_sol(1) % xsmall(ii) = parre( solve_sol(1) % lbig(kk) )
     end do
  end if
  !
  ! MPI all gather
  !
#ifdef MPI_OFF
#else
  CALL MPI_ALLGATHERV(&
       solve_sol(1) % xsmall(1:nsmall), nsmal4, MPI_DOUBLE_PRECISION,&
       solve_sol(1) % xbig(1:nbig), solve_sol(1) % lcou4(1:npart_par+1) ,&
       solve_sol(1) % disp4(1:npart_par+1), MPI_DOUBLE_PRECISION,&
       PAR_COMM_MY_CODE4, istat)
#endif
  !
  ! Gather to local array
  !
  if( ISLAVE ) then
     do ii = 1,size(parre)
        parre(ii) = 0.0_rp    
     end do
     do ii = 1,solve_sol(1) % nbig
        kk = solve_sol(1) % lbig(ii)
        parre(kk) = parre(kk) + solve_sol(1) % xbig(ii)
     end do
  end if

end subroutine par_gatsen
