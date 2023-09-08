subroutine soleig(amatr,eigva,eigen,bmatr,iter)
  !-----------------------------------------------------------------------
  !****f* master/soleig
  ! NAME 
  !    solver
  ! DESCRIPTION
  !    This routine calls the solvers
  ! USES
  !    memchk
  !    mediso
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  rp,ip
  use def_solver, only     :  eigen_sol,cpu_eigen
  use def_domain, only     :  npoin,r_dom,c_dom
  implicit none
  real(rp),    intent(in)  :: eigva(*),eigen(*)
  real(rp),    intent(in)  :: amatr(*),bmatr(*) 
  integer(ip), intent(in)  :: iter
  integer(ip)              :: izdom,ipoin,jpoin
  real(rp)                 :: time1,time2

  call cputim(time1)
  !
  ! Output  of A and B in MatrixMarket format
  !
!!$  write(90,'(a)') '%%MatrixMarket matrix coordinate real general'
!!$  write(91,'(a)') '%%MatrixMarket matrix array real general'
!!$  write(90,'(3(1x,i12))') npoin,npoin,eigen_sol(1) % nzmat
!!$  write(91,'(3(1x,i12))') npoin,1_ip
!!$  do ipoin = 1,npoin
!!$     write(91,'(1(1x,i9),1x,e12.6)') ipoin,bmatr(ipoin)
!!$     do izdom = r_dom(ipoin),r_dom(ipoin+1)-1
!!$        jpoin = c_dom(izdom)
!!$        write(90,'(2(1x,i9),1x,e12.6)') ipoin,jpoin,amatr(izdom)
!!$     end do
!!$  end do
!!$  write(92,'(3(1x,i12))') npoin
!!$  write(93,'(3(1x,i12))') npoin,eigen_sol(1) % nzmat
!!$  do ipoin = 1,npoin+1
!!$     write(92,*) r_dom(ipoin)
!!$  end do
!!$  do izdom = 1,eigen_sol(1) % nzmat
!!$     write(93,*) c_dom(izdom)
!!$  end do
!!$  flush(90_ip)
!!$  flush(91_ip)
!!$  flush(92_ip)
!!$  flush(93_ip)
!!$  call runend('MATRIX POSTPROCESSED')
  !
  ! Eigen solver
  !
  if( eigen_sol(1)%kfl_algso == 0 ) then

     call eigdir(eigen,eigva,amatr,bmatr,npoin,eigen_sol(1)%neiva)

  else if( eigen_sol(1)%kfl_algso == 1 ) then

     call eigit1(amatr,eigva,eigen,bmatr,iter)

  else if( eigen_sol(1)%kfl_algso == 2 ) then

     call eigit2(amatr,eigva,eigen,bmatr,iter)

  else if( eigen_sol(1)%kfl_algso == 3 ) then

     call eigit3(amatr,eigva,eigen,bmatr,iter)

  end if

  call cputim(time2)
  cpu_eigen = cpu_eigen + (time2-time1)

end subroutine soleig
