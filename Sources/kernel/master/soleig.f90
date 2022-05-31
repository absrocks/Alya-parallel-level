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
  
  call runend('SOLEIG: NO MORE EIGEN SOLVER')

  call cputim(time2)
  cpu_eigen = cpu_eigen + (time2-time1)

end subroutine soleig
