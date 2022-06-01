program unitt_maths_conjugate_gradient
  !
  ! Solve Ax=b with cg donde A comes from
  ! D u = 1 in [0,1], u=0 at x=0,1 
  !
  use def_kintyp_basic, only : ip, rp
  use mod_maths_solver, only : maths_conjugate_gradient
  use mod_maths_solver, only : maths_direct
  implicit none
  real(rp)    :: A(5,5),b(5),x1(5),x2(5),resid
  integer(ip) :: ierr

  A  = 0.0_rp
  b  = 1.0_rp
  x1 = 0.0_rp
  x2 = 0.0_rp
  !
  ! Boundary conditions
  !
  b(1)   = 0.0_rp
  b(5)   = 0.0_rp
  A(1,1) = 1.0_rp
  A(5,5) = 1.0_rp

  A(2,1) = -1.0_rp
  A(2,2) =  2.0_rp
  A(2,3) = -1.0_rp
  
  A(3,2) = -1.0_rp
  A(3,3) =  2.0_rp
  A(3,4) = -1.0_rp
  
  A(4,3) = -1.0_rp
  A(4,4) =  2.0_rp
  A(4,5) = -1.0_rp

  call maths_conjugate_gradient(5_ip,A,b,x1,TOLERANCE=1.0e-10_rp,ITERATIONS=100_ip,IERR=ierr)
  call maths_direct(5_ip,A,b,x2,ierr)

  resid = sqrt( dot_product(x1-x2,x1-x2) )
  print*,'RESIDUAL= ',resid
  if( resid > 1.0e-08_rp ) stop 1
  
end program unitt_maths_conjugate_gradient
