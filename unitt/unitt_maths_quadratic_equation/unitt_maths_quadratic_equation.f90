program unitt_maths_quadratic_equation 
  
  use def_kintyp_basic, only : ip, rp
  use mod_maths,        only : maths_quadratic_equation
  
  implicit none
  real(rp)                 :: aa,bb,cc,x1,x2,rr,rrmax,dd
  real(rp),    parameter   :: epsil=epsilon(1.0_8)
  integer(ip)              :: num_solutions,ii,iimax
  integer(ip), parameter   :: numb = 100000
  real(rp)                 :: rval(3)
  integer                  :: nums
  integer,     allocatable :: seed(:)
  integer                  :: values(8)

  call date_and_time(VALUES=values)
  call random_seed(size = nums)
  allocate(seed(nums))
  iimax = min(nums,size(values))
  seed(1:iimax) = values(1:iimax)
  call random_seed(put=seed)
  rrmax = 0.0_rp
  
  do ii = 1,numb
     
     call RANDOM_NUMBER(rval)     
     aa = rval(1)
     bb = rval(2)
     cc = rval(3)
     dd = max(abs(aa),abs(bb),abs(cc),epsil)
     call maths_quadratic_equation(aa,bb,cc,x1,x2,num_solutions)

     if( num_solutions > 0 ) then
        rr    = ( aa*x1**2 + bb*x1 + cc ) / dd
        rrmax = max(rr,rrmax)
        if( abs(rr) > 1.0e-06_rp ) then
           print*,'a,b,c,x= ',aa,bb,cc,x1
           print*,'rr= ',rr
           stop 1
        end if
        rr    = ( aa*x1**2 + bb*x1 + cc ) / dd
        rrmax = max(rr,rrmax)
        if( abs(rr) > 1.0e-06_rp ) then
           print*,'a,b,c,x= ',aa,bb,cc,x2
           print*,'rr= ',rr
           stop 1
        end if
     end if
     
  end do

  deallocate(seed)

  aa = 1.7567778243243737E-007_rp
  bb = 0.96967507176162515_rp
  cc = 0.37749187734578438_rp
  dd = max(abs(aa),abs(bb),abs(cc),epsil)
  call maths_quadratic_equation(aa,bb,cc,x1,x2,num_solutions)
  
  if( num_solutions > 0 ) then
     rr    = ( aa*x1**2 + bb*x1 + cc ) / dd
     rrmax = max(rr,rrmax)
     if( abs(rr) > 1.0e-06_rp ) then
        print*,'a,b,c,x= ',aa,bb,cc,x1
        print*,'rr= ',rr
        stop 1
     end if
     rr    = ( aa*x1**2 + bb*x1 + cc ) / dd
     rrmax = max(rr,rrmax)
     if( abs(rr) > 1.0e-06_rp ) then
        print*,'a,b,c,x= ',aa,bb,cc,x2
        print*,'rr= ',rr
        stop 1
     end if
  end if
  
  stop
  
  !-5519622.2877167519 
end program unitt_maths_quadratic_equation
