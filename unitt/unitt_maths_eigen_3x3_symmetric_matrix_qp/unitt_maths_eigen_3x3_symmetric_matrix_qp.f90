program unitt_maths_eigen_3x3_symmetric_matrix_qp

  use def_kintyp_basic, only : ip, rp
  use mod_maths,        only : maths_eigen_3x3_symmetric_matrix
  use mod_std

  implicit none
  real(rp)                 :: epsil = epsilon(1.0_rp)
  real(rp)                 :: A(3,3),eigVal(3),eigVec(3,3),rr,rrmax
  real(rp)                 :: rval(6)
  integer(ip)              :: iimax,inumb,num_error
  integer(ip), parameter   :: numb = 10000

  integer                  :: nums
  integer,     allocatable :: seed(:)
  integer                  :: values(8)

#ifndef __PGI

  call date_and_time(VALUES=values)
  call random_seed(size = nums)
  allocate(seed(nums))
  iimax = min(nums,size(values))
  seed(1:iimax) = values(1:iimax)
  call random_seed(put=seed)

  rrmax     = 0.0_rp
  num_error = 0

  do inumb = 1,numb

     call RANDOM_NUMBER(rval)
     rval   = (rval-0.5_rp)*2.0_rp
     A(1,1) = rval(1)
     A(1,2) = rval(2)
     A(1,3) = rval(3)
     A(2,1) = A(1,2)
     A(2,2) = rval(4)
     A(2,3) = rval(5)
     A(3,1) = A(1,3)
     A(3,2) = A(2,3)
     A(3,3) = rval(6)

     call maths_eigen_3x3_symmetric_matrix(A,eigVal,eigVec,QUAD_REAL=.true.) 
     call eigen_test(A,eigVal,eigVec,rr,rrmax,num_error)

  end do
  !
  ! Specific hard cases
  !
  do inumb = 1,numb     
     call RANDOM_NUMBER(rval)     
     A(1,1) = rval(1)
     A(1,2) = rval(2)*1.0e-10_rp
     A(1,3) = rval(3)*1.0e-10_rp
     A(2,1) = A(1,2)
     A(2,2) = rval(4)
     A(2,3) = rval(5)*1.0e-10_rp
     A(3,1) = A(1,3)
     A(3,2) = A(2,3)
     A(3,3) = rval(6)     
     call maths_eigen_3x3_symmetric_matrix(A,eigVal,eigVec,QUAD_REAL=.true.)  
     call eigen_test(A,eigVal,eigVec,rr,rrmax,num_error)     
  end do
  !
  ! max error
  !
  print*,'Maximum error=',rrmax
  if( num_error > 0 ) then
     print*,'Eigenvalues tested= ',numb*6
     print*,'Number failures=    ',num_error
     print*,'Confidency index=   ',real((numb*6-num_error),rp)/real(numb*6,rp)
     stop 1
  end if
  deallocate(seed)

#endif

contains

  subroutine eigen_test(A,eigVal,eigVec,rr,rrmax,num_error)

    implicit none
    real(rp),    intent(in)    :: A(3,3),eigVal(3),eigVec(3,3)
    real(rp),    intent(out)   :: rr
    real(rp),    intent(inout) :: rrmax
    integer(ip), intent(inout) :: num_error
    real(rp)                   :: vec1(3),vec2(3),dd
    integer(ip)                :: ii,jj

    do ii = 1,3
       !
       ! Check || A lambda - lambda v ||
       !
       vec1 = matmul(A,eigvec(:,ii))
       vec2 = eigVal(ii)*eigvec(:,ii)
       rr   = 0.0_rp
       dd   = 0.0_rp
       do jj = 1,3
          rr = rr + (vec1(jj)-vec2(jj))**2
          dd = dd + vec1(jj)**2
       end do
       rr = sqrt(rr)/sqrt(dd+epsil)
       rrmax = max(rr,rrmax)
       if( abs(rr) > 1.0e-3_rp ) then
          num_error = num_error + 1
          print*,' '
          print*,'EGEINVALUE ',ii
          print*,'---------- '
          print*,'A=         ',A(1,1:3),A(2,1:3),A(3,1:3)
          print*,'EigVal=    ',eigVal(ii)
          print*,'EigVec=    ',eigVec(:,ii)
          print*,'Residual=  ',rr
       end if
       !
       ! Check ordering
       !
       !if( eigVal(1) < eigval(2) .or. eigVal(1) < eigval(3) ) then
       !   print*,'Bad ordering'
       !   stop 1           
       !end if
       !if( eigVal(2) < eigval(3) ) then
       !   print*,'Bad ordering'
       !   stop 1           
       !end if       
    end do

  end subroutine eigen_test

end program unitt_maths_eigen_3x3_symmetric_matrix_qp
