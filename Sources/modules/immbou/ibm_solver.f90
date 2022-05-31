!  ******************************************************************
!  * solves the set of n linear equations a . x = b.  here a is     *
!  * input, not as the matrix a but rather as its lu decomposition, *
!  * determined by the routine ludcmp. indx is input as the permuta-*
!  * tion vector returned by ludcmp. b is input as the right-hand   *
!  * side vector b, and returns with the solution vector x. a, n and*
!  * indx are not modified by this routine and can be used for suc- *
!  * cessive calls with different right-hand sides. this routine is *
!  * also efficient for plain matrix inversion.                     *
!  ******************************************************************
subroutine ibm_solver(a,n,indx,bb)
  use def_kintyp, only :  ip,rp
  implicit none


  integer(ip), intent(in)    :: n
  integer(ip), intent(in)   :: indx(n)
  real(rp),    intent(in)    :: a(n,n)
  real(rp),    intent(inout) :: bb(n)      
  integer(ip) :: ii,i,ll,j
  real(rp)    :: sum

  ii = 0

  do i=1,n
     ll = indx(i)
     sum = bb(ll)
     bb(ll) = bb(i)
     if(ii.ne.0) then
        do j=ii,i-1
           sum = sum - a(i,j)*bb(j)
        end do ! j loop
     else if(sum.ne.0.0_rp) then
        ii = i
     end if
     bb(i) = sum
  end do ! i loop

  do i=n,1,-1
     sum = bb(i)
     if(i < n) then
        do j=i+1,n
           sum = sum - a(i,j)*bb(j)
        end do ! j loop
     end if
     bb(i) = sum / a(i,i)
  end do ! i loop

  return
end subroutine ibm_solver
