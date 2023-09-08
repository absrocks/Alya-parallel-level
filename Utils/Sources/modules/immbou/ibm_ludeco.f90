!  ***************************************************************
!  * Given an N x N matrix A, this routine replaces it by the LU *
!  * decomposition of a rowwise permutation of itself. A and N   *
!  * are input. INDX is an output vector which records the row   *
!  * permutation effected by the partial pivoting; D is output   *
!  * as -1 or 1, depending on whether the number of row inter-   *
!  * changes was even or odd, respectively. This routine is used *
!  * in combination with LUBKSB to solve linear equations or to  *
!  * invert a matrix. Return code is 1, if matrix is singular.   *
!  ***************************************************************
subroutine ibm_ludeco(a,n,indx,d,code)
  use def_kintyp, only :  ip,rp
  use def_master, only :  zeror
  implicit none

  integer(ip), intent(in)    :: n
  integer(ip), intent(out)   :: indx(n),d,code
  real(rp),    intent(inout) :: a(n,n)
  integer(ip)                :: i,j,k,imax
  real(rp)                   :: amax,dum,sum,tiny,vv(100)

  tiny = zeror
  !tiny = 1.5e-10_rp
  d=1; code=0

  do i=1,n
     amax=0.0_rp
     do j=1,n
        if (abs(a(i,j)) > amax) amax=abs(a(i,j))
     end do ! j loop
     if(amax < tiny) then
        call runend('KRIGING: THE COVARIANCE MATRIX IS SINGULAR')
     end if
     vv(i) = 1.0_rp / amax
  end do ! i loop

  do j=1,n
     do i=1,j-1
        sum = a(i,j)
        do k=1,i-1
           sum = sum - a(i,k)*a(k,j) 
        end do ! k loop
        a(i,j) = sum
     end do ! i loop
     amax = 0.0_rp
     do i=j,n
        sum = a(i,j)
        do k=1,j-1
           sum = sum - a(i,k)*a(k,j) 
        end do ! k loop
        a(i,j) = sum
        dum = vv(i)*abs(sum)
        if(dum >= amax) then
           imax = i
           amax = dum
        end if
     end do ! i loop  

     if(j /= imax) then
        do k=1,n
           dum = a(imax,k)
           a(imax,k) = a(j,k)
           a(j,k) = dum
        end do ! k loop
        d = -d
        vv(imax) = vv(j)
     end if

     indx(j) = imax
     if(abs(a(j,j)) < tiny) a(j,j) = tiny

     if(j /= n) then
        dum = 1.0_rp / a(j,j)
        do i=j+1,n
           a(i,j) = a(i,j)*dum
        end do ! i loop
     end if
  end do ! j loop

  return
end subroutine ibm_ludeco
