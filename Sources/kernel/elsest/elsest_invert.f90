subroutine elsest_invert(a,nmax,ndm)
  !----------------------------------------------------------------------
  !
  ! This routine performs the inversion of a ndm*ndm square matrix
  ! or just part of it (nmax*nmax)
  !
  !----------------------------------------------------------------------
  use def_elsest, only      : ip,rp
  implicit none
  integer(ip), intent(in)    :: ndm,nmax
  real(rp),    intent(inout) :: a(ndm,ndm)
  integer(ip)                :: n,i,j
  real(rp)                   :: d

  do n = 1,nmax
     d = a(n,n)
     do j = 1,nmax
        a(n,j) = -a(n,j)/d
     end do
     do i = 1,nmax
        if(n.ne.i) then
           do j = 1,nmax
              if(n/=j) a(i,j) = a(i,j) +a(i,n)*a(n,j)
           end do
        end if
        a(i,n) = a(i,n)/d
     end do
     a(n,n) = 1.0_rp/d
  end do

end subroutine elsest_invert
