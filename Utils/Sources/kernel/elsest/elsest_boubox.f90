subroutine elsest_boubox(ndime,npoin,coord,comin,comax)
  !
  ! Compute the bounding box of the domain
  !
  use def_elsest, only     : ip,rp
  implicit none
  integer(ip), intent(in)  :: ndime,npoin
  real(rp),    intent(in)  :: coord(ndime,npoin)
  real(rp),    intent(out) :: comin(ndime),comax(ndime)
  integer(ip)              :: idime,ipoin

  comax(1:ndime) = -huge(1.0_rp)
  comin(1:ndime) =  huge(1.0_rp)

  do ipoin = 1,npoin
     do idime = 1,ndime
        if( coord(idime,ipoin) > comax(idime) ) comax(idime) = coord(idime,ipoin)
        if( coord(idime,ipoin) < comin(idime) ) comin(idime) = coord(idime,ipoin)
     end do
  end do
  do idime = 1,ndime
     comax(idime) = comax(idime) + 0.000001_rp * ( comax(idime) - comin(idime) )
     comin(idime) = comin(idime) - 0.000001_rp * ( comax(idime) - comin(idime) )
  end do

end subroutine elsest_boubox
