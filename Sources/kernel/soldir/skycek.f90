subroutine skycek(gstup,nterm,saval)
!-----------------------------------------------------------------------
!
! This routine tests for rank
!
!-----------------------------------------------------------------------
  use def_kintyp
  implicit none
  integer(ip) :: nterm
  real(rp)    :: gstup(nterm)
  real(rp)    :: saval
  integer(ip) :: iterm
  
  saval = 0.0_rp
  do iterm = 1,nterm
     saval = saval + abs(gstup(iterm))
  end do
  
end subroutine skycek
