subroutine skydia(gstlo,gstup,gstdi,nterm,pivot)
!-----------------------------------------------------------------------
!
!     This routine diagonalizes element in triangular decomposition.
!     Input parameters:
!         gstup(nterm) - column of upper triangular part of matrix
!         gstlo(nterm) - row of lower triangular part of matrix
!         gstdi(nterm) - reciprocal of diagonals in triangular factors
!         nterm        - number of terms in vectors
!         pivot        - diagonal in matrix to be factored
!
!     Output parameter
!         pivot        - reduced diagonal of factor
!     
!-----------------------------------------------------------------------
  use def_kintyp
  implicit none
  integer(ip) :: nterm
  real(rp)    :: gstlo(nterm),gstup(nterm),gstdi(nterm)
  real(rp)    :: pivot
  integer(ip) :: iterm
  
  do iterm= 1,nterm
     pivot=pivot-gstlo(iterm)*gstup(iterm)*gstdi(iterm)
     gstup(iterm)=gstup(iterm)*gstdi(iterm)
     gstlo(iterm)=gstlo(iterm)*gstdi(iterm)
  end do
  
end subroutine skydia
