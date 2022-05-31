!------------------------------------------------------------------------
!> @addtogroup Mathematics
!> @{
!> @file    spcdec.f90
!> @author  Herbert Owen
!> @brief   Computes eigenval. and eigenvec. of real, symmetrix matrix A (3x3). Sorts them in accending order.
!> @details Output is a vector D
!!          containing the eigenvalues in ascending order, and a matrix V whose
!!          columns contain the corresponding eigenvectors.
!> @}
!------------------------------------------------------------------------
subroutine spcdec(A,D,V,NROT,kfl_wivec,callersub)
  use def_kintyp
  use mod_maths, only : maths_eigen_3x3_symmetric_matrix

  implicit none
  real(rp),    intent(in)  :: A(3,3)
  real(rp),    intent(out) :: D(3), V(3,3)
  integer(ip), intent(out) :: NROT
  integer(ip), intent(in)  :: kfl_wivec   ! also obtain eigenvectors

  real(rp)                 :: E(3,3)
  real(rp)                 :: daux(3,10)!,dauxi(3),error(3)
  integer(ip)              :: idime,jdime,imeth

  character(*)             :: callersub  ! who is calling spcdec


  !V = 1.0_rp
  !D = 1.0_rp
  !return
  
  call maths_eigen_3x3_symmetric_matrix(A,D,V)  

  
  return
  

end subroutine spcdec

