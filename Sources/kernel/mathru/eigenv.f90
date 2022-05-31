!------------------------------------------------------------------------
!> @addtogroup Mathematics
!> @{
!> @file    eigenv.f90
!> @author  Herbert Owen
!> @brief   computes all the eigenvalues and eigenvectors of real, symmetrix matrix A (3x3)
!> @details The eigenvalues go into matrix D, while the
!!          eigenvectors go into matrix V. NROT is the number of Jacobi rotations
!!          which were required.
!!          Beware A is altered.
!!          *This subroutine is taken from "Numerical Recipes", page 346
!> @}
!------------------------------------------------------------------------

subroutine eigenv(A,D,V,NROT)

  use def_kintyp
  implicit none
  real(rp),intent(inout)   :: A(3,3)
  real(rp),intent(out)     :: D(3), V(3,3)
  integer(ip),intent(out)  :: NROT

  real(rp)                 :: B(3), Z(3), SM, THRESH
  real(rp)                 :: C, S, G, H, T, TAU, THETA, tolei, tolej, toletole, toler
  integer(ip)              :: J, idime, jdime, i

  call runend('EIGENV: YOU SHOULD CODE THIS, THIS SHOULD NOT BE COMPLICATED...!')
  
end subroutine eigenv

