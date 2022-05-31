!------------------------------------------------------------------------
!> @addtogroup Mathematics
!> @{
!> @file    jacobi_nre.f90
!> @author  Herbert Owen
!> @brief   computes all the eigenvalues and eigenvectors of real, symmetrix matrix A (3x3)
!> @details The eigenvalues go into matrix D, while the
!!          eigenvectors go into matrix V. NROT is the number of Jacobi rotations
!!          which were required.
!!          Beware A is altered.
!!          *This subroutine is taken from "Numerical Recipes F90", page 1225
!> @}
!------------------------------------------------------------------------

SUBROUTINE jacobi_nre(a,d,v,nrot,callersub)
  USE def_kintyp
  IMPLICIT NONE

  INTEGER(ip), INTENT(OUT) :: nrot
  !  REAL(rp), DIMENSION(:), INTENT(OUT) :: d
  !  REAL(rp), DIMENSION(:,:), INTENT(INOUT) :: a
  !  REAL(rp), DIMENSION(:,:), INTENT(OUT) :: v
  REAL(rp), DIMENSION(3), INTENT(OUT) :: d     ! for the moment I restrict it to dim 3 because asset_eq gives me trouble
  REAL(rp), DIMENSION(3,3), INTENT(INOUT) :: a
  REAL(rp), DIMENSION(3,3), INTENT(OUT) :: v
  character(*)             :: callersub  ! who is calling spcdec

  call runend('JACOBI_NRE: YOU SHOULD CODE THIS.. NOT THAT COMPLEX')

END SUBROUTINE jacobi_nre
