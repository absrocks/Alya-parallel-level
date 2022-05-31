!-----------------------------------------------------------------------
!> @addtogroup Repart
!> @{
!> @file    Repart.f90
!> @author  houzeaux
!> @date    2019-06-12
!> @brief   Repartitioning
!> @details Repartitioning of the mesh. The redistribution of data is
!>          carried out by Redist
!> @} 
!-----------------------------------------------------------------------

subroutine Repart()
  
  use mod_repartitioning, only : repartitioning
  implicit none

  call repartitioning()

end subroutine Repart
