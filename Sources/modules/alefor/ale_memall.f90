!-----------------------------------------------------------------------
!> @addtogroup Alefor
!> @{
!> @file    ale_memall.f90
!> @author  Guillaume Houzeaux
!> @date    16/11/1966
!> @brief   Arrays allocation subroutine
!> @details Arrays allocation subroutine
!> @} 
!-----------------------------------------------------------------------
subroutine ale_memall()
  use      mod_ale_arrays
  implicit none

  call ale_arrays('ALLOCATE')

end subroutine ale_memall
      
