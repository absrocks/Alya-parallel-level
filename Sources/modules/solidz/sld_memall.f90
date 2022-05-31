!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_memall.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   General memory allocation
!> @details General memory allocation
!> @}
!-----------------------------------------------------------------------

subroutine sld_memall()

  use mod_sld_arrays, only : sld_arrays

  implicit none

  !----------------------------------------------------------------------
  !
  ! Primary variables
  !
  !----------------------------------------------------------------------
  call sld_arrays('ALLOCATE')

end subroutine sld_memall

