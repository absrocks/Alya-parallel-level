!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    def_material.f90
!> @author  Guido Giuntoli
!> @brief   Material Number definition
!> @details
!> @}
!------------------------------------------------------------------------
module def_material

  use def_kintyp

  integer(ip), parameter :: MAT_MICRO_NO_COUPLING  = 666
  integer(ip), parameter :: MAT_MICRO_ONE_WAY      = 667
  integer(ip), parameter :: MAT_MICRO_FULL         = 668

end module def_material
