!-----------------------------------------------------------------------
!
!> @defgroup Maths_Toolbox
!> Toolbox for mathematical functions and subroutines
!> @{
!> @name    ToolBox for mathematics operations
!> @file    def_maths.f90
!> @author  Guillaume Houzeaux
!> @brief   Variables
!> @details Variables fot mod_maths.f90
!
!-----------------------------------------------------------------------
module def_maths
  
  use def_kintyp,      only : ip,rp,lg,i1p
  use mod_memory,      only : memory_alloca
  use mod_memory,      only : memory_deallo
  use mod_memory,      only : memory_resize
  use mod_memory,      only : memory_size
  use mod_memory,      only : memory_copy
  use mod_std                                  ! defintion of qp and count
  implicit none
  
  integer(8)              :: memor(2)
#ifndef __PGI
  integer,     parameter  :: qp = 16
#endif
  real(rp),    parameter  :: epsil = epsilon(1.0_rp)
  real(rp),    parameter  :: pi    = 3.141592653589793238462643383279502884197_rp
  
end module def_maths
!> @}

