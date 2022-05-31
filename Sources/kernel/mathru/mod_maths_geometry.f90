!-----------------------------------------------------------------------
!
!> @defgroup Maths_Toolbox
!> Toolbox for mathematical solvers
!> @{
!> @name    ToolBox for mathematics operations
!> @file    mod_maths.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for basic geometry
!> @details ToolBox for basic geometry
!
!-----------------------------------------------------------------------

module mod_maths_geometry

  use def_kintyp_basic
  implicit none
  private

  public :: maths_point_plane_distance
  public :: maths_point_line_distance
  
contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-03
  !> @brief   Distance point to plane
  !> @details Compute the distance of a point to a plane
  !>          Handles point to line distance as well
  !> 
  !-----------------------------------------------------------------------

  pure function maths_point_plane_distance(x,a,b,c,d) result(dista)

    real(rp), intent(in) :: a
    real(rp), intent(in) :: b
    real(rp), intent(in) :: c
    real(rp), intent(in) :: d
    real(rp), intent(in) :: x(3)
    real(rp)             :: dista

    dista = (a*x(1)+b*x(2)+c*x(3)+d)/sqrt(a*a+b*b+c*c)
    
  end function maths_point_plane_distance
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-03
  !> @brief   Distance point to plane
  !> @details Compute the distance of a point to a plane
  !> 
  !-----------------------------------------------------------------------

  pure function maths_point_line_distance(x,a,b,c) result(dista)

    real(rp), intent(in) :: a
    real(rp), intent(in) :: b
    real(rp), intent(in) :: c
    real(rp), intent(in) :: x(2)
    real(rp)             :: dista

    dista = (a*x(1)+b*x(2)+c)/sqrt(a*a+b*b)

  end function maths_point_line_distance
  
end module mod_maths_geometry
!> @}
