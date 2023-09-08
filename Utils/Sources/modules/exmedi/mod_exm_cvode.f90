!-----------------------------------------------------------------------
!> @addtogroup ExmediCVode
!> @ingroup    Exmedi
!> @ingroup    ExmediElmoperations
!> @{
!> @file    mod_exm_cvode.f90
!> @author  Mariano Vazquez
!> @date    01/12/2016
!> @brief   Module bridge to Sundials CVode 
!> @details Module bridge to Sundials CVode 
!
!-----------------------------------------------------------------------
module mod_exm_cvode
  use      def_master
  use      def_domain
  use      def_elmtyp
  use      def_exmedi
  implicit none

  private

  integer(ip), parameter :: VECTOR = 16

  !
  ! Public stuff
  !
  public :: exm_cvode_initia


  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  
contains
  !-----------------------------------------------------------------------
  !> @author  Mariano Vazquez
  !> @date    22/11/2016
  !> @brief   Compute elemental matrix and RHS  
  !> @details Compute elemental matrix and RHS
  !-----------------------------------------------------------------------
  subroutine exm_cvode_initia
    
  
  end subroutine exm_cvode_initia


end module mod_exm_cvode
!> @}
