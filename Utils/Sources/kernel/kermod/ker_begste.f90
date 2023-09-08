!-----------------------------------------------------------------------
!> @addtogroup Kermod
!> @{
!> @file    nsi_begste.f90
!> @author  Guillaume Houzeaux
!> @date    05/06/2013
!> @brief   Begin a time step
!> @details Begin a time step
!> @} 
!-----------------------------------------------------------------------
subroutine ker_begste()
  use def_master, only : ITASK_BEGSTE
  implicit none
  !
  ! Assign a velocity function
  !
  call ker_velfun(ITASK_BEGSTE)
  
end subroutine ker_begste
