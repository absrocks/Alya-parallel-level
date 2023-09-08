!-----------------------------------------------------------------------
!> @addtogroup Kermod
!> @{
!> @file    nsi_endste.f90
!> @author  Guillaume Houzeaux
!> @date    05/06/2013
!> @brief   End a time step
!> @details End a time step
!> @} 
!-----------------------------------------------------------------------
subroutine ker_endste()
  use def_master,        only : ITASK_ENDSTE
  use mod_ker_detection, only : ker_detection_boundaries
  implicit none
  !
  ! Velocity function
  !
  call ker_velfun(ITASK_ENDSTE)
  !
  ! Feature detection
  !
  !call ker_detection_boundaries()

end subroutine ker_endste
