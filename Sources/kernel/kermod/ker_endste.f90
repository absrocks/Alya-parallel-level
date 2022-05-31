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
  use def_kintyp,        only : ip 
  use def_master,        only : ITASK_ENDSTE
  use mod_ker_detection, only : ker_detection_boundaries
  implicit none
  !
  ! Feature detection
  !
  !call ker_detection_boundaries()
  !
  ! Velocity, temperature, concentration and displacement functions
  ! (:,3) <= (:,1)
  !
  call ker_velfun(ITASK_ENDSTE)
  call ker_temfun(ITASK_ENDSTE)
  call ker_confun(ITASK_ENDSTE)
  call ker_disfun(ITASK_ENDSTE)

end subroutine ker_endste
