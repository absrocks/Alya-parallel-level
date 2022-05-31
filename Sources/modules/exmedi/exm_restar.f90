!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_restar.f90
!> @date    04(03/2020
!> @author  Guillaume Houzeaux
!> @brief   Restart subroutine
!> @details Restart subroutine
!> @}
!------------------------------------------------------------------------
subroutine exm_restar(itask)
  
  use def_kintyp,     only : ip
  use def_master,     only : ITASK_READ_RESTART
  use def_master,     only : ITASK_WRITE_RESTART
  use mod_messages,   only : livinf
  use mod_exm_arrays, only : exm_arrays
  implicit none
  integer(ip), intent(in) :: itask 

  !----------------------------------------------------------------------
  !
  ! Primary arrays
  !
  !----------------------------------------------------------------------

  if( itask == ITASK_READ_RESTART ) then
     call exm_arrays('READ RESTART')
  else if( itask == ITASK_WRITE_RESTART ) then
     call exm_arrays('WRITE RESTART')
  end if

end subroutine exm_restar
 
