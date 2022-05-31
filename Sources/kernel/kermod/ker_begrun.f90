!-----------------------------------------------------------------------
!> @addtogroup Nastin
!> @{
!> @file    ker_begrun.f90
!> @date    14/06/2019
!> @author  Guillaume Houzeaux
!> @brief   Beginning the run... 
!> @details Beginning the run... we can compute matrices!
!> @}
!-----------------------------------------------------------------------

subroutine ker_begrun()
  
  use mod_ker_noslwa, only : ker_noslwa
  implicit none
  !
  ! Wall exchange strategy
  ! 
  call ker_waexlo()
  !
  ! NO Slip Wall law - wall law adding extra viscosity 
  !
  call ker_noslwa() 

end subroutine ker_begrun
