 !------------------------------------------------------------------------
!> @addtogroup Neutro 
!> @{
!> @file    neu_begste.f90
!> @date    31/03/2016
!> @author  Guillaume Houzeaux
!> @brief   Begin time step
!> @details begin time step
!> @} 
!------------------------------------------------------------------------

subroutine neu_begste()

  use def_kintyp
  use def_neutro
  implicit none
  !
  ! Initial guess fo the velocity: u(n,0,*) <-- u(n-1,*,*).
  !
  call neu_updunk(1_ip)

end subroutine neu_begste
