!-----------------------------------------------------------------------
!> @addtogroup Timste
!> @{
!> @file    Timste.f90
!> @author  Guillaume Houzeaux
!> @brief   Compute time step
!> @details Compute time step:
!>          - Compute the time step for each module
!>          - Compute the global time step 
!>
!> @}
!-----------------------------------------------------------------------
subroutine Timste()

  use def_master
  use def_parame
  use mod_ker_subdomain, only : ker_subdomain_update_coordinates
  use mod_messages, only : livinf
  implicit none
  !
  ! Live information
  !
  call livinf(4_ip,' ',zero)
  !
  ! Initializations
  !
  call iniste(1_ip)
  !
  ! Compute the time step for each module
  !
  call moduls(ITASK_TIMSTE)
  !
  ! Initializations
  !
  call iniste(2_ip)
  !
  ! Computes the global time step
  !
  call setgts(2_ip)
  !
  ! Live information
  !
  call livinf(18_ip,' ',zero)
  !
  ! Update subdomains
  !
  call ker_subdomain_update_coordinates()

end subroutine Timste
