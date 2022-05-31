!------------------------------------------------------------------------
!> @addtogroup Neutro 
!> @{
!> @file    neu_turnon.f90
!> @author  Guillaume Houzeaux
!> @brief   Turn on module
!> @details Read data and allocate memory
!> @} 
!------------------------------------------------------------------------
subroutine neu_turnon()

  use def_kintyp, only : ip
  implicit none
  !
  ! Read the physical problem
  !
  call neu_reaphy()
  !
  ! Service: Parall
  !
  call neu_parall(1_ip)
  !
  ! Initial variables
  !
  call neu_inivar(0_ip)
  !
  ! Read the numerical treatment
  !
  call neu_reanut()
  !
  ! Read the output strategy
  !
  call neu_reaous()
  !
  ! Read the boundary conditions
  !
  call neu_reabcs()
  !
  ! Service: Parall
  !
  call neu_parall(2_ip)
  !
  ! Modify boundary conditions
  !
  call neu_inibcs()
  !
  ! Initial variables
  !
  call neu_inivar(1_ip)
  !
  ! Allocate memory
  !
  call neu_memall()
  !
  ! Compute directions
  !
  call neu_directions()
  !
  ! Compute scattering
  !
  call neu_scattering()

  ! Warnings and errors
  !
  !  call neu_outerr(1_ip)

end subroutine neu_turnon

