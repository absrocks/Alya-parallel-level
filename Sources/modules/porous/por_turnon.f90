!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_turnon.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Turn on Porous module
!> @details Read data, allocate memory and check errors
!> @} 
!------------------------------------------------------------------------
subroutine por_turnon()
  use def_parame
  use def_domain
  use def_master
  use def_porous
  implicit none
  !
  ! Initial variables
  !
  call por_inivar(zero)
  !
  ! Read the physical problem
  !
  call por_reaphy()
  !
  ! Read the numerical treatment
  !
  call por_reanut()
  !
  ! Read the output strategy
  !
  call por_reaous()
  !
  ! Read the boundary conditions
  !
  call por_reabcs()   
  !
  ! Parall service
  !
  call por_parall(1_ip)
  !
  ! Initial variables
  !
  call por_inibcs()
  !
  ! Initial variables  hhhpor For the moment I guess I dont need it
  !
  call por_inivar(1_ip)
  !
  ! Write info
  !
  call por_outinf()
  !
  ! Warnings and errors
  !
  call por_outerr()
  !
  ! Allocate memory
  !
  call por_memall()
  !
  ! Open additional files
  !
  call por_openfi(two)

end subroutine por_turnon
