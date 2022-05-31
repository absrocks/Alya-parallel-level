!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_turnon.f90
!> @author  Solidz
!> @date
!> @brief   Turn on Solidz module
!> @details
!>          \verbatim
!>          This routine performs the following tasks:
!>           - Initialize variables
!>           - Read data for the solid mechanics.
!>           - Allocate memory
!>           - Write some info
!>           - Check if there are warnings or errors
!>          \endverbatim
!>
!> @}
!------------------------------------------------------------------------

subroutine sld_turnon()

  use def_kintyp,   only : ip
  use def_master,   only : ITASK_TURNON
  use def_solidz,   only : kfl_rigid_sld
  use mod_sld_rbo,  only : sld_rbo_memall
  use mod_sld_fe2,  only : fe2_micropp_create

  implicit none
  !
  ! Initial variables
  !
  call sld_inivar(0_ip)
  !
  ! Read the physical problem
  !
  call sld_reaphy()
  !
  ! Read the numerical treatment
  !
  call sld_reanut()
  !
  ! Read the output strategy
  !
  call sld_reaous()
  !
  ! Read the boundary conditions
  !
  call sld_reabcs()
  !
  ! Service: Parall
  !
  call sld_parall(1_ip)
  !
  ! Initialize boundary conditions
  !
  call sld_inibcs()
  call sld_updbcs(ITASK_TURNON)
  !
  ! Initial variables
  !
  call sld_inivar(1_ip)
  !
  ! Allocate memory
  !
  if( kfl_rigid_sld == 0 ) then
     call sld_memall()
  else
     call sld_rbo_memall()
  end if
  !
  ! Warnings and Errors
  !
  call sld_outerr()
  !
  ! Open additional files
  !
  call sld_openfi(1_ip)
  call sld_outinf(1_ip)
  !
  ! Activate MicroPP for FE2
  !
  call fe2_micropp_create( )

end subroutine sld_turnon
