!------------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_restar.f90
!> @author  Solidz Team
!> @date    September, 2017
!>          Adds state dependent variables for any model
!> @brief   This routine reads/writes values required for a restart
!> @details Displacement is  always read and written. Velocity and
!>          and acceleration only for dynamic problems.
!>
!>          \verbatim
!>          ITASK = 1 ... Reads the initial values from the restart file
!>                  2 ... Writes restart file
!>          \endverbatim
!> @}
!------------------------------------------------------------------------

subroutine sld_restar(itask)

  use def_kintyp,             only : ip
  use def_master,             only : ITASK_READ_RESTART,ITASK_WRITE_RESTART
  use def_master,             only : coupling
  use mod_sld_arrays,         only : sld_arrays
  use mod_sld_cardiac_cycle,  only : sld_cardiac_cycle_manage_restart
  use mod_exm_sld_eccoupling, only : exm_sld_ecc_manage_restart 
  use mod_exm_sld_eccoupling, only : has_exmsld_coupling

  implicit none

  integer(ip), intent(in)    :: itask                    !< What to do

  !----------------------------------------------------------------------
  !
  ! Primary arrays
  !
  !----------------------------------------------------------------------

  if(      itask == ITASK_READ_RESTART ) then
     call sld_arrays('READ RESTART')
  else if( itask == ITASK_WRITE_RESTART ) then
     call sld_arrays('WRITE RESTART')
  end if

  !----------------------------------------------------------------------
  !
  ! Variables
  !
  !----------------------------------------------------------------------

  call sld_cardiac_cycle_manage_restart(itask)
  
  if( has_exmsld_coupling() )then
      call exm_sld_ecc_manage_restart(itask)
  endif

end subroutine sld_restar
