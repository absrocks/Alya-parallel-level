!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_output.f90
!> @author  Solidz Team
!> @date    August, 2006
!>          - Subroutine written
!> @brief   Solidz Output centralization for post-process of results
!>
!> @details
!>
!>          \verbatim
!>          Field Output:
!>          -------------
!>          Variables defined for post-processs visualization.
!>
!>          History Output:
!>          ---------------
!>          Variables defined on sets or witness points.
!>
!>          \endverbatim
!>
!> @}
!-----------------------------------------------------------------------

subroutine sld_output()

  use def_kintyp,  only : ip, rp
  use def_master,  only : ITASK_INITIA, ITASK_ENDTIM, ITASK_ENDRUN
  use def_master,  only : ittyp, nvarp
  use def_solidz,  only : kfote_sld

  implicit none

  integer(ip)          :: ivari,ivarp

  !
  ! Initialize flag fote (only after each time step)
  !
  if ( ittyp == ITASK_INITIA ) kfote_sld = 1_ip
  if ( ittyp == ITASK_ENDTIM ) kfote_sld = 0_ip
  if ( ittyp == ITASK_ENDRUN ) kfote_sld = 1_ip
  !
  ! Post-process Field Output variables
  !
  do ivarp = 1,nvarp
     ivari = ivarp
     call posdef(11_ip,ivari)
     call sld_outvar(ivari)
  end do
  !
  ! Post-process History Output variables
  !
  if ( ittyp == ITASK_INITIA .or. ittyp == ITASK_ENDTIM ) then
     !
     ! Post-process on sets
     !
     call sld_outset()
     !
     ! Post-process on witness points
     !
     call sld_outwit()
     !
     ! FEM errors (End of the run)
     !
     call sld_exaerr(2_ip)

  end if

  if ( ittyp == ITASK_ENDRUN ) then

  end if

end subroutine sld_output
