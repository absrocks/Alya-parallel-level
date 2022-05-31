!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_tistep.f90
!> @author  Solidz Team
!> @date
!> @brief   This routine sets the time step
!> @details This routine sets the time step
!>
!> @}
!-----------------------------------------------------------------------

subroutine sld_tistep

  use def_kintyp, only : ip, rp
  use def_master, only : kfl_timco, dtinv
  use def_master, only : routp, ioutp, lun_outpu
  use mod_outfor, only : outfor
  use def_solidz, only : kfl_timei_sld, kfl_stead_sld
  use def_solidz, only : dtcri_sld, dtinv_sld
  use def_solidz, only : SLD_STATIC_PROBLEM

  implicit none

  if ( kfl_timco /= 2_ip ) then

     dtinv_sld = dtinv
     if ( kfl_timei_sld == SLD_STATIC_PROBLEM ) dtinv_sld = 0.0_rp
     if ( kfl_stead_sld == 1_ip ) dtinv_sld = 0.0_rp

  end if

  !
  ! Actualize time integration parameters
  !
  routp(1) = dtcri_sld
  routp(2) = 0.0_rp
  routp(3) = 0.0_rp

  ioutp(1) = kfl_timei_sld
  ioutp(2) = kfl_stead_sld


   call outfor(8_ip,lun_outpu,' ')

end subroutine sld_tistep
