!-----------------------------------------------------------------------
!> @addtogroup Temper
!> @{
!> @file    tem_begrun.f90
!> @date    14/06/2019
!> @author  Guillaume Houzeaux
!> @brief   Beginning the run... 
!> @details Beginning the run... we can compute matrices!
!> @}
!-----------------------------------------------------------------------

subroutine tem_begrun()

  use def_master
  use def_temper
  implicit none
  !
  ! Check if some variables have been allocated
  !
  kfl_lookg_tem=0
  if (associated(sphec_gp)) then
     kfl_lookg_tem = 1
  end if
  
  kfl_diven_tem=0
  if (associated(div_enthalpy_transport)) then
     kfl_diven_tem = 1
  end if
  !
  ! Velocity subgrid scale
  !
  if( associated(vesgs) .and. kfl_advec_tem == 1 ) then
     kfl_sgsve_tem = 1
  else
     kfl_sgsve_tem = 0
  end if
  !
  ! Latex output format
  !
  call tem_outlat(1_ip)

end subroutine tem_begrun
 
