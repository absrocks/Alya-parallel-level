!------------------------------------------------------------------------
!> @addtogroup Partis 
!> @{
!> @file    pts_plugin.f90
!> @date    29/10/2014
!> @author  Guillaume Houzeaux
!> @brief   Plugin for coupling
!> @details Plugin for coupling
!> @}
!------------------------------------------------------------------------
subroutine pts_plugin(icoup)

  use def_kintyp,    only :  ip,rp
  use def_domain,    only :  ndime
  use def_kermod,    only :  kfl_vefun
  use def_coupli,    only :  coupling_type
  use mod_couplings, only :  COU_INTERPOLATE_NODAL_VALUES
  !
  ! Possible variables => 
  ! 
  use def_master,    only :  advec
  implicit none
  !
  ! <= end coupling variables
  !
  integer(ip), intent(in) :: icoup
  character(5)            :: variable
  real(rp)                :: dummr(ndime,1)

  variable = coupling_type(icoup) % variable 
  !
  ! Velocity/Advection
  ! 
  if( variable == 'VELOC' .or. variable == 'ADVEC' ) then  
     !
     ! Simple check, just in case
     !
     if( kfl_vefun /= 99 ) call runend('PTS_PLUGIN: VELOCITY FUNCITON SHOULD BE 99 TO USE WITH PLUGIN')
     call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,advec,dummr)
  end if
 
end subroutine pts_plugin
