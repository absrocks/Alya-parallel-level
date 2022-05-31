!------------------------------------------------------------------------
!> @addtogroup Insitu 
!> @{
!> @file    ins_plugin.f90
!> @date    20/04/2017
!> @author  Vishal Mehta
!> @brief   Plugin for coupling
!> @details Plugin for coupling
!> @}
!------------------------------------------------------------------------
subroutine ins_plugin(icoup)
  use mod_communications, only : PAR_GATHERV
  use def_kintyp,    only :  ip,rp
  use def_domain,    only :  ndime,npoin
  use def_kermod,    only :  kfl_vefun
  use def_coupli,    only :  coupling_type
  use mod_couplings, only :  COU_INTERPOLATE_NODAL_VALUES
  use def_master,    only :  kfl_paral,INOTMASTER,npoi1,npoi2,npoi3
  !
  ! Possible variables => 
  !
  use def_insitu
  implicit none
  real(rp),pointer          :: dummmr(:)
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
  !  if( variable == 'VELOC' .or. variable == 'UNKNO' ) then
    if( variable == 'VELOC' ) then  
     !
     ! Simple check, just in case
     !
     call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,advec,dummr)
     
  end if

  if( variable == 'ELMAG' ) then  
     !
     ! Simple check, just in case
     !
     call COU_INTERPOLATE_NODAL_VALUES(icoup,1,elmag,dummr)

     if(INOTMASTER) then
        elmag_buf(1:npoi1) = elmag(1:npoi1)
        elmag_buf(npoi2:npoi3) = elmag(npoi2:npoi3)
        npoin_all = 0
        call PAR_GATHERV(elmag_buf,dummmr,npoin_all,'IN MY CODE')
     else
        call PAR_GATHERV(dummmr,elmag,npoin_all,'IN MY CODE')
     end if

#ifdef INVIZ
     call updateframe() 
#endif
     
  end if
end subroutine ins_plugin
