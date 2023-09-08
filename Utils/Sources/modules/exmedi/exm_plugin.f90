!------------------------------------------------------------------------
!> @addtogroup Exmedi
!> @{
!> @file    exm_plugin.f90
!> @date    20/04/2017
!> @author  Vishal Mehta
!> @brief   Plugin for coupling
!> @details Plugin for coupling
!> @}
!------------------------------------------------------------------------
subroutine exm_plugin(icoup)

  use def_kintyp,    only :  ip,rp
  use def_domain,    only :  ndime,npoin
  use def_kermod,    only :  kfl_vefun
  use def_coupli,    only :  coupling_type
  use mod_couplings, only :  COU_INTERPOLATE_NODAL_VALUES
  use def_master,    only :  kfl_paral,mem_modul,modul
  use def_master,    only :  elmag,displ,vconc
  use mod_memory,        only :  memory_alloca
  !
  ! Possible variables => 
  !

  implicit none
  !
  ! <= end coupling variables
  !
  integer(ip), intent(in) :: icoup
  character(5)            :: variable
  real(rp)                :: dummr(ndime,1)

  real(rp),    pointer    :: xvalu(:,:)
  real(rp),    pointer    :: svalu(:,:)
  
  nullify(xvalu) 
  nullify(svalu)

  variable = coupling_type(icoup) % variable 
  !
  ! Electric impluse value
  ! 
  if( variable == 'ELMAG' .or. variable == 'UNKNO' ) then  
     !
     ! Simple check, just in case
     !
     call COU_INTERPOLATE_NODAL_VALUES(icoup,1_ip,dummr,elmag)

  else  if( variable == 'DISPL' ) then   
     !
     ! Displacement
     ! 
!!!     call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,displ_exm,displ)

  else if( variable == 'CALCIU' ) then   
     !
     ! Displacement
     ! 
     call memory_alloca(mem_modul(1:2,modul),'SVALU','sld_plugin',svalu,1_ip,npoin)     
     call memory_alloca(mem_modul(1:2,modul),'XVALU','sld_plugin',xvalu,1_ip,npoin)          
     
     xvalu(1_ip,:)= vconc(1_ip,:,1_ip)

     call COU_INTERPOLATE_NODAL_VALUES(icoup,1_ip,svalu,xvalu)


  end if

  if( associated(svalu) ) deallocate( svalu )
  if( associated(xvalu) ) deallocate( xvalu )

end subroutine exm_plugin
