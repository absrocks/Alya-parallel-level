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

  use def_kintyp,         only :  ip,rp, lg
  use def_domain,         only :  ndime,npoin, nmate
  use def_kermod,         only :  kfl_vefun
  use def_coupli,         only :  coupling_type
  use mod_couplings,      only :  COU_INTERPOLATE_NODAL_VALUES
  use mod_couplings,      only :  COU_SET_FIXITY_ON_TARGET
  use mod_communications, only :  PAR_SUM,PAR_BARRIER
  use def_master,         only :  INOTMASTER
  use mod_parall,         only :  PAR_MY_CODE_RANK
  use def_master,         only :  kfl_paral,mem_modul,modul,INOTMASTER
  use def_master,         only :  elmag,displ,vconc
  use def_master,         only :   kfl_exm_max_nmaterials
  use def_master,         only :  kfl_eccty, kfl_cellmod
  use mod_memory,         only :  memory_alloca
  use mod_memory,         only :  memory_deallo
  use def_exmedi,         only :  kfl_fixno_exm
  use mod_sld_cou,        only :  mod_sld_cou_initexchange
  use mod_sld_cou,        only :  mod_sld_cou_physics_initialised
  use mod_exm_sld_eccoupling, only: calcium_ecc, state_ecc, troponin_ecc
  use mod_exm_sld_eccoupling, only: EXMSLD_EMEC_LAND, EXMSLD_EMEC_LAND_BIDIR

  implicit none

  integer(ip), intent(in) :: icoup
  character(5)            :: variable
  real(rp)                :: dummr(ndime,1)
  real(rp)                :: dummr2(1,1)

  real(rp),    pointer    :: xvalu(:,:)
  real(rp),    pointer    :: svalu(:,:)

  integer(ip), save       :: ipass=0
  integer(ip)             :: i, imate
  logical(lg)             :: has_land


  nullify(xvalu) 
  nullify(svalu)

  variable = coupling_type(icoup) % variable 


  if( variable == 'CALCI' ) then   
     !
     ! Calcium concentration (vconc(1,ipoin,1) )
     ! 
     if( ipass == 0 ) then
        ipass = 1
        call COU_SET_FIXITY_ON_TARGET('CALCI',modul,kfl_fixno_exm)
        call mod_sld_cou_physics_initialised()
     end if

     has_land=.False.
     do imate = 1,nmate
        if( kfl_eccty(imate) == EXMSLD_EMEC_LAND .or. kfl_eccty(imate) == EXMSLD_EMEC_LAND_BIDIR ) then
            has_land=.True.
        endif
     end do


     if (INOTMASTER) then
        call memory_alloca(mem_modul(1:2,modul),'XVALU','sld_plugin',xvalu,1_ip,npoin)          
        call memory_alloca(mem_modul(1:2,modul),'SVALU','sld_plugin',svalu,1_ip,npoin)          
        xvalu(1_ip,:)= vconc(1_ip,:,1_ip)
     else
        call memory_alloca(mem_modul(1:2,modul),'XVALU','sld_plugin',xvalu,1_ip,1_ip)          
        call memory_alloca(mem_modul(1:2,modul),'SVALU','sld_plugin',svalu,1_ip,1_ip)          
        xvalu(1_ip,:)= 0.0_rp
     end if
     
    ! SEND CALCIUM
     call COU_INTERPOLATE_NODAL_VALUES(icoup,1_ip,svalu,xvalu)

    if(has_land) then
       xvalu(1_ip,:)= 0.0_rp

      ! TROPONIN PREV IS NOT SENT BUT COMPUTED IN SOLIDZ

      ! SEND TROPONIN
       if (INOTMASTER)  xvalu(1_ip,:)= troponin_ecc(:,1)
       call COU_INTERPOLATE_NODAL_VALUES(icoup,1_ip,svalu,xvalu)

      ! SEND STATELAND
      call runend('EXM_PLUGIN: LAND COUPLING NOT POSSIBLE IN VOLUMETRIC COUPLING. READ CODE')
      ! Sending stateland from exmedi to solidz is stupid for 2 reasons:
      !     1) It's element-wise (ielem vector)
      !     2) It's done nothing with it in exmedi, just a gather.

    endif

     

  else  if( variable == 'DISPL' ) then   
     !
     ! Displacement
     ! 
!!!     call COU_INTERPOLATE_NODAL_VALUES(icoup,ndime,displ_exm,displ)
     call runend('EXM_PLUGIN: not ready yet to couple displacement')


  end if

  if( associated(xvalu) ) call memory_deallo(mem_modul(1:2,modul),'XVALU','exm_plugin',xvalu)
  if( associated(svalu) ) call memory_deallo(mem_modul(1:2,modul),'SVALU','exm_plugin',svalu)

end subroutine exm_plugin
