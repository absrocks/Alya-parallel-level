!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_begrun.f90
!> @date    14/06/2019
!> @author  Guillaume Houzeaux
!> @brief   Beginning the run...
!> @details Beginning the run... we can compute matrices!
!> @}
!-----------------------------------------------------------------------

subroutine sld_begrun()

  use def_kintyp,    only : ip, rp
  use def_master,    only : INOTMASTER, fiber
  use def_domain,    only : lmate, xfiel, npoin, ndime
  use def_solidz,    only : kfl_local_sld, kfl_fiber_sld
  use def_solidz,    only : kfl_vofor_sld, vofor_sld, kfl_volca_sld
  use def_solidz,    only : lmate_sld
  use def_solidz,    only : kfl_conta_stent
  use def_solidz,    only : fibts_sld, fibtn_sld
  use mod_sld_csys,  only : sld_csys_build_jacrot
  use mod_sld_csys,  only : sld_csys_assign_material_axes
  use mod_sld_stent, only : sld_set_boundaries

  implicit none

  integer(ip)             :: ipoin
  real(rp)                :: dummr(3),vauxi(3),vmodu
 ! 
 ! Set boundaries for stent case 
 !
   if (INOTMASTER) then
      if (kfl_conta_stent /= 0_ip) then
         call sld_set_boundaries()
      end if
   end if
  !
  ! Build rotation matrix for local axes prescription
  !
  if( INOTMASTER ) then
     if ( kfl_local_sld /= 0_ip ) call sld_csys_build_jacrot()
  end if

  !
  ! Material axes at node level (fibers)
  !
  if (kfl_fiber_sld == 3 .and. INOTMASTER) then
     dummr(:) = 0.0_rp
     do ipoin = 1,npoin

        vauxi   = 0.0_rp
        vauxi(1)= 1.0_rp
        call vecpro(fiber(1,ipoin),vauxi,fibts_sld(1,ipoin),ndime)
        dummr(1:ndime) = fibts_sld(1:ndime,ipoin)
        vmodu= sqrt(dummr(1)*dummr(1) + dummr(2)*dummr(2) + dummr(3)*dummr(3))
        if (vmodu < 1.0e-10_rp) then
           vauxi   = 0.0_rp
           vauxi(2)= 1.0_rp
           call vecpro(fiber(1,ipoin),vauxi,fibts_sld(1,ipoin),ndime)
           dummr(1:ndime) = fibts_sld(1:ndime,ipoin)
           vmodu= sqrt(dummr(1)*dummr(1) + dummr(2)*dummr(2) + dummr(3)*dummr(3))
        end if

        fibts_sld(1:ndime,ipoin) = fibts_sld(1:ndime,ipoin)/vmodu

        call vecpro(fiber(1,ipoin),fibts_sld(1,ipoin),fibtn_sld(1,ipoin),ndime)
        dummr(1:ndime) = fibtn_sld(1:ndime,ipoin)
        vmodu= sqrt(dummr(1)*dummr(1) + dummr(2)*dummr(2) + dummr(3)*dummr(3))

        fibtn_sld(1:ndime,ipoin) = fibtn_sld(1:ndime,ipoin)/vmodu

     end do

  end if
  !
  ! Material axes at element level (fibers)
  !
  if( kfl_fiber_sld > 3 ) call sld_csys_assign_material_axes()
  !
  ! Material vector
  !
  lmate_sld => lmate
  !
  ! Concentrated loads using fields
  !
  if( kfl_vofor_sld > 0 ) then
     vofor_sld => xfiel(kfl_vofor_sld) % a(:,:,1)
  end if

end subroutine sld_begrun
