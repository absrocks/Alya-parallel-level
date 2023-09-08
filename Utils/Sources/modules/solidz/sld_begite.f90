!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_begite.f90
!> @author  Mariano Vazquez
!> @date    16/11/1966
!> @brief   This routine initializes the linearization loop
!> @details This routine initializes the linearization loop
!> @}
!-----------------------------------------------------------------------

subroutine sld_begite()

  use def_kintyp,     only : ip, rp
  use def_master,     only : itcou, itinn, modul, ITASK_BEGITE
  use def_domain,     only : kfl_elcoh
  use def_solidz,     only : kfl_goite_sld
  use def_solidz,     only : kfl_xfeme_sld
  use def_solidz,     only : volum_sld, tactv_sld
  use mod_sld_energy, only : sld_updene
  use mod_messages,   only : livinf

  implicit none

  !
  ! Initializations
  !
  kfl_goite_sld = 1
  itinn(modul)  = 0
  if ( itcou == 1 ) call sld_tistep()
  call livinf(15_ip,' ',modul)
  !
  ! Update boundary conditions
  !
  call sld_updbcs(2_ip)
  !
  ! Initialization (rhsid and amatr set to zero)
  !
  call inisol()
  !
  ! Couplings
  !
  call sld_coupli(ITASK_BEGITE)
  !
  ! Obtain the initial guess for inner iterations: unkno --> last iteration unknown
  !
  call sld_updunk(2_ip) ! Update unknowns
  call sld_updene(2_ip) ! Update energies
  !
  ! Clear counts
  !
  volum_sld = 0.0_rp
  tactv_sld = 0.0_rp
  !
  ! Update cohesive laws
  !
  if (kfl_xfeme_sld > 0 .or. kfl_elcoh > 0)  call sld_updcoh(1_ip)

end subroutine sld_begite
