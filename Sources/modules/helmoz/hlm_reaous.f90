subroutine hlm_reaous()
  !------------------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_reaous.f90
  ! NAME 
  !    hlm_reaous
  ! DESCRIPTION
  !    This routine reads the postprocess.
  ! USES
  ! USED BY
  !    hlm_turnon
  !------------------------------------------------------------------------

  use def_parame
  use def_inpout
  use def_master
  use def_helmoz
  use def_domain
  use mod_ecoute, only :  ecoute

  implicit none
  integer(ip) :: ii,dummi

  call reaous()

end subroutine hlm_reaous
 
