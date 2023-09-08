subroutine tur_turnof()
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_turnof
  ! NAME 
  !    tur_turnof
  ! DESCRIPTION
  !    This routine closes the run for the temperature equation
  ! USES
  !    tur_outcpu
  !    tur_output
  ! USED BY
  !    Turbul
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_solver
  use def_turbul
  implicit none

  !-------------------------------------------------------------------
  !
  ! Interpolation from coarse to fine mesh
  !
  !-------------------------------------------------------------------
  !
  if(kfl_meshi_tur /= 0_ip) call tur_coarfine(2_ip)

end subroutine tur_turnof
