real(rp) function rad_gamma(Absorption, Scattering, Anisotropy)  result(Gamma)
  !
  ! Compute the Gamma factor for the equation of P1 model
  !
  use def_kintyp, only     :  rp
  implicit none
  real(rp), intent(IN) :: absorption, scattering, Anisotropy

  Gamma = 1.0_rp / ( 3.0_rp * (Absorption+Scattering) - Anisotropy )

end function rad_gamma
