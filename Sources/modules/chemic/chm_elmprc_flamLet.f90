subroutine chm_elmprc_flamLet(&
  iclas,pgaus,gpden,gpmas,gphco,gpsph,gptur,gpdis,gpprd,gprrt, &
  gpdif,gpgrd,gprhs)

  !-----------------------------------------------------------------------
  !****f* chemic/chm_elmprc
  ! NAME 
  !    chm_elmprc
  ! DESCRIPTION
  !    Compute terms for each species ADR equation 
  ! USES
  ! USED BY
  !    chm_element_operations
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      :  ip,rp
  use def_domain, only      :  ndime,mgaus
  use def_chemic, only      :  nclas_chm
  use def_chemic, only      :  kfl_cotur_chm
  use def_chemic, only      :  diffu_chm

  implicit none

  integer(ip),  intent(in)  :: iclas
  integer(ip),  intent(in)  :: pgaus
  real(rp),     intent(in)  :: gpden(pgaus)
  real(rp),     intent(in)  :: gpmas(pgaus,nclas_chm)                ! Species source term (reaction rate)
  real(rp),     intent(in)  :: gphco(pgaus)                          ! heat conductivity 
  real(rp),     intent(in)  :: gpsph(pgaus)                          ! specific heat capacity
  real(rp),     intent(in)  :: gptur(pgaus)                          ! turbulent viscosity
  real(rp),     intent(in)  :: gpdis(pgaus,nclas_chm)                ! Dissipation rate term in the Flamelet model
  real(rp),     intent(in)  :: gpprd(pgaus,nclas_chm)                ! Production term in the Flamelet model
  real(rp),     intent(in)  :: gprrt(pgaus,nclas_chm)                ! Transport of reaction rate fluctuations for LES
  
  real(rp),     intent(out) :: gpdif(pgaus,nclas_chm)
  real(rp),     intent(out) :: gpgrd(ndime,pgaus)
  real(rp),     intent(out) :: gprhs(pgaus)

  integer(ip)               :: igaus,idime
 
  do igaus = 1,pgaus
     !
     ! RHS terms: Source + Production term + Transport reaction rate fluctuations + Dissipation term
     !
     gprhs(igaus) = gprhs(igaus) + gpmas(igaus,iclas) + gpprd(igaus,iclas) + gprrt(igaus,iclas) - gpdis(igaus,iclas)

     !
     ! Diffusion coefficient
     !
     gpdif(igaus,iclas) = gphco(igaus) / gpsph(igaus)

     !
     ! Adding turbulent part
     !
     if (kfl_cotur_chm > 0_ip) then
        !
        ! RANS
        !
        gpdif(igaus,iclas) = gpdif(igaus,iclas) + gptur(igaus) / diffu_chm(1,1)
     else if (kfl_cotur_chm < 0_ip) then
        !
        ! LES:
        !
        gpdif(igaus,iclas) = gpdif(igaus,iclas) + gptur(igaus) * gpden(igaus) / diffu_chm(1,1)
     end if

 enddo

end subroutine chm_elmprc_flamLet
