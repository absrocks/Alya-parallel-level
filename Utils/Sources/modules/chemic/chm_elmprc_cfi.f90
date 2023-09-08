subroutine chm_elmprc_cfi(&
  iclas,pgaus,gpcon,gpden,gpgDk,gpmas,gphco,gpsph,gptur,gpdis,gpprd,gprrt, &
  gpdif,gpgrd,gprea,gprhs)

  !-----------------------------------------------------------------------
  !****f* chemic/chm_elmprc
  ! NAME 
  !    chm_elmprc
  ! DESCRIPTION
  !    Compute terms for each species ADR equation 
  ! USES
  ! USED BY
  !    chm_elmcfi
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      :  ip,rp
  use def_domain, only      :  ndime,mgaus
  use def_chemic, only      :  nspec_chm
  use def_chemic, only      :  kfl_react_chm
  use def_chemic, only      :  kfl_cotur_chm
  use def_chemic, only      :  diffu_chm
  use def_chemic, only      :  kfl_tucfi_chm
  use def_chemic, only      :  kfl_uncfi_chm
  use def_master, only      :  kfl_paral

  implicit none

  integer(ip),  intent(in)  :: iclas
  integer(ip),  intent(in)  :: pgaus
  real(rp),     intent(in)  :: gpcon(pgaus,nspec_chm,*)
  real(rp),     intent(in)  :: gpden(pgaus)
  real(rp),     intent(in)  :: gpgDk(ndime,pgaus,nspec_chm)          ! Grad of diffusion coeff
  real(rp),     intent(in)  :: gpmas(pgaus,nspec_chm)                ! Species source term (reaction rate)
  real(rp),     intent(in)  :: gphco(pgaus)                          ! heat conductivity 
  real(rp),     intent(in)  :: gpsph(pgaus)                          ! specific heat capacity
  real(rp),     intent(in)  :: gptur(pgaus)                          ! turbulent viscosity
  real(rp),     intent(in)  :: gpdis(pgaus)                          ! Dissipation rate term in the CFI model
  real(rp),     intent(in)  :: gpprd(pgaus,nspec_chm)                ! Production term in the CFI model
  real(rp),     intent(in)  :: gprrt(pgaus)                          ! Transport of reaction rate fluctuations for LES
  
  real(rp),     intent(out) :: gpdif(pgaus)
  real(rp),     intent(out) :: gpgrd(ndime,pgaus)
  real(rp),     intent(out) :: gprea(pgaus)
  real(rp),     intent(out) :: gprhs(pgaus)

  integer(ip)               :: igaus,idime
  
  do igaus = 1,pgaus
     !
     ! Source and production terms on the RHS
     !
     gprhs(igaus) = gpmas(igaus,iclas) + gpprd(igaus,iclas)  
     !!DMM gprhs(igaus) = gpmas(igaus,iclas)
     
     ! Transport reaction rate fluctuations c * w_c or Y_c * w_c 
     !
     if (iclas == 2_ip) gprhs(igaus) = gprhs(igaus) + gprrt(igaus)
     !
     ! Diffusion coefficient
     !
     gpdif(igaus) = gphco(igaus) / gpsph(igaus)

     if (kfl_cotur_chm > 0_ip) then
       gpdif(igaus) = gpdif(igaus) + gptur(igaus) / diffu_chm(1,1)
     else if (kfl_cotur_chm < 0_ip) then
       gpdif(igaus) = gpdif(igaus) + gptur(igaus) * gpden(igaus) / diffu_chm(1,1)
     end if
     !
     ! Diffusion gradient
     !
     do idime = 1,ndime
       gpgrd(idime,igaus) = 0.0_rp 
       !!gpgrd(idime,igaus) = gpgDk(idime,igaus,iclas)
     enddo
     !
     ! Dissipation term taken implicit
     !
     if ( (iclas == 2_ip .and. (kfl_tucfi_chm == 1_ip .or. kfl_tucfi_chm == 3_ip)) .or. & 
            (iclas == 4_ip .and. (kfl_tucfi_chm == 2_ip .or. kfl_tucfi_chm == 3_ip)) ) then
     !!DMM if ( (iclas == 4_ip .and. (kfl_tucfi_chm == 2_ip .or. kfl_tucfi_chm == 3_ip)) ) then

        gprea(igaus) = gpden(igaus) * gpdis(igaus)

     endif

 enddo

end subroutine chm_elmprc_cfi
