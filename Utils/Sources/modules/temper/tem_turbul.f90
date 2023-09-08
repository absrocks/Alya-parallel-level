subroutine tem_turbul(&
     ielem,pnode,pgaus,igaui,igauf,gpsha,gpcon,gpsph, &
     gpdif,gpgrd,gpden,gptur)   
  !-----------------------------------------------------------------------
  !****f* Temper/tem_turbul
  ! NAME 
  !    tem_turbul
  ! DESCRIPTION
  !  Coupling with TURBUL
  !    Compute effective diffusion coefficient rho*(D+D_t)
  ! USES
  ! USED BY
  !    tem_elmope_new
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      : ip,rp
  use def_domain, only      : ndime
  use def_temper, only      : prtur_tem,kfl_cotur_tem,kfl_grdif_tem,kfl_tfles_tem, &
                              kfl_regim_tem
  use def_master, only      : tfles_factor,tfles_sensor,tfles_sgseff

  implicit none
  integer(ip), intent(in)   :: ielem,pnode,pgaus,igaui,igauf
  real(rp),    intent(in)   :: gpsha(pnode,pgaus)
  real(rp),    intent(in)   :: gpcon(pgaus)
  real(rp),    intent(in)   :: gpsph(pgaus)
  real(rp),    intent(in)   :: gpden(pgaus)
  real(rp),    intent(inout):: gptur(pgaus)
  real(rp),    intent(inout):: gpgrd(ndime,pgaus)
  real(rp),    intent(out)  :: gpdif(pgaus)

  integer(ip)               :: igaus,inode,idime

  !
  ! ENTHALPY EQUATION
  ! 
  if (kfl_regim_tem == 4) then
    !
    ! Laminar contribution (rho*D)
    !
    gpdif(igaui:igauf) = gpcon(igaui:igauf) / gpsph(igaui:igauf)
    ! 
    ! Turbulent contribution (rho*D_t) for RANS & LES
    !
    if(kfl_cotur_tem /= 0) then
      !
      ! Compute mu_t for LES
      !
      if(kfl_cotur_tem < 0) gptur(igaui:igauf) = gptur(igaui:igauf) * gpden(igaui:igauf)
      !
      ! Effective diffusion coefficient: rho*(D+D_t)
      !
      gpdif(igaui:igauf) = gpdif(igaui:igauf) + gptur(igaui:igauf) / prtur_tem
 
    end if
  !
  ! TEMPERATURE EQUATION 
  !
  else
    !
    ! Laminar contribution (k)
    !
    gpdif(igaui:igauf) = gpcon(igaui:igauf)
    ! 
    ! Turbulent contribution (k_t) for RANS & LES
    !
    if(kfl_cotur_tem /= 0) then
      !
      ! Compute mu_t for LES
      !
      if(kfl_cotur_tem < 0) gptur(igaui:igauf) = gptur(igaui:igauf) * gpden(igaui:igauf)
      !
      ! Effective diffusion coefficient: (k + k_t)
      !
      gpdif(igaui:igauf) = gpdif(igaui:igauf) + gpsph(igaui:igauf)*gptur(igaui:igauf)/prtur_tem
      !
      ! TFLES (to be revised)
      !
      if (kfl_tfles_tem == 1_ip) then
         do igaus=igaui,igauf
            gpdif(igaus) = gpdif(igaus)*tfles_factor(ielem)%a(igaus)*tfles_sgseff(ielem)%a(igaus) &
                         + (1.0_rp - tfles_sensor(ielem)%a(igaus))*gpsph(igaus)*gptur(igaus)/prtur_tem
         end do
      endif

    end if

  end if

  if(kfl_grdif_tem/=0) then
     !
     ! GPGRD=grad(k+kt)
     !
     if(kfl_cotur_tem /= 0) call runend('TEM_TURBUL: GRADIENT OF TURBULENT VISCOSITY NOT CODED')

  endif

end subroutine tem_turbul
