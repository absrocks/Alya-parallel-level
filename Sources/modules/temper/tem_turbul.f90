subroutine tem_turbul(&
     pnode,pgaus,igaui,igauf,gpsha,gpcon,gpsph, &
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
  use def_kermod, only      : turmu_ker
  use def_temper, only      : prtur_tem,kfl_grdif_tem, &
                              kfl_regim_tem, kfl_diven_tem, kfl_rhs_scal_tem
  use def_master, only      : div_enthalpy_transport,kfl_coupl,ID_TEMPER,ID_CHEMIC, &
                              kfl_htran

  implicit none
  integer(ip), intent(in)   :: pnode,pgaus,igaui,igauf
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
    if (kfl_coupl(ID_TEMPER,ID_CHEMIC) >= 1 .and. associated(div_enthalpy_transport) .and. kfl_htran == 1) then
       gpdif(igaui:igauf) = 0.0_rp   ! In this case div_enthalpy_transport is computed in Chemic and no additional diffusion term is required for the energy eq.
    else
       gpdif(igaui:igauf) = gpcon(igaui:igauf) / gpsph(igaui:igauf)
    end if
    ! 
    ! Turbulent contribution (rho*D_t) for RANS & LES
    !
    if( turmu_ker % kfl_exist /= 0_ip ) then
      !
      ! Compute mu_t Compute mu_t -- now all are multiplied by gpden outside
      !
      gptur(igaui:igauf) = gptur(igaui:igauf) * gpden(igaui:igauf)
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
    if( turmu_ker % kfl_exist /= 0_ip ) then
      !
      ! Compute mu_t -- now all are multiplied by gpden outside
      !
      gptur(igaui:igauf) = gptur(igaui:igauf) * gpden(igaui:igauf)
      !
      ! Effective diffusion coefficient: (k + k_t)
      !
      gpdif(igaui:igauf) = gpdif(igaui:igauf) + gpsph(igaui:igauf)*gptur(igaui:igauf)/prtur_tem

    end if


  end if

  if ( kfl_rhs_scal_tem > 0 ) gpdif(igaui:igauf) = gpdif(igaui:igauf) / gpden(igaui:igauf)  

  if(kfl_grdif_tem/=0) then
     !
     ! GPGRD=grad(k+kt)
     !
     if( turmu_ker % kfl_exist /= 0_ip ) call runend('TEM_TURBUL: GRADIENT OF TURBULENT VISCOSITY NOT CODED')

  endif

end subroutine tem_turbul
