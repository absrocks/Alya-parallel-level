subroutine tur_nodpro(ipoin,rho,mu)
  !-----------------------------------------------------------------------
  !****f* Turbul/tur_nodpro
  ! NAME 
  !    tur_nodpro
  ! DESCRIPTION
  !    This routine updates the nodal density and viscosity
  ! OUTPUT
  !    RHO ... Density
  !    MU .... Viscosity
  ! USED BY
  !    tur_updunk
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_master, only     :  densi,visco
  use def_turbul, only     :  kfl_colev_tur,lawde_tur,lawvi_tur,&
       &                      detur_tur,vitur_tur,densi_tur,visco_tur
  implicit none
  integer(ip), intent(in)  :: ipoin
  real(rp),    intent(out) :: rho,mu

  if( kfl_colev_tur /= 0 ) then
     !
     ! Projected properties
     !
     rho = detur_tur(ipoin)
     mu  = vitur_tur(ipoin)

  else
     !
     ! Properties with laws
     !
     if( lawde_tur == 0 ) then
        rho = densi_tur(1)
     else if( lawde_tur == 1 ) then
        rho = densi(ipoin,1)
     end if
     if( lawvi_tur == 0 ) then
        mu = visco_tur(1)
     else if( lawvi_tur == 1 ) then
        mu = visco(ipoin,1)
     end if

  end if

end subroutine tur_nodpro
