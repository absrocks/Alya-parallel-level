subroutine chm_nastal(&
     pgaus, enthalpy, grad_enthalpy, gpmas, gpdik, gp_grad_Dik, &
     gpcon, gpgco,gp_lap_con, enthalpy_transport, div_enthalpy_transport,chemical_heat)
 
  !-----------------------------------------------------------------------
  !****f* partis/chm_heasum
  ! NAME 
  !    chm_elmpre
  ! DESCRIPTION
  !    Sum quantities at Gauss points
  ! USES
  !-----------------------------------------------------------------------
  use def_kintyp, only      :  ip,rp
  use def_master, only      :  speci
  use def_domain, only      :  ndime,mnode,mgaus
  use def_chemic, only      :  nspec_chm,entha_chm
  implicit none
  integer(ip),  intent(in)  :: pgaus
  real(rp),     intent(in) :: enthalpy(pgaus,nspec_chm),grad_enthalpy(ndime,pgaus,nspec_chm) ! Enthalpy per species
  real(rp),     intent(in) :: gpmas(pgaus,nspec_chm)    !Mass source term
  real(rp),     intent(in) :: gpdik(pgaus,nspec_chm),gp_grad_Dik(ndime,pgaus,nspec_chm)          ! Diffusion coefficient
  real(rp),     intent(in) :: gpcon(pgaus,nspec_chm)  ! Conce gradient
  real(rp),     intent(in) :: gpgco(ndime,pgaus,nspec_chm)  ! Conce gradient
  real(rp),     intent(in) :: gp_lap_con(pgaus,nspec_chm)  ! Conce laplacian
  real(rp),     intent(out):: enthalpy_transport(ndime,pgaus)          !
  real(rp),     intent(out):: div_enthalpy_transport(pgaus)          !
  real(rp),     intent(out):: chemical_heat(pgaus)          !
  real(rp)                 :: rtemp_h

  integer(ip)               :: igaus,ispec,idime

  ! FINALLY
  ! We put the heat source into the global variable used by other modules
  enthalpy_transport=0.0_rp
  do ispec = 1,nspec_chm
     do igaus = 1,pgaus
        rtemp_h =  enthalpy(igaus,ispec) *  gpcon(igaus, ispec) *  gpDik(igaus,ispec)
        div_enthalpy_transport(igaus) = &
             enthalpy(igaus, ispec) * gpcon(igaus, ispec) * gpDik(igaus,ispec)* gp_lap_con(igaus,ispec)   
        do idime = 1,ndime
           enthalpy_transport(idime,igaus) = &
                enthalpy_transport(idime,igaus) + rtemp_h * gpgco(idime,igaus,ispec)
           div_enthalpy_transport(igaus) = &
                div_enthalpy_transport(igaus) + &
                grad_enthalpy(idime, igaus, ispec) * gpcon(igaus, ispec) *  gpDik(igaus,ispec)* gpgco(idime,igaus,ispec) +&
                enthalpy(igaus, ispec) * gpgco(idime,igaus,ispec) *  gpDik(igaus,ispec)* gpgco(idime,igaus,ispec) +&
                enthalpy(igaus, ispec) * gpcon(igaus, ispec) *  gp_grad_Dik(idime, igaus,ispec)* gpgco(idime,igaus,ispec) +&
                enthalpy(igaus, ispec) * gpcon(igaus, ispec) *  gpDik(igaus,ispec)* gp_lap_con(igaus,ispec)                
        enddo
     enddo
  enddo

  chemical_heat=0.0_rp
  do ispec = 1,nspec_chm
     do igaus = 1,pgaus
        chemical_heat(igaus) = &
             chemical_heat(igaus) - 1000.0_rp* speci(ispec)%entha(1) * gpmas(igaus,ispec)
     enddo
  enddo
  
end subroutine chm_nastal
