subroutine chm_temper(&
     pgaus, enthalpy, gpmas, gpden, gpgrt, cp, gpdik, gpcon, gpgco, &
     enthalpyTransport, chemicalHeat)
  !-----------------------------------------------------------------------
  !****f* chemic/chm_temper
  ! NAME 
  !    chm_temper
  ! DESCRIPTION
  !    Sum quantities at Gauss points
  ! USES
  !-----------------------------------------------------------------------
  use def_kintyp, only      :  ip,rp
  use def_domain, only      :  ndime,mnode,mgaus
  use def_chemic, only      :  nspec_chm,entha_chm,kfl_stagg_chm

  implicit none
  integer(ip),  intent(in)  :: pgaus
  real(rp),     intent(in)  :: enthalpy(pgaus,nspec_chm)     ! Enthalpy per species
  real(rp),     intent(in)  :: gpmas(pgaus,nspec_chm)        ! Mass source term
  real(rp),     intent(in)  :: gpden(pgaus)                  ! Density
  real(rp),     intent(in)  :: gpgrt(ndime,pgaus)            ! Temperature gradient
  real(rp),     intent(in)  :: Cp(pgaus,nspec_chm)           ! Specific heat
  real(rp),     intent(in)  :: gpdik(pgaus,nspec_chm)        ! Diffusion coefficient
  real(rp),     intent(in)  :: gpcon(pgaus,nspec_chm)        ! Mass fraction
  real(rp),     intent(in)  :: gpgco(ndime,pgaus,nspec_chm)  ! Mass fraction gradient
  real(rp),     intent(out) :: enthalpyTransport(pgaus),chemicalHeat(pgaus) !Source terms
  real(rp)                  :: temp_sum

  integer(ip)               :: igaus,ispec,idime

  ! FINALLY
  ! We put the heat source into the global variable used by other modules
  !
  enthalpyTransport=0.0_rp
  do ispec = 1,nspec_chm
     do igaus = 1,pgaus
        temp_sum=0.0_rp
        do idime = 1,ndime
           temp_sum = temp_sum + gpgrt(idime,igaus) *  gpgco(idime,igaus,ispec)
        enddo
        enthalpyTransport(igaus) = enthalpyTransport(igaus) + &
             gpden(igaus) * Cp(igaus,ispec) * gpDik(igaus,ispec) * temp_sum
     enddo
  enddo

  if (kfl_stagg_chm /= 2) then  ! Staggered level 2 changes temperature on its own
     do igaus = 1,pgaus
        chemicalHeat(igaus) = 0.0_rp
        do ispec = 1,nspec_chm
           chemicalHeat(igaus) = &
                chemicalHeat(igaus) - enthalpy(igaus,ispec) * gpmas(igaus,ispec)
        enddo
     enddo
  endif
  
end subroutine chm_temper
