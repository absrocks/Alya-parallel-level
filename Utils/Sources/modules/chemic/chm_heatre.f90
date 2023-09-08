subroutine chm_heat_of_reactions(T,densi,conce,Q) 
  !-----------------------------------------------------------------------
  ! DESCRIPTION
  !    Compute the heat of all reactions
  !-----------------------------------------------------------------------
  use def_kintyp, only      : ip,rp
  use def_chemic, only      : react_chm, nreac_chm, lreac_chm, effic_chm,nspec_chm,order_chm, &
                              stofu_chm, kfl_arreh_chm, strat_chm 
  use def_kermod, only      : gasco
  use def_master, only      : speci,ittim,kfl_paral
  use mod_ker_proper 

  implicit none

  real(rp),intent(out)      :: Q(nreac_chm)
  real(rp),intent(in)       ::  T, densi
  real(rp),intent(inout)    :: conce(nspec_chm)
  integer(ip)               :: ireac,ispec,icoef
  real(rp)                  :: kforw,kback,xconF,xconB,kpres,chape, Fcent, Ntroe,Ctroe,Ftroe
  real(rp)                  :: equivalence_ratio,factor

  !
  ! Compute heat release for each reaction
  !
  do ireac = 1,nreac_chm

     !
     ! Arrhenius coefficients k = A T^beta exp(-E/RT)
     !
     kforw = react_chm(1,ireac) * (T**react_chm(2,ireac)) * exp(-react_chm(3,ireac)/(gasco*T))
     kback = react_chm(4,ireac) * (T**react_chm(5,ireac)) * exp(-react_chm(6,ireac)/(gasco*T))

     if (kfl_arreh_chm == 1) then
        if (conce(stofu_chm(2)) /= 0.0_rp ) then
           equivalence_ratio = strat_chm * conce(stofu_chm(1))/conce(stofu_chm(2))
        else
           equivalence_ratio = -0.5_rp
        endif

        !
        ! For this correction formula see
        ! Lacaze, Richardson, and Poinsot, Combustion and Flame v156 (2009) p1993
        !
        factor=(0.5_rp*(1.0+tanh( (1.25_rp - equivalence_ratio)/0.2_rp)))**2
        kforw=kforw*factor
        kback=kback*factor
     endif
     if (lreac_chm(ireac)%l(5) == 1_ip .or. lreac_chm(ireac)%l(5) == 2_ip) then  ! Troe or Lindemann forms
        
        !
        ! The following is taken from CHEMKIN documentation
        ! http://www.me.berkeley.edu/gri-mech/data/k_form.html
        kpres = kback * chape / kforw
        Fcent = &
             (1.0_rp-react_chm(7,ireac))*exp(-T/react_chm(8,ireac))+&
             react_chm(7,ireac)*exp(-T/react_chm(9,ireac))+exp(-react_chm(10,ireac)/T)
        NTroe = 0.75_rp - 1.27_rp * log(Fcent)
        CTroe = -0.4_rp - 0.67_rp * log(Fcent)
        FTroe = exp( log(Fcent) / &
             & (1.0_rp + ((log(kpres)+Ctroe)/(NTroe-0.14*(log(kpres)+CTroe)))**2 ) )
        kforw = FTroe * kforw * kpres /(1.0_rp + kpres)
        kback = 0.0_rp ! Unidirectional reaction
     endif
     
     !
     ! Now we compute concentration factor products
     !   Third body reaction with chaperone
     if (lreac_chm(ireac)%l(4) == 1_ip) then 
        chape = 0.0
        do ispec = 1,nspec_chm
           chape = chape + effic_chm(ispec,ireac)*conce(ispec) ! Chaperone concentration depends on efficiency
        enddo
     else
        chape = 1.0_rp
     endif

     xconB = 1.0_rp
     xconF = chape

     do icoef=6, 5+lreac_chm(ireac)%l(1) !Forward
        ispec = lreac_chm(ireac)%l(  icoef )
        if (order_chm(ispec,ireac,1) /= 0.0_rp) then
           if ( conce(ispec) < 0.0_rp )  conce(ispec) = 0.0_rp
           xconF = xconF * (densi * conce(ispec)/speci(ispec)%weigh)**(order_chm(ispec,ireac,1))
        endif
     enddo
     do icoef=6 +lreac_chm(ireac)%l(1),5 +lreac_chm(ireac)%l(1)+lreac_chm(ireac)%l(2) !Backward
        ispec = lreac_chm(ireac)%l(icoef )
        if (order_chm(ispec,ireac,2) /= 0.0_rp) & 
             xconB = xconB * (densi * & 
             & conce(ispec)/speci(ispec)%weigh)**(order_chm(ispec,ireac,2))
     enddo
     Q(ireac) = kforw * xconF - kback * xconB 

  enddo
end subroutine chm_heat_of_reactions
