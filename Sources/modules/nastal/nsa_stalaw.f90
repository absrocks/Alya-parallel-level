subroutine nsa_stalaw(itask,imode,xdens,xpres,xtemp,xcoor,xmixv,mwcou,xhecp)
  !-----------------------------------------------------------------------
  !****f* Nastal/nsa_stalaw
  ! NAME 
  !    nsa_stalaw
  ! DESCRIPTION
  !    This routine computes thermodynamic variables from the state law
  ! USED BY
  !    nsa_stalaw
  !***
  !
  ! lawst_nsa = 1 --> ideal gas
  ! lawst_nsa = 2 --> isoentropic
  ! lawst_nsa = 3 --> isoentropic (wrf)
  ! lawst_nsa = 4 --> potential temperature ideal gas
  ! 
  !-----------------------------------------------------------------------
  
  use def_master
  use def_domain
  use def_nastal
  use def_kermod
  use mod_ker_proper

  implicit none

  integer(ip), intent(in):: itask,imode
  real(rp),    intent(in):: xcoor,xmixv,mwcou,xhecp                !xmixv: water vapor mixing ratio r_v, S.M.
  real(rp),    intent(inout):: xdens,xpres,xtemp

  real(rp)    :: xcosa,rgacp,xinve,xfact,vari1,vari2,xfac2,xauxi,bruva,xmowe,rgasc,adgam
  real(rp)    :: epsi     !Ratio R/Rv S.M.
  real(rp)    :: Lv       !Latent heat of vaporization: Lv = Lv0 - (Cpl - Cpv)*(T - T0)

  ! Initialize bruva:
  bruva= 0.0_rp

  !
  ! Check if the molecular weight is given as argument
  !
  if (imode /= 0 ) then
     xmowe = mwcou                   ! For computations
     rgasc = runiv_nsa / xmowe  
  else
     xmowe = mowei_nsa               ! For initilization
     rgasc = runiv_nsa / xmowe
  endif

  adgam = xhecp / (xhecp - rgasc)

  if (itask == 0) then          ! compute law constants (only by nsa_reaphy)

     afact_nsa= pbaro_nsa ** ((adgam - 1.0_rp) / adgam)
     afact_nsa= rgasc / afact_nsa

     if (lawst_nsa == 4) then
        afact_nsa= pbaro_nsa ** (adgam - 1.0_rp)
        afact_nsa= (rgasc ** adgam) / afact_nsa    ! C_o in ahmad-lindeman
     end if

  else if (itask == 1 .or. itask == -1) then        ! compute density

     if (lawst_nsa == 1) then
        xdens = xpres / rgasc / xtemp
        if (kfl_isent_nsa==1) then                                                   ! isentropic
           if (itask == -1) then
              xdens = (xpres / entro_nsa) ** (1.0_rp / adgam)
           else if (itask == 1) then
              xdens= (xtemp*rgasc/entro_nsa) ** (1.0_rp / (adgam - 1.0_rp))
           end if
        end if
     else if (lawst_nsa == 3) then
        xdens = xpres ** (1.0_rp / adgam) 
        xdens = xdens / afact_nsa / xtemp 
     else if (lawst_nsa == 4) then  !lawst_nsa == 4 means that the EQNs. are in POTENTIAL TEMPERATURE

        !with moisture:
        if (kfl_coupl_nsa == 1) then
           !Ratio R/Rv S.M.
           epsi =  rgasc/rgcou_nsa(2)  
           xdens= (xpres/afact_nsa)**(1.0_rp/adgam)
           xdens= xdens/(xtemp * ( 1.0_rp + xmixv/epsi ))
        else !dry
           xdens= (xpres/afact_nsa)**(1.0_rp/adgam)
           xdens= xdens/xtemp
        end if
        
     end if

  else if (itask == 2) then   ! compute pressure

     if (lawst_nsa == 1)  then
        xpres = rgasc * xdens * xtemp
        if (kfl_isent_nsa == 1) xpres = entro_nsa *  xdens ** adgam              ! isentropic
     else if (lawst_nsa == 3) then
        xpres = (xdens * afact_nsa * xtemp) ** adgam
     else if (lawst_nsa == 4) then
            !with moisture:
        if (kfl_coupl_nsa == 1) then
           !Ratio R/Rv S.M.
           epsi =  rgasc/rgcou_nsa(2)  
           xpres= (xdens * xtemp)**adgam 
           xpres= afact_nsa * xpres * ( 1.0_rp + xmixv/epsi )
        else !dry
           xpres= (xdens * xtemp)**adgam 
           xpres= afact_nsa * xpres
        end if
     end if

  else if (itask == 3) then   ! compute temperature

     if (lawst_nsa == 1)  then
        xauxi  = xpres / rgasc / xdens
        if (xauxi > 1.0e-8) xtemp= xauxi
        if (kfl_isent_nsa == 1)  then                                               ! isentropic
           xtemp  = xdens ** (adgam - 1.0_rp)
           xtemp  = xtemp * entro_nsa / rgasc
        end if
     else if (lawst_nsa == 3) then
        xtemp = xpres ** (1.0_rp / adgam) 
        xtemp = xtemp / xdens / afact_nsa
     else if (lawst_nsa == 4) then

        !with moisture:
        if (kfl_coupl_nsa == 1) then
           !Ratio R/Rv S.M.
           epsi =  rgasc/rgcou_nsa(2)  
           xtemp= (xpres/afact_nsa)**(1.0_rp/adgam)
           xtemp= xtemp/xdens/( 1.0_rp + xmixv/epsi )
        else !dry
           xtemp= (xpres/afact_nsa)**(1.0_rp/adgam)
           xtemp= xtemp/xdens        
        end if
     end if
     
  else if (itask == 20) then   ! compute hydrostatic pressure for initial states
  
     xfact      = - grnor_nsa / rgasc / tempe_nsa     
     vari1      = grnor_nsa * ((adgam - 1.0_rp) / adgam) / afact_nsa 
     vari2      = (adgam / (adgam - 1.0_rp))     
     xfac2      = vari1 ** vari2
     
     if (lawst_nsa == 3) then             ! isentropic
        xpres = pbaro_nsa - xfac2 * (xcoor / xtemp) ** vari2

     else if (lawst_nsa == 4) then        ! potential temperature

        rgacp = rgasc / xhecp
        xinve = 1.0_rp / rgacp
        if (kfl_infun_nsa == 3) then
           ! brunt
           bruva = stapa_nsa(1)
           xcosa = pbaro_nsa * (1.0_rp + grnor_nsa * grnor_nsa * (exp(-bruva*bruva*xcoor/grnor_nsa) - 1.0_rp) &
                / (xhecp * tempe_nsa * bruva * bruva)) ** xinve
        else
           ! isothermal
           xcosa = ((pbaro_nsa ** rgacp) / rgacp ) - (grnor_nsa * xcoor / (tempe_nsa * afact_nsa ** (1.0 / adgam)))
           xpres = (xcosa * rgacp) ** xinve 
        end if
        
     else        
        xpres = pbaro_nsa * exp(xfact * xcoor)
     end if

  else if (itask == 25) then   ! compute the equivalent potential temperature for postprocess only
  
     !
     ! Physical properties and reference values for coupled problems (CHEMIC)
     ! 1. moist air (total)
     ! 2. water vapor
     ! 3. cloud water
     ! 4. liquid water
     !
     !       mixrt_nsa,   &               ! Total mixing ration (constant for now)
     !       rgcou_nsa(5),&               ! R gas constant for species (COUPLING WITH CHEMIC)
     !       cpcou_nsa(5),&               ! Cp constant for species (COUPLING WITH CHEMIC)
     !       cvcou_nsa(5),&               ! Cv constant for species (COUPLING WITH CHEMIC)
     !       lacou_nsa(5),&               ! Reference latent heats (COUPLING WITH CHEMIC)
     !       pacou_nsa(10)                ! Coupling parameters (COUPLING WITH CHEMIC)
     
     Lv =  lacou_nsa(1) - (cpcou_nsa(4) - cpcou_nsa(2))*(xtemp - tempe_nsa)
     
     xauxi = exp(Lv * xmixv / ((rgasc + cpcou_nsa(4)*mixrt_nsa)*xtemp) )
     xfact = -rgasc / (xhecp + cpcou_nsa(4)*mixrt_nsa)
     xfact = (xpres / pbaro_nsa)**xfact
     
     xtemp = xtemp * xfact * xauxi
     

  else if (itask == 30) then   ! compute hydrostatic temperature distribution

     if (lawst_nsa == 4) then        ! potential temperature
        if (kfl_infun_nsa == 3) then
           xfact = tempe_nsa - xcoor*grnor_nsa / xhecp
           xtemp = xfact * (pbaro_nsa / xpres) ** (rgasc / xhecp)
        else
           xtemp = tempe_nsa * exp(bruva*bruva*xcoor/grnor_nsa)
        end if
     end if
     
     
  end if
  
  
end subroutine nsa_stalaw
