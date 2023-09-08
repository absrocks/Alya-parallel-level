subroutine nsa_compute_hyd_der(itask,ddens,dpres,dtemp,xcoor)
  !-----------------------------------------------------------------------
  !****f* Nastal/nsa_stalaw
  ! NAME 
  !    nsa_stalaw
  ! DESCRIPTION
  !    This routine computes the derivatives in the z-direction of 
  !    thermodynamic variables from the state law
  !
  ! USED BY
  !    nsa_gamete
  !***
  ! Simone Marras (SM)
  !-----------------------------------------------------------------------
  
  use def_master
  use def_domain
  use def_nastal
  use mod_postpr

  implicit none

  integer(ip) :: itask,ielem,ipoin,inode,pelty,pnode,pgaus
  real(rp)    :: ddens,dpres,dtemp,xcoor,xcosa,rgacp,xinve,xfact,vari1,vari2,xfac2,xauxi,bruva
  real(rp)    :: xmixr    !Water vapor mixing ratio r_v, S.M.

  !Local used for meteo only:
  real(rp)    :: rgas, cp, cv, gamma, rho_baro, a,b,c,d,e,brunt,brunt2
  real(rp)    :: pi_hyd, press_hyd, dpidz, dpdz, dthdz, xtemp

  if (itask == 1 .or. itask == -1) then        
     !
     ! Derivative of hydrostatic DENSITY with respect to z
     !
     if (lawst_nsa == 4) then 
        !
        ! It enters here for meteo:
        !
        if( kfl_benme_nsa == 1 .or. kfl_benme_nsa == 11 .or.  &
            kfl_benme_nsa == 2 .or. kfl_benme_nsa == 12) then
           !
           ! Warm bubble
           !
           pi_hyd = 1.0_rp - grnor_nsa/(cpcoe_nsa*tempe_nsa)*xcoor
           ddens = - pbaro_nsa*cvcoe_nsa*grnor_nsa/ &
                (rgasc_nsa*tempe_nsa*rgasc_nsa*cpcoe_nsa*tempe_nsa) * &
                pi_hyd**(cvcoe_nsa/rgasc_nsa - 1.0_rp)
           
        else if( kfl_benme_nsa == 3 .or. kfl_benme_nsa == 4 .or.  &
             kfl_benme_nsa == 203 .or. kfl_benme_nsa == 204) then
           !
           ! HS Linear mountain
           ! 
           pi_hyd    = exp(-grnor_nsa/(cpcoe_nsa*tempe_nsa)*xcoor)

           dpidz = -grnor_nsa*exp(-grnor_nsa*xcoor/(cpcoe_nsa*tempe_nsa))/(cpcoe_nsa*tempe_nsa)
           dpres = pbaro_nsa*cpcoe_nsa*dpidz*pi_hyd**(cpcoe_nsa/rgasc_nsa - 1.0_rp)/rgasc_nsa
           ddens = dpres/(rgasc_nsa*tempe_nsa)

        else if( kfl_benme_nsa == 5 .or. kfl_benme_nsa == 15 .or. &
                 kfl_benme_nsa == 210 .or. kfl_benme_nsa == 205) then

           !
           ! NH Linear mountain
           !
           brunt  = 0.01_rp
           brunt2 = brunt*brunt
           pi_hyd = 1.0_rp + (grnor_nsa * grnor_nsa / cpcoe_nsa / tempe_nsa / brunt2) &
                * (exp(- brunt2 * xcoor / grnor_nsa) - 1.0_rp)
           xtemp = tempe_nsa * exp(brunt2 * xcoor / grnor_nsa)
           dpidz = - grnor_nsa * exp(- brunt2 * xcoor / grnor_nsa) / cpcoe_nsa / tempe_nsa
           dtemp = tempe_nsa * brunt2 * exp(brunt2 * xcoor / grnor_nsa) / grnor_nsa
           ddens = pbaro_nsa / rgasc_nsa * (cvcoe_nsa * pi_hyd ** (cvcoe_nsa/rgasc_nsa - 1.0_rp) &
                * dpidz / rgasc_nsa / xtemp - pi_hyd ** (cvcoe_nsa/rgasc_nsa) * dtemp / xtemp / xtemp)
        end if
     end if

  else if (itask == 2) then  
     !
     ! Derivative of hydrostatic PRESSURE with respect to z
     !
     if (lawst_nsa == 4) then
        !
        ! meteo enters here:
        ! 
        if( kfl_benme_nsa == 1 .or. kfl_benme_nsa == 11 .or. &
            kfl_benme_nsa == 2 .or. kfl_benme_nsa == 12) then
           
          !! dpres = pbaro_nsa*adgam_nsa*(rgasc_nsa/pbaro_nsa)**adgam_nsa * &
          !!      (thyd(igaus)*xdhyd(igaus))**(adgam_nsa-1) * &
          !!   (xthyd(igaus)*gdhyd(ndime,igaus) + xdhyd(igaus)*gthyd(ndime,igaus))
           
        else if( kfl_benme_nsa == 3 .or. kfl_benme_nsa == 4 .or.  &
             kfl_benme_nsa == 203 .or. kfl_benme_nsa == 204) then
           !
           ! HS Linear mountain
           !
           pi_hyd    = exp(-grnor_nsa/(cpcoe_nsa*tempe_nsa)*xcoor)

           dpidz = -grnor_nsa*exp(-grnor_nsa*xcoor/(cpcoe_nsa*tempe_nsa))/(cpcoe_nsa*tempe_nsa)
           dpres = pbaro_nsa*cpcoe_nsa*dpidz*pi_hyd**(cpcoe_nsa/rgasc_nsa - 1.0_rp)/rgasc_nsa
           
        else if( kfl_benme_nsa == 5 .or. kfl_benme_nsa == 15 .or. &
                 kfl_benme_nsa == 210 .or. kfl_benme_nsa == 205) then

           brunt  = 0.01_rp
           brunt2 = brunt*brunt
           pi_hyd = 1.0_rp + (grnor_nsa * grnor_nsa / cpcoe_nsa / tempe_nsa / brunt2) &
                * (exp(- brunt2 * xcoor / grnor_nsa) - 1.0_rp)
           dpidz = - grnor_nsa * exp(- brunt2 * xcoor / grnor_nsa) / cpcoe_nsa / tempe_nsa
           dpres = pbaro_nsa * cpcoe_nsa * pi_hyd ** (cpcoe_nsa / rgasc_nsa - 1.0_rp) &
                * dpidz / rgasc_nsa
        end if
     end if

  else if (itask == 3) then   
     !
     ! Derivative of hydrostatic THETA with respect to z
     !
     if (lawst_nsa == 4) then
        !
        ! Meteo enters here
        !
        if( kfl_benme_nsa == 1 .or. kfl_benme_nsa == 11 .or. &
            kfl_benme_nsa == 2 .or. kfl_benme_nsa == 12) then
           !
           ! Warm bubble
           !
           dtemp = 0.0_rp
           
        else if( kfl_benme_nsa == 3 .or. kfl_benme_nsa == 4 .or.  &
             kfl_benme_nsa == 203 .or. kfl_benme_nsa == 204) then
           !
           ! HS Linear mountain
           !
           pi_hyd    = exp(-grnor_nsa/(cpcoe_nsa*tempe_nsa)*xcoor)
           
           dpidz = -grnor_nsa*exp(-grnor_nsa*xcoor/(cpcoe_nsa*tempe_nsa))/(cpcoe_nsa*tempe_nsa)
           dtemp = -tempe_nsa*dpidz/(pi_hyd*pi_hyd)

        else if( kfl_benme_nsa == 5 .or. kfl_benme_nsa == 15  .or. &
                 kfl_benme_nsa == 210 .or. kfl_benme_nsa == 205) then

           brunt  = 0.01_rp
           brunt2 = brunt*brunt
           dtemp = tempe_nsa * brunt2 * exp(brunt2 * xcoor / grnor_nsa) / grnor_nsa
        end if
     end if
     
  end if

end subroutine nsa_compute_hyd_der
