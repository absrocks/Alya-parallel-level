subroutine nsa_compute_hyd(itask,xdens,xpres,xtemp,xcoor)
  !-----------------------------------------------------------------------
  !****f* Nastal/nsa_stalaw
  ! NAME 
  !    nsa_stalaw
  ! DESCRIPTION
  !    This routine computes the hydrostatic background fields
  ! USED BY
  !    nsa_stalaw
  !***
  !
  !-----------------------------------------------------------------------
  
  use def_master
  use def_domain
  use def_nastal
  use mod_postpr

  implicit none

  integer(ip) :: itask,ielem,ipoin,inode,pelty,pnode,pgaus
  real(rp)    :: xdens,xpres,xtemp,xcoor,xcosa,rgacp,xinve,xfact,vari1,vari2,xfac2,xauxi,bruva
  real(rp)    :: xmixr    !Water vapor mixing ratio r_v, S.M.
  
  !Local used for meteo only:
  real(rp)    :: rgas, cp, cv, gamma, rho_baro, a,b,c,d,e,brunt,brunt2
  real(rp)    :: pi_hyd, press_hyd, densi_hyd !some auxiliary variables

  if (itask == 1 .or. itask == -1) then        
     !
     ! compute density
     !
     if (lawst_nsa == 4) then 
        !
        ! It enters here for meteo:
        !
        if( kfl_benme_nsa == 1 .or. kfl_benme_nsa == 11 .or. &
            kfl_benme_nsa == 2 .or. kfl_benme_nsa == 12) then
           !
           !Warm bubble, density current
           !
           pi_hyd = 1.0_rp - grnor_nsa/(cpcoe_nsa*tempe_nsa)*xcoor
           xdens  = pbaro_nsa/(rgasc_nsa*xtemp)*(pi_hyd)**(cvcoe_nsa/rgasc_nsa)
           !OK
        
        else if( kfl_benme_nsa == 3 .or. kfl_benme_nsa == 4 .or. &
             kfl_benme_nsa == 203 .or. kfl_benme_nsa == 204) then
           !
           ! HS Linear mountain
           !
           pi_hyd    = exp(-grnor_nsa/(cpcoe_nsa*tempe_nsa)*xcoor);
           press_hyd = pbaro_nsa*pi_hyd**(cpcoe_nsa/rgasc_nsa);
           xdens     = 1.0_rp/(rgasc_nsa*tempe_nsa)*press_hyd;

        else if( kfl_benme_nsa == 5 .or. kfl_benme_nsa == 15 .or. &
                 kfl_benme_nsa == 210 .or. kfl_benme_nsa == 205) then

           !
           ! NH Linear mountain or inertia gravity wave, or Grabowski JAS 2007
           !
           brunt  = 0.01_rp
           brunt2 = brunt*brunt
           pi_hyd = 1.0_rp + (grnor_nsa * grnor_nsa / cpcoe_nsa / tempe_nsa / brunt2) &
                * (exp(- brunt2 * xcoor / grnor_nsa) - 1.0_rp)
           xtemp = tempe_nsa * exp(brunt2 * xcoor / grnor_nsa)
           xdens = pbaro_nsa * pi_hyd ** (cvcoe_nsa / rgasc_nsa) / rgasc_nsa / xtemp
        end if

     end if

  else if (itask == 2) then  
     !
     ! compute pressure
     !
     if (lawst_nsa == 4) then
        !
        ! meteo enters here:
        !
        if( kfl_benme_nsa == 1 .or. kfl_benme_nsa == 11 .or.  &
            kfl_benme_nsa == 2 .or. kfl_benme_nsa == 12) then
           !
           ! Warm bubble and density current
           !
           xpres = pbaro_nsa*(xdens*rgasc_nsa*xtemp/pbaro_nsa)**(cpcoe_nsa/cvcoe_nsa) 
           !OK (consistent with  that computed in nsa_initial_conditions.f90
           
        else if( kfl_benme_nsa == 3 .or. kfl_benme_nsa == 4 .or.  &
             kfl_benme_nsa == 203 .or. kfl_benme_nsa == 204) then
           !
           ! HS Linear mountain
           !
           pi_hyd = exp(-grnor_nsa/(cpcoe_nsa*tempe_nsa)*xcoor)
           xpres  = pbaro_nsa*pi_hyd**(cpcoe_nsa/rgasc_nsa)
           
        else if( kfl_benme_nsa == 5  .or. kfl_benme_nsa == 15 .or. &
                 kfl_benme_nsa == 210 .or. kfl_benme_nsa == 205) then
           !
           ! NH Linear mountain
           !
           brunt  = 0.01_rp
           brunt2 = brunt*brunt
           pi_hyd = 1.0_rp + (grnor_nsa * grnor_nsa / cpcoe_nsa / tempe_nsa / brunt2) &
                * (exp(- brunt2 * xcoor / grnor_nsa) - 1.0_rp)
           xpres  = pbaro_nsa * pi_hyd ** (cpcoe_nsa / rgasc_nsa)

        end if
     
     end if

  else if (itask == 3) then
     !
     ! compute temperature
     !
     if (lawst_nsa == 4) then
        !
        ! Meteo enters here
        !
        if( kfl_benme_nsa == 1 .or. kfl_benme_nsa == 11 .or.  &
            kfl_benme_nsa == 2 .or. kfl_benme_nsa == 12) then
           !
           ! Warm bubble and density current
           !
           xtemp = tempe_nsa
           
        else if( kfl_benme_nsa == 3 .or. kfl_benme_nsa == 4 .or. &
             kfl_benme_nsa == 203 .or. kfl_benme_nsa == 204) then
           !
           ! HS Linear mountain
           !
           pi_hyd   = exp(-grnor_nsa/(cpcoe_nsa*tempe_nsa)*xcoor);
           xtemp  = tempe_nsa/pi_hyd;

        else if( kfl_benme_nsa == 5 .or. kfl_benme_nsa == 15 .or. &
                 kfl_benme_nsa == 210 .or. kfl_benme_nsa == 205) then
           !
           ! NH Linear mountain
           !
           brunt  = 0.01_rp
           brunt2 = brunt*brunt
           xtemp = tempe_nsa * exp(brunt2 * xcoor / grnor_nsa)
        end if

     end if
     
  end if

end subroutine nsa_compute_hyd
