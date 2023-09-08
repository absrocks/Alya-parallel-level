subroutine chm_stalaw(itask,xdens,xpres,xtemp,xcoor,xmixv)
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
  use def_chemic
  use def_kermod, only     : gasco
  implicit none

  integer(ip) :: itask,ielem,ipoin,inode,pelty,pnode,pgaus
  real(rp)    :: xdens,xpres,xtemp,xcoor,xcosa,rgacp,xinve,xfact,vari1,vari2,xfac2,xauxi,bruva
  real(rp)    :: xmixv    !Water vapor mixing ratio r_v, S.M.
  real(rp)    :: epsi     !Ratio R/Rv S.M.
  real(rp)    :: Lv       !Latent heat of vaporization: Lv = Lv0 - (Cpl - Cpv)*(T - T0)

  ! Initialize bruva:
  bruva= 0.0_rp
  
  afact_chm= pbaro_chm ** ((adgam_chm - 1.0_rp) / adgam_chm)
  afact_chm= gasco / afact_chm

  afact_chm= pbaro_chm ** (adgam_chm - 1.0_rp)
  afact_chm= (gasco ** adgam_chm) / afact_chm    ! C_o in ahmad-lindeman

  if (itask == 1 .or. itask == -1) then        ! compute density
     
     xdens= (xpres/afact_chm)**(1.0_rp/adgam_chm)
     xdens= xdens/xtemp

  else if (itask == 2) then   ! compute pressure

     xpres= (xdens * xtemp)**adgam_chm 
     xpres= afact_chm * xpres

  else if (itask == 3) then   ! compute temperature

     xtemp= (xpres/afact_chm)**(1.0_rp/adgam_chm)
     xtemp= xtemp/xdens   

  end if
  
end subroutine chm_stalaw
