!-----------------------------------------------------------------------
!> @addtogroup Solidz
!> @{
!> @file    sld_velsnd.f90
!> @author  Solidz Team
!> @date
!> @brief   This routine computes the velocity of sound for each material
!> @details
!>          The speed of sound is calculated by:
!>             c = \sqrt(E/\rho)
!>
!>          References:\n
!>          G. Woan. Cambridge handbook of physics formulas
!>
!> @}
!-----------------------------------------------------------------------

subroutine sld_velsnd(pmate)

  use def_kintyp,  only : ip, rp
  use def_master,  only : ittim, kfl_rstar
  use def_domain,  only : ndime
  use def_solidz,  only : parco_sld, densi_sld
  use def_solidz,  only : lawst_sld, lawco_sld, lawpl_sld
  use def_solidz,  only : velas_sld
  use def_solidz,  only : kfl_plane_sld
  use def_solidz,  only : stiff0_sld

  implicit none

  integer(ip), intent(in)    :: pmate     !< Material code
  real(rp)                   :: young,poiss,lame1,lame2
  real(rp)                   :: lambda0,mu0,K0

  if (ittim == 0_ip .or. kfl_rstar /= 0 ) then

     !
     ! Speed of sound for each material model
     !
     if (lawst_sld(pmate)==100) then
        !
        ! sm100
        !
        young = parco_sld(1,pmate)
        poiss = parco_sld(2,pmate)
        lame1 = parco_sld(3,pmate) ! lambda
        lame2 = parco_sld(4,pmate) ! mu
        if ( ndime == 2_ip .and. kfl_plane_sld == 1_ip ) then  ! 2-d plane stress
           ! E = young/(1-poiss^2)
           velas_sld(1,pmate) = sqrt(young/(1.0_rp-poiss**2)/densi_sld(1,pmate))
        else                                                   ! 3-d and 2-d plane strain
           ! E = lambda + 2*mu
           velas_sld(1,pmate) = sqrt((lame1+2.0_rp*lame2)/densi_sld(1,pmate))
        end if

     else if (lawst_sld(pmate)==101) then
        lambda0 = parco_sld(1,pmate)
        mu0     = parco_sld(2,pmate)
        velas_sld(1,pmate) = sqrt((lambda0+2.0_rp*mu0)/densi_sld(1,pmate))
     else if (lawst_sld(pmate)==102) then
        mu0  = parco_sld(1,pmate)
        K0   = parco_sld(2,pmate)
        velas_sld(1,pmate) = sqrt((K0+4.0_rp*mu0/3.0_rp)/densi_sld(1,pmate))
     else if (lawst_sld(pmate)==103) then
        velas_sld(1,pmate) = sqrt(parco_sld(1,pmate)/densi_sld(1,pmate))
     else if (lawst_sld(pmate)==104) then
        velas_sld(1,pmate) = sqrt(parco_sld(1,pmate)/densi_sld(1,pmate))
     else if (lawst_sld(pmate)==105) then
        lambda0 = parco_sld(1,pmate)
        velas_sld(1,pmate) = sqrt(lambda0/densi_sld(1,pmate))
     else if (lawst_sld(pmate)==133) then
        velas_sld(1,pmate) = sqrt(parco_sld(1,pmate)/densi_sld(1,pmate))
     else if (lawst_sld(pmate)==134) then
        velas_sld(1,pmate) = sqrt(parco_sld(1,pmate)/densi_sld(1,pmate))
     else if (lawst_sld(pmate)==135) then
        velas_sld(1,pmate) = sqrt(parco_sld(1,pmate)/densi_sld(1,pmate))
     else if (lawst_sld(pmate)==136) then
        velas_sld(1,pmate) = sqrt(parco_sld(1,pmate)/densi_sld(1,pmate))
     else if (lawst_sld(pmate)==137) then
        velas_sld(1,pmate) = sqrt(parco_sld(1,pmate)/densi_sld(1,pmate))
     else if (lawst_sld(pmate)==151) then
        !
        ! sm151
        ! E = max(Cij) (largest value of the constitutive tensor)
        velas_sld(1,pmate) = sqrt(maxval(stiff0_sld(:,:,pmate))/densi_sld(1,pmate))

     else if (lawst_sld(pmate)==152) then
        !
        ! sm152
        ! E = max(Cij) (largest value of the constitutive tensor)
        velas_sld(1,pmate) = sqrt(maxval(stiff0_sld(:,:,pmate))/densi_sld(1,pmate))

     else if (lawst_sld(pmate)==153) then
        !
        ! sm153
        !
        young = parco_sld(1,pmate)
        poiss = parco_sld(2,pmate)
        lame1 = poiss*young/((1.0_rp+poiss)*(1.0_rp-2.0_rp*poiss))
        lame2 = young/(2.0_rp*(1.0_rp+poiss))
        ! E = lambda + 2*mu
        velas_sld(1,pmate) = sqrt((lame1+2.0_rp*lame2)/densi_sld(1,pmate))

     else if (lawst_sld(pmate)==154) then
        !
        ! sm154
        ! E = max(Cij) (largest value of the constitutive tensor)
        velas_sld(1,pmate) = sqrt(maxval(stiff0_sld(:,:,pmate))/densi_sld(1,pmate))

     else if (lawst_sld(pmate)==200) then
        !
        ! sm200
        !
        if (      lawco_sld(pmate)==1 ) then ! Isotropic
           young = parco_sld(1,pmate)
           poiss = parco_sld(2,pmate)
           lame1 = poiss*young/((1.0_rp+poiss)*(1.0_rp-2.0_rp*poiss))
           lame2 = young/(2.0_rp*(1.0_rp+poiss))
           ! E = lambda + 2*mu
           velas_sld(1,pmate) = sqrt((lame1+2.0_rp*lame2)/densi_sld(1,pmate))
        else if ( lawco_sld(pmate)==2 ) then   ! Orthotropic
           ! E = max(Cij) (largest value of the constitutive tensor)
           velas_sld(1,pmate) = sqrt(maxval(stiff0_sld(:,:,pmate))/densi_sld(1,pmate))
        else if ( lawco_sld(pmate)==3 ) then   ! Orthotropic damage
           ! E = max(Cij) (largest value of the constitutive tensor)
           velas_sld(1,pmate) = sqrt(maxval(stiff0_sld(:,:,pmate))/densi_sld(1,pmate))
        end if

     else if (lawpl_sld(pmate)==1 .and. lawst_sld(pmate)==0) then
         velas_sld(1,pmate) = sqrt(parco_sld(1,pmate)/densi_sld(1,pmate))

     else if (lawst_sld(pmate)==666) then
        velas_sld(1,pmate) = sqrt(parco_sld(1,pmate)/densi_sld(1,pmate))

     else
        call runend('SLD_VELSND: STRESS MODEL NOT CONSIDERED')
     end if

  end if

end subroutine sld_velsnd
