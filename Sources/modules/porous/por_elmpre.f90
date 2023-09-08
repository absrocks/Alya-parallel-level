!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_elmpre.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Compute some Gauss values to be used by por_elmpma (or in the future by elmadr)
!> @details Compute some Gauss values
!> OUTPUT 
!>    GPDEN .......... cop *B0/Bw + cwp -------- Despite it is not actually a Density but we still call it density
!>    GPDIF .......... Diffusion term ---------- (1/Bw)(k_ro/mu_o+k_rw/mu_w)K(idime)
!>    GPGRA .......... Gravity term
!> @} 
!------------------------------------------------------------------------

subroutine por_elmpre(&
     pnode,pgaus,gpsha,elpre,elswa,elcod,&
     gpcod,gpden,gpdif,gpgra,gppre,ielem)
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,mgaus,lmate
  use def_porous, only       :  kfl_timei_por,comro_por,comwa_por
  use def_porous, only       :  comoi_por,bwref_por,boref_por,prref_por
  use def_porous, only       :  muwat_por,muoil_por,grnor_por,gravi_por
  use def_porous, only       :  poro0_por,perme_por,denoi_por,denwa_por
  use def_porous, only       :  ncomp_por,kprsa_por,denhy_por
  implicit none 
  integer(ip), intent(in)    :: pnode,pgaus                  !> 
  integer(ip), intent(in)    :: ielem                        !> 
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)           !> gauss point shape function
  real(rp),    intent(in)    :: elpre(pnode,*)               !> pressure
  real(rp),    intent(in)    :: elswa(pnode,*)               !> Water Saturation
  real(rp),    intent(in)    :: elcod(ndime,pnode)           !> Coordinates
  real(rp),    intent(out)   :: gpcod(ndime,pgaus)           !> Coordinates
  real(rp),    intent(out)   :: gpden(pgaus)                 !> Term that goes in the density place
  real(rp),    intent(out)   :: gpdif(ndime,pgaus)           !> Term that goes in the diffusion place
  real(rp),    intent(out)   :: gpgra(ndime,pgaus)           !> Terms related to gravity
  real(rp),    intent(out)   :: gppre(mgaus,ncomp_por)       !> pressure
  integer(ip)                :: idime,inode,igaus
  real(rp)                   :: gpswa(mgaus,ncomp_por)       ! Water Saturation - we only need it inside por_elmpre
  real(rp)                   :: perme(ndime),auxii,deltp
  real(rp)                   :: xx_o,xx_w,xx_r,Bo,Bw,poro
  real(rp)                   :: kro,krw,dinvbo_dp,dinvbw_dp
  real(rp)                   :: c_wp,c_op,dummr

  if( kfl_timei_por(kprsa_por) /= 0 ) then
     do igaus = 1,pgaus
        gppre(igaus,2) = 0.0_rp
        do inode = 1,pnode
           gppre(igaus,2) = gppre(igaus,2) + gpsha(inode,igaus) * elpre(inode,2)
        end do
     end do
  end if

  do igaus = 1,pgaus
     !
     ! Pressure: GPPRE & Water Saturation: GPSWA & Coordinates
     !
     gppre(igaus,1) = 0.0_rp
     gpswa(igaus,1) = 0.0_rp
     do idime = 1,ndime
        gpcod(idime,igaus) = 0.0_rp
     end do
     do inode = 1,pnode
        gppre(igaus,1) = gppre(igaus,1) + gpsha(inode,igaus) * elpre(inode,1)
        gpswa(igaus,1) = gpswa(igaus,1) + gpsha(inode,igaus) * elswa(inode,1)
        do idime = 1,ndime
           gpcod(idime,igaus) = gpcod(idime,igaus)&
                + gpsha(inode,igaus) * elcod(idime,inode)
        end do
     end do
    !
     ! Preliminaries Bo, Bw, kro, krw, phi
     !
     deltp     = gppre(igaus,1) - prref_por
     xx_o      = comoi_por * deltp
     xx_w      = comwa_por * deltp
     xx_r      = comro_por * deltp
     Bo        = boref_por / (1.0_rp + xx_o + (0.5_rp*xx_o*xx_o) )
     Bw        = bwref_por / (1.0_rp + xx_w + (0.5_rp*xx_w*xx_w) )
     poro      = poro0_por(ielem) * (1.0_rp + xx_r + (0.5_rp*xx_r*xx_r) )
     call por_gtable(gpswa(igaus,1),kro,krw,dummr,lmate(ielem),0_ip) ! obtain kro & krw
     dinvBw_dp = ( 1.0_rp + xx_w ) * comwa_por / bwref_por 
     dinvBo_dp = ( 1.0_rp + xx_o ) * comoi_por / boref_por
     c_wp      = gpswa(igaus,1) * ( (poro0_por(ielem) * comro_por / Bw) + (poro * dinvBw_dp) )
     c_op      = ( 1.0_rp - gpswa(igaus,1) ) * ( (poro0_por(ielem) * comro_por / Bo) + (poro * dinvBo_dp) )
     perme     = perme_por(:,ielem)
     !
     ! Gravity term 
     !
     auxii = ( grnor_por / Bw) * ( ( kro * denoi_por / ( muoil_por * Bo ) ) + ( krw * denwa_por / ( muwat_por * Bw ) ) -   &
                ( denhy_por  *  ( ( kro / muoil_por ) + ( krw / muwat_por ) ) ) )
     do idime = 1,ndime
        gpgra(idime,igaus) = perme(idime) * gravi_por(idime) * auxii
     end do
     !
     ! Diffusion terms -----------  (1/Bw)(k_ro/mu_o+k_rw/mu_w)K(idime)
     !
     auxii =  ( (kro / muoil_por) + (krw / muwat_por) ) / Bw
     do idime = 1,ndime
        gpdif(idime,igaus) =  perme(idime) * auxii
     end do
     !
     ! Transient terms -------- cop *B0/Bw + cwp -------- Despite it is not actually a Density but we still call it density
     !
     gpden(igaus) = c_wp + (c_op * Bo / Bw) 
  end do

end subroutine por_elmpre
