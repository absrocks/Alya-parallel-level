!------------------------------------------------------------------------
!> @addtogroup Porous 
!> @{
!> @file    por_elmprs.f90
!> @date    13/05/2013
!> @author  Herbert Owen
!> @brief   Compute some Gauss values to be used by elmadr to obtain the Saturation equation
!> @details Compute some Gauss values
!! OUTPUT 
!!    GPADV .......... the 'velocity' to be used in the advective term
!!    GPRHS .......... the terms corresponding to -f_w div(rho_w * u_star) except the one that goes multiplied    
!!    by grad(S_w) that I put inside GPADV
!> @} 
!------------------------------------------------------------------------
subroutine por_elmprs(&
     pnode,pgaus,&
     gpsha,gpcar,elpre,elswa,elvel,elcod,&
     gpcod,gpadv,gprhs,ielem)
!
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  mnode,ndime,lmate
  use def_porous, only       :  comro_por,comwa_por,comoi_por,bwref_por,boref_por,prref_por,muwat_por, &
                                muoil_por,grnor_por,gravi_por,poro0_por,perme_por,denoi_por,denwa_por,kfl_malta_por
  implicit none 
  integer(ip), intent(in)    :: pnode                        !> number of nodes
  integer(ip), intent(in)    :: pgaus                        !> number of gauss points
  integer(ip), intent(in)    :: ielem                        !> element number 
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)           !> gauss point shape function
  real(rp),    intent(in)    :: gpcar(ndime,mnode,pgaus)     !> gauss point shape function derivatives
  real(rp),    intent(in)    :: elpre(pnode,*)               !> pressure
  real(rp),    intent(in)    :: elswa(pnode,*)               !> water saturation
  real(rp),    intent(in)    :: elvel(ndime,pnode,*)         !> velocity
  real(rp),    intent(in)    :: elcod(ndime,pnode)           !> Coordinates
  real(rp),    intent(out)   :: gpcod(ndime,pgaus)           !> Coordinates
  real(rp),    intent(out)   :: gpadv(ndime,pgaus)           !> Advection velocity
  real(rp),    intent(out)   :: gprhs(pgaus)                 !> RHS

  real(rp)    :: gpswa(pgaus,1)                              ! Sw - for the moment I see no need to take it out of this sub
  real(rp)    :: gppre(pgaus,1)                              ! pressure - for the moment I see no need to take it out of this sub
  real(rp)    :: gpgrp(ndime,pgaus)                 ! pressure gradient - for the moment I see no need to take it out of this sub
  real(rp)    :: gpvel(ndime,pgaus)                   
  real(rp)    :: gpdiv(pgaus)                   
  real(rp)    :: auxi1,auxi2,auxi3,auxi4          

  integer(ip)                :: idime,inode,igaus
  real(rp)                   :: perme(ndime),deltp,xx_o,xx_w,xx_r,Bo,Bw,poro,kro,krw,dinvbo_dp,dinvbw_dp
  real(rp)                   :: lambd_o,lambd_w,f_w,f_o,gpvst(ndime),df_w,dkro

  do igaus = 1,pgaus
     !
     ! Pressure, div(velo), Coordinates, Pressure gradient & Velocity
     !
     gpswa(igaus,1) = 0.0_rp
     gppre(igaus,1) = 0.0_rp
     gpdiv(igaus) = 0.0_rp
     do idime = 1,ndime
        gpcod(idime,igaus) = 0.0_rp
        gpgrp(idime,igaus) = 0.0_rp
        gpvel(idime,igaus) = 0.0_rp
     end do
     do inode = 1,pnode
        gpswa(igaus,1) = gpswa(igaus,1) + gpsha(inode,igaus) * elswa(inode,1)
        gppre(igaus,1) = gppre(igaus,1) + gpsha(inode,igaus) * elpre(inode,1)
        do idime = 1,ndime
           gpcod(idime,igaus) = gpcod(idime,igaus)&
                + gpsha(inode,igaus) * elcod(idime,inode)
           gpgrp(idime,igaus) = gpgrp(idime,igaus)&
                + gpcar(idime,inode,igaus) * elpre(inode,1)
           gpvel(idime,igaus) = gpvel(idime,igaus)&
                + gpsha(inode,igaus) * elvel(idime,inode,1)
           gpdiv(igaus) = gpdiv(igaus)&
                + gpcar(idime,inode,igaus) * elvel(idime,inode,1)
        end do
     end do
     !
     ! Preliminaries Bo, Bw, kro, krw, dkro, phi , u_star , f_w , lambd_o , lambd_w
     !
     deltp     = gppre(igaus,1) - prref_por
     xx_o      = comoi_por * deltp
     xx_w      = comwa_por * deltp
     xx_r      = comro_por * deltp
     Bo        = boref_por / (1.0_rp + xx_o + (0.5_rp*xx_o*xx_o) )
     Bw        = bwref_por / (1.0_rp + xx_w + (0.5_rp*xx_w*xx_w) )
     poro      = poro0_por(ielem) * (1.0_rp + xx_r + (0.5_rp*xx_r*xx_r) )
     call por_gtable(gpswa(igaus,1),kro,krw,dkro,lmate(ielem),1_ip) ! obtain kro & krw & dkro
     dinvBw_dp = ( 1.0_rp + xx_w ) * comwa_por / bwref_por 
     dinvBo_dp = ( 1.0_rp + xx_o ) * comoi_por / boref_por
     lambd_o   = kro / muoil_por
     lambd_w   = krw / muwat_por
     f_w       = lambd_w  / (lambd_w + lambd_o)
     f_o       = lambd_o  / (lambd_w + lambd_o)
     !
     ! GPVST = u* = u + k * lambda0 * ( rho_w  - rho_o) g
     !
     perme(1:ndime) = perme_por(1:ndime,ielem)
     auxi4 = ( denwa_por/Bw ) - ( denoi_por/Bo )          ! rho_w - rho_o
     do idime = 1,ndime
        gpvst(idime) = gpvel(idime,igaus) + perme(idime) * lambd_o * gravi_por(idime) * grnor_por * auxi4
     end do
     !
     ! AUX1 = u* . grad(p)
     ! AUX2 = k g . grad(p)
     ! AUX3 = rho_w * dBw^-1 / dp - rho_o * dBo^-1 / dp  
     !      
     auxi1 = 0.0_rp
     auxi2 = 0.0_rp
     do idime = 1,ndime
        auxi1 = auxi1 + gpvst(idime) * gpgrp(idime,igaus) 
        auxi2 = auxi2 + gravi_por(idime) * gpgrp(idime,igaus) * perme(idime) 
     end do
     auxi2 = auxi2 * grnor_por      
     !
     ! RHS
     !      
     gprhs(igaus) = 0.0_rp

     if ( kfl_malta_por == 0 ) then 

        auxi3 = ( (denwa_por * dinvBw_dp) - (denoi_por * dinvBo_dp) ) / Bw
        gprhs(igaus) = gprhs(igaus)    &
             - f_w * denwa_por         &            ! -f_w*rho_w * [
             * (  (gpdiv(igaus) / Bw ) &            ! 1/Bw*div(u)
             +    (dinvBw_dp * auxi1)  &            ! dBw^-1 / dp * AUX1
             +    ( lambd_o * auxi3 * auxi2 ) )     ! lambda_o * AUX2 * AUX3 ]

     else if ( kfl_malta_por == 1 ) then 
        !
        ! Beware this case has only been prepared for g=0 for the moment
        !
        auxi3 = ( ( (2.0_rp * denwa_por * dinvBw_dp) - (denoi_por * dinvBo_dp) ) / Bw )  &
             -  (denoi_por * dinvBw_dp)  / Bo
                   
        gprhs(igaus) = gprhs(igaus) - f_w * denwa_por * ( f_o * (    &
             ( dinvBw_dp - ( dinvBo_dp * Bo / Bw ) ) * auxi1 ) + ( lambd_o * auxi3 * auxi2 ) )

     end if
     !
     ! Advection velocity for convective term
     !   rho_w / Bw * ( - df_w/dS_w * u* 
     ! + g * f_w * k * ( rho_w - rho_o ) * dkr_o/mu_o
     !
     df_w = - dkro / ( muoil_por * ( ( lambd_o + lambd_w )**2_ip ) )
     do idime = 1,ndime
        gpadv(idime,igaus) =  ( denwa_por / Bw ) &
             * ( ( gpvst(idime) * df_w ) &
             +   ( gravi_por(idime) * grnor_por * f_w * perme(idime) &
             *     auxi4 * dkro / muoil_por ) )
     end do
  end do

end subroutine por_elmprs
