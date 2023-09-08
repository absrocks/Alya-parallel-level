subroutine chm_hdiffu(&
     iclas,pnode,pelty,pgaus,gpcar,elvel,elcod,gpden,gpdif)
  !-----------------------------------------------------------------------
  !****f* partis/chm_hdiffu
  ! NAME 
  !    chm_hdiffu
  ! DESCRIPTION
  !    Compute horizontal diffusion
  ! USES
  ! USED BY
  !    chm_elmpro
  !***
  ! NOTE: diffu_chm contains the diffusion coefficient (m^2/2). In the 
  !       conservative form we solve using densi*diffu (kg/m.s). So, it
  !       in necessary to multiply be density...   
  !
  !-----------------------------------------------------------------------
  use def_kintyp, only      :  ip,rp
  use def_domain, only      :  ndime,mgaus,mnode,hnatu,elmar
  use def_chemic, only      :  lawdi_chm,diffu_chm
  implicit none
  integer(ip),  intent(in)  :: iclas,pnode,pelty,pgaus
  real(rp),     intent(in)  :: elvel(ndime,pnode)
  real(rp),     intent(in)  :: elcod(ndime,pnode)
  real(rp),     intent(in)  :: gpcar(ndime,mnode,pgaus)
  real(rp),     intent(in)  :: gpden(pgaus)
  real(rp),     intent(out) :: gpdif(ndime,pgaus)
  integer(ip)               :: igaus,idime,inode,ipoin
  real(rp)                  :: tragl(9),hleng(3),delta,csh,km,khn0,khn,khf,kht,alfa
  real(rp)                  :: dvxdx,dvxdy,dvydx,dvydy
  
  do igaus = 1,pgaus
     do idime = 1,ndime-1
        gpdif(idime,igaus) = 0.0_rp
     end do
  end do

  if( lawdi_chm(1,iclas) == 1 ) then
     !
     ! Constant
     !
     do igaus = 1,pgaus
        do idime = 1,ndime-1
           gpdif(idime,igaus) = diffu_chm(1,iclas)*gpden(igaus)  ! multiply by densi
        end do
     end do
 
  else if( lawdi_chm(1,iclas) == 2 ) then
     !
     ! RAMS
     !
     call elmle1(&
          ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
          hnatu(pelty),hleng)
     if( ndime == 2 ) then
        delta = hleng(1)
     else
        delta = hleng(1)*hleng(2)
     end if

     csh = 0.2275_rp

     if(ndime==3) then
        km = 0.075_rp*(delta**(2.0_rp/3.0_rp))
        do igaus=1,pgaus
           dvxdx = 0.0_rp
           dvxdy = 0.0_rp
           dvydx = 0.0_rp
           dvydy = 0.0_rp
           do inode=1,pnode
              dvxdx = dvxdx + gpcar(1,inode,igaus)*elvel(1,inode)
              dvxdy = dvxdy + gpcar(2,inode,igaus)*elvel(1,inode)
              dvydx = dvydx + gpcar(1,inode,igaus)*elvel(2,inode)
              dvydy = dvydy + gpcar(2,inode,igaus)*elvel(2,inode)        
           end do
           gpdif(1,igaus) = csh*csh*delta*sqrt( (dvxdy+dvydx)**2 + 2.0_rp*(dvxdx*dvxdx + dvydy*dvydy) )
!
!          Here we assume Pr=1
!
           gpdif(1,igaus) = gpden(igaus)*(1.0_rp*max(km,gpdif(1,igaus)))  
           gpdif(2,igaus) = gpdif(1,igaus) 
        end do
     else
        call runend('2D HORIZ. DIFFUSION TO BE CODED')
     end if

  else if( lawdi_chm(1,iclas) == 3 ) then
     !
     ! CMAQ. For large grid sizes counteract the numerical diffusion
     !       It is assumed Kf=8000 m2/s for Dx=4km = 3.6d-2 deg of meridian
     !
     call elmle1(&
          ndime,pnode,elmar(pelty)%dercg,tragl,elcod,&
          hnatu(pelty),hleng)
     if( ndime == 2 ) then
        delta = hleng(1)
     else
        delta = hleng(1)*hleng(2)
     end if
    
     khn0 = 4000_rp 
     Khf  = 8000.0_rp             ! eddy diffusivity at a fixed resolution (CMAQ)
     alfa = 0.28_rp            ! Smagorinsky parameter in L.E. Kh formula
         
     if(ndime==3) then
        do igaus=1,pgaus
           dvxdx = 0.0_rp
           dvxdy = 0.0_rp
           dvydx = 0.0_rp
           dvydy = 0.0_rp
           do inode=1,pnode
              dvxdx = dvxdx + gpcar(1,inode,igaus)*elvel(1,inode)
              dvxdy = dvxdy + gpcar(2,inode,igaus)*elvel(1,inode)
              dvydx = dvydx + gpcar(1,inode,igaus)*elvel(2,inode)
              dvydy = dvydy + gpcar(2,inode,igaus)*elvel(2,inode)        
           end do
           kht = alfa*alfa*delta*sqrt( (dvxdx-dvydy)**2 + (dvydx+dvxdy)**2 )
           khn = Khf*Khn0*Khn0/delta
           gpdif(1,igaus) = (1.0_rp/Kht) + (1.0_rp/Khn)
           gpdif(1,igaus) = gpden(igaus)*max(1.0_rp/gpdif(1,igaus),1d2)
           gpdif(2,igaus) = gpdif(1,igaus) 
        end do
     else
        call runend('2D HORIZ. DIFFUSION TO BE CODED')    
     end if
     
  end if
 
end subroutine chm_hdiffu
