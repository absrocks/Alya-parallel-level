subroutine tur_elmpro(&
     pnode,pgaus,lnods,elfle,gpsha,hleng,gpfle,gpden,gpvis)
  !------------------------------------------------------------------------
  !****f* Turbul/tur_elmpro
  ! NAME 
  !    tur_elmpro
  ! DESCRIPTION
  !    Compute element properties
  ! OUTPUT
  !    GPDEN ... Density at Gauss point
  !    GPVIS ... Viscosity at Gauss point
  ! USES
  ! USED BY
  !    tur_elmope
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_master, only     :  densi,visco
  use def_kermod, only     :  thicl
  use def_turbul, only     :  lawde_tur,lawvi_tur,densi_tur,visco_tur,&
       &                      densa_tur,visca_tur,kfl_colev_tur
  implicit none
  integer(ip), intent(in)  :: pnode,pgaus
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(in)  :: elfle(pnode)
  real(rp),    intent(in)  :: gpsha(pnode,pgaus)
  real(rp),    intent(in)  :: hleng(*)
  real(rp),    intent(out) :: gpfle(*)
  real(rp),    intent(out) :: gpden(pgaus)
  real(rp),    intent(out) :: gpvis(pgaus)
  integer(ip)              :: igaus,inode,ipoin
  real(rp)                 :: dens1,visc1,hlen1
  !
  ! Two-phase flow
  !
  if( kfl_colev_tur /=0 ) then

     do igaus = 1,pgaus
        gpfle(igaus) = 0.0_rp
        do inode = 1,pnode
           gpfle(igaus) = gpfle(igaus) + elfle(inode) * gpsha(inode,igaus)
        end do
     end do
     dens1 = densi_tur(1)
     visc1 = visco_tur(1)
     hlen1 = hleng(1)
     call coleve(&
          '11',1_ip,pgaus,pgaus,gpfle,dens1,densa_tur,&
          visc1,visca_tur,thicl,hlen1,gpden,gpvis)
  else
     !
     ! GPDEN: Density rho 
     !
     if(lawde_tur==0) then

        do igaus=1,pgaus
           gpden(igaus)=densi_tur(1)
        end do

     else if(lawde_tur==1) then

        do igaus=1,pgaus
           gpden(igaus)=0.0_rp
           do inode=1,pnode
              ipoin=lnods(inode)
              gpden(igaus)=gpden(igaus)+densi(ipoin,1)*gpsha(inode,igaus)
           end do
        end do

     end if
     !
     ! GPVIS: Viscosity mu
     !
     if(lawvi_tur==0) then

        do igaus=1,pgaus
           gpvis(igaus)=visco_tur(1)
        end do

     else if(lawvi_tur==1) then

        do igaus=1,pgaus
           gpvis(igaus)=0.0_rp
           do inode=1,pnode
              ipoin=lnods(inode)
              gpvis(igaus)=gpvis(igaus)+visco(ipoin,1)*gpsha(inode,igaus)
           end do
        end do

     end if
  end if

end subroutine tur_elmpro

