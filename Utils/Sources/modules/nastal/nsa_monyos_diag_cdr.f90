subroutine nsa_monyos_diag_cdr(&
     ielem,igaus,xresi,xsube,xortp,chale,hleng,xunkn,gunkn,gpres,gvisc,dvelo,velmo,&
     xvelo,xvisc,xdith,sound,xlade,xldve,xdtix,taudi,qufac,conme,difme)
  !-----------------------------------------------------------------------
  !***** nastal/nsa_monyos
  ! NAME 
  !    nsa_monyos
  ! DESCRIPTION
  ! USES
  !    nsa_...
  ! USED BY
  !    nsa_...
  !***
  !-----------------------------------------------------------------------
  use      def_master
  use      def_domain
  use      def_parame
  use      def_nastal
  implicit none

  integer(ip)  :: ndofc,ndofe,idime,jdime,igaus,idofn,jdofn,ifreq,jfreq,kfreq,nfreq,nfrin,itisu,ntisu,ielem

  real(rp)    :: &
       taudi(ndofn_nsa),xsube(ndofn_nsa,mgaus,3),xortp(ndofn_nsa,mgaus),dsube(ndofn_nsa), xsaux(ndofn_nsa,mgaus),&
       dundt(ndofn_nsa),spare(ndofn_nsa),emean(6+5*ndime+4*ntens),&
       hconv,htran,hleti,hsoun,chale(2),seige,veige,reate,&
       xunkn(ndofn_nsa,mgaus,3), gunkn(ndofn_nsa,ndime,mgaus),&
       gpres(ndime,mgaus,ncomp_nsa),xresi(ndofn_nsa,mgaus),xdesi(ndofn_nsa,mgaus),&
       xvelo(ndime,mgaus,ncomp_nsa),gvisc(ndime,mgaus,ncomp_nsa),xlade(mgaus,ncomp_nsa),xldve(ndime,mgaus,ncomp_nsa),&
       vitot,ditot,sound(mgaus),velmo(mgaus),xvisc(mgaus),xdith(mgaus),dvelo,qufac,tauen,taupa,tauco,hleng(ndime),&
       conme(ndofn_nsa,ndofn_nsa,ndime),difme(ndofn_nsa,ndofn_nsa,ndime,ndime),&
       emene,emvel,emve2,emso2,emve3,emesv,perio(ndime),omega(ndime)

  real(rp)     :: &
       tauxi(ndofn_nsa,ndofn_nsa),xdtix(ndofn_nsa,mgaus,2),xauxi,xsute,xunte,&
       vesta(ndofn_nsa,ndime),vista(ndofn_nsa,ndime),rasta(ndofn_nsa),absve(ndofn_nsa),absvi(ndofn_nsa), &
       gseig(ndime),gdens(ndime),xggde,xggvd,tiint,tiine,tiinc,xdens,deltat(ndofn_nsa),frequ(ndime), &
       taure(ndofn_nsa,ndofn_nsa),tauim(ndofn_nsa,ndofn_nsa), &
       adtre(ndofn_nsa,ndofn_nsa),adtim(ndofn_nsa,ndofn_nsa),detre,detim, &
       detmo
  real(rp)     :: unacosa, otracosa


  complex(rp)  :: lfo1d_z(ndofn_nsa-1,ndofn_nsa-1),lfi1d_z(ndofn_nsa-1,ndofn_nsa-1),detma_z
  complex(rp)  :: lfour_z(ndofn_nsa  ,ndofn_nsa  ),lfinv_z(ndofn_nsa  ,ndofn_nsa  )

  nfrin= 1
  nfreq= 1

  reate = (1.0_rp + rgasc_nsa / cvcoe_nsa)*dvelo

  if (dvelo < 0.0_rp) then
     reate= -reate
!     reate= 0.0_rp
  end if

  ! chales are defined in nsa_elconsxy

  hconv = chale(1)
  hsoun = chale(1)
  htran = chale(2)


  hleti = hconv*htran

  ndofc= ndime+1
  ndofe= ndime+2

  xdens        = xunkn(ndofc,      igaus,1)
  gdens(    1) = gunkn(ndofc,    1,igaus)
  gdens(    2) = gunkn(ndofc,    2,igaus)
  gdens(ndime) = gunkn(ndofc,ndime,igaus)


  tauxi= 0.0_rp

  if (kfl_relat_nsa == 1) then
     tauxi(1,1) = abs((velmo(igaus)*(1.0_rp - sound(igaus)*sound(igaus)) &
          + sound(igaus)*(1.0_rp-velmo(igaus)*velmo(igaus))) &
          / (1.0_rp - velmo(igaus)*sound(igaus)*velmo(igaus)*sound(igaus)))
     tauxi(1,2) = abs((velmo(igaus)*(1.0_rp - sound(igaus)*sound(igaus)) &
          - sound(igaus)*(1.0_rp-velmo(igaus)*velmo(igaus))) &
          / (1.0_rp - velmo(igaus)*sound(igaus)*velmo(igaus)*sound(igaus)))
     sound(igaus) = tauxi(1,1)

     if (tauxi(1,2) > sound(igaus)) sound(igaus) = tauxi(1,2)
     sound(igaus)= 1.0_rp / sound(igaus)

     velmo(igaus) = 0.0_rp
     reate = 0.0_rp
  end if

  xauxi= 0.5_rp + cvcoe_nsa/adgam_nsa/rgasc_nsa
  xauxi= 1.0_rp / xauxi

  ! c = sqrt(gamma R T) = sqrt(gamma p / rho) => p = c*c rho / gamma

  if (kfl_taufa_nsa(1,2) == 1) then
     tauxi(1,1) = velmo(igaus)/hconv
     tauxi(2,1) = velmo(igaus)/hconv
  end if
  if (kfl_taufa_nsa(2,2) == 1) then
     tauxi(1,2) = reate
     tauxi(2,2) = reate
  end if
  if (kfl_taufa_nsa(3,2) == 1) then
     tauxi(1,3) = sound(igaus) / hleng(ndime)
     tauxi(2,3) = sound(igaus) / hleng(ndime)
  end if
  if (kfl_taufa_nsa(4,2) == 1) then
     tauxi(1,4) = 4.0_rp * xvisc(igaus)/hleti / xdens
     tauxi(2,4) = 4.0_rp * xdith(igaus)/hleti / xdens / cpcoe_nsa
  end if

  tiint= tauxi(1,1) + tauxi(1,2) + tauxi(1,3) + tauxi(1,4)
  tiine= tauxi(2,1) + tauxi(2,2) + tauxi(2,3) + tauxi(2,4)
!!  tiinc= tiint                                 ! this works, but it is not logical
  tiinc= tauxi(1,1) + tauxi(1,2) + tauxi(1,3)    ! this works too, and it is more logical, so it is the best option
!!  tiinc= tauxi(1,1) + tauxi(1,2)               ! this DOES NOT work!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  tiinc = 0.0_rp
  do idime=1,ndime
     tiinc = tiinc + conve_nsa(idime) * conve_nsa(idime)
  end do
  tiinc = sqrt(tiinc) / hconv
  tiinc = tiinc + react_nsa + diffu_nsa / hleng(ndime)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  tauen= 0.0_rp
  taupa= 0.0_rp
  tauco= 0.0_rp

  if (tiine > zensa) tauen= qufac / tiine
  if (tiint > zensa) taupa= qufac / tiint
  if (tiinc > zensa) tauco= qufac / tiinc

  if (kfl_stabi_nsa == 2) then
     call runend('NSA_MONYOS: CG DEPRECATED - USE MULTISCALE STABILIZATION!!')
  end if

  taudi(1:ndime) = taupa
  taudi(ndofc) = tauco
  taudi(ndofe) = tauen

  if (kfl_isent_nsa == 1) taudi(ndofe)=0.0_rp

  xsute= 0.0_rp
  xunte= 0.0_rp

  do idofn=1,ndofn_nsa
!     xsute = xsube(idofn,igaus,2) / xdtix(idofn,igaus,2)  ! it does not work with this term     
     xsube(idofn,igaus,1) = taudi(idofn) * ( xsute + xresi(idofn,igaus) - xunte) 
  end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do idofn=1,ndofn_nsa
     xsube(idofn,igaus,1) = 0.0_rp
     xresi(idofn,igaus) = 0.0_rp
  end do
  do idime=1,ndime
     !POSAR LA DIFUSSIO!!!!!!!!!!!!!!!!
     xresi(ndime+1,igaus) = xresi(ndime+1,igaus) - conve_nsa(idime) * gunkn(ndime+1,idime,igaus)          
  end do
  xresi(ndime+1,igaus) = xsute + xresi(ndime+1,igaus) + xortp(ndime+1,igaus)

  xsube(ndime+1,igaus,1) = xresi(ndime+1,igaus) * tauco

!!$print*
!!$print*,'tauco',ielem,tauco,qufac,1.0_rp/tiinc
!!$stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end subroutine nsa_monyos_diag_cdr
