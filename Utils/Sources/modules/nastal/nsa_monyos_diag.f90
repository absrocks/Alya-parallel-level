subroutine nsa_monyos_diag(&
     ielem,igaus,xresi,xsube,chale,hleng,xunkn,gunkn,gpres,dvelo,velmo,&
     xvisc,xdith,sound,sspee,taudi,qufac,xhecp,xmowe,heats)
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
  use      def_kermod
  use      mod_ker_proper

  implicit none

  integer(ip)  :: ndofc,ndofe,idime,jdime,idofn,jdofn
  integer(ip), intent(in) :: ielem,igaus

  real(rp),  intent(in) :: xhecp,xresi(ndofn_nsa,mgaus),chale(2),hleng(ndime),xunkn(ndofn_nsa,mgaus,3),&
                           gunkn(ndofn_nsa,ndime,mgaus),xmowe,dvelo,xdith(mgaus),xvisc,&
                           sound(mgaus),sspee,qufac,gpres(ndime,mgaus,ncomp_nsa),heats
  real(rp),  intent(out) :: xsube(ndofn_nsa,mgaus,3),taudi(ndofn_nsa)
  real(rp),  intent(inout) :: velmo(mgaus)
  real(rp)    ::           hconv,hleti,hsoun,hmini,hmaxi,reate, &
                           tauen,taupa,tauco,taush, & 
                           tauxi(ndofn_nsa,5),xauxi,xsute,xunte,&
                           gdens(ndime),tiint,tiine,tiinc,xdens, &
                           rgasc,prope_tmp(1),xhecv,adgam,cfl_veloc,dtinv_phys_local, &
                           xmach,xmach_ref,cfl_sound,sound_pseud,a,b,zensa_lopre
  

  rgasc = runiv_nsa / xmowe   
  xhecv = xhecp - rgasc

  reate = (1.0_rp + rgasc / xhecv)*dvelo 

  reate = reate + heats  ! add heat source term from chemical reactions
                              ! heats = 0.0 if chemistry off
  if (dvelo < 0.0_rp) then
     reate= -reate
!     reate= 0.0_rp
  end if

  ! chales are defined in nsa_elconsxy

  hmini = chale(1) ! the smallest
  hmaxi = chale(2) ! the biggest
  hconv = hmini
  hsoun = hmini
  hleti = hmini*hmaxi
  if (kfl_higha_nsa == 1) hleti= hmini*hmini

  ndofc= ndime+1
  ndofe= ndime+2

  xdens        = xunkn(ndofc,      igaus,1)
  gdens(    1) = gunkn(ndofc,    1,igaus)
  gdens(    2) = gunkn(ndofc,    2,igaus)
  gdens(ndime) = gunkn(ndofc,ndime,igaus)

  tauxi= 0.0_rp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! DEPRECATED, BUT DO NOT ERASE.
!  if (kfl_relat_nsa == 1) then
!     tauxi(1,1) = abs((velmo(igaus)*(1.0_rp - sound(igaus)*sound(igaus)) &
!          + sound(igaus)*(1.0_rp-velmo(igaus)*velmo(igaus))) &
!          / (1.0_rp - velmo(igaus)*sound(igaus)*velmo(igaus)*sound(igaus)))
!     tauxi(1,2) = abs((velmo(igaus)*(1.0_rp - sound(igaus)*sound(igaus)) &
!          - sound(igaus)*(1.0_rp-velmo(igaus)*velmo(igaus))) &
!          / (1.0_rp - velmo(igaus)*sound(igaus)*velmo(igaus)*sound(igaus)))
!     sound(igaus) = tauxi(1,1)
!     if (tauxi(1,2) > sound(igaus)) sound(igaus) = tauxi(1,2)
!     sound(igaus)= 1.0_rp / sound(igaus)
!
!     velmo(igaus) = 0.0_rp
!     reate = 0.0_rp
!  end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  adgam = xhecp / xhecv

  xauxi= 0.5_rp + xhecv/adgam / rgasc
  xauxi= 1.0_rp / xauxi

  ! c = sqrt(gamma R T) = sqrt(gamma p / rho) => p = c*c rho / gamma

  zensa_lopre= 1.0e-6
  if(kfl_lopre_nsa > 1 .and. velmo(igaus) < zensa_lopre) then           
     velmo(igaus) = 0.0_rp
  end if

  if (kfl_taufa_nsa(1,2) == 1) then
     tauxi(      1,1) = velmo(igaus)/hconv
     tauxi(ndime+1,1) = velmo(igaus)/hconv
     tauxi(ndime+2,1) = velmo(igaus)/hconv
  end if

  if (kfl_reate_nsa == 0) reate=0.0_rp
  if (kfl_taufa_nsa(2,2) == 1) then
     tauxi(      1,2) = reate
     tauxi(ndime+1,2) = reate
     tauxi(ndime+2,2) = reate
  end if

  if (kfl_taufa_nsa(3,2) == 1) then
     tauxi(      1,3) = sound(igaus) / hsoun
     tauxi(ndime+1,3) = sound(igaus) / hsoun
     tauxi(ndime+2,3) = sound(igaus) / hsoun
  end if
  if (kfl_taufa_nsa(4,2) == 1) then     
     tauxi(      1,4) = 4.0_rp * xvisc / hleti / xdens
     tauxi(ndime+2,4) = 4.0_rp * xdith(igaus) / hleti / xdens / xhecp
  end if


! no me acuerdo para que era esto del sspee... pero asi como esta no es nada

!  tauxi(      1,5) = sspee / hleng(ndime)
!  tauxi(ndime+1,5) = sspee / hleng(ndime)
!  tauxi(ndime+2,5) = sspee / hleng(ndime)

  if (kfl_lopre_nsa > 1) then
!!tauxi(      1,4) = 0.0_rp
!!tauxi(ndime+2,4) = 0.0_rp 
     if (kfl_lopre_nsa == 3) then ! CHOI & MERKLE PRECONDITIONER IS APPLIED
        if (kfl_pseud_nsa == 1 .and. kfl_taufa_nsa(1,2) == 1) then !  PSEUDO TIME IS APPLIED
           cfl_veloc = velmo(igaus) / dtinv / hconv
!           cfl_veloc = sqrt(1.0_rp + 1.0_rp / cfl_veloc / cfl_veloc) / 2.0_rp
!!!           cfl_veloc = sqrt(1.0_rp + 1.0_rp / cfl_veloc / cfl_veloc) / 4.0_rp  ! very diffusive, but ok for the explicit
           cfl_veloc = sqrt(1.0_rp + 1.0_rp / cfl_veloc / cfl_veloc)           ! very low diffusivity, but not enough perhaps for sharp gradients
           tauxi(      1,1) = tauxi(      1,1) * cfl_veloc
           tauxi(ndime+1,1) = tauxi(ndime+1,1) * cfl_veloc
           tauxi(ndime+2,1) = tauxi(ndime+2,1) * cfl_veloc
        end if
     end if
     tauxi(      1,3) = 0.0_rp
     tauxi(ndime+1,3) = 0.0_rp
     tauxi(ndime+2,3) = 0.0_rp
  end if

!  tiint= tauxi(      1,1) + tauxi(      1,2) + tauxi(      1,3) + tauxi(      1,4) + tauxi(      1,5) 
!  tiinc= tauxi(ndime+1,1) + tauxi(ndime+1,2) + tauxi(ndime+1,3)                    + tauxi(ndime+1,5)  
!  tiine= tauxi(ndime+2,1) + tauxi(ndime+2,2) + tauxi(ndime+2,3) + tauxi(ndime+2,4) + tauxi(ndime+2,5)

  tiint= tauxi(      1,1) + tauxi(      1,3) + tauxi(      1,4)
  tiinc= tauxi(ndime+1,1) + tauxi(ndime+1,3)                    
  tiine= tauxi(ndime+2,1) + tauxi(ndime+2,3) + tauxi(ndime+2,4) 


  taush= 0.0_rp
  if (sspee > sound(igaus)) then
     taush= qufac * hleng(ndime) / sspee      
  end if
  

!!  tiinc= tiint                                 ! this works, but it is not logical
!!  tiinc= tauxi(1,1) + tauxi(1,2)               ! this DOES NOT work!!!!!!!!!

  tauen= 0.0_rp
  taupa= 0.0_rp
  tauco= 0.0_rp
  
  if (tiine > zensa) tauen= qufac / tiine
  if (tiint > zensa) taupa= qufac / tiint
  if (tiinc > zensa) tauco= qufac / tiinc

!!$  ! Tau parameter as computed in López&Nigro's paper
!!$  tauen= 0.0_rp
!!$  taupa= 0.0_rp
!!$  tauco= 0.0_rp
!!$  xmach = velmo(igaus) / sound(igaus)
!!$  cfl_sound = sound(igaus) / dtinv / hsoun
!!$  xmach_ref = sqrt(xmach * xmach + 1.0_rp / cfl_sound / cfl_sound)
!!$  a = 1.0_rp + xmach_ref * xmach_ref
!!$  sound_pseud = velmo(igaus) * velmo(igaus) * a * a + &
!!$       4.0_rp * xmach_ref * xmach_ref * (sound(igaus) * sound(igaus) - velmo(igaus) * velmo(igaus))
!!$  sound_pseud = sqrt(sound_pseud)
!!$  b = (velmo(igaus) * a + sound_pseud) / hleng(ndime)
!!$  tauen = b * b + 4.0_rp * dtinv * dtinv
!!$  tauen = 1.0_rp / sqrt(tauen)
!!$  taupa = tauen
!!$  tauco = tauen

  if (kfl_stabi_nsa == 2) then
     call runend('NSA_MONYOS: CG DEPRECATED - USE MULTISCALE STABILIZATION!!')
  end if

  taudi(1:ndime) = taupa + taush
  taudi(ndofc)   = tauco + taush
  taudi(ndofe)   = tauen + taush

  if (kfl_isent_nsa == 1) taudi(ndofe)=0.0_rp

  xsute= 0.0_rp
  xunte= 0.0_rp

  do idofn=1,ndofn_nsa
!     xsute = xsube(idofn,igaus,2) / xdtix(idofn,igaus,2)  ! it does not work with this term     
     xsube(idofn,igaus,1) = taudi(idofn) * ( xsute + xresi(idofn,igaus) - xunte) 
  end do

!if (ielem == 1 .and. igaus == 1) then
!!$print*
!!$print*,'dtinv_nsa',itinn(modul),1.0_rp/dtinv_nsa,1.0_rp/dtinv
!!$print*,'xsube',xsube(ndime+1,igaus,1)
!   write(6,*) 'tololo',dtinv,velmo(igaus),sound(igaus)
!end if


end subroutine nsa_monyos_diag
