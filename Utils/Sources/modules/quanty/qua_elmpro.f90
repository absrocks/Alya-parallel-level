subroutine qua_elmpro(&
     itask,ielem,pmate,pnode,pgaus,igaui,igauf,eledd,&
     elvel,gpqua,gpsha,gpcar,gpden,gpcon,gpsph,gpdif,&
     gprea)   
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_elmpro
  ! NAME 
  !    qua_elmpro
  ! DESCRIPTION
  !    Compute properties:
  !    1. GPDEN: Density rho
  !    2. GPSPH: Specific heat Cp
  !    3. GPDIF: Total thermal conductivity k+kt
  !    4. GPRGD: Thermal conductivity gradient grad(k+kt)
  ! USES
  ! USED BY
  !    qua_elmope
  !    qua_bouset
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_domain
  use def_master
  use def_quanty

  implicit none
  integer(ip), intent(in)  :: itask,ielem,pmate,pnode,pgaus,igaui,igauf
  real(rp),    intent(in)  :: eledd(pnode),elvel(ndime,pnode)
  real(rp),    intent(in)  :: gpqua(pgaus)
  real(rp),    intent(in)  :: gpsha(pnode,pgaus),gpcar(ndime,mnode,pgaus)
  real(rp),    intent(out) :: gpcon(pgaus),gpden(pgaus),gpsph(pgaus)
  real(rp),    intent(out) :: gpdif(pgaus),gprea(pgaus)
  integer(ip)              :: igaus,idime,inode,ipoin
  real(rp)                 :: gpmut,gpdiv,dummr
  !
  ! coeficiente GPDIF=GPCON=k
  !
!  if(lawtc_qua(pmate)==0) then

     do igaus=igaui,igauf
        gpcon(igaus)=-0.5
     end do

!  else if(lawtc_qua(pmate)==10) then
!     do igaus=igaui,igauf
!        call linint(&
!             nkint_qua(pmate),gptem(igaus),coefk_qua(1,1,pmate),&
!             gpcon(igaus),dummr)
!     end do
!  end if
  do igaus=igaui,igauf
     gpdif(igaus)=gpcon(igaus)
  end do

  if(itask/=3) then
     !
     ! Specific heat GPSPH=cp
     !
!     if(lawsp_tem(pmate)==0) then
!        do igaus=igaui,igauf
!           gpsph(igaus)=sphea_tem(1,pmate)
!        end do
!     else if(lawsp_tem(pmate)==10) then
!        do igaus=igaui,igauf
!           call linint(&
!                ncint_tem(pmate),gptem(igaus),coefc_tem(1,1,pmate),&
!                gpsph(igaus),dummr)
!        end do
!     end if
     !
     ! Total thermal diffusion GPDIF=k+kt; kt=mut*Cp/Prt
     !
!     if(kfl_cotur_tem/=0) then
!        do igaus=igaui,igauf        
!           gpmut=0.0_rp
!           do inode=1,pnode
!              gpmut=gpmut+gpsha(inode,igaus)*eledd(inode)
!           end do
!           gpdif(igaus)=gpdif(igaus)+gpsph(igaus)*gpmut/prtur_tem
!        end do
!     end if
  end if

  if(itask==2) then
     !
     ! Density GPDEN=rho
     !
!     if(kfl_regim_qua==1.or.kfl_regim_qua==2) then
        !
        ! Compressible regime
        !
        !do igaus=igaui,igauf
        !   gpden(igaus)=0.0_rp
        !end do
        !do inode=1,pnode
        !   ipoin=lnods(inode,ielem) 
        !   do igaus=igaui,igauf
        !      gpden(igaus)=gpden(igaus)+gpsha(inode,igaus)*densi(ipoin,1)
        !   end do
        !end do

!     else
        !
        ! Density law
        !
!        if(lawde_tem(pmate)==0) then
!           do igaus=igaui,igauf
!              gpden(igaus)=densi_tem(1,pmate)
!           end do
!        end if

!     end if
     !
     ! Reaction: GPREA=r
     !
!     do igaus=igaui,igauf
!        gprea(igaus)=react_tem
!     end do
     !
     ! Compressibility term: GPREA=r + rho*R*div(u)
     !
!     if(kfl_regim_tem==1.or.kfl_regim_tem==2) then
!        do igaus=igaui,igauf
!           gpdiv=0.0_rp
!           do inode=1,pnode
!              do idime=1,ndime
!                 gpdiv=gpdiv&
!                      +elvel(idime,inode)*gpcar(idime,inode,igaus)
!              end do
!           end do
!           gprea(igaus)=gprea(igaus)+gpden(igaus)*gasco_tem*gpdiv
!        end do
!     end if
  end if

end subroutine qua_elmpro
