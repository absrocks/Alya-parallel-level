subroutine got_elmtss(&
     pnode,porde,elcod,elvel,elvdr,elcdr,eldif,shacg,&
     dercg,weicg,dtcri)
  !-----------------------------------------------------------------------
  !****f* Gotita/got_elmtss
  ! NAME 
  !    got_elmtss
  ! DESCRIPTION
  !    This routine computes the element time step of the momentum
  !    and continuity equations:
  !              1                   1 
  !    dtm = --------- , dtc = -------------- 
  !          2*u               2*u   
  !          --- + sig         --- + |div(u)|  
  !           h                 h
  ! USED BY
  !    got_updtss
  !    got_outvar
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  ndime,mnode
  use def_gotita, only       :  deair_got,muair_got,ddrop_got,kfact_got,&
       &                        kfl_ellen_got,kfl_probl_got
  implicit none
  integer(ip), intent(in)    :: pnode,porde
  real(rp),    intent(in)    :: elcod(ndime,pnode),elvel(ndime,pnode)
  real(rp),    intent(in)    :: elvdr(ndime,pnode),elcdr(pnode)
  real(rp),    intent(in)    :: eldif(pnode)
  real(rp),    intent(in)    :: shacg(pnode),dercg(ndime,pnode),weicg
  real(rp),    intent(inout) :: dtcri
  integer(ip)                :: idime,inode
  real(rp)                   :: gpcar(ndime,mnode),xjaci(3),xjacm(3) 
  real(rp)                   :: gpvdr(3),gpvel(3),dtcrm,dtcrc,dummr(2)
  real(rp)                   :: rdinv,velno,gpdet,gplen,gppor,gpdiv
  real(rp)                   :: chale(2),gpdif,gpvno,gplev,gpcdr
  !
  ! Element volume
  !
  dtcri = 0.0_rp
  rdinv = 1.0_rp/real(ndime)
  call elmder(& 
       pnode,ndime,dercg,elcod,gpcar,gpdet,xjacm,xjaci)
  gplev=(weicg*gpdet)**rdinv
  !
  ! Element length
  !
  if(ndime==2) then
     call elmchl(&
          dummr,dummr,elcod,elvdr,dummr,chale,pnode,&
          porde,1.0_rp,1.0_ip,5_ip)
     gplen=chale(1)
  else
     gplen=gplev
  end if
  !
  ! Droplet and air velocity
  !
  gpdiv=0.0_rp
  gpcdr=0.0_rp
  do idime=1,ndime
     gpvdr(idime)=0.0_rp
     gpvel(idime)=0.0_rp
  end do
  do inode=1,pnode
     do idime=1,ndime
        gpvdr(idime) = gpvdr(idime)+shacg(inode)*elvdr(idime,inode)
        gpvel(idime) = gpvel(idime)+shacg(inode)*elvel(idime,inode)
        gpcdr        = gpcdr       +shacg(inode)*elcdr(inode)
        gpdiv        = gpdiv +gpcar(idime,inode)*elvdr(idime,inode)
     end do
  end do
  call vecnor(&
       gpvdr,ndime,gpvno,2_ip)
  call got_elmpro(&
       pnode,1_ip,1_ip,1_ip,eldif,shacg,gpvdr,gpvel,&
       gpvno,gpcdr,chale,gppor,gpdif)
  call vecnor(&
       gpvdr,ndime,velno,2_ip)
  !
  ! Momentum:   DTRCM=1/(4k/(alpha*h^2) + 2u/h + sig)
  !
  if(kfl_probl_got/=3) then
     !dtcrm=gpcdr/(4.0_rp*gpdif/(gplen*gplen)+gpcdr*(2.0_rp*velno/gplen+gppor))
     dtcrm=1.0_rp/(2.0_rp*velno/gplen+gppor)
  else
     dtcrm=1.0e16
  end if
  !
  ! Continuity: DTCRC=1/(2u/h + |div(u)|)
  !
  if(kfl_probl_got/=2) then
     dtcrc=1.0_rp/(2.0_rp*velno/gplen+abs(gpdiv))
  else
     dtcrc=1.0e16
  end if

  dtcri=min(dtcrm,dtcrc)

end subroutine got_elmtss
