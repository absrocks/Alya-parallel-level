subroutine qua_elmtss(&
     ielem,pelty,pnode,pmate,elcod,elvel,eledd,&
     eltem,gpcar,chale,dtcri)
  !-----------------------------------------------------------------------
  !****f* Quanty/qua_elmtss
  ! NAME 
  !    qua_elmtss
  ! DESCRIPTION
  !    This routine computes the element time step
  ! USED BY
  !    qua_updtss
  !    qua_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  use def_domain, only       :  mnode,ndime,elmar
!  use def_quanty, only       :  kfl_cotur_qua,kfl_advec_qua,&
!       &                        lawde_qua,lawtc_qua,lawsp_qua,&
!       &                        densi_qua,tcond_qua,sphea_qua,&
!       &                        kfl_sgsno_qua,kfl_sgsti_qua,&
!       &                        kfl_taust_qua,staco_qua,kfl_regim_qua
  implicit none
  integer(ip), intent(in)    :: ielem,pelty,pnode,pmate
  real(rp),    intent(in)    :: elcod(ndime,pnode)
  real(rp),    intent(in)    :: elvel(ndime,pnode),eledd(pnode)
  real(rp),    intent(in)    :: eltem(pnode),chale(2)
  real(rp),    intent(inout) :: dtcri
  real(rp),    intent(out)   :: gpcar(ndime,mnode)
  real(rp)                   :: xjaci(ndime,ndime) 
  real(rp)                   :: xjacm(ndime,ndime)
  integer(ip)                :: idime,inode
  real(rp)                   :: gpcon,gpden,gpsph,gprcp,gprea,gpdif
  real(rp)                   :: gptem,rnode,gpvno,adv,dif,rea
  real(rp)                   :: gpvol,gphes,gpvel(3)
  !
  ! Initialization
  !
  gptem = 0.0_rp
  gpvno = 0.0_rp
  gpvel = 0.0_rp
  dtcri = 0.0_rp
  rnode = 1.0_rp/real(pnode)
  !
  ! GPTEM: mean value solution 
  !
  do inode=1,pnode
     gptem=gptem+eltem(inode)
  end do
  gptem=rnode*gptem
  !
  ! GPDEN, GPSPH, GPDIF, GPREA: Properties
  !
!  call qua_elmpro(& 
!       2_ip,ielem,pmate,pnode,1_ip,1_ip,1_ip,eledd,elvel,&
!       gptem,elmar(pelty)%shacg,gpcar,gpden,gpcon,gpsph,&
!       gpdif,gprea)          
  
  gprcp=gpden*gpsph
  !
  ! DTCRI: critical time step
  !
!  adv = gprcp*gpvno
!  dif = gpdif
!  rea = gprea
!  call tauadr(kfl_taust_qua,staco_qua,adv,dif,rea,chale(1),chale(2),dtcri)
!  dtcri = gprcp*dtcri

end subroutine qua_elmtss
