subroutine nsi_elmsol(&
     pnode,pgaus,lnods,gpcar,gpvol,elmal)
  !----------------------------------------------------------------------
  !****f* Nastin/nsi_elmsol
  ! NAME 
  !    nsi_elmsol
  ! DESCRIPTION
  !    Compute the Laplacian matrix
  ! USES
  ! USED BY
  !    nsi_elmope
  !***
  !----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp 
  use def_nastin, only     :  nodpr_nsi,lperp_nsi,kfl_perip_nsi
  use def_domain, only     :  mnode,nperi,nbopo,ndime,coord,lpoty
  implicit none
  integer(ip), intent(in)  :: pnode,pgaus
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(in)  :: gpcar(ndime,mnode,pgaus),gpvol(pgaus)
  real(rp),    intent(out) :: elmal(pnode,pnode)
  integer(ip)              :: inode,jnode,kdime,igaus,ipoin,ibopo
  integer(ip)              :: ipres,iperi
  real(rp)                 :: fact1

  do inode=1,pnode
     ipoin=lnods(inode)
     !if(lmatn(ipoin)==-1) then
        elmal(inode,jnode)=0.0_rp
     !end if
  end do

end subroutine nsi_elmsol
