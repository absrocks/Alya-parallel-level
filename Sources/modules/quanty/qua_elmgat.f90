subroutine qua_elmgat(&
     pnode,lnods,elpro,elcod)
  !------------------------------------------------------------------------
  !****f* Quanty/qua_elmgat
  ! NAME 
  !    qua_elmgat
  ! DESCRIPTION
  !    This routine performs the gather operations
  ! USES
  ! USED BY
  !    qua_elmope
  !***
  !------------------------------------------------------------------------ 
  use def_kintyp, only     :  ip,rp 
  use def_domain, only     :  ndime,mnode,npoin,coord
  use def_master, only     :  phion,rhoon
  use def_quanty, only     :  kfl_timei_qua,&
                              kfl_tiacc_qua,kfl_tisch_qua
  implicit none
  integer(ip), intent(in)  :: pnode
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(out) :: elpro(mnode)
  real(rp),    intent(out) :: elcod(ndime,mnode)
  integer(ip)              :: inode,ipoin,idime,itime
  !
  ! Current function and coordinates
  !
  do inode=1,pnode
     ipoin=lnods(inode)
     elpro(inode)=rhoon(ipoin,1)
     do idime=1,ndime
        elcod(idime,inode)=coord(idime,ipoin)
     end do  
  end do
  
end subroutine qua_elmgat
