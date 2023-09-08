subroutine rad_elmgat(pnode,lnods,elrad,elcod)
  !------------------------------------------------------------------------
  !****f* Radiat/rad_elmgat
  ! NAME 
  !    rad_elmgat
  ! DESCRIPTION
  !    This routine performs the gather operations
  ! USES
  ! USED BY
  !    rad_elmope
  !    rad_bouope
  !***
  !------------------------------------------------------------------------ 
  use def_kintyp, only     :  ip,rp 
  use def_domain, only     :  ndime,mnode,npoin,coord
  use def_radiat, only      :  radav_rad

  implicit none
  integer(ip), intent(in)  :: pnode
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(out) :: elrad(pnode,*)
  real(rp),    intent(out) :: elcod(ndime,mnode)
  integer(ip)              :: inode,ipoin,idime,itime

  !
  ! Just coordinates
  !
  do inode=1,pnode
     ipoin=lnods(inode)
     elrad(inode,1)=radav_rad(ipoin,1)
     do idime=1,ndime
        elcod(idime,inode)=coord(idime,ipoin)
     end do  
  end do


end subroutine rad_elmgat
