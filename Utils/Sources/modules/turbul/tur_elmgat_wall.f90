subroutine tur_elmgat_wall(&
     pnode,lnods,elwao,elwai_coo)

  !-----------------------------------------------------------------------
  !****f* Turbul/tur_elmgat_wall
  ! NAME
  !   tur_elmgat_wall
  ! DESCRIPTION
  !    Gather operations for a generic turbulence model on the wall
  ! USES
  ! USED BY
  !    tur_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime,ntens,mnode,coord,wallcoor,wallo
  implicit none
  integer(ip), intent(in)  :: pnode
  integer(ip), intent(in)  :: lnods(pnode)
  integer(ip), intent(out) :: elwao(pnode)
  real(rp),    intent(out) :: elwai_coo(ndime,pnode)
  
  integer(ip)              :: inode,ipoin,iturb,idime,itime
  !
  ! Distance to the wall
  !
  do inode = 1,pnode
     ipoin = lnods(inode)
     elwao(inode) = wallo(ipoin)
     do idime = 1,ndime
       elwai_coo(idime,inode) = wallcoor(idime,ipoin)
     enddo
  end do

end subroutine tur_elmgat_wall
