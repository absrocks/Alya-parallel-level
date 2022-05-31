subroutine elsest_inelem(&
     mnode,ndime,nnode,lnods,ltype,ltopo,coord,xcoor,rpara,&
     ifoun,shapt,derit,coloc)
  !-----------------------------------------------------------------------
  !****f* Elsest
  ! NAME
  !    Inelem
  ! DESCRIPTION
  !    Identify shape function and derivatives in element
  ! INPUT
  ! OUTPUT
  ! USES
  ! USED BY
  !    Alya
  !***
  !-----------------------------------------------------------------------
  use def_elsest, only     :  ip,rp,iunit,nthre,kfl_memor,lmini,lmaxi
  implicit none
  integer(ip), intent(in)  :: mnode,ndime
  integer(ip), intent(in)  :: nnode(*)
  integer(ip), intent(in)  :: lnods(mnode,*),ltype(*)
  integer(ip), intent(in)  :: ltopo(*)
  real(rp),    intent(in)  :: coord(ndime,*),xcoor(*)
  real(rp),    intent(in)  :: rpara(*)
  integer(ip), intent(out) :: ifoun
  real(rp),    intent(out) :: shapt(*),derit(*),coloc(*)
  integer(ip)              :: pelty,ptopo,pnode,inode,idime,ipoin
  real(rp)                 :: elcod(ndime,mnode)

  lmini = -rpara(1)
  lmaxi = 1.0_rp + rpara(1)
  pelty = ltype(ifoun)
  ptopo = ltopo(pelty)
  pnode = nnode(pelty)
  do inode = 1,pnode
     ipoin = lnods(inode,ifoun)
     do idime = 1,ndime
        elcod(idime,inode) = coord(idime,ipoin)
     end do
  end do
  call elsest_chkelm(&
       ndime,ptopo,pnode,elcod,shapt,derit,&
       xcoor,coloc,ifoun,lmini,lmaxi)

end subroutine elsest_inelem
