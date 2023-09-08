subroutine rad_elmchl(&
     ielem,pelty,pnode,plapl,pmate,lnods,elcod,eltem,eledd,&
     elvel,gpcar,gphes,chale)
  !------------------------------------------------------------------------
  !****f* Radiat/rad_elmchl
  ! NAME 
  !    rad_elmchl
  ! DESCRIPTION
  ! USES
  ! USED BY
  !    rad_matrix
  !------------------------------------------------------------------------
  use def_kintyp, only     :  ip,rp
  use def_domain, only     :  ndime,mnode,ntens,elmar
  implicit none
  integer(ip), intent(in)  :: ielem,pelty,pnode,plapl,pmate
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(in)  :: elcod(ndime,pnode),eltem(pnode)
  real(rp),    intent(in)  :: eledd(pnode),elvel(ndime,pnode)
  real(rp),    intent(in)  :: gpcar(ndime,pnode),gphes(ntens,pnode)
  real(rp),    intent(out) :: chale(2)
  integer(ip)              :: inode,knode,jnode,idime,ipoin,klapl
  real(rp)                 :: gpdir(3),gpcod(3),gpgrt(3),rnode,gpgr2(3)
  real(rp)                 :: gpvol,gptem,gpden,gpcon,gpdif
  real(rp)                 :: gpsph,hleng(3),gprea
  real(rp)                 :: vx1,vy1,vx2,vy2,m,x1,y1,x2,y2
  real(rp)                 :: xmin1,xmin2,ymin1,ymin2,det,dinor

  call velchl(pnode,elcod,elvel,chale,hleng) 
  return

end subroutine rad_elmchl
