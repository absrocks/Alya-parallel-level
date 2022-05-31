subroutine elsest_octpos(&
     ithre,ndime,mnode,npoin,nelem,lnods,ltype,coord,nnode,ipara)
  !
  ! Post-process the Oct/Quad-tree mesh
  !
  use def_elsest
  use mod_elsest
  implicit none
  integer(ip), intent(in)  :: ithre,ndime,mnode,npoin,nelem
  integer(ip), intent(in)  :: lnods(mnode,nelem),ltype(nelem),nnode(*)
  integer(ip), intent(in)  :: ipara(*)
  real(rp),    intent(in)  :: coord(ndime,npoin)
  integer(ip)              :: ipoin,ielem,idime,kpoin

  iunit(1) = ipara( 7)
  iunit(2) = ipara(12)
  iunit(3) = ipara(13)

  if(iunit(2)>0) then
     !
     ! Coordinates
     !
     current(ithre)%o       => tree_root
     ipoin=0
     ielem=0
     if(ndime==2) then
        write(iunit(2),*) 'MESH  ELSEST_QUAD dimension 2 Elemtype Quadrilateral Nnode 4'
     else
        write(iunit(2),*) 'MESH  ELSEST_OCT  dimension 3 Elemtype Hexahedra Nnode 8'
     end if
     write(iunit(2),*) 'coordinates'
     call elsest_octdep(1_ip,ithre,ndime,ipoin,ielem)
     do kpoin=1,npoin
        write(iunit(2),*) ipoin+kpoin,(coord(idime,kpoin),idime=1,ndime)
     end do
     write(iunit(2),*)  'end coordinates'

     current(ithre)%o       => tree_root
     ipoin=0
     ielem=0
     write(iunit(2),*) 'elements'
     call elsest_octdep(2_ip,ithre,ndime,ipoin,ielem)
     write(iunit(2),*) 'end elements'
     call elsest_geogid(ndime,mnode,nelem,ielem,ipoin,nnode,lnods,ltype)
  end if

  if(iunit(3)>0) then

     write(iunit(3),*) 'GiD Post Results File 1.0'
     write(iunit(3),*) 'GaussPoints GP_QUAD4 Elemtype Quadrilateral'
     write(iunit(3),*) 'Number of Gauss Points: 1'
     write(iunit(3),*) 'Natural Coordinates: Internal'
     write(iunit(3),*) 'End GaussPoints'

     current(ithre)%o       => tree_root
     ielem=0
     write(iunit(3),*) 'Result NODE_NUMBER ELSEST 0 Scalar OnGaussPoints GP_QUAD4'
     write(iunit(3),*) 'ComponentNames NODE_NUMBER'
     write(iunit(3),*) 'values'
     call elsest_octdep(3_ip,ithre,ndime,ipoin,ielem)
     write(iunit(3),*) 'end values'

     current(ithre)%o       => tree_root
     ielem=0
     write(iunit(3),*) 'Result ELEMENT_NUMBER ELSEST 0 Scalar OnGaussPoints GP_QUAD4'
     write(iunit(3),*) 'ComponentNames ELEMENT_NUMBER'
     write(iunit(3),*) 'values'
     call elsest_octdep(4_ip,ithre,ndime,ipoin,ielem)
     write(iunit(3),*) 'end values'

     current(ithre)%o       => tree_root
     ielem=0
     write(iunit(3),*) 'Result LEVEL ELSEST 0 Scalar OnGaussPoints GP_QUAD4'
     write(iunit(3),*) 'ComponentNames LEVEL'
     write(iunit(3),*) 'values'
     call elsest_octdep(5_ip,ithre,ndime,ipoin,ielem)
     write(iunit(3),*) 'end values'

  end if

end subroutine elsest_octpos
