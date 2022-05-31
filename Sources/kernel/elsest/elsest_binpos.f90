subroutine elsest_binpos(&
     ithre,ndime,mnode,npoin,nelem,lnods,ltype,coord,nnode,ipara)
  !
  ! Post-process the bin mesh
  !
  use def_elsest
  use mod_elsest
  implicit none
  integer(ip), intent(in)  :: ithre,ndime,mnode,npoin,nelem
  integer(ip), intent(in)  :: lnods(mnode,nelem),ltype(nelem),nnode(*)
  integer(ip), intent(in)  :: ipara(*)
  real(rp),    intent(in)  :: coord(ndime,npoin)
  integer(ip)              :: ipoin,ielem,idime,kpoin
  integer(ip)              :: box_coord(3),ii,jj,iboxe
 
  iunit(1) = ipara( 7)
  iunit(2) = ipara(12)
  iunit(3) = ipara(13)
  
  if(iunit(2)>0) then
     !
     ! Coordinates
     !
     ipoin=0
     ielem=0
     if(ndime==2) then
        write(iunit(2),*) 'MESH ELSEST_BIN dimension 2 Elemtype Quadrilateral Nnode 4'
        write(iunit(2),*) 'coordinates'
        do iboxe=1,nboxe
           call elsest_boxcoo(ndime,iboxe,box_coord)
           ii = box_coord(1)
           jj = box_coord(2)
           ipoin=ipoin+1
           write(iunit(2),*) ipoin,comin(1)+real(ii-1,rp)*delta(1),comin(2)+real(jj-1,rp)*delta(2)
           ipoin=ipoin+1
           write(iunit(2),*) ipoin,comin(1)+real(ii,rp)*delta(1),  comin(2)+real(jj-1,rp)*delta(2)
           ipoin=ipoin+1
           write(iunit(2),*) ipoin,comin(1)+real(ii,rp)*delta(1),  comin(2)+real(jj,rp)*delta(2)
           ipoin=ipoin+1
           write(iunit(2),*) ipoin,comin(1)+real(ii-1,rp)*delta(1),comin(2)+real(jj,rp)*delta(2)
        end do
     end if
     do kpoin=1,npoin
        write(iunit(2),*) ipoin+kpoin,(coord(idime,kpoin),idime=1,ndime)
     end do
     write(iunit(2),*)  'end coordinates'
     !
     ! Elements
     !
     ipoin=1
     ielem=0
     write(iunit(2),*) 'elements'
     if(ndime==2) then
        do iboxe=1,nboxe
           write(iunit(2),*) iboxe,ipoin,ipoin+1,ipoin+2,ipoin+3,1
           ipoin=ipoin+4
        end do
     end if
     write(iunit(2),*) 'end elements'
     ipoin=nboxe*4
     call elsest_geogid(ndime,mnode,nelem,nboxe,ipoin,nnode,lnods,ltype)
  end if

  if(iunit(3)>0) then
     write(iunit(3),*) 'GiD Post Results File 1.0'
     write(iunit(3),*) ' '
     write(iunit(3),*) 'GaussPoints GP_QUAD4 Elemtype Quadrilateral'
     write(iunit(3),*) 'Number of Gauss Points: 1'
     write(iunit(3),*) 'Natural Coordinates: Internal'
     write(iunit(3),*) 'End GaussPoints'

     write(iunit(3),*) 'Result ELEMENT_NUMBER ELSEST 0 Scalar OnGaussPoints GP_QUAD4'
     write(iunit(3),*) 'ComponentNames ELEMENT_NUMBER'
     write(iunit(3),*) 'values'
     if(dataf==0) then
        do iboxe=1,nboxe
           if( associated(tboel(iboxe)%l) ) then
              write(iunit(3),*) iboxe,size(tboel(iboxe)%l)
           else
              write(iunit(3),*) iboxe,0
           end if
        end do
     else
        do iboxe=1,nboxe
           write(iunit(3),*) iboxe,pboel(iboxe+1)-pboel(iboxe)
        end do
     end if
     write(iunit(3),*) 'end values'
  end if

end subroutine elsest_binpos
