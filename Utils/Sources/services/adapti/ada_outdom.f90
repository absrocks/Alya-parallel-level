subroutine ada_outdom
!-----------------------------------------------------------------------
!****f* adapti/ada_outdom
! NAME 
!    ada_outdom
! DESCRIPTION
!    This routine writes the dom.dat file
! USES
!
! USED BY
!
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_adapti
  implicit none
  integer(ip)   :: &
       ielem,ipoin,idime,inode
    
  write(lun_oudom_ada,*) '$--------------------------------------------------'
  write(lun_oudom_ada,*) 'DIMENSIONS'
  write(lun_oudom_ada,*) '  NODAL_POINTS=         '  , npoin 
  write(lun_oudom_ada,*) '  ELEMENTS=             '  , nelem 
  write(lun_oudom_ada,*) '  SPACE_DIMENSIONS=     '  , ndime
  write(lun_oudom_ada,*) '  TYPES_OF_ELEMS=       '  , peltb_ada 
  write(lun_oudom_ada,*) '  BOUNDARIES=           0 ' 
  write(lun_oudom_ada,*) '  SKEW_SYSTEMS=      0  '
  write(lun_oudom_ada,*) '  SLAVES=            0  '
  write(lun_oudom_ada,*) 'END_DIMENSIONS'
  write(lun_oudom_ada,*) '$--------------------------------------------------'
  write(lun_oudom_ada,*) 'STRATEGY'
  write(lun_oudom_ada,*) '  DOMAIN_INTEGRATION_POINTS:    1  ' 
  write(lun_oudom_ada,*) '  WRITE_EXTERIOR_NORMAL:        NO '
  write(lun_oudom_ada,*) '  OUTPUT_MESH_DATA:             YES'
  write(lun_oudom_ada,*) 'END_STRATEGY'
  write(lun_oudom_ada,*) '$--------------------------------------------------'
  write(lun_oudom_ada,*) 'GEOMETRY, GID'
  write(lun_oudom_ada,*) 'TYPES'
  do ielem=1,nelem
     write(lun_oudom_ada,*) ielem,ltype(ielem)
  end do
  write(lun_oudom_ada,*) 'END_TYPES'
  write(lun_oudom_ada,*) 'ELEMENTS'
  do ielem=1,nelem
     write(lun_oudom_ada,300) ielem,(lnods(inode,ielem),inode=1,mnode)
  end do
  write(lun_oudom_ada,*) 'END_ELEMENTS'
  
300 format(10(i12,2x))
  
  write(lun_oudom_ada,*) 'COORDINATES'
  do ipoin=1,npoin
     write(lun_oudom_ada,400) ipoin,(coord(idime,ipoin),idime=1,ndime)     
  end do
  write(lun_oudom_ada,*) 'END_COORDINATES'
  
400 format(i12,3(2x,e25.15))  
  
  write(lun_oudom_ada,*) 'END_GEOMETRY'
  write(lun_oudom_ada,*) 'SETS'
  write(lun_oudom_ada,*) 'END_SETS'


end subroutine ada_outdom
