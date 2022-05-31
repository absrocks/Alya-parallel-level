subroutine elsest_geogid(ndime,mnode,nelem,kelem,kpoin,nnode,lnods,ltype)
  use def_kintyp
  use def_elsest
  implicit none
  integer(ip), intent(in)  :: ndime,mnode,nelem,kelem,kpoin
  integer(ip), intent(in)  :: nnode(*),lnods(mnode,nelem),ltype(nelem)
  integer(ip)              :: ielem,istat,iesta_dom,iesto_dom
  integer(ip)              :: ielty,inode
  character(20)            :: cetop(50)
  integer(ip), parameter :: BAR02 =   2 ! 1D
  integer(ip), parameter :: BAR03 =   3 ! 1D 
  integer(ip), parameter :: BAR04 =   4 ! 1D 
  integer(ip), parameter :: TRI03 =  10 ! 2D 
  integer(ip), parameter :: TRI06 =  11 ! 2D 
  integer(ip), parameter :: QUA04 =  12 ! 2D 
  integer(ip), parameter :: QUA08 =  13 ! 2D 
  integer(ip), parameter :: QUA09 =  14 ! 2D 
  integer(ip), parameter :: QUA16 =  15 ! 2D 
  integer(ip), parameter :: TET04 =  30 ! 3D 
  integer(ip), parameter :: TET10 =  31 ! 3D 
  integer(ip), parameter :: PYR05 =  32 ! 3D 
  integer(ip), parameter :: PYR14 =  33 ! 3D 
  integer(ip), parameter :: PEN06 =  34 ! 3D 
  integer(ip), parameter :: PEN15 =  35 ! 3D 
  integer(ip), parameter :: PEN18 =  36 ! 3D 
  integer(ip), parameter :: HEX08 =  37 ! 3D 
  !integer(ip), parameter :: HEX20 =  38 ! 3D 
  integer(ip), parameter :: HEX27 =  39 ! 3D 
  integer(ip), parameter :: HEX64 =  40 ! 3D 

  integer(ip), allocatable :: lexis(:)

  cetop(BAR02)    =  'Linear'
  cetop(BAR03)    =  'Linear'
  cetop(TRI03)    =  'Triangle'
  cetop(TRI06)    =  'Triangle'
  cetop(QUA04)   =  'Quadrilateral'
  cetop(QUA08)   =  'Quadrilateral'
  cetop(QUA09)   =  'Quadrilateral'
  cetop(TET04)  =  'Tetrahedra'
  cetop(TET10) =  'Tetrahedra'
  cetop(PYR05)   =  'Pyramid'
  cetop(PYR14)  =  'Pyramid'
  cetop(PEN06) =  'Prism'
  cetop(PEN15) =  'Prism'
  cetop(PEN18) =  'Prism'
  cetop(HEX08)   =  'Hexahedra'
  !cetop(HEX20)  =  'Hexahedra'
  cetop(HEX27)  =  'Hexahedra'

  allocate(lexis(nelem),stat=istat)
  do ielem=1,nelem
     lexis(abs(ltype(ielem)))=1
  end do

  if(ndime==1) then
     iesta_dom=BAR02
     iesto_dom=BAR03
  else if(ndime==2) then
     iesta_dom=TRI03
     iesto_dom=QUA09
  else
     iesta_dom=TET04
     iesto_dom=HEX27
  end if

  do ielty=iesta_dom,iesto_dom
     if(lexis(ielty)/=0) then
        !
        ! Header
        !
        if(ielty<10) then
           write(iunit(2),10)&
                'ELSEST_BACKGROUND',ielty,max(2_ip,ndime),&
                adjustl(trim(cetop(ielty))),nnode(ielty)
        else
           write(iunit(2),11)&
                'ELSEST_BACKGROUND',ielty,max(2_ip,ndime),&
                adjustl(trim(cetop(ielty))),nnode(ielty)
        end if
        !
        ! Connectivity
        !
        write(iunit(2),*) 'elements'
        do ielem=1,nelem
           if(abs(ltype(ielem))==ielty) then
              write(iunit(2),4) ielem+kelem,&
                   (lnods(inode,ielem)+kpoin,inode=1,nnode(ielty)),0_ip
           end if
        end do
        write(iunit(2),*) 'end elements'
     end if
  end do
  deallocate(lexis,stat=istat)

10 format('MESH ',a,i1,' dimension ',i1,' Elemtype ',a,' Nnode ',i2)
11 format('MESH ',a,i2,' dimension ',i1,' Elemtype ',a,' Nnode ',i2)
2 format(a)
3 format(i7, 3(1x,e16.8e3))
4 format(i7,50(1x,i7))

end subroutine elsest_geogid
