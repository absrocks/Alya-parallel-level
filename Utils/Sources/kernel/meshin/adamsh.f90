subroutine adamsh(itask)
  !-----------------------------------------------------------------------
  !****f* meshin/adamsh
  ! NAME
  !    adamsh
  ! DESCRIPTION
  !
  !    LMARK =  1: refine
  !          =  0: do nothing
  !          = -1: coarsen
  !    NLAY ................... Number of layer to be refined
  !    NTREE .................. 12
  !    NFTREE ................. 8
  !    HTREE(NTREE,NELEM) ..... Integer tree (=0 at beginning)
  !    HFTREE(NFTREE,NBOUN) ... Integer tree (=0 at beginning)
  !    NLAYMAX ................ Max. number of refinements
  !    GEVEC(NVAR,NPOIN) ...... Unknowns to be interpolated
  ! OUTPUT
  ! USED BY
  !    domain
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_elmtyp
  use def_master
  use def_kermod
  use def_domain
  use mod_memchk
  use mod_href
  use mod_messages, only : livinf
  implicit none
  integer(ip), intent(in) :: itask 
  integer(ip)             :: ntree,nftree,nlay,nlaymax,nvar
  integer(ip)             :: ielem,iboun,ielty
  integer(4)              :: istat
  integer(ip), pointer    :: htree(:,:)
  integer(ip), pointer    :: hftree(:,:)
  integer(ip), pointer    :: lmark(:)
  type(i1p),   pointer    :: lboib(:),lpoel(:)
  integer(ip), pointer    :: lnoel(:)
  integer(ip)             :: ipoin,inodb,pnode,pnodb,jelem,kelem
  integer(ip)             :: jpoin,icoin,kpoin,inode,icrit
  integer(ip)             :: ipoi1,ipoi2,ipoi3,ipoi4
  real(rp)                :: dista(6),dimin

  return

  select case(itask)

  case(1_ip)

     if( ISEQUEN .or. ( IMASTER .and. kfl_ptask /= 2 ) ) then

        !-------------------------------------------------------------------
        !
        ! Adapt mesh
        !
        !-------------------------------------------------------------------
        !
        ! Check errors
        !
        call livinf(73_ip,' ',zero)
        do ielty = iesta_dom,iesto_dom
           if( lexis(ielty) == 1 .and. ielty /= TET04 ) &
                call runend('ADAMSH: CAN ONLY ADAPT TET04 ELEMENTS')
        end do
        !
        ! Allocate memory
        !
        nlay      = 1
        nlaymax   = 10 
        ntree     = 12
        nftree    = 8
        nvar      = 0
        kfl_bouel = 0
        allocate(lmark(nelem),stat=istat) 
        call memchk(zero,istat,memor_dom,'LMARK','adamsh',lmark)
        allocate(htree(ntree,nelem),stat=istat) 
        call memchk(zero,istat,memor_dom,'HTREE','adamsh',htree)
        allocate(hftree(nftree,nboun),stat=istat) 
        call memchk(zero,istat,memor_dom,'HFTREE','adamsh',hftree)
        call memgen(0_ip,1_ip,1_ip) ! GEVEC
        !
        ! Mark elements
        !
        call livinf(89_ip,' ',zero)
        icrit = 2
        if( icrit == 1 ) then
           !
           ! Refine all the mesh
           !
           do ielem = 1,nelem
              lmark(ielem) = 1
           end do
        else if( icrit == 2 ) then
           do ielem = 1,nelem
              lmark(ielem) = 0
              do inode = 1,nnode(ltype(ielem))
                 ipoin = lnods(inode,ielem)
                 if( abs(coord(3,ipoin)) < 0.2_rp ) lmark(ielem) = 1
              end do
              if( lmark(ielem) == 1 ) then
                 ipoi1 = lnods(1,ielem)
                 ipoi2 = lnods(2,ielem)
                 ipoi3 = lnods(3,ielem)
                 ipoi4 = lnods(4,ielem)
                 dista(1) =   ( coord(1,ipoi1) - coord(1,ipoi2) ) ** 2 &
                      &     + ( coord(2,ipoi1) - coord(2,ipoi2) ) ** 2 &
                      &     + ( coord(3,ipoi1) - coord(3,ipoi2) ) ** 2 
                 dista(2) =   ( coord(1,ipoi1) - coord(1,ipoi3) ) ** 2 &
                      &     + ( coord(2,ipoi1) - coord(2,ipoi3) ) ** 2 &
                      &     + ( coord(3,ipoi1) - coord(3,ipoi3) ) ** 2 
                 dista(3) =   ( coord(1,ipoi1) - coord(1,ipoi4) ) ** 2 &
                      &     + ( coord(2,ipoi1) - coord(2,ipoi4) ) ** 2 &
                      &     + ( coord(3,ipoi1) - coord(3,ipoi4) ) ** 2 
                 dista(4) =   ( coord(1,ipoi2) - coord(1,ipoi3) ) ** 2 &
                      &     + ( coord(2,ipoi2) - coord(2,ipoi3) ) ** 2 &
                      &     + ( coord(3,ipoi2) - coord(3,ipoi3) ) ** 2 
                 dista(5) =   ( coord(1,ipoi2) - coord(1,ipoi4) ) ** 2 &
                      &     + ( coord(2,ipoi2) - coord(2,ipoi4) ) ** 2 &
                      &     + ( coord(3,ipoi2) - coord(3,ipoi4) ) ** 2 
                 dista(6) =   ( coord(1,ipoi3) - coord(1,ipoi4) ) ** 2 &
                      &     + ( coord(2,ipoi3) - coord(2,ipoi4) ) ** 2 &
                      &     + ( coord(3,ipoi3) - coord(3,ipoi4) ) ** 2 
                 dimin = minval( dista )
                 dimin = sqrt(dimin)
                 if( dimin < 0.15_rp ) lmark(ielem) = 0
              end if
           end do
        end if
        !
        ! Divide mesh and deallocate memory
        !
        call livinf(83_ip,' ',zero)
        call href(&
             4_ip,3_ip,6_ip,nelem,npoin,nboun,lmark,nlay,ndime,coord,&
             htree,lnods,gevec,ntree,nvar,lnodb,kfl_codbo,hftree,nftree,&
             nlaymax)
        call memgen(2_ip,1,1)
        call memchk(two,istat,memor_dom,'HFTREE','adamsh',hftree)
        deallocate(hftree,stat=istat) 
        if(istat/=0) call memerr(two,'HFTREE','adamsh',0_ip)
        call memchk(two,istat,memor_dom,'HTREE','adamsh',htree)
        deallocate(htree,stat=istat) 
        if(istat/=0) call memerr(two,'HTREE','adamsh',0_ip)
        call memchk(two,istat,memor_dom,'LMARK','adamsh',lmark)
        deallocate(lmark,stat=istat) 
        if(istat/=0) call memerr(two,'LMARK','adamsh',0_ip)
        call livinf(86_ip,' ',zero)
        !
        ! Reallocate memory
        !
        call memgeo(-1_ip)
        call memgeo( 1_ip)
        !
        ! Recompute arrays
        !
        kfl_bouel = 1
        do ielem = 1,nelem
           ltype(ielem) = TET04
        end do
        do iboun = 1,nboun
           ltypb(iboun) = TRI03
        end do
        !
        ! LBOEL
        !
        call livinf(84_ip,' ',zero)
        allocate(lboib(nboun),stat=istat)
        allocate(lpoel(npoin),stat=istat)
        allocate(lnoel(npoin),stat=istat)

        do ipoin = 1,npoin
           lnoel(ipoin) = 0
        end do
        do ielem = 1,nelem
           pnode = nnode(ltype(ielem))
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              lnoel(ipoin) = lnoel(ipoin) + 1
           end do
        end do
        do ipoin = 1,npoin
           if( lnoel(ipoin) > 0 ) then
              kelem = lnoel(ipoin)
              allocate( lpoel(ipoin)%l(kelem), stat=istat )
              do ielem = 1,kelem
                 lpoel(ipoin)%l(ielem) = 0
              end do
              lnoel(ipoin) = 0
           end if
        end do
        do ielem = 1,nelem
           pnode = nnode(ltype(ielem))
           do inode = 1,pnode
              ipoin = lnods(inode,ielem)
              lnoel(ipoin) = lnoel(ipoin) + 1
              lpoel(ipoin)%l(lnoel(ipoin)) = ielem
           end do
        end do
        do iboun = 1,nboun
           kpoin = 0
           do inodb = 1,nnode(ltypb(iboun))
              ipoin = lnodb(inodb,iboun)
              kpoin = kpoin + lnoel(ipoin)
           end do
           allocate( lboib(iboun)%l(kpoin), stat=istat )
           kpoin = 0
           do inodb = 1,nnode(ltypb(iboun))
              ipoin = lnodb(inodb,iboun)
              do kelem = 1,lnoel(ipoin)
                 kpoin = kpoin + 1
                 lboib(iboun)%l(kpoin) = lpoel(ipoin)%l(kelem)
              end do
           end do

        end do

        do iboun = 1,nboun
           pnodb = nnode(ltypb(iboun))
           kelem = 0
           do while( kelem < size(lboib(iboun)%l) )
              kelem = kelem + 1
              icoin = 0
              ielem = lboib(iboun)%l(kelem)    
              if( ielem > 0 ) then
                 pnode = nnode(ltype(ielem))              
                 do inodb = 1,pnodb
                    ipoin = lnodb(inodb,iboun)
                    inode = 0
                    do while( inode < pnode )
                       inode = inode + 1
                       jpoin = lnods(inode,ielem)
                       if( ipoin == jpoin ) then
                          icoin = icoin + 1
                          inode = pnode
                       end if
                    end do
                 end do
                 if( icoin == pnodb ) then
                    lelbo(iboun) = ielem
                    kelem = size(lboib(iboun)%l)
                 else
                    do jelem = ielem,size(lboib(iboun)%l)
                       if( ielem == lboib(iboun)%l(jelem) ) then
                          lboib(iboun)%l(jelem) = 0
                       end if
                    end do
                 end if
              end if
           end do

        end do

        do iboun = 1,nboun
           deallocate( lboib(iboun)%l, stat=istat )
        end do
        do ipoin = 1,npoin
           deallocate( lpoel(ipoin)%l, stat=istat )
        end do
        deallocate(lboib,stat=istat)
        deallocate(lpoel,stat=istat)
        deallocate(lnoel,stat=istat)
        !
        ! Output mesh
        !
        call outnew()

        call livinf(85_ip,' ',zero)

     end if

  end select

end subroutine adamsh

subroutine outnew()
  use def_kintyp
  use def_elmtyp
  use def_master
  use def_kermod
  use def_domain
  use mod_messages, only : livinf
  implicit none
  integer(ip)   :: ielem,iboun
  integer(ip)   :: ipoin,inodb,pnode,pnodb
  integer(ip)   :: inode,idime,icond,ipara
  character(12) :: wcond(10)

  call livinf(88_ip,' ',0_ip)
        
  open(unit=80,file=adjustl(trim(namda))//'-new.geo.dat',status='unknown')
  open(unit=81,file=adjustl(trim(namda))//'-new.fix.dat',status='unknown')
  open(unit=82,file=adjustl(trim(namda))//'-new.dom.dat',status='unknown')
  !
  ! Geometry file
  !
  write(80,'(a)') 'TYPES, ALL=TET04'
  write(80,'(a)') 'END_TYPES'
  write(80,'(a)') 'ELEMENTS'
  do ielem = 1,nelem
     pnode = nnode(ltype(ielem))
     write(80,'(5(1x,i9))') ielem,(lnods(inode,ielem),inode=1,pnode)        
  end do
  write(80,'(a)') 'END_ELEMENTS'
  write(80,'(a)') 'COORDINATES'
  do ipoin = 1,npoin
     write(80,'(1x,i9,3(1x,e16.8E3))') ipoin,(coord(idime,ipoin),idime=1,ndime)    
  end do
  write(80,'(a)') 'END_COORDINATES'
  write(80,'(a)') 'BOUNDARIES, ELEMENTS'
  do iboun = 1,nboun
     pnodb = nnode(ltypb(iboun))
     write(80,'(6(1x,i9))') iboun,(lnodb(inodb,iboun),inodb=1,pnodb),lelbo(iboun)
  end do
  write(80,'(a)') 'END_BOUNDARIES'
  !
  ! Boundary set and fixity file
  !
  do iboun = 1,nboun
     write(81,'(2(1x,i9))') iboun,kfl_codbo(iboun)
  end do
  !
  ! Domain file
  !
  write(82,1) '$--------------------------------------------------'
  write(82,1) 'DIMENSIONS'
  write(82,2) '  NODAL_POINTS=      ',npoin
  write(82,2) '  ELEMENTS=          ',nelem
  write(82,2) '  SPACE_DIMENSIONS=  ',ndime
  write(82,1) '  TYPES_OF_ELEMENTS=     TET04'
  write(82,2) '  BOUNDARIES=        ',nboun
  write(82,1) 'END_DIMENSIONS'
  write(82,1) '$--------------------------------------------------'
  write(82,1) 'STRATEGY'
  write(82,1) '  DOMAIN_INTEGRATION_POINTS: 0'
  write(82,1) '  OUTPUT:                    MESH'
  write(82,1) 'END_STRATEGY'
  write(82,1) '$--------------------------------------------------'
  write(82,6) 'GEOMETRY, WALL_DISTANCE= ',delta_dom,', ROUGHNESS= ',rough_dom
  write(82,1) '  INCLUDE mesh-new.geo.dat'
  write(82,1) 'END_GEOMETRY'
  write(82,1) '$--------------------------------------------------'
  write(82,1) 'SETS'
  write(82,1) '  BOUNDARIES'
  write(82,1) '    INCLUDE mesh-new.fix.dat'
  write(82,1) '  END_BOUNDARIES'
  write(82,1) 'END_SETS'
  write(82,1) '$--------------------------------------------------'
  write(82,1) 'BOUNDARY_CONDITIONS'
  write(82,1) '  ON_BOUNDARIES'
  write(82,1) '    INCLUDE mesh-new.fix.dat'
  write(82,1) '  END_ON_BOUNDARIES'
  write(82,1) '  GEOMETRICAL_CONDITIONS'

  wcond(1) = 'INFLOW=     '
  wcond(2) = 'FREESTREAM= '
  wcond(3) = 'WALL_LAW=   '
  wcond(4) = 'SYMMETRY=   '
  wcond(5) = 'WALL=       '
  wcond(6) = 'OUTFLOW=    '

  write(82,3) 'ANGLE=      ',geoan
  if( kfl_convx == 0 ) then
     write(82,7) 'CONVEX=     FREE'
  else if( kfl_convx == 1 ) then
     write(82,7) 'CONVEX=     EXNOR'
  else if( kfl_convx == 2 ) then
     write(82,7) 'CONVEX=     FIXED'
  else if( kfl_convx == 3 ) then
     write(82,7) 'CONVEX=     GEOMETRIC'
  end if
  do icond = 1,size(npbcs)
     if( npbcs(icond) /= 0 ) then
        ipara = 1
        write(82,4,advance='no') wcond(icond)
        do while( ipara < size(lsbcs,icond) .and. lsbcs(ipara,icond) /= 0 )
           write(82,5,advance='no') lsbcs(ipara,icond)           
           ipara = ipara + 1
        end do
        write(82,1) ' '
     end if
  end do
  write(82,1) '  END_GEOMETRICAL_CONDITIONS'
  write(82,1) 'END_BOUNDARY_CONDITIONS'
  write(82,1) '$--------------------------------------------------'
  !
  ! Close units
  !
  close(unit=80)
  close(unit=81)
  close(unit=82)

1 format(a)
2 format(a,i9)
3 format('    ',a12,e12.6)
4 format('    ',a12)
5 format(1x,i4)
6 format(a,e12.6,a,e12.6)
7 format('    ',a)

end subroutine outnew
