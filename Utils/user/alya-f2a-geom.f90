!-----------------------------------------------------------------------
!    
! alyaf2alya: translator from fensap (coming from icem) to alya
!
!-----------------------------------------------------------------------

program alyaf2alya
  !------------------
  !-- main
  !------------------

  !------------------
  !--    variables definition
  !----------------------------

  implicit none
  !integer, parameter  :: ip = 8             ! 8-byte integer
  integer, parameter  :: ip = 4             ! 4-byte integer
  integer, parameter  :: rp = kind(1.0d0)   ! double precision
  integer(4)          :: one4=1, two4=2, three4=3, four4=4, five4=5
  character(150)      :: &
       figrd,ficfg,chtyp,chout,tutyp,fil_ddat,fil_dgeo,fil_dbin,& 
       fil_lfix,fil_tfix,fil_bset,cfiel

  integer(ip)         :: &
       lun_ddat,lun_dgeo,lun_dbin,lun_lfix,lun_tfix,lun_figrd,lun_ficfg,lun_bset,&
       one=1,dummi
  integer(ip)         :: &
       nword,npoin,nfabc,nboun,ntotc,ndime,nelts,nelem, & 
       netet,nepyr,nehex,nepri,nnode,mnode,mnodb,nnodb,itype,itybi,knode(60)
  integer(ip)         :: & 
       iword,iauxi,ipoin,kpoin,jpoin,ifabc,iboun,ielem,idime, &
       ielts,ielty,nelty,ietet,iehex,iepyr,iepri,kelem,jelem, &
       inode,inodb,ifbou,ifpoi,kface,ieloc,neloc
  integer(ip)         :: icods(50,2),eltyp(50,2),eltyo(50),istat,nffix
  character(20)       :: wcods(50,2),welty,wngau,iffix,wbina
  real(rp)            :: pctge

  real(rp),       pointer     :: &
       coord(:,:)                     ! Coordinates
  integer(ip),    pointer     :: &
       ifico(:),&                     ! Fixity code
       lfano(:,:),&                   ! Face number of nodes and fixity code 
       lnods(:,:),&
       lnodb(:,:)
    
  !------------------
  !--    run
  !------------------


  lun_ddat  = 10
  lun_dgeo  = 11
  lun_dbin  = lun_dgeo
  lun_lfix  = 12
  lun_tfix  = 13
  lun_figrd = 14
  lun_ficfg = 15
  lun_bset  = 16
  knode     = 0

  call GETARG(one4  ,figrd)  ! *.grd
  !call GETARG(two4  ,ficfg) ! *.grd.cfg
  !call GETARG(three4,chtyp)  ! faces
  chtyp='faces'
  !call GETARG(four4 ,chout) ! bin file
  chout='asc'
  !call GETARG(five4 ,tutyp)

  write(6,*) '--|'
  write(6,*) '--| alya-fens2alya |-- '
  write(6,*) '--|'
  write(6,*) '--| Convert Fensap input files to Alya input files.'
  write(6,*) '--|'
  
  if( len(trim(figrd)) == 0 )then
     write(6,*) &
          '--| Usage: alya-f2alya "file"'
     write(6,*) &
          '--|        Files "file.grd" and "file.grd.cfg" must be located on current directory'
     write(6,*) '--|'
     write(6,*) '--|'
     write(6,*) '--| Try again.'
     write(6,*) '--|'
     stop
  end if

  fil_ddat = adjustl(trim(figrd))//'.dom.dat'
  fil_dgeo = adjustl(trim(figrd))//'.dom.geo'
  fil_dbin = adjustl(trim(figrd))//'.dom.bin'
  fil_lfix = adjustl(trim(figrd))//'.dom.fix'
  fil_tfix = adjustl(trim(figrd))//'-turbul-alya.fix'
  fil_bset = adjustl(trim(figrd))//'-alya.set' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  figrd=adjustl(trim(figrd))//'.grd'
  ficfg=adjustl(trim(figrd))//'.cfg'

  !
  !
  itybi=0
  iauxi=len(trim(chout))
  if (chout(1:iauxi)=='bin') itybi=1

  open(lun_ficfg,file=adjustl(trim(ficfg)),status='old')
  write(6,*) &
       '--| Reading cfg boundary conditions codes file...'
  read(lun_ficfg,*) nword
  write(6,*) '--|   Codes found  -->  Codes assigned:'
  write(6,*) '--|   fensap code   fensap word  -->   alya code    alya word'
!!$  print*,' nword ',nword
  do iword=1,nword
     read(lun_ficfg,*) icods(iword,1),iauxi,wcods(iword,1)     
     if (icods(iword,1)==11) then
        icods(iword,2)=1
        wcods(iword,2)='FAR_BOUNDARY'
     else if (icods(iword,1)==31) then
        icods(iword,2)=0
        wcods(iword,2)='EXIT'
     else if ((icods(iword,1)==21).or.(icods(iword,1)==22).or.(icods(iword,1)==32)) then
        icods(iword,2)=3
        wcods(iword,2)='SOLID_WALL'
     else if (icods(iword,1)==41) then
        icods(iword,2)=4
        wcods(iword,2)='SYMETRY_WALL'
     else if (icods(iword,1)==42) then
        icods(iword,2)=5
        wcods(iword,2)='SYMETRY_WALL2'
     else if (icods(iword,1)==45) then
        icods(iword,2)=2
        wcods(iword,2)='SPECIAL_BC_1'
     else if (icods(iword,1)==51) then
        icods(iword,2)=11 
        wcods(iword,2)='SPECIAL_BC_11'
     else if (icods(iword,1)==52) then
        icods(iword,2)=12 
        wcods(iword,2)='SPECIAL_BC_12'
     else if (icods(iword,1)==53) then
        icods(iword,2)=13 
        wcods(iword,2)='SPECIAL_BC_13'
     else if (icods(iword,1)==54) then
        icods(iword,2)=14 
        wcods(iword,2)='SPECIAL_BC_14' 
     end if
     if (icods(iword,1).gt.0) then
        write(6,100) ' --|    ',icods(iword,1),wcods(iword,1),' --> ',icods(iword,2),wcods(iword,2)
     end if
  end do
  write(6,*) '--| Done.'

100 format (a8,i3,10x,a10,2x,a5,2x,i3,10x,a20)

  close(lun_ficfg)

  open(lun_figrd,file=adjustl(trim(figrd)),status='old')
  write(6,*) &
       '--| Reading grd file...'
  read(lun_figrd,*) npoin,nfabc
  nboun=nfabc
  ndime=3
  ntotc=npoin*ndime
  allocate(coord(ndime,npoin),stat=istat)
  allocate(ifico(npoin)      ,stat=istat)

  read(lun_figrd,*) nelts
  write(6,*)'--|   Total nodes                   : ',npoin
  write(6,*)'--|   Total eltyps                  : ',nelts
  eltyp = 0
  netet = 0
  nepyr = 0
  nehex = 0
  nepri = 0
  eltyo = 0
  do ielts=1,nelts
     read(lun_figrd,*) ielty,nelty
     eltyp(ielty,1) = ielts
     eltyp(ielty,2) = nelty
     eltyo(ielts)   = ielty
     if (ielty==2) then
        write(6,*) '--|      Tetras  (fens=2 , alya=10) : ',nelty
        netet = nelty
     else if (ielty==5) then
        write(6,*) '--|      Pyramids(fens=5 , alya=12) : ',nelty
        nepyr = nelty
     else if (ielty==3) then
        write(6,*) '--|      Prisms  (fens=3 , alya=14) : ',nelty
        nepri = nelty
     else if (ielty==1) then
        write(6,*) '--|      Hexas   (fens=1 , alya=18) : ',nelty
        nehex = nelty
     else 
        write (6,*)'--|      Type  (',ielty,') is an unknown type of element!!'		
     end if
  end do
  read(lun_figrd,*) 
  read(lun_figrd,*) 

  nelem = netet+nepri+nepyr+nehex
  write(6,*)'--|   Total elements                : ',nelem
  write(6,*)'--|   '
  write(6,*)'--|   Size of coord vector          : ',ntotc
  write(6,*)'--|   Faces with boundary conditions: ',nfabc
  welty = ''
  wngau = ''
  mnode = 0  ! maximal nodes per element
  mnodb = 0  ! maximal nodes per boundary element
  if (netet > 0) then
!     welty = adjustl(trim(welty))//'10,'
     welty = adjustl(trim(welty))//'TET04,'
     wngau = adjustl(trim(wngau))//'1,'
     if (mnode<4) mnode = 4
     if (mnodb<3) mnodb = 3
  end if
  if (nepri > 0) then
!     welty = adjustl(trim(welty))//'14,'
     welty = adjustl(trim(welty))//'PEN06,'
     wngau = adjustl(trim(wngau))//'6,'
     if (mnode<6) mnode = 6
     if (mnodb<4) mnodb = 4
  end if
  if (nepyr > 0) then
!     welty = adjustl(trim(welty))//'12,'
     welty = adjustl(trim(welty))//'PYR05,'
     wngau = adjustl(trim(wngau))//'5,'
     if (mnode<5) mnode = 5
     if (mnodb<4) mnodb = 4
  end if
  if (nehex > 0) then
!     welty = adjustl(trim(welty))//'17,' 
     welty = adjustl(trim(welty))//'HEX08,'
     wngau = adjustl(trim(wngau))//'8,'
     if (mnode<8) mnode = 8
     if (mnodb<4) mnodb = 4
  end if
  iauxi=len(trim(welty))
  welty=adjustl(trim(welty(1:iauxi-1)))
  iauxi=len(trim(wngau))
  wngau=adjustl(trim(wngau(1:iauxi-1)))
  wngau='0'

  allocate(lnods(mnode+1,nelem),stat=istat)
  
  write(6,*)'--|   Writing dom.dat file...'

  open(lun_ddat,file=adjustl(trim(fil_ddat)),status='unknown')
  write(lun_ddat,800) '$--------------------------------------------------'
  write(lun_ddat,800) 'DIMENSIONS'
  write(lun_ddat,803) '  NODAL_POINTS=       '  , npoin 
  write(lun_ddat,803) '  ELEMENTS=           '  , nelem 
  write(lun_ddat,803) '  SPACE_DIMENSIONS=   '  , 3
  write(lun_ddat,802) '  TYPES_OF_ELEMENTS=  '  , trim(welty) 
  write(lun_ddat,803) '  BOUNDARIES=         '  , nboun 
  write(lun_ddat,800) 'END_DIMENSIONS'
  write(lun_ddat,800) '$--------------------------------------------------'
  write(lun_ddat,800) 'STRATEGY'
  write(lun_ddat,802) '  DOMAIN_INTEGRATION_POINTS: ' , trim(wngau) 
  write(lun_ddat,800) '  OUTPUT:                    MESH'
  write(lun_ddat,800) 'END_STRATEGY'
  write(lun_ddat,800) '$--------------------------------------------------'
  if (itybi == 0) then
     write(lun_ddat,800) 'GEOMETRY, WALL_DISTANCE=0.0'
     write(lun_ddat,802) '  INCLUDE  ' , adjustl(trim(fil_dgeo))
  else
     write(lun_ddat,800) 'GEOMETRY, READ_BINARY'
  end if
  write(lun_ddat,800) 'END_GEOMETRY'
  write(lun_ddat,800) '$--------------------------------------------------'
  write(lun_ddat,800) 'SETS'
  write(lun_ddat,800) '  BOUNDARIES'
  write(lun_ddat,802) '    INCLUDE ',adjustl(trim(fil_lfix))
  write(lun_ddat,800) '  END_BOUNDARIES'
  write(lun_ddat,800) 'END_SETS'
  write(lun_ddat,800) '$--------------------------------------------------'
  write(lun_ddat,800) 'BOUNDARY_CONDITIONS'
  write(lun_ddat,800) '  ON_BOUNDARIES'
  write(lun_ddat,802) '    INCLUDE ',adjustl(trim(fil_lfix))
  write(lun_ddat,800) '  END_ON_BOUNDARIES'
  write(lun_ddat,800) '   GEOMETRICAL_CONDITIONS'
  write(lun_ddat,800) '     FREESTREAM:   ???'
  write(lun_ddat,800) '     WALL_LAW:     ???'
  write(lun_ddat,800) '     SYMMETRY:     ???'
  write(lun_ddat,800) '     INFLOW:       ???'
  write(lun_ddat,800) '   END_GEOMETRICAL_CONDITIONS'
  write(lun_ddat,800) 'END_BOUNDARY_CONDITIONS'
  write(lun_ddat,800) '$--------------------------------------------------'
  close(lun_ddat)


800 format(a)
802 format(a,a)
803 format(a,i11)

  write(6,*)'--|   Done.'
  write(6,*)'--|   '

  write(6,*)'--|   Reading coordinates... '
  iauxi= 0
  kpoin= npoin/10
  jpoin= kpoin
  do ipoin=1,npoin
     if (ipoin == kpoin) then
        pctge = real(ipoin)/real(npoin) * 100.0;
        write(6,200)' --|     ',pctge,'% done... '
        kpoin = kpoin + jpoin
     end if
     read(lun_figrd,*) (coord(idime,ipoin),idime=1,ndime),ifico(ipoin)     
  end do
  write(6,*)'--|   '
  write(6,*)'--|',ipoin,' nodes were read. Done.'
  write(6,*)'--|   '
  write(6,*)'--|   Reading elements... '

  ietet=0
  iehex=0
  iepyr=0
  iepri=0
  lnods=0

  iauxi=0
  kelem= nelem/10
  jelem= kelem

  ielem= 0
  do ielts= 1,nelts
     ielty = eltyo(ielts)
     nelty = eltyp(ielty,2)     

     if (ielty==2) then
        write(6,*)'--|     Reading  tetras alya_type=30 nnode=4... '
        itype= 30
        nnode=  4
     else if (ielty==5) then
        write(6,*)'--|     Reading  pyramids alya_type=32 nnode=5... '
        itype= 32
        nnode=  5
     else if (ielty==3) then
        write(6,*)'--|     Reading  prisms alya_type=34 nnode=6... '
        itype= 34
        nnode=  6
     else if (ielty==1) then
        write(6,*)'--|     Reading  hexas alya_type=37 nnode=8... '
        itype= 37
        nnode=  8
     else 
        write(6,*)'--|     Type  (',ielty,') is an unknown type of element!!'		
        stop
     end if

     do ieloc=1,nelty
        ielem=ielem+1
        if (ielem == kelem) then
           pctge = real(ielem)/real(nelem) * 100.0;
           write(6,200)' --|     ',pctge,'% done... '
           kelem = kelem + jelem
        end if
        
        lnods(mnode+1,ielem) = itype
        read(lun_figrd,*) (lnods(inode,ielem),inode=1,nnode)
     end do

  end do

  write(6,*)'--|   '
  write(6,*)'--|',ielem,' elements were read. Done.'
  write(6,*)'--|   '

200 format(a9,2x,f10.0,a10)

  if (itybi == 0) then
     write(6,*) &
          '--|   Writing geometry ASCII file ',adjustl(trim(fil_dgeo))
     open(lun_dgeo,file=adjustl(trim(fil_dgeo)),status='unknown')
  else
     write(6,*) &
          '--|   Writing geometry BINARY file ',adjustl(trim(fil_dbin))
     open(lun_dbin,file=adjustl(trim(fil_dbin)),form='unformatted')
  end if

  knode(10) = 3
  knode(12) = 4
  knode(30) = 4
  knode(32) = 5
  knode(34) = 6
  knode(37) = 8

  if (itybi==0) then
     
     write(6,*) &
          '--|      Types... '
     write(lun_dgeo,*) 'TYPES'
     do ielem=1,nelem
        write(lun_dgeo,*) ielem,lnods(mnode+1,ielem)
     end do
     write(lun_dgeo,*) 'END_TYPES'
     
     write(lun_dgeo,*) 'ELEMENTS'
     write(6,*) &
          '--|      Elements... '
     do ielem=1,nelem
        nnode = knode(lnods(mnode+1,ielem))
        write(lun_dgeo,300) ielem,(lnods(inode,ielem),inode=1,nnode)
     end do
     
300  format(10(i12,2x))
     
     write(lun_dgeo,*) 'END_ELEMENTS'
     write(6,*) &
          '--|      Coordinates... '
     write(lun_dgeo,*) 'COORDINATES'
     do ipoin=1,npoin
        write(lun_dgeo,400) ipoin,(coord(idime,ipoin),idime=1,ndime)     
     end do
     write(lun_dgeo,*) 'END_COORDINATES'
     
400  format(i12,3(2x,e25.15))
     
     deallocate(coord)
     
  else if (itybi == 1) then

     dummi = 123456
     write(lun_dbin) dummi
     cfiel='LTYPE'
     write(lun_dbin) cfiel(1:20),one
     write(lun_dbin) (lnods(mnode+1,ielem),ielem=1,nelem)
     cfiel='LNODS'
     write(lun_dbin) cfiel(1:20),one
     write(lun_dbin) ((lnods(inode,ielem),inode=1,knode(lnods(mnode+1,ielem))),ielem=1,nelem)
     cfiel='COORD'
     write(lun_dbin) cfiel(1:20),one
     write(lun_dbin) ((coord(idime,ipoin),idime=1,ndime),ipoin=1,npoin)
     
  end if
  
  
  write(6,*) '--|   Done. '
  write(6,*) '--|    '
  
  write(6,*) &
       '--|   Read and build boundary conditions over boundary faces... '
  write(6,*) '--|    '
  
  
  !	
  !   Face number Node Ordering:
  !
  !   4-node tetrahedron
  !   ielgeom = 2
  !   1 2-3-4
  !   2 1-4-3
  !   3 1-2-4
  !   4 1-3-2
  ! 
  !   5-node pyramid
  !   ielgeom = 5
  !   1 1-4-3-2
  !   2 1-2-5
  !   3 2-3-5
  !   4 3-4-5
  !   5 4-1-5
  ! 
  !   6-node prism
  !   ielgeom = 3
  !   1 2-3-6-5
  !   2 3-1-4-6
  !   3 1-2-5-4
  !   4 1-3-2
  !   5 4-5-6
  ! 
  !   8-node brick
  !   ielgeom = 1
  !   1 1-4-3-2
  !   2 1-5-8-4
  !   3 4-8-7-3
  !   4 2-3-7-6
  !   5 1-2-6-5
  !   6 5-6-7-8
  !	
  
  allocate(lnodb(mnodb+2,nboun),stat=istat)
  allocate(lfano(2,nboun)      ,stat=istat)
  lfano = 0
  

  do iboun=1,nboun
     read(lun_figrd,*) ifbou,kface,kelem
     
     ielty= lnods(mnode+1,kelem)
     lnodb(mnodb+2,iboun) = kelem
     lfano(2,iboun)       = ifbou
     if (ielty == 30) then
        lnodb(mnodb+1,iboun) = 10   ! triangle
        if (kface == 1) then
           lfano(1,iboun)=3
           lnodb(1,iboun)=lnods(2,kelem)
           lnodb(2,iboun)=lnods(3,kelem)
           lnodb(3,iboun)=lnods(4,kelem)
        else if (kface == 2) then
           lfano(1,iboun)=3
           lnodb(1,iboun)=lnods(1,kelem)
           lnodb(2,iboun)=lnods(3,kelem)
           lnodb(3,iboun)=lnods(4,kelem)
        else if (kface == 3) then
           lfano(1,iboun)=3
           lnodb(1,iboun)=lnods(1,kelem)
           lnodb(2,iboun)=lnods(2,kelem)
           lnodb(3,iboun)=lnods(4,kelem)
        else if (kface == 4) then
           lfano(1,iboun)=3
           lnodb(1,iboun)=lnods(1,kelem)
           lnodb(2,iboun)=lnods(2,kelem)
           lnodb(3,iboun)=lnods(3,kelem)
        end if

     else if (ielty == 32) then
        if (kface == 1) then
           lnodb(mnodb+1,iboun) = 12   ! quad 4
           lfano(1,iboun)=4
           lnodb(1,iboun)=lnods(1,kelem)
           lnodb(2,iboun)=lnods(4,kelem)
           lnodb(3,iboun)=lnods(3,kelem)
           lnodb(4,iboun)=lnods(2,kelem)
        else if (kface == 2) then
           lnodb(mnodb+1,iboun) = 10   ! triangle
           lfano(1,iboun)=3
           lnodb(1,iboun)=lnods(1,kelem)
           lnodb(2,iboun)=lnods(2,kelem)
           lnodb(3,iboun)=lnods(5,kelem)
        else if (kface == 3) then
           lnodb(mnodb+1,iboun) = 10   ! triangle
           lfano(1,iboun)=3
           lnodb(1,iboun)=lnods(2,kelem)
           lnodb(2,iboun)=lnods(3,kelem)
           lnodb(3,iboun)=lnods(5,kelem)
        else if (kface == 4) then
           lnodb(mnodb+1,iboun) = 10   ! triangle
           lfano(1,iboun)=3
           lnodb(1,iboun)=lnods(3,kelem)
           lnodb(2,iboun)=lnods(4,kelem)
           lnodb(3,iboun)=lnods(5,kelem)
        else if (kface == 5) then
           lnodb(mnodb+1,iboun) = 10   ! triangle
           lfano(1,iboun)=3
           lnodb(1,iboun)=lnods(4,kelem)
           lnodb(2,iboun)=lnods(1,kelem)
           lnodb(3,iboun)=lnods(5,kelem)
        end if

     else if (ielty == 34) then
        if (kface == 1) then
           lnodb(mnodb+1,iboun) = 12   ! quad 4
           lfano(1,iboun)=4
           lnodb(1,iboun)=lnods(2,kelem)
           lnodb(2,iboun)=lnods(3,kelem)
           lnodb(3,iboun)=lnods(6,kelem)
           lnodb(4,iboun)=lnods(5,kelem)
        else if (kface == 2) then
           lnodb(mnodb+1,iboun) = 12   ! quad 4
           lfano(1,iboun)=4
           lnodb(1,iboun)=lnods(3,kelem)
           lnodb(2,iboun)=lnods(1,kelem)
           lnodb(3,iboun)=lnods(4,kelem)
           lnodb(4,iboun)=lnods(6,kelem)
        else if (kface == 3) then
           lnodb(mnodb+1,iboun) = 12   ! quad 4
           lfano(1,iboun)=4
           lnodb(1,iboun)=lnods(1,kelem)
           lnodb(2,iboun)=lnods(2,kelem)
           lnodb(3,iboun)=lnods(5,kelem)
           lnodb(4,iboun)=lnods(4,kelem)
        else if (kface == 4) then
           lnodb(mnodb+1,iboun) = 10   ! triangle
           lfano(1,iboun)=3
           lnodb(1,iboun)=lnods(1,kelem)
           lnodb(2,iboun)=lnods(3,kelem)
           lnodb(3,iboun)=lnods(2,kelem)
        else if (kface == 5) then
           lnodb(mnodb+1,iboun) = 10   ! triangle
           lfano(1,iboun)=3
           lnodb(1,iboun)=lnods(4,kelem)
           lnodb(2,iboun)=lnods(5,kelem)
           lnodb(3,iboun)=lnods(6,kelem)
        end if

     else if (ielty == 37) then
        if (kface == 1) then
           lnodb(mnodb+1,iboun) = 12   ! quad 4
           lfano(1,iboun)=4
           lnodb(1,iboun)=lnods(1,kelem)
           lnodb(2,iboun)=lnods(4,kelem)
           lnodb(3,iboun)=lnods(3,kelem)
           lnodb(4,iboun)=lnods(2,kelem)
        else if (kface == 2) then
           lnodb(mnodb+1,iboun) = 12   ! quad 4
           lfano(1,iboun)=4
           lnodb(1,iboun)=lnods(1,kelem)
           lnodb(2,iboun)=lnods(5,kelem)
           lnodb(3,iboun)=lnods(8,kelem)
           lnodb(4,iboun)=lnods(4,kelem)

        else if (kface == 3) then
           lnodb(mnodb+1,iboun) = 12   ! quad 4
           lfano(1,iboun)=4
           lnodb(1,iboun)=lnods(4,kelem)
           lnodb(2,iboun)=lnods(8,kelem)
           lnodb(3,iboun)=lnods(7,kelem)
           lnodb(4,iboun)=lnods(3,kelem)

        else if (kface == 4) then
           lnodb(mnodb+1,iboun) = 12   ! quad 4
           lfano(1,iboun)=4
           lnodb(1,iboun)=lnods(2,kelem)
           lnodb(2,iboun)=lnods(3,kelem)
           lnodb(3,iboun)=lnods(7,kelem)
           lnodb(4,iboun)=lnods(6,kelem)

        else if (kface == 5) then
           lnodb(mnodb+1,iboun) = 12   ! quad 4
           lfano(1,iboun)=4
           lnodb(1,iboun)=lnods(1,kelem)
           lnodb(2,iboun)=lnods(2,kelem)
           lnodb(3,iboun)=lnods(6,kelem)
           lnodb(4,iboun)=lnods(5,kelem)

        else if (kface == 6) then
           lnodb(mnodb+1,iboun) = 12   ! quad 4
           lfano(1,iboun)=4
           lnodb(1,iboun)=lnods(5,kelem)
           lnodb(2,iboun)=lnods(6,kelem)
           lnodb(3,iboun)=lnods(7,kelem)
           lnodb(4,iboun)=lnods(8,kelem)
        end if


     end if

  end do

  write(6,*) &
       '--|   Done.'
  write(6,*) & 
       '--|    '
  write(6,*) &
       '--|   Writing boundary elements and sets (add by LAFO)... '
  write(6,*) & 
       '--|    '

  !open(lun_bset,file=adjustl(trim(fil_bset)),status='unknown')

  if (itybi==0) then

     write(lun_dgeo,*) 'BOUNDARIES, ELEMENTS'
     !write(lun_bset,*) 'BOUNDARIES'
     do iboun=1,nboun
        nnodb=lfano(1,iboun)
        write(lun_dgeo,500) &
             iboun,(lnodb(inodb,iboun),inodb=1,nnodb),lnodb(mnodb+2,iboun)  
        !write(lun_bset,500) &
        !     iboun,lfano(2,iboun)       
     end do
     write(lun_dgeo,*) 'END_BOUNDARIES'
     !write(lun_bset,*) 'END_BOUNDARIES'
     
500  format(6(i12,2x))
     
     write(lun_dgeo,*) 'SKEW_SYSTEMS'
     write(lun_dgeo,*) 'END_SKEW_SYSTEMS'
     
     close(lun_dgeo)

  else if (itybi==1) then
     
     if(nboun/=0) then
        cfiel='LTYPB'
        write(lun_dbin) cfiel(1:20),nboun
        write(lun_dbin) (lnodb(mnodb+1,iboun),iboun=1,nboun)
        cfiel='LNODB'
        write(lun_dbin) cfiel(1:20),nboun
        write(lun_dbin) ((lnodb(inodb,iboun),inodb=1,knode(lnodb(mnodb+1,iboun))),iboun=1,nboun)
        cfiel='LBOEL'
        write(lun_dbin) cfiel(1:20),one
        write(lun_dbin) (lnodb(mnodb+2,iboun),iboun=1,nboun)
     end if
     cfiel='END_FILE'
     write(lun_dbin) cfiel(1:20),one    

     close(lun_dbin)

     !Write the set file (if binary geo file is chose (Lafo) )
     !write(lun_bset,*) 'BOUNDARIES'
     !do iboun=1,nboun
     !   nnodb=lfano(1,iboun)
     !   write(lun_bset,500) &
     !        iboun,lfano(2,iboun)       
     !end do
     !write(lun_bset,*) 'END_BOUNDARIES'

  end if

  close(lun_figrd)

  write(6,*) &
       '--|   Done.'
  write(6,*) & 
       '--|    '
  write(6,*) &
       '--|   Writing laminar boundary conditions in file:'
  write(6,*) &
       '--|      --> ',adjustl(trim(fil_lfix))

  open(lun_lfix,file=adjustl(trim(fil_lfix)),status='unknown')
  write(6,*) & 
       '--|   As required, boundary conditions set on: '
  write(6,*) & 
       '--|      --> ',adjustl(trim(chtyp))
  iauxi=len(trim(chtyp))
  if (chtyp(1:iauxi)=='faces') then
     !write(lun_lfix,*) 'ON_BOUNDARIES, CODED'
     do iboun=1,nboun
        iffix= "nofix"
        ifbou= lfano(2,iboun)



       if (ifbou >= 11 .AND. ifbou <= 19) then
           nffix= ifbou
        else if (ifbou >= 21 .AND. ifbou <= 29) then
           nffix= ifbou
        else if (ifbou >= 31 .AND. ifbou <= 39) then
           nffix= ifbou
        else if (ifbou >= 41 .AND. ifbou <= 44) then
           nffix= ifbou
        else if (ifbou >= 45 .AND. ifbou <= 49) then
           nffix= ifbou
	end if


       ! iauxi=len(trim(iffix))
       ! if (.not.(iffix(1:iauxi)=='nofix')) &
       !      write(lun_lfix,700) iboun,iffix
       !      !write(lun_lfix,700) iboun,nffix 

       write(lun_lfix,*) iboun,nffix 

     end do
     !write(lun_lfix,*) 'END_ON_BOUNDARIES'

700  format(i12,2x,a6)

  else if (chtyp(1:iauxi)=='nodes') then
     write(lun_lfix,*) 'ON_NODES, CODED'
     do ipoin=1,npoin
        iffix= "nofix"
        ifpoi= ifico(ipoin)
	if (ifpoi >= 11 .AND. ifpoi <= 19) then
           iffix= "1"
        else if (ifpoi >= 21 .AND. ifpoi <= 29) then
           iffix= "3"
        else if (ifpoi >= 31 .AND. ifpoi <= 39) then
           iffix= "0"
        else if (ifpoi >= 41 .AND. ifpoi <= 44) then
           iffix= "4"
        else if (ifpoi >= 45 .AND. ifpoi <= 49) then
           iffix= "5"
        end if
        iauxi=len(trim(iffix))
        if (.not.(iffix(1:iauxi)=='nofix')) &
             write(lun_lfix,700) ipoin,iffix             
     end do
     write(lun_lfix,*) 'END_ON_NODES'     
  else 
     write(6,*) '--|   Fail to set the boundary conditions character.. '
     write(6,*) '--|   You have entered: '
     write(6,*) '--|  ',chtyp
     write(6,*) '--|   Neither "nodes" nor "faces".'
     write(6,*) '--|'
     write(6,*) '--| Bye!'
     write(6,*) '--|'     
     stop
  end if


600 format(a7,2x,a7)

  close(lun_lfix)
  !close(lun_bset)

  !------------------
  !--    say good bye
  !------------------
  write(6,*) '--|'
  write(6,*) '--| Bye!'
  write(6,*) '--|'
  
end program alyaf2alya


  
