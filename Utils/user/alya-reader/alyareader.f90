!-----------------------------------------------------------------------
!    
! alya-reader: alya input files reader. use it as a header for reading
!              input files
!
!-----------------------------------------------------------------------

program alyareader

  implicit none
  integer, parameter  :: ip = 8             ! 8-byte integer
  !integer, parameter  :: ip = 4             ! 4-byte integer
  integer, parameter  :: rp = kind(1.0d0)   ! double precision
  integer, parameter  :: ip4 = 4
  
  character(150) :: fil_geo, fil_fix, fil_dat, fil_bou,fil_fix_disk
  character(50)   :: wopos

  integer(ip4) :: istat
  integer(ip)  :: npoin,nelem,ndime,nboun,mnode,mnodb,lauxi(10)
  integer(ip)  :: ipoin,ielem,idime,iboun,inode,inodb,iauxi,icode,ndisco,kread
  real(rp)              :: xcoord,ycoord,zcoord,xvelo,yvelo,zvelo


  integer(ip) , pointer :: lnods(:,:),lnodb(:,:), ltype(:), lboco(:),kfl_condi(:)
  real(rp)    , pointer :: coord(:,:),condi(:,:)

  nullify(lnods)
  nullify(lnodb)
  nullify(ltype)
  nullify(lboco)
  nullify(coord)

  !
  ! file names and parameters
  !

  fil_geo= 'fensap.dom.geo'
  fil_fix= 'fensap.fix.dat'
  fil_fix_disk= 'fensap.fix-disk.dat'

  ndime =  3         ! space dimension
  npoin =  862361    ! number of nodes
  nelem = 2947998    ! number of elements
  nboun =  139778    ! number of boundary elements

  mnode = 6       ! max nnode for volume elements
  mnodb = 3       ! max nnode for boundary elements

  allocate(coord(ndime  ,npoin),stat=istat)  ! vector of nodal coordinates
  allocate(lnods(mnode  ,nelem),stat=istat)  ! volume elements connectivities
  allocate(lnodb(mnodb+1,nboun),stat=istat)  ! boundary elements connectivities
  allocate(ltype(nelem)      ,stat=istat)    ! element type
  allocate(lboco(nboun)      ,stat=istat)    ! boundary condition for boundary element
  allocate(condi(ndime  ,npoin),stat=istat)  ! boundary condition for the disk
  allocate(kfl_condi(npoin),stat=istat)      ! boundary condition for the disk (1 or 0)

  lnods= -1
  lboco= -1
  kfl_condi = 0

  open(25,file=adjustl(trim(fil_geo)),status='old')
  open(26,file=adjustl(trim(fil_fix)),status='old')

  ! read dom.geo file: four fields TYPES, ELEMENTS, COORDINATES, BOUNDARIES
  kread= 0
  do while (kread==0)
     read(25,'(a10)',IOSTAT=kread) wopos  
     wopos= adjustl(wopos)
     if (wopos(1:5) == 'TYPES') then
        write(6,*) 'Reading TYPES: ',nelem
        do ielem= 1,nelem
           read(25,*) iauxi,ltype(ielem)
        end do
        read(25,'(a10)',err=5) wopos
     else if (wopos(1:5) == 'ELEME') then
        write(6,*) 'Reading ELEMENTS: ',nelem
        do ielem= 1,nelem
           if (ltype(ielem) == 30) then
              mnode= 4   ! tetras
           else if (ltype(ielem) == 34) then
              mnode= 6   ! prisms
           else
              write (6,*) 'ELEMENT NOT PROGRAMMED= ',ltype(ielem)
              stop
           end if
           read(25,*) iauxi,lnods(1:mnode,ielem)
        end do
        read(25,'(a10)',err=5) wopos
     else if (wopos(1:5) == 'COORD') then
        write(6,*) 'Reading COORDINATES: ',npoin
        do ipoin= 1,npoin
           read(25,*) iauxi,coord(1:ndime,ipoin)
        end do
        read(25,'(a10)',err=5) wopos
     else if (wopos(1:5) == 'BOUND') then
        ! TO CORRECT: DIFFERENT MNODBs ARE NOT CONSIDERED!!
        write(6,*) 'Reading BOUNDARIES: ',nboun
        do iboun= 1,nboun
           read(25,*) iauxi,lnodb(1:mnodb+1,iboun)           
        end do
        read(25,'(a10)',err=5) wopos
     end if

  end do

5 continue

  ! read fix.dat file: boundary conditions codes
  ! TO CORRECT: FOR NOW, ONLY ON BOUNDS
  write(6,*) 'Reading conditions '
  kread= 0
  do while (kread==0)
     read(26,*,IOSTAT=kread) iboun,icode
     if (kread==0) lboco(iboun)= icode
  end do

  open(27,file=adjustl(trim(fil_fix_disk)),status='old')

  ndisco = 6  ! por ejemplo               
  do iboun= 1,nboun
     if (lboco(iboun) == ndisco) then
        do inodb=1,mnodb
           ipoin=lnodb(inodb,iboun)
           xcoord= coord(1,ipoin)
           ycoord= coord(2,ipoin)
           zcoord= coord(3,ipoin)

           ! operaciones
           ! .... determinar xvelo, yvelo, zvelo

           ! escribir vectores de condicion de contorno
           
           kfl_condi(ipoin) = 1
           condi(1,ipoin)= xvelo
           condi(2,ipoin)= yvelo
           condi(3,ipoin)= zvelo
           
        end do
     end if
  end do
  
  ! escribir fichero de condiciones nuevas

  do ipoin=1,npoin
     if (kfl_condi(ipoin) > 0) then
        write(27,*) ipoin, '111', condi(1:ndime,ipoin)
     end if
  end do

end program alyareader
