program alya2pos
  use def_kintyp
  implicit none
  real(rp),    pointer     :: gesca(:),gesc2(:)
  real(rp),    pointer     :: gevec(:,:),geve2(:,:)
  integer(ip), pointer     :: gisca(:),gisc2(:)
  integer(ip), pointer     :: givec(:,:)
  integer(ip)              :: ielem,ipoin,idime,inode,mnode,kpoin,jpoin
  integer(ip)              :: nelem,npoin,ndime,ii,ieles,ipois,npart
  integer(ip)              :: npoin_total,nelem_total,ipart,dummi,pdime
  integer(ip)              :: npoin_2,nelem_2
  integer(4)               :: ihead
  real(rp),    pointer     :: coord_loc(:,:)
  real(rp)                 :: dummr

  integer(ip)              :: kfl_markm
  integer(ip)              :: kfl_elimi
  integer(ip)              :: kfl_conve

  integer(ip)              :: npart_par
  real(rp)                 :: cutim
  character(5)             :: wwwww(10)
  character(5)             :: wwww8(10)
  integer(4)               :: iiiii(10)
  real(8)                  :: rrrrr(10)
  character(150)           :: namda,filna
  integer(ip), pointer     :: npoin_par(:)
  integer(ip), pointer     :: nelem_par(:)
  integer(ip), pointer     :: nboun_par(:)
  integer(ip), pointer     :: lsubd(:)

  integer(ip), pointer     :: markm(:)
  integer(ip), pointer     :: lnods(:,:)
  integer(ip), pointer     :: ltype(:)
  integer(ip), pointer     :: lelch(:)
  real(rp),    pointer     :: coord(:,:)
  integer(ip), pointer     :: lninv(:)

  !----------------------------------------------------------------------
  !
  ! Initialization
  !
  !----------------------------------------------------------------------

  call elmtyp()
  call GETARG(1_4,namda)
  if( trim(namda) == '' ) then
     print*,'Wrong probem name'
  end if

  kfl_markm = 3 ! 3=subdomain 0=element type, 4=element chacteristic
  kfl_elimi = 0 ! =1 eliminate boundary nodes

  !----------------------------------------------------------------------
  !
  ! Check ENDIAN
  !
  !----------------------------------------------------------------------

  kfl_conve = 1

  !----------------------------------------------------------------------
  !
  ! Run data: subdomain dimensions
  !
  !----------------------------------------------------------------------

  open(unit=8,file=trim(namda)//'.post.alyalog',form='formatted') 
  read(8,*) npart_par
  allocate( nelem_par(npart_par) )
  allocate( npoin_par(npart_par) )
  allocate( nboun_par(npart_par) )
  allocate( lsubd(npart_par) )
  do ipart = 1,npart_par
     read(8,*) dummi,nelem_par(ipart),npoin_par(ipart),nboun_par(ipart)
     lsubd(ipart) = 1
  end do
  close(8)
  lsubd = 0

lsubd(1 )=1      
lsubd(2 )=1      
lsubd(3 )=1      
lsubd(7 )=1      
lsubd(8)=1       
lsubd(9 )=1      
lsubd(27 )=1     
lsubd(30)=1      
lsubd(35)=1      
lsubd(37)=1      
lsubd(38)=1      
lsubd(42)=1      
lsubd(45 )=1     
lsubd(46)=1      
lsubd(47)=1      
lsubd(48 )=1     
lsubd(56)=1      
lsubd(61)=1      
 !lsubd(2)=1
 !lsubd(32)=1;lsubd(33)=1;lsubd(35)=1;lsubd(47)=1;lsubd(65)=1
 !lsubd(36)=1 
 !lsubd(55)=1 
  !lsubd( 5)=1 
  !lsubd( 3)=1 
  !lsubd( 2)=1 
  !lsubd(54)=1 
  !lsubd(51)=1 
  !lsubd(38)=1 
  !lsubd(10)=1 

  !lsubd = 0
  !lsubd(7:11)=1 

  !----------------------------------------------------------------------
  !
  ! Geometry
  !
  !----------------------------------------------------------------------

  if( npart_par == 1 ) kfl_elimi = 0
  if( npart_par == 1 ) lsubd = 1
  if( kfl_conve == 0 ) then
     open(unit=15,file=trim(namda)//'-LTYPE.post.alyabin',form='unformatted') 
     open(unit=16,file=trim(namda)//'-LNODS.post.alyabin',form='unformatted') 
     open(unit=17,file=trim(namda)//'-COORD.post.alyabin',form='unformatted') 
     if( npart_par > 1 ) then
        open(unit=18,file=trim(namda)//'-LNINV.post.alyabin',form='unformatted') 
     end if
     open(unit=19,file=trim(namda)//'-LELCH.post.alyabin',form='unformatted') 
  else
     open(unit=15,file=trim(namda)//'-LTYPE.post.alyabin',form='unformatted',convert='BIG_ENDIAN') 
     open(unit=16,file=trim(namda)//'-LNODS.post.alyabin',form='unformatted',convert='BIG_ENDIAN') 
     open(unit=17,file=trim(namda)//'-COORD.post.alyabin',form='unformatted',convert='BIG_ENDIAN') 
     if( npart_par > 1 ) then
        open(unit=18,file=trim(namda)//'-LNINV.post.alyabin',form='unformatted',convert='BIG_ENDIAN') 
     end if
     open(unit=19,file=trim(namda)//'-LELCH.post.alyabin',form='unformatted',convert='BIG_ENDIAN') 
  end if

  open(unit=49,file=trim(namda)//      '.post.alyafil',form=  'formatted') 

  do ii = 15,19

     if( ii < 18 .or. ii == 19 .or. (ii == 18 .and. npart_par > 1) ) then
        read(ii) wwww8(1) ! AlyaPost
        read(ii) wwww8(2) ! Version
        read(ii) wwww8(3) ! NAME
        read(ii) wwww8(4) ! SCALA/VECTO        
        read(ii) wwww8(5) ! NELEM/NPOIN/NBOUN  
        read(ii) wwww8(6) ! INTEG/REAL         
        read(ii) wwww8(7) ! 4BYTE/8BYTE        
        read(ii) wwww8(8) ! SEQUE/PARAL        
        read(ii) wwww8(9) ! NOFIL/FILTE
        read(ii) iiiii(1) 
        read(ii) iiiii(2) 
        read(ii) iiiii(3) 
        read(ii) rrrrr(1) 

        npart = int(iiiii(3),ip)
        if( npart /= npart_par ) stop
     end if

     if( ii == 15 ) then
        !
        ! LTYPE
        !
        ieles = 0
        nelem_total = iiiii(2) 
        write(6,'(a,i8)') 'Start reading ltype; nelem=',nelem_total
        allocate(  ltype(nelem_total) )

        do ipart = 1,npart
           if( lsubd(ipart) == 1 ) then
              read(ii) nelem 
              read(ii) ( ltype(ielem),ielem=ieles+1,ieles+nelem)
              do ielem = ieles+1,ieles+nelem
                 lexis(abs(ltype(ielem))) = 1
              end do
              ieles = ieles + nelem_par(ipart)
           else
              read(ii) nelem
              read(ii) ( dummi,ielem=1,nelem)
           end if
        end do

        write(6,'(a)') 'End   reading ltype'

     else if( ii == 16 ) then
        !
        ! LNODS
        !
        write(6,'(a)') 'Start reading lnods'
        ieles = 0
        ipois = 0
        nelem_total = iiiii(2) 
        mnode = iiiii(1)
        allocate(  lnods(mnode,nelem_total) )

        do ipart = 1,npart
           if( lsubd(ipart) == 1 ) then
              read(ii) nelem 
              read(ii) ( (lnods(inode,ielem),inode=1,mnode),ielem=ieles+1,ieles+nelem)
              do ielem = ieles+1,ieles+nelem
                 do inode = 1,nnode(ltype(ielem))
                    lnods(inode,ielem) = lnods(inode,ielem) + ipois
                 end do
              end do
              ieles = ieles + nelem_par(ipart)
              ipois = ipois + npoin_par(ipart)
           else
              read(ii) nelem
              read(ii) ( (dummi,inode=1,mnode),ielem=1,nelem)
           end if

        end do

        write(6,'(a)') 'End   reading lnods'

     else if( ii == 17 ) then
        !
        ! COORD
        !
        ipois = 0
        npoin_total = iiiii(2)
        ndime = iiiii(1)
        write(6,'(a,i8)') 'Start reading coordinates; npoin=',npoin_total
        allocate(  coord_loc(ndime,npoin_total) )

        do ipart = 1,npart 
           if( lsubd(ipart) == 1 ) then
              read(ii) npoin
              read(ii) ( (coord_loc(idime,ipoin),idime=1,ndime),ipoin=ipois+1,ipois+npoin)
              ipois = ipois + npoin_par(ipart)
           else
              read(ii) npoin
              read(ii) ( (dummi,idime=1,ndime),ipoin=1,npoin)
           end if
        end do

        write(6,'(a)') 'End   reading coordinates'

     else if( ii == 18 .and. npart_par > 1 ) then
        !
        ! LNINV
        !
        ipois = 0
        npoin_total = iiiii(2)
        write(6,'(a,i8)') 'Start reading permutation',npoin_total
        allocate(  lninv(npoin_total) )

        do ipart = 1,npart 
           if( lsubd(ipart) == 1 ) then
              read(ii) npoin
              read(ii) ( lninv(ipoin),ipoin=ipois+1,ipois+npoin)
              ipois = ipois + npoin_par(ipart)
           else
              read(ii) npoin
              read(ii) ( dummi,ipoin=1,npoin)              
           end if
        end do

        write(6,'(a)') 'End   reading permutation'

    else if( ii == 19 ) then
        !
        ! LELCH
        !
        ieles = 0
        nelem_total = iiiii(2) 
        write(6,'(a,i8)') 'Start reading lelch; nelem=',nelem_total
        allocate(  lelch(nelem_total) )

        do ipart = 1,npart
           read(ii) nelem 
           read(ii) ( lelch(ielem),ielem=ieles+1,ieles+nelem)
           ieles = ieles + nelem_par(ipart)
        end do

        write(6,'(a)') 'End   reading lelch'

     end if
  end do
  close(15)
  close(16)
  close(17)
  if( npart_par > 1 ) close(18)

  !----------------------------------------------------------------------
  !
  ! Eliminate repeated nodes
  !
  !----------------------------------------------------------------------

  if( kfl_elimi == 1 .and. npart_par > 0 ) then
     allocate( gisca(npoin_total) )
     do ipoin = 1,npoin_total
        gisca(ipoin) = 0
     end do
     kpoin = 0
     do ipoin = 1,npoin_total
        jpoin = lninv(ipoin)
        if( gisca(jpoin) == 0 ) then
           kpoin = kpoin + 1
           gisca(jpoin) = lninv(ipoin)
        end if
     end do
     deallocate( gisca )
     npoin_total = kpoin

     allocate(  coord(ndime,npoin_total) )    
     kpoin = 0
     do ipart = 1,npart 
        do ipoin = 1,npoin_par(ipart)
           kpoin = kpoin + 1
           jpoin = lninv(kpoin)
           do idime = 1,ndime
              coord(idime,jpoin) = coord_loc(idime,kpoin)
           end do
        end do
     end do
     deallocate( coord_loc )
     do ielem = 1,nelem_total
        do inode = 1,nnode(abs(ltype(ielem)))
           ipoin = lnods(inode,ielem)
           lnods(inode,ielem) = lninv(ipoin)
        end do
     end do
  else
     coord => coord_loc
  end if

  !----------------------------------------------------------------------
  !
  ! Output geometry
  !
  !----------------------------------------------------------------------

  npoin_2 = 0
  nelem_2 = 0
  do ipart = 1,npart_par
     if( lsubd(ipart) == 1 ) then
        npoin_2 = npoin_2 + npoin_par(ipart)
        nelem_2 = nelem_2 + nelem_par(ipart)
     end if
  end do

  if( 1 == 1 ) then
     open(unit=100,file=trim(namda)//'.post.msh',form='formatted') 
     if( kfl_markm == 0 ) then
        markm => ltype
     else if( kfl_markm == 1 .or. kfl_markm == 3 ) then
        allocate( markm(nelem_2) )
        ieles = 0
        do ipart = 1,npart
           if( lsubd(ipart) == 1 ) then
              do ielem = 1,nelem_par(ipart)
                 ieles = ieles + 1
                 markm(ieles) = ipart
              end do
           end if
        end do
     else if( kfl_markm == 4 ) then
        markm => lelch
     end if
     ipois     = 0
     ieles     = 0

     call gidele(&
          mnode,npoin_2,nelem_2,100_ip,lexis,ltype,lnods,coord,&
          markm,ieles,ipois,ndime,kfl_markm,namda)
     if( kfl_markm > 0 ) deallocate( markm )
     close(100)
  end if

  !----------------------------------------------------------------------
  !
  ! GiD format
  !
  !----------------------------------------------------------------------

  if( 1 == 1 ) then
     open(unit=11,file=trim(namda)//'.post.res',form='formatted') 
     write(11,1) 'GiD Post Results File 1.0'
     write(11,1) ' '
  end if

  do
     read(49,'(a)',end=999) filna
     if( kfl_conve == 0 ) then
        open(unit=10,file=trim(filna),form='unformatted')
     else
        open(unit=10,file=trim(filna),form='unformatted',convert='BIG_ENDIAN')
     end if
     ii = 10
     write(6,'(a)') 'Convert result file: '//trim(filna)

     read(ii) wwww8(1)
     read(ii) wwww8(2)
     read(ii) wwww8(3)
     read(ii) wwww8(4)
     read(ii) wwww8(5)
     read(ii) wwww8(6)
     read(ii) wwww8(7)
     read(ii) wwww8(8)
     read(ii) wwww8(9)
     read(ii) iiiii(1) 
     read(ii) iiiii(2) 
     read(ii) iiiii(3) 
     read(ii) rrrrr(1) 

     wwwww(1) = wwww8(1)(1:5)
     wwwww(2) = wwww8(2)(1:5)
     wwwww(3) = wwww8(3)(1:5)
     wwwww(4) = wwww8(4)(1:5)
     wwwww(5) = wwww8(5)(1:5)
     wwwww(6) = wwww8(6)(1:5)
     wwwww(7) = wwww8(7)(1:5)
     wwwww(8) = wwww8(8)(1:5)
     wwwww(9) = wwww8(9)(1:5)

     pdime       = int(iiiii(1),ip)
     npoin_total = int(iiiii(2),ip)

     if( wwwww(4) == 'SCALA' ) then
        !
        ! Scalar
        !
        if( wwwww(6) == 'INTEG' ) then

           write(11,2) wwwww(3),'ALYA',rrrrr(1),'Scalar'
           write(11,3) wwwww(3)
           write(11,1) 'Values'
           if( npoin /= 0 ) then

              if( wwwww(9) == 'FILTE' ) then
                 print*,'not coded'
                 stop
                 ipois = 0
                 if( kfl_elimi == 0 ) then                 
                    do ipart = 1,npart_par
                       if( lsubd(ipart) == 1 ) then 
                          read(ii) npoin
                          if( npoin > 0 ) then
                             allocate( gisca(npoin) )
                             allocate( gesca(npoin) )
                             read(ii) ( gisca(ipoin), ipoin=1,npoin )
                             read(ii) ( gesca(ipoin), ipoin=1,npoin )
                             do ipoin = 1,npoin
                                write(11,4) ipois+gisca(ipoin),gesca(ipoin) 
                             end do
                             deallocate( gisca )
                             deallocate( gesca )
                          end if
                          ipois = ipois + npoin_par(ipart)
                       else
                          read(ii) npoin
                          if( npoin > 0 ) then
                             read(ii) ( dummi,ipoin=1,npoin)
                             read(ii) ( dummr,ipoin=1,pdime*npoin)
                          end if
                       end if
                    end do
                 else
                    write(6,*) 'NOT CODED'
                 end if
              else
                 allocate( gisca(npoin_2) )  
                 ipois = 0
                 if( kfl_elimi == 0 ) then
                    do ipart = 1,npart_par
                       if( lsubd(ipart) == 1 ) then
                          read(ii) npoin
                          read(ii) (gisca(ipoin),ipoin=1+ipois,npoin+ipois)
                          ipois = ipois + npoin_par(ipart)
                       else
                          read(ii) npoin
                          read(ii) (dummi,ipoin=1,npoin)
                       end if
                    end do
                 else
                    do ipart = 1,npart_par
                       if( lsubd(ipart) == 1 ) then
                          read(ii) npoin
                          allocate( gisc2(npoin_par(ipart)) )
                          read(ii) ( gisc2(ipoin),ipoin=1,npoin)
                          do ipoin = 1,npoin
                             jpoin = lninv(ipoin+ipois)
                             gisca(jpoin) = real(gisc2(ipoin),rp)
                          end do
                          deallocate( gesc2 )
                          ipois = ipois + npoin_par(ipart)
                       else
                          read(ii) npoin
                          read(ii) ( dummi,ipoin=1,npoin)                    
                       end if
                    end do
                 end if
                 do ipoin = 1,npoin_2
                    write(11,6) ipoin,gisca(ipoin) 
                 end do
                 deallocate(gisca)
              end if
           end if
        else
           write(11,2) wwwww(3),'ALYA',rrrrr(1),'Scalar'
           write(11,3) wwwww(3)
           write(11,1) 'Values'
           if( npoin /= 0 ) then

              if( wwwww(9) == 'FILTE' ) then
                 ipois = 0
                 if( kfl_elimi == 0 ) then                 
                    do ipart = 1,npart_par
                       if( lsubd(ipart) == 1 ) then 
                          read(ii) npoin
                          if( npoin > 0 ) then
                             allocate( gisca(npoin) )
                             allocate( gesca(npoin) )
                             read(ii) ( gisca(ipoin), ipoin=1,npoin )
                             read(ii) ( gesca(ipoin), ipoin=1,npoin )
                             do ipoin = 1,npoin
                                write(11,4) ipois+gisca(ipoin),gesca(ipoin) 
                             end do
                             deallocate( gisca )
                             deallocate( gesca )
                          end if
                          ipois = ipois + npoin_par(ipart)
                       else
                          read(ii) npoin
                          if( npoin > 0 ) then
                             read(ii) ( dummi,ipoin=1,npoin)
                             read(ii) ( dummr,ipoin=1,pdime*npoin)
                          end if
                       end if
                    end do
                 else
                    write(6,*) 'NOT CODED'
                 end if
              else
                 allocate( gesca(npoin_2) )  
                 ipois = 0
                 if( kfl_elimi == 0 ) then
                    do ipart = 1,npart_par
                       if( lsubd(ipart) == 1 ) then
                          read(ii) npoin
                          read(ii) (gesca(ipoin),ipoin=1+ipois,npoin+ipois)
                          ipois = ipois + npoin_par(ipart)
                       else
                          read(ii) npoin
                          read(ii) (dummr,ipoin=1,npoin)
                       end if
                    end do
                 else
                    do ipart = 1,npart_par
                       if( lsubd(ipart) == 1 ) then
                          read(ii) npoin
                          allocate( gesc2(npoin_par(ipart)) )
                          read(ii) ( gesc2(ipoin),ipoin=1,npoin)
                          do ipoin = 1,npoin
                             jpoin = lninv(ipoin+ipois)
                             gesca(jpoin) = gesc2(ipoin)
                          end do
                          deallocate( gesc2 )
                          ipois = ipois + npoin_par(ipart)
                       else
                          read(ii) npoin
                          read(ii) ( dummr,ipoin=1,npoin)                    
                       end if
                    end do
                 end if
                 do ipoin = 1,npoin_2
                    write(11,4) ipoin,gesca(ipoin) 
                 end do
                 deallocate(gesca)
              end if
           end if
        end if
        write(11,1) 'End values'

     else if( wwwww(4) == 'VECTO' ) then
        !
        ! Vector
        !
        write(11,2) wwwww(3),'ALYA',rrrrr(1),'Vector'
        write(11,3) wwwww(3)//'_X,'//wwwww(3)//'_Y,'//wwwww(3)//'_Z'
        write(11,1) 'Values'
        if( npoin /= 0 ) then
           if( wwwww(9) == 'FILTE' ) then
              ipois = 0
              if( kfl_elimi == 0 ) then                 
                 do ipart = 1,npart_par
                    if( lsubd(ipart) == 1 ) then 
                       read(ii) npoin
                       if( npoin > 0 ) then
                          allocate( gisca(npoin) )
                          allocate( gevec(pdime,npoin) )
                          read(ii) ( gisca(ipoin), ipoin=1,npoin )
                          read(ii) ( (gevec(idime,ipoin),idime=1,pdime),ipoin=1,npoin)
                          do ipoin = 1,npoin
                             write(11,4) ipois+gisca(ipoin),(gevec(idime,ipoin),idime=1,pdime) 
                          end do
                          deallocate( gisca )
                          deallocate( gevec )
                       end if
                       ipois = ipois + npoin_par(ipart)
                    else
                       read(ii) npoin
                       if( npoin > 0 ) then
                          read(ii) ( dummi,ipoin=1,npoin)
                          read(ii) ( dummr,ipoin=1,pdime*npoin)
                       end if
                    end if
                 end do
              else
                 write(6,*) 'NOT CODED'
              end if
           else
              allocate( gevec(pdime,npoin_total) )  
              ipois = 0
              if( kfl_elimi == 0 ) then
                 do ipart = 1,npart_par
                    if( lsubd(ipart) == 1 ) then 
                       read(ii) npoin
                       read(ii) ( (gevec(idime,ipoin),idime=1,pdime),ipoin=1+ipois,npoin+ipois)
                       ipois = ipois + npoin_par(ipart)
                    else
                       read(ii) npoin
                       read(ii) ( dummr,ipoin=1,pdime*npoin)
                    end if
                 end do
              else
                 do ipart = 1,npart_par
                    if( lsubd(ipart) == 1 ) then 
                       read(ii) npoin
                       allocate( geve2(ndime,npoin_par(ipart)) )
                       read(ii) ( (geve2(idime,ipoin),idime=1,pdime),ipoin=1,npoin)
                       do ipoin = 1,npoin
                          jpoin = lninv(ipoin+ipois)
                          do idime = 1,ndime
                             gevec(idime,jpoin) = geve2(idime,ipoin)
                          end do
                       end do
                       deallocate( geve2 )
                       ipois = ipois + npoin_par(ipart)
                    else
                       read(ii) npoin
                       read(ii) ( dummr,ipoin=1,pdime*npoin)
                    end if
                 end do
              end if
              do ipoin = 1,npoin_2
                 write(11,4) ipoin,(gevec(idime,ipoin),idime=1,ndime) 
              end do
              deallocate(gevec)
           end if
        end if
        write(11,1) 'End values'

     end if

     close(10)

  end do

999 continue
  !
  ! GiD formats
  !
1 format(a)
2 format('Result ',a,' ',a,' ',e14.8,' ',a,' OnNodes')
3 format('ComponentNames ',a)
4 format(i9, 3(1x,e16.8E3))
5 format('ComponentNames ',a,a,a)
6 format(i9, 3(1x,i8))
  !
  ! Femview formats
  !
10 format(1x,i4,a1,a6,e12.5,32x,i2,i5)
20 format(1x,i2,2x,a5,3x,3i5)
30 format(1x,i2,2x,a5,3x,2i5)
40 format(1x,i2,i5,e12.5)
50 format(1x,i2)
  !
  ! Ensight Gold formats
  !
100 format(a)
110 format(i10)
120 format(e16.8E3)
  !
  ! VU format
  !
205 format('TEXTE Time(" ',e12.6,'");')
200 format('FIELD<double> ',a,'("',a,'",',i7,',',i12,');')
201 format('FIELD<float>  ',a,'("',a,'",',i7,',',i12,');')
202 format('FIELD<int>    ',a,'("',a,'",',i7,',',i12,');')
210 format('SOLUTION Solution( ) =',/,'{')
220 format('   VARIABLE ',a,'( ',a,',',a,',',a,',',a,');')
221 format('   VARIABLE ',a,'( ',a,',',a,',',a,',',a,',',a,');')
230 format('};')
  !
  ! Gnuplot format
  !
600 format(a)
610 format(i10)
620 format(10(1x,e16.8E3))
  !
  ! Alya ASCII format
  !
800 format(a)
810 format(i10)
820 format(i8,10(1x,e16.8E3))

end program alya2pos




subroutine heapsorti1(itask,nrows,ivin)
  !------------------------------------------------------------------------
  !****f* mathru/heapsorti1
  ! NAME
  !    heapsorti1
  ! DESCRIPTION
  !    Quik sorting. The element in ivin are sorting in:
  !    ITASK = 1 ... Decreasing value, i.e., ivin(1) > ivin(2) > ...
  !    ITASK = 2 ... Increasing value, i.e., ivin(1) < ivin(2) < ...
  ! INPUT
  !    ITASK ... 1,2 for decreasing, increasing order
  !    NROWS ... Size of IVIN
  !    IVIN .... Array to be ordered
  ! OUTPUT
  !    IVIN .... Ordered array
  ! USED BY
  !    
  !***
  !------------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  implicit none
  integer(ip), intent(in)    :: itask,nrows
  integer(ip), intent(inout) :: ivin(*) 
  integer(ip)                :: len, ir, ii, jj, iaux

  select case(itask)

  case(1)
     !
     ! Decreasing order
     !
     if(nrows<2) then
        return
     end if

     len = (nrows/2) + 1
     ir  = nrows

100  continue

     if (len>1) then
        len = len - 1
        iaux = ivin(len)
     else
        iaux = ivin(ir)
        ivin(ir) = ivin(1)

        ir = ir - 1

        if (ir==1) then
           ivin(1) = iaux
           return
        endif
     end if

     ii = len
     jj = len + len

200  if (jj<=ir) then
        if (jj<ir) then
           if ( ivin(jj)>ivin(jj+1) ) then
              jj = jj + 1
           endif
        endif

        if (iaux>ivin(jj) ) then
           ivin(ii) = ivin(jj)

           ii = jj
           jj = jj + jj
        else
           jj = ir + 1
        endif

        goto 200
     end if

     ivin(ii) = iaux

     goto 100

  case(2)
     !
     ! Increasing order
     !
     if(nrows<2) then
        return
     end if

     len = (nrows/2) + 1
     ir  = nrows

300  continue

     if (len>1) then
        len = len - 1
        iaux = ivin(len)
     else
        iaux = ivin(ir)
        ivin(ir) = ivin(1)

        ir = ir - 1

        if (ir==1) then
           ivin(1) = iaux
           return
        endif
     end if

     ii = len
     jj = len + len

400  if (jj<=ir) then
        if (jj<ir) then
           if ( ivin(jj)<ivin(jj+1) ) then
              jj = jj + 1
           endif
        endif

        if (iaux<ivin(jj) ) then
           ivin(ii) = ivin(jj)

           ii = jj
           jj = jj + jj
        else
           jj = ir + 1
        endif

        goto 400
     end if

     ivin(ii) = iaux

     goto 300

  case(3)

     if(nrows<2) then
        return
     end if

     do jj=2,nrows
        iaux=ivin(jj)
        do ii=jj-1,1,-1
           if(ivin(ii)<=iaux) exit
           ivin(ii+1)=ivin(ii)
        end do
        ivin(ii+1)=iaux
     end do

  end select

end subroutine heapsorti1
