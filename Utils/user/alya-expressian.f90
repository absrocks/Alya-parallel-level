program expressian  

  implicit none

  integer, parameter  :: ip = selected_int_kind(9)   ! 4-byte integer
!!!  integer, parameter  :: ip = selected_int_kind(18)  ! 8-byte integer
  integer, parameter  :: rp = kind(1.0d0)            ! double precision 

  integer(ip)      ::&
       idime,ipoin,ielem,inode,iboun,nnode,nnodb,ndime,npoin,nelem,nboun,nzero,&
       ntolo,itolo,iello,nello,iboco,outbo,&
       nelex,neley,nelez,nepox,nepoy,nepoz,ielex,ieley,jeley,ielez,iepox,iepoy,iepoz,&
       lun_dgeo,lun_dbin,lun_scre,istat,inodb,one,ityme,&
       kboty,kflbo,kflfi,kflbl,lnohe(27),nboco(20),ntics,nticx,nticy,nticz,&
       ktype,nelau,itybi,knode(60)
  real(rp) ::&
       xmaxi,ymaxi,zmaxi,cabox(2,2),dxref,dzref,dyref,xauxi,ratio,ymaux,&
       dxcoo(10000),dycoo(10000),dzcoo(10000),copoi(3),cpcon,epcon,fpcon,facto
  real(rp)   ,    pointer     :: &
       coord(:,:)   
  integer(ip),    pointer     :: &
       lnods(:,:),&
       lnodb(:,:)
  character(3) ::&
       iflbo
  character(5) ::&
       cboco,chtyp
  character(20) ::&
       cfiel
  
  write(6,*) '--|'
  write(6,*) '--| Cartessian Express |-- '
  write(6,*) '--|'
  write(6,*) '--| Create cartessian meshes'
  write(6,*) '--|'
  write(6,*) '--|'
  if (ip == 4) then
     write(6,*) '--| INTEGERS DEFINED AS SHORT: INTEGER*4!!!'
  else 
     write(6,*) '--| INTEGERS DEFINED AS LONG:  INTEGER*8!!!'     
  end if
  write(6,*) '--|'

  lun_dgeo = 11
  lun_dbin = lun_dgeo
  lun_scre = 6
  one = 1

  nzero= 0_ip
  ktype= 37
  chtyp= 'HEX08'

!!  xmaxi= 250.0_rp
!!  ymaxi=  50.0_rp
!!  zmaxi=  50.0_rp
  write(6,'(a)',advance='no') ' --|  xmaxi = '
  read (5,*) xmaxi
  write(6,'(a)',advance='no') ' --|  ymaxi = '
  read (5,*) ymaxi
  write(6,'(a)',advance='no') ' --|  zmaxi = '
  read (5,*) zmaxi

!!$  xmaxi=  25.0_rp
!!$  ymaxi=  5.0_rp
!!$  zmaxi=  5.0_rp

  cabox(1,1) =          (xmaxi / 4.0_rp)
  cabox(2,1) = 2.0_rp * (xmaxi / 4.0_rp)
  cabox(1,2) = 2.0_rp * (ymaxi / 5.0_rp)
  cabox(2,2) = 3.0_rp * (ymaxi / 5.0_rp)

  kflbo =1    ! compute boundary elements
  kflfi =1    ! compute boundary conditions
  kflbl =1    ! boundary layer refinement
  ityme =1    ! Q1 elements (8-node bricks)  , ityme=2, Q2 (27-node bricks)


  write(6,*) '--|  The total number of mesh points is defined as: '
  write(6,*) '--|     NEPOX  *  NEPOY   *  NEPOZ'
  write(6,*) '--|  where'
  write(6,*) '--|     NEPOZ  =     nticz '
  write(6,*) '--|     NEPOY  =     nticy '
  write(6,*) '--|     NEPOX  =     nticx '
  write(6,*) '--|  and the reference value "nticx,nticy,nticz" is defined by the user.'
  write(6,*) '--|  When Q2 elements are chosen, the total number of elements '
  write(6,*) '--|  is proportionally reduced. '
  write(6,*) '--|  NOTE FOR Q2 ELEMENTS: ntics must be an odd number!'
  write(6,*) '--|  '
  write(6,*) '--|  So what are ntics values ? '
  write(6,'(a)',advance='no') ' --|  nticx = '
  read (5,*) nticx
  write(6,'(a)',advance='no') ' --|  nticy = '
  read (5,*) nticy
  write(6,'(a)',advance='no') ' --|  nticz = '
  read (5,*) nticz

  write(6,*) '--|  '
  write(6,'(a)',advance='no') ' --|  Output file is BINARY [1] or ASCII [2] = ? '
  read (5,*) itybi

  write(6,*) '--|  '
  write(6,'(a)',advance='no') ' --|  Q-1 8-node [1] or Q-2 27-node [2] bricks = ? '
  read (5,*) ityme

  do while (ityme == 2 .and. (mod(nticx,2)==0.or.mod(nticy,2)==0.or.mod(nticz,2)==0)) 
     write(6,*) '--|   ntics must be an ODD (3,5,101,...) number.'
     write(6,*) '--|  So what is ntics value ? '
     write(6,'(a)',advance='no') ' --|  ntics = '
     read (5,*) ntics
  end do

  write(6,*) '--|  '
  write(6,'(a)',advance='no') ' --|  Compute boundary elements YES [1] or NO [0] = ? '
  read (5,*) kflbo

  write(6,*) '--|  '
  write(6,'(a)',advance='no') ' --|  Boundary layer and Carter refinement YES [1] or NO [0] = ? '
  read (5,*) kflbl

  if (kflbl == 0) then
     write(6,*) '--|  '
     write(6,'(a)',advance='no') ' --|  Compute boundary conditions YES [1] or NO [0] = ? '
     read (5,*) kflfi
  else 
     kflfi=1
  end if
  
  open(10,file='cart-express.dom.dat',status='unknown')

  if (itybi==2) then
     write(6,*) '--|  ASCII output file. '
     open(lun_dgeo,file='cart-express.dom.geo',status='unknown')
  else
     write(6,*) '--|  BINARY output file. '
     open(lun_dbin,file='cart-express.geo.bin',form='unformatted')     
  end if


!  nepox= ntics

  nepox= nticx
  nepoy= nticy
  nepoz= nticz
  
  neley= nepoy-1
  nelez= nepoz-1
  nelex= nepox-1
  
  dxref= xmaxi / real(nelex)
  dyref= ymaxi / real(neley)
  dzref= zmaxi / real(nelez)

  dxcoo= dxref
  dycoo= dyref
  dzcoo= dzref
  
  nelex=nelex/ityme
  neley=neley/ityme
  nelez=nelez/ityme

  open(20,file='mesh-size.plo',status='unknown')
  
  if (kflbl == 1) then
     ! from xenophontos (2002)
     fpcon = 1.0_rp
     epcon = 0.25_rp
     cpcon = 1.0_rp - exp(-2.0_rp / (2.0_rp * fpcon + 1.0_rp) * epcon )
     facto = - epcon * 0.5_rp * (2.0_rp * fpcon + 1.0_rp)
     
     ! z-direction
     
     xauxi = 0.0_rp
     do ielez=1,nepoz-1
        dzcoo(ielez)= facto * log(1.0_rp - 2.0_rp * cpcon * real(ielez) / real(nepoz-1))  
        xauxi = xauxi + zmaxi * dzcoo(ielez)
     end do
     
     ratio = zmaxi / xauxi
     
     xauxi = 0.0_rp
     do ielez=1,nepoz-1
        dzcoo(ielez) = zmaxi * dzcoo(ielez) * ratio
        xauxi = xauxi + dzcoo(ielez)
        !!     write (20,*) ielez,xauxi
     end do
     
     ! y-direction
     
     fpcon = 1.0_rp
     epcon = 0.025_rp
     cpcon = 1.0_rp - exp(-2.0_rp / (2.0_rp * fpcon + 1.0_rp) * epcon )
     facto = - epcon * 0.5_rp * (2.0_rp * fpcon + 1.0_rp)
     
     xauxi = 0.0_rp
     nelau = (nepoy-1) / 2
     ymaux = ymaxi / 2.0_rp
     
     !  do ieley=nelau,neley
     !     dycoo(ieley)= facto * log(1.0_rp - 2.0_rp * cpcon * real(ieley) / real(nelau) )  
     !     xauxi = xauxi + ymaux * dycoo(ieley)
     !  end do
     !  ratio = ymaux / xauxi
     !  xauxi = 0.0_rp
     !  do ieley=nelau,neley
     !     dycoo(ieley) = ymaux * dycoo(ieley) * ratio
     !     xauxi = xauxi + dycoo(ieley)
     !     write (20,*) ieley,xauxi
     !  end do
     
     do ieley=1,nelau
        dycoo(ieley)= facto * log(1.0_rp - 2.0_rp * cpcon * real(ieley) / real(nelau) )  
        xauxi = xauxi + ymaux * dycoo(ieley)
     end do
     
     ratio = ymaux / xauxi
     
     xauxi = 0.0_rp
     do ieley=1,nelau
        dycoo(ieley) = ymaux * dycoo(ieley) * ratio
        xauxi = xauxi + dycoo(ieley)
     end do
     
     xauxi = 0.0_rp
     do ieley=nelau+1,nepoy-1
        dycoo(ieley) = dycoo(ieley-nelau)
        xauxi = xauxi + dycoo(ieley)
        write (20,*) ieley,xauxi
     end do
     
     xauxi = 0.0_rp
     jeley=nelau
     do ieley=nelau,1,-1     
        jeley = jeley+1
        dycoo(ieley) = dycoo(jeley)
        xauxi = xauxi + dycoo(ieley)
        write (20,*) ieley,xauxi
     end do


  end if

  close(20)
 
  npoin=   nepox  *  nepoy  *  nepoz
  nelem=   nelex  *  neley  *  nelez
  nboun=   2*nelex*neley + 2*nelez*neley + 2*nelex*nelez
  
  if (kflbo == 0) nboun = 0

  ndime= 3_ip
  nnode= 8_ip
  nnodb= 4_ip

  

  write(6,*) '--|'
  write(6,*) '--| Nodes             = ',npoin
  write(6,*) '--| Elements          = ',nelem
  write(6,*) '--|'
  write(6,*) '--|   X-Elements      = ',nelex
  write(6,*) '--|   Y-Elements      = ',neley
  write(6,*) '--|   Z-Elements      = ',nelez
  write(6,*) '--|'
  write(6,*) '--|   dx              = ',dxref
  write(6,*) '--|   dy              = ',dyref
  write(6,*) '--|   dz              = ',dzref

  write(6,*) '--|'
  write(6,*) '--| Compute connectivities... '

  if (itybi == 1) then
     allocate(lnods(nnode,nelem),stat=istat)    
  end if


  if (itybi == 2) then
     
     if (ityme==1) then   ! 8-node
        
        write(lun_dgeo,*) 'TYPES'
        do ielem=1,nelem
           write (lun_dgeo,*) ielem,ktype
        end do
        write(lun_dgeo,*) 'END_TYPES'
        write(lun_dgeo,*) 'ELEMENTS'
        
        ipoin= 0
        ielem= 0
        do ielex=1,nelex
           do ieley=1,neley
              do ielez=1,nelez
                 ipoin   = ipoin    + 1 
                 ielem   = ielem    + 1
                 lnohe(1)= ipoin
                 lnohe(2)= ipoin + nepoy * nepoz
                 lnohe(3)= ipoin + nepoy * nepoz + nepoz
                 lnohe(4)= ipoin                       + nepoz
                 lnohe(5)= lnohe(1) + 1
                 lnohe(6)= lnohe(2) + 1
                 lnohe(7)= lnohe(3) + 1
                 lnohe(8)= lnohe(4) + 1
                 write(lun_dgeo,200) ielem,(lnohe(inode), inode=1,nnode)
              end do
              ipoin   = ipoin    + 1 
           end do
           ipoin   = ipoin    + nepoz
        end do
        
        write(lun_dgeo,*) 'END_ELEMENTS'

     else if (ityme==2) then   ! 27-node

        ktype = 39
        chtyp = 'HEX27'
        if (ityme==2) then
           nnode= 27_ip
           nnodb= 9_ip     
           if (kflbo == 0) nboun = 0
        end if

        write(lun_dgeo,*) 'TYPES'
        do ielem=1,nelem
           write (lun_dgeo,*) ielem,ktype
        end do
        write(lun_dgeo,*) 'END_TYPES'
        write(lun_dgeo,*) 'ELEMENTS'
        
        ipoin= 1
        ielem= 1
        do ielex=1,nelex
           do ieley=1,neley
              do ielez=1,nelez
               
                 lnohe( 1)= ipoin
                 lnohe( 2)= ipoin + ityme * nepoy * nepoz                
                 lnohe( 3)= ipoin + ityme * nepoy * nepoz + ityme * nepoz
                 lnohe( 4)= ipoin + ityme * nepoz

                 lnohe( 5) = lnohe(1) + ityme
                 lnohe( 6) = lnohe(2) + ityme
                 lnohe( 7) = lnohe(3) + ityme
                 lnohe( 8) = lnohe(4) + ityme

                 lnohe( 9) = ipoin +   nepoy*nepoz
                 lnohe(10) = ipoin + ityme*nepoy*nepoz +       nepoz               
                 lnohe(11) = ipoin +       nepoy*nepoz + ityme*nepoz
                 lnohe(12) = ipoin +                           nepoz

                 lnohe(22) = lnohe( 9) + 1
                 lnohe(23) = lnohe(10) + 1
                 lnohe(24) = lnohe(11) + 1
                 lnohe(25) = lnohe(12) + 1

                 lnohe(17) = lnohe(22) + 1              
                 lnohe(18) = lnohe(23) + 1
                 lnohe(19) = lnohe(24) + 1
                 lnohe(20) = lnohe(25) + 1

                 lnohe(13) = lnohe(1) + 1
                 lnohe(14) = lnohe(2) + 1
                 lnohe(15) = lnohe(3) + 1
                 lnohe(16) = lnohe(4) + 1

                 lnohe(21) = ipoin + nepoy*nepoz + nepoz
                 lnohe(27) = lnohe(21) + 1
                 lnohe(26) = lnohe(27) + 1

                 write(lun_dgeo,301) ielem,(lnohe(inode), inode=1,nnode)

                 ipoin   = ipoin    + ityme 
                 ielem   = ielem    + 1

              end do
         
              ipoin   = ipoin + nepoz + 1
        
           end do
     
           ipoin   = ipoin + nepoz + nepoz*nepoy 
        end do
        
        write(lun_dgeo,*) 'END_ELEMENTS'
        
     end if

     write(6,*) '--|'
     write(6,*) '--| Compute coordinates... '
     
     write(lun_dgeo,*) 'COORDINATES'
     
     ipoin= 0
     copoi(1)= 0.0_rp
     copoi(2)= 0.0_rp
     copoi(3)= 0.0_rp
     
     do iepox=1,nepox
        do iepoy=1,nepoy
           do iepoz=1,nepoz
              ipoin= ipoin+1
              write(lun_dgeo,100) ipoin,copoi(1),copoi(2),copoi(3)
              copoi(3)= copoi(3)+dzcoo(iepoz)
           end do
           copoi(2)= copoi(2)+dycoo(iepoy)
           copoi(3)= 0.0_rp
        end do
        copoi(1)= copoi(1)+dxcoo(iepox)
        copoi(2)= 0.0_rp
        copoi(3)= 0.0_rp
     end do

     write(lun_dgeo,*) 'END_COORDINATES'
     
     
  else if (itybi == 1) then
     
     cfiel='LTYPE'
     write(lun_dbin) cfiel(1:20),one
     write(lun_dbin) (ktype, ielem=1,nelem)

!!     write(lun_scre,*) cfiel(1:20),one
!!     write(lun_scre,*) (ktype, ielem=1,nelem)

     ipoin= 0
     ielem= 0
     do ielex=1,nelex
        do ieley=1,neley
           do ielez=1,nelez
              ipoin   = ipoin    + 1 
              ielem   = ielem    + 1
              lnods(1,ielem)= ipoin
              lnods(2,ielem)= ipoin + nepoy*nepoz
              lnods(3,ielem)= ipoin + nepoy*nepoz + nepoz
              lnods(4,ielem)= ipoin                       + nepoz
              lnods(5,ielem)= lnods(1,ielem) + 1
              lnods(6,ielem)= lnods(2,ielem) + 1
              lnods(7,ielem)= lnods(3,ielem) + 1
              lnods(8,ielem)= lnods(4,ielem) + 1
           end do
           ipoin   = ipoin    + 1 
        end do
        ipoin   = ipoin    + nepoz 
     end do

     cfiel='LNODS'
     write(lun_dbin) cfiel(1:20),one
     write(lun_dbin) ((lnods(inode,ielem),inode=1,nnode),ielem=1,nelem)     

!!     write(lun_scre,*) cfiel(1:20),one
!!     write(lun_scre,*) ((lnods(inode,ielem),inode=1,nnode),ielem=1,nelem)     

     deallocate(lnods)
     allocate(coord(3,npoin),stat=istat)

     ipoin= 0
     copoi(1)= 0.0_rp
     copoi(2)= 0.0_rp
     copoi(3)= 0.0_rp
    
     do ielex=1,nelex+1
        do ieley=1,neley+1
           do ielez=1,nelez+1
              ipoin= ipoin+1
              coord(1,ipoin)= copoi(1)
              coord(2,ipoin)= copoi(2)
              coord(3,ipoin)= copoi(3)
              copoi(3)= copoi(3)+dzcoo(ielez)
           end do
           copoi(2)= copoi(2)+dycoo(ieley)
           copoi(3)= 0.0_rp
        end do
        copoi(1)= copoi(1)+dxcoo(ielex)
        copoi(2)= 0.0_rp
        copoi(3)= 0.0_rp
     end do

     cfiel='COORD'
     write(lun_dbin) cfiel(1:20),one
     write(lun_dbin) ((coord(idime,ipoin),idime=1,ndime),ipoin=1,npoin)

!!     write(lun_scre,*) cfiel(1:20),one
!!     write(lun_scre,*) ((coord(idime,ipoin),idime=1,ndime),ipoin=1,npoin)
     
  end if

  write(6,*) '--|'
  write(6,*) '--| Compute boundary elements... '

  
  if (itybi == 2) then
     
     
     if (kflbo == 1) then

         write(lun_dgeo,*) 'BOUNDARIES , ELEMENTS'
        
         ielem= 0
         iflbo= ''
         iboun= 0

        if (ityme==1) then   ! 8-node
          ! write(lun_dgeo,*) 'BOUNDARIES , ELEMENTS'
           ipoin= 0
           do ielex=1,nelex
              do ieley=1,neley
                 do ielez=1,nelez
                    ipoin   = ipoin    + 1 
                    ielem   = ielem    + 1
                    lnohe(1)= ipoin
                    lnohe(2)= ipoin + nepoy*nepoz
                    lnohe(3)= ipoin + nepoy*nepoz + nepoz
                    lnohe(4)= ipoin                       + nepoz
                    lnohe(5)= lnohe(1) + 1
                    lnohe(6)= lnohe(2) + 1
                    lnohe(7)= lnohe(3) + 1
                    lnohe(8)= lnohe(4) + 1
                    
                    iflbo= ''
                    if (ielez==1) then
                       iflbo='p'              ! piso
                    else if (ielez==nelez) then
                       iflbo='t'              ! techo
                    end if
                    if (ieley==1) then      
                       iflbo='d'//iflbo       ! derecha
                    else if (ieley==neley) then
                       iflbo='i'//iflbo       ! izquierda
                    end if
                    if (ielex==1) then
                       iflbo='e'//iflbo       ! entrada
                    else if (ielex==nelex) then
                       iflbo='s'//iflbo       ! salida
                    end if
                    
                    if (iflbo == '') then
                       iflbo='000'
                    end if
                    ! F1  -> N1,N4,N3,N2 piso
                    ! F2  -> N1,N2,N6,N5 izquierda
                    ! F3  -> N2,N3,N7,N6 salida
                    ! F4  -> N3,N4,N8,N7 derecha
                    ! F5  -> N1,N5,N8,N4 entrada
                    ! F6  -> N5,N6,N7,N8 techo
                    
                    if (iflbo == '000') then
                       ! nodo interior
                    else if (iflbo == 'e') then
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(1),lnohe(5),lnohe(8),lnohe(4),ielem
                    else if (iflbo == 's') then
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(2),lnohe(3),lnohe(7),lnohe(6),ielem
                    else if (iflbo == 'p') then
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(1),lnohe(4),lnohe(3),lnohe(2),ielem
                    else if (iflbo == 't') then
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(5),lnohe(6),lnohe(7),lnohe(8),ielem
                    else if (iflbo == 'i') then
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(1),lnohe(2),lnohe(6),lnohe(5),ielem
                    else if (iflbo == 'd') then
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(3),lnohe(4),lnohe(8),lnohe(7),ielem
                    else if (iflbo == 'dp') then
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(3),lnohe(4),lnohe(8),lnohe(7),ielem
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(1),lnohe(4),lnohe(3),lnohe(2),ielem
                    else if (iflbo == 'dt') then
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(3),lnohe(4),lnohe(8),lnohe(7),ielem
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(5),lnohe(6),lnohe(7),lnohe(8),ielem
                    else if (iflbo == 'ed') then
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(1),lnohe(5),lnohe(8),lnohe(4),ielem
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(3),lnohe(4),lnohe(8),lnohe(7),ielem
                    else if (iflbo == 'ep') then
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(1),lnohe(5),lnohe(8),lnohe(4),ielem
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(1),lnohe(4),lnohe(3),lnohe(2),ielem
                    else if (iflbo == 'et') then
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(1),lnohe(5),lnohe(8),lnohe(4),ielem
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(5),lnohe(6),lnohe(7),lnohe(8),ielem
                    else if (iflbo == 'ei') then
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(1),lnohe(5),lnohe(8),lnohe(4),ielem
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(1),lnohe(2),lnohe(6),lnohe(5),ielem
                    else if (iflbo == 'ip') then
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(1),lnohe(2),lnohe(6),lnohe(5),ielem
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(1),lnohe(4),lnohe(3),lnohe(2),ielem
                    else if (iflbo == 'it') then
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(1),lnohe(2),lnohe(6),lnohe(5),ielem
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(5),lnohe(6),lnohe(7),lnohe(8),ielem
                    else if (iflbo == 'sd') then
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(2),lnohe(3),lnohe(7),lnohe(6),ielem
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(3),lnohe(4),lnohe(8),lnohe(7),ielem
                    else if (iflbo == 'sp') then
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(2),lnohe(3),lnohe(7),lnohe(6),ielem
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(1),lnohe(4),lnohe(3),lnohe(2),ielem
                    else if (iflbo == 'st') then
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(2),lnohe(3),lnohe(7),lnohe(6),ielem
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(5),lnohe(6),lnohe(7),lnohe(8),ielem
                    else if (iflbo == 'si') then
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(2),lnohe(3),lnohe(7),lnohe(6),ielem
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(1),lnohe(2),lnohe(6),lnohe(5),ielem
                    else if (iflbo == 'edp') then
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(1),lnohe(5),lnohe(8),lnohe(4),ielem
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(3),lnohe(4),lnohe(8),lnohe(7),ielem
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(1),lnohe(4),lnohe(3),lnohe(2),ielem
                    else if (iflbo == 'edt') then
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(1),lnohe(5),lnohe(8),lnohe(4),ielem
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(3),lnohe(4),lnohe(8),lnohe(7),ielem
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(5),lnohe(6),lnohe(7),lnohe(8),ielem
                    else if (iflbo == 'eip') then
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(1),lnohe(5),lnohe(8),lnohe(4),ielem
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(1),lnohe(2),lnohe(6),lnohe(5),ielem
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(1),lnohe(4),lnohe(3),lnohe(2),ielem
                    else if (iflbo == 'eit') then
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(1),lnohe(5),lnohe(8),lnohe(4),ielem
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(1),lnohe(2),lnohe(6),lnohe(5),ielem
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(5),lnohe(6),lnohe(7),lnohe(8),ielem
                    else if (iflbo == 'sdp') then
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(2),lnohe(3),lnohe(7),lnohe(6),ielem
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(3),lnohe(4),lnohe(8),lnohe(7),ielem
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(1),lnohe(4),lnohe(3),lnohe(2),ielem
                    else if (iflbo == 'sdt') then
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(2),lnohe(3),lnohe(7),lnohe(6),ielem
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(3),lnohe(4),lnohe(8),lnohe(7),ielem
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(5),lnohe(6),lnohe(7),lnohe(8),ielem
                    else if (iflbo == 'sip') then
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(2),lnohe(3),lnohe(7),lnohe(6),ielem
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(1),lnohe(2),lnohe(6),lnohe(5),ielem
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(1),lnohe(4),lnohe(3),lnohe(2),ielem
                    else if (iflbo == 'sit') then
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(2),lnohe(3),lnohe(7),lnohe(6),ielem
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(1),lnohe(2),lnohe(6),lnohe(5),ielem
                       iboun=iboun+1
                       write (lun_dgeo,200) iboun,lnohe(5),lnohe(6),lnohe(7),lnohe(8),ielem
                    end if
                 end do
                 ipoin   = ipoin    + 1 
              end do
              ipoin   = ipoin    +  nepoz 
           end do


        else if (ityme==2) then   ! 27-node
          ! write(lun_dgeo,*) 'BOUNDARIES , ELEMENTS'
           ipoin= 1
           do ielex=1,nelex
              do ieley=1,neley
                 do ielez=1,nelez

                    ielem   = ielem    + 1
                 
                    lnohe(1)= ipoin
                    lnohe(2)= ipoin + (neley*ityme+1)*(nelez*ityme+1)*ityme
                    lnohe(3)= ipoin + (neley*ityme+1)*(nelez*ityme+1)*ityme &
                         + (nelez*ityme+1)  +  (nelez*ityme+1)                     
                    lnohe(4)= ipoin + ityme*(nelez*ityme+1)

                    lnohe(9)= ipoin + (neley*ityme+1)*(nelez*ityme+1)
                    lnohe(10)= ipoin + (neley*ityme+1)*(nelez*ityme+1)*ityme &
                         + (nelez*ityme+1)                                  
                    lnohe(11)= ipoin + (neley*ityme+1)*(nelez*ityme+1) + ityme*(nelez*ityme+1)
                    lnohe(12)= ipoin + (nelez*ityme+1)

                    lnohe(21)= ipoin + (neley*ityme+1)*(nelez*ityme+1) + (nelez*ityme+1)
                    
                    lnohe(13)=lnohe(1) + 1
                    lnohe(22)=lnohe(9) + 1
                    lnohe(14)=lnohe(2) + 1
                    lnohe(23)=lnohe(10) + 1
                    lnohe(15)=lnohe(3) + 1
                    lnohe(24)=lnohe(11) + 1
                    lnohe(16)=lnohe(4) + 1
                    lnohe(25)=lnohe(12) + 1
                    lnohe(27)=lnohe(21) + 1
                    
                    lnohe(5)=lnohe(13) + 1
                    lnohe(17)=lnohe(22) + 1              
                    lnohe(6)=lnohe(14) + 1
                    lnohe(18)=lnohe(23) + 1
                    lnohe(7)=lnohe(15) + 1
                    lnohe(19)=lnohe(24) + 1
                    lnohe(8)=lnohe(16) + 1
                    lnohe(20)=lnohe(25) + 1
                    lnohe(26)=lnohe(27) + 1                   
                                     
                    ipoin   = ipoin    + ityme
                    
                    iflbo= ''
                    if (ielez==1) then
                       iflbo='p'              ! piso
                    else if (ielez==nelez) then
                       iflbo='t'              ! techo
                    end if
                    if (ieley==1) then      
                       iflbo='d'//iflbo       ! derecha
                    else if (ieley==neley) then
                       iflbo='i'//iflbo       ! izquierda
                    end if
                    if (ielex==1) then
                       iflbo='e'//iflbo       ! entrada
                    else if (ielex==nelex) then
                       iflbo='s'//iflbo       ! salida
                    end if
                    
                    if (iflbo == '') then
                       iflbo='000'
                    end if
                    ! F1  -> N1,N4,N3,N2,N12,N11,N10,N9,N21 piso
                    ! F2  -> N1,N2,N6,N5,N9,N14,N17,N13,N22 derecha
                    ! F3  -> N2,N3,N7,N6,N10,N15,N18,N14,N23 salida
                    ! F4  -> N3,N4,N8,N7,N11,N16,N19,N15,N24 izquierda
                    ! F5  -> N1,N5,N8,N4,N13,N20,N16,N12,N25 entrada
                    ! F6  -> N5,N6,N7,N8,N17,N18,N19,N20,N26 techo

                    if (iflbo == '000') then
                       ! nodo interior
                    else if (iflbo == 'e') then
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(1),lnohe(5),lnohe(8),lnohe(4)
                       write (lun_dgeo,301) lnohe(13),lnohe(20),lnohe(16),lnohe(12),lnohe(25),ielem
                    else if (iflbo == 's') then
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(2),lnohe(3),lnohe(7),lnohe(6)
                       write (lun_dgeo,301) lnohe(10),lnohe(15),lnohe(18),lnohe(14),lnohe(23),ielem
                    else if (iflbo == 'p') then
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(1),lnohe(4),lnohe(3),lnohe(2)
                       write (lun_dgeo,301) lnohe(12),lnohe(11),lnohe(10),lnohe(9),lnohe(21),ielem
                    else if (iflbo == 't') then
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(5),lnohe(6),lnohe(7),lnohe(8)
                       write (lun_dgeo,301) lnohe(17),lnohe(18),lnohe(19),lnohe(20),lnohe(26),ielem
                    else if (iflbo == 'd') then
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(1),lnohe(2),lnohe(6),lnohe(5)
                       write (lun_dgeo,301) lnohe(9),lnohe(14),lnohe(17),lnohe(13),lnohe(22),ielem
                    else if (iflbo == 'i') then
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(3),lnohe(4),lnohe(8),lnohe(7)
                       write (lun_dgeo,301) lnohe(11),lnohe(16),lnohe(19),lnohe(15),lnohe(24),ielem
                    else if (iflbo == 'ip') then
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(3),lnohe(4),lnohe(8),lnohe(7)
                       write (lun_dgeo,301) lnohe(11),lnohe(16),lnohe(19),lnohe(15),lnohe(24),ielem
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(1),lnohe(4),lnohe(3),lnohe(2)
                       write (lun_dgeo,301) lnohe(12),lnohe(11),lnohe(10),lnohe(9),lnohe(21),ielem
                    else if (iflbo == 'it') then
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(3),lnohe(4),lnohe(8),lnohe(7)
                       write (lun_dgeo,301) lnohe(11),lnohe(16),lnohe(19),lnohe(15),lnohe(24),ielem
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(5),lnohe(6),lnohe(7),lnohe(8)
                       write (lun_dgeo,301) lnohe(17),lnohe(18),lnohe(19),lnohe(20),lnohe(26),ielem
                    else if (iflbo == 'ei') then
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(1),lnohe(5),lnohe(8),lnohe(4)
                       write (lun_dgeo,301) lnohe(13),lnohe(20),lnohe(16),lnohe(12),lnohe(25),ielem
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(3),lnohe(4),lnohe(8),lnohe(7)
                       write (lun_dgeo,301) lnohe(11),lnohe(16),lnohe(19),lnohe(15),lnohe(24),ielem
                    else if (iflbo == 'ep') then
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(1),lnohe(5),lnohe(8),lnohe(4)
                       write (lun_dgeo,301) lnohe(13),lnohe(20),lnohe(16),lnohe(12),lnohe(25),ielem
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(1),lnohe(4),lnohe(3),lnohe(2)
                       write (lun_dgeo,301)  lnohe(12),lnohe(11),lnohe(10),lnohe(9),lnohe(21),ielem
                    else if (iflbo == 'et') then
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(1),lnohe(5),lnohe(8),lnohe(4)
                       write (lun_dgeo,301) lnohe(13),lnohe(20),lnohe(16),lnohe(12),lnohe(25),ielem
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(5),lnohe(6),lnohe(7),lnohe(8)
                       write (lun_dgeo,301) lnohe(17),lnohe(18),lnohe(19),lnohe(20),lnohe(26),ielem
                    else if (iflbo == 'ed') then
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(1),lnohe(5),lnohe(8),lnohe(4)
                       write (lun_dgeo,301) lnohe(13),lnohe(20),lnohe(16),lnohe(12),lnohe(25),ielem
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(1),lnohe(2),lnohe(6),lnohe(5)
                       write (lun_dgeo,301) lnohe(9),lnohe(14),lnohe(17),lnohe(13),lnohe(22),ielem
                    else if (iflbo == 'dp') then
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(1),lnohe(2),lnohe(6),lnohe(5)
                       write (lun_dgeo,301) lnohe(9),lnohe(14),lnohe(17),lnohe(13),lnohe(22),ielem
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(1),lnohe(4),lnohe(3),lnohe(2)
                       write (lun_dgeo,301) lnohe(12),lnohe(11),lnohe(10),lnohe(9),lnohe(21),ielem
                    else if (iflbo == 'dt') then
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(1),lnohe(2),lnohe(6),lnohe(5)
                       write (lun_dgeo,301) lnohe(9),lnohe(14),lnohe(17),lnohe(13),lnohe(22),ielem
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(5),lnohe(6),lnohe(7),lnohe(8)
                       write (lun_dgeo,301) lnohe(17),lnohe(18),lnohe(19),lnohe(20),lnohe(26),ielem
                    else if (iflbo == 'si') then
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(2),lnohe(3),lnohe(7),lnohe(6)
                       write (lun_dgeo,301) lnohe(10),lnohe(15),lnohe(18),lnohe(14),lnohe(23),ielem
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(3),lnohe(4),lnohe(8),lnohe(7)
                       write (lun_dgeo,301) lnohe(11),lnohe(16),lnohe(19),lnohe(15),lnohe(24),ielem
                    else if (iflbo == 'sp') then
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(2),lnohe(3),lnohe(7),lnohe(6)
                       write (lun_dgeo,301) lnohe(10),lnohe(15),lnohe(18),lnohe(14),lnohe(23),ielem
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(1),lnohe(4),lnohe(3),lnohe(2)
                       write (lun_dgeo,301) lnohe(12),lnohe(11),lnohe(10),lnohe(9),lnohe(21),ielem
                    else if (iflbo == 'st') then
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(2),lnohe(3),lnohe(7),lnohe(6)
                       write (lun_dgeo,301) lnohe(10),lnohe(15),lnohe(18),lnohe(14),lnohe(23),ielem
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(5),lnohe(6),lnohe(7),lnohe(8)
                       write (lun_dgeo,301) lnohe(17),lnohe(18),lnohe(19),lnohe(20),lnohe(26),ielem
                    else if (iflbo == 'sd') then
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(2),lnohe(3),lnohe(7),lnohe(6)
                       write (lun_dgeo,301) lnohe(10),lnohe(15),lnohe(18),lnohe(14),lnohe(23),ielem
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(1),lnohe(2),lnohe(6),lnohe(5)
                       write (lun_dgeo,301) lnohe(9),lnohe(14),lnohe(17),lnohe(13),lnohe(22),ielem
                    else if (iflbo == 'eip') then
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(1),lnohe(5),lnohe(8),lnohe(4)
                       write (lun_dgeo,301) lnohe(13),lnohe(20),lnohe(16),lnohe(12),lnohe(25),ielem
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(3),lnohe(4),lnohe(8),lnohe(7)
                       write (lun_dgeo,301) lnohe(11),lnohe(16),lnohe(19),lnohe(15),lnohe(24),ielem
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(1),lnohe(4),lnohe(3),lnohe(2)
                       write (lun_dgeo,301) lnohe(12),lnohe(11),lnohe(10),lnohe(9),lnohe(21),ielem
                    else if (iflbo == 'eit') then
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(1),lnohe(5),lnohe(8),lnohe(4)
                       write (lun_dgeo,301) lnohe(13),lnohe(20),lnohe(16),lnohe(12),lnohe(25),ielem  
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(3),lnohe(4),lnohe(8),lnohe(7)
                       write (lun_dgeo,301) lnohe(11),lnohe(16),lnohe(19),lnohe(15),lnohe(24),ielem
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(5),lnohe(6),lnohe(7),lnohe(8)
                       write (lun_dgeo,301) lnohe(17),lnohe(18),lnohe(19),lnohe(20),lnohe(26),ielem
                    else if (iflbo == 'edp') then
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(1),lnohe(5),lnohe(8),lnohe(4)
                       write (lun_dgeo,301) lnohe(13),lnohe(20),lnohe(16),lnohe(12),lnohe(25),ielem
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(1),lnohe(2),lnohe(6),lnohe(5)
                       write (lun_dgeo,301) lnohe(9),lnohe(14),lnohe(17),lnohe(13),lnohe(22),ielem
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(1),lnohe(4),lnohe(3),lnohe(2)
                       write (lun_dgeo,301) lnohe(12),lnohe(11),lnohe(10),lnohe(9),lnohe(21),ielem
                    else if (iflbo == 'edt') then
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(1),lnohe(5),lnohe(8),lnohe(4)
                       write (lun_dgeo,301) lnohe(13),lnohe(20),lnohe(16),lnohe(12),lnohe(25),ielem
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(1),lnohe(2),lnohe(6),lnohe(5)
                       write (lun_dgeo,301) lnohe(9),lnohe(14),lnohe(17),lnohe(13),lnohe(22),ielem
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(5),lnohe(6),lnohe(7),lnohe(8)
                       write (lun_dgeo,301) lnohe(17),lnohe(18),lnohe(19),lnohe(20),lnohe(26),ielem
                    else if (iflbo == 'sip') then
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(2),lnohe(3),lnohe(7),lnohe(6)
                       write (lun_dgeo,301) lnohe(10),lnohe(15),lnohe(18),lnohe(14),lnohe(23),ielem
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(3),lnohe(4),lnohe(8),lnohe(7)
                       write (lun_dgeo,301) lnohe(11),lnohe(16),lnohe(19),lnohe(15),lnohe(24),ielem
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(1),lnohe(4),lnohe(3),lnohe(2)
                       write (lun_dgeo,301) lnohe(12),lnohe(11),lnohe(10),lnohe(9),lnohe(21),ielem
                    else if (iflbo == 'sit') then
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(2),lnohe(3),lnohe(7),lnohe(6)
                       write (lun_dgeo,301) lnohe(10),lnohe(15),lnohe(18),lnohe(14),lnohe(23),ielem
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(3),lnohe(4),lnohe(8),lnohe(7)
                       write (lun_dgeo,301) lnohe(11),lnohe(16),lnohe(19),lnohe(15),lnohe(24),ielem
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(5),lnohe(6),lnohe(7),lnohe(8)
                       write (lun_dgeo,301) lnohe(17),lnohe(18),lnohe(19),lnohe(20),lnohe(26),ielem
                    else if (iflbo == 'sdp') then
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(2),lnohe(3),lnohe(7),lnohe(6)
                       write (lun_dgeo,301) lnohe(10),lnohe(15),lnohe(18),lnohe(14),lnohe(23),ielem
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(1),lnohe(2),lnohe(6),lnohe(5)
                       write (lun_dgeo,301) lnohe(9),lnohe(14),lnohe(17),lnohe(13),lnohe(22),ielem
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(1),lnohe(4),lnohe(3),lnohe(2)
                       write (lun_dgeo,301) lnohe(12),lnohe(11),lnohe(10),lnohe(9),lnohe(21),ielem
                    else if (iflbo == 'sdt') then
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(2),lnohe(3),lnohe(7),lnohe(6)
                       write (lun_dgeo,301) lnohe(10),lnohe(15),lnohe(18),lnohe(14),lnohe(23),ielem
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(1),lnohe(2),lnohe(6),lnohe(5)
                       write (lun_dgeo,301) lnohe(9),lnohe(14),lnohe(17),lnohe(13),lnohe(22),ielem
                       iboun=iboun+1
                       write (lun_dgeo,300) iboun,lnohe(5),lnohe(6),lnohe(7),lnohe(8)
                       write (lun_dgeo,301) lnohe(17),lnohe(18),lnohe(19),lnohe(20),lnohe(26),ielem
                    end if
                    
                    
                 end do                
                 ipoin   = ipoin +  (nelez+1)*ityme                 
              end do              
              ipoin   = ipoin + 1   + nelez*ityme + (nelez*ityme+1)*(neley*ityme+1)
           end do
           

        end if
          
        
        if (nboun .ne. iboun) then
           write(6,300) nboun,iboun
           stop
        end if
        
        write(6,*) '--|'
        write(6,*) '--| Boundary elements = ',nboun
        write(6,*) '--|'
        
        write(lun_dgeo,*) 'END_BOUNDARIES'
        write(lun_dgeo,*) 'SKEW_SYSTEMS'
        write(lun_dgeo,*) 'END_SKEW_SYSTEMS'
        
     end if
     
  else if (itybi ==1) then
     
     allocate(lnodb(5,nboun),stat=istat)    
     
     if (kflbo == 1) then
        
        ipoin= 0
        ielem= 0
        iflbo= ''
        iboun= 0
        do ielex=1,nelex
           do ieley=1,neley
              do ielez=1,nelez
                 ipoin   = ipoin    + 1 
                 ielem   = ielem    + 1
                 lnohe(1)= ipoin
                 lnohe(2)= ipoin + nepoy*nepoz
                 lnohe(3)= ipoin + nepoy*nepoz + nepoz
                 lnohe(4)= ipoin                       + nepoz
                 lnohe(5)= lnohe(1) + 1
                 lnohe(6)= lnohe(2) + 1
                 lnohe(7)= lnohe(3) + 1
                 lnohe(8)= lnohe(4) + 1
                 
                 iflbo= ''
                 if (ielez==1) then
                    iflbo='p'              ! piso
                 else if (ielez==nelez) then
                    iflbo='t'              ! techo
                 end if
                 if (ieley==1) then      
                    iflbo='d'//iflbo       ! derecha
                 else if (ieley==neley) then
                    iflbo='i'//iflbo       ! izquierda
                 end if
                 if (ielex==1) then
                    iflbo='e'//iflbo       ! entrada
                 else if (ielex==nelex) then
                    iflbo='s'//iflbo       ! salida
                 end if
                 
                 if (iflbo == '') then
                    iflbo='000'
                 end if
                 ! F1  -> N1,N4,N3,N2 piso
                 ! F2  -> N1,N2,N6,N5 izquierda
                 ! F3  -> N2,N3,N7,N6 salida
                 ! F4  -> N3,N4,N8,N7 derecha
                 ! F5  -> N1,N5,N8,N4 entrada
                 ! F6  -> N5,N6,N7,N8 techo
                 
                 if (iflbo == '000') then
                    ! nodo interior
                 else if (iflbo == 'e') then
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(1)
                    lnodb(2,iboun) = lnohe(5)
                    lnodb(3,iboun) = lnohe(8)
                    lnodb(4,iboun) = lnohe(4)
                    lnodb(5,iboun) = ielem
                 else if (iflbo == 's') then
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(2)
                    lnodb(2,iboun) = lnohe(3)
                    lnodb(3,iboun) = lnohe(7)
                    lnodb(4,iboun) = lnohe(6)
                    lnodb(5,iboun) = ielem
                 else if (iflbo == 'p') then
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(1)
                    lnodb(2,iboun) = lnohe(4)
                    lnodb(3,iboun) = lnohe(3)
                    lnodb(4,iboun) = lnohe(2)
                    lnodb(5,iboun) = ielem
                 else if (iflbo == 't') then
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(5)
                    lnodb(2,iboun) = lnohe(6)
                    lnodb(3,iboun) = lnohe(7)
                    lnodb(4,iboun) = lnohe(8)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(5),lnohe(6),lnohe(7),lnohe(8),ielem
                 else if (iflbo == 'i') then
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(1)
                    lnodb(2,iboun) = lnohe(2)
                    lnodb(3,iboun) = lnohe(6)
                    lnodb(4,iboun) = lnohe(5)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(1),lnohe(2),lnohe(6),lnohe(5),ielem
                 else if (iflbo == 'd') then
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(3)
                    lnodb(2,iboun) = lnohe(4)
                    lnodb(3,iboun) = lnohe(8)
                    lnodb(4,iboun) = lnohe(7)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(3),lnohe(4),lnohe(8),lnohe(7),ielem
                 else if (iflbo == 'dp') then
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(3)
                    lnodb(2,iboun) = lnohe(4)
                    lnodb(3,iboun) = lnohe(8)
                    lnodb(4,iboun) = lnohe(7)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(3),lnohe(4),lnohe(8),lnohe(7),ielem
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(1)
                    lnodb(2,iboun) = lnohe(4)
                    lnodb(3,iboun) = lnohe(3)
                    lnodb(4,iboun) = lnohe(2)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(1),lnohe(4),lnohe(3),lnohe(2),ielem
                 else if (iflbo == 'dt') then
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(3)
                    lnodb(2,iboun) = lnohe(4)
                    lnodb(3,iboun) = lnohe(8)
                    lnodb(4,iboun) = lnohe(7)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(3),lnohe(4),lnohe(8),lnohe(7),ielem
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(5)
                    lnodb(2,iboun) = lnohe(6)
                    lnodb(3,iboun) = lnohe(7)
                    lnodb(4,iboun) = lnohe(8)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(5),lnohe(6),lnohe(7),lnohe(8),ielem
                 else if (iflbo == 'ed') then
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(1)
                    lnodb(2,iboun) = lnohe(5)
                    lnodb(3,iboun) = lnohe(8)
                    lnodb(4,iboun) = lnohe(4)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(1),lnohe(5),lnohe(8),lnohe(4),ielem
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(3)
                    lnodb(2,iboun) = lnohe(4)
                    lnodb(3,iboun) = lnohe(8)
                    lnodb(4,iboun) = lnohe(7)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(3),lnohe(4),lnohe(8),lnohe(7),ielem
                 else if (iflbo == 'ep') then
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(1)
                    lnodb(2,iboun) = lnohe(5)
                    lnodb(3,iboun) = lnohe(8)
                    lnodb(4,iboun) = lnohe(4)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(1),lnohe(5),lnohe(8),lnohe(4),ielem
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(1)
                    lnodb(2,iboun) = lnohe(4)
                    lnodb(3,iboun) = lnohe(3)
                    lnodb(4,iboun) = lnohe(2)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(1),lnohe(4),lnohe(3),lnohe(2),ielem
                 else if (iflbo == 'et') then
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(1)
                    lnodb(2,iboun) = lnohe(5)
                    lnodb(3,iboun) = lnohe(8)
                    lnodb(4,iboun) = lnohe(4)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(1),lnohe(5),lnohe(8),lnohe(4),ielem
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(5)
                    lnodb(2,iboun) = lnohe(6)
                    lnodb(3,iboun) = lnohe(7)
                    lnodb(4,iboun) = lnohe(8)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(5),lnohe(6),lnohe(7),lnohe(8),ielem
                 else if (iflbo == 'ei') then
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(1)
                    lnodb(2,iboun) = lnohe(5)
                    lnodb(3,iboun) = lnohe(8)
                    lnodb(4,iboun) = lnohe(4)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(1),lnohe(5),lnohe(8),lnohe(4),ielem
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(1)
                    lnodb(2,iboun) = lnohe(2)
                    lnodb(3,iboun) = lnohe(6)
                    lnodb(4,iboun) = lnohe(5)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(1),lnohe(2),lnohe(6),lnohe(5),ielem
                 else if (iflbo == 'ip') then
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(1)
                    lnodb(2,iboun) = lnohe(2)
                    lnodb(3,iboun) = lnohe(6)
                    lnodb(4,iboun) = lnohe(5)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(1),lnohe(2),lnohe(6),lnohe(5),ielem
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(1)
                    lnodb(2,iboun) = lnohe(4)
                    lnodb(3,iboun) = lnohe(3)
                    lnodb(4,iboun) = lnohe(2)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(1),lnohe(4),lnohe(3),lnohe(2),ielem
                 else if (iflbo == 'it') then
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(1)
                    lnodb(2,iboun) = lnohe(2)
                    lnodb(3,iboun) = lnohe(6)
                    lnodb(4,iboun) = lnohe(5)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(1),lnohe(2),lnohe(6),lnohe(5),ielem
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(5)
                    lnodb(2,iboun) = lnohe(6)
                    lnodb(3,iboun) = lnohe(7)
                    lnodb(4,iboun) = lnohe(8)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(5),lnohe(6),lnohe(7),lnohe(8),ielem
                 else if (iflbo == 'sd') then
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(2)
                    lnodb(2,iboun) = lnohe(3)
                    lnodb(3,iboun) = lnohe(7)
                    lnodb(4,iboun) = lnohe(6)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(2),lnohe(3),lnohe(7),lnohe(6),ielem
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(3)
                    lnodb(2,iboun) = lnohe(4)
                    lnodb(3,iboun) = lnohe(8)
                    lnodb(4,iboun) = lnohe(7)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(3),lnohe(4),lnohe(8),lnohe(7),ielem
                 else if (iflbo == 'sp') then
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(2)
                    lnodb(2,iboun) = lnohe(3)
                    lnodb(3,iboun) = lnohe(7)
                    lnodb(4,iboun) = lnohe(6)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(2),lnohe(3),lnohe(7),lnohe(6),ielem
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(1)
                    lnodb(2,iboun) = lnohe(4)
                    lnodb(3,iboun) = lnohe(3)
                    lnodb(4,iboun) = lnohe(2)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(1),lnohe(4),lnohe(3),lnohe(2),ielem
                 else if (iflbo == 'st') then
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(2)
                    lnodb(2,iboun) = lnohe(3)
                    lnodb(3,iboun) = lnohe(7)
                    lnodb(4,iboun) = lnohe(6)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(2),lnohe(3),lnohe(7),lnohe(6),ielem
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(5)
                    lnodb(2,iboun) = lnohe(6)
                    lnodb(3,iboun) = lnohe(7)
                    lnodb(4,iboun) = lnohe(8)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(5),lnohe(6),lnohe(7),lnohe(8),ielem
                 else if (iflbo == 'si') then
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(2)
                    lnodb(2,iboun) = lnohe(3)
                    lnodb(3,iboun) = lnohe(7)
                    lnodb(4,iboun) = lnohe(6)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(2),lnohe(3),lnohe(7),lnohe(6),ielem
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(1)
                    lnodb(2,iboun) = lnohe(2)
                    lnodb(3,iboun) = lnohe(6)
                    lnodb(4,iboun) = lnohe(5)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(1),lnohe(2),lnohe(6),lnohe(5),ielem
                 else if (iflbo == 'edp') then
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(1)
                    lnodb(2,iboun) = lnohe(5)
                    lnodb(3,iboun) = lnohe(8)
                    lnodb(4,iboun) = lnohe(4)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(1),lnohe(5),lnohe(8),lnohe(4),ielem
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(3)
                    lnodb(2,iboun) = lnohe(4)
                    lnodb(3,iboun) = lnohe(8)
                    lnodb(4,iboun) = lnohe(7)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(3),lnohe(4),lnohe(8),lnohe(7),ielem
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(1)
                    lnodb(2,iboun) = lnohe(4)
                    lnodb(3,iboun) = lnohe(3)
                    lnodb(4,iboun) = lnohe(2)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(1),lnohe(4),lnohe(3),lnohe(2),ielem
                 else if (iflbo == 'edt') then
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(1)
                    lnodb(2,iboun) = lnohe(5)
                    lnodb(3,iboun) = lnohe(8)
                    lnodb(4,iboun) = lnohe(4)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(1),lnohe(5),lnohe(8),lnohe(4),ielem
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(3)
                    lnodb(2,iboun) = lnohe(4)
                    lnodb(3,iboun) = lnohe(8)
                    lnodb(4,iboun) = lnohe(7)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(3),lnohe(4),lnohe(8),lnohe(7),ielem
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(5)
                    lnodb(2,iboun) = lnohe(6)
                    lnodb(3,iboun) = lnohe(7)
                    lnodb(4,iboun) = lnohe(8)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(5),lnohe(6),lnohe(7),lnohe(8),ielem
                 else if (iflbo == 'eip') then
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(1)
                    lnodb(2,iboun) = lnohe(5)
                    lnodb(3,iboun) = lnohe(8)
                    lnodb(4,iboun) = lnohe(4)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(1),lnohe(5),lnohe(8),lnohe(4),ielem
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(1)
                    lnodb(2,iboun) = lnohe(2)
                    lnodb(3,iboun) = lnohe(6)
                    lnodb(4,iboun) = lnohe(5)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(1),lnohe(2),lnohe(6),lnohe(5),ielem
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(1)
                    lnodb(2,iboun) = lnohe(4)
                    lnodb(3,iboun) = lnohe(3)
                    lnodb(4,iboun) = lnohe(2)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(1),lnohe(4),lnohe(3),lnohe(2),ielem
                 else if (iflbo == 'eit') then
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(1)
                    lnodb(2,iboun) = lnohe(5)
                    lnodb(3,iboun) = lnohe(8)
                    lnodb(4,iboun) = lnohe(4)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(1),lnohe(5),lnohe(8),lnohe(4),ielem
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(1)
                    lnodb(2,iboun) = lnohe(2)
                    lnodb(3,iboun) = lnohe(6)
                    lnodb(4,iboun) = lnohe(5)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(1),lnohe(2),lnohe(6),lnohe(5),ielem
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(5)
                    lnodb(2,iboun) = lnohe(6)
                    lnodb(3,iboun) = lnohe(7)
                    lnodb(4,iboun) = lnohe(8)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(5),lnohe(6),lnohe(7),lnohe(8),ielem
                 else if (iflbo == 'sdp') then
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(2)
                    lnodb(2,iboun) = lnohe(3)
                    lnodb(3,iboun) = lnohe(7)
                    lnodb(4,iboun) = lnohe(6)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(2),lnohe(3),lnohe(7),lnohe(6),ielem
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(3)
                    lnodb(2,iboun) = lnohe(4)
                    lnodb(3,iboun) = lnohe(8)
                    lnodb(4,iboun) = lnohe(7)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(3),lnohe(4),lnohe(8),lnohe(7),ielem
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(1)
                    lnodb(2,iboun) = lnohe(4)
                    lnodb(3,iboun) = lnohe(3)
                    lnodb(4,iboun) = lnohe(2)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(1),lnohe(4),lnohe(3),lnohe(2),ielem
                 else if (iflbo == 'sdt') then
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(2)
                    lnodb(2,iboun) = lnohe(3)
                    lnodb(3,iboun) = lnohe(7)
                    lnodb(4,iboun) = lnohe(6)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(2),lnohe(3),lnohe(7),lnohe(6),ielem
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(3)
                    lnodb(2,iboun) = lnohe(4)
                    lnodb(3,iboun) = lnohe(8)
                    lnodb(4,iboun) = lnohe(7)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(3),lnohe(4),lnohe(8),lnohe(7),ielem
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(5)
                    lnodb(2,iboun) = lnohe(6)
                    lnodb(3,iboun) = lnohe(7)
                    lnodb(4,iboun) = lnohe(8)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(5),lnohe(6),lnohe(7),lnohe(8),ielem
                 else if (iflbo == 'sip') then
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(2)
                    lnodb(2,iboun) = lnohe(3)
                    lnodb(3,iboun) = lnohe(7)
                    lnodb(4,iboun) = lnohe(6)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(2),lnohe(3),lnohe(7),lnohe(6),ielem
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(1)
                    lnodb(2,iboun) = lnohe(2)
                    lnodb(3,iboun) = lnohe(6)
                    lnodb(4,iboun) = lnohe(5)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(1),lnohe(2),lnohe(6),lnohe(5),ielem
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(1)
                    lnodb(2,iboun) = lnohe(4)
                    lnodb(3,iboun) = lnohe(3)
                    lnodb(4,iboun) = lnohe(2)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(1),lnohe(4),lnohe(3),lnohe(2),ielem
                 else if (iflbo == 'sit') then
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(2)
                    lnodb(2,iboun) = lnohe(3)
                    lnodb(3,iboun) = lnohe(7)
                    lnodb(4,iboun) = lnohe(6)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(2),lnohe(3),lnohe(7),lnohe(6),ielem
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(1)
                    lnodb(2,iboun) = lnohe(2)
                    lnodb(3,iboun) = lnohe(6)
                    lnodb(4,iboun) = lnohe(5)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(1),lnohe(2),lnohe(6),lnohe(5),ielem
                    iboun=iboun+1
                    lnodb(1,iboun) = lnohe(5)
                    lnodb(2,iboun) = lnohe(6)
                    lnodb(3,iboun) = lnohe(7)
                    lnodb(4,iboun) = lnohe(8)
                    lnodb(5,iboun) = ielem
!                    write (lun_dgeo,200) iboun,lnohe(5),lnohe(6),lnohe(7),lnohe(8),ielem
                 end if
              end do
              ipoin   = ipoin    + 1 
           end do
           ipoin   = ipoin    +  nepoz 
        end do
        
        write(6,*) '--|'
        write(6,*) '--| Boundary elements = ',nboun
        write(6,*) '--|'

        cfiel='LTYPB'
        ktype= 7
        write(lun_dbin) cfiel(1:20),one
        write(lun_dbin) (ktype,iboun=1,nboun)

!!        write(lun_scre,*) cfiel(1:20),one
!!        write(lun_scre,*) (ktype,iboun=1,nboun)

        cfiel='LNODB'
        write(lun_dbin) cfiel(1:20),one
        write(lun_dbin) ((lnodb(inodb,iboun),inodb=1,nnodb),iboun=1,nboun)

!!        write(lun_scre,*) cfiel(1:20),one
!!        write(lun_scre,*) ((lnodb(inodb,iboun),inodb=1,nnodb),iboun=1,nboun)

        cfiel='LBOEL'
        write(lun_dbin) cfiel(1:20),one
        write(lun_dbin) (lnodb(5,iboun),iboun=1,nboun)

!!        write(lun_scre,*) cfiel(1:20),one
!!        write(lun_scre,*) (lnodb(5,iboun),iboun=1,nboun)

     end if

     cfiel='END_FILE'
     write(lun_dbin) cfiel(1:20),one            
     close(lun_dbin)
        
     
  end if
  
  write(10,*) '$--------- GENERATED BY CARTESIAN EXPRESS ----------------'
  write(10,*) 'DIMENSIONS'
  write(10,*) '  NODAL_POINTS=      ',npoin
  write(10,*) '  ELEMENTS=          ',nelem
  write(10,*) '  SPACE_DIMENSIONS=  ',ndime
  write(10,*) '  TYPES_OF_ELEMENTS= ',chtyp
  write(10,*) '  NODES=             ',nnode
  write(10,*) '  BOUNDARIES=        ',nboun
  write(10,*) '  SKEW_SYSTEMS=      ',nzero
  write(10,*) '  SLAVES=            ',nzero
  write(10,*) 'END_DIMENSIONS'

  write(10,*) '$---------------------------------------------------------'
  write(10,*) 'STRATEGY'
  write(10,*) '  DOMAIN_INTEGRATION_POINTS: ',nnode
  write(10,*) '  OUTPUT_MESH_DATA:           YES'
  write(10,*) 'END_STRATEGY'

  write(10,*) '$---------------------------------------------------------'
  write(10,*) 'GEOMETRY, GID'
  if (itybi == 2) then
     write(10,*) '  INCLUDE cart-express.dom.geo'
  else
     write(10,*) '  BINARY'
  end if
  write(10,*) 'END_GEOMETRY'
  write(10,*) '$---------------------------------------------------------'
  write(10,*) 'SETS'
  write(10,*) 'END_SETS'
  write(10,*) '$---------------------------------------------------------'
  write(10,*) 'BOUNDARY_CONDITIONS'
  write(10,*) 'END_BOUNDARY_CONDITIONS'
  write(10,*) '$---------------------------------------------------------'

  if (kflfi == 1) then
     
     open(12,file='cart-express.nsa.fix',status='unknown')

     write(6,*) '--|'
     write(6,*) '--| Compute fixities... '
     
     write(12,*) 'ON_NODES, CODED'
     
     
     ipoin= 0
     copoi(1)= 0.0_rp
     copoi(2)= 0.0_rp
     copoi(3)= 0.0_rp
     
     nboco = 0
     
     kboty = 1   ! default, finite carter's plate
     if (kboty == 1) then
        do ielex=1,nepox
           do ieley=1,nepoy
              do ielez=1,nepoz
                 iboco= 0
                 ipoin= ipoin+1
                 if (ielez == 1) then
                    if ( copoi(1)>cabox(1,1) .and. copoi(1)<cabox(2,1) ) then
                       if ( copoi(2)>cabox(1,2) .and. copoi(2)<cabox(2,2) ) then
                          iboco = 3                  ! the plate
                          write(12,*) ipoin,iboco
                          nboco(1)= nboco(1)+1
                          nboco(2)= nboco(2)+1
                       end if
                    end if
                    if (iboco==0) then 
                       iboco = 5                     ! the floor
                       write(12,*) ipoin,iboco
                       nboco(1)= nboco(1)+1
                       nboco(3)= nboco(3)+1
                    end if
                 else if (ielez == nepoz) then                 
                    if (iboco==0) then 
                       iboco = 8                     ! the top
                       write(12,*) ipoin,iboco
                       nboco(1)= nboco(1)+1
                       nboco(4)= nboco(4)+1
                    end if
                 end if
                 if (ieley==1) then
                    iboco= 6                            ! right wall
                    write(12,*) ipoin,iboco
                    nboco(1)= nboco(1)+1
                    nboco(5)= nboco(5)+1
                 else if (ieley==nepoy) then
                    iboco= 7                            ! left wall
                    write(12,*) ipoin,iboco
                    nboco(1)= nboco(1)+1
                    nboco(6)= nboco(6)+1
                 end if
                 if (ielex==1) then
                    iboco= 1                            ! inflow
                    write(12,*) ipoin,iboco
                    nboco(1)= nboco(1)+1
                    nboco(7)= nboco(7)+1
                 else if (ielex==nepox) then
                    iboco= 0                            ! outflow
                    write(12,*) ipoin,iboco
                    nboco(1)= nboco(1)+1
                    nboco(8)= nboco(8)+1
                 end if
                 copoi(3)= copoi(3)+dzcoo(ielez)
              end do
              copoi(2)= copoi(2)+dycoo(ieley)
              copoi(3)= 0.0_rp
           end do
           copoi(1)= copoi(1)+dxcoo(ielex)
           copoi(2)= 0.0_rp
           copoi(3)= 0.0_rp
        end do
        
        write(6,*) '--|'
        write(6,*) '--| Boundary conditions:'
        write(6,*) '--|'
        write(6,*) '--|   Total fixed nodes   = ',nboco(1)
        write(6,*) '--|      Outlet           = ',nboco(8),'  (Code 0)'
        write(6,*) '--|      Inlet            = ',nboco(7),'  (Code 1)'
        write(6,*) '--|      Plate            = ',nboco(2),'  (Code 3)'
        write(6,*) '--|      Floor            = ',nboco(3),'  (Code 5)'
        write(6,*) '--|      Top              = ',nboco(4),'  (Code 8)'
        write(6,*) '--|      Right            = ',nboco(5),'  (Code 6)'
        write(6,*) '--|      Left             = ',nboco(6),'  (Code 7)'
        write(6,*) '--|'
        write(6,*) '--| Codes used: '
        write(6,*)  ' '
        write(6,*)  '  CODES '
        write(6,*)  '       0      00000     1.0     0.0   0.0    1.0   1.0'
        write(6,*)  '       1      11111     1.0     0.0   0.0    1.0   1.0'
        write(6,*)  '       3      11103     0.0     0.0   0.0    1.0   1.0'
        write(6,*)  '       5      00100     1.0     0.0   0.0    1.0   1.0'
        write(6,*)  '       8      11111     1.0     0.0   0.0    1.0   1.0'
        write(6,*)  '       6      11111     1.0     0.0   0.0    1.0   1.0'
        write(6,*)  '       7      11111     1.0     0.0   0.0    1.0   1.0'
        write(6,*)  '  END_CODES '
        write(6,*)  ' '
     
        
     end if
     
     
     
     write(12,*) 'END_ON_NODES'


  else ! no boundary conditions computed

     write(6,*) '--|'
     write(6,*) '--| No boundary conditions computed'
     write(6,*) '--| Only mesh and dat files written'
     write(6,*) '--|'

  end if

  close(10)
  close(lun_dgeo)
  close(12)

  !------------------
  !--    say good bye
  !------------------
  write(6,*) '--|'
  write(6,*) '--| Bye!'
  write(6,*) '--|'

100 format(2x,i12,2x,f17.7,f17.7,f17.7)
200 format(9(2x,i12))
300 format(5(2x,i12),'\\')
301 format(30(2x,i12))
400 format(9(2x,i12))

  
end program expressian


