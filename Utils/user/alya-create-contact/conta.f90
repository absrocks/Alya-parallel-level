program conta

  use def_parame
  use def_elmtyp
  use def_master
  use def_domain , only : lnods,lnodb,coord,ltype,lexis,lnnod,nface,nnodf,lface
  use def_domain , only : npoin,ndime,nelem,nelem_2,nboun,mnodb,mnode 
  use mod_memchk
  use mod_htable
  implicit none

  integer(ip)               :: &
       inodb,ipoin,jpoin,iboun,jboun,ielem,jelem,inode,idime,ielty,iface_local,iface,&
       pnode,idumy,icucu,nzone,izone,jzone,ielem_fluid, ielem_solid,&
       ncnta,nnodc,inodc,ielto,ipnew,npnew,ncont,icont,jcont,nelct,itype,ktype_element,&
       condition_solid,condition_fluid,toca_nodo,ifoun,nmateri,imate
  real(rp)               :: digit,digit2

  integer(ip)               :: lnoco(10,400000), lauxi(100), nelzo(2)
  real(rp)                  :: conew(3,400000)
  integer(ip), pointer      :: &
       lzone_read(:),poaux(:),pzone(:,:),lnodb_element(:),lnodb_putogid(:,:)
  character(120)     :: input_file,output_file
  character(5)       :: wopos
  character(30)      :: tcard,tcard2
  character(30)      :: ccard,ccard2
  character(50)     :: wline

!--------------------
  nzone= 2
  npoin=-1
  nelem= 0
  nboun= 0
  ndime= 2
  mnode= 4
  mnodb= 2
  ncnta= 0
  nnodc= 2
  nmateri= 1
!--------------------


  output_file = 'create-contact.input'
  write (6,*) 'Reading dimensions from:  create-contact.input'
  open (10 ,file=output_file)
  read (10 ,'(a)') input_file
  read (10 ,  *  ) npoin
  read (10 ,  *  ) nelem
  read (10 ,  *  ) nboun
  read (10 ,  *  ) ndime
  read (10 ,  *  ) nmateri
  close(10)

!  call getarg(1,input_file)  
!  call getarg(2,tcard)
!  read(tcard,'(i10)') npoin
!  call getarg(3,tcard)
!  read(tcard,'(i10)') nelem
!  call getarg(4,tcard)
!  read(tcard,'(i10)') nboun
!  call getarg(5,tcard)
!  read(tcard,'(i10)') ndime
!  write(6,*) '   Usage: conta input-file  npoin  nelem  nboun ndime'
  
  if (npoin == -1) then
     stop
  end if

  if (ndime==3) then
     mnodb= 3
     nnodc= 3
  end if

  allocate(coord(ndime,npoin))
  allocate(lnods(mnode  ,nelem))
  allocate(lnodb(mnodb+1,nboun))
  allocate(lnodb_putogid(mnodb+1,nboun))
  allocate(lnodb_element(nboun))
  allocate(lnnod(nelem))
  allocate(lzone_read(nelem))
  allocate(poaux(npoin))
  allocate(pzone(2,npoin))

  lnodb=0
  lnodb_putogid=0
  pzone= 0

  output_file = trim(input_file) // '.geo.dat'
  write (6,*) 'Reading:  ',trim(output_file) 
  open (10 ,file=output_file)
  output_file = trim(input_file) // '-cnt.geo.dat'
  open (606,file=output_file)
  output_file = trim(input_file) // '-zones-cnt.geo.dat'
  open (607,file=output_file)
  output_file = trim(input_file) // '-bocon-cnt.geo.dat'
  open (608,file=output_file)
  output_file = trim(input_file) // '-new.dom.dat'
  open (610,file=output_file)
  if (nmateri > 1) then
     output_file = trim(input_file) // '-materials-cnt.geo.dat'
     open (620,file=output_file)     
  end if
!  output_file = trim(input_file) // '-fibers-cnt.dat'
!  open (636,file=output_file)

  ! Esto es si gid te da las condiciones de contorno en el dom.dat

!  output_file = trim(input_file) // '.dom.dat'
!  open (20 ,file=output_file)
!  icucu= 0
!  do while (icucu == 0)
!     read(20,'(a)') wline
!     tcard= trim(adjustl(wline))
!     if (tcard(1:8) == 'BOUNDARY') then
!        icucu= 1
!     end if
!  end do
!  write(608,*) wline
!  read(20,'(a)') wline
!  write(608,*) wline

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Esto es si gid te da las condiciones de contorno en el fix.dat

  output_file = trim(input_file) // '.fix.dat'
  open (20 ,file=output_file)
  icucu= 0
  do while (icucu == 0)
     read(20,'(a)') wline
     tcard= trim(adjustl(wline))
     if (tcard(1:8) == 'ON_BOUND') then
        icucu= 1
     end if
  end do
  do iboun=1,nboun
     read(20,*) icucu,idumy,lnodb_putogid(1:mnodb+1,iboun)     
  end do
  write(608,*) wline


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ktype_element=0
  read(10,'(a)') wline
  tcard= trim(adjustl(wline))
  if (tcard(1:3) == 'TYP') then
     ktype_element= 1
  end if  
!  POINT =    1 ! 0D
!  BAR02 =    2 ! 1D
!  BAR03 =    3 ! 1D 
!  BAR04 =    4 ! 1D 
!  TRI03 =   10 ! 2D 
!  TRI06 =   11 ! 2D 
!  QUA04 =   12 ! 2D 
!  QUA08 =   13 ! 2D 
!  QUA09 =   14 ! 2D 
!  QUA16 =   15 ! 2D 
!  TET04 =   30 ! 3D 
!  TET10 =   31 ! 3D 
!  PYR05 =   32 ! 3D 
!  PYR14 =   33 ! 3D 
!  PEN06 =   34 ! 3D 
!  PEN15 =   35 ! 3D 
!  PEN18 =   36 ! 3D 
!  HEX08 =   37 ! 3D 
!  HEX20 =   38 ! 3D 
!  HEX27 =   39 ! 3D 
!  HEX64 =   40 ! 3D 
  do ielem=1,nelem
     read(10,*) icucu,lnnod(ielem)

     if (ktype_element == 1) then
        if (lnnod(ielem) == 30) lnnod(ielem)= 4
     else if (ktype_element == 1) then
        if (lnnod(ielem) == 30) lnnod(ielem)= 4
     end if

  end do

  read(10,'(a)') 
  read(10,'(a)') 
  do ielem=1,nelem
     pnode= lnnod(ielem)
     read(10,*) idumy,lnods(1:pnode,ielem)
  end do
  read(10,'(a)') 
  read(10,'(a)') 
  do ipoin=1,npoin
     read(10,*) idumy,coord(1:ndime,ipoin)
  end do
  read(10,'(a)') 
  read(10,'(a)') 
  do iboun=1,nboun
     read(10,*) idumy,lnodb(1:mnodb,iboun),lnodb_element(iboun)
  end do
  read(10,'(a)') 

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! reordenar la mierda de gid
  do iboun=1,nboun
     do jboun=1,nboun
        ifoun= 0
        do inode=1,mnodb
           if (lnodb(inode,iboun) == lnodb_putogid(inode,jboun)) ifoun= ifoun+1
        end do
        if (ifoun == mnodb) exit
     end do
     lnodb(mnodb+1,iboun)= lnodb_putogid(mnodb+1,jboun)
  end do
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (nmateri > 1) then
     ! hay materiales: leerlos
     read(10,'(a)') wline
     write(620,'(a)') wline
     write(6,*) 'Reading/writting materials...  '
     write(6,'(a)') wline
     icucu= 0
     do while (icucu==0)
        read(10,'(a)') ccard
        write(620,'(a)') ccard     
        tcard= trim(adjustl(ccard))
        if (tcard(1:3) == 'END') then
           icucu= 1
        end if
     end do
  end if

  write(6,*) 'Reading/writting zones...  '
  icucu= 0
  do while (icucu == 0)
     read(10,'(a)') wline
     write(607,*) wline
     wopos= trim(adjustl(wline))
     if (wopos == 'ZONES') then
        icucu= 1
     end if
  end do

  lzone_read= 0
  nelzo= 0
  do izone=1,nzone
     read(10,'(a)') wline
     write(607,*) wline     
     icucu= 0
     do while (icucu == 0)
        read(10,'(a)') ccard
        tcard= trim(adjustl(ccard))
        if (tcard(1:3) == 'END') then
           icucu= 1
        else
           read(tcard,*) digit
           ielem = int(digit)
           if (lzone_read(ielem) == 0) then
              lzone_read(ielem) = izone
              nelzo(izone) = nelzo(izone) + 1
              do inode= 1,lnnod(ielem)
                 ipoin= lnods(inode,ielem)                 
                 if (izone==1) then
                    pzone(1,ipoin) = 1
                 end if
              end do
              write(607,*) ccard     
           end if
        end if
     end do
     write(607,*) ccard     
     write(6,*) '   Elements in zone  ',izone, nelzo(izone)

  end do
  ! aqui no cierro la seccion nueva de zonas porque voy a agregar una zona mas
  ! con los contactos

  close (10)

  
  !
  ! bridge a alya que te calcula las faces
  !
  call alyafaces
  !
  !    LFACG(4,NFACG) and LELFA(1:NELEM) % L(1:NFACE(LTYPE(IELEM)))
  !
  !    LFACG(1,IFACG) = IELEM
  !    LFACG(2,IFACG) = JELEM/0 ... 0 if boundary face
  !    LFACG(3,IFACG) = IFACE ..... Local IELEM face
  !    LFACG(4,IFACG) = JFACE/0 ... Local JELEM face
  !    Caution: 0 < IELEM < NELEM_2 so IELEM can be a neighbor's element
  !
  !    LELFA(IELEM) % L(1:NFACE) ... Global face connectivity for IELEM
  !

  ncont= 0
  poaux= 0
  ipnew= 0
  do ielem= 1,nelem
     izone= lzone_read(ielem)
     inodc= 0
     lauxi= 0
     toca_nodo= 0
     
     if (izone == 1) then                    ! solo mirar los de la zona 1
        
        ielty= abs(ltype(ielem))
        do iface_local= 1,nface(ielty)
           iface= lelfa(ielem)%l(iface_local)
           ! ielem cae en la zona 1 (fluido)
           ! jelem es el elemento del otro lado de la cara 
           jelem= lfacg(1,iface)
           if (ielem == jelem) jelem= lfacg(2,iface)
           if (jelem > 0) then
              if (izone .ne. lzone_read(jelem)) then   
                 ! contacto identificado, porque jelem y ielem no estan en la misma zona
                 ncont= ncont+1        ! numero de contactos
                 if (ncont .gt. 400000) then
                    write(6,*) 'aumentar el coso ese 1!'
                    stop
                 else
                    nnodc = nnodf(ielty) % l(iface_local)
                    do inodc =1, nnodc
                       inode = lface(ielty)%l(inodc,iface_local) 
                       ! los nodos viejos del contacto se quedan del lado del solido, o sea jelem
                       ! y los nuevos iran en poaux y son los del lado del fluido, o sea ielem
                       ipoin = lnods(inode,ielem)
                       jpoin = ipoin
                       if (poaux(ipoin) == 0) then
                          ipnew= ipnew+1
                          poaux(ipoin) = ipnew + npoin  ! la numeracion nueva para el nodo nuevo
                       end if
!ojj!!!                       lnoco(inodc      ,ncont)= poaux(ipoin) ! lado del fluido
!ojj!!!                       lnoco(inodc+nnodc,ncont)= jpoin        ! lado del solido
                       lnoco(inodc      ,ncont)= jpoin        ! lado del solido
                       lnoco(inodc+nnodc,ncont)= poaux(ipoin) ! lado del fluido
                    end do
                    ! el elemento que esta del lado fluido y que usara nodo nuevo
!ojj!!                    lnoco(nnodc*2+1,ncont) = ielem   
                    ! el elemento que esta del lado solido 
!oj!!                    lnoco(nnodc*2+2,ncont) = jelem   

                    ! el elemento que esta del lado solido 
                    lnoco(nnodc*2+1,ncont) = jelem   
                    ! el elemento que esta del lado fluido y que usara nodo nuevo
                    lnoco(nnodc*2+2,ncont) = ielem   
                    
!                    if (ielem == 131078) then
!                       write (6,*) 'totooooo',ielem,jelem,ncont
!                       write (6,*) lnoco(      1:  nnodc,ncont)
!                       write (6,*) lnoco(nnodc+1:2*nnodc,ncont)                       
!                       write (6,*) 'totooooo'
!                    end if

                    ! en 2d, la parte del contacto con nodos nuevos (fluido) va en orden inverso
                    if (ndime==2) then
                       do inodc= 1,nnodc
                          lauxi(nnodc - inodc + 1)= lnoco(inodc,ncont)  
                       end do
                       lnoco(1:nnodc,ncont)  = lauxi(1:nnodc) 
                    end if
                 end if
              end if
           end if
        end do
     end if
  end do

  npnew= ipnew
  write(6,*) 'Contactos:     ', ncont
  write(6,*) 'Puntos nuevos: ', npnew

  ! poner las coordenadas de los nuevos puntos igual a las del correspondiente en el contacto
  do ipoin= 1,npoin
     if (poaux(ipoin) > 0) then
        ipnew= poaux(ipoin) - npoin
        conew(1:ndime,ipnew) = coord(1:ndime,ipoin)
     end if
  end do
  
  ! corregir lnods 
  do ielem= 1,nelem
     izone= lzone_read(ielem)
     if (izone == 1) then                 ! solo mira los de la zona 1, que son los que va a corregir
        ielty= abs(ltype(ielem))
        pnode= lnnod(ielem)
        do inode=1,pnode
           ipoin= lnods(inode,ielem)
           if (ipoin .le. npoin) then     ! cuando el ipoin es mayor que npoin es que ya fue corregido 
              if (poaux(ipoin) > 0) then
                 lnods(inode,ielem)= poaux(ipoin)
              end if
           end if
        end do
     end if
  end do
  
  ! corregir lnodb
  do iboun= 1,nboun
     ielem= lnodb_element(iboun)
     izone= lzone_read(ielem)
     if (izone == 1) then                 ! solo mira los de la zona 1, que son los que va a corregir
        do inodb= 1,mnodb
           ipoin= lnodb(inodb,iboun)
           if ((ipoin .le. npoin).and.(ipoin.gt.0)) then     ! misma cosa que arriba
              if (poaux(ipoin) > 0) then
                 lnodb(inodb,iboun)= poaux(ipoin)
              end if
           end if
        end do
     end if
  end do

  ! escribir de nuevo el asunto

  write(606,'(a)') 'NODES_PER_ELEMENT'
  do ielem= 1,nelem
     write(606,100) ielem, lnnod(ielem)
  end do
  do icont= 1,ncont
     write(606,100) icont+nelem, nnodc*2
  end do
  write(606,'(a)') 'END_NODES_PER_ELEMENT'

  write(606,'(a)') 'ELEMENTS'
  do ielem= 1,nelem
     write(606,100) ielem,lnods(1:lnnod(ielem),ielem)     
  end do
  do icont= 1,ncont
     write(606,100) icont+nelem, lnoco(1:nnodc*2,icont)
  end do
  write(606,'(a)') 'END_ELEMENTS'

  write(6,*) '  Modify dom.dat:'
  write(6,*) '     old nelem=  ',nelem
  write(6,*) '     new nelem=  ',nelem+ncont

  write(610,*) '  Modify dom.dat:'
  write(610,*) '     old nelem=  ',nelem
  write(610,*) '     new nelem=  ',nelem+ncont


  write(606,'(a)') 'COORDINATES'
  do ipoin= 1,npoin
     write(606,200) ipoin,coord(1:ndime,ipoin)
!     write(636,*) ipoin,'  0.0  0.0  1.0'   ! this is to create a bogus fiber field
  end do
  do ipnew= 1,npnew
     write(606,200) ipnew+npoin,conew(1:ndime,ipnew)
!     write(636,*) ipnew+npoin,'  0.0  0.0  1.0'   ! this is to create a bogus fiber field
  end do
  write(606,'(a)') 'END_COORDINATES'

  write(6,*) '     old npoin=  ',npoin
  write(6,*) '     new npoin=  ',npoin+npnew
  write(610,*) '     old npoin=  ',npoin
  write(610,*) '     new npoin=  ',npoin+npnew

  write(606,'(a)') 'BOUNDARIES, ELEMENT'
  ! boundaries viejas en el dom.dat nuevo, con el lnodb corregido
  do iboun= 1,nboun
     write(606,100) iboun,lnodb(1:mnodb,iboun),lnodb_element(iboun)
  end do

  ! condiciones de contorno en las boundaries viejas en el fix.dat nuevo y con lnodb corregido
  do iboun=1,nboun
     ! hay que poner lnodb porque ya esta corregido para agregar los nodos nuevos
     write(608,100) iboun,mnodb,lnodb(1:mnodb+1,iboun)

  end do

  ! boundaries nuevas, que agregan las que salen de los contactos en 606 (dom.dat)
  ! tambien se agregan en 608 (fix.dat) las condiciones 
  ! sobre boundaries nuevas, que vienen de los contactos
  ! sobre las del fluido va un 30 y sobre las del solido un 40

  condition_fluid= 30
  condition_solid= 40

  jcont= 0
  do icont= 1,ncont
!     ielem_fluid = lnoco(nnodc*2+1,icont)
!     ielem_solid = lnoco(nnodc*2+2,icont)
!     jcont= jcont+1
!     write(606,100) jcont+nboun,        lnoco(nnodc+1:nnodc*2,icont), ielem_solid
!     write(608,100) jcont+nboun, mnodb, lnoco(nnodc+1:nnodc*2,icont), condition_solid
!     jcont= jcont+1
!     write(606,100) jcont+nboun,        lnoco(      1:nnodc  ,icont), ielem_fluid
!     write(608,100) jcont+nboun, mnodb, lnoco(      1:nnodc  ,icont), condition_fluid

     ielem_solid = lnoco(nnodc*2+1,icont)
     ielem_fluid = lnoco(nnodc*2+2,icont)
     jcont= jcont+1
     write(606,100) jcont+nboun,        lnoco(      1:nnodc  ,icont), ielem_solid
     write(608,100) jcont+nboun, mnodb, lnoco(      1:nnodc  ,icont), condition_solid
     jcont= jcont+1
     write(606,100) jcont+nboun,        lnoco(nnodc+1:nnodc*2,icont), ielem_fluid
     write(608,100) jcont+nboun, mnodb, lnoco(nnodc+1:nnodc*2,icont), condition_fluid

  end do

  write(606,'(a)') 'END_BOUNDARIES'
  write(608,'(a)') 'END_ON_BOUNDARIES'
!!!  write(608,'(a)') 'END_BOUNDARY_CONDITIONS'

  write(6,*) '     old nboun=  ',nboun
  write(6,*) '     new nboun=  ',nboun+jcont
  write(6,*) '     old nzone=              2 '
  write(6,*) '     new nzone=              3 '

  write(610,*) '     old nboun=  ',nboun
  write(610,*) '     new nboun=  ',nboun+jcont
  write(610,*) '     old nzone=              2 '
  write(610,*) '     new nzone=              3 '
    
  ! agrego los contactos en las zonas y le asigno la zona numero nzone+1
  write(607,*) '   ZONE, NUMBER=',nzone+1
  do icont= 1,ncont
     write(607,100) icont+nelem
  end do
  write(607,*) '   END_ZONE'
  write(607,'(a)') 'END_ZONES'

  itype = 3
  write(607,'(a)') 'CHARACTERISTICS'
  do icont= 1,ncont
     write(607,100) icont+nelem, itype
  end do
  write(607,'(a)') 'END_CHARACTERISTICS'
  write(607,*) 

  write(6,*) '  VERY IMPORTANT:'
  write(6,*) '     DO NOT FORGET TO CORRECT DOM.DAT, INCLUDING:  '
  write(6,*) 
  write(6,*) '     1. FOR 3D PROBLEMS, ADD ELEMENT 34 TO THE TYPES OF ELEMENTS    '
  write(6,*) '        TYPES_OF_ELEMENTS= 0, 0, 0, 0, 0, 30, 0, 0, 0, 34, 0 '
  write(6,*) '        FOR 2D PROBLEMS, ADD ELEMENT 12 (when needed) TO THE TYPES OF ELEMENTS    '
  write(6,*) '        TYPES_OF_ELEMENTS= 10, 0, 12, 0, 0, 0, 0, 0, 0, 0, 0 '
  write(6,*) 
  write(6,*) '     2. ADD EXTRAPOLATE TO THE BOUNDARY CONDITIONS'
  write(6,*) '        BOUNDARY_CONDITIONS, EXTRAPOLATE  '
  write(6,*) 
  write(6,*) '     3. ADD THE NEW BOUNDARY CONDITION TO THE FLUID AND TO ALEFOR:'
  write(6,*) '        *.NSI.DAT:'
  write(6,*) '          30 --> SOLID WALL OR SLIP CONDITION'
  write(6,*) '        *.ALE.DAT:'
  write(6,*) '          30 --> DIRICHLET CONDITION TO 0.0'
  write(6,*) 
  write(6,*) '  Done!  '

  write(610,*) '  VERY IMPORTANT:'
  write(610,*) '     DO NOT FORGET TO CORRECT DOM.DAT, INCLUDING:  '
  write(610,*) 
  write(610,*) '     1. FOR 3D PROBLEMS, ADD ELEMENT 34 TO THE TYPES OF ELEMENTS    '
  write(610,*) '        TYPES_OF_ELEMENTS= 0, 0, 0, 0, 0, 30, 0, 0, 0, 34, 0 '
  write(610,*) '        FOR 2D PROBLEMS, ADD ELEMENT 12 (when needed) TO THE TYPES OF ELEMENTS    '
  write(610,*) '        TYPES_OF_ELEMENTS= 10, 0, 12, 0, 0, 0, 0, 0, 0, 0, 0 '
  write(610,*) 
  write(610,*) '     2. ADD EXTRAPOLATE TO THE BOUNDARY CONDITIONS'
  write(610,*) '        BOUNDARY_CONDITIONS, EXTRAPOLATE  '
  write(610,*) 
  write(610,*) '     3. ADD THE NEW BOUNDARY CONDITION TO THE FLUID AND TO ALEFOR:'
  write(610,*) '        *.NSI.DAT:'
  write(610,*) '          30 --> SOLID WALL OR SLIP CONDITION'
  write(610,*) '        *.ALE.DAT:'
  write(610,*) '          30 --> DIRICHLET CONDITION TO 0.0'
  write(610,*) 

  100 format(10(2x,i6))
  200 format(i8,3(2x,e12.5))
  
end program conta

