program conta

  implicit none

  integer(4)         :: istat
  integer(8)               :: &
       npoin,nelem,nboun,ndime,mnode,inodb,mnodb,ipoin,iboun,ielem,jelem,inode,idime,&
       pnode,idumy,icucu,nzone,izone,jzone,ielem_fluid, ielem_solid,&
       ncnta,nnodc,inodc,ielct,ielto,ipnew,npnew,ncont,icont,nelct,itype,&
       condicion_solido,condicion_fluido,toca_nodo
  real(8)               :: digit

  integer(8)            :: lnoco(10,400000), lnozo(2,400000), lauxi(10), nelzo(2)
  real(8)                  :: conew(3,400000)
  integer(8), pointer      :: lnods(:,:),lnope(:),lzone(:),poaux(:),pzone(:,:),lnodb(:,:)
  real(8), pointer         :: coord(:,:)
  character(120)     :: input_file,output_file
  character(5)       :: wopos
  character(30)      :: tcard
  character(30)      :: ccard
  character(50)     :: wline

!--------------------
  nzone= 2
  npoin= -1
  nelem= 0
  nboun= 0
  ndime= 2
  mnode= 4
  mnodb= 2
  ncnta= 0
  nnodc= 2
!--------------------


  call getarg(1,input_file)  
  call getarg(2,tcard)
  read(tcard,'(i10)') npoin
  call getarg(3,tcard)
  read(tcard,'(i10)') nelem
  call getarg(4,tcard)
  read(tcard,'(i10)') nboun
  call getarg(5,tcard)
  read(tcard,'(i10)') ndime

  write(6,*) '   Usage: conta input-file  npoin  nelem  nboun ndime'
  
  if (npoin == -1) then
     stop
  end if

  if (ndime==3) then
     mnodb= 3
     nnodc= 3
  end if

  allocate(coord(ndime,npoin),stat=istat)
  allocate(lnods(mnode  ,nelem),stat=istat)
  allocate(lnodb(mnodb+1,nboun),stat=istat)
  allocate(lnope(nelem),stat=istat)
  allocate(lzone(nelem),stat=istat)
  allocate(poaux(npoin),stat=istat)
  allocate(pzone(2,npoin),stat=istat)
  pzone= 0
  
  output_file = trim(input_file) // '.geo.dat'
  write (6,*) '   Reading:  ',output_file 
  open (10 ,file=output_file)
  output_file = trim(input_file) // '-cnt.geo.dat'
  open (606,file=output_file)
  output_file = trim(input_file) // '-zones-cnt.geo.dat'
  open (607,file=output_file)
  output_file = trim(input_file) // '-bocon-cnt.geo.dat'
  open (608,file=output_file)
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
  write(608,*) wline

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  read(10,'(a)') 
  do ielem=1,nelem
     read(10,*) icucu,lnope(ielem)
  end do
  read(10,'(a)') 
  read(10,'(a)') 
  do ielem=1,nelem
     pnode= lnope(ielem)
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
     read(10,*) idumy,lnodb(1:mnodb+1,iboun)
  end do
  read(10,'(a)') 

  icucu= 0
  do while (icucu == 0)
     read(10,'(a)') wline
     write(607,*) wline
     wopos= trim(adjustl(wline))
     if (wopos == 'ZONES') then
        icucu= 1
     end if
  end do

  lzone= 0
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
           if (lzone(ielem) == 0) then
              lzone(ielem) = izone
              nelzo(izone) = nelzo(izone) + 1
              do inode= 1,lnope(ielem)
                 ipoin= lnods(inode,ielem)
                 pzone(izone,ipoin) = izone
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


  ielct= 0
  ncont= 0
  do ielem= 1,nelem
     izone= lzone(ielem)
     inodc= 0
     lauxi= 0
     toca_nodo= 0
     do inode= 1,lnope(ielem)
        ipoin= lnods(inode,ielem)
        jzone= pzone(1,ipoin)
        if (izone .ne. jzone) then
           inodc = inodc + 1
           lauxi(inodc)= ipoin
           toca_nodo = 1                  ! al menos un nodo toca otra zona
        end if
     end do

     if (ielem==22277) then
        write (6,*) 'cucuuuuu', ielem, lzone(ielem), inodc
     end if

     if (inodc == nnodc) then                               
        ! este elemento de volumen tiene nnodc nodos en una zona diferente de la suya
        jelem = 0
        if ((ielem==22277)) then           
           do while (jelem <= nelem)
              ! este bucle corrige casos patologicos de tipo "mariposa"
              jelem = jelem + 1
!              if (izone .ne. lzone(jelem)) then
                 ielto= 0              
                 pnode= lnope(jelem)
                 do inode= 1,pnode
                    if (lauxi(inode) == lnods(inode,jelem)) ielto = ielto + 1
                 end do
                 if (ielto>0) then
                    write (6,*) 'cucuuuuu', ielem, jelem, lzone(ielem), lzone(jelem)
                 end if
                 if ((ielto == nnodc))then                    
                    if (ielem /= jelem) then
                       ! los nodos de contacto pertenecen a un mismo elemento de la otra zona: es bien.
                       toca_nodo = -1
                       ! salimos del bucle de nelems
                       jelem = nelem+10
                    end if
                 end if
!              end if
           end do
              
        else
           toca_nodo = -1
        end if
        if (toca_nodo == -1) then 
           ! estos son elementos que tienen un contacto: tienen nnodc nodos entre dos zonas
           ! van a ser los que reemplacen esos nodos por nuevos 
           !        
           ielct= ielct+1                                      
           ncont= ncont+1                                      
           if (ielct .gt. 400000) then
              write(6,*) 'aumentar el coso ese 1!'
              stop
           else
              do inodc= 1,nnodc
                 ! los primeros nodos en el contacto son los que no se cambiaran
                 lnoco(inodc,ielct)= lauxi(inodc)              
              end do
              lnoco(nnodc*2+2,ielct) = ielem                   ! el elemento que esta del otro lado y que usara nodo nuevo
              lnozo(2,ielct) = izone                           ! y esta es la zona de ese elemento
              !           write(6,100) lnoco(1:6,ielct)
           end if
        end if
     else if(inodc > 0) then                               
        ! en cambio, este elemento de volumen tiene algunos (aunque menos que nnodc) nodos en una zona diferente de la suya

        toca_nodo = 1

     end if
          
     if (toca_nodo > 0) then                               
        ! estos son elementos que tienen algun nodo de contacto: tienen menos de nndoc o estan "mariposeados"
        ! hay que marcarlos para reemplazar esos nodos por los nuevos. suman a ielct pero no a ncont
        ielct= ielct+1
        if (ielct .gt. 400000) then
           write(6,*) 'aumentar el coso ese 2!'
           stop
        else
           do inodc= 1,nnodc
              lnoco(inodc,ielct)= lauxi(inodc)              ! el nodo que hay que cambiar
           end do
           lnoco(nnodc*2+2,ielct) = - ielem                 ! el elemento que esta de un lado y que usara nodo nuevo
           !           write(6,100) lnoco(1:6,ielct)
        end if
     end if
  end do

  stop

  write(6,*) 'Contactos:     ', ncont
  nelct= ielct
  
  poaux= 0
  ipnew= 0
  do ielct= 1,nelct
     ielem = lnoco(nnodc*2+2,ielct)
     if (ielem > 0) then
        ielto= 0
        do inodc= 1,nnodc
           ipoin= lnoco(inodc,ielct)
           if (poaux(ipoin) == 0) then
              ielto= ielto + 1
              ipnew= ipnew+1
              poaux(ipoin)= ipnew + npoin
              do idime= 1,ndime
                 conew(idime,ipnew)= coord(idime,ipoin)
              end do
           end if
           ! agrega el nodo nuevo al elemento de contacto y así lo completa en el orden correcto (inverso en 2d)
           if (ndime == 2) then
              lnoco(2*nnodc - inodc + 1,ielct)= poaux(ipoin)  
           else if (ndime == 3) then
              lnoco(nnodc + inodc,ielct)= poaux(ipoin)                
           end if
           if (ielto == nnodc) then
              ! esto es para identificar el elemento que esta del lado del solido en el contacto
              

              
           end if
        end do
!        write(6,100) lnoco(1:6,ielct)
     end if
  end do

  npnew= ipnew
  write(6,*) 'Puntos nuevos: ', npnew


!  write(6,*) 'ooo',poaux(2457)
!  stop

  ! corregir lnods y lnodb
  do ielct= 1,nelct
     ielem = lnoco(nnodc*2+2,ielct)     
     if (ielem == 0) then
        ! no modificar elemento, pero detectar si 

     else
        if (ielem < 0) then 
           ielem = - ielem
        end if
        ! agrega el nodo nuevo al elemento de volumen
        pnode= lnope(ielem)
!        write(6,*) 'viejo:  ',lnods(1:pnode,ielem), lnoco(nnodc*2+1,ielct)
        do inode=1,pnode
           ipoin= lnods(inode,ielem)
           if (poaux(ipoin) > 0) then
              lnods(inode,ielem) = poaux(ipoin)
           end if
        end do
!        write(6,*) 'nuevo:  ',lnods(1:pnode,ielem)
     end if

     do iboun= 1,nboun
        jelem= lnodb(mnodb+1,iboun)
        if (ielem == 0) then
           ! no modificar el contorno
        else
           if (ielem == jelem) then
              do inodb=1,mnodb
                 ipoin= lnodb(inodb,iboun)
                 if (ipoin <= npoin ) then
                    if (poaux(ipoin) > 0) then
                       lnodb(inodb,iboun) = poaux(ipoin)
                    end if
                 end if
              end do
           end if
        end if
     end do
  end do

  ! check the element number on the other side of the contact
  do ielct= 1,nelct
     ielem = lnoco(nnodc*2+2,ielct)
     if (ielem > 0) then
        jelem = 0
        do while (jelem <= nelem)
           jelem = jelem + 1
           ielto= 0
           pnode= lnope(jelem)
           do inode= 1,pnode
              do inodc= 1,nnodc
                 if (lnoco(inodc,ielct) == lnods(inode,jelem)) ielto = ielto + 1
              end do
           end do
           if (ielto == nnodc) then
              if (ielem /= jelem) then
                 ! este es el elemento que esta del otro lado
                 lnoco(nnodc*2+1,ielct) = jelem
                 lnozo(1,ielct) = lzone(jelem)
                 ! salimos del bucle de nelems
                 jelem = nelem+10
              end if
           end if
        end do
     end if
  end do

  ! escribir de nuevo el asunto

  write(606,'(a)') 'NODES_PER_ELEMENT'
  do ielem= 1,nelem
     write(606,100) ielem, lnope(ielem)
  end do
  icont= 0
  do ielct= 1,nelct
     ielem = lnoco(nnodc*2+2,ielct)
     if (ielem > 0) then
        icont= icont+1
        write(606,100) icont+nelem, nnodc*2
     end if
  end do
  write(606,'(a)') 'END_NODES_PER_ELEMENT'
  write(606,'(a)') 'ELEMENTS'
  do ielem= 1,nelem
     pnode= lnope(ielem)
     write(606,100) ielem,lnods(1:pnode,ielem)     
  end do
  icont= 0
  do ielct= 1,nelct
     ielem = lnoco(nnodc*2+2,ielct)
     if (ielem > 0) then
        icont= icont+1
        write(606,100) icont+nelem, lnoco(1:nnodc*2,ielct)
        write(306,100) lnoco(1:nnodc*2+2,ielct), lnozo(1:2,ielct)
     end if
  end do
  write(606,'(a)') 'END_ELEMENTS'

  write(6,*) '  Modify dom.dat:'
  write(6,*) '     old nelem=  ',nelem
  write(6,*) '     new nelem=  ',nelem+icont

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

  write(606,'(a)') 'BOUNDARIES, ELEMENT'
  ! boundaries viejas  
  do iboun= 1,nboun
     write(606,100) iboun,lnodb(1:mnodb+1,iboun)
  end do
  ! boundaries nuevas, que agregan las que salen de los contactos
  ! tambien se agregan en 608 las condiciones sobre boundaries nuevas, que vienen de los contactos
  ! sobre las del fluido va un 30 y sobre las del solido un 40

  condicion_fluido= 30
  condicion_solido= 40

  do iboun=1,nboun
     read(20,*) icucu,idumy,lauxi(1:mnodb+1)     
! usar esto por si hay que corregir alguna condicion de contorno
!!     if (pzone(1,lnodb(1,iboun))==2) then
!!        if (lauxi(mnodb+1)==1) then
!!           lauxi(mnodb+1)  = 5  
!!           write(6,100) iboun,idumy,lnodb(1:mnodb,iboun),lauxi(mnodb+1)  
!!        end if
!!     end if

     ! hay que poner lnodb porque ya esta corregido
     write(608,100) iboun,idumy,lnodb(1:mnodb,iboun),lauxi(mnodb+1)  
  end do

  icont= 0
  do ielct= 1,nelct
     ielem_solid = lnoco(nnodc*2+1,ielct)
     ielem_fluid = lnoco(nnodc*2+2,ielct)
     if (ielem_fluid > 0) then
        !solido, zona 2
        icont= icont+1
        write(606,100) icont+nboun, lnoco(      1:nnodc  ,ielct),ielem_solid
!! la condicion del solido no hace falta agregarla en el bocon
        write(608,100) icont+nboun, mnodb, lnoco(      1:nnodc  ,ielct), condicion_solido
        !fluido, zona 1
        icont= icont+1
        write(606,100) icont+nboun, lnoco(nnodc+1:nnodc*2,ielct),ielem_fluid
        write(608,100) icont+nboun, mnodb, lnoco(nnodc+1:nnodc*2,ielct), condicion_fluido
     end if
  end do
  write(606,'(a)') 'END_BOUNDARIES'
  write(608,'(a)') 'END_ON_BOUNDARIES'
!!!  write(608,'(a)') 'END_BOUNDARY_CONDITIONS'

  write(6,*) '     old nboun=  ',nboun
  write(6,*) '     new nboun=  ',nboun+icont
  write(6,*) '     old nzone=              2 '
  write(6,*) '     new nzone=              3 '
    

  ! agrego los contactos en las zonas y le asigno la zona numero nzone+1
  write(607,*) '   ZONE, NUMBER=',nzone+1
  do icont= 1,ncont
     write(607,100) icont+nelem
  end do
  write(607,*) '   END_ZONE'
  write(607,'(a)') 'END_ZONES'

  itype = 3
  write(607,'(a)') 'CHARACTERISTICS'
  icont= 0
  do ielct= 1,nelct
     ielem = lnoco(nnodc*2+2,ielct)
     if (ielem > 0) then
        icont= icont+1
        write(607,100) icont+nelem, itype
     end if
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

  100 format(10(2x,i6))
  200 format(i8,3(2x,e12.5))
  
end program conta

