program ensight2alya 
  ! compile with: ifort -c -traceback *.f90  && ifort -traceback *.o

  use,intrinsic :: ISO_C_BINDING
  use kdtree2_module

  implicit none


  !------------------------------------------------------------------------
  !
  ! ensipart
  !
  !------------------------------------------------------------------------
  integer, parameter  :: ip = 4             ! 4-byte integer
  integer, parameter  :: rp = 8             ! Double precision 

  integer, parameter  :: maxpart = 10000      ! Maximum number of parts allowed in ensi geo

  type lntyp
     integer(ip), pointer    :: lnods(:,:)  
  end type lntyp

  type ensipart
     type(lntyp), pointer    :: ln(:)
     real(4),     pointer    :: coor(:,:)    ! Coordinates   ! ojo de tipo 4
     real(rp)                :: exten(6)
     integer(ip)             :: ltyel(10)     ! using Alya's convention  
     integer(ip)             :: np  
     integer(ip)             :: ne(10)
     integer(ip), pointer    :: lnvol(:)  
  end type ensipart

  type(kdtree2), pointer     :: tree, tree2, tree3
  type(kdtree2_result)       :: results(1)
  real(rp)                   :: t0,t1
  real(4)                    :: dists(1)
  integer(ip)                :: inds(1)


  type(ensipart), pointer    :: enspa(:)
  character(80)              :: charr
  integer(ip)                :: ii,jj,ipoin,nnode,in,kgrel,ie,iost

  integer(ip)                :: kount,iexel(40),nelem,nboun,jp,id
  real(rp)                   :: dista,auxil

  character(150)             :: namda
  integer(ip)                :: lnump(maxpart),kk
  integer(ip)                :: kfl_nidgi,kfl_eidgi,iiaux,ielem

  kfl_nidgi = 0_ip
  kfl_eidgi = 0_ip

  nullify(tree)
  nullify(tree2)
  nullify(tree3)
  nullify(enspa)

  call GETARG(1_4,namda)
  if( trim(namda) == '' ) print*,'Wrong problem name'

  open (unit=10,file=trim(namda)//'00000.geo',&
       &     form='unformatted'                   ,&
       &     recordtype='stream'                  ,&
       &     action='read'                       ,&
       &     convert='LITTLE_ENDIAN'                 ,&
       !       &     convert='BIG_ENDIAN'                 ,&
       &     access='sequential'                   )

  ! para saber que es little endian http://en.wikipedia.org/wiki/Endianness
  ! y despues jugando con los valores que salian usando la calculator del mac en programmer mode 
  ! escribimos 00 01 6b 9b  y sale 93083 que yo sabia era el primer avlor que tenia que salir
  ! Litter endian es lo que usan las maquinas intel  -- big endia marenostrum
  ! para entender un poc todo esto use un hex editor
  ! lo que es caracteres lo pone casa uno con un numero  ej 20 corresponde a space (http://www.ascii.cl/)
  ! 0A es el line finish (LF)


  allocate(enspa(maxpart))

  do ii=1,3_ip
     read(10) charr
     print*,'charr',charr
  end do

  read(10) charr
  if (index(charr,'given') /= 0_ip )   kfl_nidgi = 1_ip
  print*,'charr',charr

  read(10) charr
  if (index(charr,'given') /= 0_ip )   kfl_eidgi = 1_ip
  print*,'charr',charr

  do ii=7_ip,7_ip
     read(10) charr
     print*,'charr',charr   ! part
  end do

  jj = 0_ip
  parts: do 

     jj = jj + 1_ip

     nullify(enspa(jj) % ln)
     nullify(enspa(jj) % coor)
     nullify(enspa(jj) % lnvol)

     read(10) kk
     print*,'number of part:jj',jj,'kk',kk

     if( jj > maxpart ) then
        write(*,*) 'Stopped because the mximum number of parts allowed in ensi geo is',maxpart 
        stop
     end if
     lnump(jj) = kk

     do ii=8_ip,9_ip
        !     read(10) chara(ii)
        read(10) charr
        print*,'charr',charr
     end do

     kgrel = 0_ip
     read(10) enspa(jj) % np
     print*,'np',enspa(jj) % np

     allocate (enspa(jj) % coor (3,enspa(jj) % np))

     if ( kfl_nidgi == 1_ip ) then
        do ipoin=1,enspa(jj) % np  ! por ahora no lo quiero para nada
           read(10) iiaux
        end do
     end if


     read(10) (enspa(jj) % coor(1,ipoin),ipoin=1,enspa(jj) % np)
     read(10) (enspa(jj) % coor(2,ipoin),ipoin=1,enspa(jj) % np)
     read(10) (enspa(jj) % coor(3,ipoin),ipoin=1,enspa(jj) % np)

     print*,'first coor',enspa(jj) % coor(1:3,1_ip)
     print*,'second coor',enspa(jj) % coor(1:3,2_ip)
     print*,'last coor',enspa(jj) % coor(1:3,enspa(jj) % np)

     allocate (enspa(jj) % ln(10)) 

     do iiaux=1,10
        nullify(enspa(jj) % ln(iiaux) % lnods)
     end do

     group_elements: do 

        do ii=10_ip,10_ip
           !     read(10) chara(ii)
           read(10,iostat=iost) charr
           if (iost/=0) exit parts
           print*,'charr',charr
           print*,'index(tetra4,charr)',index(charr,'tetra4')
        end do


        if (index(charr,'tria3') /= 0) then
           nnode=3
           enspa(jj) % ltyel(kgrel + 1_ip) = 10
        else if (index(charr,'quad4') /= 0) then
           nnode=4
           enspa(jj) % ltyel(kgrel + 1_ip) = 12
        else if (index(charr,'tetra4') /= 0) then
           nnode=4
           enspa(jj) % ltyel(kgrel + 1_ip) = 30
        else if (index(charr,'pyramid5') /= 0) then
           nnode=5
           enspa(jj) % ltyel(kgrel + 1_ip) = 32
        else if (index(charr,'penta6') /= 0) then
           nnode=6
           enspa(jj) % ltyel(kgrel + 1_ip) = 34
        else if (index(charr,'hexa8') /= 0) then
           nnode=8
           enspa(jj) % ltyel(kgrel + 1_ip) = 37
        else
           print*,'exit_group_elements,charr',charr
           exit group_elements
        end if

        print*,'enspa(jj) % ltyel(kgrel + 1_ip)',enspa(jj) % ltyel(kgrel + 1_ip)


        kgrel = kgrel + 1_ip
        read(10) enspa(jj) % ne(kgrel)

        print*,'DDDDDDDDDDDDDDDDD: jj,kgrel,enspa(jj) % ne(kgrel)',jj,kgrel,enspa(jj) % ne(kgrel)


        if (kgrel>10) stop 'max 10 groups of elem per region'
        !        allocate (enspa(jj) % ln(10)) 
        allocate (enspa(jj) % ln(kgrel) % lnods (nnode,enspa(jj) % ne(kgrel)))
        if ( associated (enspa(jj) % ln(kgrel) % lnods) ) print*,'BBBBBBBBBB: asosciated kgrel',kgrel
        if ( kfl_eidgi == 1_ip ) then
           do ielem=1,enspa(jj) % ne(kgrel)  ! por ahora no lo quiero para nada
              read(10) iiaux
           end do
        end if
        read(10) ((enspa(jj) % ln(kgrel) % lnods(in,ie),in=1,nnode),ie=1,enspa(jj) % ne(kgrel))

        print*, (enspa(jj) % ln(kgrel) % lnods(in,enspa(jj) % ne(kgrel)),in=1,nnode)

     end do group_elements

  end do parts

  !
  ! process
  ! aca para las coord de las boundaries hay que encontrar sus respectivos nodos en el volumen
  ! supongo qe el volumen es la primer parte sino habrái que complicarlo un poco
  !
  jj = 1_ip   ! el volumen
  kount = 0_ip
  iexel = 0_ip
  nelem = 0_ip
  do kgrel =1,10 
     if (associated(enspa(jj) % ln(kgrel) % lnods) ) then
        kount = kount + 1
        if ( enspa(jj) % ltyel(kgrel) == 30 ) iexel(30_ip) = 1_ip
        if ( enspa(jj) % ltyel(kgrel) == 34 ) iexel(34_ip) = 1_ip
        if ( enspa(jj) % ltyel(kgrel) == 37 ) iexel(37_ip) = 1_ip
        nelem = nelem + enspa(jj) % ne(kgrel)
     end if
  end do
  !
  !   Boun
  !

  ! probar un poco y luego elegir el que mejor vaya
  call cpu_time(t0)
  tree => kdtree2_create( enspa(1) % coor ,sort=.false.,rearrange=.false.)  ! this is how you create a tree. 
  call cpu_time(t1)
  write (*,*) real( size(enspa(jj) % coor(3,:)) )/real(t1-t0), ' points per second built for non-rearranged tree.'

  call cpu_time(t0)
  tree2 => kdtree2_create( enspa(1) % coor,sort=.false.,rearrange=.true.)  ! this is how you create a tree. 
  call cpu_time(t1)
  write (*,*) real( size(enspa(jj) % coor(3,:)) )/real(t1-t0), ' points per second built for rearranged tree.'

  call cpu_time(t0)
  !  tree3 => kdtree2_create(cooau2,sort=.true.,rearrange=.true.)  ! this is how you create a tree. 
  tree3 => kdtree2_create(enspa(1) % coor ,sort=.true.,rearrange=.true.)  ! this is how you create a tree. 
  call cpu_time(t1)
  write (*,*) real( size(enspa(jj) % coor(3,:)) )/real(t1-t0), ' points per second built for 2nd rearranged tree.'


  nboun = 0_ip
  do jj = 2_ip,20_ip
     do kgrel =1,10 

        if ( associated(enspa(jj) % ln) ) then
           if ( enspa(jj) % ne(kgrel) > 0_ip) then
              nboun = nboun + enspa(jj) % ne(kgrel)
              print*,'CCCCCCCCCCCC: nboun',nboun
              print*,'EEEEEEEEEEEEEEEE: jj,kgrel,enspa(jj) % ne(kgrel)',jj,kgrel,enspa(jj) % ne(kgrel)

              allocate(enspa(jj) % lnvol(enspa(jj) % np))
              do ipoin = 1,enspa(jj) % np
                 ! find one nearest neighbors to.
                 inds = -666 ! just inicialize to a negative number - at least this is whay i guessed
                 call kdtree2_n_nearest(tp=tree3,qv=enspa(jj) % coor(:,ipoin) , nn=1_ip, results=results)
                 inds(1) = results(1)%idx
                 enspa(jj) % lnvol(ipoin) = results(1)%idx
                 !                  points_volume :do jp = 1,enspa(1_ip) % np  ! very unefficient implem   !replaced with kdtree
                 !                     dista = 0.0_rp
                 !                     do ii =1,3
                 !                        auxil = enspa(1_ip) % coor(ii,jp) - enspa(jj) % coor(ii,ipoin) 
                 !                        dista = dista + auxil * auxil 
                 !                     end do
                 !                     dista = sqrt (dista)
                 !                     if (dista < 1.0d-15) then !  usar mejor mult por una dimension del dominio
                 !                        enspa(jj) % lnvol(ipoin) = jp
                 !                        print*,'jp,ipoin,enspa(1_ip) % coor(:,jp),enspa(jj) % coor(:,ipoin)',jp,ipoin,enspa(1_ip) % coor(:,jp),enspa(jj) % coor(:,ipoin)
                 !                        exit points_volume
                 !                     end if
                 !                  end do points_volume
                 !                  print*,'jj,ipoin,enspa(jj) % lnvol(ipoin),inds,dists',jj,ipoin,enspa(jj) % lnvol(ipoin),inds,dists
              end do
           end if
        end if
     end do
  end do

  open (unit=11,file='1.dom.dat')
  open (unit=12,file='1.geo.dat')
  open (unit=13,file='1.fix')

  charr=' '
  if (iexel(30_ip) == 1_ip) charr = trim(charr) // '  TET04'
  if (iexel(32_ip) == 1_ip) charr = trim(charr) // '  PYR05'
  if (iexel(34_ip) == 1_ip) charr = trim(charr) // '  PEN06'
  if (iexel(37_ip) == 1_ip) charr = trim(charr) // '  HEX08'

  write(11,'(a)')   '$--------------------------------------------------'
  write(11,'(a)')   'DIMENSIONS'

  write(11,'(a,1x,i10)')   'NODAL_POINTS= ',enspa(1_ip) % np
  write(11,'(a,1x,i10)')   'ELEMENTS== ',nelem
  write(11,'(a)')   'SPACE_DIMENSIONS=         3'
  write(11,'(a)')   'TYPES_OF_ELEMS= '//charr
  write(11,'(a,1x,i10)')   'BOUNDARIES= ',nboun


  write(11,'(a)')   'END_DIMENSIONS'
  write(11,'(a)')   '$--------------------------------------------------'
  write(11,'(a)')   'STRATEGY'

  charr=' '
  if (iexel(30_ip) == 1_ip) charr = trim(charr) // '  1'
  if (iexel(32_ip) == 1_ip) charr = trim(charr) // '  5'
  if (iexel(34_ip) == 1_ip) charr = trim(charr) // '  6'
  if (iexel(37_ip) == 1_ip) charr = trim(charr) // '  8'


  write(11,'(a)')   '  DOMAIN_INTEGRATION_POINTS: '//charr
  write(11,'(a)')   '$  SCALE: XSCAL=0.001,YSCAL=0.001,ZSCAL=0.001'
  write(11,'(a)')   'END_STRATEGY'
  write(11,'(a)')   '$--------------------------------------------------'
  write(11,'(a)')   'GEOMETRY'
  write(11,'(a)')   '  GROUPS=650'
  write(11,'(a)')   '  INCLUDE ../mesh/1.geo.dat'
  write(11,'(a)')   'END_GEOMETRY'
  write(11,'(a)')   '$--------------------------------------------------'
  write(11,'(a)')   'SETS'
  write(11,'(a)')   '  BOUNDARIES'
  write(11,'(a)')   '    INCLUDE ../mesh/1.fix'
  write(11,'(a)')   '  END_BOUNDARIES'
  write(11,'(a)')   'END_SETS'
  write(11,'(a)')   '$--------------------------------------------------'
  write(11,'(a)')   'BOUNDARY_CONDITIONS, EXTRAPOLATE'
  write(11,'(a)')   '  ON_BOUNDARIES'
  write(11,'(a)')   '    INCLUDE ../mesh/1.fix'
  write(11,'(a)')   '  END_ON_BOUNDARIES'
  write(11,'(a)')   '   GEOMETRICAL_CONDITIONS'
  write(11,'(a)')   '     WALL_LAW:     3,6,7'
  write(11,'(a)')   '     INFLOW:       1,2'
  write(11,'(a)')   '     OUTFLOW:      5'
  write(11,'(a)')   '     ANGLE,   GEOAN=170.0'
  write(11,'(a)')   '   END_GEOMETRICAL_CONDITIONS'
  write(11,'(a)')   'END_BOUNDARY_CONDITIONS'
  write(11,'(a)')   '$--------------------------------------------------'
  !
  ! 1.dom.geo & 1.fix
  !
  write(12,'(a)')   'TYPES'
  kount = 0_ip
  jj = 1_ip
  do kgrel =1,10 
     if (associated(enspa(jj) % ln(kgrel) % lnods) ) then
        do ie = 1,enspa(jj) % ne(kgrel)
           kount = kount + 1
           write(12,*) kount,(enspa(jj) % ltyel(kgrel) )
        end do
     end if
  end do
  write(12,'(a)')   'END_TYPES'


  write(12,'(a)')   'ELEMENTS'
  kount = 0_ip
  jj = 1_ip
  do kgrel =1,10 
     if ( associated (enspa(jj) % ln(kgrel) % lnods) ) then
        print*,'asosciated kgrel',kgrel
        do ie = 1,enspa(jj) % ne(kgrel)
           kount = kount + 1
           if (enspa(jj) % ltyel(kgrel) == 37) nnode = 8
           if (enspa(jj) % ltyel(kgrel) == 34) nnode = 6
           if (enspa(jj) % ltyel(kgrel) == 30) nnode = 4
           if (enspa(jj) % ltyel(kgrel) == 32) nnode = 5
           if (enspa(jj) % ltyel(kgrel) == 34) then   ! trucheda porque sino alya da negative volume
              write(12,'(10(i10,1x))') kount,(enspa(jj) % ln(kgrel) % lnods(in,ie),in=4,6), &
                   &           (enspa(jj) % ln(kgrel) % lnods(in,ie),in=1,3)
           else
              write(12,'(10(i10,1x))') kount,(enspa(jj) % ln(kgrel) % lnods(in,ie),in=1,nnode )
           end if
        end do
     end if
  end do
  write(12,'(a)')   'END_ELEMENTS'

  write(12,'(a)')   'COORDINATES'
  kount = 0_ip
  jj = 1_ip
  do ipoin = 1,enspa(jj) % np
     kount = kount + 1
     write(12,'(i10,3(1x,e14.7))') kount,(enspa(jj) % coor(id,ipoin),id=1,3_ip )
  end do
  write(12,'(a)')   'END_COORDINATES'

  write(12,'(a)')   'BOUNDARIES'
  kount = 0_ip
  do jj = 2_ip,20_ip
     do kgrel =1,10 
        if ( associated(enspa(jj) % coor ) ) then
           if (associated(enspa(jj) % ln(kgrel) % lnods) ) then
              do ie = 1,enspa(jj) % ne(kgrel)
                 kount = kount + 1
                 if (enspa(jj) % ltyel(kgrel) == 10) nnode = 3
                 if (enspa(jj) % ltyel(kgrel) == 12) nnode = 4
                 write(12,'(10(i10,1x))') kount,( enspa(jj) % lnvol ( enspa(jj) % ln(kgrel) % lnods(in,ie)),in=1,nnode)
                 !                 ! invertí para evitar ALE_INIRBO: NEGATIVE VOLUME FOR PARTICLE 1. CHECK BOUNDARY ORIENTATION   
                 !                 ! Luego alya volvio a preferir lo orig ??
                 !                 write(12,'(10(i10,1x))') kount,( enspa(jj) % lnvol ( enspa(jj) % ln(kgrel) % lnods(in,ie)),in = nnode,1,-1)
                 write(13,'(2(i10,1x))') kount,jj
              end do
           end if
        end if
     end do
  end do
  write(12,'(a)')   'END_BOUNDARIES'

end program ensight2alya
