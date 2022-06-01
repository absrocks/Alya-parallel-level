program iris2alya


  use def_parame
  use def_elmtyp
  use def_master
  use def_domain , only : lnods,lnodb,coord,ltype,lexis,lnnod,nface,nnodf,lface
  use def_domain , only : npoin,ndime,nelem,nelem_2,nboun,mnodb,mnode
  use mod_memchk
  use mod_htable
  implicit none

  integer(ip)               :: &
       nelem_exterior, nelem_interior, npoin_total, ipoin_total,jelem,ielem_interior,&
       idumy, npoin_surfa, istack, nstack,&
       tetra_nodes(4), tetra_nodes_origi(4), tria_nodes(3), inode, ipoin, icucu, kboun, nboun_all, nboun_final, ielem,&
       ielty, iface_local, iface, nnodc,inodc, ielem_seed(3),iwrite,kbocondi(3)
  integer(ip)               :: &
       kfile,iargus,argus,ifmaterzones,iffibers,ifbocon,nelem_read,nmoh,imoh
  real(rp)               :: coord_nodes(3),fiber_nodes(3)

  real(rp)          ::  &
       angle, new_z(3), old_z(3), rot_ax(3),vecaux(3),rotvec(3), &
       dotpro,&
       modve, sinang, cosang,v(3)

  real(rp)    :: points_irregu(3,5316), points_regu(3,13442)

  integer(ip)               :: &
       nelem_material(50)

  integer(ip) , pointer  :: kpermu_nodes(:),lboun_global(:,:), lstack(:), markel(:),lnods_material(:)

!!  character(120)     :: input_file
  character(5)       :: wopos
  character(30)      :: tcard
  character(30)      :: ccard
  character(50)     :: wline
  character(1000)   :: buffi

!!  call getarg(2,tcard)
!!  read(tcard,'(i10)') npoin


!!!!! ESTAS LINEAS COMENTARLAS CUANDO NO HAY CSV CON CONDICIONES QUE SALEN DE PARAVIEW

!!  open (10 ,file='irregu0.csv')
!!  read(10,'(a)')
!!  do ipoin= 1,5316
!!     read(10,'(a)') buffi
!!     call parsecsv(buffi,v)
!!     points_irregu(1:3,ipoin) = v(1:3)
!!  end do
!!  close(10)
!!
!!  open (10 ,file='regu0.csv')
!!  read(10,'(a)')
!!  do ipoin= 1,13442
!!     read(10,'(a)') buffi
!!     call parsecsv(buffi,v)
!!     points_regu(1:3,ipoin) = v(1:3)
!!  end do
!!
!!  close(10)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !   iris2alya.input se crea haciendo 'wc -l tet* > iris2alya.input'


  ! defaults
  ifmaterzones= 0  ! materials: yes
  iffibers= 1      ! write fibers: yes
  ifbocon = 1      ! write bocon: yes

  argus= command_argument_count()
  do iargus= 1,argus
     call getarg(iargus,wline)  
     if (trim(adjustl(wline))=='materials') then
        write (6,*) 'Regions are materials...' 
        ifmaterzones=1
     end if
     if (trim(adjustl(wline))=='nofibers' ) then
        write (6,*) 'No fiber information required...' 
        iffibers=0
     end if
     if (trim(adjustl(wline))=='nobocon' )  then
        write (6,*) 'No boundary conditions computed...' 
        ifbocon=0          
     end if
  end do

  open (10 ,file='iris2alya.input')
  read(10,*) npoin_total, wline
  open (22 ,file=trim(wline))
  read(10,*) nelem_exterior, wline      ! external
  open (20 ,file=trim(wline))           ! internal regions ...
  read(10,*) nelem_interior, wline
  open (21 ,file=trim(wline))
  write(6,*) 'Total nodes read: ', npoin_total, trim(adjustl(wline)) 

  nmoh=1

  if (ifmaterzones == 0) then
     ! eliminating regions and creating contacts or boundaries
     read(10,*) nelem_exterior, wline      ! external 
     open (20 ,file=trim(wline))
     read(10,*) nelem_interior, wline      ! internal
     open (21 ,file=trim(wline))

     nelem= nelem_interior

  else if (ifmaterzones == 1) then
     ! regions are materials, do not eliminate anything          
     wopos='total'
     do
        read(10,*,end=1001) nelem_material(nmoh), wline    
        if (trim(adjustl(wline)) .ne. trim(adjustl(wopos))) then
           open (50+nmoh,file=trim(wline))                    
           nelem= nelem + nelem_material(nmoh)
           write(6,*) 'Material: ', nmoh, nelem_material(nmoh), nelem, trim(adjustl(wline)) 
           nmoh= nmoh+1
        end if
        if (nmoh == 51) then
           write(6,*) 'ERROR: More than 50 materials... correct dim in source file. ' 
           stop
        end if
     end do

1001 nmoh= nmoh-1
     write(6,*) 'Processing materials=  ',nmoh 

     
  end if

!!  nelem= nelem_exterior
  write(6,*) 'Active elements: ',nelem 

  open (31 ,file='eleme.dom.geo')

  ! 
  ! uso estas coordenadas para definir el nuevo eje z NECESARIO PARA EL CASO DE ROTAR
  new_z(1)= -166.13_rp  + 127.554_rp
  new_z(2)= -68.9136_rp - 110.858_rp
  new_z(3)= -146.659_rp + 8.4514_rp

  modve = sqrt(new_z(1)*new_z(1) + new_z(2)*new_z(2) + new_z(3)*new_z(3))
  new_z = - new_z/modve
  
  ! este es el viejo eje z 
  old_z(1) = 0.0_rp
  old_z(2) = 0.0_rp
  old_z(3) = 1.0_rp

  call cropro(new_z,old_z,rot_ax)
  cosang = sqrt(dotpro(new_z,new_z)) * sqrt(dotpro(old_z,old_z))
  cosang = dotpro(new_z,old_z) / cosang
  angle  = acos(angle)
  sinang = sin(angle)

  allocate(kpermu_nodes(npoin_total))
  kpermu_nodes= 0

  mnodb= 3
  nnodc= 3
  ndime= 3
  mnode= 4 
  allocate(lnods(mnode  ,nelem))
  allocate(lnnod(nelem))
  allocate(lnods_material(nelem))

  lnnod = 4 ! all tetra
  npoin = 0
  write(31,*) 'ELEMENTS'
!!  write(23,*) 'ELEMENTS'

  jelem= 0
  do imoh= 1,nmoh
     if (ifmaterzones == 0) then
        kfile= 21
        nelem_read= nelem
     else if (ifmaterzones == 1) then
        kfile= 50+imoh 
        nelem_read= nelem_material(imoh)
     end if
     write(6,*) 'Reading material: ',imoh,nelem_read 
     do ielem= 1,nelem_read
        jelem= jelem+1
        !!  do ielem= 1,50
        read(kfile,*) idumy, tetra_nodes(1:4)
        do inode= 1,4
           tetra_nodes_origi(inode)= tetra_nodes(inode)
           if (kpermu_nodes(tetra_nodes(inode)) == 0) then
              npoin= npoin+1
              kpermu_nodes(tetra_nodes(inode)) = npoin
           end if
           tetra_nodes(inode) = kpermu_nodes(tetra_nodes(inode))
        end do
        ! esta permutacion es necesaria por como define las normales de los tetras moh
        icucu= tetra_nodes(3)
        tetra_nodes(3) = tetra_nodes(4)
        tetra_nodes(4) = icucu
        write(31,100) jelem, tetra_nodes(1:4)
        !!     write(23,150) ielem, tetra_nodes(1:4), idumy, tetra_nodes_origi(1:4)  !! descomentar esto para chequear con moh's
        lnods(1:4,jelem)= tetra_nodes(1:4)
        lnods_material(jelem) = imoh
     end do

  end do

  write(31,*) 'END_ELEMENTS'
  close(31)

  write(6,*) 'Total nodes read:',npoin_total 
  write(6,*) 'Active nodes:    ',npoin 
  write(6,*) 'Active nodes:    ',nelem,jelem 
  allocate(coord(ndime,npoin))

  open (30 ,file='coord.dom.geo')  
  if (iffibers == 1) open (33 ,file='fakefiber.dom.geo')
  fiber_nodes=0.0_rp
  fiber_nodes(3)=1.0_rp
  write(30,*) 'COORDINATES, NOTSORTED'
  do ipoin_total= 1,npoin_total
     read(22,*) idumy, coord_nodes(1:3)
     ipoin = kpermu_nodes(ipoin_total)
     if (ipoin > 0) then
        write(30,200) ipoin,coord_nodes(1:3)
        if (iffibers == 1) write(33,200) ipoin,fiber_nodes(1:3)
        coord(1:ndime,ipoin) = coord_nodes(1:3)
     end if
  end do
  write(30,*) 'END_COORDINATES'
  close(30)
  if (iffibers == 1) close(33)


  nboun_final= 0

  if (ifbocon == 1) then  !! write boundary conditions when required
     
     open (32 ,file='bocon.dom.geo')

     !
     ! bridge a alya que te calcula las faces
     !
     write(6,*) 'Compute faces... '
     call alyafaces
     write(6,*) 'Compute faces... done.'
     
     !
     ! compute nboun_all
     !
     kboun= 0
     do ielem= 1,nelem
        ielty= abs(ltype(ielem))
        do iface_local= 1,nface(ielty)
           iface= lelfa(ielem)%l(iface_local)        
           if (lfacg(2,iface) == 0) then
              ! boundary face identified
              kboun = kboun + 1
           end if
        end do
     end do
     nboun_all = kboun
     
     !
     ! compute lboun_global(1:4,kboun)
     !
     
     allocate(lboun_global(3+2,nboun_all))  ! tiene dos posiciones mas para el elemento de la boundary y para su label
     kboun= 0
     do ielem= 1,nelem
        ielty= abs(ltype(ielem))
        do iface_local= 1,nface(ielty)
           iface= lelfa(ielem)%l(iface_local)        
           if (lfacg(2,iface) == 0) then
              ! boundary face identified
              kboun = kboun + 1
              nnodc = nnodf(ielty) % l(iface_local)
              do inodc= 1,nnodc            
                 inode = lface(ielty)%l(inodc,iface_local) 
                 ipoin = lnods(inode,ielem)
                 lboun_global(inodc,kboun) = ipoin              
              end do
              lboun_global(nnodc+1,kboun)  = ielem              
           end if
        end do
     end do
     write(6,*) 'Compute lboun... done.'
     
  !
  ! escribe las coordenadas rotadas: DESCOMENTAR ESTO SI QUERES ROTAR!!!
  !

!!  open (30 ,file='rot-coord.dom.geo')
!!  write(30,*) 'COORDINATES'
!!  do ipoin=1,npoin
!!     vecaux(1:ndime) = coord(1:ndime,ipoin)  
!!     call cropro(rot_ax,vecaux,rotvec)  
!!     ! la formula de olinde rodriguez
!!     rotvec = vecaux * cosang + rotvec * sinang + rot_ax * dotpro(rot_ax,vecaux) * (1.0_rp - cosang)
!!     coord(1:ndime,ipoin) = rotvec(1:ndime)
!!     write(30,200) ipoin,coord(1:ndime,ipoin)
!!  end do
!!  write(30,*) 'END_COORDINATES'
!!  close(30)
  

  ! ahora, escribe las boundaries
  ! tambien escribe lconditions pero para el caso "csv"
     
     write(6,*) 'Write boundaries and compute boundary conditions...'
     
     nboun_final= 0
     open (30 ,file='bounda.dom.geo')
     write(30,*) 'BOUNDARIES, ELEMENTS'
     
     write(32,*) 'BOUNDARY_CONDITIONS, EXTRAPOLATE' 
     write(32,*) 'ON_BOUNDARIES, UNKNOWN'
     
     kbocondi= 0
     do kboun= 1,nboun_all
        iwrite= 1
        do inode=1,3
           ipoin= lboun_global(inode,kboun)
           
!!! DESCOMENTAR SI QUERES PONER UNA COORDENADA PARA QUE ESCRIBA LAS FRONTERAS POR DEBAJO
           !!        if (coord(3,ipoin) > -45.0_rp) iwrite= 0
           
        end do
        
        if (iwrite > 0) then
           
           nboun_final = nboun_final + 1
           write(30,100) nboun_final, lboun_global(1:3,kboun), lboun_global(4,kboun) 
           
!!! DESCOMENTAR SI QUERES EL FAMOSO CASO CSV
           ! chequea la boun mirando en que lista esta el nodo 3 del kboun
           !!        do icucu= 1,5316
           !!           v(1:3) = points_irregu(1:3,icucu) - coord(1:3,ipoin)
           !!           modve= sqrt(dotpro(v,v))
           !!           if (modve .lt. 1.0e-10) then
           !!              iwrite= 2 
           !!              cycle
           !!           end if
           !!        end do
           !!        
           !!        if (iwrite == 1) then
           !!           do icucu= 1,13442
           !!              v(1:3) = points_regu(1:3,icucu) - coord(1:3,ipoin)
           !!              modve= sqrt(dotpro(v,v))
           !!              if (modve .lt. 1.0e-10) then
           !!                 iwrite= 3 
           !!                 cycle
           !!              end if
           !!           end do
           !!        end if
           
           kbocondi(iwrite) = kbocondi(iwrite) + 1
           
           write(32,*) nboun_final,iwrite
           
        end if
     end do
     
     write(32,*) 'END_ON_BOUNDARIES'
     write(32,*) 'END_BOUNDARY_CONDITIONS' 

     write(30,*) 'END_BOUNDARIES'
     
  end if

  close(30)
  close(32)

  write(6,*) 'Active boundaries  : ',nboun_final 
  write(6,*) 'Boundary conditions: ',kbocondi, kbocondi(1)+kbocondi(2)+kbocondi(3)


!  nboun_final= 0
!  open (30 ,file='bounda.dom.geo')
!  write(30,*) 'BOUNDARIES, ELEMENTS'
!  do kboun= 1,nboun_all
!     iwrite= 0
!     do inode=1,3
!        ipoin= lboun_global(inode,kboun)
!        if (coord(1,ipoin) > 100.0_rp) iwrite= 1
!     end do
!     if (iwrite ==1) then
!        nboun_final = nboun_final + 1
!        write(30,100) nboun_final, lboun_global(1:3,kboun), lboun_global(4,kboun) 
!     end if
!  end do
!  write(30,*) 'END_BOUNDARIES'
!  close(30)
!  write(6,*) 'Active boundaries:',nboun_final 


  !
  ! escribe el dom.dat
  ! 
  open (40 ,file='iris2alya.dom.dat')
  write (40,*) '$------------------------------------------------------------'
  write (40,*) 'DIMENSIONS'
  write (40,*) '  NODAL_POINTS= '     , npoin
  write (40,*) '  ELEMENTS=     '     , nelem
  write (40,*) '  SPACE_DIMENSIONS=               3'
  write (40,*) '  NODES=                          4'
  write (40,*) '  BOUNDARIES=        ', nboun_final
  write (40,*) '  SKEW_SYSTEMS=            0 '
  write (40,*) '  SLAVES=                  0 '
  write (40,*) '  NO_SETS '
  write (40,*) 'END_DIMENSIONS '
  write (40,*) '$------------------------------------------------------------'
  write (40,*) 'STRATEGY'
  write (40,*) '  INTEGRATION_RULE:          CLOSED'
  write (40,*) '  DOMAIN_INTEGRATION_POINTS: 1'
  write (40,*) 'END_STRATEGY'
  write (40,*) '$------------------------------------------------------------'
  write (40,*) 'GEOMETRY'
  if (ifmaterzones == 1) then
     write (40,*) 'MATERIALS, NUMBER=',nmoh
     write (40,*) '   INCLUDE materials.dom.geo'
     write (40,*) 'END_MATERIALS'
  end if
  if (iffibers == 1) then
     write (40,*) '   FIELDS, NUMBER = 1'
     write (40,*) '        FIELD=1, DIMENSION=3, NODES '
     write (40,*) '             INCLUDE  fakefiber.dom.geo'
     write (40,*) '        END_FIELD'
     write (40,*) '   END_FIELDS '
  end if
  write (40,*) '   INCLUDE  eleme.dom.geo'
  write (40,*) '   INCLUDE  coord.dom.geo'
  write (40,*) '   INCLUDE  bounda.dom.geo'
!  write (40,*) '   BOUNDARIES'
!  write (40,*) '   END_BOUNDARIES'
  write (40,*) '   SKEW_SYSTEMS '
  write (40,*) '   END_SKEW_SYSTEMS'  
  write (40,*) 'END_GEOMETRY'
  write (40,*) 'SETS'
  write (40,*) 'END_SETS'
  write (40,*) 'BOUNDARIES'
  write (40,*) 'END_BOUNDARIES'  
  close(40)

  ! finalmente, escribe un ficherito con los materiales

!!$!  open (30 ,file='materials.dom.geo')
!!$!  jelem= 0
!!$!  do ielem= 1,nelem
!!$!     iwrite= 1
!!$!     do inode=1,4
!!$!        ipoin= lnods(inode,ielem)
!!$!        if (coord(3,ipoin) < -45.0_rp) iwrite= 2
!!$!     end do
!!$!     write(30,100) ielem, iwrite
!!$!  end do
!!$!  close(30)


  if (ifmaterzones == 1) then
     
     open (30 ,file='materials.dom.geo')

     do ielem=1,nelem
        write (30,*) ielem,lnods_material(ielem)
     end do
     
     close(30)

  end if
     
  write(6,*) 'Finished OK.'



!-------------------------------------
  stop
  !---------------------------------------------
  !
  ! Esto que esta aca abajo no sirve, pero lo dejo por las dudas...
  !
  !
  ! calcular lboun_global pero con numeracion de nodos local a solo la superficie
  ! lo guardo en lnods para luego llamar a alyafaces, total no necesito mas a lnods
  !
  deallocate(lnods)
  nelem= nboun_all
  mnode= 3
  allocate(lnods(mnode  ,nelem))
  lnnod= 3      ! ahora son todos triangulos porque hara la relacion de vecindad para la superficie

  kpermu_nodes= 0
  npoin_surfa= 0
  do kboun= 1,nboun_all
     do inodc= 1, 3
        if (kpermu_nodes(lboun_global(inodc,kboun))==0) then
           npoin_surfa= npoin_surfa + 1
           kpermu_nodes(lboun_global(inodc,kboun)) = npoin_surfa
        end if
        lnods(inodc,kboun)= kpermu_nodes(lboun_global(inodc,kboun))
     end do     
  end do

  npoin= npoin_surfa
  mnodb= 2
  nnodc= 2
  ndime= 2

  !
  ! llamas de nuevo alyafaces pero ahora para los triangulitos de la superficie
  !
  call alyafaces

  ielem_seed(1)= 4

  allocate(markel(nelem))
  allocate(lstack(nelem))

  markel = 0
  lstack= 0

  lstack(1)= ielem
  nstack= 1
  markel(ielem) = 1

  istack = 0

  do 
     if( istack == nstack ) exit
     istack = istack+1   
     ielem  = lstack(istack)     
     ielty= abs(ltype(ielem))

     do iface_local= 1,nface(ielty)
        iface= lelfa(ielem)%l(iface_local)        
        jelem = lfacg(2,iface)

        if (jelem .ne. 0) then
           if (markel(jelem) == 0) then
              markel(jelem) = 1
              nstack = nstack + 1
              lstack(nstack)= jelem
           end if
        end if
     end do

  end do



  100 format(10(2x,i12))
  150 format(30(2x,i12))
  200 format(i8,3(2x,e15.8))
  
end program iris2alya

!--------------------------------------------

subroutine cropro(a,b,c)

  use def_parame
  implicit none

  real(rp)  :: c(3),a(3),b(3)

  c(1)=   a(2)*b(3) - a(3)*b(2) 
  c(2)= - a(1)*b(3) + a(3)*b(1) 
  c(3)=   a(1)*b(2) - a(2)*b(1) 

end subroutine  cropro

function dotpro(a,b)
  use def_parame
  implicit none

  real(rp)  :: dotpro,a(3),b(3)

  dotpro= a(1) * b(1) + a(2) * b(2) + a(3) * b(3) 

end function dotpro

subroutine parsecsv(buffi,v)

  use def_parame
  use def_elmtyp
  use def_master

  integer(ip) :: icount,iline,icoma(10),j,i,k
  real(rp) :: xdumy,v(3)
  character(1000) :: buffi
  
  j=0
  do i=1,100
     if (buffi(i:i) == ",") then
        j=j+1                  
        icoma(j)= i
     end if
  end do
  
  read(buffi(icoma(1)+1:icoma(2)-1),*) v(1)
  read(buffi(icoma(2)+1:icoma(3)-1),*) v(2)
  read(buffi(icoma(3)+1:icoma(3)+10),*) v(3)
!!!  write(6,*) v(1:3)
  

end subroutine parsecsv
