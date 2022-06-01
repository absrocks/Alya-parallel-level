!  cd /gpfs/projects/bsc21/WORK-HERBERT/repsol/progs/reagrid/ && ./com.sh

!  cd /Users/howen/repsol/progs/reagrid/ && ./com.sh

subroutine reagrid(coord,zcorn,ndi,ndime,nzerom,nz1,nz2)
  use def_kintyp
  implicit none

  integer(ip),intent(in)    :: ndi(3),ndime
  real(rp) ,intent(out)     :: coord(ndime,ndi(1)+1_ip,ndi(2)+1_ip,2)
  real(rp) ,intent(out)     :: zcorn(ndi(1),ndi(2),ndi(3),8)
  integer(ip),intent(in)    :: nzerom
  integer(ip),intent(in)    :: nz1,nz2

  character(120)            :: linea
  character(10)             :: wstr,wstr2
  integer(ip)               :: i,j,k,l,m,k_aux,j_aux,i_aux,l_aux,idime,ksten,i_bot_top,in_st,inode
  integer(ip)               :: ipoin,kadit,ielem,npoin,nelem,kdiff,ifoun,kauxi,kount,npnew
  integer(ip)               :: maxp,maxi,wtype,curwn,curwi,curwj,curwk,curwk2
  real(rp)                  :: zdiff,zfrac,toler,zzz(8),tol_ntg,tol_perm,tol_poro
  real(rp),allocatable      :: co_alya(:,:)
  real(rp),allocatable      :: co_aux(:,:,:,:),zcor2(:)
  real(rp),allocatable      :: perm(:,:,:,:),ntg(:,:,:),poro(:,:,:),swat(:,:,:)
  integer(ip),allocatable   :: welli(:,:,:)
  integer(ip),allocatable   :: kfl_activ(:),lnods(:,:)
  integer(ip),allocatable   :: pnumb(:,:,:,:),lpopo(:)
  integer(ip)               :: kfl_addit(8,3)
  integer(ip)               :: kfl_shift(8)
  integer(ip)               :: lauxi(8)

  integer(ip)               :: nelemtomod,kprism(2,6)
  integer(ip),allocatable   :: kzeromn(:)
  integer(ip),allocatable   :: kelemtomod(:)
  integer(ip),allocatable   :: lnadit(:,:)
  integer(ip),allocatable   :: elem_ijk(:,:), poin_ijk(:,:)
  integer(ip),allocatable   :: iel_orig(:)

  toler=1e-2

  kfl_addit(1,:)=(/0,0,0/)
  kfl_addit(2,:)=(/1,0,0/)
  kfl_addit(3,:)=(/0,1,0/)
  kfl_addit(4,:)=(/1,1,0/)
  kfl_addit(5,:)=(/0,0,1/)
  kfl_addit(6,:)=(/1,0,1/)
  kfl_addit(7,:)=(/0,1,1/)
  kfl_addit(8,:)=(/1,1,1/)

  print*,'nzerom',nzerom
  i=max(1,nzerom)
  print*,'i',i
  allocate(kzeromn(i))

  allocate(lnods(8,ndi(1)*ndi(2)*ndi(3)))
  allocate(kelemtomod(ndi(1)*ndi(2)*ndi(3)))
  allocate(kfl_activ(ndi(1)*ndi(2)*ndi(3)))

  allocate(perm(ndi(1),ndi(2),ndi(3),3))
  allocate(ntg(ndi(1),ndi(2),ndi(3)))
  allocate(poro(ndi(1),ndi(2),ndi(3)))
  allocate(swat(ndi(1),ndi(2),ndi(3)))
  allocate(welli(ndi(1)+1,ndi(2)+1,ndi(3)+1))

  allocate(elem_ijk(ndi(1)*ndi(2)*ndi(3),3))
  allocate(poin_ijk(2*(ndi(1)+1)*(ndi(2)+1)*(ndi(3)+1),3)) ! Twice as more since it is unknown yet how many points will be added

  k_aux = ndi(1)/4_ip               ! Number of full lines with 4 numbers
  j_aux = ndi(1)-k_aux*4_ip         ! Last line of the block can have less than 4 numbers


  !
  ! Porosity file
  !
  open(20,file='poro.dat')

  read(20,'(a)')linea
  do k=1,ndi(3)
     do j=1,ndi(2)
        do i_aux=1,k_aux
           read(20,*) (poro(l_aux,j,k), l_aux = (i_aux-1_ip)*4_ip + 1_ip, (i_aux-1_ip)*4_ip + 4_ip)
        end do
        read(20,*) (poro(i_aux,j,k), i_aux = k_aux*4_ip + 1_ip, k_aux*4_ip + j_aux)
     end do
  end do

  close(20)


  !
  ! Net-to-gross file, this should be multiplied by porosity
  !
  open(20,file='ntg.dat')

  read(20,'(a)')linea
  do k=1,ndi(3)
     do j=1,ndi(2)
        do i_aux=1,k_aux
           read(20,*) (ntg(l_aux,j,k), l_aux = (i_aux-1_ip)*4_ip + 1_ip, (i_aux-1_ip)*4_ip + 4_ip)
        end do
        read(20,*) (ntg(i_aux,j,k), i_aux = k_aux*4_ip + 1_ip, k_aux*4_ip + j_aux)
     end do
  end do

  close(20)


  !
  ! Water saturation file
  !
  open(20,file='swat.dat')

  read(20,'(a)')linea
  do k=1,ndi(3)
     do j=1,ndi(2)
        do i_aux=1,k_aux
           read(20,*) (swat(l_aux,j,k), l_aux = (i_aux-1_ip)*4_ip + 1_ip, (i_aux-1_ip)*4_ip + 4_ip)
        end do
        read(20,*) (swat(i_aux,j,k), i_aux = k_aux*4_ip + 1_ip, k_aux*4_ip + j_aux)
     end do
  end do

  close(20)


  !
  ! Permeability file, has 3 components - X, Y and Z
  !
  open(20,file='perm.dat')

  do idime=1,3
     readline: do
        read(20,'(a)')linea
        l_aux=index(linea,'PERM')
!        print*,'l_aux,linea',l_aux,linea
        if ((l_aux>0).and.(l_aux<5)) exit readline      ! When linea is 'PERMX', 'PERMY' or 'PERMZ'
     end do readline
     do k=1,ndi(3)
        do j=1,ndi(2)
           do i_aux=1,k_aux
              read(20,*) (perm(l_aux,j,k,idime), l_aux = (i_aux-1_ip)*4_ip + 1_ip, (i_aux-1_ip)*4_ip + 4_ip)
           end do
           read(20,*) (perm(i_aux,j,k,idime), i_aux = k_aux*4_ip + 1_ip, k_aux*4_ip + j_aux)
        end do
     end do
  end do

  close(20)


!  print*,'perm(ndi(1)-2,ndi(2),ndi(3),:)',perm(ndi(1)-2,ndi(2),ndi(3),:)


  !
  ! Binary (0-1) file
  !
  open(12,file='actnum.dat')

  k_aux=ndi(1)*ndi(2)/25_ip                   ! Now by blocks of 25 numbers
  j_aux=(ndi(1)*ndi(2))-(k_aux*25_ip)

  read(12,'(a)')linea
  do k=1,ndi(3)
     kadit = (k-1)*(ndi(1)*ndi(2))
     do i_aux=1,k_aux
        read(12,*) (kfl_activ(l_aux),l_aux = kadit+(i_aux-1_ip)*25_ip+1_ip, kadit+(i_aux-1_ip)*25_ip+25_ip)

!!!! This is probably for debugging? File 13 is not open here
!        do l_aux = kadit+(i_aux-1_ip)*25_ip+1_ip, kadit+(i_aux-1_ip)*25_ip+25_ip
!           write(13,*)'l_aux,kfl_activ(l_aux)',l_aux,kfl_activ(l_aux)
!        end do

     end do
     read(12,*) (kfl_activ(i_aux), i_aux = k_aux*25_ip+kadit+1_ip, k_aux*25_ip+kadit+j_aux)
  end do

  close(12)


  ! Elimination of very small values
  tol_ntg = 1.0e-10_rp
  tol_poro = 1.0e-10_rp
  tol_perm = 1.0e-10_rp

  do k=1,ndi(3)
     do j=1,ndi(2)
        do i=1,ndi(1)
           if ( abs(ntg(i,j,k)) < tol_ntg ) kfl_activ(i+(j-1)*ndi(1)+(k-1)*ndi(1)*ndi(2)) = 0
           if ( abs(poro(i,j,k)) < tol_poro ) kfl_activ(i+(j-1)*ndi(1)+(k-1)*ndi(1)*ndi(2)) = 0
           if ( (abs(perm(i,j,k,1)) < tol_perm) .and.(abs(perm(i,j,k,2)) < tol_perm) .and.(abs(perm(i,j,k,3)) < tol_perm) ) &
                                                                                   kfl_activ(i+(j-1)*ndi(1)+(k-1)*ndi(1)*ndi(2)) = 0
        end do
     end do
  end do
  
  
  !
  ! Grid file
  !
  open(10,file='grid.dat')
  do i=1,23
     read(10,'(a)')linea
  end do
 
  do j=1,ndi(2)+1_ip
     do i=1,ndi(1)+1_ip
        read(10,*) ((coord(idime,i,j,ksten),idime=1,ndime),ksten=1,2)
     end do
!     read(10,*)linea ! no entiendo bien porque pero si pongo esta linea va mal
  end do

  do i=1,2
     read(10,*)linea
  end do

  k_aux=ndi(1)/6_ip
  j_aux=ndi(1)-(k_aux*6_ip)

  do k=1,ndi(3)
     do i_bot_top=1,2
        do j=1,ndi(2)

           in_st=1+((i_bot_top-1)*4)
           do i_aux=1,k_aux
              read(10,*) ((zcorn(6*(i_aux-1_ip)+l_aux,j,k,inode),inode=in_st,in_st+1),l_aux=1,6_ip)
           end do
           read(10,*) ((zcorn(6*k_aux+i_aux,j,k,inode),inode=in_st,in_st+1),i_aux=1,j_aux)
           
!           read(10,*)linea

           in_st=3+((i_bot_top-1)*4)
           do i_aux=1,k_aux
              read(10,*) ((zcorn(6*(i_aux-1_ip)+l_aux,j,k,inode),inode=in_st,in_st+1),l_aux=1,6_ip)
           end do
           read(10,*) ((zcorn(6*k_aux+i_aux,j,k,inode),inode=in_st,in_st+1),i_aux=1,j_aux)
!           read(10,*)linea
        end do
     end do
  end do

  nelem = 0
  npoin = 0

  kfl_addit(1,:)=(/0,0,0/)
  kfl_addit(2,:)=(/1,0,0/)
  kfl_addit(3,:)=(/0,1,0/)
  kfl_addit(4,:)=(/1,1,0/)
  kfl_addit(5,:)=(/0,0,1/)
  kfl_addit(6,:)=(/1,0,1/)
  kfl_addit(7,:)=(/0,1,1/)
  kfl_addit(8,:)=(/1,1,1/)

  allocate(co_alya(ndime,ndi(1)*ndi(2)*ndi(3)*8))
  allocate(zcor2(ndi(1)*ndi(2)*ndi(3)*8))
  allocate(pnumb(ndi(1)+1,ndi(2)+1,ndi(3)+1,9))
  allocate(lpopo(ndi(1)*ndi(2)*ndi(3)*8))
  pnumb = 0

  do k=1,ndi(3)
     do j=1,ndi(2)
        do i=1,ndi(1)
           i_aux = i + ndi(1) * (j -1) + ndi(1) * ndi(2) * (k -1)
           if (kfl_activ(i_aux)==1_ip) then
              nelem = nelem + 1
              elem_ijk(nelem,1) = i
              elem_ijk(nelem,2) = j
              elem_ijk(nelem,3) = k
              do inode=1,8
                 npoin = npoin + 1
                 pnumb(i+kfl_addit(inode,1),j+kfl_addit(inode,2),k+kfl_addit(inode,3),9) = &
                      pnumb(i+kfl_addit(inode,1),j+kfl_addit(inode,2),k+kfl_addit(inode,3),9) + 1
                 kauxi = pnumb(i+kfl_addit(inode,1),j+kfl_addit(inode,2),k+kfl_addit(inode,3),9)
                 pnumb(i+kfl_addit(inode,1),j+kfl_addit(inode,2),k+kfl_addit(inode,3),kauxi) = npoin
                 zcor2(npoin) = zcorn(i,j,k,inode)
                 lnods(inode,nelem) = npoin
              end do
           end if
        end do
     end do
  end do
  !
  ! eliminate repeated
  !
  npnew=0
  do k=1,ndi(3)+1
     do j=1,ndi(2)+1
        do i=1,ndi(1)+1
           npnew = npnew + 1
           poin_ijk(npnew,1) = i
           poin_ijk(npnew,2) = j
           poin_ijk(npnew,3) = k
           zzz(1) = zcor2(pnumb(i,j,k,1))
           co_alya(3,npnew) = zzz(1)
           zfrac = (co_alya(3,npnew)-coord(3,i,j,1)) / (coord(3,i,j,2)-coord(3,i,j,1))
           do idime=1,2
              co_alya(idime,npnew) = coord(idime,i,j,1) +  &
                   zfrac * ( coord(idime,i,j,2) - coord(idime,i,j,1) )
           end do

           lpopo( pnumb(i,j,k,1) ) = npnew
           kount = 1
           lauxi(1)=1
           do l = 2,pnumb(i,j,k,9)
              ifoun = 0
              do m = 1,l-1
                 if ( zzz(m) == zcor2(pnumb(i,j,k,l)) ) then   ! that  z already exists
                    lpopo( pnumb(i,j,k,l) ) = lpopo( pnumb(i,j,k,lauxi(m)) )
                    ifoun = 1
                 end if
              end do
              if (ifoun==0) then
                 npnew = npnew + 1
                 poin_ijk(npnew,1) = i
                 poin_ijk(npnew,2) = j
                 poin_ijk(npnew,3) = k
                 kount = kount + 1
                 lauxi(kount)=l
                 zzz(kount) = zcor2(pnumb(i,j,k,l))
                 lpopo( pnumb(i,j,k,l) ) = npnew
                 co_alya(3,npnew) = zzz(kount)
                 zfrac = (co_alya(3,npnew)-coord(3,i,j,1)) / (coord(3,i,j,2)-coord(3,i,j,1))
                 do idime=1,2
                    co_alya(idime,npnew) = coord(idime,i,j,1) +  &
                         zfrac * ( coord(idime,i,j,2) - coord(idime,i,j,1) )
                 end do
              end if
           end do
        end do
     end do
  end do
     
  open(11,file='alya.geo.dat')


  write(11,'(a)') 'COORDINATES'
  do ipoin=1,npnew
     write(11,'(i9,3(2x,e17.10))'),ipoin,(co_alya(idime,ipoin),idime=1,3)
  end do
  write(11,'(a)') 'END_COORDINATES'

  kfl_shift(1)=1
  kfl_shift(2)=2
  kfl_shift(3)=4
  kfl_shift(4)=3
  kfl_shift(5)=5
  kfl_shift(6)=6
  kfl_shift(7)=8
  kfl_shift(8)=7


  !
  ! Zero mass file
  !
  open(15,file='zeromassn')

  do i=1,nz1
     read(15,'(a)')linea
     read(linea(40:),*)kzeromn(i)
  end do
  
! ojo no todos los nods los agarraba como zero mass -- solo los de un lado
! los otros los tuve que agarra a mano con paravie los pongo en xx
!Lo s nodos  de zero mass son los que lueg indican que hay que partir el hexa en prismas
  close(15)


  !
  ! xx(??) file
  !
  open(15,file='xx')

  do i=1,nz2
     read(15,'(a)')linea
     k=index(linea,',')
     read(linea(:k-1),*)j
     kzeromn(i+nz1)=j+1
     print*,'i+nz1,kzeromn(i+nz1)',i+nz1,kzeromn(i+nz1)
  end do
  
  close(15)


  kelemtomod = 0
  do ielem=1,nelem
     do inode=1,8
        do i=1,nzerom
           if(lpopo(lnods(kfl_shift(inode),ielem)) == kzeromn(i)) then
              print*,'i,kzeromn(i),ielem,inode',i,kzeromn(i),ielem,inode
              kelemtomod(ielem) = 1
           end if
        end do
     end do
  end do

! la forma de cojer los nodos raros ha sido un tanto trucha pero solo valida para esta malla para
! el caso gral hay que pensarlo un poco 


  kprism(1,:)=(/1,2,4,5,6,8/)
  kprism(2,:)=(/2,3,4,6,7,8/)

  allocate(lnadit(6,max(1,nzerom)))
  allocate(iel_orig(max(1,nzerom)))
  nelemtomod = 0
  write(11,'(a)') 'ELEMENTS'
  do ielem=1,nelem
     if (kelemtomod(ielem) == 0) then
        write(11,'(9(i9,2x))'),ielem,(lpopo(lnods(kfl_shift(inode),ielem)),inode=1,8)
     else
        write(11,'(9(i7,2x))'),ielem,(lpopo(lnods(kfl_shift(kprism(1,inode)),ielem)),inode=1,6)
        nelemtomod = nelemtomod + 1
        iel_orig(nelemtomod) = ielem
        do inode=1,6
           lnadit(inode,nelemtomod) = lpopo(lnods ( kfl_shift( kprism(2,inode) ) ,ielem) )
        end do
     end if
  end do
  do ielem=1,nelemtomod
     write(11,'(9(i7,2x))'),ielem+nelem,(lnadit(inode,ielem),inode=1,6)
  end do
  write(11,'(a)') 'END_ELEMENTS'

  open(15,file='types.dat')

  write(15,'(a)') 'TYPES'
  do ielem=1,nelem
     if (kelemtomod(ielem) == 0) then
        write(15,'(2(i8,2x))'),ielem,37
     else
        write(15,'(2(i8,2x))'),ielem,34
     end if
  end do
  do ielem=1,nelemtomod
     write(15,'(2(i8,2x))'),ielem+nelem,34
  end do
  write(15,'(a)') 'END_TYPES'

  close(11)
  close(15)

  print*,'npnew',npnew
  print*,'nelem+nelemtomod',nelem+nelemtomod


  ! Write porosity, permeability and water saturation into Alya-style txt files
  open(14,file='poro.txt')

! write(14,'(a)') 'POROSITY'
  do ielem=1,nelem
     write(14,'(i9,(2x,e17.10))'), ielem, poro(elem_ijk(ielem,1),elem_ijk(ielem,2),elem_ijk(ielem,3)) * &
                                            ntg(elem_ijk(ielem,1),elem_ijk(ielem,2),elem_ijk(ielem,3))
  end do
  do ielem=1,nelemtomod
     write(14,'(i9,(2x,e17.10))'), nelem+ielem, poro(elem_ijk(iel_orig(ielem),1),elem_ijk(iel_orig(ielem),2), & 
                                   elem_ijk(iel_orig(ielem),3)) * &
                                   ntg(elem_ijk(iel_orig(ielem),1),elem_ijk(iel_orig(ielem),2),elem_ijk(iel_orig(ielem),3))
  end do
! write(14,'(a)') 'END_POROSITY'
  
  close(14)


  open(14,file='perme.txt')

! write(14,'(a)') 'PERMEABILITY'
  do ielem=1,nelem
     write(14,'(i9,3(2x,e17.10))'), ielem, perm(elem_ijk(ielem,1),elem_ijk(ielem,2),elem_ijk(ielem,3),1:ndime)
  end do
  do ielem=1,nelemtomod
     write(14,'(i9,3(2x,e17.10))'), nelem+ielem, perm(elem_ijk(iel_orig(ielem),1),elem_ijk(iel_orig(ielem),2), &
                                    elem_ijk(iel_orig(ielem),3),1:ndime)
  end do
! write(14,'(a)') 'END_PERMEABILITY'
  
  close(14)


  open(14,file='swat.txt')

! write(14,'(a)') 'WATER_SATURATION'
  do ielem=1,nelem
     write(14,'(i9,(2x,e17.10))'), ielem, swat(elem_ijk(ielem,1),elem_ijk(ielem,2),elem_ijk(ielem,3))
  end do
  do ielem=1,nelemtomod
     write(14,'(i9,(2x,e17.10))'), nelem+ielem, swat(elem_ijk(iel_orig(ielem),1),elem_ijk(iel_orig(ielem),2), & 
                                   elem_ijk(iel_orig(ielem),3))
  end do
! write(14,'(a)') 'END_WATER_SATURATION'
  
  close(14)


  ! Read a well file
  open(25,file='brugge.wel')
  welli(:,:,:) = 0

  ! Reads until 'COMPDAT' is found, then scans well lines
  ! First variable: P## is a producer and I## is an injector
  ! Then I, J and K1 are used to determine the node
  maxp = 0
  maxi = 0
  do i=1,2
     readdummy: do
        read(25,'(a)') linea
        l_aux = index(linea,'COMPDAT')
        if ((l_aux>0_ip) .and. (l_aux<5_ip)) exit readdummy
     end do readdummy

     read(25,'(a)') linea
     read(25,'(a)') linea

     readwells: do
        read(25,'(a)') linea
        read(linea,*) wstr,curwi,curwj,curwk,curwk2
        l_aux = index(linea,'/')
        if ((l_aux>0_ip) .and. (l_aux<2_ip)) exit readwells    ! Last line starts with '/'
        
        wstr2 = wstr(2:len(wstr))
        read(wstr2,*) curwn

        if (wstr(1:1) == 'P') then
          wtype = 1
          if (curwn > maxp) maxp = curwn
!            welli(curwi,curwj,curwk) = curwn
          do k=curwk,curwk2
             welli(curwi,curwj,k) = curwn
          end do
        elseif (wstr(1:1) == 'I') then
          wtype = 2
          if (curwn > maxi) maxi = curwn
!             welli(curwi,curwj,curwk) = -(maxp+curwn)            ! We suppose that at this moment maxp is already known
          do k=curwk,curwk2
              welli(curwi,curwj,k) = -(maxp+curwn)            ! We suppose that at this moment maxp is already known
          end do
        else
          write(*,*) 'Unknown type!', wstr
        endif

!        write(*,*) wtype,curwn,curwi,curwj,curwk

                
     end do readwells
  end do

  close(25)

  write(*,*) 'Max # of the productors:', maxp, ', max # of the injectors:', maxi


  ! Write a well file
  open(26,file='well.txt')

! write(26,'(a)') 'WELLS'
  do ipoin=1,npnew
     write(26,'(i9,i9)'), ipoin, iabs(welli(poin_ijk(ipoin,1),poin_ijk(ipoin,2),poin_ijk(ipoin,3)) )
  end do
! write(26,'(a)') 'END_WELLS'

  close(26)

  
end subroutine reagrid
