!-----------------------------------------------------------------------
!> @addtogroup Dodeme
!> @{
!> @file    dod_merge_new_elements.f90
!> @author  Guillaume Houzeaux
!> @date    18/09/2012
!> @brief   Create 2D extensions
!> @details Create extension elements for 2D meshes
!>          What do we have for ISUBD:
!>          - SUBDOMAIN(ISUBD) % LBOCH(IBOUN)       =  BOEXT on interface
!>          - SUBDOMAIN(ISUBD) % LSUBD_NPOIN(IPOIN) =  JSUBD on interface
!>          - SUBDOMAIN(ISUBD) % LSUBD_NPOIN(IPOIN) = -JSUBD in holes
!>          - SUBDOMAIN(ISUBD) % LSUBD_NELEM(IELEM) = -JSUBD in holes
!> @} 
!-----------------------------------------------------------------------
subroutine dod_merge_new_elements()
  use def_kintyp
  use def_elmtyp
  use def_parame
  use def_master
  use def_elmtyp
  use def_domain
  use def_dodeme
  use mod_memory
  use mod_dod_extens
  use mod_outfor, only : outfor
  use mod_messages, only : livinf
  implicit none  
  integer(ip)          :: ielem,pnode,inode,kpoin,pelty,ipoin,isubd
  integer(ip)          :: izone,kelem,iboun,jpoin,jelem
  integer(ip)          :: ipoin_global,ielem_global,iboun_global
  integer(ip), pointer :: lmate_cpy(:)
  integer(ip), pointer :: lesub_cpy(:)
  integer(ip)          :: opener,inodb,idime,istat
  integer(ip)          :: ielem_max_kap,ielem_max_q,ielem_min_kap,ielem_min_q,ielem_min_asr,ielem_max_asr,kk,sign
  real(rp)             :: kappa_avg ,kappa_max ,kappa_min ,q_max ,q_avg ,q_min,asrad_max,asrad_min,asrad_avg,kappa,asrad,q,dummi(3)
!integer(ip)          :: permu(35432),ii
  call livinf(0_ip,'MERGE NEW MESH WITH EXTENSIONS WITH ORIGINAL ONE',0_ip)

  !
  ! Nullify
  !
  nullify(lmate_cpy)
  nullify(lesub_cpy)
  !
  ! New mesh dimensions
  !
  nelem_old = nelem
  nelem     = nelem + nelem_dod
  !
  ! Resize arrays
  !  
  call memory_resize(mem_servi(1:2,servi),'LTYPE','dod_merge_new_elements',ltype,nelem)
  call memory_resize(mem_servi(1:2,servi),'LNNOD','dod_merge_new_elements',lnnod,nelem)
  call memory_resize(mem_servi(1:2,servi),'LELCH','dod_merge_new_elements',lelch,nelem)
  call memory_resize(mem_servi(1:2,servi),'LNODS','dod_merge_new_elements',lnods,mnode,nelem) 
  call memory_resize(mem_servi(1:2,servi),'LESUB','dod_merge_new_elements',lesub,nelem) 
  !if( nmate > 1 ) call memory_resize(mem_servi(1:2,servi),'LMATE','dod_merge_new_elements',lmate,nelem)

  if( nmate > 0 ) then
     call memory_alloca( mem_servi(1:2,servi),'LMATE_CPY','dod_merge_new_elements',lmate_cpy,nelem_old)
     do ielem = 1,nelem_old
        lmate_cpy(ielem) = lmate(ielem)
     end do
     call memory_deallo (mem_servi(1:2,servi),'LMATE','dod_merge_new_elements', lmate )
     call memory_alloca( mem_servi(1:2,servi),'LMATE_CPY','dod_merge_new_elements',lmate,nelem)   
     do ielem = 1,nelem_old
        lmate(ielem) = lmate_cpy(ielem)
     end do
  end if
  if( neset > 0 ) call memory_resize(mem_servi(1:2,servi),'LESET','dod_merge_new_elements',leset,nelem)
  !
  ! Reallocate zones
  !
  call memory_copy(  mem_servi(1:2,servi),'LESUB_CPY','dod_merge_new_elements',lesub,lesub_cpy)
  call memory_alloca(mem_servi(1:2,servi),'LESUB',    'dod_merge_new_elements',lesub,nelem) 
  ielem = nelem_old
  !
  ! Merge new mesh
  !
  ielem = nelem_old
  if( associated(lpext) ) then
     do kpoin = 1,size(lpext)
        do kelem = 1,size( lpext(kpoin) % ltype )
           ielem = ielem + 1
           pelty = lpext(kpoin) % ltype(kelem)
           pnode = nnode(pelty)
           do inode = 1,pnode
              lnods(inode,ielem) = lpext(kpoin) % lnods(inode,kelem)
           end do
           ltype(ielem) = pelty
           lnnod(ielem) = nnode(pelty)
           lelch(ielem) = ELEXT
           if( nmate > 0 ) then
              lmate(ielem) = lpext(kpoin) % lmate(kelem)
           end if
           if( neset > 0 ) then
              leset(ielem) = 0
           end if
           lesub(ielem) = lpext(kpoin) % lesub(kelem)
        end do
     end do
  end if
  !
  ! Hole nodes
  !
  do isubd = 1,nsubd
     do ipoin = 1,subdomain(isubd) % npoin
        if( subdomain(isubd) % lsubd_npoin(ipoin) < 0 ) then
           ipoin_global = subdomain(isubd) % lnper(ipoin)
           lnoch(ipoin_global) = NOHOL
        end if
     end do
  end do
  !
  ! Hole elements
  !
  do isubd = 1,nsubd
     do ielem = 1,subdomain(isubd) % nelem
        if( subdomain(isubd) % lsubd_nelem(ielem) < 0 ) then
           ielem_global = subdomain(isubd) % leper(ielem)
           lelch(ielem_global) = ELHOL
           ltype(ielem_global) = -abs(ltype(ielem_global))
        end if
     end do
  end do

  !----------------------------------------------------------------------
  !
  ! Boundaries: new boundaries created by holes are not included
  !
  !----------------------------------------------------------------------
  !
  ! Extension boundaries
  ! IBOUN_GLOBAL = 0 for newly created holes (after hole cutting)
  !
  do isubd = 1,nsubd
     do iboun = 1,subdomain(isubd) % nboun
        iboun_global = subdomain(isubd) % lbper(iboun)
        if( iboun_global > 0 ) lboch(iboun_global) = subdomain(isubd) % lboch(iboun)
     end do
  end do

  !----------------------------------------------------------------------
  !
  ! Output info
  !
  !----------------------------------------------------------------------

  coutp(1) = 'NUMBER INTERFACE NODES'
  ioutp(1) = number_fringe_nodes
  call outfor(61_ip,lun_outpu_dod,' ')
  coutp(1) = 'NUMBER EXTENSION ELEMENTS'
  ioutp(1) = nelem_dod
  call outfor(61_ip,lun_outpu_dod,' ')
  coutp(1) = 'PERCENTAGE EXTENSION ELEMENTS'
  routp(1) = 100.0_rp*real(nelem_dod,rp)/real(nelem,rp)
  call outfor(62_ip,lun_outpu_dod,' ')
  !
  ! Count number of equal elements
  !
  if( associated(lpext) ) then
     do kpoin = 1,size(lpext)
        do kelem = 1,size( lpext(kpoin) % ltype )
           pnode = nnode(lpext(kpoin) % ltype(kelem))
           call sortin(pnode,lpext(kpoin) % lnods(1,kelem))
        end do
     end do
     kelem = 0
     kpoin = size(lpext)
     do ipoin = 1,kpoin
        do ielem = 1,size( lpext(ipoin) % ltype )
           pnode = nnode(lpext(ipoin) % ltype(ielem))
           do jpoin = ipoin+1,kpoin
              do jelem = 1,size( lpext(jpoin) % ltype )
                 inode = 0
                 do while( inode < pnode )
                    inode = inode + 1
                    if( lpext(ipoin) % lnods(inode,ielem) /= lpext(jpoin) % lnods(inode,jelem) ) then
                       inode = pnode + 1
                    end if
                 end do
                 if( inode == pnode ) kelem = kelem + 1
              end do
           end do
        end do
     end do
  end if
  coutp(1) = 'NUMBER EQUAL ELEMENTS'
  ioutp(1) = kelem
  call outfor(61_ip,lun_outpu_dod,' ')
  !
  ! Write new geometry
  !
if(1==2)then
  open(unit = 100, file ='name.geo.dat2',  action= 'write', status = 'new' , iostat = opener)
  write(100,'(a)') 'NODES_PER_ELEMENTS'
  do ielem=1,nelem
     write(100,'(i7,1(1x,i2))')ielem,nnode(ltype(ielem))
  end do
  write(100,'(a)') 'END_NODES_PER_ELEMENTS'
 write(100,'(a)') 'COORDINATES'
 do ipoin = 1,npoin
    write(100,'(i7,3(1x,e15.8))') ipoin,(coord(idime,ipoin),idime=1,ndime)
 end do
 write(100,'(a)') 'END_COORDINATES'
 write(100,'(a)') 'ELEMENTS'
 do ielem = 1,nelem
    write(100,'(i7,20(1x,i7))') ielem,(lnods(inode,ielem),inode=1,nnode(ltype(ielem)))
 end do
 write(100,'(a)') 'END_ELEMENTS'
 write(100,'(a)')  'CHARACTERISTIC ELEMENTS'
 do ielem=1,nelem
    write(100,'(i7,1(1x,i1))') ielem,lelch(ielem)
 end do
 write(100,'(a)')  'END_CHARACTERISTIC_LELCH'
 write(100,'(a)') 'BOUNDARIES'
 do iboun = 1,nboun
    write(100,'(i7,20(1x,i7))') iboun,(lnodb(inodb,iboun),inodb=1,nnode(ltypb(iboun)))
 end do
 write(100,'(a)') 'END_BOUNDARIES'
 write(100,'(a)')  'CHARACTERISTIC BOUNDARIES'
 do iboun=1,nboun
    write(100,'(i7,1(1x,i1))') iboun,lboch(iboun)
 end do
 write(100,'(a)')  'END_CHARACTERISTIC_BOUNDARIES'
 !if( nmate > 0 ) then
!  !write(100,'(a)')  'MATERIALS_NUMBER=2_DEFAULT=1'
!  !  do ielem=1,nelem
!  !      if(lmate(ielem)>1)write(100,'(i7,1(1x,i1))')ielem,lmate(ielem)
!  !  end do
  !end if
  !write(100,'(a)')  'END_MATERIALS'
  !  close(100)
end if
!permu = 0
!if(1==1)then
!  open(unit = 100, file ='exts.geo.dat',  action= 'write', status = 'new' , iostat = opener)
!  write(100,'(a)') 'NODES_PER_ELEMENTS'
!  ii = 0
!  do ielem=1,nelem
!     if(lelch(ielem)==ELEXT)then
!        ii = ii + 1
!        write(100,'(i7,1(1x,i2))')ii,nnode(ltype(ielem))
!        do inode=1,nnode(ltype(ielem))
!           permu(lnods(inode,ielem))=1
!        end do
!     end if
!  end do
!  ii = 0
!  do ipoin=1,npoin
!     if(permu(ipoin)==1)then
!        ii = ii + 1
!        permu(ipoin) = ii 
!     end if
!  end do

!  write(100,'(a)') 'END_NODES_PER_ELEMENTS'
! write(100,'(a)') 'COORDINATES'
! do ipoin = 1,npoin
!    if(permu(ipoin)/=0)then
!    write(100,'(i7,3(1x,e15.8))') permu(ipoin),(coord(idime,ipoin),idime=1,ndime)
! end if
! end do
! write(100,'(a)') 'END_COORDINATES'
! write(100,'(a)') 'ELEMENTS'
!ii = 0
! do ielem = 1,nelem
!    if(lelch(ielem)==ELEXT)then
!       ii = ii + 1
!       write(100,'(i7,20(1x,i7))') ii,(permu(lnods(inode,ielem)),inode=1,nnode(ltype(ielem)))
!    end if
! end do
! write(100,'(a)') 'END_ELEMENTS'
!!stop
!end if








if(1==2)then
     kappa_avg = 0.0_rp
     kappa_max = -100.0_rp
     kappa_min = 1.0e15_rp
     q_max = -100.0_rp
     q_avg = 0.0_rp
     q_min = 1.0e15_rp
     asrad_max = -100.0_rp
     asrad_avg = 0.0_rp
     asrad_min = 1.0e15_rp
     ielem_max_kap = 0
     ielem_min_kap = 0
     ielem_max_q = 0
     ielem_min_q = 0 
     ielem_max_asr = 0
     ielem_min_asr = 0
     do ielem=1,nelem
        if(ltype(ielem)>0 .and. lelch(ielem)==ELFEM)then
           kk = kk + 1
           call dod_extens_qual3d(2_ip,lnods(1,ielem),lnods(2,ielem),lnods(3,ielem),lnods(4,ielem),kappa,asrad,q,sign,dummi)    
           kappa_avg = kappa_avg +  kappa
           q_avg = q_avg +  q
           asrad_avg = asrad_avg + asrad
           if(kappa>kappa_max)then
              kappa_max = kappa
              ielem_max_kap = ielem
           end if
           if(q>q_max)then
              q_max = q
              ielem_max_q = ielem
           end if
           if(asrad>asrad_max)then
              asrad_max = asrad
              ielem_max_asr = ielem
           end if
           if(kappa<kappa_min)then
              kappa_min = kappa
              ielem_min_kap = ielem
           end if
           if(q<q_min)then
              q_min = q
              ielem_min_q = ielem
           end if
           if(asrad<asrad_min)then
              asrad_min = asrad
              ielem_min_asr = ielem
           end if
        end if
     end do
     kappa_avg = kappa_avg / kk
     q_avg = q_avg/kk
     asrad_avg = asrad_avg/kk
print*,'estaditica-ELFEM-merge'

     print*,'maximo-kappa',kappa_max,ielem_max_kap,'max-q',q_max,ielem_max_q,'max-asrad',asrad_max,ielem_max_asr
     print*,'promedio-kappa',kappa_avg,'average-q',q_avg,'average-asrad',asrad_avg,'number',kk
     print*,'minimo-kappa',kappa_min,ielem_min_kap,'min-q',q_min,ielem_min_q,'min-asrad',asrad_min,ielem_min_asr
     kappa_avg = 0.0_rp
     kappa_max = -100.0_rp
     kappa_min = 1.0e15_rp
     asrad_avg = 0.0_rp
     asrad_max = -100.0_rp
     asrad_min = 1.0e15_rp
     q_max = -100.0_rp
     q_avg = 0.0_rp
     q_min = 1.0e15_rp
     ielem_max_kap = 0
     ielem_max_q = 0
     ielem_max_asr = 0
     ielem_min_kap = 0
     ielem_min_q = 0
     ielem_min_asr = 0
     kk = 0
     do ielem=1,nelem
        if(ltype(ielem)>0 .and. lelch(ielem)==ELEXT)then
           kk = kk + 1
           call dod_extens_qual3d(2_ip,lnods(1,ielem),lnods(2,ielem),lnods(3,ielem),lnods(4,ielem),kappa,asrad,q,sign,dummi)           
           kappa_avg = kappa_avg +  kappa
           asrad_avg = asrad_avg +  asrad
           q_avg = q_avg +  q
           if(asrad>asrad_max)then
              asrad_max = asrad
              ielem_max_asr = ielem
           end if
           if(kappa>kappa_max)then
              kappa_max = kappa
              ielem_max_kap = ielem
           end if
           if(q>q_max)then
              q_max = q
              ielem_max_q = ielem
           end if
           if(kappa<kappa_min)then
              kappa_min = kappa
              ielem_min_kap = ielem
           end if
           if(q<q_min)then
              q_min = q
              ielem_min_q = ielem
           end if
           if(asrad<asrad_min)then
              asrad_min = asrad
              ielem_min_asr = ielem
           end if
        end if
     end do

     kappa_avg = kappa_avg / kk
     q_avg = q_avg/kk
     asrad_avg = asrad_avg/kk
print*,'estaditica-ELEXT-merge'

     print*,'maximo-kappa',kappa_max,ielem_max_kap
     print*,'promedio-kappa',kappa_avg,kk
     print*,'minimo-kappa',kappa_min,ielem_min_kap
     print*,'maximo-q',q_max,ielem_max_q
     print*,'promedio-q',q_avg,kk
     print*,'minimo-q',q_min,ielem_min_q
    print*,'maximo-asrad',asrad_max,ielem_max_asr
     print*,'promedio-asrad',asrad_avg,kk
     print*,'minimo-asrad',asrad_min,ielem_min_asr
  end if
end subroutine dod_merge_new_elements
