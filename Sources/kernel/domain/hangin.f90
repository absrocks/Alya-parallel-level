!-----------------------------------------------------------------------
!> @addtogroup Turnon
!> @{
!> @file    hangin.f90
!> @author  Guillaume Houzeaux
!> @brief   Hanging nodes
!> @details Detect hanging nodes automatically
!!
!> @} 
!-----------------------------------------------------------------------
subroutine hangin()
  use def_kintyp
  use def_domain
  use def_elmtyp
  use def_master
  use def_kermod
  use mod_elmgeo
  use mod_memory 
  implicit none
  integer(ip)           :: ipoin,ielpo,ielem,jelem,ielel,jnode
  integer(ip)           :: ptopo,pnode,kelem,ihang,jpoin,idime
  integer(ip)           :: nsize,isize,knode,kpoin,pelty
  integer(ip)           :: llist(20)
  integer(ip)           :: lelis(mnode)
  real(rp)              :: elcod(ndime,mnode)
  character(20)         :: winou,wwher
  integer(ip),  pointer :: pelem(:)     => null()
  type(i1p),    pointer :: lhang_tmp(:) => null()
  logical(lg),  pointer :: touch(:)     => null()

  if( INOTMASTER .and. kfl_hangi == 1 ) then

     call memory_alloca(memor_dom,'PELEM',    'hangin',pelem,    nelem)
     call memory_alloca(memor_dom,'LHANG_TMP','hangin',lhang_tmp,npoin)
     call memory_alloca(memor_dom,'TOUCH',    'hangin',touch,    nelem,'DO_NOT_INITIALIZE')
     do ielem = 1,nelem
        touch(ielem) = .true.
     end do
     nhang = 0

     do ipoin = 1,npoin
        if( lpoty(ipoin) /= 0 ) then
           kelem = 0
           do ielpo = pelpo_2(ipoin),pelpo_2(ipoin+1)-1
              ielem = lelpo_2(ielpo)
              do ielel = pelel_2(ielem),pelel_2(ielem+1)-1
                 jelem = lelel_2(ielel)
                 if( touch(jelem) ) then
                    touch(jelem) = .false.
                    kelem        = kelem + 1
                    pelem(kelem) = jelem
                    pelty        = abs(ltype(jelem))
                    ptopo        = ltopo(pelty)
                    pnode        = lnnod(jelem)
                    do jnode = 1,pnode
                       jpoin = lnods(jnode,jelem)
                       do idime = 1,ndime
                          elcod(idime,jnode) = coord(idime,jpoin)
                       end do
                    end do
                    call elmgeo_where_is(&
                         ndime,pnode,ptopo,pelty,elcod(1:ndime,1:pnode),&
                         coord(1:ndime,ipoin),winou,wwher,llist)                    
                    if( trim(winou) == 'INSIDE' .and. ipoin == 6 ) then
                       if( trim(wwher) == 'VERTEX' ) then
                          !print*,jelem,lnods(llist(1),jelem)
                       else if( trim(wwher) == 'MID_EDGE' ) then
                          nhang = nhang + 1
                          call memory_alloca(memor_dom,'LHANG_TMP','hangin',lhang_tmp(nhang) % l,3_ip)
                          lhang_tmp(nhang) % l(1)   = ipoin
                          lhang_tmp(nhang) % l(2:3) = lnods(llist(1:2),jelem)                    
                          print*,'MID_EDGE=',jelem,lnods(llist(1:2),jelem)                            
                       else if( trim(wwher) == 'MID_FACE' ) then
                          call runend('HAGING: NOT CODED')
                       end if
                    end if
                 end if
              end do
           end do
           do ielem = 1,kelem
              touch(pelem(ielem)) = .true.
           end do
        end if
     end do
     !
     ! Copy temporary LHANG_TMP to permanent LHANG
     !
     call memory_alloca(memor_dom,'LHANG','hangin',lhang,nhang)
     do ihang = 1,nhang
        nsize = size(lhang_tmp(ihang) % l,KIND=ip)
        call memory_alloca(memor_dom,'LHANG','hangin',lhang(ihang) % l,nsize,'DO_NOT_INITIALIZE')
        do isize = 1,nsize
           lhang(ihang) % l(isize) = lhang_tmp(ihang) % l(isize)
        end do
     end do
     !
     ! Elements can have various hanging nodes o. 
     ! Example with (e):
     ! +-------+
     ! |       |
     ! |  h1   |
     ! +---o---1------+
     ! |   | e |      |
     ! +---2---o h2   |
     ! |   |   |      |
     ! +---+---+------+
     !
     ! LELCH(IELEM) = ELHAN if element contains a hanging node
     ! NEHAN(IELEM) = Number of hanging nodes in the element
     ! LEHAN(IELEM) % l(1) % l(1:2) = h1,1 
     ! LEHAN(IELEM) % l(2) % l(1:2) = h2,1 
     ! In the example, h1 and h2 are 2 hanging nodes (NEHAN(IELEM)=2) linked to:
     ! For h1: node 1
     ! For h2: node 1
     !
     call memory_alloca(memor_dom,'NEHAN','hangin',nehan,nelem,'DO_NOT_INITIALIZE')
     do ihang = 1,nhang
        ipoin = lhang_tmp(ihang) % l(1)
        do ielpo = pelpo_2(ipoin),pelpo_2(ipoin+1)-1
           ielem = lelpo_2(ielpo)
           lelch(ielem) = ELHAN
           if( ielem <= nelem ) nehan(ielem) = nehan(ielem) + 1
        end do
     end do
     !
     ! Fill in hanging element list
     !
     call memory_alloca(memor_dom,'LEHAN','hangin',lehan,nelem,'DO_NOT_INITIALIZE')
     do ielem = 1,nelem
        if( lelch(ielem) == ELHAN ) then
           nsize = nehan(ielem)
           call memory_alloca(memor_dom,'LEHAN','hangin',lehan(ielem) % l,nsize,'DO_NOT_INITIALIZE')
           nehan(ielem) = 0
        end if
     end do
     do ihang = 1,nhang
        ipoin = lhang_tmp(ihang) % l(1)
        do ielpo = pelpo_2(ipoin),pelpo_2(ipoin+1)-1
           ielem = lelpo_2(ielpo)
           if( ielem <= nelem ) then
              nehan(ielem) = nehan(ielem) + 1
              nsize = 0 
              do knode = 1,size(lhang(ihang) % l,KIND=ip)
                 kpoin = lhang(ihang) % l(knode)
                 do jnode = 1,pnode
                    jpoin = lnods(jnode,ielem)
                    if( jpoin == kpoin ) then
                       nsize = nsize + 1
                       lelis(nsize) = jnode
                    end if
                 end do
              end do
              call memory_alloca(memor_dom,'LEHAN','hangin',lehan(ielem) % l(nehan(ielem)) % l,nsize,'DO_NOT_INITIALIZE')
              do isize = 1,nsize
                 lehan(ielem) % l(nehan(ielem)) % l(isize) = lelis(isize)
              end do
           end if
        end do
     end do

     call memory_deallo(memor_dom,'TOUCH',    'hangin',touch)
     call memory_deallo(memor_dom,'LHANG_TMP','hangin',lhang_tmp)
     call memory_deallo(memor_dom,'PELEM',    'hangin',pelem)

  end if

end subroutine hangin
