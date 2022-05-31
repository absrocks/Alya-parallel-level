subroutine elsest_octsub(&
     ithre,ndime,mnode,npoin,nelem,lnods,ltype,coord,nnode)
  !------------------------------------------------------------------------
  !****f* elsest/elsest_octsub
  ! NAME 
  !    elsest_octsub
  ! DESCRIPTION
  !    Create the quad/oct-tree structure
  ! USES
  ! USED BY
  !    elsest_octpre
  !***
  !------------------------------------------------------------------------
  use def_elsest
  use mod_elsest
  implicit none
  integer(ip), intent(in) :: ithre
  integer(ip), intent(in) :: mnode,ndime,npoin,nelem
  integer(ip), intent(in) :: nnode(*)
  integer(ip), intent(in) :: lnods(mnode,*),ltype(*)
  real(rp),    intent(in) :: coord(ndime,npoin)
  integer(ip)             :: i,ipoin,kpoin,ielem,istat,counter
  real(rp)                :: time1,time2
  logical(lg)             :: conti
  type(octbox), pointer   :: old_pointer,tm1_pointer,tm2_pointer
  real(rp),     pointer   :: xmima(:,:,:)

  !----------------------------------------------------------------------
  !
  ! Compute element bounding boxes
  !
  !----------------------------------------------------------------------

  allocate( xmima(3,2,nelem), stat = istat )
  call elsest_elbbox(&
       mnode,ndime,npoin,nelem,nnode,lnods,ltype,coord,xmima)

  !----------------------------------------------------------------------
  !
  ! Allocate memory
  !
  !----------------------------------------------------------------------

  allocate(lboel(nelem),stat=istat)
  call elsest_memchk(0_ip,ithre,istat,memor(5,ithre),'LBOEL','elsest_octsub',lboel)

  divmax  = 2**ndime
  conti   = .true.
  counter = 0

  do while( conti )
     !
     ! If maximum number of points inside current box is exceeded, subdivide
     !     
     if( current(ithre) % o % npoinbox > limit ) then
        call elsest_cputim(time1)
        allocate(current(ithre) % o % children(divmax),stat=istat)
        if(istat/=0) call elsest_memerr(0_ip,'CURRENT_CHILDREN','elsest_octsub',0_ip)
        current(ithre) % o % children(1:divmax) = octbox_init
        !
        ! Give birth to my DIVMAX children
        !
        do i = 1,divmax
           counter                                     =  counter+1
           current(ithre) % o % children(i) % id       =  counter
           current(ithre) % o % children(i) % childid  =  i
           current(ithre) % o % children(i) % level    =  current(ithre) % o % level + 1 
           current(ithre) % o % children(i) % whoiam   =  0  
           current(ithre) % o % children(i) % npoinbox =  0
           current(ithre) % o % children(i) % nelembox =  0
           current(ithre) % o % children(i) % parent   => current(ithre) % o
 
           allocate(current(ithre) % o % children(i) % nodes(current(ithre) % o % npoinbox),stat=istat)
           call elsest_memchk(0_ip,ithre,istat,memor(2,ithre),'current % children(i) % nodes',&
                'elsest_octsub',current(ithre) % o % children(i) % nodes)
        end do
        !
        ! Compute the coordinates of my children
        !
        do i = 1,ndime
           current(ithre) % o % children(1)%minc(i) = current(ithre) % o % minc(i)
           current(ithre) % o % children(1)%maxc(i) = (current(ithre) % o % maxc(i) + current(ithre) % o % minc(i))*0.5_rp
        end do
        current(ithre) % o % children(2)%minc(1) = (current(ithre) % o % maxc(1) + current(ithre) % o % minc(1))*0.5_rp
        current(ithre) % o % children(2)%minc(2) = current(ithre) % o % children(1)%minc(2)
        current(ithre) % o % children(2)%maxc(1) = current(ithre) % o % maxc(1)     
        current(ithre) % o % children(2)%maxc(2) = current(ithre) % o % children(1)%maxc(2)

        current(ithre) % o % children(3)%minc(1) = current(ithre) % o % children(1)%minc(1)
        current(ithre) % o % children(3)%minc(2) = (current(ithre) % o % minc(2) + current(ithre) % o % maxc(2))*0.5_rp
        current(ithre) % o % children(3)%maxc(1) = current(ithre) % o % children(1)%maxc(1)
        current(ithre) % o % children(3)%maxc(2) = current(ithre) % o % maxc(2)

        current(ithre) % o % children(4)%minc(1) = current(ithre) % o % children(2)%minc(1)
        current(ithre) % o % children(4)%minc(2) = current(ithre) % o % children(3)%minc(2)
        current(ithre) % o % children(4)%maxc(1) = current(ithre) % o % children(2)%maxc(1)
        current(ithre) % o % children(4)%maxc(2) = current(ithre) % o % children(3)%maxc(2)

        if( ndime == 3 ) then
           current(ithre) % o % children(2)%minc(3) = current(ithre) % o % children(1)%minc(3)
           current(ithre) % o % children(2)%maxc(3) = current(ithre) % o % children(1)%maxc(3)
           current(ithre) % o % children(3)%minc(3) = current(ithre) % o % children(1)%minc(3)
           current(ithre) % o % children(3)%maxc(3) = current(ithre) % o % children(1)%maxc(3)
           current(ithre) % o % children(4)%minc(3) = current(ithre) % o % children(1)%minc(3)
           current(ithre) % o % children(4)%maxc(3) = current(ithre) % o % children(1)%maxc(3)

           current(ithre) % o % children(5)%minc(1) = current(ithre) % o % children(1)%minc(1)
           current(ithre) % o % children(5)%minc(2) = current(ithre) % o % children(1)%minc(2)
           current(ithre) % o % children(5)%minc(3) = (current(ithre) % o % minc(3) + current(ithre) % o % maxc(3))*0.5_rp
           current(ithre) % o % children(5)%maxc(1) = current(ithre) % o % children(1)%maxc(1)
           current(ithre) % o % children(5)%maxc(2) = current(ithre) % o % children(1)%maxc(2)
           current(ithre) % o % children(5)%maxc(3) = current(ithre) % o % maxc(3)

           current(ithre) % o % children(6)%minc(1) = current(ithre) % o % children(2)%minc(1)
           current(ithre) % o % children(6)%minc(2) = current(ithre) % o % children(1)%minc(2)
           current(ithre) % o % children(6)%minc(3) = current(ithre) % o % children(5)%minc(3)
           current(ithre) % o % children(6)%maxc(1) = current(ithre) % o % children(2)%maxc(1)     
           current(ithre) % o % children(6)%maxc(2) = current(ithre) % o % children(1)%maxc(2)
           current(ithre) % o % children(6)%maxc(3) = current(ithre) % o % children(5)%maxc(3)

           current(ithre) % o % children(7)%minc(1) = current(ithre) % o % children(1)%minc(1)
           current(ithre) % o % children(7)%minc(2) = current(ithre) % o % children(3)%minc(2)
           current(ithre) % o % children(7)%minc(3) = current(ithre) % o % children(5)%minc(3)
           current(ithre) % o % children(7)%maxc(1) = current(ithre) % o % children(1)%maxc(1)
           current(ithre) % o % children(7)%maxc(2) = current(ithre) % o % children(3)%maxc(2)
           current(ithre) % o % children(7)%maxc(3) = current(ithre) % o % children(5)%maxc(3)

           current(ithre) % o % children(8)%minc(1) = current(ithre) % o % children(2)%minc(1)
           current(ithre) % o % children(8)%minc(2) = current(ithre) % o % children(3)%minc(2)
           current(ithre) % o % children(8)%minc(3) = current(ithre) % o % children(5)%minc(3)
           current(ithre) % o % children(8)%maxc(1) = current(ithre) % o % children(2)%maxc(1)
           current(ithre) % o % children(8)%maxc(2) = current(ithre) % o % children(3)%maxc(2)
           current(ithre) % o % children(8)%maxc(3) = current(ithre) % o % children(5)%maxc(3)
        end if
        !
        ! Offer my nodes to my children
        !
        if( ndime == 2 ) then
           do i=1,4
              do ipoin=1,current(ithre) % o % npoinbox
                 kpoin=current(ithre) % o % nodes(ipoin)
                 if(    coord(1,kpoin) >= current(ithre) % o % children(i)%minc(1) .and. &
                      & coord(2,kpoin) >= current(ithre) % o % children(i)%minc(2) .and. &
                      & coord(1,kpoin) <  current(ithre) % o % children(i)%maxc(1) .and. &
                      & coord(2,kpoin) <  current(ithre) % o % children(i)%maxc(2) ) then
                    current(ithre) % o % children(i)%npoinbox=current(ithre) % o % children(i)%npoinbox+1
                    current(ithre) % o % children(i)%nodes(current(ithre) % o % children(i)%npoinbox)=kpoin
                 end if
              end do
           end do
        else
           do i=1,8
              do ipoin=1,current(ithre) % o % npoinbox
                 kpoin=current(ithre) % o % nodes(ipoin)
                 if(    coord(1,kpoin) >= current(ithre) % o % children(i)%minc(1) .and. &
                      & coord(2,kpoin) >= current(ithre) % o % children(i)%minc(2) .and. &
                      & coord(3,kpoin) >= current(ithre) % o % children(i)%minc(3) .and. &
                      & coord(1,kpoin) <  current(ithre) % o % children(i)%maxc(1) .and. &
                      & coord(2,kpoin) <  current(ithre) % o % children(i)%maxc(2) .and. &
                      & coord(3,kpoin) <  current(ithre) % o % children(i)%maxc(3) ) then
                    current(ithre) % o % children(i)%npoinbox=current(ithre) % o % children(i)%npoinbox+1
                    current(ithre) % o % children(i)%nodes(current(ithre) % o % children(i)%npoinbox)=kpoin
                 end if
              end do
           end do
        end if
        call elsest_memchk(2_ip,ithre,istat,memor(2,ithre),'CURRENT%NODES','elsest_octsub',current(ithre) % o % nodes)
        deallocate(current(ithre) % o % nodes,stat=istat)
        if(istat/=0) call elsest_memerr(2_ip,'CURRENT%NODES','elsest_octsub',0_ip)
        current(ithre) % o % whoiam   =  0
        current(ithre) % o % npoinbox =  0
        current(ithre) % o            => current(ithre) % o % children(1)
        call elsest_cputim(time2)
        cputi(2,ithre)=cputi(2,ithre)+time2-time1

     else if(current(ithre) % o % id == 0 .and. current(ithre) % o % npoinbox <= limit ) then
        !
        ! If the Padrino has too few elements
        !
        allocate(current(ithre) % o % elems(nelem),stat=istat)
        call elsest_memchk(0_ip,ithre,istat,memor(2,ithre),'TREE_ROOT%ELEMS','elsest_octsub',tree_root%elems)

        do ielem = 1,nelem
           current(ithre) % o % elems(ielem) = ielem
        end do
        conti = .false.
        current(ithre) % o % nelembox =  nelem
        current(ithre) % o % whoiam   =  1
        current(ithre) % o            => old_pointer

     else 
        !
        ! if limit of points inside box is not exceeded, assign elements
        !
        call elsest_cputim(time1)
        kstat(6,ithre) = kstat(6,ithre) + 1
        if(current(ithre) % o % npoinbox < kstat(1,ithre)) kstat(1,ithre) = current(ithre) % o % npoinbox
        if(current(ithre) % o % npoinbox > kstat(2,ithre)) kstat(2,ithre) = current(ithre) % o % npoinbox 
        kstat(8,ithre) = kstat(8,ithre) + current(ithre) % o % npoinbox

        if( ndime == 1 ) then
           do ielem = 1,nelem
              if( xmima(1,1,ielem) <= current(ithre) % o % maxc(1) ) then
                 if( xmima(1,2,ielem) >= current(ithre) % o % minc(1) ) then
                    current(ithre) % o % nelembox = current(ithre) % o % nelembox + 1
                    lboel(current(ithre) % o % nelembox ) = ielem
                 end if
              end if
           end do
        else if( ndime == 2 ) then
           do ielem = 1,nelem
              if( xmima(1,1,ielem) <= current(ithre) % o % maxc(1) ) then
                 if( xmima(1,2,ielem) >= current(ithre) % o % minc(1) ) then
                    if( xmima(2,1,ielem) <= current(ithre) % o % maxc(2) ) then
                       if( xmima(2,2,ielem) >= current(ithre) % o % minc(2) ) then
                          current(ithre) % o % nelembox = current(ithre) % o % nelembox + 1
                          lboel(current(ithre) % o % nelembox ) = ielem
                       end if
                    end if
                 end if
              end if
           end do
        else if( ndime == 3 ) then
           do ielem = 1,nelem
              if( xmima(1,1,ielem) <= current(ithre) % o % maxc(1) ) then
                 if( xmima(1,2,ielem) >= current(ithre) % o % minc(1) ) then
                    if( xmima(2,1,ielem) <= current(ithre) % o % maxc(2) ) then
                       if( xmima(2,2,ielem) >= current(ithre) % o % minc(2) ) then
                          if( xmima(3,1,ielem) <= current(ithre) % o % maxc(3) ) then
                             if( xmima(3,2,ielem) >= current(ithre) % o % minc(3) ) then
                                current(ithre) % o % nelembox = current(ithre) % o % nelembox + 1
                                lboel(current(ithre) % o % nelembox ) = ielem
                             end if
                          end if
                       end if
                    end if
                 end if
              end if
           end do
        end if
        !
        ! Put elements connected to the nodes
        !
        !do ipoin=1,current(ithre) % o % npoinbox
        !   do ielpo=pelpo_els(current(ithre) % o % nodes(ipoin)),pelpo_els(current(ithre) % o % nodes(ipoin)+1)-1
        !      element_already_added =.false.
        !      if(current(ithre) % o % nelembox>0) then
        !         loop_over_elements_of_lboel: do kk=1,current(ithre) % o % nelembox
        !            if(lelpo_els(ielpo)==lboel(kk)) then 
        !               element_already_added =.true.
        !               exit loop_over_elements_of_lboel
        !            end if
        !         end do loop_over_elements_of_lboel
        !      end if
        !      if(.not.element_already_added) then
        !         current(ithre) % o % nelembox=current(ithre) % o % nelembox+1
        !         lboel(current(ithre) % o % nelembox)=lelpo_els(ielpo)
        !      end if
        !   end do
        !end do
        !
        ! Look for elements crossing the Quad
        !
        if(current(ithre) % o % nelembox < kstat(3,ithre)) kstat(3,ithre) = current(ithre) % o % nelembox
        if(current(ithre) % o % nelembox > kstat(4,ithre)) kstat(4,ithre) = current(ithre) % o % nelembox
        kstat(7,ithre) = kstat(7,ithre)+current(ithre) % o % nelembox

        call elsest_memchk(2_ip,ithre,istat,memor(2,ithre),'CURRENT(ITHRE) % O % NODES','elsest_octsub',current(ithre) % o % nodes)
        deallocate(current(ithre) % o % nodes,stat=istat)
        if(istat/=0) call elsest_memerr(2_ip,'CURRENT%NODES','elsest_octsub',0_ip)     
        !
        ! Here we assign elements
        !
        if( current(ithre) % o % nelembox /= 0 ) then
           current(ithre) % o % whoiam  = 1
           allocate( current(ithre) % o % elems(current(ithre) % o % nelembox),stat=istat)
           call elsest_memchk(0_ip,ithre,istat,memor(2,ithre),'CURRENT%ELEMS','elsest_octsub',current(ithre) % o % elems)
           do ielem = 1,current(ithre) % o % nelembox
              current(ithre) % o %elems(ielem) = lboel(ielem)
           end do
        else
           current(ithre) % o % whoiam  = 2
        end if
        call elsest_cputim(time2)
        cputi(3,ithre)=cputi(3,ithre)+time2-time1

        if( current(ithre) % o % childid < divmax .and. current(ithre) % o % id /= 0 ) then
           !
           ! Go to next children
           !
           tm1_pointer      => current(ithre)%o
           tm2_pointer      => tm1_pointer%parent%children(tm1_pointer%childid+1)
           current(ithre)%o => tm2_pointer
           goto 10

        else if(current(ithre) % o % childid == 0 ) then  
           !
           ! Padrino
           !
           goto 10

        else if(current(ithre) % o % childid == divmax ) then
           !
           ! Last children
           !
           noparent: do while( current(ithre) % o % id > 0 )
              if(current(ithre) % o % parent%id == 0) then
                 conti=.false.
                 exit noparent
              else
                 if(current(ithre) % o % parent%childid /=divmax) then 
                    tm1_pointer       => current(ithre)%o
                    tm2_pointer       => tm1_pointer%parent%parent%children(tm1_pointer%parent%childid+1)
                    current(ithre)%o  => tm2_pointer
                    exit
                 else 
                    current(ithre)%o => current(ithre) % o % parent
                 end if
              end if
           end do noparent

        else 
           !
           ! Wrong child ID
           !
           call elsest_runend('WRONG CHILD ID: '//trim(elsest_intost(current(ithre) % o % childid)))    
        end if

     end if

10   continue
     old_pointer => current(ithre) % o

  end do
  !
  ! Deallocate memory
  !
  call elsest_memchk(2_ip,ithre,istat,memor(5,ithre),'LBOEL','elsest_octsub',lboel)
  deallocate(lboel,stat=istat)
  if(istat/=0) call elsest_memerr(2_ip,'LBOEL','elsest_octsub',0_ip)
  deallocate( xmima, stat = istat )

end subroutine elsest_octsub

