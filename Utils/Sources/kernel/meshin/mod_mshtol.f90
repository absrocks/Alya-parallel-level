module mod_mshtol
  !
  !     Use openmp  $OMP1
  ! 

contains

  subroutine ptoelm( lface, nface , npoin, nnofa,ptoel1,ptoel2 )
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only          : memor_msh
    implicit none
    integer(ip), intent(in)            :: nface,npoin,nnofa
    integer(ip), intent(in)            :: lface(nnofa,nface)
    integer(ip), pointer :: ptoel1(:),ptoel2(:)
    integer(ip)                        :: iface,inode,iplace,ipoin,nptoel
    integer(4)                 :: istat
    !
    !     This subroutine gets the elements surrounding points
    !
    !     On input: the connectivity lface 
    !        
    !     On output: the pointer to elements ptoel2,ptoel1
    ! 
    if(.not.associated(ptoel2))then 
       allocate(ptoel2(npoin+1),stat=istat)
       call memchk(zero,istat,memor_msh,'PTOEL2','ptoelm',ptoel2)
    else
       call memrea(npoin+1_ip,memor_msh,'PTOEL2','ptoelm',ptoel2)
       do ipoin=1,npoin+1
          ptoel2(ipoin)=0_ip
       enddo
    endif
    !
    !    Store ahead
    !
    do iface =1,nface
       do inode=1,nnofa
          ipoin=lface(inode,iface)+1
          ptoel2(ipoin)=ptoel2(ipoin)+1
       enddo
    enddo
    !
    !     Sum up
    !
    ptoel2(1)=1
    do ipoin=2,npoin+1 
       ptoel2(ipoin)=ptoel2(ipoin)+ptoel2(ipoin-1)
    enddo
    !
    !      Allocate ptoel1
    !
    nptoel=ptoel2(npoin+1)-1

    if(.not.associated(ptoel1))then
       allocate(ptoel1(nptoel),stat=istat)
       call memchk(zero,istat,memor_msh,'PTOEL1','ptoelm',ptoel1)
    else
       call memrea(nptoel,memor_msh,'PTOEL1','ptoelm',ptoel1)
    endif
    !
    !     Store in ptoel1
    !
    do iface =1,nface
       do inode=1,nnofa
          ipoin=lface(inode,iface)
          iplace=ptoel2(ipoin)
          ptoel1(iplace)=iface
          ptoel2(ipoin)=iplace+1
       enddo
    enddo
    !
    !     Clean up 
    !  
    do ipoin=npoin+1,2,-1
       ptoel2(ipoin)=ptoel2(ipoin-1)
    enddo

    ptoel2(1)=1

  end subroutine ptoelm

  subroutine ptoelmb( lface, nface , npoin, nnofa,ptoel1,ptoel2 )
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only          : memor_msh
    implicit none
    integer(ip), intent(in)            :: nface,npoin,nnofa
    integer(ip), intent(in)            :: lface(nnofa,nface)
    integer(ip), pointer :: ptoel1(:),ptoel2(:)
    integer(ip)                        :: iface,inode,iplace,ipoin,nptoel
    integer(4)                 :: istat
    !
    !     This subroutine gets the elements surrounding points
    !     This is a special version if lface contains deleted elements
    !     marked with 0 in the first position
    !
    !     On input: the connectivity lface 
    !        
    !     On output: the pointer to elements ptoel2,ptoel1
    !

    if(.not.associated(ptoel2))then 
       allocate(ptoel2(npoin+1),stat=istat)
       call memchk(zero,istat,memor_msh,'PTOEL2','ptoelm',ptoel2)
    else
       call memrea(npoin+1_ip,memor_msh,'PTOEL2','ptoelm',ptoel2)
       do ipoin=1,npoin+1
          ptoel2(ipoin)=0_ip
       enddo
    endif
    !
    !    Store ahead
    !
    do iface =1,nface
       if(lface(1,iface)==0)cycle
       do inode=1,nnofa
          ipoin=lface(inode,iface)+1
          ptoel2(ipoin)=ptoel2(ipoin)+1
       enddo
    enddo
    !
    !     Sum up
    !
    ptoel2(1)=1
    do ipoin=2,npoin+1 
       ptoel2(ipoin)=ptoel2(ipoin)+ptoel2(ipoin-1)
    enddo
    !
    !      Allocate ptoel1
    !
    nptoel=ptoel2(npoin+1)-1

    if(.not.associated(ptoel1))then
       allocate(ptoel1(nptoel),stat=istat)
       call memchk(zero,istat,memor_msh,'PTOEL1','ptoelm',ptoel1)
    else
       call memrea(nptoel,memor_msh,'PTOEL1','ptoelm',ptoel1)
    endif
    !
    !     Store in ptoel1
    !
    do iface =1,nface
       if(lface(1,iface)==0)cycle
       do inode=1,nnofa
          ipoin=lface(inode,iface)
          iplace=ptoel2(ipoin)
          ptoel1(iplace)=iface
          ptoel2(ipoin)=iplace+1
       enddo
    enddo
    !
    !     Clean up 
    !  
    do ipoin=npoin+1,2,-1
       ptoel2(ipoin)=ptoel2(ipoin-1)
    enddo

    ptoel2(1)=1

  end subroutine ptoelmb

  subroutine ptopoi( elem, nelem , npoin, nnode,ptoel1,ptoel2,ptopo1,ptopo2 )
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only          : memor_msh
    implicit none
    integer(ip), intent(in) :: nelem,npoin,nnode
    integer(ip), intent(in) :: elem(nnode,nelem)
    integer(ip)             :: ptoel1(*),ptoel2(*)
    integer(ip), pointer    :: ptopo1(:),ptopo2(:)
    integer(ip), pointer    :: lpoin(:)
    integer(ip)             :: inode,ipoin,jpoin
    integer(ip)             :: ncont,ie,ielem
    integer(4)              :: istat

    if(.not.associated(ptopo2))then 
       allocate(ptopo2(npoin+1),stat=istat)
       call memchk(zero,istat,memor_msh,'PTOPO2','ptopoi',ptopo2)
    else
       call memrea(npoin+1_ip,memor_msh,'PTOPO2','ptopoi',ptopo2)
       do ipoin=1,npoin+1
          ptopo2(ipoin)=0_ip
       enddo
    endif
    !
    !     Allocate lpoin       
    !
    allocate(lpoin(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LPOIN','ptopoi',lpoin)
    !
    !     Count the neighbors
    !
    ncont=0_ip   
    do ipoin=1,npoin
       lpoin(ipoin)=ipoin
       do ie=ptoel2(ipoin),ptoel2(ipoin+1)-1
          ielem=ptoel1(ie)
          do inode=1,nnode 
             jpoin=elem(inode,ielem)
             if(lpoin(jpoin)/=ipoin)then
                lpoin(jpoin)=ipoin  
                ncont=ncont+1
             endif
          enddo
       enddo
    enddo
    !
    !      Allocate ptoel1
    !
    if(.not.associated(ptopo1))then
       allocate(ptopo1(ncont),stat=istat)
       call memchk(zero,istat,memor_msh,'PTOPO1','ptopoi',ptopo1)
    else
       call memrea(ncont,memor_msh,'PTOPO1','ptopoi',ptopo1)
    endif
    !
    !     Clean up lpoin
    !
    do ipoin=1,npoin
       lpoin(ipoin)=0_ip
    enddo
    !
    !     Store in ptop1
    ! 
    ncont=0_ip  
    ptopo2(1)=1_ip 
    do ipoin=1,npoin
       lpoin(ipoin)=ipoin
       do ie=ptoel2(ipoin),ptoel2(ipoin+1)-1
          ielem=ptoel1(ie)
          do inode=1,nnode 
             jpoin=elem(inode,ielem)
             if(lpoin(jpoin)/=ipoin)then
                lpoin(jpoin)=ipoin  
                ncont=ncont+1
                ptopo1(ncont)=jpoin 
             endif
          enddo
       enddo
       ptopo2(ipoin+1)=ncont+1_ip
    enddo

    call memchk(2_ip,istat,memor_msh,'LPOIN','ptopoi',lpoin)
    deallocate(lpoin,stat=istat)
    if(istat/=0) call memerr(2_ip,'LPOIN','ptopoi',0_ip)

  end subroutine ptopoi

  subroutine eltolm(lface,nface,ptoel1,ptoel2,eltoel,nnode,npoin)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       :memor_msh
    implicit none
    integer(ip), intent(in)    :: nface,nnode,npoin
    integer(ip), intent(in)    :: lface(nnode,nface)
    integer(ip), intent(in)    :: ptoel1(*),ptoel2(npoin+1)
    integer(ip), pointer       :: eltoel(:,:)
    integer(4)                 :: istat
    integer(ip)                :: ltab(2,3)
    logical(lg), pointer       :: lmark(:)
    integer(ip)                :: iface,ipoin,ipmin,inode,inodd,jptoel,nptoel,nptoelt
    integer(ip)                :: jnode,jnodd,icount,jpoin,jface  

    ltab(1,1)=2
    ltab(2,1)=3
    ltab(1,2)=3
    ltab(2,2)=1
    ltab(1,3)=1
    ltab(2,3)=2

    !
    !     Allocate  eltoel
    !   
    allocate(eltoel(nnode,nface),stat=istat)
    call memchk(zero,istat,memor_msh,'ELTOEL','eltolm',eltoel)
    !
    !     Allocate help array
    !
    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','eltolm',lmark)

    do iface=1,nface
       do inode=1,nnode 

          if(eltoel(inode,iface)==0_ip)then 
             ! 
             !     Mark lmark 
             !        
             ipoin=lface(ltab(1,inode),iface)
             lmark(ipoin)=.true.
             nptoel=ptoel2(ipoin+1)-ptoel2(ipoin)

             do inodd=2,nnode
                ipoin=lface(ltab(inodd,inode),iface)
                lmark(ipoin)=.true.
                !
                !     Take the point with the least number of surrounding elements
                !
                nptoelt=ptoel2(ipoin+1)-ptoel2(ipoin)

                if(nptoelt<nptoel)then 
                   ipmin=ipoin
                   nptoel=nptoelt
                endif
             enddo
             !
             !     Loop on the elements surrounding ipoin
             !
             do jptoel=ptoel2(ipmin),ptoel2(ipmin+1)-1
                jface=ptoel1(jptoel)
                if(jface/=iface)then
                   do jnode=1,nnode

                      icount=0
                      do jnodd=1,nnode
                         jpoin=lface(ltab(jnodd,jnode),jface)
                         if(lmark(jpoin))then
                            icount=icount+1
                         endif
                      enddo

                      if(icount==nnode)then

                         eltoel(inode,iface)=jface
                         eltoel(jnode,jface)=iface
                         goto 10

                      endif
                   enddo
                endif
             enddo
          endif
10        continue 
       enddo
    enddo

    call memchk(2_ip,istat,memor_msh,'LMARK','eltolm',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','eltolm',0_ip)

  end subroutine eltolm

  subroutine trtotr(lface,nnofa,nface,ptoel1,ptoel2,npoin,eltoel)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip),intent(in)     :: nnofa,nface,npoin 
    integer(ip),intent(in)     :: lface(nnofa,nface)
    integer(ip),intent(in)     :: ptoel1(*),ptoel2(npoin+1) 
    integer(ip), pointer       :: eltoel(:,:) 
    integer(4)                 :: istat
    integer(ip)                :: ltab(2,3),iface,jface,jp1,jp2,jnode,inode,ip1,ip2,ipmin,inofa
    integer(ip)                :: nptoelt1,nptoelt2,jptoel,jp3
    integer(ip),pointer        :: lmark(:)

    ltab(1,1)=2
    ltab(2,1)=3
    ltab(1,2)=3
    ltab(2,2)=1
    ltab(1,3)=1
    ltab(2,3)=2

    !
    !     Allocate  eltoel
    !   
    if(.not.associated(eltoel))then 
       allocate(eltoel(nnofa,nface),stat=istat)
       call memchk(zero,istat,memor_msh,'ELTOEL','trtotr',eltoel)
    else
       call memrea(nface,memor_msh,'ELTOEL','trtotr',eltoel)
       do iface=1,nface
          do inofa=1,nnofa 
             eltoel(inofa,iface)=0_ip
          enddo
       enddo 
    endif
    !
    !     Allocate help array
    !
    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','seltoel',lmark)

    do iface=1,nface
       do inode=1,3 
          if(eltoel(inode,iface)==0_ip)then 
             ! 
             !     Mark lmark 
             !        
             ip1=lface(ltab(1,inode),iface)
             lmark(ip1)=1_ip
             ip2=lface(ltab(2,inode),iface)
             lmark(ip2)=1_ip
             !
             !     Take the point with the least number of surrounding elements
             !
             nptoelt1=ptoel2(ip1+1)-ptoel2(ip1)
             nptoelt2=ptoel2(ip2+1)-ptoel2(ip2)

             if(nptoelt1<nptoelt2)then 
                ipmin=ip1
             else
                ipmin=ip2
             endif
             !
             !     Loop on the elements surrounding ipoin
             !
             do jptoel=ptoel2(ipmin),ptoel2(ipmin+1)-1
                jface=ptoel1(jptoel)
                if(jface/=iface)then
                   jp1=lface(1,jface)
                   jp2=lface(2,jface)
                   jp3=lface(3,jface)

                   if(lmark(jp1)+lmark(jp2)+lmark(jp3)==2_ip)then                   
                      do jnode=1,3
                         jp1=lface(ltab(1,jnode),jface)
                         jp2=lface(ltab(2,jnode),jface)

                         if(lmark(jp1)==1 .and. lmark(jp2)==1) then

                            eltoel(inode,iface)=jface
                            eltoel(jnode,jface)=iface
                            goto 100

                         endif
                      enddo
                   endif
                endif
             enddo
100          continue 
             !
             !     Clean up lmark
             !       
             lmark(ip1)=0
             lmark(ip2)=0 

          endif
       enddo
    enddo

    call memchk(2_ip,istat,memor_msh,'LMARK','trtotr',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','trtotr',0_ip)

  end subroutine trtotr

  subroutine trtotr2(lface,nnofa,nface,ptoel1,ptoel2,npoin,eltoel,rnofa,ndim)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip),intent(in)     :: nnofa,nface,npoin,ndim 
    integer(ip),intent(in)     :: lface(nnofa,nface)
    real(rp),intent(in)        :: rnofa(ndim,nface)
    integer(ip),intent(in)     :: ptoel1(*),ptoel2(npoin+1) 
    integer(ip), pointer       :: eltoel(:,:) 
    integer(4)                 :: istat
    integer(ip)                :: ltab(2,3),iface,jface,jp1,jp2,jnode,inode,ip1,ip2,ipmin,inofa
    integer(ip)                :: nptoelt1,nptoelt2,jptoel,jp3
    integer(ip),pointer        :: lmark(:)
    real(rp)                   :: csca,c00 
    !
    !     Compared to trtotr, this subroutine takes into account the normal of the faces
    !
    ltab(1,1)=2
    ltab(2,1)=3
    ltab(1,2)=3
    ltab(2,2)=1
    ltab(1,3)=1
    ltab(2,3)=2
    c00=0.0d+00
    !
    !     Allocate  eltoel
    !   
    if(.not.associated(eltoel))then 
       allocate(eltoel(nnofa,nface),stat=istat)
       call memchk(zero,istat,memor_msh,'ELTOEL','trtotr',eltoel)
    else
       call memrea(nface,memor_msh,'ELTOEL','trtotr',eltoel)
       do iface=1,nface
          do inofa=1,nnofa 
             eltoel(inofa,iface)=0_ip
          enddo 
       enddo 
    endif
    !
    !     Allocate help array
    !
    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','seltoel',lmark)

    do iface=1,nface
       do inode=1,3 
          if(eltoel(inode,iface)==0_ip)then 
             ! 
             !     Mark lmark 
             !        
             ip1=lface(ltab(1,inode),iface)
             lmark(ip1)=1_ip
             ip2=lface(ltab(2,inode),iface)
             lmark(ip2)=1_ip
             !
             !     Take the point with the least number of surrounding elements
             !
             nptoelt1=ptoel2(ip1+1)-ptoel2(ip1)
             nptoelt2=ptoel2(ip2+1)-ptoel2(ip2)

             if(nptoelt1<nptoelt2)then 
                ipmin=ip1
             else
                ipmin=ip2
             endif
             !
             !     Loop on the elements surrounding ipoin
             !
             do jptoel=ptoel2(ipmin),ptoel2(ipmin+1)-1
                jface=ptoel1(jptoel)
                if(jface/=iface)then
                   jp1=lface(1,jface)
                   jp2=lface(2,jface)
                   jp3=lface(3,jface)

                   if(lmark(jp1)+lmark(jp2)+lmark(jp3)==2_ip)then                   
                      do jnode=1,3
                         jp1=lface(ltab(1,jnode),jface)
                         jp2=lface(ltab(2,jnode),jface)

                         if(lmark(jp1)==1 .and. lmark(jp2)==1) then
                            !
                            !     The scalar product should be positive
                            !
                            csca=rnofa(1,iface)*rnofa(1,jface)+&
                                 rnofa(2,iface)*rnofa(2,jface)+&
                                 rnofa(3,iface)*rnofa(3,jface)

                            if(csca>c00)then 

                            eltoel(inode,iface)=jface
                            eltoel(jnode,jface)=iface
                            goto 100
                            endif

                         endif
                      enddo
                   endif
                endif
             enddo
100          continue 
             !
             !     Clean up lmark
             !       
             lmark(ip1)=0
             lmark(ip2)=0 

          endif
       enddo
    enddo

    call memchk(2_ip,istat,memor_msh,'LMARK','trtotr',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','trtotr',0_ip)

  end subroutine trtotr2

  subroutine tetote(elem,nnode,nelem,ptoel1,ptoel2,npoin,eltoel)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip),intent(in)     :: nnode,nelem,npoin 
    integer(ip),intent(in)     :: elem(nnode,nelem)
    integer(ip),intent(in)     :: ptoel1(*),ptoel2(npoin+1) 
    integer(ip),pointer        :: eltoel(:,:)
    integer(4)                 :: istat
    integer(ip)                :: ltab(3,4),ielem,jelem,jp1,jp2,jp3,jnode,inode,ip1,ip2,ip3,ipmin
    integer(ip)                :: nptoelt1,nptoelt2,nptoelt3,jptoel,npmin,jp4
    integer(ip),pointer        :: lmark(:)

    ltab(1,1)=2
    ltab(2,1)=3
    ltab(3,1)=4

    ltab(1,2)=1
    ltab(2,2)=4
    ltab(3,2)=3

    ltab(1,3)=1
    ltab(2,3)=2
    ltab(3,3)=4

    ltab(1,4)=1
    ltab(2,4)=3
    ltab(3,4)=2

    !
    !     Allocate  eltoel
    !   
    if(.not.associated(eltoel))then 
       allocate(eltoel(nnode,nelem),stat=istat)
       call memchk(zero,istat,memor_msh,'ELTOEL','tetote',eltoel)
    else
       call memrea(nelem,memor_msh,'ELTOEL','tetote',eltoel)
       do ielem=1,nelem
          do inode=1,nnode 
             eltoel(inode,ielem)=0_ip
          enddo 
       enddo 
    endif
    !
    !     Allocate help array
    !
    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','tetote',lmark)

    do ielem=1,nelem
       do inode=1,4
          if(eltoel(inode,ielem)==0_ip)then 
             ! 
             !     Mark lmark 
             !        
             ip1=elem(ltab(1,inode),ielem)
             lmark(ip1)=1_ip
             ip2=elem(ltab(2,inode),ielem)
             lmark(ip2)=1_ip
             ip3=elem(ltab(3,inode),ielem)
             lmark(ip3)=1_ip
             !
             !     Take the point with the least number of surrounding elements
             !
             nptoelt1=ptoel2(ip1+1)-ptoel2(ip1)
             nptoelt2=ptoel2(ip2+1)-ptoel2(ip2)
             nptoelt3=ptoel2(ip3+1)-ptoel2(ip3)

             if(nptoelt1<nptoelt2)then 
                ipmin=ip1
                npmin=nptoelt1
             else
                ipmin=ip2
                npmin=nptoelt2
             endif

             if(nptoelt3<npmin)then
                ipmin=ip3
             endif


             !
             !     Loop on the elements surrounding ipoin
             !
             do jptoel=ptoel2(ipmin),ptoel2(ipmin+1)-1
                jelem=ptoel1(jptoel)
                if(jelem/=ielem)then
                   jp1=elem(1,jelem) 
                   jp2=elem(2,jelem) 
                   jp3=elem(3,jelem) 
                   jp4=elem(4,jelem)
                   if(lmark(jp1)+lmark(jp2)+lmark(jp3)+lmark(jp4)==3)then
                      do jnode=1,4
                         jp1=elem(ltab(1,jnode),jelem)
                         jp2=elem(ltab(2,jnode),jelem)
                         jp3=elem(ltab(3,jnode),jelem)

                         if(lmark(jp1)==1 .and. lmark(jp2)==1 .and. lmark(jp3)==1)then

                            eltoel(inode,ielem)=jelem
                            eltoel(jnode,jelem)=ielem
                            goto 100

                         endif
                      enddo
                   endif
                endif
             enddo
100          continue          

             !
             !     Clean up lmark
             !       
             lmark(ip1)=0_ip
             lmark(ip2)=0_ip
             lmark(ip3)=0_ip

          endif
       enddo
    enddo

    call memchk(2_ip,istat,memor_msh,'LMARK','tetote',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','tetote',0_ip)

  end subroutine tetote  

  subroutine tetoteb(elem,nnode,nelem,ptoel1,ptoel2,npoin,eltoel)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip),intent(in)     :: nnode,nelem,npoin 
    integer(ip),intent(in)     :: elem(nnode,nelem)
    integer(ip),intent(in)     :: ptoel1(*),ptoel2(npoin+1) 
    integer(ip),pointer        :: eltoel(:,:)
    integer(4)                 :: istat
    integer(ip)                :: ltab(3,4),ielem,jelem,jp1,jp2,jp3,jnode,inode,ip1,ip2,ip3,ipmin
    integer(ip)                :: nptoelt1,nptoelt2,nptoelt3,jptoel,npmin,jp4
    integer(ip),pointer        :: lmark(:)

    ltab(1,1)=2
    ltab(2,1)=3
    ltab(3,1)=4

    ltab(1,2)=1
    ltab(2,2)=4
    ltab(3,2)=3

    ltab(1,3)=1
    ltab(2,3)=2
    ltab(3,3)=4

    ltab(1,4)=1
    ltab(2,4)=3
    ltab(3,4)=2

    !
    !     Allocate  eltoel
    !   
    if(.not.associated(eltoel))then 
       allocate(eltoel(nnode,nelem),stat=istat)
       call memchk(zero,istat,memor_msh,'ELTOEL','tetote',eltoel)
    else
       call memrea(nelem,memor_msh,'ELTOEL','tetote',eltoel)
       do ielem=1,nelem
          do inode=1,nnode
             eltoel(inode,ielem)=0_ip
          enddo
       enddo 
    endif
    !
    !     Allocate help array
    !
    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','tetote',lmark)

    do ielem=1,nelem
       if(elem(1,ielem)==0)cycle
       do inode=1,4
          if(eltoel(inode,ielem)==0_ip)then 
             ! 
             !     Mark lmark 
             !        
             ip1=elem(ltab(1,inode),ielem)
             lmark(ip1)=1_ip
             ip2=elem(ltab(2,inode),ielem)
             lmark(ip2)=1_ip
             ip3=elem(ltab(3,inode),ielem)
             lmark(ip3)=1_ip
             !
             !     Take the point with the least number of surrounding elements
             !
             nptoelt1=ptoel2(ip1+1)-ptoel2(ip1)
             nptoelt2=ptoel2(ip2+1)-ptoel2(ip2)
             nptoelt3=ptoel2(ip3+1)-ptoel2(ip3)

             if(nptoelt1<nptoelt2)then 
                ipmin=ip1
                npmin=nptoelt1
             else
                ipmin=ip2
                npmin=nptoelt2
             endif

             if(nptoelt3<npmin)then
                ipmin=ip3
             endif
             !
             !     Loop on the elements surrounding ipoin
             !
             do jptoel=ptoel2(ipmin),ptoel2(ipmin+1)-1
                jelem=ptoel1(jptoel)
                if(jelem/=ielem)then
                   jp1=elem(1,jelem) 
                   jp2=elem(2,jelem) 
                   jp3=elem(3,jelem) 
                   jp4=elem(4,jelem)
                   if(lmark(jp1)+lmark(jp2)+lmark(jp3)+lmark(jp4)==3)then
                      do jnode=1,4
                         jp1=elem(ltab(1,jnode),jelem)
                         jp2=elem(ltab(2,jnode),jelem)
                         jp3=elem(ltab(3,jnode),jelem)

                         if(lmark(jp1)==1 .and. lmark(jp2)==1 .and. lmark(jp3)==1)then

                            eltoel(inode,ielem)=jelem
                            eltoel(jnode,jelem)=ielem
                            goto 100

                         endif
                      enddo
                   endif
                endif
             enddo
100          continue          

             !
             !     Clean up lmark
             !       
             lmark(ip1)=0_ip
             lmark(ip2)=0_ip
             lmark(ip3)=0_ip

          endif
       enddo
    enddo

    call memchk(2_ip,istat,memor_msh,'LMARK','tetote',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','tetote',0_ip)

  end subroutine tetoteb

  subroutine gtbfa(nface,nnofa,nelem,elem,nnode,eltoel)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       :memor_msh,lface,rnofa
    implicit none
    integer(4)                 :: istat
    integer(ip),intent(in)     :: nnofa,nelem,nnode
    integer(ip),intent(in)     :: elem(nnode,nelem),eltoel(nnode,nelem)
    integer(ip),intent(inout)  :: nface
    integer(ip)                :: ltab(3,4),ielem,inode,ip1,ip2,ip3,iface,inofa

    ltab(1,1)=2
    ltab(2,1)=3
    ltab(3,1)=4

    ltab(1,2)=1
    ltab(2,2)=4
    ltab(3,2)=3

    ltab(1,3)=1
    ltab(2,3)=2
    ltab(3,3)=4

    ltab(1,4)=1
    ltab(2,4)=3
    ltab(3,4)=2

    !
    !     Count external faces
    !
    nface=0_ip
    do ielem=1,nelem
       do inode=1,nnode 
          if(eltoel(inode,ielem)==0)then
             nface=nface+1
          endif
       enddo
    enddo
    !
    !     Allocate lface && rnofa
    !
    if(.not.associated(lface))then 
       allocate(lface(nnofa,nface),stat=istat)
       call memchk(zero,istat,memor_msh,'LFACE','gtbfa',lface)
    else
       call memrea(nface,memor_msh,'LFACE','gtbfa',lface)
       do iface=1,nface
          do inofa=1,nnofa 
             lface(inofa,iface)=0_ip
          enddo 
       enddo 
    endif

    if(.not.associated(rnofa))then 
       allocate(rnofa(nnofa,nface),stat=istat)
       call memchk(zero,istat,memor_msh,'RFACE','gtbfa',rnofa)
    else
       call memrea(nface,memor_msh,'RFACE','gtbfa',rnofa)
       do iface=1,nface
          do inofa=1,nnofa 
             rnofa(inofa,iface)=0_rp
          enddo 
       enddo 
    endif

    nface=0_ip
    do ielem=1,nelem
       do inode=1,nnode 
          if(eltoel(inode,ielem)==0)then

             ip1=elem(ltab(inode,1),ielem)
             ip2=elem(ltab(inode,2),ielem)
             ip3=elem(ltab(inode,3),ielem)

             nface=nface+1
             lface(1,nface)=ip1
             lface(2,nface)=ip2
             lface(3,nface)=ip3

          endif
       enddo
    enddo


  end subroutine gtbfa

  subroutine swapfa(lface,nface,nnofa,coor,ndim,npoin,eltoel,rface)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       :memor_msh
    implicit none
    integer(4)                 :: istat
    integer(ip),intent(in)     :: nface,nnofa,npoin,ndim
    integer(ip),intent(inout)     :: lface(nnofa,nface)
    integer(ip),intent(inout)     :: eltoel(nnofa,nface)
    real(rp),intent(in)        :: coor(ndim,npoin)
    real(rp),pointer::rface(:,:)
    integer(ip)                :: ltab(2,3),ip1,ip2,ip3,ip4,inode,iface
    integer(ip)                :: ineigh1,ineigh2,jneigh1,jneigh2,jface,jnode,iopt 
    real(rp)                   :: dx12,dy12,dz12,dx13,dy13,dz13,rface1(3),rface2(3)
    real(rp)                   :: cscal1,cscal2,rl 
    integer(ip),pointer        :: lmark(:,:)


    ltab(1,1)=2
    ltab(2,1)=3
    ltab(1,2)=3
    ltab(2,2)=1
    ltab(1,3)=1
    ltab(2,3)=2



    !
    !     First compute the normal for all the faces
    !
    allocate(rface(ndim,nface),stat=istat)
    call memchk(zero,istat,memor_msh,'RFACE','swapfa',rface)

    do iface=1,nface

       ip1=lface(1,iface)
       ip2=lface(2,iface)
       ip3=lface(3,iface)

       dx12=coor(1,ip2)-coor(1,ip1)
       dy12=coor(2,ip2)-coor(2,ip1)
       dz12=coor(3,ip2)-coor(3,ip1)

       dx13=coor(1,ip3)-coor(1,ip1)
       dy13=coor(2,ip3)-coor(2,ip1)
       dz13=coor(3,ip3)-coor(3,ip1)

       rface1(1)= dy12*dz13-dz12*dy13
       rface1(2)=-dx12*dz13+dz12*dx13
       rface1(3)= dx12*dy13-dy12*dx13

       rl=sqrt(rface1(1)*rface1(1)+rface1(2)*rface1(2)+rface1(3)*rface1(3))
       rl=1_rp/rl

       rface(1,iface)=rface1(1)*rl   
       rface(2,iface)=rface1(2)*rl   
       rface(3,iface)=rface1(3)*rl   


    enddo
    !
    !     Allocate lmark
    !
    allocate(lmark(nnofa,nface),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','swapfa',lmark)


    !
    !     Loop on mesh until nothing happens
    !
    open_loop:do

       iopt=0_ip

       do iface=1,nface
          do inode=1,3
             !
             !     Is this edge swappable
             !
             if(lmark(inode,iface)==0_ip)then  
                !
                !     Neighboring element
                !               
                jface=eltoel(inode,iface)

                if(eltoel(1,jface)==iface)then 
                   jnode=1
                else if (eltoel(2,jface)==iface)then
                   jnode=2
                else
                   jnode=3
                endif
                !
                !     Scalar product
                !
                cscal1=rface(1,iface)*rface(1,jface)+rface(2,iface)*rface(2,jface)+rface(3,iface)*rface(3,jface) 
                !
                !     Compute alternate configuration
                !                 
                ip1=lface(inode,iface) 
                ip2=lface(ltab(1,inode),iface) 
                ip3=lface(ltab(2,inode),iface) 
                ip4=lface(jnode,jface) 

                dx12=coor(1,ip2)-coor(1,ip1)
                dy12=coor(2,ip2)-coor(2,ip1)
                dz12=coor(3,ip2)-coor(3,ip1)

                dx13=coor(1,ip4)-coor(1,ip1)
                dy13=coor(2,ip4)-coor(2,ip1)
                dz13=coor(3,ip4)-coor(3,ip1)

                rface1(1)= dy12*dz13-dz12*dy13
                rface1(2)=-dx12*dz13+dz12*dx13
                rface1(3)= dx12*dy13-dy12*dx13

                rl=sqrt(rface1(1)*rface1(1)+rface1(2)*rface1(2)+rface1(3)*rface1(3))
                rl=1_rp/rl
                rface1(1)=rface1(1)*rl   
                rface1(2)=rface1(2)*rl   
                rface1(3)=rface1(3)*rl   



                dx12=dx13
                dy12=dy13
                dz12=dz13

                dx13=coor(1,ip3)-coor(1,ip1)
                dy13=coor(2,ip3)-coor(2,ip1)
                dz13=coor(3,ip3)-coor(3,ip1)

                rface2(1)= dy12*dz13-dz12*dy13
                rface2(2)=-dx12*dz13+dz12*dx13
                rface2(3)= dx12*dy13-dy12*dx13 

                rl=sqrt(rface2(1)*rface2(1)+rface2(2)*rface2(2)+rface2(3)*rface2(3))
                rl=1_rp/rl
                rface2(1)=rface2(1)*rl   
                rface2(2)=rface2(2)*rl   
                rface2(3)=rface2(3)*rl 

                cscal2=rface1(1)*rface2(1)+rface1(2)*rface2(2)+rface1(3)*rface2(3) 
                !
                !     Test          
                !     
                if(cscal2>cscal1)then
                   !
                   !     Remember iopt  
                   !
                   iopt=1_ip
                   !
                   !     Optimization successfull
                   ! 
                   ineigh1=eltoel(ltab(1,inode),iface) 
                   ineigh2=eltoel(ltab(2,inode),iface) 
                   jneigh1=eltoel(ltab(1,jnode),jface) 
                   jneigh2=eltoel(ltab(2,jnode),jface) 

                   lface(1,iface)=ip1
                   lface(2,iface)=ip2
                   lface(3,iface)=ip4

                   lface(1,jface)=ip1
                   lface(2,jface)=ip4
                   lface(3,jface)=ip3

                   eltoel(1,iface)=jneigh1
                   eltoel(2,iface)=jface
                   eltoel(3,iface)=ineigh2

                   eltoel(1,jface)=jneigh2
                   eltoel(2,jface)=ineigh1
                   eltoel(3,jface)=iface
                   !
                   !     Update surrounding elements 
                   !
                   if(eltoel(1,ineigh1)==iface)then
                      eltoel(1,ineigh1)=jface
                      lmark(1,ineigh1)=0_ip
                   else if(eltoel(2,ineigh1)==iface)then
                      eltoel(2,ineigh1)=jface
                      lmark(2,ineigh1)=0_ip
                   else
                      eltoel(3,ineigh1)=jface
                      lmark(3,ineigh1)=0_ip
                   endif

                   if(eltoel(1,jneigh1)==jface)then
                      eltoel(1,jneigh1)=iface
                      lmark(1,jneigh1)=0_ip
                   else if(eltoel(2,jneigh1)==jface)then
                      eltoel(2,jneigh1)=iface
                      lmark(2,jneigh1)=0_ip
                   else
                      eltoel(3,jneigh1)=iface
                      lmark(3,jneigh1)=0_ip
                   endif

                   !
                   !     Mark all surrounding elements
                   !                 
                   lmark(1,iface)=0_ip
                   lmark(2,iface)=0_ip
                   lmark(3,iface)=0_ip
                   lmark(1,jface)=0_ip
                   lmark(2,jface)=0_ip
                   lmark(3,jface)=0_ip

                   if(eltoel(1,jneigh2)==jface)then
                      lmark(1,jneigh2)=0_ip
                   else if(eltoel(2,jneigh2)==jface)then
                      lmark(2,jneigh2)=0_ip
                   else
                      lmark(3,jneigh2)=0_ip
                   endif

                   if(eltoel(1,ineigh2)==iface)then
                      lmark(1,ineigh2)=0_ip
                   else if(eltoel(2,ineigh2)==iface)then
                      lmark(2,ineigh2)=0_ip
                   else
                      lmark(3,ineigh2)=0_ip
                   endif

                else    
                   !
                   !     Nothing to do, mark lmark
                   ! 
                   lmark(inode,iface)=1  
                   lmark(jnode,jface)=1 

                endif

             endif
          enddo
       enddo

       if(iopt==0_ip)exit open_loop 


    enddo open_loop

    call memchk(2_ip,istat,memor_msh,'LMARK','swapfa',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','swapfa',0_ip)


  end subroutine swapfa

  subroutine fatoed(nedge,nface,lface,nnofa,eltoel,ledge,lftoed)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only           : memor_msh 
    implicit none
    integer(ip),intent(in)    ::nface,nnofa
    integer(ip),intent(in)    ::lface(nnofa,nface)
    integer(ip),intent(inout) ::nedge,eltoel(nnofa,nface)
    integer(ip),pointer       :: ledge(:,:),lftoed(:,:)
    integer(ip)               ::iface,inode,jface,ltab(2,3),iedge
    integer(4)                :: istat

    ltab(1,1)=2
    ltab(2,1)=3
    ltab(1,2)=3
    ltab(2,2)=1
    ltab(1,3)=1
    ltab(2,3)=2
    !
    !     Count the number of edges
    !
    nedge=0_ip

    do iface=1,nface
       do inode=1,3 
          jface=eltoel(inode,iface)
          if(jface>0)then

             eltoel(inode,iface)=-jface   

             if(eltoel(1,jface)==iface)then
                eltoel(1,jface)=-iface
             else if(eltoel(2,jface)==iface)then
                eltoel(2,jface)=-iface
             else
                eltoel(3,jface)=-iface
             endif

             nedge=nedge+1_ip
          else if(jface==0)then
             nedge=nedge+1_ip
          endif
       enddo
    enddo
    !
    !     Allocate the edge array and the edge to-face-array
    !
    if(.not.associated(ledge))then 
       allocate(ledge(2,nedge),stat=istat)
       call memchk(zero,istat,memor_msh,'LEDGE','fatoed',ledge)
    else
       call memrea(nedge,memor_msh,'LEDGE','fatoed',ledge)
       do iedge=1,nedge
          ledge(1,iedge)=0_ip
          ledge(2,iedge)=0_ip
       enddo
    endif

    if(.not.associated(lftoed))then 
       allocate(lftoed(nnofa,nface),stat=istat)
       call memchk(zero,istat,memor_msh,'LFTOED','fatoed',lftoed)
    else
       call memrea(nface,memor_msh,'LFTOED','fatoed',lftoed)
       do iface=1,nface
          lftoed(1,iface)=0_ip
          lftoed(2,iface)=0_ip
          lftoed(3,iface)=0_ip
       enddo
    endif
    !
    !     Fill ledge
    !
    nedge=0_ip
    do iface=1,nface

       do inode=1,3 
          jface=-eltoel(inode,iface)
          if(jface>0)then

             eltoel(inode,iface)=jface   

             nedge=nedge+1_ip
             ledge(1,nedge)=lface(ltab(1,inode),iface)
             ledge(2,nedge)=lface(ltab(2,inode),iface)
             lftoed(inode,iface)=nedge 

             if(eltoel(1,jface)==-iface)then
                eltoel(1,jface)=iface
                lftoed(1,jface)=nedge 
             else if(eltoel(2,jface)==-iface)then
                eltoel(2,jface)=iface
                lftoed(2,jface)=nedge 
             else
                eltoel(3,jface)=iface
                lftoed(3,jface)=nedge 
             endif

          else if(jface==0)then

             nedge=nedge+1_ip
             ledge(1,nedge)=lface(ltab(1,inode),iface)
             ledge(2,nedge)=lface(ltab(2,inode),iface)
             lftoed(inode,iface)=nedge 

          endif
       enddo
    enddo

  end subroutine fatoed

  subroutine quatri(nx,ny,ndim,npoin,nface,nnofa,coor,lface)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       :memor_msh
    implicit none
    integer(ip),intent(in)   :: nx,ny,ndim,npoin,nnofa
    real(rp), intent(in)     :: coor(ndim,npoin) 
    integer(ip),intent(inout):: nface
    integer(ip),pointer      :: lface(:,:)
    integer(ip)              :: ip1,ip2,ip3,ip4,ix,iy,nfacnew,iface,inofa
    real(rp)                 :: dx12,dy12,dz12,dx13,dy13,dz13,rface1(3),rface2(3) 
    real(rp)                 :: cscal1,cscal2,rl,c10
    integer(4)                 :: istat
    !
    !     This sub converts the quads from the terrain to triangles
    !
    c10=1.0d+00
    !
    !     Compute the number of faces
    !
    nfacnew=(nx-1)*(ny-1)*2

    if(.not.associated(lface))then 
       allocate(lface(nnofa,nfacnew),stat=istat)
       call memchk(zero,istat,memor_msh,'LFACE','quatri',lface)
    else
       call memrea(nfacnew,memor_msh,'LFACE','quatri',lface)
       do iface=1,nface
          do inofa=1,nnofa
             lface(inofa,iface)=0_ip
          enddo
       enddo
    endif

    nface=0_ip
    do iy=1,ny-1

       do ix=1,nx-1
          !
          !     Gather points
          !
          ip1=ix+(iy-1)*nx
          ip2=ip1+1
          ip3=ip1+nx
          ip4=ip3+1 
          !
          !     Check both configurations
          ! 
          dx12=coor(1,ip2)-coor(1,ip1)
          dy12=coor(2,ip2)-coor(2,ip1)
          dz12=coor(3,ip2)-coor(3,ip1)

          dx13=coor(1,ip4)-coor(1,ip1)
          dy13=coor(2,ip4)-coor(2,ip1)
          dz13=coor(3,ip4)-coor(3,ip1)

          rface1(1)= dy12*dz13-dz12*dy13
          rface1(2)=-dx12*dz13+dz12*dx13
          rface1(3)= dx12*dy13-dy12*dx13

          rl=sqrt(rface1(1)*rface1(1)+rface1(2)*rface1(2)+rface1(3)*rface1(3))
          rl=1_rp/rl
          rface1(1)=rface1(1)*rl   
          rface1(2)=rface1(2)*rl   
          rface1(3)=rface1(3)*rl   

          dx12=dx13
          dy12=dy13
          dz12=dz13

          dx13=coor(1,ip3)-coor(1,ip1)
          dy13=coor(2,ip3)-coor(2,ip1)
          dz13=coor(3,ip3)-coor(3,ip1)

          rface2(1)= dy12*dz13-dz12*dy13
          rface2(2)=-dx12*dz13+dz12*dx13
          rface2(3)= dx12*dy13-dy12*dx13 

          rl=sqrt(rface2(1)*rface2(1)+rface2(2)*rface2(2)+rface2(3)*rface2(3))
          rl=c10/rl
          rface2(1)=rface2(1)*rl   
          rface2(2)=rface2(2)*rl   
          rface2(3)=rface2(3)*rl 

          cscal1=rface1(1)*rface2(1)+rface1(2)*rface2(2)+rface1(3)*rface2(3) 
          !
          !     Second configuration
          !
          dx12=coor(1,ip1)-coor(1,ip3)
          dy12=coor(2,ip1)-coor(2,ip3)
          dz12=coor(3,ip1)-coor(3,ip3)

          dx13=coor(1,ip2)-coor(1,ip3)
          dy13=coor(2,ip2)-coor(2,ip3)
          dz13=coor(3,ip2)-coor(3,ip3)

          rface1(1)= dy12*dz13-dz12*dy13
          rface1(2)=-dx12*dz13+dz12*dx13
          rface1(3)= dx12*dy13-dy12*dx13

          rl=sqrt(rface1(1)*rface1(1)+rface1(2)*rface1(2)+rface1(3)*rface1(3))
          rl=c10/rl
          rface1(1)=rface1(1)*rl   
          rface1(2)=rface1(2)*rl   
          rface1(3)=rface1(3)*rl   

          dx12=dx13
          dy12=dy13
          dz12=dz13

          dx13=coor(1,ip4)-coor(1,ip3)
          dy13=coor(2,ip4)-coor(2,ip3)
          dz13=coor(3,ip4)-coor(3,ip3)

          rface2(1)= dy12*dz13-dz12*dy13
          rface2(2)=-dx12*dz13+dz12*dx13
          rface2(3)= dx12*dy13-dy12*dx13 

          rl=sqrt(rface2(1)*rface2(1)+rface2(2)*rface2(2)+rface2(3)*rface2(3))
          rl=1_rp/rl
          rface2(1)=rface2(1)*rl   
          rface2(2)=rface2(2)*rl   
          rface2(3)=rface2(3)*rl 

          cscal2=rface1(1)*rface2(1)+rface1(2)*rface2(2)+rface1(3)*rface2(3) 
          !
          !      And the test        
          !
          if(cscal1>cscal2)then

             nface=nface+1
             lface(1,nface)=ip1
             lface(2,nface)=ip2
             lface(3,nface)=ip4
             nface=nface+1
             lface(1,nface)=ip1
             lface(2,nface)=ip4
             lface(3,nface)=ip3

          else

             nface=nface+1
             lface(1,nface)=ip1
             lface(2,nface)=ip2
             lface(3,nface)=ip3
             nface=nface+1
             lface(1,nface)=ip3
             lface(2,nface)=ip2
             lface(3,nface)=ip4

          endif

       enddo
    enddo

    call outfac(nnofa,nface,npoin,ndim,lface,coor)

  end subroutine quatri

  subroutine gtfnrl(lface,nface,nnofa,ndim,coor,npoin,rnofa)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       :memor_msh
    implicit none
    integer(ip),intent(in)         :: nface,nnofa,ndim,npoin 
    integer(ip),intent(in)         :: lface(nnofa,nface)
    real(rp),intent(in)            :: coor(ndim,npoin) 
    real(rp),intent(out)           :: rnofa(ndim,nface)
    integer(ip)                    :: iface,ip1,ip2,ip3
    real(rp)                       :: dx12,dy12,dz12,dx13,dy13,dz13
    real(rp)                       :: rnx,rny,rnz,rnl,c10   

    c10=1.0d+00
    !
    !     Compute face normal 
    !
    do iface=1,nface

       ip1=lface(1,iface) 
       ip2=lface(2,iface) 
       ip3=lface(3,iface) 

       dx12=coor(1,ip2)-coor(1,ip1)
       dy12=coor(2,ip2)-coor(2,ip1)
       dz12=coor(3,ip2)-coor(3,ip1)

       dx13=coor(1,ip3)-coor(1,ip1)
       dy13=coor(2,ip3)-coor(2,ip1)
       dz13=coor(3,ip3)-coor(3,ip1)

       rnx= dy12*dz13-dz12*dy13
       rny=-dx12*dz13+dz12*dx13
       rnz= dx12*dy13-dy12*dx13

       rnl=sqrt(rnx*rnx+rny*rny+rnz*rnz)
       rnl=c10/rnl

       rnofa(1,iface)=rnx*rnl
       rnofa(2,iface)=rny*rnl
       rnofa(3,iface)=rnz*rnl

    enddo

  end subroutine gtfnrl

  subroutine gtfnr2(lface,nface,nnofa,ndim,coor,npoin,rnofa,iface,rnl)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       :memor_msh
    implicit none
    integer(ip),intent(in)  :: nface,nnofa,ndim,npoin,iface 
    integer(ip),intent(in)  :: lface(nnofa,nface)
    real(rp),intent(in)     :: coor(ndim,npoin)
    real(rp),intent(inout)  :: rnl,rnofa(ndim)
    integer(ip)             :: ip1,ip2,ip3
    real(rp)                :: dx12,dy12,dz12,dx13,dy13,dz13,epsil,c00
    real(rp)                :: rnx,rny,rnz,c10,rnl2,rn1,rn2,tol   

    c10=1.0d+00
    c00=0.0d+00
    epsil=1.0d-12
    !
    !     Compute face normal 
    !

    ip1=lface(1,iface) 
    ip2=lface(2,iface) 
    ip3=lface(3,iface) 

    dx12=coor(1,ip2)-coor(1,ip1)
    dy12=coor(2,ip2)-coor(2,ip1)
    dz12=coor(3,ip2)-coor(3,ip1)

    dx13=coor(1,ip3)-coor(1,ip1)
    dy13=coor(2,ip3)-coor(2,ip1)
    dz13=coor(3,ip3)-coor(3,ip1)

    rnx= dy12*dz13-dz12*dy13
    rny=-dx12*dz13+dz12*dx13
    rnz= dx12*dy13-dy12*dx13

    rn1=abs(dx12)+abs(dy12)+abs(dz12)
    rn2=abs(dx13)+abs(dy13)+abs(dz13)
    tol=epsil*max(rn1,rn2) 

    rnl=sqrt(rnx*rnx+rny*rny+rnz*rnz)

    if(rnl<tol)then

       rnofa(1)=c00
       rnofa(2)=c00
       rnofa(3)=c00
       rnl=c00

    else

       rnl2=c10/rnl
       rnofa(1)=rnx*rnl2
       rnofa(2)=rny*rnl2
       rnofa(3)=rnz*rnl2

    endif

  end subroutine gtfnr2

  subroutine gtfnrl3(lface,nface,nnofa,ndim,coor,npoin,rnofa)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       :memor_msh
    implicit none
    integer(ip),intent(in)         :: nface,nnofa,ndim,npoin 
    integer(ip),intent(in)      :: lface(nnofa,nface)
    real(rp),intent(in)            :: coor(ndim,npoin) 
    integer(ip)                    :: iface,ip1,ip2,ip3
    real(rp)                       :: rnofa(ndim,nface),dx12,dy12,dz12,dx13,dy13,dz13
    real(rp)                       :: rnx,rny,rnz   
    !
    !     Compute face normal 
    !
    do iface=1,nface

       ip1=lface(1,iface) 
       ip2=lface(2,iface) 
       ip3=lface(3,iface) 

       dx12=coor(1,ip2)-coor(1,ip1)
       dy12=coor(2,ip2)-coor(2,ip1)
       dz12=coor(3,ip2)-coor(3,ip1)

       dx13=coor(1,ip3)-coor(1,ip1)
       dy13=coor(2,ip3)-coor(2,ip1)
       dz13=coor(3,ip3)-coor(3,ip1)

       rnx= dy12*dz13-dz12*dy13
       rny=-dx12*dz13+dz12*dx13
       rnz= dx12*dy13-dy12*dx13

       rnofa(1,iface)=rnx
       rnofa(2,iface)=rny
       rnofa(3,iface)=rnz

    enddo

  end subroutine gtfnrl3

  subroutine gtfnrl4(lface,nface,nnofa,ndim,coor,npoin,rnofa,nthread,mythread)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       :memor_msh
    implicit none
    integer(ip),intent(in)         :: nface,nnofa,ndim,npoin,nthread,mythread 
    integer(ip),intent(in)      :: lface(nnofa,nface)
    real(rp),intent(in)            :: coor(ndim,npoin) 
    integer(ip)                    :: iface,ip1,ip2,ip3,iface0,iface1,nchunk
    real(rp)                       :: rnofa(ndim,nface),dx12,dy12,dz12,dx13,dy13,dz13
    real(rp)                       :: rnx,rny,rnz   
    !
    !     Compute face normal 
    !
    nchunk=(nface+nthread-1_ip)/nthread
    iface0=(mythread-1)*nchunk+1_ip
    iface1=min(nface,mythread*nchunk) 

    do iface=iface0,iface1

       ip1=lface(1,iface) 
       ip2=lface(2,iface) 
       ip3=lface(3,iface) 

       dx12=coor(1,ip2)-coor(1,ip1)
       dy12=coor(2,ip2)-coor(2,ip1)
       dz12=coor(3,ip2)-coor(3,ip1)

       dx13=coor(1,ip3)-coor(1,ip1)
       dy13=coor(2,ip3)-coor(2,ip1)
       dz13=coor(3,ip3)-coor(3,ip1)

       rnx= dy12*dz13-dz12*dy13
       rny=-dx12*dz13+dz12*dx13
       rnz= dx12*dy13-dy12*dx13

       rnofa(1,iface)=rnx
       rnofa(2,iface)=rny
       rnofa(3,iface)=rnz

    enddo

  end subroutine gtfnrl4

  subroutine gtpnrl(nface,rnopo,npoin,ndim,ptoel1,ptoel2,rnofa)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       :memor_msh
    implicit none
    integer(ip),intent(in)         :: nface,ndim,npoin 
    integer(ip),intent(in)         :: ptoel1(*),ptoel2(*)
    real(rp),intent(in)            :: rnofa(ndim,nface)   
    real(rp),intent(out)           :: rnopo(ndim,npoin)   
    integer(ip)                    :: ifac,iface,ipoin        
    real(rp)                       :: rnx,rny,rnz,rnl,c00,c10  

    c00=0.0d+00
    c10=1.0d+00
    !
    !     Compute point normal from averaging face normal
    !
    do ipoin=1,npoin

       rnx=c00 
       rny=c00 
       rnz=c00 

       do ifac=ptoel2(ipoin),ptoel2(ipoin+1)-1
          iface=ptoel1(ifac)
          rnx=rnx+rnofa(1,iface)
          rny=rny+rnofa(2,iface)
          rnz=rnz+rnofa(3,iface)
       enddo

       rnl=sqrt(rnx*rnx+rny*rny+rnz*rnz)
       rnl=c10/rnl

       rnopo(1,ipoin)=rnx*rnl
       rnopo(2,ipoin)=rny*rnl
       rnopo(3,ipoin)=rnz*rnl

    enddo

  end subroutine gtpnrl

  subroutine orient3D(p1,p2,p3,rvol,ndim)
    use def_kintyp, only       :  ip,rp,lg
    integer(ip),intent(in)         :: ndim
    real(rp),intent(in)        :: p1(ndim),p2(ndim),p3(ndim)
    real(rp),intent(inout)     :: rvol

    rvol= p1(1)*(p2(2)*p3(3)-p2(3)*p3(2)) &
         &   - p1(2)*(p2(1)*p3(3)-p2(3)*p3(1)) &
         &   + p1(3)*(p2(1)*p3(2)-p2(2)*p3(1))

  end subroutine orient3D

  subroutine sitosa(lside,nside,npoin,nnosi,ptosi1,ptosi2,sitosi)
    use def_meshin, only       :  ID_CORNER,ID_RIDGE,ID_CUSP,ID_SMOOTH
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       :  memor_msh
    implicit none
    integer(ip),intent(in)     :: nnosi,nside,npoin 
    integer(ip),intent(in)     :: lside(nnosi,nside)
    integer(ip),intent(in)     :: ptosi1(*),ptosi2(npoin+1) 
    integer(ip),pointer        :: sitosi(:,:)
    integer(4)                 :: istat
    integer(ip)                :: iside,jside,ip1,inosi,jnosi,isid1,isid2
    integer(ip)                :: ipa,ipb,nsid

    !
    !     Allocate  sitosi
    !   
    if(.not.associated(sitosi)) then 
       allocate(sitosi(nnosi,nside),stat=istat)
       call memchk(zero,istat,memor_msh,'SITOSI','sitos',sitosi)
    else
       call memrea(nside,memor_msh,'SITOSI','sitos',sitosi)
       do iside=1,nside 
          sitosi(1,iside)=0_ip
          sitosi(2,iside)=0_ip
       enddo
    end if
    !
    !     Loop on the sides
    !  
    do iside=1,nside
       do inosi=1,nnosi
          if(sitosi(inosi,iside)==0_ip)then 
             ip1=lside(inosi,iside)
             !
             !     Look for the neighbor 
             !   
             nsid=ptosi2(ip1+1)-ptosi2(ip1)
             !
             !     Decide with the ridge number in order to 
             !     take into account the corners with two ridges
             !     with a large angle 
             !
             if(nsid/=2)then
                cycle
             endif
             isid1=ptosi1(ptosi2(ip1))
             isid2=ptosi1(ptosi2(ip1)+1)
             if(isid1/=iside)then
                jside=isid1
             else
                jside=isid2
             endif

             ipa=lside(1,jside)
             ipb=lside(2,jside)

             if(ipa==ip1)then
                jnosi=1_ip
             else
                jnosi=2_ip
             endif
             !
             !     Fill sitosi
             !
             sitosi(inosi,iside)=jside
             sitosi(jnosi,jside)=iside

          endif
       enddo
    enddo

  end subroutine sitosa

  subroutine sitofa(lside,nside,nnosi,lface,nnofa,nface,npoin,ptoel1,ptoel2,lstof)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only       : memor_msh
    implicit none
    integer(ip),intent(in)     :: nnosi,nside,npoin,nface,nnofa 
    integer(ip),intent(in)     :: lside(nnosi,nside) 
    integer(ip),intent(in)     :: lface(nnofa,nface)
    integer(ip),pointer        :: lstof(:) 
    integer(ip),intent(in)     :: ptoel1(*),ptoel2(npoin+1) 
    integer(4)                 :: istat
    integer(ip)                :: iside,ip1,ip2,ip3,ip4,ie
    integer(ip)                :: ipa,nsid1,nsid2,iface,inofa

    !
    !     Allocate  sitosi
    !   
    if(.not.associated(lstof))then 
       allocate(lstof(nside),stat=istat)
       call memchk(zero,istat,memor_msh,'LSTOF','sitofa',lstof)
    else
       call memrea(nside,memor_msh,'LSTOF','sitofa',lstof)
       do iside=1,nside 
          lstof(iside)=0_ip
       enddo
    endif

    do iside=1,nside

       ip1=lside(1,iside)
       ip2=lside(2,iside)

       nsid1=ptoel2(ip1+1)-ptoel2(ip1)
       nsid2=ptoel2(ip2+1)-ptoel2(ip2)

       if(nsid1<nsid2)then
          ip3=ip1
          ip4=ip2
       else
          ip3=ip2
          ip4=ip1
       endif

       do ie=ptoel2(ip3),ptoel2(ip3+1)-1
          iface=ptoel1(ie)
          do inofa=1,nnofa
             ipa=lface(inofa,iface)
             if(ipa==ip4)goto 10
          enddo
       enddo

       write(*,*)'Error in sitofa, face not found for point:',ip3
       stop 

10     continue

       lstof(iside)=iface

    enddo

  end subroutine sitofa

  subroutine edglem(lelem,nnode,nelem,ledglm,nsid,ptoed2,ledge,nedge)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)     :: nelem,nnode 
    integer(ip),intent(in)     :: nsid,nedge  
    integer(ip),intent(in)     :: lelem(nnode,nelem),ptoed2(*) 
    integer(ip),intent(in)     :: ledge(2,nedge)  
    integer(ip),intent(inout)  :: ledglm(nsid,nelem)  
    integer(ip)                :: lposi(2,6),ip1,ip2,ip3
    integer(ip)                :: ielem,ipmin,ipmax,iside,isid

    lposi(1,1)=1
    lposi(2,1)=2
    lposi(1,2)=1
    lposi(2,2)=3
    lposi(1,3)=1
    lposi(2,3)=4
    lposi(1,4)=2
    lposi(2,4)=3
    lposi(1,5)=2
    lposi(2,5)=4
    lposi(1,6)=3
    lposi(2,6)=4
    !
    !     Get the edge element pointer for tetrahedra
    !
    do ielem=1,nelem
       do iside=1,nsid  
          ip1=lelem(lposi(1,iside),ielem)
          ip2=lelem(lposi(2,iside),ielem)

          ipmin=min(ip1,ip2)
          ipmax=max(ip1,ip2)

          do isid=ptoed2(ipmin),ptoed2(ipmin+1)-1
             ip3=ledge(2,isid)
             if(ip3==ipmax)then
                ledglm(iside,ielem)=isid
                exit
             endif
          enddo
          if(isid==ptoed2(ipmin+1))then
             write(*,*)'Error in edglem, edge not found '
             write(*,*)ielem,iside
             stop
          endif

       enddo
    enddo

  end subroutine edglem

  subroutine edglem2(lelem,nnode,nelem,ledglm,nsid,ptoed2,ledge,nedge)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)     :: nelem,nnode 
    integer(ip),intent(in)     :: nsid,nedge  
    integer(ip),intent(in)     :: lelem(nnode,nelem),ptoed2(*) 
    integer(ip),intent(in)     :: ledge(2,nedge)  
    integer(ip),intent(inout)  :: ledglm(nsid,nelem)  
    integer(ip)                :: lposi(2,3),ip1,ip2,ip3
    integer(ip)                :: ielem,ipmin,ipmax,iside,isid 

    lposi(1,1)=2
    lposi(2,1)=3
    lposi(1,2)=3
    lposi(2,2)=1
    lposi(1,3)=1
    lposi(2,3)=2
    !
    !     Get the edge element pointer for triangles
    !
    do ielem=1,nelem
       do iside=1,nsid  
          ip1=lelem(lposi(1,iside),ielem)
          ip2=lelem(lposi(2,iside),ielem)

          ipmin=min(ip1,ip2)
          ipmax=max(ip1,ip2)

          do isid=ptoed2(ipmin),ptoed2(ipmin+1)-1
             ip3=ledge(2,isid)
             if(ip3==ipmax)then
                ledglm(iside,ielem)=isid
                exit
             endif
          enddo

          if(isid==ptoed2(ipmin+1))then
             write(*,*)'Error in edglem2, edge not found ' 
             stop
          endif
       enddo
    enddo

  end subroutine edglem2

  subroutine edglemsha(lelem,nnode,nelem,ledglm,nsid,ptoed2,ledge,nedge,&
       ptopn2,ptopn1,npoin)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)     :: nelem,nnode,npoin
    integer(ip),intent(in)     :: nsid,nedge  
    integer(ip),intent(in)     :: lelem(nnode,nelem),ptoed2(*) 
    integer(ip),intent(in)     :: ledge(2,nedge),ptopn2(npoin+1),ptopn1(*)  
    integer(ip),intent(inout)  :: ledglm(nsid,nelem)  
    integer(ip)                :: lposi(2,6),ip1,ip2,ip3
    integer(ip)                :: ielem,ipmin,ipmax,iside,isid
    integer(ip)                :: isto,jsto,jp1,jp2 
    logical(lg)                :: icheck

    lposi(1,1)=1
    lposi(2,1)=2
    lposi(1,2)=1
    lposi(2,2)=3
    lposi(1,3)=1
    lposi(2,3)=4
    lposi(1,4)=2
    lposi(2,4)=3
    lposi(1,5)=2
    lposi(2,5)=4
    lposi(1,6)=3
    lposi(2,6)=4
    !
    !     Get the edge element pointer for tetrahedra
    !
    do ielem=1,nelem
       do iside=1,nsid  
          ip1=lelem(lposi(1,iside),ielem)
          ip2=lelem(lposi(2,iside),ielem)

          ipmin=min(ip1,ip2)
          ipmax=max(ip1,ip2)

          do isid=ptoed2(ipmin),ptoed2(ipmin+1)-1
             ip3=ledge(2,isid)
             if(ip3==ipmax)then
                ledglm(iside,ielem)=isid
                exit
             endif
          enddo
          if(isid==ptoed2(ipmin+1))then

             icheck=.false.
             !
             !     Find the virtual points
             !
             jp1=ip1

             do jsto=ptopn2(ip2),ptopn2(ip2+1)-1_ip
                jp2=ptopn1(jsto)

                ipmin=min(jp1,jp2)
                ipmax=max(jp1,jp2)

                do isid=ptoed2(ipmin),ptoed2(ipmin+1)-1
                   ip3=ledge(2,isid)
                   if(ip3==ipmax)then
                      ledglm(iside,ielem)=isid
                      icheck=.true.
                      goto 100 
                   endif
                enddo
             enddo

             jp2=ip2
             do isto=ptopn2(ip1),ptopn2(ip1+1)-1_ip
                jp1=ptopn1(isto)

                ipmin=min(jp1,jp2)
                ipmax=max(jp1,jp2)

                do isid=ptoed2(ipmin),ptoed2(ipmin+1)-1
                   ip3=ledge(2,isid)
                   if(ip3==ipmax)then
                      ledglm(iside,ielem)=isid
                      icheck=.true.
                      goto 100 
                   endif
                enddo
             enddo

             do isto=ptopn2(ip1),ptopn2(ip1+1)-1_ip
                jp1=ptopn1(isto)
                do jsto=ptopn2(ip2),ptopn2(ip2+1)-1_ip
                   jp2=ptopn1(jsto)

                   ipmin=min(jp1,jp2)
                   ipmax=max(jp1,jp2)

                   do isid=ptoed2(ipmin),ptoed2(ipmin+1)-1
                      ip3=ledge(2,isid)
                      if(ip3==ipmax)then
                         ledglm(iside,ielem)=isid
                         icheck=.true.
                         goto 100 
                      endif
                   enddo
                enddo
             enddo

100          continue

             if(icheck.eqv. .false.)then 
                write(*,*)'Error in edglemsha, edge not found ' 
                stop
             endif
          endif

       enddo
    enddo

  end subroutine edglemsha

  subroutine outgidvolume(npoin,ndimn,coor,intmat,nelem,nnode)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)  :: nnode,nelem,ndimn,npoin    
    integer(ip),intent(in)  :: intmat(nnode,nelem)
    real(rp),intent(in)     ::  coor(ndimn,npoin)
    integer(ip)             :: ipoin,ip1,ip2,ip3,ip4,icont,iele
    real(rp)                :: rx,ry,rz

    open(unit=50,file='volume.msh',status='unknown')
    rewind 50
1   format('MESH dimension 3 ElemType Tetrahedra Nnode 4')
2   format('Coordinates')
3   format('#node number   coor_x   coor_y  coor_z')
100 format(i10,3e20.10)
200 format(5i10)
4   format('end coordinates')
5   format('Elements')
6   format('end elements')
    write(50,1)
    write(50,2)
    write(50,3)
    do ipoin=1,npoin
       rx=coor(1,ipoin)
       ry=coor(2,ipoin)
       rz=coor(3,ipoin)
       write(50,100)ipoin,rx,ry,rz
    end do
    write(50,4)
    write(50,5)
    icont=0
    do iele=1,nelem
       ip1=intmat(1,iele)
       ip2=intmat(2,iele)
       ip3=intmat(3,iele)
       ip4=intmat(4,iele)
       icont=icont+1
       write(50,200)icont,ip1,ip2,ip3,ip4
    end do
    write(50,6)
    close(50)
  end subroutine outgidvolume

  subroutine outgidvolume2(npoin,ndimn,coor,intmat,nelem,nnode,coold,npold,intmold,neold)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)  :: nnode,nelem,ndimn,npoin,npold,neold    
    integer(ip),intent(in)  :: intmat(nnode,nelem),intmold(nnode,neold)
    real(rp),intent(in)     ::  coor(ndimn,npoin),coold(ndimn,npold)
    integer(ip)             :: ipoin,ip1,ip2,ip3,ip4,icont,iele
    real(rp)                :: rx,ry,rz

    open(unit=50,file='volume2.msh',status='unknown')
    rewind 50
1   format('MESH dimension 3 ElemType Tetrahedra Nnode 4')
2   format('Coordinates')
3   format('#node number   coor_x   coor_y  coor_z')
100 format(i10,3e20.10)
200 format(6i10)
4   format('end coordinates')
5   format('Elements')
6   format('end elements')
    write(50,1)
    write(50,2)
    write(50,3)
    do ipoin=1,npoin
       rx=coor(1,ipoin)
       ry=coor(2,ipoin)
       rz=coor(3,ipoin)
       write(50,100)ipoin,rx,ry,rz
    end do
    do ipoin=1,npold
       rx=coold(1,ipoin)
       ry=coold(2,ipoin)
       rz=coold(3,ipoin)
       write(50,100)ipoin+npoin,rx,ry,rz
    end do
    write(50,4)
    write(50,5)
    icont=0
    do iele=1,nelem
       ip1=intmat(1,iele)
       ip2=intmat(2,iele)
       ip3=intmat(3,iele)
       ip4=intmat(4,iele)
       icont=icont+1
       write(50,200)icont,ip1,ip2,ip3,ip4,1_ip
    end do
    do iele=1,neold
       ip1=intmold(1,iele)+npoin
       ip2=intmold(2,iele)+npoin
       ip3=intmold(3,iele)+npoin
       ip4=intmold(4,iele)+npoin
       icont=icont+1
       write(50,200)icont,ip1,ip2,ip3,ip4,2_ip
    end do
    write(50,6)
    close(50)
  end subroutine outgidvolume2

  subroutine divmsh(nelem,npoin,nnode,nface,nnofa)
    use def_kintyp, only       :  ip,rp,lg
    use def_meshin, only          : coor,elem,ptoel1,ptoel2,memor_msh,lface,lsurf
    use mod_memchk
    implicit none
    integer(ip),intent(in)     :: nnode,nnofa
    integer(ip),intent(inout)  :: nelem,npoin,nface 
    integer(ip), pointer       :: lcoun(:),lremem(:),lmark(:)
    integer(ip)                :: ip1,ip2,ip3,ip4,ipoin,jpoin,lpedg(6)
    integer(ip)                :: npoin0,iel,ielem,nedg,tab(2,6),ipo
    integer(ip)                :: inode,nelemn,iedg,ipa,ipb,ipmax,ipmin,jpmax
    integer(ip)                :: iface,tab2(2,3),nfacen,isurf
    integer(4)                 :: istat
    real(rp)                   :: c05
    c05=0.5d+00  

    tab(1,1)=1  
    tab(2,1)=2  
    tab(1,2)=1  
    tab(2,2)=3  
    tab(1,3)=1  
    tab(2,3)=4  
    tab(1,4)=2  
    tab(2,4)=3  
    tab(1,5)=2  
    tab(2,5)=4  
    tab(1,6)=3  
    tab(2,6)=4  

    tab2(1,1)=1
    tab2(2,1)=2
    tab2(1,2)=2
    tab2(2,2)=3
    tab2(1,3)=3
    tab2(2,3)=1
    !
    !     Allocate help array
    !
    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','divmsh',lmark)
    allocate(lcoun(npoin+1),stat=istat)
    call memchk(zero,istat,memor_msh,'LCOUN','divmsh',lcoun)
    !
    !     First get the elements surrounding the points
    !
    call ptoelm(elem,nelem,npoin,nnode,ptoel1,ptoel2 )
    !
    !     Count the edges
    !
    npoin0=npoin
    do ipoin=1,npoin0
       lmark(ipoin)=1_ip
       do iel=ptoel2(ipoin),ptoel2(ipoin+1)-1
          ielem=ptoel1(iel)
          do inode=1,nnode
             jpoin=elem(inode,ielem)
             if(lmark(jpoin)==0)then
                lmark(jpoin)=1_ip
                if(jpoin>ipoin)then
                   npoin=npoin+1
                endif
             endif
          enddo
       enddo
    enddo
    !
    !     Reallocate the points
    !
    call memrea(npoin,memor_msh,'COOR','divmsh',coor)
    !
    !     Allocate lremem
    ! 
    nedg=npoin-npoin0 
    allocate(lremem(nedg),stat=istat)
    call memchk(zero,istat,memor_msh,'LREMEM','divmsh',lremem)
    !
    !     Store the new points
    !
    do ipoin=1,npoin
       lmark(ipoin)=0_ip
    enddo
    lcoun(1)=1_ip

    nedg=0_ip 
    do ipoin=1,npoin0
       lmark(ipoin)=1_ip
       do iel=ptoel2(ipoin),ptoel2(ipoin+1)-1
          ielem=ptoel1(iel)
          do inode=1,nnode
             jpoin=elem(inode,ielem)
             if(lmark(jpoin)==0)then
                lmark(jpoin)=1_ip
                if(jpoin>ipoin)then
                   npoin=npoin+1
                   coor(1,npoin)=(coor(1,ipoin)+coor(1,jpoin))*c05
                   coor(2,npoin)=(coor(2,ipoin)+coor(2,jpoin))*c05
                   coor(3,npoin)=(coor(3,ipoin)+coor(3,jpoin))*c05
                   nedg=nedg+1
                   lremem(nedg)=jpoin   
                endif
             endif
          enddo
       enddo
       lcoun(ipoin+1)=nedg+1
    enddo
    !
    !     Reallocate the elements
    !
    nelemn=8*nelem
    call memrea(nelemn,memor_msh,'ELEM','divmsh',elem)
    !
    !     Create the new elements
    !
    do ielem=1,nelem

       ip1=elem(1,ielem)
       ip2=elem(2,ielem)
       ip3=elem(3,ielem)
       ip4=elem(4,ielem)
       !
       !     Get the new points
       !
       do iedg=1,6
          ipa=elem(tab(1,iedg),ielem) 
          ipb=elem(tab(2,iedg),ielem) 
          if(ipa<ipb)then
             ipmin=ipa
             ipmax=ipb
          else
             ipmin=ipb
             ipmax=ipa 
          endif
          !
          !     Loop on the points created by ipmin
          !
          do ipo=lcoun(ipmin),lcoun(ipmin+1)-1
             jpmax=lremem(ipo)
             if(jpmax==ipmax)exit
          enddo
          !
          !     Did we succeed?
          !
          if(ipo>lcoun(ipmin+1)-1)then
             write(*,*)'Error in divmsh, point not found'
             stop
          endif

          lpedg(iedg)=npoin0+ip 

       enddo
       !
       !     Now use the pattern
       ! 
       elem(1,ielem)= ip1
       elem(2,ielem)= lpedg(1)
       elem(3,ielem)= lpedg(2)
       elem(4,ielem)= lpedg(3)

       nelem=nelem+1
       elem(1,nelem)=lpedg(1)
       elem(2,nelem)=ip2
       elem(3,nelem)=lpedg(4)
       elem(4,nelem)=lpedg(5)

       nelem=nelem+1
       elem(1,nelem)=lpedg(2)
       elem(2,nelem)=lpedg(4)
       elem(3,nelem)=ip3
       elem(4,nelem)=lpedg(6)

       nelem=nelem+1
       elem(1,nelem)=lpedg(3)
       elem(2,nelem)=lpedg(5)
       elem(3,nelem)=lpedg(6)
       elem(4,nelem)=ip4

       nelem=nelem+1
       elem(1,nelem)=lpedg(1)
       elem(2,nelem)=lpedg(5)
       elem(3,nelem)=lpedg(3)
       elem(4,nelem)=lpedg(2)

       nelem=nelem+1
       elem(1,nelem)=lpedg(1)
       elem(2,nelem)=lpedg(4)
       elem(3,nelem)=lpedg(5)
       elem(4,nelem)=lpedg(2)

       nelem=nelem+1
       elem(1,nelem)=lpedg(6)
       elem(2,nelem)=lpedg(2)
       elem(3,nelem)=lpedg(4)
       elem(4,nelem)=lpedg(5)

       nelem=nelem+1
       elem(1,nelem)=lpedg(6)
       elem(2,nelem)=lpedg(3)
       elem(3,nelem)=lpedg(2)
       elem(4,nelem)=lpedg(5)

    enddo
    !
    !     Create the boundary faces
    !

    !
    !     Reallocate lface
    !
    nfacen=4*nface
    call memrea(nfacen,memor_msh,'LFACE','divmsh',lface)
    call memrea(nfacen,memor_msh,'LSURF','divmsh',lsurf)
    !
    !     First get the faces surrounding the points
    !
    call ptoelm(lface,nface,npoin,nnofa,ptoel1,ptoel2)
    !
    !     Create the new elements
    !
    do iface=1,nface

       isurf=lsurf(iface)
       ip1=lface(1,ielem)
       ip2=lface(2,ielem)
       ip3=lface(3,ielem)

       !
       !     Get the new points
       !
       do iedg=1,3
          ipa=lface(tab2(1,iedg),iface) 
          ipb=lface(tab2(2,iedg),iface) 
          if(ipa<ipb)then
             ipmin=ipa
             ipmax=ipb
          else
             ipmin=ipb
             ipmax=ipa 
          endif

          !
          !     Loop on the points created by ipmin
          !
          do ipo=lcoun(ipmin),lcoun(ipmin+1)-1
             jpmax=lremem(ipo)
             if(jpmax==ipmax)exit
          enddo
          !
          !     Did we succeed?
          !
          if(ipo>lcoun(ipmin+1)-1)then
             write(*,*)'Error in divmsh, point not found'
             stop
          endif

          lpedg(iedg)=npoin0+ip 

       enddo
       !
       !     Now use the pattern
       ! 
       lface(1,iface)= ip1
       lface(2,iface)= lpedg(1)
       lface(3,iface)= lpedg(3)

       nface=nface+1
       lface(1,nface)=lpedg(1)
       lface(2,nface)=ip2
       lface(3,nface)=lpedg(2)
       lsurf(nface)=isurf

       nface=nface+1
       lface(1,nface)=lpedg(3)
       lface(2,nface)=lpedg(2)
       lface(3,nface)=ip3
       lsurf(nface)=isurf

       nface=nface+1
       lface(1,nface)=lpedg(1)
       lface(2,nface)=lpedg(2)
       lface(3,nface)=lpedg(3)
       lsurf(nface)=isurf

    enddo

    call memchk(2_ip,istat,memor_msh,'LREMEM','divmsh',lremem)
    deallocate(lremem,stat=istat)
    if(istat/=0) call memerr(2_ip,'LREMEM','divmsh',0_ip)
    call memchk(2_ip,istat,memor_msh,'LCOUN','divmsh',lcoun)
    deallocate(lcoun,stat=istat)
    if(istat/=0) call memerr(2_ip,'LCOUN','divmsh',0_ip)
    call memchk(2_ip,istat,memor_msh,'LMARK','divmsh',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','divmsh',0_ip)

  end subroutine divmsh

  subroutine outfac(nnofa,nface,npoin,ndim,lface,coor)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)      :: nnofa,nface,npoin,ndim,lface(nnofa,nface)
    real(rp), intent(in)         :: coor(ndim,npoin)
    integer(ip)                  :: i,icont

    open(unit=50,file='outface.msh',status='unknown')
    rewind 50

    write(50,1)
    write(50,2)
    write(50,3)

    do i=1,npoin
       write(50,100)i,coor(1,i),coor(2,i),coor(3,i)
    enddo

    write(50,4)
    write(50,5)

    icont=0_ip
    do i=1, nface
       icont=icont+1
       write(50,300)icont,lface(1,i),lface(2,i),lface(3,i)
    enddo

    write(50,6)
    close(50)


1   format('MESH dimension 3 ElemType  Triangle Nnode 3')
2   format('Coordinates')
3   format('#node number   coor_x   coor_y  coor_z')
100 format(i10,3e20.10)
300 format(4i10)
4   format('end coordinates')
5   format('Elements')
6   format('end elements')


  end subroutine outfac

  subroutine outfac2(nnofa,nface,npoin,ndim,lface,coor,lfold,nfold,coold,npold)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)      :: nnofa,nface,npoin,ndim,npold,nfold
    integer(ip), intent(in)      :: lface(nnofa,nface),lfold(nnofa,nfold)
    real(rp), intent(in)         :: coor(ndim,npoin),coold(ndim,npold)
    integer(ip)                  :: i,icont


    open(unit=50,file='outface2.msh',status='unknown')
    rewind 50

    write(50,1)
    write(50,2)
    write(50,3)

    do i=1,npoin
       write(50,100)i,coor(1,i),coor(2,i),coor(3,i)
    enddo
    do i=1,npold
       write(50,100)npoin+i,coold(1,i),coold(2,i),coold(3,i)
    enddo

    write(50,4)
    write(50,5)

    icont=0_ip
    do i=1, nface
       icont=icont+1
       write(50,300)icont,lface(1,i),lface(2,i),lface(3,i),1_ip
    enddo
    do i=1,nfold
       icont=icont+1
       write(50,300)icont,lfold(1,i)+npoin,lfold(2,i)+npoin,lfold(3,i)+npoin,2_ip
    enddo

    write(50,6)
    close(50)

1   format('MESH dimension 3 ElemType  Triangle Nnode 3')
2   format('Coordinates')
3   format('#node number   coor_x   coor_y  coor_z')
100 format(i10,3e20.10)
300 format(5i10)
4   format('end coordinates')
5   format('Elements')
6   format('end elements')


  end subroutine outfac2

  subroutine findelem(iguess,elem,nelem,npoin,coor,pnew,ndim,ihost,nnode,eltoel,d1,d2,d3,d4)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip), intent(in)    :: ndim,nelem,npoin,nnode,iguess
    integer(ip),intent(in)     :: elem(nnode,nelem),eltoel(nnode,nelem)
    real(rp),intent(in)        :: coor(ndim,npoin),pnew(ndim)
    real(rp),intent(inout)     :: d1,d2,d3,d4
    integer(ip),intent(inout)  :: ihost
    integer(ip)                :: iter,maxiter,ip1,ip2,ip3,ip4,l
    real(rp)                   :: c40,p1(ndim),p2(ndim),p3(ndim)
    real(rp)                   :: dtot,rtol,rl
    !
    !     This subroutine finds the host element by a random walk through
    !     the mesh
    ! 
    c40=4.0d+00  
    maxiter=50_ip
    rtol=1.0d-06

    if(iguess==0)then
       ihost=0_ip 
       return
    endif
    if(iguess<0)then
       write(*,*)'Error in findelem, iguess=',iguess
       stop
    endif

    ihost=iguess

    do iter=1,maxiter 

       if(ihost==0)return 

       if(ihost<0 .or. ihost>nelem)then
          write(*,*)'Error findelem, ihost=',ihost
          stop
       endif

       call random_number(rl)
       l=floor(c40*rl)+1_ip


       ip1=elem(1,ihost)
       ip2=elem(2,ihost)
       ip3=elem(3,ihost)
       ip4=elem(4,ihost)
       !
       !     Compute volume for tolerances
       !
       p1(1)=coor(1,ip2)-coor(1,ip1)
       p1(2)=coor(2,ip2)-coor(2,ip1)
       p1(3)=coor(3,ip2)-coor(3,ip1)
       p2(1)=coor(1,ip3)-coor(1,ip1)
       p2(2)=coor(2,ip3)-coor(2,ip1)
       p2(3)=coor(3,ip3)-coor(3,ip1)
       p3(1)=coor(1,ip4)-coor(1,ip1)
       p3(2)=coor(2,ip4)-coor(2,ip1)
       p3(3)=coor(3,ip4)-coor(3,ip1)

       call orient3D(p1,p2,p3,dtot,ndim)
       !
       !     First face
       !
       if(l==1)then

          p1(1)=coor(1,ip2)-pnew(1)
          p1(2)=coor(2,ip2)-pnew(2)
          p1(3)=coor(3,ip2)-pnew(3)
          p2(1)=coor(1,ip3)-pnew(1)
          p2(2)=coor(2,ip3)-pnew(2)
          p2(3)=coor(3,ip3)-pnew(3)
          p3(1)=coor(1,ip4)-pnew(1)
          p3(2)=coor(2,ip4)-pnew(2)
          p3(3)=coor(3,ip4)-pnew(3)

          call orient3D(p1,p2,p3,d1,ndim)
          d1=d1/dtot 

          if(d1>rtol)then

             p1(1)=coor(1,ip1)-pnew(1)
             p1(2)=coor(2,ip1)-pnew(2)
             p1(3)=coor(3,ip1)-pnew(3)
             p2(1)=coor(1,ip4)-pnew(1)
             p2(2)=coor(2,ip4)-pnew(2)
             p2(3)=coor(3,ip4)-pnew(3)
             p3(1)=coor(1,ip3)-pnew(1)
             p3(2)=coor(2,ip3)-pnew(2)
             p3(3)=coor(3,ip3)-pnew(3)

             call orient3D(p1,p2,p3,d2,ndim)
             d2=d2/dtot 

             if(d2>rtol)then

                p1(1)=coor(1,ip1)-pnew(1)
                p1(2)=coor(2,ip1)-pnew(2)
                p1(3)=coor(3,ip1)-pnew(3)
                p2(1)=coor(1,ip2)-pnew(1)
                p2(2)=coor(2,ip2)-pnew(2)
                p2(3)=coor(3,ip2)-pnew(3)
                p3(1)=coor(1,ip4)-pnew(1)
                p3(2)=coor(2,ip4)-pnew(2)
                p3(3)=coor(3,ip4)-pnew(3)

                call orient3D(p1,p2,p3,d3,ndim)
                d3=d3/dtot 

                if(d3>rtol)then

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip2)-pnew(1)
                   p3(2)=coor(2,ip2)-pnew(2)
                   p3(3)=coor(3,ip2)-pnew(3)

                   call orient3D(p1,p2,p3,d4,ndim)
                   d4=d4/dtot 

                   if(d4>rtol)then
                      return
                   else if(d4<-rtol)then
                      ihost=eltoel(4,ihost) 
                   else
                      return
                   endif


                else if(d3<-rtol)then

                   ihost=eltoel(3,ihost)

                else

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip2)-pnew(1)
                   p3(2)=coor(2,ip2)-pnew(2)
                   p3(3)=coor(3,ip2)-pnew(3)

                   call orient3D(p1,p2,p3,d4,ndim)
                   d4=d4/dtot 

                   if(d4>rtol)then
                      return
                   else if(d4<-rtol)then
                      ihost=eltoel(4,ihost) 
                   else
                      return
                   endif

                endif

             else if(d2<-rtol)then

                ihost=eltoel(2,ihost)

             else

                p1(1)=coor(1,ip1)-pnew(1)
                p1(2)=coor(2,ip1)-pnew(2)
                p1(3)=coor(3,ip1)-pnew(3)
                p2(1)=coor(1,ip2)-pnew(1)
                p2(2)=coor(2,ip2)-pnew(2)
                p2(3)=coor(3,ip2)-pnew(3)
                p3(1)=coor(1,ip4)-pnew(1)
                p3(2)=coor(2,ip4)-pnew(2)
                p3(3)=coor(3,ip4)-pnew(3)

                call orient3D(p1,p2,p3,d3,ndim)
                d3=d3/dtot 

                if(d3>rtol)then

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip2)-pnew(1)
                   p3(2)=coor(2,ip2)-pnew(2)
                   p3(3)=coor(3,ip2)-pnew(3)

                   call orient3D(p1,p2,p3,d4,ndim)
                   d4=d4/dtot 

                   if(d4>rtol)then
                      return
                   else if(d4<-rtol)then
                      ihost=eltoel(4,ihost) 
                   else
                      return
                   endif

                else if(d3<-rtol)then

                   ihost=eltoel(3,ihost)

                else

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip2)-pnew(1)
                   p3(2)=coor(2,ip2)-pnew(2)
                   p3(3)=coor(3,ip2)-pnew(3)

                   call orient3D(p1,p2,p3,d4,ndim)
                   d4=d4/dtot 

                   if(d4>rtol)then
                      return
                   else if(d4<-rtol)then
                      ihost=eltoel(4,ihost) 
                   else
                      return
                   endif

                endif

             endif

          else if(d1<-rtol)then

             ihost=eltoel(1,ihost)

          else      

             p1(1)=coor(1,ip1)-pnew(1)
             p1(2)=coor(2,ip1)-pnew(2)
             p1(3)=coor(3,ip1)-pnew(3)
             p2(1)=coor(1,ip4)-pnew(1)
             p2(2)=coor(2,ip4)-pnew(2)
             p2(3)=coor(3,ip4)-pnew(3)
             p3(1)=coor(1,ip3)-pnew(1)
             p3(2)=coor(2,ip3)-pnew(2)
             p3(3)=coor(3,ip3)-pnew(3)

             call orient3D(p1,p2,p3,d2,ndim)
             d2=d2/dtot 

             if(d2>rtol)then

                p1(1)=coor(1,ip1)-pnew(1)
                p1(2)=coor(2,ip1)-pnew(2)
                p1(3)=coor(3,ip1)-pnew(3)
                p2(1)=coor(1,ip2)-pnew(1)
                p2(2)=coor(2,ip2)-pnew(2)
                p2(3)=coor(3,ip2)-pnew(3)
                p3(1)=coor(1,ip4)-pnew(1)
                p3(2)=coor(2,ip4)-pnew(2)
                p3(3)=coor(3,ip4)-pnew(3)

                call orient3D(p1,p2,p3,d3,ndim)
                d3=d3/dtot 

                if(d3>rtol)then

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip2)-pnew(1)
                   p3(2)=coor(2,ip2)-pnew(2)
                   p3(3)=coor(3,ip2)-pnew(3)

                   call orient3D(p1,p2,p3,d4,ndim)
                   d4=d4/dtot 

                   if(d4>rtol)then
                      return
                   else if(d4<-rtol)then
                      ihost=eltoel(4,ihost) 
                   else
                      return
                   endif

                else if(d3<-rtol)then

                   ihost=eltoel(3,ihost)

                else

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip2)-pnew(1)
                   p3(2)=coor(2,ip2)-pnew(2)
                   p3(3)=coor(3,ip2)-pnew(3)

                   call orient3D(p1,p2,p3,d4,ndim)
                   d4=d4/dtot 

                   if(d4>rtol)then
                      return
                   else if(d4<-rtol)then
                      ihost=eltoel(4,ihost) 
                   else
                      return
                   endif

                endif

             else if(d2<-rtol)then

                ihost=eltoel(2,ihost)

             else

                p1(1)=coor(1,ip1)-pnew(1)
                p1(2)=coor(2,ip1)-pnew(2)
                p1(3)=coor(3,ip1)-pnew(3)
                p2(1)=coor(1,ip2)-pnew(1)
                p2(2)=coor(2,ip2)-pnew(2)
                p2(3)=coor(3,ip2)-pnew(3)
                p3(1)=coor(1,ip4)-pnew(1)
                p3(2)=coor(2,ip4)-pnew(2)
                p3(3)=coor(3,ip4)-pnew(3)

                call orient3D(p1,p2,p3,d3,ndim)
                d3=d3/dtot 

                if(d3>rtol)then

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip2)-pnew(1)
                   p3(2)=coor(2,ip2)-pnew(2)
                   p3(3)=coor(3,ip2)-pnew(3)

                   call orient3D(p1,p2,p3,d4,ndim)
                   d4=d4/dtot 

                   if(d4>rtol)then
                      return
                   else if(d4<-rtol)then
                      ihost=eltoel(4,ihost) 
                   else
                      return
                   endif

                else if(d3<-rtol)then

                   ihost=eltoel(3,ihost)

                else

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip2)-pnew(1)
                   p3(2)=coor(2,ip2)-pnew(2)
                   p3(3)=coor(3,ip2)-pnew(3)

                   call orient3D(p1,p2,p3,d4,ndim)
                   d4=d4/dtot 

                   if(d4>rtol)then
                      return
                   else if(d4<-rtol)then
                      ihost=eltoel(4,ihost) 
                   else

                      write(*,*)'Error in findelem, degenerate case'
                      stop
                   endif

                endif

             endif

          endif
          !
          !     Second face
          !
       else if(l==2)then

          p1(1)=coor(1,ip1)-pnew(1)
          p1(2)=coor(2,ip1)-pnew(2)
          p1(3)=coor(3,ip1)-pnew(3)
          p2(1)=coor(1,ip4)-pnew(1)
          p2(2)=coor(2,ip4)-pnew(2)
          p2(3)=coor(3,ip4)-pnew(3)
          p3(1)=coor(1,ip3)-pnew(1)
          p3(2)=coor(2,ip3)-pnew(2)
          p3(3)=coor(3,ip3)-pnew(3)

          call orient3D(p1,p2,p3,d2,ndim)
          d2=d2/dtot 

          if(d2>rtol)then

             p1(1)=coor(1,ip1)-pnew(1)
             p1(2)=coor(2,ip1)-pnew(2)
             p1(3)=coor(3,ip1)-pnew(3)
             p2(1)=coor(1,ip2)-pnew(1)
             p2(2)=coor(2,ip2)-pnew(2)
             p2(3)=coor(3,ip2)-pnew(3)
             p3(1)=coor(1,ip4)-pnew(1)
             p3(2)=coor(2,ip4)-pnew(2)
             p3(3)=coor(3,ip4)-pnew(3)

             call orient3D(p1,p2,p3,d3,ndim)
             d3=d3/dtot 

             if(d3>rtol)then

                p1(1)=coor(1,ip1)-pnew(1)
                p1(2)=coor(2,ip1)-pnew(2)
                p1(3)=coor(3,ip1)-pnew(3)
                p2(1)=coor(1,ip3)-pnew(1)
                p2(2)=coor(2,ip3)-pnew(2)
                p2(3)=coor(3,ip3)-pnew(3)
                p3(1)=coor(1,ip2)-pnew(1)
                p3(2)=coor(2,ip2)-pnew(2)
                p3(3)=coor(3,ip2)-pnew(3)

                call orient3D(p1,p2,p3,d4,ndim)
                d4=d4/dtot 

                if(d4>rtol)then

                   p1(1)=coor(1,ip2)-pnew(1)
                   p1(2)=coor(2,ip2)-pnew(2)
                   p1(3)=coor(3,ip2)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d1,ndim)
                   d1=d1/dtot 

                   if(d1>rtol)then
                      return
                   else if(d1<-rtol)then
                      ihost=eltoel(1,ihost) 
                   else
                      return
                   endif


                else if(d4<-rtol)then

                   ihost=eltoel(4,ihost)

                else

                   p1(1)=coor(1,ip2)-pnew(1)
                   p1(2)=coor(2,ip2)-pnew(2)
                   p1(3)=coor(3,ip2)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d1,ndim)
                   d1=d1/dtot 

                   if(d1>rtol)then
                      return
                   else if(d1<-rtol)then
                      ihost=eltoel(1,ihost) 
                   else
                      return
                   endif

                endif

             else if(d3<-rtol)then

                ihost=eltoel(3,ihost)

             else

                p1(1)=coor(1,ip1)-pnew(1)
                p1(2)=coor(2,ip1)-pnew(2)
                p1(3)=coor(3,ip1)-pnew(3)
                p2(1)=coor(1,ip3)-pnew(1)
                p2(2)=coor(2,ip3)-pnew(2)
                p2(3)=coor(3,ip3)-pnew(3)
                p3(1)=coor(1,ip2)-pnew(1)
                p3(2)=coor(2,ip2)-pnew(2)
                p3(3)=coor(3,ip2)-pnew(3)

                call orient3D(p1,p2,p3,d4,ndim)
                d4=d4/dtot 

                if(d4>rtol)then

                   p1(1)=coor(1,ip2)-pnew(1)
                   p1(2)=coor(2,ip2)-pnew(2)
                   p1(3)=coor(3,ip2)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d1,ndim)
                   d1=d1/dtot 

                   if(d1>rtol)then
                      return
                   else if(d1<-rtol)then
                      ihost=eltoel(1,ihost) 
                   else
                      return
                   endif

                else if(d4<-rtol)then

                   ihost=eltoel(4,ihost)

                else

                   p1(1)=coor(1,ip2)-pnew(1)
                   p1(2)=coor(2,ip2)-pnew(2)
                   p1(3)=coor(3,ip2)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d1,ndim)
                   d1=d1/dtot 

                   if(d1>rtol)then
                      return
                   else if(d1<-rtol)then
                      ihost=eltoel(1,ihost) 
                   else
                      return
                   endif

                endif

             endif

          else if(d2<-rtol)then

             ihost=eltoel(2,ihost)

          else      

             p1(1)=coor(1,ip1)-pnew(1)
             p1(2)=coor(2,ip1)-pnew(2)
             p1(3)=coor(3,ip1)-pnew(3)
             p2(1)=coor(1,ip2)-pnew(1)
             p2(2)=coor(2,ip2)-pnew(2)
             p2(3)=coor(3,ip2)-pnew(3)
             p3(1)=coor(1,ip4)-pnew(1)
             p3(2)=coor(2,ip4)-pnew(2)
             p3(3)=coor(3,ip4)-pnew(3)

             call orient3D(p1,p2,p3,d3,ndim)
             d3=d3/dtot 

             if(d3>rtol)then

                p1(1)=coor(1,ip1)-pnew(1)
                p1(2)=coor(2,ip1)-pnew(2)
                p1(3)=coor(3,ip1)-pnew(3)
                p2(1)=coor(1,ip3)-pnew(1)
                p2(2)=coor(2,ip3)-pnew(2)
                p2(3)=coor(3,ip3)-pnew(3)
                p3(1)=coor(1,ip2)-pnew(1)
                p3(2)=coor(2,ip2)-pnew(2)
                p3(3)=coor(3,ip2)-pnew(3)

                call orient3D(p1,p2,p3,d4,ndim)
                d4=d4/dtot 

                if(d4>rtol)then

                   p1(1)=coor(1,ip2)-pnew(1)
                   p1(2)=coor(2,ip2)-pnew(2)
                   p1(3)=coor(3,ip2)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d1,ndim)
                   d1=d1/dtot 

                   if(d1>rtol)then
                      return
                   else if(d1<-rtol)then
                      ihost=eltoel(1,ihost) 
                   else
                      return
                   endif


                else if(d4<-rtol)then

                   ihost=eltoel(4,ihost)

                else

                   p1(1)=coor(1,ip2)-pnew(1)
                   p1(2)=coor(2,ip2)-pnew(2)
                   p1(3)=coor(3,ip2)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d1,ndim)
                   d1=d1/dtot 

                   if(d1>rtol)then
                      return
                   else if(d1<-rtol)then
                      ihost=eltoel(1,ihost) 
                   else
                      return
                   endif

                endif

             else if(d3<-rtol)then

                ihost=eltoel(3,ihost)

             else

                p1(1)=coor(1,ip1)-pnew(1)
                p1(2)=coor(2,ip1)-pnew(2)
                p1(3)=coor(3,ip1)-pnew(3)
                p2(1)=coor(1,ip3)-pnew(1)
                p2(2)=coor(2,ip3)-pnew(2)
                p2(3)=coor(3,ip3)-pnew(3)
                p3(1)=coor(1,ip2)-pnew(1)
                p3(2)=coor(2,ip2)-pnew(2)
                p3(3)=coor(3,ip2)-pnew(3)

                call orient3D(p1,p2,p3,d4,ndim)
                d4=d4/dtot 

                if(d4>rtol)then

                   p1(1)=coor(1,ip2)-pnew(1)
                   p1(2)=coor(2,ip2)-pnew(2)
                   p1(3)=coor(3,ip2)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d1,ndim)
                   d1=d1/dtot 

                   if(d1>rtol)then
                      return
                   else if(d1<-rtol)then
                      ihost=eltoel(1,ihost) 
                   else
                      return
                   endif

                else if(d4<-rtol)then

                   ihost=eltoel(4,ihost)

                else

                   p1(1)=coor(1,ip2)-pnew(1)
                   p1(2)=coor(2,ip2)-pnew(2)
                   p1(3)=coor(3,ip2)-pnew(3)
                   p2(1)=coor(1,ip3)-pnew(1)
                   p2(2)=coor(2,ip3)-pnew(2)
                   p2(3)=coor(3,ip3)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d1,ndim)
                   d1=d1/dtot 

                   if(d1>rtol)then
                      return
                   else if(d1<-rtol)then
                      ihost=eltoel(4,ihost) 
                   else

                      write(*,*)'Error in findelem, degenerate case'
                      stop
                   endif

                endif

             endif

          endif
          !
          !     Third face
          !
       else if(l==3)then

          p1(1)=coor(1,ip1)-pnew(1)
          p1(2)=coor(2,ip1)-pnew(2)
          p1(3)=coor(3,ip1)-pnew(3)
          p2(1)=coor(1,ip2)-pnew(1)
          p2(2)=coor(2,ip2)-pnew(2)
          p2(3)=coor(3,ip2)-pnew(3)
          p3(1)=coor(1,ip4)-pnew(1)
          p3(2)=coor(2,ip4)-pnew(2)
          p3(3)=coor(3,ip4)-pnew(3)

          call orient3D(p1,p2,p3,d3,ndim)
          d3=d3/dtot 

          if(d3>rtol)then

             p1(1)=coor(1,ip1)-pnew(1)
             p1(2)=coor(2,ip1)-pnew(2)
             p1(3)=coor(3,ip1)-pnew(3)
             p2(1)=coor(1,ip3)-pnew(1)
             p2(2)=coor(2,ip3)-pnew(2)
             p2(3)=coor(3,ip3)-pnew(3)
             p3(1)=coor(1,ip2)-pnew(1)
             p3(2)=coor(2,ip2)-pnew(2)
             p3(3)=coor(3,ip2)-pnew(3)

             call orient3D(p1,p2,p3,d4,ndim)
             d4=d4/dtot 

             if(d4>rtol)then

                p1(1)=coor(1,ip2)-pnew(1)
                p1(2)=coor(2,ip2)-pnew(2)
                p1(3)=coor(3,ip2)-pnew(3)
                p2(1)=coor(1,ip3)-pnew(1)
                p2(2)=coor(2,ip3)-pnew(2)
                p2(3)=coor(3,ip3)-pnew(3)
                p3(1)=coor(1,ip4)-pnew(1)
                p3(2)=coor(2,ip4)-pnew(2)
                p3(3)=coor(3,ip4)-pnew(3)

                call orient3D(p1,p2,p3,d1,ndim)
                d1=d1/dtot 

                if(d1>rtol)then

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip4)-pnew(1)
                   p2(2)=coor(2,ip4)-pnew(2)
                   p2(3)=coor(3,ip4)-pnew(3)
                   p3(1)=coor(1,ip3)-pnew(1)
                   p3(2)=coor(2,ip3)-pnew(2)
                   p3(3)=coor(3,ip3)-pnew(3)

                   call orient3D(p1,p2,p3,d2,ndim)
                   d2=d2/dtot 

                   if(d2>rtol)then
                      return
                   else if(d2<-rtol)then
                      ihost=eltoel(2,ihost) 
                   else
                      return
                   endif


                else if(d1<-rtol)then

                   ihost=eltoel(1,ihost)

                else

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip4)-pnew(1)
                   p2(2)=coor(2,ip4)-pnew(2)
                   p2(3)=coor(3,ip4)-pnew(3)
                   p3(1)=coor(1,ip3)-pnew(1)
                   p3(2)=coor(2,ip3)-pnew(2)
                   p3(3)=coor(3,ip3)-pnew(3)

                   call orient3D(p1,p2,p3,d2,ndim)
                   d2=d2/dtot 

                   if(d2>rtol)then
                      return
                   else if(d2<-rtol)then
                      ihost=eltoel(2,ihost) 
                   else
                      return
                   endif

                endif

             else if(d4<-rtol)then

                ihost=eltoel(4,ihost)

             else

                p1(1)=coor(1,ip2)-pnew(1)
                p1(2)=coor(2,ip2)-pnew(2)
                p1(3)=coor(3,ip2)-pnew(3)
                p2(1)=coor(1,ip3)-pnew(1)
                p2(2)=coor(2,ip3)-pnew(2)
                p2(3)=coor(3,ip3)-pnew(3)
                p3(1)=coor(1,ip4)-pnew(1)
                p3(2)=coor(2,ip4)-pnew(2)
                p3(3)=coor(3,ip4)-pnew(3)

                call orient3D(p1,p2,p3,d1,ndim)
                d1=d1/dtot 

                if(d1>rtol)then

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip4)-pnew(1)
                   p2(2)=coor(2,ip4)-pnew(2)
                   p2(3)=coor(3,ip4)-pnew(3)
                   p3(1)=coor(1,ip3)-pnew(1)
                   p3(2)=coor(2,ip3)-pnew(2)
                   p3(3)=coor(3,ip3)-pnew(3)

                   call orient3D(p1,p2,p3,d2,ndim)
                   d2=d2/dtot 

                   if(d2>rtol)then
                      return
                   else if(d2<-rtol)then
                      ihost=eltoel(2,ihost) 
                   else
                      return
                   endif

                else if(d1<-rtol)then

                   ihost=eltoel(1,ihost)

                else

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip4)-pnew(1)
                   p2(2)=coor(2,ip4)-pnew(2)
                   p2(3)=coor(3,ip4)-pnew(3)
                   p3(1)=coor(1,ip3)-pnew(1)
                   p3(2)=coor(2,ip3)-pnew(2)
                   p3(3)=coor(3,ip3)-pnew(3)

                   call orient3D(p1,p2,p3,d2,ndim)
                   d2=d2/dtot 

                   if(d2>rtol)then
                      return
                   else if(d2<-rtol)then
                      ihost=eltoel(2,ihost) 
                   else
                      return
                   endif

                endif

             endif

          else if(d3<-rtol)then

             ihost=eltoel(3,ihost)

          else      

             p1(1)=coor(1,ip1)-pnew(1)
             p1(2)=coor(2,ip1)-pnew(2)
             p1(3)=coor(3,ip1)-pnew(3)
             p2(1)=coor(1,ip3)-pnew(1)
             p2(2)=coor(2,ip3)-pnew(2)
             p2(3)=coor(3,ip3)-pnew(3)
             p3(1)=coor(1,ip2)-pnew(1)
             p3(2)=coor(2,ip2)-pnew(2)
             p3(3)=coor(3,ip2)-pnew(3)

             call orient3D(p1,p2,p3,d4,ndim)
             d4=d4/dtot 

             if(d4>rtol)then

                p1(1)=coor(1,ip2)-pnew(1)
                p1(2)=coor(2,ip2)-pnew(2)
                p1(3)=coor(3,ip2)-pnew(3)
                p2(1)=coor(1,ip3)-pnew(1)
                p2(2)=coor(2,ip3)-pnew(2)
                p2(3)=coor(3,ip3)-pnew(3)
                p3(1)=coor(1,ip4)-pnew(1)
                p3(2)=coor(2,ip4)-pnew(2)
                p3(3)=coor(3,ip4)-pnew(3)

                call orient3D(p1,p2,p3,d1,ndim)
                d1=d1/dtot 

                if(d1>rtol)then

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip4)-pnew(1)
                   p2(2)=coor(2,ip4)-pnew(2)
                   p2(3)=coor(3,ip4)-pnew(3)
                   p3(1)=coor(1,ip3)-pnew(1)
                   p3(2)=coor(2,ip3)-pnew(2)
                   p3(3)=coor(3,ip3)-pnew(3)

                   call orient3D(p1,p2,p3,d2,ndim)
                   d2=d2/dtot 

                   if(d2>rtol)then
                      return
                   else if(d2<-rtol)then
                      ihost=eltoel(2,ihost) 
                   else
                      return
                   endif


                else if(d1<-rtol)then

                   ihost=eltoel(1,ihost)

                else

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip4)-pnew(1)
                   p2(2)=coor(2,ip4)-pnew(2)
                   p2(3)=coor(3,ip4)-pnew(3)
                   p3(1)=coor(1,ip3)-pnew(1)
                   p3(2)=coor(2,ip3)-pnew(2)
                   p3(3)=coor(3,ip3)-pnew(3)

                   call orient3D(p1,p2,p3,d2,ndim)
                   d2=d2/dtot 

                   if(d2>rtol)then
                      return
                   else if(d2<-rtol)then
                      ihost=eltoel(2,ihost) 
                   else
                      return
                   endif

                endif

             else if(d4<-rtol)then

                ihost=eltoel(4,ihost)

             else

                p1(1)=coor(1,ip2)-pnew(1)
                p1(2)=coor(2,ip2)-pnew(2)
                p1(3)=coor(3,ip2)-pnew(3)
                p2(1)=coor(1,ip3)-pnew(1)
                p2(2)=coor(2,ip3)-pnew(2)
                p2(3)=coor(3,ip3)-pnew(3)
                p3(1)=coor(1,ip4)-pnew(1)
                p3(2)=coor(2,ip4)-pnew(2)
                p3(3)=coor(3,ip4)-pnew(3)

                call orient3D(p1,p2,p3,d1,ndim)
                d1=d1/dtot 

                if(d1>rtol)then

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip4)-pnew(1)
                   p2(2)=coor(2,ip4)-pnew(2)
                   p2(3)=coor(3,ip4)-pnew(3)
                   p3(1)=coor(1,ip3)-pnew(1)
                   p3(2)=coor(2,ip3)-pnew(2)
                   p3(3)=coor(3,ip3)-pnew(3)

                   call orient3D(p1,p2,p3,d2,ndim)
                   d2=d2/dtot 

                   if(d2>rtol)then
                      return
                   else if(d2<-rtol)then
                      ihost=eltoel(2,ihost) 
                   else
                      return
                   endif

                else if(d1<-rtol)then

                   ihost=eltoel(1,ihost)

                else

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip4)-pnew(1)
                   p2(2)=coor(2,ip4)-pnew(2)
                   p2(3)=coor(3,ip4)-pnew(3)
                   p3(1)=coor(1,ip3)-pnew(1)
                   p3(2)=coor(2,ip3)-pnew(2)
                   p3(3)=coor(3,ip3)-pnew(3)

                   call orient3D(p1,p2,p3,d2,ndim)
                   d2=d2/dtot 

                   if(d2>rtol)then
                      return
                   else if(d2<-rtol)then
                      ihost=eltoel(2,ihost) 
                   else

                      write(*,*)'Error in findelem, degenerate case'
                      stop
                   endif

                endif

             endif

          endif
          !
          !     Fourth face
          !
       else

          p1(1)=coor(1,ip1)-pnew(1)
          p1(2)=coor(2,ip1)-pnew(2)
          p1(3)=coor(3,ip1)-pnew(3)
          p2(1)=coor(1,ip3)-pnew(1)
          p2(2)=coor(2,ip3)-pnew(2)
          p2(3)=coor(3,ip3)-pnew(3)
          p3(1)=coor(1,ip2)-pnew(1)
          p3(2)=coor(2,ip2)-pnew(2)
          p3(3)=coor(3,ip2)-pnew(3)

          call orient3D(p1,p2,p3,d4,ndim)
          d4=d4/dtot 

          if(d4>rtol)then

             p1(1)=coor(1,ip2)-pnew(1)
             p1(2)=coor(2,ip2)-pnew(2)
             p1(3)=coor(3,ip2)-pnew(3)
             p2(1)=coor(1,ip3)-pnew(1)
             p2(2)=coor(2,ip3)-pnew(2)
             p2(3)=coor(3,ip3)-pnew(3)
             p3(1)=coor(1,ip4)-pnew(1)
             p3(2)=coor(2,ip4)-pnew(2)
             p3(3)=coor(3,ip4)-pnew(3)

             call orient3D(p1,p2,p3,d1,ndim)
             d1=d1/dtot 

             if(d1>rtol)then

                p1(1)=coor(1,ip1)-pnew(1)
                p1(2)=coor(2,ip1)-pnew(2)
                p1(3)=coor(3,ip1)-pnew(3)
                p2(1)=coor(1,ip4)-pnew(1)
                p2(2)=coor(2,ip4)-pnew(2)
                p2(3)=coor(3,ip4)-pnew(3)
                p3(1)=coor(1,ip3)-pnew(1)
                p3(2)=coor(2,ip3)-pnew(2)
                p3(3)=coor(3,ip3)-pnew(3)

                call orient3D(p1,p2,p3,d2,ndim)
                d2=d2/dtot 

                if(d2>rtol)then

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip2)-pnew(1)
                   p2(2)=coor(2,ip2)-pnew(2)
                   p2(3)=coor(3,ip2)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d3,ndim)
                   d3=d3/dtot 

                   if(d3>rtol)then
                      return
                   else if(d3<-rtol)then
                      ihost=eltoel(3,ihost) 
                   else
                      return
                   endif


                else if(d2<-rtol)then

                   ihost=eltoel(2,ihost)

                else

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip2)-pnew(1)
                   p2(2)=coor(2,ip2)-pnew(2)
                   p2(3)=coor(3,ip2)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d3,ndim)
                   d3=d3/dtot 

                   if(d3>rtol)then
                      return
                   else if(d3<-rtol)then
                      ihost=eltoel(3,ihost) 
                   else
                      return
                   endif

                endif

             else if(d1<-rtol)then

                ihost=eltoel(1,ihost)

             else

                p1(1)=coor(1,ip1)-pnew(1)
                p1(2)=coor(2,ip1)-pnew(2)
                p1(3)=coor(3,ip1)-pnew(3)
                p2(1)=coor(1,ip4)-pnew(1)
                p2(2)=coor(2,ip4)-pnew(2)
                p2(3)=coor(3,ip4)-pnew(3)
                p3(1)=coor(1,ip3)-pnew(1)
                p3(2)=coor(2,ip3)-pnew(2)
                p3(3)=coor(3,ip3)-pnew(3)

                call orient3D(p1,p2,p3,d2,ndim)
                d2=d2/dtot 

                if(d2>rtol)then

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip2)-pnew(1)
                   p2(2)=coor(2,ip2)-pnew(2)
                   p2(3)=coor(3,ip2)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d3,ndim)
                   d3=d3/dtot 

                   if(d3>rtol)then
                      return
                   else if(d3<-rtol)then
                      ihost=eltoel(3,ihost) 
                   else
                      return
                   endif

                else if(d2<-rtol)then

                   ihost=eltoel(2,ihost)

                else

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip2)-pnew(1)
                   p2(2)=coor(2,ip2)-pnew(2)
                   p2(3)=coor(3,ip2)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d3,ndim)
                   d3=d3/dtot 

                   if(d3>rtol)then
                      return
                   else if(d3<-rtol)then
                      ihost=eltoel(3,ihost) 
                   else
                      return
                   endif

                endif

             endif

          else if(d4<-rtol)then

             ihost=eltoel(4,ihost)

          else      

             p1(1)=coor(1,ip2)-pnew(1)
             p1(2)=coor(2,ip2)-pnew(2)
             p1(3)=coor(3,ip2)-pnew(3)
             p2(1)=coor(1,ip3)-pnew(1)
             p2(2)=coor(2,ip3)-pnew(2)
             p2(3)=coor(3,ip3)-pnew(3)
             p3(1)=coor(1,ip4)-pnew(1)
             p3(2)=coor(2,ip4)-pnew(2)
             p3(3)=coor(3,ip4)-pnew(3)

             call orient3D(p1,p2,p3,d1,ndim)
             d1=d1/dtot 

             if(d1>rtol)then

                p1(1)=coor(1,ip1)-pnew(1)
                p1(2)=coor(2,ip1)-pnew(2)
                p1(3)=coor(3,ip1)-pnew(3)
                p2(1)=coor(1,ip4)-pnew(1)
                p2(2)=coor(2,ip4)-pnew(2)
                p2(3)=coor(3,ip4)-pnew(3)
                p3(1)=coor(1,ip3)-pnew(1)
                p3(2)=coor(2,ip3)-pnew(2)
                p3(3)=coor(3,ip3)-pnew(3)

                call orient3D(p1,p2,p3,d2,ndim)
                d2=d2/dtot 

                if(d2>rtol)then

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip2)-pnew(1)
                   p2(2)=coor(2,ip2)-pnew(2)
                   p2(3)=coor(3,ip2)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d3,ndim)
                   d3=d3/dtot 

                   if(d3>rtol)then
                      return
                   else if(d3<-rtol)then
                      ihost=eltoel(3,ihost) 
                   else
                      return
                   endif


                else if(d2<-rtol)then

                   ihost=eltoel(2,ihost)

                else

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip2)-pnew(1)
                   p2(2)=coor(2,ip2)-pnew(2)
                   p2(3)=coor(3,ip2)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d3,ndim)
                   d3=d3/dtot 

                   if(d3>rtol)then
                      return
                   else if(d3<-rtol)then
                      ihost=eltoel(3,ihost) 
                   else
                      return
                   endif

                endif

             else if(d1<-rtol)then

                ihost=eltoel(1,ihost)

             else

                p1(1)=coor(1,ip1)-pnew(1)
                p1(2)=coor(2,ip1)-pnew(2)
                p1(3)=coor(3,ip1)-pnew(3)
                p2(1)=coor(1,ip4)-pnew(1)
                p2(2)=coor(2,ip4)-pnew(2)
                p2(3)=coor(3,ip4)-pnew(3)
                p3(1)=coor(1,ip3)-pnew(1)
                p3(2)=coor(2,ip3)-pnew(2)
                p3(3)=coor(3,ip3)-pnew(3)

                call orient3D(p1,p2,p3,d2,ndim)
                d2=d2/dtot 

                if(d2>rtol)then

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip2)-pnew(1)
                   p2(2)=coor(2,ip2)-pnew(2)
                   p2(3)=coor(3,ip2)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d3,ndim)
                   d3=d3/dtot 

                   if(d3>rtol)then
                      return
                   else if(d3<-rtol)then
                      ihost=eltoel(3,ihost) 
                   else
                      return
                   endif

                else if(d2<-rtol)then

                   ihost=eltoel(2,ihost)

                else

                   p1(1)=coor(1,ip1)-pnew(1)
                   p1(2)=coor(2,ip1)-pnew(2)
                   p1(3)=coor(3,ip1)-pnew(3)
                   p2(1)=coor(1,ip2)-pnew(1)
                   p2(2)=coor(2,ip2)-pnew(2)
                   p2(3)=coor(3,ip2)-pnew(3)
                   p3(1)=coor(1,ip4)-pnew(1)
                   p3(2)=coor(2,ip4)-pnew(2)
                   p3(3)=coor(3,ip4)-pnew(3)

                   call orient3D(p1,p2,p3,d3,ndim)
                   d3=d3/dtot 

                   if(d3>rtol)then
                      return
                   else if(d3<-rtol)then
                      ihost=eltoel(3,ihost) 
                   else

                      write(*,*)'Error in findelem, degenerate case'
                      stop
                   endif

                endif

             endif

          endif

       endif

    enddo
    !
    !     Host element not found
    !
    ihost=0_ip

  end subroutine findelem

  subroutine length(ip1,ip2,rsize,coor,rlen)
    use def_kintyp, only       :  ip,rp,lg
    implicit none
    integer(ip),intent(in)          :: ip1,ip2
    real(rp),intent(in)             :: rsize(*),coor(3,*)
    real(rp),intent(inout)          :: rlen  
    real(rp)                        :: rx,ry,rz,rlen1,c05

    rx=coor(1,ip2)-coor(1,ip1)
    ry=coor(2,ip2)-coor(2,ip1)
    rz=coor(3,ip2)-coor(3,ip1)
    c05=0.5d+00

    rlen1=sqrt(rx*rx+ry*ry+rz*rz)
    rlen=rlen1*((1.0_rp/rsize(ip1))+(1.0_rp/rsize(ip2)))*c05

  end subroutine length

  subroutine ptoedg(lface,nface,npoin,nnofa,ptoel1,ptoel2,&
       nedge,ledge)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only          : memor_msh
    implicit none
    integer(ip), intent(in)            :: nface,npoin,nnofa
    integer(ip), intent(in)            :: ptoel1(*),ptoel2(npoin+1)
    integer(ip), intent(in)            :: lface(nnofa,nface)
    integer(ip), intent(inout)         :: nedge
    integer(ip),pointer                :: ledge(:,:)
    integer(ip)                        :: iface,j,ipoin,ip1,iel
    integer(ip), pointer               :: lmark(:)
    integer(4)                 :: istat
    !
    !     This routine gets the edges from elements lface
    !
    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','ptoedg',lmark)
    !
    !     First count the edges
    !
    !
    !     Loop on the points 
    !
    nedge=0_ip
    do ipoin=1,npoin
       !
       !     Loop on the elements surrounding the point
       !
       do iel=ptoel2(ipoin),ptoel2(ipoin+1)-1
          iface=ptoel1(iel)
          do j=1,nnofa
             ip1=lface(j,iface)
             if(ipoin<ip1)then
                if(lmark(ip1)/=ipoin)then
                   lmark(ip1)=ipoin 
                   nedge=nedge+1  
                endif
             endif
          enddo
       enddo
    enddo
    !
    !     Allocate the edges
    !
    if(.not.associated(ledge))then
       allocate(ledge(2,nedge),stat=istat)
       call memchk(zero,istat,memor_msh,'LEDGE','ptoedg',ledge)
    else
       call memrea(nedge,memor_msh,'LEDGE','ptoedg',ledge)
    endif
    !
    !     Reset lmark
    !
    do ipoin=1,npoin
       lmark(ipoin)=0_ip
    enddo
    !
    !     Fill the edges
    !
    !
    !     Loop on the points 
    !
    nedge=0_ip
    do ipoin=1,npoin
       !
       !     Loop on the elements surrounding the point
       !
       do iel=ptoel2(ipoin),ptoel2(ipoin+1)-1
          iface=ptoel1(iel)
          do j=1,nnofa
             ip1=lface(j,iface)
             if(ipoin<ip1)then
                if(lmark(ip1)/=ipoin)then
                   lmark(ip1)=ipoin  
                   nedge=nedge+1  
                   ledge(1,nedge)=ipoin
                   ledge(2,nedge)=ip1
                endif
             endif
          enddo
       enddo
    enddo

    call memchk(2_ip,istat,memor_msh,'LMARK','ptoedg',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','ptoedg',0_ip)

  end subroutine ptoedg
 
  subroutine ptoedgtet(lface,nface,npoin,nnofa,ptoel1,ptoel2,&
       nedge,ledge,ledglm)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only          : memor_msh
    implicit none
    integer(ip), intent(in)            :: nface,npoin,nnofa
    integer(ip), intent(in)            :: ptoel1(*),ptoel2(npoin+1)
    integer(ip), intent(in)            :: lface(nnofa,nface)
    integer(ip), intent(inout)         :: nedge
    integer(ip),pointer                :: ledge(:,:),ledglm(:,:)
    integer(ip)                        :: iface,j,ipoin,ip1,iel,i,iedg,ipedg,iedge
    integer(ip), pointer               :: lmark(:)
    integer(4)                 :: istat
    integer(ip)               :: ltab(4,4)=RESHAPE((/0,1,2,3,&
         1,0,4,5,&
         2,4,0,6,&
         3,5,6,0/),(/4,4/)) 
    !
    !     This routine gets the edges from elements lface
    !
    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','ptoedgtet',lmark)
    if(.not.associated(ledglm))then
       allocate(ledglm(nnofa,nface),stat=istat)
       call memchk(zero,istat,memor_msh,'LEDGLM','ptoedgtet',ledglm)
    else
       call memrea(nface,memor_msh,'LEDGLM','ptoedgtet',ledglm)
    endif
    !
    !     First count the edges
    !
    !
    !     Loop on the points 
    !
    nedge=0_ip
    do ipoin=1,npoin
       !
       !     Loop on the elements surrounding the point
       !
       do iel=ptoel2(ipoin),ptoel2(ipoin+1)-1
          iface=ptoel1(iel)
          do j=1,4
             ip1=lface(j,iface)
             if(ipoin<ip1)then
                if(lmark(ip1)/=ipoin)then
                   lmark(ip1)=ipoin 
                   nedge=nedge+1  
                endif
             endif
          enddo
       enddo
    enddo
    !
    !     Allocate the edges
    !
    if(.not.associated(ledge))then
       allocate(ledge(2,nedge),stat=istat)
       call memchk(zero,istat,memor_msh,'LEDGE','ptoedgtet',ledge)
    else
       call memrea(nedge,memor_msh,'LEDGE','ptoedgtet',ledge)
    endif
    !
    !     Reset lmark
    !
    do ipoin=1,npoin
       lmark(ipoin)=0_ip
    enddo
    !
    !     Fill the edges
    !
    !
    !     Loop on the points 
    !
    nedge=0_ip
    ipedg=1_ip
    do ipoin=1,npoin
       !
       !     Loop on the elements surrounding the point
       !
       do iel=ptoel2(ipoin),ptoel2(ipoin+1)-1
          iface=ptoel1(iel)
          do j=1,4
             ip1=lface(j,iface)
             if(ipoin<ip1)then
                !
                !     Is this a new edge?
                !
                if(lmark(ip1)/=ipoin)then
                   lmark(ip1)=ipoin  
                   nedge=nedge+1  
                   ledge(1,nedge)=ipoin
                   ledge(2,nedge)=ip1
                   !
                   !     Find the position of ipoin
                   !
                   do i=1,4
                      if(lface(i,iface)==ipoin)exit
                   enddo
                   if(i==5)then
                      write(*,*)'Error in ptoedgtet'
                      stop
                   endif
                   iedg=ltab(i,j)
                   ledglm(iedg,iface)=nedge
                else
                   !
                   !     Look in the old edges from ipedg
                   ! 
                   do iedge=ipedg,nedge
                      if(ledge(2,iedge)==ip1)exit 
                   enddo
                 
                   if(iedge==nedge+1)then
                      write(*,*)'Error in ptoedgtet 2'
                   endif
                   !
                   !     Find the position of ipoin
                   !
                   do i=1,4
                      if(lface(i,iface)==ipoin)exit
                   enddo
                   if(i==5)then
                      write(*,*)'Error in ptoedgtet 3'
                      stop
                   endif
                   iedg=ltab(i,j)
                   ledglm(iedg,iface)=iedge
 
                endif
             endif
          enddo
       enddo
       ipedg=nedge+1_ip
    enddo

    call memchk(2_ip,istat,memor_msh,'LMARK','ptoedgtet',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','ptoedgtet',0_ip)

  end subroutine ptoedgtet

  subroutine ptoedgtri(lface,nface,npoin,nnofa,ptoel1,ptoel2,&
       nedge,ledge,ledglm)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only          : memor_msh
    implicit none
    integer(ip), intent(in)            :: nface,npoin,nnofa
    integer(ip), intent(in)            :: ptoel1(*),ptoel2(npoin+1)
    integer(ip), intent(in)            :: lface(nnofa,nface)
    integer(ip), intent(inout)         :: nedge
    integer(ip),pointer                :: ledge(:,:),ledglm(:,:)
    integer(ip)                        :: iface,j,ipoin,ip1,iel,i,iedg,ipedg,iedge
    integer(ip), pointer               :: lmark(:)
    integer(4)                 :: istat
    integer(ip)               :: ltab(3,3)=RESHAPE((/0,3,2,&
                                                     3,0,1,&
                                                     2,1,0/),(/3,3/)) 
    !
    !     This routine gets the edges from elements lface
    !
    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','ptoedgtri',lmark)
    if(.not.associated(ledglm))then
       allocate(ledglm(nnofa,nface),stat=istat)
       call memchk(zero,istat,memor_msh,'LEDGLM','ptoedgtri',ledglm)
    else
       call memrea(nface,memor_msh,'LEDGLM','ptoedgtri',ledglm)
    endif
    !
    !     First count the edges
    !
    !
    !     Loop on the points 
    !
    nedge=0_ip
    do ipoin=1,npoin
       !
       !     Loop on the elements surrounding the point
       !
       do iel=ptoel2(ipoin),ptoel2(ipoin+1)-1
          iface=ptoel1(iel)
          do j=1,3
             ip1=lface(j,iface)
             if(ipoin<ip1)then
                if(lmark(ip1)/=ipoin)then
                   lmark(ip1)=ipoin 
                   nedge=nedge+1  
                endif
             endif
          enddo
       enddo
    enddo
    !
    !     Allocate the edges
    !
    if(.not.associated(ledge))then
       allocate(ledge(2,nedge),stat=istat)
       call memchk(zero,istat,memor_msh,'LEDGE','ptoedgtri',ledge)
    else
       call memrea(nedge,memor_msh,'LEDGE','ptoedgtri',ledge)
    endif
    !
    !     Reset lmark
    !
    do ipoin=1,npoin
       lmark(ipoin)=0_ip
    enddo
    !
    !     Fill the edges
    !
    !
    !     Loop on the points 
    !
    ipedg=1_ip
    nedge=0_ip
    do ipoin=1,npoin
       !
       !     Loop on the elements surrounding the point
       !
       do iel=ptoel2(ipoin),ptoel2(ipoin+1)-1
          iface=ptoel1(iel)
          do j=1,3
             ip1=lface(j,iface)
             if(ipoin<ip1)then
                !
                !     Is this a new edge
                !
                if(lmark(ip1)/=ipoin)then
                   lmark(ip1)=ipoin  
                   nedge=nedge+1  
                   ledge(1,nedge)=ipoin
                   ledge(2,nedge)=ip1
                   !
                   !     Find the position of ipoin
                   !
                   do i=1,3
                      if(lface(i,iface)==ipoin)exit
                   enddo
                   if(i==4)then
                      write(*,*)'Error in ptoedgtri'
                      stop
                   endif

                   iedg=ltab(i,j)
                   ledglm(iedg,iface)=nedge

                else   
                   !
                   !     Look in the old edges from ipedg
                   ! 
                   do iedge=ipedg,nedge
                      if(ledge(2,iedge)==ip1)exit 
                   enddo
                 
                   if(iedge==nedge+1)then
                      write(*,*)'Error in ptoedgtri 2'
                   endif
                   !
                   !     Find the position of ipoin
                   !
                   do i=1,3
                      if(lface(i,iface)==ipoin)exit
                   enddo
                   if(i==4)then
                      write(*,*)'Error in ptoedgtri 3'
                      stop
                   endif
                   iedg=ltab(i,j)
                   ledglm(iedg,iface)=iedge
               endif
             endif
          enddo
       enddo
       ipedg=nedge+1_ip
    enddo

    call memchk(2_ip,istat,memor_msh,'LMARK','ptoedgtri',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','ptoedgtri',0_ip)

  end subroutine ptoedgtri

  subroutine ptoedg2(lface, nface , npoin, nnofa,ptoel1,ptoel2,&
       nedge,ptoed2,ledge)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only          : memor_msh
    implicit none
    integer(ip), intent(in)            :: nface,npoin,nnofa,ptoel1(*)
    integer(ip), intent(in)            :: lface(nnofa,nface),ptoel2(*)
    integer(ip), intent(inout)         :: nedge
    integer(ip),pointer                :: ptoed2(:),ledge(:,:)
    integer(ip)                        :: iface,j,ipoin,ip1,iel
    integer(ip), pointer               :: lmark(:)
    integer(4)                 :: istat

    allocate(lmark(npoin),stat=istat)
    call memchk(zero,istat,memor_msh,'LMARK','ptoedg',lmark)
    !
    !     First count the edges
    !
    !
    !     Loop on the points
    !
    nedge=0_ip
    do ipoin=1,npoin
       !
       !     Loop on the elements surrounding the point
       !
       do iel=ptoel2(ipoin),ptoel2(ipoin+1)-1
          iface=ptoel1(iel)
          do j=1,nnofa
             ip1=lface(j,iface)
             if(ipoin<ip1)then
                if(lmark(ip1)/=ipoin)then
                lmark(ip1)=ipoin
                   nedge=nedge+1
                endif
             endif
          enddo
       enddo
    enddo
    !
    !     Allocate the edges
    !
    if(.not.associated(ledge))then
       allocate(ledge(2,nedge),stat=istat)
       call memchk(zero,istat,memor_msh,'LEDGE','ptoedg',ledge)
    else
       call memrea(nedge,memor_msh,'LEDGE','ptoedg',ledge)
    endif
    if(.not.associated(ptoed2))then
       allocate(ptoed2(npoin+1),stat=istat)
       call memchk(zero,istat,memor_msh,'PTOED2','ptoedg',ptoed2)
    else
       call memrea(npoin+1_ip,memor_msh,'PTOED2','ptoedg',ptoed2)
    endif
    !
    !     Reset lmark
    !
    do ipoin=1,npoin
       lmark(ipoin)=0_ip
    enddo
    !
    !     Fill the edges
    !
    !
    !     Loop on the points
    !
    nedge=0_ip
    ptoed2(1)=1_ip
    do ipoin=1,npoin
       !
       !     Loop on the elements surrounding the point
       !
       do iel=ptoel2(ipoin),ptoel2(ipoin+1)-1
          iface=ptoel1(iel)
          do j=1,nnofa
             ip1=lface(j,iface)
             if(ipoin<ip1)then
                if(lmark(ip1)/=ipoin)then
                   lmark(ip1)=ipoin
                   nedge=nedge+1
                   ledge(1,nedge)=ipoin
                   ledge(2,nedge)=ip1
                endif
             endif
          enddo
       enddo
       ptoed2(ipoin+1)=nedge+1
    enddo

    call memchk(2_ip,istat,memor_msh,'LMARK','ptoedg',lmark)
    deallocate(lmark,stat=istat)
    if(istat/=0) call memerr(2_ip,'LMARK','ptoedg',0_ip)

  end subroutine ptoedg2

  subroutine mstnrmal(rnrmal,nnrmal,ndim,mstnrm,ipoin,ierr)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only          : memor_msh
    implicit none
    integer(ip), intent(in)            :: nnrmal,ndim,ipoin
    integer(ip), intent(inout)         :: ierr
    real(rp),intent(in)                :: rnrmal(ndim,nnrmal)
    real(rp),intent(inout)             :: mstnrm(ndim)
    integer(ip)                       :: inrml,jnrml,knrml,lnrml
    real(rp)                          :: c05,rnx,rny,rnz,rnl,c10,csca,c00
    real(rp)                          :: cscamin,cscaloc,tolsca,tolcircle
    real(rp)                          :: cscaloc1,cscamin1,denom,rxc,ryc,rzc
    real(rp)                          :: scal1,scal2,scal3,epsil,epsil2 
    real(rp)                          :: rx12,ry12,rz12,rx23,ry23,rz23,rx31,ry31,rz31
    logical(lg)                       :: ifound

    c00=0.0d+00
    c05=0.5d+00
    c10=1.0d+00
    tolsca=0.9_rp
    cscamin=1.1d+00
    tolcircle=1.0d-06
    epsil=1.0d-12
    epsil2=1.0d-06
    ierr=0_ip
    !
    !     Try 'usual' method
    !
    rnx=c00
    rny=c00
    rnz=c00
    do inrml=1,nnrmal
       rnx=rnx+rnrmal(1,inrml) 
       rny=rny+rnrmal(2,inrml) 
       rnz=rnz+rnrmal(3,inrml) 
    enddo
    rnl=sqrt(rnx*rnx+rny*rny+rnz*rnz)
    rnl=c10/rnl
    rnx=rnx*rnl
    rny=rny*rnl
    rnz=rnz*rnl
    !
    !     Get max angle
    !
    do inrml=1,nnrmal
       csca=rnx*rnrmal(1,inrml)+rny*rnrmal(2,inrml)+rnz*rnrmal(3,inrml)
       if(csca<cscamin)cscamin=csca
    enddo
    !
    !     The test
    !
    if(cscamin>tolsca)then
       mstnrm(1)=rnx 
       mstnrm(2)=rny 
       mstnrm(3)=rnz 
       return
    endif
    !
    !     Second part: the real stuff 
    ! 
    !
    !     Initialize to not found
    !
    ifound=.false.
    !
    !     Initialize global angle (scalar product) to minimize
    !
    cscamin1=-1.01d+00
    cscamin=-1.01d+00 
    !
    !     First find the circles through two points
    !
    do inrml=1,nnrmal-1
       do jnrml=inrml+1,nnrmal
          !
          !     Take bisection (center in the plane)
          !
          rxc=c05*(rnrmal(1,inrml)+rnrmal(1,jnrml))
          ryc=c05*(rnrmal(2,inrml)+rnrmal(2,jnrml))
          rzc=c05*(rnrmal(3,inrml)+rnrmal(3,jnrml))

          rnl=sqrt(rxc*rxc+ryc*ryc+rzc*rzc)
          rnl=c10/rnl
          rxc=rxc*rnl
          ryc=ryc*rnl
          rzc=rzc*rnl
          !
          !     Get local angle (radius in the plane)
          !
          cscaloc=rxc*rnrmal(1,inrml)+ryc*rnrmal(2,inrml)+rzc*rnrmal(3,inrml)
          !
          !     Is the local circle smaller than the best candidate until now?
          !     If not, it is not worth following  
          !  
          if(cscaloc>cscamin1)then
             !
             !     Set tolerance
             !
             cscaloc1=cscaloc-tolcircle 
             !
             !     Loop on all the points
             !
             do knrml=1,nnrmal 
                if(knrml==inrml .or. knrml==jnrml)cycle
                !
                !     Get angle
                ! 
                csca=rxc*rnrmal(1,knrml)+ryc*rnrmal(2,knrml)+rzc*rnrmal(3,knrml)
                !
                !     Is knrml inside the local circle?   
                ! 
                if(csca<cscaloc1)exit

             enddo
             !
             !     Are all the points inside the circle or did we exit too soon?
             !
             if(knrml==nnrmal+1)then     
                !
                !     This is a candidate for the smallest circumscribed circle
                !
                ifound=.true.
                cscamin=cscaloc
                cscamin1=cscaloc+tolcircle
                mstnrm(1)=rxc 
                mstnrm(2)=ryc 
                mstnrm(3)=rzc 

             endif

          endif
       enddo
    enddo
    !
    !     Did we find already a solution?
    !     If yes, this is the solution
    !
    if(ifound)then
       return
    endif
    !
    !     Now try with three points
    ! 
    do inrml=1,nnrmal-2
       do jnrml=inrml+1,nnrmal-1
          do knrml=jnrml+1,nnrmal
             !
             !     Get relevant quantities
             !     Check input for equal normals  
             !
             rx12=rnrmal(1,inrml)-rnrmal(1,jnrml)
             ry12=rnrmal(2,inrml)-rnrmal(2,jnrml)
             rz12=rnrmal(3,inrml)-rnrmal(3,jnrml)
             rnl=abs(rx12)+abs(ry12)+abs(rz12)
             if(rnl<epsil2)cycle

             rx23=rnrmal(1,jnrml)-rnrmal(1,knrml)
             ry23=rnrmal(2,jnrml)-rnrmal(2,knrml)
             rz23=rnrmal(3,jnrml)-rnrmal(3,knrml)
             rnl=abs(rx23)+abs(ry23)+abs(rz23)
             if(rnl<epsil2)cycle

             rx31=rnrmal(1,knrml)-rnrmal(1,inrml)
             ry31=rnrmal(2,knrml)-rnrmal(2,inrml)
             rz31=rnrmal(3,knrml)-rnrmal(3,inrml)
             rnl=abs(rx31)+abs(ry31)+abs(rz31)
             if(rnl<epsil2)cycle
             !
             !     Compute center of circle
             !  
             !
             !      Three rotating denom by three cases
             !
             !
             ! First number permutation
             ! 
             denom=rx23*ry12-rx12*ry23
             if(abs(denom)>epsil)then
                denom=c10/denom
                rxc=-(rz23*ry12-rz12*ry23)*denom
                ryc=(rz23*rx12-rz12*rx23)*denom
                goto 7000
             endif

             denom=rx23*rz12-rx12*rz23
             if(abs(denom)>epsil)then
                denom=c10/denom
                rxc=-(ry23*rz12-ry12*rz23)*denom
                rzc=(ry23*rx12-ry12*rx23)*denom
                goto 7100
             endif

             denom=rz23*ry12-rz12*ry23
             if(abs(denom)>epsil)then
                denom=c10/denom
                rzc=-(rx23*ry12-rx12*ry23)*denom
                ryc=(rx23*rz12-rx12*rz23)*denom
                goto 7200
             endif
             !   
             !     Second number permutation
             ! 
             denom=rx31*ry23-rx23*ry31
             if(abs(denom)>epsil)then
                denom=c10/denom
                rxc=-(rz31*ry23-rz23*ry31)*denom
                ryc=(rz31*rx23-rz23*rx31)*denom
                goto 7000
             endif

             denom=rx31*rz23-rx23*rz31
             if(abs(denom)>epsil)then
                denom=c10/denom
                rxc=-(ry31*rz23-ry23*rz31)*denom
                rzc=(ry31*rx23-ry23*rx31)*denom
                goto 7100
             endif

             denom=rz31*ry23-rz23*ry31
             if(abs(denom)>epsil)then
                denom=c10/denom
                rzc=-(rx31*ry23-rx23*ry31)*denom
                ryc=(rx31*rz23-rx23*rz31)*denom
                goto 7200
             endif
             !
             !     Third number permutation
             !
             denom=rx12*ry31-rx31*ry12
             if(abs(denom)>epsil)then
                denom=c10/denom
                rxc=-(rz12*ry31-rz31*ry12)*denom
                ryc=(rz12*rx31-rz31*rx12)*denom
                goto 7000
             endif

             denom=rx12*rz31-rx31*rz12
             if(abs(denom)>epsil)then
                denom=c10/denom
                rxc=-(ry12*rz31-ry31*rz12)*denom
                rzc=(ry12*rx31-ry31*rx12)*denom
                goto 7100
             endif

             denom=rz12*ry31-rz31*ry12
             if(abs(denom)>epsil)then
                denom=c10/denom
                rzc=-(rx12*ry31-rx31*ry12)*denom
                ryc=(rx12*rz31-rx31*rz12)*denom
                goto 7200
             endif
             !
             !     At this point, the three normals must be aligned with a huge radius
             !         go home
             cycle

7000         continue
             rnl=sqrt(1+rxc*rxc+ryc*ryc)
             rzc=c10/(rnl)
             rxc=rxc*rzc
             ryc=ryc*rzc
             goto 7500

7100         continue
             rnl=sqrt(1+rxc*rxc+rzc*rzc)
             ryc=c10/(rnl)
             rxc=rxc*ryc
             rzc=rzc*ryc
             goto 7500

7200         continue
             rnl=sqrt(1+rzc*rzc+ryc*ryc)
             rxc=c10/(rnl)
             rzc=rzc*rxc
             ryc=ryc*rxc

7500         continue
             !  
             !      Compute the radius & orientation
             !
             scal1=rxc*rnrmal(1,knrml)+ryc*rnrmal(2,knrml)+rzc*rnrmal(3,knrml)
             if(scal1<0)then
                rxc=-rxc
                ryc=-ryc
                rzc=-rzc
                scal1=-scal1
             endif

             scal2=rxc*rnrmal(1,jnrml)+ryc*rnrmal(2,jnrml)+rzc*rnrmal(3,jnrml)
             scal3=rxc*rnrmal(1,inrml)+ryc*rnrmal(2,inrml)+rzc*rnrmal(3,inrml)
             cscaloc=scal1
             !
             !     Check that this radius is smaller that the current smallest radius
             !
             if(cscaloc>cscamin1)then
                !
                !     Set the tolerance
                !
                cscaloc1=cscaloc-tolcircle

                do lnrml=1,nnrmal
                   if(lnrml==inrml .or. lnrml==jnrml .or. lnrml==knrml)cycle
                   csca=rxc*rnrmal(1,lnrml)+ryc*rnrmal(2,lnrml)+rzc*rnrmal(3,lnrml)
                   if(csca<cscaloc1)exit
                enddo
                if(lnrml==nnrmal+1)then
                   !
                   !     Have found a circumcercle, store it
                   !
                   ifound=.true.
                   cscamin=cscaloc
                   cscamin1=cscaloc+tolcircle
                   mstnrm(1)=rxc
                   mstnrm(2)=ryc
                   mstnrm(3)=rzc
                endif
             endif
          enddo
       enddo
    enddo
    !
    !     Error 
    !
    if(.not.ifound)then
       write(*,*)'Error mstnrmal, most normal normal not found'
       write(*,*)'ipoin=',ipoin
       write(*,*)'nnrmal=',nnrmal
       do inrml=1,nnrmal
          write(*,*)rnrmal(1,inrml),rnrmal(2,inrml),rnrmal(3,inrml)
       enddo
       ierr=1_ip
    endif

  end subroutine mstnrmal

  subroutine hxtotet(ndim,nnode,coor,npoin,nelem,elem)
    use def_kintyp, only       :  ip,rp,lg
    use mod_memchk
    use def_meshin, only          : memor_msh
    implicit none
    integer(ip), intent(in)            :: ndim,nnode,npoin
    integer(ip), intent(inout)         :: nelem
    real(rp),intent(in)                :: coor(ndim,npoin)
    integer(ip),pointer                :: elem(:,:)
    integer(ip),pointer                :: elemn(:,:)
    integer(ip)                        :: nelem0,ielem,ip1,ip2,ip3,ip4,ip5,ip6,ip7,ip8
    integer(4)                 :: istat
    !
    !     This sub converts an hexahedral mesh to a tetra mesh
    !
    allocate(elemn(nnode,6*nelem),stat=istat)
    call memchk(zero,istat,memor_msh,'ELEMN','hextotet',elemn)

    nelem0=nelem
    nelem=0_ip
    do ielem=1,nelem0
       ip1=elem(1,ielem)
       ip2=elem(2,ielem)
       ip3=elem(3,ielem)
       ip4=elem(4,ielem)
       ip5=elem(5,ielem)
       ip6=elem(6,ielem)
       ip7=elem(7,ielem)
       ip8=elem(8,ielem)

       nelem=nelem+1
       elemn(1,nelem)=ip1
       elemn(2,nelem)=ip2
       elemn(3,nelem)=ip3
       elemn(4,nelem)=ip5

       nelem=nelem+1
       elemn(1,nelem)=ip5
       elemn(2,nelem)=ip2
       elemn(3,nelem)=ip3
       elemn(4,nelem)=ip7

       nelem=nelem+1
       elemn(1,nelem)=ip2
       elemn(2,nelem)=ip7
       elemn(3,nelem)=ip5
       elemn(4,nelem)=ip6

       nelem=nelem+1
       elemn(1,nelem)=ip1
       elemn(2,nelem)=ip3
       elemn(3,nelem)=ip4
       elemn(4,nelem)=ip8

       nelem=nelem+1
       elemn(1,nelem)=ip1
       elemn(2,nelem)=ip3
       elemn(3,nelem)=ip8
       elemn(4,nelem)=ip5

       nelem=nelem+1
       elemn(1,nelem)=ip5
       elemn(2,nelem)=ip3
       elemn(3,nelem)=ip8
       elemn(4,nelem)=ip7

    enddo

    call memchk(2_ip,istat,memor_msh,'ELEM','ptoedg',elem)
    deallocate(elem,stat=istat)
    if(istat/=0) call memerr(2_ip,'ELEM','ptoedg',0_ip)

    elem=>elemn

  end subroutine hxtotet

  subroutine outmsh(nnode,nelem,elem,ndim,npoin,coor)
    use def_kintyp, only       : ip,rp,lg
    implicit none
    integer(ip),intent(in)      :: nnode,nelem,npoin,ndim
    integer(ip),intent(in)      :: elem(nnode,nelem)
    real(rp), intent(in)        :: coor(ndim,npoin)
    real(rp)                    :: rx,ry,rz
    integer(ip)                 :: ipoin,ielem 
    !
    !     Print mesh
    !
    open(unit=70,file='outhexmsh.msh',status='unknown')
    rewind 70

1   format('MESH dimension 3 ElemType Tetrahedra Nnode 4')
2   format('Coordinates')
3   format('#node number   coord_x   coord_y  coord_z')
100 format(i10,3e20.10e3)
200 format(5i10)
4   format('end coordinates')
5   format('Elements')
6   format('end elements')
    write(70,1)
    write(70,2)
    write(70,3)
    do  ipoin=1,npoin
       rx=coor(1,ipoin)
       ry=coor(2,ipoin)
       rz=coor(3,ipoin)
       write(70,100)ipoin,rx,ry,rz
    enddo

    write(70,4)
    write(70,5)
    do  ielem=1,nelem
       write(70,200)ielem,elem(1,ielem),elem(2,ielem),elem(3,ielem),elem(4,ielem)
    enddo
    write(70,6)
    close(70)

  end subroutine outmsh

    
end module mod_mshtol
