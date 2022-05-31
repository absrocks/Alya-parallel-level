subroutine cargeo(ncell,lcell)
  use def_kintyp, only          :  ip,rp,lg,cell
  use mod_memchk
  use mod_cart
  use def_domain, only          :  lnods,ltype,coord,mnode,npoin,nelem,lhang,mnodb,nhang,nboun
  use def_domain, only          :  lboel,ltypb,lnodb,kfl_codbo,kfl_codno,mcono
  use def_elmtyp
  use def_meshin, only          :  memor_msh
  use mod_domain, only   :  domain_memory_allocate
  implicit none
  integer(ip),intent(inout)     :: ncell
  integer(ip)                   :: icont,iconte,icont1,icont2,icont3,icont4,icont5,icont6,icont7
  integer(ip)                   :: icont8,i,icontt,index1,ipnew,npold,iplace,jplace,nboup,iboup,p1
  integer(ip)                   ::  neigh,level,pt1,pt2,j,k,ptnew,icontpp,iter,ncell8,jpos,jcont,tab4(3)
  integer(ip)                   :: tab(4,6),tab2(4,6),icell,ilevel,ineigh,jlevel,jneigh,ipoin,inofa,nnofam1
  integer(ip)                   :: tab3(6),tab5(4,6),tab6(4,6),ipo,npnew,ncell0,ncont,nnofa,nface,iface,tab7(6),npbound
  integer(ip),pointer           :: renum(:),renum2(:),lmark(:),lsize(:),lmarkc(:),lface(:,:),lboup(:)
  real(rp), pointer             :: rsize(:)                   
  integer(4)                    :: istat
  type(cell)                    :: lcell(*)

  nnofa=5_ip
  nnofam1=nnofa-1_ip

  tab7(1)=1
  tab7(2)=2
  tab7(3)=3
  tab7(4)=4
  tab7(5)=5
  tab7(6)=6

  tab(1,1)=1
  tab(2,1)=5
  tab(3,1)=8
  tab(4,1)=4
  tab(1,2)=2
  tab(2,2)=3
  tab(3,2)=7
  tab(4,2)=6
  tab(1,3)=1
  tab(2,3)=2
  tab(3,3)=6
  tab(4,3)=5
  tab(1,4)=3
  tab(2,4)=4
  tab(3,4)=8
  tab(4,4)=7
  tab(1,5)=4
  tab(2,5)=3
  tab(3,5)=2
  tab(4,5)=1
  tab(1,6)=5
  tab(2,6)=6
  tab(3,6)=7
  tab(4,6)=8


  tab2(1,1)=2
  tab2(2,1)=6
  tab2(3,1)=7
  tab2(4,1)=3
  tab2(1,2)=1
  tab2(2,2)=4
  tab2(3,2)=8
  tab2(4,2)=5
  tab2(1,3)=4
  tab2(2,3)=3
  tab2(3,3)=7
  tab2(4,3)=8
  tab2(1,4)=2
  tab2(2,4)=1
  tab2(3,4)=5
  tab2(4,4)=6
  tab2(1,5)=8
  tab2(2,5)=7
  tab2(3,5)=6
  tab2(4,5)=5
  tab2(1,6)=1
  tab2(2,6)=2
  tab2(3,6)=3
  tab2(4,6)=4

  tab3(1)=2
  tab3(2)=1
  tab3(3)=4
  tab3(4)=3
  tab3(5)=6
  tab3(6)=5

  tab5(1,1)=1
  tab5(2,1)=5
  tab5(3,1)=4
  tab5(4,1)=8
  tab5(1,2)=2
  tab5(2,2)=6
  tab5(3,2)=3
  tab5(4,2)=7
  tab5(1,3)=1
  tab5(2,3)=5
  tab5(3,3)=2
  tab5(4,3)=6
  tab5(1,4)=4
  tab5(2,4)=8
  tab5(3,4)=3
  tab5(4,4)=7
  tab5(1,5)=1
  tab5(2,5)=4
  tab5(3,5)=2
  tab5(4,5)=3
  tab5(1,6)=5
  tab5(2,6)=8
  tab5(3,6)=6
  tab5(4,6)=7

  tab6(1,1)=2
  tab6(2,1)=6
  tab6(3,1)=3
  tab6(4,1)=7
  tab6(1,2)=1
  tab6(2,2)=5
  tab6(3,2)=4
  tab6(4,2)=8
  tab6(1,3)=4
  tab6(2,3)=8
  tab6(3,3)=3
  tab6(4,3)=7
  tab6(1,4)=1
  tab6(2,4)=5
  tab6(3,4)=2
  tab6(4,4)=6
  tab6(1,5)=5
  tab6(2,5)=8
  tab6(3,5)=6
  tab6(4,5)=7
  tab6(1,6)=1
  tab6(2,6)=4
  tab6(3,6)=2
  tab6(4,6)=3

  !
  !     Allocate renum and renum2
  !
  ncell8=ncell*8 
  allocate(renum(ncell8),stat=istat)
  call memchk(zero,istat,memor_msh,'RENUM','writeGidconf',renum)
  allocate(renum2(ncell8),stat=istat)
  call memchk(zero,istat,memor_msh,'RENUM2','writeGidconf',renum2)
    !
  !     Initialize renum
  !
  do i=1,ncell8
     renum(i)=i
  enddo
  !
  !     Conformize the cells
  !
  do
     !
     !     Initialize the "something done" check 
     !
     iter=0
     !
     !     Loop on the cells
     !
     do icell=1,ncell

        ilevel=lcell(icell)%level
        !
        !     Loop on the neighbors
        !
        do i=1,6

           ineigh=lcell(icell)%neigh(i)

           if(ineigh==0)cycle

           jlevel=lcell(ineigh)%level

           if(jlevel==ilevel)then
              !
              !     Loop on the points of the faces
              !
              do k=1,4

                 pt1=(icell-1)*8+tab(k,i) 
                 pt2=(ineigh-1)*8+tab2(k,i)
                 if(renum(pt1)<renum(pt2))then
                    renum(pt2)=renum(pt1)
                    iter=1
                 else if(renum(pt2)<renum(pt1))then   
                    renum(pt1)=renum(pt2)
                    iter=1
                 endif

              enddo

           else if(ilevel<jlevel)then

              jpos=ineigh+1_ip
              jcont=1_ip
              do j=1,8
                 if(lcell(jpos)%neigh(tab3(i))==icell)then
                    tab4(jcont)=jpos
                    jcont=jcont+1_ip
                    if(jcont==4)exit 
                 endif
                 jpos=jpos+1  
              enddo

              pt1=(icell-1)*8+tab5(1,i)
              pt2=(ineigh-1)*8+tab6(1,i)
              if(renum(pt1)<renum(pt2))then
                 renum(pt2)=renum(pt1)
                 iter=1
              else if(renum(pt2)<renum(pt1))then   
                 renum(pt1)=renum(pt2)
                 iter=1
              endif
              
              pt1=(icell-1)*8+tab5(2,i)
              pt2=(tab4(1)-1)*8+tab6(2,i)
              if(renum(pt1)<renum(pt2))then
                 renum(pt2)=renum(pt1)
                 iter=1
              else if(renum(pt2)<renum(pt1))then   
                 renum(pt1)=renum(pt2)
                 iter=1
              endif

              pt1=(icell-1)*8+tab5(3,i)
              pt2=(tab4(2)-1)*8+tab6(3,i)
              if(renum(pt1)<renum(pt2))then
                 renum(pt2)=renum(pt1)
                 iter=1
              else if(renum(pt2)<renum(pt1))then   
                 renum(pt1)=renum(pt2)
                 iter=1
              endif
              
              pt1=(icell-1)*8+tab5(4,i)
              pt2=(tab4(3)-1)*8+tab6(4,i)
              if(renum(pt1)<renum(pt2))then
                 renum(pt2)=renum(pt1)
                 iter=1
              else if(renum(pt2)<renum(pt1))then   
                 renum(pt1)=renum(pt2)
                 iter=1
              endif

           endif
        enddo
     enddo

     do i=1,ncell8
        renum(i)=renum(renum(i))
     enddo
  !
  !     DBG
  !
  !call chkarraycart(lcell,ncell,renum)

     if(iter==0)exit
     exit  

  enddo
  !call chkarraycart(lcell,ncell,renum)
  !
  !     Get the number and the new renumbering of the conformed points
  !
  icont=0
  do i=1,ncell8
     if(renum(i)==i)then
        icont=icont+1
        renum2(i)=icont
     endif
  enddo
  !
  !     Remember the conforming point number
  ! 
  npnew=icont
  !
  !     Put the new renumbering in renum
  !
  do i=1,ncell8
     renum(i)=renum2(renum(i))
  enddo
  !
  !     Get the original boundary faces
  !
  nface=0_ip
  do icell=1,ncell
     do j=1,6
        ineigh=lcell(icell)%neigh(j)
        if(ineigh==0)then
           nface=nface+1
        endif
     enddo
  enddo
  !
  !     Allocate lface
  !   
  allocate(lface(10,nface),stat=istat)
  call memchk(zero,istat,memor_msh,'LFACE','writeGidconf',lface)
  !
  !     Store in lface
  !
  nface=0_ip
  do icell=1,ncell
     do j=1,6
        ineigh=lcell(icell)%neigh(j)
        if(ineigh==0)then
           nface=nface+1
           lface(1,nface)=renum(8*(icell-1)+tab(1,j))
           lface(2,nface)=renum(8*(icell-1)+tab(2,j))
           lface(3,nface)=renum(8*(icell-1)+tab(3,j))
           lface(4,nface)=renum(8*(icell-1)+tab(4,j))
           lface(5,nface)=tab(1,j)
           lface(6,nface)=tab(2,j)
           lface(7,nface)=tab(3,j)
           lface(8,nface)=tab(4,j)
           lface(9,nface)=icell
           lface(10,nface)=tab7(j) 
        endif
     enddo
  enddo
  !
  !     Allocate lmark for the hanging nodes
  !
  nhang=npnew
  allocate(lmark(nhang),stat=istat)
  call memchk(zero,istat,memor_msh,'LMARK','writeGidconf',lmark)
  !
  !     Get the hanging nodes 
  !
  do icell=1,ncell

     ilevel=lcell(icell)%level
     !
     !     In X -
     !  
     ineigh=lcell(icell)%neigh(1)
     if(ineigh/=0)then
        jlevel=lcell(ineigh)%level
        if(ilevel<jlevel)then
           jneigh=lcell(ineigh)%neigh(2)
           if(icell/=jneigh)then 
              write(*,*)'Error conformity in hanging nodes 1'
              stop
           endif

           jpos=ineigh+1_ip
           jcont=1_ip
           do j=1,8
              if(lcell(jpos)%neigh(2)==icell)then
                 tab4(jcont)=jpos
                 jcont=jcont+1_ip
                 if(jcont==4)exit 
              endif
              jpos=jpos+1  
           enddo
           !
           !     Node at the center of the face
           !     
           ipoin=renum(8*(ineigh-1)+7_ip) 
           if(lmark(ipoin)==0)then
              lmark(ipoin)=1_ip
             ! HN lhang(0,ipoin)=4
             ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+2_ip) 
             ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+6_ip) 
             ! HN lhang(3,ipoin)=renum(8*(tab4(2)-1)+3_ip) 
             ! HN lhang(4,ipoin)=renum(8*(tab4(3)-1)+7_ip) 
           endif
           !
           !     Nodes at the edges
           !
           ipoin=renum(8*(ineigh-1)+6_ip)
           if(lmark(ipoin)==0)then
              lmark(ipoin)=1_ip
             ! HN lhang(0,ipoin)=2
             ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+2_ip) 
             ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+6_ip) 
           endif
           ipoin=renum(8*(ineigh-1)+3_ip)
           if(lmark(ipoin)==0)then
              lmark(ipoin)=1_ip
             ! HN lhang(0,ipoin)=2
             ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+2_ip) 
             ! HN lhang(2,ipoin)=renum(8*(tab4(2)-1)+3_ip) 
           endif
           ipoin=renum(8*(tab4(1)-1)+7_ip)
           if(lmark(ipoin)==0)then
              lmark(ipoin)=1_ip
             ! HN lhang(0,ipoin)=2
             ! HN lhang(1,ipoin)=renum(8*(tab4(1)-1)+6_ip) 
             ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+7_ip) 
           endif
           ipoin=renum(8*(tab4(2)-1)+7_ip)
           if(lmark(ipoin)==0)then
              lmark(ipoin)=1_ip
             ! HN lhang(0,ipoin)=2
             ! HN lhang(1,ipoin)=renum(8*(tab4(2)-1)+3_ip) 
             ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+7_ip) 
           endif
        endif
     endif
     !
     !     In X +
     !  
     ineigh=lcell(icell)%neigh(2)
     if(ineigh/=0)then  
        jlevel=lcell(ineigh)%level
        if(ilevel<jlevel)then
           jneigh=lcell(ineigh)%neigh(1)
           if(icell/=jneigh)then 
              write(*,*)'Error conformity in hanging nodes 2'
              stop
           endif

           jpos=ineigh+1_ip
           jcont=1_ip
           do j=1,8
              if(lcell(jpos)%neigh(1)==icell)then
                 tab4(jcont)=jpos
                 jcont=jcont+1_ip
                 if(jcont==4)exit 
              endif
              jpos=jpos+1  
           enddo
           !
           !     Node at the center of the face
           !     
           ipoin=renum(8*(ineigh-1)+8_ip) 
           if(lmark(ipoin)==0)then
              lmark(ipoin)=1_ip
             ! HN lhang(0,ipoin)=4
             ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+1_ip) 
             ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+5_ip) 
             ! HN lhang(3,ipoin)=renum(8*(tab4(2)-1)+4_ip) 
             ! HN lhang(4,ipoin)=renum(8*(tab4(3)-1)+8_ip) 
           endif
           !
           !     Nodes at the edges
           !
           ipoin=renum(8*(ineigh-1)+5_ip)
           if(lmark(ipoin)==0)then
              lmark(ipoin)=1_ip
             ! HN lhang(0,ipoin)=2
             ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+1_ip) 
             ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+5_ip) 
           endif
           ipoin=renum(8*(ineigh-1)+4_ip)
           if(lmark(ipoin)==0)then
              lmark(ipoin)=1_ip
             ! HN lhang(0,ipoin)=2
             ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+1_ip) 
             ! HN lhang(2,ipoin)=renum(8*(tab4(2)-1)+4_ip) 
           endif
           ipoin=renum(8*(tab4(1)-1)+8_ip)
           if(lmark(ipoin)==0)then
              lmark(ipoin)=1_ip
             ! HN lhang(0,ipoin)=2
             ! HN lhang(1,ipoin)=renum(8*(tab4(1)-1)+5_ip) 
             ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+8_ip) 
           endif
           ipoin=renum(8*(tab4(2)-1)+8_ip)
           if(lmark(ipoin)==0)then
              lmark(ipoin)=1_ip
             ! HN lhang(0,ipoin)=2
             ! HN lhang(1,ipoin)=renum(8*(tab4(2)-1)+4_ip) 
             ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+8_ip) 
           endif
        endif
     endif
     !
     !     In Y -
     !  
     ineigh=lcell(icell)%neigh(3)
     if(ineigh/=0)then  
        jlevel=lcell(ineigh)%level
        if(ilevel<jlevel)then
           jneigh=lcell(ineigh)%neigh(4)
           if(icell/=jneigh)then 
              write(*,*)'Error conformity in hanging nodes 3'
              stop
           endif

           jpos=ineigh+1_ip
           jcont=1_ip
           do j=1,8
              if(lcell(jpos)%neigh(4)==icell)then
                 tab4(jcont)=jpos
                 jcont=jcont+1_ip
                 if(jcont==4)exit 
              endif
              jpos=jpos+1  
           enddo
           !
           !     Node at the center of the face
           !     
           ipoin=renum(8*(ineigh-1)+7_ip) 
           if(lmark(ipoin)==0)then
              lmark(ipoin)=1_ip
             ! HN lhang(0,ipoin)=4
             ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+4_ip) 
             ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+8_ip) 
             ! HN lhang(3,ipoin)=renum(8*(tab4(2)-1)+3_ip) 
             ! HN lhang(4,ipoin)=renum(8*(tab4(3)-1)+7_ip) 
           endif
           !
           !     Nodes at the edges
           !
           ipoin=renum(8*(ineigh-1)+8_ip)
           if(lmark(ipoin)==0)then
              lmark(ipoin)=1_ip
             ! HN lhang(0,ipoin)=2
             ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+4_ip) 
             ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+8_ip) 
           endif
           ipoin=renum(8*(ineigh-1)+3_ip)
           if(lmark(ipoin)==0)then
              lmark(ipoin)=1_ip
             ! HN lhang(0,ipoin)=2
             ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+4_ip) 
             ! HN lhang(2,ipoin)=renum(8*(tab4(2)-1)+3_ip) 
           endif
           ipoin=renum(8*(tab4(1)-1)+7_ip)
           if(lmark(ipoin)==0)then
              lmark(ipoin)=1_ip
             ! HN lhang(0,ipoin)=2
             ! HN lhang(1,ipoin)=renum(8*(tab4(1)-1)+8_ip) 
             ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+7_ip) 
           endif
           ipoin=renum(8*(tab4(2)-1)+7_ip)
           if(lmark(ipoin)==0)then
              lmark(ipoin)=1_ip
             ! HN lhang(0,ipoin)=2
             ! HN lhang(1,ipoin)=renum(8*(tab4(2)-1)+3_ip) 
             ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+7_ip) 
           endif
        endif
     endif
     !
     !     In Y +
     !  
     ineigh=lcell(icell)%neigh(4)
     if(ineigh/=0)then  
        jlevel=lcell(ineigh)%level
        if(ilevel<jlevel)then
           jneigh=lcell(ineigh)%neigh(3)
           if(icell/=jneigh)then 
              write(*,*)'Error conformity in hanging nodes 4'
              stop
           endif

           jpos=ineigh+1_ip
           jcont=1_ip
           do j=1,8
              if(lcell(jpos)%neigh(3)==icell)then
                 tab4(jcont)=jpos
                 jcont=jcont+1_ip
                 if(jcont==4)exit 
              endif
              jpos=jpos+1  
           enddo
           !
           !     Node at the center of the face
           !     
           ipoin=renum(8*(ineigh-1)+6_ip) 
           if(lmark(ipoin)==0)then
              lmark(ipoin)=1_ip
             ! HN lhang(0,ipoin)=4
             ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+1_ip) 
             ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+5_ip) 
             ! HN lhang(3,ipoin)=renum(8*(tab4(2)-1)+2_ip) 
             ! HN lhang(4,ipoin)=renum(8*(tab4(3)-1)+6_ip) 
           endif
           !
           !     Nodes at the edges
           !
           ipoin=renum(8*(ineigh-1)+5_ip)
           if(lmark(ipoin)==0)then
              lmark(ipoin)=1_ip
             ! HN lhang(0,ipoin)=2
             ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+1_ip) 
             ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+5_ip) 
           endif
           ipoin=renum(8*(ineigh-1)+2_ip)
           if(lmark(ipoin)==0)then
              lmark(ipoin)=1_ip
             ! HN lhang(0,ipoin)=2
             ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+1_ip) 
             ! HN lhang(2,ipoin)=renum(8*(tab4(2)-1)+2_ip) 
           endif
           ipoin=renum(8*(tab4(1)-1)+6_ip)
           if(lmark(ipoin)==0)then
              lmark(ipoin)=1_ip
             ! HN lhang(0,ipoin)=2
             ! HN lhang(1,ipoin)=renum(8*(tab4(1)-1)+5_ip) 
             ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+6_ip) 
           endif
           ipoin=renum(8*(tab4(2)-1)+6_ip)
           if(lmark(ipoin)==0)then
              lmark(ipoin)=1_ip
             ! HN lhang(0,ipoin)=2
             ! HN lhang(1,ipoin)=renum(8*(tab4(2)-1)+2_ip) 
             ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+6_ip) 
           endif
        endif
     endif
     !
     !     In Z -
     !  
     ineigh=lcell(icell)%neigh(5)
     if(ineigh/=0)then  
        jlevel=lcell(ineigh)%level
        if(ilevel<jlevel)then
           jneigh=lcell(ineigh)%neigh(6)
           if(icell/=jneigh)then 
              write(*,*)'Error conformity in hanging nodes 5'
              stop
           endif

           jpos=ineigh+1_ip
           jcont=1_ip
           do j=1,8
              if(lcell(jpos)%neigh(6)==icell)then
                 tab4(jcont)=jpos
                 jcont=jcont+1_ip
                 if(jcont==4)exit 
              endif
              jpos=jpos+1  
           enddo
           !
           !     Node at the center of the face
           !     
           ipoin=renum(8*(ineigh-1)+7_ip) 
           if(lmark(ipoin)==0)then
              lmark(ipoin)=1_ip
             ! HN lhang(0,ipoin)=4
             ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+5_ip) 
             ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+8_ip) 
             ! HN lhang(3,ipoin)=renum(8*(tab4(2)-1)+6_ip) 
             ! HN lhang(4,ipoin)=renum(8*(tab4(3)-1)+7_ip) 
           endif
           !
           !     Nodes at the edges
           !
           ipoin=renum(8*(ineigh-1)+8_ip)
           if(lmark(ipoin)==0)then
              lmark(ipoin)=1_ip
             ! HN lhang(0,ipoin)=2
             ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+5_ip) 
             ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+8_ip) 
           endif
           ipoin=renum(8*(ineigh-1)+6_ip)
           if(lmark(ipoin)==0)then
              lmark(ipoin)=1_ip
             ! HN lhang(0,ipoin)=2
             ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+5_ip) 
             ! HN lhang(2,ipoin)=renum(8*(tab4(2)-1)+6_ip) 
           endif
           ipoin=renum(8*(tab4(1)-1)+7_ip)
           if(lmark(ipoin)==0)then
              lmark(ipoin)=1_ip
             ! HN lhang(0,ipoin)=2
             ! HN lhang(1,ipoin)=renum(8*(tab4(1)-1)+8_ip) 
             ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+7_ip) 
           endif
           ipoin=renum(8*(tab4(2)-1)+7_ip)
           if(lmark(ipoin)==0)then
              lmark(ipoin)=1_ip
             ! HN lhang(0,ipoin)=2
             ! HN lhang(1,ipoin)=renum(8*(tab4(2)-1)+6_ip) 
             ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+7_ip) 
           endif
        endif
     endif
     !
     !     In Z +
     !  
     ineigh=lcell(icell)%neigh(6)
     if(ineigh/=0)then  
        jlevel=lcell(ineigh)%level
        if(ilevel<jlevel)then
           jneigh=lcell(ineigh)%neigh(5)
           if(icell/=jneigh)then 
              write(*,*)'Error conformity in hanging nodes 6'
              stop
           endif

           jpos=ineigh+1_ip
           jcont=1_ip
           do j=1,8
              if(lcell(jpos)%neigh(5)==icell)then
                 tab4(jcont)=jpos
                 jcont=jcont+1_ip
                 if(jcont==4)exit 
              endif
              jpos=jpos+1  
           enddo
           !
           !     Node at the center of the face
           !     
           ipoin=renum(8*(ineigh-1)+3_ip) 
           if(lmark(ipoin)==0)then
              lmark(ipoin)=1_ip
             ! HN lhang(0,ipoin)=4
             ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+1_ip) 
             ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+4_ip) 
             ! HN lhang(3,ipoin)=renum(8*(tab4(2)-1)+2_ip) 
             ! HN lhang(4,ipoin)=renum(8*(tab4(3)-1)+3_ip) 
           endif
           !
           !     Nodes at the edges
           !
           ipoin=renum(8*(ineigh-1)+4_ip)
           if(lmark(ipoin)==0)then
              lmark(ipoin)=1_ip
             ! HN lhang(0,ipoin)=2
             ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+1_ip) 
             ! HN lhang(2,ipoin)=renum(8*(tab4(1)-1)+4_ip) 
           endif
           ipoin=renum(8*(ineigh-1)+2_ip)
           if(lmark(ipoin)==0)then
              lmark(ipoin)=1_ip
             ! HN lhang(0,ipoin)=2
             ! HN lhang(1,ipoin)=renum(8*(ineigh-1)+1_ip) 
             ! HN lhang(2,ipoin)=renum(8*(tab4(2)-1)+2_ip) 
           endif
           ipoin=renum(8*(tab4(1)-1)+3_ip)
           if(lmark(ipoin)==0)then
              lmark(ipoin)=1_ip
             ! HN lhang(0,ipoin)=2
             ! HN lhang(1,ipoin)=renum(8*(tab4(1)-1)+4_ip) 
             ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+3_ip) 
           endif
           ipoin=renum(8*(tab4(2)-1)+3_ip)
           if(lmark(ipoin)==0)then
              lmark(ipoin)=1_ip
             ! HN lhang(0,ipoin)=2
             ! HN lhang(1,ipoin)=renum(8*(tab4(2)-1)+2_ip) 
             ! HN lhang(2,ipoin)=renum(8*(tab4(3)-1)+3_ip) 
           endif
        endif
     endif
  enddo
  !
  !     Allocate lmark
  !  
  allocate(lmarkc(ncell),stat=istat)
  call memchk(zero,istat,memor_msh,'LMARKC','takeout',lmarkc)
  !
  !     Get the inner cells 
  !
  call carout(ncell,lcell,lmarkc)
  !
  !     DBG
  !
  !call chkarraycart(lcell,ncell,lmarkc)
  !
  !     Clean up renum2     
  ! 
  do ipoin=1,npnew
     renum2(ipoin)=0_ip
  enddo
  !
  !     Allocate lboup
  !
  allocate(lboup(npnew),stat=istat)
  call memchk(zero,istat,memor_msh,'LBOUP','writeGidconf',lboup)
  !
  !     Mark the exterior points
  !
  do icell=1,ncell 
     if(lmarkc(icell)==1)then
        iplace=8*(icell-1)+1 
        ipoin=renum(iplace)
        renum2(ipoin)=1
        iplace=iplace+1 
        ipoin=renum(iplace)
        renum2(ipoin)=1
        iplace=iplace+1 
        ipoin=renum(iplace)
        renum2(ipoin)=1
        iplace=iplace+1 
        ipoin=renum(iplace)
        renum2(ipoin)=1
        iplace=iplace+1 
        ipoin=renum(iplace)
        renum2(ipoin)=1
        iplace=iplace+1 
        ipoin=renum(iplace)
        renum2(ipoin)=1
        iplace=iplace+1 
        ipoin=renum(iplace)
        renum2(ipoin)=1
        iplace=iplace+1 
        ipoin=renum(iplace)
        renum2(ipoin)=1
     endif
  enddo 
  !
  !     Mark the interior points already marked
  !
  do icell=1,ncell 
     if(lmarkc(icell)/=1)then
        iplace=8*(icell-1)+1 
        ipoin=renum(iplace)
        if(renum2(ipoin)==1)then
           renum2(ipoin)=2
        endif 
        iplace=iplace+1 
        ipoin=renum(iplace)
        if(renum2(ipoin)==1)then
           renum2(ipoin)=2
        endif  
        ipoin=renum(iplace)
        if(renum2(ipoin)==1)then
           renum2(ipoin)=2
        endif  
        ipoin=renum(iplace)
        if(renum2(ipoin)==1)then
           renum2(ipoin)=2
        endif  
        ipoin=renum(iplace)
        if(renum2(ipoin)==1)then
           renum2(ipoin)=2
        endif  
        ipoin=renum(iplace)
        if(renum2(ipoin)==1)then
           renum2(ipoin)=2
        endif  
        ipoin=renum(iplace)
        if(renum2(ipoin)==1)then
           renum2(ipoin)=2
        endif  
        ipoin=renum(iplace)
        if(renum2(ipoin)==1)then
           renum2(ipoin)=2
        endif  
     endif
  enddo
  nboup=0 
  !
  !     Get the boundary points
  !
  do ipoin=1,npnew
     if(renum2(ipoin)==2)then
        nboup=nboup+1
        lboup(nboup)=ipoin
     endif
  enddo
  !
  !     Get the inner points to be deleted
  !
  npbound=0
  do icell=1,ncell 
    if(lmarkc(icell)/=1)then     
       do j=1,8
          iplace=8*(icell-1)+j
          ipoin=renum(iplace)
          !if(renum2(ipoin)/=2)then   !NOT NECESSARY AS MARKED BY OUTER  
             renum(iplace)=-ipoin     ! ELEMENTS
          !endif
        enddo
     endif
   enddo      
  
  !call chkarraycart(lcell,ncell,renum)
  !
  !     Clean up renum2     
  ! 
  do ipoin=1,npnew
     renum2(ipoin)=0_ip
  enddo
  !
  !     Transfer the valid points from cells to points 
  !
  do ipoin=1,ncell8 
     ipnew=renum(ipoin)
     if(ipnew>0)then
        renum2(ipnew)=1_ip
     endif
  enddo
  !
  !     Add the faces created by the holes 
  !
  !call outbou(ncell,nnofam1,npnew,nface,lcell,lmarkc,lface,renum,renum2)
  !
  !     Renumber the points
  !
  icont=0_ip
  do ipoin=1,npnew
     if(renum2(ipoin)==1)then
        icont=icont+1_ip 
        renum2(ipoin)=icont
     endif
  enddo
  npold=npnew
  npnew=icont
  !
  !     Renumber the hanging nodes
  !
  icont=0_ip
  do ipoin=1,npold
     !
     !     Should the point be kept?
     !
     if(renum2(ipoin)>0)then
        icont=icont+1_ip

        if(lmark(ipoin)==1)then
           !
           !     Count if all the master nodes are still there
           !
           ! HN
           !ncont=0_ip
           !do jcont=1,lhang(0,ipoin)
           !   if(lhang(jcont,ipoin)/=0)then
           !      if(renum2(lhang(jcont,ipoin))/=0)then
           !         ncont=ncont+1_ip
           !      endif
           !   endif
           !enddo
           ! HN
           !
           !     If yes, this is still a valid hanging node
           !  
           ! HN if(ncont==lhang(0,ipoin))then
           if(ncont==ncont)then

              lmark(icont)=1_ip
             ! HN lhang(0,icont)=lhang(0,ipoin)
             ! HN lhang(1,icont)=renum2(lhang(1,ipoin))
             ! HN lhang(2,icont)=renum2(lhang(2,ipoin))
              ! HN if(lhang(3,ipoin)/=0)then
             if(0/=0)then
                ! HN lhang(3,icont)=renum2(lhang(3,ipoin))
              else
                ! HN lhang(3,icont)=0_ip
              endif
              ! HN if(lhang(4,ipoin)/=0)then
              if(0/=0)then
                ! HN lhang(4,icont)=renum2(lhang(4,ipoin))
              else
                ! HN lhang(4,icont)=0_ip
              endif
             ! HN lhang(5,icont)=0_ip 
           
           else 
           
              lmark(icont)=0_ip
             ! HN lhang(0,icont)=0_ip
             ! HN lhang(1,icont)=0_ip
             ! HN lhang(2,icont)=0_ip
             ! HN lhang(3,icont)=0_ip
             ! HN lhang(4,icont)=0_ip
             ! HN lhang(5,icont)=0_ip

           endif

        else

           lmark(icont)=0_ip
          ! HN lhang(0,icont)=0_ip
          ! HN lhang(1,icont)=0_ip
          ! HN lhang(2,icont)=0_ip
          ! HN lhang(3,icont)=0_ip
          ! HN lhang(4,icont)=0_ip
          ! HN lhang(5,icont)=0_ip

        endif
     endif
  enddo
  !     
  !     Renumber the cells
  !
  icont=0_ip
  do icell=1,ncell
     if(lmarkc(icell)==1)then
        icont=icont+1_ip
        lmarkc(icell)=icont
     endif
  enddo
  !
  !     Compact the lcell array
  !

  ncont=0_ip
  ncell0=ncell
  do icell=1,ncell0
     if(lmarkc(icell)>0)then
        ncont=ncont+1
        lcell(ncont)%marked    = lcell(icell)%marked
        lcell(ncont)%level     = lcell(icell)%level
        lcell(ncont)%coor(1,1) = lcell(icell)%coor(1,1)
        lcell(ncont)%coor(2,1) = lcell(icell)%coor(2,1)
        lcell(ncont)%coor(3,1) = lcell(icell)%coor(3,1)
        lcell(ncont)%coor(1,2) = lcell(icell)%coor(1,2)
        lcell(ncont)%coor(2,2) = lcell(icell)%coor(2,2)
        lcell(ncont)%coor(3,2) = lcell(icell)%coor(3,2)
        lcell(ncont)%rsize     = lcell(icell)%rsize
        !
        !     Renumber the neighbors
        !
        ineigh=lcell(icell)%neigh(1)
        if(ineigh==0)then
           lcell(ncont)%neigh(1)=0 
        else 
           lcell(ncont)%neigh(1)=lmarkc(ineigh) 
        endif

        ineigh=lcell(icell)%neigh(2)
        if(ineigh==0)then
           lcell(ncont)%neigh(2)=0 
        else 
           lcell(ncont)%neigh(2)=lmarkc(ineigh)
        endif

        ineigh=lcell(icell)%neigh(3)
        if(ineigh==0)then
           lcell(ncont)%neigh(3)=0 
        else 
           lcell(ncont)%neigh(3)=lmarkc(ineigh)
        endif

        ineigh=lcell(icell)%neigh(4)
        if(ineigh==0)then
           lcell(ncont)%neigh(4)=0 
        else 
           lcell(ncont)%neigh(4)=lmarkc(ineigh)
        endif

        ineigh=lcell(icell)%neigh(5)
        if(ineigh==0)then
           lcell(ncont)%neigh(5)=0 
        else 
           lcell(ncont)%neigh(5)=lmarkc(ineigh)
        endif

        ineigh=lcell(icell)%neigh(6)
        if(ineigh==0)then
           lcell(ncont)%neigh(6)=0 
        else 
           lcell(ncont)%neigh(6)=lmarkc(ineigh)
        endif
        !
        !     Compact renum
        !
        iplace=8*(ncont-1_ip)+1_ip
        jplace=8*(icell-1_ip)+1_ip
        renum(iplace)=renum2(renum(jplace))
        iplace=iplace+1_ip
        jplace=jplace+1_ip 
        renum(iplace)=renum2(renum(jplace))
        iplace=iplace+1_ip
        jplace=jplace+1_ip 
        renum(iplace)=renum2(renum(jplace))
        iplace=iplace+1_ip
        jplace=jplace+1_ip 
        renum(iplace)=renum2(renum(jplace))
        iplace=iplace+1_ip
        jplace=jplace+1_ip 
        renum(iplace)=renum2(renum(jplace))
        iplace=iplace+1_ip
        jplace=jplace+1_ip 
        renum(iplace)=renum2(renum(jplace))
        iplace=iplace+1_ip
        jplace=jplace+1_ip 
        renum(iplace)=renum2(renum(jplace))
        iplace=iplace+1_ip
        jplace=jplace+1_ip 
        renum(iplace)=renum2(renum(jplace))
     endif
  enddo
  ncell=ncont 
  !
  !     Renumber the faces 
  ! 
  do iface=1,nface 
     do inofa=1,nnofam1
        lface(inofa,iface)=renum2(lface(inofa,iface)) 
     enddo
     lface(9,iface)=lmarkc(lface(9,iface))  
  enddo  
  !
  !     Renumber boundary points
  !
  do iboup=1,nboup
     lboup(iboup)=renum2(lboup(iboup))
  enddo
   call extbcs()
  !
  !     Transfer to hexahedra
  !
  nelem = ncell
  npoin = npnew
  mnode = 8
  mnodb = 4
  nboun = nface 
  !ncodn = 1_ip
  !ncodb = 6_ip
  call domain_memory_allocate('KFL_CODNO')
  call domain_memory_allocate('KFL_CODBO')

  icont=nhang
  nhang=0
  call domain_memory_allocate('GEOMETRY')
  nhang=icont
  do i=1,nelem
     ltype(i)=HEX08
  end do
  do i=1,nboun
     ltypb(i)=QUA04
  end do

  icont=0_ip
  do i=1,ncell

     icont=icont+1_ip
     ipoin=renum(icont)
     coord(1,ipoin)=lcell(i)%coor(1,1) 
     coord(2,ipoin)=lcell(i)%coor(2,1) 
     coord(3,ipoin)=lcell(i)%coor(3,1) 

     icont=icont+1_ip
     ipoin=renum(icont)
     coord(1,ipoin)=lcell(i)%coor(1,2) 
     coord(2,ipoin)=lcell(i)%coor(2,1) 
     coord(3,ipoin)=lcell(i)%coor(3,1) 

     icont=icont+1_ip
     ipoin=renum(icont)
     coord(1,ipoin)=lcell(i)%coor(1,2) 
     coord(2,ipoin)=lcell(i)%coor(2,2) 
     coord(3,ipoin)=lcell(i)%coor(3,1) 

     icont=icont+1_ip
     ipoin=renum(icont)
     coord(1,ipoin)=lcell(i)%coor(1,1) 
     coord(2,ipoin)=lcell(i)%coor(2,2) 
     coord(3,ipoin)=lcell(i)%coor(3,1) 

     icont=icont+1_ip
     ipoin=renum(icont)
     coord(1,ipoin)=lcell(i)%coor(1,1) 
     coord(2,ipoin)=lcell(i)%coor(2,1) 
     coord(3,ipoin)=lcell(i)%coor(3,2) 

     icont=icont+1_ip
     ipoin=renum(icont)
     coord(1,ipoin)=lcell(i)%coor(1,2) 
     coord(2,ipoin)=lcell(i)%coor(2,1) 
     coord(3,ipoin)=lcell(i)%coor(3,2) 

     icont=icont+1_ip
     ipoin=renum(icont)
     coord(1,ipoin)=lcell(i)%coor(1,2) 
     coord(2,ipoin)=lcell(i)%coor(2,2) 
     coord(3,ipoin)=lcell(i)%coor(3,2) 

     icont=icont+1_ip
     ipoin=renum(icont)
     coord(1,ipoin)=lcell(i)%coor(1,1) 
     coord(2,ipoin)=lcell(i)%coor(2,2) 
     coord(3,ipoin)=lcell(i)%coor(3,2) 

  enddo
  !
  !     Get the connectivities
  !
  icont=0
  icontpp=0
  do i=1,ncell 
     icont=icont+1
     icontpp=icontpp+1 
     icont1=renum(icontpp)
     icontpp=icontpp+1 
     icont2=renum(icontpp)  
     icontpp=icontpp+1 
     icont3=renum(icontpp)  
     icontpp=icontpp+1 
     icont4=renum(icontpp)  
     icontpp=icontpp+1 
     icont5=renum(icontpp)  
     icontpp=icontpp+1 
     icont6=renum(icontpp)  
     icontpp=icontpp+1 
     icont7=renum(icontpp)  
     icontpp=icontpp+1 
     icont8=renum(icontpp)  
     !index1=-lcell(i)%marked
     !index1=lcell(i)%level 
     lnods(1,icont)=icont1
     lnods(2,icont)=icont2
     lnods(3,icont)=icont3
     lnods(4,icont)=icont4
     lnods(5,icont)=icont5
     lnods(6,icont)=icont6
     lnods(7,icont)=icont7
     lnods(8,icont)=icont8
  enddo
  !
  ! Compact! HN LHANG
  !

  nhang=0
  do ipoin=1,npoin
     ! HN if(lhang(0,ipoin)/=0) then
     if(0/=0) then
        nhang=nhang+1
       ! HN lhang(0,nhang)=lhang(0,ipoin)+1
       ! HN lhang(5,nhang)=lhang(4,ipoin)
       ! HN lhang(4,nhang)=lhang(3,ipoin)
       ! HN lhang(3,nhang)=lhang(2,ipoin)
       ! HN lhang(2,nhang)=lhang(1,ipoin)
       ! HN lhang(1,nhang)=ipoin
     end if
  end do

  ! OJO: RETRECIR! HN LHANG
  
  !
  !     Transfer the faces
  !
  do iface=1,nface 
     lnodb(1,iface)=lface(1,iface)
     lnodb(2,iface)=lface(2,iface)
     lnodb(3,iface)=lface(3,iface)
     lnodb(4,iface)=lface(4,iface)
     lboel(1,iface)=lface(5,iface)
     lboel(2,iface)=lface(6,iface)
     lboel(3,iface)=lface(7,iface)
     lboel(4,iface)=lface(8,iface)
     lboel(5,iface)=lface(9,iface)
     !kfl_codbo(iface)=lface(10,iface)
  enddo
  !ncodn=1
  !call extbcs()
  !
  !     Transfer the boundary points
  !
  do iboup=1,nboup
    ! kfl_codno(1,lboup(iboup))=10_ip
  enddo
  !call coubcs(1_ip)
  mcono = 3
  !
  !     Write the marked points
  !
  !do ipoin=1,npnew
  !   if(lmark(ipoin)==1)then
  !      write(60,400)ipoin,lhang(1,ipoin),lhang(2,ipoin),lhang(3,ipoin),lhang(4,ipoin)
  !   endif
  !enddo
  !close(60)

  !
  !     Output the size
  !

  !open(unit=70,file='cartGid.res',status='unknown')
  !rewind 70

  !allocate(rsize(npnew),stat=istat)
  !call memchk(zero,istat,memor_msh,'RSIZE','writeGidconf',rsize)
  !allocate(lsize(npnew),stat=istat)
  !call memchk(zero,istat,memor_msh,'LSIZE','writeGidconf',lsize)

  !do icell=1,ncell
  !   do ipo=1,8
  !      ipoin=renum(8*(icell-1)+ipo)
  !      rsize(ipoin)=rsize(ipoin)+lcell(icell)%rsize
  !      lsize(ipoin)=lsize(ipoin)+1
  !   enddo
  !enddo

  !do  ipoin=1,npnew
  !   if(lsize(ipoin)==0)then
  !      write(*,*)'Error in writeGidconf, lsize=0 at ipoin:',ipoin
  !      stop 
  !   endif
  !   rsize(ipoin)=rsize(ipoin)/float(lsize(ipoin))
  !enddo


  !write(70,10)
  !write(70,11)
  !write(70,13)
  !do  ipoin=1,npnew
  !   write(70,500)ipoin,rsize(ipoin)
  !enddo
  !write(70,14)
  !write(70,15)
  !close(70)

  call memchk(2_ip,istat,memor_msh,'LBOUP','writeGidconf',lboup)
  deallocate(lboup,stat=istat)
  if(istat/=0) call memerr(2_ip,'LBOUP','writeGidconf',0_ip)
  !call memchk(2_ip,istat,memor_msh,'LSIZE','writeGidconf',lsize)
  !deallocate(lsize,stat=istat)
  !if(istat/=0) call memerr(2_ip,'LSIZE','writeGidconf',0_ip)
  !call memchk(2_ip,istat,memor_msh,'RSIZE','writeGidconf',rsize)
  !deallocate(rsize,stat=istat)
  !if(istat/=0) call memerr(2_ip,'RSIZE','writeGidconf',0_ip)
  call memchk(2_ip,istat,memor_msh,'LMARKC','writeGidconf',lmarkc)
  deallocate(lmarkc,stat=istat)
  if(istat/=0) call memerr(2_ip,'LMARKC','writeGidconf',0_ip)
  call memchk(2_ip,istat,memor_msh,'LMARK','writeGidconf',lmark)
  deallocate(lmark,stat=istat)
  if(istat/=0) call memerr(2_ip,'LMARK','writeGidconf',0_ip)
  call memchk(2_ip,istat,memor_msh,'RENUM','writeGidconf',renum)
  deallocate(renum,stat=istat)
  if(istat/=0) call memerr(2_ip,'RENUM','writeGidconf',0_ip)
  call memchk(2_ip,istat,memor_msh,'RENUM2','writeGidconf',renum2)
  deallocate(renum2,stat=istat)
  if(istat/=0) call memerr(2_ip,'RENUM2','writeGidconf',0_ip)

end subroutine cargeo

subroutine chkarraycart(lcell,ncell,renum)
 use def_kintyp, only          :  ip,rp,lg,cell
  implicit none
  integer(ip),intent(in)     :: ncell
  type(cell),intent(in)      :: lcell(ncell)
  integer(ip),intent(in)     :: renum(ncell*8)



end subroutine chkarraycart
