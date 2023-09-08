subroutine ada_sahole
!-----------------------------------------------------------------------
!****f* adapti/ada_sahole
! NAME 
!    ada_kpatch
! DESCRIPTION
!    This routine "sausages" a hole by eliminating nodes and elements
!    and by creating the hole surface mesh boundary (HSMB). this will
!    be filled by a patch mesh
!
!    It creates the following vectors:
!    
!    lolne_ada  -->  new numbering for the HSMB:
!       lolne_ada(1,ipoin) = jpoin
!       ipoin <-- original global numbering           
!       jpoin <-- new global numbering  
!    finally, lolne_ada is used as a permutation vector to move to the          
!    first places the HSMB nodes
!       
!    lnosu_ada  -->  boundary connectivity for the HSMB, including 
!                    the new points of the immersed boundary
!    conew_ada  -->  new points, using the new numbering     
!             
! 
! USES
!
! USED BY
!
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_adapti
  implicit none
  integer(ip)   :: &
       idime,ielem,iimmo,kneli,ifnei,iface,jface,jelem,kelem,ipoin,kpoin,jpoin,pelty,&
       pnode,inode,jnode,knode,nenew,npnew,ipoii,ielei,ielsu,jelsu,kelsu,inau1,inau2,&
       inodb,inau3,nptmp,netmp,ienum,ichck,nfacn,iffil,ielad,iolne,jolne,ifono,kfono
  integer(ip)   :: &
       lnofa(mnode),lfaux(nnodb_ada+1,mnode)
  real(rp)      :: &
       sauxi, saux2
  !
  ! Eliminate marked elements and its kneli neighbors
  !
  
  !
  !  1. mark first neighbors using element-point graph
  !
  kneli = 1

  do iimmo = 1,nimmo_ada            
     ielem = lelim_ada(1,iimmo)
     ifnei= 1
     do while (ifnei == 1)
        pelty=ltype(ielem)
        pnode=nnode(pelty)
        inode= 1
        do while (inode < pnode)
           ipoin= lnods(inode,ielem)
           do ienum= pelpo(ipoin), pelpo(ipoin+1)-1
              jelem= lelpo(ienum)
              lenew_ada(jelem)%neseg = -1                         
           end do
           inode = inode + 1
        end do
        kelem = lenew_ada(ielem)%knext         
        if (kelem == -1) then
           ifnei = 0           
           inode = pnode
        else
           ielem = kelem
        end if
     end do
  end do
  

  !
  !  2. mark kneli next neighbors (only when kneli > 1) 
  !
  
  !            to be done
  
  !
  !  3. extract the nodes that will be used to generate the new part of the mesh
  !
  
  do ielem = 1,nelem
     if (lenew_ada(ielem)%neseg == -1) then
        do iface= pelel(ielem) , pelel(ielem+1)-1             
           jelem= lelel(iface)                            
           if (lenew_ada(jelem)%neseg == 0) then
              lenew_ada(ielem)%neseg = -2            ! this ielem has a hole boundary
           end if
        end do                
     end if
     lenew_ada(ielem)%kfaie = 0
  end do  

  call ada_memall(2)
  
  !
  !     3.1 get the hole surface mesh boundary (HSMB) in lnosu_ada
  !
  kpoin= 0
  ipoin= 0
  ielsu= 0
  lnofa= 0
  lnosu_ada= 0
  lolne_ada= 0
  do ielem = 1,nelem
     if (lenew_ada(ielem)%neseg == -2) then
        do iface= pelel(ielem) , pelel(ielem+1)-1             
           nfacn= 1 + (pelel(ielem+1)- 1 - pelel(ielem)) 
           jelem= lelel(iface)                            
           if (lenew_ada(jelem)%neseg == 0 .or. lenew_ada(jelem)%neseg == -3) then
              do jface= pelel(jelem) , pelel(jelem+1)-1             
                 kelem= lelel(jface)
                 if (lenew_ada(kelem)%neseg == 0 .or. lenew_ada(kelem)%neseg == -3) then
                    inau1 = 0
                    do jnode= 1,nnode(ltype(jelem))                       
                       jpoin= lnods(jnode,jelem)
                       do inau2= 1,nnode(ltype(kelem))
                          kpoin= lnods(inau2,kelem)
                          if (jpoin == kpoin) then
                             inau1 = inau1+1
                             lnofa(inau1) = jpoin
                          end if
                       end do
                    end do
                    if (ndime == 3) then
                       !
                       ! acaaaaaaaaa by the moment only TRIANGLES and TETRAS!!!!!                       
                       !                       
                    end if                    
                    if (inau1 == nnodb_ada) then
                       ielsu= ielsu+1                       
                       lnosu_ada(1:nnodb_ada,ielsu) = lnofa(1:nnodb_ada)
                       do inau2= 1,nnodb_ada
                          kpoin= lnofa(inau2) 
                          if (lolne_ada(1,kpoin) == 0) then
                             ipoin = ipoin + 1
                             lolne_ada(1,kpoin)= ipoin                              
                          end if
                       end do
                       lenew_ada(jelem)%neseg = -3
                       lnofa = 0
                    end if
                 end if
              end do
           end if
        end do
     end if
  end do
  nptmp = ipoin
  netmp = ielsu
  
  !
  !  3.2 add those HSMB elements lying in the
  !      domain boundary (when needed, for instance when the immersed body cuts
  !      the domain boundary)
  !
!  kpoin= nptmp                             ! preparing to add new points
  ielsu= netmp                             ! preparing to add new boundary elements
  
  do ielem = 1,nelem
     if (lenew_ada(ielem)%neseg < 0) then
        nenew = nenew - 1
        nfacn= 1 + (pelel(ielem+1)- 1 - pelel(ielem)) 
        if (nfacn <= 2) then                   ! this element has a boundary face
           do inode= 1,nnode(ltype(ielem))
              jnode= mod(inode+1,nnode(ltype(ielem)))
              knode= mod(inode+2,nnode(ltype(ielem)))
              lfaux(1,inode)= lnods(inode,ielem)
              lfaux(2,inode)= lnods(jnode,ielem)
              if (ndime == 3) then
                 lfaux(3,inode)= lnods(knode,ielem)       ! recall: only for tetras in 3D!!!!
              end if              
              lfaux(nnodb_ada+1,inode) = 1
           end do

           inau2= 0
           do iface= pelel(ielem) , pelel(ielem+1)-1             
              inau2= inau2+1
              jelem= lelel(iface)                            
              inau1= 0
              do jnode= 1,nnode(ltype(jelem))
                 jpoin= lnods(inode,jelem)
                 do inode= 1,nnodb_ada
                    ipoin= lfaux(inode,inau2)
                    if (ipoin == jpoin) inau1=inau1+1
                 end do
              end do
              if (inau1 < nnodb_ada) lfaux(nnodb_ada+1,inau2) = 0
           end do

           do iface= 1,nnode(ltype(ielem))                  ! found
              if (lfaux(nnodb_ada+1,iface)==0) then
                 ielsu = ielsu + 1
                 do inau2= 1, nnodb_ada
                    lnosu_ada(inau2,ielsu) = lfaux(inau2,iface)
                 end do
              end if
           end do
           
        end if
     end if
  end do          
  netmp = ielsu                             ! updating the total number of hsmb elements

  !
  !     3.3.1 check that the HSMB is nice: repeated faces and non-convexity
  !  
  do ichck= 1,2
     do ielsu= 1,netmp
        ifnei= 0
        jelsu= 0
        do inode= 1,nnodb_ada
           lnofa(inode)= abs(lnosu_ada(inode,ielsu))
        end do
        inau1= 0
        do while (ifnei == 0)
           jelsu = jelsu + 1
           if (jelsu /= ielsu) then
              do inode= 1,nnodb_ada
                 do jnode= 1,nnodb_ada
                    if (abs(lnosu_ada(jnode,jelsu)) == abs(lnofa(inode))) then
                       if (lnofa(inode) > 0) then
                          lnofa(inode) = lnofa(inode) * (2-ichck)      ! set to 0 in the second pass
                          inau1=inau1 + 1
                       end if
                    end if
                 end do
              end do
           end if
           if (inau1 == nnodb_ada) then
              ifnei = 1 
              do inode= 1,nnodb_ada
                 do jnode= 1,nnodb_ada
                    if (abs(lnosu_ada(jnode,jelsu)) == abs(lnofa(inode))) then
                       if (lnofa(inode) > 0) inau1=inau1 + 1
                    end if
                 end do
              end do
              if (inau1 == 2*nnodb_ada) then
                 ! label repeated faces with a 0
                 lnosu_ada(1:nnodb_ada,jelsu) = 0          
              end if
           else if (jelsu == netmp) then
              ifnei = 2
           end if
        end do
        if (ifnei == 2 .and. ichck==2) then
           ! label non-convex faces with negative values (because 0 = -0)
           lnosu_ada(1:nnodb_ada,ielsu) = -lnosu_ada(1:nnodb_ada,ielsu)  
        end if
     end do
  end do
  !
  !     3.3.2 eliminate repeated faces and check non-convexity
  !  
  ielsu = 1
  do while (ielsu < netmp)
     if (lnosu_ada(1,ielsu) <= 0) then
        do jelsu = ielsu,netmp
           lnosu_ada(1:nnodb_ada,jelsu) = lnosu_ada(1:nnodb_ada,jelsu+1)  
        end do
        netmp = netmp - 1
        ielsu = ielsu - 1
     end if
     ielsu = ielsu + 1
  end do  
  nelsu_ada = ielsu 

  !
  !    3.4 store the nodes on the HSMB, now passing all to the new numbering
  !  

  ipoin= 0
  lolne_ada= 0
  do ielsu = 1,nelsu_ada
     do inau1= 1,nnodb_ada
        kpoin= lnosu_ada(inau1,ielsu)
        if (lolne_ada(1,kpoin) == 0) then
           ipoin = ipoin + 1
           lolne_ada(1,kpoin)= ipoin                              
        end if
     end do
  end do
  nposu_ada = ipoin

  jolne= 0
  do ielsu = 1,nelsu_ada
     do inau1= 1,nnodb_ada
        jpoin= lnosu_ada(inau1,ielsu)
        kpoin= lolne_ada(1,jpoin)
        !
        ! lolne_ada(3,*) and lolne_ada(4,*) keep the correspondence between new and original meshes
        ! which is used just for postprocess
        !
        if (lolne_ada(4,kpoin)==0) then
           jolne= jolne+1
           lolne_ada(3,jpoin) = kpoin
           lolne_ada(4,kpoin) = jpoin
        end if
        conew_ada(1:ndime,kpoin) = coord(1:ndime,jpoin)          
        lnosu_ada(inau1,ielsu) = kpoin
     end do
  end do

  !
  !     3.5 add the immersed object's boundary elements and its nodes
  !         and compute reference sizes
  !  
  ielsu    = nelsu_ada
  ipoin    = nposu_ada
  nelpa_ada= nelsu_ada
  npopa_ada= nposu_ada
  sizem_ada= 0.0_rp
  do iimmo = 1,nimmo_ada                   ! add the immersed object's boundary elements
     do ielei= 1,nelei_ada(iimmo)
        ielsu= ielsu + 1
        do inau1= 1,nnodb_ada
           lnosu_ada(inau1,ielsu) = lnodi_ada(inau1,ielei,iimmo) + nposu_ada
           ipoii = lnodi_ada(inau1,ielei,iimmo) 
           kpoin = lnosu_ada(inau1,ielsu)  ! and its nodes
           conew_ada(1:ndime,kpoin) = coori_ada(1:ndime,ipoii,iimmo)        
        end do

        sauxi= 0.0_rp
        saux2= 0.0_rp
        if (ndime == 2) then
           !
           ! ONLY 2D IS PROGRAMMED!!!!!!!!!
           !
           ipoin = lnodi_ada(1,ielei,iimmo)
           jpoin = lnodi_ada(2,ielei,iimmo)
           do idime= 1,ndime
              sauxi = (coori_ada(idime,ipoin,iimmo) - coori_ada(idime,jpoin,iimmo))
              saux2 = saux2 + sauxi*sauxi
           end do
           sauxi = sqrt(saux2)
           if (ielei == 1) then
              sizem_ada(1,iimmo) = sauxi
              if (iimmo == 1) sizet_ada = sauxi                 
           else
              if (sauxi < sizem_ada(1,iimmo)) sizem_ada(1,iimmo)= sauxi              
              if (sauxi < sizet_ada         ) sizet_ada         = sauxi              
           end if
        else
           !
           ! ONLY 2D IS PROGRAMMED!!!!!!!!!
           !
           
        end if

     end do
     nposu_ada = nposu_ada + npoii_ada(iimmo)
  end do
  nelsu_ada = ielsu


  !
  !  4. reorganize lnosu_ada in lnoad_ada
  !
  !
  call ada_memall(3)

  if (nnodb_ada == 2) then
     !
     ! ONLY 2D IS PROGRAMMED!!!!!!!!!
     !

     lnoad_ada = 0     
     do inodb= 1,nnodb_ada
        lnoad_ada(inodb,1) = lnosu_ada(inodb,1)
     end do
     lnosu_ada(nnodb_ada+1,1) = 0
     
     do ielad = 2,nelpa_ada
        ielsu = 2
        ipoin = lnoad_ada(2,ielad-1)
        ifono = 0
        do while (ifono==0) 
           if (lnosu_ada(nnodb_ada+1,ielsu) == 0) then
              inode = 0
              kfono = 0
              do while (kfono==0) 
                 inode= inode+1
                 if (lnosu_ada(inode,ielsu) == ipoin) then
                    ifono= 1
                    lnosu_ada(3,ielsu)= 1
                    if (inode==1) then
                       lnoad_ada(1,ielad)= lnosu_ada(1,ielsu)
                       lnoad_ada(2,ielad)= lnosu_ada(2,ielsu)
                    else
                       lnoad_ada(1,ielad)= lnosu_ada(2,ielsu)
                       lnoad_ada(2,ielad)= lnosu_ada(1,ielsu)
                    end if
                 end if
                 if (inode == nnodb_ada) kfono=1
              end do
           end if
           ielsu = ielsu+1
        end do
     end do
  else
     !
     ! ONLY 2D IS PROGRAMMED!!!!!!!!!
     !
  end if


  !
  !  5. permute node numbering to place HSMB first and discarded nodes last
  !
  !
  
  lolne_ada(1:2,1:npoin)= 0  ! reset lolne_ada, except for columns 3 and 4
  
  ! discard sausaged nodes
  do ielem = 1,nelem
     if (lenew_ada(ielem)%neseg < 0) then
        do inode = 1,mnode
           ipoin = lnods(inode,ielem)
           lolne_ada(1,ipoin)= -1
           if (lolne_ada(3,ipoin)==0) lolne_ada(3,ipoin)= -1
        end do
     end if
  end do

  ! include complement-sausage nodes
  nelco_ada = 0         ! complement-sausage number of elements
  jolne= npopa_ada
  do ielem = 1,nelem
     if (lenew_ada(ielem)%neseg == 0) then
        nelco_ada = nelco_ada + 1         ! complement-sausage number of elements
        do inode = 1,mnode
           ipoin = lnods(inode,ielem)
           if (lolne_ada(3,ipoin)==0) then              
              jolne= jolne + 1
              lolne_ada(3,ipoin)=jolne
              lolne_ada(4,jolne)=ipoin
           end if
       end do
     end if
  end do
  npoco_ada = jolne         ! complement-sausage number of nodes
 
  ! HSMB nodes:
  iolne= 0
  do ielad= 1,nelpa_ada
     do inode = 1,nnodb_ada
        ipoin = lnoad_ada(inode,ielad)
        if (lolne_ada(1,ipoin)<=0) then
           iolne= iolne + 1
           lolne_ada(1,ipoin)=iolne        ! lolne_ada(1,*) --> old to new
           lolne_ada(2,iolne)=ipoin        ! lolne_ada(2,*) --> new to old
        end if
     end do
  end do

  ! interior nodes:
  jolne= iolne
  do ipoin= 1,npoin
     if (lolne_ada(1,ipoin)==0) then
        iolne= iolne + 1
        lolne_ada(1,ipoin)=iolne
        lolne_ada(2,iolne)=ipoin
     end if
  end do


!  do ipoin=1,npoin
!     write(6,*) ipoin, lolne_ada(1,ipoin), lolne_ada(2,ipoin), lolne_ada(3,ipoin), lolne_ada(4,ipoin)
!  end do


end subroutine ada_sahole
