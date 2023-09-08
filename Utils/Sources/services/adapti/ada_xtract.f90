subroutine ada_xtract
!-----------------------------------------------------------------------
!****f* adapti/ada_segmsh
! NAME 
!    ada_modmsh
! DESCRIPTION
!    This routine segmentates the background mesh
! USES
!
! USED BY
!
!***
!-----------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_elmtyp
  use      mod_memchk

  use      def_adapti

  implicit none

  integer(ip)   :: &
       ipoii,ipoin,jpoin,kpoin,iimmo,ielei,nelei,inoii,ipoi1,ipoi2,ipoi3,istat,ieseg,&
       pelty,pnode,pface,inode,jnode,knode,ifvax,ifcut,nonew,ionew,nenew,ienew,&
       npnew,ipnew,ielem,jelem,kelem,iele1,ienum,&
       iface,jface,kface,ifnei,ifoun,idime,ifiim,&
       inau1,inau2,knau2,iveri,ivere,kpara,ksety,ifloo,idble,ifuno,kneli
  integer(ip)   :: &
       lnoie(nnodi_ada),lelie(nnodi_ada),lnofa(mnode),&
       loloc(3,10),lnaux(10,10),ifacu(20)
  real(rp)      :: &
       asegi,bsegi,asege,bsege,xinte,yinte,zinte,xdeni,xvali,dimod,semod

  real(rp)      :: &
       elcoi(3,nnodi_ada), elcor(3,10), coloc(3,10),coaux(3,10),xangi(2),xange(4,10),&
       basin(3,3),xscal(2,10)
  

  loloc= 0
  coloc(1:ndime,1:mnode)= 0.0_rp
  do ielem=1,nelem
     lnseg_ada(ielem)       = 0
     lenew_ada(ielem)%neseg = 0
     lenew_ada(ielem)%kfaie = 0
     lenew_ada(ielem)%knext = 0
  end do

  do iimmo = 1,nimmo_ada          
     nelei= nelei_ada(iimmo)
     ielei= 1
     ifiim= 0
     ifnei= 1
     ifacu= 0
     do while (ielei <= nelei)
        lnoie(1:nnodi_ada) = lnodi_ada(1:nnodi_ada,ielei,iimmo)        ! node number in the patch
        !
        ! Get the coordinates of the immersed (the patch) element nodes
        !
        basin= 0.0_rp
        do inoii=1,nnodi_ada
           ipoii               = lnoie(inoii)
           lelie(inoii)        = lelim_ada(ipoii,iimmo)
           elcoi(1:ndime,inoii)= coori_ada(1:ndime,ipoii,iimmo)           
        end do

        ifloo = 1
        if (nnodi_ada == 2) then
           if (lelie(1) == lelie(2)) then
              ielei = ielei + 1
              ifloo = 0
           end if
           ! compute the orientation vector of the immersed edge
           basin(1,1)= coori_ada(1,lnoie(2),iimmo) -  coori_ada(1,lnoie(1),iimmo)
           basin(1,2)= coori_ada(2,lnoie(2),iimmo) -  coori_ada(2,lnoie(1),iimmo)
           xdeni     = basin(1,1) * basin(1,1) + basin(1,2) * basin(1,2)
           xdeni     = sqrt(xdeni)
           basin(1,1)= basin(1,1) / xdeni
           basin(1,2)= basin(1,2) / xdeni
        else if (nnodi_ada == 3) then
           ! ---------------------------------
           ! 3D --> to be programmed...
           ! ---------------------------------
           !
        end if
        !
        ! Loop over element faces:
        !        
        ifcut = 0
        if (nnodi_ada == 2 .and. ifloo == 1) then
           ! ---------------------------------
           ! 2D --> faces are edges
           ! ---------------------------------
           !
           
           if (ifnei == 1) then
              !  Use the first edge node ielem's element as reference
              ielem= lelie(1)
              kface= 0
           else if (ifnei == 0) then
              !  Use the previously found shared cut edge 
              !  between two neighboring elements 
              ielem= kelem              
           end if
           
           pelty=ltype(ielem)
           pnode=nnode(pelty)
           pface=nface(pelty)
              
           !  Compute the patch element (edge) equation
           !  set iveri = 1 when the edge is vertical 
           xdeni= elcoi(1,1) - elcoi(1,2)
           iveri= 0
           if (xdeni <  0.0_rp) xdeni = - xdeni
           asegi = 0.0_rp
           if (xdeni > 1.0e-10) then
              asegi= (elcoi(2,1) - elcoi(2,2)) / (elcoi(1,1) - elcoi(1,2))
              bsegi= elcoi(2,2) - asegi * elcoi(1,2)              
              xangi(1)= sqrt(1.0_rp + asegi*asegi)
              xangi(1)= 1.0_rp / xangi(1)                                   ! cos(i)
              xangi(2)= asegi  * xangi(1)                                   ! sin(e)
           else
              iveri= 1
              xangi(1)= 1.0_rp
              xangi(2)= 0.0_rp
           end if

    
           !  Gather the background nodes coordinates
           do inode= 1,pnode
              ipoin= lnods(inode,ielem)
              do idime=1,ndime
                 elcor(idime,inode) = coord(idime,ipoin)  
              end do
           end do

           inode = 1
           jnode = 1
           knode = 2
           ifoun = 0
           !
           !  Check the element segments looking for an intersection
           !
           do while (ifoun == 0)
           
              jnode = jnode + 1
              knode = knode + 1
              if (jnode == pnode+1) jnode = 1
              if (knode == pnode+1) knode = 1

              idble= 0
              ipoin = lnods(inode,ielem)
              jpoin = lnods(jnode,ielem)
              if (ifacu(1) == ipoin .or. ifacu(1) == jpoin) then
                 if (ifacu(2) == ipoin .or. ifacu(2) == jpoin) then
                    !
                    ! When ifnei was 0, check that this face is the same where
                    ! the last cut was found. 
                    !
                    idble= 1
                 end if
              end if
              
              !  For each of the edges of the background element, 
              !  compute the edge equation
              !  set ivere = 1 when the edge is vertical 

              xdeni = elcor(1,inode) - elcor(1,jnode) 
              ivere = 0
              asege = 0.0_rp
              if (xdeni <  0.0_rp) xdeni = - xdeni
              if (xdeni > 1.0e-10) then
                 asege = (elcor(2,inode) - elcor(2,jnode)) / (elcor(1,inode) - elcor(1,jnode))
                 bsege = elcor(2,inode) - asege * elcor(1,inode)              
              else
                 ivere = 1
              end if

              !  Compute and verify the intersection
              !  when one of the edges is vertical, set xinte to its x-coordinate 
              !  Set kpara = 1 when both edges are co-linear              
              kpara = 0
              if (iveri+ivere == 2) then       ! both are vertical, yet co-linear
                 kpara = 1
              else if (iveri == 1) then        ! immersed is vertical
                 xinte = elcoi(1,inode)
              else if (ivere == 1) then        ! background is vertical
                 xinte = elcor(1,inode)
              else                             ! no one is vertical, then compute xinte
                 xdeni = asegi - asege
                 if (xdeni <  0.0_rp) xdeni = - xdeni
                 if (xdeni > 1.0e-10) then
                    xinte = (bsege - bsegi) / (asegi - asege)
                 else                          ! both are co-linear
                    kpara = 1
                 end if
              end if

              ifvax = 0

              if (kpara == 0) then             ! check only when they are not co-linear    
                 !    Check x 
                 if (ivere == 0) then
                    ! non-vertical edge: range for xinte is set by background edge
                    if (elcor(1,inode) >= elcor(1,jnode)) then
                       if (xinte <= elcor(1,inode) ) then
                          if (xinte >= elcor(1,jnode) ) then
                             ifvax = 1                       
                          end if
                       end if
                    else
                       if (xinte <= elcor(1,jnode) ) then
                          if (xinte >= elcor(1,inode) ) then
                             ifvax = 1                       
                          end if
                       end if
                    end if
                 else
                    ! vertical edge: range for xinte is set by immersed edge
                    if (elcoi(1,1) >= elcoi(1,2)) then
                       if (xinte <= elcoi(1,1) ) then
                          if (xinte >= elcoi(1,2) ) then
                             ifvax = 1                       
                          end if
                       end if
                    else
                       if (xinte <= elcoi(1,2) ) then
                          if (xinte >= elcoi(1,1) ) then
                             ifvax = 1                       
                          end if
                       end if
                    end if                    
                 end if
                 !    Check y
                 if (ifvax == 1) then
                    yinte = asege * xinte + bsege
                    if (ivere == 1) yinte = asegi * xinte + bsegi
                    if (elcor(2,inode) >= elcor(2,jnode)) then
                       if (yinte <= elcor(2,inode) ) then
                          if (yinte >= elcor(2,jnode) ) then
                             ifvax = 2
                          end if
                       end if
                    else
                       if (yinte <= elcor(2,jnode) ) then
                          if (yinte >= elcor(2,inode) ) then
                             ifvax = 2                       
                          end if
                       end if
                    end if
                 end if                 
              end if              
              !
              ! Check that, if idble==2, intersection point is not the same as the last one found
              !              
              if ((idble == 1) .and. (ifvax == 2)) then
                 dimod = (xinte-coloc(1,2))*(xinte-coloc(1,2))
                 dimod = dimod + (yinte-coloc(2,2))*(yinte-coloc(2,2))
                 if (dimod < 1.0e-10) then
                    idble = 0
                    ifvax = 0                    
                 end if
              end if

              if ((ifvax == 0) .or. (ifvax == 1)) then
                 !
                 ! Intersection not found
                 !
                 inode = inode + 1              
                 if (inode == pnode + 1) ifoun = 2
                            
              else if (ifvax == 2) then
                 !
                 ! Check that the intersection point is not too close to a background node
                 !
                 ifuno = 0
                 semod = (elcor(1,inode)-elcor(1,jnode))*(elcor(1,inode)-elcor(1,jnode))
                 semod = semod + (elcor(2,inode)-elcor(2,jnode))*(elcor(2,inode)-elcor(2,jnode))
                 dimod = (xinte-elcor(1,inode))*(xinte-elcor(1,inode))
                 dimod = dimod + (yinte-elcor(2,inode))*(yinte-elcor(2,inode))
                 xdeni = dimod / semod
                 xdeni = sqrt(xdeni)
                 
                 if (xdeni < 0.3_rp) then
                    ifuno = inode
                    ifvax = 3
                 else                    
                    dimod = (xinte-elcor(1,jnode))*(xinte-elcor(1,jnode))
                    dimod = dimod + (yinte-elcor(2,jnode))*(yinte-elcor(2,jnode))
                    xdeni = dimod / semod
                    if (xdeni < 0.05_rp) then
                       ifuno = jnode
                       ifvax = 3
                    end if
                 end if                                  
              
                 if (ifuno == 0) then
                    !
                    ! Intersection ok, go ahead
                    !                    
                    coloc(1,1) = xinte     ! new node coordinate x
                    coloc(2,1) = yinte     ! new node coordinate y
                    ifcut = 2
                    if (ielei == 1 .and. ifiim == 0) then                       
                       !             
                       !     jnode              
                       !       |\               
                       !       | + \            
                       !       |  +   \         
                       !       |   +     \      
                       !       |    +       \      
                       !       |      @········X············@······      
                       !       |    +     +       \      
                       !       |  +             +    \   
                       !     knode ----------------  inode
                       
                       ipoin= lnods(knode,ielem)                                            
                       coloc(1,2) = elcoi(1,1)         ! immersed node coordinate
                       coloc(2,2) = elcoi(2,1)         ! immersed node coordinate
                       ifcut = 1                             
                       ifiim = 1                       
                    end if
                    
                 else if (ifuno /= 0) then
                    
                    !
                    ! Intersection found, but it is too close to ifuno node of the
                    ! background mesh. 
                    ! coord is modified to merge both points (pray for no major distortions...)
                    ! and ielem is divided in two triangles
                    ! 
                    ! 
                    ipoin = lnods(ifuno,ielem)
                    coloc(1,1) = xinte     
                    coloc(2,1) = yinte     
                    
                    ifoun = 1
                    ifcut = 4
                    kpoin = ipoin
                    
                 end if

                 ifoun = 1
                 if (idble == 1) ifcut = 3    ! one edge is cut twice                 
                 
              end if

           end do


           if (ifoun == 1) then
              
              ! use the ielem neighbors to check if
              ! any of them has the second node inside
              ! also, do nnodi_ada loops to identify the cut face
              
              ! when ifcut = 4, use the element-node graph, otherwise
              ! use element-element graph.

              ifnei= 0
              kelem= 0

              jface= 0
              kface= 0
              iface= pelel(ielem) - 1 
              do while (kface == 0)
                 iface= iface + 1
                 if (iface == (pelel(ielem+1)-1)) kface = -1
                 jface= jface+1
                 jelem= lelel(iface)                 
                 if (jelem == lelie(2)) then
                    ifnei= 1
                 end if
                 if (kelem == 0) then
                    inau1 = 0
                    do inau2= 1,nnode(ltype(jelem))
                       ipoin= lnods(inode,ielem)
                       jpoin= lnods(inau2,jelem)
                       if (ipoin == jpoin) then
                          inau1 = inau1+1
                       end if
                    end do
                    do inau2= 1,nnode(ltype(jelem))
                       ipoin= lnods(jnode,ielem)
                       jpoin= lnods(inau2,jelem)
                       if (ipoin == jpoin) then
                          inau1 = inau1+1
                       end if
                    end do
                    if (inau1 == nnodi_ada) then
                       lnofa        = 0
                       lnofa(inode) = 1
                       lnofa(jnode) = 1
                       kface    = jface
                       kelem    = jelem
                       ifacu(1) = lnods(inode,ielem)
                       ifacu(2) = lnods(jnode,ielem)                       
                       do inau2= 1,nnode(ltype(jelem))
                          if (lnofa(inau2) == 0) knau2 = inau2
                       end do
                    end if
                 end if
              end do

              if (kface == -1) kface = 0

              ! store the next element to be segmented (knext) and the
              ! segmented face (kfape)
              ! store also, for knext, the previous element

              lenew_ada(ielem)%kfaie = kface
              lenew_ada(ielem)%knext = kelem
              lenew_ada(kelem)%kprev = ielem
              lenew_ada(kelem)%knoie = knau2

              coloc(1,2) = coloc(1,1) 
              coloc(2,2) = coloc(2,1) 

              ! allocate the new points and elements in a temporal array
              !   4-node contact element:

              if (ifnei == 1) then
                 ielei = ielei + 1
                 ifnei = 0
                 if (ielei > nelei) then                     ! last ielei reach
                    ifvax = 0
                    ielei = nelei + 3                        ! to get out of the ielei loop
                    lenew_ada(ielem)%knext = kelem
                    lenew_ada(kelem)%knext = -1
                 end if
              end if

              ! however, if kelem == 0, the cut face is not share with any neighboring
              ! element. this means that we reach the border of the background mesh.
              ! also, reposition the last node at this face intersection
              
              if (kelem == 0) then
                 ifvax = 0
                 ielei = nelei + 3                        ! to get out of the ielei loop
                 lenew_ada(ielem)%knext = -1
                 coori_ada(1,lnoie(2),iimmo) = xinte
                 coori_ada(2,lnoie(2),iimmo) = yinte
              end if

              !    - assign the first loloc and coloc to the second for the next search
              
           
           else if (ifoun == 2) then

              !
              ! No intersection was found:
              ! This means that both edge nodes are inside the same element
              ! Go the next element but keeping the previous loloc and coloc
              !

              ielei = ielei + 1
              

           end if

        else if (nnodi_ada == 3) then
           !
           ! 3D --> faces are faces
           !        


        end if
     end do
  end do


!
! By the moment, only triangles (2D) and tetras (3D)
!
  nnode_ada = 4
  nnodb_ada = 3
  nnoco_ada = 6
  peltb_ada = PEN06
  pelte_ada = TET04
  if (ndime == 2) then
     nnode_ada = 3
     nnodb_ada = 2
     nnoco_ada = 4
     peltb_ada = BAR02
     pelte_ada = TRI03
  end if

end subroutine ada_xtract
