subroutine order_boundaries
  use def_parame
  use def_master
  use def_domain
  use def_elmtyp
  implicit none
  integer(ip)                        :: ielem,pelty,iboun,pblty,pface,ipoin
  integer(ip)                        :: inode,inodb,knode,kpoin,pnodb,jnodb,cont
  integer(ip)                        :: pnode,iface,inodf,jface,jnodf,knodf
  integer(ip)                        :: nncorner, nnedg, nnface ! Num. nodes on corners, edges and face of the elm. face.
  logical(lg)                        :: change_order,select_face

  ! nboun :: number of boundaries
  ! lnodb(:,iboun) :: list of nodes composing element face
  ! nface(pelty) :: Faces on element dep. on elem. type
  ! pelty :: Element num. ID
  ! lboel(:,iboun) :: list ? on boundary i
  ! lnodb(:,iboun) :: list of nodes composing boundary
  ! nnode(pelty) :: Nodes per element type
  ! ltype(:) :: List of element types
  ! nnodf(pelty) % l(iface) :: Nodes on this face for element type 'x'
  ! lnods(:,ielem) :: List of element nodes
  ! pnodb :: nodes on boundar
  ! lelbo(iboun) :: ?

  if(INOTMASTER)then
     ! Cycle through boundaries
     do iboun= 1,nboun

        pnodb = nnode(abs(ltypb(iboun))) ! This is how it's being used originally
        ! pnodb = nnode(ltypb(iboun)) ! Nodes on boundary | how it was in Alya 1
        ielem = lelbo(iboun) ! MPIO variable
        ! ielem = lboel(pnodb+1,iboun) ! Element where face belongs | Alya 1
        pelty = abs(ltype(ielem)) ! Alya 2
        ! pelty = ltype(ielem) ! Element type
        pnode = nnode(pelty) ! Nodes on element

        if (ndime == 3) then ! 3D elements

           if (size(lnodb(:,iboun),1) == 3) then ! TET04
              nncorner = 3
              nnedg = 0
              nnface = 0
           else if (size(lnodb(:,iboun),1) == 4) then ! Could be HEX08 or PEN06
              if (lnodb(4,iboun) == 0) then ! PEN06: a triangular face with 3 nodes
                 nncorner = 3
                 nnedg = 0
                 nnface = 0
              else ! It's a quadrangular face
                 nncorner = 4
                 nnedg = 0
                 nnface = 0
              end if
           else if (size(lnodb(:,iboun),1) == 6) then ! TET10
              nncorner = 3
              nnedg = 3
              nnface = 0
           else if (size(lnodb(:,iboun),1) == 9) then ! 2nd order element; could be PEN18 or HEX27
              if (pnode == 18) then ! PEN18 is a VERY FUCKING SPECIAL CASE
                 if (lnodb(7,iboun) == 0) then ! Detects a triangular face according to boundary defs. in *.geo.dat file
                    nncorner = 3
                    nnedg = 3
                    nnface = 0
                 else ! It's a quadrangular face
                    nncorner = 4
                    nnedg = 4
                    nnface = 1
                 end if
              else ! If not PEN18 but 9 nodes on boundary, then it's a HEX27 face
                 nncorner = 4
                 nnedg = 4
                 nnface = 1
              end if        
           else ! Implies weird or non-coded elements
              print*, 'Wrong element type, please revise...'
           end if
        else if (ndime == 2) then ! 2D elements
           if (pnode == 3) then ! TRI03
              nncorner = 2
              nnedg = 0
              nnface = 0
           else if (pnode == 4) then ! QUA04
              nncorner = 2
              nnedg = 0
              nnface = 0
           else if (pnode == 6) then ! TRI06
              nncorner = 2
              nnedg = 0
              nnface = 1
           else if (pnode == 9) then ! QUA09
              nncorner = 2
              nnedg = 0
              nnface = 1
           else
              print*, '--| ELEMENT ', pelty ,' NOT SUPPORTED!'
           end if
        else
           print*, '--| 1D ELEMENTS NOT CODED!'
        end if 
        pface = nface(pelty) ! Faces on element
        inodb = 1 ! Generic boundary node counter
        inode = lboel(inodb,iboun) ! Should it be like this?
        ipoin = lnodb(inodb,iboun) ! 1st boundary node
        jface = 0 ! Element face counter
        select_face = .false. ! Flag to select element face
        ! Cycle over element faces
        do while (jface < pface)
           jface = jface + 1
           iface = jface ! Element face index
           jnodf = 0 ! Element face node counter
           ! Cycle over element face nodes
           do while (jnodf < nnodf(pelty) % l(iface))
              jnodf = jnodf + 1
              inodf = jnodf ! Elem. face node index
              ! Start operating on face when ipoin has a match on element nodes list
              ! Is this the criteria? How would work on generalized HO?
              if(ipoin == lnods( lface(pelty) % l(inodf,iface),ielem) )then ! Check if bnode exists on elem. iface
                 cont = 1
                 ! Search for other element face nodes, and increment node counter
                 do inodb = 2,pnodb
                    do jnodb = 1, nnodf(pelty) % l(iface)
                       ! Criterion for reordering
                       if(lnodb(inodb,iboun)==lnods( lface(pelty) % l(jnodb,iface),ielem) )then
                          cont = cont + 1
                       end if
                    end do
                 end do
                 ! Select face once all face nodes are found
                 if(cont==pnodb)then
                    jface = pface
                    jnodf = nnodf(pelty) % l(iface)
                    select_face = .true.
                 end if
                 change_order = .true.
                 ! Start of reorder ops.
                 if(select_face) then
                    cont = 1
                    ! Cycle from next to last elem. face node, keeping 1st node untouched
                    do jnodb=2,pnodb
                       ! Tested only for TET4-10 and HEX8-27 (HEX20 not included), as well as PEN06-18.
                       if(jnodb==2)then ! Inside corner block
                          if(inodf == nncorner)then ! == pnodb
                             knodf = 1
                          else
                             knodf = inodf+1
                          end if
                       else if(jnodb==nncorner)then ! original is == pnodb; end of corner block
                          if(inodf==1)then
                             knodf = nncorner ! pnodb
                          else
                             knodf = inodf - 1
                          end if
                       else if (nnedg /=0 .and. jnodb == nncorner+1) then ! starts edge block, if edges have nodes
                          knodf = nncorner+inodf
                       else if (nnedg /=0 .and. jnodb == nncorner+2) then
                          if (inodf == 1) then
                             knodf = jnodb
                          else if (inodf == nncorner) then
                             knodf = nncorner+1
                          else
                             knodf = nncorner+inodf+1
                          end if
                       else if (nnedg /=0 .and. jnodb == nncorner+nnedg) then ! end of edge block
                          if (inodf == 1) then
                             knodf = nncorner+nnedg
                          else
                             knodf = nncorner+inodf-1
                          end if
                       else if (nnface /= 0 .and. jnodb == pnodb) then ! If the face has a node, count it
                          knodf = pnodb
                       else ! All other nodes
                          if (jnodb > 2 .and. jnodb < nncorner) then
                             if (inodf == 1) then
                                knodf = 3
                             else if (inodf == 2) then
                                knodf = 4
                             else if (inodf == 3) then
                                knodf = 1
                             else if (inodf == 4) then
                                knodf = 2
                             end if
                          else if (jnodb > nncorner+2 .and. jnodb < nncorner+nnedg) then
                             if (inodf == 1) then
                                knodf = jnodb
                             else if (inodf == 2) then
                                knodf = jnodb+1
                             else if (inodf == 3) then
                                knodf = jnodb-2
                             else if (inodf == 4) then
                                knodf = jnodb-1
                             end if
                          end if
                       end if
                       if(lnodb(jnodb,iboun) == lnods(lface(pelty) % l(knodf,iface),ielem) )cont = cont + 1
                    end do
                    if (cont == pnodb) change_order = .false. 
                    if ( change_order ) then
                       ! Change master nodes 2 and nncorner:
                       kpoin = lnodb(nncorner,iboun) ! pnodb instead of nncorner
                       lnodb(nncorner,iboun) = lnodb(2,iboun)
                       lnodb(2,iboun) = kpoin
                       ! Reverse edge nodes:
                       do jnodb = nncorner+1,nncorner+nnedg
                          kpoin = lnods(lface(pelty) % l(jnodb,iface),ielem)
                          lnodb(jnodb,iboun) = kpoin
                       end do
                       ! Dealing with lboel...
                       knode = lboel(nncorner,iboun) ! pnodb
                       lboel(nncorner,iboun) = lboel(2,iboun)
                       lboel(2,iboun) = knode
                       if (nnedg /= 0) then
                          if (nnedg == 3) then
                             knode = lboel(nncorner+nnedg,iboun)
                             lboel(nncorner+nnedg,iboun) = lboel(nncorner+1,iboun)
                             lboel(nncorner+1,iboun) = knode
                          else if (nnedg == 4) then
                             knode = lboel(nncorner+nnedg,iboun)
                             lboel(nncorner+nnedg,iboun) = lboel(nncorner+1,iboun)
                             lboel(nncorner+1,iboun) = knode
                             knode = lboel(nncorner+nnedg-1,iboun)
                             lboel(nncorner+nnedg-1,iboun) = lboel(nncorner+2,iboun)
                             lboel(nncorner+2,iboun) = knode
                          end if
                       end if              
                    end if
                 end if
              end if
           end do
        end do
     end do

  end if


end subroutine order_boundaries
