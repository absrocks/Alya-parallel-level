!-----------------------------------------------------------------------
!> @addtogroup Dodeme
!> @{
!> @file    dod_holcut_markelements.f90
!> @author  Guillaume Houzeaux
!> @date    18/09/2012
!> @brief   Finish hole cutting 
!> @details Starting with a hole element, mark elements recursively
!>          going from face to face. The process is the following:
!>          - Find the nearest hole elements to the c.o.g of the!
!>            complete hole
!>          - Mark recursively the elements using the face-based
!>            element-element graph 
!>
!> @} 
!-----------------------------------------------------------------------

subroutine dod_holcut_markelements()
  use def_parame
  use def_master
  use def_domain
  use def_dodeme
  use mod_memchk
  use mod_memory
  use mod_postpr
  use mod_graphs
  use mod_messages, only : livinf
  implicit none
  integer(ip)            :: isubd,intij,jsubd,kelem,pnode,ielty
  integer(ip)            :: inode,ipoin,ielem,jelem,idime
  integer(ip)            :: istack,nstack,ielel,iface,nnodg_jelem
  integer(ip)            :: nfacg_subd,ifacg,jfacg,kface,jpoin
  integer(ip)            :: marked,jface,jnode,jnodg,ksubd
  integer(ip)            :: pelty_ielem,i_have_a_hole
  integer(ip)            :: pelty_jelem
  integer(ip)            :: nelem_subd
  integer(ip)            :: npoin_subd
  real(rp)               :: dista,xcog(3),xfact
  logical(lg)            :: mark_this_element

  integer(ip),   pointer :: lnnod_subd(:)
  integer(ip),   pointer :: ltype_subd(:)
  integer(ip),   pointer :: lnods_subd(:,:)
  integer(ip),   pointer :: lstack(:)
  integer(ip),   pointer :: cfael(:,:,:)
  integer(ip),   pointer :: nnodg(:,:)
  integer(ip),   pointer :: lface_subd(:,:)
  real(rp),      pointer :: center_gravity(:,:)
  type(i1p),     pointer :: lelfa_subd(:)
  integer(ip),   pointer :: marked_nodes(:)
  logical(lg),   pointer :: marked_faces(:)
  logical(lg),   pointer :: lchek(:)
!!$  integer(ip) :: ppoin,pelem
!!$  real(rp) :: xcolor

  call livinf(0_ip,'MARK HOLE ELEMENTS RECURSIVELY',0_ip)

!!$  write(100,*) 'GiD Post Results File 1.0'
!!$  write(100,*) ' '
!!$  write(100,*) 'GaussPoints GP Elemtype Triangle'
!!$  write(100,*) 'Number of Gauss Points:   1'
!!$  write(100,*) 'Natural Coordinates: Internal'
!!$  write(100,*) 'End GaussPoints'  
!!$  igene = 0
!!$  write(100,*) 'Result INITIAL ALYA ',igene,' Scalar OnGaussPoints GP'
!!$  write(100,*) 'ComponentNames INITIAL'
!!$  write(100,*) 'Values'
!!$  isubd = 1
!!$  current_subdomain => subdomain(isubd)
!!$  do kelem = 1,current_subdomain % nelem
!!$     if( current_subdomain % lsubd_nelem(kelem) < 0 ) write(100,*) current_subdomain % leper(kelem),1
!!$  end do
!!$  write(100,*) 'End values'
!!$  
!!$  open(unit=100,file='mesh.msh',status='unknown')
!!$  ppoin = 161
!!$  pelem = 0
!!$  xcolor = 0.0_rp

  nullify(lstack)
  nullify(cfael)
  nullify(nnodg)
  nullify(lface_subd)
  nullify(center_gravity)
  nullify(lelfa_subd)
  nullify(marked_nodes)
  nullify(marked_faces)
  nullify(lchek)

  do isubd = 1,nsubd

     if( ihole_dod(isubd) == 1 ) then

        current_subdomain => subdomain(isubd)
        lnods_subd        => current_subdomain % lnods
        lnnod_subd        => current_subdomain % lnnod
        ltype_subd        => current_subdomain % ltype
        nelem_subd        =  current_subdomain % nelem
        npoin_subd        =  current_subdomain % npoin
        call memory_alloca(mem_servi(1:2,servi),'CENTER_GRAVITY','dod_holcut_markelements',center_gravity,ndime,nelem_subd)
        call memory_alloca(mem_servi(1:2,servi),'MARKED_NODES'  ,'dod_holcut_markelements',marked_nodes,npoin_subd)
        !
        ! Compute list faces for hole elements
        !
        call memory_alloca(mem_servi(1:2,servi),'CFAEL','dod_holcut_markelements',cfael,mnodb,mface,nelty)
        call memory_alloca(mem_servi(1:2,servi),'NNODG','dod_holcut_markelements',nnodg,mface,nelty)     
        call memory_alloca(mem_servi(1:2,servi),'LCHEK','dod_holcut_markelements',lchek,nelem_subd)
        do ielty = iesta_dom,iesto_dom
           !if( lexis(ielty) /= 0 ) &
           !     call domfa2(&                     
           !     mnodb,nnodg(1,ielty),nface(ielty),ielty,cfael(1,1,ielty)) 
           call runend('DOD_HOLCUT_MARKELEMENTS: RECODE')
        end do
        !
        ! If a hole holed-patch exist, then compute the entire element graph
        ! If only a Chimera exists, then the recursive can be performed only
        ! on the hole elements
        !
        ksubd = 0
        do jsubd = 1,nsubd
           intij = intyp_dod(isubd,jsubd)
           if( intij > 0 ) then
              if( ictop_dod(intij) == DOD_HOLED_PATCH ) then
                 ksubd = 1
              end if
           end if
        end do
        if( ksubd == 0 ) then
           do ielem = 1,nelem_subd
              if( current_subdomain % lsubd_nelem(ielem) < 0 .or. current_subdomain % lsubd_nelem(ielem) == DOD_SOLID ) then
                 lchek(ielem) = .true.
              end if
           end do
        else
           do ielem = 1,nelem_subd
              lchek(ielem) = .true.
           end do
        end if

        call graphs_list_faces(& 
             current_subdomain % nelem,mnode,mnodb,nelty,mface,current_subdomain % lnods,&
             current_subdomain % lnnod,current_subdomain % ltype,nnodg,cfael,nface,& 
             current_subdomain % pelpo,current_subdomain % lelpo,nfacg_subd,lface_subd,&
             lelfa_subd,lchek)
        call memory_alloca(mem_servi(1:2,servi),'MARKED_FACES','dod_holcut_markelements',marked_faces,nfacg_subd)
        !
        ! Loop over neighbors
        !
        do jsubd = 1,nsubd
           intij = intyp_dod(isubd,jsubd)
           if( intij > 0 ) then
              if( ictop_dod(intij) == DOD_CHIMERA .or. ictop_dod(intij) == DOD_HOLED_PATCH ) then

                 i_have_a_hole = 0                 
                 do ielem = 1,nelem_subd
                    if( current_subdomain % lsubd_nelem(ielem) == -jsubd ) i_have_a_hole = i_have_a_hole + 1
                 end do

                 if( i_have_a_hole > 0 ) then
     
                    if( ictop_dod(intij) == DOD_CHIMERA ) then
                       !
                       ! Compute center of gravity of hole XCOG(1:NDIME)
                       ! Do not take into account solid elements not to take the
                       ! solid of another subdomain
                       !        
                       kelem = 0
                       xcog  = 0.0_rp
                       do ielem = 1,nelem_subd
                          if( current_subdomain % lsubd_nelem(ielem) == -jsubd ) then
                             pnode = current_subdomain % lnnod(ielem)
                             kelem = kelem + 1
                             center_gravity(1:ndime,ielem) = 0.0_rp
                             do inode = 1,pnode
                                ipoin = current_subdomain % lnods(inode,ielem)
                                do idime = 1,ndime
                                   center_gravity(idime,ielem) = center_gravity(idime,ielem) + current_subdomain % coord(idime,ipoin)
                                end do
                             end do
                             xfact = 1.0_rp / real(pnode,rp)
                             do idime = 1,ndime
                                center_gravity(idime,ielem) = xfact * center_gravity(idime,ielem)
                                xcog(idime) = xcog(idime) + center_gravity(idime,ielem)
                             end do
                          end if
                       end do
                       if( kelem == 0 ) call runend('DOD_HOLCUT_MARKELEMENTS: DID NOT FIND AN ELEMENT')
                       do idime = 1,ndime
                          xcog(idime) = xcog(idime) / real(kelem,rp)
                       end do
                       !
                       ! Look for nearest element KELEM from center of gravity
                       !
                       dista = 1e9_rp
                       do ielem = 1,nelem_subd
                          if( current_subdomain % lsubd_nelem(ielem) == -jsubd ) then
                             xfact = 0.0_rp
                             do idime = 1,ndime
                                xfact = xfact + ( xcog(idime) - center_gravity(idime,ielem) )**2
                             end do
                             if( xfact < dista ) then
                                dista = xfact 
                                kelem = ielem
                             end if
                          end if
                       end do
                       ksubd = -jsubd

                    else                  
                       !
                       ! For holed patch, the recursive algorithm should be applied to free elements
                       !
                       ielem = 0
                       do while( ielem < nelem_subd )
                          ielem = ielem + 1
                          if( current_subdomain % lsubd_nelem(ielem) == 0 ) then
                             kelem = ielem
                             ielem = nelem_subd
                          end if
                       end do
                       ksubd = 0
                    end if
                    !
                    ! Touch element with common faces recursively starting with KELEM
                    !
                    call memory_alloca(mem_servi(1:2,servi),'LSTACK','dod_holcut_markelements',lstack,nelem_subd)
                    call memgen(1_ip,nelem_subd,0_ip)

                    ielem        = kelem
                    nstack       = 1
                    lstack(1)    = ielem
                    gisca(ielem) = 1
                    istack       = 0                
                    pelty_ielem  = abs(ltype_subd(ielem))
                    do inode = 1,lnnod_subd(ielem)
                       ipoin = lnods_subd(inode,ielem)
                       marked_nodes(ipoin) = 1
                    end do
                    do iface = 1,nface(pelty_ielem)
                       ifacg = lelfa_subd(ielem) % l(iface)
                       marked_faces(ifacg) = .true.
                    end do

                    do 
                       if( istack == nstack ) exit
                       istack = istack+1   
                       ielem  = lstack(istack)
                       do ielel = current_subdomain % pelel(ielem), current_subdomain % pelel(ielem+1)-1
                          jelem = current_subdomain % lelel(ielel)

                          if( gisca(jelem) == 0 ) then
                             if( current_subdomain % lsubd_nelem(jelem) == ksubd .or. current_subdomain % lsubd_nelem(jelem) == DOD_SOLID ) then
                                !
                                ! Check if new nodes (o) are marked or not. Marked nodes are x
                                !    
                                !     o
                                !    / \  jelem
                                !  /     \
                                ! x-------x <= kface
                                ! |       |
                                ! | ielem |
                                ! |       |
                                ! x-------x
                                !
                                !
                                ! KFACE= face which connects IELEM and JELEM 
                                !
                                kface = 0
                                pelty_jelem = abs(ltype_subd(jelem))
                                faces_loop: do iface = 1,nface(pelty_ielem)
                                   do jface = 1,nface(pelty_jelem)
                                      if( lelfa_subd(ielem) % l(iface) == lelfa_subd(jelem) % l(jface) ) then
                                         kface = iface
                                         exit faces_loop                                     
                                      end if
                                   end do
                                end do faces_loop

                                if( kface == 0 ) then
                                   print*,'a=',current_subdomain % leper(ielem),current_subdomain % leper(jelem)
                                   print*,'b=',lelfa_subd(ielem) % l(:)
                                   print*,'c=',lelfa_subd(jelem) % l(:)                             
                                   call runend('DOD_HOLCUT_MARKELEMENTS: WE ARE IN TROUBLE')
                                end if
                                !
                                ! Check if marked nodes are only on the face
                                !
                                nnodg_jelem = nnodg(kface,pelty_jelem)
                                marked = 0
                                do jnode = 1,lnnod_subd(jelem)
                                   jpoin  = lnods_subd(jnode,jelem)
                                   marked = marked + marked_nodes(jpoin)
                                end do
                                mark_this_element = .true.

                                if( marked == nnodg_jelem ) then
                                   continue
                                else
                                   !
                                   !        no
                                   !    o-------x -1     o-------x  1     
                                   !    |       |        |       |
                                   ! no |       | si  no |       | no
                                   !    |       |        |       |
                                   !    x-------x        x-------x 
                                   !   -1   si -1       -1   si -1
                                   ! 
                                   do jface = 1,nface(pelty_jelem)
                                      jfacg = lelfa_subd(jelem) % l(jface)
                                      if( marked_faces(jfacg) ) then
                                         do jnodg = 1,nnodg(jface,pelty_jelem)
                                            jpoin = lnods_subd(cfael(jnodg,jface,pelty_jelem),jelem)
                                            marked_nodes(jpoin) = -1
                                         end do
                                      end if
                                   end do
                                   do jnode = 1,lnnod_subd(jelem)
                                      jpoin  = lnods_subd(jnode,jelem)
                                      if( marked_nodes(jpoin) > 0 ) then
                                         mark_this_element = .false.                    !  1
                                      else
                                         marked_nodes(jpoin) = abs(marked_nodes(jpoin)) ! -1,0 => 1,0
                                      end if
                                   end do
                                end if

                                if( mark_this_element ) then
                                   !
                                   ! Mark element JELEM, all its nodes and faces
                                   !
                                   gisca(jelem)   = 1
                                   nstack         = nstack + 1
                                   lstack(nstack) = jelem
                                   do jnode = 1,lnnod_subd(jelem)
                                      jpoin = lnods_subd(jnode,jelem)
                                      marked_nodes(jpoin) = 1
                                   end do
                                   do jface = 1,nface(pelty_jelem)
                                      jfacg = lelfa_subd(jelem) % l(jface)
                                      if( jfacg == 0 ) call runend('WE ARE IN VERY BIG TROUBLE')
                                      marked_faces(jfacg) = .true.
                                   end do
                                end if

!!$                             if( mark_this_element ) then
!!$                                write(100,'(a)') '# encoding utf-8'
!!$                                write(100,'(a)') 'MESH "Hole element" dimension 3 ElemType Triangle Nnode 3'
!!$                                write(100,'(a,3(1x,f12.9))') '# Color',xcolor,xcolor,xcolor
!!$                                xcolor = xcolor + 0.0333_rp
!!$                                write(100,'(a)') 'coordinates'
!!$                                do jnode = 1,lnnod_subd(jelem)
!!$                                   jpoin  = lnods_subd(jnode,jelem)
!!$                                   ppoin = ppoin + 1
!!$                                   write(100,*) ppoin,current_subdomain % coord(1,jpoin),current_subdomain % coord(2,jpoin),0.31_rp-xcolor*0.3_rp
!!$                                end do
!!$                                write(100,'(a)') 'end coordinates'
!!$                                write(100,'(a)') 'elements'
!!$                                pelem = pelem + 1
!!$                                write(100,'(10(1x,i6))') pelem,ppoin-2,ppoin-1,ppoin
!!$                                write(100,'(a)') 'end elements'
!!$                             end if

!!$                             if( mark_this_element ) then
!!$                                write(100,*) 'Result CACAC ALYA ',igene,' Scalar OnGaussPoints GP'
!!$                                write(100,*) 'ComponentNames CACAC'
!!$                                write(100,*) 'Values'
!!$                                igene = igene + 1                                                          
!!$                                do kelem = 1,nelem_subd
!!$                                   if( gisca(kelem) > 0.0_rp ) write(100,*) current_subdomain % leper(kelem),1
!!$                                end do
!!$                                write(100,*) 'End values'
!!$                             end if

                             end if
                          end if
                       end do
                    end do
                    !
                    ! Redefine hole elements using touched element (GISCA(IELEM) = 1)
                    ! Mark nodes 
                    !
                    if( ictop_dod(intij) == DOD_CHIMERA ) then
                       !
                       ! Marked elements are hole elements
                       !
                       do ielem = 1,nelem_subd
                          if( gisca(ielem) == 1 ) then
                             current_subdomain % lsubd_nelem(ielem) = -jsubd
                             pnode = current_subdomain % lnnod(ielem)
                             do inode = 1,pnode
                                ipoin = current_subdomain % lnods(inode,ielem)
                                current_subdomain % lsubd_npoin(ipoin) = -jsubd
                             end do
                          else if( current_subdomain % lsubd_nelem(ielem) == -jsubd ) then
                             current_subdomain % lsubd_nelem(ielem) = 0
                          end if
                       end do
                    else                       
                       !
                       ! Marked elements are free elements
                       !
                       do ielem = 1,nelem_subd
                          if( gisca(ielem) == 0 ) then
                             current_subdomain % lsubd_nelem(ielem) = -jsubd
                             pnode = current_subdomain % lnnod(ielem)
                             do inode = 1,pnode
                                ipoin = current_subdomain % lnods(inode,ielem)
                                current_subdomain % lsubd_npoin(ipoin) = -jsubd
                             end do
                          else 
                             current_subdomain % lsubd_nelem(ielem) = 0
                          end if
                       end do
                    end if
                    call memgen(3_ip,nelem_subd,0_ip)
                    call memory_deallo(mem_servi(1:2,servi),'LSTACK','dod_holcut_markelements',lstack)
                 end if

              end if
           end if
        end do
        !OJOOO!!!
        call memory_deallo(mem_servi(1:2,servi),'CFAEL','dod_holcut_markelements',cfael)
        call memory_deallo(mem_servi(1:2,servi),'NNODG','dod_holcut_markelements',nnodg)     
        !call memory_deallo(mem_servi(1:2,servi),'LCHEK','dod_holcut_markelements',lchek)

        call graphs_deallocate_list_faces(lface_subd,lelfa_subd,lchek)
        call memory_deallo(mem_servi(1:2,servi),'MARKED_FACES',  'dod_holcut_markelements',marked_faces)
        call memory_deallo(mem_servi(1:2,servi),'LCHEK',         'dod_holcut_markelements',lchek)
        call memory_deallo(mem_servi(1:2,servi),'MARKED_NODES'  ,'dod_holcut_markelements',marked_nodes)
        call memory_deallo(mem_servi(1:2,servi),'CENTER_GRAVITY','dod_holcut_markelements',center_gravity)       

     end if

  end do

end subroutine dod_holcut_markelements
