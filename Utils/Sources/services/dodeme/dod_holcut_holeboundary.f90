!-----------------------------------------------------------------------
!> @addtogroup Dodeme
!> @{
!> @file    dod_holcut_holeboundary.f90
!> @author  Guillaume Houzeaux
!> @date    18/09/2012
!> @brief   Contruct hole boundary
!> @details Process is the following:
!>          1. Compute faces for hole and their neighbors:
!>             A hole has: current_subdomain % lsubd_nelem(ielem) < 0
!>          
!> @} 
!-----------------------------------------------------------------------
subroutine dod_holcut_holeboundary()
  use def_parame
  use def_master
  use def_elmtyp
  use def_domain
  use def_dodeme
  use mod_memory
  use mod_kdtree
  use mod_elmgeo
  implicit none 
  integer(ip)                   :: isubd,ielty,ielem,iboun,iface,inodb,ilist
  integer(ip)                   :: inode,jelem,jface,jelty,jnode,jblty,pnodb
  integer(ip)                   :: nelel,ielel,ii,ninve,kboun
  integer(ip)                   :: nelem_subd,nboun_subd,nboun_new,nboun_ori
  integer(ip)                   :: jsubd,iboun_global,nboun_hole
  integer(ip)                   :: lnodb_aux(mnodb),lboel_aux(mnodb)
  logical(lg)                   :: equal_faces
  type(i2p),            pointer :: faces(:)
  integer(ip),          pointer :: cfael(:,:,:)
  integer(ip),          pointer :: nnodg(:,:)
  integer(ip),          pointer :: lnods_subd(:,:)
  integer(ip),          pointer :: ltype_subd(:)
  integer(ip),          pointer :: pelel_subd(:)
  integer(ip),          pointer :: lelel_subd(:)
  integer(ip),          pointer :: lboch_subd(:) 
  integer(ip),          pointer :: lnodb_subd(:,:) 
  integer(ip),          pointer :: ltypb_subd(:) 
  integer(ip),          pointer :: lboel_subd(:,:) 
  integer(ip),          pointer :: lbper_subd(:) 
  type(typ_subdomain2), target  :: copy_subdomain(1)
  logical(lg),          pointer :: lmark(:)
  real(rp)                      :: dummr(3)

 nullify(cfael)
 nullify(faces)
 nullify(nnodg)
 nullify(lnods_subd)
 nullify(ltype_subd)
 nullify(pelel_subd)
 nullify(lelel_subd)
 nullify(lboch_subd) 
 nullify(lnodb_subd) 
 nullify(ltypb_subd) 
 nullify(lboel_subd) 
 nullify(lbper_subd)
 nullify(lmark)

  do isubd = 1,nsubd

     if( ihole_dod(isubd) /= 0 ) then

        current_subdomain => subdomain(isubd) 
        nelem_subd        =  subdomain(isubd) % nelem
        lnods_subd        => subdomain(isubd) % lnods
        ltype_subd        => subdomain(isubd) % ltype
        pelel_subd        => subdomain(isubd) % pelel
        lelel_subd        => subdomain(isubd) % lelel
        !
        ! Allocate memory for FACES, CFAEL AND NNODG
        !
        call memory_alloca(mem_servi(1:2,servi),'FACES','dod_holcut_holeboundary',faces,nelem_subd)
        call memory_alloca(mem_servi(1:2,servi),'CFAEL','dod_holcut_holeboundary',cfael,mnodb,mface,nelty,'DO_NOT_INITIALIZE')
        call memory_alloca(mem_servi(1:2,servi),'NNODG','dod_holcut_holeboundary',nnodg,mface,nelty,'DO_NOT_INITIALIZE')
        call memory_alloca(mem_servi(1:2,servi),'LMARK','dod_holcut_holeboundary',lmark,nelem_subd)
        !
        ! Compute CFAEL AND NNODG
        !
        do ielty = iesta_dom,iesto_dom
           !if( lexis(ielty) /= 0 ) &
           !     call domfa2(mnodb,nnodg(1,ielty),nface(ielty),ielty,cfael(1,1,ielty)) 
           call runend('DOD_HOLCUT_HOLEBOUNDARY: RECODE')
        end do
        !
        ! Mark hole elements IELEM in contact with free elements JELEM and these free elements 
        ! LMARK(IELEM) = .TRUE.
        ! LMARK(JELEM) = .TRUE.
        !
        do ielem = 1,nelem_subd          
           if( current_subdomain % lsubd_nelem(ielem) < 0 ) then     
              do ielel = pelel_subd(ielem),pelel_subd(ielem+1)-1
                 jelem = lelel_subd(ielel)
                 if( current_subdomain % lsubd_nelem(jelem) >= 0 ) then
                    lmark(ielem) = .true.
                    lmark(jelem) = .true.
                 end if
              end do
           end if
        end do
        !
        ! Construct and sort FACES for marked element with LMARK(IELEM) = .TRUE.
        !
        do ielem = 1,nelem_subd          
           if( lmark(ielem) ) then
              ielty = abs(ltype_subd(ielem))
              nullify(faces(ielem) % l) 
              call memory_alloca(mem_servi(1:2,servi),'FACES % L','dod_holcut_holeboundary',faces(ielem) % l,mnodb+2_ip,mface)
              do iface = 1,nface(ielty)
                 do inodb = 1,nnodg(iface,ielty)
                    inode = cfael(inodb,iface,ielty)
                    faces(ielem) % l(inodb,iface) = lnods_subd(inode,ielem)
                 end do
                 call sortin(nnodg(iface,ielty),faces(ielem) % l(1,iface))
              end do
           end if
        end do
        !
        ! Compute NBOUN_NEW and fill in FACES
        !
        nboun_new = 0
        do ielem = 1,nelem_subd                                       
           !
           ! Compare faces and eliminate the repeated faces of holes
           !
           if( current_subdomain % lsubd_nelem(ielem) < 0 .and. lmark(ielem) ) then
              ielty = abs(ltype_subd(ielem))          
              do iface = 1,nface(ielty)
                 ilist = 1
                 nelel = pelel_subd(ielem+1)-pelel_subd(ielem)
                 ielel = pelel_subd(ielem)-1
                 do while( ilist <= nelel )
                    ielel = ielel + 1
                    jelem = lelel_subd(ielel)
                    if( current_subdomain % lsubd_nelem(jelem) >= 0 ) then
                       jelty = abs(ltype_subd(jelem))  
                       jface = 0
                       do while( jface /= nface(jelty) )
                          jface = jface + 1
                          !if( faces(jelem) %l(1,jface) /= 0 ) then
                          equal_faces = .true.
                          inodb       = 0
                          do while( equal_faces .and. inodb /= nnodg(jface,jelty) ) 
                             inodb = inodb + 1
                             if( faces(ielem) % l(inodb,iface) /= faces(jelem) % l(inodb,jface) ) equal_faces = .false.
                          end do
                          !
                          ! An external face to the hole has been found: JFACE in JELEM
                          !
                          if( equal_faces ) then
                             faces(ielem) % l(mnodb+2,iface) = jelem        ! Save neighboring element JELEM
                             faces(ielem) % l(mnodb+1,iface) = jface        ! Save neighboring face JFACE
                             jface                           = nface(jelty) ! Go out of face loop
                             !faces(jelem) % l(1,jface)       = 0            ! Do not test any more this free-side face
                             ilist                           = nelel        ! Go out of neighbor loop
                          end if
                          !end if
                       end do
                    end if
                    ilist = ilist + 1
                 end do
                 !
                 ! This is a new boundary: reorder and copy face JFACE of JELEM to face IFACE of IELEM
                 ! Reordering is necessary as FACES was ordered at the beginning
                 !
                 if( faces(ielem) % l(mnodb+2,iface) /= 0 ) then
                    jelem     = faces(ielem) % l(mnodb+2,iface)
                    jface     = faces(ielem) % l(mnodb+1,iface)
                    jelty     = abs(ltype_subd(jelem))
                    nboun_new = nboun_new + 1
                    do inodb = 1,nnodg(jface,jelty)
                       jnode = cfael(inodb,jface,jelty)
                       faces(ielem) %l(inodb,iface) = lnods_subd(jnode,jelem)
                    end do
                 end if
              end do
           end if
        end do
        !
        ! Eliminate hole-boundaries, that is boundaries connected to a hole element
        !
        nboun_hole = 0
        do iboun = 1,current_subdomain % nboun
           pnodb = current_subdomain % ltypb(iboun)
           ielem = current_subdomain % lelbo(iboun)
           if( current_subdomain % lsubd_nelem(ielem) < 0 ) nboun_hole = nboun_hole + 1         
        end do
        if ( nboun_hole /= 0 ) then
           kboun = 0
           do iboun = 1,current_subdomain % nboun
              pnodb = current_subdomain % ltypb(iboun)
              ielem = current_subdomain % lelbo(iboun)
              if( current_subdomain % lsubd_nelem(ielem) < 0 .and. current_subdomain % lboch(iboun) == BOFEM ) then
                 current_subdomain % lboch(iboun) = BOHOL
              end if
         !        kboun                              = kboun + 1
         !        current_subdomain % lbper(  kboun) = current_subdomain % lbper(  iboun)
         !        current_subdomain % lboch(  kboun) = current_subdomain % lboch(  iboun)
         !        current_subdomain % lnodb(:,kboun) = current_subdomain % lnodb(:,iboun)
         !        current_subdomain % ltypb(  kboun) = current_subdomain % ltypb(  iboun)
         !        current_subdomain % lboel(:,kboun) = current_subdomain % lboel(:,iboun)
         !     end if              
           end do           
         !  current_subdomain % nboun = kboun
         !  nboun_subd = current_subdomain % nboun 
         !  call memory_resize(mem_servi(1:2,servi),'LBPER','dod_holcut_holeboundary',current_subdomain % lbper,nboun_subd)
         !  call memory_resize(mem_servi(1:2,servi),'LBOCH','dod_holcut_holeboundary',current_subdomain % lboch,nboun_subd)
         !  call memory_resize(mem_servi(1:2,servi),'LNODB','dod_holcut_holeboundary',current_subdomain % lnodb,mnodb,nboun_subd)
         !  call memory_resize(mem_servi(1:2,servi),'LTYPB','dod_holcut_holeboundary',current_subdomain % ltypb,nboun_subd) 
         !  call memory_resize(mem_servi(1:2,servi),'LBOEL','dod_holcut_holeboundary',current_subdomain % lboel,mnodb,nboun_subd) 
        end if
        !
        ! Reallocate memory for boundary mesh
        !
        nboun_ori                 = current_subdomain % nboun
        current_subdomain % nboun = nboun_ori + nboun_new 
        nboun_subd                = current_subdomain % nboun
        !
        ! Copy old mesh
        !  
        nullify(copy_subdomain(1) % lbper)
        nullify(copy_subdomain(1) % lboch)
        nullify(copy_subdomain(1) % lnodb)
        nullify(copy_subdomain(1) % ltypb)
        nullify(copy_subdomain(1) % lboel) 
        nullify(copy_subdomain(1) % lsubd_nboun)

        call memory_copy(  mem_servi(1:2,servi),'LBPER'      ,'dod_membou',current_subdomain % lbper,      copy_subdomain(1) % lbper)
        call memory_copy(  mem_servi(1:2,servi),'LBOCH'      ,'dod_membou',current_subdomain % lboch,      copy_subdomain(1) % lboch)
        call memory_copy(  mem_servi(1:2,servi),'LNODB'      ,'dod_membou',current_subdomain % lnodb,      copy_subdomain(1) % lnodb)
        call memory_copy(  mem_servi(1:2,servi),'LTYPB'      ,'dod_membou',current_subdomain % ltypb,      copy_subdomain(1) % ltypb)
        call memory_copy(  mem_servi(1:2,servi),'LBOEL'      ,'dod_membou',current_subdomain % lboel,      copy_subdomain(1) % lboel)
        call memory_copy(  mem_servi(1:2,servi),'LELBO'      ,'dod_membou',current_subdomain % lelbo,      copy_subdomain(1) % lelbo)
        call memory_copy(  mem_servi(1:2,servi),'LSUBD_NBOUN','dod_membou',current_subdomain % lsubd_nboun,copy_subdomain(1) % lsubd_nboun)
        !
        ! Reallocate
        !
        call memory_alloca(mem_servi(1:2,servi),'LBPER'      ,'dod_cpymsh',current_subdomain % lbper,nboun_subd,           'DO_NOT_INITIALIZE')
        call memory_alloca(mem_servi(1:2,servi),'LBOCH'      ,'dod_cpymsh',current_subdomain % lboch,nboun_subd,           'DO_NOT_INITIALIZE')
        call memory_alloca(mem_servi(1:2,servi),'LNODB'      ,'dod_cpymsh',current_subdomain % lnodb,mnodb,nboun_subd,     'DO_NOT_INITIALIZE')
        call memory_alloca(mem_servi(1:2,servi),'LTYPB'      ,'dod_cpymsh',current_subdomain % ltypb,nboun_subd,           'DO_NOT_INITIALIZE')
        call memory_alloca(mem_servi(1:2,servi),'LBOEL'      ,'dod_cpymsh',current_subdomain % lboel,mnodb,nboun_subd,     'DO_NOT_INITIALIZE')
        call memory_alloca(mem_servi(1:2,servi),'LELBO'      ,'dod_cpymsh',current_subdomain % lelbo,nboun_subd,           'DO_NOT_INITIALIZE')
        call memory_alloca(mem_servi(1:2,servi),'LSUBD_NBOUN','dod_cpymsh',current_subdomain % lsubd_nboun,nboun_subd,     'DO_NOT_INITIALIZE')
        !
        ! Copy the copy on the new one
        !
        call memory_copy(  mem_servi(1:2,servi),'LBPER',      'dod_membou',copy_subdomain(1) % lbper,      current_subdomain % lbper)
        call memory_copy(  mem_servi(1:2,servi),'LBOCH',      'dod_membou',copy_subdomain(1) % lboch,      current_subdomain % lboch)
        call memory_copy(  mem_servi(1:2,servi),'LNODB',      'dod_membou',copy_subdomain(1) % lnodb,      current_subdomain % lnodb)
        call memory_copy(  mem_servi(1:2,servi),'LTYPB',      'dod_membou',copy_subdomain(1) % ltypb,      current_subdomain % ltypb)
        call memory_copy(  mem_servi(1:2,servi),'LBOEL',      'dod_membou',copy_subdomain(1) % lboel,      current_subdomain % lboel)
        call memory_copy(  mem_servi(1:2,servi),'LELBO',      'dod_membou',copy_subdomain(1) % lelbo,      current_subdomain % lelbo)
        call memory_copy(  mem_servi(1:2,servi),'LSUBD_NBOUN','dod_membou',copy_subdomain(1) % lsubd_nboun,current_subdomain % lsubd_nboun)
        !
        ! Merge new boundary mesh
        !         
        lnods_subd => subdomain(isubd) % lnods
        ltype_subd => subdomain(isubd) % ltype
        lboch_subd => subdomain(isubd) % lboch
        lnodb_subd => subdomain(isubd) % lnodb
        ltypb_subd => subdomain(isubd) % ltypb
        lboel_subd => subdomain(isubd) % lboel
        lelbo_subd => subdomain(isubd) % lelbo
        lbper_subd => subdomain(isubd) % lbper

        iboun = nboun_ori
        do ielem = 1,nelem_subd
           if( current_subdomain % lsubd_nelem(ielem) < 0 .and. lmark(ielem) ) then
              ielty = ltype_subd(ielem)
              do iface = 1,nface(ielty)
                 if( faces(ielem) % l(mnodb+2,iface) /= 0 ) then  
                    jelem                     = faces(ielem) %l(mnodb+2,iface)
                    jface                     = faces(ielem) %l(mnodb+1,iface)
                    jelty                     = current_subdomain % ltype(jelem)
                    pnodb                     = nnodg(jface,jelty) 
                    iboun                     = iboun + 1                               ! Global boundary number
                    lelbo_subd(iboun)         = jelem                                   ! Element on the free side
                    call fintyp(ndimb,nnodg(jface,jelty),jblty)                         ! Find boundary type JBLTY
                    ltypb_subd(iboun)         = jblty                                   ! Boundary type
                    lboch_subd(iboun)         = BOEXT                                   ! This is an extension 
                    lbper_subd(iboun)         = 0                                       ! This is a new boundary 
                    jsubd                     = current_subdomain % lsubd_nelem(ielem)                    
                    current_subdomain % lsubd_nboun(iboun) = abs(jsubd)
                    ! 
                    ! Take the inverse of the face as the boundary was generated from inside the hole!
                    !
                    !do inodb = 1,pnodb     !con hexas ya esta bien como queda el sentido
                    !   lnodb_subd(inodb,iboun) = faces(ielem) %l(inodb,iface)
                    !end do
                    !Al revesssssss...aunque no me cuandra muchoooooo...
                    if( pnodb == 4 ) then
                       do inodb = 1,pnodb     !con hexas ya esta bien como queda el sentido
                          lnodb_subd(inodb,iboun) = faces(ielem) %l(inodb,iface)
                       end do
                    else if( pnodb == 3 ) then  !con triangulos-tetras giro el sentido de la numeracion de las caras para que la normal apunte hacia fuera del subdominio
                       lnodb_subd(1,iboun) = faces(ielem) %l(1,iface)
                       ii = 0 
                       do inodb=2,pnodb
                          lnodb_subd(inodb,iboun)= faces(ielem) %l(pnodb-ii,iface)
                          ii = ii + 1
                       end do
                    else  !en 2d
                       do inodb = 1,pnodb
                          lnodb_subd(inodb,iboun) = faces(ielem) %l(inodb,iface)
                       end do
                    end if
                    ninve = 0_ip
                    call runend('SHOULD USE LELBO AND NOT LBOEL')
                    call elmgeo_bounor(&
                         1_ip,1_ip,ndime,pnodb,mnode,lnodb_subd(1:pnodb,iboun),ltypb_subd(iboun:),&
                         lelbo_subd(iboun),ltype_subd,lnods_subd,nnode(1:),subdomain(isubd) % coord,&
                         ninve,dummr)
                    if( ninve /= 0_ip ) then
                       do inodb = 1,pnodb
                          lnodb_aux(inodb) = lnodb_subd(inodb,iboun)
                          lboel_aux(inodb) = lboel_subd(inodb,iboun)
                       end do
                       ii = 0
                       do inodb = 2,pnodb
                          lnodb_subd(inodb,iboun) = lnodb_aux(pnodb-ii)
                          lboel_subd(inodb,iboun) = lboel_aux(pnodb-ii)
                          ii = ii + 1
                       end do
                    end if

if( current_subdomain % lnper(current_subdomain % lnodb(1,iboun)) == 271 .or. current_subdomain % lnper(current_subdomain % lnodb(2,iboun)) == 271 ) then
print*,'MERDEUEIWDUYEOWUFOIWFUOIEWUFOIEWF'
end if
 !do isubd = 1,nsubd
 !    current_subdomain => subdomain(isubd)
 !    do iboun = 1,current_subdomain % nboun
 !       iboun_global = current_subdomain % lbper(iboun)
 !       if( current_subdomain % lboch(iboun) == BOEXT ) then  
 !          if( current_subdomain % lnper(current_subdomain % lnodb(1,iboun)) == 271 .or. current_subdomain % lnper(current_subdomain % lnodb(2,iboun)) == 271 ) then
 !             print*,'CACACACACACACACACACA=',isubd
 !          end if
 !       end if
 !    end do
 ! end do



                 end if
              end do
           end if
        end do
        !
        ! Deallocate memory
        !
        call memory_deallo(mem_servi(1:2,servi),'CFAEL','dod_holcut_holeboundary',cfael)
        call memory_deallo(mem_servi(1:2,servi),'NNODG','dod_holcut_holeboundary',nnodg)
        do ielem = 1,nelem_subd          
           if( lmark(ielem) ) then
              call memory_deallo(mem_servi(1:2,servi),'FACES % L','dod_holcut_holeboundary',faces(ielem) % l)
           end if
        end do
        call memory_deallo(mem_servi(1:2,servi),'FACES','dod_holcut_holeboundary',faces)
        call memory_deallo(mem_servi(1:2,servi),'LMARK','dod_holcut_holeboundary',lmark)

     end if

  end do

end subroutine dod_holcut_holeboundary
