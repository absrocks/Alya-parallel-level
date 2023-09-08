!-----------------------------------------------------------------------
!> @addtogroup Dodeme
!> @{
!> @file    dod_patch_boundaries.f90
!> @author  Guillaume Houzeaux
!> @date    18/09/2012
!> @brief   Detect outer boundary of patches
!> @details Find first a seed boundary, then mark the boundaries
!>          recursively
!> @} 
!-----------------------------------------------------------------------
subroutine dod_patch_boundaries()
  use def_parame
  use def_master
  use def_elmtyp
  use def_domain
  use def_dodeme
  use def_kermod
  use mod_kdtree
  use mod_graphs
  use mod_elmgeo
  implicit none
  integer(ip)          :: isubd,ipoin,inodb,pnodb,iboun,knodb
  integer(ip)          :: istack,nstack,kboun,jboun,jsubd
  integer(ip)          :: idime,ibobo,ninve,jboun_global,jnode
  integer(ip)          :: kfl_bouno,intij,intji,kfoun,jelem
  integer(ip)          :: nboun_subd,iboun_global,ksubd,dumm1(2)
  integer(ip)          :: npoin_subd
  integer(ip), pointer :: ltypb_subd(:)
  integer(ip), pointer :: lboch_subd(:)
  integer(ip), pointer :: lboel_subd(:,:)
  integer(ip), pointer :: lnodb_subd(:,:)
  integer(ip), pointer :: pbobo_subd(:)
  integer(ip), pointer :: lbobo_subd(:)
  integer(ip), pointer :: lnods_subd(:,:)
  integer(ip), pointer :: ltype_subd(:)
  real(rp),    pointer :: bouno_subd(:,:)
  integer(ip), pointer :: lstack(:)
  real(rp),    pointer :: coord_subd(:,:)
  logical(lg)          :: insid,add_nodes
  real(rp)             :: shape_tmp(64)       
  real(rp)             :: deriv_tmp(192)      

  integer(ip)          :: pnodb_neig,jnodb,ifoun,jpoin
  real(rp)             :: coor1(3),coor2(3),dista,xieta(3)
  real(rp)             :: xcoor(3)
  real(rp)             :: plapo(3),toler,bocod(ndime,mnodb)
  real(rp)             :: bouno_tmp(3)
  real(rp)             :: neighbor_subdomain_hole_box(3,2)

  nullify(lstack)

  do isubd = 1,nsubd
     do jsubd = 1,nsubd
        intij = intyp_dod(isubd,jsubd)
        intji = intyp_dod(jsubd,isubd)           
        if( intij /= 0 .and. intji /= 0 ) then
           if( ictop_dod(intij) == DOD_PATCH .and. ictop_dod(intji) == DOD_HOLED_PATCH ) then

              ictop_dod(intji) = DOD_PRESCRIBED
              ictop_dod(intij) = DOD_PRESCRIBED
              ipres_dod(isubd) = 1
              ipres_dod(jsubd) = 1
              ipatc_dod(isubd) = 0
              ipatc_dod(jsubd) = 0

              current_subdomain  => subdomain(isubd)
              neighbor_subdomain => subdomain(jsubd)
              !
              ! Bounding box of JSUBD hole with ISUBD
              !
              neighbor_subdomain_hole_box(1:3,1) =  1e9_rp
              neighbor_subdomain_hole_box(1:3,2) = -1e9_rp
              do jelem = 1,neighbor_subdomain % nelem
                 if( neighbor_subdomain % lsubd_nelem(jelem) == -isubd ) then
                    do jnode = 1,nnode(abs(neighbor_subdomain % ltype(jelem)))
                       jpoin = neighbor_subdomain % lnods(jnode,jelem)
                       neighbor_subdomain_hole_box(1:ndime,1) = min( neighbor_subdomain_hole_box(1:ndime,1) , neighbor_subdomain % coord(1:ndime,jpoin) ) 
                       neighbor_subdomain_hole_box(1:ndime,2) = max( neighbor_subdomain_hole_box(1:ndime,2) , neighbor_subdomain % coord(1:ndime,jpoin) ) 
                    end do
                 end if
              end do
              !
              ! Check points that are inside JUSBD embedding box
              !
              call memgen(1_ip,current_subdomain % npoin,0_ip)
              do ipoin = 1,current_subdomain % npoin 
                 do idime = 1,ndime
                    if(  current_subdomain % coord(idime,ipoin) - neighbor_subdomain_hole_box(idime,1) < -zeror .or. &
                         current_subdomain % coord(idime,ipoin) - neighbor_subdomain_hole_box(idime,2) >  zeror ) then
                       gisca(ipoin) = 1
                    end if
                 end do  
              end do

              do iboun = 1,current_subdomain % nboun 
                 pnodb = abs(nnode(current_subdomain % ltypb(iboun)))
                 knodb = 0 
                 do inodb = 1,pnodb
                    ipoin = current_subdomain % lnodb(inodb,iboun)
                    if( gisca(ipoin) == 0 ) then
                       call elsest(&
                            2_ip,jsubd,ielse,mnode,ndime,neighbor_subdomain % npoin,neighbor_subdomain % nelem,&
                            nnode(1:),neighbor_subdomain % lnods,neighbor_subdomain % ltype,&
                            ltopo,neighbor_subdomain % coord,current_subdomain % coord(1,ipoin),relse,jelem,&
                            shape_tmp,deriv_tmp,xieta,dumm1)
                       if( jelem > 0 ) then
                          if( neighbor_subdomain % lsubd_nelem(jelem) == -isubd ) then
                             gisca(ipoin) = 2
                             knodb = knodb + 1
                          else
                             gisca(ipoin) = 1
                          end if
                       else
                          gisca(ipoin) = 1
                       end if
                    else if( gisca(ipoin) == 2 ) then
                       knodb = knodb + 1
                    end if
                 end do
                 if( knodb > 0 ) then
                    !current_subdomain % lboch(iboun)       = BOEXT
                    iboun_global                           = current_subdomain % lbper(iboun)
                    if( current_subdomain % lnper(current_subdomain % lnodb(1,iboun)) == 271 .or. current_subdomain % lnper(current_subdomain % lnodb(2,iboun)) == 271 ) then
                       print*,'PIPIPIPI=',isubd,jsubd
                    end if
                    prescribed_boundaries(iboun_global)    = jsubd
                    current_subdomain % lsubd_nboun(iboun) = jsubd
                 end if
              end do   
              call memgen(3_ip,current_subdomain % npoin,0_ip)

           end if
        end if
     end do
  end do

  !----------------------------------------------------------------------
  !
  ! Mark exterior boundaries: SUBDOMAIN(ISUBD) % LBOCH(IBOUN) = BOEXT
  !
  !----------------------------------------------------------------------

  do isubd = 1,nsubd

     if( ipatc_dod(isubd) == 1 ) then

        current_subdomain => subdomain(isubd)
        npoin_subd        =  current_subdomain % npoin
        nboun_subd        =  current_subdomain % nboun
        ltypb_subd        => current_subdomain % ltypb
        lboch_subd        => current_subdomain % lboch
        lnodb_subd        => current_subdomain % lnodb
        lboel_subd        => current_subdomain % lboel
        bouno_subd        => current_subdomain % bouno
        pbobo_subd        => current_subdomain % pbobo
        lbobo_subd        => current_subdomain % lbobo
        lnods_subd        => current_subdomain % lnods
        ltype_subd        => current_subdomain % ltype
        coord_subd        => current_subdomain % coord

        iboun = 0
        do while ( iboun < nboun_subd )
           iboun = iboun + 1
           if( lboch_subd(iboun) /= BOHOL ) then
              pnodb = nnode(abs(ltypb_subd(iboun)))
              inodb = 0
              do while( inodb < pnodb )
                 inodb = inodb + 1
                 ipoin = lnodb_subd(inodb,iboun)
                 !
                 ! Check if ipoin is on the embedding box: this is the seed
                 !
                 do idime = 1,ndime
                    if(  abs( coord_subd(idime,ipoin) - current_subdomain % embox(idime,1) ) < zeror .or. &
                         abs( coord_subd(idime,ipoin) - current_subdomain % embox(idime,2) ) < zeror ) then
                       inodb = mnodb + 1
                       kboun = iboun
                       iboun = nboun_subd
                    end if
                 end do
              end do
           end if
        end do

        if( inodb == mnodb + 1 ) then
           !
           ! We found the seed boundary KBOUN
           !
           call memgen(1_ip,nboun_subd,0_ip)              
           call memory_alloca(mem_servi(1:2,servi),'LSTACK','dod_patch_boundaries',lstack,nboun_subd)
           !
           ! Mark boundaries recursively
           !
           iboun        = kboun
           nstack       = 1
           lstack(1)    = iboun
           gisca(iboun) = 1
           istack       = 0 

           do 
              if( istack == nstack ) exit 
              istack = istack+1   
              iboun  = lstack(istack)
              do ibobo = pbobo_subd(iboun),pbobo_subd(iboun+1)-1
                 jboun = lbobo_subd(ibobo)
                 if( gisca(jboun) == 0 .and. lboch_subd(iboun) /= BOHOL ) then
                    gisca(jboun)   = 1
                    nstack         = nstack + 1
                    lstack(nstack) = jboun
                 end if
              end do
           end do
           !
           ! Marked boundaries are of extension type
           !
           do iboun = 1,nboun_subd
              if( gisca(iboun) == 1 .and. lboch_subd(iboun) /= BOHOL ) then                 
                 current_subdomain % lboch(iboun) = BOEXT
                 pnodb = nnode(abs(ltypb_subd(iboun)))
                 do inodb = 1,pnodb
                    ipoin = lnodb_subd(inodb,iboun)
                    current_subdomain % lsubd_npoin(ipoin) = -9999
                 end do
              end if
           end do
        end if

        call memory_deallo(mem_servi(1:2,servi),'LSTACK','dod_patch_boundaries',lstack)
        call memgen(3_ip,nboun_subd,0_ip)

     end if
  end do

  !----------------------------------------------------------------------
  !
  ! PATCH: Compute boundary normals
  !
  !----------------------------------------------------------------------

  do isubd = 1,nsubd

     kfl_bouno = 0
     do jsubd = 1,nsubd
        intij = intyp_dod(isubd,jsubd)
        intji = intyp_dod(jsubd,isubd)           
        if( intij /= 0 ) then
           if( ictop_dod(intij) == DOD_PATCH .or. ictop_dod(intij) == DOD_HOLED_PATCH ) then
              kfl_bouno = 1
           end if
        end if
        if( intij /= 0 .and. intji /= 0 ) then
           if( ( ictop_dod(intji) == DOD_PATCH .or. ictop_dod(intji) == DOD_HOLED_PATCH ) .and. ictop_dod(intij) == DOD_PRESCRIBED ) then
              kfl_bouno = 1
           end if
        end if
     end do
     if( kfl_bouno == 1 ) then
        current_subdomain => subdomain(isubd)
        nboun_subd        =  current_subdomain % nboun
        ltypb_subd        => current_subdomain % ltypb
        lnodb_subd        => current_subdomain % lnodb
        lboel_subd        => current_subdomain % lboel
        lnods_subd        => current_subdomain % lnods
        ltype_subd        => current_subdomain % ltype
        coord_subd        => current_subdomain % coord
        call memory_alloca(mem_servi(1:2,servi),'BOUNO_DOD','dod_memall',current_subdomain % bouno,ndime,nboun_subd)
        bouno_subd        => current_subdomain % bouno
        call runend('SHOULD USE LELBO AND NOT LBOEL')
        call elmgeo_bounor(&
             1_ip,nboun_subd,ndime,mnodb,mnode,lnodb_subd,ltypb_subd,lboel_subd,&
             ltype_subd,lnods_subd,nnode(1:),coord_subd,ninve,bouno_subd)
     end if
  end do

  !----------------------------------------------------------------------
  !
  ! PATCH-PATCH
  ! Idea: for patch extension, take the inverse normal to the interface
  !
  !----------------------------------------------------------------------

  do isubd = 1,nsubd

     if( ipatc_dod(isubd) == 1 ) then

        current_subdomain => subdomain(isubd)
        npoin_subd        =  current_subdomain % npoin
        nboun_subd        =  current_subdomain % nboun
        ltypb_subd        => current_subdomain % ltypb
        lnodb_subd        => current_subdomain % lnodb
        lboel_subd        => current_subdomain % lboel
        bouno_subd        => current_subdomain % bouno
        lnods_subd        => current_subdomain % lnods
        ltype_subd        => current_subdomain % ltype
        coord_subd        => current_subdomain % coord
        !
        ! Determine who is my neighbor? Use patch relation and background embedding box
        !           
        do jsubd = 1,nsubd
           intij = intyp_dod(isubd,jsubd)
           intji = intyp_dod(jsubd,isubd)

           if( intij /= 0 ) then
              if( ictop_dod(intij) == DOD_PATCH .or. ictop_dod(intij) == DOD_HOLED_PATCH ) then

                 if( ictop_dod(intji) == DOD_CHIMERA ) then
                    !
                    ! PATCH-CHIMERA interface
                    ! -----------------------
                    !
                    ! JSUBD is Chimera
                    !
                    neighbor_subdomain => subdomain(jsubd)
                    do ipoin = 1,npoin_subd
                       if( current_subdomain % lsubd_npoin(ipoin) == -9999 ) then
                          insid = .true.
                          dimensions: do idime = 1,ndime
                             if(  coord_subd(idime,ipoin) < neighbor_subdomain % embox(idime,1) .or.&
                                  coord_subd(idime,ipoin) > neighbor_subdomain % embox(idime,2) ) then
                                insid = .false.
                                exit dimensions
                             end if
                          end do dimensions
                          if( insid ) current_subdomain % lsubd_npoin(ipoin) = jsubd
                       end if
                    end do

                 else
                    !
                    ! PATCH-PATCH or PATCH-PRESCRIBED interface
                    ! -----------------------------------------
                    !
                    ! JSUBD is also a patch or prescribed
                    ! Compute the normals to the boundaries
                    ! Check if they cross JUSBD
                    ! Mark them if this is the case
                    !
                    neighbor_subdomain => subdomain(jsubd)
                    !
                    ! If JSUBD has a prescribed interface with me, compute average boundary normal
                    !
                    intji = intyp_dod(jsubd,isubd)
                    if( ictop_dod(intji) == DOD_PRESCRIBED ) then
                       bouno_tmp = 0.0_rp
                       iboun     = 0
                       do jboun = 1,neighbor_subdomain % nboun
                          jboun_global = neighbor_subdomain % lbper(jboun)      
                          if( jboun_global == 0 ) then
                             ksubd = neighbor_subdomain % lsubd_nboun(jboun) ! This is a new boundary created from a hole in dod_holcut_holeboundary
                          else
                             ksubd = prescribed_boundaries(jboun_global)     ! This is an original boundary
                          end if
                          if( ksubd == isubd ) then
                             iboun = iboun + 1
                             bouno_tmp(1:ndime) = bouno_tmp(1:ndime) + neighbor_subdomain % bouno(1:ndime,jboun)
                          end if
                       end do
                       if( iboun == 0 ) then
                          call runend('DOD_PATCH_BOUNDARIES: WRONG BOUNDARY')
                       else
                          bouno_tmp(1:ndime) = - bouno_tmp(1:ndime) / real(iboun,rp)
                       end if
                    end if

                    do iboun = 1,nboun_subd
                       !
                       ! For each boundary IBOUN, check if we intersect with neighbor
                       ! in the direction of the normal
                       !
                       if( current_subdomain % lboch(iboun) == BOEXT ) then 
                          coor1 = 0.0_rp
                          coor2 = 0.0_rp
                          pnodb = nnode(abs(ltypb_subd(iboun)))
                          do inodb = 1,pnodb
                             ipoin = lnodb_subd(inodb,iboun)
                             do idime = 1,ndime
                                coor1(idime) = coor1(idime) + coord_subd(idime,ipoin)
                             end do
                          end do
                          do idime = 1,ndime
                             coor1(idime) = coor1(idime) / real(pnodb,rp)
                          end do
                          ifoun = 0
                          jboun = 0
                          do while ( jboun < neighbor_subdomain % nboun .and. ifoun == 0 ) 
                             jboun      = jboun + 1
                             pnodb_neig = abs(nnode(neighbor_subdomain % ltypb(jboun)))
                             xcoor      = 0.0_rp
                             do jnodb = 1,pnodb_neig
                                jpoin = neighbor_subdomain % lnodb(jnodb,jboun)
                                do idime = 1,ndime                             
                                   bocod(idime,jnodb) = neighbor_subdomain % coord(idime,jpoin)
                                   xcoor(idime)       = xcoor(idime) + bocod(idime,jnodb)
                                end do
                             end do
                             xcoor = xcoor / real(pnodb_neig,rp)
                             dista = 0.0_rp
                             do idime = 1,ndime
                                dista = dista + ( coor1(idime)-xcoor(idime) ) ** 2
                             end do
                             dista = sqrt(dista)
                             if( ictop_dod(intji) == DOD_PATCH .or. ictop_dod(intji) == DOD_HOLED_PATCH ) then
                                bouno_tmp(1:ndime) = bouno_subd(1:ndime,iboun)
                             end if
                             do idime =1,ndime
                                coor2(idime) = coor1(idime) + bouno_tmp(idime) * 1.5_rp * dista
                             end do
                             toler = 1.0e-2_rp
                             call elmgeo_segfac(ndime,pnodb_neig,bocod,coor1,coor2,ifoun,plapo,toler)
                             !
                             ! Check that I do not cross my own boundary
                             !
                             if( ifoun == 1 .and. ictop_dod(intji) == DOD_PRESCRIBED ) then
                                kfoun = 0
                                kboun = 0
                                coor2 = plapo
                                plapo = coor2 - coor1
                                coor2 = coor1 + 2.0_rp * plapo
                                do while ( kboun < nboun_subd .and. kfoun == 0 ) 
                                   kboun = kboun + 1
                                   if( kboun /= iboun ) then
                                      do inodb = 1,pnodb
                                         ipoin = current_subdomain % lnodb(inodb,kboun)
                                         do idime = 1,ndime                             
                                            bocod(idime,inodb) = current_subdomain % coord(idime,ipoin)
                                         end do
                                      end do
                                      call elmgeo_segfac(ndime,pnodb_neig,bocod,coor1,coor2,kfoun,plapo,toler)
                                   end if
                                end do
                                if( kfoun == 1 ) ifoun = 0
                             end if
                          end do
                          if( ifoun /= 0 ) then
                             !
                             ! Nodes of IBOUN are interface nodes with JSUBD
                             !
                             add_nodes = .false.
                             if( ictop_dod(intji) == DOD_PRESCRIBED ) then
                                !print*,neighbor_subdomain % lnper( neighbor_subdomain % lnodb(:,jboun) ) 
                                !print*,neighbor_subdomain % lboch(jboun),neighbor_subdomain % lsubd_nboun(jboun)
                                !stop 
                                if( neighbor_subdomain % lsubd_nboun(jboun) == isubd ) then
                                   add_nodes = .true.
                                end if
                             else
                                add_nodes = .true.
                             end if

                             do inodb = 1,pnodb
                                ipoin = lnodb_subd(inodb,iboun)
                                current_subdomain % lsubd_npoin(ipoin) = jsubd
                             end do
                          end if
                       end if
                    end do
                    !
                    ! IBOUN is of extension type if all its nodes are extension nodes
                    !
                    do iboun = 1,nboun_subd
                       if( current_subdomain % lboch(iboun) == BOEXT ) then                                 
                          pnodb = nnode(abs(ltypb_subd(iboun)))
                          knodb = 0
                          do inodb = 1,pnodb
                             ipoin = lnodb_subd(inodb,iboun)
                             if( current_subdomain % lsubd_npoin(ipoin) > 0 ) knodb = knodb + 1
                          end do
                          if( knodb /= pnodb ) current_subdomain % lboch(iboun) = BOFEM
                       end if
                    end do
                 end if
              end if
           end if
        end do
        do ipoin = 1,npoin_subd
           if( current_subdomain % lsubd_npoin(ipoin) == -9999 ) then
              current_subdomain % lsubd_npoin(ipoin) = 0
              !call runend('DOD_PATCH_BOUNDARIES: A PATCH BOUNDARY NODE IS OUT OF ITS BACKGROUND')
           end if
        end do
     end if
  end do
  !
  ! Deallocate BOUNO of necessary
  !
  do isubd = 1,nsubd
     call memory_deallo(mem_servi(1:2,servi),'BOUNO_DOD','dod_memall',subdomain(isubd) % bouno)
  end do

!!$  write(100,*) 'GiD Post Results File 1.0'
!!$  write(100,*) ' '
!!$  write(100,*) 'Result CACAC ALYA ',igene,' Scalar OnNodes'
!!$  write(100,*) 'ComponentNames CACAC'
!!$  write(100,*) 'Values'
!!$  do isubd = 1,nsubd
!!$     current_subdomain => subdomain(isubd)
!!$     do ipoin = 1,current_subdomain % npoin
!!$        write(100,*) current_subdomain % lnper(ipoin),current_subdomain % lsubd_npoin(ipoin)
!!$     end do
!!$  end do
!!$  write(100,*) 'End values'
!!$  print*,'dpodpqodpq'

end subroutine dod_patch_boundaries


