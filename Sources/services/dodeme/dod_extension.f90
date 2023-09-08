!-----------------------------------------------------------------------
!> @addtogroup Dodeme
!> @{
!> @file    dod_extension.f90
!> @author  Guillaume Houzeaux
!> @date    18/09/2012
!> @brief   Create 2D extensions
!> @details Create extension elements for 2D meshes
!> @} 
!-----------------------------------------------------------------------
subroutine dod_extension()
  use def_kintyp
  use def_parame
  use def_master
  use def_elmtyp
  use def_domain
  use def_dodeme
  use mod_memory
  use mod_kdtree
  use mod_graphs
  use mod_dod_extens
  implicit none
  integer(ip)          :: isubd,jsubd,ielem
  integer(ip)          :: kelem,knode,boun1
  integer(ip)          :: perm1,perm2,perm3,boun2
  integer(ip)          :: ipoin,jpoin,jnode
  integer(ip)          :: jnodb,pnodb,pblty
  integer(ip)          :: ipoi1,ipoi2,ista0,ista1  
  integer(ip)          :: ibopo,ntext,imatn,ipoiz
  integer(ip)          :: pext1,pext2,pext3,pext4
  integer(ip)          :: inodb,iboun,nedge_ipoin
  integer(ip)          :: inext,mnext,nnext,npoin_subd
  integer(ip)          :: bandw,mbpoi,knodb,ipopo
  real(rp)             :: profi   
  real(rp)             :: time1,q,q1,q2,q3,dummy(3)
  integer(ip), pointer :: pbopo_subd(:)
  integer(ip), pointer :: lbopo_subd(:)
  integer(ip), pointer :: lnper_subd(:)
  integer(ip), pointer :: lbper_subd(:)
  integer(ip), pointer :: lboch_subd(:)
  integer(ip), pointer :: ltypb_subd(:)
  integer(ip), pointer :: lnodb_subd(:,:)
  integer(ip), pointer :: lnnob_global(:)
  integer(ip), pointer :: ppopo_on_boundaries(:)
  integer(ip), pointer :: lpopo_on_boundaries(:)
  integer(ip)          :: ipoin_global
  integer(ip)          :: iboun_global
  integer(ip)          :: ibopo_global
  integer(ip)          :: nboun_subd
  integer(ip)          :: ipoin_subd
  integer(ip)          :: iedge_global,nedge_global2,ii,iedge

  !------------------------------------------------------
  !
  ! Copy local arrays to global arrays
  ! LNODB_GLOBAL <= SUBDOMAIN(ISUBD) % LNODB
  ! LTYPB_GLOBAL <= SUBDOMAIN(ISUBD) % LTYPB
  ! LBOCH_GLOBAL <= SUBDOMAIN(ISUBD) % LBOCH
  ! PBOPO_GLOBAL <= SUBDOMAIN(ISUBD) % PBOPO
  ! LBOPO_GLOBAL <= SUBDOMAIN(ISUBD) % LBOPO
  !
  !------------------------------------------------------
  nullify(lnodb_global)
  nullify(lboch_global)
  nullify(ltypb_global)
  nullify(pbopo_global)
  nullify(lboel_global)

  nullify(pbopo_subd)
  nullify(lbopo_subd)
  nullify(lnper_subd)
  nullify(lbper_subd)
  nullify(lboch_subd)
  nullify(ltypb_subd)
  nullify(lnodb_subd)


  nboun_global = 0
  npoin_global = npoin
  nbopo_global = 0

  do isubd = 1,nsubd
     nboun_global = nboun_global + subdomain(isubd) % nboun
     nbopo_global = nbopo_global + size(subdomain(isubd) % lbopo) 
  end do
  call memory_alloca(mem_servi(1:2,servi),'LNODB_GLOBAL','dod_extension',lnodb_global,mnodb,nboun_global)
  call memory_alloca(mem_servi(1:2,servi),'LBOCH_GLOBAL','dod_extension',lboch_global,nboun_global)
  call memory_alloca(mem_servi(1:2,servi),'LTYPB_GLOBAL','dod_extension',ltypb_global,nboun_global)
  call memory_alloca(mem_servi(1:2,servi),'PBOPO_GLOBAL','dod_extension',pbopo_global,npoin_global+1_ip)
  call memory_alloca(mem_servi(1:2,servi),'LBOEL_GLOBAL','dod_extension',lboel_global,nboun_global)
  !
  ! Compute LNODB_GLOBAL, LTYPB_GLOBAL, LBOCH_GLOBAL
  !
  iboun_global = 0
  do isubd = 1,nsubd
     current_subdomain => subdomain(isubd)
     nboun_subd        =  current_subdomain % nboun
     npoin_subd        =  current_subdomain % npoin
     lnper_subd        => current_subdomain % lnper
     lboch_subd        => current_subdomain % lboch
     ltypb_subd        => current_subdomain % ltypb
     lnodb_subd        => current_subdomain % lnodb
     do iboun = 1,nboun_subd
        pnodb                      = nnode(abs(current_subdomain % ltypb(iboun)))
        iboun_global               = iboun_global + 1
        lboch_global(iboun_global) = lboch_subd(iboun)
        ltypb_global(iboun_global) = ltypb_subd(iboun)
        do inodb = 1,pnodb
           ipoin                            = lnodb_subd(inodb,iboun) 
           lnodb_global(inodb,iboun_global) = lnper_subd(ipoin)
        end do
        lboel_global(iboun_global) = current_subdomain % leper(current_subdomain % lelbo(iboun) )
     end do
  end do
  !
  ! Compute PBOPO_GLOBAL
  !
  nullify(nbpoi_global)
  call memory_alloca(mem_servi(1:2,servi),'NBPOI_GLOBAL','dod_extension',nbpoi_global,npoin_global)
  nbopo_global = 0
  do isubd = 1,nsubd
     current_subdomain => subdomain(isubd)
     pbopo_subd        => current_subdomain % pbopo
     npoin_subd        =  current_subdomain % npoin
     lnper_subd        => current_subdomain % lnper
     do ipoin = 1,npoin_subd
        !if( current_subdomain % lsubd_npoin(ipoin) > 0 ) then  !!!OJOOOO
        ipoin_global               = lnper_subd(ipoin)
        nbpoi_global(ipoin_global) = pbopo_subd(ipoin+1) - pbopo_subd(ipoin) 
        nbopo_global               = nbopo_global + nbpoi_global(ipoin_global)
        !end if
     end do
  end do

  pbopo_global(1) = 1
  do ipoin_global = 1,npoin_global
     pbopo_global(ipoin_global+1) = pbopo_global(ipoin_global) + nbpoi_global(ipoin_global)
  end do
  call memory_deallo(mem_servi(1:2,servi),'NBPOI_GLOBAL','dod_extension',nbpoi_global)
  !
  ! Fill in LBOPO_GLOBAL
  !
  call memory_alloca(mem_servi(1:2,servi),'LBOPO_GLOBAL','dod_extension',lbopo_global,nbopo_global)
  nboun_global = 0
  do isubd = 1,nsubd

     current_subdomain => subdomain(isubd)
     pbopo_subd        => current_subdomain % pbopo
     lbopo_subd        => current_subdomain % lbopo
     lnper_subd        => current_subdomain % lnper
     lbper_subd        => current_subdomain % lbper
     npoin_subd        =  current_subdomain % npoin
     nboun_subd        =  current_subdomain % nboun

     do ipoin = 1,npoin_subd
        !if( current_subdomain % lsubd_npoin(ipoin) > 0 ) then
        ipoin_global = lnper_subd(ipoin) 
        ibopo_global = pbopo_global(ipoin_global)
        do ibopo = pbopo_subd(ipoin),pbopo_subd(ipoin+1)-1
           lbopo_global(ibopo_global) = lbopo_subd(ibopo) + nboun_global
           ibopo_global               = ibopo_global + 1
        end do
        !end if
     end do
     nboun_global = nboun_global + nboun_subd

  end do

  !------------------------------------------------------
  !
  ! Construct edges: ONLY FOR 3D
  !
  !------------------------------------------------------

  if( ndime == 3 ) then
     !
     ! Create boundary graphs, only the lower part so that edges are not repeated
     !
     nullify(lnnob_global)
     nullify(ppopo_on_boundaries)
     nullify(lpopo_on_boundaries)
     call memory_alloca(mem_servi(1:2,servi),'LNNOB_GLOBAL','dod_extension',lnnob_global,nboun_global)
     do iboun = 1,nboun_global
        lnnob_global(iboun) = nnode(abs(ltypb_global(iboun)))
     end do
     call graphs_poipoi(&
          npoin,nboun_global,mnodb,lnodb_global,lnnob_global,ltypb_global,ppopo_on_boundaries,&
          lpopo_on_boundaries,bandw,profi,pbopo_global,lbopo_global,mbpoi,'EDGES')
     nedge_global2 = ppopo_on_boundaries(npoin_global+1) - 1
     call memory_alloca(mem_servi(1:2,servi),'LEDBO_GLOBAL','dod_extension',ledbo_global,4_ip,nedge_global2)
     !
     ! Count number of edges per node
     !  
     call memory_alloca(mem_servi(1:2,servi),'LEDGE_NPOIN_GLOBAL','dod_extension',ledge_npoin_global,npoin_global)
     call memgen(1_ip,npoin_global,0_ip)
     do ipoin = 1,npoin_global
        do ipopo = ppopo_on_boundaries(ipoin),ppopo_on_boundaries(ipoin+1)-1
           jpoin        = lpopo_on_boundaries(ipopo)
           gisca(ipoin) = gisca(ipoin) + 1
           gisca(jpoin) = gisca(jpoin) + 1
        end do
     end do
     do ipoin = 1,npoin_global
        if( gisca(ipoin) > 0 ) then
           call memory_alloca(mem_servi(1:2,servi),'LEDGE_NPOIN_GLOBAL(IPOIN) % L','dod_extension',ledge_npoin_global(ipoin) % l,gisca(ipoin))
           gisca(ipoin) = 0
        end if
     end do
     !
     ! LEDBO_GLOBAL(1,IEDGE_GLOBAL) = IBOUN,  LEDBO_GLOBAL(4,IEDGE_GLOBAL) = JBOUN  
     ! LEDBO_GLOBAL(2,IEDGE_GLOBAL) = INODB,  LEDBO_GLOBAL(5,IEDGE_GLOBAL) = INODB  
     ! LEDBO_GLOBAL(3,IEDGE_GLOBAL) = JNODB,  LEDBO_GLOBAL(6,IEDGE_GLOBAL) = JNODB  
     !
     iedge_global = 0
     do ipoin = 1,npoin_global
        do ipopo = ppopo_on_boundaries(ipoin),ppopo_on_boundaries(ipoin+1)-1
           jpoin                                       = lpopo_on_boundaries(ipopo)
           iedge_global                                = iedge_global + 1
           gisca(ipoin)                                = gisca(ipoin) + 1
           gisca(jpoin)                                = gisca(jpoin) + 1
           ledge_npoin_global(ipoin) % l(gisca(ipoin)) = iedge_global
           ledge_npoin_global(jpoin) % l(gisca(jpoin)) = iedge_global
           do ibopo = pbopo_global(ipoin),pbopo_global(ipoin+1)-1
              iboun = lbopo_global(ibopo)
              pnodb = lnnob_global(iboun)
              inodb = 0
              jnodb = 0
              do knodb = 1,pnodb
                 if( lnodb_global(knodb,iboun) == ipoin ) inodb = knodb
                 if( lnodb_global(knodb,iboun) == jpoin ) jnodb = knodb
              end do
              if( inodb /= 0 .and. jnodb /= 0 ) then
                 if( ledbo_global(1,iedge_global) == 0 ) then
                    ledbo_global(1,iedge_global) = iboun
                    ledbo_global(2,iedge_global) = inodb
                    ledbo_global(3,iedge_global) = jnodb
                 else if( ledbo_global(4,iedge_global) == 0 ) then
                    ledbo_global(4,iedge_global) = iboun
                 else
                    ipoi1=ipoin;ipoi2=jpoin
                    print*,'error in dod_extension'
                    print*,ipoin,jpoin,ledbo_global(1,iedge_global),ledbo_global(2,iedge_global),iboun
                    stop
                 end if

              end if
           end do

        end do
     end do
     call graphs_dealep(ppopo_on_boundaries,lpopo_on_boundaries)
     call memory_deallo(mem_servi(1:2,servi),'LNNOB_GLOBAL','dod_extension',lnnob_global)
     call memgen(3_ip,npoin_global,0_ip)
  end if


!!$  call memgen(1_ip,npoin,0_ip)
!!$  write(100,*) 'MESH boundary dimension 3 Elemtype Triangle Nnode  3'
!!$  write(100,*) 'coordinates'
!!$  isubd = 1
!!$  knodb = 0
!!$  current_subdomain => subdomain(isubd)
!!$  lnodb_subd        => current_subdomain % lnodb
!!$  ltypb_subd        => current_subdomain % ltypb
!!$  lnper_subd        => current_subdomain % lnper
!!$  npoin_subd        =  current_subdomain % npoin
!!$  nboun_subd        =  current_subdomain % nboun
!!$  do ipoin = 1,npoin_subd
!!$     if( current_subdomain % lsubd_npoin(ipoin) > 0 ) then
!!$        ipoin_global        = lnper_subd(ipoin)
!!$        knodb               = knodb + 1
!!$        gisca(ipoin_global) = knodb
!!$        write(100,*) knodb,coord(1:ndime,ipoin_global)
!!$     end if
!!$  end do
!!$  write(100,*) 'end coordinates'
!!$  write(100,*) 'elements'
!!$  do iboun = 1,nboun_subd
!!$     pnodb = nnode(ltypb_subd(iboun))
!!$     if( gisca(lnper_subd(lnodb_subd(1,iboun))) /= 0 ) &
!!$          write(100,'(10(1x,i7))') iboun,(gisca(lnper_subd(lnodb_subd(inodb,iboun))),inodb=1,pnodb)
!!$  end do
!!$  write(100,*) 'end elements'
!!$  call memgen(3_ip,npoin,0_ip)

  !------------------------------------------------------
  !
  ! Create 2D or 3D extension elements
  !
  !------------------------------------------------------
  if( ndime == 2 ) then
     call dod_extension_2d()
  else

     call dod_extension_3d()
  end if

end subroutine dod_extension
