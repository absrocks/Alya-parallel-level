!-----------------------------------------------------------------------
!> @addtogroup Dodeme
!> @{
!> @file    dod_candidate_nodes.f90
!> @author  Guillaume Houzeaux
!> @date    18/09/2012
!> @brief   Look for candidate nodes for extension elements
!> @details Find the candidate nodes for the extension elements
!>          of the fringe nodes
!>          Eleminate candidates on the wrong side of the interface
!> @} 
!-----------------------------------------------------------------------
subroutine dod_candidate_nodes()
  use def_parame
  use def_master
  use def_elmtyp
  use def_domain
  use def_dodeme
  use def_kermod
  use mod_memory
  use mod_graphs
  use mod_kdtree
  implicit none 
  integer(ip)          :: ipoin,jpoin,istack,nlaye,nstack,dummi
  integer(ip)          :: jnode,jelem,isubd,jsubd,inode,jpoin_global
  integer(ip)          :: izdom,idime,kpoin,npoin_subd,ipoin_global
  integer(ip)          :: inodb,iboun,pnodb,dumm1(2),intij,intji
  real(rp)             :: dimin,dista,xieta(3),dummr(3)
  real(rp)             :: coori(3),coorj(3),chkdi
  real(rp)             :: shape_tmp(64)       
  real(rp)             :: deriv_tmp(192)      
  integer(ip), pointer :: lmark(:)
  integer(ip), pointer :: lstack(:)
  integer(ip), pointer :: ppopo_neig(:)
  integer(ip), pointer :: lpopo_neig(:)
  logical(lg)          :: lchek
  ! 
  ! Allocate memory for extension elements
  !
  call dod_memall(3_ip)

  nullify(lmark)
  nullify(lstack)
  nullify(ppopo_neig)
  nullify(lpopo_neig)
  !
  ! Initialization and memory allocation
  !
  call memory_alloca(mem_servi(1:2,servi),'LMARK' ,'dod_candidate_nodes',lmark, npoin)
  call memory_alloca(mem_servi(1:2,servi),'LSTACK','dod_candidate_nodes',lstack,npoin)

  do isubd = 1,nsubd

     current_subdomain => subdomain(isubd)
     npoin_subd        =  current_subdomain % npoin
     do ipoin = 1,npoin_subd
 
        if( current_subdomain % lsubd_npoin(ipoin) > 0 ) then
           jsubd              =  current_subdomain % lsubd_npoin(ipoin)
           neighbor_subdomain => subdomain(jsubd)
           ppopo_neig         => neighbor_subdomain % ppopo
           lpopo_neig         => neighbor_subdomain % lpopo
           nstack             =  0
           !
           ! IPOIN has a host element JELEM? It could have been computed in dod_holcut_marknodes()
           ! In this case, JELEM is in my neighbor local numbering
           !
           do idime = 1,ndime                        
              coori(idime) = current_subdomain % coord(idime,ipoin)      
           end do
           jelem = current_subdomain % host_element(ipoin)

           if( jelem <= 0 ) then        
              call elsest(&
                   2_ip,jsubd,ielse,mnode,ndime,neighbor_subdomain % npoin,neighbor_subdomain % nelem,&
                   nnode(1:),neighbor_subdomain % lnods,neighbor_subdomain % ltype,&
                   ltopo,neighbor_subdomain % coord,coori,relse,jelem,&
                   shape_tmp,deriv_tmp,xieta,dumm1)
           end if
          if( jelem > 0 ) then
              !
              ! IPOIN has a host element
              !
              do jnode = 1,neighbor_subdomain % lnnod(jelem)                 
                 jpoin          = neighbor_subdomain % lnods(jnode,jelem)
                 nstack         = nstack + 1
                 lstack(nstack) = jpoin
                 lmark(jpoin)   = 1
              end do

           else
              !
              ! IPOIN does not have a host element: look for nearest free node in JSUBD
              ! Idea: use elsest to look first for the nearest bin/quad
              !
              dimin = 1.0e12_rp
              do jpoin = 1,neighbor_subdomain % npoin
                 if( neighbor_subdomain % lsubd_npoin(jpoin) == DOD_FREE ) then
                    dista = 0.0_rp
                    do idime = 1,ndime
                       dista = dista + ( neighbor_subdomain % coord(idime,jpoin) - coori(idime) )**2
                    end do
                    if( dista < dimin ) then
                       kpoin = jpoin
                       dimin = dista
                    end if
                 end if
              end do

              do izdom = ppopo_neig(kpoin),ppopo_neig(kpoin+1)-1
                 jpoin          = lpopo_neig(izdom)
                 nstack         = nstack + 1
                 lstack(nstack) = jpoin
                 lmark(jpoin)   = 1
              end do
           end if
           !
           ! Mark nodes recursively:
           ! LMARK(IPOIN) ...... Layer number of IPOIN
           ! LSTACK(1:NSTACK) .. List of node with non-zero layer          
           !
           nlaye = 4_ip
           call graphs_recurl(&
                neighbor_subdomain % npoin,nlaye,neighbor_subdomain % ppopo,&
                neighbor_subdomain % lpopo,nstack,lmark,lstack)
           !
           ! Compute number of candidates: 
           ! CURRENT_SUBDOMAIN % EXTENSION(IPOIN) % NUMBER_CANDIDATES = JNODE
           ! Only Include free nodes in the list
           ! If we only want to consider free nodes as candidates, uncomment the double !!
           !
           jnode = 0
           do istack = 1,nstack
              jpoin = lstack(istack)
              !! if( neighbor_subdomain % lsubd_npoin(jpoin) == DOD_FREE ) jnode = jnode + 1
              jnode = jnode + 1
           end do
           current_subdomain % extension(ipoin) % number_candidates = jnode

           if( jnode == 0_ip ) then
              print*,'NO TENGO CANDIDATOS-dod_candidate_nodes',ipoin
              stop
           end if

           !
           ! Save candidate list: current_subdomain % extension(ipoin) % candidates 
           !
           call memory_alloca(mem_servi(1:2,servi),'GRADI_DOD(IPOIN)%LNODE','dod_candidate_nodes',&
                current_subdomain % extension(ipoin) % candidates,jnode)
           jnode = 0
           do istack = 1,nstack
              jpoin = lstack(istack)
              !!if( neighbor_subdomain % lsubd_npoin(jpoin) == DOD_FREE ) then
                 jnode = jnode + 1
                 current_subdomain % extension(ipoin) % candidates(jnode) = lstack(istack)

                 !!end if
              !if(current_subdomain % lnper(ipoin)==945 .and. neighbor_subdomain % lsubd_npoin(jpoin)/=DOD_FREE)
              lstack(istack) = 0
              lmark(jpoin)   = 0
           end do

        end if

     end do

  end do
  !
  ! Deallocate memory
  !
  call memory_deallo(mem_servi(1:2,servi),'LMARK', 'dod_candidate_nodes', lmark)
  call memory_deallo(mem_servi(1:2,servi),'LSTACK','dod_candidate_nodes',lstack)

  !----------------------------------------------------------------------
  !
  ! Detect nodes on boundaries of the domain, connected at least to 
  ! a BOFEM type boundary
  ! Then if both subdomains have prescribed boundaries, only connect
  ! a BOFEM-type node to another BOFEM-type node
  !
  !----------------------------------------------------------------------

!!$  call memgen(1_ip,npoin,0_ip)
!!$
!!$  do isubd = 1,nsubd   
!!$     current_subdomain => subdomain(isubd)
!!$     do iboun = 1,current_subdomain % nboun 
!!$       if( current_subdomain % lboch(iboun) == BOFEM ) then                                 
!!$           pnodb = nnode(abs(current_subdomain % ltypb(iboun)))
!!$           do inodb = 1,pnodb
!!$              ipoin               = current_subdomain % lnodb(inodb,iboun)
!!$              ipoin_global        = current_subdomain % lnper(ipoin)
!!$              gisca(ipoin_global) = 1
!!$           end do
!!$        end if
!!$     end do
!!$  end do
!!$
!!$  do isubd = 1,nsubd
!!$     current_subdomain => subdomain(isubd)
!!$     npoin_subd        =  current_subdomain % npoin
!!$     do ipoin = 1,npoin_subd
!!$        ipoin_global = current_subdomain % lnper(ipoin)
!!$        if( current_subdomain % lsubd_npoin(ipoin) > 0 .and. gisca(ipoin_global) /= 0 ) then
!!$           jsubd =  current_subdomain % lsubd_npoin(ipoin)        
!!$           intij = intyp_dod(isubd,jsubd)
!!$           intji = intyp_dod(jsubd,isubd)
!!$           lchek = .false.
!!$           if( intij /= 0 .and. intji /= 0 ) then
!!$              if( ictop_dod(intij) == DOD_PRESCRIBED .and. ictop_dod(intji) == DOD_PRESCRIBED ) then
!!$                 lchek = .true.
!!$              end if
!!$           end if
!!$           neighbor_subdomain => subdomain(jsubd)
!!$           inode              =  0
!!$           do jnode = 1,current_subdomain % extension(ipoin) % number_candidates
!!$              jpoin        = current_subdomain % extension(ipoin) % candidates(jnode)
!!$              jpoin_global = neighbor_subdomain % lnper(jpoin)
!!$              if( gisca(jpoin_global) == 0 .and. lchek ) then
!!$                 continue
!!$              else 
!!$                 inode = inode + 1 
!!$                 current_subdomain % extension(ipoin) % candidates(inode) = jpoin
!!$              end if
!!$           end do   
!!$           current_subdomain % extension(ipoin) % number_candidates = inode
!!$        end if
!!$     end do
!!$  end do
!!$
!!$  call memgen(3_ip,npoin,0_ip)

  !----------------------------------------------------------------------
  !
  ! Remove candidates on the wrong side of the interface. Recompute
  ! current_subdomain % extension(ipoin) % number_candidates 
  ! current_subdomain % extension(ipoin) % candidates(:)
  !
  !----------------------------------------------------------------------

  chkdi = 1.0e9_rp
  do isubd = 1,nsubd
     call memgen(1_ip,npoin,0_ip)
     current_subdomain => subdomain(isubd)
     npoin_subd        =  current_subdomain % npoin
     do ipoin = 1,npoin_subd
        if( current_subdomain % lsubd_npoin(ipoin) > 0 ) then
           jsubd              =  current_subdomain % lsubd_npoin(ipoin)
           neighbor_subdomain => subdomain(jsubd)
           inode              =  0

           do jnode = 1,current_subdomain % extension(ipoin) % number_candidates
              jpoin        = current_subdomain % extension(ipoin) % candidates(jnode)
              jpoin_global = neighbor_subdomain % lnper(jpoin)

              if( gisca(jpoin_global) == 0 ) then
                 do idime = 1,ndime
                    coorj(idime) = neighbor_subdomain % coord(idime,jpoin)
                 end do

                 call dpopar(&
                      1_ip,coorj,current_subdomain % npoin_bou,mnodb,&
                      current_subdomain % nboun,chkdi,current_subdomain % ltypb,current_subdomain % lnodb_bou,&
                      current_subdomain % coord_bou,dista,dummr,dummr,dummi,&
                      current_subdomain % fabox,current_subdomain % sabox,current_subdomain % blink,&
                      current_subdomain % stru2,current_subdomain % ldist,current_subdomain % lnele)
                 
                 if( dista >= 0.0_rp ) then 
                    gisca(jpoin_global) =  1 
                 else 
                    gisca(jpoin_global) = -1 ! Remove node from list
                 end if
              end if

              if( gisca(jpoin_global) == 1 ) then
                 inode = inode + 1 
                 current_subdomain % extension(ipoin) % candidates(inode) = jpoin
              end if
           end do
           current_subdomain % extension(ipoin) % number_candidates = inode
          
           if(inode ==0_ip)then
              print*,'NO TENGO CANDIDATOS-dod_candidate_nodes-222',current_subdomain % lnper(ipoin)
              stop
           end if
        end if     
     end do

     call memgen(3_ip,npoin,0_ip)

  end do

!!$     isubd=1
!!$     if( isubd == 1 ) then
!!$        do ipoin = 1,subdomain(isubd) % npoin
!!$           if( subdomain(isubd) % lnper(ipoin) == 293 ) print*,'d=',subdomain(2) % lnper(subdomain(isubd) % extension(ipoin) % candidates( 1:subdomain(isubd) % extension(ipoin) % number_candidates ))
!!$        end do
!!$     end if

end subroutine dod_candidate_nodes
