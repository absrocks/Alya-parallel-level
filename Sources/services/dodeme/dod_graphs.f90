!-----------------------------------------------------------------------
!> @addtogroup Dodeme
!> @{
!> @file    dod_graphs.f90
!> @author  Guillaume Houzeaux
!> @date    27/02/2013
!> @brief   Compute some graphs for subdomains
!> @details Compute the following graphs:
!>          ITASK = 1 ... At the beginning
!>                = 2 ... After hole cutting, recompute boundary graph
!> @} 
!-----------------------------------------------------------------------

subroutine dod_graphs(itask)

  use def_parame
  use def_master
  use def_elmtyp
  use def_domain
  use def_dodeme 
  use mod_graphs
  use mod_messages, only : livinf
  implicit none
  integer(ip),  intent(in) :: itask
  integer(ip)              :: isubd,mepoi_tmp,ielty
  integer(ip)              :: bandw,iboun,mbpoi
  real(rp)                 :: profi
  integer(ip),  pointer    :: cfael(:,:,:)
  integer(ip),  pointer    :: nnodg(:,:)

  nullify(cfael)
  nullify(nnodg)
  
  select case ( itask )

  case ( 1_ip ) 

     !-------------------------------------------------------------------
     !
     ! Element, node and boundary graphs
     !
     !-------------------------------------------------------------------

     if( INOTSLAVE ) then

        call livinf(0_ip,'CONSTRUCT BOUNDARY GRAPHS',0_ip)
        !
        ! Compute CFAEL AND NNODG
        !
        call memory_alloca(mem_servi(1:2,servi),'CFAEL','dod_graphs',cfael,mnodb,mface,nelty)
        call memory_alloca(mem_servi(1:2,servi),'NNODG','dod_graphs',nnodg,mface,nelty)     
        do ielty = iesta_dom,iesto_dom
           !if( lexis(ielty) /= 0 ) &
           !     call domfa2(&                     
           !     mnodb,nnodg(1,ielty),nface(ielty),ielty,cfael(1,1,ielty)) 
           call runend('DOD_GRAPHS: RECODE')
        end do

        do isubd = 1,nsubd
           current_subdomain => subdomain(isubd)
           !
           ! PELPO, LELPO 
           ! PPOPO, LPOPO
           !
           call graphs_poipoi(&
                current_subdomain % npoin,current_subdomain % nelem,mnode,current_subdomain % lnods,&
                current_subdomain % lnnod,current_subdomain % ltype,current_subdomain % ppopo,&
                current_subdomain % lpopo,bandw,profi,current_subdomain % pelpo,&
                current_subdomain % lelpo,mepoi_tmp)
           !
           ! PELEL, LELEL: only needed for hole cutting
           ! 
           if( ihole_dod(isubd) /= 0 ) then
              call runend('CALL TO graphs_eleele_faces HAS CHANGED')
              !call graphs_eleele_faces(&
              !     current_subdomain % nelem,mnode,mnodb,nelty,mepoi,mface,current_subdomain % lnods,&
              !     current_subdomain % lnnod,current_subdomain % ltype,nnodg,cfael,nface,&
              !     current_subdomain % pelpo,current_subdomain % lelpo,nedge,medge,&
              !     current_subdomain % pelel,current_subdomain % lelel)
           end if
           call memgen(1_ip,current_subdomain % nboun,0_ip)
           do iboun = 1,current_subdomain % nboun
              gisca(iboun) = abs(nnode(current_subdomain % ltypb(iboun)))
           end do
           !
           ! PBOPO, LBOPO: use to create extensions
           !
           call graphs_elepoi(&
                current_subdomain % npoin,current_subdomain % nboun,mnodb,current_subdomain % lnodb,&
                gisca,mbpoi,current_subdomain % pbopo,current_subdomain % lbopo)
           !
           ! PBOBO, LBOBO: only needed for patch automatic boundary
           !
           if( ipatc_dod(isubd) /= 0 ) then
              call graphs_eleele(&
                   current_subdomain % nboun,current_subdomain % npoin,mnodb,mbpoi,&
                   current_subdomain % lnodb,gisca,current_subdomain % pbopo,current_subdomain % lbopo,&
                   nedge,medge,current_subdomain % pbobo,&
                   current_subdomain % lbobo)
           end if
           call memgen(3_ip,current_subdomain % nboun,0_ip)
        end do

        call memory_deallo(mem_servi(1:2,servi),'CFAEL','dod_graphs',cfael)
        call memory_deallo(mem_servi(1:2,servi),'NNODG','par_elmgra',nnodg)     

     end if

  case ( 2_ip ) 

     !-------------------------------------------------------------------
     !
     ! Boundary graphs after hole cutting: PBOPO and LBOPO
     !
     !-------------------------------------------------------------------

     if( INOTSLAVE ) then

        call livinf(0_ip,'CONSTRUCT BOUNDARY GRAPHS',0_ip)
        !
        ! Compute CFAEL AND NNODG
        !
        mface = maxval(nface)
        call memory_alloca(mem_servi(1:2,servi),'CFAEL','dod_graphs',cfael,mnodb,mface,nelty)
        call memory_alloca(mem_servi(1:2,servi),'NNODG','par_elmgra',nnodg,mface,nelty)     
        do ielty = iesta_dom,iesto_dom
           if( lexis(ielty) /= 0 ) &
                call domfa2(&                     
                mnodb,nnodg(1,ielty),nface(ielty),ielty,cfael(1,1,ielty)) 
        end do

        do isubd = 1,nsubd

           if( ihole_dod(isubd) /= 0 ) then
              !
              ! PBOPO, LBOPO: use to create extensions
              !
              current_subdomain => subdomain(isubd)
              call memgen(1_ip,current_subdomain % nboun,0_ip)
              do iboun = 1,current_subdomain % nboun
                 gisca(iboun) = abs(nnode(current_subdomain % ltypb(iboun)))
              end do 
              call graphs_elepoi(&
                   current_subdomain % npoin,current_subdomain % nboun,mnodb,current_subdomain % lnodb,&
                   gisca,mbpoi,current_subdomain % pbopo,current_subdomain % lbopo)
              !
              ! This is necessary because an automatic patch could have been holed
              !
              if( ipatc_dod(isubd) /= 0 ) then
                 call graphs_eleele(&
                      current_subdomain % nboun,current_subdomain % npoin,mnodb,mbpoi,&
                      current_subdomain % lnodb,gisca,current_subdomain % pbopo,current_subdomain % lbopo,&
                      nedge,medge,current_subdomain % pbobo,&
                      current_subdomain % lbobo)
              end if
              call memgen(3_ip,current_subdomain % nboun,0_ip)
           end if

        end do

        call memory_deallo(mem_servi(1:2,servi),'CFAEL','dod_graphs',cfael)
        call memory_deallo(mem_servi(1:2,servi),'NNODG','par_elmgra',nnodg)     

     end if

  end select

end subroutine dod_graphs

