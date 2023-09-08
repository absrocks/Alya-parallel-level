!-----------------------------------------------------------------------
!> @addtogroup Dodeme
!> @{
!> @file    dod_memall.f90
!> @author  Guillaume Houzeaux
!> @date    18/09/2012
!> @brief   Allocate memory
!> @details Allocate memory
!> @} 
!-----------------------------------------------------------------------
subroutine dod_memall(itask)
  use def_parame
  use def_master
  use def_domain
  use def_kermod
  use def_dodeme
  use mod_memory
  use mod_kdtree
  implicit none
  integer(ip), intent(in) :: itask
  integer(ip)             :: isubd,npoin_subd,nelem_subd,nboun_subd
  integer(ip)             :: ipoin,dummi
  real(rp)                :: dummr

  select case ( itask ) 

  case ( 1_ip ) 
     !
     ! nullify pointers
     !
     nullify(ihole_dod)
     nullify(ipatc_dod)
     nullify(ipres_dod)
     nullify(lsubd_npoin)
     nullify(lsubd_nelem)
     nullify(lsubd_nboun)
     nullify(linvp_npoin)
     nullify(linvp_nelem)
     nullify(linvp_nboun)
     nullify(lmatn_dod)
     nullify(lpoiz_dod)
     nullify(prescribed_boundaries)
     !
     ! Subdomain arrays
     !
     call memory_alloca(mem_servi(1:2,servi),'IHOLE_DOD','dod_memall',ihole_dod,nsubd)
     call memory_alloca(mem_servi(1:2,servi),'IPATC_DOD','dod_memall',ipatc_dod,nsubd)
     call memory_alloca(mem_servi(1:2,servi),'IPRES_DOD','dod_memall',ipres_dod,nsubd)
     !
     ! Mesh arrays
     !
     call memory_alloca(mem_servi(1:2,servi),'LSUBD_NPOIN'          ,'dod_memall',lsubd_npoin,npoin)
     call memory_alloca(mem_servi(1:2,servi),'LSUBD_NELEM'          ,'dod_memall',lsubd_nelem,nelem)
     call memory_alloca(mem_servi(1:2,servi),'LSUBD_NBOUN'          ,'dod_memall',lsubd_nboun,nboun)
     call memory_alloca(mem_servi(1:2,servi),'LINVP_NPOIN'          ,'dod_memall',linvp_npoin,npoin)
     call memory_alloca(mem_servi(1:2,servi),'LINVP_NELEM'          ,'dod_memall',linvp_nelem,nelem)
     call memory_alloca(mem_servi(1:2,servi),'LINVP_NBOUN'          ,'dod_memall',linvp_nboun,nboun)
     call memory_alloca(mem_servi(1:2,servi),'LMATN_DOD'            ,'dod_memall',lmatn_dod,npoin)
     call memory_alloca(mem_servi(1:2,servi),'LPOIZ_DOD'            ,'dod_memall',lpoiz_dod,npoin)
     call memory_alloca(mem_servi(1:2,servi),'PRESCRIBED_BOUNDARIES','dod_memall',prescribed_boundaries,nboun)
     !
     ! Subdomain arrays
     !
     allocate( subdomain(nsubd) )
     do isubd = 1,nsubd
        current_subdomain => subdomain(isubd)
        nullify( current_subdomain % lsubd_npoin )
        nullify( current_subdomain % lsubd_nelem )
        nullify( current_subdomain % lsubd_nboun )
        nullify( current_subdomain % host_element )
        nullify( current_subdomain % lbper )
        nullify( current_subdomain % lnper )
        nullify( current_subdomain % leper )
        nullify( current_subdomain % coord )
        nullify( current_subdomain % lnods )
        nullify( current_subdomain % lelch )
        nullify( current_subdomain % lnnod )
        nullify( current_subdomain % ltype )
        nullify( current_subdomain % lboch )
        nullify( current_subdomain % lnodb )
        nullify( current_subdomain % ltypb )
        nullify( current_subdomain % lboel )
        nullify( current_subdomain % lelbo )
        nullify( current_subdomain % bouno )

        nullify( current_subdomain % pelel )     ! Graphs
        nullify( current_subdomain % lelel )     ! 
        nullify( current_subdomain % pbobo )     ! 
        nullify( current_subdomain % lbobo )     ! 
        nullify( current_subdomain % pelpo )     ! ...
        nullify( current_subdomain % lelpo )     ! 
        nullify( current_subdomain % ppopo )     ! 
        nullify( current_subdomain % lpopo )     ! 
        nullify( current_subdomain % pbopo )     ! 
        nullify( current_subdomain % lbopo )     ! End graphs

        nullify( current_subdomain % fabox )     ! Allocated in dod_kdtree
        nullify( current_subdomain % sabox )     ! ...
        nullify( current_subdomain % blink )     ! ...
        nullify( current_subdomain % stru2 )     ! ...
        nullify( current_subdomain % ldist )     ! ...
        nullify( current_subdomain % lnele )     ! ...
        nullify( current_subdomain % coord_bou ) ! ...
        nullify( current_subdomain % lnodb_bou ) ! ...

        nullify( intyp_dod )
     end do

     call memory_alloca(mem_servi(1:2,servi),'INTYP_DOD','dod_memall',intyp_dod,nsubd,nsubd)

  case ( 2_ip )

     do isubd = 1,nsubd
        current_subdomain => subdomain(isubd)     
        !
        ! Dimensions
        !
        npoin_subd = current_subdomain % npoin
        nelem_subd = current_subdomain % nelem
        nboun_subd = current_subdomain % nboun
        !
        ! Transmission conditions
        !
        call memory_alloca(mem_servi(1:2,servi),'LSUBD_NPOIN','dod_memall',current_subdomain % lsubd_npoin,npoin_subd)
        call memory_alloca(mem_servi(1:2,servi),'LSUBD_NELEM','dod_memall',current_subdomain % lsubd_nelem,nelem_subd)
        call memory_alloca(mem_servi(1:2,servi),'LSUBD_NBOUN','dod_memall',current_subdomain % lsubd_nboun,nboun_subd)
        !
        ! Permutation arrays
        !
        call memory_alloca(mem_servi(1:2,servi),'LNPER','dod_memall',current_subdomain % lnper,npoin_subd,'DO_NOT_INITIALIZE')
        call memory_alloca(mem_servi(1:2,servi),'LEPER','dod_memall',current_subdomain % leper,nelem_subd,'DO_NOT_INITIALIZE')
        call memory_alloca(mem_servi(1:2,servi),'LBPER','dod_memall',current_subdomain % lbper,nboun_subd,'DO_NOT_INITIALIZE')
        !
        ! Others
        !
        call memory_alloca(mem_servi(1:2,servi),'HOST_ELEMENTS','dod_memall',current_subdomain % host_element,npoin_subd)
        !
        ! Geometrical arrays
        !
        call memory_alloca(mem_servi(1:2,servi),'COORD','dod_memall',current_subdomain % coord,ndime,npoin_subd,     'DO_NOT_INITIALIZE')
        call memory_alloca(mem_servi(1:2,servi),'LNODS','dod_memall',current_subdomain % lnods,mnode,nelem_subd,     'DO_NOT_INITIALIZE')
        call memory_alloca(mem_servi(1:2,servi),'LELCH','dod_memall',current_subdomain % lelch,nelem_subd,           'DO_NOT_INITIALIZE')
        call memory_alloca(mem_servi(1:2,servi),'LNNOD','dod_memall',current_subdomain % lnnod,nelem_subd,           'DO_NOT_INITIALIZE')
        call memory_alloca(mem_servi(1:2,servi),'LTYPE','dod_memall',current_subdomain % ltype,nelem_subd,           'DO_NOT_INITIALIZE')
        call memory_alloca(mem_servi(1:2,servi),'LBOCH','dod_memall',current_subdomain % lboch,nboun_subd,           'DO_NOT_INITIALIZE')
        call memory_alloca(mem_servi(1:2,servi),'LNODB','dod_memall',current_subdomain % lnodb,mnodb,nboun_subd,     'DO_NOT_INITIALIZE')
        call memory_alloca(mem_servi(1:2,servi),'LTYPB','dod_memall',current_subdomain % ltypb,nboun_subd,           'DO_NOT_INITIALIZE')
        call memory_alloca(mem_servi(1:2,servi),'LBOEL','dod_memall',current_subdomain % lboel,mnodb,nboun_subd,     'DO_NOT_INITIALIZE')
        call memory_alloca(mem_servi(1:2,servi),'LELBO','dod_memall',current_subdomain % lelbo,nboun_subd,           'DO_NOT_INITIALIZE')

     end do

  case (-2_ip )

     do isubd = 1,nsubd
        current_subdomain => subdomain(isubd)     
        !
        ! Elsest
        !
        call elsest(&
             3_ip,isubd,ielse,mnode,ndime,current_subdomain % npoin,&
             current_subdomain % nelem,nnode(1:),current_subdomain % lnods,&
             current_subdomain % ltype,ltopo,current_subdomain % coord,dummr,&
             relse,dummi,shape_dod,deriv_dod,dummr,dummi)
        !
        ! Transmission conditions
        !
        call memory_deallo(mem_servi(1:2,servi),'LSUBD_NPOIN','dod_memall',current_subdomain % lsubd_npoin)
        call memory_deallo(mem_servi(1:2,servi),'LSUBD_NELEM','dod_memall',current_subdomain % lsubd_nelem)
        call memory_deallo(mem_servi(1:2,servi),'LSUBD_NBOUN','dod_memall',current_subdomain % lsubd_nboun)
        !
        ! Permutation arrays
        !
        call memory_deallo(mem_servi(1:2,servi),'LNPER','dod_memall',current_subdomain % lnper)
        call memory_deallo(mem_servi(1:2,servi),'LEPER','dod_memall',current_subdomain % leper)
        call memory_deallo(mem_servi(1:2,servi),'LBPER','dod_memall',current_subdomain % lbper)
        !
        ! Others
        !
        call memory_deallo(mem_servi(1:2,servi),'HOST_ELEMENTS','dod_memall',current_subdomain % host_element)
        !
        ! Geometrical arrays
        !
        call memory_deallo(mem_servi(1:2,servi),'COORD','dod_memall',current_subdomain % coord)
        call memory_deallo(mem_servi(1:2,servi),'LNODS','dod_memall',current_subdomain % lnods)
        call memory_deallo(mem_servi(1:2,servi),'LELCH','dod_memall',current_subdomain % lelch)
        call memory_deallo(mem_servi(1:2,servi),'LNNOD','dod_memall',current_subdomain % lnnod)
        call memory_deallo(mem_servi(1:2,servi),'LTYPE','dod_memall',current_subdomain % ltype)
        call memory_deallo(mem_servi(1:2,servi),'LBOCH','dod_memall',current_subdomain % lboch)
        call memory_deallo(mem_servi(1:2,servi),'LNODB','dod_memall',current_subdomain % lnodb)
        call memory_deallo(mem_servi(1:2,servi),'LTYPB','dod_memall',current_subdomain % ltypb)
        call memory_deallo(mem_servi(1:2,servi),'LBOEL','dod_memall',current_subdomain % lboel)        
        call memory_deallo(mem_servi(1:2,servi),'LELBO','dod_memall',current_subdomain % lelbo)        
        !
        ! Graphs
        !
        call memory_deallo(mem_servi(1:2,servi),'pelel','dod_memall',current_subdomain % pelel)
        call memory_deallo(mem_servi(1:2,servi),'lelel','dod_memall',current_subdomain % lelel)
        call memory_deallo(mem_servi(1:2,servi),'pelpo','dod_memall',current_subdomain % pelpo)
        call memory_deallo(mem_servi(1:2,servi),'lelpo','dod_memall',current_subdomain % lelpo)
        call memory_deallo(mem_servi(1:2,servi),'pbobo','dod_memall',current_subdomain % pbobo)
        call memory_deallo(mem_servi(1:2,servi),'lbobo','dod_memall',current_subdomain % lbobo)
        call memory_deallo(mem_servi(1:2,servi),'ppopo','dod_memall',current_subdomain % ppopo)
        call memory_deallo(mem_servi(1:2,servi),'lpopo','dod_memall',current_subdomain % lpopo)
        call memory_deallo(mem_servi(1:2,servi),'pbopo','dod_memall',current_subdomain % pbopo)
        call memory_deallo(mem_servi(1:2,servi),'lbopo','dod_memall',current_subdomain % lbopo)
        !
        ! Fringe node type LPEXT
        !
        if( associated(lpext) ) then
           do ipoin = 1,number_fringe_nodes
              call memory_deallo(mem_servi(1:2,servi),'LPEXT(IPOIN) % LNODS',     'dod_memall',lpext(ipoin) % lnods)
              call memory_deallo(mem_servi(1:2,servi),'LPEXT(IPOIN) % LTYPE',     'dod_memall',lpext(ipoin) % ltype)
              call memory_deallo(mem_servi(1:2,servi),'LPEXT(IPOIN) % CANDIDATES','dod_memall',lpext(ipoin) % candidates)
              call memory_deallo(mem_servi(1:2,servi),'LPEXT(IPOIN) % LMATE',     'dod_memall',lpext(ipoin) % lmate)
              call memory_deallo(mem_servi(1:2,servi),'LPEXT(IPOIN) % LELEZ'     ,'dod_memall',lpext(ipoin) % lelez)
              call memory_deallo(mem_servi(1:2,servi),'LPEXT(IPOIN) % LESUB'     ,'dod_memall',lpext(ipoin) % lesub)
           end do
        end if
        !
        ! Global boundary
        !
        call memory_deallo(mem_servi(1:2,servi),'LNODB_GLOBAL'      ,'dod_memall',lnodb_global)
        call memory_deallo(mem_servi(1:2,servi),'LBOCH_GLOBAL'      ,'dod_memall',lboch_global)
        call memory_deallo(mem_servi(1:2,servi),'LTYPB_GLOBAL'      ,'dod_memall',ltypb_global)
        call memory_deallo(mem_servi(1:2,servi),'PBOPO_GLOBAL'      ,'dod_memall',pbopo_global)
        call memory_deallo(mem_servi(1:2,servi),'LBOPO_GLOBAL'      ,'dod_memall',lbopo_global)
        call memory_deallo(mem_servi(1:2,servi),'LEDGE_NPOIN_GLOBAL','dod_memall',ledge_npoin_global)
        call memory_deallo(mem_servi(1:2,servi),'LBOEL_GLOBAL'      ,'dod_memall',lboel_global)
        call memory_deallo(mem_servi(1:2,servi),'LELBO_GLOBAL'      ,'dod_memall',lelbo_global)

     end do
     !
     ! Fringe node structure
     !
     if( associated(lpext) ) deallocate( lpext )
     !
     ! SKD-Tree structure
     !
     call dod_kdtree(3_ip)
     !
     ! Subdomain structure
     !
     deallocate( subdomain )
      !
     ! Subdomain arrays
     !
     call memory_deallo(mem_servi(1:2,servi),'IHOLE_DOD','dod_memall',ihole_dod)
     call memory_deallo(mem_servi(1:2,servi),'IPATC_DOD','dod_memall',ipatc_dod)
     call memory_deallo(mem_servi(1:2,servi),'IPRES_DOD','dod_memall',ipres_dod)
     !
     ! Mesh arrays
     !
     call memory_deallo(mem_servi(1:2,servi),'LSUBD_NPOIN','dod_memall',lsubd_npoin)
     call memory_deallo(mem_servi(1:2,servi),'LSUBD_NELEM','dod_memall',lsubd_nelem)
     call memory_deallo(mem_servi(1:2,servi),'LSUBD_NBOUN','dod_memall',lsubd_nboun)
     call memory_deallo(mem_servi(1:2,servi),'LMATN_DOD'  ,'dod_memall',lmatn_dod)
     call memory_deallo(mem_servi(1:2,servi),'LPOIZ_DOD'  ,'dod_memall',lpoiz_dod)
     call memory_deallo(mem_servi(1:2,servi),'PRESCRIBED_BOUNDARIES','dod_memall',prescribed_boundaries)

  case ( 3_ip )
     !
     ! Extension element of fringe nodes
     !
     do isubd = 1,nsubd
        current_subdomain => subdomain(isubd) 
        npoin_subd        =  current_subdomain % npoin
        allocate( current_subdomain % extension(npoin_subd) )
        do ipoin = 1,npoin_subd
           current_subdomain % extension(ipoin) % number_candidates = 0
           nullify( current_subdomain % extension(ipoin) % lnods      )
           nullify( current_subdomain % extension(ipoin) % ltype      )
           nullify( current_subdomain % extension(ipoin) % candidates )
        end do
     end do

  case ( 4_ip ) 
     !
     ! Deallocate boundary graphs and SKD-tree structure after hole cutting
     !
     do isubd = 1,nsubd
        if( ihole_dod(isubd) /= 0 ) then
           current_subdomain => subdomain(isubd)
           call memory_deallo(mem_servi(1:2,servi),'PBOPO'    ,'dod_memall',current_subdomain % pbopo)
           call memory_deallo(mem_servi(1:2,servi),'LBOPO'    ,'dod_memall',current_subdomain % lbopo)
           call memory_deallo(mem_servi(1:2,servi),'COORD_BOU','dod_memall',current_subdomain % coord_bou)
           call memory_deallo(mem_servi(1:2,servi),'LNODB_BOU','dod_memall',current_subdomain % lnodb_bou)
           if( ipatc_dod(isubd) /= 0 ) then
              call memory_deallo(mem_servi(1:2,servi),'PBOBO'    ,'dod_memall',current_subdomain % pbobo)
              call memory_deallo(mem_servi(1:2,servi),'LBOBO'    ,'dod_memall',current_subdomain % lbobo)
           end if
           call kdtree(&
                2_ip,mnodb,&
                current_subdomain % npoin_bou , current_subdomain % nboun , current_subdomain % coord_bou, &
                current_subdomain % lnodb_bou , current_subdomain % ltypb , current_subdomain % fabox,     &
                current_subdomain % bobox     , current_subdomain % sabox , current_subdomain % blink,     &
                current_subdomain % stru2     , current_subdomain % ldist , current_subdomain % lnele      ) 
        end if
     end do

  case (-5_ip)
     !
     ! Deallocate inverse permutation arrays
     !
     call memory_deallo(mem_servi(1:2,servi),'LINVP_NPOIN','dod_memall',linvp_npoin)
     call memory_deallo(mem_servi(1:2,servi),'LINVP_NELEM','dod_memall',linvp_nelem)
     call memory_deallo(mem_servi(1:2,servi),'LINVP_NBOUN','dod_memall',linvp_nboun)

  case( 6_ip)
     !
     ! LPEXT
     !
     allocate( lpext(number_fringe_nodes) )
     do ipoin = 1,number_fringe_nodes
        nullify(lpext(ipoin) % lnods)
        nullify(lpext(ipoin) % ltype)
        nullify(lpext(ipoin) % candidates)
        nullify(lpext(ipoin) % lmate)
        nullify(lpext(ipoin) % lelez)
        nullify(lpext(ipoin) % lesub)
     end do

  end select

end subroutine dod_memall
