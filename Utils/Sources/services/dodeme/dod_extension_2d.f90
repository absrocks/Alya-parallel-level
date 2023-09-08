!-----------------------------------------------------------------------
!> @addtogroup Dodeme
!> @{
!> @file    dod_extension_2d.f90
!> @author  Guillaume Houzeaux
!> @date    18/09/2012
!> @brief   Create 2D extensions
!> @details Create extension elements for 2D meshes
!> @} 
!-----------------------------------------------------------------------
subroutine dod_extension_2d()
  use def_kintyp
  use def_parame
  use def_master
  use def_elmtyp
  use def_domain
  use def_dodeme
  use mod_memory
  use mod_kdtree
  use mod_dod_extens
  use mod_messages, only : livinf
  implicit none
  integer(ip)          :: isubd,jsubd,ielem
  integer(ip)          :: kelem,knode,boun1
  integer(ip)          :: perm1,perm2,perm3,boun2
  integer(ip)          :: ipoin,jnode,pnodb
  integer(ip)          :: ibopo,ntext,imatn,ipoiz
  integer(ip)          :: pext1,pext2,pext3,pext4
  integer(ip)          :: inodb,iboun
  integer(ip)          :: inext,mnext,nnext,npoin_subd
  real(rp)             :: time1,q,q1,q2,q3,dummy(3)
  real(rp)             :: coord_ipoin(2)
  real(rp)             :: coord_pext1(2)
  real(rp)             :: coord_pext2(2)
  real(rp)             :: coord_pext3(2)
  real(rp)             :: coord_pext4(2)
  integer(ip), pointer :: lpoex_bo1(:)
  integer(ip), pointer :: lpoex_bo2(:)
  integer(ip), pointer :: lnext(:)
  integer(ip)          :: nboun_subd
  integer(ip)          :: ipoin_subd

  call livinf(0_ip,'CREATE 2D EXTENSION ELEMENTS',0_ip)

  dummy     = 0.0_rp
  nelem_dod = 0 

  !------------------------------------------------------
  !
  ! MNEXT: maximum number of extensions elements
  !
  !------------------------------------------------------

  mnext = 0
  number_fringe_nodes = 0
  do isubd = 1,nsubd
     current_subdomain => subdomain(isubd)
     do ipoin = 1, current_subdomain % npoin
        if( current_subdomain % lsubd_npoin(ipoin) > 0 ) then
           number_fringe_nodes = number_fringe_nodes + 1
           mnext = max(mnext,current_subdomain % extension(ipoin) % number_candidates)
        end if
     end do
  end do
  nullify(lnext)
  nullify(lpoex_bo1)
  nullify(lpoex_bo2)

  call memory_alloca(mem_servi(1:2,servi),'LNEXT',    'dod_extension_elements_2d',lnext,mnext)
  call memory_alloca(mem_servi(1:2,servi),'LPOEX_BO1','dod_extension_elements_2d',lpoex_bo1,mnext)
  call memory_alloca(mem_servi(1:2,servi),'LPOEX_BO2','dod_extension_elements_2d',lpoex_bo2,mnext)

  !------------------------------------------------------
  !
  ! Allocate memory for fringe nodes type
  !
  !------------------------------------------------------

  call dod_memall(6_ip)
  
  !------------------------------------------------------
  !
  ! Para cada punto de contorno:
  ! KFL_INTOP_DOD = 0 in general case
  !               = 1 if an extension node is on the real boundary
  ! +-----oo---+---
  !    i  || j
  !       ||
  !       ++
  !
  !------------------------------------------------------

  call cputim(time1)

  ielem = 0
  number_fringe_nodes = 0

  do isubd = 1,nsubd

     current_subdomain => subdomain(isubd)
     npoin_subd        =  current_subdomain % npoin

     punto_contorno: do ipoin_subd = 1,npoin_subd

        if( current_subdomain % lsubd_npoin(ipoin_subd) > 0 ) then

           number_fringe_nodes =  number_fringe_nodes + 1
           ipoin               =  current_subdomain % lnper(ipoin_subd)
           jsubd               =  current_subdomain % lsubd_npoin(ipoin_subd)
           coord_ipoin(1)      =  current_subdomain % coord(1,ipoin_subd)
           coord_ipoin(2)      =  current_subdomain % coord(2,ipoin_subd)
           neighbor_subdomain  => subdomain(jsubd)

           nnext = current_subdomain % extension(ipoin_subd) % number_candidates 
           do jnode = 1,nnext
              lnext(jnode) = neighbor_subdomain % lnper(current_subdomain % extension(ipoin_subd) % candidates(jnode))
           end do
           do inext = 1,nnext
              lpoex_bo1(inext) = 0
              lpoex_bo2(inext) = 0
           end do
           !             pext1 o    
           !
           !      boun2  ipoin
           !   o---------x          o pext2
           ! pext4       |
           !     isubd   | boun1
           !             |
           !             o pext3
           !
           do ibopo = pbopo_global(ipoin),pbopo_global(ipoin+1)-1
              iboun = lbopo_global(ibopo)
              if( lnodb_global(2,iboun) == ipoin ) then
                 boun1 = iboun
              else
                 boun2 = iboun
              end if
           end do
           if( lboch_global(boun1) == BOFEM .or. lboch_global(boun2) == BOFEM ) then
              !
              ! Node is also on real boundary: only one extension to create
              !             
              if( lboch_global(boun1) == BOEXT ) then
                 if( lnodb_global(1,boun1) == ipoin ) then
                    pext3 = lnodb_global(2,boun1)
                 else
                    pext3 = lnodb_global(1,boun1)
                 end if
              else if( lboch_global(boun2) == BOEXT ) then
                 if( lnodb_global(1,boun2) == ipoin ) then
                    pext3 = lnodb_global(2,boun2)
                 else
                    pext3 = lnodb_global(1,boun2)
                 end if
              else
                 call runend('DOD_EXTENSION_2D: IMPOSSIBLE CASE')
              end if
print*,'holaaaaaaaaaa',ipoin,pext3,boun1,boun2,lboch_global(boun1),lboch_global(boun2)
              call dod_extension_elements_2d(&       
                ipoin,boun1,boun2,mnodb,lnodb_global,lboch_global,pbopo_global,lbopo_global,coord,&
                nnext,lnext,lpoex_bo1,lpoex_bo2,ntext,pext1,pext2)

              coord_pext1(1) = coord(1,pext1)
              coord_pext1(2) = coord(2,pext1)
              coord_pext3(1) = coord(1,pext3)
              coord_pext3(2) = coord(2,pext3)
              imatn          = lmatn_dod(pext1)
              ipoiz          = lpoiz_dod(pext1)
              
              call dod_extens_qual2d(&
                   coord_ipoin,coord_pext1,coord_pext3,pext1,pext3,q,perm1)
              print*,'holaaaaaaa-cacaaaaa',ipoin,pext1,pext3,q
              kelem = 1
              knode = 3
              call memory_alloca(mem_servi(1:2,servi),'LPEXT % LNODS','dod_extension_elements_2d',lpext(number_fringe_nodes) % lnods,knode,kelem,'DO_NOT_INITIALIZE')
              call memory_alloca(mem_servi(1:2,servi),'LPEXT % LTYPE','dod_extension_elements_2d',lpext(number_fringe_nodes) % ltype,kelem      ,'DO_NOT_INITIALIZE')
              call memory_alloca(mem_servi(1:2,servi),'LPEXT % LMATE','dod_extension_elements_2d',lpext(number_fringe_nodes) % lmate,kelem      ,'DO_NOT_INITIALIZE')
              call memory_alloca(mem_servi(1:2,servi),'LPEXT % LELEZ','dod_extension_elements_2d',lpext(number_fringe_nodes) % lelez,kelem      ,'DO_NOT_INITIALIZE')
              call memory_alloca(mem_servi(1:2,servi),'LPEXT % LESUB','dod_extension_elements_2d',lpext(number_fringe_nodes) % lesub,kelem      ,'DO_NOT_INITIALIZE')
              lpext(number_fringe_nodes) % ltype(1)   = TRI03
              lpext(number_fringe_nodes) % lnods(1,1) = ipoin
              lpext(number_fringe_nodes) % lnods(2,1) = pext1
              lpext(number_fringe_nodes) % lnods(3,1) = pext3
              lpext(number_fringe_nodes) % lmate(1)   = imatn
              lpext(number_fringe_nodes) % lelez(1)   = ipoiz
              lpext(number_fringe_nodes) % lesub(1)   = isubd
            !      if(ipoin==271) print*,lpext(number_fringe_nodes) % lnods(:,:)

           else 
              !
              ! Normal node (not belonging to real boundary)
              !              
              pext1 = 0
              pext2 = 0
   
              call dod_extension_elements_2d(&       
                   ipoin,boun1,boun2,mnodb,lnodb_global,lboch_global,pbopo_global,lbopo_global,coord,&
                   nnext,lnext,lpoex_bo1,lpoex_bo2,ntext,pext1,pext2)

              if( ntext == 3 ) then
                 !
                 ! 3 extension elements
                 !
                 pext3          = lnodb_global(1,boun1)
                 pext4          = lnodb_global(2,boun2)
                 coord_pext1(1) = coord(1,pext1)
                 coord_pext1(2) = coord(2,pext1)
                 coord_pext2(1) = coord(1,pext2)
                 coord_pext2(2) = coord(2,pext2)
                 coord_pext3(1) = coord(1,pext3)
                 coord_pext3(2) = coord(2,pext3)
                 coord_pext4(1) = coord(1,pext4)
                 coord_pext4(2) = coord(2,pext4)
                 imatn          = lmatn_dod(pext1)
                 ipoiz          = lpoiz_dod(pext1)

                 kelem = 3
                 knode = 3
                 call memory_alloca(mem_servi(1:2,servi),'LPEXT % LNODS','dod_extension_elements_2d',lpext(number_fringe_nodes) % lnods,knode,kelem,'DO_NOT_INITIALIZE')
                 call memory_alloca(mem_servi(1:2,servi),'LPEXT % LTYPE','dod_extension_elements_2d',lpext(number_fringe_nodes) % ltype,kelem      ,'DO_NOT_INITIALIZE')
                 call memory_alloca(mem_servi(1:2,servi),'LPEXT % LMATE','dod_extension_elements_2d',lpext(number_fringe_nodes) % lmate,kelem      ,'DO_NOT_INITIALIZE')
                 call memory_alloca(mem_servi(1:2,servi),'LPEXT % LELEZ','dod_extension_elements_2d',lpext(number_fringe_nodes) % lelez,kelem      ,'DO_NOT_INITIALIZE')
                 call memory_alloca(mem_servi(1:2,servi),'LPEXT % LESUB','dod_extension_elements_2d',lpext(number_fringe_nodes) % lesub,kelem      ,'DO_NOT_INITIALIZE')
                 lpext(number_fringe_nodes) % ltype(1) = TRI03
                 lpext(number_fringe_nodes) % ltype(2) = TRI03
                 lpext(number_fringe_nodes) % ltype(3) = TRI03
                 lpext(number_fringe_nodes) % lmate(1) = imatn
                 lpext(number_fringe_nodes) % lmate(2) = imatn
                 lpext(number_fringe_nodes) % lmate(3) = imatn
                 lpext(number_fringe_nodes) % lelez(1) = ipoiz
                 lpext(number_fringe_nodes) % lelez(2) = ipoiz
                 lpext(number_fringe_nodes) % lelez(3) = ipoiz
                 lpext(number_fringe_nodes) % lesub(1) = isubd
                 lpext(number_fringe_nodes) % lesub(2) = isubd
                 lpext(number_fringe_nodes) % lesub(3) = isubd
                 call dod_extens_qual2d(&
                      coord_ipoin,coord_pext3,coord_pext1,pext3,pext1,q1,perm1)
                 lpext(number_fringe_nodes) % lnods(1,1) = ipoin
                 lpext(number_fringe_nodes) % lnods(2,1) = pext3
                 lpext(number_fringe_nodes) % lnods(3,1) = pext1                 
                 if(ipoin==5333)print*,'entro1',ipoin,pext3,pext1 ,q1,perm1
                 
                 call dod_extens_qual2d(&
                      coord_ipoin,coord_pext2,coord_pext4,pext2,pext4,q2,perm2)
                 lpext(number_fringe_nodes) % lnods(1,2) = ipoin
                 lpext(number_fringe_nodes) % lnods(2,2) = pext2
                 lpext(number_fringe_nodes) % lnods(3,2) = pext4
                 if( perm1 == 1 .and. perm2 == 1 ) then
                    call dod_extens_qual2d(&
                         coord_ipoin,coord_pext3,coord_pext4,pext3,pext4,q3,perm3)
                    lpext(number_fringe_nodes) % lnods(1,3) = ipoin
                    lpext(number_fringe_nodes) % lnods(2,3) = pext3
                    lpext(number_fringe_nodes) % lnods(3,3) = pext4
                    
                 else if( perm1 == 0 .and. perm2 == 1 ) then
                    call dod_extens_qual2d(&
                         coord_ipoin,coord_pext1,coord_pext4,pext1,pext4,q3,perm3)
                    lpext(number_fringe_nodes) % lnods(1,3) = ipoin
                    lpext(number_fringe_nodes) % lnods(2,3) = pext1
                    lpext(number_fringe_nodes) % lnods(3,3) = pext4
                    
                 else if( perm1 == 1 .and. perm2 == 0 ) then
                    call dod_extens_qual2d(&
                         coord_ipoin,coord_pext3,coord_pext2,pext3,pext2,q3,perm3)
                    lpext(number_fringe_nodes) % lnods(1,3) = ipoin
                    lpext(number_fringe_nodes) % lnods(2,3) = pext3
                    lpext(number_fringe_nodes) % lnods(3,3) = pext2
                    
                 else if( perm1 == 0 .and. perm2 == 0 ) then
                    call dod_extens_qual2d(&
                         coord_ipoin,coord_pext1,coord_pext2,pext1,pext2,q3,perm3)
                    lpext(number_fringe_nodes) % lnods(1,3) = ipoin
                    lpext(number_fringe_nodes) % lnods(2,3) = pext1
                    lpext(number_fringe_nodes) % lnods(3,3) = pext2
                    
                 end if
                 
              else if( ntext == 2 ) then
                 !
                 ! 2 extension elements
                 !
                 pext2          = lnodb_global(2,boun2)
                 pext3          = lnodb_global(1,boun1)
                 coord_pext1(1) = coord(1,pext1)
                 coord_pext1(2) = coord(2,pext1)
                 coord_pext2(1) = coord(1,pext2)
                 coord_pext2(2) = coord(2,pext2)
                 coord_pext3(1) = coord(1,pext3)
                 coord_pext3(2) = coord(2,pext3)
                 imatn          = lmatn_dod(pext1)
                 ipoiz          = lpoiz_dod(pext1)
                 
                 kelem = 2
                 knode = 3
                 call memory_alloca(mem_servi(1:2,servi),'LPEXT % LNODS','dod_extension_elements_2d',lpext(number_fringe_nodes) % lnods,knode,kelem,'DO_NOT_INITIALIZE')
                 call memory_alloca(mem_servi(1:2,servi),'LPEXT % LTYPE','dod_extension_elements_2d',lpext(number_fringe_nodes) % ltype,kelem      ,'DO_NOT_INITIALIZE')
                 call memory_alloca(mem_servi(1:2,servi),'LPEXT % LMATE','dod_extension_elements_2d',lpext(number_fringe_nodes) % lmate,kelem      ,'DO_NOT_INITIALIZE')
                 call memory_alloca(mem_servi(1:2,servi),'LPEXT % LELEZ','dod_extension_elements_2d',lpext(number_fringe_nodes) % lelez,kelem      ,'DO_NOT_INITIALIZE')
                 call memory_alloca(mem_servi(1:2,servi),'LPEXT % LESUB','dod_extension_elements_2d',lpext(number_fringe_nodes) % lesub,kelem      ,'DO_NOT_INITIALIZE')
                 lpext(number_fringe_nodes) % ltype(1) = TRI03
                 lpext(number_fringe_nodes) % ltype(2) = TRI03
                 lpext(number_fringe_nodes) % lmate(1) = imatn
                 lpext(number_fringe_nodes) % lmate(2) = imatn
                 lpext(number_fringe_nodes) % lelez(1) = ipoiz
                 lpext(number_fringe_nodes) % lelez(2) = ipoiz
                 lpext(number_fringe_nodes) % lesub(1) = isubd
                 lpext(number_fringe_nodes) % lesub(2) = isubd
                 
                 call dod_extens_qual2d(&
                      coord_ipoin,coord_pext1,coord_pext2,pext1,pext2,q1,perm1)
                 lpext(number_fringe_nodes) % lnods(1,1) = ipoin
                 lpext(number_fringe_nodes) % lnods(2,1) = pext1
                 lpext(number_fringe_nodes) % lnods(3,1) = pext2
                 
                 call dod_extens_qual2d(&
                      coord_ipoin,coord_pext3,coord_pext1,pext3,pext1,q2,perm2)
                 lpext(number_fringe_nodes) % lnods(1,2) = ipoin
                 lpext(number_fringe_nodes) % lnods(2,2) = pext3
                 lpext(number_fringe_nodes) % lnods(3,2) = pext1
                 
              end if
              
          end if
           
           nelem_dod = nelem_dod + kelem
           
        end if
     
     end do punto_contorno

  end do

  call memory_deallo(mem_servi(1:2,servi),'LNEXT'       ,'dod_extension_elements_2d',lnext)
  call memory_deallo(mem_servi(1:2,servi),'LPOEX_BO1'   ,'dod_extension_elements_2d',lpoex_bo1)
  call memory_deallo(mem_servi(1:2,servi),'LPOEX_BO2'   ,'dod_extension_elements_2d',lpoex_bo2)

  call memory_deallo(mem_servi(1:2,servi),'LNODB_GLOBAL','dod_extension_elements_2d',lnodb_global)
  call memory_deallo(mem_servi(1:2,servi),'LBOCH_GLOBAL','dod_extension_elements_2d',lboch_global)
  call memory_deallo(mem_servi(1:2,servi),'LTYPB_GLOBAL','dod_extension_elements_2d',ltypb_global)
  call memory_deallo(mem_servi(1:2,servi),'PBOPO_GLOBAL','dod_extension_elements_2d',pbopo_global)
  call memory_deallo(mem_servi(1:2,servi),'LBOPO_GLOBAL','dod_extension_elements_2d',lbopo_global)

end subroutine dod_extension_2d
