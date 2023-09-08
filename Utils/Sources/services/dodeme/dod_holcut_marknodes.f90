!-----------------------------------------------------------------------
!> @addtogroup Dodeme
!> @{
!> @file    dod_holcut_marknodes.f90
!> @author  Guillaume Houzeaux
!> @date    18/09/2012
!> @brief   Treat patch subdomains
!> @details Do the following
!!          - Find host nodes
!!          - Holify free nodes when all neighbors are hole nodes
!!          - Define hole elements
!> @} 
!-----------------------------------------------------------------------
subroutine dod_holcut_marknodes()

  use def_parame
  use def_master
  use def_elmtyp
  use def_domain
  use def_dodeme
  use def_kermod
  use mod_kdtree
  use mod_graphs
  use mod_messages, only : livinf
  implicit none
  integer(ip)              :: isubd,jsubd,ipoin,idime
  integer(ip)              :: insid,intij,dummi
  integer(ip)              :: jelem,kfl_bounda
  integer(ip)              :: ipoin_global
  integer(ip)              :: dumm1(2)
  real(rp)                 :: coori(3),xieta(3),chkdi
  real(rp)                 :: dista,dummy(ndime)
  real(rp)                 :: shape_tmp(64)       
  real(rp)                 :: deriv_tmp(192)      
  integer(ip),   parameter :: LCLOS=1

  call livinf(0_ip,'MARK HOLE NODES',0_ip)

  !----------------------------------------------------------------------
  !
  ! Up to now, all nodes are free:
  ! CURRENT_SUBDOMAIN % LSUBD_NPOIN(IPOIN) = 0
  !
  ! Find hole-nodes:
  ! CURRENT_SUBDOMAIN % LSUBD_NPOIN(IPOIN) = -JSUBD
  !
  !----------------------------------------------------------------------

  kfl_bounda = 0_ip
  !
  ! ISUBD is a background and JSUBD a patch
  !
  do isubd = 1,nsubd
     do jsubd = 1,nsubd
        intij = intyp_dod(isubd,jsubd)

        if( intij > 0 ) then

           if( ictop_dod(intij) == DOD_CHIMERA .or. ictop_dod(intij) == DOD_HOLED_PATCH ) then
              !
              ! Loop over all ipoin poitns of isubd to find their host elements
              !     
              current_subdomain  => subdomain(isubd)
              neighbor_subdomain => subdomain(jsubd)
              nodes: do ipoin = 1,current_subdomain % npoin
                 !
                 ! Check if ipoin is inside the automatic hosting box of jsubd
                 !
                 insid        = 1
                 ipoin_global = current_subdomain % lnper(ipoin)
                 coori(    1) = current_subdomain % coord(    1,ipoin)
                 coori(    2) = current_subdomain % coord(    2,ipoin)
                 coori(ndime) = current_subdomain % coord(ndime,ipoin)
                 !
                 ! Check if ipoin is inside embedding box: INSID = 1
                 !
                 dimensions: do idime = 1,ndime
                    if(  coori(idime) < neighbor_subdomain % embox(idime,1) .or.&
                         coori(idime) > neighbor_subdomain % embox(idime,2) ) then
                       insid = 0
                       exit dimensions      
                    end if
                 end do dimensions
                 !
                 ! Find if there exists a host element in the background mesh
                 !
                 inside: if( insid == 1 ) then

                    if( kfl_bounda == LCLOS ) then 
                       !
                       ! Use KD-Tree
                       !
                       chkdi = 1.0e9_rp
                       call dpopar(&
                            1_ip,coori,&
                            current_subdomain % npoin_bou,mnodb,&
                            current_subdomain % nboun,chkdi,current_subdomain % ltypb,current_subdomain % lnodb_bou,&
                            current_subdomain % coord_bou,dista,dummy,dummy,dummi,&
                            current_subdomain % fabox,current_subdomain % sabox,current_subdomain % blink,&
                            current_subdomain % stru2,current_subdomain % ldist,current_subdomain % lnele)

                       if( ( dista <= 0.0_rp .and. ndime == 2 ) .or. ( dista >= 0.0_rp .and. ndime == 3 ) ) then
                          !
                          ! IPOIN is inside patch interface 
                          !
                          current_subdomain % lsubd_npoin(ipoin) = -jsubd  ! This is a hole node
                       end if

                    else
                       !
                       ! Use Elsest
                       !
                       call elsest(&
                            2_ip,jsubd,ielse,mnode,ndime,neighbor_subdomain % npoin,neighbor_subdomain % nelem,&
                            nnode(1:),neighbor_subdomain % lnods,neighbor_subdomain % ltype,&
                            ltopo,neighbor_subdomain % coord,coori,relse,jelem,&
                            shape_tmp,deriv_tmp,xieta,dumm1)
                       if( jelem > 0_ip ) then
                          !
                          ! IPOIN has host element JELEM
                          !
                          current_subdomain % lsubd_npoin(ipoin)  = -jsubd ! This is a hole node
                          current_subdomain % host_element(ipoin) =  jelem ! JELEM is in local numbering of JSUBD
                       else
                          current_subdomain % host_element(ipoin) =  0_ip
                       end if

                    end if

                 end if inside

              end do nodes

           end if

        end if

     end do

  end do

end subroutine dod_holcut_marknodes
