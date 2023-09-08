!-----------------------------------------------------------------------
!> @addtogroup Dodeme
!> @{
!> @file    dod_holcut.f90
!> @author  Guillaume Houzeaux
!> @date    18/09/2012
!> @brief   Treat patch subdomains
!> @details Do the following
!>          - Find host nodes
!>          - Holify free nodes when all neighbors are hole nodes
!>          - Define hole elements
!> @} 
!-----------------------------------------------------------------------
subroutine dod_holcut()

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
  integer(ip) :: isubd,jsubd,ipoin,inode,ielem,knode,pnode
  integer(ip) :: jpoin,ipopo,intij
  real(rp)    :: time1,time2

  call livinf(-4_ip,'HOLE CUTTING',0_ip)

  !----------------------------------------------------------------------
  !
  ! Up to now, all nodes are free:
  ! CURRENT_SUBDOMAIN % LSUBD_NPOIN(IPOIN) = 0
  !
  ! Find hole-nodes:
  ! CURRENT_SUBDOMAIN % LSUBD_NPOIN(IPOIN) = -JSUBD
  !
  !----------------------------------------------------------------------

  call cputim(time1)
  call dod_holcut_marknodes()
  call cputim(time2) ; cpu_dod_holcut_marknodes = time2-time1

  !----------------------------------------------------------------------
  !
  ! Holed patch: mark another layer of hole nodes
  !
  !----------------------------------------------------------------------

  do isubd = 1,nsubd
     current_subdomain => subdomain(isubd)
     do jsubd = 1,nsubd
        intij = intyp_dod(isubd,jsubd)        
        if( intij > 0 ) then
           if( ictop_dod(intij) == DOD_HOLED_PATCH ) then              
              do ipoin = 1,current_subdomain % npoin
                 if( current_subdomain % lsubd_npoin(ipoin) == -jsubd ) then
                    do ipopo = current_subdomain % ppopo(ipoin),current_subdomain % ppopo(ipoin+1)-1
                       jpoin = current_subdomain % lpopo(ipopo)
                       if( current_subdomain % lsubd_npoin(jpoin) /= -jsubd ) & 
                            current_subdomain % lsubd_npoin(jpoin) = -nsubd-jsubd
                    end do
                 end if
              end do
              do ipoin = 1,current_subdomain % npoin
                 if( current_subdomain % lsubd_npoin(ipoin) < -nsubd ) then
                    current_subdomain % lsubd_npoin(ipoin) = current_subdomain % lsubd_npoin(ipoin) + nsubd
                 end if
              end do
           end if
        end if 
     end do
  end do

  !----------------------------------------------------------------------
  !
  ! Define Hole element: a hole element has all hole node 
  ! CURRENT_SUBDOMAIN % LSUBD_NELEM(IELEM) = -JSUBD
  !
  !----------------------------------------------------------------------

  do isubd = 1,nsubd
     if( ihole_dod(isubd) == 1 ) then
        current_subdomain => subdomain(isubd)
        do ielem = 1,current_subdomain % nelem
           pnode = current_subdomain % lnnod(ielem)
           knode = 0 
           inode = 0
           do while( inode < pnode )
              inode = inode + 1
              ipoin = current_subdomain % lnods(inode,ielem) 
              jsubd = current_subdomain % lsubd_npoin(ipoin)
              if( jsubd < 0 ) then
                 knode = knode + 1
              else
                 inode = pnode
              end if
           end do
           if( knode == pnode ) then
              current_subdomain % lsubd_nelem(ielem) = -abs(jsubd)
           end if
        end do
     end if
  end do

  !----------------------------------------------------------------------
  !
  ! Touch free elements recursively to gather hole and solid elements
  !
  !----------------------------------------------------------------------

  call cputim(time1)
  call dod_holcut_inversehole()
  call cputim(time2) ; cpu_dod_holcut_inversehole = time2-time1

 !----------------------------------------------------------------------
  !
  ! Free hole-nodes
  ! CURRENT_SUBDOMAIN % LSUBD_NPOIN(IPOIN) = -JUSBD => JSUBD
  !
  !----------------------------------------------------------------------

  do isubd = 1,nsubd
     if( ihole_dod(isubd) == 1 ) then
        current_subdomain => subdomain(isubd)
        do ipoin = 1,current_subdomain % npoin
           current_subdomain % lsubd_npoin(ipoin) = max(0_ip,current_subdomain % lsubd_npoin(ipoin))
        end do
     end if
  end do


  !----------------------------------------------------------------------
  !
  ! From candidate elements, define recursively the hole using common
  ! faces as a criterion
  !
  !----------------------------------------------------------------------

  call cputim(time1)
  call dod_holcut_markelements()
  call cputim(time2) ; cpu_dod_holcut_markelements = time2-time1

  !----------------------------------------------------------------------
  !
  ! Construct hole boundary
  !
  !----------------------------------------------------------------------

  call cputim(time1)
  call dod_holcut_holeboundary()
  call cputim(time2) ; cpu_dod_holcut_holeboundary = time2-time1

  !----------------------------------------------------------------------
  !
  ! Detect fringe nodes
  ! Input:  current_subdomain % lsubd_npoin(ipoin)=
  !                    Free    hole  fringe
  !                      0   -jsubd  -jsubd
  ! Output: current_subdomain % lsubd_npoin(ipoin)=
  !                    Free    hole  fringe
  !                      0   -jsubd   jsubd
  !
  !----------------------------------------------------------------------
  call cputim(time1)
  call dod_holcut_fringenodes()
  call cputim(time2) ; cpu_dod_holcut_fringenodes = time2-time1

  !----------------------------------------------------------------------
  !
  ! Recompute boundary arrays
  !
  !----------------------------------------------------------------------

  call cputim(time1)
  call dod_memall(4_ip) ! Deallocate graph and SKD-Tree
  call dod_graphs(2_ip) ! Recompute boundary graphs
  call dod_kdtree(2_ip) ! Recompute SKD-Tree
  call cputim(time2) ; cpu_dod_graphs_kdtree = time2-time1

  call livinf(-5_ip,'END HOLE CUTTING',0_ip)

end subroutine dod_holcut

