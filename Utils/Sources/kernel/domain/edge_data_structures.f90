!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    edge_data_structures.f90
!> @date    04/02/2016
!> @author  Guillaume Houzeaux
!> @brief   Edge data structures
!> @details Compute edge data structures in MESHE(NDIVI):
!>          \verbatim
!>          MEDGE ....................... Maximum number of edges per elements in mesh
!>          NEDGE ....................... Tolal number of edges in the mesh
!>          LNNED(1:NELEM) .............. Element number of edges
!>          LNNEB(1:NBOUN) .............. Boundary number of edges
!>          LEDGS(1:PEDGE,1:NELEM) ...... List of global edges for all elements
!>          LEDGB(1:PEDGE,1:NBOUN) ...... List of global edges for all boundaries
!>          R_EDG(:), C_EDG(:) .......... Edge graph
!>          EDGE_TO_NODE(1:2,1:NEDGE) ... Edge nodes
!>          \endverbatim
!>
!> @} 
!-----------------------------------------------------------------------

subroutine edge_data_structures()
  use def_kintyp, only : ip
  use def_kermod, only : kfl_edge_elements
  use def_kermod, only : ndivi
  use def_domain, only : meshe
  use mod_graphs, only : graphs_edges
  use def_master, only : INOTMASTER
  use def_master, only : IPARALL
  use mod_parall, only : PAR_COMM_MY_CODE_ARRAY
  use mod_messages, only : livinf
  implicit none

  if( kfl_edge_elements == 1 ) then
     call livinf(0_ip,'COMPUTE EDGE DATA STRUCTURES',0_ip)
     !
     ! Edge arrays
     !
     if( INOTMASTER ) call graphs_edges(meshe(ndivi)) 
     !call check_edge_data_structures()
     !
     ! Communication arrays: edges are renumbered to account
     ! for interior, own and other boundary edges
     !
     if( IPARALL ) then
        call par_edge_communication_array(meshe(ndivi),PAR_COMM_MY_CODE_ARRAY)
     end if
  end if

end subroutine edge_data_structures

subroutine check_edge_data_structures()
  use def_domain
  use def_master
  use def_kermod
  implicit none
  integer(ip) :: ielem,iedge,izdom

  print*,'MEDGE=',meshe(ndivi) % medge
  print*,'NEDGE=',meshe(ndivi) % nedge
  print*,'CONNECTIVITY:'
  do ielem = 1,nelem
     print*,ielem,meshe(ndivi) % ledgs(:,ielem)
  end do
  print*,'EDGES:'
  do iedge = 1,meshe(ndivi) % nedge
     print*,iedge,meshe(ndivi) % edge_to_node(1:2,iedge)
  end do 
  print*,'NUMBER:'
  do ielem = 1,meshe(ndivi) % nelem
     print*,iedge,meshe(ndivi) % lnned(ielem)
  end do 
  print*,'GRAPHS:'
  do iedge = 1,meshe(ndivi) % nedge
     do izdom = meshe(ndivi) % r_edg(iedge),meshe(ndivi) % r_edg(iedge+1)-1
        print*,iedge,meshe(ndivi) % c_edg(izdom)
     end do
  end do 

  call runend('O.K.!')

end subroutine check_edge_data_structures
