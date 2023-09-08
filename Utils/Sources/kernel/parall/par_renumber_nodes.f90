 !-----------------------------------------------------------------------
!> @addtogroup Parall
!> @{
!> @file    par_renumber_nodes.f90
!> @date    21/03/2017
!> @author  Guillaume Houzeaux
!> @brief   Renumber nodes
!> @details Renumber nodes in this order, with npoi2 = npoi1+1.
!>          Interior nodes are renumbered using METIS.
!>
!>          1       -> npoi1: interior
!>          npoi2   -> npoi3: own boundary
!>          npoi3+1 -> npoin: oth boundary
!>
!>          Then, halos will also be renumbered later on. New and old
!>          notations are:
!>
!>          NEW RENUMBERING:
!>
!>          +---------------+------------+------------+----------+-------+
!>          |   interior    | own bound. | oth bound. |      geom halo   |
!>          +---------------+------------+------------+----------+-------+
!>
!>          ---------------->------------>------------>------------------>
!>          1          npoi1 npoi2  npoi3         npoin            npoin_2
!>
!>
!>          AFTER RENUMBERING:
!>
!>          +----------------------------+-----------------+-------------+
!>          |    Own nodes               |    comp halo    | rest halo   |
!>          +----------------------------+-----------------+-------------+
!>
!>          ----------------------------->----------------->------------->
!>          1                    npoin_own       npoin_halo        npoin_2
!>
!> @}
!-----------------------------------------------------------------------

subroutine par_renumber_nodes()

  use def_kintyp,      only : ip,rp
  use def_master,      only : ISLAVE
  use def_master,      only : INOTSLAVE
  use def_master,      only : ISEQUEN
  use def_master,      only : IMASTER
  use def_master,      only : INOTMASTER
  use def_master,      only : npoi1
  use def_master,      only : npoi2
  use def_master,      only : npoi3
  use def_kermod,      only : kfl_renumbering_npoin
  use def_kermod,      only : nsfc_renumbering_npoin
  use def_domain,      only : npoin_own
  use def_domain,      only : npoin_halo
  use def_domain,      only : npoin
  use def_domain,      only : mnode
  use def_domain,      only : nelem
  use def_domain,      only : lnods
  use def_domain,      only : ltype
  use def_domain,      only : lnnod
  use def_domain,      only : coord
  use mod_parall,      only : par_memor
  use mod_parall,      only : PAR_COMM_MY_CODE_ARRAY
  use mod_memory,      only : memory_alloca
  use mod_memory,      only : memory_deallo
  use mod_renumbering, only : renumbering_node_arrays
  use mod_renumbering, only : renumbering_nodes
  use mod_renumbering, only : renumbering_sfc_recursive
  use mod_renumbering, only : renumbering_reverse_cuthill_mckee
  use mod_graphs,      only : graphs_poipoi
  use mod_mesh_type,   only : mesh_type_update_last_mesh
  use mod_elmgeo,      only : element_type
  use mod_messages,    only : messages_live
  use mod_alya2metis,  only : alya2metis_METIS_NodeND

  use def_master
  implicit none

  integer(ip)          :: ipoin,kpoin,ielem,pelty
  integer(ip), pointer :: permR(:)
  integer(ip), pointer :: lnnod_loc(:)
  character(15)        :: my_method
  !
  ! Interior graph
  !
  integer(ip), pointer :: ia(:)
  integer(ip), pointer :: ja(:)

  nullify(permR)
  nullify(lnnod_loc)
  nullify(ia)
  nullify(ja)
  
  if( INOTSLAVE ) then
     if(      kfl_renumbering_npoin == 0 ) then
        my_method = 'NOTHING!...'
     else if( kfl_renumbering_npoin == 1 ) then
        my_method = 'METIS'
     else if( kfl_renumbering_npoin == 2 ) then
        my_method = 'SFC'
     else if( kfl_renumbering_npoin == 3 ) then
        my_method = 'CUTHILL-MCKEE'
     endif
     call messages_live('NODE RENUMBERING WITH '//trim(my_method))
  end if
   
  if( npoin > 0 ) then

     call memory_alloca(par_memor,'permR','par_renumber_nodes',permR,npoin)
     !
     ! PERMR: Identify boundary nodes. Others'=-1, Own=-2
     !
     do ipoin = npoi1+1,npoi2-1
        permR(ipoin) = -1
     end do
     do ipoin = npoi2,npoi3
        permR(ipoin) = -2
     end do
     do ipoin = npoi3+1,npoin
        permR(ipoin) = -1
     end do
     !
     ! Renumber interior nodes
     !
     kpoin = npoi1
     do ipoin = 1,npoi1
        permR(ipoin) = ipoin
     end do
     !
     ! Renumber own boundary nodes
     !
     npoi1 = kpoin
     npoi2 = kpoin + 1
     do ipoin = 1,npoin
        if( permR(ipoin) == -2 ) then
           kpoin        = kpoin + 1
           permR(ipoin) = kpoin
        end if
     end do
     npoi3 = kpoin
     !
     ! Renumber others' boundary nodes
     !
     do ipoin = 1,npoin
        if( permR(ipoin) == -1 ) then
           kpoin        = kpoin + 1
           permR(ipoin) = kpoin
        end if
     end do
     !
     ! IA, JA: graph of interior nodes without diagonal
     !
     if( kfl_renumbering_npoin == 1 .or. kfl_renumbering_npoin == 3 ) then
        if( .not. associated(lnnod) ) then
           call memory_alloca(par_memor,'LNNOD_LOC','par_renumber_nodes',lnnod_loc,nelem)
           do ielem = 1,nelem
              pelty = abs(ltype(ielem))
              lnnod_loc(ielem) = element_type(pelty) % number_nodes
           end do
        else
           lnnod_loc => lnnod
        end if
        call graphs_poipoi(&
             npoi1,nelem,mnode,lnods,lnnod_loc,ltype,ia,ja,'SQUARE REMOVE DIAGONAL')        
     end if
     !
     ! Reorder interior nodes 
     !
     if(      kfl_renumbering_npoin == 1 ) then
        !
        ! METIS renumbering
        !
        !call renumbering_nodes('METIS',npoi1,ia,ja,permR)
        call alya2metis_METIS_NodeND(npoi1,ia,ja,permr)

     else if( kfl_renumbering_npoin == 2 ) then
        !
        ! SFC renumbering
        !
        call renumbering_sfc_recursive(nsfc_renumbering_npoin,coord,permr,npoi1)

     else if( kfl_renumbering_npoin == 3 ) then
        !
        ! Cuthill-McKee renumbering
        !
        call renumbering_reverse_cuthill_mckee(npoi1,ia,ja,permr)
        
     else
        !
        ! No numbering
        !
     end if
     !
     ! Deallocate if necessary
     !
     if( .not. associated(lnnod) ) &
          call memory_deallo(par_memor,'LNNOD_LOC','par_renumber_nodes',lnnod_loc)
     !
     ! Permute nodal arrays
     !
     call renumbering_node_arrays(permR)
     !
     ! Deallocate memory
     !
     call memory_deallo(par_memor,'PERMR','par_renumber_nodes',permR)
     call memory_deallo(par_memor,'IA'   ,'par_renumber_nodes',ia)
     call memory_deallo(par_memor,'JA'   ,'par_renumber_nodes',ja)

  end if

  if( ISLAVE ) then
     !
     ! Useful parameters
     !
     npoin_own  = npoi3
     npoin_halo = npoin_own
     !
     ! Update main communicator
     !
     PAR_COMM_MY_CODE_ARRAY(1) % npoi1 = npoi1
     PAR_COMM_MY_CODE_ARRAY(1) % npoi2 = npoi2
     PAR_COMM_MY_CODE_ARRAY(1) % npoi3 = npoi3
     
  else if( ISEQUEN ) then

     npoi3      = npoin
     npoin_own  = npoin
     npoin_halo = npoin

  else if( IMASTER ) then

     npoi3      = 0
     npoin_own  = 0
     npoin_halo = 0

  end if
  !
  ! Point to mesh structure => MESHE(NDIVI)
  !
  call mesh_type_update_last_mesh()

end subroutine par_renumber_nodes
