!-----------------------------------------------------------------------
!
!> @defgroup Mesh_Type_Toolbox
!> @{
!> @name    ToolBox for mesh type
!> @file    mod_mesh_type.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for mesh type
!> @details ToolBox for mesh type
!
!-----------------------------------------------------------------------

module mod_mesh_type

  use def_kintyp,    only : ip,rp
  use def_kermod,    only : ndivi
  use def_kermod,    only : kfl_posdi

  use def_domain,    only : npoin_total
  use def_domain,    only : nelem_total
  use def_domain,    only : nboun_total
  use def_domain,    only : ndime,ntens,npoin
  use def_domain,    only : nelem,nboun,mnode,mnodb
  use def_domain,    only : mgaus,nbopo,npoin_own
  use def_domain,    only : npoin_2,nelem_2
  use def_domain,    only : nboun_2,mnodb
  use def_domain,    only : nelem,nboun
  use def_domain,    only : npoin_halo,npoin

  use def_master,    only : npoi1,npoi2,npoi3
  use def_master,    only : npoin_par,nelem_par
  use def_master,    only : nboun_par,npart
  use def_master,    only : INOTSLAVE
  use def_master,    only : INOTMASTER
  use def_master,    only : IMASTER,ISEQUEN

  use def_domain,    only : memor_dom
  use def_domain,    only : meshe
  use def_master,    only : leinv_loc
  use def_master,    only : lbinv_loc
  use def_master,    only : lninv_loc
  use def_domain,    only : lnods
  use def_domain,    only : ltype
  use def_domain,    only : lnnod
  use def_domain,    only : lelch
  use def_domain,    only : lesub
  use def_domain,    only : lmate
  use def_domain,    only : lgaus
  use def_domain,    only : lnodb
  use def_domain,    only : ltypb
  use def_domain,    only : lboch
  use def_domain,    only : lboel
  use def_domain,    only : lelbo
  use def_domain,    only : lnnob
  use def_domain,    only : coord
  use def_domain,    only : lnoch
  use def_domain,    only : lmast
  use def_domain,    only : lpoty
  use def_domain,    only : nzdom
  use def_domain,    only : r_dom
  use def_domain,    only : c_dom
  use def_domain,    only : exnor
  use def_domain,    only : vmass
  use def_domain,    only : vmasc

  use mod_memory,    only : memory_copy
  use mod_memory,    only : memory_alloca

  implicit none

  public :: mesh_type_allocate_initialize
  public :: mesh_type_save_original_mesh
  public :: mesh_type_update_last_mesh
  public :: mesh_type_allocate_minimum
  
contains

  !-----------------------------------------------------------------------
  !
  !> @date    15/04/2017
  !> @author  Guillaume Houzeaux
  !> @brief   Allocate a minimum size
  !> @details Allocate a minimum size to mesh type
  !
  !-----------------------------------------------------------------------

  subroutine mesh_type_allocate_minimum(CURRENT_MESH)

    integer(ip), intent(in) :: CURRENT_MESH
    
    call memory_alloca(memor_dom,'LNINV_LOC','mesh_type_allocate_minimum',meshe(CURRENT_MESH) % lninv_loc,2_ip)
    call memory_alloca(memor_dom,'LEINV_LOC','mesh_type_allocate_minimum',meshe(CURRENT_MESH) % leinv_loc,2_ip)
    call memory_alloca(memor_dom,'LBINV_LOC','mesh_type_allocate_minimum',meshe(CURRENT_MESH) % lbinv_loc,2_ip)
    
    call memory_alloca(memor_dom,'LNODS'    ,'mesh_type_allocate_minimum',meshe(CURRENT_MESH) % lnods,ndime,2_ip)
    call memory_alloca(memor_dom,'LTYPE'    ,'mesh_type_allocate_minimum',meshe(CURRENT_MESH) % ltype,2_ip)
    call memory_alloca(memor_dom,'LELCH'    ,'mesh_type_allocate_minimum',meshe(CURRENT_MESH) % lelch,2_ip)
    call memory_alloca(memor_dom,'LNNOD'    ,'mesh_type_allocate_minimum',meshe(CURRENT_MESH) % lnnod,2_ip)
    call memory_alloca(memor_dom,'LESUB'    ,'mesh_type_allocate_minimum',meshe(CURRENT_MESH) % lesub,2_ip)
    call memory_alloca(memor_dom,'LMATE'    ,'mesh_type_allocate_minimum',meshe(CURRENT_MESH) % lmate,2_ip)
    
    call memory_alloca(memor_dom,'LNOCH'    ,'mesh_type_allocate_minimum',meshe(CURRENT_MESH) % lnoch,2_ip)
    call memory_alloca(memor_dom,'COORD'    ,'mesh_type_allocate_minimum',meshe(CURRENT_MESH) % coord,ndime,2_ip)
    
    call memory_alloca(memor_dom,'LNODB'    ,'mesh_type_allocate_minimum',meshe(CURRENT_MESH) % lnodb,mnodb,2_ip)
    call memory_alloca(memor_dom,'LTYPB'    ,'mesh_type_allocate_minimum',meshe(CURRENT_MESH) % ltypb,2_ip)
    call memory_alloca(memor_dom,'LBOCH'    ,'mesh_type_allocate_minimum',meshe(CURRENT_MESH) % lboch,2_ip)
    call memory_alloca(memor_dom,'LELBO'    ,'mesh_type_allocate_minimum',meshe(CURRENT_MESH) % lelbo,2_ip)
    
  end subroutine mesh_type_allocate_minimum
  
  !-----------------------------------------------------------------------
  !
  !> @date    15/04/2017
  !> @author  Guillaume Houzeaux
  !> @brief   Allocate and initialize
  !> @details Allocate and initiaize mesh type
  !
  !-----------------------------------------------------------------------

  subroutine mesh_type_allocate_initialize(NUMBER_MESHES)

    integer(ip), optional, intent(in) :: NUMBER_MESHES
    integer(ip)                       :: idivi
    integer(ip)                       :: ndivi_loc

    if( present(NUMBER_MESHES) ) then
       ndivi_loc = NUMBER_MESHES
    else
       ndivi_loc = ndivi
    end if
    
    allocate( meshe(-1:ndivi_loc) )

     do idivi = -1,ndivi_loc

        meshe(idivi) % ndime      = 0
        meshe(idivi) % ntens      = 0
        meshe(idivi) % npoin      = 0
        meshe(idivi) % nelem      = 0
        meshe(idivi) % nboun      = 0
        meshe(idivi) % nedge      = 0
        meshe(idivi) % mnode      = 0
        meshe(idivi) % mnodb      = 0
        meshe(idivi) % mgaus      = 0
        meshe(idivi) % medge      = 0
        meshe(idivi) % nbopo      = 0
        meshe(idivi) % npoi1      = 0
        meshe(idivi) % npoi2      = 0
        meshe(idivi) % npoi3      = 0
        meshe(idivi) % npoin_own  = 0
        meshe(idivi) % npoin_halo = 0
        meshe(idivi) % npoin_2    = 0
        meshe(idivi) % nelem_2    = 0
        meshe(idivi) % nboun_2    = 0
        meshe(idivi) % nzdom     = 0
        meshe(idivi) % nzsym     = 0
        meshe(idivi) % nzedg     = 0
        meshe(idivi) % nzelm_2   = 0
        meshe(idivi) % nzdom_own = 0
        meshe(idivi) % nzdom_ell = 0

        nullify(meshe(idivi) % npoin_par)
        nullify(meshe(idivi) % nelem_par)
        nullify(meshe(idivi) % nboun_par)

        nullify(meshe(idivi) % leinv_loc)
        nullify(meshe(idivi) % lnods)
        nullify(meshe(idivi) % ltype)
        nullify(meshe(idivi) % lelch)
        nullify(meshe(idivi) % lnnod)
        nullify(meshe(idivi) % lesub)
        nullify(meshe(idivi) % lmate)
        nullify(meshe(idivi) % lgaus)

        nullify(meshe(idivi) % lbinv_loc)
        nullify(meshe(idivi) % lnodb)
        nullify(meshe(idivi) % lboel)
        nullify(meshe(idivi) % lelbo)
        nullify(meshe(idivi) % ltypb)
        nullify(meshe(idivi) % lboch)
        nullify(meshe(idivi) % lnnob)

        nullify(meshe(idivi) % lninv_loc)
        nullify(meshe(idivi) % coord)
        nullify(meshe(idivi) % lnoch)
        nullify(meshe(idivi) % lmast)
        nullify(meshe(idivi) % lpoty)

        nullify(meshe(idivi) % edge_to_node)
        nullify(meshe(idivi) % local_to_global_edge)
        nullify(meshe(idivi) % ledgs)
        nullify(meshe(idivi) % ledgb)
        nullify(meshe(idivi) % lnned)
        nullify(meshe(idivi) % lnneb)

        nullify(meshe(idivi) % leset)
        nullify(meshe(idivi) % lbset)
        nullify(meshe(idivi) % lnset)

        nullify(meshe(idivi) % kfl_codno)
        nullify(meshe(idivi) % kfl_codbo)

        nullify(meshe(idivi) % lgrou_dom)

        nullify(meshe(idivi) % linno)

        nullify(meshe(idivi) % c_dom)
        nullify(meshe(idivi) % r_dom)
        nullify(meshe(idivi) % coo_rows)
        nullify(meshe(idivi) % coo_cols)
        nullify(meshe(idivi) % ell_cols)
        nullify(meshe(idivi) % c_sym)
        nullify(meshe(idivi) % r_sym)
        nullify(meshe(idivi) % c_edg)
        nullify(meshe(idivi) % r_edg)
        nullify(meshe(idivi) % c_elm_2)
        nullify(meshe(idivi) % r_elm_2)
        nullify(meshe(idivi) % c_dom_2)
        nullify(meshe(idivi) % r_dom_2)
        nullify(meshe(idivi) % c_dom_own)
        nullify(meshe(idivi) % r_dom_own)
        nullify(meshe(idivi) % r_dom_end)
        nullify(meshe(idivi) % r_dom_ini)

        nullify(meshe(idivi) % exnor)
        nullify(meshe(idivi) % vmass)
        nullify(meshe(idivi) % vmasc)
     end do

     if( INOTSLAVE ) then
        do idivi = 0,ndivi_loc
           allocate( meshe(idivi) % npoin_par(npart) )
           allocate( meshe(idivi) % nelem_par(npart) )
           allocate( meshe(idivi) % nboun_par(npart) )
        end do
     end if

   end subroutine mesh_type_allocate_initialize

   !-----------------------------------------------------------------------
   !>
   !> @date    15/04/2017
   !> @author  Guillaume Houzeaux
   !> @brief   Save original mesh
   !> @details Mesh has been divided and postprocess is on original mesh
   !>          Save the original mesh in MESHE(0)
   !>
   !-----------------------------------------------------------------------

   subroutine mesh_type_save_original_mesh()

     integer(ip) :: ipart,ipoin,ielem,inodb,inode,iboun,idime
     !
     ! Save LNINV_LOC of original mesh
     ! It is useful to refer to original mesh numbering
     !   - IPOIN = lninv_loc(lnlev(JPOIN))
     !     - IPOIN:  original numbering
     !     - JPOIN:  local numbering
     !     - Check that lnlev(JPOIN) /= 0
     ! ISEQUEN does it too to avoid many if's
     ! MASTER allocate minimum memory
     ! Save also LEINV_LOC of original mesh
     !
     if( INOTMASTER ) then

        !if( ISEQUEN ) then
        !   call memory_alloca(memor_dom,'LNINV_LOC','mesh_type_save_original_mesh',lninv_loc,npoin,'IDENTITY')
        !   call memory_alloca(memor_dom,'LEINV_LOC','mesh_type_save_original_mesh',leinv_loc,nelem,'IDENTITY')
        !   call memory_alloca(memor_dom,'LBINV_LOC','mesh_type_save_original_mesh',lbinv_loc,nboun,'IDENTITY')
        !end if

        if( ndivi > 0 ) then
           call memory_alloca(memor_dom,'LNINV_LOC','mesh_type_save_original_mesh',meshe(0) % lninv_loc,npoin)
           do ipoin = 1,npoin
              meshe(0) % lninv_loc(ipoin) = lninv_loc(ipoin)
           end do
           call memory_alloca(memor_dom,'LEINV_LOC','mesh_type_save_original_mesh',meshe(0) % leinv_loc,nelem)
           do ielem = 1,nelem
              meshe(0) % leinv_loc(ielem) = leinv_loc(ielem)
           end do
           call memory_alloca(memor_dom,'LBINV_LOC','mesh_type_save_original_mesh',meshe(0) % lbinv_loc,max(1_ip,nboun))
           do iboun = 1,nboun
              meshe(0) % lbinv_loc(iboun) = lbinv_loc(iboun)
           end do
        else
           meshe(0) % lninv_loc => lninv_loc
           meshe(0) % leinv_loc => leinv_loc
           meshe(0) % lbinv_loc => lbinv_loc
        end if

     else if( IMASTER ) then

        call mesh_type_allocate_minimum(0_ip)
        
     end if

     if( ndivi > 0 .and. kfl_posdi == 0 ) then
        !
        ! Mesh has been divided and postprocess is on original mesh
        !
        if( IMASTER ) then
           !
           ! Master and sequen save global geometry parameters
           !
           meshe(0) % npoin_total = npoin_total
           meshe(0) % nelem_total = nelem_total
           meshe(0) % nboun_total = nboun_total
           do ipart = 1,npart
              meshe(0) % npoin_par(ipart) = npoin_par(ipart)
              meshe(0) % nelem_par(ipart) = nelem_par(ipart)
              meshe(0) % nboun_par(ipart) = nboun_par(ipart)
           end do

        end if

        if( INOTMASTER ) then
           !
           ! Slaves and sequen save original geometry
           !
           meshe(0) % ndime      = ndime
           meshe(0) % ntens      = ntens
           meshe(0) % npoin      = npoin
           meshe(0) % nelem      = nelem
           meshe(0) % nboun      = nboun
           meshe(0) % mnode      = mnode
           meshe(0) % mnodb      = mnodb
           meshe(0) % mgaus      = mgaus
           meshe(0) % nbopo      = nbopo
           meshe(0) % npoi1      = npoi1
           meshe(0) % npoi2      = npoi2
           meshe(0) % npoi3      = npoi3
           meshe(0) % npoin_own  = npoin_own
           meshe(0) % npoin_halo = npoin_halo
           meshe(0) % npoin_2    = npoin
           meshe(0) % nelem_2    = nelem
           meshe(0) % nboun_2    = nboun

           call memory_copy(memor_dom,'LNODS','mesh_type_save_original_mesh',lnods,meshe(0) % lnods,'DO_NOT_DEALLOCATE')
           call memory_copy(memor_dom,'LTYPE','mesh_type_save_original_mesh',ltype,meshe(0) % ltype,'DO_NOT_DEALLOCATE')
           call memory_copy(memor_dom,'LNNOD','mesh_type_save_original_mesh',lnnod,meshe(0) % lnnod,'DO_NOT_DEALLOCATE')
           call memory_copy(memor_dom,'LELCH','mesh_type_save_original_mesh',lelch,meshe(0) % lelch,'DO_NOT_DEALLOCATE')
           call memory_copy(memor_dom,'LESUB','mesh_type_save_original_mesh',lesub,meshe(0) % lesub,'DO_NOT_DEALLOCATE')
           call memory_copy(memor_dom,'LMATE','mesh_type_save_original_mesh',lmate,meshe(0) % lmate,'DO_NOT_DEALLOCATE')

           call memory_copy(memor_dom,'LNOCH','mesh_type_save_original_mesh',lnoch,meshe(0) % lnoch,'DO_NOT_DEALLOCATE')

           call memory_copy(memor_dom,'LNODB','mesh_type_save_original_mesh',lnodb,meshe(0) % lnodb,'DO_NOT_DEALLOCATE')
           call memory_copy(memor_dom,'LTYPB','mesh_type_save_original_mesh',ltypb,meshe(0) % ltypb,'DO_NOT_DEALLOCATE')
           call memory_copy(memor_dom,'LBOCH','mesh_type_save_original_mesh',lboch,meshe(0) % lboch,'DO_NOT_DEALLOCATE')
           call memory_copy(memor_dom,'LELBO','mesh_type_save_original_mesh',lelbo,meshe(0) % lelbo,'DO_NOT_DEALLOCATE')

           call memory_copy(memor_dom,'COORD','mesh_type_save_original_mesh',coord,meshe(0) % coord,'DO_NOT_DEALLOCATE')

        end if
     end if

   end subroutine mesh_type_save_original_mesh

   !-----------------------------------------------------------------------
   !>
   !> @date    15/04/2017
   !> @author  Guillaume Houzeaux
   !> @brief   Update last mesh
   !> @details Update last mesh upon changes and deallocations
   !>          e.g. repoint to last mesh after computing halo geometry
   !>          Called by mod_ghost_geometry.f90
   !>
   !-----------------------------------------------------------------------

   subroutine mesh_type_update_last_mesh(CURRENT_MESH)

     integer(ip), optional, intent(in) :: CURRENT_MESH
     integer(ip)                       :: idivi

     if( present(CURRENT_MESH) ) then
        idivi = CURRENT_MESH
     else
        idivi = ndivi
     end if
     
     if( associated(meshe) .and. size(meshe,kind=ip) >= idivi ) then       
        !
        ! Dimensions
        !
        meshe(idivi) % ndime      =  ndime
        meshe(idivi) % ntens      =  ntens
        meshe(idivi) % npoin      =  npoin
        meshe(idivi) % nelem      =  nelem
        meshe(idivi) % nboun      =  nboun
        meshe(idivi) % mnode      =  mnode
        meshe(idivi) % mnodb      =  mnodb
        meshe(idivi) % mgaus      =  mgaus
        meshe(idivi) % nbopo      =  nbopo     ! Computed in extnor
        meshe(idivi) % npoi1      =  npoi1
        meshe(idivi) % npoi2      =  npoi2
        meshe(idivi) % npoi3      =  npoi3
        meshe(idivi) % npoin_own  =  npoin_own
        meshe(idivi) % npoin_halo =  npoin_halo

        meshe(idivi) % nelem_2    =  nelem_2
        meshe(idivi) % nboun_2    =  nboun_2
        meshe(idivi) % npoin_2    =  npoin_2
        !
        ! Element arrays
        !
        meshe(idivi) % leinv_loc  => leinv_loc
        meshe(idivi) % lnods      => lnods
        meshe(idivi) % ltype      => ltype
        meshe(idivi) % lnnod      => lnnod
        meshe(idivi) % lelch      => lelch
        meshe(idivi) % lesub      => lesub
        meshe(idivi) % lmate      => lmate
        meshe(idivi) % lgaus      => lgaus
        !
        ! Boundary arrays
        !     
        meshe(idivi) % lbinv_loc  => lbinv_loc
        meshe(idivi) % lnodb      => lnodb
        meshe(idivi) % ltypb      => ltypb
        meshe(idivi) % lboch      => lboch
        meshe(idivi) % lboel      => lboel
        meshe(idivi) % lelbo      => lelbo
        meshe(idivi) % lnnob      => lnnob
        !
        ! Nodal arrays
        !     
        meshe(idivi) % lninv_loc  => lninv_loc
        meshe(idivi) % coord      => coord
        meshe(idivi) % lnoch      => lnoch
        meshe(idivi) % lmast      => lmast
        meshe(idivi) % lpoty      => lpoty
        !
        ! Graphs
        !
        meshe(idivi) % nzdom      =  nzdom
        meshe(idivi) % r_dom      => r_dom        ! Computed in domgra
        meshe(idivi) % c_dom      => c_dom        ! Computed in domgra
        !
        ! Geometrical arrays
        !
        meshe(idivi) % exnor      => exnor        ! Computed in extnor
        meshe(idivi) % vmass      => vmass        ! Computed in massma
        meshe(idivi) % vmasc      => vmasc        ! Computed in massmc

     end if

   end subroutine mesh_type_update_last_mesh

 end module mod_mesh_type
!> @}
