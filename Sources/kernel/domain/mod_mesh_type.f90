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

  use def_kintyp_basic, only : ip,rp
  use def_kintyp_mesh,  only : mesh_type_basic
  use def_kermod,       only : ndivi
  use def_kermod,       only : kfl_posdi

  use def_domain,       only : npoin_origi
  use def_domain,       only : npoin_total
  use def_domain,       only : nelem_total
  use def_domain,       only : mesh_type
  use def_domain,       only : nboun_total
  use def_domain,       only : ndime,ntens,npoin,nfiel
  use def_domain,       only : mfiel,kfl_field
  use def_domain,       only : nelem,nboun,mnode,mnodb
  use def_domain,       only : mgaus,nbopo,npoin_own
  use def_domain,       only : npoin_2,nelem_2
  use def_domain,       only : nboun_2,mnodb
  use def_domain,       only : nelem,nboun
  use def_domain,       only : npoin_halo,npoin
  use def_domain,       only : kfl_ngrou

  use def_master,       only : intost
  use def_master,       only : npoi1,npoi2,npoi3
  use def_master,       only : npoin_par,nelem_par
  use def_master,       only : nboun_par,npart
  use def_master,       only : INOTSLAVE
  use def_master,       only : IPARALL
  use def_master,       only : INOTMASTER
  use def_master,       only : IMASTER,ISEQUEN
  use def_master,       only : intost
  use def_master,       only : title

  use def_domain,       only : memor_dom
  use def_domain,       only : meshe
  use def_master,       only : leinv_loc
  use def_master,       only : lbinv_loc
  use def_master,       only : lninv_loc
  use def_domain,       only : lnods
  use def_domain,       only : ltype
  use def_domain,       only : lnnod
  use def_domain,       only : lelch
  use def_domain,       only : lesub
  use def_domain,       only : lmate
  use def_domain,       only : lgaus
  use def_domain,       only : lnodb
  use def_domain,       only : ltypb
  use def_domain,       only : lboch
  use def_domain,       only : lboel
  use def_domain,       only : lelbo
  use def_domain,       only : lnnob
  use def_domain,       only : coord
  use def_domain,       only : lnoch
  use def_domain,       only : lmast
  use def_domain,       only : lpoty
  use def_domain,       only : xfiel
  use def_domain,       only : nzdom
  use def_domain,       only : r_dom
  use def_domain,       only : c_dom
  use def_domain,       only : exnor
  use def_domain,       only : vmass
  use def_domain,       only : vmasc

  use mod_memory,       only : memory_copy
  use mod_memory,       only : memory_alloca
  use mod_memory,       only : memory_deallo

  use mod_parall,       only : PAR_INITIALIZE_COMMUNICATION_ARRAY
  use mod_parall,       only : PAR_COMM_MY_CODE
  use mod_parall,       only : PAR_MY_CODE_RANK
  use mod_parall,       only : PAR_CODE_SIZE
  use mod_parall,       only : commd

  implicit none 

  private

  integer(ip), target :: int_min_1(1),int_min_2(1,1),int_min_3(1,1,1)
  real(rp),    target :: rea_min_1(1),rea_min_2(1,1),rea_min_3(1,1,1)
  character(100)      :: vacal='mod_mesh_type'

  public :: mesh_type_allocate_initialize
  public :: mesh_type_save_original_mesh
  public :: mesh_type_update_last_mesh
  public :: mesh_type_basic_to_complete
  public :: mesh_type_allocate_minimum
  public :: mesh_type_deallocate
  public :: mesh_type_initialize
  public :: mesh_type_copy

  interface mesh_type_initialize
     module procedure mesh_type_initialize_s,&
          &           mesh_type_initialize_1
  end interface mesh_type_initialize

  interface mesh_type_deallocate
     module procedure mesh_type_deallocate_all,&
          &           mesh_type_deallocate_s
  end interface mesh_type_deallocate

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-05-15
  !> @brief   Copy an extended mesh type
  !> @details Copy an extended mesh typw 
  !> 
  !-----------------------------------------------------------------------

  subroutine mesh_type_copy(meshe_in,meshe_out) 

    type(mesh_type), intent(inout) :: meshe_in
    type(mesh_type), intent(inout) :: meshe_out

    meshe_out % name       = meshe_in % name
    meshe_out % id         = meshe_in % id
    meshe_out % ndime      = meshe_in % ndime
    meshe_out % ntens      = meshe_in % ntens
    meshe_out % npoin      = meshe_in % npoin
    meshe_out % nelem      = meshe_in % nelem
    meshe_out % nboun      = meshe_in % nboun
    meshe_out % nfiel      = meshe_in % nfiel
    meshe_out % kfl_field  = meshe_in % kfl_field
    meshe_out % mnode      = meshe_in % mnode
    meshe_out % mnodb      = meshe_in % mnodb
    meshe_out % mgaus      = meshe_in % mgaus
    meshe_out % nbopo      = meshe_in % nbopo
    meshe_out % npoi1      = meshe_in % npoi1
    meshe_out % npoi2      = meshe_in % npoi2
    meshe_out % npoi3      = meshe_in % npoi3
    meshe_out % npoin_own  = meshe_in % npoin_own
    meshe_out % npoin_halo = meshe_in % npoin_halo
    meshe_out % npoin_2    = meshe_in % npoin
    meshe_out % nelem_2    = meshe_in % nelem
    meshe_out % nboun_2    = meshe_in % nboun
    meshe_out % kfl_ngrou  = meshe_in % kfl_ngrou

    meshe_out % comm % RANK4          = meshe_in % comm % RANK4         
    meshe_out % comm % SIZE4          = meshe_in % comm % SIZE4         
    meshe_out % comm % PAR_COMM_WORLD = meshe_in % comm % PAR_COMM_WORLD

    call memory_copy(memor_dom,'MESH_OUT % LNODS'    ,'mesh_type_save_original_mesh',meshe_in % lnods,    meshe_out % lnods,    'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % LTYPE'    ,'mesh_type_save_original_mesh',meshe_in % ltype,    meshe_out % ltype,    'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % LNNOD'    ,'mesh_type_save_original_mesh',meshe_in % lnnod,    meshe_out % lnnod,    'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % LELCH'    ,'mesh_type_save_original_mesh',meshe_in % lelch,    meshe_out % lelch,    'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % LESUB'    ,'mesh_type_save_original_mesh',meshe_in % lesub,    meshe_out % lesub,    'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % LMATE'    ,'mesh_type_save_original_mesh',meshe_in % lmate,    meshe_out % lmate,    'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % LEINV_LOC','mesh_type_save_original_mesh',meshe_in % leinv_loc,meshe_out % leinv_loc,'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % PERME'    ,'mesh_type_save_original_mesh',meshe_in % perme    ,meshe_out % perme     ,'DO_NOT_DEALLOCATE')

    call memory_copy(memor_dom,'MESH_OUT % COORD'    ,'mesh_type_save_original_mesh',meshe_in % coord,    meshe_out % coord,    'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % LNOCH'    ,'mesh_type_save_original_mesh',meshe_in % lnoch,    meshe_out % lnoch,    'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % LMAST'    ,'mesh_type_save_original_mesh',meshe_in % lmast    ,meshe_out % lmast,    'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % LNINV_LOC','mesh_type_save_original_mesh',meshe_in % lninv_loc,meshe_out % lninv_loc,'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % PERMN'    ,'mesh_type_save_original_mesh',meshe_in % permn    ,meshe_out % permn     ,'DO_NOT_DEALLOCATE')

    call memory_copy(memor_dom,'MESH_OUT % LNODB'    ,'mesh_type_save_original_mesh',meshe_in % lnodb,    meshe_out % lnodb,    'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % LTYPB'    ,'mesh_type_save_original_mesh',meshe_in % ltypb,    meshe_out % ltypb,    'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % LBOCH'    ,'mesh_type_save_original_mesh',meshe_in % lboch,    meshe_out % lboch,    'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % LELBO'    ,'mesh_type_save_original_mesh',meshe_in % lelbo,    meshe_out % lelbo,    'DO_NOT_DEALLOCATE')
    call memory_copy(memor_dom,'MESH_OUT % LBINV_LOC','mesh_type_save_original_mesh',meshe_in % lbinv_loc,meshe_out % lbinv_loc,'DO_NOT_DEALLOCATE')

  end subroutine mesh_type_copy


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

    meshe(CURRENT_MESH) % lninv_loc => int_min_1
    meshe(CURRENT_MESH) % leinv_loc => int_min_1
    meshe(CURRENT_MESH) % lbinv_loc => int_min_1

    meshe(CURRENT_MESH) % lnods     => int_min_2
    meshe(CURRENT_MESH) % ltype     => int_min_1
    meshe(CURRENT_MESH) % lelch     => int_min_1
    meshe(CURRENT_MESH) % lnnod     => int_min_1
    meshe(CURRENT_MESH) % lesub     => int_min_1
    meshe(CURRENT_MESH) % lmate     => int_min_1

    meshe(CURRENT_MESH) % coord     => rea_min_2
    meshe(CURRENT_MESH) % lnoch     => int_min_1
    meshe(CURRENT_MESH) % lmast     => int_min_1

    meshe(CURRENT_MESH) % lnodb     => int_min_2
    meshe(CURRENT_MESH) % ltypb     => int_min_1
    meshe(CURRENT_MESH) % lboch     => int_min_1
    meshe(CURRENT_MESH) % lelbo     => int_min_1

  end subroutine mesh_type_allocate_minimum

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2019-06-11
  !> @brief   Deallocate mesh type
  !> @details Deallocate mesh type
  !> 
  !-----------------------------------------------------------------------

  subroutine mesh_type_deallocate_s(meshe_in,MESH_NAME)

    type(mesh_type),             intent(inout) :: meshe_in
    character(len=*),  optional, intent(in)    :: MESH_NAME
    character(20)                              :: my_mesh_name
    integer(ip)                                :: ifiel

    if( present(MESH_NAME) ) then
       my_mesh_name = trim(MESH_NAME)
    else
       my_mesh_name = 'MESH'
    end if

    call memory_deallo(memor_dom,trim(my_mesh_name)//' % NPOIN_PAR',trim(vacal),meshe_in % npoin_par)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % NELEM_PAR',trim(vacal),meshe_in % nelem_par)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % NBOUN_PAR',trim(vacal),meshe_in % nboun_par)

    call memory_deallo(memor_dom,trim(my_mesh_name)//' % LEINV_LOC',trim(vacal),meshe_in % leinv_loc)               
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lnods',trim(vacal),meshe_in % lnods)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % ltype',trim(vacal),meshe_in % ltype)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lelch',trim(vacal),meshe_in % lelch)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lnnod',trim(vacal),meshe_in % lnnod)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lesub',trim(vacal),meshe_in % lesub)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lmate',trim(vacal),meshe_in % lmate)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lgaus',trim(vacal),meshe_in % lgaus)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % perme',trim(vacal),meshe_in % perme)

    call memory_deallo(memor_dom,trim(my_mesh_name)//' % LBINV_LOC',trim(vacal),meshe_in % lbinv_loc)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lnodb',trim(vacal),meshe_in % lnodb)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lboel',trim(vacal),meshe_in % lboel)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lelbo',trim(vacal),meshe_in % lelbo)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % ltypb',trim(vacal),meshe_in % ltypb)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lboch',trim(vacal),meshe_in % lboch)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lnnob',trim(vacal),meshe_in % lnnob)

    call memory_deallo(memor_dom,trim(my_mesh_name)//' % LNINV_LOC',trim(vacal),meshe_in % lninv_loc)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % coord',trim(vacal),meshe_in % coord)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lnoch',trim(vacal),meshe_in % lnoch)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lmast',trim(vacal),meshe_in % lmast)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lpoty',trim(vacal),meshe_in % lpoty)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % permn',trim(vacal),meshe_in % permn)

    do ifiel = 1,meshe_in % nfiel
       call memory_deallo(memor_dom,trim(my_mesh_name)//' % XFIEL % A',trim(vacal),meshe_in % xfiel(ifiel) % a)
    end do

    call memory_deallo(memor_dom,trim(my_mesh_name)//' % edge_to_node',trim(vacal),meshe_in % edge_to_node)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % local_to_global_edge',trim(vacal),meshe_in % local_to_global_edge)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % ledgs',trim(vacal),meshe_in % ledgs)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % ledgb',trim(vacal),meshe_in % ledgb)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lnned',trim(vacal),meshe_in % lnned)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lnneb',trim(vacal),meshe_in % lnneb)

    call memory_deallo(memor_dom,trim(my_mesh_name)//' % leset',trim(vacal),meshe_in % leset)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lbset',trim(vacal),meshe_in % lbset)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % lnset',trim(vacal),meshe_in % lnset)

    call memory_deallo(memor_dom,trim(my_mesh_name)//' % kfl_codno',trim(vacal),meshe_in % kfl_codno)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % kfl_codbo',trim(vacal),meshe_in % kfl_codbo)

    call memory_deallo(memor_dom,trim(my_mesh_name)//' % LGROU_DOM',trim(vacal),meshe_in % lgrou_dom)

    !call memory_deallo(memor_dom,trim(my_mesh_name)//' % LINNO',trim(vacal),meshe_in % linno)

    call memory_deallo(memor_dom,trim(my_mesh_name)//' % c_dom   ',trim(vacal),meshe_in % c_dom)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % r_dom   ',trim(vacal),meshe_in % r_dom)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % coo_rows',trim(vacal),meshe_in % coo_rows)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % coo_cols',trim(vacal),meshe_in % coo_cols)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % ell_cols',trim(vacal),meshe_in % ell_cols)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % c_sym   ',trim(vacal),meshe_in % c_sym)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % r_sym   ',trim(vacal),meshe_in % r_sym)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % c_edg   ',trim(vacal),meshe_in % c_edg)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % r_edg   ',trim(vacal),meshe_in % r_edg)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % c_elm_2 ',trim(vacal),meshe_in % c_elm_2)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % r_elm_2 ',trim(vacal),meshe_in % r_elm_2)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % c_dom_2 ',trim(vacal),meshe_in % c_dom_2)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % r_dom_2 ',trim(vacal),meshe_in % r_dom_2)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % c_dom_own',trim(vacal),meshe_in % c_dom_own)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % r_dom_own',trim(vacal),meshe_in % r_dom_own)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % r_dom_end',trim(vacal),meshe_in % r_dom_end)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % r_dom_ini',trim(vacal),meshe_in % r_dom_ini)

    call memory_deallo(memor_dom,trim(my_mesh_name)//' % EXNOR',trim(vacal),meshe_in % exnor)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % VMASS',trim(vacal),meshe_in % vmass)
    call memory_deallo(memor_dom,trim(my_mesh_name)//' % VMASC',trim(vacal),meshe_in % vmasc)


  end subroutine mesh_type_deallocate_s

  subroutine mesh_type_deallocate_all(NUMBER_MESHES)

    integer(ip), optional, intent(in) :: NUMBER_MESHES
    integer(ip)                       :: idivi
    integer(ip)                       :: ndivi_loc

    if( present(NUMBER_MESHES) ) then
       ndivi_loc = NUMBER_MESHES
    else
       ndivi_loc = ndivi
    end if

    if( INOTSLAVE ) then
       do idivi = 0,ndivi_loc
          call memory_deallo(memor_dom,'MESHE('//trim(intost(idivi))//') % NPOIN_PAR','mesh_type_save_original_mesh',meshe(idivi) % npoin_par)
          call memory_deallo(memor_dom,'MESHE('//trim(intost(idivi))//') % NELEM_PAR','mesh_type_save_original_mesh',meshe(idivi) % nelem_par)
          call memory_deallo(memor_dom,'MESHE('//trim(intost(idivi))//') % NBOUN_PAR','mesh_type_save_original_mesh',meshe(idivi) % nboun_par)
       end do
    end if
    deallocate( meshe )

  end subroutine mesh_type_deallocate_all

  !-----------------------------------------------------------------------
  !
  !> @date    15/04/2017
  !> @author  Guillaume Houzeaux
  !> @brief   Allocate and initialize
  !> @details Allocate and initiaize mesh type
  !
  !-----------------------------------------------------------------------

  subroutine mesh_type_initialize_1(meshe_in)

    type(mesh_type), pointer, intent(inout) :: meshe_in(:)
    integer(ip)                    :: ii

    if( associated(meshe_in) ) then
       do ii = lbound(meshe_in,1),ubound(meshe_in,1)
          call mesh_type_initialize_s(meshe_in(ii))
       end do
    end if

  end subroutine mesh_type_initialize_1

  subroutine mesh_type_initialize_s(meshe_in)

    type(mesh_type), intent(inout) :: meshe_in
    integer(ip)                    :: ifiel

    call meshe_in % init('COMPLETE_MESH')

    return
    
!!$    meshe_in % name       = ''
!!$    meshe_in % ndime      = 0
!!$    meshe_in % ntens      = 0
!!$    meshe_in % npoin      = 0
!!$    meshe_in % nelem      = 0
!!$    meshe_in % nboun      = 0
!!$    meshe_in % nfiel      = 0
!!$    meshe_in % kfl_field  = 0
!!$    meshe_in % nedge      = 0
!!$    meshe_in % mnode      = 0
!!$    meshe_in % mnodb      = 0
!!$    meshe_in % mgaus      = 0
!!$    meshe_in % medge      = 0
!!$    meshe_in % nbopo      = 0
!!$    meshe_in % npoi1      = 0
!!$    meshe_in % npoi2      = 0
!!$    meshe_in % npoi3      = 0
!!$    meshe_in % npoin_own  = 0
!!$    meshe_in % npoin_halo = 0
!!$    meshe_in % npoin_2    = 0
!!$    meshe_in % nelem_2    = 0
!!$    meshe_in % nboun_2    = 0
!!$    meshe_in % nzdom      = 0
!!$    meshe_in % nzsym      = 0
!!$    meshe_in % nzedg      = 0
!!$    meshe_in % nzelm_2    = 0
!!$    meshe_in % nzdom_own  = 0
!!$    meshe_in % nzdom_ell  = 0
!!$    meshe_in % kfl_ngrou  = 0
!!$
!!$    nullify(meshe_in % npoin_par)
!!$    nullify(meshe_in % nelem_par)
!!$    nullify(meshe_in % nboun_par)
!!$
!!$    nullify(meshe_in % leinv_loc)
!!$    nullify(meshe_in % lnods)
!!$    nullify(meshe_in % ltype)
!!$    nullify(meshe_in % lelch)
!!$    nullify(meshe_in % lnnod)
!!$    nullify(meshe_in % lesub)
!!$    nullify(meshe_in % lmate)
!!$    nullify(meshe_in % lgaus)
!!$    nullify(meshe_in % perme)
!!$
!!$    nullify(meshe_in % lbinv_loc)
!!$    nullify(meshe_in % lnodb)
!!$    nullify(meshe_in % lboel)
!!$    nullify(meshe_in % lelbo)
!!$    nullify(meshe_in % ltypb)
!!$    nullify(meshe_in % lboch)
!!$    nullify(meshe_in % lnnob)
!!$
!!$    nullify(meshe_in % lninv_loc)
!!$    nullify(meshe_in % coord)
!!$    nullify(meshe_in % lnoch)
!!$    nullify(meshe_in % lmast)
!!$    nullify(meshe_in % lpoty)
!!$    nullify(meshe_in % permn)
!!$
!!$    do ifiel = 1,meshe_in % nfiel
!!$       nullify(meshe_in % xfiel(ifiel) % a)
!!$    end do
!!$
!!$    nullify(meshe_in % edge_to_node)
!!$    nullify(meshe_in % local_to_global_edge)
!!$    nullify(meshe_in % ledgs)
!!$    nullify(meshe_in % ledgb)
!!$    nullify(meshe_in % lnned)
!!$    nullify(meshe_in % lnneb)
!!$
!!$    nullify(meshe_in % leset)
!!$    nullify(meshe_in % lbset)
!!$    nullify(meshe_in % lnset)
!!$
!!$    nullify(meshe_in % kfl_codno)
!!$    nullify(meshe_in % kfl_codbo)
!!$
!!$    nullify(meshe_in % lgrou_dom)
!!$
!!$    nullify(meshe_in % linno)
!!$
!!$    nullify(meshe_in % c_dom)
!!$    nullify(meshe_in % r_dom)
!!$    nullify(meshe_in % coo_rows)
!!$    nullify(meshe_in % coo_cols)
!!$    nullify(meshe_in % ell_cols)
!!$    nullify(meshe_in % c_sym)
!!$    nullify(meshe_in % r_sym)
!!$    nullify(meshe_in % c_edg)
!!$    nullify(meshe_in % r_edg)
!!$    nullify(meshe_in % c_elm_2)
!!$    nullify(meshe_in % r_elm_2)
!!$    nullify(meshe_in % c_dom_2)
!!$    nullify(meshe_in % r_dom_2)
!!$    nullify(meshe_in % c_dom_own)
!!$    nullify(meshe_in % r_dom_own)
!!$    nullify(meshe_in % r_dom_end)
!!$    nullify(meshe_in % r_dom_ini)
!!$
!!$    nullify(meshe_in % exnor)
!!$    nullify(meshe_in % vmass)
!!$    nullify(meshe_in % vmasc)
!!$
!!$    if( IPARALL ) call PAR_INITIALIZE_COMMUNICATION_ARRAY(meshe_in % comm)

  end subroutine mesh_type_initialize_s

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
    integer(ip)                       :: idivi,ifiel
    integer(ip)                       :: ndivi_loc

    if( present(NUMBER_MESHES) ) then
       ndivi_loc = NUMBER_MESHES
    else
       ndivi_loc = ndivi
    end if

    allocate( meshe(-1:ndivi_loc) )

    do idivi = -1,ndivi_loc
       call mesh_type_initialize(meshe(idivi))       
       meshe(idivi) % name = trim(title)//'_MM'//trim(intost(idivi))
    end do

    if( INOTSLAVE ) then
       do idivi = 0,ndivi_loc
          call memory_alloca(memor_dom,'MESHE('//trim(intost(idivi))//') % NPOIN_PAR','mesh_type_save_original_mesh',meshe(idivi) % npoin_par,npart)
          call memory_alloca(memor_dom,'MESHE('//trim(intost(idivi))//') % NELEM_PAR','mesh_type_save_original_mesh',meshe(idivi) % nelem_par,npart)
          call memory_alloca(memor_dom,'MESHE('//trim(intost(idivi))//') % NBOUN_PAR','mesh_type_save_original_mesh',meshe(idivi) % nboun_par,npart) 
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

    integer(ip) :: ipart,ipoin,ielem,iboun
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
          meshe(0) % npoin_origi = npoin_origi
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
          meshe(0) % nfiel      = nfiel
          meshe(0) % kfl_field  = kfl_field
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
          meshe(0) % kfl_ngrou  = kfl_ngrou

          call memory_copy(memor_dom,'LNODS','mesh_type_save_original_mesh',lnods,meshe(0) % lnods,'DO_NOT_DEALLOCATE')
          call memory_copy(memor_dom,'LTYPE','mesh_type_save_original_mesh',ltype,meshe(0) % ltype,'DO_NOT_DEALLOCATE')
          call memory_copy(memor_dom,'LNNOD','mesh_type_save_original_mesh',lnnod,meshe(0) % lnnod,'DO_NOT_DEALLOCATE')
          call memory_copy(memor_dom,'LELCH','mesh_type_save_original_mesh',lelch,meshe(0) % lelch,'DO_NOT_DEALLOCATE')
          call memory_copy(memor_dom,'LESUB','mesh_type_save_original_mesh',lesub,meshe(0) % lesub,'DO_NOT_DEALLOCATE')
          call memory_copy(memor_dom,'LMATE','mesh_type_save_original_mesh',lmate,meshe(0) % lmate,'DO_NOT_DEALLOCATE')

          call memory_copy(memor_dom,'COORD','mesh_type_save_original_mesh',coord,meshe(0) % coord,'DO_NOT_DEALLOCATE')
          call memory_copy(memor_dom,'LNOCH','mesh_type_save_original_mesh',lnoch,meshe(0) % lnoch,'DO_NOT_DEALLOCATE')
          call memory_copy(memor_dom,'LMAST','mesh_type_save_original_mesh',lmast,meshe(0) % lmast,'DO_NOT_DEALLOCATE')

          call memory_copy(memor_dom,'LNODB','mesh_type_save_original_mesh',lnodb,meshe(0) % lnodb,'DO_NOT_DEALLOCATE')
          call memory_copy(memor_dom,'LTYPB','mesh_type_save_original_mesh',ltypb,meshe(0) % ltypb,'DO_NOT_DEALLOCATE')
          call memory_copy(memor_dom,'LBOCH','mesh_type_save_original_mesh',lboch,meshe(0) % lboch,'DO_NOT_DEALLOCATE')
          call memory_copy(memor_dom,'LELBO','mesh_type_save_original_mesh',lelbo,meshe(0) % lelbo,'DO_NOT_DEALLOCATE')

       end if
       !
       ! Communication
       !
       meshe(0) % comm % RANK4          = int(PAR_MY_CODE_RANK,4)        
       meshe(0) % comm % SIZE4          = int(PAR_CODE_SIZE,4) 
       meshe(0) % comm % PAR_COMM_WORLD = PAR_COMM_MY_CODE

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
    integer(ip)                       :: idivi,ipart,ifiel

    if( present(CURRENT_MESH) ) then
       idivi = CURRENT_MESH
    else
       idivi = ndivi
    end if

    if( associated(meshe) .and. size(meshe,kind=ip) >= idivi ) then       
       !
       ! Dimensions
       !
       meshe(idivi) % ndime       =  ndime
       meshe(idivi) % ntens       =  ntens
       meshe(idivi) % npoin       =  npoin
       meshe(idivi) % nelem       =  nelem
       meshe(idivi) % nboun       =  nboun
       meshe(idivi) % nfiel       =  nfiel
       meshe(idivi) % kfl_field   =  kfl_field
       meshe(idivi) % mnode       =  mnode
       meshe(idivi) % mnodb       =  mnodb
       meshe(idivi) % mgaus       =  mgaus
       meshe(idivi) % nbopo       =  nbopo     ! Computed in extnor
       meshe(idivi) % npoi1       =  npoi1
       meshe(idivi) % npoi2       =  npoi2
       meshe(idivi) % npoi3       =  npoi3
       meshe(idivi) % npoin_own   =  npoin_own
       meshe(idivi) % npoin_halo  =  npoin_halo
       meshe(idivi) % kfl_ngrou   =  kfl_ngrou

       meshe(idivi) % nelem_2     =  nelem_2
       meshe(idivi) % nboun_2     =  nboun_2
       meshe(idivi) % npoin_2     =  npoin_2
       !
       ! Dimensions of master
       !
       if( INOTSLAVE ) then
          meshe(ndivi) % npoin_origi = npoin_origi
          meshe(ndivi) % npoin_total = npoin_total
          meshe(ndivi) % nelem_total = nelem_total
          meshe(ndivi) % nboun_total = nboun_total
          do ipart = 1,npart
             meshe(ndivi) % npoin_par(ipart) = npoin_par(ipart)
             meshe(ndivi) % nelem_par(ipart) = nelem_par(ipart)
             meshe(ndivi) % nboun_par(ipart) = nboun_par(ipart)
          end do
       end if
       !
       ! Element arrays
       !
       meshe(idivi) % leinv_loc   => leinv_loc
       meshe(idivi) % lnods       => lnods
       meshe(idivi) % ltype       => ltype
       meshe(idivi) % lnnod       => lnnod
       meshe(idivi) % lelch       => lelch
       meshe(idivi) % lesub       => lesub
       meshe(idivi) % lmate       => lmate
       meshe(idivi) % lgaus       => lgaus
       !
       ! Boundary arrays
       !     
       meshe(idivi) % lbinv_loc   => lbinv_loc
       meshe(idivi) % lnodb       => lnodb
       meshe(idivi) % ltypb       => ltypb
       meshe(idivi) % lboch       => lboch
       meshe(idivi) % lboel       => lboel
       meshe(idivi) % lelbo       => lelbo
       meshe(idivi) % lnnob       => lnnob
       !
       ! Nodal arrays
       !     
       meshe(idivi) % lninv_loc   => lninv_loc
       meshe(idivi) % coord       => coord
       meshe(idivi) % lnoch       => lnoch
       meshe(idivi) % lmast       => lmast
       meshe(idivi) % lpoty       => lpoty
       !
       ! Others
       !
       do ifiel = 1,meshe(idivi) % nfiel
          meshe(idivi) % xfiel(ifiel) % a => xfiel(ifiel) % a
       end do
       !
       ! Graphs
       !
       meshe(idivi) % nzdom       =  nzdom
       meshe(idivi) % r_dom       => r_dom        ! Computed in domgra
       meshe(idivi) % c_dom       => c_dom        ! Computed in domgra
       !
       ! Geometrical arrays
       !
       meshe(idivi) % exnor       => exnor        ! Computed in extnor
       meshe(idivi) % vmass       => vmass        ! Computed in massma
       meshe(idivi) % vmasc       => vmasc        ! Computed in massmc
       !
       ! Communication
       !
       meshe(idivi) % comm % RANK4          =  int(PAR_MY_CODE_RANK,4)        
       meshe(idivi) % comm % SIZE4          =  int(PAR_CODE_SIZE,4)    
       meshe(idivi) % comm % PAR_COMM_WORLD =  PAR_COMM_MY_CODE
       if( IPARALL ) then
          meshe(idivi) % comm % nneig          =  commd % nneig
          meshe(idivi) % comm % bound_dim      =  commd % bound_dim
          meshe(idivi) % comm % neights        => commd % neights
          meshe(idivi) % comm % bound_perm     => commd % bound_perm
          meshe(idivi) % comm % bound_size     => commd % bound_size
       end if

    end if

  end subroutine mesh_type_update_last_mesh

  subroutine mesh_type_basic_to_complete(mesh_basic,mesh_complete)

    type(mesh_type_basic), intent(inout) :: mesh_basic
    type(mesh_type),       intent(inout) :: mesh_complete
    integer(ip)                          :: iboun,ii,ipoin

    mesh_complete % ndime                 =   mesh_basic % ndime
    mesh_complete % mnode                 =   mesh_basic % mnode
    mesh_complete % mnodb                 =   mnodb
    mesh_complete % npoin                 =   mesh_basic % npoin
    mesh_complete % nelem                 =   mesh_basic % nelem
    mesh_complete % lnods                 =>  mesh_basic % lnods
    mesh_complete % ltype                 =>  mesh_basic % ltype
    mesh_complete % leinv_loc             =>  mesh_basic % leinv_loc
    mesh_complete % lninv_loc             =>  mesh_basic % lninv_loc
    mesh_complete % coord                 =>  mesh_basic % coord

    mesh_complete % comm % RANK4          =   mesh_basic % comm % RANK4
    mesh_complete % comm % SIZE4          =   mesh_basic % comm % SIZE4
    mesh_complete % comm % PAR_COMM_WORLD =   mesh_basic % comm % PAR_COMM_WORLD
    mesh_complete % comm % bound_dim      =   mesh_basic % comm % bound_dim
    mesh_complete % comm % nneig          =   mesh_basic % comm % nneig
    mesh_complete % comm % neights        =>  mesh_basic % comm % neights
    mesh_complete % comm % bound_size     =>  mesh_basic % comm % bound_size
    mesh_complete % comm % bound_perm     =>  mesh_basic % comm % bound_perm

  end subroutine mesh_type_basic_to_complete

end module mod_mesh_type
!> @}
