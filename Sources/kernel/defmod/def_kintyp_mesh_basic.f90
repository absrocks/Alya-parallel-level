!-----------------------------------------------------------------------
!> @addtogroup Domain
!> @{
!> @file    def_mesh_type.f90
!> @author  houzeaux
!> @date    2020-04-24
!> @brief   Mesh type
!> @details Mesh type definitions
!>
!>          Basic type
!>          ----------
!>          This type is intended to contain all element
!>          types. Thus boundaries are not represented specifically.
!>          This type is similar to gmsh structure.
!>
!>          |-> Basic procedures are defined here. More complex ones
!>              are in the associated module: mod_mesh_type_basic.f90
!>          !-> Postprocess is in mod_postpr_mesh.f90
!>
!>          Extended type
!>          -------------
!>          Lot more info about the mesh. For example, boundaries are
!>          explicitly declared.
!>          This is the native mesh type of Alya.
!>
!>
!-----------------------------------------------------------------------

module def_kintyp_mesh_basic

  use def_kintyp_basic,  only : ip,rp,i1p,i1pp,r3p,lg
  use def_kintyp_comm,   only : comm_data_par_basic
  use mod_memory_basic,  only : memory_alloca
  use mod_memory_basic,  only : memory_deallo
  use mod_memory_basic,  only : memory_copy
  use mod_elmgeo,        only : element_type
  use mod_elmgeo,        only : HEX08_TO_TET04
  use mod_elmgeo,        only : PEN06_TO_TET04
  use mod_elmgeo,        only : PYR05_TO_TET04
  use mod_elmgeo,        only : TET04_TO_TET04
  use def_elmtyp,        only : BAR02
  use def_elmtyp,        only : TRI03
  use def_elmtyp,        only : QUA04
  use def_elmtyp,        only : HEX08
  use def_elmtyp,        only : PEN06
  use def_elmtyp,        only : PYR05
  use def_elmtyp,        only : TET04
  use def_elmtyp,        only : BAR3D
  use def_elmtyp,        only : SHELL
  use mod_maths,         only : maths_heap_sort
  use mod_htable,        only : hash_t
  use mod_htable,        only : htaini
  use mod_htable,        only : htalid
  use mod_htable,        only : htaadd
  use mod_htable,        only : htades
  use mod_htable,        only : htable_initialization
  use def_elmtyp,        only : element_num_ini
  use def_elmtyp,        only : element_num_end

  implicit none
  private

  !----------------------------------------------------------------------
  !
  ! Basic mesh type
  !
  !----------------------------------------------------------------------

  type mesh_type_basic
     character(20)                      :: name                   ! Mesh name
     integer(ip)                        :: id                     ! An identifier
     integer(ip)                        :: ndime                  ! Space dimension
     integer(ip)                        :: mnode                  ! Max number of nodes per element
     integer(ip)                        :: nelem                  ! # elements
     integer(ip)                        :: npoin                  ! # nodes
     integer(ip),               pointer :: lnods(:,:)             ! Element connectivity
     integer(ip),               pointer :: ltype(:)               ! Element type
     integer(ip),               pointer :: leinv_loc(:)           ! Element type
     integer(ip),               pointer :: lninv_loc(:)           ! Node numbering
     real(rp),                  pointer :: coord(:,:)             ! Node coordinates
     
     integer(ip),               pointer :: permn(:)               ! Permutation from a parent mesh
     integer(ip),               pointer :: perme(:)               ! Permutation from a parent mesh
     class(mesh_type_basic),    pointer :: boundary               ! Associated boundary mesh
     class(mesh_type_basic),    pointer :: parent                 ! Parent mesh mesh
     type(comm_data_par_basic)          :: comm                   ! Communication arrays
   contains
     procedure,                 pass    :: init                   ! Initialize all
     procedure,                 pass    :: alloca                 ! Allocate 
     procedure,                 pass    :: alloca_com             ! Allocate communication
     procedure,                 pass    :: deallo                 ! Deallocate
     procedure,                 pass    :: merge                  ! Merge two meshes
     procedure,                 pass    :: append                 ! Append two meshes
     procedure,                 pass    :: dims                   ! Get dimensions from pointer sizes
     procedure,                 pass    :: types                  ! Number of types and list of types
     procedure,                 pass    :: set_boundary           ! Associate a boundary mesh
     procedure,                 pass    :: boundary_mesh          ! Create a boundary mesh
     procedure,                 pass    :: output                 ! Output a mesh
     procedure,                 pass    :: extract                ! Extract a submesh using a mask
     procedure,                 pass    :: extract_boundary       ! Extract a boundary mesh
     procedure,                 pass    :: cartesian_mesh         ! Generate a Cartesian mesh
     procedure,                 pass    :: tag                    ! Put identity numbering
     procedure,                 pass    :: point                  ! Point to a mesh
     procedure,                 pass    :: copy                   ! Copy a mesh (same as equal but with optional args)
     procedure,                 pass    :: copy_comm              ! Copy communication arrays
     procedure,                 pass    :: cut_level              ! Cut a mesh given at a level field equal to zero
     procedure,                 pass    :: equal                  ! Copy a mesh
     generic,   public                  :: assignment(=) => equal
  end type mesh_type_basic

  public :: mesh_type_basic

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Initialization
  !> @details Initialization. Shoudl be called prior to any operations
  !>          on a mesh.
  !> 
  !-----------------------------------------------------------------------

  subroutine init(mesh,wname)

    class(mesh_type_basic)                      :: mesh
    character(len=*),      optional, intent(in) :: wname
    integer(ip)                                 :: wlen

    if( present(wname) ) then
       wlen        = min(len(mesh % name),len(wname))
       mesh % name = trim(wname(1:wlen))
    end if

    mesh % id                    = 0
    mesh % ndime                 = 0
    mesh % mnode                 = 0
    mesh % nelem                 = 0
    mesh % npoin                 = 0
    
    mesh % comm % offset_nelem   = 0
    mesh % comm % offset_npoin   = 0
    mesh % comm % bound_dim      = 0
    mesh % comm % nneig          = 0
    mesh % comm % RANK4          = 0_4
    mesh % comm % SIZE4          = 0_4
    mesh % comm % PAR_COMM_WORLD = 0

    nullify(mesh % lnods)
    nullify(mesh % ltype)
    nullify(mesh % perme)
    nullify(mesh % leinv_loc)
    nullify(mesh % lninv_loc)
    nullify(mesh % coord)
    nullify(mesh % permn)
    nullify(mesh % boundary)
    nullify(mesh % parent)
    
    nullify(mesh % comm % neights)
    nullify(mesh % comm % bound_size)
    nullify(mesh % comm % bound_perm)

  end subroutine init

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Allocate 
  !> @details Allocate a mesh. Dimensions can be given as arguments.
  !>          Itf not, it is assumed they were previously defined.
  !> 
  !-----------------------------------------------------------------------

  subroutine alloca(mesh,ndime,mnode,nelem,npoin,MEMORY_COUNTER)

    class(mesh_type_basic),           intent(inout) :: mesh
    integer(ip),           optional,  intent(in)    :: ndime
    integer(ip),           optional,  intent(in)    :: mnode
    integer(ip),           optional,  intent(in)    :: nelem
    integer(ip),           optional,  intent(in)    :: npoin
    integer(8),            optional,  intent(inout) :: MEMORY_COUNTER(2)
    integer(ip)                                     :: ipoin
    integer(ip)                                     :: ielem
    integer(8)                                      :: memor_loc(2)

    if( present(MEMORY_COUNTER) ) then
       memor_loc = MEMORY_COUNTER
    else
       memor_loc = 0_8
    end if
    
    if( present(ndime) ) mesh % ndime = ndime
    if( present(mnode) ) mesh % mnode = mnode
    if( present(nelem) ) mesh % nelem = nelem
    if( present(npoin) ) mesh % npoin = npoin

    call mesh % deallo(MEMORY_COUNTER)
    
    call memory_alloca(memor_loc,'MESH % LNODS'    ,'alloca',mesh % lnods    ,mesh % mnode,mesh % nelem)
    call memory_alloca(memor_loc,'MESH % LTYPE'    ,'alloca',mesh % ltype    ,mesh % nelem)             
    call memory_alloca(memor_loc,'MESH % PERME'    ,'alloca',mesh % perme    ,mesh % nelem)
    call memory_alloca(memor_loc,'MESH % LEINV_LOC','alloca',mesh % leinv_loc,mesh % nelem)             
    call memory_alloca(memor_loc,'MESH % LNINV_LOC','alloca',mesh % lninv_loc,mesh % npoin)             
    call memory_alloca(memor_loc,'MESH % COORD'    ,'alloca',mesh % coord    ,mesh % ndime,mesh % npoin)
    call memory_alloca(memor_loc,'MESH % PERMN'    ,'alloca',mesh % permn    ,mesh % npoin)

    do ielem = 1,mesh % nelem
       mesh % leinv_loc(ielem) = ielem
    end do
    do ipoin = 1,mesh % npoin
       mesh % lninv_loc(ipoin) = ipoin
    end do
 
    nullify(  mesh % comm % neights ) 
    nullify(  mesh % comm % bound_size ) 
    nullify(  mesh % comm % bound_perm )
    mesh % comm % RANK4          = 0_4
    mesh % comm % SIZE4          = 0_4
    mesh % comm % PAR_COMM_WORLD = 0
    if( present(MEMORY_COUNTER) )  MEMORY_COUNTER = memor_loc

  end subroutine alloca
  subroutine alloca_com(mesh,MEMORY_COUNTER)

    class(mesh_type_basic),           intent(inout) :: mesh
    integer(8),            optional,  intent(inout) :: MEMORY_COUNTER(2)
    integer(ip)                                     :: ipoin
    integer(ip)                                     :: ielem
    integer(8)                                      :: memor_loc(2)
    
    if( present(MEMORY_COUNTER) ) then
       memor_loc = MEMORY_COUNTER
    else
       memor_loc = 0_8
    end if
    
    call memory_alloca(memor_loc,'MESH % COMM % NEIGHTS'   ,'alloca',mesh % comm % neights   ,mesh % comm % nneig)
    call memory_alloca(memor_loc,'MESH % COMM % BOUND_SIZE','alloca',mesh % comm % bound_size,mesh % comm % nneig+1)
    call memory_alloca(memor_loc,'MESH % COMM % BOUND_PERM','alloca',mesh % comm % bound_perm,mesh % comm % bound_dim)

    if( present(MEMORY_COUNTER) )  MEMORY_COUNTER = memor_loc

  end subroutine alloca_com

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Deallocate 
  !> @details Deallocate a mesh. Dimensions are kept just in case they
  !>          would be needed.
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo(mesh,MEMORY_COUNTER,MESH_NAME)

    class(mesh_type_basic)                          :: mesh
    integer(8),            optional,  intent(inout) :: MEMORY_COUNTER(2)
    character(len=*),      optional,  intent(in)    :: MESH_NAME
    character(20)                                   :: my_mesh_name
    integer(8)                                      :: memor_loc(2)

    if( present(MESH_NAME) ) then
       my_mesh_name = trim(MESH_NAME)
    else
       my_mesh_name = 'MESH'
    end if

    if( present(MEMORY_COUNTER) ) then
       memor_loc = MEMORY_COUNTER
    else
       memor_loc = 0_8
    end if
    
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % LNODS'    ,'deallo',mesh % lnods    )
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % LTYPE'    ,'deallo',mesh % ltype    )             
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % PERME'    ,'deallo',mesh % perme    )
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % LEINV_LOC','deallo',mesh % leinv_loc)             
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % LNINV_LOC','deallo',mesh % lninv_loc)             
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % COORD'    ,'deallo',mesh % coord    )
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % PERMN'    ,'deallo',mesh % permn    )

    nullify(  mesh % boundary )
    nullify(  mesh % parent )

    call memory_deallo(memor_loc,trim(my_mesh_name)//' % COMM % NEIGHTS'    ,'deallo',mesh % comm % neights )
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % COMM % BOUND_SIZE' ,'deallo',mesh % comm % bound_size )
    call memory_deallo(memor_loc,trim(my_mesh_name)//' % COMM % BOUND_PERM' ,'deallo',mesh % comm % bound_perm )    
    
    if( present(MEMORY_COUNTER) )  MEMORY_COUNTER = memor_loc

  end subroutine deallo

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux and samaniego
  !> @date    2020-04-24
  !> @brief   Types
  !> @details Return the different element types involved in a mesh.
  !> 
  !-----------------------------------------------------------------------

  subroutine types(mesh,ntype,ltype,num_elements_per_type)

    class(mesh_type_basic), intent(in)    :: mesh
    integer(ip),            intent(out)   :: ntype
    integer(ip), pointer,   intent(inout) :: ltype(:)
    integer(ip), pointer,   intent(inout) :: num_elements_per_type(:)
    integer(ip)                           :: ktype,itype,ielem,iorde
    integer(ip),            allocatable   :: ltype_loc(:)
    integer(ip),            allocatable   :: ltype_sor(:)

    ntype = maxval(mesh % ltype)
    allocate(ltype_loc(ntype))
    allocate(ltype_sor(ntype))
    do itype = 1,ntype
       ltype_loc(itype) = 0
    end do
    
    iorde = 0_ip
    do ielem = 1,mesh % nelem
       itype = mesh % ltype(ielem)
       ! Count the number of elements for a type       
       ltype_loc(itype) = ltype_loc(itype) + 1
       ! Save the order if it is the first time that the type appears            
       if (ltype_loc(itype) == 1) then
          iorde = iorde + 1
          ltype_sor(itype) = iorde                 
       end if
    end do
   
    ntype = count(ltype_loc/=0)
    allocate(num_elements_per_type(ntype))
    allocate(ltype(ntype))

    ! Add types according the order they appear
    do itype = 1,size(ltype_loc)
       if( ltype_loc(itype) /= 0 ) then
          ktype = ltype_sor(itype)
          ltype(ktype) = itype
          num_elements_per_type(ktype) = ltype_loc(itype)
       end if
    end do
    deallocate(ltype_loc)
    deallocate(ltype_sor)

  end subroutine types

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Dimensions of a mesh
  !> @details Get the dimensions of a mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine dims(mesh,ndime,mnode,nelem,npoin)

    class(mesh_type_basic), intent(in)  :: mesh
    integer(ip),            intent(out) :: ndime
    integer(ip),            intent(out) :: mnode
    integer(ip),            intent(out) :: nelem
    integer(ip),            intent(out) :: npoin

    ndime = 0
    mnode = 0
    nelem = 0
    npoin = 0

    if(    associated(mesh % ltype)     .and. &
         & associated(mesh % lnods)     .and. &
         & associated(mesh % leinv_loc) ) then
       mnode = size(mesh % lnods,DIM=1)
       nelem = size(mesh % lnods,DIM=2)
    end if
    if(    associated(mesh % coord)     .and. &
         & associated(mesh % lninv_loc) ) then
       ndime = size(mesh % coord,DIM=1)
       npoin = size(mesh % coord,DIM=2)
    end if

  end subroutine dims

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Copy
  !> @details Copy a mesh NEW_MESH = CPY_MESH
  !> 
  !-----------------------------------------------------------------------
  
  subroutine equal(new_mesh,cpy_mesh)

    class(mesh_type_basic),           intent(in)    :: cpy_mesh
    class(mesh_type_basic),           intent(inout) :: new_mesh

    call new_mesh % copy(cpy_mesh)
    
  end subroutine equal
  
  subroutine copy(new_mesh,cpy_mesh,MEMORY_COUNTER,DO_ALLOCATE,nelem_in,npoin_in)

    class(mesh_type_basic),           intent(in)    :: cpy_mesh
    class(mesh_type_basic),           intent(inout) :: new_mesh
    integer(8),            optional,  intent(inout) :: MEMORY_COUNTER(2)
    logical(lg),           optional,  intent(in)    :: DO_ALLOCATE
    integer(ip),           optional,  intent(in)    :: nelem_in
    integer(ip),           optional,  intent(in)    :: npoin_in
    integer(ip)                                     :: ielem,ipoin
    integer(ip)                                     :: isize,ii
    integer(ip)                                     :: pelem,ppoin
    integer(ip)                                     :: pdime,pnode
    integer(ip)                                     :: nelem_use
    integer(ip)                                     :: npoin_use
    logical(lg)                                     :: if_allocate
    
    call new_mesh % dims(pdime,pnode,pelem,ppoin)

    if( present(DO_ALLOCATE) ) then
       if_allocate = DO_ALLOCATE
    else
       if_allocate = .true.
       if( ppoin >= cpy_mesh % npoin .and. pelem >= cpy_mesh % nelem ) then
          if( pdime == cpy_mesh % ndime .and. pnode >= cpy_mesh % mnode ) then 
             if_allocate = .false.
          end if
       end if
    end if
    !
    ! Nodes and elements to copy
    !
    if( present(npoin_in) ) then
       npoin_use  = npoin_in
    else
       npoin_use  = cpy_mesh % npoin
    end if
    if( present(nelem_in) ) then          
       nelem_use  = nelem_in
    else
       nelem_use  = cpy_mesh % nelem          
    end if    
    !
    ! Allocate if required
    !
    if( if_allocate ) then

       call new_mesh % deallo(MEMORY_COUNTER)

       new_mesh % id     = cpy_mesh % id
       new_mesh % name   = cpy_mesh % name
       new_mesh % ndime  = cpy_mesh % ndime
       new_mesh % mnode  = cpy_mesh % mnode
       if( present(npoin_in) ) then
          new_mesh % npoin  = npoin_in
       else
          new_mesh % npoin  = cpy_mesh % npoin
       end if
       if( present(nelem_in) ) then          
          new_mesh % nelem  = nelem_in
       else
          new_mesh % nelem  = cpy_mesh % nelem          
       end if
       
       call new_mesh % alloca(MEMORY_COUNTER=MEMORY_COUNTER)
       
    end if

    if( cpy_mesh % nelem > 0 ) pnode = min(size(cpy_mesh % lnods,1),size(new_mesh % lnods,1))
    
    do ielem = 1,nelem_use
       new_mesh % lnods(1:pnode,ielem) = cpy_mesh % lnods(1:pnode,ielem) 
       new_mesh % ltype(ielem)         = cpy_mesh % ltype(ielem)
       new_mesh % leinv_loc(ielem)     = cpy_mesh % leinv_loc(ielem)
    end do
    do ipoin = 1,npoin_use
       new_mesh % coord(:,ipoin)   = cpy_mesh % coord(:,ipoin) 
       new_mesh % lninv_loc(ipoin) = cpy_mesh % lninv_loc(ipoin)          
    end do
    
    if( associated(cpy_mesh % perme) ) then
       do ielem = 1,nelem_use
          new_mesh % perme(ielem) = cpy_mesh % perme(ielem)           
       end do
    end if
    if( associated(cpy_mesh % permn) ) then
       do ipoin = 1,npoin_use
          new_mesh % permn(ipoin) = cpy_mesh % permn(ipoin) 
       end do
    end if
    
  end subroutine copy
  
  subroutine copy_comm(new_mesh,cpy_mesh,MEMORY_COUNTER)
    
    class(mesh_type_basic),           intent(in)    :: cpy_mesh
    class(mesh_type_basic),           intent(inout) :: new_mesh
    integer(8),            optional,  intent(inout) :: MEMORY_COUNTER(2)
    integer(ip)                                     :: ii
    !
    ! Communication
    !
    new_mesh % comm % nneig          = cpy_mesh % comm % nneig
    new_mesh % comm % bound_dim      = cpy_mesh % comm % bound_dim
    new_mesh % comm % PAR_COMM_WORLD = cpy_mesh % comm % PAR_COMM_WORLD 
    new_mesh % comm % SIZE4          = cpy_mesh % comm % SIZE4 
    new_mesh % comm % RANK4          = cpy_mesh % comm % RANK4

    call alloca_com(new_mesh,MEMORY_COUNTER=MEMORY_COUNTER)
    do ii = 1,new_mesh % comm % nneig
       new_mesh % comm % neights(ii)    = cpy_mesh % comm % neights(ii)
    end do
    if( associated(cpy_mesh % comm % bound_size) ) then
       do ii = 1,new_mesh % comm % nneig + 1
          new_mesh % comm % bound_size(ii) = cpy_mesh % comm % bound_size(ii)
       end do
    end if
    do ii = 1,new_mesh % comm % bound_dim
       new_mesh % comm % bound_perm(ii) = cpy_mesh % comm % bound_perm(ii)
    end do
    
  end subroutine copy_comm

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Merge
  !> @details Merge a mesh to another. Its like an append operation
  !>          but merging duplicated nodes:
  !>          MESH = MESH  U MESH_GLU
  !>
  !>          MESH_GLU is eventually deallocated 
  !> 
  !-----------------------------------------------------------------------

  subroutine merge(mesh,mesh_glu,MEMORY_COUNTER)

    class(mesh_type_basic),            intent(inout) :: mesh         !< Output mesh
    class(mesh_type_basic),            intent(inout) :: mesh_glu     !< Input mesh (eventually deallocated)
    integer(8),             optional,  intent(inout) :: MEMORY_COUNTER(2)
    type(hash_t)                                     :: htable_3 
    integer(ip)                                      :: lid,ipoin
    integer(ip)                                      :: ielem,kelem
    integer(ip)                                      :: inode,ielty
    logical(lg)                                      :: isin
    type(mesh_type_basic)                            :: mesh_tmp

    call mesh_tmp % init()

    if( mesh_glu % nelem == 0 .and. mesh_glu % npoin == 0 ) then

       return

    else if( mesh % nelem == 0 .and. mesh % npoin == 0 ) then
       !
       ! Copy first mesh
       !
       mesh = mesh_glu
       call mesh_glu % deallo(MEMORY_COUNTER)

    else
       !
       ! Append mesh MESH_GLU to MESH_TMP
       !
       mesh_tmp % ndime = mesh % ndime
       mesh_tmp % mnode = mesh % mnode
       mesh_tmp % nelem = mesh % nelem + mesh_glu % nelem
       mesh_tmp % npoin = mesh % npoin + mesh_glu % npoin

       call htable_initialization(htable_3)
       call htades( htable_3)
       call htaini( htable_3, mesh_tmp % npoin, lidson=.true., AUTOMATIC_SIZE=.true.)
       call htaadd( htable_3, mesh % lninv_loc)
       do ipoin = 1,mesh_glu % npoin
          call htaadd(htable_3,mesh_glu % lninv_loc(ipoin),lid,isin)
       end do
       mesh_tmp % npoin = htable_3 % nelem
       call mesh_tmp % alloca(MEMORY_COUNTER=MEMORY_COUNTER)       
       !
       ! Element arrays
       !
       do ielem = 1,mesh % nelem
          mesh_tmp % lnods(:,ielem)   = mesh % lnods(:,ielem)
          mesh_tmp % ltype(ielem)     = mesh % ltype(ielem)
         mesh_tmp % perme(ielem)     = mesh % perme(ielem)
          mesh_tmp % leinv_loc(ielem) = mesh % leinv_loc(ielem)
       end do
       kelem = mesh % nelem
       do ielem = 1, mesh_glu % nelem
          kelem                       = kelem + 1
          ielty                       = mesh_glu % ltype(ielem)
          mesh_tmp % ltype(kelem)     = mesh_glu % ltype(ielem)
          mesh_tmp % perme(kelem)     = mesh_glu % perme(ielem)
          mesh_tmp % leinv_loc(kelem) = mesh_glu % leinv_loc(ielem)
          do inode = 1,element_type(ielty) % number_nodes
             ipoin                       = mesh_glu % lnods(inode,ielem)
             if( ipoin > 0 ) then
                lid                         = htalid(htable_3,mesh_glu % lninv_loc(ipoin))          
                mesh_tmp % lnods(inode,kelem) = lid
             end if
          end do
       end do
       !
       ! Node arrays
       !
       do ipoin = 1,mesh % npoin
          mesh_tmp % lninv_loc(ipoin) = mesh % lninv_loc(ipoin)
          mesh_tmp % coord(:,ipoin)   = mesh % coord(:,ipoin)
          mesh_tmp % permn(ipoin)     = mesh % permn(ipoin)
       end do
       do ipoin = 1,mesh_glu % npoin
          lid                       = htalid(htable_3,mesh_glu % lninv_loc(ipoin))
          mesh_tmp % coord(:,lid)   = mesh_glu % coord(:,ipoin)
          mesh_tmp % lninv_loc(lid) = mesh_glu % lninv_loc(ipoin)
          mesh_tmp % permn(lid)     = mesh_glu % permn(ipoin)
       end do
       !
       ! MESH = MESH_TMP
       !
       mesh = mesh_tmp
       !
       ! Deallocate
       !
       call mesh_glu % deallo(MEMORY_COUNTER)       
       call mesh_tmp % deallo(MEMORY_COUNTER)       
       call htades( htable_3 )

    end if

  end subroutine merge

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Extract a boundary mesh
  !> @details Extract a boundary mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine extract_boundary(mesh_boundary,mesh,MEMORY_COUNTER)

    class(mesh_type_basic),           intent(inout)  :: mesh_boundary       !< Output mesh (eventually deallocated)
    class(mesh_type_basic),           intent(in)     :: mesh                !< Input mesh
    integer(8),             optional, intent(inout)  :: MEMORY_COUNTER(2)   !< Memory counter
    integer(ip),            pointer                  :: list_types(:)
    integer(ip)                                      :: ntypes,itype,ii

    ntypes = element_num_end(mesh % ndime-1)-element_num_ini(mesh % ndime-1)+1
    allocate(list_types(ntypes))
    ii = 0
    do itype = element_num_ini(mesh % ndime-1),element_num_end(mesh % ndime-1)
       ii = ii + 1
       list_types(ii) = itype
    end do

    call mesh_boundary % deallo()
    call mesh_boundary % init()
    call extract(mesh_boundary,mesh,list_types=list_types,MEMORY_COUNTER=MEMORY_COUNTER)
    
    deallocate(list_types)
    
  end subroutine extract_boundary
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Extract a submesh
  !> @details Extract a submesh from MESH using a mask
  !>          Examples of use:
  !>          call mesh_out % extract(mesh,LIST_TYPES=(/TRI03/))
  !>          call mesh_out % extract(mesh,mask)
  !> 
  !-----------------------------------------------------------------------

  subroutine extract(mesh_out,mesh,mask,list_types,MEMORY_COUNTER,mesh_cmp)

    class(mesh_type_basic),                    intent(inout)  :: mesh_out            !< Output mesh
    class(mesh_type_basic),           target,  intent(in)     :: mesh                !< Input mesh
    logical(lg),            optional, pointer, intent(in)     :: mask(:)             !< Mask
    integer(ip),            optional,          intent(in)     :: list_types(:)       !< List of types to extract
    integer(8),             optional,          intent(inout)  :: MEMORY_COUNTER(2)   !< Memory counter
    class(mesh_type_basic), optional,          intent(inout)  :: mesh_cmp            !< Complementary mesh mesh
    integer(ip)                                               :: ipoin,ielem,kelem
    integer(ip)                                               :: inode,pelty,kpoin
    integer(ip)                                               :: pnode,ii,itype
    integer(ip)                                               :: kelem_cmp,kpoin_cmp
    integer(ip)                                               :: mnode_loc
    logical(lg),            pointer                           :: lmask(:)
    integer(ip),            pointer                           :: permr_nodes(:)
    integer(ip),            pointer                           :: permr_cmp(:)
    integer(8)                                                :: memor_loc(2)

    nullify(permr_nodes,permr_cmp)

    if( present(MEMORY_COUNTER) ) then
       memor_loc = MEMORY_COUNTER
    else
       memor_loc = 0_8
    end if

    if( mesh % nelem == 0 ) then
       !
       ! Just initialize structure and allocate minimum
       !
       call mesh_out % alloca(MEMORY_COUNTER=MEMORY_COUNTER)       
       if( present(mesh_cmp) ) call mesh_cmp % alloca(MEMORY_COUNTER=MEMORY_COUNTER)

    else 
       !
       ! Define mask
       !
       if( present(mask) ) then
          lmask => mask
       else if( present(list_types) ) then
          nullify(lmask)
          call memory_alloca(memor_loc,'LMASK','extract',lmask,mesh % nelem)
          do ii = 1,size(list_types)
             itype = list_types(ii)
             where( mesh % ltype == itype ) lmask = .true.
          end do
       else
          mesh_out = mesh
          if( present(mesh_cmp) ) call mesh_cmp % init()
          return
       end if
       !  
       ! Permutation
       !
       call memory_alloca(memor_loc,'PERMR_NODES','extract',permr_nodes,mesh % npoin)
       if( present(mesh_cmp) ) call memory_alloca(memor_loc,'PERMR_CMP'  ,'extract',permr_cmp  ,mesh % npoin)
       kpoin     = 0
       kelem     = 0
       kpoin_cmp = 0
       kelem_cmp = 0
       mnode_loc = 0

       do ielem = 1,mesh % nelem

          if( lmask(ielem) ) then

             kelem     = kelem + 1
             pelty     = mesh % ltype(ielem)
             pnode     = element_type(pelty) % number_nodes
             mnode_loc = max(pnode,mnode_loc)
             do inode = 1,pnode
                ipoin = mesh % lnods(inode,ielem)
                if( permr_nodes(ipoin) == 0 ) then
                   kpoin = kpoin + 1
                   permr_nodes(ipoin) = kpoin
                end if
             end do

          else if( present(mesh_cmp) ) then

             kelem_cmp = kelem_cmp + 1
             pelty     = mesh % ltype(ielem)
             pnode     = element_type(pelty) % number_nodes
             mnode_loc = max(pnode,mnode_loc)
             do inode = 1,pnode
                ipoin = mesh % lnods(inode,ielem)
                if( permr_cmp(ipoin) == 0 ) then
                   kpoin_cmp = kpoin_cmp + 1
                   permr_cmp(ipoin) = kpoin_cmp
                end if
             end do

          end if
       end do
       !
       ! New mesh
       !
       if( trim(mesh_out % name) == '' ) mesh_out % name  = trim(mesh % name) // '_extract'
       mesh_out % ndime = mesh % ndime
       mesh_out % nelem = kelem
       mesh_out % npoin = kpoin
       mesh_out % mnode = mnode_loc
       call mesh_out % alloca(MEMORY_COUNTER=MEMORY_COUNTER)

       kelem = 0
       do ielem = 1,mesh % nelem
          if( lmask(ielem) ) then
             pelty                           = mesh % ltype(ielem)
             pnode                           = element_type(pelty) % number_nodes
             kelem                           = kelem + 1
             mesh_out % perme(kelem)         = ielem
             mesh_out % ltype(kelem)         = mesh % ltype    (ielem)
             mesh_out % leinv_loc(kelem)     = mesh % leinv_loc(ielem)             
             mesh_out % lnods(1:pnode,kelem) = permr_nodes(mesh % lnods(1:pnode,ielem))
          end if
       end do

       do ipoin = 1,mesh % npoin
          kpoin = permr_nodes(ipoin)
          if( kpoin /= 0 ) mesh_out % permn(kpoin) = ipoin
       end do

       do kpoin = 1,mesh_out % npoin
          ipoin                       = mesh_out % permn(kpoin)
          mesh_out % coord(:,kpoin)   = mesh % coord(:,ipoin)  
          mesh_out % lninv_loc(kpoin) = mesh % lninv_loc(ipoin)
       end do

       call memory_deallo(memor_loc,'PERMR_NODES','extract',permr_nodes)      
       !
       ! Complementary mesh
       !
       if( present(mesh_cmp) ) then
          if( trim(mesh_cmp % name) == '' ) mesh_cmp % name  = trim(mesh % name) // '_complement'
          mesh_cmp % ndime = mesh % ndime
          mesh_cmp % nelem = kelem_cmp
          mesh_cmp % npoin = kpoin_cmp
          mesh_cmp % mnode = mnode_loc
          call mesh_cmp % alloca(MEMORY_COUNTER=MEMORY_COUNTER)

          kelem = 0
          do ielem = 1,mesh % nelem
             if( .not. lmask(ielem) ) then
                pelty                           = mesh % ltype(ielem)
                pnode                           = element_type(pelty) % number_nodes
                kelem                           = kelem + 1
                mesh_cmp % perme(kelem)         = ielem
                mesh_cmp % ltype(kelem)         = mesh % ltype    (ielem)
                mesh_cmp % leinv_loc(kelem)     = mesh % leinv_loc(ielem)             
                mesh_cmp % lnods(1:pnode,kelem) = permr_cmp(mesh % lnods(1:pnode,ielem))
             end if
          end do

          do ipoin = 1,mesh % npoin
             kpoin = permr_cmp(ipoin)
             if( kpoin /= 0 ) mesh_cmp % permn(kpoin) = ipoin
          end do

          do kpoin = 1,mesh_cmp % npoin
             ipoin                       = mesh_cmp % permn(kpoin)
             mesh_cmp % coord(:,kpoin)   = mesh % coord(:,ipoin)  
             mesh_cmp % lninv_loc(kpoin) = mesh % lninv_loc(ipoin)
          end do

          call memory_deallo(memor_loc,'PERMR_CMP','extract',permr_cmp)
          
       end if

    end if
    !
    ! Parent mesh
    !
    mesh_out % parent => mesh
    if( present(mesh_cmp) ) mesh_cmp % parent => mesh
    
    if( present(list_types) ) &
         call memory_deallo(memor_loc,'LMASK','extract',lmask)

    if( present(MEMORY_COUNTER) )  MEMORY_COUNTER = memor_loc

  end subroutine extract
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Append
  !> @details Append a mesh to another. Duplicated nodes are not
  !>          eliminated
  !>          MESH = MESH // MESH_GLU
  !> 
  !-----------------------------------------------------------------------

  subroutine append(mesh,mesh_glu)

    class(mesh_type_basic),  intent(inout) :: mesh         !< Output mesh
    class(mesh_type_basic),  intent(inout) :: mesh_glu     !< Input mesh (eventually deallocated)
    integer(ip)                            :: ipoin
    integer(ip)                            :: ielem,kelem
    integer(ip)                            :: kpoin,inode
    integer(ip)                            :: npoin_lninv
    integer(ip)                            :: nelem_leinv
    integer(ip)                            :: ielty
    type(mesh_type_basic)                  :: mesh_tmp

    if( mesh_glu % nelem == 0 ) then

       call mesh_glu % deallo
       return

    else if( mesh % nelem == 0 ) then
       !
       ! Copy first mesh
       !
       mesh = mesh_glu
       call mesh_glu % deallo

    else 
       !
       ! Append mesh MESH_GLU to MESH_TMP
       !
       call mesh_tmp % init()
       mesh_tmp     = mesh
       mesh % ndime = mesh_tmp % ndime
       mesh % mnode = max(mesh_tmp % mnode,mesh_glu % mnode)
       mesh % nelem = mesh_tmp % nelem + mesh_glu % nelem
       mesh % npoin = mesh_tmp % npoin + mesh_glu % npoin
       npoin_lninv  = maxval(mesh % lninv_loc)
       nelem_leinv  = maxval(mesh % leinv_loc)
       call mesh     % deallo()
       call mesh     % alloca()
       call mesh     % copy(mesh_tmp,DO_ALLOCATE=.false.)
       call mesh_tmp % deallo()
       !
       ! Element arrays
       !
       kelem = mesh_tmp % nelem
       kpoin = mesh_tmp % npoin
       do ielem = 1,mesh_glu % nelem
          kelem = kelem + 1
          ielty = mesh_glu % ltype(ielem)
          do inode = 1,element_type(ielty) % number_nodes
             ipoin = mesh_glu % lnods(inode,ielem)
             mesh % lnods(inode,kelem)  = mesh_glu % lnods(inode,ielem) + kpoin
          end do
          mesh % ltype(kelem)     = mesh_glu % ltype(ielem)
          mesh % leinv_loc(kelem) = mesh_glu % leinv_loc(ielem) + nelem_leinv
          mesh % perme(kelem)     = mesh_glu % perme(ielem)
       end do
       !
       ! Node arrays
       !
       do ipoin = 1,mesh_glu % npoin
          kpoin                   = kpoin + 1
          mesh % lninv_loc(kpoin) = mesh_glu % lninv_loc(ipoin) + npoin_lninv
          mesh % coord(:,kpoin)   = mesh_glu % coord(:,ipoin)
          mesh % permn(kpoin)     = mesh_glu % permn(ipoin)
       end do
       !
       ! Deallocate
       !
       call mesh_glu % deallo()       

    end if

  end subroutine append

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-25
  !> @brief   Output
  !> @details Output mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine output(mesh,lunit,filename)

    class(mesh_type_basic),           intent(inout) :: mesh                !< Mesh type
    integer(ip),            optional, intent(in)    :: lunit               !< Output unit
    character(LEN=*),       optional, intent(in)    :: filename            !< Filename
    integer(ip),  parameter                         :: melty=52
    integer(ip)                                     :: ifirs,inode
    integer(ip)                                     :: ielty,ielem
    integer(ip)                                     :: ipoin,ioerr
    integer(ip)                                     :: ipass,ndime
    character(20)                                   :: intost
    character(200)                                  :: dumml
    character(200)                                  :: title
    logical(lg)                                     :: opened
    logical(lg)                                     :: if_output
    integer(4)                                      :: iunit4

    if( mesh % nelem <= 0 ) return

    if( present(filename) ) then
       title = trim(filename)       
    else if( trim(mesh % name) /= '' ) then
       title = adjustl(trim((mesh % name)))
    else
       title = 'MESH'
    end if

    if( present(lunit) ) then
       iunit4 = int(lunit,KIND=4_ip)
    else
       do iunit4 = 90_4,1000_4
          inquire(unit=iunit4,opened=opened,iostat=ioerr)
          if( ioerr /= 0 )  cycle
          if( .not. opened ) exit
       end do
       open(iunit4,file=trim(title)//'.post.msh',status='unknown')
    end if

    ifirs = 0

    do ipass = 1,2

       do ielty = 1,melty

          if_output = .false.
          if( ipass == 1 ) then
             ndime = mesh % ndime
             if( any(abs(mesh % ltype)==ielty) ) if_output = .true.
          else if( ipass == 2 .and. associated(mesh % boundary) ) then
             ndime = mesh % boundary % ndime
             if( any(abs(mesh % boundary % ltype)==ielty) ) if_output = .true.          
          end if

          if( if_output ) then
             !
             ! Header
             !
             write(intost,*) int(ielty,4)
             dumml = trim(title)//'_'//trim(adjustl(intost))
             write(iunit4,1)&
                  adjustl(trim(dumml)),mesh % ndime,&
                  adjustl(trim(element_type(ielty) % nametopo)),element_type(ielty) % number_nodes
                  !adjustl(trim(cetop(ielty))),element_type(ielty) % number_nodes
             !
             ! Coordinates
             !
             if( ifirs == 0 ) then
                ifirs = 1
                write(iunit4,2) 'coordinates'
                do ipoin = 1,mesh % npoin
                   write(iunit4,3) ipoin,mesh % coord(1:mesh % ndime,ipoin)
                end do
                write(iunit4,2) 'end coordinates'
             end if
             !
             ! Element connectivity
             !
             write(iunit4,2) 'elements'
             if( ipass == 1 ) then
                do ielem = 1,mesh % nelem
                   if( mesh % ltype(ielem) == ielty ) then
                      write(iunit4,4) ielem,(mesh % lnods(inode,ielem),inode=1,element_type(ielty) % number_nodes)
                   end if
                end do
             else
                do ielem = 1,mesh % boundary % nelem
                   if( mesh % boundary % ltype(ielem) == ielty ) then
                      write(iunit4,4) ielem,(mesh % boundary % lnods(inode,ielem),inode=1,element_type(ielty) % number_nodes)
                   end if
                end do
             end if
             write(iunit4,2) 'end elements'
             write(iunit4,2) ''

          end if

       end do
    end do
    !
    ! Close file
    !
    if( .not. present(lunit) ) close(iunit4)

1   format('MESH ',a,' dimension ',i1,' Elemtype ',a,' Nnode ',i2)
2   format(a)
3   format(i9, 3(1x,e16.8e3))
4   format(i9,50(1x,i9))

  end subroutine output

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Boundary mesh 
  !> @details Associated a boundary mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine set_boundary(mesh,mesh_boundary)

    class(mesh_type_basic),         intent(inout) :: mesh
    class(mesh_type_basic), target, intent(in)    :: mesh_boundary
    integer(ip)                                   :: inode,ipoin,ielem
    integer(ip)                                   :: ipoin_ori

    mesh % boundary => mesh_boundary

    do ielem = 1,mesh_boundary % nelem
       do inode = 1,mesh_boundary % mnode
          ipoin = mesh_boundary % lnods(inode,ielem)
          if( ipoin > 0 ) then
             ipoin_ori = mesh_boundary % lninv_loc(ipoin)
             mesh_boundary % lnods(inode,ielem) = ipoin_ori
             !mesh % boundary % coord(idime,ipoin) = mesh % boundary % coord(idime,ipoin) = 
          end if
       end do
    end do

  end subroutine set_boundary

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2018-05-31
  !> @brief   List of boundary nodes
  !> @details List of boundary nodes of a mesh. No parallel exchange!
  !> 
  !-----------------------------------------------------------------------

  subroutine boundary_mesh(mesh_boun,mesh)

    class(mesh_type_basic), intent(inout) :: mesh_boun
    class(mesh_type_basic), intent(in)    :: mesh       
  
    integer(ip)                           :: ielty,ielem,iface,inodf
    integer(ip)                           :: inode,jelem,jface,jelty,ipoin,pnodf
    integer(ip)                           :: ielpo,mface,mnodf,nboun
    integer(ip)                           :: kpoin
    integer(8)                            :: memor_loc(2)
    logical(lg)                           :: equal_faces  

    integer(ip)                           :: mepoi,nlelp
    integer(ip), pointer                  :: lelpo(:)
    integer(ip), pointer                  :: pelpo(:)
    integer(ip), pointer                  :: nepoi(:)

    integer(ip)                           :: nelem
    integer(ip)                           :: npoin
    integer(ip)                           :: ndime
    integer(ip)                           :: mnode
    integer(ip), pointer                  :: lnods(:,:)
    integer(ip), pointer                  :: ltype(:)

    integer(ip)                           :: nfacg
    integer(ip), pointer                  :: facel(:,:,:)
    integer(ip), pointer                  :: permr(:)
    integer(ip), pointer                  :: invpr(:)

    memor_loc = 0_8

    ndime =  mesh % ndime
    npoin =  mesh % npoin
    nelem =  mesh % nelem
    lnods => mesh % lnods
    ltype => mesh % ltype

    nullify(lelpo)
    nullify(pelpo)
    nullify(nepoi)

    nullify(facel)
    nullify(permr)
    nullify(invpr)

    !----------------------------------------------------------------------
    !
    ! LELPO,PELPO: node-to-element graph
    !
    !----------------------------------------------------------------------

    mface = 0
    mnodf = 0
    do ielem = 1,nelem
       ielty = abs(ltype(ielem))
       mface = max(mface,element_type(ielty) % number_faces)
       mnodf = max(mnodf,element_type(ielty) % max_face_nodes)
    end do

    call memory_alloca(memor_loc,'NEPOI','boundary_mesh',nepoi,npoin)
    do ielem = 1,nelem
       ielty = ltype(ielem)
       do inode = 1,element_type(ielty) % number_nodes
          ipoin = lnods(inode,ielem)
          nepoi(ipoin) = nepoi(ipoin) + 1
       end do
    end do
    call memory_alloca(memor_loc,'PELPO','boundary_mesh',pelpo,npoin+1_ip)
    pelpo(1) = 1
    do ipoin = 1,npoin
       pelpo(ipoin+1) = pelpo(ipoin) + nepoi(ipoin)
    end do
    nlelp = pelpo(npoin+1)
    call memory_alloca(memor_loc,'LELPO','boundary_mesh',lelpo,nlelp)
    do ielem = 1,nelem
       ielty = ltype(ielem)
       do inode = 1,element_type(ielty) % number_nodes
          ipoin = lnods(inode,ielem)
          lelpo(pelpo(ipoin)) = ielem
          pelpo(ipoin) = pelpo(ipoin)+1
       end do
    end do
    pelpo(1) =  1
    mepoi    = -1
    do ipoin = 1,npoin
       pelpo(ipoin+1) = pelpo(ipoin) + nepoi(ipoin)
       mepoi = max(mepoi,nepoi(ipoin))
    end do
    call memory_deallo(memor_loc,'NEPOI','boundary_mesh',nepoi)

    !----------------------------------------------------------------------
    !
    ! List of global faces
    !
    !----------------------------------------------------------------------
    !
    ! Allocate memory for lelfa, FACES, CFAEL AND NNODG
    !
    call memory_alloca(memor_loc,'FACEL','boundary_mesh',facel,mnodf+1_ip,mface,nelem)
    !
    ! Construct and sort FACES
    !
    do ielem = 1,nelem                                          
       ielty = abs(ltype(ielem))
       do iface = 1,element_type(ielty) % number_faces
          pnodf = element_type(ielty) % node_faces(iface) 
          do inodf = 1,pnodf 
             inode = element_type(ielty) % list_faces(inodf,iface) 
             facel(inodf,iface,ielem) = lnods(inode,ielem)
          end do
          facel(mnodf+1,iface,ielem) = 1
          call maths_heap_sort(2_ip,pnodf,facel(:,iface,ielem))
       end do
    end do
    !
    ! Compute FACES
    !
    nfacg=0_ip
    do ielem = 1,nelem                                            ! Compare the faces and 
       ielty = abs(ltype(ielem))                                  ! eliminate the repited faces
       do iface = 1,element_type(ielty) % number_faces
          if( facel(mnodf+1,iface,ielem) > 0 ) then
             nfacg = nfacg + 1
             ipoin = facel(1,iface,ielem)
             ielpo = pelpo(ipoin)-1
             do while( ielpo < pelpo(ipoin+1)-1 )
                ielpo = ielpo + 1
                jelem = lelpo(ielpo)
                if( jelem /= ielem ) then
                   jelty = abs(ltype(jelem))                      ! eliminate the repited faces
                   jface = 0
                   do while( jface < element_type(jelty) % number_faces )
                      jface = jface + 1
                      if( facel(mnodf+1,jface,jelem) > 0 ) then
                         equal_faces = .true.
                         inodf = 0
                         do while( equal_faces .and. inodf < element_type(jelty) % node_faces(jface) )
                            inodf = inodf + 1 
                            if( facel(inodf,iface,ielem) /= facel(inodf,jface,jelem) ) equal_faces = .false.
                         end do
                         if( equal_faces ) then
                            facel(mnodf+1,iface,ielem) =  jelem                              ! Keep IELEM face
                            facel(mnodf+1,jface,jelem) = -ielem                              ! Elminate JELEM face
                            facel(      1,iface,ielem) = -jface                              ! Remember IFACE face
                            jface                      =  element_type(jelty) % number_faces ! Exit JFACE do
                            ielpo                      =  pelpo(ipoin+1)                     ! Exit JELEM do  
                         end if
                      end if
                   end do
                end if
             end do
          end if
       end do
    end do
    call memory_deallo(memor_loc,'LELPO','boundary_mesh',lelpo)
    call memory_deallo(memor_loc,'PELPO','boundary_mesh',pelpo)
    
    !----------------------------------------------------------------------
    !
    ! Count boundaries and nodes that are involved
    !
    !----------------------------------------------------------------------
    
    call memory_alloca(memor_loc,'PERMR','boundary_mesh',permr,npoin)
    call memory_alloca(memor_loc,'INVPR','boundary_mesh',invpr,npoin)
    kpoin = 0
    nboun = 0
    mnode = 0
    do ielem = 1,nelem                      
       ielty = abs(ltype(ielem))                    
       do iface = 1,element_type(ielty) % number_faces
          if( facel(mnodf+1,iface,ielem) > 0 .and. facel(1,iface,ielem) > 0 ) then
             pnodf = element_type(ielty) % node_faces(iface)
             mnode = max(mnode,pnodf)
             nboun = nboun + 1
             do inodf = 1,pnodf 
                inode = element_type(ielty) % list_faces(inodf,iface) 
                ipoin = lnods(inode,ielem)
                if( permr(ipoin) == 0 ) then
                   kpoin        = kpoin + 1
                   permr(ipoin) = kpoin
                   invpr(kpoin) = ipoin
                end if
             end do
          end if
       end do
    end do
    
    !----------------------------------------------------------------------
    !
    ! Fill in MESH_BOUN
    !
    !----------------------------------------------------------------------
 
    if( trim(mesh_boun % name) == '' ) then
       mesh_boun % name  = trim(mesh % name) // '_BOUN'
    end if
    mesh_boun % id    = mesh % id
    mesh_boun % ndime = ndime
    mesh_boun % mnode = mnode
    mesh_boun % npoin = kpoin
    mesh_boun % nelem = nboun
    call mesh_boun % alloca()
    
    do kpoin = 1,mesh_boun % npoin
       ipoin = invpr(kpoin)
       mesh_boun % coord(1:ndime,kpoin) = mesh % coord(1:ndime,ipoin)
       mesh_boun % permn(kpoin)         = ipoin
    end do

    nboun = 0
    do ielem = 1,nelem                      
       ielty = abs(ltype(ielem))
       do iface = 1,element_type(ielty) % number_faces
          if( facel(mnodf+1,iface,ielem) > 0 .and. facel(1,iface,ielem) > 0 ) then
             nboun = nboun + 1
             pnodf = element_type(ielty) % node_faces(iface)
             mesh_boun % ltype(nboun) = element_type(ielty) % type_faces(iface)
             mesh_boun % perme(nboun) = ielem
             do inodf = 1,pnodf 
                inode = element_type(ielty) % list_faces(inodf,iface) 
                ipoin = lnods(inode,ielem)
                mesh_boun % lnods(inodf,nboun) = permr(ipoin)
             end do
          end if
       end do
    end do
    !
    ! Deallocate memory
    !
    call memory_deallo(memor_loc,'FACEL','boundary_mesh',facel)
    call memory_deallo(memor_loc,'PERMR','boundary_mesh',permr)
    call memory_deallo(memor_loc,'INVPR','boundary_mesh',invpr)

  end subroutine boundary_mesh

    !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Cartesian Mesh
  !> @details Extract a submesh from MESH using a mask
  !>          Examples of use:
  !>          call mesh_out % extract(mesh,LIST_TYPES=(/TRI03/))
  !>          call mesh_out % extract(mesh,mask)
  !> 
  !-----------------------------------------------------------------------

  subroutine cartesian_mesh(mesh,ndime,boxes,comin,comax)

    class(mesh_type_basic), intent(inout)  :: mesh     !< Output mesh
    integer(ip),            intent(in)     :: ndime    !< Dimension
    integer(ip),            intent(in)     :: boxes(:) !< # elements
    real(rp),               intent(in)     :: comin(:) !< Min coordinates
    real(rp),               intent(in)     :: comax(:) !< Max coordinates
    integer(ip)                            :: ii,jj,kk
    integer(ip)                            :: ipoin
    integer(ip)                            :: kpoin
    integer(ip)                            :: ielem
    integer(ip)                            :: boxer(3)
    real(rp)                               :: delta(3)
    real(rp)                               :: xx,yy,zz

    mesh % ndime = ndime
    if( mesh % ndime == 2 ) then
       mesh % mnode = 4
       mesh % nelem = boxes(1)*boxes(2)
       mesh % npoin = (boxes(1)+1)*(boxes(2)+1)
    else if( mesh % ndime == 3 ) then
       mesh % mnode = 8
       mesh % nelem = boxes(1)*boxes(2)*boxes(3)        
       mesh % npoin = (boxes(1)+1)*(boxes(2)+1)*(boxes(3)+1)
    else
       return
    end if
    delta(1:ndime) = comax(1:ndime)-comin(1:ndime)
    boxer(1:ndime) = real(boxes(1:ndime),rp)
    
    call mesh % alloca()

    ipoin = 0
    ielem = 0
    
    if( ndime == 2 ) then
       do jj = 1,boxes(2)+1
          yy = comin(2) + real(jj-1,rp)/boxer(2) * delta(2)
          do ii = 1,boxes(1)+1
             ipoin = ipoin + 1
             mesh % coord(1,ipoin) = comin(1) + real(ii-1,rp)/boxer(1) * delta(1)
             mesh % coord(2,ipoin) = yy
          end do
       end do
       ipoin = 1
       do jj = 1,boxes(2)
          do ii = 1,boxes(1)
             ielem                 = ielem + 1
             mesh % lnods(1,ielem) = ipoin
             mesh % lnods(2,ielem) = ipoin+1
             mesh % lnods(3,ielem) = ipoin+boxes(1)+2
             mesh % lnods(4,ielem) = ipoin+boxes(1)+1
             mesh % ltype(ielem)   = QUA04
             ipoin                 = ipoin + 1
          end do
          ipoin = ipoin + 1
       end do
    else
       do kk = 1,boxes(3)+1
          zz = comin(3) + real(kk-1,rp)/boxer(3) * delta(3)
          do jj = 1,boxes(2)+1
             yy = comin(2) + real(jj-1,rp)/boxer(2) * delta(2)
             do ii = 1,boxes(1)+1
                ipoin = ipoin + 1
                mesh % coord(1,ipoin) = comin(1) + real(ii-1,rp)/boxer(1) * delta(1)
                mesh % coord(2,ipoin) = yy
                mesh % coord(3,ipoin) = zz
             end do
          end do
       end do
       ipoin = 1
       kpoin = (boxes(1)+1)*(boxes(2)+1)
       do kk = 1,boxes(3)
          do jj = 1,boxes(2)
             do ii = 1,boxes(1)
                ielem                 = ielem + 1
                
                mesh % lnods(1,ielem) = ipoin
                mesh % lnods(2,ielem) = ipoin+1
                mesh % lnods(3,ielem) = ipoin+boxes(1)+2
                mesh % lnods(4,ielem) = ipoin+boxes(1)+1
                mesh % lnods(5,ielem) = kpoin + ipoin
                mesh % lnods(6,ielem) = kpoin + ipoin+1
                mesh % lnods(7,ielem) = kpoin + ipoin+boxes(1)+2
                mesh % lnods(8,ielem) = kpoin + ipoin+boxes(1)+1
                
                mesh % ltype(ielem)   = HEX08
                ipoin                 = ipoin + 1
             end do
             ipoin = ipoin + 1
          end do
          ipoin = kpoin * kk + 1
       end do
    end if
    
  end subroutine cartesian_mesh

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-04-24
  !> @brief   Iag 
  !> @details Assign identity tag to mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine tag(mesh)

    class(mesh_type_basic), intent(inout) :: mesh
    integer(ip)                           :: ielem,ipoin
    
    do ielem = 1,mesh % nelem
       mesh % leinv_loc(ielem) = ielem
    end do
    do ipoin = 1,mesh % npoin
       mesh % lninv_loc(ipoin) = ipoin
    end do
    
  end subroutine tag

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-18
  !> @brief   Point to a mesh
  !> @details Point to a mesh
  !> 
  !-----------------------------------------------------------------------

  subroutine point(mesh,mesh_target)
    
    class(mesh_type_basic), intent(inout) :: mesh
    class(mesh_type_basic), intent(inout) :: mesh_target

    mesh % ndime                 =   mesh_target % ndime
    mesh % mnode                 =   mesh_target % mnode
    mesh % npoin                 =   mesh_target % npoin
    mesh % nelem                 =   mesh_target % nelem
    mesh % lnods                 =>  mesh_target % lnods
    mesh % ltype                 =>  mesh_target % ltype
    mesh % leinv_loc             =>  mesh_target % leinv_loc
    mesh % lninv_loc             =>  mesh_target % lninv_loc
    mesh % coord                 =>  mesh_target % coord

    mesh % comm % RANK4          =   mesh_target % comm % RANK4
    mesh % comm % SIZE4          =   mesh_target % comm % SIZE4
    mesh % comm % PAR_COMM_WORLD =   mesh_target % comm % PAR_COMM_WORLD
    mesh % comm % bound_dim      =   mesh_target % comm % bound_dim
    mesh % comm % nneig          =   mesh_target % comm % nneig
    mesh % comm % neights        =>  mesh_target % comm % neights
    mesh % comm % bound_size     =>  mesh_target % comm % bound_size
    mesh % comm % bound_perm     =>  mesh_target % comm % bound_perm

  end subroutine point

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-09-03
  !> @brief   Compute the cut mesh given a level field
  !> @details Cut a mesh when the level changes sign
  !> 
  !-----------------------------------------------------------------------
  
  subroutine cut_level(mesh_out,mesh_in,level,MEMORY_COUNTER)

    class(mesh_type_basic),          intent(inout) :: mesh_out
    class(mesh_type_basic),          intent(in)    :: mesh_in
    real(rp),               pointer, intent(in)    :: level(:)
    integer(8),   optional,          intent(inout) :: MEMORY_COUNTER(2)
    integer(ip)                                    :: ielem,pelty,pnode
    integer(ip)                                    :: inode,jnode,ipoin
    integer(ip)                                    :: ii,ij,inod1,inod2        
    integer(ip)                                    :: signn,signp,sigtn
    integer(ip)                                    :: pdime,iedge,sigtp   
    integer(ip)                                    :: num_tet,compt,compl
    integer(ip)                                    :: ledgi(10),compg
    real(rp)                                       :: elcod(mesh_in % ndime,mesh_in % mnode)
    real(rp)                                       :: ellev(mesh_in % mnode),inter(3,4)
    real(rp)                                       :: l1,lp,x1,y1,x2,y2,z1,z2
    real(rp)                                       :: ledgr(10)
    logical(lg), pointer                           :: ifcut(:)
    integer(ip), pointer                           :: XXXXX_TO_TET04(:,:)
    integer(8)                                     :: memor_loc(2)
    character(9), parameter                        :: vacal='cut_level'

    if( present(MEMORY_COUNTER) ) then
       memor_loc = MEMORY_COUNTER
    else
       memor_loc = 0_8
    end if

    nullify(ifcut)

    call memory_alloca(memor_loc,'IFCUT',vacal,ifcut,mesh_in % nelem)
    compt = 0_ip
    pdime = mesh_in % ndime
    !
    ! Count number of cut elements
    !
    do ielem = 1,mesh_in % nelem

       pelty = mesh_in % ltype(ielem)
       pnode = element_type(pelty) % number_nodes
       do inode = 1,pnode
          ipoin        = mesh_in % lnods(inode,ielem)
          ellev(inode) = level(ipoin)
       end do

       if( pdime == 2 ) then
          !
          ! 2D elements
          !
          compl = 0_ip
          do iedge = 1,element_type(pelty) % number_edges
             inode = element_type(pelty) % list_edges(1,iedge)
             jnode = element_type(pelty) % list_edges(2,iedge)
             if( ellev(inode)*ellev(jnode) <= 0.0_rp ) compl = compl + 1
          end do

          if( compl >= 2 ) then
             compt = compt + 1
             ifcut(ielem) = .true. 
          end if

       else 
          !
          ! 3D elements
          !
          signn = 0_ip
          signp = 0_ip        
          do inode = 1,pnode
             if( ellev(inode) >= 0.0_rp) then
                signp = signp+1
             else 
                signn = signn+1
             end if
          end do
          if( signp /= pnode .and. signn /= pnode ) then
             if(      pelty == HEX08 ) then
                XXXXX_TO_TET04 => HEX08_TO_TET04
             else if( pelty == PEN06 ) then
                XXXXX_TO_TET04 => PEN06_TO_TET04
             else if( pelty == PYR05 ) then
                XXXXX_TO_TET04 => PYR05_TO_TET04
             else if( pelty == TET04 ) then
                XXXXX_TO_TET04 => TET04_TO_TET04
             else
                call runend('ELEMENT NOT CODED')
             end if
             num_tet = size(XXXXX_TO_TET04,2)
             !
             ! COMPT = Number of interface triangle
             !
             do ii = 1,num_tet
                sigtn = 0_ip
                sigtp = 0_ip                 
                do ij = 1,4
                   if(ellev(XXXXX_TO_TET04(ij,ii))>=0.0_rp) then
                      sigtp=sigtp+1
                   else 
                      sigtn=sigtn+1
                   endif
                end do
                if(      sigtp == 1 .or. sigtn == 1 ) then
                   compt        = compt + 1
                   ifcut(ielem) = .true.
                else if( sigtp == 2 .or. sigtn == 2 ) then
                   compt        = compt + 2
                   ifcut(ielem) = .true.
                endif
             end do
          end if
       end if
    end do
    !
    ! Allocate mesh
    !
    call mesh_out % init()
    mesh_out % nelem = compt
    mesh_out % npoin = pdime * compt
    mesh_out % ndime = pdime
    mesh_out % mnode = pdime
    call mesh_out % alloca() 

    compt = 0_ip
    compg = 0_ip

    do ielem = 1,mesh_in % nelem

       if( ifcut(ielem) ) then

          pelty = mesh_in % ltype(ielem)
          pnode = element_type(pelty) % number_nodes
          do inode = 1,pnode
             ipoin                = mesh_in % lnods(inode,ielem)
             ellev(inode)         = level(ipoin)
             elcod(1:pdime,inode) = mesh_in % coord(1:pdime,ipoin)
          end do

          if( pdime == 2 ) then

             compl = 0_ip           
             do iedge = 1,element_type(pelty) % number_edges
                inode        = element_type(pelty) % list_edges(1,iedge)
                jnode        = element_type(pelty) % list_edges(2,iedge)
                if( ellev(inode) * ellev(jnode) <= 0.0_rp ) then
                   compl        = compl + 1
                   ledgi(compl) = iedge
                   ledgr(compl) = abs(ellev(inode) * ellev(jnode))
                end if
             end do

             call maths_heap_sort(1_ip,compl,ledgr,ledgi)
             !
             ! Compute the intersection of the elements with the surface 
             !
             compt = compt + 1
             do compl = 1,2
                compg                         = compg + 1
                iedge                         = ledgi(compl)
                inode                         = element_type(pelty) % list_edges(1,iedge)
                jnode                         = element_type(pelty) % list_edges(2,iedge)              
                l1                            = abs(ellev(inode))
                lp                            = abs(ellev(inode)-ellev(jnode))
                x1                            = elcod(1,inode)
                x2                            = elcod(1,jnode)
                y1                            = elcod(2,inode)
                y2                            = elcod(2,jnode)
                mesh_out % lnods(compl,compt) = compg
                mesh_out % coord(1,compg)     = x1 * (1.0_rp-l1/lp) + x2 * l1 / lp  
                mesh_out % coord(2,compg)     = y1 * (1.0_rp-l1/lp) + y2 * l1 / lp
             end do

          else 

             signn = 0_ip
             signp = 0_ip
             compl = 0_ip
             inod1 = 0_ip

             do inode = 1,pnode
                if( ellev(inode) >= 0.0_rp ) then
                   signp = signp+1
                else 
                   signn = signn+1
                end if
             end do
             if(      pelty == HEX08 ) then
                XXXXX_TO_TET04 => HEX08_TO_TET04
             else if( pelty == PEN06 ) then
                XXXXX_TO_TET04 => PEN06_TO_TET04
             else if( pelty == PYR05 ) then
                XXXXX_TO_TET04 => PYR05_TO_TET04
             else if( pelty == TET04 ) then
                XXXXX_TO_TET04 => TET04_TO_TET04
             else
                call runend('ELEMENT NOT CODED')
             end if

             if( pelty == TET04 ) then
                !
                ! TET04
                !
                if(signn==1) then

                   !  research of one interface triangle
                   do ij=1,pnode
                      if(ellev(ij)<0.0_rp) then
                         inod1=ij
                      endif
                   end do

                   do ij=1,pnode
                      if(ij/=inod1) then
                         compg=compg+1
                         compl=compl+1
                         if(compl==1) then
                            compt=compt+1
                         endif
                         !
                         ! Compute the intersection of the elements with the surface 
                         !
                         l1=abs(ellev(inod1))
                         lp=abs(ellev(inod1)-ellev(ij))
                         x1=elcod(1,inod1)
                         x2=elcod(1,ij)
                         y1=elcod(2,inod1)
                         y2=elcod(2,ij)
                         z1=elcod(3,inod1)
                         z2=elcod(3,ij)

                         mesh_out % lnods(compl,compt) = compg
                         mesh_out % coord(1,compg)     = x1*(1-l1/lp)+x2*l1/lp  
                         mesh_out % coord(2,compg)     = y1*(1-l1/lp)+y2*l1/lp 
                         mesh_out % coord(3,compg)     = z1*(1-l1/lp)+z2*l1/lp 

                      endif

                   end do

                else if(signp==1) then
                   !  research of one interface triangle
                   do ij=1,pnode
                      if(ellev(ij)>=0.0_rp) then
                         inod1=ij
                      endif
                   end do

                   do ij=1,pnode
                      if(ij/=inod1) then
                         compg=compg+1
                         compl=compl+1
                         if(compl==1) then
                            compt=compt+1
                         endif
                         !
                         ! Compute the intersection of the elements with the surface 
                         !
                         l1=abs(ellev(inod1))
                         lp=abs(ellev(inod1)-ellev(ij))
                         x1=elcod(1,inod1)
                         x2=elcod(1,ij)
                         y1=elcod(2,inod1)
                         y2=elcod(2,ij)
                         z1=elcod(3,inod1)
                         z2=elcod(3,ij)

                         mesh_out % lnods(compl,compt) = compg
                         mesh_out % coord(1,compg)     = x1*(1-l1/lp)+x2*l1/lp  
                         mesh_out % coord(2,compg)     = y1*(1-l1/lp)+y2*l1/lp 
                         mesh_out % coord(3,compg)     = z1*(1-l1/lp)+z2*l1/lp 

                      endif
                   end do

                else if(signp==2) then

                   !  research of two interface triangles
                   do ij=1,4
                      if(ellev(ij)<0.0_rp) then
                         if(inod1==0) then
                            inod1=ij
                         else 
                            inod2=ij
                         endif
                      endif
                   end do

                   do ij=1,4

                      if(ij/=inod1.and.ij/=inod2) then
                         !
                         ! Compute the intersection of the elements with the surface 
                         !

                         compl=compl+1
                         l1=abs(ellev(inod1))
                         lp=abs(ellev(inod1)-ellev(ij))
                         x1=elcod(1,inod1)
                         x2=elcod(1,ij)
                         y1=elcod(2,inod1)
                         y2=elcod(2,ij)
                         z1=elcod(3,inod1)
                         z2=elcod(3,ij)

                         inter(1,compl)=  x1*(1-l1/lp)+x2*l1/lp  
                         inter(2,compl)=  y1*(1-l1/lp)+y2*l1/lp 
                         inter(3,compl)=  z1*(1-l1/lp)+z2*l1/lp 


                         compl=compl+1
                         l1=abs(ellev(inod2))
                         lp=abs(ellev(inod2)-ellev(ij))
                         x1=elcod(1,inod2)
                         x2=elcod(1,ij)
                         y1=elcod(2,inod2)
                         y2=elcod(2,ij)
                         z1=elcod(3,inod2)
                         z2=elcod(3,ij)

                         inter(1,compl)=  x1*(1-l1/lp)+x2*l1/lp  
                         inter(2,compl)=  y1*(1-l1/lp)+y2*l1/lp 
                         inter(3,compl)=  z1*(1-l1/lp)+z2*l1/lp 

                      endif

                   end do

                   compg                     = compg+1
                   compt                     = compt+1
                   mesh_out % lnods(1,compt) = compg
                   mesh_out % coord(1,compg) = inter(1,1)
                   mesh_out % coord(2,compg) = inter(2,1)
                   mesh_out % coord(3,compg) = inter(3,1) 

                   compg                     = compg+1
                   mesh_out % lnods(2,compt) = compg
                   mesh_out % coord(1,compg) = inter(1,2)
                   mesh_out % coord(2,compg) = inter(2,2)
                   mesh_out % coord(3,compg) = inter(3,2) 

                   compg                     = compg+1
                   mesh_out % lnods(3,compt) = compg
                   mesh_out % coord(1,compg) = inter(1,3)
                   mesh_out % coord(2,compg) = inter(2,3)
                   mesh_out % coord(3,compg) = inter(3,3)

                   compg                     = compg+1
                   compt                     = compt+1
                   mesh_out % lnods(1,compt) = compg
                   mesh_out % coord(1,compg) = inter(1,2)
                   mesh_out % coord(2,compg) = inter(2,2)
                   mesh_out % coord(3,compg) = inter(3,2) 

                   compg                     = compg+1
                   mesh_out % lnods(2,compt) = compg
                   mesh_out % coord(1,compg) = inter(1,3)
                   mesh_out % coord(2,compg) = inter(2,3)
                   mesh_out % coord(3,compg) = inter(3,3) 

                   compg                     = compg+1
                   mesh_out % lnods(3,compt) = compg
                   mesh_out % coord(1,compg) = inter(1,4)
                   mesh_out % coord(2,compg) = inter(2,4)
                   mesh_out % coord(3,compg) = inter(3,4) 

                end if

             else if( pelty == HEX08 ) then
                !
                ! HEX08
                !
                if( signp /= pnode .and. signn /= pnode ) then

                   do ii=1,6
                      sigtn = 0_ip
                      sigtp = 0_ip
                      compl = 0_ip
                      inod1 = 0_ip

                      do ij=1,4
                         if(ellev(XXXXX_TO_TET04(ij,ii))>=0.0_rp) then
                            sigtp=sigtp+1
                         else 
                            sigtn=sigtn+1
                         endif
                      end do

                      if(sigtn==1) then
                         !  research of one interface triangle
                         do ij=1,4
                            if(ellev(XXXXX_TO_TET04(ij,ii))<0.0_rp) then
                               inod1=ij
                            endif
                         end do

                         do ij=1,4
                            if(ij/=inod1) then
                               compg=compg+1
                               compl=compl+1
                               if(compl==1) then
                                  compt=compt+1
                               endif
                               !
                               ! Compute the intersection of the elements with the surface 
                               !
                               l1=abs(ellev(XXXXX_TO_TET04(inod1,ii)))
                               lp=abs(ellev(XXXXX_TO_TET04(inod1,ii))-ellev(XXXXX_TO_TET04(ij,ii)))
                               x1=elcod(1,XXXXX_TO_TET04(inod1,ii))
                               x2=elcod(1,XXXXX_TO_TET04(ij,ii))
                               y1=elcod(2,XXXXX_TO_TET04(inod1,ii))
                               y2=elcod(2,XXXXX_TO_TET04(ij,ii))
                               z1=elcod(3,XXXXX_TO_TET04(inod1,ii))
                               z2=elcod(3,XXXXX_TO_TET04(ij,ii))

                               mesh_out % lnods(compl,compt)= compg
                               mesh_out % coord(1,compg) =  x1 * (1.0_rp-l1/lp) + x2 * l1 / lp  
                               mesh_out % coord(2,compg) =  y1 * (1.0_rp-l1/lp) + y2 * l1 / lp 
                               mesh_out % coord(3,compg) =  z1 * (1.0_rp-l1/lp) + z2 * l1 / lp 
                            endif

                         end do

                      else if(sigtp==1) then
                         !  research of one interface triangle
                         do ij=1,4
                            if(ellev(XXXXX_TO_TET04(ij,ii))>=0.0_rp) then
                               inod1=ij
                            endif
                         end do

                         do ij=1,4
                            if(ij/=inod1) then
                               compg=compg+1
                               compl=compl+1
                               if(compl==1) then
                                  compt=compt+1
                               endif
                               !
                               ! Compute the intersection of the elements with the surface 
                               !
                               l1=abs(ellev(XXXXX_TO_TET04(inod1,ii)))
                               lp=abs(ellev(XXXXX_TO_TET04(inod1,ii))-ellev(XXXXX_TO_TET04(ij,ii)))
                               x1=elcod(1,XXXXX_TO_TET04(inod1,ii))
                               x2=elcod(1,XXXXX_TO_TET04(ij,ii))
                               y1=elcod(2,XXXXX_TO_TET04(inod1,ii))
                               y2=elcod(2,XXXXX_TO_TET04(ij,ii))
                               z1=elcod(3,XXXXX_TO_TET04(inod1,ii))
                               z2=elcod(3,XXXXX_TO_TET04(ij,ii))

                               mesh_out % lnods(compl,compt)= compg
                               mesh_out % coord(1,compg)=  x1 * (1.0_rp-l1/lp) + x2 * l1 / lp  
                               mesh_out % coord(2,compg)=  y1 * (1.0_rp-l1/lp) + y2 * l1 / lp 
                               mesh_out % coord(3,compg)=  z1 * (1.0_rp-l1/lp) + z2 * l1 / lp 

                            endif
                         end do

                      else if(sigtn==2) then

                         !  research of two interface triangles
                         do ij=1,4
                            if(ellev(XXXXX_TO_TET04(ij,ii))<0.0_rp) then
                               if(inod1==0) then
                                  inod1=ij
                               else 
                                  inod2=ij
                               endif
                            endif
                         end do

                         do ij=1,4
                            if(ij/=inod1.and.ij/=inod2) then
                               !
                               ! Compute the intersection of the elements with the surface 
                               !

                               compl=compl+1
                               l1=abs(ellev(XXXXX_TO_TET04(inod1,ii)))
                               lp=abs(ellev(XXXXX_TO_TET04(inod1,ii))-ellev(XXXXX_TO_TET04(ij,ii)))
                               x1=elcod(1,XXXXX_TO_TET04(inod1,ii))
                               x2=elcod(1,XXXXX_TO_TET04(ij,ii))
                               y1=elcod(2,XXXXX_TO_TET04(inod1,ii))
                               y2=elcod(2,XXXXX_TO_TET04(ij,ii))
                               z1=elcod(3,XXXXX_TO_TET04(inod1,ii))
                               z2=elcod(3,XXXXX_TO_TET04(ij,ii))

                               inter(1,compl)=  x1*(1-l1/lp)+x2*l1/lp  
                               inter(2,compl)=  y1*(1-l1/lp)+y2*l1/lp 
                               inter(3,compl)=  z1*(1-l1/lp)+z2*l1/lp 


                               compl=compl+1
                               l1=abs(ellev(XXXXX_TO_TET04(inod2,ii)))
                               lp=abs(ellev(XXXXX_TO_TET04(inod2,ii))-ellev(XXXXX_TO_TET04(ij,ii)))
                               x1=elcod(1,XXXXX_TO_TET04(inod2,ii))
                               x2=elcod(1,XXXXX_TO_TET04(ij,ii))
                               y1=elcod(2,XXXXX_TO_TET04(inod2,ii))
                               y2=elcod(2,XXXXX_TO_TET04(ij,ii))
                               z1=elcod(3,XXXXX_TO_TET04(inod2,ii))
                               z2=elcod(3,XXXXX_TO_TET04(ij,ii))

                               inter(1,compl)=  x1*(1-l1/lp)+x2*l1/lp  
                               inter(2,compl)=  y1*(1-l1/lp)+y2*l1/lp 
                               inter(3,compl)=  z1*(1-l1/lp)+z2*l1/lp 

                            endif

                         end do

                         compg                     = compg+1
                         compt                     = compt+1
                         
                         mesh_out % lnods(1,compt) = compg
                         mesh_out % coord(1,compg) = inter(1,1)
                         mesh_out % coord(2,compg) = inter(2,1)
                         mesh_out % coord(3,compg) = inter(3,1) 

                         compg                     = compg+1
                         mesh_out % lnods(2,compt) = compg
                         mesh_out % coord(1,compg) = inter(1,2)
                         mesh_out % coord(2,compg) = inter(2,2)
                         mesh_out % coord(3,compg) = inter(3,2) 

                         compg                     = compg+1
                         mesh_out % lnods(3,compt) = compg
                         mesh_out % coord(1,compg) = inter(1,3)
                         mesh_out % coord(2,compg) = inter(2,3)
                         mesh_out % coord(3,compg) = inter(3,3)

                         compg                     = compg+1
                         compt                     = compt+1
                         
                         mesh_out % lnods(1,compt) = compg
                         mesh_out % coord(1,compg) = inter(1,2)
                         mesh_out % coord(2,compg) = inter(2,2)
                         mesh_out % coord(3,compg) = inter(3,2) 

                         compg                     = compg+1
                         mesh_out % coord(1,compg) = inter(1,3)
                         mesh_out % coord(2,compg) = inter(2,3)
                         mesh_out % coord(3,compg) = inter(3,3) 

                         compg                     = compg+1
                         mesh_out % lnods(3,compt) = compg
                         mesh_out % coord(1,compg) = inter(1,4)
                         mesh_out % coord(2,compg) = inter(2,4)
                         mesh_out % coord(3,compg) = inter(3,4) 

                      end if

                   end do

                end if

             else if( pelty == PEN06 ) then
                !
                ! PEN06
                !
                if( signp /= pnode .and. signn /= pnode ) then

                   do ii=1,3
                      sigtn=0_ip
                      sigtp=0_ip
                      compl=0_ip
                      inod1=0_ip

                      do ij=1,4
                         if(ellev(XXXXX_TO_TET04(ij,ii))>=0.0_rp) then
                            sigtp=sigtp+1
                         else 
                            sigtn=sigtn+1
                         endif
                      end do

                      if(sigtn==1) then
                         !  research of one interface triangle
                         do ij=1,4
                            if(ellev(XXXXX_TO_TET04(ij,ii))<0.0_rp) then
                               inod1=ij
                            endif
                         end do

                         do ij=1,4
                            if(ij/=inod1) then
                               compg=compg+1
                               compl=compl+1
                               if(compl==1) then
                                  compt=compt+1
                               endif
                               !
                               ! Compute the intersection of the elements with the surface 
                               !
                               l1=abs(ellev(XXXXX_TO_TET04(inod1,ii)))
                               lp=abs(ellev(XXXXX_TO_TET04(inod1,ii))-ellev(XXXXX_TO_TET04(ij,ii)))
                               x1=elcod(1,XXXXX_TO_TET04(inod1,ii))
                               x2=elcod(1,XXXXX_TO_TET04(ij,ii))
                               y1=elcod(2,XXXXX_TO_TET04(inod1,ii))
                               y2=elcod(2,XXXXX_TO_TET04(ij,ii))
                               z1=elcod(3,XXXXX_TO_TET04(inod1,ii))
                               z2=elcod(3,XXXXX_TO_TET04(ij,ii))

                               mesh_out % lnods(compl,compt) = compg
                               mesh_out % coord(1,compg)     =  x1 * (1.0_rp-l1/lp) + x2 * l1 / lp  
                               mesh_out % coord(2,compg)     =  y1 * (1.0_rp-l1/lp) + y2 * l1 / lp 
                               mesh_out % coord(3,compg)     =  z1 * (1.0_rp-l1/lp) + z2 * l1 / lp 

                            endif

                         end do

                      else if(sigtp==1) then
                         !  research of one interface triangle
                         do ij=1,4
                            if(ellev(XXXXX_TO_TET04(ij,ii))>=0.0_rp) then
                               inod1=ij
                            endif
                         end do

                         do ij=1,4
                            if(ij/=inod1) then
                               compg=compg+1
                               compl=compl+1
                               if(compl==1) then
                                  compt=compt+1
                               endif
                               !
                               ! Compute the intersection of the elements with the surface 
                               !
                               l1=abs(ellev(XXXXX_TO_TET04(inod1,ii)))
                               lp=abs(ellev(XXXXX_TO_TET04(inod1,ii))-ellev(XXXXX_TO_TET04(ij,ii)))
                               x1=elcod(1,XXXXX_TO_TET04(inod1,ii))
                               x2=elcod(1,XXXXX_TO_TET04(ij,ii))
                               y1=elcod(2,XXXXX_TO_TET04(inod1,ii))
                               y2=elcod(2,XXXXX_TO_TET04(ij,ii))
                               z1=elcod(3,XXXXX_TO_TET04(inod1,ii))
                               z2=elcod(3,XXXXX_TO_TET04(ij,ii))

                               mesh_out % lnods(compl,compt) = compg
                               mesh_out % coord(1,compg)     = x1 * (1.0_rp-l1/lp) + x2 * l1 / lp  
                               mesh_out % coord(2,compg)     = y1 * (1.0_rp-l1/lp) + y2 * l1 / lp 
                               mesh_out % coord(3,compg)     = z1 * (1.0_rp-l1/lp) + z2 * l1 / lp 

                            endif
                         end do

                      else if(sigtn==2) then

                         !  research of two interface triangles
                         do ij=1,4
                            if(ellev(XXXXX_TO_TET04(ij,ii))<0.0_rp) then
                               if(inod1==0) then
                                  inod1=ij
                               else 
                                  inod2=ij
                               endif
                            endif
                         end do

                         do ij=1,4
                            if(ij/=inod1.and.ij/=inod2) then
                               !
                               ! Compute the intersection of the elements with the surface 
                               !

                               compl=compl+1
                               l1=abs(ellev(XXXXX_TO_TET04(inod1,ii)))
                               lp=abs(ellev(XXXXX_TO_TET04(inod1,ii))-ellev(XXXXX_TO_TET04(ij,ii)))
                               x1=elcod(1,XXXXX_TO_TET04(inod1,ii))
                               x2=elcod(1,XXXXX_TO_TET04(ij,ii))
                               y1=elcod(2,XXXXX_TO_TET04(inod1,ii))
                               y2=elcod(2,XXXXX_TO_TET04(ij,ii))
                               z1=elcod(3,XXXXX_TO_TET04(inod1,ii))
                               z2=elcod(3,XXXXX_TO_TET04(ij,ii))

                               inter(1,compl)=  x1*(1-l1/lp)+x2*l1/lp  
                               inter(2,compl)=  y1*(1-l1/lp)+y2*l1/lp 
                               inter(3,compl)=  z1*(1-l1/lp)+z2*l1/lp 

                               compl=compl+1
                               l1=abs(ellev(XXXXX_TO_TET04(inod2,ii)))
                               lp=abs(ellev(XXXXX_TO_TET04(inod2,ii))-ellev(XXXXX_TO_TET04(ij,ii)))
                               x1=elcod(1,XXXXX_TO_TET04(inod2,ii))
                               x2=elcod(1,XXXXX_TO_TET04(ij,ii))
                               y1=elcod(2,XXXXX_TO_TET04(inod2,ii))
                               y2=elcod(2,XXXXX_TO_TET04(ij,ii))
                               z1=elcod(3,XXXXX_TO_TET04(inod2,ii))
                               z2=elcod(3,XXXXX_TO_TET04(ij,ii))

                               inter(1,compl)=  x1*(1-l1/lp)+x2*l1/lp  
                               inter(2,compl)=  y1*(1-l1/lp)+y2*l1/lp 
                               inter(3,compl)=  z1*(1-l1/lp)+z2*l1/lp 

                            endif

                         end do

                         compg                     = compg+1
                         compt                     = compt+1

                         mesh_out % lnods(1,compt) = compg
                         mesh_out % coord(1,compg) = inter(1,1)
                         mesh_out % coord(2,compg) = inter(2,1)
                         mesh_out % coord(3,compg) = inter(3,1) 

                         compg                     = compg+1
                         mesh_out % lnods(2,compt) = compg
                         mesh_out % coord(1,compg) = inter(1,2)
                         mesh_out % coord(2,compg) = inter(2,2)
                         mesh_out % coord(3,compg) = inter(3,2) 

                         compg                     = compg+1
                         mesh_out % lnods(3,compt) = compg
                         mesh_out % coord(1,compg) = inter(1,3)
                         mesh_out % coord(2,compg) = inter(2,3)
                         mesh_out % coord(3,compg) = inter(3,3) 

                         compg                     = compg+1
                         compt                     = compt+1

                         mesh_out % lnods(1,compt) = compg
                         mesh_out % coord(1,compg) = inter(1,2)
                         mesh_out % coord(2,compg) = inter(2,2)
                         mesh_out % coord(3,compg) = inter(3,2) 

                         compg                     = compg+1
                         mesh_out % lnods(2,compt) = compg
                         mesh_out % coord(1,compg) = inter(1,3)
                         mesh_out % coord(2,compg) = inter(2,3)
                         mesh_out % coord(3,compg) = inter(3,3) 

                         compg                     = compg+1
                         mesh_out % lnods(3,compt) = compg
                         mesh_out % coord(1,compg) = inter(1,4)
                         mesh_out % coord(2,compg) = inter(2,4)
                         mesh_out % coord(3,compg) = inter(3,4) 

                      endif


                   end do

                endif

             else if( pelty == PYR05 ) then
                !
                ! PYR05
                !
                 if( signp /= pnode .and. signn /= pnode ) then

                    do ii=1,2
                       sigtn=0_ip
                       sigtp=0_ip
                       compl=0_ip
                       inod1=0_ip

                       do ij=1,4
                          if(ellev(XXXXX_TO_TET04(ij,ii))>=0.0_rp) then
                             sigtp=sigtp+1
                          else 
                             sigtn=sigtn+1
                          endif
                       end do

                       if(sigtn==1) then
                          !  research of one interface triangle
                          do ij=1,4
                             if(ellev(XXXXX_TO_TET04(ij,ii))<0.0_rp) then
                                inod1=ij
                             endif
                          end do

                          do ij=1,4
                             if(ij/=inod1) then
                                compg=compg+1
                                compl=compl+1
                                if(compl==1) then
                                   compt=compt+1
                                endif
                                !
                                ! Compute the intersection of the elements with the surface 
                                !
                                l1=abs(ellev(XXXXX_TO_TET04(inod1,ii)))
                                lp=abs(ellev(XXXXX_TO_TET04(inod1,ii))-ellev(XXXXX_TO_TET04(ij,ii)))
                                x1=elcod(1,XXXXX_TO_TET04(inod1,ii))
                                x2=elcod(1,XXXXX_TO_TET04(ij,ii))
                                y1=elcod(2,XXXXX_TO_TET04(inod1,ii))
                                y2=elcod(2,XXXXX_TO_TET04(ij,ii))
                                z1=elcod(3,XXXXX_TO_TET04(inod1,ii))
                                z2=elcod(3,XXXXX_TO_TET04(ij,ii))

                                mesh_out % lnods(compl,compt) = compg
                                mesh_out % coord(1,compg)     =  x1 * (1.0_rp-l1/lp) + x2 * l1 / lp  
                                mesh_out % coord(2,compg)     =  y1 * (1.0_rp-l1/lp) + y2 * l1 / lp 
                                mesh_out % coord(3,compg)     =  z1 * (1.0_rp-l1/lp) + z2 * l1 / lp 

                             endif

                          end do

                       else if(sigtp==1) then
                          !  research of one interface triangle
                          do ij=1,4
                             if(ellev(XXXXX_TO_TET04(ij,ii))>=0.0_rp) then
                                inod1=ij
                             endif
                          end do

                          do ij=1,4
                             if(ij/=inod1) then
                                compg=compg+1
                                compl=compl+1
                                if(compl==1) then
                                   compt=compt+1
                                endif
                                !
                                ! Compute the intersection of the elements with the surface 
                                !
                                l1=abs(ellev(XXXXX_TO_TET04(inod1,ii)))
                                lp=abs(ellev(XXXXX_TO_TET04(inod1,ii))-ellev(XXXXX_TO_TET04(ij,ii)))
                                x1=elcod(1,XXXXX_TO_TET04(inod1,ii))
                                x2=elcod(1,XXXXX_TO_TET04(ij,ii))
                                y1=elcod(2,XXXXX_TO_TET04(inod1,ii))
                                y2=elcod(2,XXXXX_TO_TET04(ij,ii))
                                z1=elcod(3,XXXXX_TO_TET04(inod1,ii))
                                z2=elcod(3,XXXXX_TO_TET04(ij,ii))

                                mesh_out % lnods(compl,compt) = compg
                                mesh_out % coord(1,compg)     =  x1 * (1.0_rp-l1/lp) + x2 * l1 / lp  
                                mesh_out % coord(2,compg)     =  y1 * (1.0_rp-l1/lp) + y2 * l1 / lp 
                                mesh_out % coord(3,compg)     =  z1 * (1.0_rp-l1/lp) + z2 * l1 / lp 

                             endif
                          end do

                       else if(sigtn==2) then

                          !  research of two interface triangles
                          do ij=1,4
                             if(ellev(XXXXX_TO_TET04(ij,ii))<0.0_rp) then
                                if(inod1==0) then
                                   inod1=ij
                                else 
                                   inod2=ij
                                endif
                             endif
                          end do

                          do ij=1,4
                             if(ij/=inod1.and.ij/=inod2) then
                                !
                                ! Compute the intersection of the elements with the surface 
                                !

                                compl=compl+1
                                l1=abs(ellev(XXXXX_TO_TET04(inod1,ii)))
                                lp=abs(ellev(XXXXX_TO_TET04(inod1,ii))-ellev(XXXXX_TO_TET04(ij,ii)))
                                x1=elcod(1,XXXXX_TO_TET04(inod1,ii))
                                x2=elcod(1,XXXXX_TO_TET04(ij,ii))
                                y1=elcod(2,XXXXX_TO_TET04(inod1,ii))
                                y2=elcod(2,XXXXX_TO_TET04(ij,ii))
                                z1=elcod(3,XXXXX_TO_TET04(inod1,ii))
                                z2=elcod(3,XXXXX_TO_TET04(ij,ii))

                                inter(1,compl)=  x1*(1-l1/lp)+x2*l1/lp  
                                inter(2,compl)=  y1*(1-l1/lp)+y2*l1/lp 
                                inter(3,compl)=  z1*(1-l1/lp)+z2*l1/lp 


                                compl=compl+1
                                l1=abs(ellev(XXXXX_TO_TET04(inod2,ii)))
                                lp=abs(ellev(XXXXX_TO_TET04(inod2,ii))-ellev(XXXXX_TO_TET04(ij,ii)))
                                x1=elcod(1,XXXXX_TO_TET04(inod2,ii))
                                x2=elcod(1,XXXXX_TO_TET04(ij,ii))
                                y1=elcod(2,XXXXX_TO_TET04(inod2,ii))
                                y2=elcod(2,XXXXX_TO_TET04(ij,ii))
                                z1=elcod(3,XXXXX_TO_TET04(inod2,ii))
                                z2=elcod(3,XXXXX_TO_TET04(ij,ii))

                                inter(1,compl)=  x1*(1-l1/lp)+x2*l1/lp  
                                inter(2,compl)=  y1*(1-l1/lp)+y2*l1/lp 
                                inter(3,compl)=  z1*(1-l1/lp)+z2*l1/lp 

                             endif

                          end do

                          compg                     = compg+1
                          compt                     = compt+1

                          mesh_out % lnods(1,compt) = compg
                          mesh_out % coord(1,compg) = inter(1,1)
                          mesh_out % coord(2,compg) = inter(2,1)
                          mesh_out % coord(3,compg) = inter(3,1) 

                          compg                     = compg+1
                          mesh_out % lnods(2,compt) = compg
                          mesh_out % coord(1,compg) = inter(1,2)
                          mesh_out % coord(2,compg) = inter(2,2)
                          mesh_out % coord(3,compg) = inter(3,2) 

                          compg                     = compg+1
                          mesh_out % lnods(3,compt) = compg
                          mesh_out % coord(1,compg) = inter(1,3)
                          mesh_out % coord(2,compg) = inter(2,3)
                          mesh_out % coord(3,compg) = inter(3,3) 


                          compg                     = compg+1
                          compt                     = compt+1
                          mesh_out % lnods(1,compt) = compg
                          mesh_out % coord(1,compg) = inter(1,2)
                          mesh_out % coord(2,compg) = inter(2,2)
                          mesh_out % coord(3,compg) = inter(3,2) 

                          compg                     = compg+1
                          mesh_out % lnods(2,compt) = compg
                          mesh_out % coord(1,compg) = inter(1,3)
                          mesh_out % coord(2,compg) = inter(2,3)
                          mesh_out % coord(3,compg) = inter(3,3) 

                          compg                     = compg+1
                          mesh_out % lnods(3,compt) = compg
                          mesh_out % coord(1,compg) = inter(1,4)
                          mesh_out % coord(2,compg) = inter(2,4)
                          mesh_out % coord(3,compg) = inter(3,4) 


                       endif


                    end do

                 endif

              else

                 call runend('CUT_LEVEL: ELEMENT YPE NOT CODED')
                 
              end if
          end if
       end if
    end do
    !
    ! Give element types
    !
    if( pdime == 2 ) then
       do ielem = 1,mesh_out % nelem
          mesh_out % ltype(ielem) = BAR02 !BAR3D
       end do
    else
       do ielem = 1,mesh_out % nelem
          mesh_out % ltype(ielem) = TRI03 !SHELL
       end do       
    end if
    !
    ! Identity permutation
    !
    do ipoin = 1,mesh_out % npoin
       mesh_out % permn(ipoin) = ipoin
    end do
    do ielem = 1,mesh_out % nelem
       mesh_out % perme(ielem) = ielem
    end do
    !call mesh_out % output(filename='zobi')
    call memory_deallo(memor_loc,'IFCUT',vacal,ifcut)
    
    if( present(MEMORY_COUNTER) )  MEMORY_COUNTER = memor_loc

  end subroutine cut_level
      
end module def_kintyp_mesh_basic
!> @}
