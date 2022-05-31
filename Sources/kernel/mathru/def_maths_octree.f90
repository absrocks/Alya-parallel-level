!-----------------------------------------------------------------------
!
!> @defgroup Maths_Toolbox
!> Toolbox for bins
!> @{
!> @name    ToolBox for bins
!> @file    def_maths.f90
!> @author  Guillaume Houzeaux
!> @brief   Octree
!> @details Octree construction
!
!-----------------------------------------------------------------------

module def_maths_octree

  use def_kintyp_basic, only : ip,rp,lg
  use def_elmtyp,       only : QUA04
  use def_elmtyp,       only : HEX08
  use mod_memory_basic, only : memory_alloca
  use mod_memory_basic, only : memory_deallo
  use mod_memory_basic, only : memory_size

  real(rp), parameter :: epsil = epsilon(1.0_rp)
  
  type octbox
     integer(ip)               :: id          ! My global ID
     integer(ip)               :: idfilled    ! My global ID in filled bins
     integer(ip)               :: level       ! Generation
     integer(ip)               :: npoinbox    ! Number of nodes
     integer(ip)               :: childid     ! Child ID (1->4 or 1->8)
     integer(ip)               :: whoiam      ! Father or have nodes
     integer(ip),  pointer     :: nodes(:)    
     real(rp)                  :: minc(3)     ! Min coordinates
     real(rp)                  :: maxc(3)     ! Max coordinates
     type(octbox), pointer     :: parent      ! Pointer to parent
     type(octbox), pointer     :: children(:) ! Pointer to children
  end type octbox
  type(octbox)   :: octbox_init = octbox(&
       0_ip,&                                 ! id         
       0_ip,&                                 ! id filled         
       0_ip,&                                 ! level      
       0_ip,&                                 ! npoinbox   
       0_ip,&                                 ! childid    
       0_ip,&                                 ! whoiam     
       null(),&                               ! nodes(:)   
       (/0.0_rp,0.0_rp,0.0_rp/),&             ! minc(3)    
       (/0.0_rp,0.0_rp,0.0_rp/),&             ! maxc(3)    
       null(),&                               ! parent     
       null())                                ! children(:)

  type maths_octree
     type(octbox),  pointer  :: tree_root
     integer(ip)             :: divmax
     integer(ip)             :: dim
     integer(ip)             :: nleaves
     integer(ip)             :: nfilled
     real(rp)                :: comin(3)
     real(rp)                :: comax(3)
   contains
     procedure,         pass :: init           ! Initialize all
     procedure,         pass :: fill           ! Allocate 
     procedure,         pass :: deallo         ! Deallocate
     procedure,         pass :: find           ! Find a bin
     procedure,         pass :: inbox          ! If a point is inside boundaing box
     procedure,         pass :: mesh           ! Create a mesh
     procedure,         pass :: mesh_dim       ! Get a mehs dimensions
     procedure,         pass :: id             ! Get the bin ID of a point
     procedure,         pass :: weight         ! Get the weight of a bin
     procedure,         pass :: centroid       ! Centroid coordinates of the bins
  end type maths_octree

  private

  public :: maths_octree

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Octree
  !> @details Fill in a octree from a coordinate field
  !> 
  !-----------------------------------------------------------------------

  subroutine init(octree)

    class(maths_octree), intent(inout) :: octree

    nullify(octree % tree_root)
    octree % divmax  = 0
    octree % nleaves = 0
    octree % nfilled = 0
    octree % dim     = 0
    octree % comin   = huge(1.0_rp)*0.1_rp
    octree % comax   = huge(1.0_rp)*0.1_rp
    
  end subroutine init

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Deallocate
  !> @details Deallocate octree structure
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo(octree,MEMORY_COUNTER)

    class(maths_octree),            intent(inout) :: octree
    integer(8),           optional, intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                                    :: memor_loc(2)
    type(octbox), pointer                         :: current_o
    logical(lg)                                   :: conti

    if( present(MEMORY_COUNTER) ) then
       memor_loc = MEMORY_COUNTER
    else
       memor_loc = 0_8
    end if
    
    if( associated(octree % tree_root) ) then

       current_o => octree % tree_root
       conti     =  .true.

       do while ( conti )
          !
          ! First go to deepest level in first branch
          !
          do while( current_o % whoiam == 0 )
             current_o => current_o % children(1)
          end do
          !
          ! Deallocate list of elements
          !
          if( current_o % whoiam > 0 ) then
             call memory_deallo(memor_loc,'CURRENT_O % NODES','deallo',current_o % nodes)
          end if

          if( current_o % childid < octree % divmax .and. current_o % childid /= 0 ) then
             !
             ! I'm not the last child neither the Padrino
             !
             current_o => current_o % parent % children(current_o % childid+1)

          else if( current_o % childid == octree % divmax ) then
             !
             ! I'm the last child
             !
             current_o => current_o % parent 
             deallocate(current_o % children)
             current_o % whoiam = -1

          else if( current_o % id == 0 ) then
             !
             ! I'm the Padrino: end of deallocation
             !
             deallocate(current_o)
             conti = .false.

          end if

       end do

    end if

    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor_loc
    
  end subroutine deallo

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Octree
  !> @details Fill in a octree from a coordinate field
  !> 
  !-----------------------------------------------------------------------

  subroutine fill(octree,coord,limit,MEMORY_COUNTER,OFFSET,COORD_MIN,COORD_MAX)

    class(maths_octree),            intent(inout) :: octree
    real(rp),              pointer, intent(in)    :: coord(:,:)
    integer(ip),                    intent(in)    :: limit
    integer(8),  optional,          intent(inout) :: MEMORY_COUNTER(2)
    real(rp),    optional,          intent(in)    :: OFFSET
    real(rp),    optional,          intent(in)    :: COORD_MIN(:)
    real(rp),    optional,          intent(in)    :: COORD_MAX(:)
    type(octbox),          pointer                :: tree_root
    type(octbox),          pointer                :: current_o
    type(octbox),          pointer                :: old_pointer
    type(octbox),          pointer                :: tm1_pointer
    type(octbox),          pointer                :: tm2_pointer    
    integer(ip)                                   :: ii,jj,kk,idime,ndime
    integer(ip)                                   :: npoin,ipoin,kpoin
    real(rp)                                      :: delta(3)
    integer(ip)                                   :: xx(3)
    integer(8)                                    :: memor_loc(2)
    integer(ip)                                   :: counter,divmax
    logical(lg)                                   :: conti
    real(rp)                                      :: offset_loc

    if( present(OFFSET) ) then
       offset_loc = OFFSET
    else
       offset_loc = 1.0e-3_rp   
    end if
    
    allocate(octree % tree_root)
    octree % tree_root = octbox_init
    ndime              = size(coord,1)
    npoin              = size(coord,2)
    divmax             = 2**ndime
    octree % dim       = ndime
    octree % divmax    = divmax

    if( present(COORD_MIN) .and. present(COORD_MAX) ) then
       octree % comin(1:octree % dim) = COORD_MIN(1:octree % dim)
       octree % comax(1:octree % dim) = COORD_MAX(1:octree % dim)
    else
       do idime = 1,octree % dim
          octree % comin(idime) = minval(coord(idime,:))
          octree % comax(idime) = maxval(coord(idime,:))
       end do
    end if
    delta          = ( octree % comax - octree % comin ) * offset_loc 
    octree % comin = octree % comin - delta - epsil
    octree % comax = octree % comax + delta + epsil
    !
    ! Root
    !
    call memory_alloca(memor_loc,'OCTREE % TREE_ROOT % NODES','fill',octree % tree_root % nodes,npoin)
    tree_root            => octree % tree_root
    current_o            => tree_root
    current_o % npoinbox =  0
    current_o % id       =  0
    current_o % level    =  0
    current_o % childid  =  0
    current_o % whoiam   =  0
    current_o % npoinbox =  npoin

    do ipoin = 1,npoin
       current_o % nodes(ipoin) = ipoin
    end do
    current_o % minc(1:ndime) = octree % comin(1:ndime)
    current_o % maxc(1:ndime) = octree % comax(1:ndime)

    nullify( current_o % parent )

    counter = 0
    conti   = .true.
    
    do while( conti )
       !
       ! If maximum number of points inside current box is exceeded, subdivide
       !
       if( current_o % npoinbox > limit ) then
          allocate( current_o % children(divmax) )
          current_o % children(1:divmax) = octbox_init
          !
          ! Give birth to my DIVMAX children
          !
          do ii = 1,divmax
             counter                             =  counter+1
             current_o % children(ii) % id       =  counter
             current_o % children(ii) % childid  =  ii
             current_o % children(ii) % level    =  current_o % level + 1 
             current_o % children(ii) % whoiam   =  0  
             current_o % children(ii) % npoinbox =  0
             current_o % children(ii) % idfilled =  0
             current_o % children(ii) % parent   => current_o

             call memory_alloca(memor_loc,'CURRENT_O % NODES','fill',current_o % children(ii) % nodes,current_o % npoinbox)
          end do
          !
          ! Compute the coordinates of my children
          !
          do ii = 1,ndime
             current_o % children(1) % minc(ii) = current_o % minc(ii)
             current_o % children(1) % maxc(ii) = (current_o % maxc(ii) + current_o % minc(ii))*0.5_rp
          end do
          current_o % children(2) % minc(1) = (current_o % maxc(1) + current_o % minc(1))*0.5_rp
          current_o % children(2) % minc(2) = current_o % children(1) % minc(2)
          current_o % children(2) % maxc(1) = current_o % maxc(1)     
          current_o % children(2) % maxc(2) = current_o % children(1) % maxc(2)

          current_o % children(3) % minc(1) = current_o % children(1) % minc(1)
          current_o % children(3) % minc(2) = (current_o % minc(2) + current_o % maxc(2))*0.5_rp
          current_o % children(3) % maxc(1) = current_o % children(1) % maxc(1)
          current_o % children(3) % maxc(2) = current_o % maxc(2)

          current_o % children(4) % minc(1) = current_o % children(2) % minc(1)
          current_o % children(4) % minc(2) = current_o % children(3) % minc(2)
          current_o % children(4) % maxc(1) = current_o % children(2) % maxc(1)
          current_o % children(4) % maxc(2) = current_o % children(3) % maxc(2)

          if( ndime == 3 ) then
             current_o % children(2) % minc(3) = current_o % children(1) % minc(3)
             current_o % children(2) % maxc(3) = current_o % children(1) % maxc(3)
             current_o % children(3) % minc(3) = current_o % children(1) % minc(3)
             current_o % children(3) % maxc(3) = current_o % children(1) % maxc(3)
             current_o % children(4) % minc(3) = current_o % children(1) % minc(3)
             current_o % children(4) % maxc(3) = current_o % children(1) % maxc(3)

             current_o % children(5) % minc(1) = current_o % children(1) % minc(1)
             current_o % children(5) % minc(2) = current_o % children(1) % minc(2)
             current_o % children(5) % minc(3) = (current_o % minc(3) + current_o % maxc(3))*0.5_rp
             current_o % children(5) % maxc(1) = current_o % children(1) % maxc(1)
             current_o % children(5) % maxc(2) = current_o % children(1) % maxc(2)
             current_o % children(5) % maxc(3) = current_o % maxc(3)

             current_o % children(6) % minc(1) = current_o % children(2) % minc(1)
             current_o % children(6) % minc(2) = current_o % children(1) % minc(2)
             current_o % children(6) % minc(3) = current_o % children(5) % minc(3)
             current_o % children(6) % maxc(1) = current_o % children(2) % maxc(1)     
             current_o % children(6) % maxc(2) = current_o % children(1) % maxc(2)
             current_o % children(6) % maxc(3) = current_o % children(5) % maxc(3)

             current_o % children(7) % minc(1) = current_o % children(1) % minc(1)
             current_o % children(7) % minc(2) = current_o % children(3) % minc(2)
             current_o % children(7) % minc(3) = current_o % children(5) % minc(3)
             current_o % children(7) % maxc(1) = current_o % children(1) % maxc(1)
             current_o % children(7) % maxc(2) = current_o % children(3) % maxc(2)
             current_o % children(7) % maxc(3) = current_o % children(5) % maxc(3)

             current_o % children(8) % minc(1) = current_o % children(2) % minc(1)
             current_o % children(8) % minc(2) = current_o % children(3) % minc(2)
             current_o % children(8) % minc(3) = current_o % children(5) % minc(3)
             current_o % children(8) % maxc(1) = current_o % children(2) % maxc(1)
             current_o % children(8) % maxc(2) = current_o % children(3) % maxc(2)
             current_o % children(8) % maxc(3) = current_o % children(5) % maxc(3)
          end if
          !
          ! Offer my nodes to my children
          !
          if( ndime == 2 ) then
             do ii = 1,4
                do ipoin = 1,current_o % npoinbox
                   kpoin = current_o % nodes(ipoin)
                   if(    coord(1,kpoin) >= current_o % children(ii) % minc(1) .and. &
                        & coord(2,kpoin) >= current_o % children(ii) % minc(2) .and. &
                        & coord(1,kpoin) <  current_o % children(ii) % maxc(1) .and. &
                        & coord(2,kpoin) <  current_o % children(ii) % maxc(2) ) then
                      current_o % children(ii) % npoinbox = current_o % children(ii) % npoinbox + 1
                      current_o % children(ii) % nodes(current_o % children(ii) % npoinbox) = kpoin
                   end if
                end do
             end do
          else
             do ii = 1,8
                do ipoin = 1,current_o % npoinbox
                   kpoin = current_o % nodes(ipoin)
                   if(                coord(1,kpoin) >= current_o % children(ii) % minc(1) ) then
                      if(             coord(2,kpoin) >= current_o % children(ii) % minc(2) ) then
                         if(          coord(3,kpoin) >= current_o % children(ii) % minc(3) ) then 
                            if(       coord(1,kpoin) <  current_o % children(ii) % maxc(1) ) then 
                               if(    coord(2,kpoin) <  current_o % children(ii) % maxc(2) ) then 
                                  if( coord(3,kpoin) <  current_o % children(ii) % maxc(3) ) then
                                     current_o % children(ii) % npoinbox = current_o % children(ii) % npoinbox + 1
                                     current_o % children(ii) % nodes(current_o % children(ii) % npoinbox) = kpoin
                                  end if
                               end if
                            end if
                         end if
                      end if
                   end if
                end do
             end do
          end if

          call memory_deallo(memor_loc,'CURRENT_O % NODES','fill',current_o % nodes)

          current_o % whoiam   =  0
          current_o % npoinbox =  0
          current_o            => current_o % children(1)

       else if(current_o % id == 0 .and. current_o % npoinbox <= limit ) then
          !
          ! If the Padrino has too few elements
          !
          conti              =  .false.
          octree % nleaves   =  octree % nleaves + 1
          current_o % whoiam =  octree % nleaves
          if( current_o % npoinbox > 0 ) then
             octree % nfilled     =  octree % nfilled + 1
             current_o % idfilled =  octree % nfilled
          end if
          current_o          => old_pointer
          
       else 
          !
          ! if limit of points inside box is not exceeded, assign nodes
          !
          !call memory_deallo(memor_loc,'CURRENT_O % NODES','fill',current_o % nodes)
          octree % nleaves   = octree % nleaves + 1
          current_o % whoiam = octree % nleaves
          if( current_o % npoinbox > 0 ) then
             octree % nfilled     =  octree % nfilled + 1
             current_o % idfilled =  octree % nfilled
          end if
           
          if( current_o % childid < divmax .and. current_o % id /= 0 ) then
             !
             ! Go to next children
             !
             tm1_pointer => current_o
             tm2_pointer => tm1_pointer % parent % children(tm1_pointer%childid+1)
             current_o   => tm2_pointer
             goto 10

          else if(current_o % childid == 0 ) then  
             !
             ! Padrino
             !
             goto 10

          else if(current_o % childid == divmax ) then
             !
             ! Last children
             !
             noparent: do while( current_o % id > 0 )
                if(current_o % parent % id == 0) then
                   conti = .false.
                   exit noparent
                else
                   if( current_o % parent % childid /= divmax ) then 
                      tm1_pointer => current_o
                      tm2_pointer => tm1_pointer % parent % parent % children(tm1_pointer%parent%childid+1)
                      current_o   => tm2_pointer
                      exit
                   else 
                      current_o   => current_o % parent
                   end if
                end if
             end do noparent

          else 
             !
             ! Wrong child ID
             !
             stop
             
          end if

       end if

10     continue
       
       old_pointer => current_o

    end do
    
 end subroutine fill
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Location in a octree
  !> @details Find the location of a point in a octree
  !> 
  !-----------------------------------------------------------------------

  pure subroutine find(octree,coord,xx) 

    class(maths_octree), intent(in)  :: octree
    real(rp),            intent(in)  :: coord(:) !< Coordinate of the test point
    integer(ip),         intent(out) :: xx(3)    !< Box in x direction

  end subroutine find
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Location in a octree
  !> @details Find the location of a point in a octree
  !> 
  !-----------------------------------------------------------------------

  pure logical(lg) function inbox(octree,coord) 

   class(maths_octree), intent(in)  :: octree
   real(rp),            intent(in)  :: coord(:) !< Coordinate of the test point
   end function inbox 

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-09
  !> @brief   Get mesh dimensions
  !> @details Get the mesh dimensions
  !> 
  !-----------------------------------------------------------------------
  
  subroutine mesh_dim(octree,ndime,mnode,nelem,npoin,ENABLE_EMPTY)

    class(maths_octree),                    intent(in)    :: octree       !< Octree
    integer(ip),                            intent(out)   :: ndime        !< Space dimension
    integer(ip),                            intent(out)   :: mnode        !< Max number of nodes per element
    integer(ip),                            intent(out)   :: nelem        !< Number of elements
    integer(ip),                            intent(out)   :: npoin        !< Number of nodes
    logical(lg),         optional,          intent(in)    :: ENABLE_EMPTY !< If empty bins should be considered
    logical(lg)                                           :: if_empty
    
    if( present(ENABLE_EMPTY) ) then
       if_empty = ENABLE_EMPTY
    else
       if_empty = .true.
    end if
    
    ndime     =  octree % dim
    if( if_empty ) then
       nelem     = octree % nleaves
    else
       nelem     = octree % nfilled
    end if
    mnode     =  (ndime-1)*4
    npoin     =  nelem * mnode

  end subroutine mesh_dim
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-09
  !> @brief   Create a mesh
  !> @details Create a mesh from an octree
  !> 
  !-----------------------------------------------------------------------
  
  subroutine mesh(octree,ndime,mnode,nelem,npoin,lnods,ltype,coord,MEMORY_COUNTER,OFFSET_IELEM,OFFSET_IPOIN,ENABLE_EMPTY,CENTROID)

    class(maths_octree),                    intent(in)    :: octree
    integer(ip),                            intent(out)   :: ndime
    integer(ip),                            intent(out)   :: mnode
    integer(ip),                            intent(out)   :: nelem
    integer(ip),                            intent(out)   :: npoin
    integer(ip),                   pointer, intent(inout) :: lnods(:,:)
    integer(ip),                   pointer, intent(inout) :: ltype(:)
    real(rp),                      pointer, intent(inout) :: coord(:,:)
    integer(8),          optional,          intent(inout) :: MEMORY_COUNTER(2)
    integer(ip),         optional,          intent(in)    :: OFFSET_IELEM
    integer(ip),         optional,          intent(in)    :: OFFSET_IPOIN
    logical(lg),         optional,          intent(in)    :: ENABLE_EMPTY      !< If empty bins should be considered
    real(rp),            optional, pointer, intent(inout) :: CENTROID(:,:)
    type(octbox),                  pointer                :: current_o
    integer(ip)                                           :: ipoin,ielem,idime
    integer(ip)                                           :: inode,divmax
    logical(lg)                                           :: conti
    integer(8)                                            :: memor_loc(2)
    integer(ip)                                           :: offset_ielem_loc
    integer(ip)                                           :: offset_ipoin_loc
    logical(lg)                                           :: if_empty
    
    if( present(MEMORY_COUNTER) ) then
       memor_loc = MEMORY_COUNTER
    else
       memor_loc = 0_8
    end if
    if( present(OFFSET_IELEM) ) then
       offset_ielem_loc = OFFSET_IELEM
    else
       offset_ielem_loc = 0
    end if
    if( present(OFFSET_IPOIN) ) then
       offset_ipoin_loc = OFFSET_IPOIN
    else
       offset_ipoin_loc = 0
    end if
    if( present(ENABLE_EMPTY) ) then
       if_empty = ENABLE_EMPTY
    else
       if_empty = .true.
    end if

    current_o => octree % tree_root
    divmax    =  octree % divmax
    call octree % mesh_dim(ndime,mnode,nelem,npoin,ENABLE_EMPTY)
    !
    ! Allocate mesh
    !
    if( .not. associated(lnods) ) call memory_alloca(memor_loc,'LNODS','mesh',lnods,mnode,nelem)
    if( .not. associated(ltype) ) call memory_alloca(memor_loc,'LTYPE','mesh',ltype,nelem)
    if( .not. associated(coord) ) call memory_alloca(memor_loc,'COORD','mesh',coord,ndime,npoin)
    if( present(CENTROID) ) then
       if( .not. associated(CENTROID) ) call memory_alloca(memor_loc,'CENTROID','mesh',CENTROID,ndime,nelem)
    end if
    !
    ! Mesh arrays
    !    
    ipoin =  1 + offset_ipoin_loc
    conti =  .true.
    
    loop_conti2: do while( conti )
       !
       ! First go to deepest level in first branch
       !
       do while( current_o % whoiam == 0 )
          current_o => current_o % children(1)
       end do

       if( if_empty .and. current_o % whoiam /= 0 ) then
          ielem = current_o % whoiam + offset_ielem_loc
       else if( (.not. if_empty) .and. current_o % idfilled /= 0 ) then
          ielem = current_o % idfilled + offset_ielem_loc
       else
          ielem = 0
       end if
       
       if( ielem /= 0 ) then
          !
          ! Current bin is a leaf
          !
          if( ndime == 2 ) then
             coord(1:2,ipoin)   = (/ current_o % minc(1),current_o % minc(2) /)
             coord(1:2,ipoin+1) = (/ current_o % maxc(1),current_o % minc(2) /)
             coord(1:2,ipoin+2) = (/ current_o % maxc(1),current_o % maxc(2) /)
             coord(1:2,ipoin+3) = (/ current_o % minc(1),current_o % maxc(2) /)
             lnods(1,ielem)     = ipoin 
             lnods(2,ielem)     = ipoin + 1
             lnods(3,ielem)     = ipoin + 2
             lnods(4,ielem)     = ipoin + 3
             ltype(ielem)       = QUA04
             if( present(CENTROID) ) then
                do idime = 1,2
                   CENTROID(idime,ielem) = sum(coord(idime,ipoin:ipoin+3))/4.0_rp
                end do
             end if
             ipoin              = ipoin + 4
          else
             coord(1:3,ipoin)   = (/ current_o % minc(1),current_o % minc(2),current_o % minc(3) /)
             coord(1:3,ipoin+1) = (/ current_o % maxc(1),current_o % minc(2),current_o % minc(3) /)
             coord(1:3,ipoin+2) = (/ current_o % maxc(1),current_o % maxc(2),current_o % minc(3) /)
             coord(1:3,ipoin+3) = (/ current_o % minc(1),current_o % maxc(2),current_o % minc(3) /)
             coord(1:3,ipoin+4) = (/ current_o % minc(1),current_o % minc(2),current_o % maxc(3) /)
             coord(1:3,ipoin+5) = (/ current_o % maxc(1),current_o % minc(2),current_o % maxc(3) /)
             coord(1:3,ipoin+6) = (/ current_o % maxc(1),current_o % maxc(2),current_o % maxc(3) /)
             coord(1:3,ipoin+7) = (/ current_o % minc(1),current_o % maxc(2),current_o % maxc(3) /)
             lnods(1,ielem)     = ipoin 
             lnods(2,ielem)     = ipoin + 1
             lnods(3,ielem)     = ipoin + 2
             lnods(4,ielem)     = ipoin + 3
             lnods(5,ielem)     = ipoin + 4
             lnods(6,ielem)     = ipoin + 5
             lnods(7,ielem)     = ipoin + 6
             lnods(8,ielem)     = ipoin + 7
             ltype(ielem)       = HEX08
             if( present(CENTROID) ) then
                do idime = 1,3
                   CENTROID(idime,ielem) = sum(coord(idime,ipoin:ipoin+7))/8.0_rp
                end do
             end if
             ipoin              = ipoin + 8             
          end if
       end if

       if(current_o % childid < divmax .and. current_o % childid /=0 ) then
          !
          ! I'm not the last child neither the Padrino
          !
          current_o => current_o % parent % children(current_o % childid+1)

       else if( current_o % childid == divmax ) then
          !
          ! I'm the last child of this generation: postprocess 
          !
          do while(current_o % id > 0 )
             if(current_o % parent % id == 0) then
                conti = .false. 
                exit loop_conti2
             else
                if(current_o % parent % childid /=divmax) then
                   current_o => current_o % parent % parent % children(current_o % parent % childid+1)
                   exit
                else 
                   current_o => current_o % parent
                end if
             end if
          end do

       else if( current_o % id == 0 ) then
          !
          ! I'm the Padrino
          !
          conti = .false.
          
       end if

    end do loop_conti2
    
    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor_loc

  end subroutine mesh

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-09
  !> @brief   Return the centroid
  !> @details Compute the centroids of the leaves
  !> 
  !-----------------------------------------------------------------------
  
  subroutine centroid(octree,coorc,MEMORY_COUNTER,OFFSET_IELEM,ENABLE_EMPTY)

    class(maths_octree),                    intent(in)    :: octree
    real(rp),                      pointer, intent(inout) :: coorc(:,:)
    integer(8),          optional,          intent(inout) :: MEMORY_COUNTER(2)
    integer(ip),         optional,          intent(in)    :: OFFSET_IELEM
    logical(lg),         optional,          intent(in)    :: ENABLE_EMPTY      !< If empty bins should be considered
    type(octbox),                  pointer                :: current_o
    integer(ip)                                           :: ielem,divmax,ndime
    integer(ip)                                           :: nelem,mnode,idime
    integer(ip)                                           :: npoin
    logical(lg)                                           :: conti
    integer(8)                                            :: memor_loc(2)
    integer(ip)                                           :: offset_ielem_loc
    real(rp)                                              :: xx(3,8),rnode
    logical(lg)                                           :: if_empty
    
    if( present(MEMORY_COUNTER) ) then
       memor_loc = MEMORY_COUNTER
    else
       memor_loc = 0_8
    end if
    if( present(OFFSET_IELEM) ) then
       offset_ielem_loc = OFFSET_IELEM
    else
       offset_ielem_loc = 0
    end if
    if( present(ENABLE_EMPTY) ) then
       if_empty = ENABLE_EMPTY
    else
       if_empty = .true.
    end if

    call octree % mesh_dim(ndime,mnode,nelem,npoin,ENABLE_EMPTY)    
    divmax    =  octree % divmax
    current_o => octree % tree_root
    rnode     =  1.0_rp / real(mnode,rp)
    !
    ! Allocate mesh
    !
    if( .not. associated(coorc) ) call memory_alloca(memor_loc,'COORC','mesh',coorc,ndime,nelem)
    !
    ! Mesh arrays
    !    
    conti =  .true.
    
    loop_conti2: do while( conti )
       !
       ! First go to deepest level in first branch
       !
       do while( current_o % whoiam == 0 )
          current_o => current_o % children(1)
       end do

       if( if_empty .and. current_o % whoiam /= 0 ) then
          ielem = current_o % whoiam + offset_ielem_loc
       else if( (.not. if_empty) .and. current_o % idfilled /= 0 ) then
          ielem = current_o % idfilled + offset_ielem_loc
       else
          ielem = 0
       end if

       if( ielem /= 0 ) then
          !
          ! Current bin is a leaf
          !
          if( ndime == 2 ) then
             xx(1:2,1) = (/ current_o % minc(1),current_o % minc(2) /)
             xx(1:2,2) = (/ current_o % maxc(1),current_o % minc(2) /)
             xx(1:2,3) = (/ current_o % maxc(1),current_o % maxc(2) /)
             xx(1:2,4) = (/ current_o % minc(1),current_o % maxc(2) /)
          else
             xx(1:3,1) = (/ current_o % minc(1),current_o % minc(2),current_o % minc(3) /)
             xx(1:3,2) = (/ current_o % maxc(1),current_o % minc(2),current_o % minc(3) /)
             xx(1:3,3) = (/ current_o % maxc(1),current_o % maxc(2),current_o % minc(3) /)
             xx(1:3,4) = (/ current_o % minc(1),current_o % maxc(2),current_o % minc(3) /)
             xx(1:3,5) = (/ current_o % minc(1),current_o % minc(2),current_o % maxc(3) /)
             xx(1:3,6) = (/ current_o % maxc(1),current_o % minc(2),current_o % maxc(3) /)
             xx(1:3,7) = (/ current_o % maxc(1),current_o % maxc(2),current_o % maxc(3) /)
             xx(1:3,8) = (/ current_o % minc(1),current_o % maxc(2),current_o % maxc(3) /)
          end if
          do idime = 1,ndime
             coorc(idime,ielem) = sum(xx(idime,1:8))  
          end do
          coorc(:,ielem) = coorc(:,ielem) * rnode

       end if

       if(current_o % childid < divmax .and. current_o % childid /=0 ) then
          !
          ! I'm not the last child neither the Padrino
          !
          current_o => current_o % parent % children(current_o % childid+1)

       else if( current_o % childid == divmax ) then
          !
          ! I'm the last child of this generation: postprocess 
          !
          do while(current_o % id > 0 )
             if(current_o % parent % id == 0) then
                conti = .false. 
                exit loop_conti2
             else
                if(current_o % parent % childid /=divmax) then
                   current_o => current_o % parent % parent % children(current_o % parent % childid+1)
                   exit
                else 
                   current_o => current_o % parent
                end if
             end if
          end do

       else if( current_o % id == 0 ) then
          !
          ! I'm the Padrino
          !
          conti = .false.
          
       end if

    end do loop_conti2
    
    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor_loc

  end subroutine centroid

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-10
  !> @brief   Get the global ID of the bin
  !> @details Given a point, get the global ID of the bin it
  !>          is located in
  !> 
  !-----------------------------------------------------------------------
  
  integer(ip) function id(octree,coord) 

    class(maths_octree),      intent(in) :: octree
    real(rp),                 intent(in) :: coord(:) !< Coordinate of the test point
    type(octbox),     pointer            :: current_o
    integer(ip)                          :: ichild

    current_o  => octree % tree_root

    if( octree % dim == 3 ) then

       do while( current_o % whoiam == 0 )    
          childloop8: do ichild = 1,8           
             if(    coord(1) >= current_o % children(ichild) % minc(1) .and. &
                  & coord(1) <= current_o % children(ichild) % maxc(1) .and. &
                  & coord(2) >= current_o % children(ichild) % minc(2) .and. &
                  & coord(2) <= current_o % children(ichild) % maxc(2) .and. &
                  & coord(3) >= current_o % children(ichild) % minc(3) .and. &
                  & coord(3) <= current_o % children(ichild) % maxc(3) ) then
                current_o => current_o % children(ichild)
                exit childloop8
             end if
          end do childloop8
       end do

    else if( octree % dim == 2 ) then

       do while( current_o % whoiam == 0 )  
          childloop4: do ichild = 1,4
             if(    coord(1) >= current_o % children(ichild) % minc(1) .and. &
                  & coord(1) <= current_o % children(ichild) % maxc(1) .and. &
                  & coord(2) >= current_o % children(ichild) % minc(2) .and. &
                  & coord(2) <= current_o % children(ichild) % maxc(2)  ) then
                current_o => current_o % children(ichild)
                exit childloop4
             end if
          end do childloop4
       end do

    end if

    id = current_o % whoiam
    
  end function id
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-10
  !> @brief   Get the global ID of the bin
  !> @details Given a point, get the global ID of the bin it
  !>          is located in
  !> 
  !-----------------------------------------------------------------------
  
  integer(ip) function weight(octree,coord) 

    class(maths_octree),      intent(in) :: octree
    real(rp),                 intent(in) :: coord(:) !< Coordinate of the test point
    type(octbox),     pointer            :: current_o
    integer(ip)                          :: ichild
    
    current_o  => octree % tree_root

    if( octree % dim == 3 ) then

       do while( current_o % whoiam == 0 )    
          childloop8: do ichild = 1,8           
             if(    coord(1) >= current_o % children(ichild) % minc(1) .and. &
                  & coord(1) <= current_o % children(ichild) % maxc(1) .and. &
                  & coord(2) >= current_o % children(ichild) % minc(2) .and. &
                  & coord(2) <= current_o % children(ichild) % maxc(2) .and. &
                  & coord(3) >= current_o % children(ichild) % minc(3) .and. &
                  & coord(3) <= current_o % children(ichild) % maxc(3) ) then
                current_o => current_o % children(ichild)
                exit childloop8
             end if
          end do childloop8
       end do

    else if( octree % dim == 2 ) then

       do while( current_o % whoiam == 0 )  
          childloop4: do ichild = 1,4
             if(    coord(1) >= current_o % children(ichild) % minc(1) .and. &
                  & coord(1) <= current_o % children(ichild) % maxc(1) .and. &
                  & coord(2) >= current_o % children(ichild) % minc(2) .and. &
                  & coord(2) <= current_o % children(ichild) % maxc(2)  ) then
                current_o => current_o % children(ichild)
                exit childloop4
             end if
          end do childloop4
       end do

    end if

    weight = current_o % npoinbox
    
  end function weight
  
  integer(ip) function weight_id(octree,id) 

    class(maths_octree),      intent(in) :: octree
    integer(ip),              intent(in) :: id !< Coordinate of the test point
    type(octbox),     pointer            :: current_o
    integer(ip)                          :: ichild
    
  end function weight_id

  !integer(ip) function weight(octree,that)    
  !  class(maths_octree), intent(in) :: octree
  !  class(*),            intent(in) :: that
  !  select type(that)       
  !  type is (integer)
  !     weight = weight_id(octree,that)
  !  type is (real)
  !     weight = weight_x (octree,that)
  !  end select
  !end function weight
  
end module def_maths_octree
!> @}
