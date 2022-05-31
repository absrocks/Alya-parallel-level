!-----------------------------------------------------------------------
!
!> @defgroup Maths_Toolbox
!> Toolbox for octbins
!> @{
!> @name    ToolBox for octbins
!> @file    def_maths.f90
!> @author  Guillaume Houzeaux
!> @brief   Variables
!> @details Variables fot mod_maths.f90
!
!-----------------------------------------------------------------------

module def_maths_octbin

  use def_kintyp_basic, only : ip,rp,lg
  use mod_memory_basic, only : memory_alloca
  use mod_memory_basic, only : memory_deallo
  use mod_memory_basic, only : memory_size
  use def_elmtyp,       only : HEX08
  use def_elmtyp,       only : QUA04
  use def_maths_bin
  use def_maths_octree

  real(rp), parameter :: epsil = epsilon(1.0_rp)

  type test
     type(maths_octree) :: octree
  end type test
  
  type maths_octbin
     integer(ip)                     :: num_octree
     type(maths_bin)                 :: bin
     type(maths_octree), allocatable :: octree(:,:,:)
   contains
     procedure,         pass :: init           ! Initialize all
     procedure,         pass :: fill           ! Fill the octbin
     procedure,         pass :: mesh           ! Mesh
     procedure,         pass :: mesh_dim       ! Mesh dimensions
     procedure,         pass :: deallo         ! Deallocate
     procedure,         pass :: centroid       ! Get centroid
  end type maths_octbin

  private

  public :: maths_octbin

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Octbin
  !> @details Fill in a octbin from a coordinate field
  !> 
  !-----------------------------------------------------------------------

  subroutine init(octbin)

    class(maths_octbin), intent(inout) :: octbin

    call octbin % bin % init()
    octbin % num_octree = 0
    
  end subroutine init

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Deallocate
  !> @details Deallocate and octebin
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo(octbin)

    class(maths_octbin), intent(inout) :: octbin
    integer(ip)                        :: ii,jj,kk
    
    do kk = 1,octbin % bin % boxip(3)
       do jj = 1,octbin % bin % boxip(2)
          do ii = 1,octbin % bin % boxip(1)
             call octbin % octree(ii,jj,kk) % deallo()
          end do
       end do
    end do
    deallocate(octbin % octree)
    call octbin % bin % deallo()
    
  end subroutine deallo

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Octbin
  !> @details Fill in a octbin from a coordinate field
  !> 
  !-----------------------------------------------------------------------

  subroutine fill(octbin,coord,boxes,limit,MEMORY_COUNTER,OFFSET,ENABLE_GAP)

    class(maths_octbin),            intent(inout) :: octbin
    real(rp),              pointer, intent(in)    :: coord(:,:)
    integer(ip),                    intent(in)    :: boxes(:)
    integer(ip),                    intent(in)    :: limit
    integer(8),  optional,          intent(inout) :: MEMORY_COUNTER(2)
    real(rp),    optional,          intent(in)    :: OFFSET
    logical(lg), optional,          intent(in)    :: ENABLE_GAP
    integer(ip)                                   :: ii,jj,kk,idime
    integer(ip)                                   :: npoin,ll,ndime,nsize
    integer(ip)                                   :: ipoin,kpoin
    real(rp)                                      :: delta(3)
    real(rp)                                      :: comin(3)
    real(rp)                                      :: comax(3)
    real(rp),              pointer                :: xx(:,:)
    integer(8)                                    :: memor_loc(2)
    logical(lg)                                   :: if_gap

    if( present(MEMORY_COUNTER) ) then
       memor_loc = MEMORY_COUNTER
    else
       memor_loc = 0_8
    end if
    if( present(ENABLE_GAP) ) then
       if_gap = ENABLE_GAP
    else
       if_gap = .true.
    end if
    
    nullify(xx)

    call octbin % bin % init()
    call octbin % bin % fill(coord,boxes,MEMORY_COUNTER,OFFSET)
    ndime = octbin % bin % dim

    allocate(octbin % octree(octbin % bin % boxip(1),octbin % bin % boxip(2),octbin % bin % boxip(3)))

    do kk = 1,octbin % bin % boxip(3)
       do jj = 1,octbin % bin % boxip(2)
          do ii = 1,octbin % bin % boxip(1)
             call octbin % octree(ii,jj,kk) % init()
             if( associated(octbin % bin % list(ii,jj,kk) % l) ) then
                !
                ! Construct the octree of this bin
                !
                npoin = size( octbin % bin % list(ii,jj,kk) % l)
                allocate(xx(ndime,npoin))
                do ll = 1,npoin
                   ipoin          = octbin % bin % list(ii,jj,kk) % l(ll)
                   xx(1:ndime,ll) = coord(1:ndime,ipoin)
                end do
                if( if_gap ) then
                   call octbin % octree(ii,jj,kk) % fill(xx,limit,MEMORY_COUNTER,OFFSET)
                else
                   call octbin % bin % coord((/ii,jj,kk/),comin,comax)
                   call octbin % octree(ii,jj,kk) % fill(xx,limit,MEMORY_COUNTER,OFFSET,COORD_MIN=comin,COORD_MAX=comax)                   
                end if
                octbin % num_octree = octbin % num_octree + 1
                deallocate(xx)
             end if
          end do
       end do
    end do

    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor_loc

  end subroutine fill

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-09
  !> @brief   Mesh dimension
  !> @details Find mesh dimension
  !> 
  !-----------------------------------------------------------------------
  
  subroutine mesh_dim(octbin,ndime,mnode,nelem,npoin,ENABLE_EMPTY)

    class(maths_octbin),                    intent(inout) :: octbin
    integer(ip),                            intent(out)   :: ndime
    integer(ip),                            intent(out)   :: mnode
    integer(ip),                            intent(out)   :: nelem
    integer(ip),                            intent(out)   :: npoin
    logical(lg),         optional,          intent(in)    :: ENABLE_EMPTY      !< If empty bins should be considered
    integer(ip)                                           :: ielem,ipoin 
    integer(ip)                                           :: ii,jj,kk
    logical(lg)                                           :: if_empty
    
    if( present(ENABLE_EMPTY) ) then
       if_empty = ENABLE_EMPTY
    else
       if_empty = .true.
    end if
    
    nelem = 0
    npoin = 0
    mnode = (octbin % bin % dim-1)*4
    ndime = octbin % bin % dim
    
    do kk = 1,octbin % bin % boxip(3)
       do jj = 1,octbin % bin % boxip(2)
          do ii = 1,octbin % bin % boxip(1)
             if( associated(octbin % bin % list(ii,jj,kk) % l) ) then
                call octbin % octree(ii,jj,kk) % mesh_dim(ndime,mnode,ielem,ipoin,ENABLE_EMPTY)
                nelem = nelem + ielem
                npoin = npoin + ipoin
             else if( if_empty ) then
                nelem = nelem + 1
                npoin = npoin + mnode
             end if
          end do
       end do
    end do
    
  end subroutine mesh_dim
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-09
  !> @brief   Create a mesh
  !> @details Create a mesh from an octbin
  !> 
  !-----------------------------------------------------------------------
  
  subroutine mesh(octbin,ndime,mnode,nelem,npoin,lnods,ltype,coord,MEMORY_COUNTER,ENABLE_EMPTY,CENTROID)

    class(maths_octbin),                    intent(inout) :: octbin
    integer(ip),                            intent(out)   :: ndime
    integer(ip),                            intent(out)   :: mnode
    integer(ip),                            intent(out)   :: nelem
    integer(ip),                            intent(out)   :: npoin
    integer(ip),                   pointer, intent(inout) :: lnods(:,:)
    integer(ip),                   pointer, intent(inout) :: ltype(:)
    real(rp),                      pointer, intent(inout) :: coord(:,:)
    integer(8),          optional,          intent(inout) :: MEMORY_COUNTER(2)
    logical(lg),         optional,          intent(in)    :: ENABLE_EMPTY      !< If empty bins should be considered
    real(rp),            optional, pointer, intent(inout) :: CENTROID(:,:)
    integer(ip)                                           :: ielem,ipoin,pelty
    integer(ip)                                           :: ii,jj,kk
    integer(ip)                                           :: idime,inode
    integer(8)                                            :: memor_loc(2)
    logical(lg)                                           :: if_empty
    real(rp)                                              :: coord_loc(3,8),rnode
    
    if( present(ENABLE_EMPTY) ) then
       if_empty = ENABLE_EMPTY
    else
       if_empty = .true.
    end if
    
    if( present(MEMORY_COUNTER) ) then
       memor_loc = MEMORY_COUNTER
    else
       memor_loc = 0_8
    end if
    !
    ! Mesh dimensions
    !
    call octbin % mesh_dim(ndime,mnode,nelem,npoin,ENABLE_EMPTY)
    rnode = 1.0_rp / real(mnode,rp)
    if( ndime == 2 ) then
       pelty = QUA04
    else
       pelty = HEX08
    end if
    !
    ! Allocate mesh
    !
    if( .not. associated(lnods)    ) call memory_alloca(memor_loc,'LNODS','mesh',lnods,mnode,nelem)
    if( .not. associated(ltype)    ) call memory_alloca(memor_loc,'LTYPE','mesh',ltype,nelem)
    if( .not. associated(coord)    ) call memory_alloca(memor_loc,'COORD','mesh',coord,ndime,npoin)

    if( present(CENTROID) ) then
       if( .not. associated(CENTROID) ) call memory_alloca(memor_loc,'CENTROID','mesh',CENTROID,ndime,nelem)
    end if
    !
    ! Load mesh
    !
    nelem = 0
    npoin = 0
    do kk = 1,octbin % bin % boxip(3)
       do jj = 1,octbin % bin % boxip(2)
          do ii = 1,octbin % bin % boxip(1)
             if( associated(octbin % bin % list(ii,jj,kk) % l) ) then
                
                call octbin % octree(ii,jj,kk) % mesh(&
                     ndime,mnode,ielem,ipoin,&
                     lnods,ltype,coord,&
                     MEMORY_COUNTER,&
                     OFFSET_IELEM=nelem,&
                     OFFSET_IPOIN=npoin,&
                     ENABLE_EMPTY=ENABLE_EMPTY,&
                     CENTROID=CENTROID)
                nelem = nelem + ielem
                npoin = npoin + ipoin
                
             else if( if_empty ) then
                
                call octbin % bin % mesh_coord((/ii,jj,kk/),coord_loc)
                ipoin = npoin+1
                ielem = nelem+1
                do idime = 1,ndime
                   coord(idime,ipoin:ipoin+mnode-1) = coord_loc(idime,1:mnode)                   
                end do
                if( present(CENTROID) ) then
                   do idime = 1,ndime
                      CENTROID(idime,ielem) = sum(coord_loc(idime,1:mnode)) * rnode
                   end do
                end if
                ipoin = npoin+1
                ltype(ielem) = pelty
                do inode = 1,mnode
                   lnods(inode,ielem) = ipoin
                   ipoin = ipoin + 1
                end do
                nelem = nelem + 1
                npoin = npoin + mnode
                
             end if
          end do
       end do
    end do
    
    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor_loc
    
  end subroutine mesh

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-09
  !> @brief   Centroid
  !> @details Get teh centroid of each box
  !> 
  !-----------------------------------------------------------------------
  
  subroutine centroid(octbin,coorc,MEMORY_COUNTER)

    class(maths_octbin),                    intent(inout) :: octbin
    real(rp),                      pointer, intent(inout) :: coorc(:,:)
    integer(8),          optional,          intent(inout) :: MEMORY_COUNTER(2)
    integer(ip)                                           :: nelem,ndime
    integer(ip)                                           :: ii,jj,kk
    integer(8)                                            :: memor_loc(2)

    nelem = 0
    ndime = octbin % bin % dim
    
    do kk = 1,octbin % bin % boxip(3)
       do jj = 1,octbin % bin % boxip(2)
          do ii = 1,octbin % bin % boxip(1)
             if( associated(octbin % bin % list(ii,jj,kk) % l) ) then
                nelem = nelem + octbin % octree(ii,jj,kk) % nleaves
             else
                nelem = nelem + 1
             end if
          end do
       end do
    end do

    if( .not. associated(coorc) ) call memory_alloca(memor_loc,'COORC','mesh',coorc,ndime,nelem)

    nelem = 0
    do kk = 1,octbin % bin % boxip(3)
       do jj = 1,octbin % bin % boxip(2)
          do ii = 1,octbin % bin % boxip(1)
             if( associated(octbin % bin % list(ii,jj,kk) % l) ) then
                call octbin % octree(ii,jj,kk) % centroid(coorc,MEMORY_COUNTER,OFFSET_IELEM=nelem)
                nelem = nelem + octbin % octree(ii,jj,kk) % nleaves
             else
                
             end if
          end do
       end do
    end do
    
  end subroutine centroid
  
end module def_maths_octbin
!> @}
