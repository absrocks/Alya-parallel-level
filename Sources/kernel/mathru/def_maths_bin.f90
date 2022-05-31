!-----------------------------------------------------------------------
!
!> @defgroup Maths_Toolbox
!> Toolbox for bins
!> @{
!> @name    ToolBox for bins
!> @file    def_maths.f90
!> @author  Guillaume Houzeaux
!> @brief   Variables
!> @details Variables fot mod_maths.f90
!
!-----------------------------------------------------------------------

module def_maths_bin

  use def_kintyp_basic, only : ip,rp,lg,i1p
  use def_elmtyp,       only : QUA04
  use def_elmtyp,       only : HEX08
  use mod_memory_basic, only : memory_alloca
  use mod_memory_basic, only : memory_deallo
  use mod_memory_basic, only : memory_size

  real(rp), parameter :: epsil = epsilon(1.0_rp)
  
  type maths_bin
     type(i1p), pointer      :: list(:,:,:)    ! List of elements
     integer(ip)             :: dim            ! Space dimensions
     integer(ip)             :: nfilled        ! Number of filled bins
     real(rp)                :: comin(3)       ! Minimum coordinates
     real(rp)                :: comax(3)       ! Maximum coordinates
     integer(ip)             :: boxip(3)       ! Number of boxes in integer(ip)
     real(rp)                :: boxrp(3)       ! Number of boxes in real(rp)
   contains
     procedure,         pass :: init           ! Initialize all
     procedure,         pass :: fill           ! Allocate 
     procedure,         pass :: deallo         ! Deallocate
     procedure,         pass :: find           ! Find a bin
     procedure,         pass :: inbox          ! If a point is inside boundaing box
     procedure,         pass :: coord          ! Min max coordinates of a bin
     procedure,         pass :: mesh           ! Create a bin mesh
     procedure,         pass :: mesh_dim       ! Dimension of the bin mesh
     procedure,         pass :: mesh_node      ! Node number of the mesh
     procedure,         pass :: mesh_coord     ! Coordinates of bin nodes
     procedure,         pass :: centroid       ! Centroid
  end type maths_bin

  private

  public :: maths_bin

contains

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Bin
  !> @details Fill in a bin from a coordinate field
  !> 
  !-----------------------------------------------------------------------

  subroutine init(bin)

    class(maths_bin), intent(inout) :: bin

    bin % dim     =  0
    bin % nfilled =  0
    bin % comin   =  huge(1.0_rp) * 0.1_rp
    bin % comax   = -huge(1.0_rp) * 0.1_rp
    bin % boxip   =  1_ip
    bin % boxrp   =  1.0_rp
    nullify(bin % list)

  end subroutine init

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Bin
  !> @details Fill in a bin from a coordinate field
  !> 
  !-----------------------------------------------------------------------

  subroutine deallo(bin,MEMORY_COUNTER)

    class(maths_bin),               intent(inout) :: bin
    integer(8),           optional, intent(inout) :: MEMORY_COUNTER(2)
    integer(8)                                    :: memor_loc(2)

    if( present(MEMORY_COUNTER) ) then
       memor_loc = MEMORY_COUNTER
    else
       memor_loc = 0_8
    end if
    call memory_deallo(memor_loc,'BIN % LIST % L','maths_bin',bin % list)
    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor_loc

  end subroutine deallo

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-09
  !> @brief   Return the centroid
  !> @details Compute the centroids of the leaves
  !> 
  !-----------------------------------------------------------------------
  
  subroutine centroid(bin,coorc,MEMORY_COUNTER,OFFSET_IELEM,ENABLE_EMPTY)

    class(maths_bin),                       intent(inout) :: bin   
    real(rp),                      pointer, intent(inout) :: coorc(:,:)
    integer(8),          optional,          intent(inout) :: MEMORY_COUNTER(2)
    integer(ip),         optional,          intent(in)    :: OFFSET_IELEM
    logical(lg),         optional,          intent(in)    :: ENABLE_EMPTY      !< If empty bins should be considered
    logical(lg)                                           :: if_empty
    integer(ip)                                           :: mnode,ndime
    integer(ip)                                           :: nelem,inode
    integer(ip)                                           :: npoin,ielem
    integer(ip)                                           :: ii,jj,kk
    real(rp)                                              :: rnode,xx(3)
    integer(8)                                            :: memor_loc(2)
    real(rp)                                              :: coord(3,8)
    
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
    ! Dimensions
    !
    call bin % mesh_dim(ndime,mnode,nelem,npoin,ENABLE_EMPTY)
    !
    ! Allocate mesh
    !
    if( .not. associated(coorc) ) call memory_alloca(memor_loc,'COORC','mesh',coorc,ndime,nelem)
    !
    ! Centroid
    !
    rnode = 1.0_rp / real(mnode,rp)
    ielem = 0
    if( if_empty ) then
       do kk = 1,bin % boxip(3)
          do jj = 1,bin % boxip(2)
             do ii = 1,bin % boxip(1)
                call bin % mesh_coord((/ii,jj,kk/),coord)
                xx = 0.0_rp
                do inode = 1,mnode
                   xx(1:ndime) = xx(1:ndime) + coord(1:ndime,inode)
                end do
                ielem = ielem + 1
                coorc(1:ndime,ielem) = xx(1:ndime) * rnode
             end do
          end do
       end do
    else
       do kk = 1,bin % boxip(3)
          do jj = 1,bin % boxip(2)
             do ii = 1,bin % boxip(1)
                if( associated(bin % list(ii,jj,kk) % l) ) then
                   call bin % mesh_coord((/ii,jj,kk/),coord)
                   xx = 0.0_rp
                   do inode = 1,mnode
                      xx(1:ndime) = xx(1:ndime) + coord(1:ndime,inode)
                   end do
                   ielem = ielem + 1
                   coorc(1:ndime,ielem) = xx(1:ndime) * rnode
                end if
             end do
          end do
       end do       
    end if
       
  end subroutine centroid
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Bin
  !> @details Fill in a bin from a coordinate field
  !> 
  !-----------------------------------------------------------------------

  subroutine fill(bin,coord,boxes,MEMORY_COUNTER,OFFSET)

    class(maths_bin),               intent(inout) :: bin
    real(rp),              pointer, intent(in)    :: coord(:,:)
    integer(ip),                    intent(in)    :: boxes(:)
    integer(8),  optional,          intent(inout) :: MEMORY_COUNTER(2)
    real(rp),    optional,          intent(in)    :: OFFSET
    integer(ip)                                   :: ii,jj,kk,idime
    integer(ip)                                   :: npoin
    integer(ip)                                   :: ipoin
    real(rp)                                      :: delta(3)
    integer(ip)                                   :: xx(3)
    integer(ip),           pointer                :: num_bin(:,:,:)
    integer(8)                                    :: memor_loc(2)
    real(rp)                                      :: offset_loc

    if( present(OFFSET) ) then
       offset_loc = OFFSET
    else
       offset_loc = 1.0e-3_rp   
    end if

    if( present(MEMORY_COUNTER) ) then
       memor_loc = MEMORY_COUNTER
    else
       memor_loc = 0_8
    end if

    if( associated(coord) ) then

       bin % dim   = size(coord,1)
       npoin       = size(coord,2)
       bin % boxip = 1_ip
       
       if( bin % dim == 1 ) then
          bin % boxip(1)   = boxes(1)
       else if( bin % dim == 2 ) then
          bin % boxip(1:2) = boxes(1:2)
       else
          bin % boxip(1:3) = boxes(1:3)
       end if
       bin % boxrp(1:3) = real(boxes,rp)

       if( .not. associated(bin % list) ) then
          !
          ! Min and max coordinates
          !
          do idime = 1,bin % dim
             bin % comin(idime) = minval(coord(idime,:))
             bin % comax(idime) = maxval(coord(idime,:))
          end do
          delta       = ( bin % comax - bin % comin ) *  offset_loc  
          bin % comin = bin % comin - delta - epsil
          bin % comax = bin % comax + delta + epsil
          !
          ! Count bin elements
          !
          nullify(num_bin)
          call memory_alloca(memor_loc,'NUM_BIN'   ,'maths_bin',num_bin   ,bin % boxip(1),bin % boxip(2),bin % boxip(3))
          call memory_alloca(memor_loc,'BIN % LIST','maths_bin',bin % list,bin % boxip(1),bin % boxip(2),bin % boxip(3))
          do ipoin = 1, npoin
             call bin % find(coord(:,ipoin),xx)
             num_bin(xx(1),xx(2),xx(3)) = num_bin(xx(1),xx(2),xx(3)) + 1
          end do
          !
          ! Allocate bin
          !
          bin % nfilled = 0
          do kk = 1,bin % boxip(3)
             do jj = 1,bin % boxip(2)
                do ii = 1,bin % boxip(1)
                   if( num_bin(ii,jj,kk) /= 0 ) bin % nfilled = bin % nfilled + 1
                   call memory_alloca(memor_loc,'BIN % LIST % L','maths_bin',bin % list(ii,jj,kk) % l,num_bin(ii,jj,kk))
                   num_bin(ii,jj,kk) = 0
                end do
             end do
          end do
          !
          ! Fill in bin
          !
          do ipoin = 1, npoin
             call bin % find(coord(:,ipoin),xx) 
             num_bin(xx(1),xx(2),xx(3)) = num_bin(xx(1),xx(2),xx(3)) + 1
             bin % list(xx(1),xx(2),xx(3)) % l(num_bin(xx(1),xx(2),xx(3))) = ipoin
          end do
          call memory_deallo(memor_loc,'NUM_BIN','maths_bin',num_bin)

       end if

    end if

    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor_loc

  end subroutine fill

  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-09
  !> @brief   Create a mesh
  !> @details Create a mesh from an octbin
  !> 
  !-----------------------------------------------------------------------
  
  subroutine mesh_dim(bin,ndime,mnode,nelem,npoin,ENABLE_EMPTY,permn)

    class(maths_bin),                       intent(inout) :: bin
    integer(ip),                            intent(out)   :: ndime
    integer(ip),                            intent(out)   :: mnode
    integer(ip),                            intent(out)   :: nelem
    integer(ip),                            intent(out)   :: npoin
    logical(lg),         optional,          intent(in)    :: ENABLE_EMPTY !< If empty bins should be considered
    integer(ip),         optional, pointer, intent(inout) :: permn(:)
    logical(lg)                                           :: if_empty
    integer(ip)                                           :: idime,inode,ipoin
    integer(ip)                                           :: ii,jj,kk,kpoin
    integer(ip)                                           :: lnods(8)
    integer(ip),                   pointer                :: permn_loc(:)

    nullify(permn_loc)
    
    if( present(ENABLE_EMPTY) ) then
       if_empty = ENABLE_EMPTY
    else
       if_empty = .true.
    end if
    
    ndime = bin % dim
    mnode = (ndime-1)*4

    if( if_empty ) then
       nelem = 1
       npoin = 1
       do idime = 1,ndime
          nelem = nelem *   bin % boxip(idime)
          npoin = npoin * ( bin % boxip(idime) + 1 )
       end do
    else
       nelem = bin % nfilled
       npoin = 1
       do idime = 1,ndime
          npoin = npoin * ( bin % boxip(idime) + 1 )
       end do
       if( present(permn) ) then
          allocate(permn(npoin))
          permn_loc => permn
       else
          allocate(permn_loc(npoin))
       end if
       do ipoin = 1,npoin
          permn_loc(ipoin) = 0
       end do
       kpoin = 0
       do kk = 1,bin % boxip(3)
          do jj = 1,bin % boxip(2)
             do ii = 1,bin % boxip(1)
                if( associated(bin % list(ii,jj,kk) % l) ) then 
                   call bin % mesh_node((/ii,jj,kk/),lnods)
                   do inode = 1,mnode
                      ipoin = lnods(inode)
                      if( permn_loc(ipoin) == 0 ) then
                         kpoin = kpoin + 1
                         permn_loc(ipoin) = kpoin
                      end if
                   end do
                end if
             end do
          end do
       end do
       npoin = kpoin
       if( .not. present(permn) ) deallocate(permn_loc)
    end if
    
  end subroutine mesh_dim
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-09
  !> @brief   Create a mesh
  !> @details Create a mesh from an octbin
  !> 
  !-----------------------------------------------------------------------
  
  pure subroutine mesh_coord(bin,xx,coord)

    class(maths_bin),  intent(inout) :: bin
    integer(ip),       intent(in)    :: xx(:)
    real(rp),          intent(out)   :: coord(:,:)
    integer(ip)                      :: ipoin,kpoin,ndime
    real(rp)                         :: comin(3),comax(3)
    real(rp)                         :: delta(3)
    
    ndime          = bin % dim
    delta(1:ndime) = ( bin % comax(1:ndime) - bin % comin(1:ndime) ) / bin % boxrp(1:ndime)
    comin(1:ndime) = bin % comin(1:ndime) + real(xx(1:ndime) - 1,rp) * delta(1:ndime)
    comax(1:ndime) = bin % comin(1:ndime) + real(xx(1:ndime)    ,rp) * delta(1:ndime)
    
    if( bin % dim == 2 ) then

       coord(1:ndime,1) = (/ comin(1) , comin(2) /)
       coord(1:ndime,2) = (/ comax(1) , comin(2) /)
       coord(1:ndime,3) = (/ comax(1) , comax(2) /)
       coord(1:ndime,4) = (/ comin(1) , comax(2) /)
              
    else if( bin % dim == 3 ) then

       coord(1:ndime,1) = (/ comin(1) , comin(2) , comin(3) /)
       coord(1:ndime,2) = (/ comax(1) , comin(2) , comin(3) /)
       coord(1:ndime,3) = (/ comax(1) , comax(2) , comin(3) /)
       coord(1:ndime,4) = (/ comin(1) , comax(2) , comin(3) /)
       coord(1:ndime,5) = (/ comin(1) , comin(2) , comax(3) /)
       coord(1:ndime,6) = (/ comax(1) , comin(2) , comax(3) /)
       coord(1:ndime,7) = (/ comax(1) , comax(2) , comax(3) /)
       coord(1:ndime,8) = (/ comin(1) , comax(2) , comax(3) /)
 
    end if
    
  end subroutine mesh_coord
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-09
  !> @brief   Create a mesh
  !> @details Create a mesh from an octbin
  !> 
  !-----------------------------------------------------------------------
  
  pure subroutine mesh_node(bin,xx,lnods)

    class(maths_bin),  intent(inout) :: bin
    integer(ip),       intent(in)    :: xx(:)
    integer(ip),       intent(out)   :: lnods(:)
    integer(ip)                      :: ipoin,kpoin
    
    if( bin % dim == 2 ) then

       ipoin    = xx(1) + (bin % boxip(1)+1)*(xx(2)-1) 
       lnods(1) = ipoin
       lnods(2) = ipoin+1
       lnods(3) = ipoin+bin % boxip(1)+2
       lnods(4) = ipoin+bin % boxip(1)+1
       
    else if( bin % dim == 3 ) then

       kpoin    = (bin % boxip(1)+1)*(bin % boxip(2)+1)
       ipoin    = xx(1) + (bin % boxip(1)+1)*(xx(2)-1) + (bin % boxip(1)+1)*(bin % boxip(2)+1)*(xx(3)-1)
       lnods(1) = ipoin
       lnods(2) = ipoin + 1
       lnods(3) = ipoin + bin % boxip(1)+2
       lnods(4) = ipoin + bin % boxip(1)+1
       lnods(5) = ipoin                    + kpoin
       lnods(6) = ipoin + 1                + kpoin
       lnods(7) = ipoin + bin % boxip(1)+2 + kpoin
       lnods(8) = ipoin + bin % boxip(1)+1 + kpoin              
       
    end if
    
  end subroutine mesh_node
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-09
  !> @brief   Create a mesh
  !> @details Create a mesh from an octbin
  !> 
  !-----------------------------------------------------------------------
  
  subroutine mesh(bin,ndime,mnode,nelem,npoin,lnods,ltype,coord,MEMORY_COUNTER,ENABLE_EMPTY,CENTROID)

    class(maths_bin),                       intent(inout) :: bin
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
    integer(ip)                                           :: ielem,ipoin,idime
    integer(ip)                                           :: ii,jj,kk,kpoin
    integer(8)                                            :: memor_loc(2)
    real(rp)                                              :: comin(3)
    real(rp)                                              :: comax(3)
    real(rp)                                              :: delta(3)
    real(rp)                                              :: yy,zz
    integer(ip),                   pointer                :: permn(:)
    integer(ip)                                           :: lnods_loc(8)
    logical(lg)                                           :: if_empty
    
    nullify(permn)

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

    call bin % mesh_dim(ndime,mnode,nelem,npoin,ENABLE_EMPTY,permn)

    comin(1:ndime) = bin % comin(1:ndime)
    comax(1:ndime) = bin % comax(1:ndime)
    delta(1:ndime) = comax(1:ndime) - comin(1:ndime)
    !
    ! Allocate mesh
    !
    if( .not. associated(lnods) ) call memory_alloca(memor_loc,'LNODS','mesh',lnods,mnode,nelem)
    if( .not. associated(ltype) ) call memory_alloca(memor_loc,'LTYPE','mesh',ltype,nelem)
    if( .not. associated(coord) ) call memory_alloca(memor_loc,'COORD','mesh',coord,ndime,npoin)
    !
    ! Mesh arrays
    !
    ipoin = 0
    ielem = 0
    if( if_empty ) then
       if( ndime == 2 ) then
          do jj = 1,bin % boxip(2)+1
             yy = comin(2) + real(jj-1,rp)/bin % boxrp(2) * delta(2)
             do ii = 1,bin % boxip(1)+1
                ipoin = ipoin + 1
                coord(1,ipoin) = comin(1) + real(ii-1,rp)/bin % boxrp(1) * delta(1)
                coord(2,ipoin) = yy
             end do
          end do
          ipoin = 1
          do jj = 1,bin % boxip(2)
             do ii = 1,bin % boxip(1)
                ielem          = ielem + 1
                lnods(1,ielem) = ipoin
                lnods(2,ielem) = ipoin+1
                lnods(3,ielem) = ipoin+bin % boxip(1)+2
                lnods(4,ielem) = ipoin+bin % boxip(1)+1
                ltype(ielem)   = QUA04
                ipoin          = ipoin + 1
             end do
             ipoin = ipoin + 1
          end do
       else
          do kk = 1,bin % boxip(3)+1
             zz = comin(3) + real(kk-1,rp)/bin % boxrp(3) * delta(3)
             do jj = 1,bin % boxip(2)+1
                yy = comin(2) + real(jj-1,rp)/bin % boxrp(2) * delta(2)
                do ii = 1,bin % boxip(1)+1
                   ipoin = ipoin + 1
                   coord(1,ipoin) = comin(1) + real(ii-1,rp)/bin % boxrp(1) * delta(1)
                   coord(2,ipoin) = yy
                   coord(3,ipoin) = zz
                end do
             end do
          end do
          ipoin = 1
          kpoin = (bin % boxip(1)+1)*(bin % boxip(2)+1)
          do kk = 1,bin % boxip(3)
             do jj = 1,bin % boxip(2)
                do ii = 1,bin % boxip(1)
                   ielem          = ielem + 1                
                   lnods(1,ielem) = ipoin
                   lnods(2,ielem) = ipoin+1
                   lnods(3,ielem) = ipoin+bin % boxip(1)+2
                   lnods(4,ielem) = ipoin+bin % boxip(1)+1
                   lnods(5,ielem) = kpoin + ipoin
                   lnods(6,ielem) = kpoin + ipoin+1
                   lnods(7,ielem) = kpoin + ipoin+bin % boxip(1)+2
                   lnods(8,ielem) = kpoin + ipoin+bin % boxip(1)+1                
                   ltype(ielem)   = HEX08
                   ipoin          = ipoin + 1
                end do
                ipoin = ipoin + 1
             end do
             ipoin = kpoin * kk + 1
          end do
       end if
       
    else

       if( ndime == 2 ) then
          do jj = 1,bin % boxip(2)+1
             yy = comin(2) + real(jj-1,rp)/bin % boxrp(2) * delta(2)
             do ii = 1,bin % boxip(1)+1
                ipoin          = ipoin + 1
                kpoin          = permn(ipoin)
                if( kpoin > 0 ) then
                   coord(1,kpoin) = comin(1) + real(ii-1,rp)/bin % boxrp(1) * delta(1)
                   coord(2,kpoin) = yy
                end if
             end do
          end do
          do jj = 1,bin % boxip(2)
             do ii = 1,bin % boxip(1)
                if( associated(bin % list(ii,jj,1_ip) % l) ) then
                   ielem                = ielem + 1
                   call bin % mesh_node((/ii,jj,1_ip/),lnods_loc)
                   lnods(1:mnode,ielem) = permn(lnods_loc(1:mnode))
                   ltype(ielem)         = QUA04
                end if
             end do
          end do
       else
          do kk = 1,bin % boxip(3)+1
             zz = comin(3) + real(kk-1,rp)/bin % boxrp(3) * delta(3)
             do jj = 1,bin % boxip(2)+1
                yy = comin(2) + real(jj-1,rp)/bin % boxrp(2) * delta(2)
                do ii = 1,bin % boxip(1)+1
                   ipoin          = ipoin + 1
                   kpoin          = permn(ipoin)
                   if( kpoin > 0 ) then
                      coord(1,kpoin) = comin(1) + real(ii-1,rp)/bin % boxrp(1) * delta(1)
                      coord(2,kpoin) = yy
                      coord(3,kpoin) = zz
                   end if
                end do
             end do
          end do
          do kk = 1,bin % boxip(3)
             do jj = 1,bin % boxip(2)
                do ii = 1,bin % boxip(1)
                   if( associated(bin % list(ii,jj,kk) % l) ) then
                      ielem                = ielem + 1
                      call bin % mesh_node((/ii,jj,kk/),lnods_loc)
                      lnods(1:mnode,ielem) = permn(lnods_loc(1:mnode))
                      ltype(ielem)         = HEX08
                   end if
                end do
             end do
          end do
       end if
    end if

    if( associated(permn) ) deallocate(permn)
    if( present(MEMORY_COUNTER) ) MEMORY_COUNTER = memor_loc
    
    if( present(CENTROID) ) then
       call bin % centroid(CENTROID,MEMORY_COUNTER=MEMORY_COUNTER,ENABLE_EMPTY=ENABLE_EMPTY)
    end if

  end subroutine mesh
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Coordinates of a bin
  !> @details Coordinates of a bin
  !> 
  !-----------------------------------------------------------------------

  pure subroutine coord(bin,xx,comin,comax) 

    class(maths_bin), intent(in)  :: bin
    integer(ip),      intent(in)  :: xx(3)    !< Box in x direction
    real(rp),         intent(out) :: comin(:) !< Coordinate of the test point
    real(rp),         intent(out) :: comax(:) !< Coordinate of the test point
    integer(ip)                   :: ndime
    real(rp)                      :: delta(3)
    
    comin          =  huge(1.0_rp) * 0.1_rp
    comax          = -huge(1.0_rp) * 0.1_rp
    ndime          =  bin % dim
    delta(1:ndime) = ( bin % comax(1:ndime) - bin % comin(1:ndime) ) / bin % boxrp(1:ndime)
    
    comin(1:ndime) = bin % comin(1:ndime) + real(xx(1:ndime) - 1,rp) * delta(1:ndime)
    comax(1:ndime) = bin % comin(1:ndime) + real(xx(1:ndime)    ,rp) * delta(1:ndime)
    
  end subroutine coord
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Location in a bin
  !> @details Find the location of a point in a bin
  !> 
  !-----------------------------------------------------------------------

  pure subroutine find(bin,coord,xx) 

    class(maths_bin), intent(in)  :: bin
    real(rp),         intent(in)  :: coord(:) !< Coordinate of the test point
    integer(ip),      intent(out) :: xx(3)    !< Box in x direction

    select case ( bin % dim ) 

    case ( 1_ip )

       xx(1) = int( ( (coord(1)-bin % comin(1)-epsil) / (bin % comax(1)-bin % comin(1)) ) * bin % boxrp(1), ip ) + 1
       xx(1) = min(max(1_ip,xx(1)),bin % boxip(1))
       xx(2) = 1
       xx(3) = 1

    case ( 2_ip )

       xx(1:2) = int( ( (coord(1:2)-bin % comin(1:2)-epsil) / (bin % comax(1:2)-bin % comin(1:2)) ) * bin % boxrp(1:2), ip ) + 1
       xx(1:2) = min(max(1_ip,xx(1:2)),bin % boxip(1:2))
       xx(3)   = 1

    case ( 3_ip )

       xx(1:3) = int( ( (coord(1:3)-bin % comin(1:3)-epsil) / (bin % comax(1:3)-bin % comin(1:3)) ) * bin % boxrp(1:3), ip ) + 1
       xx(1:3) = min(max(1_ip,xx(1:3)),bin % boxip(1:3))

    end select

  end subroutine find
  
  !-----------------------------------------------------------------------
  !> 
  !> @author  houzeaux
  !> @date    2020-06-03
  !> @brief   Location in a bin
  !> @details Find the location of a point in a bin
  !> 
  !-----------------------------------------------------------------------

  pure logical(lg) function inbox(bin,coord) 

    class(maths_bin), intent(in)  :: bin
    real(rp),         intent(in)  :: coord(:) !< Coordinate of the test point

    inbox = .true.

    select case (  bin % dim )

    case ( 1_ip )
       
       if( coord(1) > bin % comax(1) .or.  coord(1) < bin % comin(1) ) inbox = .false.
       
    case ( 2_ip )
       
       if(     coord(1) > bin % comax(1) ) then
          inbox = .false. ; return
       else if( coord(1) < bin % comin(1) ) then
          inbox = .false. ; return
       else if( coord(2) > bin % comax(2) ) then
          inbox = .false. ; return
       else if( coord(2) < bin % comin(2) ) then
          inbox = .false. ; return
       end if
       
     case ( 3_ip )
      
       if(     coord(1) > bin % comax(1) ) then
          inbox = .false. ; return
       else if( coord(1) < bin % comin(1) ) then
          inbox = .false. ; return
       else if( coord(2) > bin % comax(2) ) then
          inbox = .false. ; return
       else if( coord(2) < bin % comin(2) ) then
          inbox = .false. ; return
       else if( coord(3) > bin % comax(3) ) then
          inbox = .false. ; return
       else if( coord(3) < bin % comin(3) ) then
          inbox = .false. ; return
       end if
       
    end select

  end function inbox
  
end module def_maths_bin
!> @}
