!------------------------------------------------------------------------
!>
!> @defgroup Memory_Toolbox
!> @{
!> @name    ToolBox for memory management
!> @file    mod_memory.f90
!> @author  Guillaume Houzeaux
!> @brief   ToolBox for memory management
!> @details ToolBox for memory management: allocate, deallocate
!>
!------------------------------------------------------------------------

module mod_memory

  use def_kintyp, only : ip,rp,lg,i1p,i2p,i3p,r1p,r2p,r3p,r4p,i1pi1p,i1pp,spmat
  use mod_std

  implicit none

  private

  real(rp),      parameter :: Kbytes           = 1024.0_rp
  real(rp),      parameter :: Mbytes           = 1024.0_rp*1024.0_rp
  real(rp),      parameter :: Gbytes           = 1024.0_rp*1024.0_rp*1024.0_rp

  integer(8),    parameter :: Kbytes_i8        = 1024_8
  integer(8),    parameter :: Mbytes_i8        = 1024_8*1024_8
  integer(8),    parameter :: Gbytes_i8        = 1024_8*1024_8*1024_8

  integer(8),    parameter :: size_i1p_type    = 8
  integer(8),    parameter :: size_i1pp_type   = 8
  integer(8),    parameter :: size_i1pi1p_type = 8
  integer(8),    parameter :: size_i2p_type    = 8
  integer(8),    parameter :: size_i3p_type    = 8
  integer(8),    parameter :: size_r1p_type    = 8
  integer(8),    parameter :: size_r2p_type    = 8
  integer(8),    parameter :: size_r3p_type    = 8
  integer(8),    parameter :: size_r4p_type    = 8

  character(26), parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  character(26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'
  
  integer(ip),   parameter :: kfl_debug = 0                 ! Debugging flag
  integer(8)               :: lbytm                         ! Allocated bytes
  integer(8)               :: mem_curre                     ! Current memory allocation in bytes
  integer(8)               :: mem_maxim                     ! Maximum memory allocation in bytes
  integer(ip)              :: mem_alloc                     ! Number of allocations
  integer(ip),   parameter :: kfl_paral = 0                 ! Select the one in def_master instead to debug
  integer(ip)              :: lun_memor                     ! Memory file unit
  integer(ip)              :: kfl_memor                     ! Memory output
  integer(ip)              :: lun_varcount                  ! Variable memory counter file unit
  integer(ip)              :: kfl_varcount                  ! If memory counter is activated
  integer(ip)              :: varcount_number               ! Number of memory-tracked variables
  integer(ip),   parameter :: varcount_max = 1000           ! Number of variables to track memory 
  integer(8)               :: varcount_memory(varcount_max) ! Memory counter of variables
  character(50)            :: varcount_name(varcount_max)   ! Variable names
  !
  ! MEMORY_ALLOC: allocate
  !
  interface memory_alloca
     module procedure &
          &           memory_allr81   , memory_allr82   , memory_allr83   , memory_allr84   , & ! Pointer:     Real(rp)(:), Real(rp)(:,:), Real(rp)(:,:,:), Real(rp)(:,:,:,:)
          &           memory_allr81_8 ,                                                       & !
          &           memory_allr41   , memory_allr42   , memory_allr43   , memory_allr44   , & ! Pointer:     Real(rp)(:), Real(rp)(:,:), Real(rp)(:,:,:), Real(rp)(:,:,:,:)
          &           memory_allr41_8 ,                                                       &
          &           memory_allr161   , memory_allr162   , memory_allr163 , memory_allr164 , & ! Pointer:     Real(rp)(:), Real(rp)(:,:), Real(rp)(:,:,:), Real(rp)(:,:,:,:)
          &           memory_allr161_8 ,                                                      &
          &           memory_alli41   , memory_alli42   , memory_alli43   ,                   & ! Pointer:     Inte(4) (:), Inte(4) (:,:), Inte(4) (:,:,:)
          &           memory_alli418  ,                                                       & ! Pointer:     Inte(4) (:) with special argument
          &           memory_alli81   , memory_alli82   , memory_alli83   ,                   & ! Pointer:     Inte(8) (:), Inte(8) (:,:), Inte(8) (:,:,:)
          &           memory_alli1p_1 , memory_alli2p_1 , memory_alli3p_1 ,                   &
          &           memory_alli1p_2 , memory_alli2p_2 ,                                     &
          &           memory_alli1p_3 ,                                                       &
          &           memory_allr1p_1 , memory_allr2p_1 , memory_allr3p_1 , memory_allr4p_1 , &
          &                             memory_allr2p_2 , memory_allr3p_2 ,                   &
          &           memory_alllg1   , memory_alllg2   ,                                     &
          &           memory_alli1pp_1,                                                       &
          &           memory_alli1pi1p,                                                       &
          &           memory_allcha_4 , memory_allcha_8,  memory_allcha_0,                    &
          &           memory_allxp1,    memory_allxp2,                                        &   ! Pointer:     complex(rp)(:), complex(rp)(:,:)
          &           memory_allspmat , memory_allspmat_1
  end interface memory_alloca

  interface memory_size
     module procedure &
          &           memory_sizerp1     , memory_sizerp2   , memory_sizerp3 ,                  &
          &           memory_sizeip1_4   , memory_sizeip2_4 , memory_sizeip3_4,                 &
          &           memory_sizeip1_8   , memory_sizeip2_8 , memory_sizeip3_8,                 &
          &           memory_sizelg
  end interface memory_size

  interface memory_alloca_min
     module procedure &
          &           memory_alloca_min_rp1,memory_alloca_min_rp2,memory_alloca_min_rp3,memory_alloca_min_rp4,&
          &           memory_alloca_min_ip1,memory_alloca_min_ip2,memory_alloca_min_ip3,&
          &           memory_alloca_min_lg1,memory_alloca_min_lg2,memory_alloca_min_lg3
  end interface memory_alloca_min

  interface memory_deallo
     module procedure &
          &           memory_dear81   , memory_dear82   , memory_dear83   , memory_dear84   , &
          &           memory_dear41   , memory_dear42   , memory_dear43   , memory_dear44   , &
          &           memory_dear161  , memory_dear162  , memory_dear163  , memory_dear164  , &
          &           memory_deai41   , memory_deai42   , memory_deai43   ,                   &
          &           memory_deai81   , memory_deai82   , memory_deai83   ,                   &
          &           memory_deai1p_1 , memory_deai2p_1 , memory_deai3p_1 ,                   &
          &           memory_deai1p_2 , memory_deai2p_2 ,                                     &
          &           memory_deai1p_3 ,                                                       &
          &           memory_deai1pp_1,                                                       &
          &           memory_dealg1,    memory_dealg2,                                        &
          &           memory_dear1p_1 , memory_dear2p_1 , memory_dear3p_1 , memory_dear4p_1 , &
          &                             memory_dear2p_2 , memory_dear3p_2 ,                   &
          &           memory_deacha,    memory_deacha_0,                                      &
          &           memory_deaxp1   , memory_deaxp2,                                        &
          &           memory_deaspmat , memory_deaspmat_1
  end interface memory_deallo

  interface memory_copy
     module procedure &
          &           memory_cpyi41   , memory_cpyi42   , memory_cpyi43   , & ! Inte(4)   (:),  Inte(4) (:,:),  Inte(4) (:,:,:)
          &           memory_cpyi81   , memory_cpyi82   , memory_cpyi83   , & ! Inte(8)   (:),  Inte(8) (:,:),  Inte(8) (:,:,:)
          &           memory_cpyi1p_1 ,                                     & ! Type(i1p) (:)
          &           memory_cpyrp1   , memory_cpyrp2   , memory_cpyrp3       ! Real(rp)  (:),  Real(rp)(:,:),  Real(rp)(:,:,:)
  end interface memory_copy

  interface memory_resize
     module procedure &
          &           memory_newi41   , memory_newi42   ,  memory_newi43  , & ! Inte(4)   (:),  Inte(4) (:,:),  Inte(4)  (:,:,:)
          &           memory_newi81   , memory_newi82   ,  memory_newi83  , & ! Inte(8)   (:),  Inte(8) (:,:),  Inte(8)  (:,:,:)
          &                                                memory_newr83  , & !                                 Real(rp) (:,:,:)
          &           memory_newi1p_1                                         ! Type(i1p) (:)
  end interface memory_resize

  interface memory_renumber
     module procedure &
          &           memory_renrp1   , memory_renrp2   ,  memory_renrp3  , & ! Real(rp)  (:),  Real(rp)(:,:)
          &           memory_reni41   , memory_reni81                         ! Inte(4)   (:),  Inte(8) (:)
  end interface memory_renumber

  public :: memory_initialization             ! Initializaiton: call at the beginning
  public :: memory_alloca                     ! Allocate memory
  public :: memory_size                       ! Gives the size of a pointer (=0 if not associated)
  public :: memory_alloca_min                 ! Allocate a minimum memory for arrays
  public :: memory_deallo                     ! Deallocate memory
  public :: memory_copy                       ! Copy an array
  public :: memory_renumber                   ! Renumber an array
  public :: memory_resize                     ! Resize an array
  public :: memory_unit                       ! Memory scaling and units
  public :: memory_output_info                ! Info about memory
  public :: lbytm                             ! Allocated bytes
  public :: mem_curre                         ! Current memory allocation in bytes
  public :: mem_maxim                         ! Current memory allocation in bytes
  public :: lun_memor                         ! Output unit
  public :: lun_varcount                      ! Variable memory counter unit
  public :: kfl_memor                         ! Output
  public :: kfl_varcount                      ! Output
  public :: mem_alloc                         ! Number of allocation
  public :: Kbytes                            ! Memory units
  public :: Mbytes                            ! Memory units
  public :: Gbytes                            ! Memory units
  public :: memory_add_to_memory_counter      ! Add bytes to memory counter
  public :: memory_remove_from_memory_counter ! Add bytes to memory counter
  public :: memory_output_variable_counter    ! Output variable counter
  
contains

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Variables counter
  !> @details Enable variable counter to track memory
  !>
  !----------------------------------------------------------------------

  integer(ip) function memory_enable_variable_counter()
    
    kfl_varcount = 1
    memory_enable_variable_counter = 0
    
  end function memory_enable_variable_counter

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Output counter
  !> @details Output variable memory counter
  !>
  !----------------------------------------------------------------------

  subroutine memory_output_variable_counter(nunit,OUTPUT_FORMAT)
    
    character(len=*), optional, intent(in) :: OUTPUT_FORMAT
    integer(ip),      optional, intent(in) :: nunit
    integer(4)                             :: nunit4
    integer(ip)                            :: ii

    if( kfl_varcount == 1 ) then

       if( present(nunit) ) then
          nunit4 = int(nunit,4)
       else
          nunit4 = int(lun_varcount,4)
       end if

       if( present(OUTPUT_FORMAT) ) then
          do ii = 1,varcount_number
             write(nunit4,OUTPUT_FORMAT) varcount_name(ii),real(varcount_memory(ii),rp)
          end do
       else
          do ii = 1,varcount_number
             write(nunit4,*) varcount_name(ii),varcount_memory(ii)
          end do
       end if
       
    end if

  end subroutine memory_output_variable_counter
  
  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Variables counter
  !> @details Enable variable counter to track memory
  !>
  !----------------------------------------------------------------------

  subroutine memory_variable_counter(vanam,lbyts)

    character(len=*), intent(in) :: vanam         !< Variable name
    integer(8),       intent(in) :: lbyts         !< Number of bytes
    character(len_trim(vanam))   :: vanam_cap
    integer(ip)                  :: ii(1)
    
    vanam_cap = memory_to_upper_case(trim(vanam))

    ii = maxloc(index(varcount_name,trim(vanam_cap)),MASK=index(varcount_name,trim(vanam_cap))/=0)

    if( ii(1) <= 0 ) then
       varcount_number = varcount_number + 1
       if( varcount_number > varcount_max ) varcount_number = varcount_max
       ii(1) = varcount_number
       varcount_name(varcount_number) = trim(vanam_cap)
    end if

    varcount_memory(ii(1)) = varcount_memory(ii(1)) + lbyts
    
  end subroutine memory_variable_counter

  function memory_to_upper_case (str) Result (string)

    character(*),       intent(In) :: str
    character(len(str))            :: string
    integer(ip)                    :: ic, i

!   Capitalize each letter if it is lowecase
    string = str
    do i = 1, len_trim(str)
        ic = index(low, str(i:i))
        if (ic > 0) string(i:i) = cap(ic:ic)
    end do

  end Function memory_to_upper_case

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   Initialize module variables
  !> @details Initialize all the vairables of this module
  !>
  !----------------------------------------------------------------------

  subroutine memory_initialization()

    lbytm           = 0          ! Allocated bytes
    mem_curre       = 0          ! Current memory allocation in bytes
    mem_maxim       = 0          ! Maximum memory allocation in bytes
    mem_alloc       = 0          ! Number of allocations
    lun_memor       = 0          ! Memory file unit
    kfl_memor       = 0          ! Memory output
    lun_varcount    = 0          ! Variable memory counter unit
    kfl_varcount    = 0          ! Variable emory counter not activated
    varcount_number = 0          ! Number of memory-t<racked variables
    varcount_memory = 0          ! Memory counter of variables

  end subroutine memory_initialization

  !----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    30/01/2018
  !> @brief   This routine writes some info on memory
  !> @details This routine writes some info on memory
  !>
  !----------------------------------------------------------------------

  subroutine memory_output_info(memor,vanam,vacal,vatyp)

    character(*), intent(in)    :: vanam         !< Variable name
    character(*), intent(in)    :: vacal         !< Variable name
    integer(8),   intent(inout) :: memor(2)      !< Memory counter
    character(*), intent(in)    :: vatyp
    integer(4)                  :: lun_memor4
    real(rp)                    :: rbyte,rbyt2
    character(6)                :: lbyte,lbyt2

    if( kfl_memor == 1 ) then

       !call memory_unit(mem_curre ,lbyte,rbyte)
       !call memory_unit(abs(lbytm),lbyt2,rbyt2)
       rbyte = 1.0_rp
       rbyt2 = 1.0_rp

       lun_memor4 = int(lun_memor,4)
       mem_alloc  = mem_alloc + 1
       if( mem_alloc == 1 ) write(lun_memor4,1)
       !
       ! Write info
       !
       write(lun_memor4,2) &
            real(mem_curre,rp),real(lbytm,rp),trim(vanam),&
            trim(vatyp),trim(vacal)
       flush(lun_memor4)

    end if

1   format('# Information on memory in bytes',/,&
         & '# --|  Columns displayed:',/,&
         & '#      1. Current memory      2. Variable memory    3. Variable name                  ',/,&
         & '#      4. Variable type       5. Calling subroutine ',/,&
         & '# ')
2   format(e16.8E3,1x,e16.8E3,1x,a20,1x,a10,1x,a)

  end subroutine memory_output_info

  subroutine memory_allr81(memor,vanam,vacal,varia,ndim1,wzero)
    !
    ! Real(rp)(:)
    !
    implicit none
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(4),   intent(in)              :: ndim1
    real(8),     intent(inout), pointer   :: varia(:)
    character(*), intent(in),    optional :: wzero
    integer(4)                            :: istat
    integer(4)                            :: idim1
    logical(lg)                           :: lzero

    if( ndim1 > 0 ) then

       if(kfl_debug == 1_ip) then
          write(kfl_paral+5000,*) 'trim(vanam),trim(vacal) ',trim(vanam),trim(vacal)
          flush(kfl_paral+5000)
       end if

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       allocate( varia(ndim1) , stat = istat )
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim1 = 1,ndim1
                varia(idim1) = 0.0_8
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'real')

    else

       nullify(varia)

    end if

  end subroutine memory_allr81

  subroutine memory_allr81_8(memor,vanam,vacal,varia,ndim1,wzero)
    !
    ! Real(rp)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(8),   intent(in)            :: ndim1
    real(8),     intent(inout), pointer :: varia(:)
    character(*), intent(in),  optional :: wzero
    integer(4)                          :: istat
    integer(8)                          :: idim1
    logical(lg)                         :: lzero

    if( ndim1 > 0 ) then

       if(kfl_debug == 1_ip) then
          write(kfl_paral+5000,*) 'trim(vanam),trim(vacal) ',trim(vanam),trim(vacal)
          flush(kfl_paral+5000)
       end if

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       allocate( varia(ndim1) , stat = istat )
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim1 = 1,ndim1
                varia(idim1) = 0.0_8
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'real')

    else

       nullify(varia)

    end if

  end subroutine memory_allr81_8

  subroutine memory_allr82(memor,vanam,vacal,varia,ndim1,ndim2,wzero,lboun1,lboun2)
    !
    ! Real(rp)(:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    integer(ip),  intent(in)            :: ndim2
    real(8),     intent(inout), pointer  :: varia(:,:)
    character(*), intent(in),  optional :: wzero
    integer(ip),  intent(in),  optional :: lboun1
    integer(ip),  intent(in),  optional :: lboun2
    integer(4)                          :: istat
    integer(ip)                         :: idim1
    integer(ip)                         :: idim2
    logical(lg)                         :: lzero

    if( ndim1*ndim2 > 0 ) then

       if(kfl_debug == 1_ip) then
          write(kfl_paral+5000,*) 'trim(vanam),trim(vacal) ',trim(vanam),trim(vacal)
          flush(kfl_paral+5000)
       end if

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( present(lboun1) .and. present(lboun2) ) then
          allocate( varia(lboun1:lboun1+ndim1-1,lboun2:lboun2+ndim2-1) , stat = istat )
       else
          allocate( varia(ndim1,ndim2) , stat = istat )
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim2 = lbound(varia,2,KIND=ip),ubound(varia,2,KIND=ip)
                do idim1 = lbound(varia,1,KIND=ip),ubound(varia,1,KIND=ip)
                   varia(idim1,idim2) = 0.0_8
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'real')

    else

       nullify(varia)

    end if

  end subroutine memory_allr82

  subroutine memory_allr83(memor,vanam,vacal,varia,ndim1,ndim2,ndim3,wzero,lboun1,lboun2,lboun3)
    !
    ! Real(rp)(:,:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    integer(ip),  intent(in)            :: ndim2
    integer(ip),  intent(in)            :: ndim3
    real(8),     intent(inout), pointer  :: varia(:,:,:)
    character(*), intent(in),  optional :: wzero
    integer(ip),  intent(in),  optional :: lboun1
    integer(ip),  intent(in),  optional :: lboun2
    integer(ip),  intent(in),  optional :: lboun3
    integer(4)                          :: istat
    integer(ip)                         :: idim1
    integer(ip)                         :: idim2
    integer(ip)                         :: idim3
    logical(lg)                         :: lzero

    if( ndim1*ndim2*ndim3 > 0 ) then

       if(kfl_debug == 1_ip) then
          write(kfl_paral+5000,*) 'trim(vanam),trim(vacal) ',trim(vanam),trim(vacal)
          flush(kfl_paral+5000)
       end if

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( present(lboun1) .and. present(lboun2) .and. present(lboun3) ) then
          allocate( varia(lboun1:lboun1+ndim1-1,lboun2:lboun2+ndim2-1,lboun3:lboun3+ndim3-1) , stat = istat )
       else
          allocate( varia(ndim1,ndim2,ndim3) , stat = istat )
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim3 = lbound(varia,3,KIND=ip),ubound(varia,3,KIND=ip)
                do idim2 = lbound(varia,2,KIND=ip),ubound(varia,2,KIND=ip)
                   do idim1 = lbound(varia,1,KIND=ip),ubound(varia,1,KIND=ip)
                      varia(idim1,idim2,idim3) = 0_8
                   end do
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'real')

    else

       nullify(varia)

    end if
  end subroutine memory_allr83

  subroutine memory_allr84(memor,vanam,vacal,varia,ndim1,ndim2,ndim3,ndim4,wzero,lboun1,lboun2,lboun3,lboun4)
    !
    ! Real(rp)(:,:,:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    integer(ip),  intent(in)            :: ndim2
    integer(ip),  intent(in)            :: ndim3
    integer(ip),  intent(in)            :: ndim4
    real(8),     intent(inout), pointer :: varia(:,:,:,:)
    character(*), intent(in),  optional :: wzero
    integer(ip),  intent(in),  optional :: lboun1
    integer(ip),  intent(in),  optional :: lboun2
    integer(ip),  intent(in),  optional :: lboun3
    integer(ip),  intent(in),  optional :: lboun4
    integer(4)                          :: istat
    logical(lg)                         :: lzero

    if( ndim1*ndim2*ndim3*ndim4 > 0 ) then

       if(kfl_debug == 1_ip) then
          write(kfl_paral+5000,*) 'trim(vanam),trim(vacal) ',trim(vanam),trim(vacal)
          flush(kfl_paral+5000)
       end if

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( present(lboun1) .and. present(lboun2) .and. present(lboun3) .and. present(lboun4) ) then
          allocate( varia(lboun1:lboun1+ndim1-1,lboun2:lboun2+ndim2-1,lboun3:lboun3+ndim3-1,lboun4:lboun4+ndim4-1) , stat = istat )
       else
          allocate( varia(ndim1,ndim2,ndim3,ndim4) , stat = istat )
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) varia = 0.0_8
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'real')

    else

       nullify(varia)

    end if

  end subroutine memory_allr84

  subroutine memory_allr41(memor,vanam,vacal,varia,ndim1,wzero)
    !
    ! Real(rp)(:)
    !
    implicit none
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(4),   intent(in)              :: ndim1
    real(4),     intent(inout), pointer   :: varia(:)
    character(*), intent(in),    optional :: wzero
    integer(4)                            :: istat
    integer(4)                            :: idim1
    logical(lg)                           :: lzero

    if( ndim1 > 0 ) then

       if(kfl_debug == 1_ip) then
          write(kfl_paral+5000,*) 'trim(vanam),trim(vacal) ',trim(vanam),trim(vacal)
          flush(kfl_paral+5000)
       end if

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       allocate( varia(ndim1) , stat = istat )
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim1 = 1,ndim1
                varia(idim1) = 0.0_4
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'real')

    else

       nullify(varia)

    end if

  end subroutine memory_allr41

  subroutine memory_allr41_8(memor,vanam,vacal,varia,ndim1,wzero)
    !
    ! Real(rp)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(8),   intent(in)            :: ndim1
    real(4),     intent(inout), pointer :: varia(:)
    character(*), intent(in),  optional :: wzero
    integer(4)                          :: istat
    integer(8)                          :: idim1
    logical(lg)                         :: lzero

    if( ndim1 > 0 ) then

       if(kfl_debug == 1_ip) then
          write(kfl_paral+5000,*) 'trim(vanam),trim(vacal) ',trim(vanam),trim(vacal)
          flush(kfl_paral+5000)
       end if

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       allocate( varia(ndim1) , stat = istat )
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim1 = 1,ndim1
                varia(idim1) = 0.0_4
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'real')

    else

       nullify(varia)

    end if

  end subroutine memory_allr41_8

  subroutine memory_allr42(memor,vanam,vacal,varia,ndim1,ndim2,wzero,lboun1,lboun2)
    !
    ! Real(rp)(:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    integer(ip),  intent(in)            :: ndim2
    real(4),     intent(inout), pointer  :: varia(:,:)
    character(*), intent(in),  optional :: wzero
    integer(ip),  intent(in),  optional :: lboun1
    integer(ip),  intent(in),  optional :: lboun2
    integer(4)                          :: istat
    integer(ip)                         :: idim1
    integer(ip)                         :: idim2
    logical(lg)                         :: lzero

    if( ndim1*ndim2 > 0 ) then

       if(kfl_debug == 1_ip) then
          write(kfl_paral+5000,*) 'trim(vanam),trim(vacal) ',trim(vanam),trim(vacal)
          flush(kfl_paral+5000)
       end if

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( present(lboun1) .and. present(lboun2) ) then
          allocate( varia(lboun1:lboun1+ndim1-1,lboun2:lboun2+ndim2-1) , stat = istat )
       else
          allocate( varia(ndim1,ndim2) , stat = istat )
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim2 = lbound(varia,2,KIND=ip),ubound(varia,2,KIND=ip)
                do idim1 = lbound(varia,1,KIND=ip),ubound(varia,1,KIND=ip)
                   varia(idim1,idim2) = 0.0_4
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'real')

    else

       nullify(varia)

    end if

  end subroutine memory_allr42

  subroutine memory_allr43(memor,vanam,vacal,varia,ndim1,ndim2,ndim3,wzero,lboun1,lboun2,lboun3)
    !
    ! Real(rp)(:,:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    integer(ip),  intent(in)            :: ndim2
    integer(ip),  intent(in)            :: ndim3
    real(4),     intent(inout), pointer  :: varia(:,:,:)
    character(*), intent(in),  optional :: wzero
    integer(ip),  intent(in),  optional :: lboun1
    integer(ip),  intent(in),  optional :: lboun2
    integer(ip),  intent(in),  optional :: lboun3
    integer(4)                          :: istat
    integer(ip)                         :: idim1
    integer(ip)                         :: idim2
    integer(ip)                         :: idim3
    logical(lg)                         :: lzero

    if( ndim1*ndim2*ndim3 > 0 ) then

       if(kfl_debug == 1_ip) then
          write(kfl_paral+5000,*) 'trim(vanam),trim(vacal) ',trim(vanam),trim(vacal)
          flush(kfl_paral+5000)
       end if

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( present(lboun1) .and. present(lboun2) .and. present(lboun3) ) then
          allocate( varia(lboun1:lboun1+ndim1-1,lboun2:lboun2+ndim2-1,lboun3:lboun3+ndim3-1) , stat = istat )
       else
          allocate( varia(ndim1,ndim2,ndim3) , stat = istat )
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim3 = lbound(varia,3,KIND=ip),ubound(varia,3,KIND=ip)
                do idim2 = lbound(varia,2,KIND=ip),ubound(varia,2,KIND=ip)
                   do idim1 = lbound(varia,1,KIND=ip),ubound(varia,1,KIND=ip)
                      varia(idim1,idim2,idim3) = 0_4
                   end do
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'real')

    else

       nullify(varia)

    end if
  end subroutine memory_allr43

  subroutine memory_allr44(memor,vanam,vacal,varia,ndim1,ndim2,ndim3,ndim4,wzero,lboun1,lboun2,lboun3,lboun4)
    !
    ! Real(rp)(:,:,:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    integer(ip),  intent(in)            :: ndim2
    integer(ip),  intent(in)            :: ndim3
    integer(ip),  intent(in)            :: ndim4
    real(4),     intent(inout), pointer :: varia(:,:,:,:)
    character(*), intent(in),  optional :: wzero
    integer(ip),  intent(in),  optional :: lboun1
    integer(ip),  intent(in),  optional :: lboun2
    integer(ip),  intent(in),  optional :: lboun3
    integer(ip),  intent(in),  optional :: lboun4
    integer(4)                          :: istat
    logical(lg)                         :: lzero

    if( ndim1*ndim2*ndim3*ndim4 > 0 ) then

       if(kfl_debug == 1_ip) then
          write(kfl_paral+5000,*) 'trim(vanam),trim(vacal) ',trim(vanam),trim(vacal)
          flush(kfl_paral+5000)
       end if

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( present(lboun1) .and. present(lboun2) .and. present(lboun3) .and. present(lboun4) ) then
          allocate( varia(lboun1:lboun1+ndim1-1,lboun2:lboun2+ndim2-1,lboun3:lboun3+ndim3-1,lboun4:lboun4+ndim4-1) , stat = istat )
       else
          allocate( varia(ndim1,ndim2,ndim3,ndim4) , stat = istat )
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) varia = 0.0_4
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'real')

    else

       nullify(varia)

    end if

  end subroutine memory_allr44

#ifdef __PGI

  subroutine memory_allr161(dummy)
    implicit none
    type(pgi_dummy1), intent(in)          :: dummy
  end subroutine memory_allr161

  subroutine memory_allr161_8(dummy)
    implicit none
    type(pgi_dummy2), intent(in)          :: dummy
  end subroutine memory_allr161_8

  subroutine memory_allr162(dummy)
    implicit none
    type(pgi_dummy3), intent(in)          :: dummy
  end subroutine memory_allr162

  subroutine memory_allr163(dummy)
    implicit none
    type(pgi_dummy4), intent(in)          :: dummy
  end subroutine memory_allr163

  subroutine memory_allr164(dummy)
    implicit none
    type(pgi_dummy5), intent(in)          :: dummy
  end subroutine memory_allr164

#else

  subroutine memory_allr161(memor,vanam,vacal,varia,ndim1,wzero)
    !
    ! Real(rp)(:)
    !
    implicit none
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(4),   intent(in)              :: ndim1
    real(16),     intent(inout), pointer   :: varia(:)
    character(*), intent(in),    optional :: wzero
    integer(4)                            :: istat
    integer(4)                            :: idim1
    logical(lg)                           :: lzero

    if( ndim1 > 0 ) then

       if(kfl_debug == 1_ip) then
          write(kfl_paral+5000,*) 'trim(vanam),trim(vacal) ',trim(vanam),trim(vacal)
          flush(kfl_paral+5000)
       end if

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       allocate( varia(ndim1) , stat = istat )
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim1 = 1,ndim1
                varia(idim1) = 0.0_16
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'real')

    else

       nullify(varia)

    end if

  end subroutine memory_allr161

  subroutine memory_allr161_8(memor,vanam,vacal,varia,ndim1,wzero)
    !
    ! Real(rp)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(8),   intent(in)            :: ndim1
    real(16),     intent(inout), pointer :: varia(:)
    character(*), intent(in),  optional :: wzero
    integer(4)                          :: istat
    integer(8)                          :: idim1
    logical(lg)                         :: lzero

    if( ndim1 > 0 ) then

       if(kfl_debug == 1_ip) then
          write(kfl_paral+5000,*) 'trim(vanam),trim(vacal) ',trim(vanam),trim(vacal)
          flush(kfl_paral+5000)
       end if

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       allocate( varia(ndim1) , stat = istat )
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim1 = 1,ndim1
                varia(idim1) = 0.0_16
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'real')

    else

       nullify(varia)

    end if

  end subroutine memory_allr161_8

  subroutine memory_allr162(memor,vanam,vacal,varia,ndim1,ndim2,wzero,lboun1,lboun2)
    !
    ! Real(rp)(:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    integer(ip),  intent(in)            :: ndim2
    real(16),     intent(inout), pointer  :: varia(:,:)
    character(*), intent(in),  optional :: wzero
    integer(ip),  intent(in),  optional :: lboun1
    integer(ip),  intent(in),  optional :: lboun2
    integer(4)                          :: istat
    integer(ip)                         :: idim1
    integer(ip)                         :: idim2
    logical(lg)                         :: lzero

    if( ndim1*ndim2 > 0 ) then

       if(kfl_debug == 1_ip) then
          write(kfl_paral+5000,*) 'trim(vanam),trim(vacal) ',trim(vanam),trim(vacal)
          flush(kfl_paral+5000)
       end if

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( present(lboun1) .and. present(lboun2) ) then
          allocate( varia(lboun1:lboun1+ndim1-1,lboun2:lboun2+ndim2-1) , stat = istat )
       else
          allocate( varia(ndim1,ndim2) , stat = istat )
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim2 = lbound(varia,2,KIND=ip),ubound(varia,2,KIND=ip)
                do idim1 = lbound(varia,1,KIND=ip),ubound(varia,1,KIND=ip)
                   varia(idim1,idim2) = 0.0_16
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'real')

    else

       nullify(varia)

    end if

  end subroutine memory_allr162

  subroutine memory_allr163(memor,vanam,vacal,varia,ndim1,ndim2,ndim3,wzero,lboun1,lboun2,lboun3)
    !
    ! Real(rp)(:,:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    integer(ip),  intent(in)            :: ndim2
    integer(ip),  intent(in)            :: ndim3
    real(16),     intent(inout), pointer  :: varia(:,:,:)
    character(*), intent(in),  optional :: wzero
    integer(ip),  intent(in),  optional :: lboun1
    integer(ip),  intent(in),  optional :: lboun2
    integer(ip),  intent(in),  optional :: lboun3
    integer(4)                          :: istat
    integer(ip)                         :: idim1
    integer(ip)                         :: idim2
    integer(ip)                         :: idim3
    logical(lg)                         :: lzero

    if( ndim1*ndim2*ndim3 > 0 ) then

       if(kfl_debug == 1_ip) then
          write(kfl_paral+5000,*) 'trim(vanam),trim(vacal) ',trim(vanam),trim(vacal)
          flush(kfl_paral+5000)
       end if

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( present(lboun1) .and. present(lboun2) .and. present(lboun3) ) then
          allocate( varia(lboun1:lboun1+ndim1-1,lboun2:lboun2+ndim2-1,lboun3:lboun3+ndim3-1) , stat = istat )
       else
          allocate( varia(ndim1,ndim2,ndim3) , stat = istat )
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim3 = lbound(varia,3,KIND=ip),ubound(varia,3,KIND=ip)
                do idim2 = lbound(varia,2,KIND=ip),ubound(varia,2,KIND=ip)
                   do idim1 = lbound(varia,1,KIND=ip),ubound(varia,1,KIND=ip)
                      varia(idim1,idim2,idim3) = 0.0_16
                   end do
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'real')

    else

       nullify(varia)

    end if
  end subroutine memory_allr163

  subroutine memory_allr164(memor,vanam,vacal,varia,ndim1,ndim2,ndim3,ndim4,wzero,lboun1,lboun2,lboun3,lboun4)
    !
    ! Real(rp)(:,:,:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    integer(ip),  intent(in)            :: ndim2
    integer(ip),  intent(in)            :: ndim3
    integer(ip),  intent(in)            :: ndim4
    real(16),     intent(inout), pointer :: varia(:,:,:,:)
    character(*), intent(in),  optional :: wzero
    integer(ip),  intent(in),  optional :: lboun1
    integer(ip),  intent(in),  optional :: lboun2
    integer(ip),  intent(in),  optional :: lboun3
    integer(ip),  intent(in),  optional :: lboun4
    integer(4)                          :: istat
    logical(lg)                         :: lzero

    if( ndim1*ndim2*ndim3*ndim4 > 0 ) then

       if(kfl_debug == 1_ip) then
          write(kfl_paral+5000,*) 'trim(vanam),trim(vacal) ',trim(vanam),trim(vacal)
          flush(kfl_paral+5000)
       end if

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( present(lboun1) .and. present(lboun2) .and. present(lboun3) .and. present(lboun4) ) then
          allocate( varia(lboun1:lboun1+ndim1-1,lboun2:lboun2+ndim2-1,lboun3:lboun3+ndim3-1,lboun4:lboun4+ndim4-1) , stat = istat )
       else
          allocate( varia(ndim1,ndim2,ndim3,ndim4) , stat = istat )
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) varia = 0.0_16
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'real')

    else

       nullify(varia)

    end if

  end subroutine memory_allr164

#endif


  subroutine memory_alli41(memor,vanam,vacal,varia,ndim1,wzero,lboun)
    !
    ! Integer(4)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(4),   intent(in)            :: ndim1
    integer(4),   intent(inout), pointer  :: varia(:)
    character(*), intent(in),  optional :: wzero
    integer(4),   intent(in),  optional :: lboun
    integer(4)                          :: istat
    integer(ip)                         :: idim1
    logical(lg)                         :: lzero
    logical(lg)                         :: linde
    integer(4)                          :: xvalu

    if( ndim1 > 0 ) then

       if( present(wzero) ) then
          if( wzero == 'REALLOCATE') then
             call memory_deallo(memor,vanam,vacal,varia)
          end if
       end if

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       lzero = .true.
       linde = .false.
       xvalu = 0_4
       if( present(wzero) ) then
          if(      trim(wzero) == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          else if( trim(wzero) == 'HUGE') then
             xvalu = huge(0_4)
          else if( trim(wzero) == 'IDENTITY') then
             lzero = .false.
             linde = .true.
          end if
       end if

       if( present(lboun) ) then
          allocate( varia(lboun:lboun+ndim1-1) , stat = istat )
       else
          allocate( varia(ndim1) , stat = istat )
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim1 = lbound(varia,1,KIND=ip),ubound(varia,1,KIND=ip)
                varia(idim1) = xvalu
             end do
          else if( linde ) then
             do idim1 = lbound(varia,1,KIND=ip),ubound(varia,1,KIND=ip)
                varia(idim1) = int(idim1,4)
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'integer(4)')

    else

       nullify(varia)

    end if

  end subroutine memory_alli41



  subroutine memory_alli418(memor,vanam,vacal,varia,ndim1,wzero)
    !
    ! Integer(4)(:): wrt memory_alli41, the dimension is integer(8)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(8),   intent(in)            :: ndim1
    integer(4),   intent(inout), pointer  :: varia(:)
    character(*), intent(in),  optional :: wzero
    integer(4)                          :: istat
    integer(ip)                         :: idim1
    logical(lg)                         :: lzero
    logical(lg)                         :: linde
    integer(4)                          :: xvalu

    if( ndim1 > 0 ) then

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       lzero = .true.
       linde = .false.
       xvalu = 0_4
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          else if( wzero == 'HUGE') then
             xvalu = huge(0_4)
          else if( wzero == 'IDENTITY') then
             lzero = .false.
             linde = .true.
          end if
       end if

       allocate( varia(ndim1) , stat = istat )

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim1 = 1,ndim1
                varia(idim1) = xvalu
             end do
          else if( linde ) then
             do idim1 = lbound(varia,1,KIND=ip),ubound(varia,1,KIND=ip)
                varia(idim1) = int(idim1,4)
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'integer(4)')

    else

       nullify(varia)

    end if

  end subroutine memory_alli418

  subroutine memory_alli42(memor,vanam,vacal,varia,ndim1,ndim2,wzero,lboun1,lboun2)
    !
    ! Integer(4)(:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(4),   intent(in)            :: ndim1
    integer(4),   intent(in)            :: ndim2
    integer(4),   intent(inout), pointer  :: varia(:,:)
    character(*), intent(in),  optional :: wzero
    integer(4),   intent(in),  optional :: lboun1
    integer(4),   intent(in),  optional :: lboun2
    integer(4)                          :: istat
    integer(4)                          :: idim1
    integer(4)                          :: idim2
    logical(lg)                         :: lzero

    if( ndim1*ndim2 > 0 ) then

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( present(lboun1) .and. present(lboun2) ) then
          allocate( varia(lboun1:lboun1+ndim1-1,lboun2:lboun2+ndim2-1) , stat = istat )
       else
          allocate( varia(ndim1,ndim2) , stat = istat )
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim2 = lbound(varia,2),ubound(varia,2)
                do idim1 = lbound(varia,1),ubound(varia,1)
                   varia(idim1,idim2) = 0_ip
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'integer(4)')

    else

       nullify(varia)

    end if

  end subroutine memory_alli42

  subroutine memory_alli43(memor,vanam,vacal,varia,ndim1,ndim2,ndim3,wzero,lboun1,lboun2,lboun3)
    !
    ! Integer(4)(:,:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(4),   intent(in)            :: ndim1
    integer(4),   intent(in)            :: ndim2
    integer(4),   intent(in)            :: ndim3
    integer(4),   intent(inout), pointer  :: varia(:,:,:)
    character(*), intent(in),  optional :: wzero
    integer(4),   intent(in),  optional :: lboun1
    integer(4),   intent(in),  optional :: lboun2
    integer(4),   intent(in) , optional :: lboun3
    integer(4)                          :: istat
    integer(8)                          :: idim1
    integer(8)                          :: idim2
    integer(8)                          :: idim3
    logical(lg)                         :: lzero

    if( ndim1*ndim2*ndim3 > 0 ) then

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( present(lboun1) .and. present(lboun2) .and. present(lboun3) ) then
          allocate( varia(lboun1:lboun1+ndim1-1,lboun2:lboun2+ndim2-1,lboun3:lboun3+ndim3-1) , stat = istat )
       else
          allocate( varia(ndim1,ndim2,ndim3) , stat = istat )
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim3 = lbound(varia,3,KIND=8),ubound(varia,3,KIND=8)
                do idim2 = lbound(varia,2,KIND=8),ubound(varia,2,KIND=8)
                   do idim1 = lbound(varia,1,KIND=8),ubound(varia,1,KIND=8)
                      varia(idim1,idim2,idim3) = 0_4
                   end do
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'integer(4)')

    else

       nullify(varia)

    end if

  end subroutine memory_alli43

  subroutine memory_alllg1(memor,vanam,vacal,varia,ndim1,wzero,lboun)
    !
    ! Logical(lg)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    logical(lg),  intent(inout), pointer  :: varia(:)
    character(*), intent(in),  optional :: wzero
    integer(ip),  intent(in),  optional :: lboun
    integer(4)                          :: istat
    integer(8)                          :: idim1
    logical(lg)                         :: lzero

    if( ndim1 > 0 ) then

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( present(lboun) ) then
          allocate( varia(lboun:lboun+ndim1-1) , stat = istat )
       else
          allocate( varia(ndim1) , stat = istat )
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim1 = lbound(varia,1,KIND=8),ubound(varia,1,KIND=8)
                varia(idim1) = .false.
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'logical')

    else

       nullify(varia)

    end if

  end subroutine memory_alllg1

  subroutine memory_alllg2(memor,vanam,vacal,varia,ndim1,ndim2,wzero,lboun1,lboun2)
    !
    ! Logical(lg)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    integer(ip),  intent(in)            :: ndim2
    logical(lg),  intent(inout), pointer  :: varia(:,:)
    character(*), intent(in),  optional :: wzero
    integer(ip),  intent(in),  optional :: lboun1
    integer(ip),  intent(in),  optional :: lboun2
    integer(4)                          :: istat
    logical(lg)                         :: lzero

    if( ndim1*ndim2 > 0 ) then

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( present(lboun1) .and. present(lboun2) ) then
          allocate( varia(lboun1:lboun1+ndim1-1,lboun2:lboun2+ndim2-1) , stat = istat )
       else
          allocate( varia(ndim1,ndim2) , stat = istat )
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             varia = .false.
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'logical')

    else

       nullify(varia)

    end if

  end subroutine memory_alllg2

  subroutine memory_allr1p_1(memor,vanam,vacal,varia,ndim1,wzero,lboun1)
    !
    ! Type(r1p)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    type(r1p),    intent(inout), pointer  :: varia(:)
    character(*), intent(in),  optional :: wzero
    integer(ip),  intent(in),  optional :: lboun1
    integer(ip)                         :: idim1
    integer(4)                          :: istat
    logical(lg)                         :: lzero

    if( ndim1 > 0 ) then

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       if( present(lboun1) ) then
          allocate( varia(lboun1:lboun1+ndim1-1) , stat = istat )
       else
          allocate( varia(ndim1) , stat = istat )
       end if
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*size_r1p_type
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim1 = lbound(varia,1),ubound(varia,1)
                nullify(varia(idim1) % a)
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type(r1p)')

    else

       nullify(varia)

    end if

  end subroutine memory_allr1p_1

  subroutine memory_allr2p_1(memor,vanam,vacal,varia,ndim1,wzero,lboun1)
    !
    ! Type(r2p)(:)
    !
    implicit none
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(ip),  intent(in)              :: ndim1
    type(r2p),    intent(inout), pointer  :: varia(:)
    character(*), intent(in),    optional :: wzero
    integer(ip),  intent(in),    optional :: lboun1
    integer(ip)                           :: idim1
    integer(4)                            :: istat
    logical(lg)                           :: lzero

    if( ndim1 > 0 ) then

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       if( present(lboun1) ) then
          allocate( varia(lboun1:lboun1+ndim1-1) , stat = istat )
       else
          allocate( varia(ndim1) , stat = istat )
       end if

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*size_r2p_type
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim1 = lbound(varia,1),ubound(varia,1)
                nullify(varia(idim1) % a)
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type(r2p)')

    else

       nullify(varia)

    end if

  end subroutine memory_allr2p_1

  subroutine memory_allr2p_2(memor,vanam,vacal,varia,ndim1,ndim2,wzero)
    !
    ! Type(r2p)(:,:)
    !
    implicit none
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(ip),  intent(in)              :: ndim1
    integer(ip),  intent(in)              :: ndim2
    type(r2p),    intent(inout), pointer  :: varia(:,:)
    character(*), intent(in),    optional :: wzero
    integer(ip)                           :: idim1
    integer(ip)                           :: idim2
    integer(4)                            :: istat
    logical(lg)                           :: lzero

    if( ndim1*ndim2 > 0 ) then

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       allocate( varia(ndim1,ndim2) , stat = istat )
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*size_r2p_type
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim2 = 1,ndim2
                do idim1 = 1,ndim1
                   nullify(varia(idim1,idim2) % a)
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type(r2p)')

    else

       nullify(varia)

    end if

  end subroutine memory_allr2p_2

  subroutine memory_allr3p_1(memor,vanam,vacal,varia,ndim1,wzero)
    !
    ! Type(r3p)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    type(r3p),    intent(inout), pointer  :: varia(:)
    character(*), intent(in),  optional :: wzero
    integer(ip)                         :: idim1
    integer(4)                          :: istat
    logical(lg)                         :: lzero

    if( ndim1 > 0 ) then

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       allocate( varia(ndim1) , stat = istat )
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*size_r3p_type
          !lbytm = int(sizeof(varia),8) ! ndim1*ip
          if( lzero ) then
             do idim1 = 1,ndim1
                nullify(varia(idim1) % a)
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type(r3p)')

    else

       nullify(varia)

    end if

  end subroutine memory_allr3p_1

  subroutine memory_allr4p_1(memor,vanam,vacal,varia,ndim1,wzero)
    !
    ! Type(r4p)(:)
    !
    implicit none
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(ip),  intent(in)              :: ndim1
    type(r4p),    intent(inout), pointer  :: varia(:)
    character(*), intent(in),    optional :: wzero
    integer(ip)                           :: idim1
    integer(4)                            :: istat
    logical(lg)                           :: lzero

    if( ndim1 > 0 ) then

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       allocate( varia(ndim1) , stat = istat )
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*size_r4p_type
          !lbytm = int(sizeof(varia),8) ! ndim1*ip
          if( lzero ) then
             do idim1 = 1,ndim1
                nullify(varia(idim1) % a)
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type(r4p)')

    else

       nullify(varia)

    end if

  end subroutine memory_allr4p_1

  subroutine memory_allspmat_1(memor,vanam,vacal,varia,ndim1,wzero)
     !
     ! Type(spmat)
     !
     implicit none
     integer(8),   intent(inout)           :: memor(2)      !< Memory counter
     character(*), intent(in)              :: vanam         !< Variable name
     character(*), intent(in)              :: vacal         !< Calling subroutine name
     type(spmat),  intent(inout), pointer  :: varia(:)
     integer(ip),  intent(in)              :: ndim1
     character(*), intent(in),    optional :: wzero
     integer(ip)                           :: idim1
     integer(4)                            :: istat
     logical(lg)                           :: lzero

     if( ndim1 > 0 ) then

        if( associated(varia)) then
           write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
           call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
        end if

        istat=0
        allocate( varia(ndim1) , stat = istat )

        lzero = .true.
        if(present(wzero)) then
           if( wzero == 'DO_NOT_INITIALIZE') then
              lzero = .false.
           end if
        end if

        if( istat == 0 ) then
           lbytm = int(ndim1*ip,KIND=8)
           if( lzero ) then
              do idim1 = 1,ndim1
                 varia(idim1) % ndof  = 0
                 varia(idim1) % nrows = 0
                 varia(idim1) % ncols = 0
                 nullify(varia(idim1) % iA)
                 nullify(varia(idim1) % jA)
                 nullify(varia(idim1) % vA)
              end do
           end if
        else
           call memory_error(0_ip,vanam,vacal,istat)
        end if

        call memory_info(memor,vanam,vacal,'type(spmat)')

     else

        nullify(varia)

     end if

  end subroutine memory_allspmat_1

  subroutine memory_allspmat(memor,vanam,vacal,varia,ndof,ndim1,wzero)
     !
     ! Type(spmat)
     !
     implicit none
     integer(8),   intent(inout)         :: memor(2)      !< Memory counter
     character(*), intent(in)            :: vanam         !< Variable name
     character(*), intent(in)            :: vacal         !< Calling subroutine name
     type(spmat), intent(out)            :: varia
     integer(ip),  intent(in)            :: ndim1
     integer(ip),  intent(in)            :: ndof
     character(*), intent(in),  optional :: wzero
     integer(4)                          :: istat
     integer(4)                          :: aux_istat
     logical(lg)                         :: lzero

     if( ndim1 > 0 ) then

        if( associated(varia % iA) .or. associated(varia % jA) .or. associated(varia % vA) ) then
           write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
           call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
        end if

        istat=0
        allocate( varia % iA(ndim1) , stat = aux_istat )
        istat = istat + aux_istat
        allocate( varia % jA(ndim1) , stat = aux_istat )
        istat = istat + aux_istat
        allocate( varia % vA(ndof, ndof, ndim1) , stat = aux_istat )
        istat = istat + aux_istat
        lzero = .true.
        if(present(wzero)) then
           if( wzero == 'DO_NOT_INITIALIZE') then
              lzero = .false.
           end if
        end if

        if( istat == 0 ) then
           lbytm = int(ndim1 * (ip * 2 + rp * ndof * ndof),KIND=8)
           if( lzero ) then
              varia % ndof  = ndof
              varia % nrows = 0
              varia % ncols = 0
              varia % iA(1:ndim1) = 0_ip
              varia % jA(1:ndim1) = 0_ip
              varia % vA(1:ndof,1:ndof,1:ndim1) = 0_ip
           end if
        else
           call memory_error(0_ip,vanam,vacal,istat)
        end if

        call memory_info(memor,vanam,vacal,'type(spmat)')

     else

        varia % ndof  = 0
        varia % nrows = 0
        varia % ncols = 0
        nullify(varia % iA)
        nullify(varia % jA)
        nullify(varia % vA)

     end if

  end subroutine memory_allspmat

  subroutine memory_allr3p_2(memor,vanam,vacal,varia,ndim1,ndim2,wzero)
    !
    ! Type(r3p)(:,:)
    !
    implicit none
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(ip),  intent(in)              :: ndim1
    integer(ip),  intent(in)              :: ndim2
    type(r3p),    intent(inout), pointer  :: varia(:,:)
    character(*), intent(in),    optional :: wzero
    integer(ip)                           :: idim1,idim2
    integer(4)                            :: istat
    logical(lg)                           :: lzero

    if( ndim1*ndim2 > 0 ) then

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       allocate( varia(ndim1,ndim2) , stat = istat )
       if( istat == 0 ) then
          lbytm = int(ndim1,KIND=8)*int(ndim2,KIND=8)*int(ip,KIND=8)
          if( lzero ) then
             do idim2 = 1,ndim2
                do idim1 = 1,ndim1
                   nullify(varia(idim1,idim2) % a)
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type(r3p)')

    else

       nullify(varia)

    end if

  end subroutine memory_allr3p_2

  subroutine memory_alli1p_2(memor,vanam,vacal,varia,ndim1,ndim2,wzero,lboun1,lboun2)
    !
    ! Type(i1p)(:,:)
    !
    implicit none
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    integer(ip),  intent(in)              :: ndim1
    integer(ip),  intent(in)              :: ndim2
    type(i1p),    intent(inout), pointer  :: varia(:,:)
    character(*), intent(in),    optional :: wzero
    integer(ip),  intent(in),    optional :: lboun1
    integer(ip),  intent(in),    optional :: lboun2
    integer(ip)                           :: idim1
    integer(ip)                           :: idim2
    integer(4)                            :: istat
    logical(lg)                           :: lzero

    if( ndim1*ndim2 > 0 ) then

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( present(lboun1) .and. present(lboun2) ) then
          allocate( varia(lboun1:lboun1+ndim1-1,lboun2:lboun2+ndim2-1) , stat = istat )
       else
          allocate( varia(ndim1,ndim2) , stat = istat )
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)
          if( lzero ) then
             do idim2 = lbound(varia,2),ubound(varia,2)
                do idim1 = lbound(varia,1),ubound(varia,1)
                   nullify(varia(idim1,idim2) % l)
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type(i1p)')

    else

       nullify(varia)

    end if

  end subroutine memory_alli1p_2

  subroutine memory_alli1p_3(memor,vanam,vacal,varia,ndim1,ndim2,ndim3,wzero,lboun1,lboun2,lboun3)
    !
    ! Type(i1p)(:,:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    integer(ip),  intent(in)            :: ndim2
    integer(ip),  intent(in)            :: ndim3
    type(i1p),    intent(inout), pointer  :: varia(:,:,:)
    character(*), intent(in),  optional :: wzero
    integer(ip),  intent(in),  optional :: lboun1
    integer(ip),  intent(in),  optional :: lboun2
    integer(ip),  intent(in),  optional :: lboun3
    integer(ip)                         :: idim1
    integer(ip)                         :: idim2
    integer(ip)                         :: idim3
    integer(4)                          :: istat
    logical(lg)                         :: lzero

    if( ndim1*ndim2*ndim3 > 0 ) then

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       if( present(lboun1) .and. present(lboun2) .and. present(lboun3) ) then
          allocate( varia(lboun1:lboun1+ndim1-1,lboun2:lboun2+ndim2-1,lboun3:lboun3+ndim3-1) , stat = istat )
       else
          allocate( varia(ndim1,ndim2,ndim3) , stat = istat )
       end if
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*size_i1p_type
          if( lzero ) then
             do idim3 = lbound(varia,3),ubound(varia,3)
                do idim2 = lbound(varia,2),ubound(varia,2)
                   do idim1 = lbound(varia,1),ubound(varia,1)
                      nullify(varia(idim1,idim2,idim3) % l)
                   end do
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type(i1p)')

    else

       nullify(varia)

    end if

  end subroutine memory_alli1p_3

  subroutine memory_alli1p_1(memor,vanam,vacal,varia,ndim1,wzero,lboun)
    !
    ! Type(i1p)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    type(i1p),    intent(inout), pointer  :: varia(:)
    character(*), intent(in),  optional :: wzero
    integer(ip),  intent(in),  optional :: lboun
    integer(ip)                         :: idim1
    integer(4)                          :: istat
    logical(lg)                         :: lzero

    if( ndim1 > 0 ) then

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       if( present(lboun) ) then
          allocate( varia(lboun:lboun+ndim1-1) , stat = istat )
       else
          allocate( varia(ndim1) , stat = istat )
       end if
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*size_i1p_type
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim1 = lbound(varia,1),ubound(varia,1)
                nullify(varia(idim1) % l)
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type(i1p)')

    else

       nullify(varia)

    end if

  end subroutine memory_alli1p_1

  subroutine memory_alli1pp_1(memor,vanam,vacal,varia,ndim1,wzero)
    !
    ! Type(i1p)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    type(i1pp),   intent(inout), pointer  :: varia(:)
    character(*), intent(in),  optional :: wzero
    integer(ip)                         :: idim1
    integer(4)                          :: istat
    logical(lg)                         :: lzero

    if( ndim1 > 0 ) then

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       allocate( varia(ndim1) , stat = istat )
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*size_i1pp_type
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim1 = 1,ndim1
                varia(idim1) % n = 0
                nullify(varia(idim1) % l)
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type(i1p)')

    else

       nullify(varia)

    end if

  end subroutine memory_alli1pp_1

  subroutine memory_alli2p_1(memor,vanam,vacal,varia,ndim1,wzero,lboun)
    !
    ! Type(i2p)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    type(i2p),    intent(inout), pointer  :: varia(:)
    character(*), intent(in),  optional :: wzero
    integer(ip),  intent(in),  optional :: lboun
    integer(ip)                         :: idim1
    integer(4)                          :: istat
    logical(lg)                         :: lzero

    if( ndim1 > 0 ) then

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( present(lboun) ) then
          allocate( varia(lboun:lboun+ndim1-1) , stat = istat )
       else
          allocate( varia(ndim1) , stat = istat )
       end if

      if( istat == 0 ) then
          lbytm = size(varia,kind=8)*size_i2p_type
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim1 = lbound(varia,1),ubound(varia,1)
                nullify(varia(idim1) % l)
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type(i2p)')

    else

       nullify(varia)

    end if

  end subroutine memory_alli2p_1

  subroutine memory_alli2p_2(memor,vanam,vacal,varia,ndim1,ndim2,wzero)
    !
    ! Type(i2p)(:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    integer(ip),  intent(in)            :: ndim2
    type(i2p),    intent(inout), pointer  :: varia(:,:)
    character(*), intent(in),  optional :: wzero
    integer(ip)                         :: idim1
    integer(ip)                         :: idim2
    integer(4)                          :: istat
    logical(lg)                         :: lzero

    if( ndim1*ndim2 > 0 ) then

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       allocate( varia(ndim1,ndim2) , stat = istat )
       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*size_i2p_type
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim1 = 1,ndim1
                do idim2 = 1,ndim2
                   nullify(varia(idim1,idim2) % l)
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type(i2p)')

    else

       nullify(varia)

    end if

  end subroutine memory_alli2p_2

  subroutine memory_alli3p_1(memor,vanam,vacal,varia,ndim1,wzero)
    !
    ! Type(i3p)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    type(i3p),    intent(inout), pointer  :: varia(:)
    character(*), intent(in),  optional :: wzero
    integer(ip)                         :: idim1
    integer(4)                          :: istat
    logical(lg)                         :: lzero

    if( ndim1 > 0 ) then

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       allocate( varia(ndim1) , stat = istat )
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*size_i3p_type
          !lbytm = int(sizeof(varia),8) ! ndim1*ip
          if( lzero ) then
             do idim1 = 1,ndim1
                nullify(varia(idim1) % l)
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type(i3p)')

    else

       nullify(varia)

    end if

  end subroutine memory_alli3p_1

  subroutine memory_alli1pi1p(memor,vanam,vacal,varia,ndim1,wzero)
    !
    ! Type(i1p)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    type(i1pi1p), intent(inout), pointer  :: varia(:)
    character(*), intent(in),  optional :: wzero
    integer(ip)                         :: idim1
    integer(4)                          :: istat
    logical(lg)                         :: lzero

    if( ndim1 > 0 ) then

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       allocate( varia(ndim1) , stat = istat )
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*size_i1pi1p_type
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim1 = 1,ndim1
                nullify(varia(idim1) % l)
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'type(i1p)')

    else

       nullify(varia)

    end if

  end subroutine memory_alli1pi1p

  subroutine memory_deai41(memor,vanam,vacal,varia)
    !
    ! Integer(4)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(4),   pointer               :: varia(:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8)

       deallocate( varia , stat = istat )

       nullify(varia)

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)

       end if
       call memory_info(memor,vanam,vacal,'integer(4)')

    end if

  end subroutine memory_deai41

  subroutine memory_deai42(memor,vanam,vacal,varia)
    !
    ! Integer(4)(:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(4),   pointer               :: varia(:,:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'integer(4)')

    end if

  end subroutine memory_deai42

  subroutine memory_deai43(memor,vanam,vacal,varia)
    !
    ! Integer(4)(:,:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(4),   pointer               :: varia(:,:,:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8) !-size(varia,kind=8)*int(kind(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'integer(4)')

    end if

  end subroutine memory_deai43

  subroutine memory_dear81(memor,vanam,vacal,varia)
    !
    ! Real(rp)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    real(8),     pointer                :: varia(:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8) ! -size(varia,kind=8)*int(kind(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'real')

    end if

  end subroutine memory_dear81

  subroutine memory_dear82(memor,vanam,vacal,varia)
    !
    ! Real(rp)(:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    real(8),     pointer                :: varia(:,:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8) ! -size(varia,kind=8)*int(kind(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'real')

    end if

  end subroutine memory_dear82

  subroutine memory_dear83(memor,vanam,vacal,varia)
    !
    ! Real(rp)(:,:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    real(8),     pointer                :: varia(:,:,:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'real')
    end if

  end subroutine memory_dear83

  subroutine memory_dear84(memor,vanam,vacal,varia)
    !
    ! Real(rp)(:,:,:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    real(8),     pointer                :: varia(:,:,:,:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8) ! -size(varia,kind=8)*int(kind(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'real')
    end if

  end subroutine memory_dear84

    subroutine memory_dear41(memor,vanam,vacal,varia)
    !
    ! Real(rp)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    real(4),     pointer                :: varia(:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8) ! -size(varia,kind=8)*int(kind(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'real')

    end if

  end subroutine memory_dear41

  subroutine memory_dear42(memor,vanam,vacal,varia)
    !
    ! Real(rp)(:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    real(4),     pointer                :: varia(:,:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8) ! -size(varia,kind=8)*int(kind(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'real')

    end if

  end subroutine memory_dear42

  subroutine memory_dear43(memor,vanam,vacal,varia)
    !
    ! Real(rp)(:,:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    real(4),     pointer                :: varia(:,:,:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'real')
    end if

  end subroutine memory_dear43

  subroutine memory_dear44(memor,vanam,vacal,varia)
    !
    ! Real(rp)(:,:,:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    real(4),     pointer                :: varia(:,:,:,:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8) ! -size(varia,kind=8)*int(kind(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'real')
    end if

  end subroutine memory_dear44

#ifdef __PGI

  subroutine memory_dear161(dummy)
    implicit none
    type(pgi_dummy1), intent(in)          :: dummy
  end subroutine memory_dear161

  subroutine memory_dear162(dummy)
    implicit none
    type(pgi_dummy2), intent(in)          :: dummy
  end subroutine memory_dear162

  subroutine memory_dear163(dummy)
    implicit none
    type(pgi_dummy3), intent(in)          :: dummy
  end subroutine memory_dear163

  subroutine memory_dear164(dummy)
    implicit none
    type(pgi_dummy4), intent(in)          :: dummy
  end subroutine memory_dear164

#else

  subroutine memory_dear161(memor,vanam,vacal,varia)
    !
    ! Real(rp)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    real(16),     pointer               :: varia(:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8) ! -size(varia,kind=8)*int(kind(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'real')

    end if

  end subroutine memory_dear161

  subroutine memory_dear162(memor,vanam,vacal,varia)
    !
    ! Real(rp)(:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    real(16),     pointer               :: varia(:,:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8) ! -size(varia,kind=8)*int(kind(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'real')

    end if

  end subroutine memory_dear162

  subroutine memory_dear163(memor,vanam,vacal,varia)
    !
    ! Real(rp)(:,:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    real(16),     pointer               :: varia(:,:,:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'real')
    end if

  end subroutine memory_dear163

  subroutine memory_dear164(memor,vanam,vacal,varia)
    !
    ! Real(rp)(:,:,:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    real(16),     pointer               :: varia(:,:,:,:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8) ! -size(varia,kind=8)*int(kind(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'real')
    end if

  end subroutine memory_dear164

#endif

  subroutine memory_deaspmat_1(memor,vanam,vacal,varia)
     !
     ! type(spmat)
     !
     implicit none
     character(*), intent(in)            :: vanam         !< Variable name
     character(*), intent(in)            :: vacal         !< Calling subroutine name
     integer(8),   intent(inout)         :: memor(2)      !< Memory counter
     type(spmat), pointer                :: varia(:)
     integer(4)                          :: istat
     integer(ip)                         :: idim1
     if( associated(varia) ) then

        do idim1 = lbound(varia,1),ubound(varia,1)
           call memory_deaspmat( memor, vanam, vacal, varia(idim1) )
        end do
        lbytm = -size(varia,kind=8)*ip
        deallocate( varia , stat = istat )

        if( istat /= 0 ) then
           call memory_error(2_ip,vanam,vacal,istat)
        else
           nullify (varia)
        end if
        call memory_info(memor,vanam,vacal,'real')
     end if

  end subroutine memory_deaspmat_1

  subroutine memory_deaspmat(memor,vanam,vacal,varia)
     !
     ! type(spmat)
     !
     implicit none
     character(*), intent(in)            :: vanam         !< Variable name
     character(*), intent(in)            :: vacal         !< Calling subroutine name
     integer(8),   intent(inout)         :: memor(2)      !< Memory counter
     type(spmat)                         :: varia
     integer(4)                          :: istat
     integer(4)                          :: istat_aux

     istat = 0
     varia % ndof  = 0
     varia % nrows = 0
     varia % ncols = 0

     if( associated(varia % iA) .or.  associated(varia % jA) .or. associated(varia % vA) ) then

        lbytm = - size(varia % iA,kind=8) * (ip * 2 + rp * size(varia % vA,2,kind=8)  * size(varia % vA,2,kind=8))
        if(associated(varia % iA) ) then
           deallocate(varia % iA, stat=istat_aux)
           istat = istat + istat_aux
        endif
        if(associated(varia % jA) ) then
           deallocate(varia % jA, stat=istat_aux)
           istat = istat + istat_aux
        endif
        if(associated(varia % vA) ) then
           deallocate(varia % vA, stat=istat_aux)
           istat = istat + istat_aux
        endif

        if( istat /= 0 ) then
           call memory_error(2_ip,vanam,vacal,istat)
        else
           nullify (varia % iA)
           nullify (varia % jA)
           nullify (varia % vA)
        end if
        call memory_info(memor,vanam,vacal,'real')
     end if

  end subroutine memory_deaspmat

  subroutine memory_dealg1(memor,vanam,vacal,varia)
    !
    ! Logical(lg)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    logical(lg),  pointer               :: varia(:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8) !-size(varia,kind=8)*int(kind(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'real')

    end if

  end subroutine memory_dealg1

  subroutine memory_dealg2(memor,vanam,vacal,varia)
    !
    ! Logical(lg)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    logical(lg),  pointer               :: varia(:,:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8) ! -size(varia,kind=8)*int(kind(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'real')

    end if

  end subroutine memory_dealg2

  subroutine memory_dear1p_1(memor,vanam,vacal,varia)
    !
    ! Type(r1p)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(r1p),    pointer               :: varia(:)
    integer(4)                          :: istat
    integer(ip)                         :: idim1

    if( associated(varia) ) then

       lbytm = 0
       do idim1 = lbound(varia,1),ubound(varia,1)
          if( associated(varia(idim1) % a) ) then
             lbytm = lbytm-size(varia(idim1) % a,1,kind=8)*int(kind(varia(idim1) % a),8)
             !lbytm = lbytm-int(sizeof(varia(idim1) % a),8)
             deallocate( varia(idim1) % a )
          end if
       end do
       lbytm = lbytm-size(varia,kind=8)*size_r1p_type
       !lbytm = lbytm-int(sizeof(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type(r1p)')
    end if

  end subroutine memory_dear1p_1

  subroutine memory_dear2p_2(memor,vanam,vacal,varia)
    !
    ! Type(r2p)(:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(r2p),    pointer               :: varia(:,:)
    integer(4)                          :: istat
    integer(ip)                         :: idim1,idim2

    if( associated(varia) ) then

       lbytm = 0
       do idim1 = lbound(varia,1),ubound(varia,1)
          do idim2 = lbound(varia,2),ubound(varia,2)
             if( associated(varia(idim1,idim2) % a) ) then
                lbytm = lbytm - int(size(varia(idim1,idim2) % a),8)*int(kind(varia(idim1,idim2) % a),8)
                !lbytm = lbytm - int(sizeof(varia(idim1,idim2) % a),8)
                deallocate( varia(idim1,idim2) % a )
             end if
          end do
       end do
       lbytm = lbytm-size(varia,kind=8)*size_r2p_type
       !lbytm = lbytm - int(sizeof(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type(r2p)')
    end if

  end subroutine memory_dear2p_2

  subroutine memory_dear2p_1(memor,vanam,vacal,varia)
    !
    ! Type(r2p)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(r2p),    pointer               :: varia(:)
    integer(4)                          :: istat
    integer(ip)                         :: idim1

    if( associated(varia) ) then

       lbytm = 0
       do idim1 = lbound(varia,1),ubound(varia,1)
          if( associated(varia(idim1) % a) ) then
             lbytm = lbytm - size(varia(idim1) % a,kind=8)*int(kind(varia(idim1) % a),8)
             !lbytm = lbytm - int(sizeof(varia(idim1) % a),8)
             deallocate( varia(idim1) % a )
          end if
       end do

       lbytm = lbytm-size(varia,kind=8)*size_r2p_type
       !lbytm = lbytm-int(sizeof(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type(r2p)')
    end if

  end subroutine memory_dear2p_1

  subroutine memory_dear3p_1(memor,vanam,vacal,varia)
    !
    ! Type(r3p)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(r3p),    pointer               :: varia(:)
    integer(4)                          :: istat
    integer(ip)                         :: idim1

    if( associated(varia) ) then

       lbytm = 0
       do idim1 = lbound(varia,1),ubound(varia,1)
          if( associated(varia(idim1) % a) ) then
             lbytm = lbytm - int(size(varia(idim1) % a),8)*int(kind(varia(idim1) % a),8)
             !lbytm = lbytm - int(sizeof(varia(idim1) % a),8)
             deallocate( varia(idim1) % a )
          end if
       end do
       lbytm = lbytm-size(varia,kind=8)*size_r3p_type
       !lbytm = lbytm-int(sizeof(varia),8) !int(size(varia,1),8)*ip
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type(r3p)')
    end if

  end subroutine memory_dear3p_1

  subroutine memory_dear4p_1(memor,vanam,vacal,varia)
    !
    ! Type(r3p)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(r4p),    pointer               :: varia(:)
    integer(4)                          :: istat
    integer(ip)                         :: idim1

    if( associated(varia) ) then

       lbytm = 0
       do idim1 = lbound(varia,1),ubound(varia,1)
          if( associated(varia(idim1) % a) ) then
             lbytm = lbytm - size(varia(idim1) % a,kind=8)*int(kind(varia(idim1) % a),8)
             !lbytm = lbytm - int(sizeof(varia(idim1) % a),8)
             deallocate( varia(idim1) % a )
          end if
       end do
       lbytm = lbytm-size(varia,kind=8)*size_r4p_type
       !lbytm = lbytm-int(sizeof(varia),8)  !-int(size(varia,1),8)*ip
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type(r4p)')
    end if

  end subroutine memory_dear4p_1

  subroutine memory_dear3p_2(memor,vanam,vacal,varia)
    !
    ! Type(r3p)(:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(r3p),    pointer               :: varia(:,:)
    integer(4)                          :: istat
    integer(ip)                         :: idim1
    integer(ip)                         :: idim2

    if( associated(varia) ) then

       lbytm = 0
       do idim2 = lbound(varia,2),ubound(varia,2)
          do idim1 = lbound(varia,1),ubound(varia,1)
             if( associated(varia(idim1,idim2) % a) ) then
                lbytm = lbytm - size(varia(idim1,idim2) % a,kind=8)*int(kind(varia(idim1,idim2) % a),8)
                !lbytm = lbytm - int(sizeof(varia(idim1,idim2) % a),8)
                deallocate( varia(idim1,idim2) % a )
             end if
          end do
       end do

       lbytm = lbytm-size(varia,kind=8)*size_r3p_type
       !lbytm = lbytm-int(sizeof(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type(r3p)')
    end if

  end subroutine memory_dear3p_2

  subroutine memory_dear3p_3(memor,vanam,vacal,varia)
    !
    ! Type(r3p)(:,:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(r3p),    pointer               :: varia(:,:,:)
    integer(4)                          :: istat
    integer(ip)                         :: idim1
    integer(ip)                         :: idim2
    integer(ip)                         :: idim3

    if( associated(varia) ) then

       lbytm = 0
       do idim3 = lbound(varia,3),ubound(varia,3)
          do idim2 = lbound(varia,2),ubound(varia,2)
             do idim1 = lbound(varia,1),ubound(varia,1)
                if( associated(varia(idim1,idim2,idim3) % a) ) then
                   lbytm = lbytm - size(varia(idim1,idim2,idim3) % a,kind=8)*int(kind(varia(idim1,idim2,idim3) % a),8)
                   !lbytm = lbytm - int(sizeof(varia(idim1,idim2,idim3) % a),8)
                   deallocate( varia(idim1,idim2,idim3) % a )
                end if
             end do
          end do
       end do

       lbytm = lbytm-size(varia,kind=8)*size_r3p_type
       !lbytm = lbytm-int(sizeof(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type(r3p)')
    end if

  end subroutine memory_dear3p_3

  subroutine memory_deai1p_1(memor,vanam,vacal,varia)
    !
    ! Type(i1p)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(i1p),    pointer               :: varia(:)
    integer(4)                          :: istat
    integer(ip)                         :: idim1

    if( associated(varia) ) then

       lbytm = 0
       do idim1 = lbound(varia,1),ubound(varia,1)
          if( associated(varia(idim1) % l) ) then
             lbytm = lbytm-size(varia(idim1) % l,1,kind=8)*int(kind(varia(idim1) % l),8)
             !lbytm = lbytm-int(sizeof(varia(idim1) % l),8)
             deallocate( varia(idim1) % l )
          end if
       end do
       lbytm = lbytm-size(varia,kind=8)*size_i1p_type
       !lbytm = lbytm-int(sizeof(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type(i1p)')

    end if

  end subroutine memory_deai1p_1

  subroutine memory_deai1p_2(memor,vanam,vacal,varia)
    !
    ! Type(i1p)(:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(i1p),    pointer               :: varia(:,:)
    integer(4)                          :: istat
    integer(ip)                         :: idim1
    integer(ip)                         :: idim2

    if( associated(varia) ) then

       lbytm = 0
       do idim2 = lbound(varia,2),ubound(varia,2)
          do idim1 = lbound(varia,1),ubound(varia,1)
             if( associated(varia(idim1,idim2) % l) ) then
                lbytm = lbytm-size(varia(idim1,idim2) % l,kind=8)*int(kind(varia(idim1,idim2) % l),8)
                !lbytm = lbytm-int(sizeof(varia(idim1,idim2) % l),8)
                deallocate( varia(idim1,idim2) % l )
             end if
          end do
       end do

       lbytm = lbytm-size(varia,kind=8)*size_i1p_type
       !lbytm = lbytm-int(sizeof(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type(i1p)')

    end if

  end subroutine memory_deai1p_2

  subroutine memory_deai1p_3(memor,vanam,vacal,varia)
    !
    ! Type(i1p)(:,:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(i1p),    pointer               :: varia(:,:,:)
    integer(4)                          :: istat
    integer(ip)                         :: idim1
    integer(ip)                         :: idim2
    integer(ip)                         :: idim3

    if( associated(varia) ) then

       lbytm = 0
       do idim3 = lbound(varia,3),ubound(varia,3)
          do idim2 = lbound(varia,2),ubound(varia,2)
             do idim1 = lbound(varia,1),ubound(varia,1)
                if( associated(varia(idim1,idim2,idim3) % l) ) then
                   lbytm = lbytm-size(varia(idim1,idim2,idim3) % l,kind=8)*int(kind(varia(idim1,idim2,idim3) % l),8)
                   !lbytm = lbytm-int(sizeof(varia(idim1,idim2,idim3) % l),8)
                   deallocate( varia(idim1,idim2,idim3) % l )
                end if
             end do
          end do
       end do
       lbytm = lbytm-size(varia,kind=8)*size_i1p_type
       !lbytm = lbytm-int(sizeof(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type(i1p)')

    end if

  end subroutine memory_deai1p_3

  subroutine memory_deai1pp_1(memor,vanam,vacal,varia)
    !
    ! Type(i1pp)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(i1pp),   pointer               :: varia(:)
    integer(4)                          :: istat
    integer(ip)                         :: idim1

    if( associated(varia) ) then

       lbytm = 0
       do idim1 = lbound(varia,1),ubound(varia,1)
          if( associated(varia(idim1) % l) ) then
             lbytm = lbytm-size(varia(idim1) % l,kind=8)*int(kind(varia(idim1) % l),8)
             !lbytm = lbytm-int(sizeof(varia(idim1) % l),8)
             deallocate( varia(idim1) % l , stat = istat )
          end if
       end do
       lbytm = lbytm-size(varia,kind=8)*size_i1pp_type
       !lbytm = lbytm-int(sizeof(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type(i1p)')

    end if

  end subroutine memory_deai1pp_1

  subroutine memory_deai2p_2(memor,vanam,vacal,varia)
    !
    ! Type(i2p)(:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(i2p),    pointer               :: varia(:,:)
    integer(4)                          :: istat
    integer(ip)                         :: idim1,idim2

    if( associated(varia) ) then

       lbytm = 0
       do idim1 = 1,size(varia,1,kind=ip)
          do idim2 = 1,size(varia,2,kind=ip)
             if( associated(varia(idim1,idim2) % l) ) then
                lbytm = lbytm-size(varia(idim1,idim2) % l,kind=8)*int(kind(varia(idim1,idim2) % l),8)
                !lbytm = lbytm-int(sizeof(varia(idim1,idim2) % l),8)
                deallocate( varia(idim1,idim2) % l )
             end if
          end do
       end do

       lbytm = lbytm-size(varia,kind=8)*size_i2p_type
       !lbytm = lbytm-int(sizeof(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type(i2p)')

    end if

  end subroutine memory_deai2p_2

  subroutine memory_deai2p_1(memor,vanam,vacal,varia)
    !
    ! Type(i2p)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(i2p),    pointer               :: varia(:)
    integer(4)                          :: istat
    integer(ip)                         :: idim1

    if( associated(varia) ) then

       lbytm = 0
       do idim1 = lbound(varia,1),ubound(varia,1)
          if( associated(varia(idim1) % l) ) then
             lbytm = lbytm-size(varia(idim1) % l,kind=8)*int(kind(varia(idim1) % l),8)
             !lbytm = lbytm-int(sizeof(varia(idim1) % l),8)
             deallocate( varia(idim1) % l )
          end if
       end do

       lbytm = lbytm-size(varia,kind=8)*size_i2p_type
       !lbytm = lbytm-int(sizeof(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type(i2p)')

    end if

  end subroutine memory_deai2p_1

  subroutine memory_deai3p_1(memor,vanam,vacal,varia)
    !
    ! Type(i3p)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(i3p),    pointer               :: varia(:,:,:)
    integer(4)                          :: istat
    integer(ip)                         :: idim1
    integer(ip)                         :: idim2
    integer(ip)                         :: idim3

    if( associated(varia) ) then

       lbytm = 0
       do idim3 = lbound(varia,3),ubound(varia,3)
          do idim2 = lbound(varia,2),ubound(varia,2)
             do idim1 = lbound(varia,1),ubound(varia,1)
                if( associated(varia(idim1,idim2,idim3) % l) ) then
                   lbytm = lbytm-size(varia(idim1,idim2,idim3) % l,kind=8)*int(kind(varia(idim1,idim2,idim3) % l),8)
                   !lbytm = lbytm-int(sizeof(varia(idim1,idim2,idim3) % l),8)
                   deallocate( varia(idim1,idim2,idim3) % l )
                end if
             end do
          end do
       end do

       lbytm = lbytm-size(varia,kind=8)*size_i3p_type
       !lbytm = lbytm-int(sizeof(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'type(i3p)')

    end if

  end subroutine memory_deai3p_1

  subroutine memory_deaxp1(memor,vanam,vacal,varia)
    !
    ! complex(rp)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    complex(rp),     pointer            :: varia(:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)*2
       !lbytm = -int(sizeof(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'real')

    end if

  end subroutine memory_deaxp1

  subroutine memory_deaxp2(memor,vanam,vacal,varia)
    !
    ! complex(rp)(:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    complex(rp),  pointer               :: varia(:,:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)*2
       !lbytm = -int(sizeof(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'real')

    end if

  end subroutine memory_deaxp2

  subroutine memory_alli81(memor,vanam,vacal,varia,ndim1,wzero,lboun)
    !
    ! Integer(8)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(8),   intent(in)            :: ndim1
    integer(8),   intent(inout), pointer  :: varia(:)
    character(*), intent(in),  optional :: wzero
    integer(ip),  intent(in),  optional :: lboun
    integer(4)                          :: istat
    integer(8)                          :: idim1
    logical(lg)                         :: lzero
    logical(lg)                         :: linde
    integer(8)                          :: xvalu

    if( ndim1 > 0 ) then

       if( present(wzero) ) then
          if( wzero == 'REALLOCATE') then
             call memory_deallo(memor,vanam,vacal,varia)
          end if
       end if

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       if( present(lboun) ) then
          allocate( varia(lboun:lboun+int(ndim1,KIND=ip)-1_ip) , stat = istat )
       else
          allocate( varia(ndim1) , stat = istat )
       end if
       lzero = .true.
       linde = .false.
       xvalu = 0_8
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          else if( wzero == 'HUGE') then
             xvalu = huge(0_8)
          else if( wzero == 'IDENTITY') then
             lzero = .false.
             linde = .true.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8) !int(ndim1,8)*int(kind(varia),8)
          if( lzero ) then
             do idim1 = lbound(varia,1,KIND=8),ubound(varia,1,KIND=8)
                varia(idim1) = xvalu
             end do
          else if( linde ) then
             do idim1 = lbound(varia,1,KIND=8),ubound(varia,1,KIND=8)
                varia(idim1) = idim1
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'integer(4)')

    else

       nullify(varia)

    end if

  end subroutine memory_alli81

  subroutine memory_alli82(memor,vanam,vacal,varia,ndim1,ndim2,wzero,lboun1,lboun2)
    !
    ! Integer(8)(:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(8),   intent(in)            :: ndim1
    integer(8),   intent(in)            :: ndim2
    integer(8),   intent(inout), pointer  :: varia(:,:)
    character(*), intent(in),  optional :: wzero
    integer(ip),  intent(in),  optional :: lboun1
    integer(ip),  intent(in),  optional :: lboun2
    integer(4)                          :: istat
    integer(ip)                          :: idim1
    integer(ip)                          :: idim2
    logical(lg)                         :: lzero

    if( ndim1*ndim2 > 0 ) then

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       if( present(lboun1) .and. present(lboun2) ) then
          allocate( varia(lboun1:lboun1+int(ndim1,KIND=ip)-1_ip,lboun2:lboun2+int(ndim2,KIND=ip)-1_ip) , stat = istat )
       else
          allocate( varia(ndim1,ndim2) , stat = istat )
       end if
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim2 = lbound(varia,2,KIND=ip),ubound(varia,2,KIND=ip)
                do idim1 = lbound(varia,1,KIND=ip),ubound(varia,1,KIND=ip)
                   varia(idim1,idim2) = 0_ip
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'integer(4)')

    else

       nullify(varia)

    end if

  end subroutine memory_alli82

  subroutine memory_alli83(memor,vanam,vacal,varia,ndim1,ndim2,ndim3,wzero,lboun1,lboun2,lboun3)
    !
    ! Integer(8)(:,:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(8),   intent(in)            :: ndim1
    integer(8),   intent(in)            :: ndim2
    integer(8),   intent(in)            :: ndim3
    integer(8),   intent(inout), pointer  :: varia(:,:,:)
    character(*), intent(in),  optional :: wzero
    integer(8),   intent(in),  optional :: lboun1
    integer(8),   intent(in),  optional :: lboun2
    integer(8),   intent(in),  optional :: lboun3
    integer(4)                          :: istat
    integer(ip)                         :: idim1
    integer(ip)                         :: idim2
    integer(ip)                         :: idim3
    logical(lg)                         :: lzero

    if( ndim1*ndim2*ndim3 > 0 ) then

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       if( present(lboun1) .and. present(lboun2) .and. present(lboun3) ) then
          allocate( varia(lboun1:lboun1+ndim1-1,lboun2:lboun2+ndim2-1,lboun3:lboun3+ndim3-1) , stat = istat )
       else
          allocate( varia(ndim1,ndim2,ndim3) , stat = istat )
       end if
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = size(varia,kind=8)*int(kind(varia),8)
          !lbytm = int(sizeof(varia),8)
          if( lzero ) then
             do idim3 = lbound(varia,3,KIND=ip),ubound(varia,3,KIND=ip)
                do idim2 = lbound(varia,2,KIND=ip),ubound(varia,2,KIND=ip)
                   do idim1 = lbound(varia,1,KIND=ip),ubound(varia,1,KIND=ip)
                      varia(idim1,idim2,idim3) = 0_8
                   end do
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'integer(4)')

    else

       nullify(varia)

    end if

  end subroutine memory_alli83

  subroutine memory_deai81(memor,vanam,vacal,varia)
    !
    ! Integer(8)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(8),   pointer               :: varia(:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'integer(4)')

    end if

  end subroutine memory_deai81

  subroutine memory_deai82(memor,vanam,vacal,varia)
    !
    ! Integer(8)(:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(8),   pointer               :: varia(:,:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'integer(4)')

    end if

  end subroutine memory_deai82

  subroutine memory_deai83(memor,vanam,vacal,varia)
    !
    ! Integer(8)(:,:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(8),   pointer               :: varia(:,:,:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'integer(4)')

    end if

  end subroutine memory_deai83

  subroutine memory_allcha_0(memor,vanam,vacal,varia,ndim1,wzero)
    !
    ! Character(:)
    !
    integer(8),       intent(inout)           :: memor(2)      !< Memory counter
    character(*),     intent(in)              :: vanam         !< Variable name
    character(*),     intent(in)              :: vacal         !< Calling subroutine name
    integer(ip),      intent(in)              :: ndim1
    character(len=:), intent(inout), pointer  :: varia
    character(*),     intent(in),    optional :: wzero
    integer(4)                                :: istat
    integer(ip)                               :: idim1
    logical(lg)                               :: lzero

    if( ndim1 > 0 ) then

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       allocate( character(ndim1) :: varia , stat = istat )
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = -int(ndim1,kind=8)*int(kind(varia),8)
          if( lzero ) then
             varia = ''
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'character')

    else

       nullify(varia)

    end if

  end subroutine memory_allcha_0
  
  subroutine memory_allcha_4(memor,vanam,vacal,varia,ndim1,wzero)
    !
    ! Character(*)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(4),   intent(in)            :: ndim1
    character(*), intent(inout), pointer  :: varia(:)
    character(*), intent(in),  optional :: wzero
    integer(4)                          :: istat
    integer(4)                          :: idim1
    logical(lg)                         :: lzero

    if( ndim1 > 0 ) then

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       allocate( varia(ndim1) , stat = istat )
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = -size(varia,kind=8)*int(kind(varia),8)
          !lbytm = -int(sizeof(varia),8)
          if( lzero ) then
             do idim1 = 1,ndim1
                varia(idim1) = ''
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'character')

    else

       nullify(varia)

    end if

  end subroutine memory_allcha_4

  subroutine memory_allcha_8(memor,vanam,vacal,varia,ndim1,wzero)
    !
    ! Character(*)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(8),   intent(in)            :: ndim1
    character(*), intent(inout), pointer  :: varia(:)
    character(*), intent(in),  optional :: wzero
    integer(4)                          :: istat
    integer(8)                          :: idim1
    logical(lg)                         :: lzero

    if( ndim1 > 0 ) then

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       allocate( varia(ndim1) , stat = istat )
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = -size(varia,kind=8)*int(kind(varia),8)
          !lbytm = -int(sizeof(varia),8)
          if( lzero ) then
             do idim1 = 1,ndim1
                varia(idim1) = ''
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'character')

    else

       nullify(varia)

    end if

  end subroutine memory_allcha_8

  subroutine memory_deacha(memor,vanam,vacal,varia)
    !
    ! Character(*)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    character(*), pointer               :: varia(:)
    integer(4)                          :: istat

    if( associated(varia) ) then

       lbytm = -size(varia,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'character')

    end if

  end subroutine memory_deacha

  subroutine memory_deacha_0(memor,vanam,vacal,varia)
    !
    ! Character(*)(:)
    !
    implicit none
    character(*),     intent(in)            :: vanam         !< Variable name
    character(*),     intent(in)            :: vacal         !< Calling subroutine name
    integer(8),       intent(inout)         :: memor(2)      !< Memory counter
    character(len=:), pointer               :: varia
    integer(4)                              :: istat
    integer(ip)                             :: ndim1
    
    if( associated(varia) ) then

       ndim1 = len(varia)
       lbytm = -int(ndim1,kind=8)*int(kind(varia),8)
       !lbytm = -int(sizeof(varia),8)
       deallocate( varia , stat = istat )

       if( istat /= 0 ) then
          call memory_error(2_ip,vanam,vacal,istat)
       else
          nullify (varia)
       end if
       call memory_info(memor,vanam,vacal,'character')

    end if

  end subroutine memory_deacha_0

  !----------------------------------------------------------------------
  !
  ! Copy arrays
  !
  ! 1. Allocate varia if it is null
  ! 2. varia <= vacpy
  ! 3. Deallocate vacpy if 'DO_NOT_DEALLOCATE' not present
  !
  !----------------------------------------------------------------------

  subroutine memory_cpyi41(memor,vanam,vacal,vacpy,varia,wzero)
    !
    ! Integer(4)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(4),   intent(inout), pointer  :: varia(:)
    integer(4),   intent(inout), pointer  :: vacpy(:)
    character(*), intent(in),  optional :: wzero
    integer(ip)                         :: idim1
    logical(lg)                         :: lzero
    integer(ip)                         :: ndim1

    if( .not. associated(vacpy) ) then
       ndim1 = 0
    else
       ndim1 = size(vacpy,1_4,kind=ip)
    end if

    if( ndim1 > 0 ) then

       if( .not. associated(varia) ) then
          call memory_alloca(memor,vanam,vacal,varia,ndim1,'DO_NOT_INITIALIZE')
       else
          ndim1 = min(ndim1,size(varia,kind=ip))
       end if

       do idim1 = 1,ndim1
          varia(idim1) = vacpy(idim1)
       end do

    else

       nullify(varia)

    end if

    lzero = .true.
    if(present(wzero)) then
       if( wzero == 'DO_NOT_DEALLOCATE') lzero = .false.
    end if
    if( lzero ) call memory_deallo(memor,vanam,vacal,vacpy)

  end subroutine memory_cpyi41

  subroutine memory_cpyi81(memor,vanam,vacal,vacpy,varia,wzero)
    !
    ! Integer(8)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(8),   intent(inout), pointer  :: varia(:)
    integer(8),   intent(inout), pointer  :: vacpy(:)
    character(*), intent(in),  optional :: wzero
    integer(8)                          :: idim1
    logical(lg)                         :: lzero
    integer(8)                          :: ndim1

    if( .not. associated(vacpy) ) then
       ndim1 = 0
    else
       ndim1 = size(vacpy,kind=8)
    end if

    if( ndim1 > 0 ) then

       if( .not. associated(varia) ) then
          call memory_alloca(memor,vanam,vacal,varia,ndim1,'DO_NOT_INITIALIZE')
       else
          ndim1 = min(ndim1,size(varia,kind=8))
       end if

       do idim1 = 1,ndim1
          varia(idim1) = vacpy(idim1)
       end do

    else

       nullify(varia)

    end if

    lzero = .true.
    if(present(wzero)) then
       if( wzero == 'DO_NOT_DEALLOCATE') lzero = .false.
    end if
    if( lzero ) call memory_deallo(memor,vanam,vacal,vacpy)

  end subroutine memory_cpyi81

  subroutine memory_cpyi42(memor,vanam,vacal,vacpy,varia,wzero)
    !
    ! Integer(4)(:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(4),   intent(inout), pointer  :: varia(:,:)
    integer(4),   intent(inout), pointer  :: vacpy(:,:)
    character(*), intent(in),  optional :: wzero
    integer(4)                          :: idim1,idim2
    logical(lg)                         :: lzero
    integer(4)                          :: ndim1,ndim2

    if( .not. associated(vacpy) ) then
       ndim1 = 0
       ndim2 = 0
    else
       ndim1 = size(vacpy,1_ip,kind=4)
       ndim2 = size(vacpy,2_ip,kind=4)
    end if

    if( ndim1 > 0 .and. ndim2 > 0 ) then

       if( .not. associated(varia) ) then
          call memory_alloca(memor,vanam,vacal,varia,ndim1,ndim2,'DO_NOT_INITIALIZE')
       else
          ndim1 = min(ndim1,size(varia,1,kind=4))
          ndim2 = min(ndim2,size(varia,2,kind=4))
       end if

       do idim2 = 1,ndim2
          do idim1 = 1,ndim1
             varia(idim1,idim2) = vacpy(idim1,idim2)
          end do
       end do

    else

       nullify(varia)

    end if

    lzero = .true.
    if(present(wzero)) then
       if( wzero == 'DO_NOT_DEALLOCATE') lzero = .false.
    end if
    if( lzero ) call memory_deallo(memor,vanam,vacal,vacpy)

  end subroutine memory_cpyi42

  subroutine memory_cpyi43(memor,vanam,vacal,vacpy,varia,wzero)
    !
    ! Integer(4)(:,:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(4),   intent(inout), pointer  :: varia(:,:,:)
    integer(4),   intent(inout), pointer  :: vacpy(:,:,:)
    character(*), intent(in),  optional :: wzero
    integer(4)                          :: idim1,idim2,idim3
    logical(lg)                         :: lzero
    integer(4)                          :: ndim1,ndim2,ndim3

    if( .not. associated(vacpy) ) then
       ndim1 = 0
       ndim2 = 0
       ndim3 = 0
    else
       ndim1 = size(vacpy,1_ip,kind=4)
       ndim2 = size(vacpy,2_ip,kind=4)
       ndim3 = size(vacpy,3_ip,kind=4)
    end if

    if( ndim1 > 0 .and. ndim2 > 0 .and. ndim3 > 0 ) then

       if( .not. associated(varia) ) then
          call memory_alloca(memor,vanam,vacal,varia,ndim1,ndim2,ndim3,'DO_NOT_INITIALIZE')
       else
          ndim1 = min(ndim1,size(varia,1,kind=4))
          ndim2 = min(ndim2,size(varia,2,kind=4))
          ndim3 = min(ndim3,size(varia,3,kind=4))
       end if

       do idim3 = 1,ndim3
          do idim2 = 1,ndim2
             do idim1 = 1,ndim1
                varia(idim1,idim2,idim3) = vacpy(idim1,idim2,idim3)
             end do
          end do
       end do

    else

       nullify(varia)

    end if

    lzero = .true.
    if(present(wzero)) then
       if( wzero == 'DO_NOT_DEALLOCATE') lzero = .false.
    end if
    if( lzero ) call memory_deallo(memor,vanam,vacal,vacpy)

  end subroutine memory_cpyi43

  subroutine memory_cpyi82(memor,vanam,vacal,vacpy,varia,wzero)
    !
    ! Integer(8)(:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(8),   intent(inout), pointer  :: varia(:,:)
    integer(8),   intent(inout), pointer  :: vacpy(:,:)
    character(*), intent(in),  optional :: wzero
    integer(8)                          :: idim1,idim2
    logical(lg)                         :: lzero
    integer(8)                          :: ndim1,ndim2

    if( .not. associated(vacpy) ) then
       ndim1 = 0
       ndim2 = 0
    else
       ndim1 = size(vacpy,1_ip,kind=8)
       ndim2 = size(vacpy,2_ip,kind=8)
    end if

    if( ndim1 > 0 .and. ndim2 > 0 ) then

       if( .not. associated(varia) ) then
          call memory_alloca(memor,vanam,vacal,varia,ndim1,ndim2,'DO_NOT_INITIALIZE')
       else
          ndim1 = min(ndim1,size(varia,1,kind=8))
          ndim2 = min(ndim2,size(varia,2,kind=8))
       end if

       do idim2 = 1,ndim2
          do idim1 = 1,ndim1
             varia(idim1,idim2) = vacpy(idim1,idim2)
          end do
       end do

    else

       nullify(varia)

    end if

    lzero = .true.
    if(present(wzero)) then
       if( wzero == 'DO_NOT_DEALLOCATE') lzero = .false.
    end if
    if( lzero ) call memory_deallo(memor,vanam,vacal,vacpy)

  end subroutine memory_cpyi82

  subroutine memory_cpyi83(memor,vanam,vacal,vacpy,varia,wzero)
    !
    ! Integer(8)(:,:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(8),   intent(inout), pointer  :: varia(:,:,:)
    integer(8),   intent(inout), pointer  :: vacpy(:,:,:)
    character(*), intent(in),  optional :: wzero
    integer(8)                          :: idim1,idim2,idim3
    logical(lg)                         :: lzero
    integer(8)                          :: ndim1,ndim2,ndim3

    if( .not. associated(vacpy) ) then
       ndim1 = 0
       ndim2 = 0
       ndim3 = 0
    else
       ndim1 = size(vacpy,1_ip,kind=8)
       ndim2 = size(vacpy,2_ip,kind=8)
       ndim3 = size(vacpy,3_ip,kind=8)
    end if

    if( ndim1 > 0 .and. ndim2 > 0 .and. ndim3 > 0 ) then

       if( .not. associated(varia) ) then
          call memory_alloca(memor,vanam,vacal,varia,ndim1,ndim2,ndim3,'DO_NOT_INITIALIZE')
       else
          ndim1 = min(ndim1,size(varia,1,kind=8))
          ndim2 = min(ndim2,size(varia,2,kind=8))
          ndim3 = min(ndim3,size(varia,3,kind=8))
       end if

       do idim3 = 1,ndim3
          do idim2 = 1,ndim2
             do idim1 = 1,ndim1
                varia(idim1,idim2,idim3) = vacpy(idim1,idim2,idim3)
             end do
          end do
       end do

    else

       nullify(varia)

    end if

    lzero = .true.
    if(present(wzero)) then
       if( wzero == 'DO_NOT_DEALLOCATE') lzero = .false.
    end if
    if( lzero ) call memory_deallo(memor,vanam,vacal,vacpy)

  end subroutine memory_cpyi83

  subroutine memory_cpyi1p_1(memor,vanam,vacal,vacpy,varia,wzero)
    !
    ! Type(i1p)(4)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(i1p),    intent(inout), pointer  :: varia(:)
    type(i1p),    intent(inout), pointer  :: vacpy(:)
    character(*), intent(in),  optional :: wzero
    integer(ip)                         :: idim1,idim2
    logical(lg)                         :: lzero
    integer(ip)                         :: ndim1
    integer(ip)                         :: ndim2

    if( .not. associated(vacpy) ) then
       ndim1 = 0
    else
       ndim1 = size(vacpy,1_ip,kind=ip)
    end if

    if( ndim1 > 0 ) then

       if( .not. associated(varia) ) then
          call memory_alloca(memor,vanam,vacal,varia,ndim1,'DO_NOT_INITIALIZE')
          do idim1 = 1,ndim1
             ndim2 = size(vacpy(idim1) % l,1_ip,kind=ip)
             nullify(varia(idim1) % l)
             call memory_alloca(memor,vanam,vacal,varia(idim1) % l,ndim2,'DO_NOT_INITIALIZE')
          end do
       else
          ndim1 = min(ndim1,size(varia,kind=ip))
       end if

       do idim1 = 1,ndim1
          ndim2 = min(size(vacpy(idim1) % l,1_ip,kind=ip),size(varia(idim1) % l,1_ip,kind=ip))
          do idim2 = 1,ndim2
             varia(idim1) % l(idim2) = vacpy(idim1) % l(idim2)
          end do
       end do

    else

       nullify(varia)

    end if

    lzero = .true.
    if(present(wzero)) then
       if( wzero == 'DO_NOT_DEALLOCATE') lzero = .false.
    end if
    if( lzero ) call memory_deallo(memor,vanam,vacal,vacpy)

  end subroutine memory_cpyi1p_1

  subroutine memory_cpyrp1(memor,vanam,vacal,vacpy,varia,wzero)
    !
    ! Real(rp)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    real(rp),     intent(inout), pointer  :: varia(:)
    real(rp),     intent(inout), pointer  :: vacpy(:)
    character(*), intent(in),  optional :: wzero
    integer(ip)                         :: idim1
    logical(lg)                         :: lzero
    integer(ip)                         :: ndim1

    ndim1 = memory_size(vacpy)

    if( ndim1 > 0 ) then

       if( .not. associated(varia) ) then
          call memory_alloca(memor,vanam,vacal,varia,ndim1,'DO_NOT_INITIALIZE')
       else
          ndim1 = min(ndim1,memory_size(varia))
       end if

       do idim1 = 1,ndim1
          varia(idim1) = vacpy(idim1)
       end do

    else

       nullify(varia)

    end if

    lzero = .true.
    if(present(wzero)) then
       if( wzero == 'DO_NOT_DEALLOCATE') lzero = .false.
    end if
    if( lzero ) call memory_deallo(memor,vanam,vacal,vacpy)

  end subroutine memory_cpyrp1

  subroutine memory_cpyrp2(memor,vanam,vacal,vacpy,varia,wzero)
    !
    ! Real(rp)(:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    real(rp),     intent(inout), pointer  :: varia(:,:)
    real(rp),     intent(inout), pointer  :: vacpy(:,:)
    character(*), intent(in),  optional :: wzero
    integer(ip)                         :: idim1,idim2
    logical(lg)                         :: lzero
    integer(ip)                         :: ndim1,ndim2

    if( .not. associated(vacpy) ) then
       ndim1 = 0
       ndim2 = 0
    else
       ndim1 = size(vacpy,1_ip,kind=ip)
       ndim2 = size(vacpy,2_ip,kind=ip)
    end if

    if( ndim1 > 0 .and. ndim2 > 0 ) then

       if( .not. associated(varia) ) then
          call memory_alloca(memor,vanam,vacal,varia,ndim1,ndim2,'DO_NOT_INITIALIZE')
       else
          ndim1 = min(ndim1,size(varia,1,kind=ip))
          ndim2 = min(ndim2,size(varia,2,kind=ip))
       end if

       do idim2 = 1,ndim2
          do idim1 = 1,ndim1
             varia(idim1,idim2) = vacpy(idim1,idim2)
          end do
       end do

    else

       nullify(varia)

    end if

    lzero = .true.
    if(present(wzero)) then
       if( wzero == 'DO_NOT_DEALLOCATE') lzero = .false.
    end if
    if( lzero ) call memory_deallo(memor,vanam,vacal,vacpy)

  end subroutine memory_cpyrp2

 subroutine memory_cpyrp3(memor,vanam,vacal,vacpy,varia,wzero)
    !
    ! Real(rp)(:,:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    real(rp),     intent(inout), pointer  :: varia(:,:,:)
    real(rp),     intent(inout), pointer  :: vacpy(:,:,:)
    character(*), intent(in),  optional :: wzero
    integer(ip)                         :: idim1,idim2,idim3
    logical(lg)                         :: lzero
    integer(ip)                         :: ndim1,ndim2,ndim3

    if( .not. associated(vacpy) ) then
       ndim1 = 0
       ndim2 = 0
       ndim3 = 0
    else
       ndim1 = size(vacpy,1_ip,kind=ip)
       ndim2 = size(vacpy,2_ip,kind=ip)
       ndim3 = size(vacpy,3_ip,kind=ip)
    end if

    if( ndim1 > 0 .and. ndim2 > 0 .and. ndim3 > 0 ) then

       if( .not. associated(varia) ) then
          call memory_alloca(memor,vanam,vacal,varia,ndim1,ndim2,ndim3,'DO_NOT_INITIALIZE')
       else
          ndim1 = min(ndim1,size(varia,1,kind=ip))
          ndim2 = min(ndim2,size(varia,2,kind=ip))
          ndim3 = min(ndim3,size(varia,3,kind=ip))
       end if

       do idim3 = 1,ndim3
          do idim2 = 1,ndim2
             do idim1 = 1,ndim1
                varia(idim1,idim2,idim3) = vacpy(idim1,idim2,idim3)
             end do
          end do
       end do

    else

       nullify(varia)

    end if

    lzero = .true.
    if(present(wzero)) then
       if( wzero == 'DO_NOT_DEALLOCATE') lzero = .false.
    end if
    if( lzero ) call memory_deallo(memor,vanam,vacal,vacpy)

  end subroutine memory_cpyrp3

  !----------------------------------------------------------------------
  !
  ! Resize arrays
  !
  ! 1. Nullify(vacpy)
  ! 2. vacpy <= varia, deallocate varia
  ! 3. Allocate varia with new dimension
  ! 4. varia <= vacpy, deallocate vacpy
  !
  !----------------------------------------------------------------------

  subroutine memory_newi41(memor,vanam,vacal,varia,ndim1)
    !
    ! Integer(4)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(4),   intent(inout), pointer  :: varia(:)
    integer(4),   intent(in)            :: ndim1
    integer(4),                pointer  :: vacpy(:)
    logical(lg)                         :: varia_associated

    if( ndim1 == 0_4 ) then
       call memory_deallo(memor,vanam,vacal,varia)
    else
       nullify(vacpy)

       if( associated(varia) ) then
          varia_associated = .true.
       else
          varia_associated = .false.
       end if

       if( varia_associated ) then
          call memory_copy(  memor,vanam,vacal,varia,vacpy)
       else
          call memory_alloca(memor,vanam,vacal,vacpy,ndim1)
       end if

       call memory_alloca(memor,vanam,vacal,varia,ndim1)
       if( varia_associated ) call memory_copy(  memor,vanam,vacal,vacpy,varia)
    end if

  end subroutine memory_newi41

  subroutine memory_newi42(memor,vanam,vacal,varia,ndim1,ndim2)
    !
    ! Integer(4)(:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(4),   intent(inout), pointer  :: varia(:,:)
    integer(4),   intent(in)            :: ndim1
    integer(4),   intent(in)            :: ndim2
    integer(4),                pointer  :: vacpy(:,:)
    logical(lg)                         :: varia_associated

    if( ndim1 * ndim2 == 0_4 ) then
       call memory_deallo(memor,vanam,vacal,varia)
    else
       nullify(vacpy)

       if( associated(varia) ) then
          varia_associated = .true.
       else
          varia_associated = .false.
       end if

       if( varia_associated ) then
          call memory_copy(  memor,vanam,vacal,varia,vacpy)
       else
          call memory_alloca(memor,vanam,vacal,vacpy,ndim1,ndim2)
       end if

       call memory_alloca(memor,vanam,vacal,varia,ndim1,ndim2)
       if( varia_associated ) call memory_copy(  memor,vanam,vacal,vacpy,varia)
    end if

  end subroutine memory_newi42

  subroutine memory_newi43(memor,vanam,vacal,varia,ndim1,ndim2,ndim3)
    !
    ! Integer(4)(:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(4),   intent(inout), pointer  :: varia(:,:,:)
    integer(4),   intent(in)            :: ndim1
    integer(4),   intent(in)            :: ndim2
    integer(4),   intent(in)            :: ndim3
    integer(4),                pointer  :: vacpy(:,:,:)
    logical(lg)                         :: varia_associated

    if( ndim1 * ndim2 * ndim3 == 0_4 ) then
       call memory_deallo(memor,vanam,vacal,varia)
    else
       nullify(vacpy)

       if( associated(varia) ) then
          varia_associated = .true.
       else
          varia_associated = .false.
       end if

       if( varia_associated ) then
          call memory_copy(  memor,vanam,vacal,varia,vacpy)
       else
          call memory_alloca(memor,vanam,vacal,vacpy,ndim1,ndim2,ndim3)
       end if

       call memory_alloca(memor,vanam,vacal,varia,ndim1,ndim2,ndim3)
       if( varia_associated ) call memory_copy(  memor,vanam,vacal,vacpy,varia)
    end if

  end subroutine memory_newi43

  subroutine memory_newi81(memor,vanam,vacal,varia,ndim1)
    !
    ! Integer(8)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(8),   intent(inout), pointer  :: varia(:)
    integer(8),   intent(in)            :: ndim1
    integer(8),                pointer  :: vacpy(:)
    logical(lg)                         :: varia_associated

    if( ndim1 == 0_8 ) then
       call memory_deallo(memor,vanam,vacal,varia)
    else
       nullify(vacpy)

       if( associated(varia) ) then
          varia_associated = .true.
       else
          varia_associated = .false.
       end if

       if( varia_associated ) then
          call memory_copy(  memor,vanam,vacal,varia,vacpy)
       else
          call memory_alloca(memor,vanam,vacal,vacpy,ndim1)
       end if

       call memory_alloca(memor,vanam,vacal,varia,ndim1)
       if( varia_associated ) call memory_copy(  memor,vanam,vacal,vacpy,varia)
    end if

  end subroutine memory_newi81

  subroutine memory_newi82(memor,vanam,vacal,varia,ndim1,ndim2)
    !
    ! Integer(8)(:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(8),   intent(inout), pointer  :: varia(:,:)
    integer(8),   intent(in)            :: ndim1
    integer(8),   intent(in)            :: ndim2
    integer(8),                pointer  :: vacpy(:,:)
    logical(lg)                         :: varia_associated

    if( ndim1 * ndim2 == 0_8 ) then
       call memory_deallo(memor,vanam,vacal,varia)
    else
       nullify(vacpy)

       if( associated(varia) ) then
          varia_associated = .true.
       else
          varia_associated = .false.
       end if

       if( varia_associated ) then
          call memory_copy(  memor,vanam,vacal,varia,vacpy)
       else
          call memory_alloca(memor,vanam,vacal,vacpy,ndim1,ndim2)
       end if

       !call memory_deallo(memor,vanam,vacal,varia)
       !varia => vacpy

       call memory_alloca(memor,vanam,vacal,varia,ndim1,ndim2)
       if( varia_associated ) call memory_copy(  memor,vanam,vacal,vacpy,varia)
    end if

  end subroutine memory_newi82

  subroutine memory_newi83(memor,vanam,vacal,varia,ndim1,ndim2,ndim3)
    !
    ! Integer(8)(:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(8),   intent(inout), pointer  :: varia(:,:,:)
    integer(8),   intent(in)            :: ndim1
    integer(8),   intent(in)            :: ndim2
    integer(8),   intent(in)            :: ndim3
    integer(8),                pointer  :: vacpy(:,:,:)
    logical(lg)                         :: varia_associated

    if( ndim1 * ndim2 * ndim3 == 0_8 ) then
       call memory_deallo(memor,vanam,vacal,varia)
    else
       nullify(vacpy)

       if( associated(varia) ) then
          varia_associated = .true.
       else
          varia_associated = .false.
       end if

       if( varia_associated ) then
          call memory_copy(  memor,vanam,vacal,varia,vacpy)
       else
          call memory_alloca(memor,vanam,vacal,vacpy,ndim1,ndim2,ndim3)
       end if

       !call memory_deallo(memor,vanam,vacal,varia)
       !varia => vacpy

       call memory_alloca(memor,vanam,vacal,varia,ndim1,ndim2,ndim3)
       if( varia_associated ) call memory_copy(  memor,vanam,vacal,vacpy,varia)
    end if

  end subroutine memory_newi83

  subroutine memory_newr83(memor,vanam,vacal,varia,ndim1,ndim2,ndim3)
    !
    ! Real(rp)(:,:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    real(rp),     intent(inout), pointer  :: varia(:,:,:)
    integer(ip),  intent(in)            :: ndim1
    integer(ip),  intent(in)            :: ndim2
    integer(ip),  intent(in)            :: ndim3
    real(rp),                  pointer  :: vacpy(:,:,:)
    logical(lg)                         :: varia_associated

    if( ndim1 *ndim2 * ndim3 == 0_ip ) then
       call memory_deallo(memor,vanam,vacal,varia)
    else
       nullify(vacpy)

       if( associated(varia) ) then
          varia_associated = .true.
       else
          varia_associated = .false.
       end if

       if( varia_associated ) then
          call memory_copy(  memor,vanam,vacal,varia,vacpy)
       else
          call memory_alloca(memor,vanam,vacal,vacpy,ndim1,ndim2,ndim3)
       end if

       call memory_alloca(memor,vanam,vacal,varia,ndim1,ndim2,ndim3)
       if( varia_associated ) call memory_copy(  memor,vanam,vacal,vacpy,varia)
    end if

  end subroutine memory_newr83

  subroutine memory_newi1p_1(memor,vanam,vacal,varia,ndim1)
    !
    ! Type(I1P)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    type(i1p),    intent(inout), pointer  :: varia(:)
    integer(ip),  intent(in)            :: ndim1
    type(i1p),                 pointer  :: vacpy(:)
    logical(lg)                         :: varia_associated

    if( ndim1 == 0_ip ) then
       call memory_deallo(memor,vanam,vacal,varia)
    else
       nullify(vacpy)

       if( associated(varia) ) then
          varia_associated = .true.
       else
          varia_associated = .false.
       end if

       if( varia_associated ) then
          call memory_copy(  memor,vanam,vacal,varia,vacpy)
       else
          call memory_alloca(memor,vanam,vacal,vacpy,ndim1)
       end if

       !call memory_deallo(memor,vanam,vacal,varia)
       !varia => vacpy
       call memory_alloca(memor,vanam,vacal,varia,ndim1)
       if( varia_associated ) call memory_copy(  memor,vanam,vacal,vacpy,varia)
    end if

  end subroutine memory_newi1p_1

  subroutine memory_allxp1(memor,vanam,vacal,varia,ndim1,wzero)
    !
    ! Complex(rp)(:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    complex(rp),  intent(inout), pointer  :: varia(:)
    character(*), intent(in),  optional :: wzero
    integer(4)                          :: istat
    integer(ip)                         :: idim1
    logical(lg)                         :: lzero

    if( ndim1 > 0 ) then

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       allocate( varia(ndim1) , stat = istat )
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = -size(varia,kind=8)*int(kind(varia),8)*2
          !lbytm = -int(sizeof(varia),8)
          if( lzero ) then
             do idim1 = 1,ndim1
                varia(idim1) = CMPLX(0.0_rp,0.0_rp,kind=rp)
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if

       call memory_info(memor,vanam,vacal,'real')

    else

       nullify(varia)

    end if

  end subroutine memory_allxp1

  subroutine memory_allxp2(memor,vanam,vacal,varia,ndim1,ndim2,wzero)
    !
    ! Complex(rp)(:,:)
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    integer(ip),  intent(in)            :: ndim1
    integer(ip),  intent(in)            :: ndim2
    complex(rp),  intent(inout), pointer  :: varia(:,:)
    character(*), intent(in),  optional :: wzero
    integer(4)                          :: istat
    integer(ip)                         :: idim1,idim2
    logical(lg)                         :: lzero

    if( ndim1*ndim2 > 0 ) then

       if( associated(varia) ) then
          write(*,*) 'POINTER ALREADY ASSOCIATED: ',kfl_paral,trim(vanam),trim(vacal)
          call runend('MOD_MEMORY: POINTER ALREADY ASSOCIATED')
       end if

       allocate( varia(ndim1,ndim2) , stat = istat )
       lzero = .true.
       if(present(wzero)) then
          if( wzero == 'DO_NOT_INITIALIZE') then
             lzero = .false.
          end if
       end if

       if( istat == 0 ) then
          lbytm = -size(varia,kind=8)*int(kind(varia),8)*2
          !lbytm = -int(sizeof(varia),8)
          if( lzero ) then
             do idim2 = 1,ndim2
                do idim1 = 1,ndim1
                   varia(idim1,idim2) = CMPLX(0.0_rp,0.0_rp,kind=rp)
                end do
             end do
          end if
       else
          call memory_error(0_ip,vanam,vacal,istat)
       end if
       call memory_info(memor,vanam,vacal,'real')

    else

       nullify(varia)

    end if

  end subroutine memory_allxp2

  !----------------------------------------------------------------------
  !
  ! Renumber arrays
  !
  ! 1. Nullify vacpy
  ! 2. Copy vacpy <= varia, deallocate varia
  ! 3. Reallocate varia
  ! 4. Renumber varia using vacpy
  ! 5. Deallocate vacpy
  !
  !----------------------------------------------------------------------

  subroutine memory_reni41(memor,vanam,vacal,varia,lrenu)
    !
    ! Integer(4)(:)
    !
    implicit none
    character(*), intent(in)           :: vanam         !< Variable name
    character(*), intent(in)           :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)        :: memor(2)      !< Memory counter
    integer(4),   intent(inout), pointer :: varia(:)
    integer(ip),  intent(in),  pointer :: lrenu(:)
    integer(4),                pointer :: vacpy(:)
    integer(ip)                        :: idim1_new,idim1_old
    integer(ip)                        :: ndim1_new,ndim1_old

    nullify(vacpy)

    ndim1_old = memory_size(varia)
    call memory_copy(memor,vanam,vacal,varia,vacpy)
    ndim1_new = 0
    do idim1_old = 1,ndim1_old
       if( lrenu(idim1_old) /= 0 ) ndim1_new = ndim1_new + 1
    end do
    call memory_alloca(memor,vanam,vacal,varia,ndim1_new)
    do idim1_old = 1,ndim1_old
       if( lrenu(idim1_old) /= 0 ) then
          idim1_new        = lrenu(idim1_old)
          varia(idim1_new) = vacpy(idim1_old)
       end if
    end do
    call memory_deallo(memor,vanam,vacal,vacpy)

  end subroutine memory_reni41

  subroutine memory_reni81(memor,vanam,vacal,varia,lrenu)
    !
    ! Integer(8)(:)
    !
    implicit none
    character(*), intent(in)           :: vanam         !< Variable name
    character(*), intent(in)           :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)        :: memor(2)      !< Memory counter
    integer(8),   intent(inout), pointer :: varia(:)
    integer(8),   intent(in),  pointer :: lrenu(:)
    integer(8),                pointer :: vacpy(:)
    integer(8)                         :: idim1_new,idim1_old
    integer(8)                         :: ndim1_new,ndim1_old

    nullify(vacpy)

    ndim1_old = int(memory_size(varia),8)
    call memory_copy(memor,vanam,vacal,varia,vacpy)
    ndim1_new = 0
    do idim1_old = 1,ndim1_old
       if( lrenu(idim1_old) /= 0 ) ndim1_new = ndim1_new + 1
    end do
    call memory_alloca(memor,vanam,vacal,varia,ndim1_new)
    do idim1_old = 1,ndim1_old
       if( lrenu(idim1_old) /= 0 ) then
          idim1_new        = lrenu(idim1_old)
          varia(idim1_new) = vacpy(idim1_old)
       end if
    end do
    call memory_deallo(memor,vanam,vacal,vacpy)

  end subroutine memory_reni81

  subroutine memory_renrp1(memor,vanam,vacal,varia,lrenu)
    !
    ! Real(rp)(:)
    !
    implicit none
    character(*), intent(in)           :: vanam         !< Variable name
    character(*), intent(in)           :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)        :: memor(2)      !< Memory counter
    real(rp),     intent(inout), pointer :: varia(:)
    integer(ip),  intent(in),  pointer :: lrenu(:)
    real(rp),                  pointer :: vacpy(:)
    integer(ip)                        :: idim1_new,idim1_old
    integer(ip)                        :: ndim1_new,ndim1_old

    nullify(vacpy)

    ndim1_old = memory_size(varia)
    call memory_copy(memor,vanam,vacal,varia,vacpy)
    ndim1_new = 0
    do idim1_old = 1,ndim1_old
       if( lrenu(idim1_old) /= 0 ) ndim1_new = ndim1_new + 1
    end do
    call memory_alloca(memor,vanam,vacal,varia,ndim1_new)
    do idim1_old = 1,ndim1_old
       if( lrenu(idim1_old) /= 0 ) then
          idim1_new        = lrenu(idim1_old)
          varia(idim1_new) = vacpy(idim1_old)
       end if
    end do
    call memory_deallo(memor,vanam,vacal,vacpy)

  end subroutine memory_renrp1

  subroutine memory_renrp2(memor,vanam,vacal,varia,lrenu,idime)
    !
    ! Real(rp)(:,:)
    ! If renumbering LRENU is present, it assumes by default that
    ! the second dimension is the one that should be renumbered
    !
    implicit none
    character(*), intent(in)            :: vanam         !< Variable name
    character(*), intent(in)            :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)         :: memor(2)      !< Memory counter
    real(rp),     intent(inout), pointer  :: varia(:,:)
    integer(ip),  intent(in),  pointer  :: lrenu(:)
    real(rp),                  pointer  :: vacpy(:,:)
    integer(ip),  intent(in),  optional :: idime
    integer(ip)                         :: idim1_new,idim1_old
    integer(ip)                         :: ndim1_new,ndim1_old
    integer(ip)                         :: ndim2,kdime_1,kdime_2

    nullify(vacpy)
    !
    ! If IDIME is not present: KDIME_1=2, KDIME_2=1
    !
    if( present(idime) ) then
       kdime_1 = idime
       call runend('MEMORY_RENRP2: NOT CODED')
    else
       kdime_1 = 2
    end if
    if( kdime_1 == 1 ) then
       kdime_2 = 2
    else if( kdime_1 == 2 ) then
       kdime_2 = 1
    else
       call runend('MEMORY_RENRP2: WRONG DIMENSION')
    end if
    !
    ! VARIA(KDIME_2,KDIME_1)
    ! NDIM1_OLD <= size(VARIA,2)
    !
    ndim1_old = memory_size(varia,kdime_1)
    ndim2     = memory_size(varia,kdime_2)
    call memory_copy(memor,vanam,vacal,varia,vacpy)
    ndim1_new = 0
    do idim1_old = 1,ndim1_old
       if( lrenu(idim1_old) /= 0 ) ndim1_new = ndim1_new + 1
    end do
    call memory_alloca(memor,vanam,vacal,varia,ndim2,ndim1_new)
    do idim1_old = 1,ndim1_old
       if( lrenu(idim1_old) /= 0 ) then
          idim1_new = lrenu(idim1_old)
          varia(1:ndim2,idim1_new) = vacpy(1:ndim2,idim1_old)
       end if
    end do
    call memory_deallo(memor,vanam,vacal,vacpy)

  end subroutine memory_renrp2

  subroutine memory_renrp3(memor,vanam,vacal,varia,lrenu,idime)
    !
    ! Real(rp)(:,:)
    ! If renumbering LRENU is present, it assumes by default that
    ! the second dimension is the one that should be renumbered
    !
    implicit none
    character(*), intent(in)              :: vanam         !< Variable name
    character(*), intent(in)              :: vacal         !< Calling subroutine name
    integer(8),   intent(inout)           :: memor(2)      !< Memory counter
    real(rp),     intent(inout), pointer  :: varia(:,:,:)
    integer(ip),  intent(in),    pointer  :: lrenu(:)
    real(rp),                    pointer  :: vacpy(:,:,:)
    integer(ip),  intent(in),    optional :: idime
    integer(ip)                           :: idim1_new,idim1_old
    integer(ip)                           :: ndim1_new,ndim1_old
    integer(ip)                           :: ndim2,kdime_1,kdime_2
    integer(ip)                           :: ndim3

    nullify(vacpy)
    !
    ! If IDIME is not present: KDIME_1=2, KDIME_2=1
    !
    if( present(idime) ) then
       kdime_1 = idime
       call runend('RENPR3 NOT CODED')
    else
       kdime_1 = 2
    end if
    if( kdime_1 == 1 ) then
       kdime_2 = 2
    else if( kdime_1 == 2 ) then
       kdime_2 = 1
    else
       call runend('MEMORY_RENRP2: WRONG DIMENSION')
    end if
    !
    ! VARIA(KDIME_2,KDIME_1)
    ! NDIM1_OLD <= size(VARIA,2)
    !
    ndim1_old = size(varia,kdime_1,kind=ip)
    ndim2     = size(varia,kdime_2,kind=ip)
    ndim3     = size(varia,3,kind=ip)
    call memory_copy(memor,vanam,vacal,varia,vacpy)
    ndim1_new = 0
    do idim1_old = 1,ndim1_old
       if( lrenu(idim1_old) /= 0 ) ndim1_new = ndim1_new + 1
    end do
    call memory_alloca(memor,vanam,vacal,varia,ndim2,ndim1_new,ndim3)
    do idim1_old = 1,ndim1_old
       if( lrenu(idim1_old) /= 0 ) then
          idim1_new = lrenu(idim1_old)
          varia(1:ndim2,idim1_new,1:ndim3) = vacpy(1:ndim2,idim1_old,1:ndim3)
       end if
    end do
    call memory_deallo(memor,vanam,vacal,vacpy)

  end subroutine memory_renrp3

  !----------------------------------------------------------------------
  !
  ! Allocate minimum size if pointers are not associated
  !
  !----------------------------------------------------------------------
  !
  ! Real(rp)(:)
  !
  subroutine memory_alloca_min_rp1(varia)
    real(rp), intent(inout), pointer :: varia(:)
    if( .not. associated(varia) ) then
       allocate( varia(1) )
       varia= 0.0_rp
    end if
  end subroutine memory_alloca_min_rp1
  !
  ! Real(rp)(:,:)
  !
  subroutine memory_alloca_min_rp2(varia)
    real(rp), intent(inout), pointer :: varia(:,:)
    if( .not. associated(varia) ) then
       allocate( varia(1,1) )
       varia = 0.0_rp
    end if
  end subroutine memory_alloca_min_rp2
  !
  ! Real(rp)(:,:,:)
  !
  subroutine memory_alloca_min_rp3(varia)
    real(rp), intent(inout), pointer :: varia(:,:,:)
    if( .not. associated(varia) ) then
       allocate( varia(1,1,1) )
       varia = 0.0_rp
    end if
  end subroutine memory_alloca_min_rp3
  !
  ! Real(rp)(:,:,:,:)
  !
  subroutine memory_alloca_min_rp4(varia)
    real(rp), intent(inout), pointer :: varia(:,:,:,:)
    if( .not. associated(varia) ) then
       allocate( varia(1,1,1,1) )
       varia = 0.0_rp
    end if
  end subroutine memory_alloca_min_rp4
  !
  ! integer(ip)(:)
  !
  subroutine memory_alloca_min_ip1(varia)
    integer(ip), intent(inout), pointer :: varia(:)
    if( .not. associated(varia) ) then
       allocate( varia(1) )
       varia = 0_ip
    end if
  end subroutine memory_alloca_min_ip1
  !
  ! integer(ip)(:,:)
  !
  subroutine memory_alloca_min_ip2(varia)
    integer(ip), intent(inout), pointer :: varia(:,:)
    if( .not. associated(varia) ) then
       allocate( varia(1,1) )
       varia = 0_ip
    end if
  end subroutine memory_alloca_min_ip2
  !
  ! integer(ip)(:,:,:)
  !
  subroutine memory_alloca_min_ip3(varia)
    integer(ip), intent(inout), pointer :: varia(:,:,:)
    if( .not. associated(varia) ) then
       allocate( varia(1,1,1) )
       varia = 0_ip
    end if
  end subroutine memory_alloca_min_ip3
  !
  ! logical(lg)(:)
  !
  subroutine memory_alloca_min_lg1(varia)
    logical(lg), intent(inout), pointer :: varia(:)
    if( .not. associated(varia) ) then
       allocate( varia(1) )
       varia = .false.
    end if
  end subroutine memory_alloca_min_lg1
  !
  ! logical(lg)(:,:)
  !
  subroutine memory_alloca_min_lg2(varia)
    logical(lg), intent(inout), pointer :: varia(:,:)
    if( .not. associated(varia) ) then
       allocate( varia(1,1) )
       varia = .false.
    end if
  end subroutine memory_alloca_min_lg2
  !
  ! logical(lg)(:,:,:)
  !
  subroutine memory_alloca_min_lg3(varia)
    logical(lg), intent(inout), pointer :: varia(:,:,:)
    if( .not. associated(varia) ) then
       allocate( varia(1,1,1) )
       varia = .false.
    end if
  end subroutine memory_alloca_min_lg3

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    19/11/2015
  !> @brief   Errro message
  !> @details Write an error message when memory could not be allocated
  !>          or deallocated
  !>
  !-----------------------------------------------------------------------

  subroutine memory_error(itask,vanam,vacal,istat)
    implicit none
    integer(ip),   intent(in) :: itask
    integer(4),    intent(in) :: istat
    integer(ip)               :: ibyte
    real(rp)                  :: rbyte
    character(6)              :: lbyte
    integer(ip)               :: ibyte_tot
    real(rp)                  :: rbyte_tot
    character(6)              :: lbyte_tot
    character*(*), intent(in) :: vanam         !< Variable name
    character*(*), intent(in) :: vacal         !< Variable name
    character(200)            :: wmess
    character(20)             :: wmes2,wmes3,wmes4

    if( itask == 0 ) then
       !
       ! Allocation
       !
       if(      lbytm >= Gbytes_i8 ) then
          rbyte = Gbytes
          lbyte = 'GBYTES'
       else if( lbytm >= Mbytes_i8 ) then
          rbyte = Mbytes
          lbyte = 'MBYTES'
       else if( lbytm >= Kbytes_i8 ) then
          rbyte = Kbytes
          lbyte = 'KBYTES'
       else
          rbyte = 1.0_rp
          lbyte = 'BYTES'
       end if
       ibyte = int(real(lbytm,rp)/rbyte,KIND=ip)
       wmes2 = memory_intost(ibyte)
       wmes3 = memory_intost(int(istat,ip))

       if(      mem_curre >= Gbytes_i8 ) then
          rbyte_tot = Gbytes
          lbyte_tot = 'GBYTES'
       else if( mem_curre >= Mbytes_i8 ) then
          rbyte_tot = Mbytes
          lbyte_tot = 'MBYTES'
       else if( mem_curre >= Kbytes_i8 ) then
          rbyte_tot = Kbytes
          lbyte_tot = 'KBYTES'
       else
          rbyte_tot = 1.0_rp
          lbyte_tot = 'BYTES'
       end if
       ibyte_tot = int(real(mem_curre,rp)/rbyte_tot)
       wmes4     = memory_intost(ibyte_tot)

       wmess = trim(vacal)&
            //': MEMORY FOR '//trim(vanam)&
            //' COULD NOT BE ALLOCATED.'&
            //' RUN TIME ERROR= '//trim(wmes3)&
            //' WHEN ALLOCATING '//trim(wmes2)//' '//trim(lbyte)&
            //' OUT OF '//trim(wmes4)//' '//trim(lbyte_tot)
       call runend(trim(wmess))

    else if( itask == 1 ) then
       !
       ! Reallocation
       !
       call runend(trim(vacal)//': MEMORY FOR '//trim(vanam)//' COULD NOT BE REALLOCATED')

    else
       !
       ! Deallocation
       !
       call runend(trim(vacal)//': MEMORY FOR '//trim(vanam)//' COULD NOT BE DEALLOCATED')

    end if

  end subroutine memory_error

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    03/03/2016
  !> @brief   Memory info
  !> @details Returns memory scaling info
  !>
  !-----------------------------------------------------------------------

  subroutine memory_unit(memor_value,memor_char,memor_factor)
    implicit none
    integer(8),    intent(in)  :: memor_value  !< Memory in bytes
    character(6),  intent(out) :: memor_char   !< Memory unit character
    real(rp),      intent(out) :: memor_factor !< Memory scaling

    if( memor_value >= Gbytes_i8 ) then
       memor_factor = 1.0_rp / Gbytes
       memor_char   = 'Gbytes'
    else if( memor_value >= Mbytes_i8 ) then
       memor_factor = 1.0_rp / Mbytes
       memor_char   = 'Mbytes'
    else if( memor_value >= Kbytes_i8 ) then
       memor_factor = 1.0_rp / Kbytes
       memor_char   = 'kbytes'
    else
       memor_factor = 1.0_rp
       memor_char   = 'bytes '
    end if

  end subroutine memory_unit

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    19/11/2015
  !> @brief   Integer to string
  !> @details Convert an integer(ip) to a string
  !>
  !-----------------------------------------------------------------------

  function memory_intost(integ)

    integer(ip)   :: integ
    integer(4)    :: integ4
    character(20) :: memory_intost
    character(20) :: intaux

    integ4 = int(integ,4)
    write(intaux,*) integ4
    memory_intost = adjustl(intaux)

  end function memory_intost

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    19/11/2015
  !> @brief   Memory info
  !> @details This routine computes some info on memory
  !>
  !-----------------------------------------------------------------------

  subroutine memory_info(memor,vanam,vacal,vatyp)

    character(*), intent(in)    :: vanam         !< Variable name
    character(*), intent(in)    :: vacal         !< Calling subroutine name
    integer(8),   intent(inout) :: memor(2)      !< Memory counter
    character(*), intent(in)    :: vatyp
    !
    ! Update input memory counter
    !
    memor(1) = memor(1)+lbytm
    memor(2) = max(memor(2),memor(1))
    !
    ! Update this module memory counters
    !
    mem_curre = mem_curre+lbytm
    mem_maxim = max(mem_curre,mem_maxim)
    !
    ! Write memory info
    !
    call memory_output_info(memor,vanam,vacal,vatyp)
    !
    ! Memory counter of variables
    !
    if( kfl_varcount == 1 ) call memory_variable_counter(vanam,lbytm)
    
  end subroutine memory_info

  !-----------------------------------------------------------------------
  !>
  !> @author  Guillaume Houzeaux
  !> @date    19/11/2015
  !> @brief   Size of some arrays
  !> @details This routine computes the total size of arrays.
  !>          If the pointer is nullified, it returns 0.
  !>
  !-----------------------------------------------------------------------

  function memory_sizerp1(varia)
    !
    ! Real(rp)(:)
    !
    implicit none
    real(rp),     intent(in),    pointer  :: varia(:)
    integer(ip)                           :: memory_sizerp1

    if( associated(varia) ) then
       memory_sizerp1 = size(varia,kind=ip)
    else
       memory_sizerp1 = 0_ip
    end if

  end function memory_sizerp1

  function memory_sizerp2(varia,ndim1)
    !
    ! Real(rp)(:,:)
    !
    implicit none
    real(rp),     intent(in),    pointer  :: varia(:,:)
    integer(ip),  intent(in),    optional :: ndim1
    integer(ip)                           :: memory_sizerp2

    if( associated(varia) ) then
       if( present(ndim1) ) then
          if( ndim1 >= 1 .and. ndim1 <= 2 ) then
             memory_sizerp2 = size(varia,ndim1,kind=ip)
          else
             call runend('MEMORY_SIZE: WRONG SIZE')
          end if
       else
          memory_sizerp2 = size(varia,kind=ip)
       end if
    else
       memory_sizerp2 = 0_ip
    end if

  end function memory_sizerp2

  function memory_sizerp3(varia,ndim1)
    !
    ! Real(rp)(:,:,:)
    !
    implicit none
    real(rp),     intent(in),    pointer  :: varia(:,:,:)
    integer(ip),  intent(in),    optional :: ndim1
    integer(ip)                           :: memory_sizerp3

    if( associated(varia) ) then
       if( present(ndim1) ) then
          if( ndim1 >= 1 .and. ndim1 <= 3 ) then
             memory_sizerp3 = size(varia,ndim1,kind=ip)
          else
             call runend('MEMORY_SIZE: WRONG SIZE')
          end if
       else
          memory_sizerp3 = size(varia,kind=ip)
       end if
    else
       memory_sizerp3 = 0_ip
    end if

  end function memory_sizerp3

  function memory_sizeip1_4(varia)
    !
    ! Int(ip)(:)
    !
    implicit none
    integer(4),   intent(in), pointer :: varia(:)
    integer(ip)                       :: memory_sizeip1_4

    if( associated(varia) ) then
       memory_sizeip1_4 = size(varia,kind=ip)
    else
       memory_sizeip1_4 = 0_ip
    end if

  end function memory_sizeip1_4

  function memory_sizeip1_8(varia)
    !
    ! Int(ip)(:)
    !
    implicit none
    integer(8),   intent(in), pointer :: varia(:)
    integer(ip)                       :: memory_sizeip1_8

    if( associated(varia) ) then
       memory_sizeip1_8 = size(varia,kind=ip)
    else
       memory_sizeip1_8 = 0_ip
    end if

  end function memory_sizeip1_8

  function memory_sizeip2_4(varia,ndim1)
    !
    ! Int(ip)(:,:)
    !
    implicit none
    integer(4),   intent(in), pointer  :: varia(:,:)
    integer(ip),  intent(in), optional :: ndim1
    integer(ip)                        :: memory_sizeip2_4

    if( associated(varia) ) then
       if( present(ndim1) ) then
          if( ndim1 >= 1 .and. ndim1 <= 2 ) then
             memory_sizeip2_4 = size(varia,ndim1,kind=ip)
          else
             call runend('MEMORY_SIZE: WRONG SIZE')
          end if
       else
          memory_sizeip2_4 = size(varia,kind=ip)
       end if
    else
       memory_sizeip2_4 = 0_ip
    end if

  end function memory_sizeip2_4

  function memory_sizeip2_8(varia,ndim1)
    !
    ! Int(ip)(:,:)
    !
    implicit none
    integer(8),   intent(in), pointer  :: varia(:,:)
    integer(ip),  intent(in), optional :: ndim1
    integer(ip)                        :: memory_sizeip2_8

    if( associated(varia) ) then
       if( present(ndim1) ) then
          if( ndim1 >= 1 .and. ndim1 <= 2 ) then
             memory_sizeip2_8 = size(varia,ndim1,kind=ip)
          else
             call runend('MEMORY_SIZE: WRONG SIZE')
          end if
       else
          memory_sizeip2_8 = size(varia,kind=ip)
       end if
    else
       memory_sizeip2_8 = 0_ip
    end if

  end function memory_sizeip2_8

  function memory_sizeip3_4(varia,ndim1)
    !
    ! Int(ip)(:,:,:)
    !
    implicit none
    integer(4),   intent(in), pointer  :: varia(:,:,:)
    integer(ip),  intent(in), optional :: ndim1
    integer(ip)                        :: memory_sizeip3_4

    if( associated(varia) ) then
       if( present(ndim1) ) then
          if( ndim1 >= 1 .and. ndim1 <= 3 ) then
             memory_sizeip3_4 = size(varia,ndim1,kind=ip)
          else
             call runend('MEMORY_SIZE: WRONG SIZE')
          end if
       else
          memory_sizeip3_4 = size(varia,kind=ip)
       end if
    else
       memory_sizeip3_4 = 0_ip
    end if

  end function memory_sizeip3_4

  function memory_sizeip3_8(varia,ndim1)
    !
    ! Int(ip)(:,:,:)
    !
    implicit none
    integer(8),   intent(in), pointer  :: varia(:,:,:)
    integer(ip),  intent(in), optional :: ndim1
    integer(ip)                        :: memory_sizeip3_8

    if( associated(varia) ) then
       if( present(ndim1) ) then
          if( ndim1 >= 1 .and. ndim1 <= 3 ) then
             memory_sizeip3_8 = size(varia,ndim1,kind=ip)
          else
             call runend('MEMORY_SIZE: WRONG SIZE')
          end if
       else
          memory_sizeip3_8 = size(varia,kind=ip)
       end if
    else
       memory_sizeip3_8 = 0_ip
    end if

  end function memory_sizeip3_8

  function memory_sizelg(varia)
    !
    ! Int(ip)(:)
    !
    implicit none
    logical(lg),  intent(in), pointer :: varia(:)
    integer(ip)                       :: memory_sizelg

    if( associated(varia) ) then
       memory_sizelg = int(size(varia,1,ip),ip)
    else
       memory_sizelg = 0_ip
    end if

  end function memory_sizelg

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-09-19
  !> @brief   Add memory to memory counter
  !> @details Add memory to the memory counter of this module
  !>
  !-----------------------------------------------------------------------

  subroutine memory_add_to_memory_counter(lbytm_in,MEMORY_COUNTER,VARIABLE_NAME,CALLING_SUBROUTINE)

    integer(8),   intent(in)              :: lbytm_in           !< Memory in bytes
    integer(8),   intent(inout), optional :: MEMORY_COUNTER(2)  !< Memory counter
    character(*), intent(in),    optional :: VARIABLE_NAME      !< Name of variable
    character(*), intent(in),    optional :: CALLING_SUBROUTINE !< The calling subroutine
    integer(8)                            :: memor_loc(2)
    character(200)                        :: vanam_loc
    character(200)                        :: vacal_loc

    lbytm     = lbytm_in
    memor_loc = 0_8
    vanam_loc = 'UNKNOWN'
    vacal_loc = 'unknown'
    if( present(MEMORY_COUNTER)     ) memor_loc = MEMORY_COUNTER
    if( present(VARIABLE_NAME)      ) vanam_loc = trim(VARIABLE_NAME)
    if( present(CALLING_SUBROUTINE) ) vacal_loc = trim(CALLING_SUBROUTINE)

    call memory_info(memor_loc,trim(vanam_loc),trim(vacal_loc),'unknown')

  end subroutine memory_add_to_memory_counter

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2018-09-19
  !> @brief   Add memory to memory counter
  !> @details Add memory to the memory counter of this module
  !>
  !-----------------------------------------------------------------------

  subroutine memory_remove_from_memory_counter(lbytm_in,MEMORY_COUNTER,VARIABLE_NAME,CALLING_SUBROUTINE)

    integer(8),   intent(in)              :: lbytm_in           !< Memory in bytes
    integer(8),   intent(inout), optional :: MEMORY_COUNTER(2)  !< Memory counter
    character(*), intent(in),    optional :: VARIABLE_NAME      !< Name of variable
    character(*), intent(in),    optional :: CALLING_SUBROUTINE !< The calling subroutine
    integer(8)                            :: memor_loc(2)
    character(200)                        :: vanam_loc
    character(200)                        :: vacal_loc

    lbytm     = -lbytm_in
    memor_loc = 0_8
    vanam_loc = 'UNKNOWN'
    vacal_loc = 'unknown'
    if( present(MEMORY_COUNTER)     ) memor_loc = MEMORY_COUNTER
    if( present(VARIABLE_NAME)      ) vanam_loc = trim(VARIABLE_NAME)
    if( present(CALLING_SUBROUTINE) ) vacal_loc = trim(CALLING_SUBROUTINE)

    call memory_info(memor_loc,trim(vanam_loc),trim(vacal_loc),'unknown')

  end subroutine memory_remove_from_memory_counter

end module mod_memory
!> @}

