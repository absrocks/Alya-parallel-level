!-----------------------------------------------------------------------
!> @addtogroup Kernel
!> @{
!> @file    mod_std.f90
!> @author  houzeaux
!> @date    2019-03-06
!> @brief   Standard correction
!> @details Some compilers do not respect fortran2008 standard
!>          This subroutine tries to solver this...
!>
!-----------------------------------------------------------------------

module mod_std

  use def_kintyp, only : ip,rp,lg
  implicit none
  integer,     parameter :: qp = 8 ! Quadric<< precision for PGI
  real(rp),    pointer   :: std_real_1(:)
  integer(ip), pointer   :: std_int_1(:)

  type pgi_dummy1
     logical(lg) :: dummy
  end type pgi_dummy1

  type pgi_dummy2
     logical(lg) :: dummy
  end type pgi_dummy2

  type pgi_dummy3
     logical(lg) :: dummy
  end type pgi_dummy3

  type pgi_dummy4
     logical(lg) :: dummy
  end type pgi_dummy4

  type pgi_dummy5
     logical(lg) :: dummy
  end type pgi_dummy5

  private

#ifdef __PGI
  interface count
     module procedure alya_count_1_4, &
          &           alya_count_1_8
  end interface count

  interface lbound
     module procedure alya_lbound_41, alya_lbound_42, alya_lbound_43,  &
          &           alya_lbound_81, alya_lbound_82, alya_lbound_83,  &
          &           alya_lbound_r41,alya_lbound_r42,alya_lbound_r43, &
          &           alya_lbound_r81,alya_lbound_r82,alya_lbound_r83, &
          &           alya_lbound_lg1
  end interface lbound
  interface ubound
     module procedure alya_ubound_41, alya_ubound_42, alya_ubound_43,  &
          &           alya_ubound_81, alya_ubound_82, alya_ubound_83,  &
          &           alya_ubound_r41,alya_ubound_r42,alya_ubound_r43, &
          &           alya_ubound_r81,alya_ubound_r82,alya_ubound_r83, &
          &           alya_ubound_lg1
  end interface ubound
  public :: count
  public :: qp
  public :: lbound
  public :: ubound
  public :: pgi_dummy1, pgi_dummy2, pgi_dummy3, pgi_dummy4, pgi_dummy5
  public :: std_real_1
  public :: std_int_1
#endif

  public :: std_initialization
  
contains

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2019-03-06
  !> @brief   Initialization
  !> @details Initialization, mainly for PGI
  !>
  !-----------------------------------------------------------------------

  subroutine std_initialization()

    nullify(std_real_1)
    
  end subroutine std_initialization

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2019-03-06
  !> @brief   Intrinsic count for PGI
  !> @details For PGI, count does not include the KIND argument
  !>
  !-----------------------------------------------------------------------

  integer(ip) function alya_count_1_4(MASK,DIM,KIND)

    logical(lg), intent(in)           :: MASK(:)
    integer(4),  intent(in), optional :: DIM
    integer(4),  intent(in)           :: KIND
    integer(ip)                       :: lboun1,uboun1
    integer(ip)                       :: ii

    lboun1 = lbound(MASK,1)
    uboun1 = ubound(MASK,1)

    alya_count_1_4 = 0_ip
    do ii = lboun1,uboun1
       if( MASK(ii) ) alya_count_1_4 = alya_count_1_4 + 1_ip
    end do
    if( int(KIND,ip) /= ip ) call runend('ALYA_COUNT: WRONG TYPE')

  end function alya_count_1_4

  integer(ip) function alya_count_1_8(MASK,DIM,KIND)

    logical(lg), intent(in)           :: MASK(:)
    integer(8),  intent(in), optional :: DIM
    integer(8),  intent(in)           :: KIND
    integer(ip)                       :: lboun1,uboun1
    integer(ip)                       :: ii

    lboun1 = lbound(MASK,1)
    uboun1 = ubound(MASK,1)

    alya_count_1_8 = 0_ip
    do ii = lboun1,uboun1
       if( MASK(ii) ) alya_count_1_8 = alya_count_1_8 + 1_ip
    end do
    if( int(KIND,ip) /= ip ) call runend('ALYA_COUNT: WRONG TYPE')

  end function alya_count_1_8

  !-----------------------------------------------------------------------
  !>
  !> @author  houzeaux
  !> @date    2019-03-06
  !> @brief   Intrinsic count for PGI
  !> @details For PGI, LBOUND does not include the KIND argument when
  !>          the array is a pointer
  !>
  !-----------------------------------------------------------------------

  integer(ip) function alya_lbound_41(xx,DIM,KIND) result(my_lbound)
    integer(4),  pointer, intent(in) :: xx(:)
    integer(ip),          intent(in) :: DIM
    integer(ip),          intent(in) :: KIND
    if( associated(xx) ) then
       my_lbound = lbound(xx,DIM=DIM)
    else
       my_lbound = 0_ip
    end if
  end function alya_lbound_41

  integer(ip) function alya_lbound_42(xx,DIM,KIND) result(my_lbound)
    integer(4),  pointer, intent(in) :: xx(:,:)
    integer(ip),          intent(in) :: DIM
    integer(ip),          intent(in) :: KIND
    if( associated(xx) ) then
       my_lbound = lbound(xx,DIM=DIM)
    else
       my_lbound = 0_ip
    end if
  end function alya_lbound_42

  integer(ip) function alya_lbound_43(xx,DIM,KIND) result(my_lbound)
    integer(4),  pointer, intent(in) :: xx(:,:,:)
    integer(ip),          intent(in) :: DIM
    integer(ip),          intent(in) :: KIND
    if( associated(xx) ) then
       my_lbound = lbound(xx,DIM=DIM)
    else
       my_lbound = 0_ip
    end if
  end function alya_lbound_43

  integer(ip) function alya_lbound_81(xx,DIM,KIND) result(my_lbound)
    integer(8),  pointer, intent(in) :: xx(:)
    integer(ip),          intent(in) :: DIM
    integer(ip),          intent(in) :: KIND
    if( associated(xx) ) then
       my_lbound = lbound(xx,DIM=DIM)
    else
       my_lbound = 0_ip
    end if
  end function alya_lbound_81

  integer(ip) function alya_lbound_82(xx,DIM,KIND) result(my_lbound)
    integer(8),  pointer, intent(in) :: xx(:,:)
    integer(ip),          intent(in) :: DIM
    integer(ip),          intent(in) :: KIND
    if( associated(xx) ) then
       my_lbound = lbound(xx,DIM=DIM)
    else
       my_lbound = 0_ip
    end if
  end function alya_lbound_82

  integer(ip) function alya_lbound_83(xx,DIM,KIND) result(my_lbound)
    integer(8),  pointer, intent(in) :: xx(:,:,:)
    integer(ip),          intent(in) :: DIM
    integer(ip),          intent(in) :: KIND
    if( associated(xx) ) then
       my_lbound = lbound(xx,DIM=DIM)
    else
       my_lbound = 0_ip
    end if
  end function alya_lbound_83

  integer(ip) function alya_lbound_r41(xx,DIM,KIND) result(my_lbound)
    real(4),    pointer, intent(in) :: xx(:)
    integer(ip),          intent(in) :: DIM
    integer(ip),          intent(in) :: KIND
    if( associated(xx) ) then
       my_lbound = lbound(xx,DIM=DIM)
    else
       my_lbound = 0_ip
    end if
  end function alya_lbound_r41

  integer(ip) function alya_lbound_r42(xx,DIM,KIND) result(my_lbound)
    real(4),    pointer, intent(in) :: xx(:,:)
    integer(ip),          intent(in) :: DIM
    integer(ip),          intent(in) :: KIND
    if( associated(xx) ) then
       my_lbound = lbound(xx,DIM=DIM)
    else
       my_lbound = 0_ip
    end if
  end function alya_lbound_r42

  integer(ip) function alya_lbound_r43(xx,DIM,KIND) result(my_lbound)
    real(4),    pointer, intent(in) :: xx(:,:,:)
    integer(ip),          intent(in) :: DIM
    integer(ip),          intent(in) :: KIND
    if( associated(xx) ) then
       my_lbound = lbound(xx,DIM=DIM)
    else
       my_lbound = 0_ip
    end if
  end function alya_lbound_r43

  integer(ip) function alya_lbound_r81(xx,DIM,KIND) result(my_lbound)
    real(8),    pointer, intent(in) :: xx(:)
    integer(ip),          intent(in) :: DIM
    integer(ip),          intent(in) :: KIND
    if( associated(xx) ) then
       my_lbound = lbound(xx,DIM=DIM)
    else
       my_lbound = 0_ip
    end if
  end function alya_lbound_r81

  integer(ip) function alya_lbound_r82(xx,DIM,KIND) result(my_lbound)
    real(8),    pointer, intent(in) :: xx(:,:)
    integer(ip),          intent(in) :: DIM
    integer(ip),          intent(in) :: KIND
    if( associated(xx) ) then
       my_lbound = lbound(xx,DIM=DIM)
    else
       my_lbound = 0_ip
    end if
  end function alya_lbound_r82

  integer(ip) function alya_lbound_r83(xx,DIM,KIND) result(my_lbound)
    real(8),    pointer, intent(in) :: xx(:,:,:)
    integer(ip),          intent(in) :: DIM
    integer(ip),          intent(in) :: KIND
    if( associated(xx) ) then
       my_lbound = lbound(xx,DIM=DIM)
    else
       my_lbound = 0_ip
    end if
  end function alya_lbound_r83

  integer(ip) function alya_lbound_lg1(xx,DIM,KIND) result(my_lbound)
    logical(lg), pointer, intent(in) :: xx(:)
    integer(ip),          intent(in) :: DIM
    integer(ip),          intent(in) :: KIND
    if( associated(xx) ) then
       my_lbound = lbound(xx,DIM=DIM)
    else
       my_lbound = 0_ip
    end if
  end function alya_lbound_lg1

  integer(ip) function alya_ubound_41(xx,DIM,KIND) result(my_ubound)
    integer(4),  pointer, intent(in) :: xx(:)
    integer(ip),          intent(in) :: DIM
    integer(ip),          intent(in) :: KIND
    if( associated(xx) ) then
       my_ubound = ubound(xx,DIM=DIM)
    else
       my_ubound = 0_4
    end if
  end function alya_ubound_41

  integer(ip) function alya_ubound_42(xx,DIM,KIND) result(my_ubound)
    integer(4),  pointer, intent(in) :: xx(:,:)
    integer(ip),          intent(in) :: DIM
    integer(ip),          intent(in) :: KIND
    if( associated(xx) ) then
       my_ubound = ubound(xx,DIM=DIM)
    else
       my_ubound = 0_4
    end if
  end function alya_ubound_42

  integer(ip) function alya_ubound_43(xx,DIM,KIND) result(my_ubound)
    integer(4),  pointer, intent(in) :: xx(:,:,:)
    integer(ip),          intent(in) :: DIM
    integer(ip),          intent(in) :: KIND
    if( associated(xx) ) then
       my_ubound = ubound(xx,DIM=DIM)
    else
       my_ubound = 0_4
    end if
  end function alya_ubound_43

  integer(ip) function alya_ubound_81(xx,DIM,KIND) result(my_ubound)
    integer(8),  pointer, intent(in) :: xx(:)
    integer(ip),          intent(in) :: DIM
    integer(ip),          intent(in) :: KIND
    if( associated(xx) ) then
       my_ubound = ubound(xx,DIM=DIM)
    else
       my_ubound = 0_4
    end if
  end function alya_ubound_81

  integer(ip) function alya_ubound_82(xx,DIM,KIND) result(my_ubound)
    integer(8),  pointer, intent(in) :: xx(:,:)
    integer(ip),          intent(in) :: DIM
    integer(ip),          intent(in) :: KIND
    if( associated(xx) ) then
       my_ubound = ubound(xx,DIM=DIM)
    else
       my_ubound = 0_4
    end if
  end function alya_ubound_82

  integer(ip) function alya_ubound_83(xx,DIM,KIND) result(my_ubound)
    integer(8),  pointer, intent(in) :: xx(:,:,:)
    integer(ip),          intent(in) :: DIM
    integer(ip),          intent(in) :: KIND
    if( associated(xx) ) then
       my_ubound = ubound(xx,DIM=DIM)
    else
       my_ubound = 0_4
    end if
  end function alya_ubound_83

  integer(ip) function alya_ubound_r41(xx,DIM,KIND) result(my_ubound)
    real(4),    pointer, intent(in) :: xx(:)
    integer(ip),          intent(in) :: DIM
    integer(ip),          intent(in) :: KIND
    if( associated(xx) ) then
       my_ubound = ubound(xx,DIM=DIM)
    else
       my_ubound = 0_4
    end if
  end function alya_ubound_r41

  integer(ip) function alya_ubound_r42(xx,DIM,KIND) result(my_ubound)
    real(4),    pointer, intent(in) :: xx(:,:)
    integer(ip),          intent(in) :: DIM
    integer(ip),          intent(in) :: KIND
    if( associated(xx) ) then
       my_ubound = ubound(xx,DIM=DIM)
    else
       my_ubound = 0_4
    end if
  end function alya_ubound_r42

  integer(ip) function alya_ubound_r43(xx,DIM,KIND) result(my_ubound)
    real(4),    pointer, intent(in) :: xx(:,:,:)
    integer(ip),          intent(in) :: DIM
    integer(ip),          intent(in) :: KIND
    if( associated(xx) ) then
       my_ubound = ubound(xx,DIM=DIM)
    else
       my_ubound = 0_4
    end if
  end function alya_ubound_r43

  integer(ip) function alya_ubound_r81(xx,DIM,KIND) result(my_ubound)
    real(8),    pointer, intent(in) :: xx(:)
    integer(ip),          intent(in) :: DIM
    integer(ip),          intent(in) :: KIND
    if( associated(xx) ) then
       my_ubound = ubound(xx,DIM=DIM)
    else
       my_ubound = 0_4
    end if
  end function alya_ubound_r81

  integer(ip) function alya_ubound_r82(xx,DIM,KIND) result(my_ubound)
    real(8),    pointer, intent(in) :: xx(:,:)
    integer(ip),          intent(in) :: DIM
    integer(ip),          intent(in) :: KIND
    if( associated(xx) ) then
       my_ubound = ubound(xx,DIM=DIM)
    else
       my_ubound = 0_4
    end if
  end function alya_ubound_r82

  integer(ip) function alya_ubound_r83(xx,DIM,KIND) result(my_ubound)
    real(8),    pointer, intent(in) :: xx(:,:,:)
    integer(ip),          intent(in) :: DIM
    integer(ip),          intent(in) :: KIND
    if( associated(xx) ) then
       my_ubound = ubound(xx,DIM=DIM)
    else
       my_ubound = 0_4
    end if
  end function alya_ubound_r83

  integer(ip) function alya_ubound_lg1(xx,DIM,KIND) result(my_ubound)
    logical(lg), pointer, intent(in) :: xx(:)
    integer(ip),          intent(in) :: DIM
    integer(ip),          intent(in) :: KIND
    if( associated(xx) ) then
       my_ubound = ubound(xx,DIM=DIM)
    else
       my_ubound = 0_4
    end if
  end function alya_ubound_lg1

end module mod_std
!> @}
