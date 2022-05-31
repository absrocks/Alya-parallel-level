!-----------------------------------------------------------------------
!> @defgroup Kinds_and_types
!> Kinds ands types of Alya
!> @{
!> @file    def_kintyp.f90
!> @author  houzeaux
!> @date    2018-12-28
!> @brief   Definition of kinds and types.
!> @details Definition of kinds and types.
!>          "The range of the default integers is not specified in the language
!>          but on a computer with a word size of n bits, is often from
!>          -2^{n-1} to +2^{n-1}-1. Thus on a 32-bit computer the range is
!>          often -2.14*10^9 to +2.14*10^9."
!>          M. Metclaf and J. Reid, FORTRAN 90/95 explained, 2nd edition.
!>
!>          Defaults are:
!>          Integers: 4-bytes
!>          Reals:    8-bytes
!>
!-----------------------------------------------------------------------

module def_kintyp_basic

  !----------------------------------------------------------------------
  !
  ! Symbolc names for integers, reals and logicals
  !
  !----------------------------------------------------------------------
  !
  ! Symbolic names for integers
  !
#ifdef I8
  integer, parameter  :: ip = 8             ! 8-byte integer
#else
  integer, parameter  :: ip = 4             ! 4-byte integer
#endif
  !
  ! Symbolic names for reals
  !
#ifdef R4
  integer, parameter  :: rp = 4             ! Simple precision
#elif R16
  integer, parameter  :: rp = 16            ! Very high precision - Beware does not work with MPI
#else
  integer, parameter  :: rp = 8             ! Double precision
#endif
  !
  ! Symbolic name for kind type of default logical
  !
  integer, parameter  :: lg = kind(.true.)

  !----------------------------------------------------------------------
  !
  ! General types
  !
  !----------------------------------------------------------------------

  type i1p
     integer(ip), pointer :: l(:)
  end type i1p
  type i2p
     integer(ip), pointer :: l(:,:)
  end type i2p
  type i3p
     integer(ip), pointer :: l(:,:,:)
  end type i3p
  type r1p
     real(rp),    pointer :: a(:)
  end type r1p
  type r2p
     real(rp),    pointer :: a(:,:)
  end type r2p
  type r3p
     real(rp),    pointer :: a(:,:,:)
  end type r3p
  type r4p
     real(rp),    pointer :: a(:,:,:,:)
  end type r4p
  
  type, extends(i1p)      :: i1pp
     integer(ip)          :: n
  end type i1pp
  
  !----------------------------------------------------------------------
  !
  ! Sparse matrix 
  !
  !----------------------------------------------------------------------
  
  type spmat
     integer(ip)          :: ndof
     integer(ip)          :: nrows
     integer(ip)          :: ncols
     integer(ip), pointer :: iA(:)
     integer(ip), pointer :: jA(:)
     real(rp),    pointer :: vA(:,:,:)
  end type spmat
  
  !----------------------------------------------------------------------
  !
  ! Element interpolation
  !
  !----------------------------------------------------------------------
  
  type typ_interp
     integer(ip), pointer  :: lelem(:)   ! List of element
     real(rp),    pointer  :: shapf(:,:) ! List of shape functions
  end type typ_interp

end module def_kintyp_basic
!> @}
