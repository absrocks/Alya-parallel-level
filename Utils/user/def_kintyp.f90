module def_kintyp
  !-----------------------------------------------------------------------
  !****f* defmod/def_kintyp
  ! NAME
  !   def_kintyp
  ! DESCRIPTION
  !   Definition of kinds and types. 
  !   "The range of the default integers is not specified in the language
  !   but on a computer with a word size of n bits, is often from 
  !   -2^{n-1} to +2^{n-1}-1. Thus on a 32-bit computer the range is
  !   often -2.14*10^9 to +2.14*10^9."
  !   M. Metclaf and J. Reid, FORTRAN 90/95 explained, 2nd edition.
  !
  !   Defaults are:
  !   Integers: 4-bytes 
  !   Reals:    8-bytes
  !  
  !***
  !-----------------------------------------------------------------------

  !----------------------------------------------------------------------
  !
  ! Symbolc names for integers, reals and logicals
  !
  !----------------------------------------------------------------------
  !
  ! Symbolic names for integers
  !
  integer, parameter  :: ip = 4             ! 4-byte integer
  integer, parameter  :: rp = 8             ! Double precision 
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
  type i1pp
     integer(ip)          :: n
     integer(ip), pointer :: l(:)
  end type i1pp
  type i1pi1p
     type(i1p),   pointer :: l(:)
  end type i1pi1p

end module def_kintyp
