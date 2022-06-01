module def_elmtyp
  !-----------------------------------------------------------------------
  !****f* defmod/def_elmtyp
  ! NAME 
  !    def_elmtyp
  ! DESCRIPTION
  !    Elements available in Alya
  !    The original numbering followed CGNS standrad:
  !    http://www.grc.nasa.gov/WWW/cgns/sids/conv.html
  !    As new elements have been introduced, CGNS standard is no 
  !    longer used.
  !    The numbering is:
  !    1D elements:  2 ->  9
  !    2D elements: 10 -> 29
  !    3D elements: 30 -> 50
  ! USES
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only   :  ip

  integer(ip), parameter :: POINT =    1 ! 0D
  integer(ip), parameter :: BAR02 =    2 ! 1D
  integer(ip), parameter :: BAR03 =    3 ! 1D 
  integer(ip), parameter :: BAR04 =    4 ! 1D 
  integer(ip), parameter :: TRI03 =   10 ! 2D 
  integer(ip), parameter :: TRI06 =   11 ! 2D 
  integer(ip), parameter :: QUA04 =   12 ! 2D 
  integer(ip), parameter :: QUA08 =   13 ! 2D 
  integer(ip), parameter :: QUA09 =   14 ! 2D 
  integer(ip), parameter :: QUA16 =   15 ! 2D 
  integer(ip), parameter :: TET04 =   30 ! 3D 
  integer(ip), parameter :: TET10 =   31 ! 3D 
  integer(ip), parameter :: PYR05 =   32 ! 3D 
  integer(ip), parameter :: PYR14 =   33 ! 3D 
  integer(ip), parameter :: PEN06 =   34 ! 3D 
  integer(ip), parameter :: PEN15 =   35 ! 3D 
  integer(ip), parameter :: PEN18 =   36 ! 3D 
  integer(ip), parameter :: HEX08 =   37 ! 3D 
  integer(ip), parameter :: HEX20 =   38 ! 3D 
  integer(ip), parameter :: HEX27 =   39 ! 3D 
  integer(ip), parameter :: HEX64 =   40 ! 3D 
  integer(ip), parameter :: SHELL =   51 ! 3D 
  integer(ip), parameter :: BAR3D =   52 ! 3D 

  integer(ip), parameter :: DDESS = -101 ! DD 
  integer(ip), parameter :: DDNAT = -102 ! DD
  integer(ip), parameter :: DDEXT = -103 ! DD
  !
  ! Element characteristics
  !
  integer(ip), parameter :: ELFEM =    0 ! Finite element
  integer(ip), parameter :: ELEXT =    1 ! Extension finite element
  integer(ip), parameter :: ELHOL =    2 ! Hole
  integer(ip), parameter :: ELCNT =    3 ! Contact element
  !
  ! Node characteristics
  !
  integer(ip), parameter :: NOFEM =    0 ! Normal node
  integer(ip), parameter :: NOEXT =    1 ! Extension  node
  integer(ip), parameter :: NOHOL =    2 ! Hole node
  integer(ip), parameter :: NOCNT =    3 ! Contact node


end module def_elmtyp

