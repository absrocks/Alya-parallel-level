subroutine elmtyp()
  !-----------------------------------------------------------------------
  !****f* Domain/elmtyp
  ! NAME
  !    cderda
  ! DESCRIPTION
  !    This routine defines the different types of elements:
  !
  !    NNODE ... number fo nodes
  !
  !                 1D           2D           3D
  !                 ----------------------------------
  !    LTOPO ... -1 Lines          -
  !               0    -    Quadrilateral  Hexahedra
  !               1    -    Triangle       Tetrahedra
  !               2    -           -       Pentahedra
  !               3    -           -       Pyramid
  !
  !    LDIME ... 1,2,3: Dimension of the element
  !    LLAPL ... 1 if Hessian matrix should be considered
  !              0 otherwise
  !    CENAM ... Element name (3 letters + 2 figures=# nodes)
  !    CENAL ... Alternative element name (for Ensight format)
  !
  ! USED BY
  !    Domain
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_elmtyp
  implicit none
  integer(ip) :: ielty

  nnode = 0
  !
  ! 0D elements
  !
  nnode(POINT)  =  1        !   o
  ltopo(POINT)  = -2 
  ldime(POINT)  =  0
  llapl(POINT)  =  0
  lorde(POINT)  =  1
  nface(BAR02)  =  0
  cenam(POINT)  =  'POINT'
  cenal(POINT)  =  'poin'
  cetop(POINT)  =  'Point'
  cepos(BAR02)  =  'LagrLigne02'
  !
  ! 1D elements
  !
  nnode(BAR02)  =  2        !   o-----o
  ltopo(BAR02)  = -1 
  ldime(BAR02)  =  1
  llapl(BAR02)  =  0
  lorde(BAR02)  =  1
  nface(BAR02)  =  2
  cenam(BAR02)  =  'BAR02'
  cenal(BAR02)  =  'bar2'
  cetop(BAR02)  =  'Linear'

  nnode(BAR03)  =  3      
  ltopo(BAR03)  = -1 
  ldime(BAR03)  =  1 
  llapl(BAR03)  =  1
  lorde(BAR03)  =  2
  nface(BAR03)  =  2
  cenam(BAR03)  =  'BAR03'
  cenal(BAR03)  =  'bar3'
  cetop(BAR03)  =  'Linear'
  cepos(BAR03)  =  'LagrLigne03'
 
  nnode(BAR04)  =  4      
  ltopo(BAR04)  = -1 
  ldime(BAR04)  =  1 
  llapl(BAR04)  =  1
  lorde(BAR04)  =  3
  nface(BAR04)  =  2
  cenam(BAR04)  =  'BAR04'
  cenal(BAR04)  =  'bar4'
  cetop(BAR04)  =  'Linear'
  !
  ! 2D elements 
  !
  nnode(TRI03)  =  3        !   o-----o
  ltopo(TRI03)  =  1        !    \   /
  ldime(TRI03)  =  2        !     \ /
  llapl(TRI03)  =  0        !      o
  lorde(TRI03)  =  1 
  nface(TRI03)  =  3
  cenam(TRI03)  =  'TRI03'
  cenal(TRI03)  =  'tria3'
  cetop(TRI03)  =  'Triangle'
  cepos(TRI03)  =  'LagrTrian03'

  nnode(TRI06)  =  6     
  ltopo(TRI06)  =  1     
  ldime(TRI06)  =  2     
  llapl(TRI06)  =  1     
  lorde(TRI06)  =  2
  nface(TRI06)  =  3
  cenam(TRI06)  =  'TRI06'  
  cenal(TRI06)  =  'tria6'  
  cetop(TRI06)  =  'Triangle'
  cepos(TRI06)  =  'LagrTrian06'

  nnode(QUA04)  =  4        !   o-----o
  ltopo(QUA04)  =  0        !   |     |
  ldime(QUA04)  =  2        !   |     |
  llapl(QUA04)  =  1        !   o-----o
  lorde(QUA04)  =  1
  nface(QUA04)  =  4
  cenam(QUA04)  =  'QUA04'
  cenal(QUA04)  =  'quad4'
  cetop(QUA04)  =  'Quadrilateral'
  cepos(QUA04)  =  'LagrQuadr04'

  nnode(QUA08)  =  8
  ltopo(QUA08)  =  0
  ldime(QUA08)  =  2
  llapl(QUA08)  =  1
  lorde(QUA08)  =  2
  nface(QUA08)  =  4
  cenam(QUA08)  =  'QUA08'
  cenal(QUA08)  =  'quad8'
  cetop(QUA08)  =  'Quadrilateral'
  cepos(QUA08)  =  'LagrQuadr08'
       
  nnode(QUA09)  =  9
  ltopo(QUA09)  =  0
  ldime(QUA09)  =  2
  llapl(QUA09)  =  1
  lorde(QUA09)  =  2
  nface(QUA09)  =  4
  cenam(QUA09)  =  'QUA09'
  cenal(QUA09)  =  'quad9'
  cetop(QUA09)  =  'Quadrilateral'
  cepos(QUA09)  =  'LagrQuadr09'

  nnode(QUA16)  =  16
  ltopo(QUA16)  =  0
  ldime(QUA16)  =  2
  llapl(QUA16)  =  1
  lorde(QUA16)  =  2
  nface(QUA16)  =  4
  cenam(QUA16)  =  'QUA16'
  cenal(QUA16)  =  'quad16'
  cetop(QUA16)  =  'Quadrilateral'
  !
  ! 3D elements
  !
  nnode(TET04)  =  4        !    o
  ltopo(TET04)  =  1        !   / \
  ldime(TET04)  =  3        !  o-|-o
  llapl(TET04)  =  0        !    o
  lorde(TET04)  =  1
  nface(TET04)  =  4
  cenam(TET04)  =  'TET04'
  cenal(TET04)  =  'tetra4'
  cetop(TET04)  =  'Tetrahedra'
  cepos(TET04)  =  'LagrTetra04'
  needg(TET04)  =  6
  nnode(TET10)  = 10 
  ltopo(TET10)  =  1 
  ldime(TET10)  =  3 
  llapl(TET10)  =  1
  lorde(TET10)  =  2
  nface(TET10)  =  4
  cenam(TET10)  =  'TET10'
  cenal(TET10)  =  'tetra10'
  cetop(TET10)  =  'Tetrahedra'
  cepos(TET10)  =  'LagrTetra10'

  nnode(PYR05)  =  5        !     o
  ltopo(PYR05)  =  3        !    / \
  ldime(PYR05)  =  3        !   o---o
  llapl(PYR05)  =  1        !   o---o
  lorde(PYR05)  =  1
  nface(PYR05)  =  5
  cenam(PYR05)  =  'PYR05'
  cenal(PYR05)  =  'pyra5'
  cetop(PYR05)  =  'Pyramid'
  cepos(PYR05)  =  'LagrPyram05'

  nnode(PYR14)  = 14 
  ltopo(PYR14)  =  3 
  ldime(PYR14)  =  3 
  llapl(PYR14)  =  1
  lorde(PYR14)  =  2
  nface(PYR14)  =  5
  cenam(PYR14)  =  'PYR14'
  cenal(PYR14)  =  'pyra14'
  cetop(PYR14)  =  'Pyramid'
  cepos(PYR14)  =  'LagrPyram14'

  nnode(PEN06)  =  6        !   o---o  
  ltopo(PEN06)  =  2        !     o
  ldime(PEN06)  =  3        !   o-|-o
  llapl(PEN06)  =  0        !     o
  lorde(PEN06)  =  1
  nface(PEN06)  =  5
  cenam(PEN06)  =  'PEN06'
  cenal(PEN06)  =  'penta6'
  cetop(PEN06)  =  'Prism'
  cepos(PEN06)  =  'LagrPrism06'

  nnode(PEN15)  = 15 
  ltopo(PEN15)  =  2 
  ldime(PEN15)  =  3 
  llapl(PEN15)  =  1
  lorde(PEN15)  =  2
  nface(PEN15)  =  6
  cenam(PEN15)  =  'PEN15'
  cenal(PEN15)  =  'penta15'
  cetop(PEN15)  =  'Prism'
  cepos(PEN15)  =  'LagrPrism15'

  nnode(PEN18)  = 18 
  ltopo(PEN18)  =  2 
  ldime(PEN18)  =  3 
  llapl(PEN18)  =  1
  lorde(PEN18)  =  2
  nface(PEN18)  =  6
  cenam(PEN18)  =  'PEN18'
  cenal(PEN18)  =  'penta18'
  cetop(PEN18)  =  'Prism'
  cepos(PEN18)  =  'LagrPrism18'

  nnode(HEX08)  =  8       !     o---o  
  ltopo(HEX08)  =  0       !   o---o |  
  ldime(HEX08)  =  3       !   | o-|-o  
  llapl(HEX08)  =  1       !   o---o   
  lorde(HEX08)  =  1
  nface(HEX08)  =  6
  cenam(HEX08)  =  'HEX08'
  cenal(HEX08)  =  'hexa8'
  cetop(HEX08)  =  'Hexahedra'
  cepos(HEX08)  =  'LagrHexae08'

  nnode(HEX20)  = 20  
  ltopo(HEX20)  =  0  
  ldime(HEX20)  =  3  
  llapl(HEX20)  =  1 
  lorde(HEX20)  =  2
  nface(HEX20)  =  6
  cenam(HEX20)  =  'HEX20'
  cenal(HEX20)  =  'hexa20'
  cetop(HEX20)  =  'Hexahedra'
  cepos(HEX20)  =  'LagrHexae20'

  nnode(HEX27)  = 27 
  ltopo(HEX27)  =  0
  ldime(HEX27)  =  3
  llapl(HEX27)  =  1
  lorde(HEX27)  =  2
  nface(HEX27)  =  6
  cenam(HEX27)  =  'HEX27'
  cenal(HEX27)  =  'hexa27'
  cetop(HEX27)  =  'Hexahedra'
  cepos(HEX27)  =  'LagrHexae27'

  nnode(SHELL)  =  3
  ltopo(SHELL)  =  1 
  ldime(SHELL)  =  3
  llapl(SHELL)  =  0
  lorde(SHELL)  =  2
  nface(SHELL)  =  3
  cenam(SHELL)  =  'SHELL'
  cenal(SHELL)  =  'shell'
  cetop(SHELL)  =  'Triangle'
  cepos(SHELL)  =  'Triangle'

  nnode(BAR3D)  =  2      
  ltopo(BAR3D)  = -1 
  ldime(BAR3D)  =  1 
  llapl(BAR3D)  =  0
  lorde(BAR3D)  =  1
  nface(BAR3D)  =  2
  cenam(BAR3D)  =  'BAR3D'
  cenal(BAR3D)  =  'bar3'
  cetop(BAR3D)  =  'Linear'
  cepos(BAR3D)  =  'LagrLigne03'

  do ielty = 1,nelty
     nnode(-ielty) = nnode(ielty)
  end do

end subroutine elmtyp
