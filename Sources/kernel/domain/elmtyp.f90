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
  use def_elmtyp
  use def_domain
  use def_postpr
!  use mod_elmgeo, only     :  list_edges_HEX27,list_edges_QUA09,element_type ! this taht would have been nice does not work due to the order in which things are called

  implicit none
  integer(ip) :: pelty,inode,pnode

  !---------------------------------------------------------------------
  !
  ! 0D elements
  !
  !---------------------------------------------------------------------

  nnode(POINT)  =  1    
  ltopo(POINT)  = -2 
  ldime(POINT)  =  0
  llapl(POINT)  =  0
  lorde(POINT)  =  1
  nface(POINT)  =  0
  cenam(POINT)  =  'POINT'
  cenal(POINT)  =  'poin'
  cetop(POINT)  =  'Linear'
  cepos(BAR02)  =  'LagrLigne02'

  !---------------------------------------------------------------------
  !
  ! 1D elements
  !
  !---------------------------------------------------------------------
  !
  ! BAR02
  !
  nnode(BAR02)     =  2   
  ltopo(BAR02)     = -1 
  ldime(BAR02)     =  1
  llapl(BAR02)     =  0
  lorde(BAR02)     =  1
  nface(BAR02)     =  2
  cenam(BAR02)     =  'BAR02'
  cenal(BAR02)     =  'bar2'
  cetop(BAR02)     =  'Linear'
  needg(BAR02)     =  1
  leedg(1,1,BAR02) = 1
  leedg(2,1,BAR02) = 2
  !
  ! BAR03
  !
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
  leedg(1,1,BAR03) = 1
  leedg(2,1,BAR03) = 2
  ! 
  ! BAR04
  !
  nnode(BAR04)  =  4      
  ltopo(BAR04)  = -1 
  ldime(BAR04)  =  1 
  llapl(BAR04)  =  1
  lorde(BAR04)  =  3
  nface(BAR04)  =  2
  cenam(BAR04)  =  'BAR04'
  cenal(BAR04)  =  'bar4'
  cetop(BAR04)  =  'Linear'
  leedg(1,1,BAR04) = 1
  leedg(2,1,BAR04) = 2

  !---------------------------------------------------------------------
  !
  ! 2D elements 
  !
  !---------------------------------------------------------------------
  !
  ! TRI03
  !
  nnode(TRI03)  =  3      
  ltopo(TRI03)  =  1     
  ldime(TRI03)  =  2   
  llapl(TRI03)  =  0     
  lorde(TRI03)  =  1 
  nface(TRI03)  =  3
  cenam(TRI03)  =  'TRI03'
  cenal(TRI03)  =  'tria3'
  cetop(TRI03)  =  'Triangle'
  cepos(TRI03)  =  'LagrTrian03'
  needg(TRI03)  =  3
  leedg(1,1,TRI03) = 1
  leedg(2,1,TRI03) = 2
  leedg(1,2,TRI03) = 2
  leedg(2,2,TRI03) = 3
  leedg(1,3,TRI03) = 3
  leedg(2,3,TRI03) = 1
  !
  ! TRI06
  !
  nnode(TRI06)  =  6     
  ltopo(TRI06)  =  1     
  ldime(TRI06)  =  2     
  llapl(TRI06)  =  1     
  lorde(TRI06)  =  2
  nface(TRI06)  =  3
  needg(TRI06)  =  3
  cenam(TRI06)  =  'TRI06'  
  cenal(TRI06)  =  'tria6'  
  cetop(TRI06)  =  'Triangle'
  cepos(TRI06)  =  'LagrTrian06'
  leedg(1,1,TRI06) = 1
  leedg(2,1,TRI06) = 2
  leedg(1,2,TRI06) = 2
  leedg(2,2,TRI06) = 3
  leedg(1,3,TRI06) = 3
  leedg(2,3,TRI06) = 1
  !
  ! TRI10
  !
  nnode(TRI10)  =  10
  ltopo(TRI10)  =  1     
  ldime(TRI10)  =  2     
  llapl(TRI10)  =  0    
  lorde(TRI10)  =  3
  nface(TRI10)  =  3
  needg(TRI10)  =  3
  cenam(TRI10)  =  'TRI10'  
  cenal(TRI10)  =  'tria10'  
  cetop(TRI10)  =  'Triangle'
  cepos(TRI10)  =  'LagrTrian10'
  leedg(1,1,TRI10) = 1
  leedg(2,1,TRI10) = 2
  leedg(1,2,TRI10) = 2
  leedg(2,2,TRI10) = 3
  leedg(1,3,TRI10) = 3
  leedg(2,3,TRI10) = 1
  !
  ! QUA04
  !
  nnode(QUA04)  =  4       
  ltopo(QUA04)  =  0       
  ldime(QUA04)  =  2       
  llapl(QUA04)  =  1       
  lorde(QUA04)  =  1
  nface(QUA04)  =  4
  cenam(QUA04)  =  'QUA04'
  cenal(QUA04)  =  'quad4'
  cetop(QUA04)  =  'Quadrilateral'
  cepos(QUA04)  =  'LagrQuadr04'
  needg(QUA04)  =  4
  leedg(1,1,QUA04) = 1
  leedg(2,1,QUA04) = 2
  leedg(1,2,QUA04) = 2
  leedg(2,2,QUA04) = 3
  leedg(1,3,QUA04) = 3
  leedg(2,3,QUA04) = 4
  leedg(1,4,QUA04) = 4
  leedg(2,4,QUA04) = 1
  !
  ! QUA08
  !
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
  needg(QUA08)  =  4
  leedg(1,1,QUA08) = 1
  leedg(2,1,QUA08) = 2
  leedg(1,2,QUA08) = 2
  leedg(2,2,QUA08) = 3
  leedg(1,3,QUA08) = 3
  leedg(2,3,QUA08) = 4
  leedg(1,4,QUA08) = 4
  leedg(2,4,QUA08) = 1
  !
  ! QUA09
  !
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
  needg(QUA09)  =  4 ! element_type(QUA09) % number_edges - THIS that looked nice cna not be done because when this subroutine is called the value has not been defined
  leedg(1:2,1:4,QUA09) = reshape ( (/ 1,2,     2,3,     3,4,     4,1                             /), (/2,4 /) )
  !
  ! QUA16
  !
  nnode(QUA16)  =  16
  ltopo(QUA16)  =  0
  ldime(QUA16)  =  2
  llapl(QUA16)  =  1
  lorde(QUA16)  =  3
  nface(QUA16)  =  4
  cenam(QUA16)  =  'QUA16'
  cenal(QUA16)  =  'quad16'
  cetop(QUA16)  =  'Quadrilateral'

  !---------------------------------------------------------------------
  !
  ! 3D elements
  !
  !---------------------------------------------------------------------
  !
  ! TET04
  !
  nnode(TET04)  =  4      
  ltopo(TET04)  =  1      
  ldime(TET04)  =  3      
  llapl(TET04)  =  0      
  lorde(TET04)  =  1
  nface(TET04)  =  4
  cenam(TET04)  =  'TET04'
  cenal(TET04)  =  'tetra4'
  cetop(TET04)  =  'Tetrahedra'
  cepos(TET04)  =  'LagrTetra04'
  needg(TET04)  =  6
  leedg(1,1,TET04) = 1
  leedg(2,1,TET04) = 2
  leedg(1,2,TET04) = 1
  leedg(2,2,TET04) = 3
  leedg(1,3,TET04) = 1
  leedg(2,3,TET04) = 4
  leedg(1,4,TET04) = 2
  leedg(2,4,TET04) = 3
  leedg(1,5,TET04) = 4
  leedg(2,5,TET04) = 2
  leedg(1,6,TET04) = 3
  leedg(2,6,TET04) = 4
  !
  ! TET10
  !
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
  needg(TET10)  =  6
  leedg(1,1,TET10) = 1
  leedg(2,1,TET10) = 2
  leedg(1,2,TET10) = 1
  leedg(2,2,TET10) = 3
  leedg(1,3,TET10) = 1
  leedg(2,3,TET10) = 4
  leedg(1,4,TET10) = 2
  leedg(2,4,TET10) = 3
  leedg(1,5,TET10) = 4
  leedg(2,5,TET10) = 2
  leedg(1,6,TET10) = 3
  leedg(2,6,TET10) = 4
  !
  ! TET20
  !
  nnode(TET20)  = 20 
  ltopo(TET20)  =  1 
  ldime(TET20)  =  3 
  llapl(TET20)  =  0
  lorde(TET20)  =  3
  nface(TET20)  =  4
  cenam(TET20)  =  'TET20'
  cenal(TET20)  =  'tetra20'
  cetop(TET20)  =  'Tetrahedra'
  cepos(TET20)  =  'LagrTetra20'
  needg(TET20)  =  6
  leedg(1,1,TET20) = 1
  leedg(2,1,TET20) = 2
  leedg(1,2,TET20) = 1
  leedg(2,2,TET20) = 3
  leedg(1,3,TET20) = 1
  leedg(2,3,TET20) = 4
  leedg(1,4,TET20) = 2
  leedg(2,4,TET20) = 3
  leedg(1,5,TET20) = 4
  leedg(2,5,TET20) = 2
  leedg(1,6,TET20) = 3
  leedg(2,6,TET20) = 4
  !
  ! PYR05
  !
  nnode(PYR05)  =  5       
  ltopo(PYR05)  =  3       
  ldime(PYR05)  =  3       
  llapl(PYR05)  =  1       
  lorde(PYR05)  =  1
  nface(PYR05)  =  5
  cenam(PYR05)  =  'PYR05'
  cenal(PYR05)  =  'pyra5'
  cetop(PYR05)  =  'Pyramid'
  cepos(PYR05)  =  'LagrPyram05'
  needg(PYR05)  =  8
  leedg(1,1,PYR05) = 1
  leedg(2,1,PYR05) = 2
  leedg(1,2,PYR05) = 2
  leedg(2,2,PYR05) = 3
  leedg(1,3,PYR05) = 3
  leedg(2,3,PYR05) = 4
  leedg(1,4,PYR05) = 4
  leedg(2,4,PYR05) = 1
  leedg(1,5,PYR05) = 1
  leedg(2,5,PYR05) = 5
  leedg(1,6,PYR05) = 2
  leedg(2,6,PYR05) = 5
  leedg(1,7,PYR05) = 3
  leedg(2,7,PYR05) = 5
  leedg(1,8,PYR05) = 4
  leedg(2,8,PYR05) = 5
  !
  ! PYR14
  !
  nnode(PYR14)  = 14 
  ltopo(PYR14)  =  3 
  ldime(PYR14)  =  3 
  llapl(PYR14)  =  1
  lorde(PYR14)  =  2
  nface(PYR14)  =  0
  cenam(PYR14)  =  'PYR14'
  cenal(PYR14)  =  'pyra14'
  cetop(PYR14)  =  'Pyramid'
  cepos(PYR14)  =  'LagrPyram14'
  !
  ! PEN06
  !
  nnode(PEN06)  =  6       
  ltopo(PEN06)  =  2       
  ldime(PEN06)  =  3       
  llapl(PEN06)  =  0       
  lorde(PEN06)  =  1
  nface(PEN06)  =  5
  cenam(PEN06)  =  'PEN06'
  cenal(PEN06)  =  'penta6'
  cetop(PEN06)  =  'Prism'
  cepos(PEN06)  =  'LagrPrism06'
  needg(PEN06)  =  9
  leedg(1,1,PEN06) = 1
  leedg(2,1,PEN06) = 2
  leedg(1,2,PEN06) = 2
  leedg(2,2,PEN06) = 3
  leedg(1,3,PEN06) = 3
  leedg(2,3,PEN06) = 1
  leedg(1,4,PEN06) = 4
  leedg(2,4,PEN06) = 5
  leedg(1,5,PEN06) = 5
  leedg(2,5,PEN06) = 6
  leedg(1,6,PEN06) = 6
  leedg(2,6,PEN06) = 4
  leedg(1,7,PEN06) = 1
  leedg(2,7,PEN06) = 4
  leedg(1,8,PEN06) = 3
  leedg(2,8,PEN06) = 6
  leedg(1,9,PEN06) = 2
  leedg(2,9,PEN06) = 5
  !
  ! PEN15
  !
  nnode(PEN15)  = 15
  ltopo(PEN15)  =  2
  ldime(PEN15)  =  3
  llapl(PEN15)  =  1
  lorde(PEN15)  =  2
  nface(PEN15)  =  5
  cenam(PEN15)  =  'PEN15'
  cenal(PEN15)  =  'penta15'
  cetop(PEN15)  =  'Prism'
  cepos(PEN15)  =  'LagrPrism15'
  needg(PEN15)  =  9
  leedg(1,1,PEN15) = 1
  leedg(2,1,PEN15) = 2
  leedg(1,2,PEN15) = 2
  leedg(2,2,PEN15) = 3
  leedg(1,3,PEN15) = 3
  leedg(2,3,PEN15) = 1
  leedg(1,4,PEN15) = 4
  leedg(2,4,PEN15) = 5
  leedg(1,5,PEN15) = 5
  leedg(2,5,PEN15) = 6
  leedg(1,6,PEN15) = 6
  leedg(2,6,PEN15) = 4
  leedg(1,7,PEN15) = 1
  leedg(2,7,PEN15) = 4
  leedg(1,8,PEN15) = 3
  leedg(2,8,PEN15) = 6
  leedg(1,9,PEN15) = 2
  leedg(2,9,PEN15) = 5
  !
  ! PEN18
  !
  nnode(PEN18)  = 18
  ltopo(PEN18)  =  2
  ldime(PEN18)  =  3
  llapl(PEN18)  =  1
  lorde(PEN18)  =  2
  nface(PEN18)  =  5
  cenam(PEN18)  =  'PEN18'
  cenal(PEN18)  =  'penta18'
  cetop(PEN18)  =  'Prism'
  cepos(PEN18)  =  'LagrPrism18'
  needg(PEN18)  =  9
  leedg(1,1,PEN18) = 1
  leedg(2,1,PEN18) = 2
  leedg(1,2,PEN18) = 2
  leedg(2,2,PEN18) = 3
  leedg(1,3,PEN18) = 3
  leedg(2,3,PEN18) = 1
  leedg(1,4,PEN18) = 4
  leedg(2,4,PEN18) = 5
  leedg(1,5,PEN18) = 5
  leedg(2,5,PEN18) = 6
  leedg(1,6,PEN18) = 6
  leedg(2,6,PEN18) = 4
  leedg(1,7,PEN18) = 1
  leedg(2,7,PEN18) = 4
  leedg(1,8,PEN18) = 3
  leedg(2,8,PEN18) = 6
  leedg(1,9,PEN18) = 2
  leedg(2,9,PEN18) = 5
  !
  ! HEX08
  !
  nnode(HEX08)  =  8     
  ltopo(HEX08)  =  0     
  ldime(HEX08)  =  3     
  llapl(HEX08)  =  1     
  lorde(HEX08)  =  1
  nface(HEX08)  =  6
  cenam(HEX08)  =  'HEX08'
  cenal(HEX08)  =  'hexa8'
  cetop(HEX08)  =  'Hexahedra'
  cepos(HEX08)  =  'LagrHexae08'
  needg(HEX08)  =  12
  leedg(1, 1,HEX08) = 1
  leedg(2, 1,HEX08) = 2
  leedg(1, 2,HEX08) = 1
  leedg(2, 2,HEX08) = 4
  leedg(1, 3,HEX08) = 1
  leedg(2, 3,HEX08) = 5
  leedg(1, 4,HEX08) = 2
  leedg(2, 4,HEX08) = 3
  leedg(1, 5,HEX08) = 2
  leedg(2, 5,HEX08) = 6
  leedg(1, 6,HEX08) = 3
  leedg(2, 6,HEX08) = 4
  leedg(1, 7,HEX08) = 3
  leedg(2, 7,HEX08) = 7
  leedg(1, 8,HEX08) = 4
  leedg(2, 8,HEX08) = 8
  leedg(1, 9,HEX08) = 5
  leedg(2, 9,HEX08) = 6
  leedg(1,10,HEX08) = 5
  leedg(2,10,HEX08) = 8
  leedg(1,11,HEX08) = 6
  leedg(2,11,HEX08) = 7
  leedg(1,12,HEX08) = 7
  leedg(2,12,HEX08) = 8
  !!
  !! HEX20
  !!
  !nnode(HEX20)  = 20  
  !ltopo(HEX20)  =  0  
  !ldime(HEX20)  =  3  
  !llapl(HEX20)  =  1 
  !lorde(HEX20)  =  2
  !nface(HEX20)  =  0
  !cenam(HEX20)  =  'HEX20'
  !cenal(HEX20)  =  'hexa20'
  !cetop(HEX20)  =  'Hexahedra'
  !cepos(HEX20)  =  'LagrHexae20'
  !
  ! HEX27
  !
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
  needg(HEX27)  =  12 ! element_type(HEX27) % number_edges 
  leedg(1:2,1:12,HEX27) = reshape ( (/ 1,2, 1,4, 1,5, 2,3, 2,6, 3,4, 3,7, 4,8, 5,6, 5,8, 6,7, 7,8 /), (/2,12/) )


  !---------------------------------------------------------------------
  !
  ! 3D 1D/2D-elements
  !
  !---------------------------------------------------------------------
  !
  ! SHELL
  !
  nnode(SHELL)     = 3    
  ltopo(SHELL)     = 1    
  ldime(SHELL)     = 2    
  llapl(SHELL)     = 0    
  lorde(SHELL)     = 1 
  nface(SHELL)     = 3
  cenam(SHELL)     = 'SHELL'
  cenal(SHELL)     = 'shell'
  cetop(SHELL)     = 'Triangle'
  cepos(SHELL)     = 'LagrTrian03'
  needg(SHELL)     = 3
  leedg(1,1,SHELL) = 1
  leedg(2,1,SHELL) = 2
  leedg(1,2,SHELL) = 2
  leedg(2,2,SHELL) = 3
  leedg(1,3,SHELL) = 3
  leedg(2,3,SHELL) = 1
  !
  ! BAR03
  !
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
  needg(BAR3D)     = 1
  leedg(1,1,BAR3D) = 1
  leedg(2,1,BAR3D) = 2

  !---------------------------------------------------------------------
  !
  ! Derived parameters
  !
  !---------------------------------------------------------------------
  !
  ! Next element node
  !
  do pelty = 1,nelty
     pnode = nnode(pelty)
     if( pnode > 0 ) then
        do inode = 1,pnode-1
           lenex(inode,pelty) = inode + 1
        end do
        lenex(pnode,pelty) = 1
     end if
  end do
  !
  ! Mirror NNODE
  !
  do pelty = 1,nelty
     nnode(-pelty) = nnode(pelty)
  end do
  mface = maxval(nface)

end subroutine elmtyp
