*if(ndime==2)
*set elems(Triangle)  
*add elems(Quadrilateral)
*else
*set elems(Tetrahedra)  
*add elems(Hexahedra)
*add elems(Prisma)
*endif
  ELEMENTS
*if(ndime==2)
*set var ielem=0
*set Cond DOMAIN_element_2D       *elems
*loop elems *onlyincond *NoCanRepeat
*set var entit=cond(1,int)
*set var ielem=ielem+1
*set var kelem=elemsNum
   *ielem *entit
*end
*else
*set var ielem=0
*loop elems
*set var ielem=ielem+1
*set var kelem=elemsNum
   *ielem 1
*end
*endif
  END_ELEMENTS
  BOUNDARIES
*set var iboun=0
*if(ndime==2)
*set Cond DOMAIN_boundary_2D       *elems
*loop elems *onlyincond *NoCanRepeat
*set var entit=cond(1,int)
*set var iboun=iboun+1 
*iboun *entit
*end elems
*else
*set Cond DOMAIN_boundary_3D       *elems
*loop elems *onlyincond *NoCanRepeat
*set var entit=cond(1,int)
*set var iboun=iboun+1
*iboun *entit 
*end elems
*endif
  END_BOUNDARIES

  NODES
  END_NODES
