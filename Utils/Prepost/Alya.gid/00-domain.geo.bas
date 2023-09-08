*if(ndime==2)
*set elems(Triangle)  
*add elems(Quadrilateral)
*else
*set elems(Tetrahedra)  
*add elems(Hexahedra)
*add elems(Prism)
*add elems(Pyramid)
*endif
*set var ielem=0
*\------------------------------------------------------------------------
*\
*\ NODES PER ELEMENTS
*\ 
*\------------------------------------------------------------------------
  NODES_PER_ELEMENT
*loop elems
   *ElemsNum *ElemsNnode
*end
  END_NODES_PER_ELEMENT
*\------------------------------------------------------------------------
*\
*\ ELEMENTS
*\ 
*\------------------------------------------------------------------------
  ELEMENTS
*loop elems
  *ElemsNum *elemsConec  
*end
  END_ELEMENTS  
*\------------------------------------------------------------------------
*\
*\ COODRINATES
*\ 
*\------------------------------------------------------------------------
*realformat "%14.6e"
  COORDINATES  
*loop nodes
   *NodesNum *NodesCoord   
*end
  END_COORDINATES
*\------------------------------------------------------------------------
*\
*\ BOUNDARIES
*\ 
*\------------------------------------------------------------------------
  LELBO
*set var iboun=0
*if(ndime==2)
*set Cond DOMAIN_boundary_2D       *elems
*loop elems *onlyincond *canrepeat
*set var iboun=iboun+1
*format "%i %i "
*iboun *ElemsNum 
*end elems
*else
*set Cond DOMAIN_boundary_3D       *elems
*loop elems *onlyincond *canrepeat
*set var iboun=iboun+1
*format "%i %i "
*iboun *ElemsNum 
*end elems
*endif   
  END_LELBO
  BOUNDARIES
*set var iboun=0
*if(ndime==2)
*set Cond DOMAIN_boundary_2D       *elems
*loop elems *onlyincond *canrepeat
*set var entit=cond(1,int)
*set var iboun=iboun+1
*format "%i %i %i %i %i %i %i %i %i %i %i %i "
*iboun *globalnodes 
*end elems
*else
*set Cond DOMAIN_boundary_3D       *elems
*loop elems *onlyincond *canrepeat
*set var entit=cond(1,int)
*set var iboun=iboun+1
*format "%i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i %i "
*iboun *globalnodes
*end elems
*endif
  END_BOUNDARIES
*\------------------------------------------------------------------------
*\
*\ MATERIALS
*\ 
*\------------------------------------------------------------------------
*if(nmats!=0) 
  MATERIALS
*loop elems 
   *elemsNum *ElemsMat
*end
  END_MATERIALS 
*endif 
*\------------------------------------------------------------------------
*\
*\ SUBDOMAINS
*\ 
*\------------------------------------------------------------------------
  SUBDOMAINS
*if(ndime==2)
*set Cond SUBDOMAINS_2D       *elems
*else
*set Cond SUBDOMAINS_3D       *elems
*endif
*loop elems *onlyincond *NoCanRepeat
*set var entit=cond(1,int)
   *ElemsNum *entit
*end
  END_SUBDOMAINS
*\------------------------------------------------------------------------
*\
*\ CHARACTERISTICS
*\ 
*\------------------------------------------------------------------------
  CHARACTERISTICS
*set var kflcontact=3
*set var epsilon=0.000000001
*set var ielem=0
*if(ndime==2)
*#
*#  2D problem, check only for quadris
*#
*set elems(Quadrilateral)
*loop elems
*set var ielem=ielem+1
*set var kelem=elemsNum 
*if((NodesCoord(1,2)>NodesCoord(4,2) - epsilon) && (NodesCoord(1,2)<NodesCoord(4,2) + epsilon))
*if((NodesCoord(1,1)>NodesCoord(4,1) - epsilon) && (NodesCoord(1,1)<NodesCoord(4,1) + epsilon))
*#  contact found
  *kelem *kflcontact 
*endif
*endif
*end
*else
*#
*#  2D problem, check only for hexas and prisms
*#
*set elems(Hexahedra)
*loop elems
*set var ielem=ielem+1
*set var kelem=elemsNum 
*if((NodesCoord(1,3)>NodesCoord(5,3) - epsilon) && (NodesCoord(1,3)<NodesCoord(5,3) + epsilon))
*if((NodesCoord(1,2)>NodesCoord(5,2) - epsilon) && (NodesCoord(1,2)<NodesCoord(5,2) + epsilon))
*if((NodesCoord(1,1)>NodesCoord(5,1) - epsilon) && (NodesCoord(1,1)<NodesCoord(5,1) + epsilon))
*#  contact found
  *kelem *kflcontact 
*endif
*endif
*endif
*end
*set elems(Prism)
*loop elems
*set var ielem=ielem+1
*set var kelem=elemsNum 
*if((NodesCoord(1,3)>NodesCoord(4,3) - epsilon) && (NodesCoord(1,3)<NodesCoord(4,3) + epsilon))
*if((NodesCoord(1,2)>NodesCoord(4,2) - epsilon) && (NodesCoord(1,2)<NodesCoord(4,2) + epsilon))
*if((NodesCoord(1,1)>NodesCoord(4,1) - epsilon) && (NodesCoord(1,1)<NodesCoord(4,1) + epsilon))
*#  contact found
  *kelem *kflcontact 
*endif
*endif
*endif
*end
*endif
  END_CHARACTERISTICS
