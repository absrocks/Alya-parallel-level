*realformat "%16.6f"
*intformat "%6i"
*set var nboun=0
*if(ndime==2)
*set Cond DOMAIN_boundary_2D       *elems
*loop elems *onlyincond *canrepeat
*set var entit=cond(1,int)
*set var nboun=nboun+1
*end elems
*else
*set Cond DOMAIN_boundary_3D       *elems
*loop elems *onlyincond *canrepeat
*set var entit=cond(1,int)
*set var nboun=nboun+1
*end elems
*endif
*set var bar02=0
*set var bar03=0
*set var tri03=0
*set var tri06=0
*set var qua04=0
*set var qua08=0
*set var qua09=0
*set var tet04=0
*set var tet10=0
*set var hex08=0
*set var hex20=0
*set var hex27=0
*set var pri06=0
*if(ndime==1)
*set elems(Linear)  
*loop elems
*if(ElemsNnode==2)
*set var bar02=bar02+1
*else
*set var bar03=bar03+1
*end
*endif
*elseif(ndime==2)
*set elems(Triangle)  
*loop elems
*if(ElemsNnode==3)
*set var tri03=tri03+1
*else
*set var tri06=tri06+1
*endif
*end
*set elems(Quadrilateral)  
*loop elems
*if(ElemsNnode==4)
*set var qua04=qua04+1
*elseif(ElemsNnode==8)
*set var qua08=qua08+1
*else
*set var qua09=qua09+1
*endif
*end
*else
*set elems(Tetrahedra)
*loop elems
*if(ElemsNnode==4)
*set var tet04=tet04+1
*else
*set var tet10=tet10+1
*endif
*end
*set elems(Hexahedra)  
*loop elems
*if(ElemsNnode==27)
*set var hex27=hex27+1
*elseif(ElemsNnode==20)
*set var hex20=hex20+1
*else
*set var hex08=hex08+1
*endif
*end
*set elems(Prisma)  
*loop elems
*set var pri06=pri06+1
*end
*endif
*set var nskew=0
*set var nslav=0
*set var ngive=0
$------------------------------------------------------------
*if(ndime==2)
*set elems(Triangle)  
*add elems(Quadrilateral)
*else
*set elems(Tetrahedra)
*add elems(Hexahedra)
*add elems(Prisma)
*endif 
DIMENSIONS
  NODAL_POINTS=       *npoin
*set var kelem=nelem
*if(tri03>0)
*set var tri03=10
*endif
*if(tri06>0)
*set var tri06=11
*endif
*if(qua04>0)
*set var qua04=12
*endif
*if(qua08>0)
*set var qua08=13
*endif
*if(qua09>0)
*set var qua09=14
*endif
*if(tet04>0)
*set var tet04=30
*endif
*if(tet10>0)
*set var tet10=31
*endif
*if(pri06>0)
*set var pri06=34
*endif
*if(hex08>0)
*set var hex08=37
*endif
*if(hex20>0)
*set var hex20=38
*endif
*if(hex27>0)
*set var hex27=39
*endif
*set var nperi=0
*if(ndime==2)
*set Cond DOMAIN_Periodic_Master_2D *nodes
*else
*set Cond DOMAIN_Periodic_Master_3D *nodes
*endif
*loop nodes *OnlyInCond
*format " %i "
*set var nperi=nperi+1
*end
  ELEMENTS=           *kelem
  SPACE_DIMENSIONS=   *ndime
*if(ndime==2)
  TYPES_OF_ELEMENTS=  *tri03, *tri06, *qua04, *qua08, *qua09
*else  
  TYPES_OF_ELEMENTS=  *tet04, *tet10, *pri06, *hex08, *hex20, *hex27
*endif
  BOUNDARIES=         *nboun
*set var field=0
*set var field1=0
*set var field2=0
*set var field3=0
*if(ndime==2)
*Set Cond FIELD_ELEMENT_1_2D *elems
*else
*Set Cond FIELD_ELEMENT_1_3D *elems
*endif
*loop elems *onlyincond
*set var field1=1
*end
*if(ndime==2)
*Set Cond FIELD_ELEMENT_2_2D *elems
*else
*Set Cond FIELD_ELEMENT_2_3D *elems
*endif
*loop elems *onlyincond
*set var field2=1
*end
*if(ndime==2)
*Set Cond FIELD_ELEMENT_2_2D *elems
*else
*Set Cond FIELD_ELEMENT_2_3D *elems
*endif
*loop elems *onlyincond
*set var field3=1
*end
*set var field=field1+field2+field3
*if(field>0)
  FIELDS=         *field
*if(field1>0)
*if(ndime==2)
*Set Cond FIELD_ELEMENT_1_2D *elems
*else
*Set Cond FIELD_ELEMENT_1_3D *elems
*endif
    FIELD=1, DIMENSION=1, Elements
*endif
*if(field2>0)
*if(ndime==2)
*Set Cond FIELD_ELEMENT_2_2D *elems
*else
*Set Cond FIELD_ELEMENT_2_3D *elems
*endif
    FIELD=2, DIMENSION=1, Elements
*endif
*if(field3>0)
*if(ndime==2)
*Set Cond FIELD_ELEMENT_3_2D *elems
*else
*Set Cond FIELD_ELEMENT_3_3D *elems
*endif
    FIELD=3, DIMENSION=1, Elements
*endif
  END_FIELDS
*endif
END_DIMENSIONS
$------------------------------------------------------------
STRATEGY
  INTEGRATION_RULE:                *GenData(Integration_rule)  
  DOMAIN_INTEGRATION_POINTS:       *GenData(Integration_points)  
*if(strcasecmp(GenData(Scaling_of_Geometry),"On")==0)
  SCALE_FACTORS:                   X_SCALE: *GenData(X_Scale), Y_SCALE: *GenData(Y_Scale), Z_SCALE: *GenData(Z_Scale)
*endif
*if(strcasecmp(GenData(Extrapolate_from_boundaries_to_nodes),"On")==0)
  EXTRAPOLATE_BOUNDARY_CONDITIONS: On
*endif
  BOUNDARY_ELEMENT:                On
  GROUPS=                          *GenData(Number_groups,int), *GenData(Groups_strategy)  
END_STRATEGY
$-------------------------------------------------------------
GEOMETRY
  INCLUDE  *tcl(AlyaProblemSolo).geo.dat
END_GEOMETRY  
$-------------------------------------------------------------
SETS
  INCLUDE  *tcl(AlyaProblemSolo).set.dat  
END_SETS
$-------------------------------------------------------------
BOUNDARY_CONDITIONS
  INCLUDE  *tcl(AlyaProblemSolo).fix.dat  
END_BOUNDARY_CONDITIONS
$-------------------------------------------------------------
FIELDS
*if(field>0)
  INCLUDE *tcl(AlyaProblemSolo).fie.dat  
*endif
END_FIELDS
$-------------------------------------------------------------
