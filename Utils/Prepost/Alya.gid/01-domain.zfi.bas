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
*if(field1>0)
  FIELD: 1 
*loop elems *onlyincond *canrepeat
    *elemsNum *cond(1,int)
*end elems
  END_FIELD
*endif
*if(ndime==2)
*Set Cond FIELD_ELEMENT_2_2D *elems
*else
*Set Cond FIELD_ELEMENT_2_3D *elems 
*endif
*loop elems *onlyincond
*set var field2=1
*end
*if(field2>0)
  FIELD: 2 
*loop elems *onlyincond *canrepeat
    *elemsNum *cond(1,int)
*end elems
  END_FIELD
*endif
*if(ndime==2)
*Set Cond FIELD_ELEMENT_3_2D *elems
*else
*Set Cond FIELD_ELEMENT_3_3D *elems
*endif
*loop elems *onlyincond
*set var field3=1
*end
*if(field3>0)
  FIELD: 3 
*loop elems *onlyincond *canrepeat
    *elemsNum *cond(1,int)
*end elems
  END_FIELD
*endif
