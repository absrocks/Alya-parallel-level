$-------------------------------------------------------------------
PHYSICAL_PROBLEM
  PROPERTIES
*if(nmats!=0)     
*loop materials
    MATERIAL: *MatNum,     NAME=*MatProp(0) 
      DENSITY=       CONSTANT, VALUE=*MatProp(Density)
      VISCOSITY=     CONSTANT, VALUE=*MatProp(Viscosity)
      SPECIFIC_HEAT= CONSTANT, VALUE=*MatProp(Specific_heat)
      CONDUCTIVITY=  CONSTANT, VALUE=*MatProp(Conductivity)
    END_MATERIAL
*end materials
*else
    MATERIAL: 1,     NAME=GENERIC
      DENSITY=       CONSTANT, VALUE=1.0
      VISCOSITY=     CONSTANT, VALUE=1.0
      SPECIFIC_HEAT= CONSTANT, VALUE=1.0
      CONDUCTIVITY=  CONSTANT, VALUE=1.0
    END_MATERIAL
*endif
  END_PROPERTIES
END_PHYSICAL_PROBLEM
$-------------------------------------------------------------------
NUMERICAL_TREATMENT 
  ELSEST
    STRATEGY:    *GenData(KER_ELSEST_Strategy)
*if(strcasecmp(GenData(KER_ELSEST_Strategy),"Bin")==0)
    NUMBER_BINS= *GenData(KER_Number_bins)
*else
    MAXIMUM=     *GenData(KER_Max_nodes)
*endif
  END_ELSEST
  MESH
    MULTIPLICATION= *GenData(KER_Multiplication)
  END_MESH
END_NUMERICAL_TREATMENT  
$-------------------------------------------------------------------
OUTPUT_&_POST_PROCESS  
  *GenData(KER_Postprocess)
*if( GenData(KER_Post_frequency,int) != -2 )
   STEPS = *GenData(KER_Post_frequency,int)
*endif
END_OUTPUT_&_POST_PROCESS  
$-------------------------------------------------------------------
