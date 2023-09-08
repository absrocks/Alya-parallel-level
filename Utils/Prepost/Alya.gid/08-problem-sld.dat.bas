*if(strcasecmp(GenData(SOLIDZ_Module),"On")==0)
*realformat "%16.6f"
*intformat "%8i"
$------------------------------------------------------------
PHYSICAL_PROBLEM
$------------------------------------------------------------ 
PROBLEM_DEFINITION       
 TEMPORAL_DERIVATIVES:      *GenData(SLD_Temporal_derivatives)  
 MODEL:                     *GenData(SLD_Model)  
END_PROBLEM_DEFINITION  
$------------------------------------------------------------
PROPERTIES
  DENSITY=                  *GenData(SLD_Density)   
  LAW_DENSITY:              *GenData(SLD_Law_density)  
  CONSTITUTIVE_LAW:         *GenData(SLD_Constitutive_law), PARAMETERS=*GenData(SLD_Parameters_Constitutive_law)
END_PROPERTIES  
$------------------------------------------------------------
END_PHYSICAL_PROBLEM  
$------------------------------------------------------------
NUMERICAL_TREATMENT 
  TIME_TREATMENT:           *GenData(SLD_Time_treatment)
  TIME_INTEGRATION:         *GenData(SLD_Integration_scheme)
  DAMPING:                  *GenData(SLD_Numerical_Damping), PARAMETERS=*GenData(SLD_Value_damping)
END_NUMERICAL_TREATMENT  
$------------------------------------------------------------
OUTPUT_&_POST_PROCESS  
  START_POSTPROCES_AT       STEP =*GenData(SLD_Start_postprocess_at_step)
*if(strcasecmp(GenData(SLD_Postprocess_displacement),"0")!=0)
  POSTPROCESS DISPLACEMENT, STEPS=*GenData(SLD_Postprocess_displacement)
*endif
*if(strcasecmp(GenData(SLD_Postprocess_dis_velocity),"0")!=0)
  POSTPROCESS DIS_VELOC,    STEPS=*GenData(SLD_Postprocess_dis_velocity)
*endif
*if(strcasecmp(GenData(SLD_Postprocess_dis_acceleration),"0")!=0)
  POSTPROCESS DIS_ACCEL,    STEPS=*GenData(SLD_Postprocess_dis_acceleration)
*endif
*if(strcasecmp(GenData(SLD_Postprocess_initial_condition),"No")!=0)
  POSTPROCESS INITIAL
*endif
END_OUTPUT_&_POST_PROCESS  
$------------------------------------------------------------
BOUNDARY_CONDITIONS, \
*if(strcasecmp(GenData(NSI_fix_pressure_mode),"Yes")==0)
      FIX_PRESSURE=*GenData(NSI_Fix_pressure_on_node), \
*elseif(strcasecmp(GenData(NSI_fix_pressure_mode),"No")==0)
      DONT_FIX_PRESSURE, \
*endif
*if(strcasecmp(GenData(NSI_Initial_stokes),"On")==0)
      STOKES \
*endif
      RELAXATION=*GenData(NSI_Relaxation_boundary_condition) \
      INITIAL=*GenData(NSI_Initial_step_boundary_condition) \
      UNKNOWN_BOUNDARIES
  INCLUDE *tcl(AlyaProblemName).nsi.fix
END_BOUNDARY_CONDITIONS  
$------------------------------------------------------------
*endif
