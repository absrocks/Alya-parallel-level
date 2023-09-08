*if(strcasecmp(GenData(TEMPER_Module),"on")==0)
*realformat "%16.6f"
*intformat  "%8i"
$------------------------------------------------------------------------
PHYSICAL_PROBLEM
PROBLEM_DEFINITION
  TEMPORAL_DERIVATIVES:   *GenData(TEM_Temporal_derivatives) 
*if(strcasecmp(GenData(TEM_Convective_velocity),"VELOC")==0) 
  CONVECTIVE_TERM:        *GenData(TEM_Convective_term) , VELOC: NASTIN
*else
  CONVECTIVE_TERM:        *GenData(TEM_Convective_term) , VELOC: FUNCTION=*GenData(TEM_Convective_velocity)
*endif 
  CONDUCTION_TERM:        *GenData(TEM_Conductivity_term)
  SOURCE_TERMS:           *GenData(TEM_Source_terms), VALUE: *GenData(TEM_Source_value)  
*if(strcasecmp(GenData(TEM_Turbulence_model),"Off")==0) 
  TURBULENCE_MODEL:       Off 
*elseif(strcasecmp(GenData(TEM_Turbulence_model),"From_TURBUL")==0)
  TURBULENCE_MODEL:       RANS_DIFFERENTIAL
*elseif(strcasecmp(GenData(TEM_Turbulence_model),"Smagorinsky")==0)
  TURBULENCE_MODEL:       LES_ALGORITHM,   SMAGORINSKY,   PARAMETER: *GenData(TEM_Para_turbulence) 
*elseif(strcasecmp(GenData(TEM_Turbulence_model),"Mixing_length")==0)
  TURBULENCE_MODEL:       RANS_ALGEBRAIC,  MIXING_LENGTH, PARAMETER: *GenData(TEM_Para_turbulence) 
*elseif(strcasecmp(GenData(TEM_Turbulence_model),"Constant")==0)
  TURBULENCE_MODEL:       RANS_ALGEBRAIC,  C0NSTANT,      PARAMETER: *GenData(TEM_Para_turbulence) 
*endif   
END_PROBLEM_DEFINITION
$------------------------------------------------------------------------
*if(nmats==0)
PROPERTIES:
  DENSITY=                *GenData(TEM_Density)
  SPECIFIC_HEAT=          *GenData(TEM_Specific_heat)  
  THERMAL_COND=           *GenData(TEM_Thermal_conductivity)  
  VISCOSITY=              *GenData(TEM_Viscosity)  
  TURBULENT_PRANDTL=      *GenData(TEM_Turbulent_Prandtl_number)  
  REACTION=               *GenData(TEM_Reaction)
  LAW_DENSITY=            *GenData(TEM_Law_density)
  LAW_SPECIFIC_HEAT=      *GenData(TEM_Law_specific_heat)  
  LAW_THERMAL_COND=       *GenData(TEM_Law_thermal_conductivity)  
  LAW_VISCOSITY=          *GenData(TEM_Law_viscosity)  
END_PROPERTIES
*else
PROPERTIES
  DENSITY               
*loop materials
  *MatNum  *MatProp(TEM_Density) 
*end    
  END_DENSITY
  SPECIFIC_HEAT                
*loop materials
  *MatNum  *MatProp(TEM_Specific_heat) 
*end    
  END_SPECIFIC_HEAT
  THERMAL_COND              
*loop materials
  *MatNum  *MatProp(TEM_Thermal_conductivity) 
*end    
  END_THERMAL_COND
  REACTION=*GenData(TEM_Reaction)
  TURBULENT_PRANDTL=*GenData(TEM_Turbulent_Prandtl_number)  
  LAW_DENSITY                
*loop materials
  *MatNum  *MatProp(TEM_Law_density) 
*end    
  END_LAW_DENSITY
  LAW_SPECIFIC_HEAT                
*loop materials
  *MatNum  *MatProp(TEM_Law_specific_heat) 
*end    
  END_LAW_SPECIFIC_HEAT
  LAW_THERMAL_COND              
*loop materials
  *MatNum  *MatProp(TEM_Law_thermal_conductivity) 
*end    
  END_LAW_THERMAL_COND
END_PROPERTIES
*endif
$------------------------------------------------------------------------
END_PHYSICAL_PROBLEM
$------------------------------------------------------------------------
NUMERICAL_TREATMENT
*if(strcasecmp(GenData(TEM_Stability_constants),"Default")==0) 
  STABILITY_CONSTANTS:      1.0,1.0,1.0
*else
  STABILITY_CONSTANTS:      *GenData(TEM_Tau_1) *GenData(TEM_Tau_2) *GenData(TEM_Tau_3)
*endif
  TAU_STRATEGY:            *GenData(TEM_Tau_strategy)
  ELEMENT_LENGTH:          *GenData(TEM_Element_characteristic_length)
  SHOCK_CAPTURING=         *GenData(TEM_Shock_capturing), VALUE: *GenData(TEM_Value_sh)  
*if(strcasecmp(GenData(TEM_Subscale_time_tracking),"On")==0)
  TRACKING:                 TIME
*endif  
  TIME_INTEGRATION:        *GenData(tem_Time_integration_type), ORDER: *GenData(TEM_Time_integration_order), EULER=*GenData(TEM_Euler_iterations)
  SAFETY_FACTOR=           *GenData(TEM_Safety_factor)  
  STEADY_STATE_TOLERANCE=  *GenData(TEM_Steady_state_tolerance)  
  NORM_OF_CONVERGENCE:     *GenData(TEM_Norm_of_convergence)   
  MAXIMUM_NUMBER_OF_ITER=  *GenData(TEM_Number_of_internal_iterations)  
  RELAXATION_FACTOR=       *GenData(TEM_Relaxation_factor)
  CONVERGENCE_TOLERANCE=   *GenData(TEM_Convergence_tolerance)  
*if(strcasecmp(GenData(TEM_Algebraic_solver),"EXPLI")==0) 
  ALGEBRAIC_SOLVER:        *Gendata(TEM_Algebraic_solver)  
*elseif(strcasecmp(GenData(TEM_Algebraic_solver),"GMRES")==0)  
  ALGEBRAIC_SOLVER:        *Gendata(TEM_Algebraic_solver)  
  SOLVER_ITERATIONS=       *GenData(TEM_Solver_iterations)  
  TOLERANCE_TO_CONVER:     *GenData(TEM_Solver_tolerance)  
  KRYLOV_DIMENSION=        *GenData(TEM_Krylov_dimension)   
*else
  ALGEBRAIC_SOLVER:        *Gendata(TEM_Direct_solver)  
*endif
END_NUMERICAL_TREATMENT
$------------------------------------------------------------------------
OUTPUT_&_POST_PROCESS
  POSTPROCESS INITIAL
  START_POSTPROCES_AT      STEP =*GenData(TEM_Start_postprocess_at_step)
  POSTPROCESS TEMPER,      STEPS=*GenData(TEM_Postprocess_temperature)
  POSTPROCESS HEAT_FLUX,   STEPS=*GenData(TEM_Postprocess_heat_flux)
  POSTPROCESS RESIDUAL,    STEPS=*GenData(TEM_Postprocess_residual)
  POSTPROCESS VELOCITY,    STEPS=*GenData(TEM_Postprocess_velocity)
*if(strcasecmp(GenData(TEM_Exact_solution),"No")!=0)  
  OUTPUT      ERROR,       SOLUTION=*GenData(TEM_Exact_solution)
*endif
  BOUNDARY_SET
*if(strcasecmp(GenData(TEM_Boundary_set_mean_temperature),"Yes")==0)  
    MEAN_TEMPERATURE
*endif
*if(strcasecmp(GenData(TEM_Boundary_set_heat_flux),"Yes")==0)  
    HEAT_FLUX
*endif
  END_BOUNDARY_SET
  ELEMENT_SET
*if(strcasecmp(GenData(TEM_Element_set_mean_convective_temperature),"Yes")==0) 
*if(strcasecmp(GenData(TEM_Element_set_axes_mean_convective_temperature),"x")==0) 
    MEAN_CONVECTIVE_TEMPERATURE, AXES=1.0,0.0,0.0
*elseif(strcasecmp(GenData(TEM_Element_set_axes_mean_convective_temperature),"y")==0) 
    MEAN_CONVECTIVE_TEMPERATURE, AXES=0.0,1.0,0.0
*else
    MEAN_CONVECTIVE_TEMPERATURE, AXES=0.0,0.0,1.0
*endif
*endif
  END_ELEMENT_SET
END_OUTPUT_&_POST_PROCESS
$------------------------------------------------------------------------
*if(strcasecmp(GenData(TEM_First_iteration),"Off")==0) 
BOUNDARY_CONDITIONS UNKNOWN_BOUNDARIES, NON_CONSTANT
*else
BOUNDARY_CONDITIONS UNKNOWN_BOUNDARIES, NON_CONSTANT, INITIAL: DIFFUSION 
*endif
*if(strcasecmp(GenData(TEM_Node_codes),"On")==0) 
  CODES, NODES
    *GenData(TEM_Node_value_code_1) 
    *GenData(TEM_Node_value_code_2) 
    *GenData(TEM_Node_value_code_3) 
    *GenData(TEM_Node_value_code_4) 
    *GenData(TEM_Node_value_code_5) 
  END_CODES
*endif
  INCLUDE *tcl(AlyaProblemSolo).tem.fix
END_BOUNDARY_CONDITIONS
$------------------------------------------------------------------------
*endif
