*if(strcasecmp(GenData(TURBUL_Module),"on")==0)
*realformat "%16.6f"
*intformat "%8i" 
$------------------------------------------------------------
PHYSICAL_PROBLEM
$------------------------------------------------------------
PROBLEM_DEFINITION
  MODEL:                     *GenData(TUR_model)
  TEMPORAL_DERIVATIVES:      *GenData(TUR_temporal_derivatives) 
*if(strcasecmp(GenData(TUR_Boussinesq_coupling),"Off")==0) 
  TEMPERATURE_COUPLING:      Off
*else    
  TEMPERATURE_COUPLING:      BOUSSINESQ, BETA: *GenData(TUR_Beta), T_R: *GenData(TUR_T_reference)
  GRAVITY:                   NORM: *GenData(TUR_Norm_gravity), \
                             GX:   *GenData(TUR_X_gravity), \
                             GY:   *GenData(TUR_Y_gravity), \
                             GZ:   *GenData(TUR_Z_gravity)  
*endif
END_PROBLEM_DEFINITION
$------------------------------------------------------------
PROPERTIES:
  DENSITY=                   *GenData(TUR_Density)
  VISCOSITY=                 *GenData(TUR_Viscosity)  
  LAW_DENSITY=               *GenData(TUR_Law_density)
  LAW_VISCOSITY=             *GenData(TUR_Law_viscosity)
  TURBULENT_PRANDTL=         *GenData(TUR_Turbulent_Prandtl_number)  
*if(strcasecmp(GenData(TUR_Parameters),"Default")==0)
  REAL_PARAMETERS=           AUTOMATIC
  INTEGER_PARAMETERS=        AUTOMATIC
*else
  REAL_PARAMETERS=           *GenData(TUR_Real_parameters)
  INTEGER_PARAMETERS=        *GenData(TUR_Integer_parameters)
*endif
END_PROPERTIES
$------------------------------------------------------------
END_PHYSICAL_PROBLEM
$------------------------------------------------------------
NUMERICAL_TREATMENT
*if(strcasecmp(GenData(TUR_Stability_constants),"Default")==0)
  STABILITY_CONSTANTS=       1.0,1.0,1.0
*else
  STABILITY_CONSTANTS=       *GenData(TUR_Tau_1) *GenData(TUR_Tau_2) *GenData(TUR_Tau_3)
*endif
  ELEMENT_LENGTH:            *GenData(TUR_Element_characteristic_length)
  SHOCK_CAPTURING=           *GenData(TUR_Shock_capturing), VALUE: *GenData(TUR_Value_sh)  
  TEMPORAL_TERM_WEIGHT=      All
  TIME_INTEGRATION:          *GenData(TUR_Time_integration_type), ORDER: *GenData(TUR_Time_integration_order), EULER=*GenData(TUR_Euler_iterations)
  SAFETY_FACTOR=             *GenData(TUR_Safety_factor)  
  STEADY_STATE_TOLER=        *GenData(TUR_Steady_state_tolerance)  
  NORM_OF_CONVERGENCE:       *GenData(TUR_Norm_of_convergence)   
  MAXIMUM_NUMBER_OF_ITER=    *GenData(TUR_Number_of_internal_iterations)  
  RELAXATION_FACTOR=         *GenData(TUR_Relaxation_factor)
  CONVERGENCE_TOLERANCE=     *GenData(TUR_Convergence_tolerance)  
*if(strcasecmp(GenData(TUR_Algebraic_solver),"Explicit")==0)  
  ALGEBRAIC_SOLVER:          Explicit
*elseif(strcasecmp(GenData(TUR_Algebraic_solver),"GMRES")==0)  
  ALGEBRAIC_SOLVER:          *Gendata(TUR_Algebraic_solver)  
  SOLVER_ITERATIONS=         *GenData(TUR_Solver_iterations)  
  TOLERANCE_TO_CONVER:       *GenData(TUR_Solver_tolerance)  
  KRYLOV_DIMENSION=          *GenData(TUR_Krylov_dimension)   
*else
*if(strcasecmp(GenData(TUR_Direct_solver),"LDU")==0)  
  ALGEBRAIC_SOLVER:          Direct
*else
  ALGEBRAIC_SOLVER:          MUMPS
*endif
*endif
END_NUMERICAL_TREATMENT
$------------------------------------------------------------
OUTPUT_&_POST_PROCESS
  START_POSTPROCES_AT        STEP =*GenData(TUR_Start_postprocess_at_step)
  POSTPROCESS KEY,           STEPS=*GenData(TUR_Postprocess_key) 
  POSTPROCESS EPSILON,       STEPS=*GenData(TUR_Postprocess_epsilon) 
  POSTPROCESS OMEGA,         STEPS=*GenData(TUR_Postprocess_omega) 
  POSTPROCESS NU_TILDE,      STEPS=*GenData(TUR_Postprocess_nu_tilde) 
  POSTPROCESS LENGTH_SCALE,  STEPS=*GenData(TUR_Postprocess_length_scale)
  POSTPROCESS TUR_VISCOSITY, STEPS=*GenData(TUR_Postprocess_turbulent_viscosity) 
  POSTPROCESS WALL_DISTANCE, STEPS=*GenData(TUR_Postprocess_wall_distance) 
  POSTPROCESS YPLUS,         STEPS=*GenData(TUR_Postprocess_yplus) 
  POSTPROCESS YSTAR,         STEPS=*GenData(TUR_Postprocess_ystar)
  POSTPROCESS YVSTAR,        STEPS=*GenData(TUR_Postprocess_yvstar)
END_OUTPUT_&_POST_PROCESS
$------------------------------------------------------------
BOUNDARY_CONDITIONS
 PARAMETERS
   INLET_TYPE:           *GenData(TUR_Inlet_type)
   INTENSITY_TURBULENCE= *GenData(TUR_Intensity_of_turbulence)
   HYDRAULIC_DIAMETER=   *GenData(TUR_Hydraulic_diameter)
   RELAXATION=           *GenData(TUR_Relaxation_boundary_conditions)
   OMEGA_CONDITION=      *GenData(TUR_Omega_boundary_condition)
 END_PARAMETERS
 INCLUDE *tcl(AlyaProblemName).tur.fix
END_BOUNDARY_CONDITIONS
$------------------------------------------------------------
*endif
