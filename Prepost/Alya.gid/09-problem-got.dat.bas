*if(strcasecmp(GenData(GOTITA_Module),"on")==0)
*realformat "%16.6f"
*intformat  "%8i"
$------------------------------------------------------------------------
PHYSICAL_PROBLEM
PROBLEM_DEFINITION
  TEMPORAL_DERIVATIVES:    *GenData(GOT_Temporal_derivatives) 
*if(strcasecmp(GenData(GOT_Air_velocity),"VELOC")==0) 
  AIR_VELOCITY:            VELOC
*else 
  AIR_VELOCITY:            FUNCTION=*GenData(GOT_Air_velocity)
*endif 
END_PROBLEM_DEFINITION
$------------------------------------------------------------------------
PROPERTIES:
  DENSITY=                 *GenData(GOT_Density) 
  AIR_DENSITY=             *GenData(GOT_Air_density)      
  AIR_VISCOSITY=           *GenData(GOT_Air_viscosity)    
  VELOCITY=                *GenData(GOT_Velocity)   
  LENGTH=                  *GenData(GOT_Length)   
  DROPLET_DIAMETER=        *GenData(GOT_Droplet_diameter)   
END_PROPERTIES
$------------------------------------------------------------------------
END_PHYSICAL_PROBLEM
$------------------------------------------------------------------------
NUMERICAL_TREATMENT
  STABILITY_CONSTANTS=     *GenData(GOT_Tau_1) *GenData(GOT_Tau_2)
  ELEMENT_LENGTH:          *GenData(GOT_Element_characteristic_length)
  SAFETY_FACTOR=           *GenData(GOT_Safety_factor)  
  STEADY_STATE_TOLERANCE=  *GenData(GOT_Steady_state_tolerance)  
  NORM_OF_CONVERGENCE:     *GenData(GOT_Norm_of_convergence)   
  MAXIMUM_NUMBER_OF_ITER=  *GenData(GOT_Number_of_internal_iterations)  
  RELAXATION_FACTOR=       *GenData(GOT_Relaxation_factor)
  CONVERGENCE_TOLERANCE=   *GenData(GOT_Convergence_tolerance)  
*if(strcasecmp(GenData(GOT_Algebraic_solver),"EXPLI")==0) 
  ALGEBRAIC_SOLVER:        *Gendata(GOT_Algebraic_solver)  
*elseif(strcasecmp(GenData(GOT_Algebraic_solver),"GMRES")==0)  
  ALGEBRAIC_SOLVER:        *Gendata(GOT_Algebraic_solver)  
  SOLVER_ITERATIONS=       *GenData(GOT_Solver_iterations)  
  TOLERANCE_TO_CONVER:     *GenData(GOT_Solver_tolerance)  
  KRYLOV_DIMENSION=        *GenData(GOT_Krylov_dimension)   
*else
  ALGEBRAIC_SOLVER:        *Gendata(GOT_Direct_solver)  
*endif
END_NUMERICAL_TREATMENT
$------------------------------------------------------------------------
OUTPUT_&_POST_PROCESS
  START_POSTPROCES_AT      STEP =*GenData(GOT_Start_postprocess_at_step)
  POSTPROCESS CDROP,       STEPS=*GenData(GOT_Postprocess_water_volume_fraction)
  POSTPROCESS VDROP,       STEPS=*GenData(GOT_Postprocess_droplet_velocity)
  POSTPROCESS MDROP,       STEPS=*GenData(GOT_Postprocess_droplet_momentum)
*if(strcasecmp(GenData(GOT_Exact_solution),"No")!=0)  
  OUTPUT      ERROR,       SOLUTION=*GenData(GOT_Exact_solution)
*endif
END_OUTPUT_&_POST_PROCESS
$------------------------------------------------------------------------
BOUNDARY_CONDITIONS UNKNOWN_BOUNDARIES
  INCLUDE *tcl(AlyaProblemName).got.fix
END_BOUNDARY_CONDITIONS
$------------------------------------------------------------------------
*endif
