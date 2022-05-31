*if(strcasecmp(GenData(WAVEQU_Module),"on")==0)
*realformat "%16.6f"
*intformat  "%8i"
$------------------------------------------------------------------------
PHYSICAL_PROBLEM
PROBLEM_DEFINITION
*if(strcasecmp(GenData(WAV_Source_terms),"Function")==0)
  SOURCE_TERMS:           FUNCTION=*GenData(WAV_Source_function), PARAMETERS=*GenData(WAV_Source_value) 
*else
  SOURCE_TERMS:           *GenData(WAV_Source_terms), VALUE: *GenData(WAV_Source_value)  
*endif
END_PROBLEM_DEFINITION
$------------------------------------------------------------------------
*if(strcasecmp(GenData(WAV_Materials),"One")==0)
PROPERTIES:
  DENSITY=                *GenData(WAV_Density)
  KAPPA=                  *GenData(WAV_Kappa) 
  LAW_DENSITY=            *GenData(WAV_Law_density)
  LAW_SPECIFIC_HEAT=      *GenData(WAV_Law_kappa)  
END_PROPERTIES
*else
PROPERTIES,    MATERIALS: *nmats
*if(ndime==2)
*set elems(Triangle)  
*add elems(Quadrilateral)
*else
*set elems(Tetrahedra)
*add elems(Hexahedra)
*add elems(Prisma)
*endif
  DENSITY=                
*loop materials
  *MatNum  *MatProp(WAV_Density) 
*end    
  END_DENSITY
  KAPPA_COND=               
*loop materials
  *MatNum  *MatProp(WAV_Kappa) 
*end    
  END_KAPPA_COND
  MATERIALS
*loop elems
   *elemsNum *ElemsMat
*end
  END_MATERIALS
END_PROPERTIES
*endif
$------------------------------------------------------------------------
END_PHYSICAL_PROBLEM
$------------------------------------------------------------------------
NUMERICAL_TREATMENT
  SUBGRID_SCALE:           *GenData(WAV_Subgrid_scale)
  TIME_TREATMENT:          *GenData(WAV_Time_treatment)
  TIME_INTEGRATION:        *GenData(WAV_Time_integration_type), \
                           BETA=*GenData(WAV_Newmark_beta), \
                           GAMMA=*GenData(WAV_Newmark_gamma), \ 
                           *GenData(WAV_Newmark_mass_matrix)
  SAFETY_FACTOR=           *GenData(WAV_Safety_factor)  
  STEADY_STATE_TOLERANCE=  *GenData(WAV_Steady_state_tolerance)
  MAXIMUM_NUMBER_OF_ITER=  *GenData(WAV_Number_of_internal_iterations)  
  CONVERGENCE_TOLERANCE=   *GenData(WAV_Convergence_tolerance)    
*if(strcasecmp(GenData(WAV_Algebraic_solver),"GMRES")==0)  
  ALGEBRAIC_SOLVER:        *Gendata(WAV_Algebraic_solver)  
  SOLVER_ITERATIONS=       *GenData(WAV_Solver_iterations)  
  TOLERANCE_TO_CONVER:     *GenData(WAV_Solver_tolerance)  
  KRYLOV_DIMENSION=        *GenData(WAV_Krylov_dimension)   
*else
  ALGEBRAIC_SOLVER:        *Gendata(WAV_Direct_solver)  
*endif
END_NUMERICAL_TREATMENT
$------------------------------------------------------------------------
OUTPUT_&_POST_PROCESS
  START_POSTPROCES_AT         STEP =*GenData(WAV_Start_postprocess_at_step)
  POSTPROCESS WAVE_AMPLITUDE, STEPS=*GenData(WAV_Postprocess_wav_amplitude)
  NODE_SET
*if(strcasecmp(GenData(WAV_Node_set_wave_amplitude),"Yes")==0)  
    WAVE_AMPLITUDE
*endif
  END_NODE_SET
END_OUTPUT_&_POST_PROCESS
$------------------------------------------------------------------------
BOUNDARY_CONDITIONS UNKNOWN_BOUNDARIES
  INCLUDE *tcl(AlyaProblemName).wav.fix
END_BOUNDARY_CONDITIONS
$------------------------------------------------------------------------
*endif
