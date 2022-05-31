*if(strcasecmp(GenData(LEVELS_Module),"on")==0)
*realformat "%16.6f"
*intformat  "%8i"
$------------------------------------------------------------------------
PHYSICAL_PROBLEM
PROBLEM_DEFINITION
*if(strcasecmp(GenData(LEV_Source_terms),"Function")==0)
  SOURCE_TERMS:           FUNCTION=*GenData(LEV_Source_function), PARAMETERS=*GenData(LEV_Source_value) 
*else
  SOURCE_TERMS:           *GenData(LEV_Source_terms), VALUE: *GenData(LEV_Source_value)  
*endif
END_PROBLEM_DEFINITION
$------------------------------------------------------------------------
*if(strcasecmp(GenData(LEV_Materials),"One")==0)
PROPERTIES:
  DENSITY=                *GenData(LEV_Density)
  KAPPA=                  *GenData(LEV_Kappa) 
  LAW_DENSITY=            *GenData(LEV_Law_density)
  LAW_SPECIFIC_HEAT=      *GenData(LEV_Law_kappa)  
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
  *MatNum  *MatProp(LEV_Density) 
*end    
  END_DENSITY
  KAPPA_COND=               
*loop materials
  *MatNum  *MatProp(LEV_Kappa) 
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
  SUBGRID_SCALE:           *GenData(LEV_Subgrid_scale)
  TIME_TREATMENT:          *GenData(LEV_Time_treatment)
  TIME_INTEGRATION:        *GenData(LEV_Time_integration_type), \
                           BETA=*GenData(LEV_Newmark_beta), \
                           GAMMA=*GenData(LEV_Newmark_gamma), \ 
                           *GenData(LEV_Newmark_mass_matrix)
  SAFETY_FACTOR=           *GenData(LEV_Safety_factor)  
  STEADY_STATE_TOLERANCE=  *GenData(LEV_Steady_state_tolerance)
  MAXIMUM_NUMBER_OF_ITER=  *GenData(LEV_Number_of_internal_iterations)  
  CONVERGENCE_TOLERANCE=   *GenData(LEV_Convergence_tolerance)    
*if(strcasecmp(GenData(LEV_Algebraic_solver),"GMRES")==0)  
  ALGEBRAIC_SOLVER:        *Gendata(LEV_Algebraic_solver)  
  SOLVER_ITERATIONS=       *GenData(LEV_Solver_iterations)  
  TOLERANCE_TO_CONVER:     *GenData(LEV_Solver_tolerance)  
  KRYLOV_DIMENSION=        *GenData(LEV_Krylov_dimension)   
*else
  ALGEBRAIC_SOLVER:        *Gendata(LEV_Direct_solver)  
*endif
END_NUMERICAL_TREATMENT
$------------------------------------------------------------------------
OUTPUT_&_POST_PROCESS
  START_POSTPROCES_AT         STEP =*GenData(LEV_Start_postprocess_at_step)
  POSTPROCESS LEVEL_SET, STEPS=*GenData(LEV_Postprocess_level_set)
  NODE_SET
*if(strcasecmp(GenData(LEV_Node_set_level_set),"Yes")==0)  
    WAVE_AMPLITUDE
*endif
  END_NODE_SET
END_OUTPUT_&_POST_PROCESS
$------------------------------------------------------------------------
BOUNDARY_CONDITIONS UNKNOWN_BOUNDARIES
  INCLUDE *tcl(AlyaProblemName).lev.fix
END_BOUNDARY_CONDITIONS
$------------------------------------------------------------------------
*endif
