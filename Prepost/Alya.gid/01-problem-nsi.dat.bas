*if(strcasecmp(GenData(NASTIN_Module),"Off")!=0)
*realformat "%16.6f"
*intformat "%8i"
$------------------------------------------------------------
PHYSICAL_PROBLEM
$------------------------------------------------------------ 
PROBLEM_DEFINITION       
  TEMPORAL_DERIVATIVES:     *GenData(NSI_temporal_derivatives)  
  CONVECTIVE_TERM:          *GenData(NSI_convective_term)  
  VISCOUS_TERM:             *GenData(NSI_viscous_term)  
*if(strcasecmp(GenData(NSI_Gravity_acceleration),"Off")!=0) 
  GRAVITY:                  NORM: *GenData(NSI_Norm_gravity), \
                            GX:   *GenData(NSI_X_gravity), \
                            GY:   *GenData(NSI_Y_gravity), \
                            GZ:   *GenData(NSI_Z_gravity)  
*endif
*if(strcasecmp(GenData(NSI_Axes_rotation),"Off")!=0)
  AXES_ROTATION:            NORM: *GenData(NSI_Omega_norm), \
                            OX:   *GenData(NSI_X_omega), \
                            OY:   *GenData(NSI_Y_omega), \
                            OZ:   *GenData(NSI_Z_omega)  
  CENTER_ROTATION:           *GenData(NSI_Rotation_center)  
*endif
*if(strcasecmp(GenData(NSI_Axes_velocity),"Off")!=0)
  AXES_VELOCITY:            NORM: *GenData(NSI_Velocity_norm), \
                            VX:   *GenData(NSI_X_velocity), \
                            VY:   *GenData(NSI_Y_velocity), \
                            VZ:   *GenData(NSI_Z_velocity) 
*endif
*if(strcasecmp(GenData(NSI_Turbulence_model),"Off")==0) 
  TURBULENCE_MODEL:         Off 
*elseif(strcasecmp(GenData(NSI_Turbulence_model),"From_TURBUL")==0)
*if(strcasecmp(GenData(NSI_Tur_Pressure),"On")==0) 
  TURBULENCE_MODEL:         RANS_DIFFERENTIAL, TUR_PRESSURE
*else
  TURBULENCE_MODEL:         RANS_DIFFERENTIAL
*endif
*elseif(strcasecmp(GenData(NSI_Turbulence_model),"Smagorinsky")==0)
  TURBULENCE_MODEL:         LES_ALGORITHM,   SMAGORINSKY,   PARAMETER: *GenData(NSI_Para_turbulence) 
*elseif(strcasecmp(GenData(NSI_Turbulence_model),"Mixing_length")==0)
  TURBULENCE_MODEL:         RANS_ALGEBRAIC,  MIXING_LENGTH, PARAMETER: *GenData(NSI_Para_turbulence) 
*elseif(strcasecmp(GenData(NSI_Turbulence_model),"Constant")==0)
  TURBULENCE_MODEL:         RANS_ALGEBRAIC,  C0NSTANT,      PARAMETER: *GenData(NSI_Para_turbulence) 
*endif   
*if(strcasecmp(GenData(NSI_Boussinesq_coupling),"Off")==0) 
  BOUSSINESQ_COUPLING:      Off    
*else 
  BOUSSINESQ_COUPLING:      On, BETA: *GenData(NSI_Beta), T_R: *GenData(NSI_T_reference), G: *GenData(NSI_Modulus_of_G)   
*endif
END_PROBLEM_DEFINITION  
$------------------------------------------------------------
END_PHYSICAL_PROBLEM  
$------------------------------------------------------------
NUMERICAL_TREATMENT 
  TAU_STRATEGY:             *GenData(NSI_Tau_strategy)
*if(strcasecmp(GenData(NSI_Stability_constants),"Default")==0) 
  STABILITY_CONSTANTS:      1.0,1.0,1.0,1.0
*else
  STABILITY_CONSTANTS:      *GenData(NSI_Tau_1) *GenData(NSI_Tau_2) *GenData(NSI_Tau_3) *GenData(NSI_Tau_4)
*endif
  ELEMENT_LENGTH:           *GenData(NSI_Element_characteristic_length)
  TRACKING:                 \
*if(strcasecmp(GenData(NSI_Subscale_convection_tracking),"On")==0) 
                            CONVECTION, ITERA=*GenData(NSI_Subscale_iterations), TOLER=*GenData(NSI_Subscale_tolerance), RELAX=*GenData(NSI_Subscale_relaxation) \
*endif
*if(strcasecmp(GenData(NSI_Subscale_Time_tracking),"On")==0) 
                            TIME \
*endif 
                            CONTINUE                  
  SHOCK_CAPTURING:          *GenData(NSI_Shock_capturing), VALUE: *GenData(NSI_Value_sh) 
*if(strcasecmp(GenData(NSI_Penalty_strategy),"No")==0) 
  PENALTY_METHOD:           Off  
*else
*if(strcasecmp(GenData(NSI_Penalty_strategy),"Automatic")==0) 
  PENALTY_METHOD:           Iterative, VALUE: *GenData(NSI_Epsilon_value), AUTOMATIC
*else
  PENALTY_METHOD:           *GenData(NSI_Penalty_strategy), VALUE: *GenData(NSI_Epsilon_value)  
*endif
*endif
  TIME_INTEGRATION:         *GenData(NSI_Time_integration_type), ORDER: *GenData(NSI_Time_integration_order), EULER=*GenData(NSI_Euler_iterations)
  SAFETY_FACTOR:            *GenData(NSI_Safety_factor)  
  STEADY_STATE_TOLER:       *GenData(NSI_Steady_state_tolerance)   
  NORM_OF_CONVERGENCE:      *GenData(NSI_Norm_of_convergence)   
*if(strcasecmp(GenData(NSI_Linearization_method),"Picard")==0) 
  LINEARIZATION_METHOD:     Picard   
*else
  LINEARIZATION_METHOD:     *GenData(NSI_Linearization_method), PICARD_ITERATIONS: *GenData(NSI_Picard_Iterations)  
*endif
  MAXIMUM_NUMBER_OF_IT:     *GenData(NSI_Number_of_internal_iterations)
  CONVERGENCE_TOLERANCE:    *GenData(NSI_Convergence_tolerance)
  ASSEMBLY:                 *GenData(NSI_Assembly)

*if(strcasecmp(GenData(NSI_Scheme),"Orthomin_momentum")==0)
  ALGORITHM:                SCHUR
   SOLVER:                  ORTHOMIN, MOMENTUM_PRESERVING
    PRECONDITIONER:         TAU
    ELEMENT_LENGTH:         Minimum
    TAU_STRATEGY:           codina
    CORRECTION:             Open
  END_ALGORITHM
*elseif(strcasecmp(GenData(NSI_Scheme),"Orthomin_continuity")==0)
  ALGORITHM:                SCHUR
   SOLVER:                  ORTHOMIN, CONTINUITY_PRESERVING
    PRECONDITIONER:         TAU
    ELEMENT_LENGTH:         Minimum
    TAU_STRATEGY:           codina
    CORRECTION:             Open
  END_ALGORITHM
*endif

*if(strcasecmp(GenData(NSI_Scheme),"Monolithic")!=0)
  MOMENTUM
*else
  ALL_EQUATIONS
*endif
    ALGEBRAIC_SOLVER:       *GenData(NSI_Algebraic_solver), ITERA= *GenData(NSI_Solver_iterations), \
                            TOLER= *GenData(NSI_Solver_tolerance), KRYLOV= *GenData(NSI_Krylov_dimension)     
    PRECONDITIONING:        DIAGONAL
*if(strcasecmp(GenData(NSI_Scheme),"Monolithic")!=0)
  END_MOMENTUM
*else
  ALL_EQUATIONS
*endif

*if(strcasecmp(GenData(NSI_Scheme),"Monolithic")!=0)
  CONTINUITY
    ALGEBRAIC_SOLVER:       *GenData(NSI_Continuity_algebraic_solver), ITERA= *GenData(NSI_Continuity_solver_iterations), \  
                            TOLER= *GenData(NSI_Continuity_solver_tolerance), \
                            KRYLOV= *GenData(NSI_Continuity_krylov_dimension)
    PRECONDITIONING:        *GenData(NSI_Continuity_Preconditioner), TOLERANCE=*GenData(NSI_Continuity_Linelet_tolerance)
  END_CONTINUITY
*endif

*if(strcasecmp(GenData(NSI_Scheme),"Block_Gauss_Seidel")==0)
*if(strcasecmp(GenData(NSI_Continuity_Pressure_Laplacian),"Classical")==0)
  PRESSURE_LAPLACIAN:       CLASSICAL
*elseif(strcasecmp(GenData(NSI_Continuity_Pressure_Laplacian),"Iterative")==0)
  PRESSURE_LAPLACIAN:       ITERATIVE, GAMMA=1.0
*elseif(strcasecmp(GenData(NSI_Continuity_Pressure_Laplacian),"Staggered")==0)
  PRESSURE_LAPLACIAN:       ITERATIVE, GAMMA=0.0
*endif
*endif

END_NUMERICAL_TREATMENT  
$------------------------------------------------------------
OUTPUT_&_POST_PROCESS  
  START_POSTPROCES_AT       STEP =*GenData(NSI_Start_postprocess_at_step)
*if(strcasecmp(GenData(NSI_Postprocess_velocity),"0")!=0)
  POSTPROCESS VELOCITY,     STEPS=*GenData(NSI_Postprocess_velocity)
*endif
*if(strcasecmp(GenData(NSI_Postprocess_pressure),"0")!=0)
  POSTPROCESS PRESSURE,     STEPS=*GenData(NSI_Postprocess_pressure)
*endif
*if(strcasecmp(GenData(NSI_Postprocess_streamlines),"0")!=0)
  POSTPROCESS STREAMLINES,  STEPS=*GenData(NSI_Postprocess_streamlines)
*endif
*if(strcasecmp(GenData(NSI_Postprocess_viscosity),"0")!=0)
  POSTPROCESS VISCOSITY,    STEPS=*GenData(NSI_Postprocess_viscosity)    
*endif
*if(strcasecmp(GenData(NSI_Postprocess_permeability),"0")!=0)
  POSTPROCESS PERMEABILITY, STEPS=*GenData(NSI_Postprocess_permeability)    
*endif
*if(strcasecmp(GenData(NSI_Postprocess_stress),"0")!=0)
  POSTPROCESS STRESS,       STEPS=*GenData(NSI_Postprocess_stress)    
*endif
*if(strcasecmp(GenData(NSI_Postprocess_tangential_traction),"0")!=0)
  POSTPROCESS TANGENTIAL,   STEPS=*GenData(NSI_Postprocess_tangential_traction) 
*endif
*if(strcasecmp(GenData(NSI_Boundary_conditions),"Yes")==0)  
  POSTPROCESS BOUNDARY_COND
*endif
  ELEMENT_SET
*if(strcasecmp(GenData(NSI_Element_set_mean_velocity_module),"Yes")==0)  
    VELOCITY_MEAN_MODULE
*endif
*if(strcasecmp(GenData(NSI_Element_set_mean_vorticity_module),"Yes")==0)  
    VORTICITY_MEAN_MODULE
*endif
*if(strcasecmp(GenData(NSI_Element_set_kinetic_energy),"Yes")==0)  
    KINETIC_ENERGY
*endif
  END_ELEMENT_SET
  BOUNDARY_SET
*if(strcasecmp(GenData(NSI_Boundary_set_mean_pressure),"Yes")==0)  
    MEAN_PRESSURE
*endif
*if(strcasecmp(GenData(NSI_Boundary_set_mass),"Yes")==0)  
    MASS
*endif
*if(strcasecmp(GenData(NSI_Boundary_set_force),"Yes")==0)  
    FORCE
*endif
*if(strcasecmp(GenData(NSI_Boundary_set_torque),"Yes")==0)  
    TORQUE, CENTER=*GenData(NSI_Boundary_set_center_torque)
*endif
*if(strcasecmp(GenData(NSI_Boundary_set_mean_yplus),"Yes")==0)  
    MEAN_YPLUS
*endif
*if(strcasecmp(GenData(NSI_Boundary_set_mean_velocity),"Yes")==0)  
    MEAN_VELOCITY
*endif
  END_BOUNDARY_SET
  NODE_SET
*if(strcasecmp(GenData(NSI_Node_set_pressure),"Yes")==0)  
    PRESSURE
*endif
*if(strcasecmp(GenData(NSI_Node_set_velocity),"Yes")==0)  
    VELOCITY
*endif
*if(strcasecmp(GenData(NSI_Node_set_yplus),"Yes")==0)  
    YPLUS
*endif
  END_NODE_SET
END_OUTPUT_&_POST_PROCESS  
$------------------------------------------------------------
BOUNDARY_CONDITIONS, NON_CONSTANT UNKNOWN_BOUNDARIES
  PARAMETERS
    INITIAL_SOLUTION:               *GenData(NSI_Initial_solution_type), VALUES=*GenData(NSI_Initial_solution_values)
    FIX_PRESSURE=                   *GenData(NSI_Fix_pressure_mode), \
                                    ON_NODE:*GenData(NSI_Fix_pressure_on_node),\
                                    VALUE=*GenData(NSI_Fix_pressure_value)
    RELAXATION_BOUNDARY_CONDITIONS= *GenData(NSI_Relaxation_boundary_condition)
    START_BOUNDARY_CONDITIONS=      *GenData(NSI_Initial_step_boundary_condition)
    WIND_ANGLE=                     *GenData(NSI_Wind_angle)
    WIND_FRICTION_VELOCITY=         *GenData(NSI_Wind_friction_velocity)
    WIND_ZERO_PLANE_DISPLACEMENT=   *GenData(NSI_Wind_zero_plane_displacement)
    WIND_SURFACE_ROUGHNESS=         *GenData(NSI_Wind_surface_roughness)
  END_PARAMETERS
*if(strcasecmp(GenData(NSI_Node_codes),"Off")!=0) 
  CODES, NODES
*if(strcasecmp(GenData(NSI_Node_value_code_01),"Off")!=0) 
$    *GenData(NSI_Node_name_code_01) 
     *GenData(NSI_Node_value_code_01) 
*endif
*if(strcasecmp(GenData(NSI_Node_value_code_02),"Off")!=0) 
$    *GenData(NSI_Node_name_code_02) 
     *GenData(NSI_Node_value_code_02) 
*endif
*if(strcasecmp(GenData(NSI_Node_value_code_03),"Off")!=0) 
$    *GenData(NSI_Node_name_code_03) 
     *GenData(NSI_Node_value_code_03) 
*endif
*if(strcasecmp(GenData(NSI_Node_value_code_04),"Off")!=0) 
$    *GenData(NSI_Node_name_code_04) 
     *GenData(NSI_Node_value_code_04) 
*endif
*if(strcasecmp(GenData(NSI_Node_value_code_05),"Off")!=0) 
$    *GenData(NSI_Node_name_code_05) 
     *GenData(NSI_Node_value_code_05) 
*endif
*if(strcasecmp(GenData(NSI_Node_value_code_06),"Off")!=0) 
$    *GenData(NSI_Node_name_code_06) 
     *GenData(NSI_Node_value_code_06) 
*endif
*if(strcasecmp(GenData(NSI_Node_value_code_07),"Off")!=0) 
$    *GenData(NSI_Node_name_code_07) 
     *GenData(NSI_Node_value_code_07) 
*endif
*if(strcasecmp(GenData(NSI_Node_value_code_08),"Off")!=0) 
$    *GenData(NSI_Node_name_code_08) 
     *GenData(NSI_Node_value_code_08) 
*endif
*if(strcasecmp(GenData(NSI_Node_value_code_09),"Off")!=0) 
$    *GenData(NSI_Node_name_code_09) 
     *GenData(NSI_Node_value_code_09) 
*endif
*if(strcasecmp(GenData(NSI_Node_value_code_10),"Off")!=0) 
$    *GenData(NSI_Node_name_code_10) 
     *GenData(NSI_Node_value_code_10) 
*endif
  END_CODES
*endif
*if(strcasecmp(GenData(NSI_Pressure_Node_codes),"Off")!=0) 
  CODES, NODES, PRESSURE
*if(strcasecmp(GenData(NSI_Pressure_Node_value_code_01),"Off")!=0) 
$    *GenData(NSI_Pressure_Node_name_code_01) 
     *GenData(NSI_Pressure_Node_value_code_01) 
*endif
*if(strcasecmp(GenData(NSI_Pressure_Node_value_code_02),"Off")!=0) 
$    *GenData(NSI_Pressure_Node_name_code_02) 
     *GenData(NSI_Pressure_Node_value_code_02) 
*endif
*if(strcasecmp(GenData(NSI_Pressure_Node_value_code_03),"Off")!=0) 
$    *GenData(NSI_Pressure_Node_name_code_03) 
     *GenData(NSI_Pressure_Node_value_code_03) 
*endif
*if(strcasecmp(GenData(NSI_Pressure_Node_value_code_04),"Off")!=0) 
$    *GenData(NSI_Pressure_Node_name_code_04) 
     *GenData(NSI_Pressure_Node_value_code_04) 
*endif
*if(strcasecmp(GenData(NSI_Pressure_Node_value_code_05),"Off")!=0) 
$    *GenData(NSI_Pressure_Node_name_code_05) 
     *GenData(NSI_Pressure_Node_value_code_05) 
*endif
*if(strcasecmp(GenData(NSI_Pressure_Node_value_code_06),"Off")!=0) 
$    *GenData(NSI_Pressure_Node_name_code_06) 
     *GenData(NSI_Pressure_Node_value_code_06) 
*endif
*if(strcasecmp(GenData(NSI_Pressure_Node_value_code_07),"Off")!=0) 
$    *GenData(NSI_Pressure_Node_name_code_07) 
     *GenData(NSI_Pressure_Node_value_code_07) 
*endif
*if(strcasecmp(GenData(NSI_Pressure_Node_value_code_08),"Off")!=0) 
$    *GenData(NSI_Pressure_Node_name_code_08) 
     *GenData(NSI_Pressure_Node_value_code_08) 
*endif
*if(strcasecmp(GenData(NSI_Pressure_Node_value_code_09),"Off")!=0) 
$    *GenData(NSI_Pressure_Node_name_code_09) 
     *GenData(NSI_Pressure_Node_value_code_09) 
*endif
*if(strcasecmp(GenData(NSI_Pressure_Node_value_code_10),"Off")!=0) 
$    *GenData(NSI_Pressure_Node_name_code_10) 
     *GenData(NSI_Pressure_Node_value_code_10) 
*endif
  END_CODES
*endif
*if(strcasecmp(GenData(NSI_Boundary_codes),"Off")!=0) 
  CODES, BOUNDARIES
*if(strcasecmp(GenData(NSI_Boundary_value_code_01),"Off")!=0) 
$    *GenData(NSI_Boundary_name_code_01) 
     *GenData(NSI_Boundary_value_code_01) 
*endif
*if(strcasecmp(GenData(NSI_Boundary_value_code_02),"Off")!=0) 
$    *GenData(NSI_Boundary_name_code_02) 
     *GenData(NSI_Boundary_value_code_02) 
*endif
*if(strcasecmp(GenData(NSI_Boundary_value_code_03),"Off")!=0) 
$    *GenData(NSI_Boundary_name_code_03) 
     *GenData(NSI_Boundary_value_code_03) 
*endif
*if(strcasecmp(GenData(NSI_Boundary_value_code_04),"Off")!=0) 
$    *GenData(NSI_Boundary_name_code_04) 
     *GenData(NSI_Boundary_value_code_04) 
*endif
*if(strcasecmp(GenData(NSI_Boundary_value_code_05),"Off")!=0) 
$    *GenData(NSI_Boundary_name_code_05) 
     *GenData(NSI_Boundary_value_code_05) 
*endif
*if(strcasecmp(GenData(NSI_Boundary_value_code_06),"Off")!=0) 
$    *GenData(NSI_Boundary_name_code_06) 
     *GenData(NSI_Boundary_value_code_06) 
*endif 
*if(strcasecmp(GenData(NSI_Boundary_value_code_07),"Off")!=0) 
$    *GenData(NSI_Boundary_name_code_07) 
     *GenData(NSI_Boundary_value_code_07) 
*endif
*if(strcasecmp(GenData(NSI_Boundary_value_code_08),"Off")!=0) 
$    *GenData(NSI_Boundary_name_code_08) 
     *GenData(NSI_Boundary_value_code_08) 
*endif
*if(strcasecmp(GenData(NSI_Boundary_value_code_09),"Off")!=0) 
$    *GenData(NSI_Boundary_name_code_09) 
     *GenData(NSI_Boundary_value_code_09) 
*endif
*if(strcasecmp(GenData(NSI_Boundary_value_code_10),"Off")!=0) 
$    *GenData(NSI_Boundary_name_code_10) 
     *GenData(NSI_Boundary_value_code_10)  
*endif
  END_CODES
*endif
  INCLUDE *tcl(AlyaProblemSolo).nsi.fix
END_BOUNDARY_CONDITIONS  
$------------------------------------------------------------
*endif
