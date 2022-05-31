*if(strcasecmp(GenData(Codire_Module),"Off")!=0)
*realformat "%16.6f"
*intformat "%8i"
$------------------------------------------------------------
PHYSICAL_PROBLEM
 TEMPORAL_DERIVATIVES: *GenData(CDR_temporal_derivatives)
*if(strcasecmp(GenData(CDR_problem),"General")==0)
  MASS_MATRIX:               "constant"
                             *Gendata(CDR_Mass_matrix)
  RES_DIFFUSION_MATRICES:    "constant"
                             *Gendata(CDR_Residual_diffusion_matrices)
  RES_CONVECTION_MATRICES:   "constant" 
                             *Gendata(CDR_Residual_convection_matrices)
  RES_FLUX_MATRICES:         "constant" 
                             *Gendata(CDR_Residual_flux_matrices)
  RES_SOURCE_MATRICES:       "constant" 
                             *Gendata(CDR_Residual_source_matrices)
  TES_DIFFUSION_MATRICES:    "constant"
                             *Gendata(CDR_Test_diffusion_matrices)
  TES_CONVECTION_MATRICES:   "constant" 
                             *Gendata(CDR_Test_convection_matrices)
  TES_SOURCE_MATRICES:       "constant" 
                             *Gendata(CDR_Test_source_matrices)
  FORCE_VECTOR:              "constant"
                             *Gendata(CDR_force_vector)
  DtN_MATRICES:              "constant" 
                             *Gendata(CDR_dtn_matrices)
*elseif(strcasecmp(GenData(CDR_problem),"Navier_Stokes")==0)
$------------------------------------------------------------
  SECOND_ORDER_TERM:         *Gendata(CDR_Second_order_term)
  LAW_VISCOSITY:             *Gendata(CDR_Law_viscosity)
*if(strcasecmp(GenData(CDR_Law_viscosity),"Constant")==0) 
  VISCOSITY:                 *Gendata(CDR_Viscosity_01)
*else
  VISCOSITY:                 *Gendata(CDR_Viscosity_01) *Gendata(CDR_Viscosity_02) *Gendata(CDR_Viscosity_03) *Gendata(CDR_Viscosity_04) *Gendata(CDR_Viscosity_05) *Gendata(CDR_Viscosity_06) *Gendata(CDR_Viscosity_07) *Gendata(CDR_Viscosity_08) *Gendata(CDR_Viscosity_09) *Gendata(CDR_Viscosity_10)
*endif
  DENSITY:                   *Gendata(CDR_Density)
  PERMEABILITY:              *Gendata(CDR_Permeability)
*if(strcasecmp(GenData(CDR_Angular_velocity),"Off")==0) 
  ANGULAR_VELOCITY:          NORM: 0.0, GX: 0.0, GY: 0.0,  GZ: 0.0  
*else
  ANGULAR_VELOCITY:          NORM: *GenData(CDR_Norm_angular_velocity), GX: *GenData(CDR_X_angular_velocity), GY: *GenData(CDR_Y_angular_velocity), GZ: *GenData(CDR_Z_angular_velocity)  
*endif
  PENALTY_PARAMETER:         *Gendata(CDR_Penalty_parameter)
  SPEED_OF_SOUND:            *Gendata(CDR_speed_of_sound)
*elseif(strcasecmp(GenData(CDR_problem),"Plates")==0)
$------------------------------------------------------------
  THICKNESS:                 *Gendata(CDR_Thickness)
  YOUNG_MODULUS:             *Gendata(CDR_Young_modulus)
  POISSON_COEFFICIENT:       *Gendata(CDR_Poisson_coefficient)
  SHEAR_CORRECTION:          *Gendata(CDR_Shear_correction)
  NORMALIZED_LOAD:           *Gendata(CDR_Normalized_load)
*elseif(strcasecmp(GenData(CDR_problem),"Helmholtz")==0)
$------------------------------------------------------------
  WAVE_NUMBER:               *Gendata(CDR_Wave_number)
  CONSTANT_FORCE:            *Gendata(CDR_Constant_force)
*elseif(strcasecmp(GenData(CDR_problem),"Boussinesq")==0)
$------------------------------------------------------------
  SECOND_ORDER_TERM:         *Gendata(CDR_Second_order_term)
  LAW_VISCOSITY:             *Gendata(CDR_Law_viscosity)
*if(strcasecmp(GenData(CDR_Law_viscosity),"Constant")==0) 
  VISCOSITY:                 *Gendata(CDR_Viscosity_01)
*else
  VISCOSITY:                 *Gendata(CDR_Viscosity_01) *Gendata(CDR_Viscosity_02) *Gendata(CDR_Viscosity_03) *Gendata(CDR_Viscosity_04) *Gendata(CDR_Viscosity_05) *Gendata(CDR_Viscosity_06) *Gendata(CDR_Viscosity_07) *Gendata(CDR_Viscosity_08) *Gendata(CDR_Viscosity_09) *Gendata(CDR_Viscosity_10)
*endif
  DENSITY:                   *Gendata(CDR_Density)
  PERMEABILITY:              *Gendata(CDR_Permeability)
*if(strcasecmp(GenData(CDR_Angular_velocity),"Off")==0) 
  ANGULAR_VELOCITY:          NORM: 0.0, GX: 0.0, GY: 0.0,  GZ: 0.0  
*else
  ANGULAR_VELOCITY:          NORM: *GenData(CDR_Norm_angular_velocity), OX: *GenData(CDR_X_angular_velocity), OY: *GenData(CDR_Y_angular_velocity), OZ: *GenData(CDR_Z_angular_velocity)  
*endif
  PENALTY_PARAMETER:         *Gendata(CDR_Penalty_parameter)
  SPEED_OF_SOUND:            *Gendata(CDR_speed_of_sound)
*if(strcasecmp(GenData(CDR_Gravity_acceleration),"Off")==0) 
  GRAVITY:                   NORM: 0.0, GX: 0.0, GY: 0.0,  GZ: 0.0  
*else
  GRAVITY:                   NORM: *GenData(CDR_Norm_gravity), GX: *GenData(CDR_X_gravity), GY: *GenData(CDR_Y_gravity), GZ: *GenData(CDR_Z_gravity)  
*endif
  SPECIFIC_HEAT:             *GenData(CDR_Specific_heat)
  THERMAL_CONDUCTIVITY       *GenData(CDR_Thermal_conductivity)
  DILATATION_COEFFICIENT:    *GenData(CDR_Dilatation_coefficient)
  REFERENCE_TEMPERATURE:     *GenData(CDR_Reference_temperature)
*elseif(strcasecmp(GenData(CDR_problem),"LowMach")==0)
$------------------------------------------------------------
*if(strcasecmp(GenData(CDR_Low_Mach_system),"System1")==0) 
  SYSTEM:                    1
*else
  SYSTEM:                    2
*endif
  LAW_VISCOSITY:             *Gendata(CDR_Law_viscosity)
*if(strcasecmp(GenData(CDR_Law_viscosity),"Constant")==0) 
  VISCOSITY:                 *Gendata(CDR_Viscosity_01)
*else
  VISCOSITY:                 *Gendata(CDR_Viscosity_01) *Gendata(CDR_Viscosity_02) *Gendata(CDR_Viscosity_03) *Gendata(CDR_Viscosity_04) *Gendata(CDR_Viscosity_05) *Gendata(CDR_Viscosity_06) *Gendata(CDR_Viscosity_07) *Gendata(CDR_Viscosity_08) *Gendata(CDR_Viscosity_09) *Gendata(CDR_Viscosity_10)
*endif
  DENSITY:                   *Gendata(CDR_Density)
  PERMEABILITY:              *Gendata(CDR_Permeability)
*if(strcasecmp(GenData(CDR_Angular_velocity),"Off")==0) 
  ANGULAR_VELOCITY:          NORM: 0.0, OX: 0.0, OY: 0.0,  OZ: 0.0  
*else
  ANGULAR_VELOCITY:          NORM: *GenData(CDR_Norm_angular_velocity), OX: *GenData(CDR_X_angular_velocity), OY: *GenData(CDR_Y_angular_velocity), OZ: *GenData(CDR_Z_angular_velocity)  
*endif
  PENALTY_PARAMETER:         *Gendata(CDR_Penalty_parameter)
  SPEED_OF_SOUND:            *Gendata(CDR_speed_of_sound)
*if(strcasecmp(GenData(CDR_Gravity_acceleration),"Off")==0) 
  GRAVITY:                   NORM: 0.0, GX: 0.0, GY: 0.0,  GZ: 0.0  
*else
  GRAVITY:                   NORM: *GenData(CDR_Norm_gravity), GX: *GenData(CDR_X_gravity), GY: *GenData(CDR_Y_gravity), GZ: *GenData(CDR_Z_gravity)  
*endif
  SPECIFIC_HEAT:             *GenData(CDR_Specific_heat)
  THERMAL_CONDUCTIVITY       *GenData(CDR_Thermal_conductivity)
  GAMMA:                     *GenData(CDR_Gamma)
  TH_PRE_EQUATION:           *GenData(CDR_Th_pre_equation)
  INITIAL_THPRESS:           *GenData(CDR_Initial_th_pressure)
*endif
END_PHYSICAL_PROBLEM 
$------------------------------------------------------------
NUMERICAL_TREATMENT    
 STABILITY_CONSTANTS:     *GenData(CDR_Tau_1) *GenData(CDR_Tau_2) *GenData(CDR_Tau_3) *GenData(CDR_Tau_4)  *GenData(CDR_Tau_5)  *GenData(CDR_Tau_6)  *GenData(CDR_Tau_7)  *GenData(CDR_Tau_8)
 TYPE_OF_STABILIZATION:   *GenData(CDR_Type_of_stabilization)  
 SHOCK_CAPTURING:         *GenData(CDR_Shock_capturing), VALUE: *GenData(CDR_Value_sh) 
 TEMPORAL_TERM_WEIGHTING: *GenData(CDR_Temporal_Term_weighting)
*if(strcasecmp(GenData(CDR_Time_integration_order),"First")==0) 
 TIME_ACCURACY:           1
*else
 TIME_ACCURACY:           2
*endif
*if(strcasecmp(GenData(CDR_Time_integration_type),"Gear")==0)
 TIME_INTEGRATION:        Gear
*else
 TIME_INTEGRATION:        Trapezoidal
*endif
 SAFETY_FACTOR:           *GenData(CDR_Safety_factor)  
 TRACKING_OF_SUBSCALES:   *GenData(CDR_Tracking_of_subscales)
 STEADY_STATE_TOLERANCE:  *GenData(CDR_Steady_state_tolerance)   
 NORM_OF_CONVERGENCE:     *GenData(CDR_Norm_of_convergence)
*if(strcasecmp(GenData(CDR_Linearization_method),"Picard")==0) 
 LINEARIZATION_METHOD: Picard   
*else
 LINEARIZATION_METHOD:    *GenData(CDR_Linearization_method), PICARD_ITERATIONS: *GenData(CDR_Picard_Iterations)  \
                          PARAMETERS: *GenData(CDR_Lin_par_01) *GenData(CDR_Lin_par_02) *GenData(CDR_Lin_par_03) *GenData(CDR_Lin_par_04) *GenData(CDR_Lin_par_05) *GenData(CDR_Lin_par_06) *GenData(CDR_Lin_par_07) *GenData(CDR_Lin_par_08) *GenData(CDR_Lin_par_09) *GenData(CDR_Lin_par_10)
*endif
*if(strcasecmp(GenData(CDR_Line_search),"Test_convergence")==0) 
 LINEARIZATION_METHOD:    Test_convergence
*elseif(strcasecmp(GenData(CDR_Line_search),"Armijo")==0) 
 LINEARIZATION_METHOD:    *GenData(CDR_Line_search), PARAMETER: *GenData(CDR_LSparameter)  \
                          ITERATIONS: *GenData(CDR_LSiterations)
*endif
 MAXIMUM_NUMBER_OF_ITER:  *GenData(CDR_Number_of_internal_iterations)  
 CONVERGENCE_TOLERANCE:   *GenData(CDR_Convergence_tolerance)   
*if(strcasecmp(GenData(CDR_Algebraic_solver),"GMRES")==0)
 ALGEBRAIC_SOLVER:        *GenData(CDR_Algebraic_solver)  
 SOLVER_ITERATIONS:       *GenData(CDR_Solver_iterations)  
 TOLERANCE_OF_CONVER:     *GenData(CDR_Solver_tolerance)  
 KRYLOV_DIMENSION:        *GenData(CDR_Krylov_dimension)     
*else
 ALGEBRAIC_SOLVER:        *GenData(CDR_Algebraic_solver)  
*endif
END_NUMERICAL_TREATMENT  
$------------------------------------------------------------
OUTPUT_&_POST-PROCESS  
 START_POSTPROCES_AT      STEP : *GenData(CDR_Start_Postprocess_at_step)
*if(GenData(CDR_Postprocess_unknowns,int)>0)
 POSTPROCESS UNKNOWNS,    STEPS: *GenData(CDR_Postprocess_unknowns)
*endif
*if(GenData(CDR_Postprocess_streamlines,int)>0)
 POSTPROCESS STREAMLINES, STEPS: *GenData(CDR_Postprocess_streamlines)
*endif
*if(strcasecmp(GenData(CDR_Tracking_of_points),"On")==0)
 TRACKING
*if(GenData(CDR_Nodal_Points_tracked,int)>0)
   NODAL_POINTS:          *GenData(CDR_Nodal_Points_tracked)
*endif
*if(GenData(CDR_Points_tracked,int)==1)
   POINTS:                1 *GenData(CDR_1_Pt_X) *GenData(CDR_1_Pt_Y) *GenData(CDR_1_Pt_Z)
*elseif(GenData(CDR_Points_tracked,int)==2)
   POINTS:                2 \ 
                          *GenData(CDR_1_Pt_X) *GenData(CDR_1_Pt_Y) *GenData(CDR_1_Pt_Z) \
                          *GenData(CDR_2_Pt_X) *GenData(CDR_2_Pt_Y) *GenData(CDR_2_Pt_Z)
*elseif(GenData(CDR_Points_tracked,int)==3)
   POINTS:                3 \ 
                          *GenData(CDR_1_Pt_X) *GenData(CDR_1_Pt_Y) *GenData(CDR_1_Pt_Z) \
                          *GenData(CDR_2_Pt_X) *GenData(CDR_2_Pt_Y) *GenData(CDR_2_Pt_Z) \
                          *GenData(CDR_3_Pt_X) *GenData(CDR_3_Pt_Y) *GenData(CDR_3_Pt_Z)
*elseif(GenData(CDR_Points_tracked,int)==4)
   POINTS:                4 \ 
                          *GenData(CDR_1_Pt_X) *GenData(CDR_1_Pt_Y) *GenData(CDR_1_Pt_Z) \
                          *GenData(CDR_2_Pt_X) *GenData(CDR_2_Pt_Y) *GenData(CDR_2_Pt_Z) \
                          *GenData(CDR_3_Pt_X) *GenData(CDR_3_Pt_Y) *GenData(CDR_3_Pt_Z) \
                          *GenData(CDR_4_Pt_X) *GenData(CDR_4_Pt_Y) *GenData(CDR_4_Pt_Z)
*elseif(GenData(CDR_Points_tracked,int)==5)
   POINTS:                5 \ 
                          *GenData(CDR_1_Pt_X) *GenData(CDR_1_Pt_Y) *GenData(CDR_1_Pt_Z) \
                          *GenData(CDR_2_Pt_X) *GenData(CDR_2_Pt_Y) *GenData(CDR_2_Pt_Z) \
                          *GenData(CDR_3_Pt_X) *GenData(CDR_3_Pt_Y) *GenData(CDR_3_Pt_Z) \
                          *GenData(CDR_4_Pt_X) *GenData(CDR_4_Pt_Y) *GenData(CDR_4_Pt_Z) \
                          *GenData(CDR_5_Pt_X) *GenData(CDR_5_Pt_Y) *GenData(CDR_5_Pt_Z)
*endif
 END_TRACKING
*endif
*if(strcasecmp(GenData(CDR_Plotting_of_sections),"On")==0)
 PLOT_SECTIONS, STEPS:   *GenData(CDR_Section_every_steps)
*endif
*if(strcasecmp(GenData(CDR_Compute_forces),"On")==0)
 FORCES_AND_MOMENTS
  STEPS: *GenData(CDR_Forces_every_steps)
*if(strcasecmp(GenData(CDR_Output_fluxes),"On")==0)
  OUTPUT_FLUXES
*endif
  FORCE_COEFFICIENTS:  X_FACTOR: *GenData(CDR_X_Force_coefficient) Y_FACTOR: *GenData(CDR_Y_Force_coefficient) Z_FACTOR: *GenData(CDR_Z_Force_coefficient)
  MOMENT_COEFFICIENTS: X_FACTOR: *GenData(CDR_X_Moment_coefficient) Y_FACTOR: *GenData(CDR_Y_Moment_coefficient) Z_FACTOR: *GenData(CDR_Z_Moment_coefficient)
  APPLICATION_POINT:   X_COORD:  *GenData(CDR_X_Moment_application_point) Y_COORD:  *GenData(CDR_Y_Moment_application_point) Z_COORD:  *GenData(CDR_Z_Moment_application_point)
  BODY_DEFINITION
    INCLUDE *tcl(AlyaProblemName).bod.dat
  END_BODY_DEFINITION
 END_FORCES_AND_MOMENTS
*endif
END_OUTPUT_&_POST_PROCESS  
$------------------------------------------------------------
BOUNDARY_CONDITIONS
 INCLUDE *tcl(AlyaProblemName).cdr.fix
END_BOUNDARY_CONDITIONS  
$------------------------------------------------------------
*endif
