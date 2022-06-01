$PHYSICAL_PROBLEM                                                           
  PROBLEM_DEFINITION                                                         
    TEMPORAL_DERIVATIVES:   On                                               
    CONVECTIVE_TERM:        On , VELOC: NASTIN                               
    CONDUCTION_TERM:        On                                               
    TURBULENCE_MODEL:       from_turbul                                      
  END_PROBLEM_DEFINITION                                                     
  PROPERTIES
    GAS_CONSTANT = 287.0                                                     
    DENSITY=                    EXTERNAL                                     
    TURBULENT_PRANDTL=          0.9                                          
    LAW_SPECIFIC_HEAT=          EXTERNAL                                     
    LAW_THERMAL_COND=           EXTERNAL                                     
    LAW_VISCOSITY=              EXTERNAL                                     
  END_PROPERTIES                                                             
$END_PHYSICAL_PROBLEM
