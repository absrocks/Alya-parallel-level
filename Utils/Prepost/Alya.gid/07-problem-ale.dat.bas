*if(strcasecmp(GenData(ALEFOR_Module),"on")==0)
*realformat "%16.6f"
*intformat "%8i"
$------------------------------------------------------------
PHYSICAL_PROBLEM
$------------------------------------------------------------
PROBLEM_DEFINITION

END_PROBLEM_DEFINITION
$------------------------------------------------------------
END_PHYSICAL_PROBLEM
$------------------------------------------------------------
NUMERICAL_TREATMENT
END_NUMERICAL_TREATMENT
$------------------------------------------------------------
OUTPUT_&_POST_PROCESS
  POSTPROCESS MESH_DISPLACEMENT, STEPS: *GenData(ALE_Postprocess_mesh_displacement)
  POSTPROCESS MESH_VELOCITY,     STEPS: *GenData(ALE_Postprocess_mesh_velocity)
END_OUTPUT_&_POST-PROCESS
$------------------------------------------------------------
BOUNDARY_CONDITIONS
   INCLUDE *tcl(AlyaProblemName).ale.fix
END_BOUNDARY_CONDITIONS
$------------------------------------------------------------
*endif
