*realformat "%16.6f"
*intformat "%8i"
$-------------------------------------------------------------------
RUN_DATA
  ALYA:                   *GenData(Title)
*if(strcasecmp(GenData(Preliminary_run),"On")==0)
  RUN_TYPE:               Preliminary, Frequency=*GenData(Preliminary_output_frequency) \
*else
  RUN_TYPE:               \
*endif
*if(strcasecmp(GenData(Restart_run),"Continue")==0)
                          Continue_restart
*elseif(strcasecmp(GenData(Restart_run),"Initial")==0)
                          Initial_restart
*elseif(strcasecmp(GenData(Restart_run),"Interpolate_initial")==0)
                          Interpolate_initial
*else
                          No_restart
*endif
END_RUN_DATA
$-------------------------------------------------------------------
PROBLEM_DATA
*if(strcasecmp(GenData(Time_coupling),"Global_from_critical")==0)
  TIME_COUPLING:          Global, From_critical
*elseif(strcasecmp(GenData(Time_coupling),"Global_prescribed")==0)
  TIME_COUPLING:          Global, Prescribed
*elseif(strcasecmp(GenData(Time_coupling),"Local")==0)
  TIME_COUPLING:          Local
*endif
  TIME_INTERVAL=          *GenData(Time_initial), *GenData(Time_final)  
  TIME_STEP_SIZE=         *GenData(Time_step_size)  
  NUMBER_OF_STEPS=        *GenData(Maximum_number_of_steps)   
*if(strcasecmp(GenData(Block_iterations),"Off")!=0)
  MAXIMUM_NUMBER_GLOBAL=  *GenData(Maximum_block_coupling_iterations)
  BLOCK_ITERATION: *GenData(Block_iterations)
    1 *GenData(Block_1)
*if(strcasecmp(GenData(Block_iterations),"1")!=0)
    2 *GenData(Block_2)
*if(strcasecmp(GenData(Block_iterations),"2")!=0)
    3 *GenData(Block_3)
*endif
*endif
  END_BLOCK_ITERATION
*else
  MAXIMUM_NUMBER_GLOBAL=  *GenData(Maximum_coupling_iterations)
*endif
*\
*\ NASTIN ----------------------------------------------------------
*\ 
*if(strcasecmp(GenData(NASTIN_Module),"Off")!=0)
  NASTIN_MODULE:          *GenData(NASTIN_Module)
  DELAY STEPS:            *GenData(NASTIN_Delay_Steps)
*if(strcasecmp(GenData(NASTIN_Convergence),"Required")==0)
    CONVERGENCE_REQUIRED: Yes
*else
    CONVERGENCE_REQUIRED: No
*endif
  END_NASTIN_MODULE
*endif
*\
*\ TEMPER ----------------------------------------------------------
*\ 
*if(strcasecmp(GenData(TEMPER_Module),"Off")!=0)
  TEMPER_MODULE:          *GenData(TEMPER_Module)
  DELAY STEPS:            *GenData(TEMPER_Delay_Steps)
*if(strcasecmp(GenData(TEMPER_Convergence),"Required")==0)
    CONVERGENCE_REQUIRED: Yes
*else
    CONVERGENCE_REQUIRED: No
*endif
  END_TEMPER_MODULE
*endif
*\
*\ TURBUL ----------------------------------------------------------
*\ 
*if(strcasecmp(GenData(TURBUL_Module),"Off")!=0)
  TURBUL_MODULE:          *GenData(TURBUL_Module)
  DELAY STEPS:            *GenData(TURBUL_Delay_Steps)
*if(strcasecmp(GenData(TURBUL_Convergence),"Required")==0)
    CONVERGENCE_REQUIRED: Yes
*else
    CONVERGENCE_REQUIRED: No
*endif
  END_TURBUL_MODULE
*endif
*\
*\ PARALL ----------------------------------------------------------
*\ 
*if(strcasecmp(GenData(PARALL_Service),"Off")!=0)
  PARALL_SERVICE:         *GenData(PARALL_Service)
    OUTPUT:               *GenData(PARALL_Output)
    POSTPROCESS:          *GenData(PARALL_Postprocess)
    TASK:                 *GenData(PARALL_Task), SUBDOMAINS=*GenData(PARALL_Subdomains), *GenData(PARALL_File_format)
    PARTITION:            *GenData(PARALL_Partition)
    FILE_HIERARCHY:       *GenData(PARALL_File_hierarchy) 
  END_PARALL_SERVICE
*endif
*\
END_PROBLEM_DATA
$-------------------------------------------------------------------
MPI_IO: Off
  GEOMETRY:    Off
  RESTART:     Off
  POSTPROCESS: Off
END_MPI_IO
$-------------------------------------------------------------------
