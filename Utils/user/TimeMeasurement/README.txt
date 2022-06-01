A) NO ACOPLADO 

0.0) Sources/kernel/coupli/mod_coupling_timer.f90 

0.1) NO acoplado/acoplado   
     COU_STATISTICS = 0_ip -> COU_STATISTICS = 1_ip  


1.0) EXEC   

1.1) run (RUNNER.sh+ALONEi)  

1.2) ls *.npart *.tms  *.ncou

1.3) ANALISE  
     module load python
     python checkTimeMeasuraments01.py -F "." 
     ls diric_*  neuma_* quartile50th*  


B) NO ACOPLADO

0.2) acoplado 

     NUMERICAL_TREATMENT
       ... 
       STATISTICS
     END_NUMERICAL_TREATMENT



