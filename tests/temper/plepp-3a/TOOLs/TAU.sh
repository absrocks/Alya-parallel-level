#  
export TAU_PATH=/home/bsc21/bsc21704/z2016/REPOSITORY/TAU/TAU2251/EXEC2/x86_64
export TAU_MAKEFILE=$TAU_PATH/lib/Makefile.tau-icpc-mpi-pdt
export PATH=$PATH:$TAU_PATH/bin


### Tracing with Jumpshot (export TAU_TRACE=1)  
# 1)
tau_treemerge.pl
# 2)
tau2slog2 tau.trc tau.edf -o tau.slog2
# 3)
tar -czvf tau.tar tau.slog2



Exec01()
{
  rm F$1xS$2.ppk
  paraprof --pack F$1xS$2.ppk
  rm profile.*
  paraprof -v F$1xS$2.ppk

  pprof -t -f profile.Mean            > F$1xS$2.tau   
  pprof -t -f profile.Mean,AllThreads > F$1xS$2.tau2
 #paraprof F$1xS$2.ppk  
}

ALYA_PATH=/home/bsc21/bsc21704/z2016/REPOSITORY/ALYAs/ALYA_2016DEC23/Executables/plepp03

ALYAi=256  
ALYAj=32  

#Exec01 $ALYAi $ALYAj  


