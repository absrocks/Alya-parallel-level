
# /home/bsc21/bsc21704/z2016/REPOSITORY/HPCTOOLKIT/PAPI543/build_gcc/Exec/bin
# ./papi_avail 
#
# PAPI_FP_OPS  0x80000066  Yes   Yes  Floating point operations
# PAPI_TOT_CYC 0x8000003b  Yes   No   Total cycles
# PAPI_FP_INS  0x80000034  Yes   Yes  Floating point instructions
#

Exec01()
{
  mpirun -np 4 hpcprof-mpi  \
  -S $ALYA_PATH/Alya.x.hpcstruct \
  -I $ALYA_PATH/../../Sources/+ \
  -o $1$2_database \
     $1$2_measurements
}

ALYA_PATH=/home/bsc21/bsc21704/z2016/REPOSITORY/ALYAs/ALYA_2016/Executables/plepp03

ALYAi=1 
ALYAj=1  

export PATH=$PATH:/home/bsc21/bsc21704/z2016/REPOSITORY/HPCTOOLKIT/HPCTK/hpctoolkit/build_gcc_papi/Execs/bin

Exec01 'D' $ALYAi 
Exec01 'N' $ALYAj 
#Exec01 'P' $((ALYAi+ALYAj)) 
 
rm hpc.tar
tar -cvf HPC.tar *_database *_measurements

