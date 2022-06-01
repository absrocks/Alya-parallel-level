#BSUB -n 17    
#BSUB -R"span[ptile=16]"
#BSUB -o OUT00.txt
#BSUB -e ERR00.txt
#BSUB -J 2016MAR29  
#BSUB -W 00:29    
# #BSUB -q bsc_case  

module purge
module load PYTHON/2.7.3
module load  gcc/4.9.1 intel/15.0.2 openmpi/1.8.1

export ALYAi=3  
export ALYAj=3    

export ALYA_PATH=/home/bsc21/bsc21704/z2016/REPOSITORY/ALYAs/ALYA_2016AUG24
time mpirun -np $ALYAi $ALYA_PATH/Executables/mui00/Alya.x Interior01 --name NEUMA : -np $ALYAj $ALYA_PATH/Executables/mui00/Alya.x vortex2D --name DIRIC


#
# vortex2D.dat 
#  CODE:               2
#  ALYA:               DIRIC 
# 
# Interior01.dat  
#  CODE:               1
#  ALYA:               NEUMA
# 
