
#BSUB -n 288  
#BSUB -R"span[ptile=16]"
#BSUB -o out.hpc   
#BSUB -e err.hpc
#BSUB -J HPC   
#BSUB -W 00:29   
# #BSUB -q bsc_case  

export PATH=$PATH:/home/bsc21/bsc21704/z2016/REPOSITORY/HPCTOOLKIT/HPCTK/hpctoolkit/build_gcc_papi/Execs/bin

ALYA_PATH=/home/bsc21/bsc21704/z2016/REPOSITORY/ALYAs/ALYA_2016/Executables/plepp03
#ALYA_PATH=/home/bsc21/bsc21704/z2016/REPOSITORY/ALYAs/ALYA_2016DEC23/Executables/plepp03
#ALYA_PATH=/home/bsc21/bsc21704/z2017/REPOSITORY/ALYA_2017JAN05/Executables/plepp03
ALYA_PATH=/home/bsc21/bsc21704/z2017/REPOSITORY/ALYA_2017JAN07/Executables/plepp03b 
ALYA_PATH=/home/bsc21/bsc21704/z2017/REPOSITORY/ALYA_2017JAN10/Executables/plepp-3
ALYA_PATH=/home/bsc21/bsc21704/z2017/REPOSITORY/ALYA_2017MAY09_MNT/Executables/plepp-3


ALYAi=3  
ALYAj=3 
ALYAij=$((ALYAi+ALYAj))

time mpirun \
-np $ALYAi  \
 $ALYA_PATH/Alya.x fluid --name DIRIC  \
 : \
-np $ALYAj \
 $ALYA_PATH/Alya.x solid --name NEUMA

