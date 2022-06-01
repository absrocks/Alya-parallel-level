#BSUB -n 17    
#BSUB -R"span[ptile=16]"
#BSUB -o OUT00.txt
#BSUB -e ERR00.txt
#BSUB -J 2016MAR29  
#BSUB -W 00:29    
# #BSUB -q bsc_case  

ALYA_PATH=/home/bsc21/bsc21704/z2016/REPOSITORY/ALYAs/ALYA_2016Mar29
ALYA_PATH=/home/bsc21/bsc21704/z2016/REPOSITORY/ALYAs/ALYA_2016MAY27/
ALYA_PATH=/home/bsc21/bsc21704/z2016/REPOSITORY/ALYAs/ALYA_2016AUG24/
ALYA_PATH=/home/bsc21/bsc21704/z2016/REPOSITORY/ALYAs/ALYA_2016DEC23/
ALYA_PATH=/home/bsc21/bsc21704/z2017/REPOSITORY/ALYA_2017JAN09/


ALYAi=4 
ALYAj=1    

#time mpirun -np $ALYAi $ALYA_PATH/Executables/unix/Alya.x vortex2D
#time mpirun -np $ALYAi $ALYA_PATH/Executables/unix/Alya.x Interior01 
time mpirun -np $ALYAj $ALYA_PATH/Executables/plepp03/Alya.x Interior01 --name DIRIC : -np $ALYAi $ALYA_PATH/Executables/plepp03/Alya.x vortex2D --name NEUMA

## DEBUG 
#time mpirun \
#-np $ALYAi gdb -ex=r --args $ALYA_PATH/Executables/plepp03/Alya.g Interior01 --name DIRIC : \
#-np $ALYAj gdb -ex=r --args $ALYA_PATH/Executables/plepp03/Alya.g vortex2D --name NEUMA 

## TEST
# Alya2pos.x Interior01  
# vim -d base/1p/Interior01.ensi.GRATE-000002 Interior01.ensi.GRATE-000002 

#
#  ~ 0.3 segs / 1900stp / 30 min / 3 + 12 cores  
# 
