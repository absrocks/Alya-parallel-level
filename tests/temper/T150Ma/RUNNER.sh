#BSUB -n 64  
#BSUB -R"span[ptile=16]"
#BSUB -o out00.run
#BSUB -e err00.run 
#BSUB -J R0040 
#BSUB -W 00:59  
# #BSUB -q bsc_case  

ALYA_PATH=/home/bsc21/bsc21704/z2016/REPOSITORY/ALYAs/ALYA_2016Mar12/
ALYA_PATH=/home/bsc21/bsc21704/z2016/REPOSITORY/ALYAs/ALYA_2016JUL18/
ALYA_PATH=/home/bsc21/bsc21704/z2017/REPOSITORY/ALYA_2017JAN10/Executables/unix
ALYA_PATH=/home/bsc21/bsc21704/z2017/REPOSITORY/ALYA_2017MAY09_MNT/Executables/unix

#-------------------------------------------------------------------------||--# 
#
ALYAi=3
ALYAj=3 
#
#-------------------------------------------------------------------------||--# 
find . -type l -delete
#
CASEi=fluid
CASEj=solid
#
ln -s cou.dat  $CASEi.cou.dat
ln -s cou.dat  $CASEj.cou.dat
#
#ALYAij=64 
#ALYAj=4 
#ALYAi=$((ALYAij-ALYAj))
#
time mpirun \
-np $ALYAi  \
 $ALYA_PATH/Alya.x $CASEi  \
 : \
-np $ALYAj \
 $ALYA_PATH/Alya.x $CASEj 
#
find . -type l -delete
#-------------------------------------------------------------------------||--# 
#
#NOTA:
python TOOLs/analizeSets02.py -F fluid 
#
