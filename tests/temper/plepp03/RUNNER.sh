#BSUB -n 32  
#BSUB -R"span[ptile=16]"
#BSUB -o OUT00.txt
#BSUB -e ERR00.txt
#BSUB -J 2016MAY10  
#BSUB -W 00:29  
# #BSUB -q bsc_case  

ALYA_PATH=/home/bsc21/bsc21704/z2016/REPOSITORY/ALYAs/ALYA_2016Mar12/

ALYAi=4  
ALYAj=1  

#ALYAi=8
#ALYAj=24

ln -s cou.dat  Interior01.cou.dat
ln -s cou.dat    vortex2D.cou.dat
#
time mpirun -np $ALYAj $ALYA_PATH/Executables/unix/Alya.x Interior01 : -np $ALYAi $ALYA_PATH/Executables/unix/Alya.x vortex2D 
#
rm *.cou.dat 

#
#  2016MAY02: CPLNG01_C01 <- CPLNG01_01   
#    
