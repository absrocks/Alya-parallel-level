#!/bin/bash
##BSUB -J IQN01 
##BSUB -W 00:05 
##BSUB -n 17 
# ##BSUB -q sequential
##BSUB -eo out01.err
##BSUB -oo out01.out

ALYA_PATH=/home/bsc21/bsc21704/z2017/REPOSITORY/ALYA_2017ABR05_MNT/Executables/unix


ALYAij=6  
ALYAi=3    

#ALYAij=16   
#ALYAi=8  
ALYAj=$((ALYAij-ALYAi)) 

CASEi=esbelta_fld
CASEj=esbelta_sld


COU()
{ 
  rm *.cou.dat *.npart *.tms  *.ncou
  # 
  ln -s cou.alya $CASEi.cou.dat
  ln -s cou.alya $CASEj.cou.dat
  #
  time mpirun \
  -np $ALYAi \
   $ALYA_PATH/Alya.x $CASEi \
   : \
  -np $ALYAj \
   $ALYA_PATH/Alya.x $CASEj   
  #
  rm *.cou.dat 
  #
  ls *.npart *.tms  *.ncou
  #
}


ALONEi()
{
  rm *.cou.dat *.npart *.tms  *.ncou  
  #
  time mpirun \
  -np $ALYAi \
   $ALYA_PATH/Alya.x $CASEi 

  ls *.npart *.tms  *.ncou  
}


#
COU
#ALONEi

# module load python  
# python /home/bsc21/bsc21704/z2017/RUNNER/HPC02/DLR2D_02/CHT_PLEPP04/checkTimeMeasuraments01.py -F "."
 
