#BSUB -n 17    
#BSUB -R"span[ptile=16]"
#BSUB -o OUT00.txt
#BSUB -e ERR00.txt
#BSUB -J 2016MAY03  
#BSUB -W 00:10 
# #BSUB -q bsc_case  

ALYAi=1
ALYAj=1

ALYA_PATH=/home/bsc21/bsc21704/z2016/REPOSITORY/ALYAs/ALYA_2016May03/
ALYA_PATH=/home/bsc21/bsc21704/z2017/REPOSITORY/ALYA_2017JAN09/
ALYA_PATH=../../Alya/

DUMMY_PATH=$ALYA_PATH/Thirdparties/libple/PLEPP/Wrappers/Cpp/

time mpirun -np $ALYAi $ALYA_PATH/Executables/unix/Alya_libple_four.x vortex2D --name NEUMA : -np  $ALYAj $DUMMY_PATH/wloader.x DIRIC

##time mpirun -np $ALYAi $ALYA_PATH/Executables/plepp04/Alya.x Interior01 --name DIRIC : -np $ALYAj $ALYA_PATH/Executables/plepp04/Alya.x vortex2D --name NEUMA

#
#
#          3
#    -----------
#   |    _      |
# 1 |   (_)4    | 2
#   |           |
#    -----------
#
#
#    y
#    |
#    |___ x
#   /
#  /
# z
# 
# 
