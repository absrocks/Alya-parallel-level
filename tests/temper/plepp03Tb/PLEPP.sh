#BSUB -n 32  
#BSUB -R"span[ptile=16]"
#BSUB -o out00.run
#BSUB -e err00.run 
#BSUB -J SerialMPMD  
#BSUB -W 00:29  
# #BSUB -q bsc_case  

#-------------------------------------------------------------------------||--# 
# 
ALYA_PATH=/home/bsc21/bsc21704/z2017/REPOSITORY/ALYA_2017JAN10/Executables/plepp03
# 
ALYAi=4
ALYAj=4 
#
#-------------------------------------------------------------------------||--# 
echo "RUNNIG:" $(date +'%Y-%m-%d') $(date +'%T')
#-------------------------------------------------------------------------||--# 
#
CASEi=fluid
CASEj=solid 
#
ALYAij=32
#ALYAj=4 
#ALYAi=$((ALYAij-ALYAj))
#
time mpirun \
-np $ALYAi  \
 $ALYA_PATH/Alya.x $CASEi --name NEUMA  \
 : \
-np $ALYAj \
 $ALYA_PATH/Alya.x $CASEj --name DIRIC  

#-------------------------------------------------------------------------||--# 
echo "RUNNIG:" $(date +'%Y-%m-%d') $(date +'%T')
#-------------------------------------------------------------------------||--#
# 
# NOTAS:
#      a) LOS TIEMPOS DEBEN DE SER LOS MISMOS!!   
#      b) 
#                i       j
#     CASE:   fluid   solid 
#     CODE:      1       2
#     ALYA:   NEUMA   DIRIC   <- OJO CON EL ORDEN (POSIBLE ERROR EN LA IMPLEMENTACION, PERO FUNCION BIEN TODO) 
#     COUP:      6       3    
#      c) 
#     EN AMBOS CASOS OPTIONS: FIXITY 
#      d) 
#     parallel MPMD 
#
