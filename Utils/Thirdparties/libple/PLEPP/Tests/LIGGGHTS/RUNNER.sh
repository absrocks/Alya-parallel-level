# 
# 2015feb05, BSC, Barcelona, Spain 
# 2017ABR24, 
#           la compilacion funciona, pero deje las pruebas por falta de tiempo. 
# 
##  0
#compile libple and libplepp -> COMMDOM_PATH  

# git clone https://github.com/CFDEMproject/LIGGGHTS-PUBLIC.git

## 1 
# cd LIGGGHTS/src
# vim MAKE/Makefile.openmpi
# 
## 2
#make clean 


## 3.1
make openmpi -j 4 CC=mpic++ LINK=mpic++ LMP_INC="-DLAMMPS_GZIP -DQUIP_GFORTRAN" FFT_INC='' FFT_PATH='' FFT_LIB='' COMMDOM_PATH=/home/bsc21/bsc21704/z2016/REPOSITORY/ALYAs/ALYA_2016Feb23/Thirdparties/libple/PLEPP/

## >> size ../lmp_openmpi
## >>    text	   data	    bss	    dec	    hex	filename
## >> 17094044	 539432	  21344	17654820	10d6424	../lmp_openmpi

## 3.2
##make openmpi -j 4 CC=mpic++ LINK=mpic++ LMP_INC="-DLAMMPS_GZIP -DQUIP_GFORTRAN" FFT_INC='' FFT_PATH='' FFT_LIB='' SATURNE=/home/bsc21/bsc21704/z2016/REPOSITORY/SATURNE334/ COMMDOM_PATH=/home/bsc21/bsc21704/z2016/REPOSITORY/ALYAs/ALYA_2016Feb23/Thirdparties/libple/PLEPP/

## 3.3 
# ln -s lmp_openmpi Ligghts.x 


## 4.0 
cd /LIGGGHTS/examples/LIGGGHTS/Tutorials_public/cohesion


#
#python -i /home/bsc21/bsc21704/z2016/REPOSITORY/LIGGGHTS/LPP_2016FEB05/src/lpp.py --cpunum 4 --chunksize 1
