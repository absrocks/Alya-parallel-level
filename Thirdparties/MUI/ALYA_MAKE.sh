module purge
mdule load PYTHON/2.7.3
module load  gcc/4.9.1 intel/15.0.2 openmpi/1.8.1

make -j4
#
# config.in 
#   F90      = mpif90 -DMUI
#   EXTRALIB = /home/bsc21/bsc21704/z2016/REPOSITORY/MUI/MUI_2016AUG22/wrapper_f3/mui_3df.o -lmpi_cxx -lstdc++
# 
