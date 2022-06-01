#
# http://ix.cs.uoregon.edu/~khuck/tau-faq/
# https://wiki.mpich.org/mpich/index.php/TAU_by_example
# https://www.alcf.anl.gov/user-guides/tuning-and-analysis-utilities-tau   
# https://www.cs.uoregon.edu/research/tau/home.php
# http://nic.uoregon.edu/pipermail/tau-users/2015-February/000961.html
# 

## 1)
#./configure -ICPC -prefix=/home/bsc21/bsc21704/z2016/REPOSITORY/TAU/PDT322/EXEC
# make 
# make install 

## 2) 
./configure -mpi \
-cc=mpicc \
-c++=mpic++ \
-pdt=/home/bsc21/bsc21704/z2016/REPOSITORY/TAU/PDT322/EXEC/ \
-prefix=/home/bsc21/bsc21704/z2016/REPOSITORY/TAU/TAU2251/EXEC2 
#-bfd=download \ ??

#
#-DISABLESHARED 
#NOTE: Not building TAU's DSOs
#

# http://nic.uoregon.edu/pipermail/tau-users/2015-February/000961.html
#
# with out -fortran
#Setting F90 compiler based on requested: intel
#Default Fortran compiler will be Intel ifort

# using: 
#-fortran=mpif90 \
# Setting F90 compiler based on requested: mpif90

## 3) 
#make install 


## 4)
#cd EXEC 
#../tau_validate --html x86_64 &> results.html 

## 5) 
#export TAU_PATH:=/home/bsc21/bsc21704/z2016/REPOSITORY/TAU/TAU2251/EXEC2/x86_64
#export TAU_MAKEFILE:=$(TAU_PATH)/lib/Makefile.tau-icpc-mpi-pdt
#export PATH:=$(PATH):$(TAU_PATH)/bin
#pprof
#paraprof profile.* 

## RING 
#export TAU_COMM_MATRIX=1

## Tracing with Jumpshot 
#export TAU_TRACE=1
#tau_treemerge.pl 
#tau2slog2 tau.trc tau.edf -o tau.slog2 

# https://www.cs.uoregon.edu/research/tau/docs/newguide/bk01ch04s03.html

 
