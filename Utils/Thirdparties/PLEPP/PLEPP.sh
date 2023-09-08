##
##  http://bsccase02.bsc.es/TestSuite/NewTS/
##

module purge
module load  mkl/2017.4
module load  gcc/7.1.0    
module load  openmpi/1.10.7    

make -f Makefile.plepp ple 
make -f Makefile.plepp plepp


## 2017SEP18 @ marenostrum 4 
## OK 
## 1) intel/2017.4   2) impi/2017.4   3) mkl/2017.4   4) bsc/1.0   5) python/2.7.13
##
## gcc/4.9.4 mkl/2017.4 openmpi/1.10.7|openmpi/1.10.4  
##
## NO 
## mkl/2017.4 gcc/7.1.0 openmpi/1.10.7|openmpi/1.10.4     
##
## gcc/5.4.0  mkl/2017.4 openmpi/1.10.7|openmpi/1.10.4  
##
## gcc/7.1.0  mkl/2017.4 openmpi/1.10.7|openmpi/1.10.4  
##/usr/bin/ld: warning: libgfortran.so.3, needed by /usr/mpi/gcc/openmpi-1.10.4-hfi/lib64/libmpi_usempi.so, may conflict with libgfortran.so.4
## 
 
