##
## http://pythonhosted.org/mpi4py/
## http://pythonhosted.org/mpi4py/mpi4py.pdf 
##
## 1) 
#    python setup.py build --mpicc=/apps/OPENMPI/1.8.1-mellanox/MULTITHREAD/bin/mpicc
## 2) 
#    python setup.py install --home=/home/bsc21/bsc21704/z2015/REPOSITORY/MPI4PY131/Execs --prefix=''  
## 3)  
##   mpirun -np 16 python helloworld.py  
##   import sys
##   sys.path.append("/home/bsc21/bsc21704/z2015/REPOSITORY/MPI4PY131/Execs/lib64/python/") 
## 
