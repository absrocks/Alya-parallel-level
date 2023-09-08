#
# /opt/cuda/5.0/bin/nvcc -m64  -gencode arch=compute_10,code=sm_10 -gencode arch=compute_20,code=sm_20 -gencode arch=compute_30,code=sm_30  -o simpleMPI.o -c simpleMPI.cu
#
# nvcc cuda-mpi4.cu -I$CUDA_SDK/shared/inc -I$MPI_HOME/include -L$MPI_HOME/lib64/ -lmpi -arch sm_13  -o cuda-mpi.exe
#
# /gpfs/projects/bsc21/bsc21704/RUNNER_2015/REPOSITORY/ALYA_2015May25/Thirdparties/libple/PLEPP/Wrappers/Cpp/
#
#
module unload cuda/4.1 && module load cuda/5.0


/opt/cuda/5.0/bin/nvcc cuda_mpi.cu  -o cuda_mpi.exe -lmpi -L/apps/OPENMPI/1.8.1/lib/ -lmpi -I/apps/OPENMPI/1.8.1/include/ /gpfs/projects/bsc21/bsc21704/RUNNER_2015/REPOSITORY/ALYA_2015May25/Thirdparties/libple/PLEPP/Wrappers/Cpp/libcommdom.a


mpirun -np 2 ./cuda_mpi.exe 
