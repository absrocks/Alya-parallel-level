#
# 2015Jun24 created  
# 
include ../../Makefile.in

SHARED = $(FCOMPILER) -shared
STATIC = ar     cuda rcs
FLAGS =

ROOT_CUDA    =  /Developer/NVIDIA/CUDA-6.5/
LIB_CUDA     = -L$(ROOT_CUDA)/lib      -lcudart         
FLAGS_CUDA   =  -Xlinker -rpath -Xlinker $(ROOT_CUDA)/lib -stdlib=libstdc++ 
CUDACOMPILER =  $(ROOT_CUDA)/bin/nvcc 



INCLUDE_MPI  = -I/Users/aerozita/z2015/Openmpi185/include
ROOT_MPI    := /opt/local/
INCLUDE_MPI := -I$(ROOT_MPI)/include/mpich-mp/

###
## ERROR: Undefined symbols for architecture x86_64 
##       -arch x86_64  -stdlib=libstdc++ 
##
## Minotauro
## @module unload cuda/4.1 && module load cuda/5.0
##

all:
	@echo
	#====================================================#
	@echo "--> make [static|cuda|run]"
	#====================================================#
	@echo


cuda:   static  
	$(CUDACOMPILER) -c           cuda_mpi.cu $(INCLUDE_PLE) $(INCLUDE_MPI)
	$(CPPCOMPILER)  -o wloader.x cuda_mpi.o  $(FLAGS_CUDA)  $(LIB_CUDA)    $(ROOT_PLE)/lib/libple.a ./libcommdom.a
	#$(CPPCOMPILER)  -o wloader.x cuda_mpi.o -Xlinker -rpath -Xlinker /Developer/NVIDIA/CUDA-6.5/lib  -stdlib=libstdc++ -L/Developer/NVIDIA/CUDA-6.5/lib -lcudart $(ROOT_PLE)/lib/libple.a ./libcommdom.a 
	#====================================================#
	#@mpirun-mpich-mp  -np 2 ./wloader.x FLUID : -np 2 ./wloader.x SOLID
	#====================================================#


static:
	$(CCOMPILER)    -c $(INCLUDES)/build_octree.c      $(INCLUDE_PLE)          -arch x86_64  
	$(CPPCOMPILER)  -c $(INCLUDES)/commdom.cpp         $(INCLUDE_PLE) $(FLAGS) -arch x86_64  -stdlib=libstdc++
	$(CPPCOMPILER)  -c $(WRAPPERS)/commdom_wrapper.cpp $(INCLUDE_PLE)          -arch x86_64  -stdlib=libstdc++ 
	$(STATIC) libcommdom.a build_octree.o commdom.o $(ROOT_PLE)/lib/libple.a 
	#$(CPPCOMPILER)  -o wloader.x cuda_mpi.o build_octree.o commdom.o commdom_wrapper.o -Xlinker -rpath -Xlinker /Developer/NVIDIA/CUDA-6.5/lib  -stdlib=libstdc++ -L/Developer/NVIDIA/CUDA-6.5/lib -lcudart $(ROOT_PLE)/lib/libple.a 


run:    cuda 
	#====================================================#
	#@mpirun-mpich-mp  -np 2 ./wloader.x "AAA" : -np 3 ./wloader.x "BBB" :  -np 4 ./wloader.x "CCC"
	#@mpirun-mpich-mp  -np 2 ./wloader.x "CUDA" : -np 1 ../Cpp//wloader.x CPP  
	#@mpirun-mpich-mp  -np 2 ./wloader.x "CUDA" : -np 1 ../Fortran//wloader.x Fortran  
	@mpirun-mpich-mp  -np 2 ./wloader.x FLUID : -np 1 /Users/aerozita/z2015/REPOSITORY/ALYA_2015JUN04/Executables/unix/Alya.x the_case --name SOLID  
	#====================================================#


clean:
	rm -rf *~ *.o *.x *.mod *.so *.a


