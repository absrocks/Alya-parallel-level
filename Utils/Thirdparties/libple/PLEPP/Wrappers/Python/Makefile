include ../../Makefile.in

all:
	@echo
	#====================================================# 
	@echo "--> make [Commdomm/Read_file/Loader_alya]"
	@echo
	@echo "--> In order to use PLE++ and Paraview together:" 
	@echo "      CXX = mpic++.mpich2"
	@echo "      ./configure CC=mpicc.mpich2 --prefix=PREFIX_PATH"
	@echo "      export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:PREFIX_PATH"
	@echo "      mpirun.mpich2 -np 1 pvpython XXX.x "
	@echo 
	@echo "--> To use mpi4py:"
	@echo "      from mpi4py import MPI"
	@echo "      sys.path.append('ParaView_PATH/lib/paraview-4.1/site-packages/')"
	#====================================================#
	@echo 

## make libcommdom PCOMMDOM=1
## PCOMMDOM = 0 # PRINT COMMDOM
## CXX = mpic++.mpich2

Commdomm: swig_commdom 
	$(CXX) -c -fPIC $(INCLUDES)/build_octree.c  $(INCLUDE_PLE)
	$(CXX) -c -fPIC $(INCLUDES)/commdom.cpp     $(INCLUDE_PLE) 
	$(CXX) -c -fPIC ./commdom_wrap.cxx          $(INCLUDE_PLE) $(LIBS_PLE) $(INCLUDE_PY)
	$(CXX) -shared  build_octree.o commdom.o commdom_wrap.o -o _Commdomm.so $(LIBS_PLE) $(LIBS_PY) #-lc -lstdc++  
	#====================================================#
	#@mpirun -np 1 python $(HOME)/Tests/Locator/Python/main_mesh02.py : -np 1 python $(HOME)/Tests/Locator/Python/main_mesh01.py
	#====================================================#
	@mpirun -np 3 python  $(TESTS)/Python/main.py CCCC : -np 2 python  $(TESTS)/Python/main.py AAAA
	#====================================================#
	#@mpirun -np 1 python  $(TESTS)/Python/main.py CCCC : -np 1 python  $(TESTS)/Python/main.py AAAA


swig_commdom: 
	#swig -python -c++ commdom.i
	#====================================================#
	@$(ROOT_SWIG)/swig -python -c++ $(INCLUDE_SWIG) commdom.i 
	#====================================================#


Read_file: swig_read_file
	$(CXX) -c -fPIC $(SRC)/read_file.cpp $(INCLUDE)
	$(CXX) -c -fPIC ./read_file_wrap.cxx $(INCLUDE_PY) 
	$(CXX) -shared  read_file.o read_file_wrap.o -o _Read_file.so #-lc -lstdc++
	#====================================================#
	#====================================================#

swig_read_file: 
	swig -python -c++ read_file.i 


Loader_alya: swig_loader_alya
	$(CXX) -c -fPIC $(SRC)/read_file.cpp $(INCLUDE) 
	$(CXX) -c -fPIC $(INCLUDES)/loader_alya.cpp $(INCLUDE)
	$(CXX) -c -fPIC ./loader_alya_wrap.cxx $(INCLUDE_PY) 
	$(CXX) -shared  loader_alya.o loader_alya_wrap.o read_file.o -o _Loader_alya.so -lc -lstdc++
	#====================================================#
	#====================================================#

swig_loader_alya:
	swig -python -c++ loader_alya.i


#syrthes: Commdomm
#	#====================================================#
#	@mpirun -np 1 python $(HOME)/Tests/Python/ 
#	#====================================================#	


clean:
	rm -rf *~ *.o *.x *.mod *.so *.a *.pyc *_wrap.cxx *.py  
#
# apt-get install python-dev swig python-mpi4py
#
