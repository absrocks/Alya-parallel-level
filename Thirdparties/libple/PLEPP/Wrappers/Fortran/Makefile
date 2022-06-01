include ../../Makefile.in 
#
# incluir en el configure de alya 
# print MAKEFILE "commdom:\n\t\@\$(MAKE) -C /home/jmake/z2014/Cplng/Mpi4python/All_Beta/Wrappers/Fortran \n\n";
#


SHARED = $(FCOMPILER) -shared
STATIC = ar rcs 

## make libcommdom PCOMMDOM=1

FLAGS =

PCOMMDOM=0 # PRINT COMMDOM
ifneq ($(PCOMMDOM), 0)
FLAGS += -D__PRINT__=$(PCOMMDOM)
endif

SATURNE=0
ifneq ($(SATURNE), 0)
FLAGS += -DSATURNE
endif

#LIBS_SATURNE

all:
	@echo
	#====================================================# 
	@echo "--> make [test01/test02/libcommdom/shared]"
	#====================================================#
	@echo 


libcommdom: 
	$(CPPCOMPILER) -c -g -fPIC $(INCLUDES)/build_octree.c $(INCLUDE_PLE) 
	$(CPPCOMPILER) -c -g -fPIC $(INCLUDES)/commdom.cpp $(INCLUDE_PLE) $(FLAGS)
	$(CPPCOMPILER) -c -g -fPIC $(WRAPPERS)/commdom_wrapper.cpp $(INCLUDE_PLE) 


static: 
	$(CPPCOMPILER) -c $(INCLUDES)/build_octree.c $(INCLUDE_PLE) 
	$(CPPCOMPILER) -c $(INCLUDES)/commdom.cpp $(INCLUDE_PLE) $(FLAGS) 
	$(CPPCOMPILER) -c $(WRAPPERS)/commdom_wrapper.cpp $(INCLUDE_PLE) $(FLAGS) 
	$(STATIC) libcommdom.a build_octree.o commdom.o commdom_wrapper.o $(ROOT_PLE)/lib/libple.a


shared: libcommdom
	$(SHARED) -o libcommdom.so build_octree.o commdom.o commdom_wrapper.o $(LIBS_PLE) -lstdc++  $(LIBS_SATURNE) -lmpi -lmpi_cxx 


test01: shared  
	$(FCOMPILER)   -c $(TESTS)/Fortran/main_loader.f90 
	$(FCOMPILER)   -c $(TESTS)/Fortran/main_loader02.f90
	$(CPPCOMPILER) -c $(SRC)/read_file.cpp $(INCLUDE_PLE)
	$(CPPCOMPILER) -c $(INCLUDES)/loader_alya.cpp $(INCLUDE_PLE) 
	$(CPPCOMPILER) -c $(WRAPPERS)/loader_alya_wrapper.cpp $(INCLUDE_PLE)
	$(CPPCOMPILER) -c $(WRAPPERS)/loader_alya_fortran.c   $(INCLUDE_PLE) 
	$(FCOMPILER)   -o wloader.x main_loader.o loader_alya_fortran.o loader_alya_wrapper.o loader_alya.o read_file.o $(LIBS_PLE) -lstdc++ -L. -lcommdom -lmpi $(LIBS_SATURNE) #-lmpi_cxx
	$(FCOMPILER)   -o wloader02.x main_loader02.o loader_alya_fortran.o loader_alya_wrapper.o loader_alya.o read_file.o $(LIBS_PLE) -lstdc++ -L. -lcommdom -lmpi $(LIBS_SATURNE) #-lmpi_cxx 
	#====================================================#
	@mpirun -np 2 ./wloader.x CCCC : -np 2 ./wloader.x BBBB : -np 3 ./wloader.x AAAA
	#====================================================#
	@mpirun -np 2 ./wloader.x CCCC : -np 2 ./wloader02.x DDDD : -np 1 ./wloader.x COSA --name AAAA 
	#====================================================#


test02: static  
	$(FCOMPILER)   -c $(TESTS)/Fortran/main_loader.f90 
	$(FCOMPILER)   -c $(TESTS)/Fortran/main_loader02.f90
	$(CPPCOMPILER) -c $(SRC)/read_file.cpp $(INCLUDE_PLE)
	$(CPPCOMPILER) -c $(INCLUDES)/loader_alya.cpp $(INCLUDE_PLE) 
	$(CPPCOMPILER) -c $(WRAPPERS)/loader_alya_wrapper.cpp $(INCLUDE_PLE)
	$(CPPCOMPILER) -c $(WRAPPERS)/loader_alya_fortran.c   $(INCLUDE_PLE) 
	$(FCOMPILER)   -o wloader.x main_loader.o loader_alya_fortran.o loader_alya_wrapper.o loader_alya.o read_file.o -lstdc++ -L. -lcommdom -lmpi $(ROOT_PLE)/lib/libple.a $(LIBS_SATURNE) -lmpi_cxx
	$(FCOMPILER)   -o wloader02.x main_loader02.o loader_alya_fortran.o loader_alya_wrapper.o loader_alya.o read_file.o -lstdc++ -L. -lcommdom -lmpi $(ROOT_PLE)/lib/libple.a $(LIBS_SATURNE) -lmpi_cxx 
	#====================================================#
	@mpirun -np 2 ./wloader.x CCCC : -np 2 ./wloader.x BBBB : -np 3 ./wloader.x AAAA
	#====================================================#

clean:
	rm -rf *~ *.o *.x *.mod *.so *.a
 
