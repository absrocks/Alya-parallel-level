include ../../Makefile.in
#
# incluir en el configure de alya
# print MAKEFILE "commdom:\n\t\@\$(MAKE) -C /home/jmake/z2014/Cplng/Mpi4python/All_Beta/Wrappers/Fortran \n\n";
#


## make libcommdom PCOMMDOM=1

SHARED = $(FCOMPILER) -shared
STATIC = ar rcs

FLAGS =

PCOMMDOM=0 # PRINT COMMDOM
ifneq ($(PCOMMDOM), 0)
FLAGS += -D__PRINT__=$(PCOMMDOM)
endif

SATURNE=0
ifneq ($(SATURNE), 0)
FLAGS += -DSATURNE
endif

BOOST=0
ifneq ($(BOOST), 0)
FLAGS += -I$(BOOST)
endif

all:
	@echo
	#====================================================#
	@echo "--> make [test01/test02/libcommdom/shared]"
	#====================================================#
	@echo


static:
	$(CCOMPILER) -c $(INCLUDES)/build_octree.c $(INCLUDE_PLE)
	$(CPPCOMPILER) -c $(INCLUDES)/commdom.cpp $(INCLUDE_PLE) $(FLAGS) 
	$(CPPCOMPILER) -c $(WRAPPERS)/commdom_wrapper.cpp $(INCLUDE_PLE)
	$(STATIC) libcommdom.a build_octree.o commdom.o $(ROOT_PLE)/lib/libple.a 


test02: static
	$(CPPCOMPILER) -c ./main_loader.cpp $(INCLUDE_PLE)
	$(CPPCOMPILER) -o wloader.x main_loader.o ./libcommdom.a $(ROOT_PLE)/lib/libple.a $(LIBS_SATURNE)
	#====================================================#
	@mpirun -np 2 ./wloader.x CCCC : -np 2 ./wloader.x BBBB : -np 3 ./wloader.x AAAA
	#====================================================#


test04: static
	$(CPPCOMPILER) -c ../../Src/read_file.cpp  $(INCLUDE_PLE)
	$(CPPCOMPILER) -c ./test04.cpp     $(INCLUDE_PLE) $(FLAGS)  
	$(CPPCOMPILER) -o wloader.x test04.o read_file.o ./libcommdom.a $(ROOT_PLE)/lib/libple.a $(LIBS_SATURNE)
	#====================================================#
	#@mpirun -np 1 ./wloader.x DIRIC  
	#====================================================#

test05: static
	$(CPPCOMPILER) -c ../../Src/read_file.cpp  $(INCLUDE_PLE)
	$(CPPCOMPILER) -c ./$@.cpp     $(INCLUDE_PLE) $(FLAGS)  
	$(CPPCOMPILER) -o $@.x $@.o read_file.o ./libcommdom.a $(ROOT_PLE)/lib/libple.a $(LIBS_SATURNE) 
	#====================================================#
	@mpirun -np 1 ./$@.x 
	#====================================================#


plepp2.0: static
	$(CPPCOMPILER) -c ../Src/read_file.cpp  $(INCLUDE_PLE)
	$(CPPCOMPILER) -c ./multi_plepp.cpp     $(INCLUDE_PLE)
	$(CPPCOMPILER) -o wloader.x multi_plepp.o read_file.o ./libcommdom.a $(ROOT_PLE)/lib/libple.a $(LIBS_SATURNE)

	#====================================================#
	#
	#1) ple_locator.c
	#  comment "_ple_locator_t"
	#	
	#2) ple_locator.h  
	#  before "typedef struct _ple_locator_t ple_locator_t;" 
	#  add
	#  "
	#  struct _ple_locator_t {
	#  ...
	#  } 
	#  "
	#:)
	#  enjoy 
	#
	#@mpirun-mpich-mp -np 1 ./wloader.x Mesh03 : -np 1 ./wloader.x Mesh02 : -np 1 ./wloader.x Mesh01
	#
	#====================================================#

clean:
	rm -rf *~ *.o *.x *.mod *.so *.a


