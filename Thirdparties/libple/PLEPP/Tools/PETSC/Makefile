ALL: exec 
CFLAGS   = -I./Blitz010/ -I./IO/
CPPFLAGS = -I./Blitz010/ -I./IO/
#FFLAGS =
#FPPFLAGS =
CLEANFILES = *.dat *.mtx *.vtk ./Bin/exec 


#PETSC_DIR=../Pectsc
PETSC_DIR= ../../../Developet23Oct/Cpp/Petsc3.3/
include ${PETSC_DIR}/conf/variables
include ${PETSC_DIR}/conf/rules


#sourse files. *.c* 
VPATH  = ./Src
VPATH += ./IO/


exec: main.o mmio.o  
	g++ -o ./Execs/$@ main.o ./Execs/mmio.o ${PETSC_LIB} 
	${RM} main.o

mmio.o: mmio.cpp
	g++ -c $^ -o ./Execs/$@ $(LIBS) $(INCLUDE)  
