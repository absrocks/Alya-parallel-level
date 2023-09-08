# Which compiler to use                                                                                                                                        
 FC = mpif90
 NVCC = nvcc
 CC = mpicc
 CPP = mpic++
# What optimization level to use                                                                                                                               
OPTFLAGS = -O3 --restrict -gencode arch=compute_35,code=compute_35

# What archiving to use                                                                                                                                        
AR = ar rv

# What to use for indexing the archive                                                                                                                         
RANLIB = ranlib

INCLUDES = -I../../Executables/unix/Objects_x/ -fPIC
CUDAINCLUDES = -I${CUDA_HOME}/include/ -Xcompiler=-fPIC
LD = -lcublas_static -lcusparse_static -lculibos -lcudart -lpthread -lstdc++ -ldl -lcusolver
LD_LIB = -L${CUDA_HOME}lib64/
CFLAGS = $(OPTFLAGS) $(INCLUDES) 

#lightspmv = LightSpMV-1.0/objs/*.o

OBJS = launcher.o passcomm.o

all: libinviz.a

libinviz.a: $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@

.c.o:
	$(CC) $(CFLAGS) -c $*.c

clean:
	rm -f *.o *~ *.mod
	rm -rf nvindexsrc/objs/

realclean:
	make clean
	rm -f libinviz.a libviewer.a

%.o:%.f90
	$(FC) -c -o $@ $< $(INCLUDES)

%.o:%.cu
	$(NVCC) -c -o $@ $< $(OPTFLAGS) $(CUDAINCLUDES)

%.o:%.c
	$(CC) -c -o $@ $< $(INCLUDES)

%.o:%.cpp
	$(CPP) -c -o $@ $< $(INCLUDES)
