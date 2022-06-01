ifeq ($(strip $(NINJA_CC)),)
    CC = mpicc
else
    CC = $(NINJA_CC)
endif
ifeq ($(strip $(NINJA_CPP)),)
    CPP = mpic++
else
    CPP = $(NINJA_CPP)
endif
ifeq ($(strip $(NINJA_FC)),)
    FC = mpif90
else
    FC = $(NINJA_FC)
endif
ifeq ($(strip $(NINJA_NVCC)),)
    NVCC = nvcc
else
    NVCC = $(NINJA_NVCC)
endif

# What optimization level to use                                                                                                                               
OPTFLAGS = -O3 --restrict -gencode arch=compute_70,code=compute_70

# What archiving to use                                                                                                                                        
AR = ar rv

# What to use for indexing the archive                                                                                                                         
RANLIB = ranlib

INCLUDES = -fPIC -I../../Executables/${ALDIR}/Objects_x
CUDAINCLUDES = -I${CUDA_HOME}/include/ -Xcompiler=-fPIC
LD = -lcublas_static -lcusparse_static -lculibos -lcudart -lpthread -lstdc++ -ldl -lcusolver
LD_LIB = -L${CUDA_HOME}lib64/
CFLAGS = $(OPTFLAGS) $(INCLUDES) 

#lightspmv = LightSpMV-1.0/objs/*.o

#OBJS = diapregen.o gpuvar.o gpumemmanager.o gpublas.o gpusparse.o precon.o spmv.o codecg.o codecgdef.o codecgpipe.o codecgpipedef.o codegmres.o codeqr.o codeamgx.o  amgxcomm.o cudaapi.o cudakernel.o cudakernelapi.o cudasolverapi.o cudacublasapi.o cudacusparseapi.o cudaamgxapi.o
#in PIZ DAINT use
OBJS = diapregen.o gpuvar.o gpumemmanager.o gpublas.o gpusparse.o precon.o spmv.o codecg.o codecgnoprec.o codecgdef.o codecgpipe.o codegmres.o cudaapi.o cudakernel.o cudakernelapi.o cudacublasapi.o cudacusparseapi.o lasterrorninja.o

libninja.a: $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@

.c.o:
	$(CC) $(CFLAGS) -c $*.c

clean:
	rm -f *.o *~ *.mod

realclean:
	make clean
	rm -f libninja.a

diapregen.o:diapregen.f90
	$(FC) -c -o $@ $< $(INCLUDES) 

%.o:%.f90
	$(FC) -c -o $@ $< $(INCLUDES)

%.o:%.cu
	$(NVCC) -c -o $@ $< $(OPTFLAGS) $(CUDAINCLUDES)

%.o:%.c
	$(CC) -c -o $@ $< $(INCLUDES)

%.o:%.cc
	$(CPP) -c -o $@ $< $(INCLUDES)
