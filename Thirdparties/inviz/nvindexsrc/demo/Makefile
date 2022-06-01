#*****************************************************************************
# Copyright 1986, 2016 NVIDIA Corporation. All rights reserved.
#*****************************************************************************
HAS_CUDA_SRC := yes
NVCC := nvcc -ccbin=gcc-4.9

CXX      = mpic++
CC       = mpicc

CPPFLAGS ?= \
        -DBIT64 -D_GNU_SOURCE -D_THREAD_SAFE -D_REENTRANT \
        -D__USE_XOPEN2K8 -pthread -DLINUX     \
        -march=nocona -m64 -fPIC -fno-strict-aliasing     \
         -I${NVINDEX_ROOT}/include -I.. -I. -I../../ -I${CUDA_HOME}/include                

LDFLAGS  ?=  
LIBS     := -ldl -lz -lm -lHalf -llog4cplus -lrt   -pthread


ifdef LLIBS
    LIBS += $(LLIBS)
endif	

EXEC_PATH := ../objs
COMMONLIB := $(EXEC_PATH)/libcommon.a

LIBS      += -lcudart

# Build viewer as MPI application
CPPFLAGS  += -DNVINDEX_HAS_MPI_IPC_COMPUTE

# optimize build options
CPPFLAGS  += -O3 -DNDEBUG
# debug build options
# CPPFLAGS += -g -DDEBUG

NVCC_FLAGS := -m64 -use_fast_math -O3 -gencode arch=compute_35,code=compute_35

# Defines $(SRC)
include sources.mk

# Add source files for OpenGL version
SRC += $(SRC_NO_GL) $(SRC_MPI_IPC_COMPUTE)

CUDASRC := $(SRC_MPI_IPC_COMPUTE_CUDA)

OBJS := $(SRC:.cpp=.o)
OBJS := $(OBJS:.c=.o)
OBJS := $(addprefix $(EXEC_PATH)/,$(OBJS))
CUDAOBJ := $(CUDASRC:.cu=.o)
CUDAOBJ := $(addprefix $(EXEC_PATH)/,$(CUDAOBJ))

EXEC := libviewer.a
all:   $(EXEC_PATH)/$(EXEC)

$(EXEC_PATH)/%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) -o $@ $<

$(EXEC_PATH)/%.o: %.c
	$(CC)  -c $(CPPFLAGS) -o $@ $<

$(EXEC_PATH)/%.o: %.cu $(HEADERS)
	$(NVCC) -c $(NVCC_FLAGS) --compiler-options  " $(CPPFLAGS) "  -o $@ $<

$(EXEC_PATH)/$(EXEC) : $(OBJS) $(CUDAOBJ) $(COMMONLIB) $(SHAREDLIB) 



$(COMMONLIB):
	$(MAKE) -C ../common

$(EXEC_PATH):
	mkdir -p $@
Makefile: $(EXEC_PATH)

clean:
	rm -rf $(OBJS) $(EXEC_PATH)/$(EXEC) $(CUDAOBJ)

# Local Variables: ***
# mode: Makefile ***
# End: ***
