#*****************************************************************************
# Copyright 1986, 2016 NVIDIA Corporation. All rights reserved.
#*****************************************************************************

CXX      = g++
CC       = gcc
CPPFLAGS ?= \
            -DBIT64 -D_GNU_SOURCE -D_THREAD_SAFE -D_REENTRANT \
	    -D__USE_XOPEN2K8 -pthread -DLINUX                 \
            -march=nocona -m64 -fPIC -fno-strict-aliasing     \
            -I. -std=c++11 -I${NVINDEX_ROOT}/include/

EXEC_PATH = ../objs

# optimize build options
CPPFLAGS += -O3 -DNDEBUG
# debug build options
# CPPFLAGS += -g -DDEBUG

include sources.mk

OBJS   = $(SRC:%.cpp=$(EXEC_PATH)/%.o)

all: $(OBJS)

$(EXEC_PATH)/%.o: %.cpp
	$(CXX) -c $(CPPFLAGS) -o $@ $<

$(EXEC_PATH):
	mkdir -p $@
Makefile: $(EXEC_PATH)

clean:
	rm -rf $(OBJS) $(EXEC_PATH)/$(COMMON_LIB) $(COMMON_LIB)
	rm -f *~

# Local Variables: ***
# mode: Makefile ***
# End: ***
