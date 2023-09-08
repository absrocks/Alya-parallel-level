#*****************************************************************************
# Copyright 1986, 2016 NVIDIA Corporation. All rights reserved.
#*****************************************************************************

#
# src/Makefile
#
# Makefile to build all NVIDIA IndeX demo and example programs.
#

CC  = mpicc
CXX = mpic++
libraries := common
programs := demo
alya := alya
staticlib := libviewer.a
subdirs := $(libraries) $(programs) $(alya)

.PHONY: all $(subdirs)

all:$(staticlib)

$(staticlib): allobjs
	ar rv $@ objs/*.o
	ranlib $@
	cp $(staticlib) ../
allobjs: $(subdirs)


$(libraries):
	$(MAKE) CC=$(CC) CXX=$(CXX) --directory=$@
$(programs): 
	$(MAKE) CC=$(CC) CXX=$(CXX) --directory=$@
$(alya):
	$(MAKE) CC=$(CC) CXX=$(CXX) --directory=$@	

clean:
	@for subdir in $(subdirs); do           \
	   $(MAKE) --directory=$$subdir clean;  \
	done

#----------------------------------------------------------------------
# Local Variables:
# mode: Makefile
# End:
#----------------------------------------------------------------------
