# Used by Eclipse 
ROOT       = $(shell pwd)
ALYA4CPLNG = $(ROOT)
LIBPLE     = $(ROOT)/../


all: subsystem

subsystem:
	@$(MAKE) -C $(ALYA4CPLNG)/Src/ -f Makefile.free
	@mv $(ALYA4CPLNG)/Src/libalya4cplng.a $(ALYA4CPLNG)/Lib 
	@echo OK!!
	@echo  

libple:
	@echo 
	@echo ./configure CC=mpicc --prefix=$(LIBPLE)/Execs
	@echo make 
	@echo make install 
	@echo 
clean:
	@$(MAKE) -C $(ALYA4CPLNG)/Src/ -f Makefile.free clean
	@rm $(ALYA4CPLNG)/Lib/libalya4cplng.a  
