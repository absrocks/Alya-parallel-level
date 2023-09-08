
#FC=mpif90
#CXX=mpicxx
#ROOT_BOOST=/apps/BOOST/1.64.0/INTEL/IMPI
#ROOT_CANTERA=/home/bsc21/bsc21704/z2018_1/REPOSITORY/CANTERA/Cantera230/ExecsINTEL

INCLUDES =-I$(ROOT_CANTERA)/include
INCLUDES+=-I$(ROOT_BOOST)
INCLUDES+=-std=c++0x 

STATIC=ar rcs

$(eval export PYTHONPATH=$(PYTHONPATH):$(ROOT_CANTERA)/lib/python2.7/site-packages)


canteratest :  
	@$(FC) -c $@.f90 -I$(ROOT_CANTERA)/include/cantera 
	@$(FC) -o $@ $@.o $(ROOT_CANTERA)/lib/libcantera_fortran.a $(ROOT_CANTERA)/lib/libcantera.a -lstdc++
	@./$@  > outputCanteratest.txt    
	@echo CANTERA TEST !!  

static :  
	@$(CXX) -c cantera_wrapper.cpp $(INCLUDES) 
	@$(STATIC) libcantera4alya.a cantera_wrapper.o  


alyatest : static 
	@$(FC) -o $@ $@.f90 libcantera4alya.a -I$(ROOT_CANTERA)/include/cantera $(ROOT_CANTERA)/lib/libcantera_fortran.a $(ROOT_CANTERA)/lib/libcantera.a -lstdc++  
	@./$@ > outputAlyatest.txt  
	@echo ALYA TEST  !!  

clean :
	@rm *.o *.a alyatest  canteratest  output.txt   
