
#CXX=mpicxx
#ROOT_BOOST=/apps/BOOST/1.64.0/INTEL/IMPI
#ROOT_CANTERA=/home/bsc21/bsc21704/z2018_1/REPOSITORY/CANTERA/Cantera230/ExecsINTEL

INCLUDES =-I$(ROOT_CANTERA)/include
INCLUDES+=-I$(ROOT_BOOST)
INCLUDES+=-std=c++0x 

$(eval export PYTHONPATH=$(PYTHONPATH):$(ROOT_CANTERA)/lib/python2.7/site-packages)


canteratest : 
	@$(CXX) -c $@.cpp $(INCLUDES) 
	@$(CXX) -o $@ $@.o $(ROOT_CANTERA)/lib/libcantera.a -pthread -O3         
	@./$@   > output.txt  
	@echo CPP TEST !!  

clean :
	@rm *.o canteratest output.txt   
