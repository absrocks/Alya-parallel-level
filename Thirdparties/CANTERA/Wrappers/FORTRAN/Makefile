INCLUDES =-I$(ROOT_CANTERA)/include
INCLUDES+=-I$(ROOT_BOOST)
INCLUDES+=-std=c++0x 

STATIC=ar rcs

static :  
	@$(CXX) -c cantera_wrapper.cpp $(INCLUDES) 
	@$(STATIC) libcantera4alya.a cantera_wrapper.o  

clean :
	@rm *.o *.a alyatest  canteratest  output.txt   
