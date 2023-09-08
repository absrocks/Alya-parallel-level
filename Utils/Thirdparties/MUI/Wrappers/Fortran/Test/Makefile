#!/bin/bash

CC	= mpic++
CFLAGS	= -std=c++11 -O3 -w -fpermissive -c  

SRCS=$(wildcard *.cpp) 
OBJS=$(SRCS:.cpp=.o )
EXEC=a.x

default: $(EXEC)
#all:$(OBJS)

%.o: %.cpp  
	@$(CC) $(CFLAGS) -c $< -o $@


$(EXEC): $(OBJS) 
	@$(CC) -o $@ $(LFLAGS) $(OBJS)  
	@#mpirun -np 1 gdb -ex=r --args ./$(EXEC) MESH02 : -np 1 gdb -ex=r --args ./$(EXEC) MESH01 
	mpirun -np 1 ./$(EXEC) MESH02 : -np 1 ./$(EXEC) MESH01 


clean:
	rm -f $(EXEC) *.txt *.o *.x 
