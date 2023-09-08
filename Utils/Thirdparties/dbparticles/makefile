#CC=gcc -O3
#CC=icc -O3 -axMIC-AVX512,xCORE-AVX2

#-----------CFLAGS---------------------------
CFLAGS=-std=c++11 -I$(HECUBA_ROOT)/include -L$(HECUBA_ROOT)/lib -I$(HFETCH_ROOT) -I/apps/INTEL/2017.4/impi/2017.3.196/include64
NOLINK = -c

#------------LIBS------------------
#For cassandra interface
LIBS = -lhfetch -lcassandra


AR = ar rv
# What to use for indexing the archive                                                                                                                         
RANLIB = ranlib

all: libdbparticles.so 

static: libdbparticles.a


# Read binary output and write to Cassandra
AlyaMock.x: ParticlesToCassandraAPI.h BinaryToCassandra.cpp libdbparticles.so
	$(CXX) $(CFLAGS) $(NOLINK) BinaryToCassandra.cpp
	$(CXX) $(CFLAGS) -L. -ldbparticles $(LIBS) BinaryToCassandra.o  -o AlyaMock.x 



# Output to Cassandra, object for the static lib
ParticlesToCassandraConnector.o: ParticlesToCassandraAPI.h ParticlesToCassandraConnector.cpp ParticlesToCassandraConnector.h
	$(CXX) $(CFLAGS) $(NOLINK) ParticlesToCassandraConnector.cpp


# Output to Cassandra, object for the dynamic lib
ParticlesToCassandraConnector_sh.o: ParticlesToCassandraAPI.h ParticlesToCassandraConnector.cpp ParticlesToCassandraConnector.h
	$(CXX) $(CFLAGS) $(NOLINK) -fPIC ParticlesToCassandraConnector.cpp -o ParticlesToCassandraConnector_sh.o


# Output to binary files
BinaryConnector.o: ParticlesToCassandraAPI.h ConnectorBinaryFile.cpp ConnectorBinaryFile.h 
	$(CXX) $(CFLAGS) $(NOLINK) -fPIC ConnectorBinaryFile.cpp -o BinaryConnector.o



# Static DBParticles
libdbparticles.a: ParticlesToCassandraConnector.o
	$(AR) $@ ParticlesToCassandraConnector.o
	$(RANLIB) $@


# Dynamic DBParticles for binary files
libdbparticlesbin.so: BinaryConnector.o
	$(CXX) $(CFLAGS) -shared -o libdbparticlesbin.so BinaryConnector.o


# Dynamic DBParticles writing to Cassandra
libdbparticles.so: ParticlesToCassandraConnector_sh.o
	$(CXX) $(CFLAGS) $(LIBS) -shared -o libdbparticles.so ParticlesToCassandraConnector_sh.o

ConnectorFlushAfterWrite.o: ParticlesToCassandraAPI.h ConnectorFlushAfterWrite.h ConnectorFlushAfterWrite.cpp
	$(CXX) $(CFLAGS) $(NOLINK) -fPIC ConnectorFlushAfterWrite.cpp -o ConnectorFlushAfterWrite.o

libdbparticlesflush.so: ConnectorFlushAfterWrite.o
	$(CXX) $(CFLAGS) $(LIBS) -shared -o libdbparticlesflush.so ConnectorFlushAfterWrite.o


clean:
	rm -rf *.o libdbparticles* AlyaMock.x
	
