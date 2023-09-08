CC=icc
FC=ifort
CFLAGS= -g 
CCCFLAGS = -std=c99
LDFLAGS = -lcrypto -lssl -lintlc -I../openssl-1.1.0e/include -L../openssl-1.1.0e -Wl,-rpath=$(PWD)/../openssl-1.1.0e

#LDFLAGS = -lcrypto -lssl -lgfortran -I../openssl-1.1.0e/include -L../openssl-1.1.0e -Wl,-rpath=$(PWD)/../openssl-1.1.0e

objects = hello.o cwrapper.o chiper.o

TARGET = hello

all:$(TARGET)

hello:  cwrapper.o chiper.o test.o
	$(FC) $(CFLAGS) $(CINCS) $+ $(LDFLAGS) -o hello

lib: $(objects)
	ar cr libburgogne.a $(objects) 

hello.o: hello.f90
	$(FC) $(CFLAGS) $(CINCS) $+ $(LDFLAGS) -c #-S -fcode-asm

test.o: test.f90
	$(FC) $(CFLAGS) $(CINCS) $+ $(LDFLAGS) -c #-S -fcode-asm

cwrapper.o: cwrapper.c
	$(CC) $(CFLAGS) $(CCCFLAGS) $(CINCS) $+ $(LDFLAGS) -c

chiper.o: chiper.c
	$(CC) $(CFLAGS) $(CCCFLAGS) $(CINCS) $+ $(LDFLAGS) -c	

generator: generateLicense.c chiper.o
	$(CC) $(CFLAGS) $(CCCFLAGS) $(CINCS) $+ $(LDFLAGS) -o generateLicense

test: hello generator
	./test.plx

clean:
	rm -r *.o hello generateLicense libburgogne.a






