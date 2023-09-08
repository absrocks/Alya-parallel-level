FCOMP = gfortran
FCFLAGS = -c -J./ -I./ -ffree-line-length-none -x f95-cpp-input -O3
LDFLAGS =
PROGRAM = test_boast.x
SRCS = def_kintyp.f90 def_parameters.f90 runend.f90 nsi_element_assembly_split_oss_default.f90 mod_nsi_assembly.f90 test_boast.f90
OBJECTS = $(SRCS:.f90=.o)

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(FCOMP) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(FCOMP) $(FCFLAGS) ./$<

.PHONY: clean veryclean

clean:
	rm -f *.o *.mod *.MOD

veryclean: clean
	rm -f *~ $(PROGRAM)
