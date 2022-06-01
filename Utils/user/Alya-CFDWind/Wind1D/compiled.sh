ifort -c -i4 -r8 -fpp  -g  def_master.f90   -o def_master.o
ifort -c -i4 -r8 -fpp  -g  Turnon.f90       -o Turnon.o
ifort -c -i4 -r8 -fpp  -g   Doiter.f90      -o Doiter.o
ifort -c -i4 -r8 -fpp  -g   Reapro.f90      -o Reapro.o
ifort -c -i4 -r8 -fpp  -g   Solite.f90      -o Solite.o
ifort -c -i4 -r8 -fpp  -g   Solvtr.f90      -o Solvtr.o
ifort -c -i4 -r8 -fpp  -g   Turnof.f90      -o Turnof.o
ifort -c -i4 -r8 -fpp  -g    Wind.f90       -o Wind.o
ifort -c -i4 -r8 -fpp  -g    solvpa.f       -o solvpa.o
ifort  -O0 -traceback -g -debug all -nogen-interfaces -ftrapuv -lmkl_core -lmkl_sequential -lmkl_intel_lp64	-lpthread  -o  Wind.g *.o 
rm -rf *.o rm *_genmod.f90 *.mod
