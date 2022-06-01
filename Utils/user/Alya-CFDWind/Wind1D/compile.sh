ifort -c -i4 -r8   -fpp -O3  def_master.f90   -o def_master.o
ifort -c  -i4 -r8  -fpp  -O3  Turnon.f90       -o Turnon.o
ifort -c -i4 -r8   -fpp  -O3  Doiter.f90       -o Doiter.o
ifort -c -i4 -r8   -fpp -O3  Reapro.f90       -o Reapro.o
ifort -c -i4 -r8   -fpp  -O3  Solite.f90       -o Solite.o
ifort -c -i4 -r8   -fpp  -O3  Solvtr.f90       -o Solvtr.o
ifort -c -i4 -r8   -fpp -O3  Turnof.f90       -o Turnof.o
ifort -c -i4 -r8   -fpp  -O3  Wind.f90         -o Wind.o
ifort -c -i4 -r8   -fpp  -O3    solvpa.f       -o solvpa.o
ifort  -O3 -lmkl_core -lmkl_sequential -lmkl_intel_lp64	-lpthread -o Wind.x *.o 
rm -rf *.o rm *_genmod.f90 *.mod


