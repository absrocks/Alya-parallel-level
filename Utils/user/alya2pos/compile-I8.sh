ifort -c -traceback   -O3 -fpp -DI8 def_kintyp.f90   -o def_kintyp.o
ifort -c -traceback   -O3  mod_maths.f90    -o mod_maths.o
ifort -c -traceback   -O3  def_elmtyp.f90   -o def_elmtyp.o
ifort -c -traceback   -O3  def_inpout.f90   -o def_inpout.o
ifort -c -traceback   -O3  mod_output.f90    -o mod_output.o
ifort -c -traceback   -O3  gidres_he.f90    -o gidres_he.o
ifort -c -traceback   -O3  gidres_gp.f90    -o gidres_gp.o
ifort -c -traceback   -O3  elmtyp.f90       -o elmtyp.o
ifort -c -traceback   -O3  ecoute.f90       -o ecoute.o
ifort -c -traceback   -O3  runend.f90       -o runend.o
ifort -c -traceback   -O3  connpo.f90       -o connpo.o
ifort -c -traceback   -O3  vu_msh.f90       -o vu_msh.o
ifort -c -traceback   -O3  vu_res.f90       -o vu_res.o
ifort -c -traceback   -O3  vu_filter.f90    -o vu_filter.o
ifort -c -traceback   -O3  ensmsh.f90       -o ensmsh.o
ifort -c -traceback   -O3  ensmsh_bin.f90       -o ensmsh_bin.o
ifort -c -traceback   -O3  ensres_bin.f90       -o ensres_bin.o
ifort -c -traceback   -O3  ensres_filter.f90       -o ensres_filter.o
ifort -c -traceback   -O3  ensmsh_filter.f90       -o ensmsh_filter.o
ifort -c -traceback   -O3  ensres.f90       -o ensres.o
ifort -c -traceback   -O3  txtres.f90       -o txtres.o
ifort -c -traceback   -O3  alya2pos.f90     -o alya2pos.o
ifort -c -traceback  -O3   wristl.f90       -o wristl.o
ifort -c -traceback  -O3   ensmsh_filter.f90  -o ensmsh_filter.o
ifort -c -traceback  -O3  reahed.f90 -o reahed.o
ifort -c -traceback  -O3  zfemres.f90 -o zfemres.o
ifort -traceback -O3 -o alya2pos.x *.o
rm -rf *.o rm *_genmod.f90 *.mod
