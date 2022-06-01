mpif90 -c -cpp   -O3  def_kintyp.f90   -o def_kintyp.o
mpif90 -c -cpp   -O3  mod_maths.f90    -o mod_maths.o
mpif90 -c -cpp   -O3  mod_output.f90    -o mod_output.o
mpif90 -c -cpp   -O3  def_elmtyp.f90   -o def_elmtyp.o
mpif90 -c -cpp   -O3  def_inpout.f90   -o def_inpout.o
mpif90 -c -cpp   -O3  gidres_he.f90    -o gidres_he.o
mpif90 -c -cpp   -O3  gidres_gp.f90    -o gidres_gp.o
mpif90 -c -cpp   -O3  elmtyp.f90       -o elmtyp.o
mpif90 -c -cpp   -O3  ecoute.f90       -o ecoute.o
mpif90 -c -cpp   -O3  runend.f90       -o runend.o
mpif90 -c -cpp   -O3  connpo.f90       -o connpo.o
mpif90 -c -cpp   -O3  vu_msh.f90       -o vu_msh.o
mpif90 -c -cpp   -O3  vu_res.f90       -o vu_res.o
mpif90 -c -cpp   -O3  vu_filter.f90    -o vu_filter.o
mpif90 -c -cpp   -O3  ensmsh.f90       -o ensmsh.o
mpif90 -c -cpp   -O3  ensmsh_bin.f90       -o ensmsh_bin.o
mpif90 -c -cpp   -O3  ensres_bin.f90       -o ensres_bin.o
mpif90 -c -cpp   -O3  ensres_filter.f90       -o ensres_filter.o
mpif90 -c -cpp   -O3  ensmsh_filter.f90       -o ensmsh_filter.o
mpif90 -c -cpp   -O3  ensres.f90       -o ensres.o
mpif90 -c -cpp   -O3  txtres.f90       -o txtres.o
mpif90 -c -cpp   -O3  alya2pos.f90     -o alya2pos.o
mpif90 -c -cpp  -O3   wristl.f90       -o wristl.o
mpif90 -c -cpp  -O3   ensmsh_filter.f90  -o ensmsh_filter.o
mpif90 -c -cpp  -O3  reahed.f90 -o reahed.o
mpif90 -c -cpp  -O3  zfemres.f90 -o zfemres.o
mpif90 -cpp -O3 -o alya2pos.x *.o
rm -rf *.o rm *_genmod.f90 *.mod
