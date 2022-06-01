ifort -c -debug -traceback -fpp def_kintyp.f90	
ifort -c -debug -traceback main.f90	
ifort -c -debug -traceback reagrid.f90
ifort -debug -traceback *.o
