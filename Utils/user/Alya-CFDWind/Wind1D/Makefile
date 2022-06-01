#! /bin/sh -f
#
#  1. Define the compiler (gfortran, ifort, f90, xlf90)
#
FC=ifort
#  
#  2. Code vesion
#
MACHINE=MN4
ver=dev

EXE=Wind1d.x
#
#  3. Location of the NetCDF library (at MN)
#
NETCDF=YES
#
#needs the flag -fpp to recognise "#ifdef _NETCDF_" inside the code
#This allows Wind1D compilation without: read_meso.o & netCDF_output.o
#DEFINES=-D_NETCDF_   
ifeq ($(NETCDF),YES)
ifeq ($(MACHINE),MN4)
  # LIB_NetCDF= $(NETCDF_LIB)/libnetcdf.a 
  # INC_NetCDF= -I$(NETCDF_INC)
   LIB_NetCDF=-L/apps/NETCDF/4.4.1.1/INTEL/IMPI/lib -lnetcdff -lnetcdf  
   INC_NetCDF= -I$(NETCDF_INC)
endif
ifeq ($(MACHINE),N3)
   LIB_NetCDF= $(NETCDF_LIB)/libnetcdff.a $(NETCDF_LIB)/libnetcdf.a
   INC_NetCDF= -I$(NETCDF_INC)
endif 
#    LIB_NetCDF= /apps/NETCDF/3.6.2/lib/libnetcdf.a
#    INC_NetCDF= -I/apps/NETCDF/3.6.2/include/
    DEFINES=-D_NETCDF_   
endif

#
#----------------------------------------
# Do not change anything beyond this line
#----------------------------------------
#
# Compiler selection
#
ifeq ($(FC),gfortran)
#  FFLAGS= -O1 -ffree-line-length-none  -fdefault-real-8
  FFLAGS= -O0 -ffree-line-length-none -fdefault-real-8  -fbounds-check -Wall -fbacktrace -ftrapv
  LINKER= $(FC)
  LFLAGS=
endif	
ifeq ($(FC),f90)
  FFLAGS= -extend_source -r8 -64
  LINKER= $(FC) -64
  LFLAGS=
endif	
ifeq ($(FC),ifort)
  FFLAGS=  -r8 -i4 -fpp -O3 -traceback  
#  FFLAGS= -132 -r8 -i4 -fpp -O0  -g -check all -traceback -debug full -warn all -init=snan   -nogen-interfaces -ftrapuv -fp-stack-check              
#   FFLAGS= -132 -r8 -i4 -fpp -g  -traceback -ftrapuv # -init=snan                 
  LINKER= $(FC)
  LFLAGS= -O1 -traceback 
#  PFLAGS= -profile-loops=all -profile-loops-report=2
endif
ifeq ($(FC),xlf90)
  LINKER= xlf90
  FOPT= -qstrict -qtune=ppc970 -qarch=ppc970 -qcache=auto -q64 -WF,-D__USE_LARGEFILE64
  FFLAGS= -qfree=f90 -qrealsize=8 -qextname=flush $(FOPT)
  LINKER= $(FC)
  LFLAGS= $(FFLAGS) 
endif
#
# Lists
#
FMODS=   def_master.o KindType.o InpOut.o def_netCDF.o mod_usenetcdf.o
ifeq ($(NETCDF),YES)
    FOBJS=  runend.o Turnon.o Doiter.o Reapro.o Solite.o Solvtr.o Turnof.o Wind.o readinp.o Begste.o Postpr.o GABLS.o read_meso.o netCDF_output.o Restart.o Endste.o read_qrad.o
else
    FOBJS=  runend.o Turnon.o Doiter.o Reapro.o Solite.o Solvtr.o Turnof.o Wind.o readinp.o Begste.o Postpr.o GABLS.o Restart.o Endste.o read_qrad.o
endif


  
.SUFFIXES:.o .f90
.f90.o:
	 $(FC) -c $(FFLAGS) $(PFLAGS) ${NetCDF_FLAG} ${INC_NetCDF} ${DEFINES} $<
#
#----------------------------------------------
#  Tasks
#----------------------------------------------
#
Wind1d: $(FMODS) $(FOBJS)
	$(LINKER) $(LFLAGS) -o $(EXE) $(FMODS) $(FOBJS) $(LIB_NetCDF)
	@chmod 770 $(EXE)
	@echo '---------------------------->>> END OF COMPILATION'
#
new:
	@rm -rf *.o *.mod core*
	@make
#
clean:
	@rm -f core core.* *~ *.o *.mod
