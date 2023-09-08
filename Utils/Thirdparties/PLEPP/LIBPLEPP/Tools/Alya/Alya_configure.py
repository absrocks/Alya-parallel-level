#!/usr/bin/python
#------------------------------------------------------------------------------#
CONFIGURE ="""
./configure -%s -f=configure.in/%s parall exmedi solidz 
"""

INTEL01 = """#
#>> configure-mpif90-ifort-01-traceback.txt
#
f90   ==   mpif90  -traceback -module $O -c
f77   ==   mpif90  -traceback -module $O -c
fpp90 ==   mpif90  -traceback -module $O -c  -fpp
fomp90==   mpif90  -traceback -module $O -c  -fpp
fpp77 ==   mpif90  -traceback -module $O -c  -fpp
#cpp  ==   mpicc   -no-multibyte-chars -c
cpp   ==   mpicc   -c
link  ==   mpif90  -traceback

libs  ==    -L../../Thirdparties/metis-4.0 -lmetis
opti  ==    -O1

fa2p  == mpif90 -module ../../Utils/user/alya2pos -c  -fpp
"""

GCC01 = """#
#>> configure-macos-mpi-gfortran-optimized.txt
#
f90==     mpif90 -I /usr/include -m64  -c -J$O -I$O  -ffree-line-length-none
fpp90==   mpif90 -I /usr/include -m64  -c -x f95-cpp-input -J$O -I$O -ffree-line-length-none
fomp90==  mpif90 -I /usr/include -m64  -c -x f95-cpp-input -J$O -I$O  -ffree-line-length-none
opti==   -O1
cpp==     gcc -c
link==    mpif90 -m64  -lc
libs==   -L../../Thirdparties/metis-4.0 -lmetis

fa2p== mpif90 -m64 -c -x f95-cpp-input -DMPI_OFF -J../../Utils/user/alya2pos -I../../Utils/user/alya2pos
"""

YOSEMITE = """#
#>> configure-macos-mpi-gfortran-optimized.txt
#
# 1.0) 
#  sudo port install mpich
#
# 2.1) /Executables/unix/makefile/configure.in/configure-macos-mpi-gfortran-optimized.txt
#     mpiXXX  -> mpiXXX-mpich-mp
#
# 3.0) /Executables/unix/makefile  
#     To change @$(MAKE) -C ../../Thirdparties/metis-4.0/  --> @$(MAKE) -C ../../Thirdparties/metis-4.0/Lib/ CC=mpicc-mpich-mp
#
# 4.0) /Executables/unix/makefile metis4        
#
# 5.0) /Executables/unix/makefile -j 4 
#
# 6.0) enjoy 'Alya.x' 
#      :) 
#
# 7.0) /Executables/unix/makefile alya2pos.x 
# 7.1) enjoy '../../Utils/user/alya2pos/alya2pos.x' 
#      :) 
#
# 8.0) /Executables/unix/makefile   
#       @cd $(LIBPLE); ./configure CC=mpicc --prefix=$(LIBPLE)/Execs -> @cd $(LIBPLE); ./configure CC=mpicc-mpich-mp --prefix=$(LIBPLE)/Execs
#       @$(MAKE) -C $(LIBPLE)/src CFLAGS_OPT='' 
#
# 9.0) /Executables/unix/makefile
#  @$(MAKE) -C $(LIBPLEPP)/Wrappers/Cpp     static HOME=$(LIBPLEPP) ROOT_PLE=$(LIBPLE)/Execs CPPCOMPILER=mpic++-mpich-mp CCOMPILER=mpicc-mpich-mp FCOMPILER=mpif90-mpich-mp
#  @$(MAKE) -C $(LIBPLEPP)/Wrappers/Fortran static HOME=$(LIBPLEPP) ROOT_PLE=$(LIBPLE)/Execs CPPCOMPILER=mpic++-mpich-mp CCOMPILER=mpicc-mpich-mp FCOMPILER=mpif90-mpich-mp
#
#10.0) /Executables/unix/makefile
#
#11.0) /Executables/unix/makefile COMMDOM=[1|2]
# 
#12.0) enjoy 'alya+plepp'
#      time mpirun-mpich-mp -np $ALYAi $ALYA_PATH/Executables/plepp/Alya.x Small01 --name FLUID  : -np $ALYAj $ALYA_PATH/Executables/plepp/Alya.x Big01 --name SOLID
#      :) 
#

f90==     mpif90-mpich-mp -I /usr/include -m64  -c -J$O -I$O  -ffree-line-length-none
fpp90==   mpif90-mpich-mp -I /usr/include -m64  -c -x f95-cpp-input -J$O -I$O -ffree-line-length-none
fomp90==  mpif90-mpich-mp -I /usr/include -m64  -c -x f95-cpp-input -J$O -I$O  -ffree-line-length-none
opti==    -O1
cpp==     gcc -c
link==    mpif90-mpich-mp -m64  -lc
libs==    -L../../Thirdparties/metis-4.0 -lmetis

fa2p== mpif90-mpich-mp -m64 -c -x f95-cpp-input -DMPI_OFF -J../../Utils/user/alya2pos -I../../Utils/user/alya2pos

"""

#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
import argparse
parser = argparse.ArgumentParser()
#
#parser.add_argument('-X', action='append', dest='XXX',
#                    default=[], type=str, nargs="+",
#                    help='gcc|intel|yosemite',
#                    )
#
parser.add_argument('-T', action='store', dest='Type',
                    default='x', type=str, nargs=1,
                    help='x|g',
                    )

parser.add_argument('-C', action='store', dest='Compiler',
                    default='xxx', type=str, nargs=1,
                    help='gcc|intel|yosemite',
                    )

Results    = parser.parse_args()
Type       = Results.Type[0]
Compiler   = Results.Compiler[0]
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
if(Compiler=="gcc"):
  configure = "configure-macos-mpi-gfortran-optimized.txt"
elif(Compiler=="intel"):
  configure = "configure-mpif90-ifort-01-traceback.txt"
elif(Compiler=="yosemite"):
  configure = "mpich_at_yosemite.txt"
  Fout = open("./configure.in/" + configure, "w") 
  print>> Fout, YOSEMITE 
  Fout.close()  
else: 
  import sys
  parser.print_help()
  print "EXIT!!\n"
  sys.exit(1) 

import os
os.system( CONFIGURE % (Type, configure) ) 

print "ok! \n"
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
