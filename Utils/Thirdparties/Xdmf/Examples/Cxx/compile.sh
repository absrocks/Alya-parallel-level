#
# Temporal script to compile XdmfFortranExample.f90 file
# Please modify this script to adjust it to your requirements
#

#!/bin/bash

#mpif90 -g XdmfFortranExample.f90  ~/local/lib/libXdmf.a ~/local/lib/libhdf5.a /apps/SZIP/2.1/lib/libsz.a -lxml2 -I ~/local/lib -L /opt/intel/composerxe/lib/intel64/ -lifcore -lstdc++ -o XdmfFortranExample
mpif90 -g XdmfFortranExample.f90  ~/local/lib/libXdmf.a ~/local/lib/libhdf5.a /apps/SZIP/2.1/lib/libsz.a -lxml2 -I ~/local/lib -I ~/local/include -L ~/local/lib -L /opt/intel/composerxe/lib/intel64/ -lifcore -lstdc++ -lhdf5_fortran -lhdf5 -o XdmfFortranExample
