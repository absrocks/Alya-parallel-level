#!/bin/bash

module load intel
module load mkl
module load impi
module load gcc
module load boost/1.64.0
module load vtk/8.0.1
module load cmake/3.9.2
cd ../../../Executables/unix
./configure -x parall nastin
make alya-mpio-tools
