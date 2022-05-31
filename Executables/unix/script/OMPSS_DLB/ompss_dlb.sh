#!/bin/bash

#Nanos configuration flags
export NX_ARGS+=" --stack-size=64M"

#Nanox configuration flags for DLB
export NX_ARGS+=" --summary --force-tie-master --enable-dlb --enable-block "

#DLB configuration flags
export DLB_ARGS=" --lewi --lewi-mpi "

export LD_PRELOAD=$DLB_HOME/lib/libdlb_mpif.so # For Fortran apps

## Run the desired program
$*
 
