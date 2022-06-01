#!/bin/bash
#SBATCH --job-name=ts
#SBATCH --output=output.out
#SBATCH --error=output.err
#SBATCH --ntasks=4
#SBATCH --qos=debug
#SBATCH --time=00:02:00
export CASE=${PWD##*/}  
export ROMIO_HINTS=./io_hints
export I_MPI_EXTRA_FILESYSTEM_LIST=gpfs
export I_MPI_EXTRA_FILESYSTEM=on
module unload intel
module load gcc
mpirun $ALYA/alya/Executables/unix-i8-TS/Alya.x $CASE
mpirun mpio2txt $CASE
mv *.post.mpio.txt base/4p/
