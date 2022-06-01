#!/bin/bash
#SBATCH --job-name=cavtri03_mpio_i4_hybrid
#SBATCH --output=output.out
#SBATCH --error=output.err
#SBATCH --ntasks=1
#SBATCH --time=01:00:00
export ROMIO_HINTS=./io_hints
export I_MPI_EXTRA_FILESYSTEM_LIST=gpfs
export I_MPI_EXTRA_FILESYSTEM=on
srun /gpfs/projects/bsc21/damien/apps/alya/alya/Executables/unix-TS/Alya.x cavtri03_mpio_i4_hybrid
