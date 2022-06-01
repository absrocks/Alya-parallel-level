#!/bin/bash
#
#  Submit jobs in MN-IV
#     sbatch < job.sh
#
#
#SBATCH --job-name=sm1521T
#SBATCH --workdir=.
#SBATCH --error=%j.err
#SBATCH --output=%j.out
#SBATCH --qos=debug
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=4
#SBATCH --tasks-per-node=8
#SBATCH --time=00:01:00

module purge
module load intel/2018.1
module load impi/2018.1

ALYAPATH="/home/bsc21/bsc21946/solidz-dev"
#ALYAPATH="/home/bsc21/bsc21946/alya"
PROBLEMNAME=sm152_cube_d1T
#
# Launches ALYA
#
srun $ALYAPATH/Executables/unix/Alya.x $PROBLEMNAME
#mpirun -np 6 $ALYAPATH/Executables/unix/Alya.x $PROBLEMNAME
