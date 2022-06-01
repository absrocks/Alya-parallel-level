#!/bin/bash
#
#  Submit jobs in MN-IV
#     sbatch < job.sh
#
#
#SBATCH --job-name=filam
#SBATCH --workdir=.
#SBATCH --error=%j.err
#SBATCH --output=%j.out
#SBATCH --qos=debug
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=4
#SBATCH --tasks-per-node=8
#SBATCH --time=00:01:00

ALYAPATH="/home/bsc21/bsc21946/alya"
PROBLEMNAME=filament_sstate
#
# Launches ALYA
#
srun $ALYAPATH/Executables/unix/Alya.x $PROBLEMNAME
#mpirun -np 4 $ALYAPATH/Executables/unix/Alya.x $PROBLEMNAME
