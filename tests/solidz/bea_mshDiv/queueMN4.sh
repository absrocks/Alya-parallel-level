#!/bin/bash
#
#  Submit jobs in MN-IV
#     sbatch < job.sh
#
#SBATCH --job-name=beamshDiv
#SBATCH --workdir=.
#SBATCH --error=%j.err
#SBATCH --output=%j.out
#SBATCH --qos=debug
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=1
#SBATCH --tasks-per-node=48
#SBATCH --time=00:03:00
#
# Paths and problem name
#
#ALYAPATH="/home/bsc21/bsc21946/alya"
ALYAPATH="/home/bsc21/bsc21946/solidz-contaExpl"
PROBLEMNAME=bea_mshDiv
#
# Launches ALYA
#
#mpirun -np 24 $ALYAPATH/Executables/unix/Alya.x $PROBLEMNAME
srun $ALYAPATH/Executables/unix/Alya.x $PROBLEMNAME
