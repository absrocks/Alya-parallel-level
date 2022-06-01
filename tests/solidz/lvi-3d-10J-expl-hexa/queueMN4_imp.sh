#!/bin/bash
#
#  Submit jobs in MN-IV
#     sbatch < job.sh
#
#SBATCH --job-name=lvi-3d-expl
#SBATCH -D .
#SBATCH --error=%j.err
#SBATCH --output=%j.out
#SBATCH --qos=debug
#SBATCH --ntasks=10
#SBATCH --tasks-per-node=48
#SBATCH --time=00:05:00

#
ALYAPATH="/home/bsc21/bsc21946/alya/master/Executables/plepp/Alya.x"
#ALYAPATH="/home/bsc21/bsc21946/alya/solidz-rbo/Executables/plepp/Alya.x"
MODEL1=block
MODEL2=impactor
CPUSM1=5
CPUSM2=5
#
# Launches ALYA
#
mpirun -np $CPUSM1 $ALYAPATH $MODEL1 --name SOLID : -np $CPUSM2 $ALYAPATH $MODEL2 --name FLUID
