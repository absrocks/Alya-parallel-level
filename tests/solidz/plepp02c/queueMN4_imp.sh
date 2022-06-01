#!/bin/bash
#
#  Submit jobs in MN-IV
#     sbatch < job.sh
#
#SBATCH --job-name=plepp02c
#SBATCH -D .
#SBATCH --error=%j.err
#SBATCH --output=%j.out
#SBATCH --qos=debug
#SBATCH --ntasks=6
#SBATCH --tasks-per-node=48
#SBATCH --time=00:10:00

#
#ALYAPATH="/home/bsc21/bsc21946/alya/master/Executables/plepp/Alya.x"
#ALYAPATH="/home/bsc21/bsc21946/alya/solidz-rbo/Executables/plepp/Alya.x"
MODEL1=block
MODEL2=identer
CPUSM1=3
CPUSM2=3
#
# Launches ALYA
#
mpirun -np $CPUSM1 $ALYAPATH $MODEL1 --name SOLID : -np $CPUSM2 $ALYAPATH $MODEL2 --name FLUID
