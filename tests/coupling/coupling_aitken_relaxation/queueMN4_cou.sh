#!/bin/bash
#
#  Submit jobs in MN-IV
#     sbatch < job.sh
#
#SBATCH --job-name=iqn
#SBATCH -D .
#SBATCH --error=%j.err
#SBATCH --output=%j.out
#SBATCH --qos=debug
#SBATCH --ntasks=6
#SBATCH --tasks-per-node=48
#SBATCH --time=00:10:00

module purge
module load gcc/8.1.0
module load openmpi/3.1.1
module load mkl/2019.2
#
#ALYAPATH="/home/bsc21/bsc21946/alya/Executables/unix/Alya.x"
ALYAPATH="/home/bsc21/bsc21946/solidz-2.6-gnu/Executables/unix/Alya.x"

MODEL1=angrybird_sld
MODEL2=angrybird_fld
CPUSM1=3
CPUSM2=3
#
# Launches ALYA
#
mpirun -np $CPUSM1 $ALYAPATH $MODEL1 : -np $CPUSM2 $ALYAPATH $MODEL2
