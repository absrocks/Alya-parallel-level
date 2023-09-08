#!/bin/csh
#
#  Submit jobs in MN-III  
#     bsub < Script
#
#BSUB -J  jube_bench_n#NODES#_N#TASKS#
#BSUB -W  2:00
#BSUB -q  bsc_case
#BSUB -n  #TASKS#
#BSUB -oo output_%J.out
#BSUB -eo output_%J.err
#BSUB -x
#
set PROBLEMNAME=#PROBLEMNAME#
#
# Launches ALYA
#

module purge
module load intel/16.0.1  
module load impi/5.1.2.150 
module load transfer/1.0
module load bsc/current

date
time mpirun ../app/Alya.x ${PROBLEMNAME}
date

