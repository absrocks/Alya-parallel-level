#!/bin/bash
#SBATCH --job-name=my_run
#SBATCH --output=alya.out
#SBATCH --error=alya.err
#SBATCH --ntasks=48
#SBATCH --time=00:10:00
#SBATCH -c 4

module load ompss
module load dlb

# Check bindings with this script
#srun --cpu_bind=cores nanox-bindings
srun --cpu_bind=cores $ALYA_DIR/Executables/unix/script/ompss_dlb.sh Alya.x my_run

