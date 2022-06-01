#!/bin/bash
#SBATCH --job-name="a2e"
#SBATCH --time=01:00:00
#SBATCH --qos=debug
#SBATCH --ntasks=48
#SBATCH --output=output_%J.out
#SBATCH --error=output_%J.err

# You can choose the parallel environment through modules
# this example is for MN-IV
module load gcc/7.2.0 impi/2018.1 mkl/2018.1
module load python/3.6.1

srun python /correct/path/to/alya2ensight-mpi.py problem_name /path/where/your/alyabin/are /path/where/your/ensi/will/be  --format alyabin
