#!/bin/bash 
#BSUB -n 16 
#BSUB -o alya.out
#BSUB -e Error_alya.txt
#BSUB -J  10
#BSUB -W 0:19
### You can choose the parallel environment through
mpirun /gpfs/projects/bsc21/WORK-HERBERT/svnmn3/Alya/Executables/unix/Alya.x c
###/gpfs/projects/bsc21/WORK-HERBERT/svnmn3/Alya/Executables/unix_deb/Alya.g c
###mpirun valgrind -v --gen-suppressions=all --suppressions=${ALYA_DIR}/Utils/user/valgrind/alya.supp  --suppressions=/apps/INTEL/2016.1.056/itac/9.1.2.024/lib/impi.supp --leak-check=full  --track-origins=yes  --log-file=valgrind.out ${ALYA_DIR}/Executables/unixm5i8_deb_nd3/Alya.g c
