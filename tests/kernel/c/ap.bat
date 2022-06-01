#!/bin/bash 
#BSUB -n 1
#BSUB -o ap.out
#BSUB -e ap.err
#BSUB -R"span[ptile=16]"
#BSUB -J ap
#BSUB -W 1:00
/gpfs/projects/bsc21/WORK-HERBERT/svnmn3/Alya/Utils/user/alya2pos/alya2pos.x c

