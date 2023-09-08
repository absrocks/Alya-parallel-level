#! /bin/bash
#BSUB -n 1
#BSUB -o alya-make.out
#BSUB -e alya-make.err
#BSUB -R"span[ptile=16]"
#BSUB -J alya-make
#BSUB -W 00:30
echo '--| Starting alya-make process at : ' `date`
time make -j 8
echo '--| End alya-make at              : ' `date`
