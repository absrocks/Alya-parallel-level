#!/bin/bash

DATE=$(date +%Y%m%d-%H%M%S)

#
# Mac
#
ALYAPATH="/Users/gguillamet/BSC/alya/master"
#ALYAPATH="/Users/gguillamet/BSC/alya/369-bug-parallel-pdn-contact"
#
# MareNostrum
#
#ALYAPATH="/home/bsc21/bsc21946/alya/master"
#ALYAPATH="/home/bsc21/bsc21946/alya/369-bug-parallel-pdn-contact"

if [ $# -eq 0 ]; then
    echo "Arguments are missing"
    echo "The argument must be the case name and the number of process for each Alya"
    echo "Example: \"./runner.sh Big01 3\" will run Big01 case using -np 3"
    echo "If you are using plepp then you can execute the multi-parallel as following:"
    echo "\"./runner Big01 FLUID Small01 SOLID 1\". This will execute Big01 as FLUID,"
    echo "Small01 as SOLID and it will use for both problems 1 processor."
    echo "Also you can use \"./runner.sh clean\" or \"./runner.sh pos CASE_NAME\""
    exit
fi
if [ $1 == 'clean' ]; then
    echo "doing alya clean..."
    $ALYAPATH/Utils/user/alya-clean
    echo "alya-clean ended succesfully"

elif [ $1 == 'pos' ]; then
    if [ $# -eq 4 ]; then
        echo -ne '\n' | $ALYAPATH/Utils/user/alya2pos/alya2pos.x $2 & $ALYAPATH/Utils/user/alya2pos/alya2pos.x $3 & $ALYAPATH/Utils/user/alya2pos/alya2pos.x $4
    else
        echo -ne '\n' | $ALYAPATH/Utils/user/alya2pos/alya2pos.x $2 & $ALYAPATH/Utils/user/alya2pos/alya2pos.x $3
    fi
elif [ $# -eq 5 ]; then
    mpirun -np $5 $ALYAPATH/Executables/plepp/Alya.x $1 --name $2 : -np $5 $ALYAPATH/Executables/plepp/Alya.x $3 --name $4 #| tee screen-$DATE.txt

elif [ $# -eq 6 ]; then
    mpirun -np $3 $ALYAPATH/Executables/plepp/Alya.x $1 --name $2 : -np $6 $ALYAPATH/Executables/plepp/Alya.x $4 --name $5 #| tee screen-$DATE.txt

elif [ $# -eq 7 ]; then
    mpirun -np $7 $ALYAPATH/Executables/plepp/Alya.x $1 --name $2 : -np $7 $ALYAPATH/Executables/plepp/Alya.x $3 --name $4 : -np $7 $ALYAPATH/Executables/plepp/Alya.x $5 --name $6 #| tee screen-$DATE.txt

else
    echo "Check number and type of arguments"
fi
