#!/bin/bash

#Environment

BASE=".."
EXP_DIR="$BASE/experiments"
SRC_DIR="$BASE"

#Parameters by default

LANG=FORTRAN

#Help

help_script()
{
    cat << EOF
Usage: $0 [options]

Script to generate BOAST split function

OPTIONS:
   -h      --help                           Show this message
   -l      --language                       Language (C/FORTRAN)
EOF
}

#Managing arguments

while [ $# -gt 0 ]
do
  if [ $1 = "-h" ] || [ $1 = "--help" ]
  then
    help_script
    exit 0
  elif [ $1 = "-l" ] || [ $1 = "--language" ]
  then
    shift
    LANG=$1
  fi
  shift
done

DATA_FOLD_DAY=`date +%Y_%m_%d`
DATA_FOLD_DAY="$EXP_DIR/$DATA_FOLD_DAY"
BKUP=`date +%H_%M_%S`
DATA_FOLD_HOST=`hostname`
DATA_FOLD_HOST="$DATA_FOLD_DAY/$DATA_FOLD_HOST"
DATA_FOLD_TIME="$DATA_FOLD_HOST/$BKUP"
INFO_FILE="$DATA_FOLD_TIME/Info.yaml"
DATA_FILE="$DATA_FOLD_TIME/Data.yaml"

mkdir -p $DATA_FOLD_DAY
mkdir -p $DATA_FOLD_HOST
mkdir -p $DATA_FOLD_TIME

echo "(1) EXECUTING ruby $SRC_DIR/Split/run_COMP.rb --data=$DATA_FILE --info=$INFO_FILE --lang=$LANG ..."
ruby $SRC_DIR/Split/run_COMP.rb --data=$DATA_FILE --info=$INFO_FILE --lang=$LANG
echo "(2) EXECUTING ruby format_data_comparison.rb $DATA_FILE test.csv ..."
ruby format_data_comparison.rb $DATA_FILE test.csv
echo "(3) GENERATING FIGURE ..."
R.sh
