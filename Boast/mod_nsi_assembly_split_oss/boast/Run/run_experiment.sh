#!/bin/bash

BASE=".."
EXP_DIR="$BASE/experiments"
SRC_DIR="$BASE"

help_script()
{
    cat << EOF
Usage: $0 [options] /path_to_data_directory

Script for to get machine information before doing the experiment

OPTIONS:
   -h      Show this message
EOF
}

while getopts "d:" opt; do
    case $opt in
	d)
	    DATADIR="$OPTARG"
	    ;;
	h)
	    help_script
	    exit 4
	    ;;
	\?)
	    echo "Invalid option: -$OPTARG"
	    help_script
	    exit 3
	    ;;
    esac
done

DATA_FOLD_DAY=`date +%Y_%m_%d`
DATA_FOLD_DAY="$EXP_DIR/$DATA_FOLD_DAY"
BKUP=`date +%H_%M_%S`
DATA_FOLD_HOST=`hostname`
DATA_FOLD_HOST="$DATA_FOLD_DAY/$DATA_FOLD_HOST"
DATA_FOLD_TIME="$DATA_FOLD_HOST/$BKUP"
INFO_FILE="$DATA_FOLD_TIME/Info.yaml"
DATA_FILE="$DATA_FOLD_TIME/Data.yaml"

# Set vector size (Beware of KERNEL_FILE name)
VSIZE=1
F_KERNEL_FILE="$DATA_FOLD_TIME/nsi_element_assembly_split_oss_boast.f90"
C_KERNEL_FILE="$DATA_FOLD_TIME/nsi_element_assembly_split_oss_boast.c"
SETLANG=FORTRAN #FORTRAN #C

if [[ "$SETLANG" == "C" ]]; then
	KERNEL_FILE=$C_KERNEL_FILE
else
	KERNEL_FILE=$F_KERNEL_FILE
fi


mkdir -p $DATA_FOLD_DAY
mkdir -p $DATA_FOLD_HOST
mkdir -p $DATA_FOLD_TIME



### Choose the script you want to run ###
# Beware of the language you want


###########################################
## (1) Verification of the Boast version ##  
###########################################
#echo "EXECUTING ruby $SRC_DIR/Split/run_VERIF.rb --data=$DATA_FILE --info=$INFO_FILE --kernel=$KERNEL_FILE"
#ruby $SRC_DIR/Split/run_VERIF.rb --data=$DATA_FILE --info=$INFO_FILE --kernel=$KERNEL_FILE

###############################################
## (2) Ref VS Boast with different options   ##
###############################################
########  (2)(a) Uncomment to test all nests
#echo "(1) EXECUTING ruby $SRC_DIR/Split/run_COMP.rb --data=$DATA_FILE --info=$INFO_FILE --lang=$SETLANG ..."
#ruby $SRC_DIR/Split/run_COMP.rb --data=$DATA_FILE --info=$INFO_FILE --lang=$SETLANG
#echo "(2) EXECUTING ruby format_data_comparison.rb $DATA_FILE test.cvs ..."
#ruby format_data_comparison.rb $DATA_FILE test.cvs
#echo "(3) PARSING ..."
#sh parse.sh
#echo "(4) GENERATING FIGURE ..."
#gnuplot results.gplot
#########  (2)(b) Uncomment to test each nest separately
#echo "EXECUTING ruby $SRC_DIR/Split/run_COMP_nests.rb --data=$DATA_FILE --info=$INFO_FILE --kernel=$KERNEL_FILE --lang=$SETLANG"
#ruby $SRC_DIR/Split/run_COMP_nests.rb --data=$DATA_FILE --info=$INFO_FILE --kernel=$KERNEL_FILE  --lang=$SETLANG


############################################################
## (3) Comparison of Boast versions (with different options)
############################################################
#echo "EXECUTING ruby $SRC_DIR/Split/run_boast_COMP.rb --data=$DATA_FILE --info=$INFO_FILE --kernel=$KERNEL_FILE --vector_size=$VSIZE --lang=$SETLANG"
#ruby $SRC_DIR/Split/run_boast_COMP.rb --data=$DATA_FILE --info=$INFO_FILE --kernel=$KERNEL_FILE --vector_size=$VSIZE --lang=$SETLANG

#######################################################
## (4) Generate the Boast kernel with specific options
#######################################################
#  Set the following variables:
OFLAGS="-O3"
UNROLL=true
INLINE=:included
echo "EXECUTING ruby $SRC_DIR/Split/run_boast.rb --vector_size=$VSIZE --unroll=$UNROLL --inline=$INLINE --oflags=$OFLAGS --lang=$SETLANG"
ruby $SRC_DIR/Split/run_boast.rb --vector_size=$VSIZE --unroll=$UNROLL --inline=$INLINE --oflags=$OFLAGS --lang=$SETLANG


############################
## (5) Test the Ref version
############################
#echo "EXECUTING ruby $SRC_DIR/Split/run_ref.rb --vector_size=$VSIZE"
#ruby $SRC_DIR/Split/run_ref.rb --vector_size=$VSIZE
