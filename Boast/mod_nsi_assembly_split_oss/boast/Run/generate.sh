#!/bin/bash

#Environment

BOASTDIR=".."
ALYADIR="../../../.."

#Parameters by default

OFLAGS="-O3"
UNROLL=false
INLINE=:inlined
VSIZE=4
EXECDIR=unix

#Help

help_script()
{
    cat << EOF
Usage: $0 [options]

Script to generate BOAST split function

OPTIONS:
   -h      --help                           Show this message
   -O2/-O3                                  Optimization level (2/3)
   -u      --unroll                         Enable unrolling
   -i      --inline [inlined/call/included] Inlining
   -v      --vsize [number]                 Vector size
   -e      --executable                     Executable directory (unix by default)
EOF
}

export_config()
{
echo '#START_SECTION_BOAST' > config_boast.txt
echo 'CSALYA := $(CSALYA) '"-DBOAST" >> config_boast.txt
echo 'CSALYA := $(CSALYA) '"-DVECTOR_SIZE=$VSIZE" >> config_boast.txt
echo "FOPT = $OFLAGS" >> config_boast.txt
echo '#END_SECTION_BOAST' >> config_boast.txt
cp $EXECDIR/config.in $EXECDIR/config.in.bak
line_start=`grep -nr 'START_SECTION_BOAST' $EXECDIR/config.in | cut -d : -f 1`
line_end=`grep -nr 'END_SECTION_BOAST' $EXECDIR/config.in | cut -d : -f 1`
sed -i "${line_start},${line_end}d" $EXECDIR/config.in
let line_start--
sed -i "${line_start}r config_boast.txt" $EXECDIR/config.in
}

#Managing arguments

while [ $# -gt 0 ]
do
  if [ $1 = "-h" ] || [ $1 = "--help" ]
  then
    help_script
    exit 0
  elif [ $1 = "-O2" ]
  then
    OFLAGS="-O2"
  elif [ $1 = "-O3" ]
  then
    OFLAGS="-O3"
  elif [ $1 = "-u" ] || [ $1 = "--unroll" ]
  then
    UNROLL=true
  elif [ $1 = "-i" ] || [ $1 = "--inline" ]
  then
    shift
    INLINE=:$1
  elif [ $1 = "-v" ] || [ $1 = "--vsize" ]
  then
    shift
    VSIZE=$1
  elif [ $1 = "-e" ] || [ $1 = "--executable" ]
  then
    shift
    EXECDIR=$1
  fi
  shift
done

rm nsi_element_assembly_split_oss_boast.f90
rm $ALYADIR/Sources/modules/nastin/nsi_element_assembly_split_oss_boast.f90
EXECDIR=$ALYADIR/Executables/$EXECDIR
export_config
ruby $BOASTDIR/Split/run_boast.rb --vector_size=$VSIZE --unroll=$UNROLL --inline=$INLINE --oflags=$OFLAGS --lang="FORTRAN"
cp nsi_element_assembly_split_oss_boast.f90 $ALYADIR/Sources/modules/nastin/
