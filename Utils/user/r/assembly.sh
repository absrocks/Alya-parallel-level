#!/bin/bash

usage_fct(){
  echo "Usage : $0 [options] application"
  echo "Options:"
  echo "-h --help: Print help"
}

platform=unknown
mpi=unknown
openmp=unknown
vec=unknown
vecgpu=unknown
cores_per_gpu=unknown
cores_per_node=unknown
processes_per_gpu=unknown
tag=unknown
while [ $# -gt 1 ]
do
  if [ $1 = "-h" ] || [ $1 = "--help" ]
  then
    usage_fct
    exit 0
  elif [ $1 = "-cg" ] || [ $1 = "--cores-per-gpu" ]
  then
    shift
    cores_per_gpu=$1
  elif [ $1 = "-cn" ] || [ $1 = "--cores-per-node" ]
  then
    shift
    cores_per_node=$1
  elif [ $1 = "-pg" ] || [ $1 = "--processes-per-gpu" ]
  then
    shift
    processes_per_gpu=$1
  elif [ $1 = "-m" ] || [ $1 = "--mpi" ]
  then
    shift
    mpi=$1
  elif [ $1 = "-o" ] || [ $1 = "--openmp" ]
  then
    shift
    openmp=$1
  elif [ $1 = "-v" ] || [ $1 = "--vector-size" ]
  then
    shift
    vec=$1
  elif [ $1 = "-u" ] || [ $1 = "--gpu-vector-size" ]
  then
    shift
    vecgpu=$1
  elif [ $1 = "-t" ] || [ $1 = "--tag" ]
  then
    shift
    tag=$1
  elif [ $1 = "-p" ] || [ $1 = "--platform" ]
  then
    shift
    platform=$1
  fi
  shift
done
FILE=assembly.csv
rm $FILE
touch $FILE
#Global
for log in *[0-9]-partition.par.post.res 
do
  tmp=tmp
  iteration=`echo $log | cut -d '-' -f 2`
  component=assembly_time_nastin
  COMP='ELEMENT\_ASSEMBLY\_TIME\_NASTIN'
  cat $log | awk "/ComponentNames\ $COMP/,/End/"> $tmp
  sed -i '1d' $tmp
  sed -i '1d' $tmp
  sed -i '$d' $tmp
  while read line; do
    num=`echo $line | tr -s ' ' | cut -d ' ' -f 1`
    tim=`echo $line | tr -s ' ' | cut -d ' ' -f 2`
    echo $platform, $mpi, $openmp, $vec, $vecgpu, $cores_per_node, $cores_per_gpu, $processes_per_gpu, $tag, $component, $iteration, $num, $tim >> $FILE 
  done < "$tmp"
  component=num_elements_partition
  COMP='NUM\_ELEMENTS'
  cat $log | awk "/ComponentNames\ $COMP/,/End/"> $tmp
  sed -i '1d' $tmp
  sed -i '1d' $tmp
  sed -i '$d' $tmp
  while read line; do
    num=`echo $line | tr -s ' ' | cut -d ' ' -f 1`
    tim=`echo $line | tr -s ' ' | cut -d ' ' -f 2`
    echo $platform, $mpi, $openmp, $vec, $vecgpu, $cores_per_node, $cores_per_gpu, $processes_per_gpu, $tag, $component, $iteration, $num, $tim >> $FILE 
  done < "$tmp"
done
