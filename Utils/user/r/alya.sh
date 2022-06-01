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
case=$1

rm global.csv
rm nsi.csv

csv=global.csv
#Global
for log in $case-[0-9].log
do
  iteration=`echo $log | cut -d '-' -f 2 | cut -d '.' -f 1`
  tmp=$log.tmp
  cat $log | awk '/TOTAL\ CPU \TIME/,/\|-\ END/'> $tmp
  sed -i '$d' $tmp
  sed -i '$d' $tmp
  sed -i '$d' $tmp
  #GLOBAL
  #cpu time
  a=`grep "CPU" $tmp | cut -d':' -f2 | tr -d ' '`
  b=`grep "STARTING" $tmp | cut -d':' -f2 | cut -d'(' -f1 | tr -d ' '`
  c=`grep "Geometry" $tmp | cut -d':' -f2 | cut -d'(' -f1 | tr -d ' '`
  d=`grep "Set" $tmp | cut -d':' -f2 | cut -d'(' -f1 | tr -d ' '`
  e=`grep "Boundary" $tmp | cut -d':' -f2 | cut -d'(' -f1 | tr -d ' '`
  f=`grep "partitioning" $tmp | cut -d':' -f2 | cut -d'(' -f1 | tr -d ' '`
  g=`grep "distribution" $tmp | cut -d':' -f2 | cut -d'(' -f1 | tr -d ' '`
  h=`grep "multiplication" $tmp | cut -d':' -f2 | cut -d'(' -f1 | tr -d ' '`
  i=`grep "construction" $tmp | cut -d':' -f2 | cut -d'(' -f1 | tr -d ' '`
  j=`grep "arrays" $tmp | cut -d':' -f2 | cut -d'(' -f1 | tr -d ' '`
  #header:
  echo "$platform, $mpi, $openmp, $vec, $vecgpu, $cores_per_node, $cores_per_gpu, $processes_per_gpu, $tag, $iteration, total, total, total, $a" >> $csv
  echo "$platform, $mpi, $openmp, $vec, $vecgpu, $cores_per_node, $cores_per_gpu, $processes_per_gpu, $tag, $iteration, start, start, total, $b" >> $csv
  echo "$platform, $mpi, $openmp, $vec, $vecgpu, $cores_per_node, $cores_per_gpu, $processes_per_gpu, $tag, $iteration, start, start, geome, $c" >> $csv
  echo "$platform, $mpi, $openmp, $vec, $vecgpu, $cores_per_node, $cores_per_gpu, $processes_per_gpu, $tag, $iteration, start, start, setre, $d" >> $csv
  echo "$platform, $mpi, $openmp, $vec, $vecgpu, $cores_per_node, $cores_per_gpu, $processes_per_gpu, $tag, $iteration, start, start, bound, $e" >> $csv
  echo "$platform, $mpi, $openmp, $vec, $vecgpu, $cores_per_node, $cores_per_gpu, $processes_per_gpu, $tag, $iteration, start, start, parti, $f" >> $csv
  echo "$platform, $mpi, $openmp, $vec, $vecgpu, $cores_per_node, $cores_per_gpu, $processes_per_gpu, $tag, $iteration, start, start, distr, $g" >> $csv
  echo "$platform, $mpi, $openmp, $vec, $vecgpu, $cores_per_node, $cores_per_gpu, $processes_per_gpu, $tag, $iteration, start, start, multi, $h" >> $csv
  echo "$platform, $mpi, $openmp, $vec, $vecgpu, $cores_per_node, $cores_per_gpu, $processes_per_gpu, $tag, $iteration, start, start, const, $i" >> $csv
  echo "$platform, $mpi, $openmp, $vec, $vecgpu, $cores_per_node, $cores_per_gpu, $processes_per_gpu, $tag, $iteration, start, start, array, $j" >> $csv
  for MOD in `grep "MODULE" $tmp| cut -d':' -f1 | tr -d ' ' | sed 's/MODULE//g'`
  do
    #echo $MOD
    mod=`echo $MOD | tr '[:upper:]' '[:lower:]'`
    line_start=`grep -nr $MOD $tmp | gawk '{print $1}' FS=":"`
    line=${line_start}
    a=`sed -n "${line}p" $tmp | cut -d':' -f2 | cut -d'(' -f1 | tr -d ' '`
    line=`expr $line + 2`
    b=`sed -n "${line}p" $tmp | cut -d':' -f2 | cut -d'(' -f1 | tr -d ' '`
    line=`expr $line + 1`
    c=`sed -n "${line}p" $tmp | cut -d':' -f2 | cut -d'(' -f1 | tr -d ' '`
    line=`expr $line + 1`
    d=`sed -n "${line}p" $tmp | cut -d':' -f2 | cut -d'(' -f1 | tr -d ' '`
    line=`expr $line + 2`
    e=`sed -n "${line}p" $tmp | cut -d':' -f2 | cut -d'(' -f1 | tr -d ' '`
    line=`expr $line + 1`
    f=`sed -n "${line}p" $tmp | cut -d':' -f2 | cut -d'(' -f1 | tr -d ' '`
    line=`expr $line + 1`
    g=`sed -n "${line}p" $tmp | cut -d':' -f2 | cut -d'(' -f1 | tr -d ' '`
    line=`expr $line + 2`
    h=`sed -n "${line}p" $tmp | cut -d':' -f2 | cut -d'(' -f1 | tr -d ' '`
    line=`expr $line + 1`
    i=`sed -n "${line}p" $tmp | cut -d':' -f2 | cut -d'(' -f1 | tr -d ' '`
    line=`expr $line + 1`
    j=`sed -n "${line}p" $tmp | cut -d':' -f2 | cut -d'(' -f1 | tr -d ' '`
    echo "$platform, $mpi, $openmp, $vec, $vecgpu, $cores_per_node, $cores_per_gpu, $processes_per_gpu, $tag, $iteration, $mod, total, total, $a" >> $csv
    echo "$platform, $mpi, $openmp, $vec, $vecgpu, $cores_per_node, $cores_per_gpu, $processes_per_gpu, $tag, $iteration, $mod, matri, avera, $b" >> $csv
    echo "$platform, $mpi, $openmp, $vec, $vecgpu, $cores_per_node, $cores_per_gpu, $processes_per_gpu, $tag, $iteration, $mod, matri, maxim, $c" >> $csv
    echo "$platform, $mpi, $openmp, $vec, $vecgpu, $cores_per_node, $cores_per_gpu, $processes_per_gpu, $tag, $iteration, $mod, matri, loadb, $d" >> $csv
    echo "$platform, $mpi, $openmp, $vec, $vecgpu, $cores_per_node, $cores_per_gpu, $processes_per_gpu, $tag, $iteration, $mod, solve, avera, $e" >> $csv
    echo "$platform, $mpi, $openmp, $vec, $vecgpu, $cores_per_node, $cores_per_gpu, $processes_per_gpu, $tag, $iteration, $mod, solve, maxim, $f" >> $csv
    echo "$platform, $mpi, $openmp, $vec, $vecgpu, $cores_per_node, $cores_per_gpu, $processes_per_gpu, $tag, $iteration, $mod, solve, loadb, $g" >> $csv
    echo "$platform, $mpi, $openmp, $vec, $vecgpu, $cores_per_node, $cores_per_gpu, $processes_per_gpu, $tag, $iteration, $mod, postp, avera, $h" >> $csv
    echo "$platform, $mpi, $openmp, $vec, $vecgpu, $cores_per_node, $cores_per_gpu, $processes_per_gpu, $tag, $iteration, $mod, postp, maxim, $i" >> $csv
    echo "$platform, $mpi, $openmp, $vec, $vecgpu, $cores_per_node, $cores_per_gpu, $processes_per_gpu, $tag, $iteration, $mod, postp, loadb, $j" >> $csv
  done
  rm $tmp
done
