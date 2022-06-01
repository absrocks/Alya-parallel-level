#!/bin/bash

usage_fct(){
  echo "Usage : $0 [options] application"
  echo "Options:"
  echo "-h --help: Print help"
}

FILE=elements.csv
rm $FILE
touch $FILE
#Global
for log in *[0-9]-partition.par.post.res 
do
  tmp=tmp
  iteration=`echo $log | cut -d '-' -f 2`
  component=tri03
  COMP='NUM\_TRI03'
  cat $log | awk "/ComponentNames\ $COMP/,/End/"> $tmp
  sed -i '1d' $tmp
  sed -i '1d' $tmp
  sed -i '$d' $tmp
  while read line; do
    num=`echo $line | tr -s ' ' | cut -d ' ' -f 1`
    tim=`echo $line | tr -s ' ' | cut -d ' ' -f 2`
    echo $iteration, $component, $num, $tim >> $FILE 
  done < "$tmp"
  component=qua04
  COMP='NUM\_QUA04'
  cat $log | awk "/ComponentNames\ $COMP/,/End/"> $tmp
  sed -i '1d' $tmp
  sed -i '1d' $tmp
  sed -i '$d' $tmp
  while read line; do
    num=`echo $line | tr -s ' ' | cut -d ' ' -f 1`
    tim=`echo $line | tr -s ' ' | cut -d ' ' -f 2`
    echo $iteration, $component, $num, $tim >> $FILE 
  done < "$tmp"
  component=tet04
  COMP='NUM\_TET04'
  cat $log | awk "/ComponentNames\ $COMP/,/End/"> $tmp
  sed -i '1d' $tmp
  sed -i '1d' $tmp
  sed -i '$d' $tmp
  while read line; do
    num=`echo $line | tr -s ' ' | cut -d ' ' -f 1`
    tim=`echo $line | tr -s ' ' | cut -d ' ' -f 2`
    echo $iteration, $component, $num, $tim >> $FILE 
  done < "$tmp"
  component=pyr05
  COMP='NUM\_PYR05'
  cat $log | awk "/ComponentNames\ $COMP/,/End/"> $tmp
  sed -i '1d' $tmp
  sed -i '1d' $tmp
  sed -i '$d' $tmp
  while read line; do
    num=`echo $line | tr -s ' ' | cut -d ' ' -f 1`
    tim=`echo $line | tr -s ' ' | cut -d ' ' -f 2`
    echo $iteration, $component, $num, $tim >> $FILE 
  done < "$tmp"
  component=pen06
  COMP='NUM\_PEN06'
  cat $log | awk "/ComponentNames\ $COMP/,/End/"> $tmp
  sed -i '1d' $tmp
  sed -i '1d' $tmp
  sed -i '$d' $tmp
  while read line; do
    num=`echo $line | tr -s ' ' | cut -d ' ' -f 1`
    tim=`echo $line | tr -s ' ' | cut -d ' ' -f 2`
    echo $iteration, $component, $num, $tim >> $FILE 
  done < "$tmp"

done
