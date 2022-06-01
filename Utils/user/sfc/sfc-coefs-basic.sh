#!/bin/bash

# chec parameters
if [ "$#" -ne 4 ]; then 
    echo "Four parameters required: [case name],[iteration],[folder tag],[npart]"
    exit 1
fi 

## inputs
case=$1
iter=$(( $2 -1 ))
flag=$3
npart=$4
phase=assembly

##
## check imputs
##

if (( $iter < 0 )); then 
  echo "   ERROR: 2th argument must be > 0 "
   exit 1 
fi
if (( $npart < 2 )); then 
   echo "   ERROR: 4th argument must be > 1 "
   exit 1 
fi

##
## prepare files (cor,imb,ent)
##

alyafile="rank-elements.dat"
filein="$case-partition.par.post.res"
itant=$(($iter-1))
dir="SFCopt_$flag"
mkdir -p $dir

if (( $itant \< 10 )); then
   filecorant="$dir/$1-$phase$part"0"$itant.cor"
   fileimbant="$dir/$1-$phase$part"0"$itant.imb"
else
   filecorant="$dir/$1-$phase$part$itant.cor"
   fileimbant="$dir/$1-$phase$part$itant.imb"
fi
if (( $iter \< 10 )); then
   fileimb="$dir/$1-$phase$part"0"$iter.imb"
   filecor="$dir/$1-$phase$part"0"$iter.cor"
   fileent="$dir/$1-$phase$part"0"$iter.ent"
else
   fileimb="$dir/$1-$phase$part$iter.imb"
   filecor="$dir/$1-$phase$part$iter.cor"
   fileent="$dir/$1-$phase$part$iter.ent"
fi 

if (( $iter == 0 )); then		
   for (( i=1; i<=$npart; i++ )); do echo "$i 1"; done > $filecor
   cp $filecor $alyafile 
   exit 0
fi

if [ ! -e "$filecorant" ]; then
   echo "   ERROR: Correction file for previous iteration not present: $filecorant "
   exit 1 
fi

##
## take data from filein
##

ratio=`cat $filein | awk '/ComponentNames\ ELEMENT_ASSEMBLY_RATIO_NASTIN/,/End\ Values/' | sed -e '1,3d' | sed -e '$d' | sed -e 's/^[ ]*//' | sed 's/  \+/ /g'`
laten=`cat $filein | awk '/ComponentNames\ ELEMENT_ASSEMBLY_TIME_NASTIN/,/End\ Values/' | sed -e '1,3d' | sed -e '$d' | sed -e 's/^[ ]*//' | sed 's/  \+/ /g'`
acoef=($(awk -F' ' '{printf "%.8f\n", $2}' $filecorant))
atime=($(cat $filein | awk '/ComponentNames\ ELEMENT_ASSEMBLY_TIME_NASTIN/,/End\ Values/' | sed -e '1,3d' | sed -e '$d'| sed -e 's/^[ ]*//' | sed 's/  \+/ /g' | awk -F' ' '{printf "%.8f\n", $2}'))
##
## get imblalance
##
echo "$ratio" | awk -F' ' 'BEGIN{t=1} {printf "%8s%14s\n", t, $2; t+=1}' > $fileimb
size=`wc -l $fileimb  | cut -f1 -d' '`

if (( $size != $npart )); then 
   echo "   ERROR: 6th argument wrong size $size $npart"
   exit 1 
fi

avgtime=`echo "$laten" | awk -F' ' 'BEGIN {SUM=0; t=0}; {SUM=SUM+$2;t+=1}; END{SUM=SUM/t; printf "%.2f", SUM}'`
maxtime=`echo "$laten" | awk -F' ' 'BEGIN {MAX=0}; {if($2>MAX){MAX=$2}}; END{printf "%.2f", MAX}'`
mintime=`echo "$laten" | awk -F' ' -v max="$maxtime" 'BEGIN {MIN=max}; {if($2<MIN){MIN=$2}}; END{printf "%.2f", MIN}'`
imbalan=`bc <<<"scale=4; 1.0 - ($avgtime/$maxtime)" | awk '{printf "%08f\n", $0}'`

##
## generate new coefficients
##
rm -f $filecor
sum=0
for (( i=0; i<$npart; i++ )); do
   c[$i]=`bc <<<"scale=8; ($avgtime/${atime[$i]})*${acoef[$i]}"`
   sum=`bc <<<"scale=8; $sum + ${c[$i]}"`
done

for (( i=0; i<$npart; i++ )); do
   c[$i]=`bc <<<"scale=8; ${c[$i]}*$npart/$sum" | awk '{printf "%08f\n", $0}'`
   j=`bc <<< "$i+1"`
   echo "$j ${c[$i]}" >> $filecor
done

##
## files copies
##

cp $filecor $alyafile 
if (( $iter \< 10 )); then
   cp $1.log $dir/$1-$phase"0"$iter.log
   cp $1-partition.par.post.res $dir/$1-"0"$iter-partition.par.post.res
else
   cp $1.log $dir/$1-$phase$iter.log
   cp $1-partition.par.post.res $dir/$1-$iter-partition.par.post.res
fi



## output info
echo "-----------------------------------------------------------------------------------"
echo "Case: $1"
echo -e "Number lines: \t\t$size"
echo -e "Sum correct : \t\t$sum"
echo -e "Average time: \t\t$avgtime"
echo -e "Maximum time: \t\t$maxtime"
echo -e "Minimum time: \t\t$mintime"
echo -e "Imbalance:    \t\t$imbalan"
echo "-----------------------------------------------------------------------------------"
