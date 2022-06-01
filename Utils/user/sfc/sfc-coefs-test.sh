#!/bin/bash

# chec parameters
if [ "$#" -ne 6 ]; then 
    echo "Parameters it must be 6: [case name],[iteration],[folder tag],[relax factor/10],[spmv/assembly][nparts]"
    exit 1
fi 

## inputs
case=$1
iter=$(( $2 -1 ))
flag=$3
relx=$4
phase=$5
npart=$6

##
## check imputs
##

if (( $iter < 0 )); then 
  echo "   ERROR: 2th argument must be > 0 "
   exit 1 
fi
if [ $phase != "assembly" ] && [ $phase != "spmv" ]; then 
   echo "   ERROR: 4th argument must be <spmv> or <assembly> "
   exit 1 
fi
if (( $npart < 2 )); then 
   echo "   ERROR: 6th argument must be > 1 "
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

if [ $phase == "spmv" ]; then
   ratio=`cat $filein | awk '/ComponentNames\ SPMV_RATIO_COMPUTATION_CONTINUITY/,/End\ Values/' | sed -e '1,3d' | sed -e '$d' | sed -e 's/^[ ]*//' | sed 's/  \+/ /g'`
   laten=`cat $filein | awk '/ComponentNames\ SPMV_TIME_COMPUTATION_CONTINUITY/,/End\ Values/' | sed -e '1,3d' | sed -e '$d' | sed -e 's/^[ ]*//' | sed 's/  \+/ /g'`
   entrs=`cat $filein | awk '/ComponentNames\ MATRIX_NUM_COEFF_CONTINUITY/,/End\ Values/' | sed -e '1,3d' | sed -e '$d' | sed -e 's/^[ ]*//' | sed 's/  \+/ /g'`
   avgent=`echo "$entrs" | awk -F' ' 'BEGIN {SUM=0; t=1}; {SUM=SUM+$2;t+=1}; END{SUM=SUM/t; printf "%.2f", SUM}'`
   echo "$entrs" | awk -F' ' -v avgent="$avgent" 'BEGIN{t=1} {printf "%8s%14f\n", t, $2/avgent; t+=1}' > $fileent
else 
   ratio=`cat $filein | awk '/ComponentNames\ ELEMENT_ASSEMBLY_RATIO_NASTIN/,/End\ Values/' | sed -e '1,3d' | sed -e '$d' | sed -e 's/^[ ]*//' | sed 's/  \+/ /g'`
   laten=`cat $filein | awk '/ComponentNames\ ELEMENT_ASSEMBLY_TIME_NASTIN/,/End\ Values/' | sed -e '1,3d' | sed -e '$d' | sed -e 's/^[ ]*//' | sed 's/  \+/ /g'`
fi

##
## get imblalance
##
echo "$ratio" | awk -F' ' 'BEGIN{t=1} {printf "%8s%14s\n", t, $2; t+=1}' > $fileimb
size=`wc -l $fileimb  | cut -f1 -d' '`

if (( $size != $npart )); then 
   echo "   ERROR: 6th argument wrong size $size $npart"
   exit 1 
fi


avgtime=`echo "$laten" | awk -F' ' 'BEGIN {SUM=0; t=0}; {SUM=SUM+$2;t+=1}; END{SUM=SUM/t; printf "%.8f", SUM}'`
maxtime=`echo "$laten" | awk -F' ' 'BEGIN {MAX=0}; {if($2>MAX){MAX=$2}}; END{printf "%.8f", MAX}'`
mintime=`echo "$laten" | awk -F' ' -v max="$maxtime" 'BEGIN {MIN=max}; {if($2<MIN){MIN=$2}}; END{printf "%.8f", MIN}'`
imbalan=`bc <<<"scale=4; ($avgtime/$maxtime)" | awk '{printf "%8f\n", $0}'`

maxbound=1.01
minbound=0.99

##
## generate new coefficients
##
rm $filecor
coefk[0]=0
sum=0
for (( i=1; i<=$npart; i++ )); do 
   sumk=`echo "$laten" | awk -F' ' -v k="$i" 'BEGIN {SUM=0; t=1}; {if(t<=k){SUM=SUM+$2};t+=1}; END{printf"%.8f", SUM}'`
   coefk[$i]=`bc <<<"scale=8; ($avgtime*$i)/$sumk"`
   j=`bc <<< "$i-1"`
   c[$i]=`bc <<<"scale=8; ($i*${coefk[$i]})-($j*${coefk[$j]})" | awk '{printf "%08f\n", $0}'`
   sum=`bc <<<"scale=8; $sum + ${c[$i]}"`
   echo "$i ${c[$i]}" >> $filecor
done

#
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
echo -e "Imbalance:    \t\t0$imbalan"
echo "-----------------------------------------------------------------------------------"
