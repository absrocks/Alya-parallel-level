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


## check imputs
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

## prepare files
alyafile="rank-elements.dat"
filein="$case-partition.par.post.res"
itant=$(($iter-1))
itaant=$(($iter-2))
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

## take data from filein
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

## get imblalance
echo "$ratio" | awk -F' ' 'BEGIN{t=1} {printf "%8s%14s\n", t, $2; t+=1}' > $fileimb
size=`wc -l $fileimb  | cut -f1 -d' '`

if (( $size != $npart )); then 
   echo "   ERROR: 6th argument wrong size $size $npart"
   exit 1 
fi
avgtime=`echo "$laten" | awk -F' ' 'BEGIN {SUM=0; t=0}; {SUM=SUM+$2;t+=1}; END{SUM=SUM/t; printf "%.2f", SUM}'`
maxtime=`echo "$laten" | awk -F' ' 'BEGIN {MAX=0}; {if($2>MAX){MAX=$2}}; END{printf "%.2f", MAX}'`
avgrati=`echo "$ratio" | awk -F' ' 'BEGIN {SUM=0; t=0}; {SUM=SUM+$2;t+=1}; END{SUM=SUM/t; printf "%.2f", SUM}'`
maxrati=`echo "$ratio" | awk -F' ' 'BEGIN {MAX=0}; {if($2>MAX){MAX=$2}}; END{printf "%.2f", MAX}'`
minrati=`echo "$ratio" | awk -F' ' 'BEGIN {MIN=1}; {if($2<MIN){MIN=$2}}; END{printf "%.2f", MIN}'`

#echo "values: $minrati $maxrati $avgrati"
#maxbound="$( bc <<<"1.0 + ($maxrati-1.0)*1.0" )"
#minbound="$( bc <<<"1.0 + ($minrati-1.0)*1.0" )"

#echo "BOUNDS: $minbound $maxbound"
#if (( $(echo "$maxbound < 1.01" | bc -l) )); then
#  maxbound=1.01
#fi
#if (( $(echo "$minbound > 0.99" | bc -l) )); then
#  minbound=0.99
#fi
#echo "BOUNDS CORRECTED: $minbound $maxbound"

maxbound=1.10
minbound=0.90

## generate new coefficients
if [ $relx != "10" ]; then

   for (( i=1; i<=$npart; i++ )); do relxa[$i]=$relx; done
   suma=`paste $fileimb $filecorant | sed -e 's/^[ ]*//' | sed 's/  \+/ /g' | awk -F' ' -v relxa="${relxa[*]}" -v maxbound="$maxbound" -v minbound="$minbound" 'BEGIN {split(relxa,arr, / /); SUM=0; t=1}{
        if( $2 < maxbound && $2 > minbound )
	   SUM=SUM+$4
	else
	   SUM=SUM+($4*(1+((1/$2-1)*arr[t])))
        t+=1} 
        END {printf "%s\n", SUM}'`
   sumd=`paste $fileimb $filecorant | sed -e 's/^[ ]*//' | sed 's/  \+/ /g' | awk -F' ' -v relxa="${relxa[*]}" -v maxbound="$maxbound" -v minbound="$minbound" 'BEGIN {split(relxa,arr, / /); SUM1=0; SUM2=0; t=1}{
        if( $2 < maxbound && $2 > minbound )
	   SUM1=SUM1+$4
	else
	   SUM2=SUM2+($4*(1+((1/$2-1)*arr[t])))
        t+=1} 
        END {printf "%s, %s\n", SUM1,SUM2}'`
   sumd1=`echo $sumd | cut -f1 -d','`
   sumd2=`echo $sumd | cut -f2 -d','`
   ##echo $suma "=" $sumd1 "+" "$sumd2"
   sized="$( bc <<<"$size - $sumd1" )"
   ##echo "hola" $sized
   paste $fileimb $filecorant | awk -F' ' -v suma="$sumd2" -v size="$sized" -v relxa="${relxa[*]}" -v maxbound="$maxbound" -v minbound="$minbound" 'BEGIN{split(relxa,arr, / /); t=1}{
        if( $2 < maxbound && $2 > minbound )
	   printf "%8s%14s\n", t, $4
	else
           printf "%8s%14s\n", t, (size/suma)*$4*(1+(1/$2-1)*arr[t]) 
        t+=1}' > $filecor
   sumb=`cat $filecor | sed -e 's/^[ ]*//' | sed 's/  \+/ /g' | awk -F' ' 'BEGIN {SUM=0}; {SUM=SUM+($2)}; END {printf "%s\n", SUM}'`

else
   if (( $itaant \< 10 )); then
      filecoraant="$dir/$1-$phase$part"0"$itaant.cor"
   else
      filecoraant="$dir/$1-$phase$part$itaant.cor"
   fi
   if (( $iter < 2 )); then
      for (( i=1; i<=$npart; i++ )); do relxa[$i]=1.0; done
      suma=`paste $fileimb $filecorant | sed -e 's/^[ ]*//' | sed 's/  \+/ /g' | awk -F' ' -v relxa="${relxa[*]}" 'BEGIN {split(relxa,arr, / /); SUM=0; t=1}
         {SUM=SUM+($4*(1+((1/$2-1)*arr[t]))); t+=1}; END {printf "%s\n", SUM}'`
      paste $fileimb $filecorant | awk -F' ' -v suma="$suma" -v size="$size" -v relxa="${relxa[*]}" 'BEGIN{split(relxa,arr, / /); t=1} 
         {printf "%8s%14s\n", t, (size/suma)*$4*(1+(1/$2-1)*arr[t]); t+=1}' > $filecor
      sumb=`cat $filecor | sed -e 's/^[ ]*//' | sed 's/  \+/ /g' | awk -F' ' 'BEGIN {SUM=0}; {SUM=SUM+($2)}; END {printf "%s\n", SUM}'`
   else

      echo "here: $fileimb $fileimbant $filecorant $filecoraant"
      relxa=($( paste $fileimb $fileimbant $filecorant $filecoraant | awk '{ 
      a1= $2-$4
      a2= $6-$8 
      if( a2 < 0.001 || a2 < 0.001 )
         c = $6
      else
	a = a1/a2 
      	b = $2-a$6
        c= (1-b)/a
      if( c < 0.5 ) c = 0.5
      if( c > 1.5 ) c = 1.5
      print c
      }' ))
      suma=`paste $fileimb $filecorant | sed -e 's/^[ ]*//' | sed 's/  \+/ /g' | awk -F' ' -v relxa="${relxa[*]}" 'BEGIN {split(relxa,arr, / /); SUM=0; t=1}
         {SUM=SUM+(arr[t]); t+=1}; END {printf "%s\n", SUM}'`
      paste $fileimb $filecorant | awk -F' ' -v suma="$suma" -v size="$size" -v relxa="${relxa[*]}" 'BEGIN{split(relxa,arr, / /); t=1} 
         {printf "%8s%14s\n", t, (size/suma)*(arr[t]); t+=1}' > $filecor
      sumb=`cat $filecor | sed -e 's/^[ ]*//' | sed 's/  \+/ /g' | awk -F' ' 'BEGIN {SUM=0}; {SUM=SUM+($2)}; END {printf "%s\n", SUM}'`
      
   fi
fi

## last file copies
cp $filecor $alyafile 
if (( $iter \< 10 )); then
   cp $1.log $dir/$1-$phase"0"$iter.log
   cp $1-partition.par.post.res $dir/$1-"0"$iter-partition.par.post.res
else
   cp $1.log $dir/$1-$phase$iter.log
   cp $1-partition.par.post.res $dir/$1-$iter-partition.par.post.res
fi


imbalan=`bc <<<"scale=4; 1.0 - ($avgtime/$maxtime)"`

## output info
echo "-----------------------------------------------------------------------------------"
echo "Case: $1"
echo -e "Number of lines: \t\t$size"
echo -e "Sum initial correction : \t$suma"
echo -e "Sum normalized correct : \t$sumb"
echo -e "Average time: \t\t\t$avgtime"
echo -e "Maximum time: \t\t\t$maxtime"
echo -e "Imbalance:    \t\t\t$imbalan"
echo "-----------------------------------------------------------------------------------"
