#!/bin/bash

FILE=test.cvs
RES_FILE=res.dat

echo "============================="
echo "=    Parsing file $FILE     ="
echo "============================="

SUM=0
AVG=0
COUNTER=0

# Number of time you launched the kernel (see the ruby file used)
MAX_COUNTER=100

COUNTER_USAGE=1
# Uncomment if you compare usage included, inlined and call
#MAX_USAGE=3
# Uncomment if you compare ref and boast kernel
MAX_USAGE=2

NEST=1
# Uncomment if you test the kernels with all nests in the main procedure
#MAX_NEST=1
# Uncomment if you test the kernels with one nest at a time in the main procedure
MAX_NEST=13


# Uncomment if you compare usage included, inlined and call
#echo "#NEST\t Included\t Inlined\t call" >> $RES_FILE
# Uncomment if you compare ref and boast kernel
echo "#NEST\t Ref\t Boast" >> $RES_FILE


/bin/echo -n "$NEST " >> $RES_FILE

echo "** NEST = $NEST **"

while read Line; do

	#DEBUG: echo "LINE = $Line"
  printf -v Line "%.7f" "$Line"
	#DEBUG: echo "NOW LINE = $Line"

	NUMBER=$( echo "$Line")
	SUM=$( echo $SUM + $NUMBER | bc -l )
	COUNTER="$(($COUNTER+1))"

	if [[ $COUNTER -eq $MAX_COUNTER ]]; then
		echo "COUNTER = $COUNTER"
		echo "-> SUM = $SUM"
		#AVG=$( echo $SUM/$MAX_COUNTER | bc -l)
		AVG=$( echo $SUM/$MAX_COUNTER | awk '{printf("%.7f\n", $1);}')
		echo "-> AVG = $AVG"
		/bin/echo -n "$AVG " >> $RES_FILE

		COUNTER=0
		#DEBUG: echo "COUNTER_USAGE = $COUNTER_USAGE"
    COUNTER_USAGE="$(($COUNTER_USAGE+1))"
		SUM=0
		AVG=0
 fi

 if [[ $COUNTER_USAGE -gt $MAX_USAGE ]]; then
	#DEBUG: echo "COUNTER_USAGE = $COUNTER_USAGE"
	echo  " " >> $RES_FILE
	COUNTER_USAGE=1
  NEST="$(($NEST+1))"
	#DEBUG: echo "* NEST = $NEST *"
	/bin/echo -n "$NEST " >> $RES_FILE
 fi

	if [[ $NEST -gt $MAX_NEST ]]; then
		#echo "ERROR !!!!"
		exit 1
	fi
	
done < $FILE

echo "Done"


