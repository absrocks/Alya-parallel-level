#!/bin/bash

FILE=test.cvs
RES_FILE=res.dat

echo "=> Parsing file $FILE"


printf "#Options\t Time\n" >> $RES_FILE

while read Line; do

	#/bin/echo -n $Line | cut -d',' -f1,2,3,4,5 >> $RES_FILE
	OPTIONS=`echo $Line | cut -d',' -f1,2,3,4,5`
	printf "\"$OPTIONS\" " >> $RES_FILE
	TIME=$(echo $Line | cut -d',' -f6) 
	printf -v TIME "%.7f" "$TIME"
	printf "\t $TIME\n" >> $RES_FILE

done < $FILE

echo "=> Done"


