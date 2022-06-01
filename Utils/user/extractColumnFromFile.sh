#!/bin/bash
# usage: ./extractColumnFromFile.sh fileName col1 col2 ...
FILE=$1
for ((i=2; i <= $#; i++)); do
	COL=${!i}
	OUT=$FILE"-col"$COL".out"
	echo "Extracting column $COL of $FILE to $OUT"
	awk "/^[^#]/ { print $"$COL" }" $FILE > $OUT
	sed -i "s/E\([+-]\)\([0-9][0-9]\)\>/E\10\2/g" $OUT
done
