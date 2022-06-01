#!/bin/bash
for arg in "$@"; do
    sed -i "s/E\([+-]\)\([0-9][0-9]\)\>/E\10\2/g" $arg
done
