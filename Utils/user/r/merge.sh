#!/bin/bash
output=$1
shift
cat $@ > $output.csv 
