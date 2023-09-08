#!/bin/bash

export EXTRAE_HOME=/apps/BSCTOOLS/extrae/3.6.1/impi_2017_4
export EXTRAE_CONFIG_FILE=./extrae.xml
export LD_PRELOAD=${EXTRAE_HOME}/lib/libmpitracef.so # For Fortran apps

## Run the desired program
$*

