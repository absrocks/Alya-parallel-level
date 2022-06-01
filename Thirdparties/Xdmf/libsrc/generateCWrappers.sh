#!/bin/bash

Xdmf_SOURCE_DIR=/home/bsc21/bsc21399/Xdmf
Xdmf_BINARY_DIR=/home/bsc21/bsc21399/Xdmf/Linux
ICE_INCLUDES="-I${Xdmf_SOURCE_DIR} -I${Xdmf_SOURCE_DIR}/libsrc -I${Xdmf_BINARY_DIR}/libsrc -I${Xdmf_BINARY_DIR}/Ice/libsrc"

echo "$ICE_INCLUDES"
swig -v -c++ -includeall -c $ICE_INCLUDES -o XdmfC.cxx Xdmf.i

