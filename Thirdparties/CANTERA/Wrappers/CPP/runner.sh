##include ../Cantera230/ExecsINTEL/include/cantera/Cantera.mak    

export PYTHONPATH=${PYTHONPATH}:/home/bsc21/bsc21704/z2018_1/REPOSITORY/CANTERA/Cantera230/ExecsINTEL/lib/python2.7/site-packages

## c++ 
icpc combustor.cpp \
-I/home/bsc21/bsc21704/z2018_1/REPOSITORY/CANTERA/Cantera230/ExecsINTEL/include \
-I/apps/BOOST/1.64.0/INTEL/IMPI/include/ \
/home/bsc21/bsc21704/z2018_1/REPOSITORY/CANTERA/Cantera230/ExecsINTEL/lib/libcantera.a \
-pthread -O3 -std=c++0x 

./a.out

