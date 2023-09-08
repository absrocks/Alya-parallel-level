

ROOT_CANTERA=/home/bsc21/bsc21704/z2018_1/REPOSITORY/CANTERA/Cantera230/ExecsINTEL

mpif90 demo.f90 \
-I$ROOT_CANTERA/include/cantera  \
  $ROOT_CANTERA/lib/libcantera_fortran.a \
  $ROOT_CANTERA/lib/libcantera.a \
-lstdc++ 

export PYTHONPATH=${PYTHONPATH}:$ROOT_CANTERA/lib/python2.7/site-packages
./a.out 

# -I/apps/BOOST/1.64.0/INTEL/IMPI/include/ \

