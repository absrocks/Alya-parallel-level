/*
  swig -python -c++ xxx.i   
*/
%module Loader_alya 
%include "std_string.i"


%include "typemaps.i"
//%apply int *OUTPUT { int *result01 };

%include "carrays.i"
%array_class(float,  farray);
%array_class(double, darray);
%array_class(int,    iarray);

%{
  #include "../../Include/loader_alya.hpp"
%}


%include mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);


//%include "numpy.i"

%include "../../Include/loader_alya.hpp"

