/*
  swig -python -c++ xxx.i   
*/
%module Read_file 
%include "std_string.i"


%include "typemaps.i"
//%apply int *OUTPUT { int *result01 };

%include "carrays.i"
%array_class(float,  farray);
%array_class(double, darray);
%array_class(int,    iarray);

%{
  #include "../../Include/read_file.hpp"
%}


%include mpi4py.i
%mpi4py_typemap(Comm, MPI_Comm);


//%include "numpy.i"

%include "../../Include/read_file.hpp"

/*
void sayhello(MPI_Comm comm);
int  alya_coupling_mpi_name_to_id(MPI_Comm, char*);
void alya_coupling_set_create(char*, char*, MPI_Comm, MPI_Comm, int, int[2], int* result); 

MPI_Comm 
alya_coupling_mpi_intracomm_create(MPI_Comm, MPI_Comm, int[2], int[2]); 

void alya_locator_create(double,  MPI_Comm, int, int); 
void alya_locator_set_mesh02(MPI_Comm, int, int, int, int, double*, int*, int, int); 
//void alya_locator_set_mesh(MPI_Comm, int, int, void*); 
//void alya_coupling_mpi_intracomm_create(MPI_Comm); 
*/
