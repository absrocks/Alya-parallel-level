#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

void insitu_communicator_init_(MPI_Fint *f_handle)
{
  MPI_Comm comm;
  comm = MPI_Comm_f2c(*f_handle);
  commtoinsitu(comm);
  return;
}
