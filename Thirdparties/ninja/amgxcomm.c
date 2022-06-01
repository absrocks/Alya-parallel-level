#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include "amgx_c.h"

typedef size_t devptr_t;

void amgx_communicator_init_(MPI_Fint *f_handle)
{
  MPI_Comm comm;
  comm = MPI_Comm_f2c(*f_handle);
  int rank;
  MPI_Comm_rank(comm,&rank);
  printf("my rank in sin master is %d\n\n",rank);  
  return;
}

void amgx_create_resource_(devptr_t *resource, devptr_t *cfg,MPI_Fint *f_handle,int *rank)
{
  AMGX_resources_handle rsrc;
  MPI_Comm comm;
  int N,lrank;
  
  if(*rank == -1)
    {
      AMGX_SAFE_CALL(AMGX_resources_create_simple(&rsrc,(AMGX_config_handle)(*cfg)));
      *resource = (devptr_t)rsrc;
    }
  else
    {
      cudaGetDeviceCount(&N);
      lrank = (*rank) % N;
      comm = MPI_Comm_f2c(*f_handle);
      AMGX_SAFE_CALL(AMGX_resources_create(&rsrc,(AMGX_config_handle)(*cfg),comm,1,&lrank));
    }
  return; 
}
