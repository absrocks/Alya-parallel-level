#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "commdom.hpp"

//#include "shrUtils.h"
// utilities and system includes
// // CUDA-C includes
#include "cuda.h"
#define BUFSIZE 256
#define TAG 0
//  
int devCount;
int myid;
int ihavecuda;
int nodes[256];
int nocuda[256];
int deviceselector=0;
//     
string IAM = "MMMESH";

int main(int argc, char *argv[])
{
  char idstr[256];
  char idstr2[256];
  char buff[BUFSIZE];
  int i;
  int numprocs, rank, namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  freopen("/dev/null", "w", stderr); /* Hide errors from nodes with no CUDA cards */

  MPI_Status stat;
  MPI_Init(NULL, NULL);
  MPI_Comm  world_comm = MPI_COMM_WORLD;

  string  namei = "";
  if(argc==2) namei = argv[1];

  // CommDom
  CommDom  CD = CommDom();
  CD.init();
  CD.set_app_type(IAM);
  CD.set_world_comm(world_comm);
  CD.set_app_name(namei);

  int  app_id = -1;
  int  n_apps = -1;
  MPI_Comm  local_comm;
  local_comm = CD.set_mpi_comms(); 
  MPI_Barrier(local_comm);

  int local_rank = -1;
  MPI_Comm_rank(local_comm, &local_rank);

/*
  MPI_Comm_size(local_comm, &numprocs);
  MPI_Comm_rank(local_comm, &rank);
  MPI_Get_processor_name(processor_name, &namelen);
  MPI_Comm_rank(local_comm, &myid);
  if (myid == 0)
  {
    printf("  We have %d processors\n", numprocs);
    printf("  Spawning from \'%s.%s\' \n", processor_name, namei.c_str());
    printf("  CUDA MPI\n");
    printf("\n");
    for(i=1; i<numprocs;i++)
    {
      buff[0]='\0';
      MPI_Send(buff, BUFSIZE, MPI_CHAR, i, TAG, local_comm);
    }
    printf("\n\n\n");
    printf("  Probing nodes...\n");
    printf("     Node        Psid  CUDA Cards (devID)\n");
    printf("     ----------- ----- ---- ----------\n");
    for(i=1; i<numprocs;i++)
    {
      MPI_Recv(buff, BUFSIZE, MPI_CHAR, i, TAG, local_comm, &stat);
      printf("%s\n", buff);
    }
    printf("\n");
    MPI_Finalize(); 
  }
  else
  {
   MPI_Recv(buff, BUFSIZE, MPI_CHAR, 0, TAG, local_comm, &stat);
   MPI_Get_processor_name(processor_name, &namelen);
   cudaGetDeviceCount(&devCount);
   buff[0]='\0';
   idstr[0]='\0';
   if (devCount == 0) 
   {
     sprintf(idstr,"- %-11s %5d %4d NONE", processor_name, rank, devCount);
     ihavecuda=0;
   }
   else
   {
     ihavecuda=1;
     if (devCount >= 1)
     {
       sprintf(idstr, "+ '%s'.%-11s %5d %4d", namei.c_str(), processor_name, rank, devCount);
       idstr2[0]='\0';
       for (int i = 0; i < devCount; ++i)
       {
         cudaDeviceProp devProp;
         cudaGetDeviceProperties(&devProp, i);
         sprintf(idstr2, " %s (%d) ", devProp.name, i);
         strncat(idstr,idstr2,BUFSIZE);
       }
     }
     else
     {
       cudaDeviceProp devProp;
       cudaGetDeviceProperties(&devProp, i);
       sprintf(idstr, "%-11s %5d %4d %s", processor_name, rank, devCount, devProp.name);
     } 
   } 
   strncat(buff, idstr, BUFSIZE);
   MPI_Send(buff, BUFSIZE, MPI_CHAR, 0, TAG, local_comm);
  }
*/

  int n_send = 1;
  int n_recv = 1;
  double send[1];
  double recv[1];

  send[0] =  1e20;
  recv[0] = -1e20;

  string  namej = "";
  if(namei=="SOLID")
  {
    namej   = "FLUID";
    send[0] = 1.0; 
  }
  if(namei=="FLUID")
  {
    namej   = "SOLID";
    send[0] = -1.0;  
  }

  MPI_Comm  commij;
  CD.get_commij(namej, &commij);
  CD.__mpi_sendrecv_real__(send, n_send, recv, n_recv, local_comm, commij); 

  MPI_Barrier(commij);
  if(local_rank==0) cout<<" ===>'"<< namei <<"': "<< recv[0] <<"<=== \n"; 

  CD.locator_create2(local_comm, commij, 1e-3); 

/*
  int        n_vertices_j = 0 
  double *vertex_coords_j = NULL; 

  int     n_vertices_i    = 0 
  int     n_elements_i    = 0
  double *vertex_coords_i = NULL; 
  int    *vertex_num_i    = NULL;  

  CD.locator_set_mesh(n_vertices_i, 
                      n_elements_i, 
                      vertex_coords_i, 
                      vertex_num_i,
                      n_vertices_j, 
                      vertex_coords_j)

  CD.save_dist_coords(local_rank); 
*/  

  MPI_Finalize();
  return 0;
} 

