// General
#include <iostream> // cin, cout, endl, cerr
#include <vector>   // vector
#include <map>
#include <mpi.h>
#include "commdom.hpp"
using namespace std;

string IAM = "MMMESH";

int main(int argc, char** argv)
{
  //------------------------------------------------------| INIT COUPLING |---//
  MPI_Init(NULL, NULL);
  MPI_Comm  world_comm = MPI_COMM_WORLD;

  int  app_id = -1;
  int  n_apps = -1;
  MPI_Comm  local_comm;
  string  namei = "";
  if(argc==2) namei = argv[1];

  // CommDom
  CommDom  CD = CommDom();
  CD.init();
  CD.set_app_type(IAM);
  CD.set_world_comm(world_comm);
  CD.set_app_name(namei);
  CD.name_to_id(&app_id, &n_apps, &local_comm);
  CD.__create_interaction__(); //MPI_Barrier(local_comm);
  CD.create_commij(local_comm);


  //------------------------------------------------------------| WORKING |---//
  MPI_Barrier(local_comm);
  int local_rank = -1;
  MPI_Comm_rank(local_comm, &local_rank);

  int commij_size = -1;
  CD.get_commij_size(&commij_size);

  //----------------------------------------------------------------| END |---//
  MPI_Barrier(local_comm);
  MPI_Finalize();

  //if(local_rank==0) cout<<"|\n+OK! \n\n";
  return 0;
} // main

//=======================================================================||===//
