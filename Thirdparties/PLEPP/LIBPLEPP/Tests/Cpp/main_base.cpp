// General 
#include <iostream> // cin, cout, endl, cerr
#include <vector>   // vector 
#include <map>
#include <mpi.h>
using namespace std;

#include "loader_alya.hpp"
#include "commdom.hpp"


double tolerance  = 1e-3;
int*   point_list = NULL;


string          IAM = "BBBASE"; 
vector<string>  FRIENDS;

int main(int argc, char** argv)
{
  //------------------------------------------------------| SET COUPLING |---// 
  FRIENDS.push_back("PLOT"); 
  FRIENDS.push_back("XXXX"); 


  //------------------------------------------------------| INIT COUPLING |---// 
  CommDom CD = CommDom(); 
  CD.set_app_type(IAM); 

  MPI_Init(NULL, NULL);  
  MPI_Comm  world_comm = MPI_COMM_WORLD;
  CD.set_world_comm(world_comm); 

  string  namei = ""; 
  if(argc==2) namei = argv[1]; 
  CD.set_app_name(namei); 

  int  app_id = -1; 
  int  n_apps = -1; 
  MPI_Comm  local_comm;
  CD.name_to_id(&app_id, &n_apps, &local_comm); 


  //---------------------------------------------------| INIT INTERACTION |---//   
  //for(int i=0; i<HOMELY.size(); i++) CD.set_applications( HOMELY[i] ); 
  for(int i=0; i<FRIENDS.size(); i++) CD.set_iteration( FRIENDS[i] ); 

  CD.__create_interaction__();
  MPI_Barrier(local_comm); 

  CD.create_commij(local_comm); 


  //------------------------------------------------------------| WORKING |---// 
  MPI_Barrier(world_comm); 
  int local_rank = -1; 
  MPI_Comm_rank(local_comm, &local_rank); 

  int commij_size = -1; 
  CD.get_commij_size(&commij_size); 

  if(commij_size>0)
  {
    //if(local_rank==0) CD.print_commij_names();
    string      appj = ""; // WhoIsSending?
    MPI_Comm  commij; 

    appj = "PLOT"; // WhoIsSending?
    CD.get_commij(appj, &commij); 
    if(commij != MPI_COMM_NULL)
    {
      CD.__print_comm__(commij); 
      //cout<<"\n";


      if(local_rank==0)
      {
        int rval = -1;
        //CD.mpi_tbicpl(NULL, 0, &rval, 1, -1, MPI_COMM_NULL, appj); 
        //if(local_rank==0) 
        //CD.mpi_Irecv(NULL, 0, &rval, 1, -1, MPI_COMM_NULL, appj); 
        cout<<"|=>["<< namei <<"]: "<< rval;  
        cout<<"\n"; 
      }
    }
  }

/*
  if(commij_size>0)
  {
    string    mesh_app = ""; // WhoIsSending?
    MPI_Comm  mesh_comm; 

    mesh_app = "PLOT"; // WhoIsSending?
    CD.get_commij(mesh_app, &mesh_comm); 
    if(mesh_comm != MPI_COMM_NULL)
    {
      CD.__print_comm__(mesh_comm); 
      cout<<"\n";

      int rval = -1;
      if(local_rank==0) CD.mpi_tbicpl(NULL, 0, &rval, 1, -1, MPI_COMM_NULL, mesh_app); 
      //cout<<"|==["<< namei <<"]: "<< rval;  
      cout<<"\n"; 
    }
    mesh_comm = MPI_COMM_NULL; 
*/
/*
    mesh_app = "YYYY"; // WhoIsSending?
    CD.get_commij(mesh_app, &mesh_comm);     
    if(mesh_comm != MPI_COMM_NULL)
    {
      CD.__print_comm__(mesh_comm); 
      //cout<<"\n";

      int rval = -1;
      if(local_rank==0) CD.mpi_tbicpl(NULL, 0, &rval, 1, -1, MPI_COMM_NULL, mesh_app); 
      //cout<<"|==["<< namei <<"]: "<< rval;  
      cout<<"\n"; 
    }
    mesh_comm = MPI_COMM_NULL; 
*/
//  }


  //----------------------------------------------------------------| END |---// 
  MPI_Barrier(local_comm);
  MPI_Finalize();
  
  //if(local_rank==0) cout<<"|\n+OK! \n\n";
  return 0; 
} // main 

//=======================================================================||===//
