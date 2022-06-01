// General
#include <iostream> // cin, cout, endl, cerr
#include <vector>   // vector
#include <map>
#include <mpi.h>
#include "commdom.hpp"
#include "read_file.hpp"

using namespace std;

string IAM = "PLEPP";

int main(int argc, char** argv)
{
  int idx = -1; 

  string  name_argv = "";
  if(argc==2) name_argv = argv[1];

  //--------------------------------------------------------------------||---//

  string         namei  = ":(";
  vector<string> Files;
  map<string,MPI_Comm>            Commij;
  map<string,MPI_Comm>::iterator  It;

  if( name_argv == "Mesh01")
  {
    namei            = name_argv; 
    Commij["Mesh02"] = MPI_COMM_NULL;
    Commij["Mesh03"] = MPI_COMM_NULL;
    Files.push_back("/Users/poderozo/z2015/Runner/TRIPLE01/MESH01/Mesh01_COORDINATES.alya"); 
    Files.push_back("/Users/poderozo/z2015/Runner/TRIPLE01/MESH01/Mesh01_ELEMENTS.alya"); 
  }

  if( name_argv == "Mesh02")
  {
    namei            = name_argv; 
    Commij["Mesh01"] = MPI_COMM_NULL;
    Commij["Mesh03"] = MPI_COMM_NULL;
    Files.push_back("/Users/poderozo/z2015/Runner/TRIPLE01/MESH02/Mesh02_COORDINATES.alya"); 
    Files.push_back("/Users/poderozo/z2015/Runner/TRIPLE01/MESH02/Mesh02_ELEMENTS.alya"); 
  }

  if( name_argv == "Mesh03")
  {
    namei            = name_argv; 
    Commij["Mesh01"] = MPI_COMM_NULL;
    Commij["Mesh02"] = MPI_COMM_NULL;
    Files.push_back("/Users/poderozo/z2015/Runner/TRIPLE01/MESH03/Mesh03_COORDINATES.alya"); 
    Files.push_back("/Users/poderozo/z2015/Runner/TRIPLE01/MESH03/Mesh03_ELEMENTS.alya"); 
  }
  /*
   --- ---
  |   | 3 |
  | 1 |---
  |   | 2 |
   --- ---
  */
    
  //--------------------------------------------------------------| PLEPP |---//
  MPI_Init(NULL, NULL);
  MPI_Comm global_comm = MPI_COMM_WORLD;

  int  app_id = -1;
  int  n_apps = -1;
  MPI_Comm  local_comm;

  // CommDom
  CommDom  CD = CommDom();
  CD.init();
  CD.set_app_type(IAM);
  CD.set_world_comm(global_comm);
  CD.set_app_name(namei);
  CD.name_to_id(&app_id, &n_apps, &local_comm);
  CD.__create_interaction__();   
  CD.create_commij(local_comm);

  int local_rank = -1;
  MPI_Comm_rank(local_comm, &local_rank);


  //----------------------------------------------------------------| PTS |---//
  MPI_Barrier(global_comm);
  vector<read_log_file> Data(2);

  idx = 0; 
  Data[idx].set_name( Files[idx] );
  Data[idx].run();
  vector< vector<double> > vcoords( Data[idx].get_vdata() );
  Data[idx].end();

  int n_vertices_i = vcoords.size(); 
  if(local_rank==0) cout<< namei <<".n_vertices = "<< n_vertices_i <<"\n";

  vector<double> vertex_coords_i;
  for(int i=0; i<n_vertices_i; i++) 
  {
    for(int j=1; j<vcoords[i].size(); j++) 
    {
      vertex_coords_i.push_back( vcoords[i][j] );
    }
  }

  // vector -> ptr
  double* vertex_coords_i_ptr = NULL;
  vertex_coords_i_ptr         = new double[ vertex_coords_i.size() ];
  memcpy(vertex_coords_i_ptr, vertex_coords_i.data(), sizeof(double) * vertex_coords_i.size() );

  vertex_coords_i.clear();   
  vcoords.clear();

  int            n_vertices_j = n_vertices_i;
  double* vertex_coords_j_ptr = vertex_coords_i_ptr;


  //--------------------------------------------------------------| CELLs |---//
  idx++; 
  Data[idx].set_name( Files[idx] );
  Data[idx].run();
  vector< vector<double> > vcells( Data[idx].get_vdata() );
  Data[idx].end();
  
  int n_elements_i = vcells.size(); 
  if(local_rank==0) cout<< namei <<".n_elements_i  = "<< n_elements_i <<"\n";

  vector<int> vertex_num_i;
  vector<int> vertex_type_i;
  for(int i=0; i<n_elements_i; i++) 
  {
    for(int j=1; j<vcells[i].size(); j++) 
    {
       vertex_num_i.push_back( (int)vcells[i][j]-0 ); // Fortran style (not C!!). 
      vertex_type_i.push_back( -1 ); 
    }
  }

  // vector -> ptr
  int* vertex_num_i_ptr = NULL;
  vertex_num_i_ptr      = new int[ vertex_num_i.size() ];
  memcpy(vertex_num_i_ptr, vertex_num_i.data(), sizeof(int) * vertex_num_i.size() );

  vertex_num_i.clear();   
  vcells.clear();

  int* vertex_num_j_ptr  = vertex_num_i_ptr;
  

  //-----------------------------------------------------------| GET COMMs|---//
  MPI_Barrier(global_comm);
  
  for(It = Commij.begin(); It != Commij.end(); ++It)
  {
    string    name_j = It->first;
    MPI_Comm  commij = It->second;
    
    if(CD.__get_friends__(name_j) == 1) 
    {
      commij         = CD.get_mpi_commij(name_j); 
      Commij[name_j] = commij; 
    }
  }

//  It = Commij.find("b");
//  if(It != Commij.end())  Commij.erase( It );


  //-----------------------------------------------------------| LOCATION |---//
  MPI_Barrier(global_comm);

  for(It = Commij.begin(); It != Commij.end(); ++It)
  {
    string    name_j = It->first;
    MPI_Comm  commij = It->second;
    
    if(commij != MPI_COMM_NULL) 
    {
      CD.locator_create2(local_comm, commij, 1e-3, name_j);

      CD.locator_set_mesh(n_vertices_i, 
                          n_elements_i, 
                          vertex_coords_i_ptr, 
                          vertex_num_i_ptr,
                          n_vertices_j, 
                          vertex_coords_j_ptr, 
                          name_j);

      CD.save_dist_coords(local_comm, name_j);

      int n_recv = CD.get_n_interior(name_j); 
      int n_send = CD.get_n_dist_points(name_j); 

      //CD.locator_destroy();
    }
  }


  //--------------------------------------------------------------------||---//
  MPI_Barrier(global_comm);




  //--------------------------------------------------------------------||---//
  //--------------------------------------------------------------------||---//


  //--------------------------------------------------------------------||---//
  MPI_Finalize();
  //if(local_rank==0) cout<<"|\n+OK! \n\n";
  return 0;
} // main

//=======================================================================||===//