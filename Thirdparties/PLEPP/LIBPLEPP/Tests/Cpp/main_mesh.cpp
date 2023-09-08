// General 
#include <iostream> // cin, cout, endl, cerr
#include <vector>   // vector 
#include <map>
#include <mpi.h>
using namespace std;

#include "loader_alya.hpp"
#include "commdom.hpp"

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

  if(commij_size>0)
  {
    //if(local_rank==0) CD.print_commij_names();
    string      namej = ""; 
    MPI_Comm  commij; 

    namej = "XXXX"; // WhoIsSending?
    CD.get_commij(namej, &commij); 
    if(commij != MPI_COMM_NULL)
    {
      CD.__print_comm__(commij); //cout<<"\n";

      int sval = 69; 
      if(local_rank==0) CD.mpi_tbicpl(&sval, 1, NULL, 0, namej); 
      //cout<<"|==["<< namei <<"]: "<< sval;
      cout<<"\n"; 
      cout<<"\n"; 
      //MPI_Barrier(commij);

      
      // GET MESH 
      int appi = 0; 
      int appj = 1; 
      vector<loader_alya> GetMesh(2);
      if(local_rank == 0)
      {
        //GetMesh[appi].set_ple_data("/home/jmake/z2014/Cplng/MeshCplng/Mesh02/mesh02", 1);
        //GetMesh[appj].set_ple_data("/home/jmake/z2014/Cplng/MeshCplng/Mesh01/mesh01", 2);
        GetMesh[appi].set_ple_data("/home/jmake/z2014/Cplng/MeshCplng/MeshCoarse/disk", 1); // n_Vrtx=47, n_elems=74
        GetMesh[appj].set_ple_data("/home/jmake/z2014/Cplng/MeshCplng/MeshCoarse/square", 1); // n_Vrtx=57, n_elems=88
        //vertex_coords_size = GetMesh[].vertex_coords.size();
        //vertex_num_size    = GetMesh[].vertex_num.size();      
      }
      MPI_Barrier(local_comm);  


      // LOCATOR
      CD.__locator_create__(namej, NULL);
      ple_locator_t  *locatorj = CD.__locator_get__(); 

      CD.__locator_init_meshes__(); 
      if(local_rank == 0)
      {
        CD.__locator_set_meshi__(GetMesh[appi].n_vertices, GetMesh[appi].n_elements, 
                                 GetMesh[appi].vertex_coords_ptr, GetMesh[appi].vertex_num_ptr); 
        CD.__locator_set_meshj__(GetMesh[appj].n_vertices, GetMesh[appj].vertex_coords_ptr); 
      }
      CD.__locator_set_mesh__(locatorj); 
      MPI_Barrier(local_comm);  

      
      // INTERPOLATION 
      MPI_Barrier(local_comm);
      double *var_ij = NULL; // to be sent
      double *var_ji = NULL; // to be received
      int stride    = 1;
      int reverse   = 0;

      int n_dist_points = ple_locator_get_n_dist_points(locatorj);
      var_ij = new double[stride*n_dist_points];
      for(int i=0; i<stride*n_dist_points; i++) var_ij[i] = (i+1.123);

      //int n_interior    = ple_locator_get_n_interior(locatorj);
      //var_ji = new double[stride*n_interior];
      //for(int i=0; i<stride*n_interior;    i++) var_ji[i] = -666.66; 

      CD.__locator_exchange_point_var_send_double__(locatorj, var_ij, stride, reverse); 


      ple_locator_destroy(locatorj);

    }
    
    MPI_Barrier(local_comm);  
  }

  //----------------------------------------------------------------| END |---// 
  MPI_Barrier(local_comm);
  MPI_Finalize();
  
  //if(local_rank==0) cout<<"|\n+OK! \n\n";
  return 0; 
} // main 

//=======================================================================||===//
