/*
  2016JAN22. BSC, BARCELONA, SPAIN 
  Migue Zavala 

  FROM:
    test04.cpp 
*/
#include <iostream> // cin, cout, endl, cerr
#include <vector>   // vector
#include <map>
#include <mpi.h>
#include <omp.h>
#include <cmath>
#include <assert.h>     /* assert */
#include <algorithm>    // std::fill
#include "../mui/mui.h"
#include "read_file.hpp"


using namespace std;

string IAM   = "PLEPP";
string NAMEi = "MESH01"; 
string NAMEj = "MESH02"; 

int main(int argc, char** argv)
{
  int idx = -1; 
  int init_col=-1; 

  string  name_argv = "";
  if(argc==2) name_argv = argv[1];

  //-----------------------------------------------------------------| MUI |---//
  MPI_Init(NULL, NULL);
  MPI_Comm PLEPP_COMM_WORLD = MPI_COMM_WORLD;
  MPI_Comm       local_comm = MPI_COMM_NULL;

  //using namespace mui;
  local_comm = mui::mpi_split_by_app();

  // CommDom
  int local_rank = -1;
  MPI_Comm_rank(local_comm, &local_rank);

  string mui_type; 
  if(name_argv == NAMEi) mui_type="mpi://brownian/ifs"; 
  if(name_argv == NAMEj) mui_type="mpi://vortex/ifs";       
  mui::uniface3d interface( mui_type );
  //--------------------------------------------------------------------||---//

  string         namei  = ":(";
  string         namej  = ":(";
  string         dicti  = ":(";

  map<string, vector<int>    >            ArrayI;
  map<string, vector<double> >            ArrayD;

  map<string,string>                      Files;
  typedef map<string,string>::iterator    FilesIt;

  map<string,MPI_Comm>                    Commij;
  typedef map<string,MPI_Comm>::iterator  It;

  if( name_argv == NAMEi )
  {
    namei         = name_argv; 
    namej         = NAMEj; 
    Commij[NAMEj] = MPI_COMM_NULL;
    // 
    string PATH = "./"; ///home/jmake/z2016/REPOSITORY/COMMDOMs/SLIDING01/MESH01/"; 
    Files["types"]        = PATH +"Mesh01_TYPES.dat";       // 8040  
    Files["connectivity"] = PATH +"Mesh01_ELEMENTS.dat";
    Files["points"]       = PATH +"Mesh01_COORDINATES.dat"; // 4184 
    init_col = 1; 
  }

  if( name_argv == NAMEj )
  {
    namei         = name_argv; 
    namej         = NAMEi; 
    Commij[NAMEi] = MPI_COMM_NULL;
    // 
    string PATH = "./"; ///home/jmake/z2016/REPOSITORY/COMMDOMs/SLIDING01/MESH01/"; 
    Files["types"]        = PATH +"A_TYPES.dat";       // 2022 
    Files["connectivity"] = PATH +"A_ELEMENTS.dat"; 
    Files["points"]       = PATH +"A_COORDINATES.dat"; // 1071  
    init_col = 1;                                            // NOTE: 1==ALYA TYPE 
  }
 
  //---------------------------------------------------------------------||---//
  MPI_Barrier(PLEPP_COMM_WORLD);

  for(FilesIt it=Files.begin(); it!=Files.end(); it++)
  {
    string name = it->first; 
    string path = it->second;   

    read_log_file DATA; 
    DATA.set_name( path  );
    DATA.run();
    
    vector< vector<double> > vdata( DATA.get_vdata() );

    if(name=="points")
    {
      for(int i=0,k=0; i<vdata.size(); i++)
      {
        mui::point3d loc( vdata[i][0], vdata[i][1], vdata[i][2] );
        interface.push( "ux", loc, 0.0);
      }
      string sent = ( interface.commit(0.0) )?("ON"):("OFF");
      cout<<"sending:"<< sent <<"\n";

      mui::sampler_gauss3d<double> gauss( 1.0, 1.0/2 );
      mui::chrono_sampler_exact1d  exact;
      for(int i=0,k=0; i<vdata.size(); i++)
      {
        mui::point3d loc( vdata[i][0], vdata[i][1], vdata[i][2] );
        auto ux = interface.fetch( "ux", loc, 0.0, gauss, exact );
      }

    }
  }

  //---------------------------------------------------------------------||---//


  cout<<" '"<< name_argv <<"' ";  
  cout<<"\n";
  //---------------------------------------------------------------------||---//


  if(local_comm != MPI_COMM_NULL)
  {
/*
    mui::point3d loc( 0.0, 0.0, 0.0);
    interface.push( "ux", loc, 0.0);
    interface.commit( 0.0 );
*/
/*
    interface.push( "u", x[i], u[i] );
    interface.commit( k );
*/ 
  } 

/*
  //-----------------------------------------------------------| LOCATION |---//
  MPI_Barrier(PLEPP_COMM_WORLD);

  MPI_Comm  commij = MPI_COMM_NULL;
    if(commij != MPI_COMM_NULL)
    {

      int n_vertices_i = points.size()/3; 
      double* vertex_coords_i_ptr = points.data(); 

      int n_elements_i = types.size();   
      int* vertex_num_i_ptr  = connectivity.data(); 
      int* vertex_type_i_ptr = types.data(); 

      int n_vertices_j = points.size()/3;
      double* vertex_coords_j_ptr = points.data(); 

    }
*/
  //--------------------------------------------------------------------||---//
  MPI_Barrier(PLEPP_COMM_WORLD);
  MPI_Finalize();

  return 0;
} // main


//=======================================================================||===//
//=======================================================================||===//
