/*
  2019FEB18. BSC, BARCELONA, SPAIN 
  Migue Zavala 
*/
#include <iostream> // cin, cout, endl, cerr
#include <vector>   // vector
#include <map>
#include <mpi.h>
//#include <omp.h>
#include <cmath>
#include <assert.h>     /* assert */
#include <algorithm>    // std::fill
#include <sstream>      // std::stringstream
#include <iomanip>
#include <fstream>      // std::ofstream
#include <limits>       // 

#include "commdom.hpp"
#include "read_file.hpp"
using namespace std;

void get_tetra_coords_j(double*  coords, int* vertices, int n_dist_j, int* elemts_i, double*      coords_j, double* tetracoords_j);
void        propi2propj(double* props_i, int* vertices, int n_dist_j, int* elemts_i, double* tetracoords_j, double*       props_j);

int n_node = -1; // tetras:4, tria:3 
int    DIM =  3;
int    DIM_PROPS; //MATMATMAT 

CommDom   CD = CommDom();
string IAM   = "PLEPP";
string NAMEi = "TETRAS"; 
string NAMEj = "HEXAS"; 

int main(int argc, char** argv)
{
  //--------------------------------------------------------------------||---//
  int idx = -1;
  int init_col=-1;
  int elemt_type=-1;  

  string  name_argv = "";
  if(argc==2) name_argv = argv[1];

  string         namei  = ":(";
  string         namej  = ":(";
  string         dicti  = ":(";

  map<string, vector<int>    >            ArrayI;
  map<string, vector<double> >            ArrayD;

  map<string,string>                      Files;
  typedef map<string,string>::iterator    FilesIt;

  map<string,MPI_Comm>                    Commij;
  typedef map<string,MPI_Comm>::iterator  It;

  string PATH = "./"; 
  if( name_argv == NAMEi )
  {
    namei         = name_argv; 
    namej         = NAMEj; 
    Commij[NAMEj] = MPI_COMM_NULL;
    // 
    PATH += "/TETRAS01_01/Tetras01_01";
    Files["types"]        = PATH +"_TYPES.dat";       
    Files["connectivity"] = PATH +"_ELEMENTS.dat";
    Files["points"]       = PATH +"_COORDINATES.dat"; 
    init_col   = 1;   
    elemt_type = 4;    
    n_node     = 4;  
  }

  if( name_argv == NAMEj )
  {
    namei         = name_argv; 
    namej         = NAMEi; 
    Commij[NAMEi] = MPI_COMM_NULL;
    // 
    PATH += "/HEXAS01_01/Hexas01_01";
    Files["types"]        = PATH +"_TYPES.dat";       // 10370 
    Files["connectivity"] = PATH +"_ELEMENTS.dat";
    Files["points"]       = PATH +"_COORDINATES.dat"; // 5441 
  //Files["metis"]        = PATH +"Mesh02_ELEMENTS.metis.epart.3";
    init_col   = 1;                                                             // <-- NOTE: 1==ALYA TYPE 
    elemt_type = 6;
    n_node     = 8;
  }

  vector<read_log_file> DATA( Files.size() );

  //--------------------------------------------------------------| PLEPP |---//
  MPI_Init(NULL, NULL);
  MPI_Comm PLEPP_COMM_WORLD = MPI_COMM_WORLD;

  int  app_id = -1;
  int  n_apps = -1;
  MPI_Comm  local_comm;

  // CommDom
  CD.init();
  CD.set_app_type(IAM);
  CD.set_world_comm(PLEPP_COMM_WORLD);
  CD.set_app_name(namei); 
  CD.name_to_id(&app_id, &n_apps, &local_comm);
  CD.__create_interaction__();
  CD.create_commij(local_comm);
  
  int local_rank = -1;
  int local_size = -1;
  MPI_Comm_rank(local_comm, &local_rank);
  MPI_Comm_size(local_comm, &local_size);

  //---------------------------------------------------------------------||---//
  //MPI_Barrier(PLEPP_COMM_WORLD);

  int n_cells = -1; 
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
      ArrayD[name] = vector<double>(0);
      for(int i=0,k=0; i<vdata.size(); i++)
      {
      //for(int j=init_col; j<vdata[i].size(); j++) ArrayD[name].push_back( vdata[i][j] );
        for(int j=init_col,l=0; j<DIM+init_col; j++) ArrayD[name].push_back( vdata[i][j] );  
      }
    }
    else if(name=="metis")
    {
      ArrayI[name] = vector<int>(0);
      for(int i=0,k=0; i<vdata.size(); i++)
      {
        for(int j=0; j<vdata[i].size(); j++) ArrayI[name].push_back( (int)vdata[i][j] );
      }
    }
    else
    {
      ArrayI[name] = vector<int>(0);  
      for(int i=0,k=0; i<vdata.size(); i++) 
      {
        for(int j=init_col; j<vdata[i].size(); j++) ArrayI[name].push_back( (int)vdata[i][j] ); 
      }
    }

    DATA.end();
  }

  vector<int>        types( ArrayI.find("types")->second );
  vector<int> connectivity;
  vector<double>    points( ArrayD.find("points")->second );

  //---------------------------------------------------------------------||---//
  if(ArrayI.find("metis") != ArrayI.end())
  {
  }
  else 
  {
    read_log_file DATA;
    DATA.set_name( Files.find("connectivity")->second );
    DATA.run();
    DATA.end();
    vector< vector<double> > aux( DATA.get_vdata() );
 
    for(int i=0; i<aux.size(); i++)
    {
      for(int j=init_col; j<aux[i].size(); j++)
      {
        connectivity.push_back( (int)aux[i][j]  );
      }
    }
    cout<<"["<< namei <<"] n_pts:"<<  points.size()/DIM <<"\n";
    cout<<"["<< namei <<"] n_cells:"<<  aux.size()  <<"\n";  

    aux.clear();
  }

  //-----------------------------------------------------------| LOCATION |---//
  MPI_Barrier(PLEPP_COMM_WORLD);

  MPI_Comm  commij = MPI_COMM_NULL;
  if( (CD.__get_app_name__() == namei)&&(CD.__get_friends__(namej) == 1) ) commij = CD.get_mpi_commij(namej);
  if( (CD.__get_app_name__() == namej)&&(CD.__get_friends__(namei) == 1) ) commij = CD.get_mpi_commij(namei); 

  int            n_vertices_i = 0;
  double* vertex_coords_i_ptr = NULL;

  int            n_vertices_j = 0;
  double* vertex_coords_j_ptr = NULL; 
  double*  vertex_props_j_ptr = NULL; //MATMATMAT 

  int            n_elements_i = 0;
  int*       vertex_num_i_ptr = NULL;
  int*      vertex_type_i_ptr = NULL;

  if(commij != MPI_COMM_NULL)
  {

    DIM_PROPS = 4; //MATMATMAT 
    int num_nodes = points.size()/DIM;
    vector<double>    propspropsi(DIM_PROPS*num_nodes); //MATMATMAT
 
    // Iterations  
    double time = 0.0;   
    for(int itime=0; itime<1; itime++) 
    { // for2  

             n_vertices_i = points.size()/DIM; 
      vertex_coords_i_ptr = points.data(); 

             n_elements_i = types.size();   
         vertex_num_i_ptr = connectivity.data(); 
        vertex_type_i_ptr = types.data(); 

      if(local_rank==0)
      {
               n_vertices_j = points.size()/DIM;
        vertex_coords_j_ptr = points.data();
         vertex_props_j_ptr = propspropsi.data(); //MATMATMAT
      }

      CD.locator_create2(local_comm, commij, 1e-3, elemt_type);
      CD.locator_set_cs_mesh(n_vertices_i,
                             n_elements_i,
                             vertex_coords_i_ptr,
                             vertex_num_i_ptr,
                             vertex_type_i_ptr, 
                             n_vertices_j,
                             vertex_coords_j_ptr, DIM, vertex_props_j_ptr, DIM_PROPS); //MATMATMAT
  
      CD.save_dist_coords(0, local_comm);
      cout<<" \n"; 
/* 
      double dt = 0.0;  
      // Time step 
      {
        double send[1]={1.0e-2}, recv[1]={1.0/HUGE_VAL};
        int n_send=1, n_recv=1; 
        CD.__mpi_sendrecv_real__( send, n_send, recv, n_recv, local_comm, commij);
        CD.__mpi_bcast_real__(                  recv, n_recv, local_comm, commij);  

        time += 1.0/recv[0]; 
        dt    = 1.0/recv[0];  
      }

      // Local properties (to be sent)  
      vector<double>  props_i(n_vertices_j*DIM, HUGE_VAL);
      vector<double>  disp(2); 
      disp[0] =  0.0; 
      disp[1] =  0.0; //-1.0e-5; 
     

      for(int i=0,j=0; i<DIM; i++)
      {
        for(int k=0; k<n_vertices_j; k++) props_i[j++] = disp[i]; //-(i+1);  
      } 
      
      // 
      int n_recv = CD.get_n_interior();
      vector<double> var_ji(n_recv*DIM, HUGE_VAL);
      vector<int> interior_list_j(n_recv, -1); 
      CD.__locator_get_interior_list__( interior_list_j.data() ); 

      // 
      int n_send = CD.get_n_dist_points();

      vector<int>    dist_locations_i(n_send,                 -1);
      vector<double>    dist_coords_j(n_send*DIM,       HUGE_VAL);
      vector<double>     dist_props_j(n_send*DIM_PROPS, HUGE_VAL); //MATMATMAT

      CD.__locator_get_dist_locations__( dist_locations_i.data() );
      CD.__locator_get_dist_coords__(       dist_coords_j.data() );
      CD.__locator_get_dist_props__(         dist_props_j.data() ); //MATIAS

      vector<double>   shapef_j( n_node*n_send, HUGE_VAL);
      get_tetra_coords_j(vertex_coords_i_ptr, vertex_num_i_ptr,
                         n_send, dist_locations_i.data(), dist_coords_j.data(), shapef_j.data() );

      vector<double> var_ij(n_send*DIM, HUGE_VAL);
      for(int i=0; i<DIM; i++)
      {
        propi2propj(props_i.data()+n_vertices_j*i, vertex_num_i_ptr, n_send, dist_locations_i.data(), shapef_j.data(), var_ij.data()+n_send*i); 
      } 

      vector<double> aux(var_ij); 
      for(int i=0,j=0; i<n_send; i++)
      {
        for(int k=0; k<DIM; k++, j++) var_ij[j] = aux[k*n_send + i]; 
      }
      aux.clear(); 

      // Exchange  
      CD.__locator_exchange_double_scalar__(var_ij.data(), var_ji.data(), DIM); 
*/ 
      // 
      CD.locator_destroy();

    } // for2 
 
  } //  commij!=MPI_COMM_NULL

  //--------------------------------------------------------------------||---//
//  MPI_Barrier(PLEPP_COMM_WORLD);
//  MPI_Abort(commij, -1);        
  MPI_Finalize();

  return 0;
} // main
//=======================================================================||===//
//=======================================================================||===//

void
propi2propj(double* props_i, int* vertices, int n_dist_j, int* elemts_i, double* tetracoords_j, double* props_j )
{
  int    ii, ielem;
  int      vertices_i[n_node];
  double vol_coords_j[n_node];
  double       prop_j[n_node];

  if(n_dist_j > 0)
  {
    for(int ii=0; ii<n_dist_j; ii++ )
    {
      ielem = elemts_i[ii]-1; //  <- Fortran style ?? 
      for(int jj=0; jj<n_node; jj++)   vertices_i[jj]  =      vertices[n_node*ielem + jj]-1; // <- Fortran style -> C style !!  
      for(int jj=0; jj<n_node; jj++) vol_coords_j[jj]  = tetracoords_j[n_node*ii    + jj];
      for(int jj=0; jj<n_node; jj++)       prop_j[jj]  =       props_i[ vertices_i[jj]  ];

      props_j[ii] = 0.0;
      for(int jj=0; jj<n_node; jj++)      props_j[ii] += prop_j[jj] * vol_coords_j[jj];
    }
  }
}


void
get_tetra_coords_j(double* coords, int* vertices, int n_dist_j, int* elemts_i, double* coords_j, double* tetracoords_j)
{
  int    ii, ielem;
  int      vertices_i[n_node];
  double vol_coords_j[n_node];
  double      point_j[DIM];

  if(n_dist_j > 0)
  {
    for(int ii=0; ii<n_dist_j; ii++ )
    {
      for(int jj=0; jj<DIM;     jj++)   point_j[jj] = coords_j[DIM*ii + jj];

      ielem = elemts_i[ii]-1; //  <- Fortran style ?? 
      for(int jj=0; jj<n_node; jj++) vertices_i[jj] = vertices[n_node*ielem + jj]; // <- Fortran style, yes!! (ple rest one) 

      CD.__simple_interpolation__( coords, vertices_i, point_j, vol_coords_j );

      for(int jj=0; jj<n_node; jj++) tetracoords_j[n_node*ii + jj] = vol_coords_j[jj];
    }
  }
}

//=======================================================================||===//
//=======================================================================||===//
