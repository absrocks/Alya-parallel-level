#include "commdom.hpp"
#include <algorithm>  // std::search 
#include <iomanip>    // setw 
#include <cmath>      // fabs
#include <fstream>
#include <sstream>     
#ifdef SATURNE==1
#include "cs_alya.h"
#include "cs_alya.c"
#endif 
//-----------------------------------------------------------------------||---//
ple_mesh_extents_t           *mesh_extents_f   = NULL;
ple_mesh_elements_contain_t  *locate_inside_f  = NULL;
ple_mesh_elements_closest_t  *locate_closest_f = NULL;
//-----------------------------------------------------------------------||---//


#ifndef __PRINT__
#define __PRINT__ 0 
#endif 
//=======================================================================||===//
//=======================================================================||===//


//=======================================================================||===//
//================================================================| COMMs |===//

//-----------------------------------------------------------------------||---//
CommDom::CommDom()
{
  if(__PRINT__) cout<<"+CommDom Create \n"; 
}


CommDom::~CommDom()
{
  if(__PRINT__) cout<<"+CommDom \'"<< IAM <<"\' End \n"; 
}


void 
CommDom::init()
{  
       IAM = "";
  app_name = "";

  quitze_comm = MPI_COMM_NULL;
 //local_comm = MPI_COMM_NULL;
   world_comm = MPI_COMM_NULL;
       commij = MPI_COMM_NULL;

   local_size = -1; 
   local_rank = -1;
   world_size = -1; 
   world_rank = -1;

  app_namei = ""; 
  app_namej = ""; 

  cplng    = NULL; 

  n_meshes =  0; 
  n_locs   =  0; 
  n_apps   = -1; 
  n_types  = -1; 
  app_id   = -1; 

  if(__PRINT__) cout<<"|_CommDom.init() \n"; 
}


void 
CommDom::end()
{  
  if(__PRINT__) cout<<"|_CommDom.end() \n"; 
}
//-----------------------------------------------------------------------||---//

//-----------------------------------------------------------------------||---//
void
CommDom::set_app_type(string last_name, int length) 
{
  //last_name.erase(std::remove_if(last_name.begin(), last_name.end(), ::isspace), last_name.end());
  if(length==-1) 
  {
    app_type = last_name; 
  }
  else
  {
  /*
    //  cout<< last_name.size() <<" \n\n";
    cout<<"(1)->\'"<< last_name <<"\' \n\n";
    app_type = string_f2c(last_name.c_str(), last_name.size());
    cout<<"(2)->\'"<< app_type <<"\' \n\n";
  */
  }
  
  if(__PRINT__) cout<<"|_CommDom.set_app_type \n"; 
}


void
CommDom::set_app_name(string first_name)
{
  app_name = first_name; 
  IAM      = first_name; 

  if(__PRINT__) cout<<"|_CommDom.set_app_name \n"; 
}


void 
CommDom::__genesis__() 
{
//
//  |_ ./src/syrthes-kernel/src/syrthes.c
//    |_ syrthes_initmpi
//
//  if exit "--name"
//  cfd.do_coupling = 1; // mainsyrthes.c 
//

  MPI_Comm  Gucumatz = MPI_COMM_WORLD; 

  int w_rank, w_size;
  MPI_Comm_rank(Gucumatz, &w_rank);
  MPI_Comm_size(Gucumatz, &w_size);

  int id_type = ple_coupling_mpi_name_to_id(Gucumatz, app_type.c_str());
/*
int id_name = ple_coupling_mpi_name_to_id(Gucumatz, app_name.c_str());
*/
  quitze_comm = MPI_COMM_NULL;
  if(id_type > -1) 
  {

    MPI_Comm_split(Gucumatz, id_type, w_rank, &quitze_comm);

    int sync_flag = _cs_coupling_sync_flag;
    ple_coupling_mpi_set_t *cplng_aux = 
                           ple_coupling_mpi_set_create(sync_flag, app_type.c_str(), app_name.c_str(), Gucumatz, quitze_comm);
    
    
    n_types = ple_coupling_mpi_set_n_apps(cplng_aux);

    cout<<"|_["; 
    //cout<<"n_types:"<< n_types; 
    cout<< w_rank+1 <<"."<< app_type; 
    //cout<<"."<< app_name;  
    cout<< "] ";
    //cout<< id_type+1 <<"/"<< id_name+1 <<" ";  
    __print_comm__(quitze_comm);
    cout<<"\n";
  }
  else
  {
    n_types = 1; 
  }

  //cout<< "<--------";
  //cout<<"\n";
}


void 
CommDom::set_world_comm(MPI_Comm  comm)
{
  world_comm = comm; 
  if(quitze_comm != MPI_COMM_NULL) world_comm = quitze_comm; // communicator given by app_type  

  MPI_Comm_size(world_comm, &world_size);
  MPI_Comm_rank(world_comm, &world_rank); 

  if(__PRINT__) cout<<"|_CommDom.set_world_comm \n"; 
}


void 
CommDom::name_to_id(int* appi, int* napps, MPI_Comm* lcomm)
{
  MPI_Comm local_comm = world_comm;
 
  app_id = ple_coupling_mpi_name_to_id(world_comm, app_name.c_str());
  if(app_id > -1) MPI_Comm_split(world_comm, app_id, world_rank, &local_comm); 

  cout<<"|_["<< world_rank+1 <<"."<< app_type <<"."<< app_name << "] ";
  __print_comm__(local_comm); cout<<"\n";

  __set_create__(app_name, local_comm, &n_apps); 
  if(app_id < 0){ app_id = 0; n_apps= 1;}

   appi[0] = app_id; 
  napps[0] = n_apps; 
  lcomm[0] = local_comm; 

  __set_get_info__(); 

  if(__PRINT__) cout<<"|_CommDom.name_to_id \n"; 
}


void
CommDom::__set_create__(string namei, MPI_Comm commi, int* napps) 
{
  int sync_flag = _cs_coupling_sync_flag; // & PLE_COUPLING_UNSTEADY 

  cplng = ple_coupling_mpi_set_create(sync_flag, app_type.c_str(), namei.c_str(), world_comm, commi);
  napps[0] = ple_coupling_mpi_set_n_apps(cplng);
}


void
CommDom::__set_get_info__()
{
  for(int i = 0; i < n_apps; i++) // n_apps <-- name_to_id() 
  { 
    ple_coupling_mpi_set_info_t  infoi = ple_coupling_mpi_set_get_info(cplng, i);
    Info.push_back( infoi ); 

    WhoIam[infoi.app_name] = i; 
    //CplngIJ.push_back( vector<int>() ); 
  }
}


void
CommDom::__create_interaction__()
{
  cout<<"|_["<< world_rank+1 <<".\'"<< app_type <<"\'.\'"<< app_name << "\'] "; 

  string namei = app_name; 

//for(Whoi = WhoIam.begin(); Whoi!=WhoIam.end(); ++Whoi) 
//{
//  namei = Whoi->first;
  cout<<" | "; 

  for(int i = 0; i < n_apps; i++) CplngIJ.push_back( vector<int>() ); 
  
  for(Whoj = WhoIam.begin(); Whoj!=WhoIam.end(); ++Whoj) 
  {
    string  namej = Whoj->first; 
    if( (namei==namej) )
    { 
      //namej = "----";
      //cout<<"["<< namej <<"]";
    }
    else
    {
      bool you_are_friendly = true; //find( friendsi.begin(), friendsi.end(), namej) != friendsi.end(); 
      if(you_are_friendly==true)
      {
        cout<<"["<< namej <<"]";
        CplngIJ[ WhoIam[namei] ].push_back( WhoIam[namej] ); 
      }
      else
      {
        //namej = "----";
        //cout<<"["<< namej <<"]";
      }
    }
  }
  cout<<" | "; 
//}


  int appi = WhoIam[app_name]; 
  //cout<<"CplngIJ="<< CplngIJ[appi].size(); 
  cout<<"\n";

  if(__PRINT__) cout<<"|_CommDom.__create_interaction__ \n"; 
}


//=======================================================================||===//
//================================================================| COMMs |===//
//-----------------------------------------------------------------------||---//
void
CommDom::create_commij(MPI_Comm lcomm) 
{
  int appi = -1, appj = -1; 
  
  appi = WhoIam[app_name];  // <-- Sooo... important!!
  //if(CplngIJ[appi].size()<1)  __error__(true, "interaction doesnt exist");

  string namei = app_name; //cout<<"|_["<< namei;
  for(int j = 0; j < CplngIJ[appi].size(); j++)
  {
    appj = CplngIJ[appi][j]; 

    string namej = Info[appj].app_name; 
    if(world_size<=4) cout<<"|_["<< namei <<"."<< Info[appi].root_rank+1 <<"]->["<< namej <<"."<< Info[appj].root_rank+1 <<"]: "; 

    vector<int> rangei(2); 
    rangei[0] = Info[appi].root_rank; 
    rangei[1] = Info[appi].root_rank + Info[appi].n_ranks; 
    
    vector<int> rangej(2); 
    rangej[0] = Info[appj].root_rank; 
    rangej[1] = Info[appj].root_rank + Info[appj].n_ranks; 

    if(world_size<=4) 
    {
      cout<<"|";
      cout<<"["<< appi+1 <<"]("<< rangei[0]+1 <<","<< rangei[1] <<")->"; 
      cout<<"["<< appj+1 <<"]("<< rangej[0]+1 <<","<< rangej[1] <<")| "; 
      //cout<<"|\n";
    }

    int commij_size = -1; 
    int commij_rank = -1;

    commij = MPI_COMM_NULL;
    ple_coupling_mpi_intracomm_create(world_comm, lcomm, rangej[0], &commij, rangei.data(), rangej.data());  
    CommIJ[namej] = commij; 

    MPI_Comm_size(commij, &commij_size);
    MPI_Comm_rank(commij, &commij_rank); 

    if(world_size<=4)
    {
      cout<< commij_rank+1 <<"/"<< world_rank+1 <<","<< commij_size <<"/"<< world_size; 
      cout<<"\n";
    }

    int rooti = -1; 
    int rootj = -1; 
    __mpi_get_roots__(lcomm, commij, &rooti, &rootj); 
    RootIJ[namej].push_back(rooti);
    RootIJ[namej].push_back(rootj);

  }
  MPI_Barrier(lcomm); 

  if(__PRINT__)
  {
    cout<<"|_CommDom.create_commij: |"; 
    //__print_comm__(lcomm); 
    //__print_comm__(world_comm); 
    print_commij_names(); 
    cout<<"|\n";
  }
}
//-----------------------------------------------------------------------||---//


//-----------------------------------------------------------------------||---//
void
CommDom::get_commij_size(int* size) 
{
  size[0] = CommIJ.size(); 
}


void
CommDom::print_commij_names() 
{
  cout<<"["<< app_type <<"."<< app_name <<"."<< world_rank+1 <<"] -> "; 
  for(CommIT = CommIJ.begin(); CommIT != CommIJ.end(); ++CommIT) 
  cout<<"["<< CommIT->first <<"."<< app_type <<"."<< MPI_Comm_c2f(CommIT->second) <<"]"; 
}


void
CommDom::get_commij(string  namej, MPI_Comm*  commij) 
{
  commij[0] = MPI_COMM_NULL; 

  CommIT = CommIJ.find(namej); 
  if(CommIT == CommIJ.end()) 
  {
    cout<<"\n\nERROR: ["; 
    cout<< app_type <<".";  
    cout<< namej;
    cout<<"] not found. INTERNSHIPS->"; 
    print_commij_names(); 
    cout<<"!!\n\n";
    exit(1);
  }
  else
  {
    if(local_size<=4) cout<<"|_["<< IAM <<"."<< world_rank+1 <<"]<->["<< namej <<"] "; 
    commij[0] = CommIJ[namej]; 
  }
}


MPI_Comm
CommDom::set_mpi_comms() // used by python 
{
  MPI_Comm lcomm = MPI_COMM_NULL;
  int appi  =-1; 
  int napps =-1; 
  
  name_to_id(&appi, &napps, &lcomm); 
  __create_interaction__(); 
  create_commij(lcomm);

  return lcomm; 
}


int
CommDom::get_mpi_commij_size() // used by python 
{
  int size;
  get_commij_size(&size);
  return size; 
}


MPI_Comm
CommDom::get_mpi_commij(string namej) // used by python 
{
  MPI_Comm  commij = MPI_COMM_NULL; 
  get_commij(namej, &commij);
  return commij; 
}


//=======================================================================||===//
//===================================================================| MPI|===//
//-----------------------------------------------------------------------||---//
int
CommDom::__mpi_get_roots__(MPI_Comm commi, MPI_Comm commij, 
                           int* rooti, int* rootj)
{
  int ranki  = -1;
  int rankij = -1;
  int sizei  = -1;
  int sizeij = -1;

  MPI_Comm_size(commi,  &sizei );
  MPI_Comm_size(commij, &sizeij);
  MPI_Comm_rank(commi,  &ranki );
  MPI_Comm_rank(commij, &rankij);

  int rootij = -666; 
  if(ranki==rankij)
  {
    rooti[0] = 0;
    rootj[0] = sizei;
    rootij = 1;
  }
  else
  {
    rooti[0] = sizeij-sizei;
    rootj[0] = 0;
  }

  //if(rankij==0) rootij = 1; 

  return rootij; 
}



//=======================================================================||===//
//==============================================================| LOCATOR |===//

//-----------------------------------------------------------------------||---//
void
CommDom::locator_create2(
                         MPI_Comm   commi,
                         MPI_Comm   commij, 
                         double     tol, 
                         string     namej)
{
  int ranki  = -1;
  int rankij = -1;
  int sizei  = -1;
  int sizeij = -1;

  MPI_Comm_size(commi,  &sizei );
  MPI_Comm_size(commij, &sizeij);
  MPI_Comm_rank(commi,  &ranki );
  MPI_Comm_rank(commij, &rankij);

  int rooti=-1, rootj=-1; 
  if(ranki==rankij)
  {
    rooti = 0;
    rootj = sizei;
  }
  else
  {
    rooti = sizeij-sizei;
    rootj = 0;
  }


  int root_rank = rootj;        //= Info[appj].n_ranks;
  int n_rank    = sizeij-sizei; //= Info[appj].root_rank;

  Locators[namej] = ple_locator_create(tol, commij, n_rank, root_rank);  
}


//-----------------------------------------------------------------------||---//
void
CommDom::locator_set_mesh(
                          int         n_vertices_i, 
                          int         n_elements_i, 
                          REAL*       vertex_coords_i, //double*     vertex_coords_i, 
                          INTEGER*    vertex_num_i,    //int*        vertex_num_i, 
                          int         n_vertices_j, 
                          REAL*       vertex_coords_j,  //double*     vertex_coords_j
                          string      namej
                          ) 
{
  int idx = 0; 
  vector<syr_cfd_mesh_t>      LocalMesh;
  LocalMesh.push_back( syr_cfd_mesh_t() ); 

  // Set Meshi
  LocalMesh[idx].dim           = __DIM__;
  LocalMesh[idx].element_type  = SYR_CFD_TETRA;
  LocalMesh[idx].n_vertices    = n_vertices_i;
  LocalMesh[idx].n_elements    = n_elements_i;
  LocalMesh[idx].vertex_coords = NULL;
  LocalMesh[idx].vertex_num    = NULL;

  if(n_vertices_i >= 0) LocalMesh[idx].vertex_coords = vertex_coords_i;
  if(n_elements_i >= 0) LocalMesh[idx].vertex_num    = vertex_num_i;


  // Set Meshj
  n_coords_dist = 0;
  coords_dist   = NULL;
  distance_dist = NULL;

  if(n_vertices_j > 0)
  {
    n_coords_dist = n_vertices_j;
      coords_dist = vertex_coords_j;
    distance_dist = new float[n_coords_dist];
  }


  // Set locator  
  int               dim = __DIM__;
  const int *point_list = NULL;

  ple_locator_t *locatorj = Locators[namej]; 

  __locator_set_meshij__(locatorj,
                        &LocalMesh[idx],
                         dim,
                         n_coords_dist,
                         point_list,
                         coords_dist,
                         distance_dist);

  LocalMesh.clear(); 
}


void
CommDom::__locator_set_meshij__(ple_locator_t      *locator, 
                                syr_cfd_mesh_t     *Mesh, 
                                int                 dim, 
                                int                 n_coords_dist, 
                                const int          *point_list, 
                                double             *coords_dist, 
                                float              *distance_dist 
                              ) 
{
      mesh_extents_f   = syr_cfd_point_location_extents;
      locate_inside_f  = syr_cfd_point_location_contain;
      locate_closest_f = NULL;               // Why dosent it work? R. 3d not implemented 

      ple_locator_set_mesh(locator,          //*ple_locator_t
                           Mesh,             // void*
                           dim,              // coord_dim!! 
                           n_coords_dist,    // ple_lnum_t 
                           point_list,       // const ple_lnum_t[] 
                           coords_dist,      // const ple_coord_t[]
                           distance_dist,    // float[]
                           mesh_extents_f,   // ple_mesh_extents_t.          cs_coupling_mesh_extents
                           locate_inside_f,  // ple_mesh_elements_contain_t. cs_coupling_point_in_mesh. 
                           locate_closest_f);// ple_mesh_elements_closest_t. 
                                             // locate_on_closest = cs_coupling_point_closest_mesh. optional 

    //ple_locator_dump(locator); 

    //syr_cfd_coupling._interpolation_init ?? 
    ple_locator_exchange_point_var(locator,
                                   NULL,
                                   distance_dist,
                                   NULL,
                                   sizeof(float),
                                   1,
                                   1);
    if(world_size<=8)
    {
        cout<<"|_["<< app_type <<"."<< app_name <<"] ";
        cout<<"n_elemets: "; 
        cout<< Mesh[0].n_elements <<", "; 
        cout<<"n_exterior: ";
        cout<< ple_locator_get_n_exterior(locator)<<", ";
        cout<<"n_interior: ";
        cout<< ple_locator_get_n_interior(locator)<<", ";
        cout<<"n_dist_points: "; 
        cout<< ple_locator_get_n_dist_points(locator); 
        cout<<"\n";
    }

}
//-----------------------------------------------------------------------||---//

//-----------------------------------------------------------------------||---//
int
CommDom::get_n_dist_points(string namej) 
{
  ple_locator_t *locatorj = Locators[namej];
  return ple_locator_get_n_dist_points(locatorj); 
}


int
CommDom::get_n_interior(string namej) 
{
  ple_locator_t *locatorj = Locators[namej];
  return ple_locator_get_n_interior(locatorj); 
}

//-----------------------------------------------------------------------||---//
void
CommDom::__print_comm__(MPI_Comm  comm)
{
  int size=-1, rank=-1;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank); 

  MPI_Fint  fcomm = MPI_Comm_c2f(comm); 

  cout<<"["<< rank+1 <<"/"<< size <<"]"; 
}
//-----------------------------------------------------------------------||---//


//-----------------------------------------------------------------------||---//
void
CommDom::__error__(bool ok, string message)
{
  if(ok==true) 
  {
    cout<<"\n\nERROR: "<< message  <<"!! \n\n"; 
    exit(1);
  }
}


void
CommDom::__strncmp__(const char* x, const char* y, int* size, int* ok)
{
  ok[0] = -1;
  if(strncmp(x, y, size[0]) == 0) 
  {
    ok[0] = 1;
  }
  else
  {
   //printf(" \'%s\', \'%s\' \n", x, y); 
   /*
    cout<<"\n\nWARNING: ";
    cout<<"\'"<< x <<"\' != \'"<< y <<"\'";   
    cout<<"!! \n\n";
    //exit(1);
   */
  } 
}


const char* 
CommDom::string_f2c(const char* fchar, int fchar_length)
{
  int i = fchar_length-1; 
  while( ((fchar[i]==' ')||(fchar[i]==0)) && (i>=0) ) i--; 
  i+=1; 
    
  std::string  cstr; 
  for(int j=0; j<i; j++) cstr.push_back(fchar[j]); 

  return cstr.c_str(); 
}


vector<string> 
CommDom::__split_string__(const string& s, const string& delim, const bool keep_empty) 
{
    vector<string> result;
    if(delim.empty()) 
    {
        result.push_back(s);
        return result;
    }
    string::const_iterator substart = s.begin(), subend;

    while(true) 
    {
        subend = search(substart, s.end(), delim.begin(), delim.end());
        string temp(substart, subend);
        if(keep_empty || !temp.empty()) 
        {
            result.push_back(temp);
        }
        
        if(subend == s.end()) 
        {
            break;
        }
        
        substart = subend + delim.size();
    }
    return result;
}
//-----------------------------------------------------------------------||---//


//-----------------------------------------------------------------------||---//
string
CommDom::__get_app_type__()
{
  return app_type; 
}


string
CommDom::__get_app_name__()
{
  return IAM; 
}


int 
CommDom::__get_n_apps__()
{
  return n_apps; 
}


int 
CommDom::__get_n_types__()
{
  return n_types; 
}


MPI_Comm 
CommDom::__get_world_comm__()
{
  return world_comm; 
}


int 
CommDom::__get_friends__(string  namej) 
{
  int ok[] = {1}; 

  CommIT = CommIJ.find(namej); 
  if(CommIT == CommIJ.end()) 
  {
    cout<<"WARNINIG: ["<< namej <<"] not found!! \n"; 
    ok[0] = -1; 
  }

  return ok[0]; 
}
//-----------------------------------------------------------------------||---//


//-----------------------------------------------------------------------||---//
void
CommDom::__to_vtk__(string filename, int npts, const double* pts) 
{
  std::ofstream  fout(filename.c_str()); 
  fout<<"# vtk DataFile Version 3.0\n";
  fout<<"chingon!!"<<"\n"; 
  fout<<"ASCII \nDATASET POLYDATA \nPOINTS "<< npts <<" float \n";
  for(int i=0, k=0; i<npts; i++) 
  {
    for(int j=0; j<3; j++, k++) fout<< pts[k] <<" "; 
    //fout<< dist_locations[i] <<" "; 
    fout<<"\n";
  }
  fout.close();   
}
//-----------------------------------------------------------------------||---//


//-----------------------------------------------------------------------||---//
void
CommDom::__tetra_interpolation__(double* coords, int* vertices, double* point, double* shapef)
{
  double _epsilon_denom = 1.e-28; 
  double isop_0, isop_1, isop_2;
  double t00, t10, t20, t01, t02, t03, t11, t12, t13, t21, t22, t23;
  double v01[3], v02[3], v03[3]; //, shapef[4];
  double vol6;
  int i, j, k;


  ple_coord_t  tetra_coords[4][3];
  for(i = 0; i < 4; i++)
  {
    ple_lnum_t   coord_idx = vertices[i] - 1; //<-NOTE: fortran numeration 
    for(j = 0; j < 3; j++) tetra_coords[i][j] = coords[(coord_idx * 3) + j];
  }

  // vol[i] = (<x_j,y_j,z_j> - <x_0,y_0,z_0>)_i = <x_j-x_0, y_j-y_0, z_j-z_0>_i = R_j - R_0 = R_0j 
  for(i = 0; i < 3; i++)
  {
    v01[i] = tetra_coords[1][i] - tetra_coords[0][i];
    v02[i] = tetra_coords[2][i] - tetra_coords[0][i];
    v03[i] = tetra_coords[3][i] - tetra_coords[0][i];
  }

  // 6*Vol = J = A.BxC = R_01 . R_02 x R_03
  vol6 = fabs(  v01[0] * (v02[1]*v03[2] - v02[2]*v03[1])
              - v02[0] * (v01[1]*v03[2] - v01[2]*v03[1])
              + v03[0] * (v01[1]*v02[2] - v01[2]*v02[1]));
  if(vol6 < _epsilon_denom) return;


  // Tetrahedral coordinates 
  t00  = - tetra_coords[0][0] + point[0];
  t10  = - tetra_coords[0][1] + point[1];
  t20  = - tetra_coords[0][2] + point[2];

  t01  = - tetra_coords[0][0] + tetra_coords[1][0];
  t02  = - tetra_coords[0][0] + tetra_coords[2][0];
  t03  = - tetra_coords[0][0] + tetra_coords[3][0];

  t11  = - tetra_coords[0][1] + tetra_coords[1][1];
  t12  = - tetra_coords[0][1] + tetra_coords[2][1];
  t13  = - tetra_coords[0][1] + tetra_coords[3][1];

  t21  = - tetra_coords[0][2] + tetra_coords[1][2];
  t22  = - tetra_coords[0][2] + tetra_coords[2][2];
  t23  = - tetra_coords[0][2] + tetra_coords[3][2];

  isop_0 = (  t00 * (t12*t23 - t13*t22)
            - t10 * (t02*t23 - t22*t03)
            + t20 * (t02*t13 - t12*t03)) / vol6;
  isop_1 = (- t00 * (t11*t23 - t13*t21)
            + t10 * (t01*t23 - t21*t03)
            - t20 * (t01*t13 - t03*t11)) / vol6;
  isop_2 = (  t00 * (t11*t22 - t21*t12)
            - t10 * (t01*t22 - t21*t02)
            + t20 * (t01*t12 - t11*t02)) / vol6;

  shapef[0] = 1.0 - isop_0 - isop_1 - isop_2;
  shapef[1] =       isop_0;
  shapef[2] =                isop_1;
  shapef[3] =                         isop_2;

}
//-----------------------------------------------------------------------||---//

//----------------------------------------------------------------| argvs |---//
void
CommDom::__set_argvs__(string dummy)
{
  argvs.push_back(dummy); 
}


void
CommDom::__get_argvs__(string& dummy)
{
  dummy = argv_found; 
}


void
CommDom::__analyse_argvs__(string token)
{  
  argv_found = ""; 

  string s;
  
  int argc = argvs.size(); 
  int numarg = 0;
  while(++numarg < argc) 
  {
    s = argvs[numarg];
    if( strncmp(s.c_str(), token.c_str(), token.size()) == 0 ) 
    {
      //cout<<"["<< token <<"] = "; 
      
      if(++numarg < argc) 
      {
        s   = argvs[numarg];
        int len = s.size();
        char* log_name = (char *)malloc(sizeof(char)*(len+1));
        strncpy(log_name, s.c_str(), len);
        log_name[len] = '\0';

        //cout<<"\'"<< log_name <<"\'"; 

        argv_found = s; 
      }

      //cout<<"\n"; 
    }
  }

}


void
CommDom::save_dist_coords(MPI_Comm comm, string namej)
{

  int id = 0; 
  int rankij=-1;
  if( (comm != MPI_COMM_NULL) || (!comm<0) )
  {
    MPI_Comm_rank(comm, &rankij);
    id = rankij+1;
  }

  int appj = WhoIam[app_name];  // <-- Sooo... important!!

  std::string        filename(app_name);
  std::stringstream  srut;
  srut<< std::setfill('0') << std::setw(3) << id;
  
  filename += "_from_" + namej;  
  filename += "_" + srut.str();
  filename += ".vtk";

  ple_locator_t     *locatorj = Locators[namej];
  int           n_dist_points = ple_locator_get_n_dist_points(locatorj);
  const double   *dist_coord  = ple_locator_get_dist_coords(locatorj);

  __to_vtk__(filename, n_dist_points, dist_coord);

}

//=======================================================================||===//
//=======================================================================||===//
