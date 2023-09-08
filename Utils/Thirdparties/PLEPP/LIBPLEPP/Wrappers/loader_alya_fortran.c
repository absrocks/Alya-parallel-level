#include "loader_alya_wrapper.h"
#include "commdom_wrapper.h"

  // pointers
  //LOADER_ALYA  *loader, *loader_array;

//LOADER_ALYA  *loader, *loader_array;
#ifdef __cplusplus
extern "C" {
#endif


void 
loader_alya_wrapper_create()
{
  //loader = loader_alya_create();
  //loader_alya_set_ple_data(loader, "/home/jmake/z2014/Cplng/MeshCplng/MeshCoarse/disk", 1);
  //loader_alya_delete(loader);

  //loader_array = c_loader_alya_create_array(2);
  //c_loader_alya_set_ple_data(loader_array, 0, "/home/jmake/z2014/Cplng/MeshCplng/MeshCoarse/disk", 1);
  //c_loader_alya_set_ple_data(loader_array, 1, "/home/jmake/z2014/Cplng/MeshCplng/MeshCoarse/square", 1);
}



void 
loader_alya_wrapper_create_array(int num)
{
  //loader_array = c_loader_alya_create_array(num);
}


void 
loader_alya_wrapper_set_ple_data_array(int idx)
{
  //c_loader_alya_set_ple_data(loader_array, idx, "/home/jmake/z2014/Cplng/MeshCplng/MeshCoarse/disk", 1);
}


void
commdom_wrapper_create(char* iam, MPI_Fint* fcomm, char* namei, MPI_Fint* lcomm, MPI_Fint* commij, char* nameij[], int* Ruts)
{

  MPI_Comm  ccomm = MPI_Comm_f2c(fcomm[0]); 

  CommDom *CD = new CommDom();
  CD->init();
  CD->set_world_comm(ccomm);
  //CD->__print_comm__(ccomm); 

  string ciam = string(iam);
  CD->set_app_type(ciam);

  string cnamei = string(namei);
  CD->set_app_name(cnamei);
  CD->print_commij_names(); 


  MPI_Comm  local_comm = MPI_COMM_NULL; 
  int app_id = -1;
  int n_apps = -1; 
  CD->name_to_id(&app_id, &n_apps, &local_comm);
  //CD->__print_comm__(local_comm); 
  lcomm[0] = MPI_Comm_c2f(local_comm);
  
  
  CD->__create_interaction__(); 
  MPI_Barrier(local_comm); 
  CD->create_commij(local_comm);
  MPI_Barrier(local_comm); 


  //CD->print_commij_names(); 
  //cout<<"\n"; 
/*
  vector<MPI_Comm> Comms; 
  CD->get_ids_commij(Comms); 

  vector<string> Names; 
  CD->get_names_commij(Names); 

  vector<int> Roots; 
  CD->get_roots_commij(Roots); 


  //commij = new MPI_Fint[ Comms.size() ];
  for(int i=0; i<Comms.size(); i++) Ruts[i] = Roots[i];
  for(int i=0; i<Comms.size(); i++) commij[i] = MPI_Comm_c2f( Comms[i] );
  for(int i=0; i<Comms.size(); i++) memcpy( &nameij[i], Names[i].c_str(), strlen(Names[i].c_str()) );
  MPI_Barrier(local_comm); 
*/


/*
  CD.name_to_id(&app_id, &n_apps, &local_comm);
  CD.__create_interaction__(); //MPI_Barrier(local_comm); 
  CD.create_commij(local_comm);
*/

  delete CD; 
}



#ifdef __cplusplus
}
#endif

