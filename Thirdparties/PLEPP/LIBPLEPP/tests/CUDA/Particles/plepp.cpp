#include "plepp.h"

PLEPP::PLEPP()
{
  n_iters = 0;
  n_pts   = 0; 

     mesh_created = false; 
  locator_created = false; 
            debug = false; 

  if(debug) cout<<"[PLEPP]"<<"\n"; 
}


PLEPP::~PLEPP()
{

  if(debug) cout<<"[~PLEPP]"<<"\n";

}


void 
PLEPP::create(int argc, char **argv)
{ 
  IAM   = "PARTICLES";
  namei = "";
  namej = "";

//  MPI_Init(NULL, NULL);
  world_comm = MPI_COMM_WORLD;
  local_comm = MPI_COMM_NULL;
  commij     = MPI_COMM_NULL; 


  if(argc==2)      namei = argv[1];
  printf("\'%s\'", namei.c_str() );

  if(namei=="SOLID") namej = "FLUID";
  if(namei=="FLUID") namej = "SOLID";

  CD = CommDom();
  CD.init();
  CD.set_app_type(IAM);
  CD.set_world_comm(world_comm);
  CD.set_app_name(namei);


  int  app_id = -1;
  int  n_apps = -1;
  CD.name_to_id(&app_id, &n_apps, &local_comm);
  CD.__create_interaction__(); //MPI_Barrier(local_comm);
  CD.create_commij(local_comm);
  MPI_Barrier(local_comm);

  int commij_size = -1;
  CD.get_commij_size(&commij_size); 
    
  if(commij_size>0)
  { 
    CD.get_commij(namej, &commij);
//    MPI_Barrier(commij);
  }

  if(debug) cout<<"[ '"<<namei<<"' PLEPP.create]"<<"\n";

}


void
PLEPP::reset_mesh() 
{

  if( (mesh_created)&&(commij != MPI_COMM_NULL) )
  {
     if( vertex_coords_j )
     {
       delete[] vertex_coords_j; 
     }

     vertex_coords_j = NULL;
        n_vertices_j = 0;

     n_vertices_i    = 0;
     n_elements_i    = 0;
     vertex_coords_i = NULL;
     vertex_num_i    = NULL;

     mesh_created = false; 
  }

  if(false) //&&locator_created)
  {
    CD.locator_destroy();  
    locator_created = false; 
  }

  if(debug) cout<<"[ '"<<namei<<"' PLEPP.reset_mesh]"<<"\n";

}


void
PLEPP::set_mesh()
{

  if( (!mesh_created)&&(commij != MPI_COMM_NULL) )
  {
       n_vertices_j = n_pts;
               dim  = 3; 
    vertex_coords_j = new double[n_pts*dim];
    for(int i=0; i<n_pts*dim; i++) vertex_coords_j[i] = 0.0;  

    mesh_created = true;

    if( !locator_created )
    {
       CD.locator_create2(local_comm, commij, 1e-3);
       locator_created = true;
    }
  }

  if(debug) cout<<"[ '"<<namei<<"' PLEPP.set_mesh]"<<"\n";

}


void
PLEPP::create_mesh(float* pts)
{

  if( (mesh_created)&&(commij != MPI_COMM_NULL) )
  {
    for(int i=0; i<=n_pts; i++) 
    {
      vertex_coords_j[i*3+0] =  pts[i*4+0];  
      vertex_coords_j[i*3+1] =  pts[i*4+1];   
      vertex_coords_j[i*3+2] =  pts[i*4+2];   
/*
cout<< vertex_coords_j[i*3+0] <<" "; 
cout<< vertex_coords_j[i*3+1] <<" "; 
cout<< vertex_coords_j[i*3+2] <<" "; 
cout<<"\n";
*/
    }

    if( locator_created )
    {
  
      CD.locator_set_mesh(n_vertices_i,
                          n_elements_i,
                          vertex_coords_i,
                          vertex_num_i,
                          n_vertices_j,
                          vertex_coords_j);

      //CD.save_dist_coords(local_rank);
    } 
  }

  if(debug) cout<<"[ '"<<namei<<"' PLEPP.create_mesh]"<<"\n";

}



void
PLEPP::count_iters()  
{
  n_iters += 1;
  if(debug) cout<<"n_iters: "<< n_iters <<" \n"; 
}


void 
PLEPP::set_n_particles( int __n_pts__ )
{
  n_pts = __n_pts__; 
}



