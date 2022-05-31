#include "commdom.hpp"


#ifndef __PLEPP_H__
#define __PLEPP_H__


class PLEPP
{
  public:
    string namei;
    string namej;
    string IAM; 

    CommDom  CD; // = CommDom();

    MPI_Comm  world_comm; 
    MPI_Comm  local_comm;
    MPI_Comm  commij;


    int        n_vertices_j; // = 0;
    double *vertex_coords_j; // = NULL;

    int     n_vertices_i   ; // = 0;
    int     n_elements_i   ; // = 0;
    double *vertex_coords_i; // = NULL;
    int    *vertex_num_i   ; // = NULL;

    void      create(int, char**); 
    void  reset_mesh(); 
    void    set_mesh();
    void create_mesh(float* dummy_ptr=NULL);

    void count_iters(); 
    void set_n_particles(int); 

    PLEPP(); 
   ~PLEPP(); 

    int n_iters; 
    int n_pts; 
    int dim; 

    bool    mesh_created; 
    bool locator_created; 
    bool           debug; 

}; 

#endif 
