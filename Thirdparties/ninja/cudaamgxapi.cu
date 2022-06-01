#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cuda_runtime.h"
#include "amgx_c.h"

typedef size_t devptr_t;
FILE *file;
/* print callback (could be customized) */
void print_callback(const char *msg, int length)
{
  int rank;
  rank =0;
  if (rank == 0) fprintf(file,"%s", msg);
}


extern "C" void amgx_init_()
{
  file = fopen("AMGX.log","w");
  AMGX_SAFE_CALL(AMGX_initialize());
  AMGX_SAFE_CALL(AMGX_initialize_plugins());
  AMGX_SAFE_CALL(AMGX_install_signal_handler());
  AMGX_SAFE_CALL(AMGX_register_print_callback(&print_callback));
  return;
}

extern "C" void amgx_create_config_(devptr_t *config,char *file)
{
  char *tf = strstr(file,".json");
  if(tf == NULL)
    {
      printf("The AMGX confing file has to be a valid .json\n\n\n");
      exit(1);
    }
  int len = (size_t)(tf) - (size_t)(file);
  len = len + 5;
  char *name = (char*)malloc(len+1);
  strncpy(name,file,len);
  name[len] = '\0';
  AMGX_config_handle cfg;
  AMGX_SAFE_CALL(AMGX_config_create_from_file(&cfg,name));
  *config = (devptr_t)cfg;
  free(name);
  return;
}

extern "C" void amgx_destroy_config_(devptr_t *obj)
{
  AMGX_SAFE_CALL(AMGX_config_destroy((AMGX_config_handle)(*obj)));
  return;
}

extern "C" void amgx_destroy_resource_(devptr_t *obj)
{
  AMGX_SAFE_CALL(AMGX_resources_destroy((AMGX_resources_handle)(*obj)));
  return;
}

extern "C" void amgx_create_matrix_(devptr_t *mat, devptr_t *res)
{
  AMGX_matrix_handle A;
  AMGX_SAFE_CALL(AMGX_matrix_create(&A,(AMGX_resources_handle)(*res),AMGX_mode_dDDI));
  *mat = (devptr_t)A;
  return;
}

extern "C" void amgx_destroy_matrix_(devptr_t *obj)
{
  AMGX_SAFE_CALL(AMGX_matrix_destroy((AMGX_matrix_handle)(*obj)));
  return;
}

extern "C" void amgx_create_vector_(devptr_t *vec, devptr_t *res)
{
  AMGX_vector_handle vv;
  AMGX_SAFE_CALL(AMGX_vector_create(&vv,(AMGX_resources_handle)(*res),AMGX_mode_dDDI));
  *vec = (devptr_t)vv;
  return;
}

extern "C" void amgx_destroy_vector_(devptr_t *obj)
{
  AMGX_SAFE_CALL(AMGX_vector_destroy((AMGX_vector_handle)(*obj)));
  return;
}

extern "C" void amgx_create_solver_(devptr_t *sol, devptr_t *res,devptr_t *cfg)
{
  AMGX_solver_handle S;
  AMGX_SAFE_CALL(AMGX_solver_create(&S,(AMGX_resources_handle)(*res),AMGX_mode_dDDI,(AMGX_config_handle)(*cfg)));
  *sol = (devptr_t)S;
  return;
}

extern "C" void amgx_destroy_solver_(devptr_t *obj)
{
  AMGX_SAFE_CALL(AMGX_solver_destroy((AMGX_solver_handle)(*obj)));
  return;
}

extern "C" void amgx_upload_matrix_(devptr_t *mat,int *nd,int *nnzb,int *bx,int *by,devptr_t *row,devptr_t *col,devptr_t *values,devptr_t *diag)
{
  AMGX_SAFE_CALL(AMGX_matrix_upload_all((AMGX_matrix_handle)(*mat),*nd,*nnzb,*bx,*by,(int*)(*row),(int*)(*col),(double*)(*values),NULL));
  return;
}

extern "C" void amgx_upload_vector_(devptr_t *vec,int *nd,int *bx,devptr_t *val)
{
  AMGX_SAFE_CALL(AMGX_vector_upload((AMGX_vector_handle)(*vec),*nd,*bx,(double*)(*val)));
  return;
}

extern "C" void amgx_bind_vector_(devptr_t *vec, devptr_t *mat)
{
  AMGX_SAFE_CALL(AMGX_vector_bind((AMGX_vector_handle)(*vec),(AMGX_matrix_handle)(*mat)));
  return;
}

extern "C" void amgx_setup_solver_(devptr_t *sol,devptr_t *mat)
{
  AMGX_SAFE_CALL(AMGX_solver_setup((AMGX_solver_handle)(*sol),(AMGX_matrix_handle)(*mat)));
  return;
}

extern "C" void amgx_solve_solver_(devptr_t *sol,devptr_t *rhs,devptr_t *unk)
{
  AMGX_SOLVE_STATUS status;
  AMGX_SAFE_CALL(AMGX_solver_solve((AMGX_solver_handle)(*sol),(AMGX_vector_handle)(*rhs),(AMGX_vector_handle)(*unk)));

  AMGX_SAFE_CALL(AMGX_solver_get_status((AMGX_solver_handle)(*sol),&status));

  if(status != AMGX_SOLVE_SUCCESS)
    {
      printf("AMGX solver not success\n\n\n");exit(1);
    }
  return;
}

extern "C" void amgx_memory_pin_(double *mem,size_t bytes)
{
  AMGX_SAFE_CALL(AMGX_pin_memory(mem,bytes));
  return;
}

extern "C" void amgx_download_vector_(devptr_t *vec,devptr_t *val)
{
  AMGX_SAFE_CALL(AMGX_vector_download((AMGX_vector_handle)(*vec),(double*)(*val)));
  return;
}

extern "C" void amgx_solver_get_iterations_number_(devptr_t *handle,int *n)
{
  int iters;
  AMGX_SAFE_CALL(AMGX_solver_get_iterations_number((AMGX_solver_handle)(*handle),&iters));
  *n = iters;
  return;
}

extern "C" void amgx_solver_get_iteration_residual_(devptr_t *handle,int *n,int *idx,double *res)
{
  double residual;
  AMGX_SAFE_CALL(AMGX_solver_get_iteration_residual((AMGX_solver_handle)(*handle),*n,*idx,&residual));
  *res = residual;
  return;
}


extern "C" void amgx_matrix_comm_from_maps_one_ring_(devptr_t *mathand,int *num_nei,int *nei,int *sendsz, int *send_maps,int *recvsz,int *recv_maps)
{
  AMGX_SAFE_CALL(AMGX_matrix_comm_from_maps_one_ring((AMGX_matrix_handle)(*mathand), 1, *num_nei, nei, sendsz, (const int **)send_maps, recvsz, (const int **)recv_maps));
  return;
}

extern "C" void amgx_config_get_rings_(devptr_t *cfg,int *nrings)
{
  int nr;
  AMGX_SAFE_CALL(AMGX_config_get_default_number_of_rings((AMGX_config_handle)(*cfg), &nr));
  *nrings = nr;
  return;
}
