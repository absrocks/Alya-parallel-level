#include <mpi.h>
#include <string.h>
/*!
* \file    alya_coupling_mpi_name_to_id.c
* \brief   PLE_WARP 
* \details taken from libple/src/ple_coupling.c 
* \autor   zavala-ake
* \autor   juan cajas 
*/

inline static void
_order_names_descend_tree(const char  *name[],
                          int          level,
                          const int    n_ents,
                          int          order[])
{
  int i_save, i1, i2, lv_cur;

  i_save = order[level];

  while (level <= (n_ents/2)) {

    lv_cur = (2*level) + 1;

    if (lv_cur < n_ents - 1) {

      i1 = order[lv_cur+1];
      i2 = order[lv_cur];

      if (strcmp(name[i1], name[i2]) > 0) lv_cur++;
    }

    if (lv_cur >= n_ents) break;

    i1 = i_save;
    i2 = order[lv_cur];

    if (strcmp(name[i1], name[i2]) >= 0) break;

    order[level] = order[lv_cur];
    level = lv_cur;
  }

  order[level] = i_save;
}


static void
_order_names(const char  *name[],
             int          order[],
             const int    n_ents)
{
  int i;
  int o_save;

  /* Initialize ordering array */

  for (i = 0 ; i < n_ents ; i++)
    order[i] = i;

  if (n_ents < 2)
    return;

  /* Create binary tree */

  i = (n_ents / 2) ;
  do {
    i--;
    _order_names_descend_tree(name, i, n_ents, order);
  } while (i > 0);

  /* Sort binary tree */

  for (i = n_ents - 1 ; i > 0 ; i--) {
    o_save   = order[0];
    order[0] = order[i];
    order[i] = o_save;
    _order_names_descend_tree(name, 0, i, order);
  }
}


void PLE_FREE_(void *ptr)
{
//  if(ptr != NULL) free(ptr);
}


int
alya_coupling_mpi_name_to_id(MPI_Comm     comm,
                            const char  *group_name)
{
  int i, eq_prev, eq_all;
  MPI_Status status;

  int l = 0, l_prev = 0;
  int rank_prev = MPI_PROC_NULL, rank_next = MPI_PROC_NULL;
  int rank_id = 0, n_ranks = 1, tag = 1;

  char *_group_name = NULL;
  char *buf = NULL, *names_buf = NULL;
  int *recv_count = NULL, *recv_displ = NULL, *app_id = NULL;

  int retval = -1;

  /* Initialization */

  if (_group_name == NULL) {
    l = strlen(group_name);
    //PLE_MALLOC(_group_name, l + 1, char); //  Allocate memory for _ni elements of type _type 
    _group_name = malloc((l + 1)*sizeof(char)); 
    strcpy(_group_name, group_name);
  }
  else {
    //PLE_MALLOC(_group_name, 1, char);
    _group_name = malloc(1*sizeof(char)); 
    _group_name[0] = '\0';
  }

  if (comm != MPI_COMM_NULL) {
    MPI_Comm_rank(comm, &rank_id);
    MPI_Comm_size(comm, &n_ranks);
    if (rank_id > 0)
      rank_prev = rank_id -1;
    if (rank_id + 1 < n_ranks)
      rank_next = rank_id + 1;
  }

  /* Check if multiple names are present using "light" algorithm */
  /*-------------------------------------------------------------*/

  if (rank_id %2 == 0) {
    MPI_Send(&l, 1, MPI_INT, rank_next, tag, comm);
    MPI_Recv(&l_prev, 1, MPI_INT, rank_prev, tag, comm, &status);
  }
  else {
    MPI_Recv(&l_prev, 1, MPI_INT, rank_prev, tag, comm, &status);
    MPI_Send(&l, 1, MPI_INT, rank_next, tag, comm);
  }

  //PLE_MALLOC(buf, l_prev + 1, char);
  buf = malloc((l_prev + 1)*sizeof(char)); 

  if (rank_id %2 == 0) {
    MPI_Send(_group_name, l, MPI_CHAR, rank_next, tag, comm);
    MPI_Recv(buf, l_prev, MPI_CHAR, rank_prev, tag, comm, &status);
  }
  else {
    MPI_Recv(buf, l_prev, MPI_CHAR, rank_prev, tag, comm, &status);
    MPI_Send(_group_name, l, MPI_CHAR, rank_next, tag, comm);
  }

  eq_prev = 1;
  if (rank_id > 0) {
    buf[l_prev] = '\0';
    if (strcmp(_group_name, buf))
      eq_prev = 0;
  }
  MPI_Allreduce(&eq_prev, &eq_all, 1, MPI_INT, MPI_MIN, comm);

  //PLE_FREE(buf);
  if(buf != NULL) free(buf);

  if (eq_all == 1) {
    //PLE_FREE(_group_name);
    if(_group_name != NULL) free(_group_name);
    return -1;
  }

  /* Now gather to rank 0 for full algorithm */
  /*-----------------------------------------*/

  if (rank_id == 0) {
    //PLE_MALLOC(recv_count, n_ranks, int);
    //PLE_MALLOC(recv_displ, n_ranks, int);
    recv_count = malloc((n_ranks)*sizeof(int)); 
    recv_displ = malloc((n_ranks)*sizeof(int)); 
  }

  MPI_Gather(&l, 1, MPI_INT, recv_count, 1, MPI_INT, 0, comm);

  if (rank_id == 0) {

    recv_displ[0] = 0;
    for (i = 1; i < n_ranks; i++)
      recv_displ[i] = recv_displ[i-1] + recv_count[i-1] + 1;
/*
    PLE_MALLOC(names_buf,
               recv_displ[n_ranks - 1] + recv_count[n_ranks - 1] + 1,
               char);
*/
    names_buf = malloc((recv_displ[n_ranks - 1] + recv_count[n_ranks - 1] + 1)*sizeof(char)); 
  }

  MPI_Gatherv(_group_name, l, MPI_CHAR,
              names_buf, recv_count, recv_displ, MPI_CHAR,
              0, comm);

  //PLE_FREE(_group_name);
  if(_group_name != NULL) free(_group_name);

  /* Order groups for rank 0 */

  if (rank_id == 0) {

    int n_apps = 1;
    int *order = NULL;
    char *name_prev = NULL;
    char **names_ptr = NULL;

    //PLE_MALLOC(names_ptr, n_ranks, char *);
    names_ptr = malloc((n_ranks)*sizeof(char*)); 
    
    for (i = 0; i < n_ranks; i++) {
      names_ptr[i] = names_buf + recv_displ[i];
      (names_ptr[i])[recv_count[i]] = '\0';
      recv_count[i] = -1;
    }

    /* Re-use arrays */
    order = recv_displ; recv_displ = NULL;
    app_id = recv_count; recv_count = NULL;

    _order_names((const char **)names_ptr, order, n_ranks);

    name_prev = names_ptr[order[0]];
    app_id[order[0]] = n_apps - 1;
    for (i = 1; i < n_ranks; i++) {
      if (strcmp(names_ptr[order[i]], name_prev)) {
        n_apps += 1;
        name_prev = names_ptr[order[i]];
      }
      app_id[order[i]] = n_apps - 1;
    }

    //PLE_FREE(names_ptr);
    //PLE_FREE(names_buf);
    //PLE_FREE(order);
     if(names_ptr != NULL) free(names_ptr);
     if(names_buf != NULL) free(names_buf);
     if(order     != NULL) free(order);
  }

  /* Now send app_id to all ranks */

  MPI_Scatter(app_id, 1, MPI_INT, &retval, 1, MPI_INT, 0, comm);

  //if (rank_id == 0) PLE_FREE(app_id);
  if (rank_id == 0) { if(app_id != NULL) free(app_id); }

  return retval;
}


void
alya_coupling_mpi_name_to_id_(MPI_Fint *ccomm, char* app_name, int* app_id)
{
  MPI_Comm fcomm;
  fcomm = MPI_Comm_f2c( *ccomm );

  app_id[0] = alya_coupling_mpi_name_to_id(fcomm, app_name);
  if(app_id[0] < 0) app_id[0] = 0;
}

