#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ple_config_defs.h"
#include "ple_defs.h"
#include "ple_coupling.h"

/* NOTAS: 
  >>> 4.12.4. Transfer of Handles
  // C wrapper  
  void MPI_X_TYPE_COMMIT( MPI_Fint *f_handle, MPI_Fint *ierr) 
  { 
    MPI_Datatype datatype; 
    datatype = MPI_Type_f2c( *f_handle ); 
    *ierr = (MPI_Fint)MPI_Type_commit( &datatype ); 
    *f_handle = MPI_Type_c2f( datatype ); 
    return; 
  } 

  >> error dereferencing pointer to incomplete type
  Si en el archivo .h aparece únicamente:
    typedef struct _mitipo mitipo;
  No va a ser posible acceder a los miembros internos de esa estructura, ya que es un tipo incompleto.  

  "/libple/src/ple_coupling.h"
  // Opaque code coupling information structure 
  typedef struct _ple_coupling_mpi_set_t  ple_coupling_mpi_set_t;

  // Info structure for code coupling
  typedef struct {
    int          status;    // Status flag for synchronization info 
    int          root_rank; // Application root rank in MPI_COMM_WORLD 
    int          n_ranks;   // Number of ranks associated with application 
    const char  *app_type;  // Application type name (may be empty) 
    const char  *app_name;  // Application instance name (may be empty) 
  } ple_coupling_mpi_set_info_t;
*/

typedef struct 
{
  int       n_apps;       // Number of distinct applications 
  int       app_id;       // Id of the local application in the application info 
  int       app_names_l;  // Length of application names array 
  int      *app_info;     // For each application, 4 integers: 
                          // associated root in base_comm, n_ranks, and indexes in app_names 
  char     *app_names;    // Array of application type names and instance names 
  int      *app_status;   // Synchronization status for each application 
  double   *app_timestep; // Current time step for each application 
  MPI_Comm  base_comm;    // Handle to base communicator 
  MPI_Comm  app_comm;     // Handle to local application communicator 
} ple_coupling_mpi_set_t2;



void 
alya_coupling_mpi_intracomm_create_(MPI_Fint  *commi,
                                    int        rangei[2], 
                                    int        rangej[2], MPI_Fint  *fcommij) 
{
  MPI_Comm fcommi = MPI_COMM_NULL;
  fcommi = MPI_Comm_f2c( *commi );

  MPI_Comm commij = MPI_COMM_NULL;  
  ple_coupling_mpi_intracomm_create(MPI_COMM_WORLD,
                                    fcommi,
                                    rangej[0],
                                   &commij,
                                    rangei,  // <- local
                                    rangej); // <- distant

  fcommij[0] = MPI_Comm_c2f( commij );  

  printf("commi, rangei, rangej: %d, [%d,%d], [%d,%d], %d \n", commi[0], 
                                                           rangei[0], rangei[1]-1, 
                                                           rangej[0], rangej[1]-1, fcommij[0]); 

}


void
alya_coupling_mpi_name_to_id_(MPI_Fint *ccomm, char* app_name, int* app_id) 
{
  MPI_Comm fcomm;
  fcomm = MPI_Comm_f2c( *ccomm );


  app_id[0] = alya_coupling_mpi_name_to_id(fcomm, app_name); 
  //app_id[0] = ple_coupling_mpi_name_to_id(fcomm, app_name);
  if(app_id[0] < 0) app_id[0] = 0;

  printf("--| libple available!! \n");
}


void
alya_coupling_mpi_set_n_apps_(ple_coupling_mpi_set_t *cplng, int* napps) 
{
  napps[0] = ple_coupling_mpi_set_n_apps(cplng);
}


void
alya_coupling_mpi_set_get_app_id_(ple_coupling_mpi_set_t *cplng, int* app_id) 
{
  app_id[0] = ple_coupling_mpi_set_get_app_id(cplng);
}  


void 
alya_coupling_mpi_get_app_comm_(ple_coupling_mpi_set_t2  *cplng, 
                   MPI_Fint                 *fcomm)
{
  fcomm[0] = MPI_Comm_c2f( cplng->app_comm );  
}


void 
alya_coupling_mpi_get_base_comm_(ple_coupling_mpi_set_t2  *cplng, 
                   MPI_Fint                 *fcomm)
{
  fcomm[0] = MPI_Comm_c2f( cplng->base_comm );  
}


void
alya_coupling_mpi_set_get_timestep_(ple_coupling_mpi_set_t   *cplng, 
                                    float                    *fretval)
{
  int i; 
  int napps = ple_coupling_mpi_set_n_apps(cplng);
  const double *cretval = ple_coupling_mpi_set_get_timestep(cplng);
  for(i=0; i<napps; i++) fretval[i] = cretval[i]; 
/*
  printf("f1: "); 
  for(i=0; i<napps; i++) printf(" %f", retval[i]); 
  printf("\n"); 

  printf("c1: "); 
  for(i=0; i<napps; i++) 
  {
    printf(" %f", cretval[i]); 
    retval[i] = cretval[i]; 
  }
  printf("\n"); 
*/
}


void
alya_coupling_mpi_set_get_info_(ple_coupling_mpi_set_t       *cplng, 
                                int                          *app_id, 
                                ple_coupling_mpi_set_info_t  *info) 
{
  info[0] = ple_coupling_mpi_set_get_info(cplng, app_id[0]);
}


void
alya_coupling_mpi_set_destroy_(ple_coupling_mpi_set_t *cplng)
{
  ple_coupling_mpi_set_destroy(&cplng); 
}


//ple_coupling_mpi_set_t2 * <- my_change
void 
alya_coupling_mpi_set_create3_(int       *sync_flag,
                            const char  *app_type,
                            const char  *app_name,
                            MPI_Fint    *fbase_comm,
                            MPI_Fint    *fapp_comm, 
                            ple_coupling_mpi_set_t2 *s)
{

  int i, j;
  MPI_Status status;

  int set_rank = -1, app_rank = -1, n_app_ranks = 0;
  int root_marker = 0;
  int info_tag = 1, name_tag = 2;

  int counter[2] = {0, 0};
  int l_rank_info[5] = {-1, -1, -1, 1, 1}; // Status (1) + rank info (4) 
  int *rank_info = NULL;
  char *app_names = NULL;

printf("--| alya_coupling_mpi_set_create \n");


  sync_flag[0] = PLE_COUPLING_NO_SYNC;            // <- my_change
//  ple_coupling_mpi_set_t2 s;                 // <- my_change

  // Initialization

//  PLE_MALLOC(s2, 1, ple_coupling_mpi_set_t); // <- my_change

  s->n_apps = 0;
  s->app_id = -1;
  s->app_names_l = 0;
  s->app_info = NULL;
  s->app_names = NULL;

  // Initialization 
  MPI_Comm base_comm;                         // <- my_change
  base_comm = MPI_Comm_f2c( *fbase_comm );    // <- my_change

  MPI_Comm app_comm;                          // <- my_change
  app_comm = MPI_Comm_f2c( *fapp_comm );      // <- my_change

  s->base_comm = base_comm;
  s->app_comm = app_comm;

  MPI_Comm_rank(base_comm, &set_rank);

  if (app_comm != MPI_COMM_NULL) {
    MPI_Comm_rank(app_comm, &app_rank);
    MPI_Comm_size(app_comm, &n_app_ranks);
  }
  else {
    app_rank = 0;
    n_app_ranks = 1;
  }

  l_rank_info[0] = sync_flag[0] | PLE_COUPLING_INIT;
  l_rank_info[1] = set_rank;
  l_rank_info[2] = n_app_ranks;
  if (app_type != NULL)
    l_rank_info[3] = strlen(app_type) + 1;
  if (app_name != NULL)
    l_rank_info[4] = strlen(app_name) + 1;

  if (app_rank == 0)
    root_marker = 1;

  // Root rank of base_comm counts applications and receives info 

  MPI_Reduce(&root_marker, &(counter[0]), 1, MPI_INT, MPI_SUM, 0,
             base_comm);

//printf("--| %d \n", counter[0]); 

  // Root of base_comm collects all info 

  if (set_rank == 0) {

    int start = 0;

    PLE_MALLOC(rank_info, counter[0]*5, int);

    if (app_rank == 0) {
      for (i = 0; i < 5; i++)
        rank_info[i] = l_rank_info[i];
      start = 1;
    }

    // Use of different tags for info and strings is important
    // here as we use MPI_ANY_SOURCE and messages could be mixed 

    for (i = start; i < counter[0]; i++)
      MPI_Recv(rank_info + i*5, 5, MPI_INT, MPI_ANY_SOURCE, info_tag,
               base_comm, &status);

    // Convert rank_info count to index values 

    for (i = 0; i < counter[0]; i++)
      counter[1] += (rank_info[i*5 + 3] + rank_info[i*5 + 4]);

    PLE_MALLOC(app_names, counter[1], char);
    memset(app_names, 0, counter[1]);

    counter[1] = 0;

    if (app_rank == 0) {
      strcpy(app_names, app_type);
      if (app_name != NULL)
        strcpy(app_names + rank_info[3], app_name);
      else
        app_names[rank_info[3]] = '\0';
      counter[1] += (rank_info[3] + rank_info[4]);
      rank_info[4] = rank_info[3];
      rank_info[3] = 0;
    }

    for (i = start; i < counter[0]; i++) {
      int app_type_size = rank_info[i*5 + 3];
      int app_name_size = rank_info[i*5 + 4];
      int msg_len = app_type_size + app_name_size;
      rank_info[i*5 + 3] = counter[1];
      rank_info[i*5 + 4] = counter[1] + app_type_size;
      MPI_Recv(app_names + counter[1], msg_len, MPI_CHAR, rank_info[i*5 +1],
               name_tag, base_comm, &status);
      counter[1] += msg_len;
    }

  }

  // Other root ranks send info to root 

  else if (app_rank == 0) { // set_rank != 0

    char *sendbuf = NULL;
    int   sendbuf_l = l_rank_info[3] + l_rank_info[4];

    PLE_MALLOC(sendbuf, sendbuf_l, char);

    if (app_type != NULL)
      strcpy(sendbuf, app_type);
    else
      sendbuf[0] = '\0';
    if (app_name != NULL)
      strcpy(sendbuf + l_rank_info[3], app_name);
    else
      sendbuf[l_rank_info[3]] = '\0';

    MPI_Send(l_rank_info, 5, MPI_INT, 0, info_tag, base_comm);
    MPI_Send(sendbuf, sendbuf_l, MPI_CHAR, 0, name_tag, base_comm);

    PLE_FREE(sendbuf);
  }

  // Now root broadcasts application info 

  MPI_Bcast(counter, 2, MPI_INT, 0, base_comm);

  if (set_rank != 0) {
    PLE_MALLOC(rank_info, counter[0]*5, int);
    PLE_MALLOC(app_names, counter[1], char);
  }

  MPI_Bcast(rank_info, counter[0]*5, MPI_INT, 0, base_comm);
  MPI_Bcast(app_names, counter[1], MPI_CHAR, 0, base_comm);

  // Set global values 

  s->n_apps = counter[0];
  s->app_names_l = counter[1];
  s->app_names = app_names;

  PLE_MALLOC(s->app_info, s->n_apps*4, int);
  PLE_MALLOC(s->app_status, s->n_apps, int);
  PLE_MALLOC(s->app_timestep, s->n_apps, double);

  for (i = 0; i < s->n_apps && s->app_id < 0; i++) {
    for (j = 0; j < 4; j++)
      s->app_info[i*4 + j] = rank_info[i*5 + j+1];
    s->app_status[i] = rank_info[i*5];
    s->app_timestep[i] = -(i+1);
  }

  PLE_FREE(rank_info);

  // Set rank set to that of the application root for matching

  //MPI_Bcast(&set_rank, 1, MPI_INT, 1, app_comm); // ESTO DA PROBLEMA Y NO TENGO IDEA DE PORQUE!!

  for (i = 0; i < s->n_apps && s->app_id < 0; i++) {
    if (s->app_info[i*4] == set_rank)
      s->app_id = i;
  }
 
 
//  return s; // <- my_change  
}



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

//printf("a.01 \n");

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

  PLE_FREE(buf);

  if (eq_all == 1) {
    PLE_FREE(_group_name);
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

  PLE_FREE(_group_name);

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


/*
void
alya_ple_xyz_(MPI_Fint *base_comm, 
              MPI_Fint *app_comm, 
              char     *app_type, 
              char     *app_name, 
              ple_coupling_mpi_set_t *cplng) 
{
  int sync_flag = PLE_COUPLING_NO_SYNC;

  MPI_Comm fbase_comm;
  fbase_comm = MPI_Comm_f2c( *base_comm );

  MPI_Comm fapp_comm;
  fapp_comm = MPI_Comm_f2c( *app_comm );


//  ple_coupling_mpi_set_t *cplng2 = malloc( sizeof(_ple_coupling_mpi_set_t) );
//  ple_coupling_mpi_set_tcplng 0 = NULL; 
//  cplng[0].; 
  
  //cplng = ple_coupling_mpi_set_create(sync_flag,
  //                           app_type,
  //                           app_name,
  //                           fbase_comm,
  //                           fapp_comm); 
  
  ple_coupling_mpi_set_t *cplng1 = malloc(sizeof(cplng)); 
  cplng = 0; 
//  memcpy(cplng0, cplng, sizeof(cplng)); 
//  printf("%d %d \n", n_apps, app_id);

  int n_apps = ple_coupling_mpi_set_n_apps(cplng); // s->n_apps 
  int app_id = ple_coupling_mpi_set_get_app_id(cplng);
//  printf("%d %d \n", n_apps, app_id);

  ple_coupling_mpi_set_info_t info; 
  info = ple_coupling_mpi_set_get_info(cplng, app_id);
//  printf("%d %d %s\n", n_apps, app_id, info.app_type);
}
*/

/*
void 
alya_coupling_set_create_(   const char  *app_type,
                             const char  *app_name,
                             MPI_Fint    *base_comm,
                             MPI_Fint    *app_comm, 
                             int         *app_id, 
                             int         *range, 
                             int         *n_apps) 
{
  int sync_flag = PLE_COUPLING_NO_SYNC;

  MPI_Comm fbase_comm;
  fbase_comm = MPI_Comm_f2c( *base_comm );

  MPI_Comm fapp_comm;
  fapp_comm = MPI_Comm_f2c( *app_comm );

  ple_coupling_mpi_set_t *cplng = NULL; 
  
  cplng = ple_coupling_mpi_set_create(sync_flag,
                                      app_type, app_name,
                                      fbase_comm, fapp_comm); 

//  n_apps[0] = ple_coupling_mpi_set_n_apps(cplng);

  ple_coupling_mpi_set_info_t info; 
  
  int appi = app_id[0]; 
  info = ple_coupling_mpi_set_get_info(cplng, appi);

  range[0] = info.root_rank;
  range[1] = range[0] + info.n_ranks;
  
//printf("%d,%s) range: [%d,%d] \n", app_id[0], app_type, range[0], range[1]-1); 
//KEEP COMMENTED: ple_coupling_mpi_set_destroy(&cplng);
}
*/
