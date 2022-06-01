#include <cstring>
#include "ParticlesToCassandraConnector.h"
#include "mpi.h"

/*
 * WRITING_MODE goes
 * - 1 all worker by itself
 * - 2 intranode share
 * - 3 full share
 *
 */
#define WRITING_MODE 1
//#define SYNCH

#if WRITING_MODE == 2 //no share

//#define TRACE
#ifdef TRACE
#include "extrae.h"
#endif

int ts=1;
void MySetup::sendparticles_internal(int recvcounts, int nvar2, int kfl_posla_pts, int ilimi, double current_time,

                                     double *recvbuf_rp) {
    #ifdef TRACE
    Extrae_event(10000,ts);
    ts+=1;
    #endif
    if (!s.particle_writer) {
        s.setup_writer(kfl_posla_pts);

    }
     //}
#ifdef COLLECT_STATS
    // MEASURE TIME
    auto start = std::chrono::system_clock::now();
    uint32_t local_writes = 0;
#endif

    
     #ifdef COLLECT_STATS
    local_writes=recvcounts/nvar2;
    #endif

    for (int row = 0; row < recvcounts; row += nvar2) {
        char *keys = nullptr;
        char *values = nullptr;
        //if (recvbuf_rp[row + 13] > ilimi) continue;
        uint32_t offset = 0;

        //Allocate Keys
        keys = (char *) malloc(s.keys_size);
        offset = 0;
        // COPY KEYS
        std::memcpy(keys, recvbuf_rp + row, sizeof(int32_t)); //partId
        offset += sizeof(int32_t);
        std::memcpy(keys + offset, &current_time, sizeof(double)); //time
        //offset += sizeof(double);

        //Allocate Values
        values = (char *) malloc(s.values_size);
        offset = 0;


        // COPY VALUES
        // x, y, z, vx, vy, vz, ax, ay, az  + particle type
        size_t delta = sizeof(double) * 9 + sizeof(int32_t);
        std::memcpy(values + offset, &recvbuf_rp[row+1],delta);
        offset += delta;

        std::memcpy(values + offset, &worker_id, sizeof(int32_t));
        offset += sizeof(int32_t);


        if (kfl_posla_pts >= 3) {
            //  Cd, Stk,vfx, vfy, vfz, adx, ady, adz, aex, aey, aez, agx, agy, agz
            std::memcpy(values + offset, &recvbuf_rp[row + 14], sizeof(double) * 15);
            //offset += sizeof(double) * 9;
        }else if (kfl_posla_pts >= 2) {
            // Cd, Stk.vfx, vfy, vfz
            std::memcpy(values + offset, &recvbuf_rp[row + 14], sizeof(double) * 6);
        }else if (kfl_posla_pts >= 1) {
            // Cd, Stk1,2
            std::memcpy(values + offset, &recvbuf_rp[row + 14], sizeof(double) * 3);
        }

        // WRITE OPERATION
        try {
            s.particle_writer->write_to_cassandra(keys, values);
            // Deletion of keys and values performed inside write_to_cassandra
        }
        catch (std::exception &e) {
            std::cerr << "Can't insert data on table " << s.problem + "_particle" << std::endl;
            std::cerr << e.what() << std::endl;
        }

    }
    #ifdef SYNCH
    s.particle_writer->flush_elements();
    #endif

    #ifdef TRACE
    Extrae_event(10000,0);
    #endif

#ifdef COLLECT_STATS
    // MEASURE TIME
    auto end = std::chrono::system_clock::now();
    auto iter_time = (end - start);

    this->update_stats(iter_time, local_writes);
#endif
}

#elif  WRITING_MODE == 2 //intranode share

int shm_rank, shm_size, global_rank;
int name_len=-1;
int ts=1;
MPI_Comm shmcomm;
char processor_name[MPI_MAX_PROCESSOR_NAME];
void MySetup::sendparticles_internal(int recvcounts, int nvar2, int kfl_posla_pts, int ilimi, double current_time,
                                     double *input_data) {
    #ifdef TRACE
    Extrae_event(10000,ts);
    ++ts;
    #endif
    if (!s.particle_writer) {
        s.setup_writer(kfl_posla_pts);

    }
    if(name_len==-1) {
        MPI_Get_processor_name(processor_name, &name_len);
        MPI_Comm_rank(MPI_COMM_WORLD, &global_rank);
        MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm);
        MPI_Comm_rank(shmcomm, &shm_rank);
        MPI_Comm_size(shmcomm, &shm_size);
    }
    //}
#ifdef COLLECT_STATS
    // MEASURE TIME
    auto start = std::chrono::system_clock::now();
    uint32_t local_writes = 0;
#endif

    //std::cout <<  processor_name<< " global=>"<< global_rank << "local" << shm_rank<< " with recvcounts: " << recvcounts <<std::endl;
    int *rows_per_worker;
    rows_per_worker = (int *) malloc(shm_size * sizeof(int));
    MPI_Allgather(&recvcounts, 1, MPI_INT, rows_per_worker, 1, MPI_INT, shmcomm);

    int totalRows = 0;
    for (int i = 0; i < shm_size; i++) {
        totalRows += rows_per_worker[i];
    }

    if (totalRows == 0) {
        if (shm_rank == 0) {
            std::cout << processor_name << ": totalRows 0 " << std::endl;
        }
        #ifdef TRACE
        Extrae_event(10000,0);
        #endif
        return;

    }

    if (shm_rank == 0) {
        std::cout << processor_name << ": allocating " << totalRows << " rows, per node : " << " ";
        for (int i = 0; i < shm_size; i++) {
            std::cout << rows_per_worker[i]/nvar2 << " ";
        }
        std::cout << std::endl;

    }
    free(rows_per_worker);
    /*  if(recvcounts>0) {
        std::cout << processor_name << "." << shm_rank << " writes " << recvcounts << " has position at: "
                  << my_position << std::endl;
    }*/
    MPI_Info win_info;
    MPI_Info_create(&win_info);
    MPI_Info_set(win_info, "alloc_shared_noncontig", "false");
    MPI_Win win;
    double *recvbuf_rp = NULL;
    int chunk_size = (totalRows / nvar2 / shm_size)*nvar2;
    int start_segment = shm_rank * chunk_size;
    int end_segment;
    if (shm_rank == shm_size - 1) {
        end_segment = totalRows;
    } else {
        end_segment = start_segment + chunk_size;
    }


    // Shared memory
    if (MPI_SUCCESS !=
        MPI_Win_allocate_shared(recvcounts * sizeof(double), 1, win_info, shmcomm, &recvbuf_rp, &win)) {
        std::cout << processor_name << "." << shm_rank << " can't allocate shared window" << std::endl;
        return;
    }
    MPI_Info_free(&win_info);


    MPI_Win_lock_all(0, win);
    if (recvcounts > 0) {
        std::memcpy(recvbuf_rp, input_data,recvcounts * sizeof(double));
    }



    MPI_Win_sync(win);     /* memory fence to sync node exchanges */
    MPI_Barrier(shmcomm);  /* time barrier to make sure all ranks have updated their info */
    MPI_Win_sync(win);     /* additional memory fence maybe needed on some platforms */
    MPI_Win_unlock_all(win);

    //MPI_Barrier(shmcomm);
    // std::cout << processor_name << "." << shm_rank <<" writes from "<<start_segment<< " to "<< end_segment<<std::endl;

    MPI_Aint size;
    int disp;
    if (MPI_SUCCESS != MPI_Win_shared_query(win, MPI_PROC_NULL, &size, &disp, &recvbuf_rp)) {
        std::cout << processor_name << "." << shm_rank << " can't query the first node " << std::endl;
        return;
    }
    #ifdef COLLECT_STATS
    local_writes=(end_segment-start_segment)/nvar2;
    #endif

    for (int row = start_segment; row < end_segment; row += nvar2) {
        char *keys = nullptr;
        char *values = nullptr;
        //if (recvbuf_rp[row + 13] > ilimi) continue;
        uint32_t offset = 0;

        //Allocate Keys
        keys = (char *) malloc(s.keys_size);
        offset = 0;
        // COPY KEYS
        std::memcpy(keys, recvbuf_rp + row, sizeof(int32_t)); //partId
        offset += sizeof(int32_t);
        std::memcpy(keys + offset, &current_time, sizeof(double)); //time
        //offset += sizeof(double);

        //Allocate Values
        values = (char *) malloc(s.values_size);
        offset = 0;


        // COPY VALUES
        // x, y, z, vx, vy, vz, ax, ay, az  + particle type
        size_t delta = sizeof(double) * 9 + sizeof(int32_t);
        std::memcpy(values + offset, &recvbuf_rp[row+1],delta);
        offset += delta;

        std::memcpy(values + offset, &worker_id, sizeof(int32_t));
        offset += sizeof(int32_t);


        if (kfl_posla_pts >= 3) {
            //  Cd, Stk,vfx, vfy, vfz, adx, ady, adz, aex, aey, aez, agx, agy, agz
            std::memcpy(values + offset, &recvbuf_rp[row + 14], sizeof(double) * 15);
            //offset += sizeof(double) * 9;
        }else if (kfl_posla_pts >= 2) {
            // Cd, Stk.vfx, vfy, vfz
            std::memcpy(values + offset, &recvbuf_rp[row + 14], sizeof(double) * 6);
        }else if (kfl_posla_pts >= 1) {
            // Cd, Stk1,2
            std::memcpy(values + offset, &recvbuf_rp[row + 14], sizeof(double) * 3);
        }

        // WRITE OPERATION
        try {
            s.particle_writer->write_to_cassandra(keys, values);
            // Deletion of keys and values performed inside write_to_cassandra
        }
        catch (std::exception &e) {
            std::cerr << "Can't insert data on table " << s.problem + "_particle" << std::endl;
            std::cerr << e.what() << std::endl;
        }

    }
    #ifdef SYNCH
    s.particle_writer->flush_elements();
    #endif

    #ifdef TRACE
    Extrae_event(10000,0);
    #endif

    MPI_Win_free(&win);



#ifdef COLLECT_STATS
    // MEASURE TIME
    auto end = std::chrono::system_clock::now();
    auto iter_time = (end - start);

    this->update_stats(iter_time, local_writes);
#endif
}

#elif WRITING_MODE == 3 // full shuffle share

int  global_size, global_rank;
void MySetup::sendparticles_internal(int recvcounts, int nvar2, int kfl_posla_pts, int ilimi, double current_time,
                                     double *input_data) {
    #ifdef TRACE
    Extrae_event(10000,1);
    #endif
    if (!s.particle_writer) {
        s.setup_writer(kfl_posla_pts);

    }

    //}
#ifdef COLLECT_STATS
    // MEASURE TIME
    auto start = std::chrono::system_clock::now();
    uint32_t local_writes = 0;
#endif

    //std::cout <<  processor_name<< " global=>"<< global_rank << "local" << shm_rank<< " with recvcounts: " << recvcounts <<std::endl;
    int *double_per_worker;
    double_per_worker = (int *) malloc(global_size * sizeof(int));
    MPI_Allgather(&recvcounts, 1, MPI_INT, double_per_worker, 1, MPI_INT, MPI_COMM_WORLD);

    int totalRows = 0;
    for (int i = 0; i < global_size; i++) {
        totalRows += double_per_worker[i];
    }

    if (totalRows == 0) {
        if (global_rank == 0) {
            std::cout << "TotalRows 0 " << std::endl;
        }
        #ifdef TRACE
        Extrae_event(10000,0);
        #endif
        return;

    }

    if (global_rank == 0) {
        std::cout << "Allocating " << totalRows << " rows, per node : " << " ";
        for (int i = 0; i < global_size; i++) {
            std::cout << double_per_worker[i]/nvar2 << " ";
        }
        std::cout << std::endl;

    }

    double *recvbuf_rp = NULL;
    recvbuf_rp = (double *)(malloc(totalRows * sizeof(double)));


    MPI_Allgather(input_data,recvcounts,MPI_DOUBLE,
            recvbuf_rp,totalRows,MPI_DOUBLE,
            MPI_COMM_WORLD);


    int chunk_size = (totalRows / nvar2 / global_size)*nvar2;
    int start_segment = global_rank * chunk_size;
    int end_segment;
    if (global_rank == global_size - 1) {
        end_segment = totalRows;
    } else {
        end_segment = start_segment + chunk_size;
    }

    free(double_per_worker);

    #ifdef COLLECT_STATS
    local_writes=(end_segment-start_segment)/nvar2;
    #endif

    for (int row = start_segment; row < end_segment; row += nvar2) {
        char *keys = nullptr;
        char *values = nullptr;
        //if (recvbuf_rp[row + 13] > ilimi) continue;
        uint32_t offset = 0;

        //Allocate Keys
        keys = (char *) malloc(s.keys_size);
        offset = 0;
        // COPY KEYS
        std::memcpy(keys, recvbuf_rp + row, sizeof(int32_t)); //partId
        offset += sizeof(int32_t);
        std::memcpy(keys + offset, &current_time, sizeof(double)); //time
        //offset += sizeof(double);

        //Allocate Values
        values = (char *) malloc(s.values_size);
        offset = 0;


        // COPY VALUES
        // x, y, z, vx, vy, vz, ax, ay, az  + particle type
        size_t delta = sizeof(double) * 9 + sizeof(int32_t);
        std::memcpy(values + offset, &recvbuf_rp[row+1],delta);
        offset += delta;

        std::memcpy(values + offset, &worker_id, sizeof(int32_t));
        offset += sizeof(int32_t);


        if (kfl_posla_pts >= 3) {
            //  Cd, Stk,vfx, vfy, vfz, adx, ady, adz, aex, aey, aez, agx, agy, agz
            std::memcpy(values + offset, &recvbuf_rp[row + 14], sizeof(double) * 15);
            //offset += sizeof(double) * 9;
        }else if (kfl_posla_pts >= 2) {
            // Cd, Stk.vfx, vfy, vfz
            std::memcpy(values + offset, &recvbuf_rp[row + 14], sizeof(double) * 6);
        }else if (kfl_posla_pts >= 1) {
            // Cd, Stk1,2
            std::memcpy(values + offset, &recvbuf_rp[row + 14], sizeof(double) * 3);
        }

        // WRITE OPERATION
        try {
            s.particle_writer->write_to_cassandra(keys, values);
            // Deletion of keys and values performed inside write_to_cassandra
        }
        catch (std::exception &e) {
            std::cerr << "Can't insert data on table " << s.problem + "_particle" << std::endl;
            std::cerr << e.what() << std::endl;
        }

    }
    #ifdef SYNCH
    s.particle_writer->flush_elements();
    #endif
    free(recvbuf_rp);
    #ifdef TRACE
    Extrae_event(10000,0);
    #endif



#ifdef COLLECT_STATS
    // MEASURE TIME
    auto end = std::chrono::system_clock::now();
    auto iter_time = (end - start);

    this->update_stats(iter_time, local_writes);
#endif
}
#endif

void MySetup::sendsurfacedepos(double cutim, int subdomain, int sum_parts) {
    int32_t kfl_posla_pts = 1;
    if (!s.surfdepos_writer) s.setup_writer(kfl_posla_pts);
    int32_t offset = 0;
    char *keys = (char *) malloc(sizeof(double) + sizeof(int32_t));
    std::memcpy(keys + offset, &cutim, sizeof(double));
    offset += sizeof(double);

    int32_t sd = (int32_t) subdomain;
    std::memcpy(keys + offset, &sd, sizeof(int32_t));


    int32_t *values = (int32_t *) malloc(sizeof(int32_t));
    *values = (int32_t) sum_parts;
    try {
        s.surfdepos_writer->write_to_cassandra(keys, values);
    }
    catch (std::exception &e) {
        std::cerr << "Can't insert data on table " << s.problem + "_particle" << std::endl;
        std::cerr << e.what() << std::endl;
    }
}


void MySetup::writedeposition(double pos_tim, int partid, int parttype, double x, double y, double z, int kfl_exist,
                              int pos_set, double stk_1, double stk_2) {
/*
    int32_t kfl_posla_pts = 1;
    if (!s.depos_writer) s.setup_writer(kfl_posla_pts);


#ifdef COLLECT_STATS
    // MEASURE TIME
    auto start = std::chrono::system_clock::now();
#endif


    int32_t offset = 0;
    char *keys = (char *) malloc(sizeof(int32_t) + sizeof(double));
    std::memcpy(keys + offset, &partid, sizeof(int32_t));
    offset += sizeof(int32_t);

    std::memcpy(keys + offset, &pos_tim, sizeof(double));


    offset = 0;
    char *values = (char *) malloc(sizeof(double) * 5 + sizeof(int32_t) * 3);
    std::memcpy(values + offset, &parttype, sizeof(int32_t));
    offset += sizeof(int32_t);

    std::memcpy(values + offset, &x, sizeof(double));
    offset += sizeof(double);

    std::memcpy(values + offset, &y, sizeof(double));
    offset += sizeof(double);

    std::memcpy(values + offset, &z, sizeof(double));
    offset += sizeof(double);

    std::memcpy(values + offset, &pos_set, sizeof(int32_t));
    offset += sizeof(int32_t);

    std::memcpy(values + offset, &kfl_exist, sizeof(int32_t));
    offset += sizeof(int32_t);

    std::memcpy(values + offset, &stk_1, sizeof(double));
    offset += sizeof(double);

    std::memcpy(values + offset, &stk_2, sizeof(double));
    offset += sizeof(double);

    try {
        s.depos_writer->write_to_cassandra(keys, values);
    }
    catch (std::exception &e) {
        std::cerr << "Can't insert data on table " << s.problem + "_depositions" << std::endl;
        std::cerr << e.what() << std::endl;
    }


#ifdef COLLECT_STATS
    // MEASURE TIME
    auto end = std::chrono::system_clock::now();
    auto iter_time = (end - start);

    update_stats(iter_time, 1);
#endif 
*/
}


void MySetup::setup_writer(uint8_t kfl_posla_pts) {
    std::vector<config_map> keysnames, colsnames;


    keysnames = {{{"name", "partid"}},
                 {{"name", "time"}}};
    keys_size = sizeof(int32_t) + sizeof(double);

    colsnames = {{{"name", "x"}},
                 {{"name", "y"}},
                 {{"name", "z"}},
                 {{"name", "vx"}},
                 {{"name", "vy"}},
                 {{"name", "vz"}},
                 {{"name", "ax"}},
                 {{"name", "ay"}},
                 {{"name", "az"}},
                 {{"name", "parttype"}},
                 {{"name", "subdomain"}}};
    values_size = sizeof(double) * 9 + sizeof(int32_t) * 2;


    if (kfl_posla_pts >= 1) {
        values_size += sizeof(double) * 3; //Cd, stk
        std::vector<config_map> problem1 = {{{"name", "Cd"}},
                                            {{"name", "Stk_inst"}},
                                            {{"name", "Stk_eff"}}};

        colsnames.insert(colsnames.end(), problem1.begin(), problem1.end());
    }
    if (kfl_posla_pts >= 2) {
        values_size += sizeof(double) * 3; //fluid ccs
        std::vector<config_map> problem2 = {{{"name", "vfx"}},
                                            {{"name", "vfy"}},
                                            {{"name", "vfz"}}};
        colsnames.insert(colsnames.end(), problem2.begin(), problem2.end());
    }
    if (kfl_posla_pts >= 3) {
        values_size += sizeof(double) * 9; //3d ccs
        std::vector<config_map> problem3 = {{{"name", "adx"}},
                                            {{"name", "ady"}},
                                            {{"name", "adz"}},
                                            {{"name", "aex"}},
                                            {{"name", "aey"}},
                                            {{"name", "aez"}},
                                            {{"name", "agx"}},
                                            {{"name", "agy"}},
                                            {{"name", "agz"}}};
        colsnames.insert(colsnames.end(), problem3.begin(), problem3.end());
    }
    if (kfl_posla_pts >= 4) {
        // "Not implemented problem type 4"
        throw std::exception();
    }

    /*** setup config ***/
    config_map config = {{"writer_par",    std::to_string(writer_parallelism)},
                         {"writer_buffer", std::to_string(writer_queue)}};

    std::string table_name = problem + "_particle";
    try {
        particle_writer = SI->make_writer(table_name.c_str(), keyspace.c_str(), keysnames, colsnames, config);
    }
    catch (std::exception &e) {
        std::cerr << "Can't create table " << table_name << " on keyspace: " << keyspace << std::endl;
        std::cerr << e.what() << std::endl;
        throw e;
    }


    /*** SURFACE DEPOSITIONS

    keysnames = {{{"name", "time"}},{{"name", "subdomain"}}};

    colsnames = {{{"name", "depositions"}}};

    table_name = problem + "_surfdepos";
    try {
        surfdepos_writer = SI->make_writer(table_name.c_str(), keyspace.c_str(), keysnames, colsnames, config);
    }
    catch (std::exception &e) {
        std::cerr << "Can't create table " << table_name << " on keyspace: " << keyspace << std::endl;
        std::cerr << e.what() << std::endl;
        throw e;
    }

    ***/



    /*** DEPOSITIONS ***/


    keysnames = {{{"name", "partid"}},
                 {{"name", "time"}}};

    colsnames = {{{"name", "parttype"}},
                 {{"name", "x"}},
                 {{"name", "y"}},
                 {{"name", "z"}},
                 {{"name", "pos_set"}},
                 {{"name", "kfl_exist"}},
                 {{"name", "stk_inst"}},
                 {{"name", "stk_eff"}}};

    table_name = problem + "_depositions";
    try {
        depos_writer = SI->make_writer(table_name.c_str(), keyspace.c_str(), keysnames, colsnames, config);
    }
    catch (std::exception &e) {
        std::cerr << "Can't create table " << table_name << " on keyspace: " << keyspace << std::endl;
        std::cerr << e.what() << std::endl;
        throw e;
    }
}
