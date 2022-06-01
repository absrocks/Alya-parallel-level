/******************************************************************************
 * Copyright 1986, 2016 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief IndeX stream viewer main

#ifdef WIN32
# include <windows.h>
#undef min // remove the ugly defines in windows.h
#undef max
#endif  // WIN32

#ifdef LINUX
#include <signal.h>
#endif // LINUX

#ifdef NVINDEX_HAS_MPI_IPC_COMPUTE
#include <cuda_runtime.h>
#include <mpi.h>
#endif // NVINDEX_HAS_MPI_IPC_COMPUTE

#include <set>
#include <memory>
#include <map>
#include <sstream>

#include "common/common_utility.h"
#include "common/forwarding_logger.h"
#include "common/string_dict.h"
#include "../alya/alyaglobal.h"
#include "geostream_viewer.h"
#include "nvindex_appdata.h"
#include "config_utility.h"
#ifdef NVINDEX_HAS_MPI_IPC_COMPUTE
#include "cuda_ipc_volume_editing.h"
#include "distributed_volume_application.h"
#include "cuda_ipc_volume_editing_kernel.cuh"
#define INVALID_BRICK 0xFFFFFFFF
#endif // NVINDEX_HAS_MPI_IPC_COMPUTE

using namespace nv::index_common;

//----------------------------------------------------------------------
/// host rendering loop: test mode
///
/// Test mode:
///  - limited number of iterations
///
/// options:
///   app::run_mode = test
///      mode name. Fixed to 'test'.
///   test::max_iteration = INT_N
///      max number of iteration of the test mode. After N iteration, quit the process.
///   test::output_image_fname = OUTPUT_IMAGE_FILE_BASENAME
///      output snapshot image file basename 
///   test::snapshot_frame_list = int_F1 int_F2 int_F3 ... int_Fn
///      snapshot frame number list. The snapshot file name is 
///      'OUTPUT_IMAGE_FILE_BASENAME_%03d.ppm', %03d is replaced with the frame number.
///   e.g., 
///     test::output_image_fname = test_demo_viewer_
///     test::snapshot_frame_list = 5 9
///  
/// \return true when succeeded.


//GLOBAL variables for easy inter - portability between C++ and Fortran with ALYA
Nvindex_rendering_context irc;
Index_application index_app;
MPI_Comm insitucomm;
mi::math::Vector<mi::Uint32, 3> VOLUME_SIZE,BRICK_SIZE;
mi::Sint32 name_len,nb_global_ranks,cur_global_rank,nb_local_ranks,cur_local_rank;
mi::Uint32 nb_collected_bricks,nb_bricks_per_rank,nb_bricks;
std::vector<Brick_locality> local_brick_list;
int rendercount;

REGISTRY_SIM registry_sim;
#ifdef LINUX

//----------------------------------------------------------------------
// signal handler: prints signal info before calling default signal handler
static void signal_handler(int signum)
{
    const char* signame;
    switch (signum)
    {
      case SIGABRT:
          signame = "SIGABRT";
          break;
      case SIGINT:
          signame = "SIGINT";
          break;
      case SIGTERM:
          signame = "SIGTERM";
          break;
      default:
          signame = "unknown";
          break;
    }

    pid_t pid = getpid();
    printf("nvindex-viewer (pid %d): Caught signal %s\n", pid, signame);

   // Call default signal handler
   signal(signum, SIG_DFL);
   raise(signum);
}

#endif // LINUX

#ifdef NVINDEX_HAS_MPI_IPC_COMPUTE

//----------------------------------------------------------------------
/// main mpi_main
int mpi_main(int argc, char **argv,int grank, int gsize, int lrank, int lsize)
{
    int result = 0;
    
    // Init MPI and get host/gpu/rank information
    // ------------------------------------------

    name_len  = -1;
    nb_global_ranks = lsize;
    cur_global_rank = lrank;
    nb_local_ranks = lsize;
    cur_local_rank = lrank;
    char name[256];
    
    MPI_Comm_size(insitucomm, &nb_global_ranks);
    MPI_Comm_rank(insitucomm, &cur_global_rank);
    MPI_Get_processor_name(name,  &name_len);
    
    // Number of local ranks must match number of GPU
    mi::Sint32 nb_gpus = cuda::get_nb_gpus();
    //nb_local_ranks = nb_gpus;
    cur_local_rank = cur_global_rank%nb_local_ranks;

    // Number of local ranks must match number of GPU
    /*if(nb_local_ranks != nb_gpus)
    {
        INFO_LOG << "Number of GPUs: " << nb_gpus << ", " << "doesn't match number of local ranks: " 
            << nb_local_ranks;
            
        MPI_Finalize();
        exit(1);
	}*/
    
    // Number of ranks/gpus per host must be the same.
    if(nb_global_ranks%nb_local_ranks != 0)
    {
        INFO_LOG << "Number of GPU/Ranks per host is not the same";
            
        MPI_Finalize();
        exit(1);
    }

    mi::Sint32 nb_hosts = nb_global_ranks/nb_local_ranks;
    mi::Sint32 cur_host = cur_global_rank/nb_local_ranks;
    mi::Sint32 base_rank = (cur_global_rank/nb_local_ranks)*nb_local_ranks;
    
    INFO_LOG << "Node name(" << cur_host << "):" << name
             <<", Rank/Worldsize: " << cur_global_rank << "/" << nb_global_ranks
             << ", LocalRank/LocalSize: " << cur_local_rank << "/" << nb_local_ranks
             << ", Base rank: " << base_rank;
             
    // Make a communicator group only with ranks where IndeX instance are running (cur_local_rank == 0).
    // We want to collect in rank 0 (IndeX viewer) the node_id (using IndeX numeration) of all hosts.
    MPI_Group MPI_GROUP_WORLD, group_index_ranks, insitu_group;
    MPI_Comm comm_index_ranks;
    MPI_Comm insitu_comm;
    
    std::vector<mi::Sint32> index_ranks(nb_hosts);
    for(mi::Sint32 i=0; i<nb_hosts; i++)
        index_ranks[i] = i*nb_local_ranks;
    
    MPI_Comm_group(insitucomm, &MPI_GROUP_WORLD);
    MPI_Group_incl(MPI_GROUP_WORLD, nb_hosts, &index_ranks[0], &group_index_ranks);
    MPI_Comm_create(insitucomm, group_index_ranks, &comm_index_ranks);
    
    // Create application data:
    // ------------------------

    String_dict *p_app_dict = Nvindex_AppData::instance()->peek_app_proj();
    
    // Parse in application related and session related data from the
    // argument list and store it in the dictionary
    if(!index_app.process_command_line_arguments(argc, argv, *p_app_dict))
    {
        // Failure in command line or project file parsing
        MPI_Finalize();
        return 1;
    }

    
#ifdef LINUX
    // Catch common signals
    signal(SIGABRT, signal_handler);  // Triggered by assertions
    signal(SIGINT,  signal_handler);  // Ctrl+C
    signal(SIGTERM, signal_handler);  // 'kill' command
#endif // LINUX

    // We have to start an IndeX instance per host. We will use all cur_local_rank == 0 for that.
    // Rendering loop
    if(cur_local_rank == 0)
    {
        INFO_LOG << "Initializing Index at rank " << cur_global_rank;
        
        // index viewer
        if(cur_global_rank == 0)
        {
            // Setup
            if (!index_app.setup_viewer(irc, false, false))
            {
                ERROR_LOG << "Failed to set up IndeX service. Exiting.";
                MPI_Finalize();
                return 1;
            }
        }
        // index remote
	else
        {
            // Setup
            if(!index_app.setup_remote(irc))
            {
                ERROR_LOG << "Failed to set up IndeX remote service. Exiting.";
                MPI_Finalize();
                return 1;
            }
        }	
    }
    MPI_Barrier(insitucomm);
    return result;
}

extern "C" void index_viewer_shutdown_()
{
 
  // Shut down servers.
  index_app.shutdown_all_server(irc);
  irc.shutdown();
  index_app.shutdown();
  
  MPI_Finalize();
  return;
}

#endif // NVINDEX_HAS_MPI_IPC_COMPUTE
//----------------------------------------------------------------------
/// main
extern "C" int index_viewer_main_(int *grank, int *gsize, int *lrank,int *lsize,double *data,int *ND, int *V)
{
  int argc=0;
  std::string cmd,secondline;
  std::ifstream infile("VIZconfig.dat",std::ifstream::in);
  if(!infile.is_open())
    {
      std::cout << " Error opening VIZconfig.dat, needed for index" << std::endl;
    }
  std::getline(infile,cmd);
  std::getline(infile,secondline);
  std::istringstream iss2(secondline);
  iss2 >> secondline;
  iss2 >> rendercount;

  std::vector<char *> args;
  std::istringstream iss(cmd);
  std::string token;
  while(iss >> token) {
    char *arg = new char[token.size() + 1];
    copy(token.begin(), token.end(), arg);
    arg[token.size()] = '\0';
    args.push_back(arg);
    argc++;
  }
  args.push_back(0);

  
#ifdef NVINDEX_HAS_MPI_IPC_COMPUTE
  return mpi_main(argc,&args[0],*grank,*gsize,*lrank,*lsize);
#endif // NVINDEX_HAS_MPI_IPC_COMPUTE
}

extern "C" void commtoinsitu(MPI_Comm comm)
{
  MPI_Comm_dup(comm, &insitucomm);
  return;
}

extern "C" void update_framedata_(int *globalnum,double *data,int *ND, int *V)
{

  return;
}


extern "C" void render_frame_()
{
  for(int frame =0;frame<rendercount;frame++)
    //for(int frame =0;;)
    {
      if(cur_local_rank == 0)
	{
	  if(cur_global_rank == 0) // index viewer only
	    {
	      // check for end of application
	      if (Nvindex_AppData::instance()->is_app_run())
		{
		  //ERROR_LOG << "START RENDERING FRAME: " << frame;
		  index_app.render_frame(irc, frame, true);
		  //ERROR_LOG << "END RENDERING FRAME: " << frame;
		}
	      else
		{
		  if(Nvindex_AppData::instance()->is_any_video_stream_enabled())
		    {
		      // Shut down servers.
		      index_app.shutdown_all_server(irc);
		    }
		  
		  irc.shutdown();
		  index_app.shutdown();
		  
		  MPI_Finalize();
		  return;
		}
	    } 
	}
    }
}


