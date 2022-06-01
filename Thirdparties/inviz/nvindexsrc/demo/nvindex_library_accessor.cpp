/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

/**
 * In this example, it will use IndeX library
 * There are two possibilities: is an argument is given, the a scene is loaded
 * if nothing is given, then a scene is created.
 */

#include "nvindex_library_accessor.h"

// common
#include "common/affinity_information.h"
#include "common/cluster_change_callback.h"
#include "common/distributed_heightfield_elevation_change.h"
#include "common/distributed_voxel_value_change.h"
#include "common/distributed_compute_techniques.h"
#include "common/forwarding_logger.h"
#include "common/scene_parameter_transfer.h"
#include "common/statistics_transfer.h"
#include "common/string_dict.h"

// triangle mesh importer
#include "common/triangle_mesh_importer.h"
#include "common/synthetic_reservoir_patch_generator.h"
#include "common/reservoir_grid_importer.h"

// height field data importer/generator
#include "common/raw_heightfield_data_importer.h"
#include "common/raw_heightfield_data_importer_ssv.h" // FIXME: obsolete
#include "common/synthetic_heightfield_generator.h"
#include "common/heightfield_skeleton_generator.h"
#include "common/ppm_heightfield_importer.h"

// volume data importer/generator
#include "common/multi_attribute_scaling_sequence_importer.h"
#include "common/raw_volume_data_importer.h"
#include "common/raw_volume_data_sequence_importer.h"
#include "common/repeat_raw_volume_data_importer.h"
#include "common/rgba_raw_volume_data_importer.h"
#include "common/single_value_raw_volume_data_importer.h"
#include "common/short_value_raw_volume_data_importer.h"
#include "common/synthetic_volume_generator.h"
#include "common/volume_replica_generation.h" // the importers/generators in <common> should go to <importers>

// irregular volume importers
#include "common/irregular_volume_importer.h"
#include "common/irregular_volume_importer_vtk.h"

#ifdef NV_IDX_USE_OPENVDB_IMPORTER
// sparse volume importers
#include "common/sparse_volume_importer.h"
#endif

#ifdef NVINDEX_HAS_MPI_IPC_COMPUTE
#include "cuda_ipc_volume_editing.h"
#endif // NVINDEX_HAS_MPI_IPC_COMPUTE

#ifdef NVINDEX_HAS_COMPUTE
#include "compute/compute_wrapper.h"
#endif

#include "geostream_viewer.h"

#include "autotracking_workflow_operation.h"
#include "bridge_server.h"
#include "config_utility.h"
#include "distributed_amplitude_map_job.h"
#include "distributed_command_execution.h"
#include "distributed_heightfield_elevation_delete.h"
#include "distributed_heightfield_seeding_updates.h"
#include "distributed_seed_point_evaluation.h"
#include "distributed_volume_bordermap_generation.h"
#include "distributed_volume_filter.h"
#include "heightfield_data_retrieval_appjob.h"
#include "irregular_volume_data_editing.h"
#include "simple_gridding_operation.h"
#include "volume_bordermap_database_element.h"
#include "volume_brick_element.h"
#include "volume_data_export.h"
#include "volume_data_retrieval_appjob.h"

#ifndef MI_PLATFORM_WINDOWS
#include <dlfcn.h>
#include <netinet/in.h>
#include <sys/time.h>
#else
#include <Winsock2.h>
#include <windows.h>
#endif  // MI_PLATFORM_WINDOWS

#ifdef LINUX
#include <unistd.h>
#endif  // LINUX

#include <time.h>
#include <stdlib.h>
#include <cassert>
#include <string>

#include <nv/index/iindex_debug_configuration.h>

//----------------------------------------------------------------------
Nvindex_library_accessor::Nvindex_library_accessor()
  : m_nvindex_interface(0)
{
    // empty
}

//----------------------------------------------------------------------
Nvindex_library_accessor::~Nvindex_library_accessor()
{
    clear();

    // Shutting down the nvindex library
    if (m_nvindex_interface.is_valid_interface()) {
        m_nvindex_interface->shutdown();
        m_nvindex_interface = NULL;
    }
}

//----------------------------------------------------------------------
bool Nvindex_library_accessor::initialize(
    const nv::index_common::String_dict& app_project,
    Nvindex_rendering_context& irc)
{
    if (is_initialized())
    {
        {
            ERROR_LOG << "Nvindex_library_accessor::initialize: double initialization.";
        }
        exit(3);
    }

    // Access nvindex library interface
    m_nvindex_interface = access_index_interface("libnvindex");
    if (m_nvindex_interface == NULL)
    {
        ERROR_LOG << "Could not acquire IIndex interface.";
        return false;
    }

    mi::base::Handle<mi::neuraylib::INeuray> dice(
        m_nvindex_interface->get_dice_interface());
    if (!dice.is_valid_interface())
    {
        ERROR_LOG << "Could not acquire DiCE interface.";
        return false;
    }

    INFO_LOG << get_interface()->get_product_name()
             << " version " <<  get_interface()->get_version()
             << ", build " << get_interface()->get_revision();

    const mi::Uint32 nvindex_api_version = get_interface()->get_interface_version();
    const mi::Uint32 expected_api_version = 2;
    if (nvindex_api_version != expected_api_version)
    {
        ERROR_LOG << "This reference viewer is based on " << get_interface()->get_product_name()
                  << " interface version " << expected_api_version << ", "
                  << "but only " << nvindex_api_version << " is available. Exiting with exit code 4.";
        exit(4);
    }

    configure_log_module(m_nvindex_interface, app_project);

    // DiCE debug configuration options
    {
        mi::base::Handle<mi::neuraylib::IDebug_configuration> debug_configuration(
            get_interface()->get_api_component<mi::neuraylib::IDebug_configuration>());
        assert(debug_configuration.is_valid_interface());

        // Configure networking.
        if (set_config_dice_network_setting(debug_configuration, app_project) < 0)
        {
            ERROR_LOG << "Some DiCE debug network configuration failed. Please check your project file.";
        }

        std::vector< std::pair< std::string, std::string > > ret_opt_vec =
            get_configuration_vector("dice::debug_configuration::", app_project, true);

        for(size_t i = 0; i < ret_opt_vec.size(); ++i)
        {
            const std::string opt = ret_opt_vec[i].first + "=" + ret_opt_vec[i].second;
            INFO_LOG << "Extra DiCE debug configuration: [" << opt << "]";
            mi::Sint32 ret = debug_configuration->set_option(opt.c_str());
            if(ret != 0)
            {
                ERROR_LOG << "Failed to set DiCE debug configuration: [" << opt << "]";
            }
        }
    }

    // IndeX debug configuration options
    {
        std::vector< std::pair< std::string, std::string > > ret_opt_vec =
            get_configuration_vector("index::debug_configuration::", app_project, true);
        
        mi::base::Handle<nv::index::IIndex_debug_configuration> index_debug_config(
            m_nvindex_interface->get_api_component<nv::index::IIndex_debug_configuration>());
        assert(index_debug_config.is_valid_interface());

        for(size_t i = 0; i < ret_opt_vec.size(); ++i)
        {
            const std::string opt = ret_opt_vec[i].first + "=" + ret_opt_vec[i].second;
            INFO_LOG << "Extra IndeX debug configuration: [" << opt << "]";
            mi::Sint32 ret = index_debug_config->set_option(opt.c_str());
            if(ret != 0)
            {
                ERROR_LOG << "Failed to set IndeX debug configuration: [" << opt << "]";
            }
        }
    }

    // Configure scheduling
    //mi::base::Handle<mi::neuraylib::IScheduling_configuration> ischeduling_configuration(
    //    get_interface()->get_api_component<mi::neuraylib::IScheduling_configuration>());
    //assert(ischeduling_configuration.is_valid_interface());
    //ischeduling_configuration->set_thread_pool_size(256);
    //ischeduling_configuration->set_max_number_of_active_threads(256);

    // Configure networking
    mi::base::Handle<mi::neuraylib::INetwork_configuration> inetwork_configuration(
        get_interface()->get_api_component<mi::neuraylib::INetwork_configuration>());
    assert(inetwork_configuration.is_valid_interface());

    if (!app_project.is_defined("dice::network::mode"))
    {
        ERROR_LOG << "NETWORK: No 'dice::network::mode' key found, cannot configure network, networking disabled.";
    }

    bool enable_networking = true;
    std::string netmode = app_project.get("dice::network::mode", "OFF");
    if (netmode == "OFF")
    {
        enable_networking = false;
    }
    else if (netmode == "UDP")
    {
        INFO_LOG << "NETWORK: Enabling UDP networking mode.";
        inetwork_configuration->set_mode(mi::neuraylib::INetwork_configuration::MODE_UDP);
    }
    else if (netmode == "TCP")
    {
        WARN_LOG << "NETWORK: TCP networking mode is no longer supported."; // backward compatibility. Can be removed.
    }
    else if (netmode == "TCP_WITH_DISCOVERY")
    {
        INFO_LOG << "NETWORK: Enabling TCP (with discovery) networking mode.";
        inetwork_configuration->set_mode(mi::neuraylib::INetwork_configuration::MODE_TCP_WITH_DISCOVERY);
    }
    else if (netmode == "UDP_WITH_DISCOVERY")
    {
        WARN_LOG << "NETWORK: UDP_WITH_DISCOVERY networking mode is no longer supported."; // backward compatibility. Can be removed.
    }
    else if (netmode == "")
    {
        WARN_LOG << "NETWORK: Networking mode is not specified. Networking disabled. "
                 << "To specify the mode, add 'dice::network::mode = mode' line in the project file.";
        enable_networking = false;
    }
    else
    {
        WARN_LOG << "NETWORK: Unknown networking mode '" << netmode << "', networking disabled."
                 << "Please check the 'dice::network::mode = mode' line in the project file.";
        enable_networking = false;
    }

    // Configure if network is on
    if (enable_networking)
    {
        if (netmode == "UDP")
        {
            // Empty of boolean (one digit), then do not set the multicast address
            std::string const multicast_addr = app_project.get("dice::network::multicast_address", "");
            if (multicast_addr.size() > 1) {
                // Mutilcast address exists, set it
                INFO_LOG << "NETWORK: UDP network is configured with multicast address ["
                         << multicast_addr << "].";
                inetwork_configuration->set_multicast_address(multicast_addr.c_str());
            }
            else
            {
                WARN_LOG << "NETWORK: UDP network is defined without multicast address. "
                         << "Using DiCE default multicast address ["
                         << inetwork_configuration->get_multicast_address()->get_c_str() << "].";
            }
        }
        else if (netmode == "TCP_WITH_DISCOVERY" )
        {
            // TCP_WITH_DISCOVERY means: communication is done by TCP,
            // but the host discovery is done by UDP multicast or
            // unicast. TCP_WITH_DISCOVERY uses dice::network::discovery_address for
            // discovery. See the dice documentation.

            std::string const discovery_addr = app_project.get("dice::network::discovery_address", "");
            // Empty of boolean (one digit), then do not set the multicast address
            if (discovery_addr.size() > 1) {
                // Discovery address exists, set it
                INFO_LOG << "NETWORK: TCP_WITH_DISCOVERY network is configured with discovery address ["
                         << discovery_addr << "].";
                inetwork_configuration->set_discovery_address(discovery_addr.c_str());
            }
            else
            {
                WARN_LOG << "NETWORK: TCP_WITH_DISCOVERY network is defined without discovery address. "
                         << "Using DiCE default discovery address.";
            }
        }

        // Compression level
        mi::Uint32 const comp_lv = nv::index_common::get_uint32(app_project.get("dice::network::compression_level", "0"));
        if (inetwork_configuration->set_compression_level(comp_lv) != 0)
        {
            ERROR_LOG << "NETWORK: Failed to set_compression_level to " << comp_lv << ".";
        }

        // Cluster interface
        std::string const cluster_if_addr = app_project.get("dice::network::cluster_interface_address", "");
        if (!cluster_if_addr.empty())
        {
            INFO_LOG << "NETWORK: Set cluster interface [" << cluster_if_addr << "]";
            const mi::Sint32 ret_stat =
                inetwork_configuration->set_cluster_interface(cluster_if_addr.c_str());
            if (ret_stat != 0)
            {
                ERROR_LOG << "NETWORK: Failed to set_cluster_interface with [" << cluster_if_addr
                          << "] code: " << ret_stat;
            }
        }

        // Redundancy level of DiCE database
        std::string const redundancy_level = app_project.get("dice::network::redundancy_level", "");
        if (!redundancy_level.empty())
        {
            if (inetwork_configuration->set_redundancy_level(nv::index_common::get_uint32(redundancy_level)) == 0)
            {
                INFO_LOG << "NETWORK: Redundancy level set to " << redundancy_level;
            }
            else
            {
                ERROR_LOG << "NETWORK: Failed to set redundancy level to " << redundancy_level
                          << ", using value " << inetwork_configuration->get_redundancy_level() << " instead "
                          << "(value should be between 1 and 4)";
            }
        }

        // Set up usage of RDMA Infiniband in the networking layer
        const std::string use_rdma = app_project.get("dice::network::use_rdma", "");
        if (!use_rdma.empty())
        {
            const bool is_use_rdma = nv::index_common::get_bool(use_rdma);
            if (inetwork_configuration->set_use_rdma(is_use_rdma) == 0)
            {
                INFO_LOG << "NETWORK: " << (is_use_rdma ? "Enabling" : "Disabling") << " RDMA";
            }
            else
            {
                ERROR_LOG << "NETWORK: Failed to set_use_rdma to " << (is_use_rdma ? "enabled" : "disabled");
            }
        }
    }

    mi::base::Handle<mi::neuraylib::IGeneral_configuration> igeneral_configuration(
        get_interface()->get_api_component<mi::neuraylib::IGeneral_configuration>());
    assert(igeneral_configuration.is_valid_interface());

    // Activate admin HTTP interface if requested
    const std::string admin_http_port = app_project.get("dice::admin_http_port", "");
    if (!admin_http_port.empty())
    {
        const std::string admin_serv_param = app_project.get("dice::admin_http_listen") + std::string(":") + admin_http_port;
        igeneral_configuration->set_admin_http_address(admin_serv_param.c_str());
    }

    // DiCE-specific configuration
    mi::base::Handle<mi::neuraylib::IDice_configuration> dice_configuration(
        get_interface()->get_api_component<mi::neuraylib::IDice_configuration>());
    assert(dice_configuration.is_valid_interface());

    // Register classes for serialization
    bool is_registered = false;
    is_registered = get_interface()->register_serializable_class<Volume_brick_element>();
    assert(is_registered);
    nv::index_common::no_unused_variable_warning_please(is_registered);

    // Register serializable classes    
    is_registered = get_interface()->register_serializable_class<nv::index_common::Affinity_information>();
    assert(is_registered);
    is_registered = get_interface()->register_serializable_class<nv::index_common::Domain_specific_subdivision>();
    assert(is_registered);

    // Register Bridge client/server classes for serialization
    is_registered = get_interface()->register_serializable_class<nv::index_common::Scene_parameter_transfer>();
    assert(is_registered);
    is_registered = get_interface()->register_serializable_class<nv::index_common::Statistics_transfer>();
    assert(is_registered);


    // Register volume retrieve job for application
    is_registered = get_interface()->register_serializable_class<Volume_data_retrieval_appjob>();
    assert(is_registered);

    // ------------------------------------------------------------------------------------------------------------
    // Register volume data importers and generators
    is_registered = get_interface()->register_serializable_class<nv::index_common::Multi_attribute_scaling_sequence_importer>();
    assert(is_registered);
    is_registered = get_interface()->register_serializable_class<nv::index_common::Raw_volume_data_importer>();
    assert(is_registered);
    is_registered = get_interface()->register_serializable_class<nv::index_common::Raw_volume_data_sequence_importer>();
    assert(is_registered);
    is_registered = get_interface()->register_serializable_class<nv::index_common::Repeat_raw_volume_data_importer>();
    assert(is_registered);
    is_registered = get_interface()->register_serializable_class<nv::index_common::Rgba_raw_volume_data_importer>();
    assert(is_registered);
    is_registered = get_interface()->register_serializable_class<nv::index_common::Single_value_raw_volume_data_importer>();
    assert(is_registered);
    is_registered = get_interface()->register_serializable_class<nv::index_common::Short_value_raw_volume_data_importer>();
    assert(is_registered);
    is_registered = get_interface()->register_serializable_class<nv::index_common::Synthetic_volume_generator>();
    assert(is_registered);
    is_registered = get_interface()->register_serializable_class<nv::index_common::Volume_replica_generation>();
    assert(is_registered);
    // ------------------------------------------------------------------------------------------------------------

    // ------------------------------------------------------------------------------------------------------------
    // Register height field data importers and generators
    is_registered = get_interface()->register_serializable_class<nv::index_common::Raw_heightfield_data_importer>();
    assert(is_registered);
    is_registered = get_interface()->register_serializable_class<nv::index_common::Raw_heightfield_data_importer_ssv>();
    assert(is_registered);
    is_registered = get_interface()->register_serializable_class<nv::index_common::Synthetic_heightfield_generator>();
    assert(is_registered);
    is_registered = get_interface()->register_serializable_class<nv::index_common::Heightfield_skeleton_generator>();
    assert(is_registered);
    is_registered = get_interface()->register_serializable_class<nv::index_common::PPM_heightfield_importer>();
    assert(is_registered);
    // ------------------------------------------------------------------------------------------------------------

    // ------------------------------------------------------------------------------------------------------------
    // Register triangle mesh data importers
    is_registered = get_interface()->register_serializable_class<nv::index_common::Triangle_mesh_importer>();
    assert(is_registered);
    // Register reservoir grid synthetic generator configuration and the respective job
    is_registered = get_interface()->register_serializable_class<nv::index_common::Synthetic_reservoir_patch_generator>();
    assert(is_registered);
    // Register reservoir grid importer configuration and the respective job
    is_registered = get_interface()->register_serializable_class<nv::index_common::Reservoir_grid_importer>();
    assert(is_registered);
    // ------------------------------------------------------------------------------------------------------------

    // ------------------------------------------------------------------------------------------------------------
    // Irregular volume data importers
    is_registered = get_interface()->register_serializable_class<nv::index_common::Irregular_volume_importer>();
    assert(is_registered);
    // Special Irregular volume data importer for VTK files
    is_registered = get_interface()->register_serializable_class<nv::index_common::Irregular_volume_importer_vtk>();
    assert(is_registered);    
    // ------------------------------------------------------------------------------------------------------------

#ifdef NV_IDX_USE_OPENVDB_IMPORTER
    // ------------------------------------------------------------------------------------------------------------
    // Sparse volume data importers
    is_registered = get_interface()->register_serializable_class<nv::index_common::Sparse_volume_importer>();
    assert(is_registered);
#endif

    // ------------------------------------------------------------------------------------------------------------
    // Register texturing techniques

    // Register distributed compute techniques
    is_registered = dice_configuration->register_serializable_class<nv::index_common::Distributed_compute_checkerboard_2d>();
    assert(is_registered);
    is_registered = dice_configuration->register_serializable_class<nv::index_common::Distributed_compute_checkerboard_3d>();
    assert(is_registered);
    is_registered = dice_configuration->register_serializable_class<nv::index_common::Distributed_compute_mandelbrot_2d>();
    assert(is_registered);
    is_registered = dice_configuration->register_serializable_class<nv::index_common::Distributed_compute_bitmap_mapping_2d>();
    assert(is_registered);
    is_registered = dice_configuration->register_serializable_class<nv::index_common::Distributed_compute_height_color_mapping_2d>();
    assert(is_registered);

    // Register heightfield access/edit job
    is_registered = dice_configuration->register_serializable_class<nv::index_common::Distributed_heightfield_elevation_change>();
    assert(is_registered);
    is_registered = dice_configuration->register_serializable_class<Heightfield_data_retrieval_appjob>();
    assert(is_registered);

    // Register volume edit job
    is_registered = dice_configuration->register_serializable_class<nv::index_common::Distributed_voxel_value_change>();
    assert(is_registered);

    is_registered = dice_configuration->register_serializable_class<Distributed_volume_bordermap_generation>();
    assert(is_registered);

    is_registered = dice_configuration->register_serializable_class<Distributed_volume_filter>();
    assert(is_registered);

    is_registered = dice_configuration->register_serializable_class<Volume_bordermap_database_element>();
    assert(is_registered);

    // Register the simple gridder algorithm
    is_registered = dice_configuration->register_serializable_class<Simple_gridding_operation>();
    assert(is_registered);
    is_registered = dice_configuration->register_serializable_class<Gridder_elevation_mean_compute>();
    assert(is_registered);

    // Register classes implementing the distributed autotracking algorithm
    is_registered = dice_configuration->register_serializable_class<Autotracking_workflow_operation>();
    assert(is_registered);
    is_registered = dice_configuration->register_serializable_class<Distributed_heightfield_elevation_delete>();
    assert(is_registered);
    is_registered = dice_configuration->register_serializable_class<Distributed_seed_point_evaluation>();
    assert(is_registered);
    is_registered = dice_configuration->register_serializable_class<Distributed_heightfield_seedings_updates>();
    assert(is_registered);
    is_registered = dice_configuration->register_serializable_class<Flow_grid>();
    assert(is_registered);

    // Register draping job
    is_registered = dice_configuration->register_serializable_class<Distributed_amplitude_map_job>();
    assert(is_registered);

    // Register draping job
    is_registered = dice_configuration->register_serializable_class<Irregular_volume_data_editing>();
    assert(is_registered);

    // Register helpers
    is_registered = get_interface()->register_serializable_class<Distributed_command_execution>();
    assert(is_registered);
        
#ifdef NVINDEX_HAS_MPI_IPC_COMPUTE    
    is_registered = get_interface()->register_serializable_class<CUDA_IPC_volume_editing>();
    assert(is_registered);
#endif // NVINDEX_HAS_MPI_IPC_COMPUTE   

#ifdef NVINDEX_HAS_COMPUTE
    Compute_wrapper::register_serializable_classes(get_interface().get());
#endif


    // We load the default video codec plugin which will be used to encode the rendered frames
    mi::base::Handle<mi::neuraylib::IPlugin_configuration> plugin_configuration(
        get_interface()->get_api_component<mi::neuraylib::IPlugin_configuration>());
    assert(plugin_configuration.is_valid_interface());

    // check video streaming: rtmp, bridge, html5
    assert(app_project.is_defined("dice::rtmp_video_streaming::enabled"));    
    assert(app_project.is_defined("dice::bridge_video_streaming::enabled"));    
    assert(app_project.is_defined("dice::html5_video_streaming::enabled"));    
    const bool is_rtmp_enabled = 
        nv::index_common::get_bool(app_project.get("dice::rtmp_video_streaming::enabled"));
    const bool is_html5_enabled = 
        nv::index_common::get_bool(app_project.get("dice::html5_video_streaming::enabled"));
    const bool is_bridge_enabled = 
        nv::index_common::get_bool(app_project.get("dice::bridge_video_streaming::enabled"));
    const bool is_one_of_video_stream_enabled = 
        is_rtmp_enabled || is_bridge_enabled || is_html5_enabled;

    if (is_one_of_video_stream_enabled)
    {
#ifndef MI_PLATFORM_WINDOWS
        // for all video streams
        if (plugin_configuration->load_plugin_library("screen_video.so") != 0)
        {
            ERROR_LOG << "Failed to load video plugin 'screen video'.";
        }
        if (plugin_configuration->load_plugin_library("x264_video.so") != 0)
        {
            WARN_LOG << "Optional video plugin 'x264_video' could not be loaded.";
        }
        if (plugin_configuration->load_plugin_library("nvenc_video.so") != 0)
        {
            WARN_LOG << "Optional video plugin 'nvenc_video' could not be loaded.";
        }

        if (is_html5_enabled)
        {
            if (plugin_configuration->load_plugin_library("mp4_encoder.so") != 0)
            {
                ERROR_LOG << "'mp4_encoder' plugin could not be loaded. "
                          << "html5 video streaming feature needs this plugin. "
                          << "Please check the plugin path configuration.";
                return false;
            }
        }
#else  // MI_PLATFORM_WINDOWS
        // for all video streams
        if (plugin_configuration->load_plugin_library("screen_video.dll") != 0)
        {
            ERROR_LOG << "Failed to load video plugin 'screen video'.";
        }
        if (plugin_configuration->load_plugin_library("x264_video.dll") != 0)
        {
            ERROR_LOG << "Failed to load video plugin 'x264_video'.";
        }

        if (is_html5_enabled)
        {
            if (plugin_configuration->load_plugin_library("mp4_encoder.dll") != 0)
            {
                ERROR_LOG << "'mp4_encoder' plugin could not be loaded. "
                          << "html5 video streaming feature needs this plugin. "
                          << "Please check the plugin path configuration.";
                return false;
            }
        }
#endif  // MI_PLATFORM_WINDOWS
    }

    // Define service mode of the cluster machine
    mi::base::Handle<nv::index::ICluster_configuration> cluster_configuration(
        get_interface()->get_api_component<nv::index::ICluster_configuration>());
    irc.m_cluster_change_callback = new nv::index_common::Cluster_change_callback();
        cluster_configuration->register_callback(irc.m_cluster_change_callback.get());
    const std::string service_mode = app_project.get("index::service", "rendering_and_compositing");
    if(cluster_configuration->set_service_mode(service_mode.c_str()))
    {
        if (enable_networking)
        {
            INFO_LOG << "Setting IndeX service mode to [" << service_mode << "].";
        }
    }

    const mi::Sint32 sub_cluster_id = nv::index_common::get_sint32(app_project.get("index::sub_cluster_id", "-1"));
    if(sub_cluster_id >= 0)
    {
        if(cluster_configuration->set_sub_cluster_id(mi::Uint32(sub_cluster_id)))
        {
            if (enable_networking)
            {
                INFO_LOG << "Using custom sub-clustering, setting sub_cluster_id to [" << sub_cluster_id << "].";
            }
        }
        else
        {
            ERROR_LOG << "Error setting sub_cluster_id to [" << sub_cluster_id << "].";
        }
        // igeneral_configuration->set_host_property("sub_cluster_id", sub_cluster_id.c_str());
    }
    else
    {
        const mi::Uint32 min_sub_cluster_size
            = nv::index_common::get_uint32(app_project.get("dice::min_sub_cluster_size", "1"));
        const mi::Uint32 max_nr_of_sub_clusters
            = nv::index_common::get_uint32(app_project.get("dice::max_nr_of_sub_clusters", "1"));

        if(!cluster_configuration->set_automatic_subclustering(min_sub_cluster_size, max_nr_of_sub_clusters))
        {
            ERROR_LOG <<"Error setting DiCE sub-cluster parameters min_sub_cluster_size[" << min_sub_cluster_size 
                      << "] and max_nr_of_sub_clusters[" << max_nr_of_sub_clusters << "]";
        }           
    }

    std::string const host_name = nv::index_common::get_host_name();
    INFO_LOG << "Starting up NVIDIA IndeX service on host '" << host_name <<  "'...";
    // Finally start the IndeX service. The license is checked now.
    if (m_nvindex_interface->start(true) != 0)
    {
        return false;
    }
    
    // Check if hostname reported by DiCE matches the one return by get_host_name()
    {
        mi::base::Handle<const mi::neuraylib::IHost_properties> hprop(
            igeneral_configuration->get_host_properties());
        assert(hprop.is_valid_interface());
        mi::base::Handle<const mi::IString> dice_host_name(hprop->get_property("host_name"));
        if (dice_host_name.is_valid_interface())
        {
            if (host_name != dice_host_name->get_c_str())
            {
                WARN_LOG << "Hostname '" << dice_host_name->get_c_str()
                         <<  "' reported by IHost_properties differs from gethostname() '"
                         << host_name << "'";
            }
        }
        else
        {
            // No network (cluster) mode
        }
    }

    return true;
}

//----------------------------------------------------------------------
nv::index::IIndex* Nvindex_library_accessor::access_index_interface(
    const std::string& nvindex_library_variant)
{
    typedef nv::index::IIndex*(Factory());

#ifdef MI_PLATFORM_WINDOWS
    // Loading the IndeX windows dynamics link library
    const std::string library_name = nvindex_library_variant + ".dll";
    void* handle = LoadLibrary(TEXT(library_name.c_str()));

    if (!handle)
    {
        ERROR_LOG << "Could not retrieve a handle to the '"<< library_name << "' library.";
        return NULL;
    }
    void* entry_point = GetProcAddress((HMODULE)handle, "nv_index_factory");
#else
    const std::string library_name = nvindex_library_variant + ".so";
    // Loading the DiCE linux shared library
    // RTLD_LAZY causes crash on exit. see Bugzilla bug
    // 11379. However, RTLD_NOW doesn't solve the problem.
    // However, Environment variable LD_BIND_NOW=1 works.
    void* handle = dlopen(library_name.c_str(), RTLD_LAZY);
    if (!handle)
    {
        ERROR_LOG << "Could not retrieve a handle to the '" << library_name << "' library: " << dlerror();
        return NULL;
    }
    void* entry_point = dlsym(handle, "nv_index_factory");
#endif // MI_PLATFORM_WINDOWS

    if (!entry_point)
    {
        ERROR_LOG << "Could not retrieve the entry point into the '" << library_name << "' library.";
        return NULL;
    }

    Factory* factory = (Factory*)entry_point;
    return factory();
}

//----------------------------------------------------------------------
void Nvindex_library_accessor::clear()
{
    if (m_nvindex_interface.is_valid_interface())
    {
        mi::base::Handle<mi::neuraylib::INeuray> dice(
            m_nvindex_interface->get_dice_interface());
        if (dice.is_valid_interface())
        {
            // Removing the logger
            mi::base::Handle<mi::neuraylib::ILogging_configuration> logging_configuration(
                get_interface()->get_api_component<mi::neuraylib::ILogging_configuration>());
            assert(logging_configuration.is_valid_interface());
            logging_configuration->set_receiving_logger(NULL);
            logging_configuration = NULL;

            // shutdown forwarding logger factory. Also need
            // delete_instance which will be later
            nv::index_common::Forwarding_logger_factory::instance()->shutdown();
        }
    }
}

//----------------------------------------------------------------------
