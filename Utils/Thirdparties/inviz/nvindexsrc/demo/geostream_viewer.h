/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief NVIDIA IndeX Technical Reference Viewer

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_GEOSTREAM_VIEWER_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_GEOSTREAM_VIEWER_H

#include "common/string_dict.h"

#ifdef NVINDEX_HAS_COMPUTE
#include "compute/compute_wrapper.h"
#endif

#include "nvindex_rendering_context.h"

#ifdef NVINDEX_HAS_MPI_IPC_COMPUTE
#include "cuda_ipc_volume_editing.h"
#endif // NVINDEX_HAS_MPI_IPC_COMPUTE


class Nvindex_AppData;
class Opengl_application_buffer;

/// Index reference viewer application state class
class Index_application
{
public:
    /// constructor
    Index_application();
    /// destructor
    ~Index_application();

    //
    // Setup
    //
    
    /// setup the IndeX viewer. This also setup the application rendering context.
    bool setup_viewer(
        Nvindex_rendering_context& irc,
        bool                       is_remote_viewer,
        bool                       update_session = true);

    /// setup the IndeX remote process.
    bool setup_remote(Nvindex_rendering_context& irc);

    //
    // Render calls
    //

    /// viewer render call
    void render_frame(
        Nvindex_rendering_context& irc,
        mi::Uint32                 frame_num,
        bool                       need_compute = false);

    /// remote render call
    void render_frame_remote(Nvindex_rendering_context& irc);

    //
    // Rendering loops
    //
    
    /// remote rendering loop
    bool run_remote(Nvindex_rendering_context& irc);

    //
    // Tools
    //

    /// Process the command line arguments to import project and config files
    ///
    /// \param[in] argc main argc
    /// \param[in] argv main argv
    /// \return true when command line option parsing succeeded
    bool process_command_line_arguments(
        mi::Sint32                     argc,
        char**                         argv,
        nv::index_common::String_dict& result_opt);

    /// save span buffer snapshot
    /// 
    /// \param[in] span_buffer the span buffer to be saved
    /// \param[in] frame_num   current frame number
    void save_span_buffer_snapshot(Span_renderer_IF* span_buffer,
                                   mi::Sint32 frame_num);

    /// Shutdown all Servers (http, rtmp, bridge, websocket)
    ///
    /// \param[in] irc index rendering context
    void shutdown_all_server(Nvindex_rendering_context& irc);

    /// shutdown the application
    ///   de-allocate application state
    void shutdown();
    
#ifdef NVINDEX_HAS_MPI_IPC_COMPUTE

    //----------------------------------------------------------------------
    /// Get the local host id
    /// \return     the host id of the local machine running IndeX
    ///
    mi::Uint32 get_local_host_id(Nvindex_rendering_context& irc) const;

    //----------------------------------------------------------------------
    /// set volume bricks affinity information
    /// Used by IndeX to distribute the bricks over the cluster.
    ///
    /// \param[in] bricks_affinity  The list brick affinities
    ///
    void set_bricks_affinity(
        Nvindex_rendering_context& irc,
        const std::vector<Brick_affinity> &bricks_affinity);
        
    //----------------------------------------------------------------------
    /// set volume bricks locality information
    /// Used during rendering to locate and retrieve volume bricks data.
    ///
    /// \param[in] bricks_locality  The list brick localities
    ///
    void set_bricks_locality(const std::vector<Brick_locality> &bricks_locality);
    
#endif // NVINDEX_HAS_MPI_IPC_COMPUTE

private:
    /// set up application, for both viewer and remote common
    bool setup_common(Nvindex_rendering_context& irc);

    /// Set up Servers: start rtmp, html, bridge, websockt servers
    bool setup_all_server(Nvindex_rendering_context& irc);

    /// Get http address from the project
    std::string get_http_address() const;

    /// Set up RTMP and HTTP server    
    bool setup_rtmp_http(Nvindex_rendering_context& irc);

    /// Set up DiCE bridge
    bool setup_dice_bridge(Nvindex_rendering_context& irc);

    /// Set up html5 video stream
    bool setup_html5_video_stream(Nvindex_rendering_context& irc);

    /// set up the camera
    void setup_camera(
        Nvindex_rendering_context&        irc,
        mi::neuraylib::IDice_transaction* dice_transaction);

    /// set up application configuration
    ///  application status setup via project file    
    void setup_application_configuration(
        nv::index_common::String_dict*   prj);

    /// set up IndeX configuration.
    ///   IndeX configuration setup via API calls
    void setup_index_configuration(
        Nvindex_rendering_context&        irc,
        nv::index_common::String_dict*    prj,
        mi::neuraylib::Tag                config_tag,
        mi::neuraylib::IDice_transaction* dice_transaction);
    
    ///----------------------------------------------------------------------

    /// render frame camera update
    ///
    /// \param[in] irc IndeX rendering context
    void render_frame_update_camera(
        Nvindex_rendering_context& irc);

    /// render frame stereo view (experimental)
    ///
    /// \param[in] irc IndeX rendering context
    /// \return current camera tag
    void render_frame_stereo_view(
        mi::neuraylib::IDice_transaction* dice_transaction);

    /// render frame check host state. Does any host leave the cluster?
    /// Shall we quite the rendering call?
    /// 
    /// \param[in] dice_transaction dice transaction
    /// \return when false, we should quit the rendering
    bool render_frame_check_cluster_state_is_quit(
        Nvindex_rendering_context& irc);

    /// render frame compute related
    /// 
    /// \param[in] dice_transaction dice transaction
    void render_frame_compute(
        Nvindex_rendering_context&        irc,
        mi::Uint32                        frame_num,
        bool                              need_compute,
        mi::Float32&                      application_compute_time,
        mi::Uint64&                       application_compute_memory,
        mi::neuraylib::IDice_transaction* dice_transaction);

    /// render frame animation handling
    ///
    /// \param[in] dice_transaction dice transaction
    void render_frame_handle_animation(
        Nvindex_rendering_context&                          irc,
        mi::Uint32                                          frame_num,
        mi::base::Handle<mi::neuraylib::IDice_transaction>& dice_transaction);

    /// render frame export session handling
    /// 
    /// \param[in] dice_transaction dice transaction
    void render_frame_handle_export_session(
        Nvindex_rendering_context&        irc,
        mi::neuraylib::IDice_transaction* dice_transaction);

    /// render frame check the quit option: case 0
    /// 
    /// \return when false, we should quit the rendering
    bool render_frame_quit_handle_case_0();

    /// render frame check the quit option: case n
    /// 
    /// \return when false, we should quit the rendering
    bool render_frame_quit_handle_case_n(mi::Uint32 frame_num);

private:
    // render frame common: before database access
    void render_frame_common_setup_before_db_access(    
        Nvindex_rendering_context& irc,        
        mi::Uint32                 frame_num);

    /// render frame common: with database access, before render call
    /// \return false needs to quit rendering
    bool render_frame_common_setup_with_db_access(    
        Nvindex_rendering_context&        irc,        
        mi::Uint32                        frame_num,
        bool                              need_compute,
        mi::Float32&                      application_compute_time,
        mi::Uint64&                       application_compute_memory,
        mi::neuraylib::IDice_transaction* dice_transaction);

    /// render frame common: post process, after render call
    /// \return false needs to quit rendering
    bool render_frame_common_postprocess(    
        Nvindex_rendering_context&                          irc,
        mi::Uint32                                          frame_num,
        mi::base::Handle<mi::neuraylib::IDice_transaction>& dice_transaction);

    /// index canvas setup before the rendering
    void render_frame_common_canvas_setup(
        Nvindex_rendering_context& irc,
        bool                       is_snapshot_on_this_frame);

    /// OpenGL canvas setup if OpenGL is enabled.
    void render_frame_common_opengl_canvas_setup(
        Nvindex_rendering_context& irc);

    /// restore canvas size if snapshot.
    void render_frame_common_restore_canvas_size(
        bool is_snapshot_on_this_frame);

    /// check OpenGL error
    void render_frame_common_check_gl_error();

private:
    // single view rendering
    
    /// rendering, the frame's performance information update, and
    /// statistics.  For single view.
    void render_frame_render_performance_statistics_single_view(
        Nvindex_rendering_context&        irc,
        mi::Uint32                        frame_num,
        mi::Float32                       application_compute_time,
        mi::Uint64                        application_compute_memory,
        mi::neuraylib::IDice_transaction* dice_transaction);

    /// render frame timestep handling (experimental)
    /// 
    void render_frame_timestep_handling(
        Nvindex_rendering_context&        irc,
        mi::neuraylib::IDice_transaction* dice_transaction);

    /// subroutine for render_frame: single view entry
    void render_frame_single_view(
        Nvindex_rendering_context& irc,
        mi::Uint32                 frame_num,
        bool                       need_compute);

    /// subroutine for render_frame: single view OpenGL integration setup
    void render_frame_single_view_opengl_integration_setup(
        Nvindex_rendering_context&        irc,
        mi::neuraylib::IDice_transaction* dice_transaction);

    /// subroutine for render_frame: single view render
    mi::base::Handle<nv::index::IFrame_results> render_frame_single_view_render_call(
        Nvindex_rendering_context&        irc,
        bool                              is_snapshot_on_this_frame,
        mi::neuraylib::IDice_transaction* dice_transaction);

    /// Additional opengl rendering (debug info, polygon delete visualization)
    void render_frame_single_view_additional_gl_rendering(
        Nvindex_rendering_context&        irc,
        mi::neuraylib::IDice_transaction* dice_transaction);
    
    /// Render slice wireframe
    /// TODO: support multiple volumes
    void render_frame_wireframe_slice(
        mi::neuraylib::IDice_transaction* dice_transaction);

    /// get OpenGL app buffer when enabled.
    Opengl_application_buffer* render_frame_get_opengl_app_buffer() const;

    /// render frame: show frame info
    void render_frame_show_frame_info(
        Frame_info_callbacks* frame_info_callbacks) const;

    /// visualize performance etc. with 2d OpenGL
    void render_frame_single_view_visualize_perf_with_2d_opengl(
        Nvindex_rendering_context&        irc,
        nv::index::IFrame_results*        frame_results,
        mi::neuraylib::IDice_transaction* dice_transaction);

private:
    /// multi-view mode
    void render_frame_multi_view(
        Nvindex_rendering_context& irc,
        mi::Uint32                 frame_num,
        bool                       need_compute);    

    /// multi-view: render call, performance, statistics
    void render_frame_render_performance_statistics_multi_view(
        Nvindex_rendering_context& irc,
        mi::Uint32                 frame_num,
        mi::Float32                application_compute_time,
        mi::Uint64                 application_compute_memory);

    /// multi-view render call
    /// \return true when rendering succeeded
    bool render_frame_multi_view_render_call(
        Nvindex_rendering_context&        irc,
        bool                              is_snapshot_on_this_frame);

    /// visualize performance etc. with 2d OpenGL
    void render_frame_multi_view_visualize_perf_with_2d_opengl(
        Nvindex_rendering_context&        irc);

private:
#ifdef NVINDEX_HAS_MPI_IPC_COMPUTE
    /// edit compute volume
    void edit_volume_data(
        mi::neuraylib::Tag                                           volume_tag,
        mi::neuraylib::Tag                                           session_tag,
        mi::neuraylib::IDice_transaction*                            dice_transaction);

#endif // NVINDEX_HAS_MPI_IPC_COMPUTE       

    // DEMO HACK/FEATURE/SOLUTION: GTC 2017
    void unstructured_volume_compute(
        mi::Uint32                                          frame_id,
        const std::string&                                  db_dataset_name,
        mi::neuraylib::Tag                                  session_tag,
        mi::neuraylib::IDice_transaction*                   dice_transaction);
    
    /// Setup for multi-view mode.
    ///
    /// \param[in] irc  index demo viewer rendering context
    /// \param[in] dice_transaction dice transaction
    void setup_multi_view_mode(
        Nvindex_rendering_context& irc,
        mi::neuraylib::IDice_transaction* dice_transaction);

private:
    Nvindex_AppData*               m_appdata;
    nv::index_common::String_dict* m_prj;
    /// number of hosts of the last frame
    mi::Uint32                     m_previous_nb_hosts;

#ifdef NVINDEX_HAS_COMPUTE
    Compute_wrapper                m_compute_wrapper;
#endif

#ifdef NVINDEX_HAS_MPI_IPC_COMPUTE
    std::vector<Brick_locality>    m_bricks_locality;       // bricks locality information.
#endif // NVINDEX_HAS_MPI_IPC_COMPUTE
};

//----------------------------------------------------------------------
#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_GEOSTREAM_VIEWER_H
