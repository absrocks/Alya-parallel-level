/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/

#ifdef WIN32
# include <windows.h>
#undef min // remove the ugly defines in windows.h
#undef max
#endif  // WIN32

#include "geostream_viewer.h"
#include "bridge_server.h"

#include <cassert>
#include <iomanip>
#include <iterator>
#include <sstream>

#include <nv/index/icamera.h>
#include <nv/index/icolormap.h>
#include <nv/index/iconfig_settings.h>
#include <nv/index/iregular_heightfield.h>
#include <nv/index/iregular_volume.h>
#include <nv/index/iirregular_volume_scene_element.h>
#include <nv/index/iscene.h>
#include <nv/index/isession.h>
#include <nv/index/islice_scene_element.h>

#include "common/affinity_information.h"
#include "common/clock_pulse_generator.h"
#include "common/canvas.h"
#include "common/forwarding_logger.h"
#include "common/ppm_io.h"
#include "common/string_dict.h"

#include "camera_utility.h"
#include "colormap_manager.h"
#include "colormap_util.h"
#include "compute_utility.h"
#include "examiner_manipulator.h"
#include "geostream_io.h"
#include "heightfield_workflow_functionality.h"
#include "html5_video_stream.h"
#include "irregular_volume_data_editing.h"
#include "nvindex_appdata.h"
#include "nvindex_library_accessor.h"
#include "opengl_appdata.h"
#include "performance_log.h"
#include "progress_callback.h"
#include "rtc_parameter_buffer_manip.h"
#include "rtmp_handler.h"
#include "scene_utility.h"
#include "span_renderer_no_gl.h"
#include "utilities.h"
#include "websocket_utility.h"

// OpenGL related includes: the macro USE_OPENGL is defined in Makefile
// But the following headers don't depend on OpenGL
#include "opengl_drawing_utilities.h"
#include "opengl_application_draw.h"
#include "opengl_application_buffer.h"

#ifdef USE_OPENGL
#  include <GL/glew.h>
#  include "span_renderer_gl.h"
#  include "opengl_offscreen_context.h"

#ifndef WIN_NT
#include <X11/Xlib.h>
#include "xwin_scoped_error_handler.h"
#endif  // WIN_NT

#endif // USE_OPENGL

#include "stereo_camera.h"      // Stereo ICamera test
#include "multiple_camera.h"    // Multiple stereo camera test
#include "camera_animator.h"    // Performance test with camera movement

#include "rtmp_call_event_handler.h"
#include "rtmp_connect_event_handler.h"

using namespace nv::index_common;

namespace {

String_dict get_default_command_line_options()
{
    String_dict default_opt;

    default_opt.insert("app::project_file",                         "");

    // dice network
    default_opt.insert("dice::http_listen",                         "0.0.0.0");
    default_opt.insert("dice::http_port",                           "8080");
    default_opt.insert("dice::admin_http_listen",                   "0.0.0.0");
    default_opt.insert("dice::admin_http_port",                     ""); // disabled by default
    default_opt.insert("dice::network::multicast_address",          "");
    default_opt.insert("dice::network::mode",                       "");
    default_opt.insert("dice::min_sub_cluster_size",                "1");
    default_opt.insert("dice::max_nr_of_sub_clusters",              "1");

    // dice rtmp video streaming
    default_opt.insert("dice::rtmp_video_streaming::enabled",             "yes");
    default_opt.insert("dice::rtmp_video_streaming::listen",              "0.0.0.0");
    default_opt.insert("dice::rtmp_video_streaming::port",                "1935");
    default_opt.insert("dice::rtmp_video_streaming::video_codec",         "screen video");
    default_opt.insert("dice::rtmp_video_streaming::video_bitrate",       "30000000");
    default_opt.insert("dice::rtmp_video_streaming::video_bitrate_error", "0.5");
    default_opt.insert("dice::rtmp_video_streaming::video_framerate",     "25");
    default_opt.insert("dice::rtmp_video_streaming::flash_client",        "");

    // dice bridge video stream
    default_opt.insert("dice::bridge_video_streaming::enabled",     "false");

    // dice html5 video streaming
    default_opt.insert("dice::html5_video_streaming::enabled",      "no");    
    default_opt.insert("dice::html5_video_streaming::port",         "3333");    
    default_opt.insert("dice::html5_video_streaming::mp4_video_encoder", "h264");

    // index
    default_opt.insert("index::service",                            "rendering_and_compositing");
    default_opt.insert("app::host_mode",                            "viewer");
    default_opt.insert("index::horizontal_spans",                   "1");
    default_opt.insert("index::result_queue",                       "1");
    default_opt.insert("index::enable_automatic_span_control",      "off");
    default_opt.insert("index::max_spans_per_machine",              "1");
    default_opt.insert("index::sub_cluster_id",                     "-1");       // default is invalid
    default_opt.insert("index::cuda_debug_checks",                  "off");      // disabled by default
    default_opt.insert("index::timesteps",                          "0");
    default_opt.insert("index::dynamic_memory_management",          "on");
    default_opt.insert("index::workload_balancing",                 "off");
    
    // demo viewer canvas
    default_opt.insert("app::canvas::background_color",             "0 0 0 1");
    default_opt.insert("app::canvas::background_color_1",           "1 1 1 1");
    default_opt.insert("index::image_file_canvas_resolution",       "1024 1024"); // default image size

    // demo application performance monitoring related default
    default_opt.insert("index::performance::set_monitor_performance_values", "no");
    default_opt.insert("app::performance::logger_type",          "csv");
    default_opt.insert("app::performance::global::enable_log",   "no");
    default_opt.insert("app::performance::global::file",         "/dev/null");
    default_opt.insert("app::performance::global::item_list",    "<all>");
    default_opt.insert("app::performance::per_host::enable_log", "no");
    default_opt.insert("app::performance::per_host::file",       "/dev/null");
    default_opt.insert("app::performance::per_host::item_list",  "<all>");
    default_opt.insert("app::performance::per_span::enable_log", "no");
    default_opt.insert("app::performance::per_span::file",       "/dev/null");
    default_opt.insert("app::performance::per_span::item_list",  "<all>");

    // Buffer resolution
    default_opt.insert("index::canvas_resolution",          "1024 1024"); // default screen size

    // examiner settings
    default_opt.insert("app::examiner::initial_rotation_center::type",  "roi_center");
    default_opt.insert("app::examiner::initial_rotation_center::world_coordinates",  "0 0 0");

    return default_opt;
}

//----------------------------------------------------------------------
/// visualize region of interest and volume extents for debugging purpose
///
/// \param[in] is_render_roi when true, show the region of interest
/// using OpenGL functions
/// \param[in] is_render_extents when true, show the extents using
/// OpenGL functions
/// \param[in] session_tag       session tag
/// \param[in] roi               the region of interest bbox
/// \param[in] dice_transaction  dice transaction
void visualize_region_of_interest(
    bool                                  is_render_roi,
    bool                                  is_render_extents,
    mi::neuraylib::Tag                    session_tag,
    const mi::math::Bbox<mi::Float32, 3>& roi,
    mi::neuraylib::IDice_transaction*     dice_transaction)
{
    if((!is_render_roi) && (!is_render_extents))
    {
        return;
    }
    assert(dice_transaction != 0);
    assert(session_tag.is_valid());
    gl_visualize_roi(is_render_roi, roi, mi::math::Color(0.6f, 0.6f, 1.f, 1.f));
            
    const mi::math::Bbox<mi::Float32, 3> scene_bbox_f32 = 
        get_scene_bounding_box(session_tag, dice_transaction);
    gl_visualize_bbox(
        is_render_extents, scene_bbox_f32, mi::math::Color(0.0f, 0.0f, 1.0f, 1.0f));

    const mi::math::Color base_color(0.8f, 0.5f, 0.5f, 0.8f);

    std::vector<mi::neuraylib::Tag> volume_tag_vec = Nvindex_AppData::instance()->get_volume_tag_vec();
    for (std::vector<mi::neuraylib::Tag>::const_iterator vi = volume_tag_vec.begin();
        vi != volume_tag_vec.end(); ++vi)
    {
        assert(vi->is_valid());
        mi::base::Handle<const nv::index::IRegular_volume> volume(
            dice_transaction->access<const nv::index::IRegular_volume>(*vi));

        gl_visualize_bbox(is_render_extents, volume->get_XYZ_bounding_box(), base_color);

        const mi::math::Matrix<mi::Float32, 4, 4> mat = volume->get_transform();
        mi::math::Bbox<mi::Float32, 3> ijk_roi = volume->get_IJK_region_of_interest();
        // Border voxels are half the size of interiour ones, so adapt for rendering wireframe box.
        ijk_roi.max -= mi::math::Vector<mi::Float32, 3>(1.0f);

        const mi::math::Color color(0.7f, 0.7f, 0.7f, 0.9f);
        gl_visualize_roi(is_render_roi, ijk_roi, color, mat);
    }

    std::vector<mi::neuraylib::Tag> height_tag_vec = Nvindex_AppData::instance()->get_heightfield_tag_vec();
    for (std::vector<mi::neuraylib::Tag>::const_iterator hi = height_tag_vec.begin();
         hi != height_tag_vec.end(); ++hi)
    {
        assert(hi->is_valid());
        mi::base::Handle<const nv::index::IRegular_heightfield> heightfield_scene_element(
            dice_transaction->access<const nv::index::IRegular_heightfield>(*hi));
        mi::math::Color color(1.0f);
        gl_visualize_bbox(is_render_extents, heightfield_scene_element->get_XYZ_clipped_bounding_box(), color);

        const mi::math::Matrix<mi::Float32, 4, 4> mat = heightfield_scene_element->get_transform();

        mi::math::Bbox<mi::Float32, 3> ijk_roi = heightfield_scene_element->get_IJK_region_of_interest();
        ijk_roi.max -= mi::math::Vector<mi::Float32, 3>(1.0f);
        gl_visualize_roi(is_render_roi, ijk_roi, color, mat);

        mi::math::Bbox<mi::Float32, 3> ijk_bbox = heightfield_scene_element->get_IJK_region_of_interest();
        ijk_bbox.max -= mi::math::Vector<mi::Float32, 3>(1.0f);
        gl_visualize_bbox(is_render_extents, ijk_bbox, color, mat);
    }
}

//----------------------------------------------------------------------
/// visualization:
///  - horizontal spans
///  - statistics by OpenGL
///  - compositing workload
///  - performance overlay
///  - color table
///
/// \param[in] irc IndeX application rendering context
/// \param[in] appdata IndeX application data
/// \param[in] colormap_tag colormap tag to visualize
/// \param[in] performance_values performance values to visualize
/// \param[in] dice_transaction db transaction
void visualize_span_and_performance(
    Nvindex_rendering_context& irc,
    Nvindex_AppData *          appdata,
    const mi::neuraylib::Tag & colormap_tag,
    nv::index::IPerformance_values * performance_values,
    mi::neuraylib::IDice_transaction * dice_transaction
    )
{
    assert(appdata != 0);
    assert(performance_values != 0);

    // visualize horizontal span
    gl_visualize_horizontal_span(appdata->is_show_horizontal_spans_overlay(),
                                 irc.m_span_buffer.get(),
                                 performance_values);

    // Draw the per host compositing workload
    std::vector<mi::math::Bbox_struct<mi::Uint32, 2> > screen_space_subdivision;
    irc.m_span_buffer->get_screen_space_subdivision(screen_space_subdivision);

    // Render some statistics into the OpenGL window contents
    if (appdata->m_show_rendering_statistics_overlay)
    {
        draw_rendering_workload(
            irc.m_icluster_configuration->get_number_of_hosts(), 0.5f, performance_values);
    }
    if (appdata->m_show_compositing_statistics_overlay)
    {
        draw_compositing_workload(screen_space_subdivision.size(), 0.5f, performance_values);
    }
    if (appdata->m_show_all_performance_statistics_overlay)
    {
        draw_performance(0.5f, performance_values);
    }
    // Render color table for editing
    gl_render_color_table(Nvindex_AppData::instance()->get_opengl_appdata()->is_show_color_table(),
                          colormap_tag, dice_transaction);
}

void show_statistics_overlay(
    bool                                     is_show_statictics,
    const nv::index::ICluster_configuration* rendering_properties_query,
    const nv::index::IPerformance_values*    performance_values,
    bool                                     show_large)
{
    if (!is_show_statictics)
    {
        return;                 // no need to show statistics
    }

    std::vector<int>   per_host_rendered_cubes;
    std::vector<float> per_host_time_rendering;
    std::vector<float> per_host_time_heightfield;
    std::vector<float> per_host_time_volume;
    std::vector<float> per_host_time_volume_and_horizons;
    std::vector<bool>  per_host_using_cpu;
    std::vector<bool>  per_host_streaming;
    for (mi::Uint32 i=0; i < rendering_properties_query->get_number_of_hosts(); ++i)
    {
        mi::Uint32 host_id = i + 1; //TODO: This assumes a continuous host id list
        per_host_rendered_cubes.push_back(static_cast<mi::Sint32>(performance_values->get("nb_subcubes_rendered", host_id)));
        per_host_time_rendering.  push_back(performance_values->get_time("time_rendering", host_id));
        per_host_time_heightfield.push_back(performance_values->get_time("time_rendering_horizon", host_id));
        per_host_time_volume.    push_back(performance_values->get_time("time_rendering_volume", host_id));
        per_host_time_volume_and_horizons.push_back(performance_values->get_time("time_rendering_volume_and_horizon", host_id));
        per_host_using_cpu.push_back(!performance_values->get("is_using_gpu", host_id));
        per_host_streaming.push_back(performance_values->get("size_volume_data_upload", host_id) > 0.f);
    }

    draw_host_workload(
        show_large,
        per_host_rendered_cubes, per_host_time_rendering, per_host_time_heightfield, per_host_time_volume,
        per_host_time_volume_and_horizons, per_host_using_cpu, per_host_streaming);

    INFO_LOG << "*************************** Performance counter values:";
    mi::Uint32 nb_types = performance_values->get_nb_type_names();
    for (mi::Uint32 i = 0; i < nb_types; ++i)
    {
        const std::string type_name = performance_values->get_type_name(i);

        if (type_name.find("time_") == 0 || type_name == "frames_per_second")
            INFO_LOG << type_name << "\t"
                     << performance_values->get_time(type_name.c_str()) << " ms";
        else if (type_name.find("size_") == 0)
            INFO_LOG << type_name << "\t"
                     << performance_values->get(type_name.c_str()) / (1024.f * 1024.f) << " MB";
        else if (type_name.find("is_") == 0)
            INFO_LOG << type_name << "\t"
                     << (performance_values->get(type_name.c_str()) == 0 ? "no" : "yes");
        else
            INFO_LOG << type_name << "\t"
                     << performance_values->get(type_name.c_str());

    }
}

//----------------------------------------------------------------------
/// Handle playback of recorded user commands from a text file
///
void handle_playback(const std::string& filename, Nvindex_rendering_context& irc)
{
    static mi::Sint32 last_line = 1;
    if (last_line <= 0)
        return;

    static Call_event_handler* handler = 0;
    if (handler == 0)
    {
        handler = new Call_event_handler(
            0,  // The recording shall only done for the standard viewing scenario
            irc);

        INFO_LOG << "Starting playback of file '" << filename << "'";
    }

    // Pause playback until performance measuring is done.
    if (Nvindex_AppData::instance()->get_remaining_performance_measuring_frames() > 0){
        return;
    }

    // Wait a few frames before starting playback, to make sure everything is loaded
    static mi::Sint32 frame = 0;
    frame++;
    if (frame >= 5)
    {
        last_line = handler->play_recording(filename, last_line);

        if (last_line == 0)
        {
            INFO_LOG << "Playback of file '" << filename << "' finished successfully.";
            handler->release();
            handler = 0;
        }
        else if (last_line == -1)
        {
            ERROR_LOG << "Playback of file '" << filename << "' failed.";
            handler->release();
            handler = 0;
        }
        else
        {
            // Bail out to the render loop to take a snapshot or start a performance measurement
        }
    }
}

//----------------------------------------------------------------------
/// Returns whether rendering should continue. Will return false if no video stream is connected and
/// "app::render_always" is set to false.
bool is_rendering_needed()
{
    if (Nvindex_AppData::instance()->m_render_always)
    {
        return true;
    }

    // Stop rendering if no video stream is connected
    mi::Sint32 nb_connections = Play_event_handler::get_nb_connections();
    if (nb_connections == 0)
    {
        if (Nvindex_AppData::instance()->m_last_nb_video_connections < -1)
        {
            // Render some initial frames even if no client has connected yet, so that data gets
            // loaded now and not when a client connects.
            Nvindex_AppData::instance()->m_last_nb_video_connections++;
        }
        else
        {
            // Nobody is watching, so don't render
            if (Nvindex_AppData::instance()->m_last_nb_video_connections == -1)
                INFO_LOG << "Waiting for first video stream connection to start rendering again...";
            else if (Nvindex_AppData::instance()->m_last_nb_video_connections > 0)
                INFO_LOG << "Last video stream was closed, stopping rendering.";

            Nvindex_AppData::instance()->m_last_nb_video_connections = 0;
            sleep(0.1f); // Wait a bit
            return false;
        }
    }
    else if (Nvindex_AppData::instance()->m_last_nb_video_connections <= 0)
    {
        INFO_LOG << "New video stream connected, restarting rendering...";
        Nvindex_AppData::instance()->m_last_nb_video_connections = nb_connections;
    }
    else
    {
        Nvindex_AppData::instance()->m_last_nb_video_connections = nb_connections;
    }

    return true;
}

//----------------------------------------------------------------------
// Example how to find out which database elements were modified since the last frame
void example_how_to_check_modification_of_db_element(
    Nvindex_rendering_context& irc)
{
    // currently disabled this function
    return;

    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        irc.get_current_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());

    {
        static std::string s_last_timestamp = ""; // Stores the time stamp when this was called last

        // Query the changes
        mi::base::Handle<const mi::neuraylib::ITag_set> changed_database_elems(
            dice_transaction->get_changed_elements_since_time_stamp(s_last_timestamp.c_str()));

        if (changed_database_elems)
        {
            for (mi::Size i = 0; i < changed_database_elems->get_length(); ++i)
            {
                mi::neuraylib::Tag tag = changed_database_elems->get_tag(i);

                mi::base::Handle<const mi::neuraylib::IElement> elem(
                    dice_transaction->access<const mi::neuraylib::IElement>(tag));
                if (elem)
                    INFO_LOG << "Changed database element " << i << ": tag " << tag.id
                             << ", class '" << elem->get_class_name() << "'";
            }
        }
        else if (s_last_timestamp.empty())
        {
            // First run, everything changed so the method returned null as expected.
        }
        else
        {
            // The method returned null and this is not the first run, so there probably were too
            // many changes since the given time stamp.
            INFO_LOG << "There were too many database changes to have a detailed log.";
        }
        s_last_timestamp = dice_transaction->get_time_stamp(); // Update time stamp
    }
    dice_transaction->commit();
}

} // anonymous namespace

//----------------------------------------------------------------------
//======================================================================
// Index_application method implementation
//======================================================================
//----------------------------------------------------------------------

Index_application::Index_application() :
    m_previous_nb_hosts(0)
{
    m_appdata = Nvindex_AppData::instance();
    assert(m_appdata != 0);
    m_prj = m_appdata->peek_app_proj();
    assert(m_prj != 0);
}

Index_application::~Index_application()
{
    // empty
}

//----------------------------------------------------------------------
bool Index_application::setup_viewer(
    Nvindex_rendering_context& irc,
    bool                       is_remote_viewer,
    bool                       update_session)
{
    if (!setup_common(irc))
        return false;

    Nvindex_library_accessor* p_nvindex_access = m_appdata->get_nvindex_library_accessor();
    {
        irc.m_database = p_nvindex_access->get_interface()->get_api_component<mi::neuraylib::IDatabase>();
        assert(irc.m_database.is_valid_interface());

        irc.m_global_scope = irc.m_database->get_global_scope();
        assert(irc.m_global_scope.is_valid_interface());

        mi::Sint32 extra_scopes = get_sint32(m_prj->get("app::local_scopes", "0"));
        for (mi::Sint32 i = 0; i < extra_scopes; ++i)
        {
            irc.m_extra_scopes.push_back(
                mi::base::make_handle(irc.m_database->create_scope(0, 0, false)));
        }

        irc.m_viewing_scenario_scopes.resize(1);
        irc.m_viewing_scenario_scopes[0] = irc.m_database->create_scope(
            irc.m_global_scope.get(),   // The parent scope, if NULL then the global scope (setting to global explicitly)
            0,                          // The privacy level of the parent scope +1
            false);                     // Indicates if the scope is temporary then when the host created the scope is
                                        // removed the scope and all data contained in the scope will be removed.
        assert(irc.m_viewing_scenario_scopes[0].is_valid_interface());
    }

    if (get_bool(m_prj->get("app::setup::enabled", "no")))
    {
        // Start up RTMP and HTTP server
        if (!setup_all_server(irc))
            return false;

        INFO_LOG << "Scene is not configured yet, waiting for user to finish scene setup...";

        while (!irc.m_scene_setup_done)
        {
            sleep(0.5f);
        }
        INFO_LOG << "Scene setup complete, continuing";
    }

    if (irc.m_http_request_handler.is_valid_interface())
    {
        // Inform the HTTP handler about the resolution, as it might have been changed from the initial
        // value by the scene setup dialog
        std::vector<mi::math::Vector<mi::Uint32, 2> > resolutions;
        resolutions.push_back(
            m_appdata->get_user_interaction(0)->get_examiner()->get_main_window_resolution_uint32());
        resolutions.push_back(mi::math::Vector<mi::Uint32, 2>(1024, 1024));
        irc.m_http_request_handler->set_resolutions(resolutions);

        // Guest authentication (after scene setup dialog)
        if (!m_prj->get("dice::http_auth_guest_password").empty())
        {
            irc.m_http_request_handler->set_guest_authentication(
                m_prj->get("dice::http_auth_user"),             // same user
                m_prj->get("dice::http_auth_guest_password"));  // different password
        }
    }

    // create a user interaction (the user_interaction has an examiner)
    {
        m_appdata->create_user_interaction();
        assert(m_appdata->get_nb_user_interaction() == 1);
        m_appdata->get_user_interaction(0)->initialize_by_dict(*m_prj);
    }

    {
        mi::base::Lock::Block block(&m_appdata->m_scene_edit_lock);

        // Begin of internal application scope ...

        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            irc.get_current_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
        assert(dice_transaction.is_valid_interface());

        VERBOSE_LOG << "Application build environment: " << get_app_build_info();
        VERBOSE_LOG << "Creating offscreen context";
        gl_initialize_offscreen_context();

        // ---------------------------------------------------------------------------------------
        // Access the IndeX configurations to create a session
        VERBOSE_LOG << "Retrieving IIndex_session";
        irc.m_iindex_session = irc.m_iindex_if->get_api_component<nv::index::IIndex_session>();
        assert(irc.m_iindex_session.is_valid_interface());

        // Access the IndeX rendering interface for rendering
        VERBOSE_LOG << "Creating an instance of the rendering interface IIndex_rendering";
        irc.m_iindex_rendering = irc.m_iindex_if->create_rendering_interface();
        assert(irc.m_iindex_rendering.is_valid_interface());

        // Access the IndeX rendering query interface for querying performance values and pick results
        VERBOSE_LOG << "Retrieving ICluster_configuration";
        irc.m_icluster_configuration = irc.m_iindex_if->get_api_component<nv::index::ICluster_configuration>();
        assert(irc.m_icluster_configuration.is_valid_interface());

        const bool is_network_mode_on = (m_prj->get("dice::network::mode", "OFF") != "OFF");
        if (is_network_mode_on)
        {
            const mi::Uint32 cluster_size = get_uint32(m_prj->get("app::cluster_size", "0"));
            if (cluster_size > 0)
            {
                mi::Uint32 old_nb_hosts = 0;
                mi::Uint32 nb_hosts = 0;
                INFO_LOG << "Waiting until cluster size reaches " << cluster_size << " hosts";
                while (nb_hosts < cluster_size)
                {
                    nb_hosts = irc.m_icluster_configuration->get_number_of_hosts();
                    if (nb_hosts > old_nb_hosts)
                    {
                        INFO_LOG << "Cluster now has " << nb_hosts << " host" << (nb_hosts == 1 ? "" : "s")
                                 << ", need " << (cluster_size - std::min(nb_hosts, cluster_size))
                                 << " more to continue";
                        old_nb_hosts = nb_hosts;
                    }
                    else if (nb_hosts < old_nb_hosts)
                    {
                        INFO_LOG << "Aborting because at least one host has left - had "
                                 << old_nb_hosts << " hosts, only " << nb_hosts << " left";
                        return false;
                    }
                    sleep(0.5f);
                }
                INFO_LOG << "Cluster size " << cluster_size << " reached, continuing...";
            }
            else
            {
                // join-latency workaround!
                mi::Uint32 old_nb = irc.m_icluster_configuration->get_number_of_hosts();
                INFO_LOG << "Waiting for remote hosts to join...";
                sleep(5.0f);

                mi::Uint32 new_nb = irc.m_icluster_configuration->get_number_of_hosts();
                if (old_nb != new_nb)
                {
                    INFO_LOG << "Number of hosts changed from " << old_nb << " to " << new_nb << " while waiting.";
                }
            }
        }

        if(!m_prj->get("app::clock_pulse::interval").empty())
        {
            const mi::math::Vector_struct<mi::Float32, 2> entry = get_vec_float32_2(m_prj->get("app::clock_pulse::interval", "0 1"));
            nv::index_common::Clock_pulse_generator* clock_pulse =
                new Clock_pulse_generator(static_cast<mi::Float64>(entry.x), static_cast<mi::Float64>(entry.y));
            // clock_pulse->start();
            irc.m_iindex_session->set_clock(clock_pulse);
        }

        // Verify if the project file defines any affinity information: 
        {
            const mi::Uint32 use_affinity_only
                = get_uint32(m_prj->get("index::domain_subdivision::use_affinity_only", "1"));
            const mi::Uint32 nb_spatial_regions
                = get_uint32(m_prj->get("index::domain_subdivision::nb_spatial_regions", "0"));
            if (nb_spatial_regions > 0)
            {
                nv::index_common::Affinity_information* affinity_information = NULL;
                nv::index_common::Domain_specific_subdivision* domain_subdivision = NULL;
                if(use_affinity_only)
                {
                    affinity_information = new nv::index_common::Affinity_information;
                }
                else
                {
                    domain_subdivision = new nv::index_common::Domain_specific_subdivision;
                }

                for (mi::Uint32 i = 0; i < nb_spatial_regions; ++i)
                {
                    std::ostringstream bbox_string;
                    bbox_string << "index::domain_subdivision::spatial_region_" << i << "::bbox";

                    const mi::math::Bbox<mi::Float32, 3> subregion_bbox = get_bbox_float32_3(m_prj->get(bbox_string.str()));
                    std::ostringstream host_string;
                    host_string << "index::domain_subdivision::affinity_information_" << i << "::host";
                    const mi::Uint32 host_id = get_uint32(m_prj->get(host_string.str(), "0"));

                    std::ostringstream gpu_string;
                    gpu_string << "index::domain_subdivision::affinity_information_" << i << "::gpu";
                    const mi::Uint32 gpu_id = get_uint32(m_prj->get(gpu_string.str(), "0"));

                    if (affinity_information)
                    {
                        affinity_information->add_affinity(subregion_bbox, host_id, gpu_id);
                    }
                    else // domain_subdivision
                    {
                        domain_subdivision->add(subregion_bbox, host_id, gpu_id);
                    }
                }

                // The library takes ownership
                if (affinity_information)
                {
                    irc.m_iindex_session->set_affinity_information(affinity_information);
                }
                else // domain_subdivision
                {
                    irc.m_iindex_session->set_affinity_information(domain_subdivision);
                }
            }
        }

        // Setup session information
        bool created_new_session = false;
        if (!is_remote_viewer)
        {
            const std::string session_name = "nv_reference_viewer_session_0";
            irc.m_session_tag = dice_transaction->name_to_tag(session_name.c_str());
            if (!irc.m_session_tag.is_valid())
            {
                INFO_LOG << "Creating new viewer session";
                created_new_session = true;
                irc.m_session_tag = irc.m_iindex_session->create_session(
                    dice_transaction.get(),
                    session_name.c_str());
            }
            else
            {
                INFO_LOG << "Using previously stored viewer session";
            }
            assert(irc.m_session_tag.is_valid());
            VERBOSE_LOG << "The viewer session tag " << irc.m_session_tag.id
                        << " has identifier name: '" << session_name << "'";
        }

        mi::base::Handle<const nv::index::ISession> session(
            dice_transaction->access<nv::index::ISession>(irc.m_session_tag));
        assert(session.is_valid_interface());

        mi::math::Bbox<mi::Float32, 3> xyz_roi_bbox;
        if (created_new_session) // Initialize a new session only.
        {
            if (!m_prj->is_defined("index::region_of_interest"))
            {
                ERROR_LOG << "'index::region_of_interest' missing from project file!";
            }

            xyz_roi_bbox = get_bbox_float32_3(m_prj->get("index::region_of_interest"));

            // ROI is specified as "value range" in project file, so convert it to the bounding box format
            // that is used internally.
            xyz_roi_bbox.max += mi::math::Vector<mi::Float32, 3>(1.f);
            set_XYZ_global_region_of_interest_bbox(xyz_roi_bbox, session.get(), dice_transaction.get());
        }
        else
        {
            xyz_roi_bbox = get_XYZ_global_region_of_interest_bbox(session.get(), dice_transaction.get());
        }
        m_appdata->m_initial_roi = xyz_roi_bbox;

        // Initialize multi-view mode if enabled
        setup_multi_view_mode(irc, dice_transaction.get());

        // Access the IndeX rendering interface for scene queries such as mouse over (pick) operations
        VERBOSE_LOG << "Retrieving IIndex_scene_query";
        irc.m_iindex_scene_query = irc.m_iindex_if->get_api_component<nv::index::IIndex_scene_query>();
        assert(irc.m_iindex_scene_query.is_valid_interface());

        m_appdata->m_background_color
            = get_color(m_prj->get("app::canvas::background_color", "1 1 1 1"));

        // For default assigned colormap configuration
        m_appdata->set_current_colormap_index(get_uint32(m_prj->get("app::colormap::startup_colormap_index", "0")));

        // Create heightfield workflow functionality
        VERBOSE_LOG << "Creating heightfield workflow functionality";
        Heightfield_workflow_functionality::init(irc.m_iindex_session, irc.m_session_tag);

        {
            VERBOSE_LOG << "Loading colormaps";
            mi::base::Handle<const nv::index::IScene> scene(
                dice_transaction->access<const nv::index::IScene>(session->get_scene()));
            assert(scene.is_valid_interface());

            if (!load_all_colormap(scene.get(), *m_prj, dice_transaction.get()))
            {
                ERROR_LOG << "Failed to load colormap files.";
                return false;
            }
            if (!check_all_colormap_valid(dice_transaction.get()))
            {
                WARN_LOG << "Some colormaps are not valid.";
            }
        }

        if (created_new_session)
        {
            VERBOSE_LOG << "Creating the scene.";

            // This can be used to disable all volume rendering initially
            m_appdata->set_render_volume(get_bool(m_prj->get("index::is_render_seismic_volume", "1")));

            // Add scene data
            add_scene_graph(dice_transaction, irc.m_session_tag, *m_prj);
        }

        update_host_map(irc, *(m_appdata->peek_host_info()));

        irc.m_client_connection_issued = false;

        // Initialize the GL
        gl_initialize_opengl();
        gl_init_display_lists();   // merely for displaying the coord system

        // A progress callback passed to the rendering call.
        irc.m_progress_callback = new Progress_callback();
        assert(irc.m_progress_callback.is_valid_interface());

        // A frame info callback passed to the rendering call.
        irc.m_frame_info_callbacks = new Frame_info_callbacks();
        assert(irc.m_frame_info_callbacks.is_valid_interface());

        // Initialize OpenGL
#ifdef USE_OPENGL
        mi::base::Handle<Span_renderer_IF> span_buffer_gl(new Span_renderer_gl());
        irc.m_span_buffer.swap(span_buffer_gl);
#else
        mi::base::Handle<Span_renderer_IF> span_buffer_no_gl(new Span_renderer_no_gl());
        irc.m_span_buffer.swap(span_buffer_no_gl);
#endif  // USE_OPENGL

        gl_application_initialize(
            Nvindex_AppData::instance()->get_opengl_appdata()->get_opengl_application_buffer_ptr(),
            m_appdata->get_user_interaction(0)->get_examiner()->get_main_window_resolution());

        const std::string span_renderer_name = irc.m_span_buffer->get_class_name();
        VERBOSE_LOG << "Span renderer type: " << span_renderer_name;
        m_prj->insert("span_renderer_name", span_renderer_name);

        if (created_new_session)
        {
            // General configuration settings
            setup_index_configuration(irc, m_prj, session->get_config(), dice_transaction.get());

            // Set up camera and its examiner
            setup_camera(irc, dice_transaction.get());
        }
        else
        {
            // Create local camera object from the camera already stored in the scene
            mi::base::Handle<nv::index::IScene> scene(
                dice_transaction->edit<nv::index::IScene>(session->get_scene()));
            mi::neuraylib::Tag cam = scene->get_camera();
            Multiple_stereo_camera* ms_camera =
                m_appdata->get_user_interaction(0)->get_multiple_stereo_camera();
            ms_camera->set_main_camera(Stereo_camera(cam, cam));
        }

        // Set scene bounding box to the manipulator:
        // Default scene bbox: region of interest bounding box: ijk space
        m_appdata->get_user_interaction(0)->get_examiner()->set_scene_bbox(
            get_XYZ_global_region_of_interest_bbox(session.get(), dice_transaction.get()));

        // Set up performance monitoring
        set_performance_monitoring_by_project(irc, *m_prj, dice_transaction.get());

        // Initialize resolution (except the http server resolution)
        {
            const mi::Sint32_2 win_res = get_vec_sint32_2(m_prj->get("index::canvas_resolution"));
            request_canvas_resolution(win_res);

            // We have six canvases/buffers, but the http canvas is
            // not ready at here. the encoding canvas is not needed to
            // update here. So we initialize three canvases here.

            // 1. examiner resolution update
            set_examiner_window_resolution(win_res);

            // 2. index canvas (span buffer) resolution update
            irc.m_span_buffer->set_buffer_resolution(win_res);

            // 3. OpenGL offscreen context when running in OpenGL mode
            Nvindex_AppData::instance()->get_opengl_appdata()->
                resize_offscreen_context(win_res.x, win_res.y);
        }

        // Committing dice transaction and store all changes ...
        dice_transaction->commit();

        //----------------------------------------------------------------------

        dice_transaction = irc.get_current_scope()->create_transaction<mi::neuraylib::IDice_transaction>();
        assert(dice_transaction.is_valid_interface());

        // Create canvas buffers
        assert(irc.m_canvas_buffers.empty());
        irc.m_canvas_buffers.resize(2);
        // ... the first buffer is used for the main view
        const mi::math::Vector<mi::Sint32, 2> main_window_resolution =
            m_appdata->get_user_interaction(0)->get_examiner()->get_main_window_resolution();

        irc.m_canvas_buffers[0] = new Canvas_buffers(main_window_resolution.x, main_window_resolution.y);
        // ... the second buffer is used for multi-view tests
        irc.m_canvas_buffers[1] = new Canvas_buffers(1024, 1024);

        // Set default camera animator
        m_appdata->set_camera_animator_data(get_default_camera_animator(m_prj));
        m_appdata->peek_camera_animator_data()->reset(false); // no animation mode

        // Examiner setup
        if (created_new_session)
        {
            // scene translation setting
            const mi::math::Vector<mi::Float32, 3> translation = get_vec_float32_3(m_prj->get("index::translation", "0 0 0"));
            // scene scaling setting
            const mi::math::Vector_struct<mi::Float32, 3> scaling = get_vec_float32_3(m_prj->get("index::scaling", "1 1 1"));

            mi::Float32 angle;
            mi::math::Matrix<mi::Float32, 4, 4> rotation_mat(1.0f); // unit matrix
            const mi::Uint32 nb_rotations = get_uint32(m_prj->get("index::nb_rotations", "0"));

            for (mi::Uint32 i=0; i< nb_rotations; i++)
            {

                std::ostringstream v;
                v << "index::rotation_" << i;

                mi::math::Vector_struct<mi::Float32, 2> entry = get_vec_float32_2(m_prj->get(v.str(), "-1 -1"));
                mi::Sint32 axis = (mi::Sint32) entry.x;
                angle = entry.y;

                if (axis == 0)
                {
                    const mi::math::Matrix<mi::Float32, 4, 4> cur_rotation = get_rotation_matrix_about_x(angle);
                    rotation_mat = cur_rotation * rotation_mat;
                }
                else if (axis == 1)
                {
                    const mi::math::Matrix<mi::Float32, 4, 4> cur_rotation = get_rotation_matrix_about_y(angle);
                    rotation_mat = cur_rotation * rotation_mat;
                }
                else if (axis == 2)
                {
                    const mi::math::Matrix<mi::Float32, 4, 4> cur_rotation = get_rotation_matrix_about_z(angle);
                    rotation_mat = cur_rotation * rotation_mat;
                }
                else
                {
                    ERROR_LOG << "Invalid or missing rotation parameters";
                }
            }

            if (m_prj->is_defined("index::global_transform"))
            {
                mi::base::Handle<nv::index::IScene> scene(
                    dice_transaction->edit<nv::index::IScene>(session->get_scene()));
                assert(scene.is_valid_interface());
                scene->set_transform_matrix(
                    get_mat_float32_4_4(m_prj->get("index::global_transform")));
            }
            else
            {
                mi::base::Handle<nv::index::IScene> scene(
                    dice_transaction->edit<nv::index::IScene>(session->get_scene()));
                assert(scene.is_valid_interface());
                scene->set_transform_matrix(translation, rotation_mat, scaling);
            }

            if (get_bool(m_prj->get("app::examiner::predefined_view::startup_update", "0")))
            {
                mi::base::Handle<const nv::index::ISession> session(
                    dice_transaction->access<const nv::index::ISession>(irc.m_session_tag));
                assert(session.is_valid_interface());
                mi::base::Handle<const nv::index::IScene> scene(
                    dice_transaction->access<const nv::index::IScene>(session->get_scene()));
                assert(scene.is_valid_interface());
                INFO_LOG << "Adjusting the camera parameter according to predefined view "
                         << "(app::examiner::predefined_view::startup_update = yes).";
                Multiple_stereo_camera* ms_camera = 
                    m_appdata->get_user_interaction(0)->get_multiple_stereo_camera();
                m_appdata->get_user_interaction(0)->get_examiner()->set_predefined_view_parameter_to_camera(
                    session.get(), scene.get(), ms_camera, dice_transaction.get());
            }

            if (m_prj->get("app::examiner::initial_rotation_center::type", "<none>") == "roi_center")
            {
                if (!m_prj->is_defined("index::region_of_interest"))
                {
                    WARN_LOG << "initial_rotation_center::type = roi_center, but no index::region_of_interest defined. "
                             << "Will use [0 0 0 1024 1024 1024].";
                }
                const mi::math::Bbox<mi::Float32, 3> roi_bbox =
                    nv::index_common::get_bbox_float32_3(
                        m_prj->get("index::region_of_interest", "0 0 0 1024 1024 1024"));
                const mi::math::Vector<mi::Float32, 3> xyz_initial_rotation_center = roi_bbox.center();
                mi::base::Handle<nv::index::IScene> scene(
                    dice_transaction->edit<nv::index::IScene>(session->get_scene()));
                assert(scene.is_valid_interface());
                const mi::math::Matrix<mi::Float32, 4, 4> mat = scene->get_transform_matrix();
                const mi::math::Vector<mi::Float32, 3> initial_rotation_center =
                    mi::math::transform_point(mat, xyz_initial_rotation_center);
                m_appdata->get_user_interaction(0)->get_examiner()->set_examiner_rotation_center(initial_rotation_center);
                INFO_LOG << "Setting the initial rotation center " << initial_rotation_center << " from the region of interest.";
            }
        }

        if (update_session && irc.m_session_tag.is_valid())
        {
            // Update the session to the initial state of the scene
            mi::base::Lock::Block block(&(m_appdata->m_scene_update_lock));
            irc.m_iindex_session->update(irc.m_session_tag, dice_transaction.get());
        }

        // Commit the transaction
        dice_transaction->commit();

        //
        // When multiple scopes are used then localize some database elements that should not be
        // shared.
        //
        if (!m_appdata->is_enable_multi_view_mode() ||
            !get_bool(m_prj->get("app::viewport::share_camera", "on")))
        {
            for (size_t i=0; i < irc.m_extra_scopes.size(); ++i)
            {
                dice_transaction = irc.m_extra_scopes[i]->create_transaction<mi::neuraylib::IDice_transaction>();
                assert(dice_transaction.is_valid_interface());

                Multiple_stereo_camera* ms_camera =
                    m_appdata->get_user_interaction(0)->get_multiple_stereo_camera();

                // Cameras should be independent
                dice_transaction->localize(
                    ms_camera->get_current_main_camera_tag(),
                    mi::neuraylib::IDice_transaction::LOCAL_SCOPE);
                dice_transaction->localize(
                    ms_camera->get_ortho_camera_tag(),
                    mi::neuraylib::IDice_transaction::LOCAL_SCOPE);

                dice_transaction->commit();
            }
        }
    } // end of scene edit lock

    setup_application_configuration(m_prj);

    irc.m_initialized = true;

    // Start up RTMP and HTTP server if not done already
    if (irc.m_rtmp_server == 0)
    {
        if (!setup_all_server(irc))
            return false;
    }

    INFO_LOG << "Initialization complete.";

    return true;
}

//----------------------------------------------------------------------
bool Index_application::setup_remote(Nvindex_rendering_context& irc)
{
    if (!setup_common(irc))
        return false;

    Nvindex_library_accessor* p_nvindex_access = m_appdata->get_nvindex_library_accessor();

    irc.m_database = p_nvindex_access->get_interface()->get_api_component<mi::neuraylib::IDatabase>();
    assert(irc.m_database.is_valid_interface());

    irc.m_global_scope = irc.m_database->get_global_scope();
    assert(irc.m_global_scope.is_valid_interface());

    return true;
}

//----------------------------------------------------------------------
void Index_application::render_frame(
    Nvindex_rendering_context& irc, 
    mi::Uint32                 frame_num,
    bool                       need_compute)
{
    if (!is_rendering_needed())
    {
        // No rendering if no video stream is connected (unless "render_always" is enabled)
        return;
    }

    // render frame: common setup
    render_frame_common_setup_before_db_access(irc, frame_num);

    if (m_appdata->is_enable_multi_view_mode())
    {
        render_frame_multi_view(irc, frame_num, need_compute);
    }
    else
    {
        render_frame_single_view(irc, frame_num, need_compute);
    }
}

//----------------------------------------------------------------------
void Index_application::render_frame_remote(
    Nvindex_rendering_context& irc)
{
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        irc.m_viewing_scenario_scopes[0]->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());

    {
        mi::base::Lock::Block block(&(m_appdata->m_scene_update_lock));
        irc.m_iindex_session->update(irc.m_session_tag, dice_transaction.get());
    }

    // Create a new canvas and render into that canvas.
    // NOTE: Creating and deleting a canvas for every frame should generally be for performance reasons.
    //       Here, it just simplifies the source code of the reference viewer!
    mi::base::Handle<nv::index::IIndex_canvas> rendering_canvas(new Span_renderer_no_gl());
    (static_cast<Span_renderer_IF*>(rendering_canvas.get()))->set_background_color(
        m_appdata->m_background_color);
    const mi::math::Vector_struct<mi::Sint32, 2> canvas_resolution = { 1024, 1024, };
    (static_cast<Span_renderer_IF*>(rendering_canvas.get()))->set_buffer_resolution(canvas_resolution);
    mi::base::Handle<Frame_info_callbacks> frame_info_callbacks(new Frame_info_callbacks());

    // Create a NVIDIA IndeX rendering interface instance.
    // NOTE: Creating and deleting a rendering interface for every frame should be avoided
    //       1) for performance reasons and
    //       2) because the rendering interface instance keeps an internal state for optimizing
    //          rendering of subsequent frames.
    //       Here, it just simplifies the source code of the reference viewer!
    mi::base::Handle<nv::index::IIndex_rendering> index_rendering_interface(
        irc.m_iindex_if->create_rendering_interface());

    // Main NVIDIA IndeX rendering call. A special fallback mode when the first viewer is down.
    mi::base::Handle<nv::index::IFrame_results> frame_results(
        index_rendering_interface->render(
            irc.m_session_tag,
            rendering_canvas.get(),
            dice_transaction.get(),
            NULL,               // Progress callback
            frame_info_callbacks.get(),
            true,               // Composite immediately
            NULL));             // OpenGL application buffer

    const mi::base::Handle<nv::index::IError_set> err_set(frame_results->get_error_set());
    if (err_set->any_errors())
    {
        std::ostringstream os;

        const mi::Uint32 nb_err = err_set->get_nb_errors();
        for (mi::Uint32 e = 0; e < nb_err; ++e)
        {
            if (e != 0) os << '\n';
            os << err_set->get_error(e)->get_error_string();
        }

        ERROR_LOG << "IIndex_rendering rendering call failed with the following error(s): " << '\n'
                  << os.str();

        return;
    }

    std::vector<Frame_info_callbacks::Dynamic_alloc_event> dyn_alloc_events =
        frame_info_callbacks->get_dynamic_allocation_events();
    if (!dyn_alloc_events.empty())
    {
        INFO_LOG << "Dynamic allocation events occurred during the IIndex_rendering rendering call:";
        for (mi::Size i = 0; i < dyn_alloc_events.size(); ++i)
        {
            INFO_LOG << " * event " << i << ":\n"
                     << "   - host id:           " << dyn_alloc_events[i].m_host_id << "\n"
                     << "   - device id:         " << dyn_alloc_events[i].m_device_id << "\n"
                     << "   - requested size:    " << dyn_alloc_events[i].m_memory_allocation_size << "\n"
                     << "   - available memory:  " << dyn_alloc_events[i].m_memory_available << "\n"
                     << "   - memory freed up:   " << dyn_alloc_events[i].m_memory_freed_up;
        }
    }

    std::vector<Frame_info_callbacks::Device_reset_event> dev_reset_events =
        frame_info_callbacks->get_device_reset_events();
    if (!dev_reset_events.empty())
    {
        INFO_LOG << "Device reset events occurred during the IIndex_rendering rendering call:";
        for (mi::Size i = 0; i < dev_reset_events.size(); ++i)
        {
            INFO_LOG << " * event " << i << ":\n"
                     << "   - host id:           " << dev_reset_events[i].m_host_id << "\n"
                     << "   - device id:         " << dev_reset_events[i].m_device_id;
        }
    }

    std::vector<mi::base::Handle<nv::index::IBalancing_operations> > balancing_ops =
        frame_info_callbacks->get_workload_balancing_operations();
    if (!balancing_ops.empty())
    {
        for (mi::Size i = 0; i < balancing_ops.size(); ++i)
        {
            const mi::Uint32 nb_operations = balancing_ops[i]->get_nb_operation();
            const char* event_desc = balancing_ops[i]->get_event_description();

            INFO_LOG << i << ".) Balancing event: '" << event_desc << "'" 
                     << " results in " << nb_operations << " balancing operations: ";
            for (mi::Uint32 j=0; j<nb_operations; ++j)
            {
                const nv::index::IBalancing_operation* op = balancing_ops[i]->get_operation(j);
                INFO_LOG << "Balance operation " << j << ": " << op->get_description();
            }
        }
    }

    // Prepare frame buffer to be passed into the video stream or for taking a snapshot image (tripple buffer)
    mi::base::Handle<Canvas> canvas(irc.m_canvas_buffers[0]->get_render_canvas());
    assert(canvas);
    
    // Copy the result pixels to the canvas (tripple buffer)
    copy_result_pixel_to_canvas(static_cast<Span_renderer_IF*>(rendering_canvas.get()), canvas.get());

    // Trigger internal buffer swapping
    irc.m_canvas_buffers[0]->rendering_finished();

    dice_transaction->commit();
}

//----------------------------------------------------------------------
bool Index_application::run_remote(Nvindex_rendering_context& irc)
{
    mi::base::Handle<mi::neuraylib::INetwork_configuration> inetwork_configuration(
        irc.m_iindex_if->get_api_component<mi::neuraylib::INetwork_configuration>());
    assert(inetwork_configuration.is_valid_interface());
    if (inetwork_configuration->get_mode() == mi::neuraylib::INetwork_configuration::MODE_OFF)
    {
        ERROR_LOG << "No network mode available, exit the remote service...";
        return false;
    }

    assert(irc.m_iindex_if.is_valid_interface());
    mi::base::Handle<nv::index::ICluster_configuration> rendering_properties_query(
        irc.m_iindex_if->get_api_component<nv::index::ICluster_configuration>());
    mi::Uint32 number_of_hosts = rendering_properties_query->get_number_of_hosts();

    INFO_LOG << "****************************************************************************";
    INFO_LOG << " Running in remote service mode, waiting for connections, "
             << "host id = " << rendering_properties_query->get_local_host_id();
    INFO_LOG << "****************************************************************************";

    const bool started_as_viewing_scenario = get_bool(m_prj->get("app::remote_viewing_capability", "no"));
    while(true)
    {
        mi::Uint32 new_number_of_hosts = rendering_properties_query->get_number_of_hosts();

        // Finish on lost host
        if (new_number_of_hosts < number_of_hosts)
        {
            INFO_LOG << "At least one host has left the cluster - had "
                     << number_of_hosts << " hosts, only " << new_number_of_hosts << " left";
            const bool fail_safety_enabled = get_bool(m_prj->get("app::enable_fail_safety", "no"));
            if (!fail_safety_enabled)
            {
                INFO_LOG << "Leaving because at least one host has left - had "
                         << number_of_hosts << " hosts, only " << new_number_of_hosts << " left";
                break;
            }
        }

        if (started_as_viewing_scenario && irc.m_client_connection_issued)
        {
            if (!irc.m_session_tag.is_valid())
            {
                INFO_LOG << "No previously stored session available. "
                         << "Will continue in remote service mode.";
            }
            else
            {
                INFO_LOG << "Finished remote loop to serve the requested rendering.";
                return false;
            }
        }

        number_of_hosts = new_number_of_hosts;
        sleep(0.1f); // seconds
    }

    INFO_LOG << "Finished remote rendering loop.";
    return true;
}

//----------------------------------------------------------------------
bool Index_application::process_command_line_arguments(
    mi::Sint32                                      argc,
    char**                                          argv,
    String_dict &                                   result_opt)
{
    // default for command line option
    String_dict opt = get_default_command_line_options();

    // command line option
    String_dict com_line_opt;
    const bool is_true_key_only = true;
    string_array_to_string_dict(argc, argv, com_line_opt, is_true_key_only);

    // check the project file. should be there
    if (!com_line_opt.is_defined("app::project_file"))
    {
        ERROR_LOG << "No project file specified. use -app::project_file option.";
        return false;
    }

    // project file option
    String_dict proj_opt;
    const std::string proj_fname = com_line_opt.get("app::project_file");
    if (!load_application_project_file(proj_fname, proj_opt))
    {
        ERROR_LOG << "Failed to load application project file [" << proj_fname << "].";
        return false;
    }

    // final option. priority is: default < app::project_file < app::project_file::0 < .. < command_line
    opt.insert_all(proj_opt);      // override the default with proj_opt

    std::vector< String_dict > option_vec;
    option_vec.push_back(proj_opt);
    option_vec.push_back(com_line_opt);
    const String_dict project_fname_opt =
        get_prefix_entry_key("app::project_file::", option_vec);
    const String_dict extra_project_file_entry =
        get_extra_project_option(project_fname_opt);
    opt.insert_all(extra_project_file_entry); // in case more project file are specified

    if (opt.is_defined("config_file"))
    {
        // obsolete the config file (config.prj)
        WARN_LOG << "The use of an external config file is obsolete and no longer in effect.";
    }

    opt.insert_all(com_line_opt);  // override again with the command line option

    const std::string service_mode = opt.get("index::service");
    if ((service_mode != "rendering_and_compositing") &&
        (service_mode != "rendering")                 &&
        (service_mode != "compositing")               &&
        (service_mode != "none"))
    {
        ERROR_LOG << "Unknown service: '" << service_mode
                  << "'. Assuming 'rendering and compositing'";
        opt.insert("index::service", "rendering_and_compositing"); // update the service mode
    }

    const std::string startup_state = opt.get("app::host_mode");
    if ((startup_state != "remote_service") && (startup_state != "viewer"))
    {
        ERROR_LOG << "unknown -mode [" << startup_state << "]. set as a viewer mode.";
    }

    INFO_LOG << "Running in " << startup_state << " mode (" << opt.get("index::service") << ")";

    result_opt = opt;

    std::ostringstream default_level;
    default_level << mi::base::MESSAGE_SEVERITY_VERBOSE;
    if (get_uint32(result_opt.get("dice::verbose", default_level.str())) >= mi::base::MESSAGE_SEVERITY_DEBUG)
    {
        std::stringstream sstr;
        result_opt.write(sstr, "info: [debug]: ");
        INFO_LOG << "[debug]: options\n" << sstr.str();
    }

    // check key has white space and if so, warn them.
    for (String_dict::const_iterator si = result_opt.begin();
        si != result_opt.end(); ++si)
    {
        const bool is_valid = nv::index_common::is_valid_key(si->first);
        if (!is_valid)
        {
            WARN_LOG << "Found an invalid key [" << si->first << "]. Please check the project file "
                     << "(the key might have invalid characters (e.g., white space).)";
        }
    }

    if (get_bool(result_opt.get("app::wait_for_cin", "0"))) // for debugging
    {
        INFO_LOG << "waiting for cin. please type.";
        std::string x;
        std::cin >> x;
    }

    return true;
}

//----------------------------------------------------------------------
void Index_application::save_span_buffer_snapshot(
    Span_renderer_IF * span_buffer,
    mi::Sint32 frame_num)
{
    assert(span_buffer != 0);

    Canvas canvas(1, 1);        // for temporary copy
    copy_result_pixel_to_canvas(span_buffer, &canvas);

    const mi::math::Vector<mi::Sint32, 2> span_buffer_size = span_buffer->get_buffer_resolution();

    const std::string fname = get_snapshot_filename_by_mode(m_appdata, frame_num);
    const mi::base::Handle<mi::neuraylib::ITile> tile(canvas.get_tile(0, 0, 0));
    const bool ret = save_pixel_buffer_to_ppm_file(span_buffer_size.x, span_buffer_size.y,
                                                   reinterpret_cast< mi::Uint8 *>(tile->get_data()),
                                                   fname);
    if (ret)
    {
        INFO_LOG << "Saved snapshot '" << fname << "'";
    }
    else
    {
        ERROR_LOG << "Failed to save snapshot '" << fname << "'";
    }
}

//----------------------------------------------------------------------
void Index_application::shutdown_all_server(Nvindex_rendering_context& irc)
{
    // Shut down servers.
    // Requests currently being handled will be served until they are done.
    if (irc.m_bridge_application)
    {
        irc.m_bridge_application->close();
    }

    if (irc.m_http_server)
    {
        irc.m_http_server->shutdown();
    }

    if (irc.m_rtmp_server)
    {
        irc.m_rtmp_server->shutdown();
    }

    if (irc.m_html5_video_stream_manager)
    {
        irc.m_html5_video_stream_manager->shutdown();
        delete irc.m_html5_video_stream_manager;
        irc.m_html5_video_stream_manager = 0;
    }
}

//----------------------------------------------------------------------
void Index_application::shutdown()
{
    // shutdown/delete before DiCE shutdown
    gl_application_shutdown();
    Heightfield_workflow_functionality::delete_instance();
    Forwarding_logger_factory::delete_instance();
    RTC_parameter_buffer_manip::delete_instance();
    nv::index_common::Forwarding_logger_factory::instance()->shutdown();        

    // Note: This shutdown the IndeX and DiCE via
    // nvindex_library_accessor deletion.
    Nvindex_AppData::delete_instance();
}

//----------------------------------------------------------------------
//======================================================================
//----------------------------------------------------------------------
bool Index_application::setup_common(Nvindex_rendering_context& irc)
{
    check_project_compatibility(*m_prj);

    // nvindex library access
    Nvindex_library_accessor* p_nvindex_access = new Nvindex_library_accessor();
    assert(p_nvindex_access != 0);
    m_appdata->set_nvindex_library_accessor(p_nvindex_access); // keep this for delete later
    if (!p_nvindex_access->initialize(*m_prj, irc))
    {
        ERROR_LOG << "Initialization of the NVIDIA IndeX library failed.";
        // here, we do nothing, but later we do shutdown.
        return false;
    }

    // For keeping track of the nvindex library interface
    irc.m_iindex_if = p_nvindex_access->get_interface();
    assert(irc.m_iindex_if.is_valid_interface());

    // Store dice version to m_prj (used by performance monitor later)
    m_prj->insert("app::info::dice_version",   p_nvindex_access->get_interface()->get_dice_version());
    m_prj->insert("app::info::index_revision", p_nvindex_access->get_interface()->get_revision());
    m_prj->insert("app::info::driver_version", p_nvindex_access->get_interface()->get_nvidia_driver_version());
    m_prj->insert("app::info::cuda_version",
                  nv::index_common::to_string(p_nvindex_access->get_interface()->get_cuda_runtime_version()));

    // project file startup message out
    if (m_prj->is_defined("app::startup_message::warn"))
    {
        const std::string startup_warn = m_prj->get("app::startup_message::warn");
        if (!(startup_warn.empty()))
        {
            WARN_LOG << startup_warn;
        }
    }
    if (m_prj->is_defined("app::startup_message::info"))
    {
        const std::string startup_info = m_prj->get("app::startup_message::info");
        if (!(startup_info.empty()))
        {
            INFO_LOG << startup_info;
        }
    }

    return true;
}

//----------------------------------------------------------------------
bool Index_application::setup_all_server(Nvindex_rendering_context& irc)
{
    const bool rtmp_ok        = setup_rtmp_http(irc);
    const bool dice_bridge_ok = setup_dice_bridge(irc);
    const bool websocket_ok   = setup_html5_video_stream(irc);

    return rtmp_ok && dice_bridge_ok && websocket_ok;
}

//----------------------------------------------------------------------
std::string Index_application::get_http_address() const
{
    const std::string http_listen  = m_prj->get("dice::http_listen");
    const std::string http_port    = m_prj->get("dice::http_port");
    const std::string http_address = http_listen + ":" + http_port;

    return http_address;
}

//----------------------------------------------------------------------
bool Index_application::setup_rtmp_http(
    Nvindex_rendering_context & irc)
{
    const bool is_rtmp_video_streaming = get_bool(m_prj->get("dice::rtmp_video_streaming::enabled", "false"));
    if (is_rtmp_video_streaming)
    {
        m_appdata->set_any_video_stream_enabled(is_rtmp_video_streaming);
    }
    else
    {
        return true;            // still valid (other video stream possible)
    }
    INFO_LOG << "Enabled RTMP video streaming.";

    // Create an RTMP server instance
    irc.m_rtmp_factory = irc.m_iindex_if->get_api_component<mi::rtmp::IFactory>();
    if (!irc.m_rtmp_factory.get())
    {
        ERROR_LOG << "No RTMP server exposed. Aborting ...";
        return false;
    }
    irc.m_rtmp_server = irc.m_rtmp_factory->create_server();
    assert(irc.m_rtmp_server.is_valid_interface());

    // Install our RTMP connect handler
    // mandatory project key list when video streaming is on.
    char const * const p_necessary_key[] = { 
        "dice::rtmp_video_streaming::video_codec", 
        "dice::rtmp_video_streaming::video_bitrate",
        0,
    };
    for (mi::Sint32 i = 0; p_necessary_key[i] != 0; ++i)
    {
        if (m_prj->get("dice::rtmp_video_streaming::video_codec", "") == "")
        {
            WARN_LOG << "Video streaming is on, but " << p_necessary_key[i] << " is not specified.";
        }
    }

    // Start RTMP server
    assert(m_prj->is_defined("dice::rtmp_video_streaming::port"));
    assert(m_prj->is_defined("dice::rtmp_video_streaming::listen"));
    const std::string rtmp_port    = m_prj->get("dice::rtmp_video_streaming::port");
    const std::string rtmp_listen  = m_prj->get("dice::rtmp_video_streaming::listen");
    const std::string rtmp_address = rtmp_listen + ":" + rtmp_port;
    if (irc.m_rtmp_server->start(rtmp_address.c_str()) != 0)
    {
        ERROR_LOG << "Could not start RTMP server on " << rtmp_address << ", aborting.";
        return false;
    }

    // Create an HTTP server instance
    assert(irc.m_iindex_if.is_valid_interface());
    irc.m_http_factory = irc.m_iindex_if->get_api_component<mi::http::IFactory>();
    if (!irc.m_http_factory.get())
    {
        ERROR_LOG << "No HTTP server exposed. Aborting ...";
        return false;
    }

    assert(irc.m_http_factory.is_valid_interface());
    // Create the http server for standard application
    irc.m_http_server = irc.m_http_factory->create_server();

    // Install the HTTP request handle for the standard application
    const std::string flash_client_file = m_prj->get("dice::rtmp_video_streaming::flash_client");
    if (flash_client_file.empty() && is_rtmp_video_streaming)
    {
        WARN_LOG << "No flash client file specified. Please check the project file "
                 << "'dice::rtmp_video_streaming::flash_client = FILENAME' line.";
    }

    std::vector<mi::math::Vector<mi::Uint32, 2> > resolutions;
    resolutions.push_back(m_appdata->get_user_interaction(0)->
                          get_examiner()->get_main_window_resolution_uint32());
    resolutions.push_back(mi::math::Vector<mi::Uint32, 2>(1024, 1024));

    irc.m_http_request_handler = new HTTP_request_handler(
        flash_client_file, resolutions, rtmp_port);
    assert(irc.m_http_request_handler.is_valid_interface());

    irc.m_http_server->install(irc.m_http_request_handler.get());

    // Activate HTTP authentication if a username is specified in the configuration
    std::string session_cookie;
    std::string session_cookie_guest;
    if (!m_prj->get("dice::http_auth_user").empty())
    {
        irc.m_http_request_handler->set_authentication(
            m_prj->get("dice::http_auth_user"),
            m_prj->get("dice::http_auth_password"));

        // Guest authentication
        if (!m_prj->get("dice::http_auth_guest_password").empty())
        {
            irc.m_http_request_handler->set_guest_authentication(
                m_prj->get("dice::http_auth_user"),             // same user
                m_prj->get("dice::http_auth_guest_password"));  // different password
        }

        // A random session cookie is used for authentication of the RTMP connection.
        // It is enabled together with HTTP authentication.
        session_cookie       = irc.m_http_request_handler->get_session_cookie();
        session_cookie_guest = irc.m_http_request_handler->get_session_cookie_guest();
    }

    // Install RTMP connect handler, passing the session cookie
    mi::base::Handle<Connect_event_handler> conhnd(
        new Connect_event_handler(irc, *m_prj, session_cookie, session_cookie_guest));
    irc.m_connect_handler = conhnd;
    irc.m_rtmp_server->install(irc.m_connect_handler.get());

    // Optionally use SSL encryption
    const bool use_ssl = get_bool(m_prj->get("dice::ssl", "no"));

    // Start HTTP(S) server for the standard application
    const std::string http_address = get_http_address();

    mi::Sint32 result;
    if (use_ssl)
    {
        const std::string ssl_password = m_prj->get("dice::ssl_password");
        result = irc.m_http_server->start_ssl(
            http_address.c_str(),
            m_prj->get("dice::ssl_cert_file").c_str(),
            m_prj->get("dice::ssl_private_key_file").c_str(),
            ssl_password.empty() ? 0 : ssl_password.c_str());
        if (result != 0)
        {
            std::string error = "unknown";
            switch (result)
            {
            case -1: error = "listen_address already in use";           break;
            case -2: error = "The method has already been called";      break;
            case -3: error = "Invalid certificate";                     break;
            case -4: error = "Invalid private key";                     break;
            case -5: error = "Mismatch of certificate and private key"; break;
            case -6: error = "Cannot load SSL module";                  break;
            }
            ERROR_LOG << "start_ssl() failed: " << error;
        }
    }
    else
    {
        result = irc.m_http_server->start(http_address.c_str());
    }

    if (result != 0)
    {
        ERROR_LOG << "Could not start HTTP" << (use_ssl ? "S" : "") << " server on " << http_address << ". "
                  << "Aborting...";
        return false;
    }
    INFO_LOG << "HTTP" << (use_ssl ? "S" : "") << " server listens on " << http_address;

    INFO_LOG << "****************************************************************************";
    std::string url = "http" + std::string(use_ssl ? "s" : "") + "://" + get_host_name();
    const std::string port = m_prj->get("dice::http_port");
    if ((!use_ssl && port != "80") || (use_ssl && port != "443"))
    {
        url += ":" + port;
    }
    INFO_LOG << " Server running at " << url;
    INFO_LOG << "****************************************************************************";

    return true;
}

//----------------------------------------------------------------------
bool Index_application::setup_dice_bridge(Nvindex_rendering_context& irc)
{
    const bool bridge_enabled = get_bool(m_prj->get("dice::bridge_video_streaming::enabled", "false"));
    if (bridge_enabled)
    {
        m_appdata->set_any_video_stream_enabled(bridge_enabled);
    }
    else
    {
        return true;            // still valid (other video stream possible)
    }
    INFO_LOG << "Enabled DiCE bridge.";

    // start DiCE Bridge Server
    mi::base::Handle<mi::neuraylib::INeuray> dice(
        irc.m_iindex_if->get_dice_interface());
    assert(dice.is_valid_interface());
            
    // Access the API component for the Bridge server and create and application.
    const std::string http_address = get_http_address();
    irc.m_bridge_server = dice->get_api_component<mi::bridge::IBridge_server>();
    INFO_LOG << "DiCE-Bridge server listens on " << http_address;
            
    irc.m_bridge_application = irc.m_bridge_server->create_application(
        m_prj->get("dice::bridge_video_streaming::app_path").c_str(), irc.m_http_server.get());
                
    if ( irc.m_bridge_application.is_valid_interface() ) 
    {
        // Configure the application.
        irc.m_bridge_application->set_disk_cache( m_prj->get("dice::bridge_video_streaming::disk_cache_path").c_str());
                
        // register Render_batch_job factory
        const bool is_output_command =
            get_bool(m_prj->get("dice::bridge_video_streaming::is_output_command", "no"));
        mi::base::Handle<mi::neuraylib::IUser_class_factory> bridge_job_serverside_factory(
            new Bridge_job_serverside_factory(&irc, is_output_command));

        // register bridge job taking care of client/server communication
        if (!irc.m_bridge_application->register_job(
                Bridge_job_serverside::IID(),
                bridge_job_serverside_factory.get()))
        {
            ERROR_LOG << "Failed to register factory for Bridge communication job";
            return 0;
        }
                    
        // Bridge_init_job_serverside
        {
            // register bridge_init_job_serverside factory
            mi::base::Handle<mi::neuraylib::IUser_class_factory> bridge_init_job_serverside_factory( 
                new Bridge_init_job_serverside_factory(&irc));

            // register server initialization job
            if (!irc.m_bridge_application->register_job(
                    Bridge_init_job_serverside::IID(), 
                    bridge_init_job_serverside_factory.get()))
            {
                ERROR_LOG << "Failed to register factory for Bridge initialization job";
                return 0;
            }
        }

        irc.m_bridge_application_session_handler = 
            new Application_session_handler(m_prj->get("dice::bridge_video_streaming::security_token", ""));
                
        if (irc.m_bridge_application->set_session_handler(irc.m_bridge_application_session_handler.get()) != 0)
        {
            ERROR_LOG << "It could not set application session handler. Aborting...";
            return false;
        }

        // Run the application
        irc.m_bridge_application->open();

    }

    return true;
}

//----------------------------------------------------------------------
bool Index_application::setup_html5_video_stream(Nvindex_rendering_context& irc)
{
    const bool html5_enabled = get_bool(m_prj->get("dice::html5_video_streaming::enabled", "false"));
    if (html5_enabled)
    {
        m_appdata->set_any_video_stream_enabled(html5_enabled);
    }
    else
    {
        return true; // No html5 video streaming, but still valid. (other video stream possible)
    }
    INFO_LOG << "Enabled html5 video streaming.";

    assert(irc.m_html5_video_stream_manager == 0);

    irc.m_html5_video_stream_manager = new Html5_video_stream_manager();
    const bool ret = irc.m_html5_video_stream_manager->init(&irc, m_prj);
    
    return ret;
}

//----------------------------------------------------------------------
void Index_application::setup_camera(
    Nvindex_rendering_context&        irc,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    // If camera already exists then just update the settings
    Multiple_stereo_camera* ms_camera = m_appdata->get_user_interaction(0)->get_multiple_stereo_camera();
    const bool camera_exists = ms_camera->get_current_main_camera_tag().is_valid();

    // TODO camera set: use Camera_parameter class
    mi::math::Vector<mi::Float32, 3> from = get_vec_float32_3(m_prj->get("index::camera::from", "0.0 0.0 -5.0"));
    mi::math::Vector<mi::Float32, 3> dir  = get_vec_float32_3(m_prj->get("index::camera::dir",  "0.0 0.0  1.0"));

    // This only for compatibility for the lookat position set.
    if (m_prj->is_defined("index::camera::to"))
    {
        INFO_LOG << "index::camera::to is deprecated, please use index::camera::dir, which is viewing direction vector. ";
        const mi::math::Vector<mi::Float32, 3> to = get_vec_float32_3(m_prj->get("index::camera::to",   "0.0 0.0 0.0"));
        dir = to - from;
        if (!dir.normalize())
        {
            ERROR_LOG << "Cannot normalize viewing direction. Use (0, 0, 1).";
            dir = mi::math::Vector<mi::Float32, 3>(0.0f, 0.0f, 1.0f);
        }
    }
    mi::math::Vector<mi::Float32, 3> up = get_vec_float32_3(m_prj->get("index::camera::up", "0.0 1.0 0.0"));
    mi::Float32 aspect   = get_float32(m_prj->get("index::camera::aspect",     "1.0"));
    mi::Float32 aperture = get_float32(m_prj->get("index::camera::aperture",   "0.033"));
    mi::Float32 focal    = get_float32(m_prj->get("index::camera::focal",      "0.03"));
    mi::Float64 clip_min = get_float32(m_prj->get("index::camera::clip_min",   "0.01"));
    mi::Float64 clip_max = get_float32(m_prj->get("index::camera::clip_max",   "100.0"));
    bool orthographic    = get_bool(m_prj->get("index::camera::orthographic", "no"));

    if (clip_max <= clip_min)
    {
        ERROR_LOG << "setup_camera(): illegal clipping plane, clip_min must be < clip_max. "
                  << "Setting to 0.01, 100.0";
        clip_min = 0.01;
        clip_max = 100.0;
    }

    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<const nv::index::ISession>(irc.m_session_tag));
    assert(session.is_valid_interface());

    if (!camera_exists)
    {
        // Create main camera
        ms_camera->set_main_camera(
            Stereo_camera(session.get()->create_camera(dice_transaction),
                          session.get()->create_camera(dice_transaction)));

        ms_camera->set_ortho_camera_tag(
            session.get()->create_camera(dice_transaction, nv::index::IOrthographic_camera::IID()));
    }

    // Set camera parameters
    {
        mi::base::Handle< nv::index::IPerspective_camera > cam(
            dice_transaction->edit< nv::index::IPerspective_camera >(
                ms_camera->get_current_main_camera_tag()));
        assert(cam.is_valid_interface());
        cam->set(from, dir, up);
        cam->set_aperture(aperture);
        cam->set_aspect(aspect);
        cam->set_focal(focal);
        cam->set_clip_min(clip_min);
        cam->set_clip_max(clip_max);

        mi::base::Handle< nv::index::IOrthographic_camera > ortho_cam(
            dice_transaction->edit< nv::index::IOrthographic_camera >(
                ms_camera->get_ortho_camera_tag()));
        assert(ortho_cam.is_valid_interface());
        ortho_cam->set(from, dir, up);
        ortho_cam->set_aperture(aperture);
        ortho_cam->set_aspect(aspect);
        ortho_cam->set_clip_min(clip_min);
        ortho_cam->set_clip_max(clip_max);

        // Cam handle is released here. We can have another
        // access/edit, otherwise crash in initialize_storage_camera_by_main_camera
    }

    ms_camera->set_use_ortho_camera(orthographic);

    if (!camera_exists)
    {
        // Add three storage cameras (GUI has 3 storage cameras)
        for (mi::Uint32 i = 0; i < 3; ++i)
        {
            ms_camera->append_camera(
                Stereo_camera(session.get()->create_camera(dice_transaction),
                              session.get()->create_camera(dice_transaction)));
        }

        // Init all storage cameras
        ms_camera->initialize_storage_camera_by_main_camera(dice_transaction);
    }
}

//----------------------------------------------------------------------
void Index_application::setup_application_configuration(
    String_dict* prj)
{
    m_appdata->set_snapshot_filename(
        prj->get("app::camera::snapshot_filename", "snap_%03d.ppm"));
    m_appdata->set_snapshot_sequence_mode(
        get_bool(prj->get("app::snapshot::is_sequence_mode", "1")));
    m_appdata->set_animation_mode(
        get_bool(m_prj->get("app::animation::enabled", "0")));

    m_appdata->m_show_rendering_statistics_overlay = get_bool(
        prj->get("app::statistics::enable_opengl_rendering_statistics", "0"));
    m_appdata->m_show_compositing_statistics_overlay = get_bool(
        prj->get("app::statistics::enable_opengl_compositing_statistics", "0"));
    m_appdata->m_show_all_performance_statistics_overlay = get_bool(
        prj->get("app::statistics::enable_opengl_performance_statistics", "0"));
    m_appdata->m_render_always = get_bool(prj->get("app::render_always", "yes"));
    m_appdata->set_pick_operation_enabled(get_bool(prj->get("app::enable_pick_operations", "yes")));

    m_appdata->set_show_bounding_boxes_overlay(
        get_bool(prj->get("overlay::show_bounding_boxes", "0")));
    m_appdata->set_show_region_of_interest_overlay(
        get_bool(prj->get("overlay::show_region_of_interest", "0")));
    m_appdata->set_show_horizontal_spans_overlay(
        get_bool(prj->get("overlay::show_horizontal_spans", "0")));
    m_appdata->set_show_coordinate_system_overlay(
        get_bool(prj->get("overlay::show_coordinate_system", "1")));

    // experimental features
    m_appdata->set_enable_opengl_z_buffer_for_rendering(
        get_bool(prj->get("app::experimental::enable_opengl_z_buffer_for_rendering", "0")));
    if (m_appdata->is_enable_opengl_z_buffer_for_rendering())
    {
        if (m_appdata->is_use_opengl())
        {
            INFO_LOG << "Enabled experimental feature OpenGL z-buffer integration.";
        }
        else
        {
            WARN_LOG << "Experimental feature OpenGL z-buffer integration requires OpenGL support, "
                     << "which is not supported in this build. You may see unexpected results. "
                     << "Maybe you forgot the --gl option?";
        }
    }

    const std::string compute_technique = prj->get("app::compute_technique");
    if (!compute_technique.empty())
    {
#ifdef NVINDEX_HAS_COMPUTE
        INFO_LOG << "Activating compute technique '" << compute_technique << "'";
        m_compute_wrapper.init(
            compute_technique,
            prj,
            m_appdata->get_nvindex_library_accessor()->get_interface().get());
#else
        ERROR_LOG << "Compute technique '" << compute_technique << "' was requested, "
                  << "but this application was compiled without NVINDEX_HAS_COMPUTE!";
#endif
    }

    // Option to quit before after a certain number of frames or before rendering the first image
    // ("-quit 0")
    m_appdata->set_quit_after_frame_number(get_sint32(m_prj->get("quit", "-1")));
}

//----------------------------------------------------------------------
void Index_application::setup_index_configuration(
    Nvindex_rendering_context&        irc,
    String_dict*                      prj,
    mi::neuraylib::Tag                config_tag,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    mi::base::Handle<nv::index::IConfig_settings> config_settings(
        dice_transaction->edit<nv::index::IConfig_settings>(config_tag));

    // Compression settings for image transfer
    nv::index::IConfig_settings::Data_transfer_config cfg = config_settings->get_data_transfer_config();
    cfg.span_compression_level = get_uint32(prj->get("dice::network::span_compression", "1"));
    cfg.tile_compression_level = get_uint32(prj->get("dice::network::tile_compression", "1"));
    cfg.span_image_encoding = get_bool(prj->get("dice::network::span_encoding", "yes"));
    cfg.tile_image_encoding = get_bool(prj->get("dice::network::tile_encoding", "yes"));
    cfg.span_alpha_channel = get_bool(prj->get("dice::network::span_alpha", "no"));
    cfg.span_compression_mode = get_uint32(prj->get("dice::network::span_compression_mode", "0"));
    cfg.tile_compression_mode = get_uint32(prj->get("dice::network::tile_compression_mode", "0"));
    cfg.span_compression_threads = get_uint32(prj->get("dice::network::span_compression_threads", "0"));
    cfg.tile_compression_threads = get_uint32(prj->get("dice::network::tile_compression_threads", "0"));
    cfg.z_buffer_compression_level = get_uint32(prj->get("dice::network::z_buffer_compression", "0"));
    cfg.z_buffer_compression_mode = get_uint32(prj->get("dice::network::z_buffer_compression_mode", "0"));
    cfg.z_buffer_compression_threads = get_uint32(prj->get("dice::network::z_buffer_compression_threads", "0"));
    config_settings->set_data_transfer_config(cfg);

    // Force data uploading
    config_settings->set_force_data_upload(get_bool(prj->get("index::force_data_uploading", "no")));

    // CPU/GPU rendering mode
    if (get_bool(prj->get("index::allow_cpu_fallback", "no")))
        config_settings->set_rendering_mode(nv::index::IConfig_settings::GPU_OR_CPU_RENDERING);

    if (get_bool(prj->get("index::force_cpu_rendering", "no")))
        config_settings->set_rendering_mode(nv::index::IConfig_settings::CPU_RENDERING);

    //
    // Horizontal spans
    //

    // This will be set after rendering the first frame
    const mi::Uint32 horizontal_spans = get_sint32(prj->get("index::horizontal_spans", "1").c_str());
    const mi::Sint32 horizontal_spans_per_machine = get_sint32(prj->get("index::horizontal_spans_per_machine", "0"));

    if (get_bool(prj->get("index::enable_automatic_span_control")))
    {
        config_settings->set_automatic_span_control(true);
        mi::Sint32 max_spans_per_machine = get_sint32(prj->get("index::max_spans_per_machine", "4"));
        if (max_spans_per_machine < 1)
        {
            WARN_LOG << "Max spans per machine: ["<<max_spans_per_machine<<"] is invalid and should not be negative, current value is ignored";
            max_spans_per_machine = 1;
        }
        config_settings->set_max_spans_per_machine(mi::Uint32(max_spans_per_machine));
    }
    else
    {
        config_settings->set_automatic_span_control(false);
        // Setting the number of spans
        if (horizontal_spans_per_machine > 0)
        {
            config_settings->set_nb_spans(
                mi::Uint32(horizontal_spans_per_machine) * irc.m_icluster_configuration->get_number_of_hosts());
        }
        else
        {
            config_settings->set_nb_spans(horizontal_spans); // at least set to "1" (see above)
        }
    }

    // Size of the result queue
    mi::Sint32 result_queue = get_sint32(prj->get("index::result_queue"));
    if (result_queue < 1)
    {
        WARN_LOG<<"Result queue ["<<result_queue<<"] should not be negative, using default value of 1";
        result_queue = 1;
    }
    config_settings->set_size_of_rendering_results_in_queue(mi::Uint32(result_queue));

    const bool enable_fb_blending = get_bool(prj->get("index::enable_frame_buffer_blending", "no"));
    config_settings->enable_frame_buffer_blending(enable_fb_blending);

    // Additional CUDA debug checks
    config_settings->set_cuda_debug_checks_enabled(get_bool(prj->get("index::cuda_debug_checks")));

    mi::Uint32 nb_timesteps = get_uint32(prj->get("index::timesteps"));
    config_settings->set_nb_timesteps(nb_timesteps);

    // Dynamic memory management
    config_settings->set_dynamic_memory_management_enabled(get_bool(prj->get("index::dynamic_memory_management")));

    // Automatic workload balancing
    config_settings->set_workload_balancing(get_bool(prj->get("index::workload_balancing")));

    // Allow access to CUDA memory for user-supplied compute tasks if requested
    config_settings->set_cuda_compute_surface_access_enabled(
        get_bool(prj->get("index::cuda_compute_surface_access", "no")));

    //
    // Initial rendering configuration: These settings must be applied before rendering
    // starts and can't be changed afterwards.
    //

    // Subcube size
    mi::Uint32_3 subcube_size = get_vec_uint32_3(prj->get("index::subcube_size", "510 510 510"));

    // Subcube border size
    mi::Uint32 subcube_border_size = get_uint32(prj->get("index::subcube_border_size", "1"));

    // Supporting continuous translation requires additional memory compared to when only
    // integer translation is supported
    bool support_continuous_translation = get_bool(prj->get("index::support_continuous_translation", "no"));

    // Subcube size will be changed internally if rotation of volume scene elements is supported
    bool support_rotation = get_bool(prj->get("index::support_rotation", "no"));

    // Same if they are scaled with a scale factor < 1, need to specify the minimum factor then.
    mi::Float32_3 min_scale = get_vec_float32_3(prj->get("index::minimal_scale", "1 1 1"));

    if (!config_settings->set_subcube_configuration(
            subcube_size, subcube_border_size, support_continuous_translation, support_rotation, min_scale))
    {
        ERROR_LOG << "Failed to configure the subcubes.";
    }

    // Sparse volumes
    const mi::Uint32_3 svol_brick_size    = get_vec_uint32_3(prj->get("index::sparse_volume_brick_size", "62 62 62"));
    const mi::Uint32   svol_brick_overlap = get_uint32(prj->get("index::sparse_volume_brick_shared_border_size", "1"));

    nv::index::IConfig_settings::Sparse_volume_config svol_cfg;
    svol_cfg.brick_dimensions         = svol_brick_size;
    svol_cfg.brick_shared_border_size = svol_brick_overlap;

    if (!config_settings->set_sparse_volume_configuration(svol_cfg)) 
    {
        ERROR_LOG << "Failed to set sparse volume configuration.";
    }

    // Texture filter mode used for volume rendering
    {
        using nv::index::IConfig_settings;

        const std::string filter_name = prj->get("index::volume_filter", "nearest");
        IConfig_settings::Volume_filtering_modes filter_mode = IConfig_settings::VOLUME_FILTER_NEAREST;
        if (filter_name == "nearest")
            filter_mode = IConfig_settings::VOLUME_FILTER_NEAREST;
        else if (filter_name == "trilinear_post_hw")
            filter_mode = IConfig_settings::VOLUME_FILTER_TRILINEAR_POST_HW;
        else if (filter_name == "trilinear_post_sw")
            filter_mode = IConfig_settings::VOLUME_FILTER_TRILINEAR_POST_SW;
        else if (filter_name == "trilinear_pre_sw")
            filter_mode = IConfig_settings::VOLUME_FILTER_TRILINEAR_PRE_SW;
        else if (filter_name == "tricubic_catmull")
            filter_mode = IConfig_settings::VOLUME_FILTER_TRICUBIC_CATMULL;
        else if (filter_name == "tricubic_catmull_post_hw")
            filter_mode = IConfig_settings::VOLUME_FILTER_TRICUBIC_CATMULL_POST_HW;
        else if (filter_name == "tricubic_bspline")
            filter_mode = IConfig_settings::VOLUME_FILTER_TRICUBIC_BSPLINE;
        else if (filter_name == "tricubic_bspline_post_hw")
            filter_mode = IConfig_settings::VOLUME_FILTER_TRICUBIC_BSPLINE_POST_HW;
        else
            ERROR_LOG << "Invalid volume filter mode specified for index::volume_filter: " << filter_name;

        config_settings->set_volume_filter(filter_mode);
    }

    // Threshold for volume picking
    config_settings->set_volume_picking_threshold(get_float32(prj->get("index::volume_picking_threshold", "-1.0")));;

    // Ray segment accumulation
    config_settings->set_ray_segment_accumulation_technique(get_bool(prj->get("index::ray_segment_accum", "no")));
    config_settings->set_rgba_volume_opacity_mapping(get_bool(prj->get("index::rgba_opacity_mapping", "no")));
}

//----------------------------------------------------------------------
void Index_application::render_frame_update_camera(
    Nvindex_rendering_context& irc)
{
    static mi::neuraylib::Tag last_camera_tag = mi::neuraylib::NULL_TAG;

    Multiple_stereo_camera* ms_camera = m_appdata->get_user_interaction(0)->get_multiple_stereo_camera();
    const mi::neuraylib::Tag current_camera_tag =
        ms_camera->get_current_camera_tag_considering_ortho();
    
    if (last_camera_tag != current_camera_tag)
    {
        // Camera changed: Edit the scene using a new transaction, but only after acquiring the
        // lock, to prevent race conditions.
        mi::base::Lock::Block block(&m_appdata->m_scene_edit_lock);
        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            irc.get_current_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
        assert(dice_transaction.is_valid_interface());
        {
            mi::base::Handle<const nv::index::ISession> session(
                dice_transaction->access<nv::index::ISession>(irc.m_session_tag));

            // Edit scene description to set the camera used
            {
                mi::base::Handle<nv::index::IScene> scene(
                    dice_transaction->edit<nv::index::IScene>(session->get_scene()));
                assert(scene.is_valid_interface());

                scene->set_camera(current_camera_tag);
                last_camera_tag = current_camera_tag;
            }
        }
        dice_transaction->commit();
    }
}

//----------------------------------------------------------------------
void Index_application::render_frame_stereo_view(
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    assert(dice_transaction != 0);
    if (m_appdata->is_stereo_mode())
    {
        Multiple_stereo_camera* ms_camera = 
            m_appdata->get_user_interaction(0)->get_multiple_stereo_camera();
        ms_camera->update_swap_main_camera_view(dice_transaction);
    }
}

//----------------------------------------------------------------------
bool Index_application::render_frame_check_cluster_state_is_quit(
    Nvindex_rendering_context& irc)
{
    // Terminate if we lost some hosts, to prevent the viewer node from taking over all the work
    // from the now missing hosts and potentially running out of memory.
    mi::Uint32 new_nb_hosts = irc.m_icluster_configuration->get_number_of_hosts();
    if (new_nb_hosts < m_previous_nb_hosts && m_previous_nb_hosts != 0)
    {
        INFO_LOG << "At least one host has left - had "
                 << m_previous_nb_hosts << " hosts, only " << new_nb_hosts << " left";
        const bool fail_safety_enabled = get_bool(m_prj->get("app::enable_fail_safety", "no"));
        if (!fail_safety_enabled)
        {
            INFO_LOG << "Leaving because at least one host has left - had "
                     << m_previous_nb_hosts << " hosts, only " << new_nb_hosts << " left";
            m_appdata->set_app_run(false);
            return false;
        }
    }
    m_previous_nb_hosts = new_nb_hosts;

    return true;
}

//----------------------------------------------------------------------
void Index_application::render_frame_compute(
    Nvindex_rendering_context&        irc,
    mi::Uint32                        frame_num,
    bool                              need_compute,
    mi::Float32&                      application_compute_time,
    mi::Uint64&                       application_compute_memory,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
#ifdef NVINDEX_HAS_COMPUTE
    if (m_compute_wrapper.is_active() && frame_num > 10) // skip the first few frames to ensure data is there
    {
        std::vector<mi::neuraylib::Tag> volume_tag_vec = m_appdata->get_volume_tag_vec();
        if (volume_tag_vec.begin() != volume_tag_vec.end())
        {
            const mi::neuraylib::Tag volume_tag = *(volume_tag_vec.begin());
            const mi::Float64 compute_start = nv::index_common::get_time();
            m_compute_wrapper.step(volume_tag, irc.m_session_tag, dice_transaction);
            const mi::Float64 compute_end   = nv::index_common::get_time();
            application_compute_time = static_cast<mi::Float32>(compute_end - compute_start);
            application_compute_memory = m_compute_wrapper.get_allocated_device_memory();
        }
    }
#endif

    const bool GTC_DEMO_2017 = true;
    if (GTC_DEMO_2017 && frame_num > 1)
    {
        unstructured_volume_compute(frame_num, "ivol_data_01", irc.m_session_tag, dice_transaction);
    }

#ifdef NVINDEX_HAS_MPI_IPC_COMPUTE
    if (frame_num > 1 && need_compute)
    {
        std::vector<mi::neuraylib::Tag> volume_tag_vec = m_appdata->get_volume_tag_vec();
        if (volume_tag_vec.begin() != volume_tag_vec.end())
        {
            const mi::neuraylib::Tag volume_tag = *(volume_tag_vec.begin());
            edit_volume_data(volume_tag, irc.m_session_tag, dice_transaction);
        }
    }
#endif // NVINDEX_HAS_MPI_IPC_COMPUTE
}

//----------------------------------------------------------------------
void Index_application::render_frame_handle_animation(
    Nvindex_rendering_context&                          irc,
    mi::Uint32                                          frame_num,
    mi::base::Handle<mi::neuraylib::IDice_transaction>& dice_transaction)
{
    mi::Sint32 animation_step = frame_num;

    mi::math::Vector<mi::Sint32, 2> range(
        get_vec_sint32_2(m_prj->get("app::animation::range", "0 0")));
    if (range.y > range.x && range.x >= 0)
    {
        // Map requested animation range to frame number
        animation_step = range.x + frame_num;
        if (animation_step >= range.y)
        {
            // Quit after last animation step
            m_appdata->set_app_run(false);
        }

        // Number used for writing snapshots (when 'is_sequence_mode' is enabled)
        m_appdata->set_snapshot_index(animation_step);

        INFO_LOG << "Animation step " << animation_step
                 << " (" << frame_num + 1 << "/" << (range.y - range.x) + 1 << ")";
    }
    else
    {
        // No range given, just use frame number as animation step and continue until no .prj file
        // is found
        INFO_LOG << "Animation step " << animation_step;
    }

    const std::string filename_format = m_prj->get("app::animation::project_files", "anim_%03d.prj");
    const int BUFSIZE = 1024;
    char buf[BUFSIZE];
    snprintf(buf, BUFSIZE, filename_format.c_str(), animation_step);
    const std::string filename(buf);

    INFO_LOG << "Loading animation project file '" << filename << "'";
    String_dict extra_dict;
    if (load_application_project_file(filename, extra_dict))
    {
        m_prj->insert_all(extra_dict);

        if (get_bool(m_prj->get("app::animation::update_scene", "1")))
        {
            // Update scene graph
            add_scene_graph(dice_transaction, irc.m_session_tag, *m_prj);
        }

        // Update camera
        setup_camera(irc, dice_transaction.get());

        // Trigger snapshot
        m_appdata->set_snapshot_on(true);
    }
    else
    {
        ERROR_LOG << "Could not open animation project file '" << filename << "', aborting.";
        m_appdata->set_app_run(false);
    }
}

//----------------------------------------------------------------------
void Index_application::render_frame_handle_export_session(
    Nvindex_rendering_context&        irc,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    // Export the scene if requested via the command line ("-export file" or "-export stdout")
    static bool export_done = false; // Only export once
    if (!export_done && m_prj->is_defined("export"))
    {
        export_session(irc, dice_transaction, m_prj->get("export"));
    }
    export_done = true;
}

//----------------------------------------------------------------------
bool Index_application::render_frame_quit_handle_case_0()
{
    // Option to quit before rendering the first image ("-quit 0")
    if (m_appdata->get_quit_after_frame_number() == 0)
    {
        INFO_LOG << "Shutting down before rendering the first frame, as requested.";
        m_appdata->set_app_run(false);
        return false;
    }

    return true;
}

//----------------------------------------------------------------------
bool Index_application::render_frame_quit_handle_case_n(
    mi::Uint32 frame_num)
{
    // Option to quit after rendering a given number of frames ("-quit <number>")
    if ((m_appdata->get_quit_after_frame_number() > 0) && 
        (frame_num + 1) >= static_cast<mi::Uint32>(m_appdata->get_quit_after_frame_number()))
    {
        INFO_LOG << "Shutting down after rendering " << (frame_num + 1) << " frame" << (frame_num > 0 ? "s" : "")
                 << ", as requested.";
        m_appdata->set_app_run(false);
        return false;
    }

    return true;
}

//----------------------------------------------------------------------
void Index_application::render_frame_common_setup_before_db_access(
    Nvindex_rendering_context& irc,        
    mi::Uint32                 frame_num)
{
    irc.m_frame_num = frame_num;

    /// Camera update 
    render_frame_update_camera(irc);

    // Example how to find out which database elements were modified since the last frame
    example_how_to_check_modification_of_db_element(irc);

    // traverse scene to get the scene status information
    // Scene_traverse_print_element print_scene_element_strategy;
    Scene_traverse_volume_heightfield_element volume_heightfiled_status;
    {
        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            irc.get_current_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
        assert(dice_transaction.is_valid_interface());
        {
            get_scene_status(&volume_heightfiled_status, irc.m_session_tag, dice_transaction.get());
        }
        dice_transaction->commit();
    }

    m_appdata->set_volume_tag_vec(volume_heightfiled_status.get_volume_tag_vec());
    m_appdata->set_sparse_volume_tag_vec(volume_heightfiled_status.get_sparse_volume_tag_vec());
    m_appdata->set_heightfield_tag_vec(volume_heightfiled_status.get_heightfield_tag_vec());
    m_appdata->set_rtc_param_buffer_tag_vec(volume_heightfiled_status.get_rtc_param_buffer_tag_vec());

    // update the self/remote host status
    update_host_map(irc, *(m_appdata->peek_host_info()));
}

//----------------------------------------------------------------------
bool Index_application::render_frame_common_setup_with_db_access(    
    Nvindex_rendering_context&        irc,        
    mi::Uint32                        frame_num,    
    bool                              need_compute,
    mi::Float32&                      application_compute_time,
    mi::Uint64&                       application_compute_memory,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    // Experimental: camera animator
    Multiple_stereo_camera* ms_camera = m_appdata->get_user_interaction(0)->get_multiple_stereo_camera();
    update_camera_by_animator(ms_camera->get_current_camera_tag_considering_ortho(),
                              m_appdata->peek_camera_animator_data(),
                              dice_transaction);

    // Experimental: stereo view    
    render_frame_stereo_view(dice_transaction);

    // Check the cluster host state
    if (!render_frame_check_cluster_state_is_quit(irc))
    {
        return false;                 // quit if false ... some host left the cluster
    }

    // Clear previous progress and frame info
    irc.m_progress_callback->set_progress(0.f);

    // Set a progress callback passed to the rendering call.
    irc.m_frame_info_callbacks = new Frame_info_callbacks();
    assert(irc.m_frame_info_callbacks.is_valid_interface());

    // compute related task
    render_frame_compute(irc, frame_num, need_compute, 
                         application_compute_time, application_compute_memory,
                         dice_transaction);

    return true;
}

//----------------------------------------------------------------------
bool Index_application::render_frame_common_postprocess(    
    Nvindex_rendering_context&                          irc,
    mi::Uint32                                          frame_num,
    mi::base::Handle<mi::neuraylib::IDice_transaction>& dice_transaction)
{
    // INFO_LOG << "DEBUG: render_frame_common_postprocess";    
    
    // For one of the video streaming enabled, copy the rendering buffer 
    if (m_appdata->is_any_video_stream_enabled())
    {
        // INFO_LOG << "DEBUG: render_frame_common_postprocess: video stream enabled.";    
        mi::base::Handle<Canvas> canvas(irc.m_canvas_buffers[0]->get_render_canvas());
        copy_result_pixel_to_canvas(irc.m_span_buffer.get(), canvas.get());

        // Tripple buffering for video-rendering swapping to accelerate the video encoding.
        // - Trigger internal buffer swapping
        irc.m_canvas_buffers[0]->rendering_finished();
    }

    // For rtmp video streaming 
    // Nothing to do for RTMP video streaming here

    // For dice bridge video streaming 
    if (irc.m_bridge_video_source.is_valid_interface())
    {
        Bridge_video_stream *video_source = (Bridge_video_stream *)irc.m_bridge_video_source.get();
        video_source->frame_ready(frame_num);
    }

    // For html5 video streaming 
    // if (irc.m_html5_video_stream_manager != 0)
    // Nothing to do for html5 video streaming here

    // Playback of recorded user commands
    if (m_prj->is_defined("playback"))
    {
        handle_playback(m_prj->get("playback"), irc);
    }

    // Option quit handling
    if (!render_frame_quit_handle_case_n(frame_num))
    {
        return false;          // quit if false ... the user requested
    }

    return true;
}


//----------------------------------------------------------------------
void Index_application::render_frame_render_performance_statistics_single_view(
    Nvindex_rendering_context&        irc,
    mi::Uint32                        frame_num,
    mi::Float32                       application_compute_time,
    mi::Uint64                        application_compute_memory,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    // Read the snapshot trigger here and not read m_appdata later in
    // this method to avoid to use the snapshot trigger while
    // rendering. Otherwise you need a lock since rtmp handler thread
    // can turn on this.
    const bool is_snapshot_on_this_frame = m_appdata->is_snapshot_on();

    // index canvas setup
    render_frame_common_canvas_setup(irc, is_snapshot_on_this_frame);

    // render single view opengl integration setup
    render_frame_single_view_opengl_integration_setup(irc, dice_transaction);

    // enable opengl copy buffer if possible
    render_frame_common_opengl_canvas_setup(irc);

    // Main render call (rendering and storing the performance values)
    mi::base::Handle<nv::index::IFrame_results> frame_results(
        render_frame_single_view_render_call(irc, is_snapshot_on_this_frame, dice_transaction));

    const mi::base::Handle<nv::index::IError_set> err_set(frame_results->get_error_set());
    if (err_set->any_errors())
    {
        // shutdown in case of errors
        INFO_LOG << "Shutting down after encountering errors during rendering.";
        m_appdata->set_app_run(false);
    }

    {
        mi::base::Lock::Block block(&m_appdata->m_performance_values_lock);
        m_appdata->m_performance_values = mi::base::make_handle(frame_results->get_performance_values());
        m_appdata->m_application_compute_time   = application_compute_time;
        m_appdata->m_application_compute_memory = application_compute_memory;
    }

    // show statistics if needed
    show_statistics_overlay(m_appdata->m_show_statistics_overlay, irc.m_icluster_configuration.get(),
                            m_appdata->m_performance_values.get(), m_appdata->m_statistics_overlay_large);

    // Update recent n frames statistics
    update_recent_n_statistics(m_appdata->m_performance_values.get());

    // experimental timestep
    render_frame_timestep_handling(irc, dice_transaction);
        
    // log performance value depends on the performance logging state
    const mi::base::Handle<nv::index::IPerformance_values> perf_value(frame_results->get_performance_values());
    log_performance_information(perf_value.get(), frame_num);

    // performance measurement frames
    if (m_appdata->get_remaining_performance_measuring_frames() > 0)
    {
        m_appdata->decrement_remaining_performance_measuring_frames();
        if (m_appdata->get_remaining_performance_measuring_frames() == 0)
        {
            turn_off_logging();            // turn off the logging when reach to 0
            INFO_LOG << "Turn off the performance logging.";
        }
    }

    // Display nice axis in the lower part of the screen
    if (m_appdata->is_show_coordinate_system_overlay())
    {
        Multiple_stereo_camera* ms_camera =
            m_appdata->get_user_interaction(0)->get_multiple_stereo_camera();
        gl_draw_axis(ms_camera->get_current_main_camera_tag(), dice_transaction);
    }

    // Snapshot
    if (is_snapshot_on_this_frame)
    {
        save_span_buffer_snapshot(static_cast<Span_renderer_IF*>(irc.m_span_buffer.get()),
                                  irc.m_frame_num);
        // But turn off can be done through m_appdata since this is
        // the only single place to do that.
        m_appdata->set_snapshot_on(false); 
    }
}

//----------------------------------------------------------------------
void Index_application::render_frame_timestep_handling(
    Nvindex_rendering_context&        irc,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    static mi::Float64 previous_time = nv::index_common::get_time();
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<nv::index::ISession>(irc.m_session_tag));
    mi::base::Handle<nv::index::IConfig_settings> configs(
        dice_transaction->edit<nv::index::IConfig_settings>(session->get_config()));
    mi::Uint32 current_step = configs->get_current_timestep();
    if (irc.m_is_tsteps_playforward)
    {
        mi::Float64 current_time = nv::index_common::get_time();
        if (current_time>=previous_time+0.04f)
        {
            current_step++;
            if (current_step >= configs->get_nb_timesteps())
            {
                current_step = 0;
            }
            configs->set_current_timestep(current_step);
            previous_time = current_time;
        }
    }
    else if (irc.m_is_tsteps_playbackward)
    {
        mi::Float64 current_time = nv::index_common::get_time();
        if (current_time>=previous_time+0.04f)
        {
            if (current_step <= 0)
            {
                current_step = configs->get_nb_timesteps();
            }
            current_step--;
            configs->set_current_timestep(current_step);
            previous_time = current_time;
        }
    }
    else if (irc.m_is_pause_tsteps)
    {
        configs->set_current_timestep(current_step);
        irc.m_is_pause_tsteps = false;
    }
    else if (irc.m_is_tsteps_stepback)
    {
        if (current_step <= 0)
        {
            current_step = configs->get_nb_timesteps();
        }
        current_step--;
        configs->set_current_timestep(current_step);
        irc.m_is_tsteps_stepback = false;
    }
    else if (irc.m_is_tsteps_stepforward)
    {
        current_step++;
        if (current_step >= configs->get_nb_timesteps())
        {
            current_step = 0;
        }
        configs->set_current_timestep(current_step);
        irc.m_is_tsteps_stepforward = false;
    }
}

//----------------------------------------------------------------------
void Index_application::render_frame_single_view(
    Nvindex_rendering_context& irc,
    mi::Uint32                 frame_num,
    bool                       need_compute)
{
    // start transaction
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        irc.get_current_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    {
        // rendering setup with db access
        mi::Float32 application_compute_time   = 0.0f;
        mi::Uint64  application_compute_memory = 0;
        {
            if (m_appdata->is_animation_mode())
            {
                // This will load a project file snippet for the current frame and set snapshot mode
                render_frame_handle_animation(irc, frame_num, dice_transaction);
            }

            const bool is_success = 
                render_frame_common_setup_with_db_access(irc, frame_num, need_compute, 
                                                         application_compute_time, application_compute_memory,
                                                         dice_transaction.get());
            if (!is_success)
            {
                dice_transaction->abort();
                return;
            }

            // Update the session, this must be called before rendering, using the same transaction
            {
                mi::base::Lock::Block block(&(m_appdata->m_scene_update_lock));
                irc.m_iindex_session->update(irc.m_session_tag, dice_transaction.get());
            }

            // Export the scene if requested via the command line ("-export file" or "-export stdout")
            render_frame_handle_export_session(irc, dice_transaction.get());

            // Check "-quit 0" option
            if (!render_frame_quit_handle_case_0())
            {
                dice_transaction->commit();
                return;
            }
        }

        // render call
        // Single view: rendering, performance update, statistics, snapshot
        render_frame_render_performance_statistics_single_view(
            irc, frame_num,
            application_compute_time, application_compute_memory,
            dice_transaction.get());

        // rendering post process
        {
            const bool is_success = render_frame_common_postprocess(irc,        
                                                                    frame_num,
                                                                    dice_transaction);
            if (!is_success)
            {
                dice_transaction->abort();
                return;
            }
        }
    }
    // Finally commit the transaction
    dice_transaction->commit();
}

//----------------------------------------------------------------------
void Index_application::render_frame_common_canvas_setup(
    Nvindex_rendering_context& irc,
    bool is_snapshot_on_this_frame)
{
    if (is_snapshot_on_this_frame)
    {
        const mi::Sint32_2 image_res =
            get_vec_sint32_2(m_appdata->peek_app_proj()->get("index::image_file_canvas_resolution"));
        request_canvas_resolution(image_res);
    }
    update_canvas_resolution(irc);

    // The following call resizes the buffer (CPU version) or clears the GL canvas (OpenGL version)
    assert(irc.m_span_buffer.is_valid_interface());
    irc.m_span_buffer->set_background_color(m_appdata->m_background_color);
    const mi::math::Vector<mi::Sint32, 2> main_window_resolution =
        m_appdata->get_user_interaction(0)->
        get_examiner()->get_main_window_resolution();
    irc.m_span_buffer->set_buffer_resolution(main_window_resolution);
}

//----------------------------------------------------------------------
void Index_application::render_frame_common_opengl_canvas_setup(
    Nvindex_rendering_context& irc)
{
    // No 3D rendering from here on anymore. Copy the tile if OpenGL is available.
    gl_push_matrix_state();

    const mi::math::Vector<mi::Sint32, 2> buffer_resolution = irc.m_span_buffer.get()->get_buffer_resolution();
    gl_setup_2d_rendering(static_cast< mi::Float32 >(buffer_resolution.x),
                          static_cast< mi::Float32 >(buffer_resolution.y));
}

//----------------------------------------------------------------------
void Index_application::render_frame_common_restore_canvas_size(
    bool is_snapshot_on_this_frame)
{
    // Restore snap shot canvas size 
    if (is_snapshot_on_this_frame)
    {
        // Reset to the original resolution after a snapshot
        const mi::Sint32_2 original_canvas_res =
            get_vec_sint32_2(m_appdata->peek_app_proj()->get("index::canvas_resolution"));
        request_canvas_resolution(original_canvas_res);
    }
}

//----------------------------------------------------------------------
void Index_application::render_frame_single_view_opengl_integration_setup(
    Nvindex_rendering_context&        irc,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    // OpenGL integration 
    Multiple_stereo_camera* ms_camera = 
        m_appdata->get_user_interaction(0)->get_multiple_stereo_camera();
    const mi::neuraylib::Tag current_camera_tag =
        ms_camera->get_current_camera_tag_considering_ortho();

    mi::base::Handle<const nv::index::ICamera > cam(
        dice_transaction->access< nv::index::ICamera >(current_camera_tag));
    assert(cam.is_valid_interface());

    // Initialize matrices
    gl_setup_projection_matrix(cam.get());

    // Access session and configuration for later use
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<const nv::index::ISession>(irc.m_session_tag));

    // Setting the modelview transformation
    mi::math::Vector<mi::Float32, 3> eye, target ,up;
    Camera_tool::get_glu_lookat_vector(cam.get(), eye, target ,up);
    gl_setup_lookat(eye, target, up);

    mi::base::Handle<const nv::index::IScene> scene(
        dice_transaction->access<const nv::index::IScene>(session->get_scene()));
    assert(scene.is_valid_interface());
    mi::math::Matrix<mi::Float32, 4, 4> const mat = scene->get_transform_matrix();
    gl_set_world_transform(mat);

    // render application primitive and get the z-buffer.
    if (m_appdata->is_enable_opengl_z_buffer_for_rendering())
    {
        const mi::math::Bbox<mi::Float32, 3> xyz_roi =
            get_XYZ_global_region_of_interest_bbox(session.get(), dice_transaction);

        // Create well geometry, unless it is already there
        gl_application_example_well_creation(xyz_roi, *m_prj);

        // Well rendering
        gl_application_draw_object(
            Nvindex_AppData::instance()->get_opengl_appdata()->get_opengl_application_buffer_ptr());
    }

    // render your primitives in world space.
    const bool is_show_coordinate_system_axes = get_bool(m_prj->get("app::is_show_coordinate_system_axes", "0"));
    gl_show_icamera_coordinate_system(is_show_coordinate_system_axes, scene->get_transform_matrix());

    // Visualize additional information for such as bounding box of the entire volume using OpenGL (if available)
    const bool show_opengl_guides = get_bool(m_prj->get("app::show_opengl_guides", "0"));
    if (show_opengl_guides)
    {
        const mi::math::Bbox<mi::Float32, 3> roi = scene->get_clipped_bounding_box();
        visualize_region_of_interest(
            m_appdata->is_show_region_of_interest_overlay(),
            m_appdata->is_show_bounding_boxes_overlay(),
            irc.m_session_tag,
            roi,
            dice_transaction);
    }

    // Render slice wireframe
    // render_frame_wireframe_slice(dice_transaction); // DISABLED


    // ------------------------------------------------------------------------------------------------
    // No 3D rendering from here on anymore ...
}


//----------------------------------------------------------------------
mi::base::Handle<nv::index::IFrame_results> Index_application::render_frame_single_view_render_call(
    Nvindex_rendering_context&        irc,
    bool                              is_snapshot_on_this_frame,
    mi::neuraylib::IDice_transaction* dice_transaction)
{

    // Render the subsurface data and display the horizontal spans using the user-defined span buffer renderer
    irc.m_frame_info_callbacks->clear_dynamic_allocation_events();

    Opengl_application_buffer* opengl_app_buffer = render_frame_get_opengl_app_buffer();
    assert(irc.m_span_buffer.is_valid_interface());

    // Main NVIDIA IndeX rendering call
    mi::base::Handle<nv::index::IFrame_results> frame_results(
        irc.m_iindex_rendering->render(
            irc.m_session_tag,
            irc.m_span_buffer.get(),
            dice_transaction,
            irc.m_progress_callback.get(),
            irc.m_frame_info_callbacks.get(),
            m_appdata->m_immediate_final_parallel_compositing,
            opengl_app_buffer));

    const mi::base::Handle<nv::index::IError_set> err_set(frame_results->get_error_set());
    if (err_set->any_errors())
    {
        std::ostringstream os;
        const mi::Uint32 nb_err = err_set->get_nb_errors();
        for (mi::Uint32 e = 0; e < nb_err; ++e)
        {
            if (e != 0) os << '\n';
            os << err_set->get_error(e)->get_error_string();
        }

        ERROR_LOG << "IIndex_rendering rendering call failed with the following error(s): " << '\n'
                  << os.str();

        return frame_results;
    }

    // Show frame info callback's data
    render_frame_show_frame_info(irc.m_frame_info_callbacks.get());

    // restore the canvas size if snapshot
    render_frame_common_restore_canvas_size(is_snapshot_on_this_frame);

    // visualize performance value with 2D OpenGL
    render_frame_single_view_visualize_perf_with_2d_opengl(irc, frame_results.get(), dice_transaction);

    // Additional gl rendering: debug, visualize heightfield delete polygon
    render_frame_single_view_additional_gl_rendering(irc, dice_transaction);

    // check GL error of this rendering
    render_frame_common_check_gl_error();

    return frame_results;
}

//----------------------------------------------------------------------
void Index_application::render_frame_single_view_additional_gl_rendering(
    Nvindex_rendering_context&        irc,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    // ------------------------------------------------------------------------------------------------    
    // Initialize matrices for 3D rendering again ... for debugging
    gl_push_matrix_state();

    Multiple_stereo_camera* ms_camera = 
        m_appdata->get_user_interaction(0)->get_multiple_stereo_camera();
    const mi::neuraylib::Tag current_camera_tag =
        ms_camera->get_current_camera_tag_considering_ortho();
    mi::base::Handle<const nv::index::ICamera > cam(
        dice_transaction->access< nv::index::ICamera >(current_camera_tag));
    assert(cam.is_valid_interface());

    gl_setup_projection_matrix(cam.get());

    if (    Heightfield_workflow_functionality::instance()->is_heightfield_workflow_delete_operation()
           || Heightfield_workflow_functionality::instance()->is_heightfield_workflow_elevation_change_operation()
           || Heightfield_workflow_functionality::instance()->is_heightfield_workflow_gridding_operation() )
    {
        // Render the heightfield delete polygon
        mi::base::Handle<const nv::index::IRegular_heightfield> heightfield_scene_element(
            dice_transaction->access<const nv::index::IRegular_heightfield>(
                Heightfield_workflow_functionality::instance()->get_heightfield()));
        const mi::math::Matrix<mi::Float32, 4, 4> mat = heightfield_scene_element->get_transform();
        gl_draw_heightfield_delete_polygon(mat);
    }

    gl_pop_matrix_state();
}

//---------------------------------------------------------------------- 
void Index_application::render_frame_common_check_gl_error()
{
    // Check for OpenGL errors
    std::string msg = gl_error();
    if (!msg.empty())
    {
        WARN_LOG << "glGetError(): " << msg << " (" << __FILE__ << ":" << __LINE__ << ")";
    }
}

//----------------------------------------------------------------------
//======================================================================
//----------------------------------------------------------------------
void Index_application::render_frame_wireframe_slice(
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    std::vector<mi::neuraylib::Tag> volume_tag_vec = m_appdata->get_volume_tag_vec();
    if (!volume_tag_vec.empty())
    {
        const mi::neuraylib::Tag volume_tag = volume_tag_vec.at(0);
        assert(volume_tag.is_valid());
        mi::base::Handle<const nv::index::IRegular_volume> volume(
            dice_transaction->access<const nv::index::IRegular_volume>(volume_tag));

        const mi::Uint32 nb_slices = volume->get_nb_slices();
        const mi::math::Matrix<mi::Float32, 4, 4> mat = volume->get_transform();
        mi::math::Bbox<mi::Float32, 3> roi = volume->get_IJK_region_of_interest();
        roi.max -= mi::math::Vector<mi::Float32, 3>(1.f); // Adapt for half voxels on boundary

        for (mi::Uint32 i=0; i<nb_slices; ++i)
        {
            const mi::neuraylib::Tag slice_tag = volume->get_slice(i);
            mi::base::Handle<const nv::index::ISlice_scene_element> slice(
                dice_transaction->access<const nv::index::ISlice_scene_element>(slice_tag));

            if (!slice->get_enabled())
                continue;

            mi::base::Handle<const nv::index::ISection_scene_element> section_slice(
                slice->get_interface<const nv::index::ISection_scene_element>());
            if (section_slice.get())
            {
                const nv::index::ISection_scene_element::Slice_orientation orientation =
                    section_slice->get_orientation();
                if (orientation==nv::index::ISection_scene_element::INLINE_SECTION)
                {
                    const mi::math::Vector<mi::Float32, 3> llf(section_slice->get_position(), roi.min.y, roi.min.z);
                    const mi::math::Vector<mi::Float32, 3> urb(section_slice->get_position(), roi.max.y, roi.max.z);
                    gl_draw_slice(llf, urb, mat);
                }
                else if (orientation==nv::index::ISection_scene_element::CROSS_LINE_SECTION)
                {
                    const mi::math::Vector<mi::Float32, 3> llf(roi.min.x, section_slice->get_position(), roi.min.z);
                    const mi::math::Vector<mi::Float32, 3> urb(roi.max.x, section_slice->get_position(), roi.max.z);

                    gl_draw_slice(llf, urb, mat);
                }
                else if (orientation==nv::index::ISection_scene_element::HORIZONTAL_SECTION)
                {
                    const mi::math::Vector<mi::Float32, 3> llf(roi.min.x, roi.min.y, section_slice->get_position());
                    const mi::math::Vector<mi::Float32, 3> urb(roi.max.x, roi.max.y, section_slice->get_position());
                    gl_draw_slice(llf, urb, mat);
                }
            }
            else
            {
                mi::base::Handle<const nv::index::IVertical_profile_scene_element> profile(
                    slice->get_interface<const nv::index::IVertical_profile_scene_element>());
                const mi::Uint32 nb_vertices = profile->get_nb_vertices();

                if (nb_vertices > 1)
                {
                    mi::math::Vector_struct<mi::Float32, 2> v0;
                    mi::math::Vector_struct<mi::Float32, 2> v1;
                    for (mi::Uint32 i=0; i < nb_vertices-1; ++i)
                    {
                        profile->get_vertex(i, v0); profile->get_vertex(i+1, v1);
                        gl_draw_profile(v0, v1, roi, mat);
                    }
                }
            }
        }
    }
}

//----------------------------------------------------------------------
Opengl_application_buffer* Index_application::render_frame_get_opengl_app_buffer() const
{
    Opengl_application_buffer* opengl_app_buffer = 0;

    if (m_appdata->is_enable_opengl_z_buffer_for_rendering())
    {
        opengl_app_buffer = 
            Nvindex_AppData::instance()->get_opengl_appdata()->get_opengl_application_buffer_ptr();
        const mi::Sint32_2 requested_res = m_appdata->get_request_canvas_resolution();
        assert(opengl_app_buffer != 0);
        opengl_app_buffer->resize_buffer(requested_res);

        static bool warning_printed = false;
        if (!m_appdata->is_use_opengl() && !warning_printed)
        {
            WARN_LOG << "Experimental feature OpenGL z-buffer integration requires OpenGL support, "
                     << "which is not supported in this build. You may see unexpected results. "
                     << "Maybe you forgot the --gl option?";
            warning_printed = true;
        }
    }

    return opengl_app_buffer;
}

//----------------------------------------------------------------------
void Index_application::render_frame_show_frame_info(
    Frame_info_callbacks* frame_info_callbacks    
    ) const
{
    std::vector<Frame_info_callbacks::Dynamic_alloc_event> dyn_alloc_events =
        frame_info_callbacks->get_dynamic_allocation_events();
    if (!dyn_alloc_events.empty())
    {
        INFO_LOG << "Dynamic allocation events occurred during the IIndex_rendering rendering call:";
        for (mi::Size i = 0; i < dyn_alloc_events.size(); ++i)
        {
            INFO_LOG << " * event " << i << ":\n"
                     << "   - host id:           " << dyn_alloc_events[i].m_host_id << "\n"
                     << "   - device id:         " << dyn_alloc_events[i].m_device_id << "\n"
                     << "   - requested size:    " << dyn_alloc_events[i].m_memory_allocation_size << "\n"
                     << "   - available memory:  " << dyn_alloc_events[i].m_memory_available << "\n"
                     << "   - memory freed up:   " << dyn_alloc_events[i].m_memory_freed_up;
        }
    }

    std::vector<Frame_info_callbacks::Device_reset_event> dev_reset_events =
        frame_info_callbacks->get_device_reset_events();
    if (!dev_reset_events.empty())
    {
        INFO_LOG << "Device reset events occurred during the IIndex_rendering rendering call:";
        for (mi::Size i = 0; i < dev_reset_events.size(); ++i)
        {
            INFO_LOG << " * event " << i << ":\n"
                     << "   - host id:           " << dev_reset_events[i].m_host_id << "\n"
                     << "   - device id:         " << dev_reset_events[i].m_device_id;
        }
    }

    std::vector<mi::base::Handle<nv::index::IBalancing_operations> > balancing_ops =
        frame_info_callbacks->get_workload_balancing_operations();
    if (!balancing_ops.empty())
    {
        for (mi::Size i = 0; i < balancing_ops.size(); ++i)
        {
            const mi::Uint32 nb_operations = balancing_ops[i]->get_nb_operation();
            const char* event_desc = balancing_ops[i]->get_event_description();

            INFO_LOG << i << ".) Balancing event: '" << event_desc << "'" 
                     << " results in " << nb_operations << " balancing operations: ";
            for (mi::Uint32 j=0; j<nb_operations; ++j)
            {
                const nv::index::IBalancing_operation* op = balancing_ops[i]->get_operation(j);
                INFO_LOG << "Balance operation " << j << ": " << op->get_description();
            }
        }
    }
}

//----------------------------------------------------------------------
void Index_application::render_frame_single_view_visualize_perf_with_2d_opengl(
    Nvindex_rendering_context& irc,
    nv::index::IFrame_results* frame_results,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    if (frame_results != 0)
    {
        // visualize horizontal span, colormap, workload, performance
        // statistics (only for OpenGL version) using 2d OpenGL rendering
        const mi::Uint32 colormap_index = m_appdata->get_current_colormap_index();
        mi::neuraylib::Tag colormap_tag;
        if (colormap_index < get_number_of_colormap())
        {
            colormap_tag = get_colormap_tag(colormap_index);
        }

        mi::base::Handle<nv::index::IPerformance_values> perf_value(frame_results->get_performance_values());
        visualize_span_and_performance(irc, m_appdata, colormap_tag,
                                       perf_value.get(), dice_transaction);

    }

    // End of 2D OpenGL rendering
    gl_pop_matrix_state();
}

//----------------------------------------------------------------------
//======================================================================
//----------------------------------------------------------------------
void Index_application::render_frame_multi_view(
    Nvindex_rendering_context& irc,
    mi::Uint32                 frame_num, 
    bool                       need_compute)
{
    mi::Float32 application_compute_time   = 0.0f;
    mi::Uint64  application_compute_memory = 0;
    {
        // start transaction
        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            irc.get_current_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
        assert(dice_transaction.is_valid_interface());
        {
            const bool is_success = 
                render_frame_common_setup_with_db_access(irc, frame_num, need_compute, 
                                                         application_compute_time, application_compute_memory,
                                                         dice_transaction.get());
            if (!is_success)
            {
                dice_transaction->abort();
                return;
            }

            // Export the scene if requested via the command line ("-export file" or "-export stdout")
            render_frame_handle_export_session(irc, dice_transaction.get());

            // Check "-quit 0" option
            if (!render_frame_quit_handle_case_0())
            {
                dice_transaction->commit();
                return;
            }
        }
        dice_transaction->commit();
    }

    // render call
    // Multi-view: rendering, performance update, statistics, snapshot
    render_frame_render_performance_statistics_multi_view(
        irc, frame_num,
        application_compute_time, application_compute_memory);

    // rendering post process
    {
        // start transaction
        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            irc.get_current_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
        assert(dice_transaction.is_valid_interface());
        {
            const bool is_success = render_frame_common_postprocess(irc,        
                                                                    frame_num,
                                                                    dice_transaction);
            if (!is_success)
            {
                dice_transaction->abort();
                return;
            }
        }
        dice_transaction->commit();
    }
}

//----------------------------------------------------------------------
void Index_application::render_frame_render_performance_statistics_multi_view(
    Nvindex_rendering_context&        irc,
    mi::Uint32                        frame_num,
    mi::Float32                       application_compute_time,
    mi::Uint64                        application_compute_memory)
{
    // Read the snapshot trigger here and not read m_appdata later in
    // this method to avoid to use the snapshot trigger while
    // rendering. Otherwise you need a lock since rtmp handler thread
    // can turn on this.
    const bool is_snapshot_on_this_frame = m_appdata->is_snapshot_on();

    // index canvas setup
    render_frame_common_canvas_setup(irc, is_snapshot_on_this_frame);

    // No opengl integration

    // enable opengl copy buffer if possible
    render_frame_common_opengl_canvas_setup(irc);

    // Main render call (rendering and storing the performance values)
    const bool is_success = render_frame_multi_view_render_call(irc, is_snapshot_on_this_frame);

    if (!is_success)
    {
        // shutdown in case of errors
        INFO_LOG << "Shutting down after encountering errors during multi-view rendering.";
        m_appdata->set_app_run(false);
    }

    // Not yet performance value

    // {
    //     mi::base::Lock::Block block(&m_appdata->m_performance_values_lock);
    //     m_appdata->m_performance_values = mi::base::make_handle(frame_results->get_performance_values());
    //     m_appdata->m_application_compute_time   = application_compute_time;
    //     m_appdata->m_application_compute_memory = application_compute_memory;
    // }

    // // show statistics if needed
    // show_statistics_overlay(m_appdata->m_show_statistics_overlay, irc.m_icluster_configuration.get(),
    //                         m_appdata->m_performance_values.get(), m_appdata->m_statistics_overlay_large);

    // // Update recent n frames statistics
    // update_recent_n_statistics(m_appdata->m_performance_values.get());

    // // experimental timestep
    // render_frame_timestep_handling(irc, dice_transaction);
        
    // log performance value depends on the performance logging state
    // const mi::base::Handle<nv::index::IPerformance_values> perf_value(frame_results->get_performance_values());
    // log_performance_information(perf_value.get(), frame_num);

    // // performance measurement frames
    // if (m_appdata->get_remaining_performance_measuring_frames() > 0)
    // {
    //     m_appdata->decrement_remaining_performance_measuring_frames();
    //     if (m_appdata->get_remaining_performance_measuring_frames() == 0)
    //     {
    //         turn_off_logging();            // turn off the logging when reach to 0
    //         INFO_LOG << "Turn off the performance logging.";
    //     }
    // }

    // Display nice axis in the lower part of the screen
    if (m_appdata->is_show_coordinate_system_overlay())
    {
        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            irc.get_current_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
        assert(dice_transaction.is_valid_interface());
        {
            Multiple_stereo_camera* ms_camera = 
                m_appdata->get_user_interaction(0)->get_multiple_stereo_camera();
            gl_draw_axis(ms_camera->get_current_main_camera_tag(), dice_transaction.get());
        }
        dice_transaction->commit();
    }

    // Snapshot
    if (is_snapshot_on_this_frame)
    {
        save_span_buffer_snapshot(static_cast<Span_renderer_IF*>(irc.m_span_buffer.get()),
                                  irc.m_frame_num);
        // But turn off can be done through m_appdata since this is
        // the only single place to do that.
        m_appdata->set_snapshot_on(false); 
    }
}

//----------------------------------------------------------------------
bool Index_application::render_frame_multi_view_render_call(
    Nvindex_rendering_context& irc,
    bool                       is_snapshot_on_this_frame)
{
    // Render the subsurface data and display the horizontal spans using the user-defined span buffer renderer
    irc.m_frame_info_callbacks->clear_dynamic_allocation_events();

    // No OpenGL integration support in this mode
    if (m_appdata->is_enable_opengl_z_buffer_for_rendering())
    {
        static bool has_warned = false;
        if (!has_warned)
        {
            WARN_LOG << "No OpenGL integration support for multi-view mode.";
            has_warned = true;
        }
    }

    // set advisory to see the multi-view rendering details
    irc.m_viewport_list->set_advisory_enabled(m_appdata->is_enable_multi_view_advisory());

    // Main NVIDIA IndeX multi-view rendering call
    mi::base::Handle<nv::index::IFrame_results_list> frame_results_list(
        irc.m_iindex_rendering->render(
            irc.m_session_tag,
            irc.m_span_buffer.get(),
            irc.m_viewport_list.get()));
    assert(frame_results_list.is_valid_interface());

    bool is_success = true;
    if (frame_results_list->size() == 0)
    {
        ERROR_LOG << "IIndex_rendering rendering call has no results.";
        is_success = false;
    }
    else
    {
        const mi::Size nb_results = frame_results_list->size();
        for (mi::Size i = 0; i < nb_results; ++i)
        {
            mi::base::Handle<nv::index::IFrame_results>
                frame_results(frame_results_list->get(i));
            assert(frame_results.is_valid_interface());

            const mi::base::Handle<nv::index::IError_set> err_set(frame_results->get_error_set());
            assert(err_set.is_valid_interface());
            if (err_set->any_errors())
            {
                std::ostringstream os;

                const mi::base::Handle<nv::index::IError_set> err_set(frame_results->get_error_set());
                const mi::Uint32 nb_err = err_set->get_nb_errors();
                for (mi::Uint32 e = 0; e < nb_err; ++e)
                {
                    if (e != 0) os << '\n';
                    const mi::base::Handle<nv::index::IError> err(err_set->get_error(e));
                    os << err->get_error_string();
                }

                ERROR_LOG << "IIndex_rendering rendering call failed with the following error(s): " << '\n'
                          << os.str();
                is_success = false;
            }
        }
    }

    if (!is_success)
    {
        return false;           // some error in the rendering
    }

    // Show frame info callback's data
    render_frame_show_frame_info(irc.m_frame_info_callbacks.get());

    // restore the canvas size if snapshot
    render_frame_common_restore_canvas_size(is_snapshot_on_this_frame);
    
    // Visualize performance value with 2D OpenGL.
    render_frame_multi_view_visualize_perf_with_2d_opengl(irc);

    // Additional gl rendering: debug, visualize heightfield delete polygon
    // render_frame_single_view_additional_gl_rendering(irc, dice_transaction);

    // check GL error of this rendering
    render_frame_common_check_gl_error();

    return true;
}

//----------------------------------------------------------------------
void Index_application::render_frame_multi_view_visualize_perf_with_2d_opengl(
    Nvindex_rendering_context&        irc)
{
    // No multi-view performance visualization for now

    // End of 2D OpenGL rendering
    gl_pop_matrix_state();
}

//----------------------------------------------------------------------
//======================================================================
//----------------------------------------------------------------------
#ifdef NVINDEX_HAS_MPI_IPC_COMPUTE

mi::Uint32 Index_application::get_local_host_id(Nvindex_rendering_context& irc) const
{
    assert(irc.m_iindex_if.is_valid_interface());
    
    mi::base::Handle<nv::index::ICluster_configuration> rendering_properties_query(
        irc.m_iindex_if->get_api_component<nv::index::ICluster_configuration>());
        
    mi::Uint32 local_host_id = rendering_properties_query->get_local_host_id();

    if (local_host_id == 0)
    {
        // Only effective after NVindex::initialize() call.
        ERROR_LOG << "Local Host id not available yet. NVindex be initialized first";
    }

    return local_host_id;
}


void Index_application::set_bricks_affinity(
    Nvindex_rendering_context& irc,
    const std::vector<Brick_affinity> &bricks_affinity)
{
    // setup affinity information
    nv::index_common::Affinity_information* affinity_information = new nv::index_common::Affinity_information;
    for (mi::Uint32 i=0; i<bricks_affinity.size(); ++i)
    {
        affinity_information->add_affinity(
            bricks_affinity[i].bbox,
            bricks_affinity[i].host_id,
            bricks_affinity[i].device_id);
    }

    // The library takes ownership
    irc.m_iindex_session->set_affinity_information(affinity_information);
    
    // DiCE database access
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        irc.m_global_scope->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());

    // Update IndeX session so that data distribution takes place
    irc.m_iindex_session->update(irc.m_session_tag, dice_transaction.get());
    dice_transaction->commit();
    
}

void Index_application::set_bricks_locality(const std::vector<Brick_locality> &bricks_locality)
{
    m_bricks_locality.clear();
    m_bricks_locality = bricks_locality;
}

void Index_application::edit_volume_data(
    mi::neuraylib::Tag                                           volume_tag,
    mi::neuraylib::Tag                                           session_tag,
    mi::neuraylib::IDice_transaction*                            dice_transaction)
{
    // Access the session
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<nv::index::ISession>(session_tag));
    assert(session.is_valid_interface());

    // Access the volume scene element
    mi::base::Handle<const nv::index::IRegular_volume> volume(
        dice_transaction->access<nv::index::IRegular_volume>(volume_tag));
    assert(volume.is_valid_interface());

    // The entire region of interest is used for computing the volume data
    const mi::math::Bbox<mi::Uint32, 3> query_bbox = volume->get_IJK_bounding_box();

    // Access the distribution scheme
    const mi::neuraylib::Tag dist_layout_tag = session->get_distribution_layout();
    assert(dist_layout_tag.is_valid());
    
    mi::base::Handle<const nv::index::IData_distribution> distribution_layout(
        dice_transaction->access<nv::index::IData_distribution>(dist_layout_tag));
    assert(distribution_layout.is_valid_interface());

    // Retrieve configuration to get the border size
    mi::base::Handle<const nv::index::IConfig_settings> config_settings(
        dice_transaction->access<nv::index::IConfig_settings>(session->get_config()));
    assert(config_settings.is_valid_interface());

    // Distribution layout
    mi::base::Handle<nv::index::IRegular_volume_data_locality> data_locality(
        distribution_layout->retrieve_data_locality(volume_tag, query_bbox, dice_transaction));

    // Determine the host ids of all hosts in the cluster which store parts of the volume data
    for (mi::Uint32 i=0; i < data_locality->get_nb_cluster_nodes(); ++i)
    {
        // Host id 0 specifies that this subregion has not yet been loaded.
        // This should not happen.
        assert(data_locality->get_cluster_node(i) != 0);
        const mi::Uint32 host_id = data_locality->get_cluster_node(i);
        // INFO_LOG << "Cluster host "
        // << data_locality->get_cluster_node(i)
        // << " is responsible for processing "
        // << data_locality->get_nb_bounding_box(host_id)
        // << " problems.";
    }

    // Set up distributed editing process.
    const mi::Float64 start_time = nv::index_common::get_time();
    mi::base::Handle<nv::index::IDistributed_compute_algorithm> task(
        new CUDA_IPC_volume_editing(
            session_tag,
            volume_tag,
            query_bbox,
            volume->get_volume_size(),
            config_settings->get_subcube_border_size(),
            data_locality.get(),
            m_bricks_locality));

    // Start the fragmented job
    dice_transaction->execute_fragmented(task.get(), task->get_nb_of_fragments());
    
    // Performance measurements
    const mi::Float64 end_time = nv::index_common::get_time();
    const mi::Float64 elapsed_time = end_time - start_time;
    // INFO_LOG << "Total elapsed time for volume edit: " << elapsed_time;
}

#endif // NVINDEX_HAS_MPI_IPC_COMPUTE

void Index_application::unstructured_volume_compute(
    mi::Uint32                                          frame_id,
    const std::string&                                  db_dataset_name,
    mi::neuraylib::Tag                                  session_tag,
    mi::neuraylib::IDice_transaction*                   dice_transaction)
{
    // Set up distributed editing process.
    const mi::Float64 start_time = nv::index_common::get_time();

    const mi::neuraylib::Tag_struct irregular_volume_tag = dice_transaction->name_to_tag(db_dataset_name.c_str());

    // Access the session
    mi::base::Handle<const nv::index::ISession> session(dice_transaction->access<nv::index::ISession>(session_tag));
    assert(session.is_valid_interface());

    // Access the irregular volume scene element
    mi::base::Handle<const nv::index::IIrregular_volume_scene_element> irregular_volume(
        dice_transaction->access<nv::index::IIrregular_volume_scene_element>(irregular_volume_tag));
    assert(irregular_volume.is_valid_interface());

    // The entire region of interest is used for computing the irregular volume data
    const mi::math::Bbox<mi::Float32, 3> query_bbox = irregular_volume->get_bounding_box();

    DEBUG_LOG << "Compute for datset " << db_dataset_name.c_str() << " with tag id " << irregular_volume_tag.id << " will be launched for bounding box: " << query_bbox;

    // Access the distribution scheme
    const mi::neuraylib::Tag dist_layout_tag = session->get_distribution_layout();
    assert(dist_layout_tag.is_valid());

    mi::base::Handle<const nv::index::IData_distribution> distribution_layout(dice_transaction->access<nv::index::IData_distribution>(dist_layout_tag));
    assert(distribution_layout.is_valid_interface());

    // Distribution layout
    mi::base::Handle<nv::index::IIrregular_volume_data_locality> data_locality(
        distribution_layout->retrieve_data_locality(irregular_volume_tag, query_bbox, dice_transaction));

    // Determine the host ids of all hosts in the cluster which store parts of the volume data
    for (mi::Uint32 i = 0; i < data_locality->get_nb_cluster_nodes(); ++i)
    {
        assert(data_locality->get_cluster_node(i) != 0);
        const mi::Uint32 host_id = data_locality->get_cluster_node(i);
        DEBUG_LOG << "Cluster host "
                << data_locality->get_cluster_node(i)
                << " is responsible for processing "
                << data_locality->get_nb_bounding_box(host_id)
                << " subsets.";
    }

   mi::base::Handle<nv::index::IDistributed_compute_algorithm> task(
        new Irregular_volume_data_editing(
            frame_id,
            session_tag,
            irregular_volume_tag,
            query_bbox,
            data_locality.get()));
    // Start the fragmented job
    dice_transaction->execute_fragmented(task.get(), task->get_nb_of_fragments());

    // Performance measurements
    const mi::Float64 end_time = nv::index_common::get_time();
    const mi::Float64 elapsed_time = end_time - start_time;
    DEBUG_LOG << "Total elapsed time for the volume update: " << elapsed_time;
}


//----------------------------------------------------------------------
void Index_application::setup_multi_view_mode(
    Nvindex_rendering_context&        irc,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    const mi::Sint32 nb_viewports =
        get_sint32(m_appdata->peek_app_proj()->get("app::viewport::count", "0"));
    m_appdata->set_enable_multi_view_mode(nb_viewports > 0);
    if (nb_viewports <= 0)
    {            
        return; // multi-view mode is off
    }

    const bool is_advisory = 
        get_bool(m_appdata->peek_app_proj()->get("app::viewport::advisory", "off"));
    m_appdata->set_enable_multi_view_advisory(is_advisory);

    assert(!irc.m_viewport_list.is_valid_interface()); // should not yet setup
    assert(dice_transaction!= 0);

    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<nv::index::ISession>(irc.m_session_tag));
    assert(session.is_valid_interface());

    irc.m_viewport_list = session->create_viewport_list();
    assert(irc.m_viewport_list.is_valid_interface());

    const std::string viewport_prefix = "app::viewport::";

    // create viewports
    for (mi::Sint32 i = 0; i < nb_viewports; ++i)
    {
        mi::base::Handle<nv::index::IViewport> viewport(session->create_viewport());
        assert(viewport.is_valid_interface());

        std::stringstream pos_sstr;
        std::stringstream size_sstr;
        // get viewport anchor_position
        {
            pos_sstr << viewport_prefix << i << "::anchor_position";
            if (!m_appdata->peek_app_proj()->is_defined(pos_sstr.str()))
            {
                ERROR_LOG << "Could not find entry [" << pos_sstr.str() << "], may fail to run multi-view mode.";
                continue;
            }
        }
        // get viewport size
        {
            size_sstr << viewport_prefix << i << "::size";
            if (!m_appdata->peek_app_proj()->is_defined(size_sstr.str()))
            {
                ERROR_LOG << "Could not find entry [" << size_sstr.str() << "], may fail to run multi-view mode.";
                continue;
            }
        }

        const mi::math::Vector<mi::Sint32, 2> anchor_position = get_vec_sint32_2(m_appdata->peek_app_proj()->get(pos_sstr.str()));
        const mi::math::Vector<mi::Sint32, 2> resolution_size = get_vec_sint32_2(m_appdata->peek_app_proj()->get(size_sstr.str()));
        if ((resolution_size.x <= 0) || (resolution_size.y <= 0))
        {
            ERROR_LOG << "Illegal viewport size [" << resolution_size << "], may fail to run multi-view mode.";
            continue;
        }

        viewport->set_position(anchor_position);
        viewport->set_size(resolution_size);

        if (i > 0 && (i - 1) < static_cast<mi::Sint32>(irc.m_extra_scopes.size()))
        {
            // Use extra scopes when available
            viewport->set_scope(irc.m_extra_scopes[i - 1].get());
        }
        else
        {
            // Use global scope
            viewport->set_scope(irc.m_global_scope.get());
        }

        irc.m_viewport_list->append(viewport.get());
    }

    for (mi::Size i = 0; i < irc.m_viewport_list->size(); ++i)
    {
        mi::base::Handle<nv::index::IViewport> viewport(irc.m_viewport_list->get(i));
        INFO_LOG << "Multi-view mode: Viewport " << i << ": "
                 << viewport->get_position() << ": " << viewport->get_size();
    }
}

//----------------------------------------------------------------------
