/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/

#include "rtmp_call_event_handler.h"

#include <iomanip>
#include <sstream>
#include <ctype.h>

#include <nv/index/icolormap.h>
#include <nv/index/iconfig_settings.h>
#include <nv/index/iindex_debug_configuration.h>
#include <nv/index/imaterial.h>

#include "common/clock_pulse_generator.h"
#include "common/colormap_io.h"
#include "common/scenegraph_utility.h"

#include "command_handler_colormap.h"
#include "command_handler_scene_element.h"
#include "command_handler_settings.h"
#include "config_utility.h"
#include "distributed_command_execution.h"
#include "examiner_manipulator.h"
#include "heightfield_seed_manipuration.h"
#include "heightfield_workflow_functionality.h"
#include "performance_log.h"
#include "rtc_parameter_buffer_manip.h"
#include "rtmp_handler.h"
#include "shapes_utility.h"
#include "span_renderer_no_gl.h"

using namespace nv::index_common;

namespace {

/// Set resolution with command line
bool command_set_resolution(
    const std::string&                  com_str,
    const Nvindex_rendering_context&    irc,
    const mi::neuraylib::Tag&           session_tag)
{
    const std::string fun_name = "command_set_resolution: ";
    const std::string com_name = "set_resolution ";
    const std::string param = (com_str.substr(com_name.size()));
    if(param.size() == 0)
    {
        ERROR_LOG << ("Empty parameter: " + com_str);
        return false;
    }

    mi::Sint32_2 canvas_res(0, 0);
    {
        mi::Sint32 width  = -1;
        mi::Sint32 height = -1;
        std::stringstream instr(param);
        instr >> width >> height;
        if(instr.fail())
        {
            ERROR_LOG << ("Invalid parameters: " + com_str);
            return false;
        }
        canvas_res.x = width;
        canvas_res.y = height;
    }

    const std::string video_codec = Nvindex_AppData::instance()->peek_app_proj()->
        get("dice::rtmp_video_streaming::video_codec");
    if((video_codec == "h264") || (video_codec == "h264_nvenc")){
        // h264 software encoder has a minimal encode buffer resolution (128, 128)
        if(canvas_res.x < 128){
            WARN_LOG << "requested x resolution is too small for h264 encoder, set to 128.";
            canvas_res.x = 128;
        }
        if(canvas_res.y < 128){
            WARN_LOG << "requested y resolution is too small for h264 encoder, set to 128.";
            canvas_res.y = 128;
        }

        // also force to even number resolution
        if((canvas_res.x % 2) != 0){
            ++canvas_res.x;
            WARN_LOG << "cannot request odd number x resolution for for h264 encoder, set to " << canvas_res.x;
        }
        if((canvas_res.y % 2) != 0){
            ++canvas_res.y;
            WARN_LOG << "cannot request odd number y resolution for for h264 encoder, set to " << canvas_res.y;
        }
    }

    INFO_LOG << fun_name << com_name << ", width = " << canvas_res.x << ", height = " << canvas_res.y;

    const bool reset_aspect = get_bool(Nvindex_AppData::instance()->peek_app_proj()->get("index::update_camera_aspect_ratio_on_resolution_change", "0"));
    if(reset_aspect)
    {
        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            irc.m_viewing_scenario_scopes[0]->create_transaction<mi::neuraylib::IDice_transaction>());
        assert(dice_transaction.is_valid_interface());
        {
            mi::base::Handle<const nv::index::ISession> session(
                dice_transaction->access<const nv::index::ISession>(session_tag));
            assert(session.is_valid_interface());

            mi::base::Handle<const nv::index::IScene> scene(
                dice_transaction->access<const nv::index::IScene>(session->get_scene()));
            assert(scene.is_valid_interface());

            mi::base::Handle<nv::index::ICamera> camera(
                dice_transaction->edit<nv::index::ICamera>(scene->get_camera()));
            assert(camera.is_valid_interface());

            mi::base::Handle<nv::index::IPerspective_camera> camera_perspective(camera->get_interface<nv::index::IPerspective_camera>());
            if(camera_perspective)
            {
                const mi::Float32 aspect = static_cast<mi::Float32>(canvas_res.x)/static_cast<mi::Float32>(canvas_res.y);
                camera_perspective->set_aspect(aspect);
            }
        }
        dice_transaction->commit();
    }

    request_canvas_resolution(canvas_res);
    // set the canvas resolution for the application
    Nvindex_AppData::instance()->peek_app_proj()->insert("index::canvas_resolution", nv::index_common::to_string(canvas_res));
    // inform resize to the frame event handlers
    Frame_event_handler::set_resized(true);
    INFO_LOG << "Resize processed. Stopped the stream. Please re-load the viewer in the browser.";

    return true;
}

//----------------------------------------------------------------------
/// Set image file resolution with command line
bool command_set_image_file_resolution(const std::string & com_str)
{
    const std::string fun_name = "command_set_image_file_resolution: ";
    const std::string com_name = "set_image_file_resolution ";
    const std::string param = (com_str.substr(com_name.size()));
    if(param.size() == 0)
    {
        ERROR_LOG << ("Empty parameter: " + com_str);
        return false;
    }

    mi::Sint32_2 image_res(0, 0);
    {
        mi::Sint32 width  = -1;
        mi::Sint32 height = -1;
        std::stringstream instr(param);
        instr >> width >> height;
        if(instr.fail())
        {
            ERROR_LOG << ("Invalid parameters: " + com_str);
            return false;
        }
        image_res.x = width;
        image_res.y = height;
    }

    INFO_LOG << fun_name << com_name
             << ", width = " << image_res.x << ", height = " << image_res.y;

    // set the image canvas resolution for the application
    set_image_canvas_resolution(image_res);

    return true;
}

//----------------------------------------------------------------------
void change_camera_for_scenario(
    const Nvindex_rendering_context& irc,
    const mi::neuraylib::Tag&        session_tag)
{
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        irc.m_viewing_scenario_scopes[0]->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    {
        mi::base::Handle<const nv::index::ISession> session(
            dice_transaction->access<const nv::index::ISession>(session_tag));
        assert(session.is_valid_interface());

        mi::base::Handle<const nv::index::IScene> scene(
            dice_transaction->access<const nv::index::IScene>(session->get_scene()));
        assert(scene.is_valid_interface());

        dice_transaction->localize(scene->get_camera(), 1);
        mi::base::Handle<nv::index::ICamera> camera(
            dice_transaction->edit<nv::index::ICamera>(scene->get_camera()));
        assert(camera.is_valid_interface());

        mi::math::Vector_struct<mi::Float32, 3> origin = camera->get_eye_point();
        origin.x = origin.x+0.1f;
        camera->set_eye_point(origin);
    }
    dice_transaction->commit();
}

void enable_volume_slice_for_scenario(
    const Nvindex_rendering_context& irc)
{
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        irc.m_viewing_scenario_scopes[0]->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    {
        const std::vector<mi::neuraylib::Tag> volume_tag_vec = Nvindex_AppData::instance()->get_volume_tag_vec();
        const mi::Uint32 volume_id = 0;
        // actually when volume_id != 0 case, assert(volume_id < volume_tag_vec.size()),
        // but currently we have only volume_id == 0. So cppcheck complain the
        // assert(volume_id < volume_tag_vec.size()).
        assert(!volume_tag_vec.empty());

        mi::base::Handle<const nv::index::IRegular_volume> volume(
            dice_transaction->access<const nv::index::IRegular_volume>(volume_tag_vec.at(volume_id)));

        const mi::Uint32 slice_id = 0;
        const mi::neuraylib::Tag slice_tag = volume->get_slice(slice_id);

        dice_transaction->localize(slice_tag, 1); // TODO: needed?
        mi::base::Handle<nv::index::ISlice_scene_element> slice(
            dice_transaction->edit<nv::index::ISlice_scene_element>(slice_tag));

        slice->set_enabled(!slice->get_enabled());
    }
    dice_transaction->commit();
}

void move_volume_slice_for_scenario(
    const Nvindex_rendering_context& irc,
    mi::Float32                      value)
{
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        irc.m_viewing_scenario_scopes[0]->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());

    {
        const std::vector<mi::neuraylib::Tag> volume_tag_vec = Nvindex_AppData::instance()->get_volume_tag_vec();
        const mi::Uint32 volume_id = 0;
        // actually when volume_id != 0 case, assert(volume_id < volume_tag_vec.size()),
        // but currently we have only volume_id == 0. So cppcheck complain the
        // assert(volume_id < volume_tag_vec.size()).
        assert(!volume_tag_vec.empty());
        mi::base::Handle<const nv::index::IRegular_volume> volume(
            dice_transaction->access<const nv::index::IRegular_volume>(volume_tag_vec.at(volume_id)));

        const mi::Uint32 slice_id = 0;
        const mi::neuraylib::Tag slice_tag = volume->get_slice(slice_id);

        dice_transaction->localize(slice_tag, 1); // FIXME: needed?
        mi::base::Handle<nv::index::ISlice_scene_element> slice(
            dice_transaction->edit<nv::index::ISlice_scene_element>(slice_tag));
        mi::base::Handle<nv::index::ISection_scene_element> section_slice(
            slice->get_interface<nv::index::ISection_scene_element>());
        if(section_slice.get())
        {
            const mi::Float32 position = section_slice->get_position();
            section_slice->set_position(position+value);
        }
    }
    dice_transaction->commit();
}

void enable_region_of_interest_annotation(
    const mi::base::Handle<mi::neuraylib::IDice_transaction>&   dice_transaction,
    const mi::neuraylib::Tag&                                   session_tag)
{
    const std::string roi_scene_group = "roi_annotations";
    const mi::neuraylib::Tag group_tag = get_scene_group(session_tag, roi_scene_group, dice_transaction.get());

    if(group_tag!=mi::neuraylib::NULL_TAG)
    {
        mi::base::Handle<const nv::index::ISession> session(
            dice_transaction->access<nv::index::ISession>(session_tag));
        assert(session.is_valid_interface());

        mi::base::Handle<const nv::index::IScene> scene(
            dice_transaction->access<const nv::index::IScene>(session->get_scene()));
        assert(scene.is_valid_interface());

        mi::base::Handle<nv::index::ITransformed_scene_group> group(
            dice_transaction->edit<nv::index::ITransformed_scene_group>(group_tag));

        { // Display the scene's region of interest
            const mi::math::Bbox<mi::Float32, 3> roi = scene->get_clipped_bounding_box();
            mi::math::Color color(0.8f, 0.8f, 0.8f, 0.8f);
            const mi::neuraylib::Tag line_set_tag = create_line_set(roi, color, dice_transaction, session_tag, "scene_roi");
            group->append(line_set_tag, dice_transaction.get());
        }

        const std::vector<mi::neuraylib::Tag> volume_tag_vec = Nvindex_AppData::instance()->get_volume_tag_vec();
        for(std::vector<mi::neuraylib::Tag>::const_iterator vi = volume_tag_vec.begin();
            vi != volume_tag_vec.end(); ++vi)
        {
            mi::base::Handle<const nv::index::IRegular_volume> volume_scene_element(
                dice_transaction->access<const nv::index::IRegular_volume>(*vi));

            { // Display the region of interest
                // TODO: const mi::math::Matrix<mi::Float32, 4, 4> mat = volume_scene_element->get_transform();
                mi::math::Bbox<mi::Float32, 3> roi = volume_scene_element->get_IJK_region_of_interest();
                // Border voxels are half the size of interiour ones, so adapt for rendering wireframe box.
                roi.max -= mi::math::Vector<mi::Float32, 3>(1.f);
                const mi::math::Color color(0.7f, 0.7f, 0.7f, 0.8f);
                std::string name = "volume_roi_" + std::string(volume_scene_element->get_name());
                const mi::neuraylib::Tag line_set_tag = create_line_set(roi, color, dice_transaction, session_tag, name);
                group->append(line_set_tag, dice_transaction.get());
            }
        }

        std::vector<mi::neuraylib::Tag> height_tag_vec = Nvindex_AppData::instance()->get_heightfield_tag_vec();
        for(std::vector<mi::neuraylib::Tag>::const_iterator hi = height_tag_vec.begin();
            hi != height_tag_vec.end(); ++hi)
        {
            assert(hi->is_valid());
            mi::base::Handle<const nv::index::IRegular_heightfield> heightfield(
                dice_transaction->access<const nv::index::IRegular_heightfield>(*hi));
            mi::math::Color color(1.0f);
            { // Display the region of interest
                // TODO: const mi::math::Matrix<mi::Float32, 4, 4> mat = heightfield->get_transform();
                mi::math::Bbox<mi::Float32, 3> roi = heightfield->get_IJK_region_of_interest();
                // Border voxels are half the size of interiour ones, so adapt for rendering wireframe box.
                roi.max -= mi::math::Vector<mi::Float32, 3>(1.0f);
                std::string name = "heightfield_roi_" + std::string(heightfield->get_name());
                const mi::neuraylib::Tag line_set_tag = create_line_set(roi, color, dice_transaction, session_tag, name);
                group->append(line_set_tag, dice_transaction.get());
            }
        }
    }
}

void disable_region_of_interest_annotation(
    const mi::base::Handle<mi::neuraylib::IDice_transaction>&   dice_transaction,
    const mi::neuraylib::Tag&                                   session_tag)
{
    const std::string roi_scene_group = "roi_annotations";
    const mi::neuraylib::Tag group_tag = get_scene_group(session_tag, roi_scene_group, dice_transaction.get());
    if(group_tag!=mi::neuraylib::NULL_TAG)
    {
        mi::base::Handle<nv::index::ITransformed_scene_group> group(
            dice_transaction->edit<nv::index::ITransformed_scene_group>(group_tag));
        group->clear(dice_transaction.get());
    }
}

void enable_bounding_box_annotation(
    const mi::base::Handle<mi::neuraylib::IDice_transaction>&   dice_transaction,
    const mi::neuraylib::Tag&                                   session_tag)
{
    const std::string roi_scene_group = "bb_annotations";
    const mi::neuraylib::Tag group_tag = get_scene_group(session_tag, roi_scene_group, dice_transaction.get());

    if(group_tag!=mi::neuraylib::NULL_TAG)
    {
        mi::base::Handle<nv::index::ITransformed_scene_group> group(
            dice_transaction->edit<nv::index::ITransformed_scene_group>(group_tag));

        {
            // Display the scene's region of interest
            const mi::math::Bbox<mi::Float32, 3> roi = 
                get_scene_bounding_box(session_tag, dice_transaction.get());

            mi::math::Color color(0.8f, 0.8f, 0.8f, 0.8f);
            const mi::neuraylib::Tag line_set_tag = create_line_set(roi, color, dice_transaction, session_tag, "scene_bb");
            group->append(line_set_tag, dice_transaction.get());
        }

        std::vector<mi::neuraylib::Tag> volume_tag_vec = Nvindex_AppData::instance()->get_volume_tag_vec();
        for(std::vector<mi::neuraylib::Tag>::const_iterator vi = volume_tag_vec.begin();
            vi != volume_tag_vec.end(); ++vi)
        {
            mi::base::Handle<const nv::index::IRegular_volume> volume_scene_element(
                dice_transaction->access<const nv::index::IRegular_volume>(*vi));

            { // Display the bounding box in global space
                const mi::math::Bbox<mi::Float32, 3> roi = volume_scene_element->get_XYZ_bounding_box();
                mi::math::Color color(0.9f, 0.5f, 0.5f, 0.8f);
                // std::string name = "volume_bb_" + std::string(volume_scene_element->get_name());
                const mi::neuraylib::Tag line_set_tag = create_line_set(roi, color, dice_transaction, session_tag, "name");
                group->append(line_set_tag, dice_transaction.get());
            }
        }

        const std::vector<mi::neuraylib::Tag> svol_tag_vec = Nvindex_AppData::instance()->get_sparse_volume_tag_vec();
        mi::Uint32 svol_idx = 0u;
        for(std::vector<mi::neuraylib::Tag>::const_iterator vi = svol_tag_vec.begin();
            vi != svol_tag_vec.end(); ++vi)
        {
            mi::base::Handle<const nv::index::ISparse_volume_scene_element> svol_scene_element(
                dice_transaction->access<const nv::index::ISparse_volume_scene_element>(*vi));

            { // Display the region of interest
                // TODO: const mi::math::Matrix<mi::Float32, 4, 4> mat = volume_scene_element->get_transform();
                mi::math::Bbox<mi::Float32, 3> roi = svol_scene_element->get_bounding_box();
                // Border voxels are half the size of interiour ones, so adapt for rendering wireframe box.
                roi.max -= mi::math::Vector<mi::Float32, 3>(1.f);
                const mi::math::Color color(1.0f, 0.7f, 0.7f, 0.8f);
                std::stringstream name;
                name << "svol_bbox_" << svol_idx++;
                //std::string name = "svol_roi_" + std::string(svol_scene_element->get_name());
                const mi::neuraylib::Tag line_set_tag = create_line_set(roi, color, dice_transaction, session_tag, name.str());
                group->append(line_set_tag, dice_transaction.get());
            }
        }

        std::vector<mi::neuraylib::Tag> height_tag_vec = Nvindex_AppData::instance()->get_heightfield_tag_vec();
        for(std::vector<mi::neuraylib::Tag>::const_iterator hi = height_tag_vec.begin();
            hi != height_tag_vec.end(); ++hi)
        {
            assert(hi->is_valid());
            mi::base::Handle<const nv::index::IRegular_heightfield> heightfield(
                dice_transaction->access<const nv::index::IRegular_heightfield>(*hi));
            mi::math::Color color(1.f);
            { // Display the region of interest
                mi::math::Bbox<mi::Float32, 3> roi = heightfield->get_XYZ_clipped_bounding_box();
                std::string name = "heightfield_bb_" + std::string(heightfield->get_name());
                const mi::neuraylib::Tag line_set_tag = create_line_set(roi, color, dice_transaction, session_tag, name);
                group->append(line_set_tag, dice_transaction.get());
            }
        }
    }
}

void disable_bounding_box_annotation(
    const mi::base::Handle<mi::neuraylib::IDice_transaction>&   dice_transaction,
    const mi::neuraylib::Tag&                                   session_tag)
{
    const std::string roi_scene_group = "bb_annotations";
    const mi::neuraylib::Tag group_tag = get_scene_group(session_tag, roi_scene_group, dice_transaction.get());
    if(group_tag!=mi::neuraylib::NULL_TAG)
    {
        mi::base::Handle<nv::index::ITransformed_scene_group> group(
            dice_transaction->edit<nv::index::ITransformed_scene_group>(group_tag));
        group->clear(dice_transaction.get());
    }
}

//----------------------------------------------------------------------
void get_quality_level(
    const std::string&                                   level,
    nv::index::IConfig_settings::Volume_filtering_modes& filter_mode,
    bool&                                                ray_segment_accum)
{
    std::string default_filter = "nearest";
    std::string default_ray_segment_accum = "no";
    if (level == "1")
    {
        // Medium quality
        default_filter = "trilinear_post_hw";
        default_ray_segment_accum = "yes";
    }
    else if (level == "2")
    {
        // High quality
        default_filter = "tricubic_bspline_post_hw";
        default_ray_segment_accum = "yes";
    }

    const String_dict* prj = Nvindex_AppData::instance()->peek_app_proj();
    const std::string prefix = "app::quality_level::" + level + "::";

    using nv::index::IConfig_settings;
    const std::string filter_name = prj->get(prefix + "volume_filter", default_filter);
    filter_mode = IConfig_settings::VOLUME_FILTER_NEAREST;
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
        ERROR_LOG << "Invalid volume filter mode specified for " << prefix << "volume_filter: " << filter_name;

    ray_segment_accum = get_bool(prj->get(prefix + "ray_segment_accum", default_ray_segment_accum));
}

} // namespace

//////////////////// Call_event_handler ////////////////////

bool Call_event_handler::s_is_recording = false;
mi::base::Lock Call_event_handler::s_recording_lock;
std::ofstream Call_event_handler::s_recording_file;

Call_event_handler::Call_event_handler(
    mi::Uint32                                       viewing_scenario,
    Nvindex_rendering_context&                       irc)
  : m_viewing_scenario(viewing_scenario),
    m_irc(irc),
    m_iindex_if(irc.m_iindex_if),
    m_iindex_session(irc.m_iindex_session),
    m_iindex_rendering(irc.m_iindex_rendering),
    m_session_tag(irc.m_session_tag),
    m_progress_callback(irc.m_progress_callback),
    m_frame_info_callbacks(irc.m_frame_info_callbacks),
    m_iindex_scene_query(irc.m_iindex_scene_query)
{
    Nvindex_AppData* appdata = Nvindex_AppData::instance();

    // Set up the command handlers
    m_command_handlers.push_back(new Command_handler_colormap(irc, appdata));
    m_command_handlers.push_back(new Command_handler_scene_element(irc, appdata));
    m_command_handlers.push_back(new Command_handler_settings(irc, appdata));

    // GUI initial state
    // copy volume
    m_gui_state_dict.insert("attrgen_copy_source_volume_index", "0");
    // filter volume
    m_gui_state_dict.insert("attrgen_filter_source_volume_index", "0");
    m_gui_state_dict.insert("attrgen_filter_name_index", "0");
    // edit volume
    m_gui_state_dict.insert("attrgen_edit_volume_index", "0");
}

Call_event_handler::~Call_event_handler()
{
    for (size_t i=0; i < m_command_handlers.size(); ++i)
        delete m_command_handlers[i];
}

bool Call_event_handler::handle( //-V553 PVS
    mi::rtmp::IConnection*      /*connection*/,
    const char*                 /*procedure_name*/,
    const mi::IData*            /*command_args*/,
    const mi::IData*            user_args,
    mi::IData**                 response_args)
{
    // Unwrap the user arguments
    mi::base::Handle<const mi::IMap> imap(user_args->get_interface<const mi::IMap>());
    if (!imap.is_valid_interface())
        return false;

    // Determine the requested command
    mi::base::Handle<const mi::IString> cmd_istr(imap->get_value<mi::IString>("cmd"));
    if (!cmd_istr.is_valid_interface())
        return false;
    const std::string cmd_str = cmd_istr->get_c_str();

    // Record the command if recording is enabled
    record_received_command(cmd_str, imap);

    // New command handler (version 2)
    mi::base::Handle<const mi::IString> cmd_version(imap->get_value<mi::IString>("version"));
    if (cmd_version.is_valid_interface() && std::string(cmd_version->get_c_str()) == "2")
    {
        // Call the command handlers
        Call_event_arguments call_args(imap);
        for (size_t i=0; i < m_command_handlers.size(); ++i)
        {
            if (m_command_handlers[i]->handle(cmd_str, call_args))
            {
                // Command was handled successfully, set the response
                *response_args = call_args.generate_response(m_iindex_if);
                return true;
            }
            call_args.clear_response(); // Clean up, just in case
        }
    }

    //
    // Old-style handling code follows below (no version field).
    // Note: New code should use the command handler mechanism above.
    //

    Arguments_get_wrapper user_arguments(imap);
    Arguments_set_wrapper response_arguments(response_args, m_iindex_if);

    const char* cmd = cmd_str.c_str();
    Nvindex_AppData* appdata = Nvindex_AppData::instance();

    mi::base::Handle<mi::neuraylib::IScope> cur_scope(get_scope());

    //
    // Run the command
    //
    if (strcmp(cmd, "mouseMove") == 0)
    {
        const char* pos_x = user_arguments.get_value("x");
        const char* pos_y = user_arguments.get_value("y");
        {
            mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
                cur_scope->create_transaction<mi::neuraylib::IDice_transaction>());
            assert(dice_transaction.is_valid_interface());
            {
                mi::base::Handle<const nv::index::ISession> session(
                    dice_transaction->access<const nv::index::ISession>(m_session_tag));
                assert(session.is_valid_interface());

                Nvindex_AppData::instance()->get_user_interaction(0)->
                    mouse_motion(Examiner_manipulator::Left_button,
                                 get_sint32(pos_x), 
                                 get_sint32(pos_y),
                                 session->get_scene(),
                                 dice_transaction.get());
            }
            dice_transaction->commit();
        }
    }
    else if (strcmp(cmd, "mouseUp") == 0)
    {
        const char* pos_x  = user_arguments.get_value("x");
        const char* pos_y  = user_arguments.get_value("y");
        {
            mi::math::Vector_struct<mi::Sint32, 2> pick_location;
            {
                mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
                    cur_scope->create_transaction<mi::neuraylib::IDice_transaction>());
                assert(dice_transaction.is_valid_interface());
                {
                    pick_location = Nvindex_AppData::instance()->get_user_interaction(0)->
                        mouse_button_action(Examiner_manipulator::Left_button,
                                            Examiner_manipulator::Button_up,
                                            get_sint32(pos_x), 
                                            get_sint32(pos_y),
                                            dice_transaction.get());
                }
                dice_transaction->commit();
            }

            if ((pick_location.x >= 0) && (pick_location.y >= 0))
            {
                const mi::math::Vector<mi::Uint32, 2> pick_location_u32(
                    static_cast<mi::Uint32>(pick_location.x), static_cast<mi::Uint32>(pick_location.y));

                if (Heightfield_workflow_functionality::instance()->is_heightfield_workflow_manual_pick())
                {
                    // Select the active volume to be used for slice picking
                    {
                        std::vector<mi::neuraylib::Tag> volume_tag_vec =
                            appdata->get_volume_tag_vec();

                        // Use the first available volume
                        mi::neuraylib::Tag volume_tag;
                        if (!(volume_tag_vec.empty()))
                        {
                            volume_tag = volume_tag_vec.at(0);
                            assert(volume_tag.is_valid());
                        }
                        Heightfield_workflow_functionality::instance()->set_volume(volume_tag);
                    }

                    // FIXME: separate viewport information && lock when access
                    const nv::index::IIndex_canvas* pick_canvas = m_irc.m_span_buffer.get();
                    heightfield_workflow_manual_picking(cur_scope, pick_location_u32, pick_canvas);

                }
                else if (Heightfield_workflow_functionality::instance()->is_heightfield_workflow_delete_operation()
                         ||  Heightfield_workflow_functionality::instance()->is_heightfield_workflow_elevation_change_operation()
                         ||  Heightfield_workflow_functionality::instance()->is_heightfield_workflow_gridding_operation() )
                {
                    // FIXME: separate viewport information && lock when access
                    const nv::index::IIndex_canvas* pick_canvas = m_irc.m_span_buffer.get();
                    heightfield_workflow_append_vertex_to_bounding_polygon(
                        cur_scope, m_iindex_scene_query, m_session_tag, pick_location_u32, pick_canvas);
                }
                else if (appdata->is_pick_operation_enabled())
                {
                    mi::neuraylib::Tag picked_tag = Nvindex_AppData::instance()->get_user_interaction(0)->
                        issue_pick_command(
                            &m_irc, 
                            Nvindex_AppData::instance()->is_enable_multi_view_mode(),
                            pick_location);

                    if (picked_tag.is_valid())
                    {
                        // Send the tag of the first picked element back
                        std::ostringstream os;
                        os << picked_tag;
                        response_arguments.set_value("pickedTag", os.str().c_str());
                    }
                }
            }
        }
    }
    else if (strcmp(cmd, "mouseDown") == 0)
    {
        const char* pos_x  = user_arguments.get_value("x");
        const char* pos_y  = user_arguments.get_value("y");
        {
            mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
                cur_scope->create_transaction<mi::neuraylib::IDice_transaction>());
            assert(dice_transaction.is_valid_interface());
            Nvindex_AppData::instance()->get_user_interaction(0)->
                mouse_button_action(Examiner_manipulator::Left_button,
                                    Examiner_manipulator::Button_down,
                                    get_sint32(pos_x), 
                                    get_sint32(pos_y),
                                    dice_transaction.get());
            dice_transaction->commit();
        }
    }
    else if (strcmp(cmd, "gesturePan") == 0)
    {
        const char* offset_x = user_arguments.get_value("offsetX");
        const char* offset_y = user_arguments.get_value("offsetY");
        INFO_LOG << "gesturePan: " << offset_x << " /" << offset_y;
    }
    else if (strcmp(cmd, "gestureZoom") == 0)
    {
        const char* scale_x = user_arguments.get_value("scaleX");
        const char* scale_y = user_arguments.get_value("scaleY");

        // Average x and y scale to get a single scale factor to work with
        mi::Float32 scale = (get_float32(scale_x) + get_float32(scale_y)) / 2.0f;

        //            INFO_LOG << "gestureZoom: " << scale_x << " /" << scale_y << " => " << scale;

        // Typical values for scale are 1.002 (for zooming in) or 0.976 (for zooming out)
        mi::Float32 zoom = -(1.f - scale) * 4.f;

        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            cur_scope->create_transaction<mi::neuraylib::IDice_transaction>());
        assert(dice_transaction.is_valid_interface());
        {
            mi::base::Handle< const nv::index::ISession >
                session(
                    dice_transaction->access< const nv::index::ISession >(
                        m_session_tag));
            assert(session.is_valid_interface());

            User_interaction* ui = Nvindex_AppData::instance()->get_user_interaction(0);
            ui->get_examiner()->camera_zoom(
                zoom,
                ui->get_multiple_stereo_camera()->get_main_base_camera_tag(),
                session->get_scene(),
                dice_transaction.get());
            ui->get_examiner()->camera_zoom(
                zoom,
                ui->get_multiple_stereo_camera()->get_ortho_camera_tag(),
                session->get_scene(),
                dice_transaction.get());
        }
        dice_transaction->commit();
    }
    else if (strcmp(cmd, "keyDown") == 0)
    {
        const std::string keycode  = user_arguments.get_value("keycode");
        const std::string shiftkey = user_arguments.get_value("shiftkey");
        const std::string ctrlkey  = user_arguments.get_value("ctrlkey");
        const std::string altkey   = user_arguments.get_value("altkey");

        const bool is_shift_on = (shiftkey == "true") ? true : false;
        const bool is_ctrl_on  = (ctrlkey  == "true") ? true : false;
        const bool is_alt_on   = (altkey   == "true") ? true : false;
        const bool is_meta_on  = false; // RTMP doesn't handle the meta key
        const mi::Sint32 keycode_s32 = get_sint32(keycode);

        bool is_key_handled = true;
        if (is_shift_on || is_ctrl_on || is_alt_on || iscntrl(keycode_s32))
        {
            // A special key is pushed.
            Nvindex_AppData::instance()->get_user_interaction(0)->
                key_down_action(&m_irc, keycode_s32, is_shift_on, is_ctrl_on, is_alt_on, is_meta_on);
        }
        else
        {
            // No special key is pushed.
            is_key_handled = Nvindex_AppData::instance()->get_user_interaction(0)->
                key_press_action(&m_irc, keycode_s32, is_shift_on, is_ctrl_on, is_alt_on, is_meta_on);
        }

        if (!is_key_handled)
        {
            bool is_key_handled_in_rtmp_handler = false;
            switch (keycode_s32)
            {
            case 'B':
            {
                if(m_viewing_scenario==1)
                {
                    INFO_LOG << "Changing camera position for viewing scenario " << m_viewing_scenario << " (multi-view, multi-user).";
                    change_camera_for_scenario(m_irc, m_session_tag);

                    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
                        m_irc.m_viewing_scenario_scopes[0]->create_transaction<mi::neuraylib::IDice_transaction>());
                    assert(dice_transaction.is_valid_interface());

                    render_viewing_scenario(
                        m_viewing_scenario,
                        dice_transaction);

                    dice_transaction->commit();
                }
                is_key_handled_in_rtmp_handler = true;
                break;
            }
            case 'N':
            {
                if(m_viewing_scenario==1)
                {
                    INFO_LOG << "En/disable slice for viewing scenario " << m_viewing_scenario << " (minor scene variation, multi-user).";
                    enable_volume_slice_for_scenario(m_irc);

                    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
                        m_irc.m_viewing_scenario_scopes[0]->create_transaction<mi::neuraylib::IDice_transaction>());
                    assert(dice_transaction.is_valid_interface());

                    render_viewing_scenario(
                        m_viewing_scenario,
                        dice_transaction);

                    dice_transaction->commit();
                }
                is_key_handled_in_rtmp_handler = true;
                break;
            }
            case 'M':
            {
                if(m_viewing_scenario==1)
                {
                    INFO_LOG << "Moving slice for viewing scenario " << m_viewing_scenario << " (minor scene variation, multi-user).";
                    if (shiftkey == "true")
                    {
                        move_volume_slice_for_scenario(m_irc, 10);
                    }
                    else
                    {
                        move_volume_slice_for_scenario(m_irc, -10);
                    }
                    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
                        m_irc.m_viewing_scenario_scopes[0]->create_transaction<mi::neuraylib::IDice_transaction>());
                    assert(dice_transaction.is_valid_interface());

                    render_viewing_scenario(
                        m_viewing_scenario,
                        dice_transaction);

                    dice_transaction->commit();
                }
                is_key_handled_in_rtmp_handler = true;
                break;
            }

            default:
                // no key handled
                break;
            }
            
            if (!is_key_handled_in_rtmp_handler)
            {
                ERROR_LOG << "The keycode " << keycode_s32 << " has not been handled.";
            }
        }
    }
    else if (strcmp(cmd, "keyUp") == 0)
    {
        const std::string keycode  = user_arguments.get_value("keycode");
        const std::string shiftkey = user_arguments.get_value("shiftkey");
        const std::string ctrlkey  = user_arguments.get_value("ctrlkey");
        const std::string altkey   = user_arguments.get_value("altkey");

        const bool is_shift_on = (shiftkey == "true") ? true : false;
        const bool is_ctrl_on  = (ctrlkey  == "true") ? true : false;
        const bool is_alt_on   = (altkey   == "true") ? true : false;
        const bool is_meta_on  = false; // RTMP doesn't handle the meta key
        const mi::Sint32 keycode_s32 = get_sint32(keycode);

        Nvindex_AppData::instance()->get_user_interaction(0)->
            key_up_action(&m_irc, keycode_s32, is_shift_on, is_ctrl_on, is_alt_on, is_meta_on);
    }
    else if (strcmp(cmd, "getoptions") == 0)
    {
        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            cur_scope->create_transaction<mi::neuraylib::IDice_transaction>());
        assert(dice_transaction.is_valid_interface());
        {
            mi::base::Handle<const nv::index::ISession> session(
                dice_transaction->access<const nv::index::ISession>(m_session_tag));

            mi::base::Handle<const nv::index::IConfig_settings> configs(
                dice_transaction->access<nv::index::IConfig_settings>(session->get_config()));

            {
                std::ostringstream os;
                os << configs->get_nb_spans();
                response_arguments.set_value("numhspans", os.str().c_str());
            }

            {
                std::ostringstream os;
                os << configs->get_size_of_rendering_results_in_queue();
                response_arguments.set_value("render_result_queue", os.str().c_str());
            }

            {
                std::ostringstream os;
                os << configs->is_automatic_span_control();
                response_arguments.set_value("autospancontrol", os.str().c_str());
            }

            if (configs->get_rendering_mode() == nv::index::IConfig_settings::CPU_RENDERING)
                response_arguments.set_value("cpu_rendering", "true");
            if (configs->get_rendering_mode() == nv::index::IConfig_settings::GPU_OR_CPU_RENDERING)
                response_arguments.set_value("cpu_fallback", "true");

            String_dict* prj = appdata->peek_app_proj();
            if (prj->is_defined("app::roi_animation_1"))
            {
                response_arguments.set_value("roi_animation_1", prj->get("app::roi_animation_1").c_str());
                response_arguments.set_value("roi_animation_2", prj->get("app::roi_animation_2").c_str());
                response_arguments.set_value("roi_animation_3", prj->get("app::roi_animation_3").c_str());
            }

#ifndef USE_OPENGL
            // Parallel span rendering is not possible (nor necessary) when using OpenGL
            response_arguments.set_value("parallelSpanRenderingSupported", "true");
#endif  // USE_OPENGL

            {
                std::ostringstream os;
                os << static_cast<int>(configs->get_volume_filter());
                response_arguments.set_value("filter", os.str().c_str());
            }

            response_arguments.set_value(
                "ray_segment_accum",
                configs->get_ray_segment_accumulation_technique() ? "on" : "off");

            // Find a visual quality level that matches the current settings
            for (int i = 0; i < 3; ++i)
            {
                std::ostringstream os;
                os << i;
                nv::index::IConfig_settings::Volume_filtering_modes filter;
                bool ray_segment_accum;
                get_quality_level(os.str(), filter, ray_segment_accum);

                if (filter == configs->get_volume_filter() &&
                    ray_segment_accum == configs->get_ray_segment_accumulation_technique())
                {
                    response_arguments.set_value("visual_quality_level", os.str().c_str());
                    break;
                }
            }

            if (prj->is_defined("index::timesteps"))
            {
                response_arguments.set_value("nb_timesteps", prj->get("index::timesteps").c_str());
            }

            {
                std::ostringstream os;
                os << m_irc.m_current_scope;
                response_arguments.set_value("current_scope", os.str().c_str());
            }

            if (!m_irc.m_extra_scopes.empty())
            {
                std::ostringstream os;
                os << (m_irc.m_extra_scopes.size() + 1);
                response_arguments.set_value("nb_scopes", os.str().c_str());
            }
        }
        dice_transaction->commit();
    }
    else if (strcmp(cmd, "setoptions") == 0)
    {
        // Create transaction in global scope to make sure bbox/ROI annotations are visible
        // everywhere
        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            m_irc.m_global_scope->create_transaction<mi::neuraylib::IDice_transaction>());
        assert(dice_transaction.is_valid_interface());

        const std::string on = "on"; // To simplify comparison to C-style string

        bool local_monitor_performance_values =
            (user_arguments.get_value("perfMonitor") && user_arguments.get_value("perfMonitor") == on);

        // update perforamnce logging state
        bool const log_perf_gui_state[Performance_logger_base::LK_COUNT] = {
            (user_arguments.get_value("logPerfGlobal") == on),
            (user_arguments.get_value("logPerfHost") == on),
            (user_arguments.get_value("logPerfSpan") == on),
        };
        Performance_logger_app_state & perf_log_state = appdata->peek_log_state();
        for(mi::Uint32 i = 0; i < Performance_logger_base::LK_COUNT; ++i)
        {
            // when log perf gui is turn off from on
            if(is_bool_state_change_to(perf_log_state.is_logging_on(i), log_perf_gui_state[i], false)){
                // Not implemented Perf_data::instance()->set_performance_log_need_close(i, true);
            }
            perf_log_state.set_logging_on(i, log_perf_gui_state[i]);
        }

        nv::index::IConfig_settings::Rendering_mode local_rendering_mode =
            static_cast<nv::index::IConfig_settings::Rendering_mode>(get_sint32(user_arguments.get_value("rendermode")));

        if (get_bool(appdata->peek_app_proj()->get("index::force_cpu_rendering", "no")))
            local_rendering_mode = nv::index::IConfig_settings::CPU_RENDERING;

        const mi::Uint32 local_cpu_thread_count = get_sint32(user_arguments.get_value("cputhreads"));
        const bool local_parallel_rendering_and_compositing = (user_arguments.get_value("parallel") == on);

        using nv::index::IConfig_settings;
        IConfig_settings::Compositing_mode local_compositing_mode = IConfig_settings::COMPOSITING_ALL;
        if (user_arguments.get_value("compositingMode") == std::string("1"))
            local_compositing_mode = IConfig_settings::COMPOSITING_LOCAL_ONLY;
        else if (user_arguments.get_value("compositingMode") == std::string("2"))
            local_compositing_mode = IConfig_settings::COMPOSITING_REMOTE_ONLY;

        appdata->m_immediate_final_parallel_compositing =
            (user_arguments.get_value("immediateFinalCompositing") == on);

        const bool parallel_span_rendering = (user_arguments.get_value("parallelSpanRendering") == on);
        m_irc.m_span_buffer->set_use_parallel_rendering(parallel_span_rendering);

        const std::string numhspans  = user_arguments.get_value("numhspans");
        mi::Uint32 local_nb_spans_per_host = (mi::Uint32)(get_float32(numhspans));

        const std::string rendering_result_queue  = user_arguments.get_value("render_result_queue");
        const mi::Uint32 local_size_of_rendering_result_queue =  (mi::Uint32)(get_float32(rendering_result_queue));

        // Enable disable bbox/ROI annotation
        bool local_visualize_horizontal_spans = (user_arguments.get_value("showhspans") == on);
        bool local_visualize_volume = (user_arguments.get_value("showvolume") == on);
        bool local_visualize_region_of_interest = (user_arguments.get_value("showvolumeofinterest") == on);
        if(local_visualize_region_of_interest && !appdata->is_show_region_of_interest_overlay())
            enable_region_of_interest_annotation(dice_transaction, m_session_tag);
        if(!local_visualize_region_of_interest && appdata->is_show_region_of_interest_overlay())
            disable_region_of_interest_annotation(dice_transaction, m_session_tag);
        appdata->set_show_region_of_interest_overlay(local_visualize_region_of_interest);

        if(local_visualize_volume && !appdata->is_show_bounding_boxes_overlay())
            enable_bounding_box_annotation(dice_transaction, m_session_tag);
        if(!local_visualize_volume && appdata->is_show_bounding_boxes_overlay())
            disable_bounding_box_annotation(dice_transaction, m_session_tag);
        appdata->set_show_bounding_boxes_overlay(local_visualize_volume);

        appdata->set_show_horizontal_spans_overlay(local_visualize_horizontal_spans);

        // Scope used for accessing and possibly editing the config settings
        {
            mi::base::Handle<const nv::index::ISession> session(
                dice_transaction->access<nv::index::ISession>(m_session_tag));
            mi::base::Handle<const nv::index::IConfig_settings> configs(
                dice_transaction->access<nv::index::IConfig_settings>(session->get_config()));

            nv::index::IConfig_settings::Volume_filtering_modes
                local_filtering_of_volume = configs->get_volume_filter();
            bool ray_segment_accum = configs->get_ray_segment_accumulation_technique();

            // Visual quality level first, if specified
            if (user_arguments.get_value("qualityLevel"))
            {
                get_quality_level(
                    user_arguments.get_value("qualityLevel"), local_filtering_of_volume, ray_segment_accum);
            }

            if (user_arguments.get_value("filter"))
            {
                const std::string filter  = user_arguments.get_value("filter");
                mi::Uint32 filter_id = get_sint32(filter);
                local_filtering_of_volume = (nv::index::IConfig_settings::Volume_filtering_modes)filter_id;
            }

            if (user_arguments.get_value("ray_segment_accum"))
            {
                ray_segment_accum = (on == user_arguments.get_value("ray_segment_accum"));
            }

            mi::Float32 local_step_size_min = get_float32(user_arguments.get_value("stepmin"));
            mi::Float32 local_step_size_max = get_float32(user_arguments.get_value("stepmax"));

            mi::Uint32 local_geometry_samples_max = 1;
            if (user_arguments.get_value("geometrySamplesMax"))
                local_geometry_samples_max = static_cast<mi::Uint32>(get_float32(user_arguments.get_value("geometrySamplesMax")));
            mi::Uint32 local_geometry_samples_min = 1;
            if (user_arguments.get_value("geometrySamplesMin"))
                local_geometry_samples_min = static_cast<mi::Uint32>(get_float32(user_arguments.get_value("geometrySamplesMin")));
            mi::Float32 local_geometry_samples_error = 0.f;
            if (user_arguments.get_value("geometrySamplesError"))
                local_geometry_samples_error = get_float32(user_arguments.get_value("geometrySamplesError"));
            mi::Uint32 local_rendering_samples_max = 1;
            if (user_arguments.get_value("renderingSamplesMax"))
                local_rendering_samples_max = static_cast<mi::Uint32>(get_float32(user_arguments.get_value("renderingSamplesMax")));

            appdata->set_current_colormap_index(get_sint32(user_arguments.get_value("colormap")));
            bool boost_geometry_colormap_opacity = (user_arguments.get_value("opacityboost") == on);

            if (
                (configs->get_rendering_mode()   != local_rendering_mode)                                       ||
                (configs->get_cpu_thread_count() != local_cpu_thread_count)                                     ||
                (configs->get_volume_filter() != local_filtering_of_volume)                                     ||
                (configs->get_step_size_min() != local_step_size_min)                                           ||
                (configs->get_step_size_max() != local_step_size_max)                                           ||
                (configs->get_ray_segment_accumulation_technique() != ray_segment_accum)                        ||
                (configs->get_geometry_samples_error() != local_geometry_samples_error)                         ||
                (configs->get_geometry_samples_min() != local_geometry_samples_min)                             ||
                (configs->get_geometry_samples_max() != local_geometry_samples_max)                             ||
                (configs->get_rendering_samples() != local_rendering_samples_max)                               ||
                (configs->get_compositing_mode() != local_compositing_mode)                                     ||
                (configs->is_parallel_rendering_and_compositing() != local_parallel_rendering_and_compositing)  ||
                (configs->get_size_of_rendering_results_in_queue() != local_size_of_rendering_result_queue)     ||
                (configs->is_monitor_performance_values() != local_monitor_performance_values)                  ||
                (configs->get_nb_spans() != local_nb_spans_per_host && local_nb_spans_per_host > 0)             ||
                (configs->is_boost_geometry_colormap_opacity() != boost_geometry_colormap_opacity)
                )
            {
                {   /// Retrieve a handle to the config settings for editing ....
                    mi::base::Handle<nv::index::IConfig_settings> edit_config_settings(
                        dice_transaction->edit<nv::index::IConfig_settings>(session->get_config()));

                    /// ... edit the config settigns
                    edit_config_settings->set_rendering_mode(local_rendering_mode);
                    edit_config_settings->set_cpu_thread_count(local_cpu_thread_count);
                    edit_config_settings->set_volume_filter(local_filtering_of_volume);
                    edit_config_settings->set_step_size_min(local_step_size_min);
                    edit_config_settings->set_step_size_max(local_step_size_max);
                    edit_config_settings->set_ray_segment_accumulation_technique(ray_segment_accum);
                    edit_config_settings->set_geometry_samples_error(local_geometry_samples_error);
                    edit_config_settings->set_geometry_samples_min(local_geometry_samples_min);
                    edit_config_settings->set_geometry_samples_max(local_geometry_samples_max);
                    edit_config_settings->set_rendering_samples(local_rendering_samples_max);
                    edit_config_settings->set_compositing_mode(local_compositing_mode);
                    edit_config_settings->set_parallel_rendering_and_compositing(local_parallel_rendering_and_compositing);
                    edit_config_settings->set_size_of_rendering_results_in_queue(local_size_of_rendering_result_queue);
                    set_monitor_performance_values(edit_config_settings, local_monitor_performance_values); // for sync state
                    if (local_nb_spans_per_host > 0)
                        edit_config_settings->set_nb_spans(local_nb_spans_per_host);
                    edit_config_settings->set_boost_geometry_colormap_opacity(boost_geometry_colormap_opacity);
                }   /// Leaving scope causes edited config settings to release and thus commit its changes to the database
            }
        }   /// Leaving scope causes accessed config setting to release

        dice_transaction->commit();
    }
    else if (strcmp(cmd, "set_current_scope") == 0)
    {
        const int scope_nb = get_uint32(user_arguments.get_value("scope"));
        m_irc.set_current_scope(scope_nb);
    }
    else if (strcmp(cmd, "timestepAnimation") == 0)
    {
        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            cur_scope->create_transaction<mi::neuraylib::IDice_transaction>());
        assert(dice_transaction.is_valid_interface());

        const std::string on = "on"; // To simplify comparison to C-style string
        mi::Uint32 local_current_timestep = 0;
        // Time step animation
        bool local_tsteps_playforward = false;
        bool local_tsteps_playbackward = false;
        bool local_tsteps_stepback = false;
        bool local_tsteps_stepforward = false;
        bool local_pause_tsteps = false;

        if (user_arguments.get_value("currenttimestep"))
        {
            local_current_timestep = get_uint32(user_arguments.get_value("currenttimestep"));
            local_tsteps_playforward = (user_arguments.get_value("playforwardtimestep") == on);
            local_tsteps_playbackward = (user_arguments.get_value("playbackwardtimestep") == on);
            local_tsteps_stepback = (user_arguments.get_value("stepbacktimestep") == on);
            local_tsteps_stepforward = (user_arguments.get_value("stepforwardtimestep") == on);
            local_pause_tsteps = (user_arguments.get_value("pausetimestep") == on);
        }
        m_irc.m_is_tsteps_playforward = local_tsteps_playforward;
        m_irc.m_is_tsteps_playbackward = local_tsteps_playbackward;
        m_irc.m_is_pause_tsteps = local_pause_tsteps;
        m_irc.m_is_tsteps_stepback = local_tsteps_stepback;
        m_irc.m_is_tsteps_stepforward = local_tsteps_stepforward;

        {
            mi::base::Handle<const nv::index::ISession> session(
                dice_transaction->access<nv::index::ISession>(m_session_tag));
            mi::base::Handle<const nv::index::IConfig_settings> configs(
                dice_transaction->access<nv::index::IConfig_settings>(session->get_config()));
            if (configs->get_current_timestep() != local_current_timestep)
            {
                /// Retrieve a handle to the config settings for editing ....
                mi::base::Handle<nv::index::IConfig_settings> edit_config_settings(
                    dice_transaction->edit<nv::index::IConfig_settings>(session->get_config()));

                /// ... edit the config settigns
                edit_config_settings->set_current_timestep(local_current_timestep);
            }
        }
        dice_transaction->commit();
    }
    else if (strcmp(cmd, "getslices") == 0)
    {
        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            cur_scope->create_transaction<mi::neuraylib::IDice_transaction>());
        assert(dice_transaction.is_valid_interface());

        const std::vector<mi::neuraylib::Tag> volume_tag_vec =
            appdata->get_volume_tag_vec();
        const mi::Sint32 volume_id = get_sint32(user_arguments.get_value("volume_id"));

        std::string names;
        if((volume_id >= 0) && (volume_id < static_cast<mi::Sint32>(volume_tag_vec.size())))
        {
            mi::base::Handle<const nv::index::IRegular_volume> volume(
                dice_transaction->access<const nv::index::IRegular_volume>(
                    volume_tag_vec.at(volume_id)));
            if (volume.is_valid_interface())
            {
                const mi::Uint32 nb_slices = volume->get_nb_slices();
                std::ostringstream os;
                for(mi::Uint32 i=0; i<nb_slices; ++i)
                {
                    const mi::neuraylib::Tag slice_tag = volume->get_slice(i);
                    mi::base::Handle<const nv::index::ISlice_scene_element> slice(
                        dice_transaction->access<const nv::index::ISlice_scene_element>(slice_tag));

                    mi::base::Handle<const nv::index::ISection_scene_element> section_slice(
                        slice->get_interface<const nv::index::ISection_scene_element>());
                    if(section_slice.get())
                    {
                        const nv::index::ISection_scene_element::Slice_orientation orientation =
                            section_slice->get_orientation();
                        if(orientation==nv::index::ISection_scene_element::INLINE_SECTION)
                            os << "Inline section";
                        else if(orientation==nv::index::ISection_scene_element::CROSS_LINE_SECTION)
                            os << "Cross-line section";
                        else if(orientation==nv::index::ISection_scene_element::HORIZONTAL_SECTION)
                            os << "Horizontal section";
                    }
                    else
                    {
                        os << "Vertical profile";
                    }

                    if (i < nb_slices - 1)
                        os << "|";

                    names = os.str();
                }
            }
        }

        // Handle interaction group
        const nv::index_common::String_dict* dict = appdata->peek_app_proj();
        const std::string interact = dict->get("app::interaction_group::group");
        if (!interact.empty())
        {
            if (!names.empty())
                names += "|";
            names += interact;
        }

        if (!names.empty())
            response_arguments.set_value("names", names.c_str());

        dice_transaction->commit();
    }
    else if (strcmp(cmd, "getslice") == 0)
    {
        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            cur_scope->create_transaction<mi::neuraylib::IDice_transaction>());
        assert(dice_transaction.is_valid_interface());

        const std::vector<mi::neuraylib::Tag> volume_tag_vec =
            appdata->get_volume_tag_vec();
        const mi::Sint32 volume_id = get_sint32(user_arguments.get_value("volume_id"));

        bool done = false;
        if ((volume_id >= 0) && (volume_id < static_cast<mi::Sint32>(volume_tag_vec.size())))
        {
            mi::base::Handle<const nv::index::IRegular_volume> volume(
                dice_transaction->access<const nv::index::IRegular_volume>(
                    volume_tag_vec.at(volume_id)));
            assert(volume.is_valid_interface());

            const mi::Sint32 slice_id = get_sint32(user_arguments.get_value("slice_id"));
            if (slice_id < static_cast<int>(volume->get_nb_slices()))
            {
                done = true;
                const mi::neuraylib::Tag slice_tag = volume->get_slice(slice_id);
                mi::base::Handle<const nv::index::ISlice_scene_element> slice(
                    dice_transaction->access<const nv::index::ISlice_scene_element>(slice_tag));
                assert(slice.is_valid_interface());

                mi::base::Handle<const nv::index::ISection_scene_element> section_slice(
                    slice->get_interface<const nv::index::ISection_scene_element>());
                if(section_slice.get())
                {
                    mi::math::Bbox<mi::Uint32, 3> box =
                        volume->get_IJK_bounding_box();
                    box.max -= mi::math::Vector<mi::Uint32, 3>(1); // Adapt for half voxels on boundary

                    response_arguments.set_value("section", "1");

                    {
                        std::ostringstream os;
                        os << section_slice->get_position();
                        response_arguments.set_value("pos", os.str().c_str());
                    }

                    const nv::index::ISection_scene_element::Slice_orientation orientation =
                        section_slice->get_orientation();
                    mi::Uint32 slice_min = 0;
                    mi::Uint32 slice_max = 0;
                    if(orientation==nv::index::ISection_scene_element::INLINE_SECTION)
                    { slice_min = mi::Uint32(box.min.x); slice_max = mi::Uint32(box.max.x);}
                    else if(orientation==nv::index::ISection_scene_element::CROSS_LINE_SECTION)
                    { slice_min = mi::Uint32(box.min.y); slice_max = mi::Uint32(box.max.y);}
                    else if(orientation==nv::index::ISection_scene_element::HORIZONTAL_SECTION)
                    { slice_min = mi::Uint32(box.min.z); slice_max = mi::Uint32(box.max.z);}
                    {
                        std::ostringstream os;
                        os << slice_min;
                        response_arguments.set_value("range_min", os.str().c_str());
                    }
                    {
                        std::ostringstream os;
                        os << slice_max;
                        response_arguments.set_value("range_max", os.str().c_str());
                    }
                }
                else
                {
                    response_arguments.set_value("section", "0");
                    response_arguments.set_value("pos", "0.0");
                }

                {
                    const mi::neuraylib::Tag colormap_tag = slice->assigned_colormap();
                    const mi::Uint32 nb_colormaps = get_number_of_colormap();
                    mi::Uint32 colormap_index = 0;
                    for(mi::Uint32 i=0; i<nb_colormaps; ++i)
                    {
                        if(get_colormap_tag(i)==colormap_tag)
                        {
                            colormap_index = i;
                        }
                    }

                    std::ostringstream os;
                    os << colormap_index;
                    response_arguments.set_value("slice_cm_id", os.str().c_str());
                }

                const bool rendered_enable = slice->get_enabled();
                if(rendered_enable)
                    response_arguments.set_value("enabled", "1");
                else
                    response_arguments.set_value("enabled", "0");
            }
        }

        // Handle interaction group
        const nv::index_common::String_dict* dict = appdata->peek_app_proj();
        const std::string interact = dict->get("app::interaction_group::group");
        if (!done && !interact.empty())
        {
            response_arguments.set_value("section", "-1");
            response_arguments.set_value(
                "range_min", dict->get("app::interaction_group::min", "-1000").c_str());
            response_arguments.set_value(
                "range_max", dict->get("app::interaction_group::max", "1000").c_str());
        }

        dice_transaction->commit();
    }
    else if (strcmp(cmd, "setslicechanges") == 0)
    {
        mi::Sint32 volume_id = get_sint32(user_arguments.get_value("volume_id"));
        mi::Sint32 slice_id  = get_sint32(user_arguments.get_value("slice_id"));
        mi::Uint32  const current_volume_idx   = volume_id;
        mi::Float32 const slice_position       = get_float32_userarg(user_arguments, "slice_pos");
        mi::Uint32  const slice_colormap_index = get_uint32_userarg(user_arguments, "slice_colormap");
        bool const is_slice_enabled            = get_sint32_userarg(user_arguments, "enabled") == 1;

        mi::neuraylib::Tag slice_colormap_tag;
        if (slice_colormap_index < get_number_of_colormap())
        {
            slice_colormap_tag = get_colormap_tag(slice_colormap_index);
        }
        set_slice_change(cur_scope.get(), current_volume_idx,
                         slice_id, slice_position, slice_colormap_tag, is_slice_enabled);
    }
    else if (strcmp(cmd, "getsteering") == 0)
    {
        std::vector< std::string > steering_name_vec;
        Nvindex_AppData::instance()->get_steering_names(steering_name_vec);
        const std::string steering_name_gui = concat_string_with_separator(steering_name_vec, "|");
        response_arguments.set_value("name", steering_name_gui.c_str());
    }
    else if (strcmp(cmd, "setsteering") == 0)
    {
        const mi::Sint32 steering_id = get_sint32(user_arguments.get_value("steering_id"));
        const std::string steering_name = user_arguments.get_value("steering_name");

        Nvindex_AppData::instance()->set_steering(steering_id, steering_name);
    }
    else if (strcmp(cmd, "getvolumes") == 0)
    {
        std::vector< std::string > volume_name_vec;
        Nvindex_AppData::instance()->get_volume_name_vec(cur_scope, volume_name_vec);
        const std::string volume_name_gui = concat_string_with_separator(volume_name_vec, "|");
        response_arguments.set_value("names", volume_name_gui.c_str());
    }
    else if (strcmp(cmd, "getvolume") == 0)
    {
        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            cur_scope->create_transaction<mi::neuraylib::IDice_transaction>());
        assert(dice_transaction.is_valid_interface());

        const std::vector<mi::neuraylib::Tag> volume_tag_vec =
            appdata->get_volume_tag_vec();
        const mi::Sint32 volume_id = get_sint32(user_arguments.get_value("volume_id"));

        if ((volume_id >= 0) && (volume_id < static_cast<mi::Sint32>(volume_tag_vec.size())))
        {
            mi::base::Handle<const nv::index::IRegular_volume> volume(
                dice_transaction->access<nv::index::IRegular_volume>(
                    volume_tag_vec.at(volume_id)));

            if (volume.is_valid_interface())
            {
                response_arguments.set_value("enabled", volume->get_enabled() ? "1" : "0");

                const mi::neuraylib::Tag colormap_tag = volume->assigned_colormap();
                mi::Uint32 colormap_index = 0;
                for (mi::Uint32 i=0; i < get_number_of_colormap(); ++i) {
                    if (get_colormap_tag(i) == colormap_tag)
                        colormap_index = i;
                }

                std::ostringstream os;
                os << colormap_index;
                response_arguments.set_value("colormap_id", os.str().c_str());
            }
        }
        dice_transaction->commit();
    }
    else if (strcmp(cmd, "setvolume") == 0)
    {
        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            cur_scope->create_transaction<mi::neuraylib::IDice_transaction>());
        assert(dice_transaction.is_valid_interface());

        const std::vector<mi::neuraylib::Tag> volume_tag_vec = appdata->get_volume_tag_vec();
        const mi::Sint32 volume_id = get_sint32(user_arguments.get_value("volume_id"));

        if ((volume_id >= 0) && (volume_id < static_cast<mi::Sint32>(volume_tag_vec.size())))
        {
            mi::base::Handle<nv::index::IRegular_volume> volume(
                dice_transaction->edit<nv::index::IRegular_volume>(
                    volume_tag_vec.at(volume_id)));

            if (volume)
            {
                const std::string enabled = user_arguments.get_value("enabled");
                volume->set_enabled(enabled == "1");

                const std::string colormap = user_arguments.get_value("colormap_id");
                if (get_sint32(colormap) < static_cast<mi::Sint32>(get_number_of_colormap()))
                {
                    mi::neuraylib::Tag colormap_tag = get_colormap_tag(get_sint32(colormap));
                    volume->assign_colormap(colormap_tag);
                }
            }
        }
        dice_transaction->commit();
    }
    else if (strcmp(cmd, "gethorizons") == 0)
    {
        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            cur_scope->create_transaction<mi::neuraylib::IDice_transaction>());
        assert(dice_transaction.is_valid_interface());

        {
            std::ostringstream os;

            const std::vector<mi::neuraylib::Tag> hf_tag_vec = appdata->get_heightfield_tag_vec();
            const mi::Uint32 nb_heightfield = hf_tag_vec.size();

            for (mi::Uint32 i = 0; i < nb_heightfield; ++i)
            {
                mi::base::Handle<const nv::index::IRegular_heightfield> heightfield(
                    dice_transaction->access<const nv::index::IRegular_heightfield>(hf_tag_vec.at(i)));
                assert(heightfield.is_valid_interface());

                os << heightfield->get_name();
                if (i < nb_heightfield-1)
                {
                    os << "|";
                }
            }

            response_arguments.set_value("names", os.str().c_str());
        }

        dice_transaction->commit();
    }
    else if (strcmp(cmd, "gethorizon") == 0)
    {
        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            cur_scope->create_transaction<mi::neuraylib::IDice_transaction>());
        assert(dice_transaction.is_valid_interface());

        const std::vector<mi::neuraylib::Tag> hf_tag_vec = appdata->get_heightfield_tag_vec();
        const mi::Sint32 heightfield_id = get_sint32(user_arguments.get_value("horizon_id"));
        if ((heightfield_id >= 0) && (heightfield_id < static_cast<mi::Sint32>(hf_tag_vec.size())))
        {
            const mi::neuraylib::Tag heightfield_tag = hf_tag_vec.at(heightfield_id);
            mi::base::Handle<const nv::index::IRegular_heightfield> heightfield(
                dice_transaction->access<const nv::index::IRegular_heightfield>(heightfield_tag));

            if (heightfield.is_valid_interface())
            {
                if(heightfield->get_enabled())
                    response_arguments.set_value("enabled", "1");
                else
                    response_arguments.set_value("enabled", "0");

                if(heightfield->get_colormap_mapping())
                    response_arguments.set_value("colormap_mapping_enabled", "1");
                else
                    response_arguments.set_value("colormap_mapping_enabled", "0");
                {
                    // colormap GUI is effective all the time
                    const mi::neuraylib::Tag colormap_tag = heightfield->assigned_colormap();
                    const mi::Uint32 nb_colormaps = get_number_of_colormap();
                    mi::Uint32 colormap_index = 0;
                    for(mi::Uint32 i=0; i<nb_colormaps; ++i){
                        if(get_colormap_tag(i) == colormap_tag)
                        {
                            colormap_index = i;
                        }
                    }

                    std::ostringstream os;
                    os << colormap_index;
                    response_arguments.set_value("horizon_colormap_id", os.str().c_str());
                }

                const char* given_name = dice_transaction->tag_to_name(heightfield_tag);
                if(given_name!=NULL)
                {
                    const std::string material_tag_name = "mat_" + std::string(given_name);
                    const mi::neuraylib::Tag material_tag = dice_transaction->name_to_tag(material_tag_name.c_str());

                    if (material_tag.is_valid())
                    {
                        mi::base::Handle<const nv::index::IPhong_gl> material(
                            dice_transaction->access<nv::index::IPhong_gl>(material_tag));
                        if(material.is_valid_interface())
                        {
                            mi::math::Color new_color;
                            new_color = material->get_ambient();
                            new_color /= 0.4f; // Normalize, assuming 0.4 ambient and diffuse term
                            new_color.a = material->get_opacity();

                            {
                                std::ostringstream os;
                                os << new_color.r;
                                response_arguments.set_value("color_r", os.str().c_str());
                            }
                            {
                                std::ostringstream os;
                                os << new_color.g;
                                response_arguments.set_value("color_g", os.str().c_str());
                            }
                            {
                                std::ostringstream os;
                                os << new_color.b;
                                response_arguments.set_value("color_b", os.str().c_str());
                            }
                            {
                                std::ostringstream os;
                                os << new_color.a;
                                response_arguments.set_value("color_a", os.str().c_str());
                            }
                        }
                    }
                }
            }
        }
        dice_transaction->commit();
    }
    else if (strcmp(cmd, "sethorizon") == 0)
    {
        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            cur_scope->create_transaction<mi::neuraylib::IDice_transaction>());
        assert(dice_transaction.is_valid_interface());

        const std::vector<mi::neuraylib::Tag> hf_tag_vec = appdata->get_heightfield_tag_vec();
        const mi::Uint32 heightfield_id = get_sint32(user_arguments.get_value("horizon_id"));

        if(heightfield_id < hf_tag_vec.size())
        {
            const mi::neuraylib::Tag heightfield_tag = hf_tag_vec.at(heightfield_id);
            assert(heightfield_tag.is_valid());

            mi::base::Handle<nv::index::IRegular_heightfield> heightfield(
                dice_transaction->edit<nv::index::IRegular_heightfield>(heightfield_tag));

            if(heightfield.is_valid_interface())
            {
                const std::string enabled = user_arguments.get_value("enabled");
                const mi::Uint32 enabled_int = get_uint32(enabled);
                if(enabled_int==1)
                {
                    heightfield->set_enabled(true);
                }
                else
                {
                    heightfield->set_enabled(false);
                }

                const std::string colormap_enabled = user_arguments.get_value("colormap_mapping_enabled");
                const mi::Uint32 colormap_enabled_int = get_uint32(colormap_enabled);
                if(colormap_enabled_int==1)
                {
                    heightfield->set_colormap_mapping(true);
                }
                else
                {
                    heightfield->set_colormap_mapping(false);
                }

                const std::string heightfield_colormap = user_arguments.get_value("horizon_colormap_id");
                const mi::Uint32 colormap_index = get_uint32(heightfield_colormap);

                if (colormap_index < get_number_of_colormap())
                {
                    const mi::neuraylib::Tag colormap_tag = get_colormap_tag(colormap_index);
                    heightfield->assign_colormap(colormap_tag);  // Assign the colormap to the heightfield
                }

                std::string color_r;
                if (user_arguments.get_value("color_r"))
                    color_r = user_arguments.get_value("color_r");
                std::string color_g;
                if (user_arguments.get_value("color_g"))
                    color_g = user_arguments.get_value("color_g");
                std::string color_b;
                if (user_arguments.get_value("color_b"))
                    color_b = user_arguments.get_value("color_b");

                const std::string color_a = user_arguments.get_value("color_a"); // always there
                const bool set_color = !color_r.empty();

                mi::math::Color new_color(0.f);
                if (set_color)
                {
                    new_color.r = get_float32(color_r);
                    new_color.g = get_float32(color_g);
                    new_color.b = get_float32(color_b);
                }
                new_color.a = get_float32(color_a);

                const char* given_name = dice_transaction->tag_to_name(heightfield_tag);
                if(given_name!=NULL)
                {
                    const std::string material_tag_name = "mat_" + std::string(given_name);
                    const mi::neuraylib::Tag material_tag = dice_transaction->name_to_tag(material_tag_name.c_str());

                    if(material_tag.is_valid())
                    {
                        mi::base::Handle<nv::index::IPhong_gl> material(
                            dice_transaction->edit<nv::index::IPhong_gl>(material_tag));
                        if(material.is_valid_interface())
                        {
                            if (set_color)
                            {
                                material->set_ambient(new_color * 0.4f); // Assume 0.4 ambient and diffuse
                                material->set_diffuse(new_color * 0.4f);
                            }
                            material->set_opacity(new_color.a);
                        }
                    }
                    else
                    {
                        WARN_LOG << "Material attribute '" << material_tag_name << "' not found for heightfield "
                                 << "'" << given_name << "', color change is not possible.";
                    }
                }
            }
        }
        dice_transaction->commit();
    }
    else if (strcmp(cmd, "loadSeedLines") == 0)
    {
        const std::string file_name  = user_arguments.get_value("filename");

        // Do strict filename checking for some added security:
        // a) File must be located in the current working directory
        // b) It must have ".asc" extension
        // c) It may only contain alphanumeric characters, '-' or '_'.
        bool valid_file_name = true;
        if ((file_name.find(".asc") != file_name.size() - 4) || file_name.size() <= 4)
        {
            valid_file_name = false;
        }
        else
        {
            for (size_t i=0; i < file_name.size() - 4; ++i)
            {
                char c = file_name[i];
                if (!isalpha(c) && !isdigit(c) && c != '-' && c != '_')
                    valid_file_name = false;
            }
        }

        if (valid_file_name)
        {
            const std::string heightfield_id_str = user_arguments.get_value("horizon_id");
            const mi::Uint32  heightfield_id = get_uint32(heightfield_id_str);

            // Import the seed lines from the given file
            import_seed_lines(cur_scope, m_session_tag, heightfield_id, file_name);
        }
        else
        {
            ERROR_LOG << "Invalid filename specified. Must only contain alphanumeric characters, "
                "'-' or '_', and must have '.asc' extension.";
        }
    }
    else if (strcmp(cmd, "removeArbirarySeedLine") == 0)
    {
        const std::string heightfield_id_str = user_arguments.get_value("horizon_id");
        const mi::Uint32  heightfield_id     = get_uint32(heightfield_id_str);
        remove_arbitrary_seed_line(cur_scope, m_session_tag, heightfield_id);
    }
    else if (strcmp(cmd, "removeArbirarySeedPoint") == 0)
    {
        const std::string heightfield_id_str = user_arguments.get_value("horizon_id");
        const mi::Uint32  heightfield_id     = get_uint32(heightfield_id_str);
        remove_arbitrary_seed_point(cur_scope, m_session_tag, heightfield_id);
    }
    else if (strcmp(cmd, "getcurrentcolormap") == 0)
    {
        get_user_current_colormap(response_arguments);
    }
    else if (strcmp(cmd, "getcolormaps") == 0)
    {
        get_user_current_colormap_all(response_arguments);
    }
    else if (strcmp(cmd, "setcolormap") == 0)
    {
        set_user_colormap(user_arguments);
    }
    else if (strcmp(cmd, "localizeColormap") == 0)
    {
        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            get_scope()->create_transaction<mi::neuraylib::IDice_transaction>());

        const std::string cm_id = user_arguments.get_value("id");
        mi::Uint32 colormap_id = nv::index_common::get_uint32(cm_id);

        dice_transaction->localize(
            get_colormap_tag(colormap_id),
            mi::neuraylib::IDice_transaction::LOCAL_SCOPE);

        dice_transaction->commit();
    }
    else if (strcmp(cmd, "getvolumeofinterest") == 0)
    {
        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            cur_scope->create_transaction<mi::neuraylib::IDice_transaction>());
        assert(dice_transaction.is_valid_interface());

        mi::math::Bbox<mi::Float32, 3> extent;
        mi::math::Bbox<mi::Float32, 3> roi;

        mi::Sint32 scene_element_id = get_sint32(user_arguments.get_value("scene_element"));
        if (scene_element_id > 0)
        {
            // User selected a volume or a heightfield
            std::vector<mi::neuraylib::Tag> volume_tag_vec = appdata->get_volume_tag_vec();
            const mi::Sint32 nb_volumes = static_cast<mi::Sint32>(volume_tag_vec.size());
            if (scene_element_id <= nb_volumes)
            {
                // Volumes come first
                mi::Sint32 volume_id = scene_element_id - 1; // the 0-th one is the entire scene.
                assert(volume_id >= 0);
                assert(volume_id < static_cast<mi::Sint32>(volume_tag_vec.size()));
                mi::base::Handle<const nv::index::IRegular_volume> volume(
                    dice_transaction->access<const nv::index::IRegular_volume>(
                        volume_tag_vec.at(volume_id)));
                assert(volume.is_valid_interface());

                const mi::math::Bbox<mi::Uint32, 3> bbox = volume->get_IJK_bounding_box();
                extent = mi::math::Bbox<mi::Float32, 3>(bbox);

                roi = volume->get_IJK_region_of_interest();
            }
            else
            {
                std::vector<mi::neuraylib::Tag> heightfield_tag_vec =
                    appdata->get_heightfield_tag_vec();
                mi::Sint32 heightfield_id = scene_element_id - nb_volumes - 1;
                if (heightfield_id < static_cast<mi::Sint32>(heightfield_tag_vec.size()))
                {
                    mi::base::Handle<const nv::index::IRegular_heightfield> heightfield(
                        dice_transaction->access<const nv::index::IRegular_heightfield>(
                            heightfield_tag_vec.at(heightfield_id)));

                    if (heightfield.is_valid_interface())
                    {
                        extent = heightfield->get_IJK_bounding_box();
                        roi = heightfield->get_IJK_region_of_interest();
                    }
                }
            }
        }
        else
        {
            // Use extent and region of interest of the entire scene
            //
            // The extent is calculated as the union of the maximum size of all scene objects
            // and the initial region of interest specified in the project file.
            mi::base::Handle<const nv::index::ISession> session(
                dice_transaction->access<const nv::index::ISession>(m_session_tag));
            assert(session.is_valid_interface());
            mi::base::Handle<const nv::index::IScene> scene(
                dice_transaction->access<const nv::index::IScene>(session->get_scene()));
            assert(scene.is_valid_interface());

            extent = get_scene_bounding_box(m_session_tag, dice_transaction.get());
            extent.push_back(appdata->m_initial_roi);

            // Limit the extent if requested in the configuration, to prevent users from setting
            // a region of interest that is too large.
            String_dict* app_proj = appdata->peek_app_proj();
            assert(app_proj != 0);
            if (app_proj->is_defined("index::region_of_interest_limit"))
            {
                mi::math::Bbox<mi::Float32, 3> limit(
                    get_bbox_float32_3(app_proj->get("index::region_of_interest_limit")));
                // ROI is specified as "value range" in project file, so convert it to the bounding box format
                // that is used internally.
                limit.max += mi::math::Vector<mi::Float32, 3>(1.f);
                extent = mi::math::clip(extent, limit); // Compute intersection
            }

            roi = get_XYZ_global_region_of_interest_bbox(session.get(), dice_transaction.get());

            // These are actually not using "value range": add 1 here so that the conversion
            // below will result in the original value again
            extent.max += mi::math::Vector<mi::Float32, 3>(1.f);
            roi.max += mi::math::Vector<mi::Float32, 3>(1.f);
        }

        // Return total extent (used for min/max values of the sliders)
        {
            // Flash viewer uses "value range" format, so we must convert the bounding box here.
            std::ostringstream os;
            os << floorf(extent.min.x) << "," << ceilf(extent.max.x) - 1 << ","
               << floorf(extent.min.y) << "," << ceilf(extent.max.y) - 1 << ","
               << floorf(extent.min.z) << "," << ceilf(extent.max.z) - 1;
            response_arguments.set_value("total_extent", os.str().c_str());
        }

        // Return region of interest
        {
            // Flash viewer uses "value range" format, so we must convert the bounding box here.
            std::ostringstream os;
            os << (mi::Sint32)(roi.min.x)     << ",";
            os << (mi::Sint32)(roi.max.x) - 1 << ",";
            os << (mi::Sint32)(roi.min.y)     << ",";
            os << (mi::Sint32)(roi.max.y) - 1 << ",";
            os << (mi::Sint32)(roi.min.z)     << ",";
            os << (mi::Sint32)(roi.max.z) - 1;
            response_arguments.set_value("interest", os.str().c_str());
        }

        dice_transaction->commit();

    }
    else if (strcmp(cmd, "setvolumeofinterest") == 0)
    {
        mi::base::Lock::Block block(&appdata->m_scene_edit_lock); // Access lock before starting transaction

        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            cur_scope->create_transaction<mi::neuraylib::IDice_transaction>());
        assert(dice_transaction.is_valid_interface());

        mi::base::Handle<const nv::index::ISession> session(
            dice_transaction->access<const nv::index::ISession>(m_session_tag));

        mi::math::Bbox<mi::Float32, 3> interest_bbox;
        {
            const std::string min_x_in = user_arguments.get_value("min_x");
            const mi::Sint32 min_x_value = get_sint32(min_x_in);

            const std::string max_x_in = user_arguments.get_value("max_x");
            const mi::Sint32 max_x_value = get_sint32(max_x_in);

            const std::string min_y_in = user_arguments.get_value("min_y");
            const mi::Sint32 min_y_value = get_sint32(min_y_in);

            const std::string max_y_in = user_arguments.get_value("max_y");
            const mi::Sint32 max_y_value = get_sint32(max_y_in);

            const std::string min_z_in = user_arguments.get_value("min_z");
            const mi::Sint32 min_z_value = get_sint32(min_z_in);

            const std::string max_z_in = user_arguments.get_value("max_z");
            const mi::Sint32 max_z_value = get_sint32(max_z_in);

            // Flash viewer uses "value range" format, so we must convert it to a bounding box here.
            interest_bbox.min.x = static_cast<mi::Float32>(min_x_value);
            interest_bbox.max.x = static_cast<mi::Float32>(max_x_value + 1);
            interest_bbox.min.y = static_cast<mi::Float32>(min_y_value);
            interest_bbox.max.y = static_cast<mi::Float32>(max_y_value + 1);
            interest_bbox.min.z = static_cast<mi::Float32>(min_z_value);
            interest_bbox.max.z = static_cast<mi::Float32>(max_z_value + 1);
        }

        const mi::Sint32 scene_element_id = get_sint32(user_arguments.get_value("scene_element"));
        if (scene_element_id > 0)
        {
            const std::vector<mi::neuraylib::Tag> volume_tag_vec = appdata->get_volume_tag_vec();
            const mi::Sint32 nb_volumes = static_cast<mi::Sint32>(volume_tag_vec.size());
            if (scene_element_id <= nb_volumes)
            {
                // Volumes come first
                mi::Sint32 volume_id = scene_element_id - 1;
                mi::base::Handle<nv::index::IRegular_volume> volume(
                    dice_transaction->edit<nv::index::IRegular_volume>(volume_tag_vec.at(volume_id)));

                if (volume.is_valid_interface())
                    volume->set_IJK_region_of_interest(interest_bbox);
            }
            else
            {
                const std::vector<mi::neuraylib::Tag> hf_tag_vec = appdata->get_heightfield_tag_vec();
                const mi::Sint32 heightfield_id = scene_element_id - nb_volumes - 1;
                if (heightfield_id < static_cast<mi::Sint32>(hf_tag_vec.size()))
                {
                    mi::base::Handle<nv::index::IRegular_heightfield> heightfield(
                        dice_transaction->edit<nv::index::IRegular_heightfield>(hf_tag_vec.at(heightfield_id)));

                    if (heightfield.is_valid_interface())
                        heightfield->set_IJK_region_of_interest(interest_bbox);
                }
            }
        }
        else
        {
            // Set global region of interest for the entire scene
            interest_bbox.max -= mi::math::Vector<mi::Float32, 3>(1.f); // Not using "value range"
            set_XYZ_global_region_of_interest_bbox(interest_bbox, session.get(), dice_transaction.get());
        }

        if(appdata->is_show_region_of_interest_overlay())
        {
            // We don't update the shapes but delete and rebuild the shapes.
            // This is not a really clever way of handling roi changes but ...
            // ... sufficient for now.
            //
            // TODO: implement an updated methods and replace the quick and dirty solution
            disable_region_of_interest_annotation(dice_transaction, m_session_tag);
            enable_region_of_interest_annotation(dice_transaction, m_session_tag);
        }

        dice_transaction->commit();
    }
    else if (strcmp(cmd, "stats") == 0)
    {
        // Generic performance values
        {
            mi::Sint32 host = get_sint32(user_arguments.get_value("host")); // 0 means show global values
            if (host < 0)
                host = 0;
            mi::Uint64 nb_active_hosts = 0;

            mi::base::Lock::Block block(&appdata->m_performance_values_lock);
            if (appdata->m_performance_values) // This is not available until the first frame is rendered.
            {
                mi::Uint32 nb_types = appdata->m_performance_values->get_nb_type_names();
                for (mi::Uint32 i = 0; i < nb_types; ++i)
                {
                    const std::string type_name = appdata->m_performance_values->get_type_name(i);

                    // For these we always want to use the global value even though per-host values
                    // would also be available
                    bool force_global = (type_name == "nb_fragments");

                    mi::Uint64 value = appdata->m_performance_values->get(
                        type_name.c_str(), force_global ? 0 : host);

                    std::ostringstream os;
                    if (type_name.find("time_") == 0 || type_name == "frames_per_second")
                        // Convert raw time into miliseconds
                        os << value / static_cast<mi::Float32>(nv::index::IPerformance_values::TIME_RESOLUTION);
                    else
                        os << value;

                    // Use type name as identifier with "perf_" prefix added
                    response_arguments.set_value(std::string("perf_" + type_name).c_str(), os.str().c_str());
                    // Send all the key value, but GUI chooses what to show.
                    // DEBUG_LOG << "DEBUG: stat string: " << std::string("perf_" + type_name) << ":" << os.str();
                }

                mi::Float32 const fps_max = appdata->get_stat_recent_n_max_fps();
                mi::Float32 const fps_avg = appdata->get_stat_recent_n_ave_fps();

                {
                    std::ostringstream os;
                    os << fps_max;
                    response_arguments.set_value("perf_frames_per_second_max", os.str().c_str());
                }
                {
                    std::ostringstream os;
                    os << fps_avg;
                    response_arguments.set_value("perf_frames_per_second_avg", os.str().c_str());
                }

                if (appdata->m_application_compute_time > 0.f)
                {
                    std::ostringstream os;
                    os << (1.f / appdata->m_application_compute_time);
                    response_arguments.set_value("perf_compute_steps_per_second", os.str().c_str());
                }

                if (appdata->m_application_compute_memory > 0)
                {
                    std::ostringstream os;
                    os << appdata->m_application_compute_memory;
                    response_arguments.set_value("perf_compute_memory", os.str().c_str());
                }

                mi::Uint32 const system_host_id = 0;
                mi::Uint32 const auto_hspan_number = static_cast<mi::Uint32>(
                    appdata->m_performance_values->get("nb_horizontal_spans", system_host_id));
                {
                    std::ostringstream os;
                    os << auto_hspan_number;
                    response_arguments.set_value("auto_hspan_number", os.str().c_str());
                }

                // get number of active hosts
                nb_active_hosts = appdata->m_performance_values->get("nb_active_hosts", system_host_id);
            }

            std::ostringstream os;
            os << nb_active_hosts;
            response_arguments.set_value("nb_hosts", os.str().c_str());

            // Hostnames
            const std::vector<std::string>& hostnames =
                appdata->peek_host_info()->get_sorted_host_name_vec();

            std::string hostname_string;
            for (size_t i = 0; i < hostnames.size(); ++i)
            {
                std::string hostname = hostnames[i];

                // Only use the part before the first dot (if any)
                std::string::size_type p = hostname.find(".");
                if (p != std::string::npos)
                    hostname = hostname.substr(0, p - 1);

                hostname_string += hostname + "|";
            }
            response_arguments.set_value("hostnames", hostname_string.c_str());
        }

        // Progress
        {
            const mi::Float32 progress = m_progress_callback->get_progress();
            std::ostringstream os;
            os << progress;
            response_arguments.set_value("progress", os.str().c_str());
        }

        // Dataset information
        {
            String_dict* p_app_proj = appdata->peek_app_proj();
            assert(p_app_proj != 0);
            const std::string dataset_related_information = p_app_proj->get("app::dataset_note", "");
            response_arguments.set_value("dataset_note", dataset_related_information.c_str());
        }

        // Timestep slider
        {
            mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
                cur_scope->create_transaction<mi::neuraylib::IDice_transaction>());
            assert(dice_transaction.is_valid_interface());
            mi::base::Handle<const nv::index::ISession> session(
                dice_transaction->access<const nv::index::ISession>(
                    m_session_tag));
            mi::base::Handle<const nv::index::IConfig_settings> configs(
                dice_transaction->access<nv::index::IConfig_settings>(
                    session->get_config()));
            std::ostringstream os;
            os << configs->get_current_timestep();
            response_arguments.set_value("current_timestep", os.str().c_str());

            dice_transaction->commit();
        }
    }
    else if (strcmp(cmd, "user_command") == 0 &&
             get_bool(appdata->peek_app_proj()->get("app::demo", "no")))
    {
        // Disallow user commands in demo mode
        const std::string msg = "User commands are disabled in demo mode (app::demo = yes).";
        response_arguments.set_value("command_output", msg.c_str());
    }
    else if (strcmp(cmd, "user_command") == 0)
    {
        std::string command = user_arguments.get_value("user_command");
        INFO_LOG << "Received user command: '" << command << "'";

        // Define some shortcuts
        if (command == "screen")
            command = "set dice::rtmp_video_streaming::video_codec = screen video";
        else if (command == "h264")
            command = "set dice::rtmp_video_streaming::video_codec = h264";
        else if (command == "nvenc")
            command = "set dice::rtmp_video_streaming::video_codec = h264_nvenc";
        else if (command.find("bit ") == 0)
        {
            // Parameter is bitrate in Mbit/s
            std::istringstream is(command.substr(4));
            mi::Float32 bitrate;
            is >> bitrate;
            std::ostringstream os;
            os << "set dice::rtmp_video_streaming::video_bitrate = " << std::setprecision(15)
               << (bitrate * 1000000); // convert to bit/s
            command = os.str();
        }

        const std::string usage(
            "Available commands:\n"
            "animate_camera - camera animation\n"
            "cancel_rendering - cancel the current rendering of a frame\n"
            "cmp - set compression parameters for cluster compositing\n"
            "copy_volume <src_volume_tag> <dst_group_tag> <x> <y> <z>\n"
            "edit_volume <volume_tag> <value> - Set given value for all voxels inside ROI of volume\n"
            "export - write scene export to stdout\n"
            "filter_volume <volume_index> <filter_type> - Apply a filter (eg. dip_median)\n"
            "help - Show this command list\n"
            "hm_map <horizon_idx> <volume_idx> - generate heightfield amplitude value map\n"
            "perf - measure performance, help with '?' arg\n"
            "perf_logger - switch performance logger, print help with '?' arg\n"
            "quit - Terminate the server application\n"
            "save_colormap <file.cmap> - write current colormap to file\n"
            "set - overwrite project file config value\n"
            "set_image_file_resolution <width> <height> - set image file canvas resolution\n"
            "set_resolution <width> <height> - set canvas resolution\n"
            "snap <comment> - take a snapshot with comment in the filename\n"
            "start(_recording) - record all user commands for later playback\n"
            "stop(_recording) - finish recording of user commands\n"
            );

        // commands: mind the trailing blank!
        const std::string CMD_EDIT_VOLUME     = "edit_volume ";
        const std::string CMD_COPY_VOLUME     = "copy_volume ";
        const std::string CMD_FILTER_VOLUME   = "filter_volume ";
        const std::string CMD_TEST_CONTENTS   = "test_contents ";
        const std::string CMD_HEIGHTFIELD_AMPMAP  = "hm_map ";
        const std::string CMD_SET_RESOLUTION  = "set_resolution ";
        const std::string CMD_SET_IMAGE_FILE_RESOLUTION  = "set_image_file_resolution ";
        const std::string CMD_COMPRESSION     = "cmp ";
        const std::string CMD_SAVE_COLORMAP   = "save_colormap ";
        const std::string CMD_CLUSTER         = "cluster ";

        if (command.find(CMD_EDIT_VOLUME) == 0)
        {
            std::string const msg = this->user_command_edit_volume(command);
            response_arguments.set_value("command_output", msg.c_str());
        }
        else if (command.find(CMD_COPY_VOLUME) == 0)
        {
            std::string const msg = this->user_command_copy_volume(command);
            response_arguments.set_value("command_output", msg.c_str());
        }
        else if (command.find(CMD_FILTER_VOLUME) == 0)
        {
            std::string const msg = this->user_command_filter_volume(command);
            response_arguments.set_value("command_output", msg.c_str());
        }
        else if (command.find(CMD_TEST_CONTENTS) == 0)
        {
            // verify the volume contents
            std::string const param = command.substr(CMD_TEST_CONTENTS.size());

            std::ostringstream msg;
            msg << "test synthetic data: " << param << ".";
            test_verify_synthetic_volume_contents(param, cur_scope, m_session_tag);

            response_arguments.set_value("command_output", msg.str().c_str());
        }
        else if (command.find(CMD_HEIGHTFIELD_AMPMAP) == 0)
        {
            // FIXME: Need test due to valuerange/bbox IF update
            bool const ret = generate_heightfield_amplitude_map(command, cur_scope,
                                                                m_iindex_session,
                                                                m_session_tag);
            std::string const retstr = CMD_HEIGHTFIELD_AMPMAP + (ret ? "... success" : "... failed");
            response_arguments.set_value("command_output", retstr.c_str());
        }
        else if (command.find(CMD_SET_RESOLUTION) == 0)
        {
            bool const ret = command_set_resolution(command, m_irc, m_session_tag);
            std::string const retstr = CMD_SET_RESOLUTION + (ret ? "... success" : "... failed");
            response_arguments.set_value("command_output", retstr.c_str());
        }
        else if (command.find(CMD_SET_IMAGE_FILE_RESOLUTION) == 0)
        {
            bool const ret = command_set_image_file_resolution(command);
            std::string const retstr = CMD_SET_IMAGE_FILE_RESOLUTION + (ret ? "... success" : "... failed");
            response_arguments.set_value("command_output", retstr.c_str());
        }
        else if (command.find(CMD_SAVE_COLORMAP) == 0)
        {
            const std::string result = user_command_save_colormap(command);
            response_arguments.set_value("command_output", result.c_str());
        }
        else if (command.find(CMD_COMPRESSION) == 0)
        {
            std::istringstream args(command.substr(CMD_COMPRESSION.size()));

            mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
                cur_scope->create_transaction<mi::neuraylib::IDice_transaction>());
            assert(dice_transaction.is_valid_interface());
            mi::base::Handle<const nv::index::ISession> session(
                dice_transaction->access<const nv::index::ISession>(
                    m_session_tag));

            bool show_usage = false;
            if (!args.eof())
            {
                mi::base::Handle<nv::index::IConfig_settings> configs_edit(
                    dice_transaction->edit<nv::index::IConfig_settings>(
                        session->get_config()));

                nv::index::IConfig_settings::Data_transfer_config cfg
                    = configs_edit->get_data_transfer_config();

                std::string subcmd;
                mi::Sint32 param = -1;
                if ((args >> subcmd) && (args >> param) && param >= 0)
                {
                    if (subcmd == "show")
                    {
                        // Do nothing
                    }
                    else if (subcmd == "span" || subcmd == "sl")
                    {
                        cfg.span_compression_level = param;
                    }
                    else if (subcmd == "tile" || subcmd == "tl")
                    {
                        cfg.tile_compression_level = param;
                    }
                    else if (subcmd == "spanenc" || subcmd == "se")
                    {
                        cfg.span_image_encoding = (param != 0);
                    }
                    else if (subcmd == "tileenc" || subcmd == "te")
                    {
                        cfg.tile_image_encoding = (param != 0);
                    }
                    else if (subcmd == "spanalpha" || subcmd == "sa")
                    {
                        cfg.span_alpha_channel = (param != 0);
                    }
                    else if (subcmd == "smode" || subcmd == "sm")
                    {
                        cfg.span_compression_mode = param;
                    }
                    else if (subcmd == "tmode" || subcmd == "tm")
                    {
                        cfg.tile_compression_mode = param;
                    }
                    else if (subcmd == "sthreads" || subcmd == "st")
                    {
                        cfg.span_compression_threads = param;
                    }
                    else if (subcmd == "tthreads" || subcmd == "tt")
                    {
                        cfg.tile_compression_threads = param;
                    }
                    else if (subcmd == "zlevel" || subcmd == "zl")
                    {
                        cfg.z_buffer_compression_level = param;
                    }
                    else if (subcmd == "zmode" || subcmd == "zm")
                    {
                        cfg.z_buffer_compression_mode = param;
                    }
                    else if (subcmd == "zthreads" || subcmd == "zt")
                    {
                        cfg.z_buffer_compression_threads = param;
                    }
                    else
                    {
                        show_usage = true;
                    }
                }
                else
                {
                    show_usage = true;
                }

                configs_edit->set_data_transfer_config(cfg);
            }

            mi::base::Handle<const nv::index::IConfig_settings> configs(
                dice_transaction->access<nv::index::IConfig_settings>(
                    session->get_config()));

            nv::index::IConfig_settings::Data_transfer_config cfg = configs->get_data_transfer_config();

            std::ostringstream msg;

            if (show_usage)
                msg << "Usage: cmp <show|span|tile|spanenc|tileenc|smode|tmode|sthreads|threads|spanalpha|zlevel|zmode|zthreads> [param]\n\n";

            msg << "Compression settings:\n"
                << "  span compression level: " << (int)cfg.span_compression_level << "\n"
                << "  tile compression level: " << (int)cfg.tile_compression_level << "\n"
                << "  span image encoding: " << (int)cfg.span_image_encoding << "\n"
                << "  tile image encoding: " << (int)cfg.tile_image_encoding << "\n"
                << "  span alpha channel: " << (int)cfg.span_alpha_channel << "\n"
                << "  span mode:          " << (int)cfg.span_compression_mode << "\n"
                << "  tile mode:          " << (int)cfg.tile_compression_mode << "\n"
                << "  span threads:       " << (int)cfg.span_compression_threads << "\n"
                << "  tile threads:       " << (int)cfg.tile_compression_threads << "\n"
                << "  z-buffer level:     " << (int)cfg.z_buffer_compression_level << "\n"
                << "  z-buffer mode:      " << (int)cfg.z_buffer_compression_mode << "\n"
                << "  z-buffer threads:   " << (int)cfg.z_buffer_compression_threads << "\n";

            response_arguments.set_value("command_output", msg.str().c_str());

            dice_transaction->commit();
        }
        else if (command == "cancel_rendering")
        {
            const bool canceled = m_iindex_rendering->cancel_rendering(m_frame_info_callbacks->get_frame_identifier());
            if(!canceled)
            {
                INFO_LOG << "Tried to cancel a frame but the rendering was already finished.";
            }
        }
        else if (command == "start_recording" || command == "start")
        {
            // Do not allow the user to specify the filename for security reasons (it could be
            // coming of an unsecure network connection)
            const std::string filename = "index.rec";
            start_recording(filename);
            response_arguments.set_value(
                "command_output", std::string("Recording of all commands to file "
                                              + filename + " was started.").c_str());
        }
        else if (command == "stop_recording" || command == "stop")
        {
            stop_recording();
            response_arguments.set_value("command_output", "Recording stopped.");
        }
        else if (command == "snap" || command.find("snap ") == 0)
        {
            std::string comment;
            if (command.length() > 5)
                comment = command.substr(5);
            INFO_LOG << "Taking snapshot with comment '" << comment << "'";

            if (!comment.empty())
            {
                // Build a safe filename out of the snap comment
                comment = comment.substr(0, 100); // Restrict length
                for (size_t i = 0; i < comment.length(); ++i)
                {
                    char& c = comment[i];
                    if (!isalnum(c) && c != '-') // replace non-alphanumeric characters
                        c = '_';
                }
                comment = "_" + comment;
            }

            std::string format = appdata->peek_app_proj()->get(
                "app::camera::snapshot_filename", "snap_%03d.ppm");

            // Insert comment before file extension
            size_t p = format.find_last_of(".");
            if (p != std::string::npos)
                format = format.substr(0, p) + comment + format.substr(p);

            appdata->set_snapshot_filename(format);
            appdata->set_snapshot_on(true);
            response_arguments.set_value("command_output", "Taking a snapshot.");

            if (s_is_recording)
            {
                // Record information about current camera parameters. This will allow keeping the
                // same view when interaction changes. The "ICamera::Set" command will still have to
                // be uncommented and moved before the "snap" to become active.
                mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
                    cur_scope->create_transaction<mi::neuraylib::IDice_transaction>());
                {
                    Multiple_stereo_camera* ms_camera = Nvindex_AppData::instance()->
                        get_user_interaction(0)->get_multiple_stereo_camera();
                    mi::neuraylib::Tag camera_tag = ms_camera->get_main_base_camera_tag();
                    if (ms_camera->get_use_ortho_camera())
                    {
                        camera_tag = ms_camera->get_ortho_camera_tag();
                    }

                    const mi::base::Handle< const nv::index::ICamera > cam(
                        dice_transaction->access< nv::index::ICamera >(camera_tag));

                    const mi::math::Vector<mi::Float32, 3> eye = cam->get_eye_point();
                    const mi::math::Vector<mi::Float32, 3> dir = cam->get_view_direction();
                    const mi::math::Vector<mi::Float32, 3> up  = cam->get_up_direction();

                    std::ostringstream os;
                    os << std::setprecision (15)
                       << "#ICamera::Set:\n"
                       << "    comment: \"Camera parameters used for previous 'snap':\"\n"
                       << "    from: "<< eye.x << " " << eye.y << " " << eye.z << "\n"
                       << "    dir: " << dir.x << " " << dir.y << " " << dir.z << "\n"
                       << "    up: "  << up.x  << " " << up.y  << " " << up.z  << "\n";

                    mi::base::Handle<const nv::index::IPerspective_camera>  perspective_camera(
                        cam->get_interface<const nv::index::IPerspective_camera>());
                    if (perspective_camera)
                    {
                        os << "    aspect: "   << perspective_camera->get_aspect()   << "\n"
                           << "    aperture: " << perspective_camera->get_aperture() << "\n"
                           << "    focal: "    << perspective_camera->get_focal()    << "\n"
                           << "    clip_min: " << perspective_camera->get_clip_min() << "\n"
                           << "    clip_max: " << perspective_camera->get_clip_max() << "\n"
                           << "    orthographic: " << "no" << "\n";
                    }
                    mi::base::Handle<const nv::index::IOrthographic_camera> ortho_camera(
                        cam->get_interface<const nv::index::IOrthographic_camera>());
                    if (ortho_camera)
                    {
                        os << "    aspect: "   << ortho_camera->get_aspect()   << "\n"
                           << "    aperture: " << ortho_camera->get_aperture() << "\n"
                           << "    clip_min: " << ortho_camera->get_clip_min() << "\n"
                           << "    clip_max: " << ortho_camera->get_clip_max() << "\n"
                           << "    orthographic: " << "yes" << "\n";
                    }

                    mi::base::Lock::Block block(&s_recording_lock);
                    s_recording_file << os.str();
                }
                dice_transaction->commit();
            }
        }
        else if (command == "perf" || command.find("perf ") == 0)
        {
            // measure the performance for specified frames
            std::string ret_str;
            log_performance_by_command_str(command, ret_str);
            response_arguments.set_value("command_output", ret_str.c_str());
        }
        else if (command.find("perf_logger") == 0)
        {
            // switching performance logger
            std::string ret_str;
            setup_performance_logger_by_command_str(command, ret_str);
            response_arguments.set_value("command_output", ret_str.c_str());
        }
        else if (command.find("animate_camera") == 0)
        {
            // camera animation
            std::string ret_str;
            animate_camera_by_command_str(command, ret_str);
            response_arguments.set_value("command_output", ret_str.c_str());
        }
        else if (command == "dump" || command == "export")
        {
            response_arguments.set_value("command_output", "Exporting internal state...");

            mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
                cur_scope->create_transaction<mi::neuraylib::IDice_transaction>());
            assert(dice_transaction.is_valid_interface());
            mi::base::Handle<const nv::index::ISession> session(
                dice_transaction->access<const nv::index::ISession>(
                    m_session_tag));


            std::cout << "Exporting internal state of the session:\n\n"
                      << "##################################################\n" << std::endl;

            export_session(m_irc, dice_transaction.get());

            std::cout << "##################################################\n" << std::endl;

            dice_transaction->commit();
        }
        else if (command.find(CMD_CLUSTER) == 0)
        {
            const std::string param = command.substr(CMD_CLUSTER.size());

            mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
                cur_scope->create_transaction<mi::neuraylib::IDice_transaction>());
            assert(dice_transaction.is_valid_interface());

            // Pass command to all cluster nodes
            Distributed_command_execution job(param, m_session_tag);
            dice_transaction->execute_fragmented(&job, job.get_nb_of_fragments());

            dice_transaction->commit();

            return true;
        }
        else if (command.find("set ") == 0)
        {
            // Overwrite project file settings.
            // Will only have an effect if the settings are read from the dictionary again after
            // they have been changed. This works for the video streaming settings, since they
            // are read each time a new client connects, i.e. when reloading the page.
            std::string s = command.substr(4);
            size_t p = s.find(" = ");
            if (p != std::string::npos)
            {
                std::string key = s.substr(0, p);
                std::string value = s.substr(p + 3);

                // Only allow specific settings for security reasons
                if ((key == "dice::rtmp_video_streaming::video_codec" &&
                     (value == "screen video" || value == "h264" || value == "h264_nvenc")) ||
                    key == "dice::rtmp_video_streaming::video_bitrate" ||
                    key == "dice::rtmp_video_streaming::video_framerate")
                {
                    appdata->peek_app_proj()->insert(key, value);

                    response_arguments.set_value(
                        "command_output",
                        std::string("Set config value '" + key + "' to '" + value + "'").c_str());
                }
                else
                {
                    response_arguments.set_value("command_output", "Denied.");
                }
            }
            else
            {
                response_arguments.set_value("command_output", "Error: Expected 'key = value'.");
            }
        }
        else if (command == "quit")
        {
            response_arguments.set_value("command_output", "End of Line.");
            appdata->set_app_run(false);
        }
        else if (command == "help")
        {
            response_arguments.set_value("command_output", usage.c_str());
        }
        else
        {
            const std::string msg = "Unknown command '" + command + "'\n\n" + usage;
            response_arguments.set_value("command_output", msg.c_str());
        }
    }
    else if (strcmp(cmd, "ICamera::predefined_view::setoption") == 0)
    {
        this->camera_predefined_view_set_option(user_arguments);
    }
    else if (strcmp(cmd, "ICamera::predefined_view::getoption") == 0)
    {
        this->camera_predefined_view_get_option(response_arguments);
    }
    else if (strcmp(cmd, "ICamera::predefined_view::action") == 0)
    {
        this->camera_predefined_view_action();
    }
    else if (strcmp(cmd, "ICamera::view_all::setoption") == 0)
    {
        ERROR_LOG << "ICamera::view_all::setoption is no longer used. Ignored."; // 2015-12
    }
    else if (strcmp(cmd, "ICamera::view_all::getoption") == 0)
    {
        ERROR_LOG << "ICamera::view_all::getoption is no longer used. Ignored."; // 2015-12
    }
    else if (strcmp(cmd, "ICamera::view_all") == 0)
    {
        ERROR_LOG << "ICamera::view_all is no longer used. Ignored."; // 2015-12
    }
    else if (strcmp(cmd, "ICamera::Camera::SetOption") == 0)
    {
        this->camera_set_option(user_arguments);
    }
    else if (strcmp(cmd, "ICamera::Camera::GetOption") == 0)
    {
        this->camera_get_option(response_arguments);
    }
    else if (strcmp(cmd, "ICamera::Set") == 0)
    {
        this->camera_set_parameter(user_arguments);
    }
    else if (strcmp(cmd, "ICamera::SetView") == 0)
    {
        this->camera_set_view_action();
    }
    else if (strcmp(cmd, "ICamera::RestoreView") == 0)
    {
        this->camera_restore_view_action();
    }
    else if (strcmp(cmd, "ICamera::print_param") == 0)
    {
        this->camera_print_param_action();
    }
    else if (strcmp(cmd, "IColormap::userdef::option") == 0)
    {
        mi::Sint32 op_type = get_sint32_userarg(user_arguments, "colormap_userdef_operation_type");
        Nvindex_AppData::instance()->get_colormap_manager()->set_colormap_userdef_operation_type(op_type);
    }
    else if (strcmp(cmd, "IColormap::userdef::action") == 0)
    {
        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            cur_scope->create_transaction<mi::neuraylib::IDice_transaction>());
        assert(dice_transaction.is_valid_interface());
        {
            bool const ret = exec_userdef_colormap_operation(
                get_colormap_tag(Nvindex_AppData::instance()->get_current_colormap_index()),
                dice_transaction.get());
            if (!ret)
            {
                WARN_LOG << "Fail to userdef colormap operation.";
            }
        }
        dice_transaction->commit();
    }
    else if(strcmp(cmd, "horizon_workflow_task") == 0)
    {
        const std::string workflow_task = user_arguments.get_value("workflow_task");
        const mi::Uint32  workflow_task_id = get_uint32(workflow_task);

        // choose the workflow
        Heightfield_workflow_functionality::instance()->set_workflow(
            Heightfield_workflow_functionality::Workflow(workflow_task_id));
        response_arguments.set_value("next_workflow_operation",
                                     Heightfield_workflow_functionality::instance()->get_next_step().c_str());

        // set the heightfield names to GUI
        std::vector< std::string > result_heightfield_name_vec;
        Nvindex_AppData::instance()->get_heightfield_name_vec(cur_scope, result_heightfield_name_vec);

        const std::string heightfield_name_gui = concat_string_with_separator(result_heightfield_name_vec, "|");
        response_arguments.set_value("names", heightfield_name_gui.c_str());
    }
    else if(strcmp(cmd, "horizon_workflow_task_horizon_id") == 0)
    {
        const mi::Sint32 selected_heightfield = get_sint32_userarg(user_arguments, "horizon_id");
        const mi::neuraylib::Tag& selected_heightfield_tag =
            Nvindex_AppData::instance()->get_heightfield_tag_from_idx(selected_heightfield);
        if(selected_heightfield_tag.is_valid())
        {
            Heightfield_workflow_functionality::instance()->set_heightfield(selected_heightfield_tag);
            response_arguments.set_value("next_workflow_operation",
                                         Heightfield_workflow_functionality::instance()->get_next_step().c_str());
        }
    }
    else if(strcmp(cmd, "horizon_workflow_task_seismic_id") == 0)
    {
        const mi::Sint32 selected_volume_id = get_sint32_userarg(user_arguments, "seismic_id");
        const mi::neuraylib::Tag& selected_volume_tag =
            Nvindex_AppData::instance()->get_volume_tag_from_idx(selected_volume_id);
        if(selected_volume_tag.is_valid())
        {
            Heightfield_workflow_functionality::instance()->set_associated_scene_element(selected_volume_tag);
            response_arguments.set_value("next_workflow_operation",
                                         Heightfield_workflow_functionality::instance()->get_next_step().c_str());
        }
    }
    else if(strcmp(cmd, "horizon_polygon_delete") == 0)
    {
        if ( Heightfield_workflow_functionality::instance()->is_heightfield_workflow_delete_operation())
        {
            const bool is_use_convex_hull = false;
            heightfield_workflow_delete_polygon(cur_scope, is_use_convex_hull);
        }
    }
    else if(strcmp(cmd, "horizon_elevation_change") == 0)
    {
        const std::string elevation_value_str = user_arguments.get_value("elevation_value");
        const mi::Float32 elevation_value = get_float32(elevation_value_str);
        if(Heightfield_workflow_functionality::instance()->is_heightfield_workflow_elevation_change_operation() )
        {
            const bool is_scale_op        = false;
            const bool is_use_convex_hull = false;
            heightfield_workflow_change_elevation_values(cur_scope, elevation_value,
                                                         is_scale_op, is_use_convex_hull);
        }
    }
    else if(strcmp(cmd, "horizon_gridding") == 0)
    {
        if (Heightfield_workflow_functionality::instance()->is_heightfield_workflow_gridding_operation())
        {
            heightfield_workflow_gridding(cur_scope, false);
        }
    }
    else if(strcmp(cmd, "horizon_autotracking") == 0)
    {
        if (Heightfield_workflow_functionality::instance()->is_autotracking_workflow_operation())
        {
            heightfield_workflow_autotracking(cur_scope);
        }
    }
    // export heightfield
    else if(strcmp(cmd, "export_heightfield") == 0)
    {
        this->export_heightfield_by_gui(user_arguments);
    }
    else if(strcmp(cmd, "get_export_heightfield_name_list") == 0)
    {
        this->get_export_heightfield_name_list(response_arguments);
    }
    else if(strcmp(cmd, "set_heightfield_index_change") == 0)
    {
        this->set_export_heightfield_index_change(user_arguments);
    }
    // export volume
    else if(strcmp(cmd, "export_volume") == 0)
    {
        this->export_volume_by_gui(user_arguments);
    }
    else if(strcmp(cmd, "get_export_volume_name_list") == 0)
    {
        this->get_export_volume_name_list(response_arguments);
    }
    else if(strcmp(cmd, "set_export_volume_index_change") == 0)
    {
        // Changed the selected volume of the combo box GUI. But nothing to do for now.
    }
    else if(strcmp(cmd, "get_attrgen_gui_section_info") == 0)
    {
        this->get_attrgen_gui_section_info(response_arguments);
    }
    else if(strcmp(cmd, "attrgen_copy_volume") == 0)
    {
        ERROR_LOG << "called attrgen_copy_volume. But not implemented yet.";
        this->attrgen_copy_volume_by_gui(user_arguments);
    }
    else if(strcmp(cmd, "attrgen_filter_volume") == 0)
    {
        ERROR_LOG << "called attrgen_filter_volume. But not implemented yet.";
        this->attrgen_filter_volume_by_gui(user_arguments);
    }
    else if(strcmp(cmd, "attrgen_edit_volume") == 0)
    {
        ERROR_LOG << "called attrgen_filter_volume. But not implemented yet.";
        this->attrgen_edit_volume_by_gui(user_arguments);
    }
    else if (strcmp(cmd, "set_rtc_params") == 0)
    {
        this->rtc_edit_param_buffer_by_gui(user_arguments);
    }
    else if(strcmp(cmd, "get_performance_logging_state") == 0)
    {
        // performance monitor panel
        this->get_performance_logging_state_by_gui(response_arguments);
    }
    else if(strcmp(cmd, "test_interactive_image_file_snapshot") == 0)
    {
        // Snap shot by the GUI. This is just a trigger. When the next
        // frame is finished, snapshot is written.
        INFO_LOG << "Image file snapshot button has been pushed.";
        Nvindex_AppData::instance()->set_snapshot_on(true);
    }
    else if(strcmp(cmd, "prepare_last_snapshot") == 0)
    {
        // image file download to the browser
        if (Nvindex_AppData::instance()->get_snapshot_index() > 0)
        {
            std::string snap = Nvindex_AppData::instance()->prepare_last_snapshot();
            if (snap.empty()){
                response_arguments.set_value("no_snapshot_yet", "1");
                WARN_LOG << "Snapshot image not ready.";
            }
            else{
                response_arguments.set_value("snapshot_filename", snap.c_str());
                INFO_LOG << "Downloading snapshot image. [" << snap << "]";
            }
        }
        else
        {
            response_arguments.set_value("no_snapshot_yet", "1");
            WARN_LOG <<"Snapshot image not found.";
        }
    }

    return true;
}

    mi::Sint32 Call_event_handler::play_recording(const std::string& filename, mi::Sint32 start_line)
    {
        if (s_is_recording)
        {
            ERROR_LOG << "Can't play while recording";
            return -1;
        }

        std::ifstream f(filename.c_str());
        if (!f.good())
        {
            ERROR_LOG << "Could not open recording file '" << filename << "' for playback";
            return -1;
        }

        mi::base::Handle<mi::neuraylib::IFactory> factory(m_iindex_if->get_api_component<mi::neuraylib::IFactory>());
        mi::base::Handle<mi::IMap> imap(factory->create<mi::IMap>("Map<Interface>"));

        std::string cmd;
        mi::Sint32 line_counter = 0;
        bool is_snap = false;
        bool is_perf = false;
        while(!f.eof())
        {
            std::string line;
            getline(f, line);
            line_counter++;

            if (line_counter < start_line)
                continue;

            size_t first_non_whitespace = line.find_first_not_of(" \t");

            // Skip lines containing only whitespace or empty lines
            if (first_non_whitespace == std::string::npos)
                continue;

            // Skip lines starting with '#'
            if (line[first_non_whitespace] == '#')
                continue;

            // Initial whitespace: Expect a parameters key/value pair
            if (first_non_whitespace > 0)
            {
                size_t colon = line.find(": ");
                if (colon == std::string::npos)
                {
                    ERROR_LOG << "Missing ': ' after key in " << filename << ":" << line_counter;
                    return -1;
                }
                else
                {
                    std::string key = line.substr(first_non_whitespace, colon - first_non_whitespace);
                    std::string value = line.substr(colon + 2);
                    if (cmd.empty())
                    {
                        ERROR_LOG << "Read parameters '" << key << "' but no previous command in "
                                  << filename << ":" << line_counter;
                        return -1;
                    }
                    else
                    {
                        // Note if the command is for taking a snapshot
                        if (cmd == "user_command" && key == "user_command" && value.find("snap") == 0)
                            is_snap = true;

                        // ...or a performance measurement
                        if (cmd == "user_command" && key == "user_command" && value.find("perf") == 0)
                            is_perf = true;

                        // Store key/value pair in the map
                        mi::base::Handle<mi::IString> str(factory->create<mi::IString>("String"));
                        str->set_c_str(value.c_str());
                        imap->insert(key.c_str(), str.get());
                    }
                }
            }
            // No initial whitespace: It must be a command
            else
            {
                if (line[line.length() - 1] != ':')
                {
                    ERROR_LOG << "Missing ':' as last character in " << filename << ":" << line_counter;
                    return -1;
                }
                else
                {
                    // Finish previous command
                    if (!imap->empty()) {
                        mi::IData* response_args = 0;
                        handle(0, 0, 0, imap.get(), &response_args);
                        if (response_args)
                            response_args->release();

                        // If the previous command was 'snap' then we exit now to let the rendering
                        // loop run and actually render and take the snapshot. After that this
                        // function will be called again starting at the current line.
                        if (is_snap || is_perf)
                            return line_counter;
                    }

                    // Store command into the map
                    cmd = line.substr(0, line.length() - 1);
                    imap->clear();
                    mi::base::Handle<mi::IString> str(factory->create<mi::IString>("String"));
                    str->set_c_str(cmd.c_str());
                    imap->insert("cmd", str.get());
                }
            }
        }

        // Handle the last command
        if (!imap->empty()) {
            mi::IData* response_args = 0;
            handle(0, 0, 0, imap.get(), &response_args);
            if (response_args)
                response_args->release();
        }

        return 0; // We are done
    }

    void Call_event_handler::start_recording(const std::string& filename)
    {
        if (s_is_recording)
            return;

        mi::base::Lock::Block block(&s_recording_lock);
        s_recording_file.open(filename.c_str());
        if (s_recording_file.good())
        {
            s_is_recording = true;
            INFO_LOG << "Recording started to file: '" << get_current_directory_name_string()
                     << "/" << filename << "'";
        }
        else
        {
            ERROR_LOG << "Could not open recording file '" << filename << "'";
            s_is_recording = false;
        }
    }

    void Call_event_handler::stop_recording()
    {
        if (!s_is_recording)
            return;

        mi::base::Lock::Block block(&s_recording_lock);
        s_recording_file.close();
        s_is_recording = false;
        INFO_LOG << "Recording stopped.";
    }

    void Call_event_handler::record_received_command(
        const std::string&               command,
        mi::base::Handle<const mi::IMap> imap)
    {
        if (!s_is_recording)
            return;

        // Skip uninteresting commands that are sent automatically be the viewer
        if (command == "stats" || command.find("get") == 0)
            return;

        mi::base::Lock::Block block(&s_recording_lock);

        // Write output in a simple YAML-format:
        // command1:
        //     key1: foo bar
        //     key2: bar baz
        // command2:
        //     key3: foobar
        s_recording_file << command << ":\n"; // Write command

        // Write parameters
        for (mi::Uint32 i = 0; i < imap->get_length(); ++i)
        {
            std::string key = imap->get_key(i);
            if (key == "cmd")
                continue; // Command was already written

            mi::base::Handle<const mi::IString> value(imap->get_value<mi::IString>(i));
            s_recording_file << "    " << key << ": " << value->get_c_str() << "\n";
        }
    }

//----------------------------------------------------------------------
    std::string Call_event_handler::user_command_edit_volume(const std::string & com)
    {
        const std::string com_usage_msg =
            "Usage: edit_volume [volume_tag] [amplitude_value]\n"
            "  - volume_tag: the volume tag to edit\n"
            "  - amplitude_value: the amplitude value to be filled in the volume.\n"
            "    When < 0, set host assign mode"
            ;

        std::vector< std::string > tokens;
        Tokenizer::parse(com, " ", tokens);
        if(tokens.size() != 3){
            const std::string mes = "Illegal command.\n";
            ERROR_LOG << (mes + com_usage_msg);
            return (mes + com_usage_msg);
        }

        // [0] command
        if(tokens[0] != "edit_volume"){
            const std::string mes = "Not edit volume command.\n";
            ERROR_LOG << (mes + com_usage_msg);
            return (mes + com_usage_msg);
        }

        // [1] volume tag
        mi::neuraylib::Tag volume_tag;
        volume_tag.id = get_uint32(tokens[1]);
        if(!volume_tag.is_valid()){
            const std::string mes = "Illegal volume tag: " + tokens[1];
            ERROR_LOG << mes;
            return mes;
        }

        // [2] amplitude value
        const mi::Sint32 amplitude_value_sint = get_sint32(tokens[2]);
        bool is_host_assign_mode = false;
        if(amplitude_value_sint < 0){
            is_host_assign_mode = true;
        }
        else if(amplitude_value_sint > 255){
            const std::string mes = "Illegal amplitude value: out of range (should be [0,255]): "
                + tokens[2];
            ERROR_LOG << mes;
            return mes;
        }

        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            get_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
        assert(dice_transaction.is_valid_interface());

        const mi::Uint8 amplitude_value = static_cast<mi::Uint8>(amplitude_value_sint);
        std::string result_str = "volume_set_amplitude_value: ok";
        const bool is_success =
            volume_set_amplitude_value(
                is_host_assign_mode,
                amplitude_value,
                m_session_tag,
                volume_tag,
                dice_transaction.get());
        if(is_success){
            dice_transaction->commit();
        }
        else{
            dice_transaction->abort();
            result_str = "volume_set_amplitude_value: failed!";
        }

        return result_str;
    }

//----------------------------------------------------------------------
/// user_command_copy_volume's subroutine.
///
/// Note: This function the system has enough memory for the copy.
///
/// \param[in]  session_tag      session tag
/// \param[in]  src_volume_tag   input copy source volume tag
/// \param[in]  dst_group_tag    destination group tag. The volume will be created under this group.
/// \param[in]  trans_vec        translation vector for the copy destination volume
/// \param[in]  dice_transaction db transaction
/// \param[out] err_mes          error message when failed
/// \return a newly created volume tag. NULL_TAG when copy failed.
static mi::neuraylib::Tag user_command_copy_volume_sub(
    const mi::neuraylib::Tag & session_tag,
    const mi::neuraylib::Tag & src_volume_tag,
    const mi::neuraylib::Tag & dst_group_tag,
    const mi::math::Vector<mi::Float32, 3> trans_vec,
    mi::neuraylib::IDice_transaction * dice_transaction,
    std::string & err_mes)
{
    assert(dice_transaction != 0);

    // check the input volume
    mi::base::Handle<const nv::index::IRegular_volume> volume(
        dice_transaction->access<const nv::index::IRegular_volume>(src_volume_tag));
    if(!(volume.is_valid_interface())){
        std::stringstream sstr;
        sstr << "Invalid source volume tag: " << src_volume_tag << ", must be exist and a IRegular_volume.";
        err_mes = sstr.str();
        return mi::neuraylib::NULL_TAG;
    }
    const std::string src_volume_name = volume->get_name();
    const std::string dst_volume_name = src_volume_name + ".copy";

    // check the destination static group
    mi::base::Handle<nv::index::IStatic_scene_group> static_group_node(
        dice_transaction->edit<nv::index::IStatic_scene_group>(dst_group_tag));
    if(!(static_group_node.is_valid_interface())){
        std::stringstream sstr;
        sstr << "Invalid destination group tag: " << dst_group_tag << ", must be exist and a IStatic_scene_group.";
        err_mes = sstr.str();
        return mi::neuraylib::NULL_TAG;
    }

    // generate a new volume by copy (with translation)
    mi::neuraylib::Tag copy_volume_tag =
        generate_volume_by_copy(
            session_tag,
            src_volume_tag,
            dst_volume_name,
            trans_vec,
            dice_transaction);
    if(!copy_volume_tag.is_valid()){
        std::stringstream sstr;
        sstr << "Failed to create a volume.";
        err_mes = sstr.str();
        return mi::neuraylib::NULL_TAG;
    }

    // add to the same group of the source volume
    static_group_node->append(copy_volume_tag, dice_transaction);
    err_mes = "";

    return copy_volume_tag;
}

//----------------------------------------------------------------------
std::string Call_event_handler::user_command_copy_volume(const std::string & com)
{
    const std::string com_usage_msg =
        "Usage: copy_volume [src_volume_tag] [dst_group_tag] [x] [y] [z]\n"
        "  - src_volume_tag: copy source volume tag. Must be a volume\n"
        "  - dst_group_tag:  copy destination group tag. Must be a Static_group_node.\n"
        "  - x y z:          translation (x, y, z)\n"
        "Note: This command assumes the system has enough memory for the new volume!";

    std::vector< std::string > tokens;
    Tokenizer::parse(com, " ", tokens);
    if(tokens.size() != 6){
        const std::string mes = "Illegal command.\n";
        ERROR_LOG << (mes + com_usage_msg);
        return (mes + com_usage_msg);
    }

    if(tokens[0] != "copy_volume"){
        const std::string mes = "Not copy volume command.\n";
        ERROR_LOG << (mes + com_usage_msg);
        return (mes + com_usage_msg);
    }

    mi::neuraylib::Tag src_volume_tag;
    mi::neuraylib::Tag dst_group_tag;
    src_volume_tag.id = get_uint32(tokens[1]);
    dst_group_tag.id  = get_uint32(tokens[2]);
    const mi::Float32 tr_x = get_float32(tokens[3]);
    const mi::Float32 tr_y = get_float32(tokens[4]);
    const mi::Float32 tr_z = get_float32(tokens[5]);
    const mi::math::Vector<mi::Float32, 3> trans_vec(tr_x, tr_y, tr_z);

    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        get_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    std::string err_mes;

    const mi::neuraylib::Tag copy_volume_tag =
        user_command_copy_volume_sub(
            m_session_tag,
            src_volume_tag,
            dst_group_tag,
            trans_vec,
            dice_transaction.get(),
            err_mes);

    if(copy_volume_tag.is_valid()){
        dice_transaction->commit();
    }
    else{
        ERROR_LOG << err_mes;
        dice_transaction->abort();
    }

    return err_mes;
}

//----------------------------------------------------------------------
std::string Call_event_handler::user_command_filter_volume(const std::string & com)
{
    const std::string com_usage_msg =
        "  filter_volume <volume_tag> <filter_name>\n"
        "ex.\n"
        "  filte_volume 42 dip_median\n";

    std::vector< std::string > tokens;
    Tokenizer::parse(com, " ", tokens);
    if(tokens.size() != 3){
        const std::string mes = "Illegal command.\n";
        ERROR_LOG << (mes + com_usage_msg);
        return (mes + com_usage_msg);
    }

    if(tokens[0] != "filter_volume"){
        const std::string mes = "Not filter_volume command.\n";
        ERROR_LOG << (mes + com_usage_msg);
        return (mes + com_usage_msg);
    }

    mi::neuraylib::Tag volume_tag;
    volume_tag.id = get_uint32(tokens[1]);

    const mi::Sint32 ftype = Volume_data_filter::get_filter_type(tokens[2]);
    if(!Volume_data_filter::is_valid_filter_type(ftype)){
        return (std::string("Unknown filter name [") + tokens[2] + "]\n" +
                com_usage_msg);
    }

    // FIXME: parallel numbers
    const std::string pcount_str = Nvindex_AppData::instance()->peek_app_proj()->
        get("app::distributed_seismic_filter_filter_job::parallel_count", "1");
    const mi::Sint32 pcount = get_sint32(pcount_str);

    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        get_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    bool is_success = false;
    {
        is_success = volume_inset_filter_amplitude_values(
            m_session_tag, volume_tag, tokens[2], pcount, dice_transaction.get());
    }
    dice_transaction->commit();

    if(is_success){
        return "volume_filter_amplitude_values ... done.";
    }

    return "volume_inset_filter_amplitude_values ... failed.";
}

//----------------------------------------------------------------------
std::string Call_event_handler::user_command_save_colormap(const std::string & com)
{
    std::string const com_err_msg = "illegal usage of save_colormap.\n"
        "  save_colormap [colormap_fname]\n"
        "ex.\n"
        "  save_colormap 00cmap.cmap\n";

    // parse command
    std::vector< std::string > tokens;
    Tokenizer::parse(com, " ", tokens);
    if(tokens.size() != 2)
    {
        return com_err_msg;
    }

    // get color table
    std::vector<mi::math::Color_struct> color_table;
    this->get_user_current_colormap_value(color_table);

    std::string const colormap_fname = tokens[1];
    std::string err_mes;
    bool const is_succeeded = write_16bit_colormap(colormap_fname, color_table, err_mes);
    if(!is_succeeded)
    {
        return com_err_msg + err_mes;
    }
    const std::string res_mes = "Saved colormap '" + colormap_fname + "'.";
    INFO_LOG << res_mes;

    return res_mes;
}

//----------------------------------------------------------------------
void Call_event_handler::camera_predefined_view_set_option(
    Arguments_get_wrapper& user_args)
{
    assert(user_args.get_value("cmd") == std::string("ICamera::predefined_view::setoption"));

    const mi::Sint32 predef_view_index = get_sint32_userarg(user_args, "predefined_view::view_index");
    assert(predef_view_index >= 0);
    Nvindex_AppData::instance()->get_user_interaction(0)->
        get_examiner()->set_predefined_view_index(static_cast<mi::Sint32>(predef_view_index));
}

//----------------------------------------------------------------------
void Call_event_handler::camera_predefined_view_get_option(Arguments_set_wrapper& response_arguments)
{
    std::vector<std::string> predefined_view_name_vec = 
        Nvindex_AppData::instance()->get_user_interaction(0)->
        get_examiner()->get_predefined_view_name_vec();
    const std::string predefined_view_name_gui = concat_string_with_separator(predefined_view_name_vec, ",");
    response_arguments.set_value("predefined_view_name_list", predefined_view_name_gui.c_str());

    std::ostringstream sstr;
    sstr << Nvindex_AppData::instance()->get_user_interaction(0)->
        get_examiner()->get_predefined_view_index();
    response_arguments.set_value("view_index", sstr.str().c_str());
}

//----------------------------------------------------------------------
void Call_event_handler::camera_predefined_view_action()
{
    // INFO_LOG << "Call_event_handler::camera_predefined_view_action.";
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        get_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    {
        mi::base::Handle<const nv::index::ISession> session(
            dice_transaction->access<const nv::index::ISession>(m_session_tag));
        assert(session.is_valid_interface());

        mi::base::Handle<const nv::index::IScene> scene(
            dice_transaction->access<nv::index::IScene>(session->get_scene()));
        assert(scene.is_valid_interface());


        User_interaction* ui              = Nvindex_AppData::instance()->get_user_interaction(0);
        Multiple_stereo_camera* ms_camera = ui->get_multiple_stereo_camera();
        ui->get_examiner()->set_predefined_view_parameter_to_camera(
            session.get(), scene.get(), ms_camera, dice_transaction.get());
    }
    dice_transaction->commit();
}

//----------------------------------------------------------------------
void Call_event_handler::camera_set_option(Arguments_get_wrapper& user_args)
{
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        get_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    {
        assert(user_args.get_value("cmd") == std::string("ICamera::Camera::SetOption"));
        Multiple_stereo_camera* ms_camera =
            Nvindex_AppData::instance()->get_user_interaction(0)->get_multiple_stereo_camera();

        // animation on/off: TODO should be a push button
        bool const cam_anim = get_bool_userarg(user_args, "icamera_animation_on");
        if (cam_anim)
        {
            Nvindex_AppData::instance()->peek_camera_animator_data()->reset(true);
        }

        // clip min/max
        mi::Float32 clip_min = get_float32_userarg(user_args, "icamera_clipmin");
        mi::Float32 clip_max = get_float32_userarg(user_args, "icamera_clipmax");
        if (clip_min >= clip_max)
        {
            ERROR_LOG << "Can not set clip_min (" << clip_min
                      << ") >= clip_max (" << clip_max << ") value.";
        }
        else
        {
            mi::base::Handle< nv::index::IPerspective_camera > perspective_cam(
                dice_transaction->edit< nv::index::IPerspective_camera >(ms_camera->get_main_base_camera_tag()));
            assert(perspective_cam.is_valid_interface());
            perspective_cam->set_clip_min(clip_min);
            perspective_cam->set_clip_max(clip_max);


            mi::base::Handle< nv::index::IOrthographic_camera > ortho_cam(
                dice_transaction->edit< nv::index::IOrthographic_camera >(ms_camera->get_ortho_camera_tag()));
            assert(ortho_cam.is_valid_interface());
            ortho_cam->set_clip_min(clip_min);
            ortho_cam->set_clip_max(clip_max);

            if (!ms_camera->get_use_ortho_camera() &&
                get_bool_userarg(user_args, "icamera_orthographic"))
            {
                // transfer parameters from perspective_cam to ortho_cam
                ortho_cam->set(perspective_cam->get_eye_point(),
                               perspective_cam->get_view_direction(),
                               perspective_cam->get_up_direction());
                ortho_cam->set_aspect(perspective_cam->get_aspect());
                ortho_cam->set_clip_min(perspective_cam->get_clip_min());
                ortho_cam->set_clip_max(perspective_cam->get_clip_max());
            }
            else if (ms_camera->get_use_ortho_camera() &&
                     !get_bool_userarg(user_args, "icamera_orthographic"))
            {
                // transfer parameters from ortho_cam to perspective_cam
                perspective_cam->set(ortho_cam->get_eye_point(),
                                     ortho_cam->get_view_direction(),
                                     ortho_cam->get_up_direction());
                perspective_cam->set_aspect(ortho_cam->get_aspect());
                perspective_cam->set_clip_min(ortho_cam->get_clip_min());
                perspective_cam->set_clip_max(ortho_cam->get_clip_max());
            }
        }

        {
            ms_camera->set_use_ortho_camera(get_bool_userarg(user_args, "icamera_orthographic"));
        }

        // selected storage camera
        Nvindex_AppData::instance()->set_storage_camera_index(
            get_sint32_userarg(user_args, "icamera_multiple_camera"));

        // stereo mode
        const bool is_stereo = get_bool_userarg(user_args, "icamera_stereo_on");
        Nvindex_AppData::instance()->set_stereo_mode(is_stereo);
        if(!is_stereo)
        {
            // set stereo camera's current to base camera (left)
            ms_camera->peek_main_camera()->set_current_camera(Stereo_camera::SC_Left);
        }
        mi::Float32 eyesep = get_float32_userarg(user_args, "icamera_stereo_eyesep");
        ms_camera->peek_main_camera()->set_eye_separation(eyesep);
        // INFO_LOG << "DEBUG: GUI selected camera index: " << m_gui_selected_cam_idx;
    }
    dice_transaction->commit();
}

void Call_event_handler::camera_get_option(Arguments_set_wrapper& response_arguments)
{
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        get_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    {
        Multiple_stereo_camera* ms_camera =
            Nvindex_AppData::instance()->get_user_interaction(0)->get_multiple_stereo_camera();
        bool const is_anim = Nvindex_AppData::instance()->peek_camera_animator_data()->is_valid();

        std::ostringstream os;
        const bool is_ortho_cam = ms_camera->get_use_ortho_camera();
        if (is_ortho_cam)
        {
            mi::base::Handle< const nv::index::IOrthographic_camera > cam(
                dice_transaction->access< nv::index::IOrthographic_camera >(
                    ms_camera->get_ortho_camera_tag()));
            assert(cam.is_valid_interface());

            os << (is_anim ? "1" : "0") << ","
               << cam->get_clip_min() << ","
               << cam->get_clip_max() << ","
               <<  Nvindex_AppData::instance()->get_storage_camera_index() << ","
               << (Nvindex_AppData::instance()->is_stereo_mode() ? "1" : "0") << ","
               << (int)(ms_camera->get_main_camera().get_eye_separation()) << ","
               << (is_ortho_cam ? "1" : "0") << ","
               << "\n";
        }
        else
        {
            mi::base::Handle< const nv::index::IPerspective_camera > cam(
                dice_transaction->access< nv::index::IPerspective_camera >(
                    ms_camera->get_main_base_camera_tag()));
            assert(cam.is_valid_interface());

            os << (is_anim ? "1" : "0") << ","
               << cam->get_clip_min() << ","
               << cam->get_clip_max() << ","
               <<  Nvindex_AppData::instance()->get_storage_camera_index() << ","
               << (Nvindex_AppData::instance()->is_stereo_mode() ? "1" : "0") << ","
               << (int)(ms_camera->get_main_camera().get_eye_separation()) << ","
               << (is_ortho_cam ? "1" : "0") << ","
               << "\n";
        }
        response_arguments.set_value("icamera_option", os.str().c_str());
    }
    dice_transaction->commit();
}

void Call_event_handler::camera_set_parameter(Arguments_get_wrapper& args)
{
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        get_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    {
        const mi::math::Vector<mi::Float32, 3> from = get_vec_float32_3(get_string_userarg(args, "from"));
        const mi::math::Vector<mi::Float32, 3> up   = get_vec_float32_3(get_string_userarg(args, "up"));
        const mi::math::Vector<mi::Float32, 3> dir  = get_vec_float32_3(get_string_userarg(args, "dir"));
        const mi::Float32 aspect   = get_float32_userarg(args, "aspect");
        const mi::Float32 aperture = get_float32_userarg(args, "aperture");
        const mi::Float64 clip_min = get_float32_userarg(args, "clip_min");
        const mi::Float64 clip_max = get_float32_userarg(args, "clip_max");
        const bool orthographic    = get_bool(get_string_userarg(args, "orthographic"));

        Multiple_stereo_camera* ms_camera = 
            Nvindex_AppData::instance()->get_user_interaction(0)->get_multiple_stereo_camera();
        
        if (orthographic)
        {
            mi::base::Handle<nv::index::IOrthographic_camera> ortho_cam(
                dice_transaction->edit<nv::index::IOrthographic_camera>(
                    ms_camera->get_ortho_camera_tag()));
            assert(ortho_cam.is_valid_interface());

            ortho_cam->set(from, dir, up);
            ortho_cam->set_clip_min(clip_min);
            ortho_cam->set_clip_max(clip_max);
            ortho_cam->set_aperture(aperture);
            ortho_cam->set_aspect(aspect);
        }
        else
        {
            mi::base::Handle<nv::index::IPerspective_camera> perspective_cam(
                dice_transaction->edit<nv::index::IPerspective_camera>(
                    ms_camera->get_main_base_camera_tag()));
            assert(perspective_cam.is_valid_interface());

            perspective_cam->set(from, dir, up);
            perspective_cam->set_clip_min(clip_min);
            perspective_cam->set_clip_max(clip_max);
            perspective_cam->set_aperture(aperture);
            perspective_cam->set_aspect(aspect);
            // "focal" is in only in perspective camera
            const mi::Float32 focal = get_float32_userarg(args, "focal");
            perspective_cam->set_focal(focal);
        }
        ms_camera->set_use_ortho_camera(orthographic);
    }
    dice_transaction->commit();
}

void Call_event_handler::camera_set_view_action()
{
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        get_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    {
        Multiple_stereo_camera* ms_camera = 
            Nvindex_AppData::instance()->get_user_interaction(0)->get_multiple_stereo_camera();
        ms_camera->save_camera_state(Nvindex_AppData::instance()->get_storage_camera_index(),
                                  dice_transaction.get());
    }
    dice_transaction->commit();

    // to show the stored camera parameter
    this->camera_print_param_action();
}

void Call_event_handler::camera_restore_view_action()
{
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        get_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    {
        Multiple_stereo_camera* ms_camera = 
            Nvindex_AppData::instance()->get_user_interaction(0)->get_multiple_stereo_camera();
        ms_camera->load_camera_state(Nvindex_AppData::instance()->get_storage_camera_index(),
                                  dice_transaction.get());
    }
    dice_transaction->commit();
}

void Call_event_handler::camera_print_param_action()
{
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        get_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    {
        Multiple_stereo_camera* ms_camera = 
            Nvindex_AppData::instance()->get_user_interaction(0)->get_multiple_stereo_camera();
        if (ms_camera->get_use_ortho_camera())
        {
            Camera_tool::print_camera_param(
                ms_camera->get_ortho_camera_tag(), dice_transaction.get());
        }
        else
        {
            Camera_tool::print_camera_param(
                ms_camera->get_main_base_camera_tag(), dice_transaction.get());
        }
    }
    dice_transaction->commit();
}

void Call_event_handler::get_user_current_colormap_value(std::vector<mi::math::Color_struct> & color_table)
{
    // response_arguments.set_value(const char* key, const char* value)
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        get_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<nv::index::ISession>(m_session_tag));

    if (Nvindex_AppData::instance()->get_current_colormap_index() < get_number_of_colormap())
    {
        mi::neuraylib::Tag colormap_tag = get_colormap_tag(
            Nvindex_AppData::instance()->get_current_colormap_index());
        mi::base::Handle<const nv::index::IColormap> colormap(
            dice_transaction->access<nv::index::IColormap>(colormap_tag));
        if (colormap)
        {
            const mi::Uint32 nb_colormap_entries = static_cast<mi::Uint32>(colormap->get_number_of_entries());
            color_table.resize(nb_colormap_entries);
            for(mi::Uint32 i = 0; i < nb_colormap_entries; ++i){
                color_table[i] = colormap->get_color(i);
            }
        }
    }

    dice_transaction->commit();
}

void Call_event_handler::get_user_current_colormap(Arguments_set_wrapper& response_arguments)
{
    std::vector<mi::math::Color_struct>  color_table;
    this->get_user_current_colormap_value(color_table);
    const mi::Uint32 nb_colormap_entries = color_table.size();
    std::ostringstream os;
    for(mi::Uint32 i = 0; i < nb_colormap_entries; ++i){
        mi::math::Color const & entry = color_table[i];
        os << entry[0] << " " << entry[1] << " " << entry[2] << " " << entry[3];
        if(i < (nb_colormap_entries - 1)){
            os << " ";
        }
    }
    response_arguments.set_value("values", os.str().c_str());

    // Retrieve and output the domain of the colormap
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        get_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<nv::index::ISession>(m_session_tag));

    if (Nvindex_AppData::instance()->get_current_colormap_index() < get_number_of_colormap())
    {
        mi::neuraylib::Tag colormap_tag = get_colormap_tag(
            Nvindex_AppData::instance()->get_current_colormap_index());
        mi::base::Handle<const nv::index::IColormap> colormap(
            dice_transaction->access<nv::index::IColormap>(colormap_tag));
        if (colormap)
        {
            mi::math::Vector<mi::Float32, 2> domain;
            colormap->get_domain(domain.x, domain.y);
            {
                std::ostringstream os;
                os << domain.x;
                response_arguments.set_value("domainFirst", os.str().c_str());
            }

            {
                std::ostringstream os;
                os << domain.y;
                response_arguments.set_value("domainLast", os.str().c_str());
            }

            {
                std::ostringstream os;
                os << static_cast<mi::Uint32>(colormap->get_domain_boundary_mode()) + 1;
                response_arguments.set_value("domainBoundaryMode", os.str().c_str());
            }

            if (   dice_transaction->get_privacy_level(colormap_tag)
                   == static_cast<mi::Sint32>(get_scope()->get_privacy_level()))
            {
                response_arguments.set_value("isLocalized", "1");
            }
        }
    }

    dice_transaction->commit();
}


void Call_event_handler::get_user_current_colormap_all(Arguments_set_wrapper& response_arguments)
{
    // response_arguments.set_value(const char* key, const char* value)
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        get_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<nv::index::ISession>(m_session_tag));
    {
        std::ostringstream os;
        const mi::Uint32 nb_colormaps = get_number_of_colormap();
        for(mi::Uint32 i=0; i<nb_colormaps; ++i)
        {
            os << i;
            if(i<nb_colormaps-1)
                os << "|";
        }
        response_arguments.set_value("names", os.str().c_str());
    }

    {
        std::ostringstream os;
        os << Nvindex_AppData::instance()->get_current_colormap_index();
        response_arguments.set_value("cm_id", os.str().c_str());
    }

    const std::vector<mi::neuraylib::Tag> volume_tag_vec = Nvindex_AppData::instance()->get_volume_tag_vec();

    if (!volume_tag_vec.empty())
    {
        mi::base::Handle<const nv::index::IRegular_volume> volume_scene_element(
            dice_transaction->access<const nv::index::IRegular_volume>(volume_tag_vec.at(0)));
        assert(volume_scene_element.is_valid_interface());

        const mi::neuraylib::Tag& volume_current_colormap_tag = volume_scene_element->assigned_colormap();
        const mi::Uint32 nb_colormaps = get_number_of_colormap();
        mi::Uint32 volume_current_colormap_id = 0;
        for(mi::Uint32 i=0; i<nb_colormaps; ++i)
        {
            if(get_colormap_tag(i) == volume_current_colormap_tag)
            {
                volume_current_colormap_id = i;
                break;          // found
            }
        }
        {
            std::ostringstream os; os << (mi::Uint32(volume_current_colormap_id));
            response_arguments.set_value("seismic_cm_id", os.str().c_str());
        }
    }

    dice_transaction->commit();
}

void Call_event_handler::set_user_colormap(Arguments_get_wrapper& user_arguments)
{
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        get_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<nv::index::ISession>(m_session_tag));
    const std::string cm_id  = user_arguments.get_value("id");
    mi::Uint32 colormap_id = nv::index_common::get_uint32(cm_id);

    mi::neuraylib::Tag colormap_tag = get_colormap_tag(colormap_id);
    {
        mi::base::Handle<nv::index::IColormap> colormap(
            dice_transaction->edit<nv::index::IColormap>(colormap_tag));
        if (colormap)
        {

            const std::string values  = user_arguments.get_value("values");
            std::istringstream is (values);

            const mi::Uint32 nb_colormap_entries = static_cast<mi::Uint32>(colormap->get_number_of_entries());
            for(mi::Uint32 i=0; i<nb_colormap_entries; ++i)
            {
                // Leave unspecified colormap entries unchanged
                if (!is)
                    break;

                mi::math::Color entry;
                is >> entry[0] >> entry[1] >> entry[2] >> entry[3];
                colormap->set_color(i, entry);
            }

            // Set the domain of the colormap
            const char* domain_first = user_arguments.get_value("domainFirst");
            const char* domain_last  = user_arguments.get_value("domainLast");
            if (domain_first != 0 && !std::string(domain_first).empty() &&
                domain_last  != 0 && !std::string(domain_last).empty())
            {
                colormap->set_domain(get_float32(domain_first), get_float32(domain_last));
            }

            const char* mode_str = user_arguments.get_value("domainBoundaryMode");
            if (mode_str != 0)
            {
                nv::index::IColormap::Domain_boundary_mode mode = nv::index::IColormap::CLAMP_TO_EDGE;
                if (std::string(mode_str) == "2")
                {
                    mode = nv::index::IColormap::CLAMP_TO_TRANSPARENT;
                }
                colormap->set_domain_boundary_mode(mode);
            }
        }
    }

    dice_transaction->commit();
}

//----------------------------------------------------------------------
void Call_event_handler::export_heightfield_by_gui(Arguments_get_wrapper& user_arguments)
{
    nv::index_common::String_dict* p_app_proj = Nvindex_AppData::instance()->peek_app_proj();
    assert(p_app_proj != 0);

    const mi::Sint32 heightfield_gui_idx = get_sint32_userarg(user_arguments, "export_heightfield_id");
    if(heightfield_gui_idx < 0){
        ERROR_LOG << "Cannot export a heightfield: invalid heightfield_gui_idx [" << heightfield_gui_idx << "].";
        return;
    }

    const std::string heightfield_name = get_string_userarg(user_arguments, "export_heightfield_name");
    if((heightfield_gui_idx == 0) && heightfield_name.empty()){
        ERROR_LOG << "no heightfield in the scene or heightfield has no name, no export.";
        return;
    }
    if(heightfield_name.empty()){
        ERROR_LOG << "invalid heightfield name: no heightfield export.";
        return;
    }

    const std::string export_dirname = p_app_proj->get("app::export_heightfield::export_dirname", ".");
    std::string export_fname = "heightfield_export_file.bin";
    if(!export_dirname.empty()){
        export_fname = export_dirname + "/" + export_fname;
    }
    if(!heightfield_name.empty()){
        export_fname = export_dirname + "/" + heightfield_name + ".bin";
    }
    export_fname = p_app_proj->get("app::export_heightfield::export_filename", export_fname);

    // GUI's heightfield index order is based on heightfield_tag_vec, look it up.
    const std::vector<mi::neuraylib::Tag> heightfield_tag_vec =
        Nvindex_AppData::instance()->get_heightfield_tag_vec();
    if(heightfield_gui_idx >= static_cast<mi::Sint32>(heightfield_tag_vec.size())){
        ERROR_LOG << "Cannot export a heightfield: heightfield_gui_idx [" << heightfield_gui_idx
                  << "] is out of range.";
        return;
    }
    const mi::neuraylib::Tag heightfield_tag = heightfield_tag_vec.at(heightfield_gui_idx);
    assert(heightfield_tag.is_valid());
    INFO_LOG << "Exporting a heightfield data [" << heightfield_name << "] to [" << export_fname << "]...";

    // export the heightfield
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        get_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    {
        // get heightfield's bounding box and use it as the exporting region
        const mi::math::Bbox<mi::Sint32, 3> ijk_roi_bbox_uint =
            get_heightfield_roi_bbox(heightfield_tag, dice_transaction.get());
        // Note: you could set the heightfield bbox manually, too.
        const mi::math::Bbox<mi::Sint32, 2> ij_bbox(ijk_roi_bbox_uint.min.x, ijk_roi_bbox_uint.min.y,
                                                    ijk_roi_bbox_uint.max.x, ijk_roi_bbox_uint.max.y);

        if(export_heightfield_data(heightfield_tag, ij_bbox,
                                   export_fname, m_session_tag, dice_transaction.get()))
        {
            INFO_LOG << "Exported a heightfield data [" << export_fname << "] with bbox: "
                     << ij_bbox;
        }
        else{
            ERROR_LOG << "Failed exporting a heightfield data [" << export_fname << "].";
        }
    }
    dice_transaction->commit();
}

//----------------------------------------------------------------------
void Call_event_handler::get_export_heightfield_name_list(Arguments_set_wrapper& response_arguments)
{
    std::vector< std::string > heightfield_name_vec;
    mi::base::Handle<mi::neuraylib::IScope> scope = get_scope();
    Nvindex_AppData::instance()->get_heightfield_name_vec(scope, heightfield_name_vec);

    // set the list
    std::string const heightfield_name_gui = concat_string_with_separator(heightfield_name_vec, "|");
    response_arguments.set_value("names", heightfield_name_gui.c_str());

    // set exporting heightfield selection
    std::string heightfield_idx_str =
        Nvindex_AppData::instance()->peek_app_proj()->get("app::gui::export_heightfield_id", "-1");
    mi::Sint32 heightfield_idx = nv::index_common::get_sint32(heightfield_idx_str);
    if(heightfield_idx >= static_cast<mi::Sint32>(heightfield_name_vec.size())){
        // heightfield index has become invalid (maybe the
        // heightfield has been removed from the scene.)
        heightfield_idx_str = "-1";
    }
    response_arguments.set_value("export_heightfield_id", heightfield_idx_str.c_str());
}

//----------------------------------------------------------------------
    void Call_event_handler::set_export_heightfield_index_change(Arguments_get_wrapper& user_arguments)
    {
        const std::string heightfiled_idx_str = get_string_userarg(user_arguments, "export_heightfield_id");
        Nvindex_AppData::instance()->peek_app_proj()->insert("app::gui::export_heightfield_id", heightfiled_idx_str);
    }

//----------------------------------------------------------------------
void Call_event_handler::export_volume_by_gui(Arguments_get_wrapper& user_arguments)
{
    nv::index_common::String_dict* p_app_proj = Nvindex_AppData::instance()->peek_app_proj();
    assert(p_app_proj != 0);

    const mi::Sint32 volume_gui_idx = get_sint32_userarg(user_arguments, "export_volume_id");
    if(volume_gui_idx < 0){
        ERROR_LOG << "Cannot export a volume: invalid volume_gui_idx [" << volume_gui_idx << "].";
        return;
    }
    const std::string volume_name = get_string_userarg(user_arguments, "export_volume_name");
    if((volume_gui_idx == 0) && volume_name.empty()){
        ERROR_LOG << "no volume in the scene or volume has no name, no export.";
        return;
    }
    if(volume_name.empty()){
        ERROR_LOG << "invalid volume name: no volume export.";
        return;
    }

    const std::string export_dirname = p_app_proj->get("app::export_volume::export_dirname", ".");
    std::string export_fname = "volume_export_file.bin";
    if(!export_dirname.empty()){
        export_fname = export_dirname + "/" + export_fname;
    }
    if(!volume_name.empty()){
        export_fname = export_dirname + "/" + volume_name + ".bin";
    }
    export_fname = p_app_proj->get("app::export_volume::export_filename", export_fname);

    // GUI's volume index order is based on volume_tag_vec, look it up.
    const std::vector<mi::neuraylib::Tag> volume_tag_vec = Nvindex_AppData::instance()->get_volume_tag_vec();
    if(volume_gui_idx >= static_cast<mi::Sint32>(volume_tag_vec.size())){
        ERROR_LOG << "Cannot export a volume: volume_gui_idx [" << volume_gui_idx << "] is out of range.";
        return;
    }
    const mi::neuraylib::Tag volume_tag = volume_tag_vec.at(volume_gui_idx);
    assert(volume_tag.is_valid());
    INFO_LOG << "Exporting a volume data [" << volume_name << "] to [" << export_fname << "]...";

    // export the volume
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        get_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    {
        // get volume's bounding box and use it as the exporting region
        const mi::math::Bbox<mi::Sint32, 3> volume_ijk_roi_bbox_sint =
            get_volume_roi_bbox(volume_tag, dice_transaction.get());

        if(export_volume_data(volume_tag, volume_ijk_roi_bbox_sint,
                              export_fname, m_session_tag, dice_transaction.get()))
        {
            INFO_LOG << "Exported a volume data [" << export_fname << "] with bbox: "
                     << volume_ijk_roi_bbox_sint;
        }
        else{
            ERROR_LOG << "Failed exporting a volume data [" << export_fname << "].";
        }
    }
    dice_transaction->commit();
}

//----------------------------------------------------------------------
void Call_event_handler::get_export_volume_name_list(Arguments_set_wrapper& response_arguments)
{
    std::vector< std::string > volume_name_vec;
    mi::base::Handle<mi::neuraylib::IScope> scope = get_scope();
    Nvindex_AppData::instance()->get_volume_name_vec(scope, volume_name_vec);

    // set the list
    std::string const volume_name_gui = concat_string_with_separator(volume_name_vec, "|");
    response_arguments.set_value("names", volume_name_gui.c_str());

    // set exporting volume selection
    // If there is no volume in the scene, no selection. Otherwise the first one.
    std::ostringstream os;
    os << (volume_name_vec.empty() ? -1 : 0);
    response_arguments.set_value("export_volume_id", os.str().c_str());
}

//----------------------------------------------------------------------
void Call_event_handler::attrgen_copy_volume_by_gui(Arguments_get_wrapper& user_args)
{
    const mi::Sint32  volume_gui_idx =
        get_sint32_userarg(user_args, "id_attrgen_copy_volume_source_index");
    const mi::Float32 x_translation =
        get_float32_userarg(user_args, "id_attrgen_copy_volume_x_translation_position");
    const mi::Float32 y_translation =
        get_float32_userarg(user_args, "id_attrgen_copy_volume_y_translation_position");
    const mi::Float32 z_translation =
        get_float32_userarg(user_args, "id_attrgen_copy_volume_z_translation_position");

    m_gui_state_dict.insert("attrgen_copy_source_volume_index", nv::index_common::to_string(volume_gui_idx));

    // get the copy source tag
    const std::vector<mi::neuraylib::Tag> volume_tag_vec = Nvindex_AppData::instance()->get_volume_tag_vec();
    if (volume_gui_idx >= static_cast<mi::Sint32>(volume_tag_vec.size()))
    {
        ERROR_LOG << "Cannot copy a volume: volume_gui_idx [" << volume_gui_idx << "] is out of range.";
        return;
    }
    const mi::neuraylib::Tag volume_tag = volume_tag_vec.at(volume_gui_idx);
    assert(volume_tag.is_valid());

    // get the copy destination group
    Scene_traverse_get_parent parent_finder(volume_tag);
    bool is_success = false;
    {
        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            m_irc.get_current_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
        assert(dice_transaction.is_valid_interface());
        {
            is_success = get_scene_status(&parent_finder, m_irc.m_session_tag, dice_transaction.get());
        }
        dice_transaction->commit();
    }

    if (!is_success)
    {
        ERROR_LOG << "Fail to traverse the scene.";
        return;
    }
    mi::neuraylib::Tag parent_tag = parent_finder.get_parent_tag();
    if (!parent_tag.is_valid())
    {
        ERROR_LOG << "Tag (" << parent_tag << ") has no parent.";
        return;

    }
    const mi::math::Vector<mi::Float32, 3> translation_vec(x_translation, y_translation, z_translation);
    INFO_LOG << "Copy a volume [" << volume_tag << "] under [" << parent_tag << "] with translation: "
             << translation_vec;

    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        get_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    std::string err_mes;
    mi::neuraylib::Tag new_copied_tag;
    {
        new_copied_tag = user_command_copy_volume_sub(
            m_session_tag,
            volume_tag,
            parent_tag,
            translation_vec,
            dice_transaction.get(),
            err_mes);
    }

    if(new_copied_tag.is_valid()){
        INFO_LOG << "created copy tag: " << new_copied_tag;
        dice_transaction->commit();
    }
    else{
        ERROR_LOG << "Failed to copy a volume: " << err_mes;
        dice_transaction->abort();
    }
}

//----------------------------------------------------------------------
void Call_event_handler::attrgen_filter_volume_by_gui(Arguments_get_wrapper& user_args)
{
    const mi::Sint32  volume_gui_idx =
        get_sint32_userarg(user_args, "id_attrgen_filter_volume_index");
    const mi::Sint32  filter_gui_idx =
        get_sint32_userarg(user_args, "id_attrgen_filter_volume_filter_index");

    ERROR_LOG << "attrgen_filter_volume_by_gui: volume_gui_idx: " << volume_gui_idx;
    ERROR_LOG << "attrgen_filter_volume_by_gui: id_attrgen_filter_volume_filter_index: " << filter_gui_idx;

    m_gui_state_dict.insert("attrgen_filter_source_volume_index", nv::index_common::to_string(volume_gui_idx));
    m_gui_state_dict.insert("attrgen_filter_name_index",          nv::index_common::to_string(filter_gui_idx));

    // get the filtering volume tag
    const std::vector<mi::neuraylib::Tag> volume_tag_vec = Nvindex_AppData::instance()->get_volume_tag_vec();
    if(volume_gui_idx >= static_cast<mi::Sint32>(volume_tag_vec.size())){
        ERROR_LOG << "Cannot filter a volume: volume_gui_idx [" << volume_gui_idx << "] is out of range.";
        return;
    }
    const mi::neuraylib::Tag volume_tag = volume_tag_vec.at(volume_gui_idx);
    assert(volume_tag.is_valid());

    // get the filter
    if(!Volume_data_filter::is_valid_filter_type(filter_gui_idx)){
        ERROR_LOG << "Unknown filter index: " << filter_gui_idx;
        return;
    }
    const std::string filter_name = Volume_data_filter::get_filter_name(filter_gui_idx);

    // FIXME: parallel numbers
    const std::string pcount_str = Nvindex_AppData::instance()->peek_app_proj()->
        get("app::distributed_seismic_filter_filter_job::parallel_count", "1");
    const mi::Sint32 pcount = get_sint32(pcount_str);

    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        get_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    bool is_success = false;
    {
        is_success = volume_inset_filter_amplitude_values(
            m_session_tag, volume_tag, filter_name, pcount, dice_transaction.get());
    }
    if(is_success){
        INFO_LOG << "Volume_filter_amplitude_values ... done.";
        dice_transaction->commit();
    }
    else{
        ERROR_LOG << "Volume_filter_amplitude_values ... failed.";
        dice_transaction->abort();
    }
}

//----------------------------------------------------------------------
void Call_event_handler::attrgen_edit_volume_by_gui(Arguments_get_wrapper& user_args)
{
    const mi::Sint32  volume_gui_idx =
        get_sint32_userarg(user_args, "id_attrgen_edit_volume_index");
    ERROR_LOG << "attrgen_edit_volume_by_gui: volume_gui_idx: " << volume_gui_idx;

    m_gui_state_dict.insert("attrgen_edit_volume_index", nv::index_common::to_string(volume_gui_idx));

    // get the volume tag to change the amplitude value
    const std::vector<mi::neuraylib::Tag> volume_tag_vec = Nvindex_AppData::instance()->get_volume_tag_vec();
    if(volume_gui_idx >= static_cast<mi::Sint32>(volume_tag_vec.size())){
        ERROR_LOG << "Cannot edit a volume: volume_gui_idx [" << volume_gui_idx << "] is out of range.";
        return;
    }
    const mi::neuraylib::Tag volume_tag = volume_tag_vec.at(volume_gui_idx);
    assert(volume_tag.is_valid());

    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        get_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());

    const mi::Uint8 amplitude_value = 128;
    const bool is_host_assign_mode = true; // show the cluster assignment in this case
    const bool is_success =
        volume_set_amplitude_value(
            is_host_assign_mode,
            amplitude_value,
            m_session_tag,
            volume_tag,
            dice_transaction.get());
    if(is_success){
        INFO_LOG << "attrgen_edit_volume_index: set host id for the assignment.";
        dice_transaction->commit();
    }
    else{
        ERROR_LOG << "attrgen_edit_volume_index: fail to set host id.";
        dice_transaction->abort();
    }
}

//----------------------------------------------------------------------
void Call_event_handler::get_attrgen_gui_section_info(Arguments_set_wrapper& response_args)
{
    std::vector< std::string > generate_type_name_vec;
    mi::base::Handle<mi::neuraylib::IScope> scope = get_scope();

    // copy volume info
    {
        std::vector< std::string > volume_name_vec;
        Nvindex_AppData::instance()->get_volume_name_vec(scope, volume_name_vec);

        std::string const volume_name_gui = concat_string_with_separator(volume_name_vec, "|");
        response_args.set_value("attrgen_copy_source_volume_names", volume_name_gui.c_str());

        // set selection
        response_args.set_value("attrgen_copy_source_volume_index",
                                m_gui_state_dict.
                                get("attrgen_copy_source_volume_index", "0").c_str());
    }

    // filter volume info
    {
        // volume
        std::vector< std::string > volume_name_vec;
        mi::base::Handle<mi::neuraylib::IScope> scope = get_scope();
        Nvindex_AppData::instance()->get_volume_name_vec(scope, volume_name_vec);

        std::string const volume_name_gui = concat_string_with_separator(volume_name_vec, "|");
        response_args.set_value("attrgen_filter_source_volume_names", volume_name_gui.c_str());
        // set selection
        response_args.set_value("attrgen_filter_source_volume_index",
                                m_gui_state_dict.
                                get("attrgen_filter_source_volume_index", "0").c_str());

        // filter
        std::vector< std::string > generate_filter_name_vec;
        for(mi::Sint32 i = 0; i < Volume_data_filter::VDFT_Count; ++i){
            generate_filter_name_vec.push_back(Volume_data_filter::get_filter_name(i));
        }
        std::string const generate_filter_name_gui =
            concat_string_with_separator(generate_filter_name_vec, "|");
        response_args.set_value("attrgen_filter_name_list", generate_filter_name_gui.c_str());
        // set selection
        response_args.set_value("attrgen_filter_name_index",
                                m_gui_state_dict.get("attrgen_filter_name_index", "0").
                                c_str());
    }

    // edit volume
    {
        // volume
        std::vector< std::string > volume_name_vec;
        Nvindex_AppData::instance()->get_volume_name_vec(scope, volume_name_vec);

        std::string const volume_name_gui = concat_string_with_separator(volume_name_vec, "|");
        response_args.set_value("attrgen_edit_volume_names", volume_name_gui.c_str());

        // set selection
        response_args.set_value("attrgen_edit_volume_index",
                                m_gui_state_dict.
                                get("attrgen_edit_volume_index", "0").c_str());
    }
}

//----------------------------------------------------------------------
void Call_event_handler::rtc_edit_param_buffer_by_gui(Arguments_get_wrapper& user_args)
{
    const std::vector<mi::neuraylib::Tag> rtc_buf_tag_vec = Nvindex_AppData::instance()->get_rtc_param_buffer_tag_vec();

    if (rtc_buf_tag_vec.size() > 0u)
    {
        const mi::Float32 id_rtc_param_plane_pos_value =
            get_float32_userarg(user_args, "id_rtc_param_plane_pos_value");

        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            get_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
        assert(dice_transaction.is_valid_interface());

        RTC_parameter_buffer_manip::instance()->edit_parameter_buffer_instance(
            rtc_buf_tag_vec[0], id_rtc_param_plane_pos_value, dice_transaction.get());
        
        dice_transaction->commit();
    }
}

//----------------------------------------------------------------------
void Call_event_handler::get_performance_logging_state_by_gui(Arguments_set_wrapper& response_args)
{
    const Performance_logger_app_state & app_perf_state = Nvindex_AppData::instance()->peek_log_state();
    const mi::Sint32 opt_count = 4;
    std::string perf_key[opt_count] = {
        "is_monitor_performance_values",
        "is_log_performance_global",
        "is_log_performance_per_host",
        "is_log_performance_per_span",
    };
    bool perf_state[opt_count] = {
        app_perf_state.is_monitoring_on(),
        app_perf_state.is_logging_on(Performance_logger_base::LK_System_global),
        app_perf_state.is_logging_on(Performance_logger_base::LK_Per_host),
        app_perf_state.is_logging_on(Performance_logger_base::LK_Per_span),
    };
    assert((sizeof(perf_state)/sizeof(bool)) == opt_count);

    for (mi::Sint32 i = 0; i < opt_count; ++i)
    {
        std::ostringstream os;
        os << (perf_state[i] ? "1" : "0");
        response_args.set_value(perf_key[i].c_str(), os.str().c_str());
    }
}

//----------------------------------------------------------------------
void Call_event_handler::render_viewing_scenario(
    mi::Uint32                                                id,
    const mi::base::Handle<mi::neuraylib::IDice_transaction>& dice_transaction)
{
    // Create a new NVIDIA IndeX canvas and render into that canvas.
    // NOTE: Creating and deleting a canvas for every frame should generally be for performance reasons.
    //       Here, it just simplifies the source code of the reference viewer!
    mi::base::Handle<nv::index::IIndex_canvas> rendering_canvas(new Span_renderer_no_gl());
    (static_cast<Span_renderer_IF*>(rendering_canvas.get()))->set_background_color(
        Nvindex_AppData::instance()->m_background_color);
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
        m_iindex_if->create_rendering_interface());

    // NVIDIA IndeX rendering call
    mi::base::Handle<nv::index::IFrame_results> frame_results(
        index_rendering_interface->render(
            m_session_tag,
            rendering_canvas.get(),
            dice_transaction.get(),
            NULL,
            frame_info_callbacks.get(),
            true,
            NULL));

    const mi::base::Handle<nv::index::IError_set> err_set(frame_results->get_error_set());
    if (err_set->any_errors())
    {
        std::ostringstream os;

        const mi::Uint32 nb_err = err_set->get_nb_errors();
        for (mi::Uint32 e = 0; e < nb_err; ++e)
        {
            if (e != 0) os << '\n';
            const mi::base::Handle<nv::index::IError> err(err_set->get_error(e));
            os << err->get_error_string();
        }

        ERROR_LOG << "IIndex_rendering rendering call failed with the following error(s): " << '\n'
                  << os.str();
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
            for(mi::Uint32 j=0; j<nb_operations; ++j)
            {
                const nv::index::IBalancing_operation* op = balancing_ops[i]->get_operation(j);
                INFO_LOG << "Balance operation " << j << ": " << op->get_description();
            }
        }
    }

    // Prepare frame buffer to be passed into the video stream or for taking a snapshot image (tripple buffer)
    mi::base::Handle<Canvas> canvas(m_irc.m_canvas_buffers[1]->get_render_canvas());
    copy_result_pixel_to_canvas(static_cast<Span_renderer_IF*>(rendering_canvas.get()), canvas.get());

    // Tripple buffering for video-rendering swapping to accelerate the video encoding.
    // - Trigger internal buffer swapping
    m_irc.m_canvas_buffers[1]->rendering_finished();
}

//////////////////// Scene_setup_call_event_handler ////////////////////

mi::base::Lock Scene_setup_call_event_handler::s_scene_setup_lock;

Scene_setup_call_event_handler::Scene_setup_call_event_handler(
    Nvindex_rendering_context& irc,
    bool                       admin)
    : Call_event_handler(0, irc),
      m_admin(admin)
{
    String_dict* prj = Nvindex_AppData::instance()->peek_app_proj();

    // Optional directory where all the project files reside
    const std::string directory = prj->get("app::setup::scene_directory");

    // Iterate over the available scenes
    std::istringstream is_scenes(prj->get("app::setup::scenes"));
    std::string scene;
    while (std::getline(is_scenes, scene, ' '))
    {
        const std::string prefix = "app::setup::scenes::" + scene + "::";

        if (!get_bool(prj->get(prefix + "enabled", "true")))
            continue;

        m_scene_names.push_back(prj->get(prefix + "name", scene));
        m_scene_cluster_size_min.push_back(get_sint32(prj->get(prefix + "cluster_size_min", "0")));

        // Add directory to each listed project file
        std::istringstream is_files(prj->get(prefix + "files"));
        std::string file;
        std::string files;
        while (std::getline(is_files, file, ' '))
        {
            if (!files.empty())
                files += " ";
            if (!directory.empty() && !file.empty() && file[0] != '/')
                files += directory + "/";
            files += file;
        }
        m_scene_files.push_back(files);
    }

    if (m_scene_names.empty())
    {
        m_scene_names.push_back("Default");
        m_scene_files.push_back("");
        m_scene_cluster_size_min.push_back(0);
    }
}

    namespace {

    std::string quote(const std::string& s)
    {
        return '"' + escape_JSON(s) + '"';
    }

    } // namespace

bool Scene_setup_call_event_handler::handle(
    mi::rtmp::IConnection*      /*connection*/,
    const char*                 /*procedure_name*/,
    const mi::IData*            /*command_args*/,
    const mi::IData*            user_args,
    mi::IData**                 response_args)
{
    mi::base::Lock::Block block(&s_scene_setup_lock);

    mi::base::Handle<const mi::IMap> imap(user_args->get_interface<const mi::IMap>());
    if (!imap.is_valid_interface())
        return false;

    Call_event_arguments args(imap);

    std::string cmd;
    if (args.has("cmd"))
        cmd = args.get("cmd");

    String_dict* prj = Nvindex_AppData::instance()->peek_app_proj();

    if (cmd == "get_is_scene_initialized")
    {
        args.set("cmd", cmd);
        args.set("initialized", (m_irc.m_initialized ? "true" : "false"));
        args.set("scene_setup_done", (m_irc.m_scene_setup_done ? "true" : "false"));
        args.set("admin", (m_admin ? "true" : "false"));
    }
    if (cmd == "get_scene_options" && m_admin)
    {
        args.set("cmd", cmd);

        std::ostringstream os;
        os << "{";

        mi::Uint32 cluster_size = get_uint32(prj->get("app::cluster_size", "0"));
        mi::Uint32 cluster_size_min = get_uint32(prj->get("app::setup::cluster_size_min", "0"));
        mi::Uint32 cluster_size_max = get_uint32(prj->get("app::setup::cluster_size_max", "0"));
        if (prj->get("dice::network::mode", "OFF") == "OFF")
        {
            cluster_size = 1;
            cluster_size_min = 1;
            cluster_size_max = 1;
        }

        os << quote("cluster_size") << ":" << cluster_size << ","
           << quote("cluster_size_min") << ":" << cluster_size_min << ","
           << quote("cluster_size_max") << ":" << cluster_size_max << ","
           << quote("has_password") << ":"
           << quote(prj->get("dice::http_auth_user").empty() ? "false" : "true") << ","
           << quote("resolutions") << ":" << quote(prj->get("app::setup::resolutions")) << ","
           << quote("resolution") << ":" << quote(prj->get("index::canvas_resolution")) << ","
           << quote("video_codec") << ":" << quote(prj->get("dice::rtmp_video_streaming::video_codec")) << ","
           << quote("video_bitrate") << ":" << prj->get("dice::rtmp_video_streaming::video_bitrate") << ",";

        os << quote("scenes") << ":[";
        for (size_t i=0; i < m_scene_names.size(); ++i)
        {
            if (i > 0)
                os << ",";
            os << "{"
               << quote("name") << ":" << quote(m_scene_names[i]) << ","
               << quote("scene_idx") << ":" << i << ","
               << quote("cluster_size_min") << ":" << m_scene_cluster_size_min[i]
               << "}";
        }
        os << "]";

        os << "}";
        args.set("data", os.str());
    }
    else if (cmd == "setup_scene" && !m_irc.m_scene_setup_done && m_admin)
    {
        INFO_LOG << "Received scene setup from user";

        mi::Uint32 scene_idx = get_uint32(args.get("scene_idx"));
        if (scene_idx < m_scene_files.size())
        {
            // Additional project files
            INFO_LOG << "  Selected scene " << scene_idx << ": " << m_scene_names[scene_idx];
            std::istringstream is(m_scene_files[scene_idx]);
            std::string file;
            while (std::getline(is, file, ' '))
            {
                INFO_LOG << "  Loading extra project file '" << file << "'";
                String_dict extra_dict;
                load_application_project_file(file, extra_dict);
                prj->insert_all(extra_dict);

                // Also read IndeX debug configuration settings from the additional project file.
                // Note that not all settings will work here, some need to be specified in the
                // initial project file.
                const std::vector<std::pair<std::string, std::string> > config =
                    get_configuration_vector("index::debug_configuration::", extra_dict, true);
                if (!config.empty())
                {
                    mi::base::Handle<nv::index::IIndex_debug_configuration> index_debug_config(
                        m_iindex_if->get_api_component<nv::index::IIndex_debug_configuration>());
                    for (size_t i = 0; i < config.size(); ++i)
                    {
                        const std::string opt = config[i].first + "=" + config[i].second;
                        INFO_LOG << "  Extra IndeX debug configuration: [" << opt << "]";
                        if (index_debug_config->set_option(opt.c_str()))
                            ERROR_LOG << "Failed to set IndeX debug configuration: [" << opt << "]";
                    }
                }
            }
        }
        else
        {
            ERROR_LOG << "Received invalid scene index " << scene_idx;
        }

        // Selected resolution
        const std::string resolution_w = args.get("resolution_w");
        const std::string resolution_h = args.get("resolution_h");
        INFO_LOG << "  Setting resolution to " << resolution_w << "x" << resolution_h;
        prj->insert("index::canvas_resolution", resolution_w + " " + resolution_h);

        // Aspect ratio
        {
            std::ostringstream os;
            os << (get_float32(resolution_w) / get_float32(resolution_h));
            prj->insert("index::camera::aspect", os.str());
        }

        // Video codec
        const std::string video_codec = args.get("video_codec");
        if (!video_codec.empty())
        {
            INFO_LOG << "  Setting rtmp video codec to '" << video_codec << "'";
            prj->insert("dice::rtmp_video_streaming::video_codec", video_codec);
        }

        // Video bitrate
        const mi::Uint32 video_bitrate = get_uint32(args.get("video_bitrate"));
        if (video_bitrate > 0)
        {
            std::ostringstream os;
            os << (video_bitrate * 1000000); // Mbit -> bit
            INFO_LOG << "  Setting video bitrate to " << os.str();
            prj->insert("dice::rtmp_video_streaming::video_bitrate", os.str());
        }

        // Cluster size
        const mi::Uint32 cluster_size = get_uint32(args.get("cluster_size"));
        if (cluster_size > 0)
        {
            INFO_LOG << "  Setting minimum cluster size to " << cluster_size;
            prj->insert("app::cluster_size", args.get("cluster_size"));
        }

        // Guest password
        const std::string guest_password = args.get("guest_password");
        if (!guest_password.empty())
        {
            INFO_LOG << "  Setting guest password";
            prj->insert("dice::http_auth_guest_password", guest_password);
        }

        m_irc.m_scene_setup_done = true;
    }

    *response_args = args.generate_response(m_iindex_if);
    return true;
}
