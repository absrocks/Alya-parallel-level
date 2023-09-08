/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "user_interaction.h"

#include <nv/index/iconfig_settings.h>

#include "common/forwarding_logger.h"
#include "common/clock_pulse_generator.h"

#include "multiple_camera.h"
#include "nvindex_rendering_context.h"
#include "nvindex_appdata.h"
#include "pick_tool.h"
#include "scene_utility.h"


//----------------------------------------------------------------------
User_interaction::User_interaction()
    :
    m_last_press_mouse_x(-1),
    m_last_press_mouse_y(-1)
{
    // empty
}

//----------------------------------------------------------------------
User_interaction::~User_interaction()
{
    // empty
}

//----------------------------------------------------------------------
void User_interaction::initialize_by_dict(const nv::index_common::String_dict& dict)
{
    // Currently only initialize the examiner.
    get_examiner()->initialize_by_dict(dict);
}

//----------------------------------------------------------------------
Examiner_manipulator* User_interaction::get_examiner()
{
    return &m_examiner;
}

//----------------------------------------------------------------------
Multiple_stereo_camera* User_interaction::get_multiple_stereo_camera()
{
    return &m_multiple_camera;
}

//----------------------------------------------------------------------
void User_interaction::mouse_motion(
    mi::Sint32                        button,
    mi::Sint32                        mouse_x,
    mi::Sint32                        mouse_y,
    const mi::neuraylib::Tag&         scene_tag,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    assert(dice_transaction != 0);

    if (get_multiple_stereo_camera()->get_use_ortho_camera())
    {
        get_examiner()->mouse_motion_event(
            button,
            mouse_x,
            mouse_y,
            get_multiple_stereo_camera()->get_ortho_camera_tag(),
            scene_tag,
            dice_transaction);
    }
    else
    {
        get_examiner()->mouse_motion_event(
            button,
            mouse_x,
            mouse_y,
            get_multiple_stereo_camera()->get_main_base_camera_tag(),
            scene_tag,
            dice_transaction);
    }
}

//----------------------------------------------------------------------
mi::math::Vector<mi::Sint32, 2> User_interaction::mouse_button_action(
    mi::Sint32                         button,
    Examiner_manipulator::Mouse_action mouse_action,
    mi::Sint32                         mouse_x,
    mi::Sint32                         mouse_y,
    mi::neuraylib::IDice_transaction*  dice_transaction)
{
    assert(dice_transaction != 0);

    // WARN_LOG << "mouse: " << mouse_x << ", " << mouse_y;

    if (get_multiple_stereo_camera()->get_use_ortho_camera())
    {
        get_examiner()->mouse_press_event(
            button, 
            mouse_action, 
            mouse_x, 
            mouse_y,
            get_multiple_stereo_camera()->get_ortho_camera_tag(),
            dice_transaction);
    }
    else
    {
        get_examiner()->mouse_press_event(
            button, 
            mouse_action,
            mouse_x, 
            mouse_y,
            get_multiple_stereo_camera()->get_main_base_camera_tag(),
            dice_transaction);
    }

    // Negative pick point means no pick
    mi::math::Vector<mi::Sint32, 2> pick_position(-1, -1);

    if(mouse_action == Examiner_manipulator::Button_down)
    {
        m_last_press_mouse_x = mouse_x;
        m_last_press_mouse_y = mouse_y;
    }
    else if(mouse_action == Examiner_manipulator::Button_up)
    {
        // When the locations of push and release are the same, we pick.
        if((mouse_x == m_last_press_mouse_x) && (mouse_y == m_last_press_mouse_y))
        {
            const mi::math::Vector<mi::Sint32,2> window_resolution =
                get_examiner()->get_main_window_resolution();

            pick_position.x = mouse_x;
            // Top left origin to bottom left origin coordinate
            pick_position.y = window_resolution.y - mouse_y;
        }

        m_last_press_mouse_x = -1;
        m_last_press_mouse_y = -1;
    }

    return pick_position;
}

//----------------------------------------------------------------------
mi::neuraylib::Tag User_interaction::issue_pick_command(
    Nvindex_rendering_context*             irc_ref,
    bool                                   is_multi_view_mode,
    const mi::math::Vector<mi::Sint32, 2>& pick_position)
{
    return Pick_tool::issue_pick_command(
        irc_ref,
        is_multi_view_mode,
        pick_position,
        get_examiner());
}

//----------------------------------------------------------------------
void User_interaction::mouse_wheel_action(
    Nvindex_rendering_context*        irc_ref,
    const mi::neuraylib::Tag&         scene_tag,
    mi::Sint32                        mouse_x,
    mi::Sint32                        mouse_y,
    mi::Sint32                        delta,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    if (get_multiple_stereo_camera()->get_use_ortho_camera())
    {
        get_examiner()->
            mouse_wheel_event(
                delta,
                get_multiple_stereo_camera()->get_ortho_camera_tag(),
                scene_tag,
                dice_transaction);
    }
    else
    {
        get_examiner()->
            mouse_wheel_event(
                delta,
                get_multiple_stereo_camera()->get_main_base_camera_tag(),
                scene_tag,
                dice_transaction);
    }
}

//----------------------------------------------------------------------
void User_interaction::key_down_action(
    Nvindex_rendering_context* irc_ref,
    const mi::Sint32           keycode, 
    const bool                 is_shift_on,
    const bool                 is_ctrl_on,
    const bool                 is_alt_on,
    const bool                 is_meta_on)
{
    // (Any key with a special key || only a special key) handling

    Nvindex_AppData* appdata = Nvindex_AppData::instance();

    if (is_shift_on)
    {
        get_examiner()->set_manipulation_mode(Examiner_manipulator::EM_Zoom);
    }

    if (is_ctrl_on)
    {
        get_examiner()->set_manipulation_mode(Examiner_manipulator::EM_Pan);
    }

    if (is_alt_on)
    {
        // Set Alt key
    }

    switch (keycode)
    {
    case 27: // Escape key
    {
        shutdown_server_action(irc_ref);
        break;
    }

    case '0':
        if (is_ctrl_on)
        {
            appdata->set_enable_opengl_z_buffer_for_rendering(
                !appdata->is_enable_opengl_z_buffer_for_rendering());
        }
        break;

    case '1':
        if (is_ctrl_on)
        {
            appdata->m_show_rendering_statistics_overlay
                = !appdata->m_show_rendering_statistics_overlay;
        }
        break;

    case '2':
        if (is_ctrl_on)
        {
            appdata->m_show_compositing_statistics_overlay
                = !appdata->m_show_compositing_statistics_overlay;
        }
        break;

    case '3':
        if (is_ctrl_on)
        {
            appdata->m_show_all_performance_statistics_overlay
                = !appdata->m_show_all_performance_statistics_overlay;
        }
        break;

    case '4':
        if(is_ctrl_on)
        {
            appdata->set_enabled_recent_n_stat(
                !appdata->is_enabled_recent_n_stat());
            bool const is_enabled_n_stat = appdata->is_enabled_recent_n_stat();
            INFO_LOG << "cumulative_fps_statistics (reset with shift when toggle to off): "
                     << (is_enabled_n_stat ? "on" : "off");
            if (!is_enabled_n_stat)
            {
                // deactivated
                appdata->set_stat_recent_n_max_fps(-1.0);
                appdata->set_stat_recent_n_ave_fps(-1.0);

                if (is_shift_on)
                {
                    // reset the statistics
                    appdata->peek_cumulative_stat()->clear();
                    INFO_LOG << "reset cumulative_fps_statistics.";
                }
            }
        }
        break;

    case '5':
        if (is_ctrl_on)    // Set rotation Center to a volume center.
        {
            set_rotation_center_to_volume_roi_center(irc_ref);
        }
        break;

    case '6':
        if (is_ctrl_on)
        {
            INFO_LOG << "No command is assigned to key 'Ctrl-6'";
        }
        break;

    case '7':
        if (is_ctrl_on)
        {
            nv::index_common::String_dict* p_app_proj = appdata->peek_app_proj();
            assert(p_app_proj != 0);

            static mi::Uint32 color_idx = 0;
            if(color_idx==1) color_idx=0;
            else if(color_idx==0) color_idx=1;

            INFO_LOG << "Changing background color to " << color_idx;

            if(color_idx==1)
            {
                Nvindex_AppData::instance()->m_background_color
                    = nv::index_common::get_color(p_app_proj->get("app::canvas::background_color_1", "1 1 1 1"));
            }
            else if(color_idx==0)
            {
                Nvindex_AppData::instance()->m_background_color
                    = nv::index_common::get_color(p_app_proj->get("app::canvas::background_color", "1 1 1 1"));
            }
        }
        break;

    case 'P':  // CLOCK STOP/RESUME
    {
        if (is_ctrl_on)
        {
            mi::base::Handle<nv::index::IClock_pulse_generator> iclock(
                irc_ref->m_iindex_session->get_clock());
            if (iclock)
            {
                mi::base::Handle<nv::index_common::Clock_pulse_generator> clock(
                    iclock->get_interface<nv::index_common::Clock_pulse_generator>());
                if (clock)
                {
                    if (clock->is_running())
                    {
                        clock->stop();
                    }
                    else
                    {
                        clock->resume();
                    }
                }
            }
        }
        break;
    }

    case 'F': // CLOCK PLAYBACK
    {
        if (is_ctrl_on)
        {
            mi::base::Handle<nv::index::IClock_pulse_generator> iclock(
                irc_ref->m_iindex_session->get_clock());
            if (iclock)
            {
                mi::base::Handle<nv::index_common::Clock_pulse_generator> clock(
                    iclock->get_interface<nv::index_common::Clock_pulse_generator>());
                if (clock)
                {
                    if(clock->get_playback() == nv::index_common::Clock_pulse_generator::REWIND)
                    {
                        clock->set_playback(nv::index_common::Clock_pulse_generator::FORWARD);
                    }
                    else
                    {
                        clock->set_playback(nv::index_common::Clock_pulse_generator::REWIND);
                    }
                }
            }
        }
        break;
    }

    case 'C': // CLOCK REPEAT  //FIXME: This collides with Ctrl+C for copying to clipboard
    {
        if (is_ctrl_on)
        {
            mi::base::Handle<nv::index::IClock_pulse_generator> iclock(
                irc_ref->m_iindex_session->get_clock());
            if (iclock)
            {
                mi::base::Handle<nv::index_common::Clock_pulse_generator> clock(
                    iclock->get_interface<nv::index_common::Clock_pulse_generator>());
                if (clock)
                {
                    clock->set_repeat(!clock->get_repeat());
                }
            }
        }
        break;
    }

    case 'S':
    {
        if (is_ctrl_on)
        {
            if (is_shift_on)
            {
                appdata->m_statistics_overlay_large
                    = !appdata->m_statistics_overlay_large;
                appdata->m_show_statistics_overlay = true;
            }
            else
            {
                appdata->m_show_statistics_overlay
                    = !appdata->m_show_statistics_overlay;
            }
        }
        break;
    }

    case 'T':
    {
        if (is_ctrl_on)
        {
            mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
                irc_ref->get_current_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
            assert(dice_transaction.is_valid_interface());
            {
                mi::base::Handle<const nv::index::ISession> session(
                    dice_transaction->access<const nv::index::ISession>(irc_ref->m_session_tag));

                mi::base::Handle<nv::index::IConfig_settings> configs(
                    dice_transaction->edit<nv::index::IConfig_settings>(session->get_config()));

                mi::Uint32 time_step          = configs->get_current_timestep();
                time_step++;
                const mi::Uint32 nb_time_step = configs->get_nb_timesteps();
                if(time_step==nb_time_step)
                {
                    time_step = 0;
                }
                INFO_LOG << "Time step (+): " << time_step;
                configs->set_current_timestep(time_step);
            }
            dice_transaction->commit();
        }
        break;
    }

    case 't':
    {
        if (is_ctrl_on)
        {
            mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
                irc_ref->get_current_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
            assert(dice_transaction.is_valid_interface());
            {
                mi::base::Handle<const nv::index::ISession> session(
                    dice_transaction->access<const nv::index::ISession>(irc_ref->m_session_tag));

                mi::base::Handle<nv::index::IConfig_settings> configs(
                    dice_transaction->edit<nv::index::IConfig_settings>(session->get_config()));

                mi::Sint32 time_step          = configs->get_current_timestep();
                time_step--;
                const mi::Uint32 nb_time_step = configs->get_nb_timesteps();
                if (time_step==-1)
                {
                    time_step = nb_time_step-1;
                }
                INFO_LOG << "Time step (-): " << time_step; 
                configs->set_current_timestep(time_step);
            }
            dice_transaction->commit();
        }
        break;
    }

    case 'R':               // reset the pan reference point
    {
        if (is_ctrl_on && is_shift_on)
        {
            INFO_LOG << "Invalidate the pan pick point.";
            get_examiner()->set_pan_reference_point_valid(false);
        }
        break;
    }

    case 'Z':               // zoom-in/zoom-out to the reference point
    {
        if (is_ctrl_on)
        {
            if (is_shift_on)
            {
                const bool is_zoom_in = true;
                camera_zoom_in_reference_point(irc_ref, is_zoom_in);
            }
            else if(is_alt_on)
            {
                const bool is_zoom_in = false;
                camera_zoom_in_reference_point(irc_ref, is_zoom_in);
            }
        }
        break;
    }

    default:
        // no key assigned
        break;
    }
}

//----------------------------------------------------------------------
bool User_interaction::key_press_action(
    Nvindex_rendering_context* irc_ref,
    const mi::Sint32           keycode, 
    const bool                 is_shift_on,
    const bool                 is_ctrl_on,
    const bool                 is_alt_on,
    const bool                 is_meta_on)
{
    // ASCII key without any special key handling
    
    switch (keycode)
    {
    case 'B':
    {
        // WARN_LOG << "No viewing scenario support in this user interaction. 'B'.";
        return false; // This action is not supported here
        break;
    }
    case 'N':
    {
        // WARN_LOG << "No viewing scenario support in this user interaction. 'N'.";
        return false; // This action is not supported here
        break;
    }
    case 'M':
    {
        // WARN_LOG << "No viewing scenario support in this user interaction. 'M'.";
        return false; // This action is not supported here
        break;
    }

    default:
        // no key assigned
        break;
    }

    // Reached here, OK.
    return true;
}

//----------------------------------------------------------------------
void User_interaction::key_up_action(
    Nvindex_rendering_context* irc_ref,
    const mi::Sint32           keycode, 
    const bool                 is_shift_on,
    const bool                 is_ctrl_on,
    const bool                 is_alt_on,
    const bool                 is_meta_on)
{ 
    if (!is_shift_on)
    {
        get_examiner()->
            set_manipulation_mode(Examiner_manipulator::EM_TrackballRotation);
    }

    if (!is_ctrl_on)
    {
        get_examiner()->
            set_manipulation_mode(Examiner_manipulator::EM_TrackballRotation);
    }

    if (!is_alt_on)
    {
        // Release the alt key
    }

    if (!is_meta_on)
    {
        // Release the meta key
    }
}

//----------------------------------------------------------------------
void User_interaction::video_resize_action(
    Nvindex_rendering_context* irc_ref,
    const mi::Sint32           video_width, 
    const mi::Sint32           video_height)
{
    INFO_LOG << "video_resize_action: " << video_width << ", " << video_height;
}

//----------------------------------------------------------------------
void User_interaction::shutdown_server_action(
    Nvindex_rendering_context* irc_ref)
{
    Nvindex_AppData* appdata = Nvindex_AppData::instance();
    if (!nv::index_common::get_bool(appdata->peek_app_proj()->get("app::demo", "no")))
    {
        appdata->set_app_run(false); // exit the application
    }
}

//----------------------------------------------------------------------
void User_interaction::set_rotation_center_to_volume_roi_center(
    Nvindex_rendering_context* irc_ref)
{
    std::vector<mi::neuraylib::Tag> volume_tag_vec =
        Nvindex_AppData::instance()->get_volume_tag_vec();
    if (volume_tag_vec.empty())
    {
        WARN_LOG << "No volume found in the scene. Cannot set the volume center. (Ctrl+5 action)";
        return;
    }

    // get the volume bbox
    const mi::neuraylib::Tag volume_tag = volume_tag_vec[0];
    assert(volume_tag.is_valid());
    mi::math::Vector<mi::Float32, 3> volume_center(0.0f, 0.0f, 0.0f);
    if (get_center_of_volume_in_roi(irc_ref, volume_tag, volume_center))
    {
        INFO_LOG << "Set rotation Center to the volume center: " << volume_center;
        get_examiner()->set_examiner_rotation_center(volume_center);
    }
    else
    {
        WARN_LOG << "Cannot find the volume center. (Ctrl+5 action).";
    }
}

//----------------------------------------------------------------------
void User_interaction::camera_zoom_in_reference_point(
    Nvindex_rendering_context* irc_ref,
    bool                       is_zoom_in)
{
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        irc_ref->get_current_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    {
        assert(irc_ref->m_session_tag.is_valid());

        mi::base::Handle<const nv::index::ISession> session(
            dice_transaction->access<const nv::index::ISession>(irc_ref->m_session_tag));
        assert(session.is_valid_interface());

        mi::base::Handle<const nv::index::IScene> scene(
            dice_transaction->access<const nv::index::IScene>(session->get_scene()));
        assert(scene.is_valid_interface());
        {
            mi::base::Handle< nv::index::ICamera > main_cam(
                dice_transaction->edit< nv::index::ICamera>(get_multiple_stereo_camera()->get_main_base_camera_tag()));
            assert(main_cam.is_valid_interface());
            get_examiner()->camera_zoom_to_pan_reference_point(is_zoom_in, scene.get(), main_cam.get());
        }
        {
            mi::base::Handle< nv::index::ICamera > ortho_cam(
                dice_transaction->edit< nv::index::ICamera>(get_multiple_stereo_camera()->get_main_base_camera_tag()));
            assert(ortho_cam.is_valid_interface());
            get_examiner()->camera_zoom_to_pan_reference_point(is_zoom_in, scene.get(), ortho_cam.get());
        }
    }
    dice_transaction->commit();
}

//----------------------------------------------------------------------
