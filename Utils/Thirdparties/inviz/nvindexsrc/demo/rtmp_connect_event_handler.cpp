/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/

#include "rtmp_connect_event_handler.h"

#include <nv/index/iscene.h>

#include "heightfield_workflow_functionality.h"
#include "multiple_camera.h"
#include "nvindex_appdata.h"
#include "rtmp_call_event_handler.h"
#include "rtmp_handler.h"

using namespace nv::index_common;

Connect_event_handler::Connect_event_handler(
    Nvindex_rendering_context&           irc,
    const nv::index_common::String_dict& application_project,
    const std::string&                   session_cookie,
    const std::string&                   session_cookie_guest)
  : m_irc(irc),
    m_application_project(application_project),
    m_session_cookie(session_cookie),
    m_session_cookie_guest(session_cookie_guest)
{
    const char* p_necessary_key_list[] = {
        "dice::rtmp_video_streaming::video_codec",
        "dice::rtmp_video_streaming::video_bitrate",
        "dice::rtmp_video_streaming::video_framerate",
        "dice::rtmp_video_streaming::video_bitrate_error",
        0,
    };
    std::vector< std::string > undef_list;
    bool const is_keys_ok = is_all_keys_defined(m_application_project, p_necessary_key_list, &undef_list);
    std::stringstream sstr;
    if(!is_keys_ok)
    {
        std::copy(undef_list.begin(), undef_list.end(), std::ostream_iterator< std::string >(sstr, " "));
        ERROR_LOG << "Undefined project entries: " << sstr.str();
    }
    assert(is_keys_ok);
}

bool Connect_event_handler::handle(
    bool                    is_create,
    mi::rtmp::IConnection*  connection,
    const mi::IData*        command_arguments,
    const mi::IData*        user_arguments)
{
    if (is_create)
    {
        // Authentication using a session cookie
        std::string received_session_cookie;

        std::string flash_player = "unknown";

        mi::Uint32 viewing_scenario_id = 0;
        if (user_arguments)
        {
            // Extract the cookie sent by the client
            mi::base::Handle<const mi::IMap> map(user_arguments->get_interface<const mi::IMap>());
            if (map && map->has_key("session_cookie"))
            {
                mi::base::Handle<const mi::IString> value(
                    map->get_value("session_cookie")->get_interface<const mi::IString>());
                if (value)
                {
                    received_session_cookie = value->get_c_str();
                }
            }

            // Get Flash Player version
            if (map && map->has_key("flash_player"))
            {
                mi::base::Handle<const mi::IString> value(
                    map->get_value("flash_player")->get_interface<const mi::IString>());
                if (value)
                {
                    flash_player = value->get_c_str();
                }
            }

            // Get Flash Player version
            if (map && map->has_key("viewing_session"))
            {
                mi::base::Handle<const mi::IString> value(
                    map->get_value("viewing_session")->get_interface<const mi::IString>());
                if (value)
                {
                    viewing_scenario_id = get_uint32(value->get_c_str());
                }
            }
        }

        INFO_LOG << "RTMP connection from " << connection->get_peer_address()
                 << ", Flash Player: " << flash_player;

        // Only allow connection if the cookies match (or both are empty)
        if (received_session_cookie != m_session_cookie && received_session_cookie != m_session_cookie_guest)
        {
            ERROR_LOG << "rtmp::Connect_event_handler [" << connection->get_peer_address() << "]: "
                      << "Client denied because of invalid session cookie!";

            return false;
        }

        // We are authenticated now ...

        // When not initialized yet then start the scene setup handler, but no video stream.
        if (!m_irc.m_initialized && received_session_cookie == m_session_cookie)
        {
            // Admin user
            mi::base::Handle<mi::rtmp::ICall_event_handler> call_event_handler(
                new Scene_setup_call_event_handler(m_irc));
            connection->register_remote_call_handler(call_event_handler.get(), "client_event");
            return true;
        }
        else if (!m_irc.m_initialized && received_session_cookie == m_session_cookie_guest)
        {
            // Guest user
            mi::base::Handle<mi::rtmp::ICall_event_handler> call_event_handler(
                new Scene_setup_call_event_handler(m_irc, false));
            connection->register_remote_call_handler(call_event_handler.get(), "client_event");
            return true;
        }

        // Handle remote viewing scenario
        String_dict* p_app_proj = Nvindex_AppData::instance()->peek_app_proj();
        const bool started_as_viewing_scenario = get_bool(p_app_proj->get("app::remote_viewing_capability", "no"));
        if(started_as_viewing_scenario && !m_irc.m_client_connection_issued)
        {
            mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
                m_irc.get_current_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
            assert(dice_transaction.is_valid_interface());
            // Setup session information
            const std::string session_name = "nv_reference_viewer_session_0";
            m_irc.m_session_tag = dice_transaction->name_to_tag(session_name.c_str());
            if(!(m_irc.m_session_tag.is_valid()))
            {
                INFO_LOG << "No previously stored session available. Simply continue with the remote service.";
            }
            else
            {
                assert(m_irc.m_session_tag.is_valid());
                INFO_LOG << "Using previously stored session. "
                         << "The session tag " << m_irc.m_session_tag.id
                         << " has identifier name: '" << session_name << "'";
                mi::base::Handle<const nv::index::ISession> session(
                    dice_transaction->access<const nv::index::ISession>(m_irc.m_session_tag));
                assert(session.is_valid_interface());
                mi::base::Handle<const nv::index::IScene> scene(
                    dice_transaction->access<const nv::index::IScene>(session->get_scene()));
                Multiple_stereo_camera* ms_cam = 
                    Nvindex_AppData::instance()->get_user_interaction(0)->get_multiple_stereo_camera();
                ms_cam->set_main_camera(Stereo_camera(scene->get_camera(), scene->get_camera()));

                const std::vector<mi::neuraylib::Tag> volume_tag_vec
                    = Nvindex_AppData::instance()->get_volume_tag_vec();
                if(!(volume_tag_vec.empty()))
                {
                    mi::base::Handle<const nv::index::IRegular_volume> volume_scene_element(
                        dice_transaction->access<const nv::index::IRegular_volume>(volume_tag_vec.at(0)));
                    assert(volume_scene_element.is_valid_interface());

                    const mi::neuraylib::Tag& current_volume_colormap_tag
                        = volume_scene_element->assigned_colormap();
                    set_colormap_tag(0, current_volume_colormap_tag);
                }

                // Create heightfield workflow functionality
                INFO_LOG << "Creating heightfield workflow functionality";
                Heightfield_workflow_functionality::init(m_irc.m_iindex_session, m_irc.m_session_tag);
            }
            dice_transaction->commit();
        }

        // Choose the appropriate canvas for the viewing scenario
        Canvas_buffers* canvas = NULL;
        if (m_irc.m_initialized)
        {
            if(viewing_scenario_id<m_irc.m_canvas_buffers.size() )
            {
                canvas = m_irc.m_canvas_buffers[viewing_scenario_id].get();
            }
            else
            {
                canvas = m_irc.m_canvas_buffers[0].get();
            }
        }

        // Install rtmp stream event handler and remote call event handler
        mi::base::Handle<mi::rtmp::IStream_event_handler> stream_event_handler(
            new Stream_event_handler(m_irc.m_global_scope, canvas));
        connection->register_stream_event_handler(stream_event_handler.get());
        mi::base::Handle<Call_event_handler> call_event_handler(
            new Call_event_handler(viewing_scenario_id, m_irc));
        connection->register_remote_call_handler(call_event_handler.get(), "client_event");

        if(viewing_scenario_id==1)
        {
            INFO_LOG << "Enabling rendering of viewing scenario " << viewing_scenario_id
                     << " using a separate video stream.";
            WARN_LOG << "Check Connect_event_handler::handle()";

            mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
                m_irc.m_viewing_scenario_scopes[0]->create_transaction<mi::neuraylib::IDice_transaction>());
            assert(dice_transaction.is_valid_interface());

            call_event_handler->render_viewing_scenario(
                viewing_scenario_id,
                dice_transaction);

            dice_transaction->commit();
        }

        m_irc.m_client_connection_issued = true;
    }
    return true;
}
