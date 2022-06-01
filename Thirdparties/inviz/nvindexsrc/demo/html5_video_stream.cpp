/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief html5 video stream mechanism implementation

#include "html5_video_stream.h"

#include "common/canvas.h"
#include "common/common_utility.h"
#include "common/forwarding_logger.h"
#include "common/string_dict.h"
#include "common/tokenizer.h"

#include "nvindex_appdata.h"
#include "nvindex_rendering_context.h"
#include "websocket_utility.h"

#include <set>


//======================================================================
void Html5_client_command::mouse_move_command(
    Nvindex_rendering_context* irc_ref,
    Json::Value&               client_arg)
{
    assert(irc_ref != 0);

    const mi::Sint32 mouse_x = client_arg["x"].asInt();
    const mi::Sint32 mouse_y = client_arg["y"].asInt();
    const mi::Sint32 button  = client_arg["button"].asInt();

    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        irc_ref->get_current_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
    {
        mi::base::Handle<const nv::index::ISession> session(
            dice_transaction->access<const nv::index::ISession>(irc_ref->m_session_tag));
        assert(session.is_valid_interface());

        Nvindex_AppData::instance()->get_user_interaction(0)->
            mouse_motion(button, 
                         mouse_x,
                         mouse_y,
                         session->get_scene(),
                         dice_transaction.get());
    }
    dice_transaction->commit();
}

//----------------------------------------------------------------------
void Html5_client_command::mouse_down_command(
    Nvindex_rendering_context* irc_ref,
    Json::Value&               client_arg)
{
    assert(irc_ref != 0);

    const mi::Sint32 mouse_x = client_arg["x"].asInt();
    const mi::Sint32 mouse_y = client_arg["y"].asInt();
    const mi::Sint32 button  = client_arg["button"].asInt();

    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        irc_ref->get_current_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
    {
        Nvindex_AppData::instance()->get_user_interaction(0)->
            mouse_button_action(button, 
                                Examiner_manipulator::Button_down,
                                mouse_x, 
                                mouse_y,
                                dice_transaction.get());
    }
    dice_transaction->commit();
}


//----------------------------------------------------------------------
void Html5_client_command::mouse_up_command(
    Nvindex_rendering_context* irc_ref,
    Json::Value&               client_arg)
{
    assert(irc_ref != 0);

    const mi::Sint32 mouse_x = client_arg["x"].asInt();
    const mi::Sint32 mouse_y = client_arg["y"].asInt();
    const mi::Sint32 button  = client_arg["button"].asInt();

    mi::math::Vector<mi::Sint32, 2> pick_position(-1, -1);
    {
        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            irc_ref->get_current_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
        {
            pick_position = Nvindex_AppData::instance()->get_user_interaction(0)->
                mouse_button_action(button,
                                    Examiner_manipulator::Button_up,
                                    mouse_x, 
                                    mouse_y,
                                    dice_transaction.get());

            // INFO_LOG << "mouse up pick: " << pick_position;
        }
        dice_transaction->commit();
    }

    if (Nvindex_AppData::instance()->is_pick_operation_enabled() && 
        (button == Examiner_manipulator::Left_button)            && // Only left button can pick
        (pick_position.x >= 0.0f) && (pick_position.y >= 0.0f))
    {
        Nvindex_AppData::instance()->get_user_interaction(0)->
            issue_pick_command(irc_ref, 
                               Nvindex_AppData::instance()->is_enable_multi_view_mode(),
                               pick_position);
    }
}

//----------------------------------------------------------------------
void Html5_client_command::mouse_wheel_command(
    Nvindex_rendering_context* irc_ref,
    Json::Value&               client_arg)
{
    const mi::Sint32 mouse_x = client_arg["mouse_x"]. asInt();
    const mi::Sint32 mouse_y = client_arg["mouse_y"]. asInt();
    const mi::Sint32 delta   = client_arg["delta"].   asInt();

    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        irc_ref->get_current_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
    {
        mi::base::Handle<const nv::index::ISession> session(
            dice_transaction->access<const nv::index::ISession>(irc_ref->m_session_tag));
        assert(session.is_valid_interface());

        Nvindex_AppData::instance()->get_user_interaction(0)->
            mouse_wheel_action(irc_ref, 
                               session->get_scene(),
                               mouse_x, 
                               mouse_y,
                               delta,
                               dice_transaction.get());
    }
    dice_transaction->commit();
}

//----------------------------------------------------------------------
void Html5_client_command::key_down_command(
    Nvindex_rendering_context* irc_ref,
    Json::Value&               client_arg)
{
    const mi::Sint32 keycode     = client_arg["keyCode"]. asInt();
    const bool       is_shift_on = client_arg["shiftKey"].asBool();
    const bool       is_ctrl_on  = client_arg["ctrlKey"]. asBool();
    const bool       is_alt_on   = client_arg["altKey"].  asBool();
    const bool       is_meta_on  = client_arg["metaKey"]. asBool();

    Nvindex_AppData::instance()->get_user_interaction(0)->
        key_down_action(
            irc_ref,
            keycode, 
            is_shift_on,
            is_ctrl_on,
            is_alt_on,
            is_meta_on);

}

//----------------------------------------------------------------------
void Html5_client_command::key_press_command(
    Nvindex_rendering_context* irc_ref,
    Json::Value&               client_arg)
{
    const mi::Sint32 keycode     = client_arg["keyCode"]. asInt();
    const bool       is_shift_on = client_arg["shiftKey"].asBool();
    const bool       is_ctrl_on  = client_arg["ctrlKey"]. asBool();
    const bool       is_alt_on   = client_arg["altKey"].  asBool();
    const bool       is_meta_on  = client_arg["metaKey"]. asBool();

    Nvindex_AppData::instance()->get_user_interaction(0)->
        key_press_action(
            irc_ref,
            keycode, 
            is_shift_on,
            is_ctrl_on,
            is_alt_on,
            is_meta_on);

}

//----------------------------------------------------------------------
void Html5_client_command::key_up_command(
    Nvindex_rendering_context* irc_ref,
    Json::Value&               client_arg)
{
    const mi::Sint32 keycode     = client_arg["keyCode"]. asInt();
    const bool       is_shift_on = client_arg["shiftKey"].asBool();
    const bool       is_ctrl_on  = client_arg["ctrlKey"]. asBool();
    const bool       is_alt_on   = client_arg["altKey"].  asBool();
    const bool       is_meta_on  = client_arg["metaKey"]. asBool();

    Nvindex_AppData::instance()->get_user_interaction(0)->
        key_up_action(
            irc_ref,
            keycode, 
            is_shift_on,
            is_ctrl_on,
            is_alt_on,
            is_meta_on);
}

//----------------------------------------------------------------------
void Html5_client_command::video_resize_command(
    Nvindex_rendering_context* irc_ref,
    Json::Value&               client_arg)
{
    mi::Sint32 video_width  = client_arg["video_width"]. asInt();
    mi::Sint32 video_height = client_arg["video_height"]. asInt();
    if (video_width <= 0)
    {
        ERROR_LOG << "zero or negative video width.";
        video_width = 1;
    }
    if (video_height <= 0)
    {
        ERROR_LOG << "zero or negative video height.";
        video_height = 1;
    }

    Nvindex_AppData::instance()->get_user_interaction(0)->
        video_resize_action(
            irc_ref,
            video_width, 
            video_height);
}

//----------------------------------------------------------------------
void Html5_client_command::shutdown_command(
    Nvindex_rendering_context* irc_ref,
    Json::Value&               client_arg)
{
    Nvindex_AppData::instance()->get_user_interaction(0)->
        shutdown_server_action(irc_ref);
}

//----------------------------------------------------------------------
//======================================================================
//----------------------------------------------------------------------
Html5_video_encoder_thread::Html5_video_encoder_thread(Html5_video_connection* connection)
    :
    m_connection_ref(connection),
    m_is_encoder_initialized(false),
    m_encode_loop_running(false),
    m_websocket_state(mi::http::IWeb_socket::WS_STATE_INIT),
    m_encoded_frame_number(-1)
{
    // empty
}

//----------------------------------------------------------------------
Html5_video_encoder_thread::~Html5_video_encoder_thread()
{
    // INFO_LOG << "Html5_video_encoder_thread dtor is called";
    shutdown();
}

//----------------------------------------------------------------------
bool Html5_video_encoder_thread::create_video_encoder()
{
    assert(!m_websocket.is_valid_interface());

    m_websocket = mi::base::Handle<mi::http::IWeb_socket>(m_connection_ref->get_websocket());
    mi::base::Handle<mi::neuraylib::IVideo_codec_factory> video_codec_factory(
        m_connection_ref->get_rendering_context_ref()->m_iindex_if->
        get_api_component<mi::neuraylib::IVideo_codec_factory>());
    mi::base::Handle<mi::neuraylib::IVideo_encoder> video_encoder(
        video_codec_factory->create_video_encoder("mp4_encoder"));
    if (!video_encoder.is_valid_interface())
    {
        ERROR_LOG << "Cannot find an mp4 video encoder.";
        return false;
    }
    m_video_encoder.swap(video_encoder);
    INFO_LOG << "An mp4 encoder is created.";

    assert(m_connection_ref != 0);
    assert(m_connection_ref->get_project());

    // mp4_encoder specific parameters
    nv::index_common::String_dict mp4_params;
    const mi::Sint32 ret = string_dict_key_prefix_filter(
        *(m_connection_ref->get_project()), "dice::html5_video_streaming::mp4_encoder", mp4_params, true);
    nv::index_common::no_unused_variable_warning_please(ret);    
    assert(ret == 0);           // currently no mp4_encoder options

    for (nv::index_common::String_dict::const_iterator si = mp4_params.begin();
         si != mp4_params.end();
         ++si)
    {
        const std::string key = si->first;
        const std::string val = si->second;

        if (m_video_encoder->set_parameter(key.c_str(), val.c_str()))
        {
            INFO_LOG << "Set mp4_encoder option: {" << key << ", " << val << "}";
        }
        else 
        {
            ERROR_LOG << "Fail to set mp4_encoder option: {" << key << ", " << val << "}";
            return false;
        }
    }

    // Setup the internal encoder
    if (!(m_connection_ref->get_project()->is_defined("dice::html5_video_streaming::mp4_video_encoder")))
    {
        ERROR_LOG << "No dice::html5_video_streaming::mp4_video_encoder option. Failed to create an encoder.";
        return false;
    }
    const std::string encoder_list_str = m_connection_ref->get_project()->get("dice::html5_video_streaming::mp4_video_encoder");
    const std::string separators = ",";
    std::vector<std::string> encoder_str_vec;
    nv::index_common::Tokenizer::parse(encoder_list_str, separators, encoder_str_vec);
    if (encoder_str_vec.empty())
    {
        ERROR_LOG << "No value of dice::html5_video_streaming::mp4_video_encoder. Failed to create an encoder.";
        return false;
    }

    // check the current supported encoders
    std::set<std::string> supported_encoder_set;
    supported_encoder_set.insert("h264");
    supported_encoder_set.insert("h264_nvenc");

    const mi::Size nb_encoder = encoder_str_vec.size();
    for (mi::Size i = 0; i < nb_encoder; ++i)
    {
        const std::string internal_encoder = encoder_str_vec[i];
        if (supported_encoder_set.find(internal_encoder) == supported_encoder_set.end())
        {
            ERROR_LOG << "Unsupported mp4_video_encoder: " << internal_encoder << ", skipped.";
            continue;
        }

        INFO_LOG << "Creating mp4_video_encoder: " << internal_encoder;
        if (m_video_encoder->set_parameter("mp4_video_encoder", internal_encoder.c_str()))
        {
            break;
        }
        else
        {
            ERROR_LOG << "Failed to create mp4_video_encoder: " << internal_encoder;
            // Try again
        }
    }

    const std::string created_internal_encoder = m_video_encoder->get_parameter("mp4_video_encoder");
    if (!created_internal_encoder.empty())
    { 
        INFO_LOG << "Succeeded to create mp4_video_encoder: " << created_internal_encoder;
        return true;
    }
        
    ERROR_LOG << "No encoder was created.";
    return false;
}

//----------------------------------------------------------------------
bool Html5_video_encoder_thread::initialize_video_encoder(mi::Sint32 width, mi::Sint32 height)
{
    assert(width  > 0);
    assert(height > 0);
    assert(m_video_encoder.is_valid_interface());

    // Set the encode parameter depends on the current internal encoder
    std::string internal_encoder_name = "";
    if (m_video_encoder->get_parameter("mp4_video_encoder") != 0)
    {
        internal_encoder_name = m_video_encoder->get_parameter("mp4_video_encoder");
    }
    if (internal_encoder_name.empty())
    {
        ERROR_LOG << "No internal encoder found.";
        return false;
    }
    // DEBUG_LOG << "mp4 encoder name: " << internal_encoder_name;

    // set internal encoder (h264, h264_nvenc) parameters here (mp4_encoder has an internal encoder)
    nv::index_common::String_dict const* src_params = m_connection_ref->get_project();
    assert(src_params != 0);

    nv::index_common::String_dict  internal_encoder_params;
    const bool is_delete_prefix = true;
    const std::string param_prefix_key = "dice::html5_video_streaming::" + internal_encoder_name + "::";
    const mi::Sint32 nb_params = nv::index_common::string_dict_key_prefix_filter(
        *src_params, param_prefix_key, internal_encoder_params, is_delete_prefix);

    if (nb_params == 0)
    {
        INFO_LOG << "no " << internal_encoder_name << " parameters found.";
    }
    else
    {
        INFO_LOG << "Setting " << internal_encoder_name << " parameters.";
        for (nv::index_common::String_dict::const_iterator si = internal_encoder_params.begin();
             si != internal_encoder_params.end(); ++si)
        {
            if (m_video_encoder->set_parameter(si->first.c_str(), si->second.c_str()))
            {
                INFO_LOG << "Set " << internal_encoder_name << " parameter: "
                         << si->first << ": " << si->second;
            }
            else
            {
                ERROR_LOG << "Fail to set: " << si->first << ": " << si->second;
            }
        }
    }
    
    // Now initialize all 
    mi::neuraylib::IVideo_data* video_data = 0;
    const bool result = (m_video_encoder->init(width, height, &video_data) == 0);
    if (!result)
    {
        ERROR_LOG << "Cannot initialize the video encoder.";
        // Now, we need to shutdown this encoder and the encoding thread.
        return false;
    }

    // If the video encoder has any initial data, send it out
    if(video_data)
    {
        // video_data contains some metadata to the decoder (or may
        // empty, depends on the video_encoder)
        Websocket_video_data_buffer video_data_buffer(video_data);
        // Note: This shoule be after web socket connection has established.
        // send the data to client
        m_websocket->write( &video_data_buffer, true);
    }

    return true;
}

//----------------------------------------------------------------------
void Html5_video_encoder_thread::set_websocket_state(mi::http::IWeb_socket::State websocket_state)
{
    {
        mi::base::Lock::Block block(&m_websocket_state_lock);
        // DEBUG_LOG << "websocket state change: " << m_websocket_state << " to "
        //           << websocket_state;
        m_websocket_state = websocket_state;
    }
}

//----------------------------------------------------------------------
mi::http::IWeb_socket::State Html5_video_encoder_thread::get_websocket_state()
{
    mi::http::IWeb_socket::State websocket_state = mi::http::IWeb_socket::WS_STATE_INIT;
    {
        mi::base::Lock::Block block(&m_websocket_state_lock);
        websocket_state = m_websocket_state;
    }
    return websocket_state;
}

//----------------------------------------------------------------------
bool Html5_video_encoder_thread::encode_frame()
{
    if (get_websocket_state() != mi::http::IWeb_socket::WS_STATE_CONNECTED)
    {
        INFO_LOG << "Connection has not yet handshaked, no encoding.";
        return true;            // still next frame can be rendered when connected
    }

    // The encoder should be initialized *after* the WebSocket has
    // been connected, since the initial frame should be sent.
    if (!m_is_encoder_initialized)
    {
        const mi::math::Vector<mi::Sint32,2> canvas_resolution = 
            nv::index_common::get_vec_sint32_2(m_connection_ref->get_project()->
                                               get("index::canvas_resolution", "800 800"));
        if (initialize_video_encoder(canvas_resolution.x, canvas_resolution.y))
        {
            m_is_encoder_initialized = true;
        }
        else
        {
            // fail to initialize the encoder, shutdown and stop encoding the frame.
            m_is_encoder_initialized = false; // just make sure that initialization was failed.
            return false;                     // encode should be terminated.
        }
    }

    mi::neuraylib::IVideo_data* video_data = 0;
    bool encode_result = false;
    {
        mi::base::Lock::Block block(&m_canvas_lock);
        if (!m_connection_ref->get_rendering_context_ref()->m_canvas_buffers[0]->prepare_encoder_canvas())
        {
            // no canvas available, wait 0.005 seconds
            nv::index_common::sleep(0.005f);
            block.release();
            return true;        // still next frame can be rendered rendered when the canvas is ready.
        }

        const mi::Sint64 rendered_frame_number =
            m_connection_ref->get_rendering_context_ref()->m_canvas_buffers[0]->get_render_frame_number();

        if (m_encoded_frame_number < rendered_frame_number)
        {
            mi::base::Handle<nv::index_common::Canvas> encoder_canvas(
                m_connection_ref->get_rendering_context_ref()->m_canvas_buffers[0]->get_encoder_canvas());
            assert(encoder_canvas.is_valid_interface());
            // encoder: [in] encoding_canvas, [out] video_data
            // If frame_number are the same, no encode, thus encode_result == false, so no send.
            encode_result = (m_video_encoder->encode_canvas(encoder_canvas.get(), &video_data)
                             == 0);
            m_encoded_frame_number = rendered_frame_number;
        }
    }

    if(encode_result && (video_data != 0))
    {
        Websocket_video_data_buffer video_data_buffer(video_data);
        m_websocket->write(&video_data_buffer, true);
    }

    return true;
}

//----------------------------------------------------------------------
void Html5_video_encoder_thread::shutdown()
{
    m_connection_ref = 0;       // not owner, just reference: no delete

    m_encode_loop_running = false;

    if (m_video_encoder.is_valid_interface())
    {
        mi::neuraylib::IVideo_data* video_data = 0;
        mi::Sint32 ret_code = m_video_encoder->close(&video_data);
        if (ret_code != 0)
        {
            WARN_LOG << "Fail to close the video encoder. error_code: " << ret_code;
        }
        m_video_encoder = 0;
    }

    // This WebSocket may have already shutdown since the server
    // should be shutdown first. However, this = 0 is fine because our
    // ref count mechanism.
    m_websocket = 0;
}

//----------------------------------------------------------------------
void Html5_video_encoder_thread::finish()
{
    m_encode_loop_running = false;
    join();
    INFO_LOG << "Terminated the thread of Html5_video_encoder_thread.";
}

//----------------------------------------------------------------------
void Html5_video_encoder_thread::run()
{
    m_encode_loop_running = true;
    while (m_encode_loop_running)
    {
        if (!this->encode_frame())
        {
            break;              // exit the encoding loop
        }
    }
    INFO_LOG << "Terminated Html5_video_encoder_thread::run().";    
    // m_connection_ref->shutdown_encodeing();
}

//----------------------------------------------------------------------
//======================================================================
//----------------------------------------------------------------------
Html5_video_connection::Html5_video_connection(
    Nvindex_rendering_context*           irc_ref,
    nv::index_common::String_dict const* prj_ref,
    Html5_video_websocket_server*        html5_video_websocket_server_ref,
    mi::Sint32                           connection_id,
    mi::http::IWeb_socket*               websocket)
    :
    m_irc_ref(irc_ref),
    m_project_ref(prj_ref),
    m_html5_video_websocket_server_ref(html5_video_websocket_server_ref),
    m_connection_id(connection_id),
    m_websocket(mi::base::make_handle_dup<mi::http::IWeb_socket>(websocket)),
    m_encoder_thread(0)
{
    assert(m_irc_ref != 0);
    assert(m_html5_video_websocket_server_ref != 0);
    assert(m_websocket.is_valid_interface());

    // INFO_LOG << "Html5_video_connection::Installing a data handler and a state handler";
    mi::base::Handle<Websocket_data_handler<Html5_video_connection> > data_handler(
        new Websocket_data_handler<Html5_video_connection>(
            this, &Html5_video_connection::data_handler));
    m_websocket->set_data_handler(data_handler.get());

    mi::base::Handle<Websocket_state_handler<Html5_video_connection> > state_handler(
        new Websocket_state_handler<Html5_video_connection>(
            this, &Html5_video_connection::state_handler));
    m_websocket->set_state_handler(state_handler.get());
}

//----------------------------------------------------------------------
Html5_video_connection::~Html5_video_connection()
{
    shutdown_encodeing();
}

//----------------------------------------------------------------------
bool Html5_video_connection::init()
{
    assert(m_websocket.is_valid_interface());

    m_encoder_thread = new Html5_video_encoder_thread(this);
    if (!m_encoder_thread->create_video_encoder())
    {
        ERROR_LOG << "Failed to create and initialize an mp4 encoder thread.";
        return false;
    }

    return true;
}

//----------------------------------------------------------------------
void Html5_video_connection::data_handler(
    mi::http::IWeb_socket*  socket,
    mi::neuraylib::IBuffer* buffer, 
    bool /* binary_frame */)
{
    std::stringstream sstr(std::string(reinterpret_cast<const char*>(
                                           buffer->get_data()), buffer->get_data_size()));
    // INFO_LOG << "data_handler: " << sstr.str();

    Json::Value  client_command;
    Json::Reader json_reader;
    json_reader.parse(sstr.str(), client_command);

    if (!client_command.isMember("command"))
    {
        ERROR_LOG << "data_handler: no command from the client.";
        return;
    }

    const std::string com  = client_command["command"].asString();
    Json::Value client_arg = client_command[com];

    // This is a command dispatcher.
    if (com == "mouse_move")
    {
        Html5_client_command::mouse_move_command(m_irc_ref, client_arg);
    }
    else if (com == "mouse_down")
    {
        Html5_client_command::mouse_down_command(m_irc_ref, client_arg);
    }
    else if (com == "mouse_up")
    {
        Html5_client_command::mouse_up_command(m_irc_ref, client_arg);
    }
    else if (com == "mouse_wheel")
    {
        Html5_client_command::mouse_wheel_command(m_irc_ref, client_arg);
    }
    else if (com == "key_down")
    {
        Html5_client_command::key_down_command(m_irc_ref, client_arg);
    }
    else if (com == "key_press")
    {
        Html5_client_command::key_press_command(m_irc_ref, client_arg);
    }
    else if (com == "key_up")
    {
        Html5_client_command::key_up_command(m_irc_ref, client_arg);
    }
    else if (com == "disconnect")
    {
        INFO_LOG << "A html5 video stream client has been disconnected.";
        // Html5_client_command::disconnect_command(m_irc_ref, client_arg);
    }
    else if (com == "video_resize")
    {
        INFO_LOG << "received video resize command.";
        Html5_client_command::video_resize_command(m_irc_ref, client_arg);
    }
    else if (com == "shutdown_server")
    {
        INFO_LOG << "received shutdown_server command.";
        Html5_client_command::shutdown_command(m_irc_ref, client_arg);
    }
    else
    {
        ERROR_LOG << "No such command: " << sstr.str();
    }
}

//----------------------------------------------------------------------
void Html5_video_connection::state_handler(mi::http::IWeb_socket* websocket)
{
    m_encoder_thread->set_websocket_state(websocket->get_state());

    switch(websocket->get_state())
    {
    case mi::http::IWeb_socket::WS_STATE_CLOSED:
    case mi::http::IWeb_socket::WS_STATE_ERROR:
    {
        INFO_LOG << "Connection to " << websocket->get_peer_address()
                 << " has been lost. Close the connection.";
        m_websocket->close();
        m_html5_video_websocket_server_ref->close_connection_by_id(m_connection_id);
    }
    break;

    case mi::http::IWeb_socket::WS_STATE_CONNECTED:
    {
        assert(m_encoder_thread != 0);
        // start encoding
        m_encoder_thread->start();
    }
    break;

    default:
        break;
    }
}

//----------------------------------------------------------------------
Nvindex_rendering_context* Html5_video_connection::get_rendering_context_ref()
{
    return m_irc_ref;
}

//----------------------------------------------------------------------
mi::Sint32 Html5_video_connection::get_connection_id() const
{
    return m_connection_id;
}

//----------------------------------------------------------------------
mi::http::IWeb_socket* Html5_video_connection::get_websocket()
{
    {
        mi::base::Lock::Block block(&m_websocket_lock);
        m_websocket->retain();
    }
    return m_websocket.get();
}

//----------------------------------------------------------------------
void Html5_video_connection::set_websocket_state(mi::http::IWeb_socket::State websocket_state)
{
    // forward to the encoder thread
    m_encoder_thread->set_websocket_state(websocket_state);
}

//----------------------------------------------------------------------
mi::http::IWeb_socket::State Html5_video_connection::get_websocket_state()
{
    // forward to the encoder thread
    return m_encoder_thread->get_websocket_state();
}

//----------------------------------------------------------------------
nv::index_common::String_dict const* Html5_video_connection::get_project() const
{
    return m_project_ref;
}

//----------------------------------------------------------------------
void Html5_video_connection::shutdown_encodeing()
{
    // Finish encoder thread
    if (m_encoder_thread != 0)
    {
        m_encoder_thread->finish();
    }

    m_irc_ref                          = 0;
    m_project_ref                      = 0;
    m_html5_video_websocket_server_ref = 0;
    // m_connection_id

    m_websocket = 0;

    if (m_encoder_thread != 0)
    {
        m_encoder_thread->shutdown();
        delete m_encoder_thread;
        m_encoder_thread = 0;
    }
}

//----------------------------------------------------------------------
//======================================================================
//----------------------------------------------------------------------
Html5_video_websocket_server::Html5_video_websocket_server(
    Nvindex_rendering_context*           irc_ref,
    nv::index_common::String_dict const* prj_ref
    )
    :
    m_irc_ref(irc_ref),
    m_project_ref(prj_ref),
    m_connection_id(0)
{
    // empty
}

//----------------------------------------------------------------------
Html5_video_websocket_server::~Html5_video_websocket_server()
{
    m_irc_ref     = 0; // de-reference only since this is not the owner
    m_project_ref = 0; // de-reference only since this is not the owner
    assert(m_connection_video_map.empty());
}

//----------------------------------------------------------------------
bool Html5_video_websocket_server::websocket_handler(mi::http::IWeb_socket* websocket)
{
    INFO_LOG << "A client has requested a WebSocket connection.";

    // Reject web socket connections not matching our fixed (hardcoded) URL.
    const char* url_path = websocket->get_url_path();
    if (url_path == 0)
    {
        ERROR_LOG << "Websocket_handler didn't get a url path.";
        return false;
    }
    const std::string index_websocket_url_path = "/index_app";
    if (std::string(url_path) != index_websocket_url_path)
    {
        ERROR_LOG << "Websocket_handler url path is not " << index_websocket_url_path;
        return false;
    }
    // INFO_LOG << "websocket_handler url: " << url;

    const mi::Sint32 connection_id = get_and_update_id();
    {
        mi::base::Lock::Block block(&m_connection_map_lock);
        assert(m_connection_video_map.find(connection_id) == m_connection_video_map.end());
        // The socket is now owned by the Html5_video_connection
        Html5_video_connection* vcon =
            new Html5_video_connection(m_irc_ref, m_project_ref, this, connection_id, websocket);
        m_connection_video_map[connection_id] = vcon;
        vcon->init();

        // set the current socket state
        vcon->set_websocket_state(websocket->get_state());
    }

    return true;
}

//----------------------------------------------------------------------
void Html5_video_websocket_server::close_connection_by_id(mi::Sint32 connection_id)
{
    // Lock the map and delete one element.
    {
        mi::base::Lock::Block block(&m_connection_map_lock);
        if (m_connection_video_map.find(connection_id) == m_connection_video_map.end())
        {
            ERROR_LOG << "No such connection ("<< connection_id << ") exists.";
            return;
        }

        Html5_video_connection* vcon = m_connection_video_map[connection_id];
        m_connection_video_map[connection_id] = 0;
        m_connection_video_map.erase(connection_id);
        delete vcon;
        vcon = 0;
    }
}

//----------------------------------------------------------------------
void Html5_video_websocket_server::shutdown()
{
    // Lock the m_connection_video_map and get all remaining connection ids.
    std::vector<mi::Sint32> con_id_vec;
    {
        mi::base::Lock::Block block(&m_connection_map_lock);
        for (std::map<mi::Sint32, Html5_video_connection*>::iterator vi = m_connection_video_map.begin();
             vi != m_connection_video_map.end();
             ++vi)
        {
            const mi::Sint32 con_id = vi->first;
            con_id_vec.push_back(con_id);
        }
    }

    // Delete each connection id. Since close_connection_by_id() can
    // be called other threads and that locks the connection map. This
    // locks the map recursively. Thus, we do this two step way.
    for (std::vector<mi::Sint32>::const_iterator ci = con_id_vec.begin(); ci != con_id_vec.end(); ++ci)
    {
        const mi::Sint32 con_id = *ci;
        DEBUG_LOG << "Closing WebSocket connection: " << con_id;
        close_connection_by_id(con_id);
    }
    assert(m_connection_video_map.empty());
}

//----------------------------------------------------------------------
mi::Sint32 Html5_video_websocket_server::get_and_update_id()
{
    mi::Sint32 cur_id = -1;
    {
        mi::base::Lock::Block block(&m_connection_id_lock);
        cur_id = m_connection_id;
        ++m_connection_id;
    }
    assert(cur_id >= 0);

    return cur_id;
}

//----------------------------------------------------------------------
//======================================================================
//----------------------------------------------------------------------
Html5_video_connection_http_server::Html5_video_connection_http_server()
{
    // empty
}

//----------------------------------------------------------------------
Html5_video_connection_http_server::~Html5_video_connection_http_server()
{
    // empty
}

//----------------------------------------------------------------------
bool Html5_video_connection_http_server::handle(mi::http::IConnection* connection)
{
    INFO_LOG << "Html5_video_connection_http_server gets a connection request.";

    // If we need to inline the connection page, we can write
    // connection->print("HTML file contents");

    // Currently, re-direction a html file is implemented due to its
    // simplicity.

    return true;
}

//----------------------------------------------------------------------
//======================================================================
//----------------------------------------------------------------------
Html5_video_stream_manager::Html5_video_stream_manager()
    :
    m_project_ref(0),
    // m_html5_http_server,
    m_html5_websocket_server(0)
{
    // empty
}

//----------------------------------------------------------------------
Html5_video_stream_manager::~Html5_video_stream_manager()
{
    m_project_ref = 0;           // de-reference (This is not the owner)

    // First WebSocket (order matters)
    if (m_html5_websocket_server != 0)
    {
        WARN_LOG << "html5 WebSocket server is still running. Please shutdown first.";
        m_html5_websocket_server->shutdown();
        delete m_html5_websocket_server;
        m_html5_websocket_server = 0;
    }

    // Second http server (order matters)
    if (m_html5_http_server.is_valid_interface())
    {
        WARN_LOG << "html5 http server is still running. Please shutdown first.";
        m_html5_http_server->shutdown();
        m_html5_http_server = 0;
    }
}

//----------------------------------------------------------------------
bool Html5_video_stream_manager::init(
    Nvindex_rendering_context*     irc_ref,
    nv::index_common::String_dict const* prj_ref)
{
    assert(irc_ref != 0);
    assert(prj_ref != 0);
    m_project_ref = prj_ref;

    // Create a web socket request handler instance
    mi::base::Handle<mi::http::IFactory> http_factory(
        irc_ref->m_iindex_if->get_api_component<mi::http::IFactory>());
    assert(http_factory.is_valid_interface());

    mi::base::Handle<mi::http::IServer> http_server(http_factory->create_server());
    m_html5_http_server.swap(http_server);
    assert(m_html5_http_server.is_valid_interface());

    if (!(m_project_ref->is_defined("dice::html5_video_streaming::port")))
    {
        ERROR_LOG << "dice::html5_video_streaming::port is not defined. "
                  << "Failed to initialize the html5 video stream.";
        return false;
    }
    const std::string port_str = m_project_ref->get("dice::html5_video_streaming::port");
    const std::string redirect_fname =
        m_project_ref->get("dice::html5_video_streaming::redirect_html_file", "/html5_video_stream.html");
    INFO_LOG << "html5 http server redirect html file path: " << redirect_fname;

    // Redirect HTTP requests to the file html5_video_stream.html which
    // contains the JavaScript for open our video stream
    mi::base::Handle<mi::http::IRequest_handler> redirect_handler_1(
        http_factory->create_redirect_handler("/", redirect_fname.c_str()));
    m_html5_http_server->install(redirect_handler_1.get());

    mi::base::Handle<mi::http::IRequest_handler> redirect_handler_2(
        http_factory->create_redirect_handler("/index.html", redirect_fname.c_str()));
    m_html5_http_server->install(redirect_handler_2.get());

    mi::base::Handle<mi::http::IRequest_handler> redirect_handler_3(
        http_factory->create_file_handler("/", ".", true));
    m_html5_http_server->install(redirect_handler_3.get());

    mi::base::Handle<mi::http::IRequest_handler> connection_handler(
        new Html5_video_connection_http_server());
    m_html5_http_server->install(connection_handler.get()); // m_html5_http_server retains the request_handler.

    // Install the websocket handler
    assert(m_html5_websocket_server == 0);
    std::string url_path = "/index_app";
    m_html5_websocket_server = new Html5_video_websocket_server(irc_ref, prj_ref);
    mi::base::Handle<Websocket_handler<Html5_video_websocket_server> > websocket_handler;
    websocket_handler = new Websocket_handler<Html5_video_websocket_server>(
        url_path.c_str(), m_html5_websocket_server, &Html5_video_websocket_server::websocket_handler);
    m_html5_http_server->install(websocket_handler.get());

    // Assemble server address
    std::string address = "0.0.0.0:";
    address += port_str;

    const mi::Uint32 ret = m_html5_http_server->start(address.c_str());
    if (ret != 0)
    {
        ERROR_LOG << "WebSocket server listen: " << address << ", failed to start.";
    }
    else
    {
        INFO_LOG << "WebSocket server listen: " << address << ", start succeeded.";
    }

    return ret == 0;
}

//----------------------------------------------------------------------
void Html5_video_stream_manager::shutdown()
{
    // First, we need to close the encoder and its websocket Since
    // when a websocket shutdowns, it wants to release the server side
    // connection data. In this moment, the server should still alive.
    // Therefore, we first shutdown the websocket side, then shutdown
    // the server.
    if (m_html5_websocket_server != 0)
    {
        m_html5_websocket_server->shutdown();
        delete m_html5_websocket_server;
        m_html5_websocket_server = 0;
    }

    // Second, shutdown the server.
    if (m_html5_http_server.is_valid_interface())
    {
        m_html5_http_server->shutdown();
        m_html5_http_server = 0;
    }
}

//----------------------------------------------------------------------
//======================================================================
//----------------------------------------------------------------------
