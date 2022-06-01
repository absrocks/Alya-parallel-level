/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief html5 video stream implementation

#ifndef NVIDIA_INDEX_HTML5_VIDEO_STREAM_H
#define NVIDIA_INDEX_HTML5_VIDEO_STREAM_H

#include <mi/dice.h>

#include <nv/index/iindex.h>

#include <queue>
#include <string>
#include <map>

#include "common/string_dict.h"

#include "json.h"
#include "thread_utility.h"

namespace mi {
namespace neuraylib {
class ICanvas;
}
}

class Nvindex_rendering_context;
class Html5_video_connection;
class Html5_video_websocket_server;

//======================================================================
/// HTML5 client command handling static methods collection
///
/// The client commands come here and manipulates the IndeX
/// visualization server.
class Html5_client_command
{
public:
    /// default constructor
    Html5_client_command();
    /// destructor
    ~Html5_client_command();

public:
    /// mouse move handling
    ///
    /// \param[in] irc_ref    IndeX rendering context reference 
    /// \param[in] client_arg client argument. a json object.
    static void mouse_move_command(
        Nvindex_rendering_context* irc_ref,
        Json::Value&               client_arg);

    /// mouse down handling
    ///
    /// \param[in] irc_ref    IndeX rendering context reference 
    /// \param[in] client_arg client argument. a json object.
    static void mouse_down_command(
        Nvindex_rendering_context* irc_ref,
        Json::Value&               client_arg);

    /// mouse up handling
    ///
    /// \param[in] irc_ref    IndeX rendering context reference 
    /// \param[in] client_arg client argument. a json object.
    static void mouse_up_command(
        Nvindex_rendering_context* irc_ref,
        Json::Value&               client_arg);

    /// mouse wheel handling
    ///
    /// \param[in] irc_ref    IndeX rendering context reference 
    /// \param[in] client_arg client argument. a json object.
    static void mouse_wheel_command(
        Nvindex_rendering_context* irc_ref,
        Json::Value&               client_arg);

    /// key down handling (push special key (e.g., shift, alt only))
    ///
    /// \param[in] irc_ref    IndeX rendering context reference 
    /// \param[in] client_arg client argument. a json object.
    static void key_down_command(
        Nvindex_rendering_context* irc_ref,
        Json::Value&               client_arg);

    /// key press handling (push ascii key)
    ///
    /// \param[in] irc_ref    IndeX rendering context reference 
    /// \param[in] client_arg client argument. a json object.
    static void key_press_command(
        Nvindex_rendering_context* irc_ref,
        Json::Value&               client_arg);

    /// key up handling
    ///
    /// \param[in] irc_ref    IndeX rendering context reference 
    /// \param[in] client_arg client argument. a json object.
    static void key_up_command(
        Nvindex_rendering_context* irc_ref,
        Json::Value&               client_arg);

    /// Video resize handling
    ///
    /// \param[in] irc_ref    IndeX rendering context reference 
    /// \param[in] client_arg client argument. a json object.
    static void video_resize_command(
        Nvindex_rendering_context* irc_ref,
        Json::Value&               client_arg);

    /// Shutdown server command
    ///
    /// \param[in] irc_ref    IndeX rendering context reference 
    /// \param[in] client_arg client argument. a json object.
    static void shutdown_command(
        Nvindex_rendering_context* irc_ref,
        Json::Value&               client_arg);

private:
    /// copy constructor. prohibit until proved useful.
    Html5_client_command(Html5_client_command const &);
    /// operator=. prohibit until proved useful.
    Html5_client_command const & operator=(Html5_client_command const &);
};



//======================================================================
/// Html5 video encoder thread
class Html5_video_encoder_thread: public Application_thread_base
{
public:
    /// Ctor
    Html5_video_encoder_thread(Html5_video_connection* connection);
    /// Dtor
    virtual ~Html5_video_encoder_thread();

    /// Create the video encoder
    bool create_video_encoder();

    /// Initialize the video encoder
    ///
    /// \param[in] width  video size width
    /// \param[in] height video size height
    /// \return true when succeeded.
    bool initialize_video_encoder(mi::Sint32 width, mi::Sint32 height);

    /// set WebSocket state
    /// \param[in] websocket_state updated WebSocket state information
    void set_websocket_state(mi::http::IWeb_socket::State websocket_state);
    /// get WebSocket state
    /// \return current WebSocket state
    mi::http::IWeb_socket::State get_websocket_state();

    /// encode one frame and sent to the client
    /// \return true when encode can be continued
    bool encode_frame();

    /// shutdown 
    void shutdown();

    /// Stop the encoder loop and this thread
    void finish();
    
    /// thread runner, implementation of Application_thread_base
    virtual void run();

private:
    /// Reference to the video connection object
    Html5_video_connection*                         m_connection_ref;
    /// Video encoder (currently this is mp4 encoder)
    mi::base::Handle<mi::neuraylib::IVideo_encoder> m_video_encoder;
    /// WebSocket communication channel
    mi::base::Handle<mi::http::IWeb_socket>         m_websocket;
    /// Flag of the encoder initialized state
    bool                                            m_is_encoder_initialized;
    /// Current encoder loop state
    bool                                            m_encode_loop_running;
    /// current web socket state (reported by the state handler)
    mi::http::IWeb_socket::State                    m_websocket_state;
    /// current encoded frame number
    mi::Sint64                                      m_encoded_frame_number;
    /// Lock for canvas access
    mi::base::Lock                                  m_canvas_lock;
    /// Lock for web socket state change
    mi::base::Lock                                  m_websocket_state_lock;
};

//======================================================================
/// html5 video stream connection
class Html5_video_connection
{
public:
    /// Ctor
    /// 
    /// \param[in] irc_ref index rendering context reference (This is
    /// not the owner)
    /// \param[in] prj_ref project const object reference (This is
    /// not the owner)
    /// \param[in] html5_video_websocket_server_ref WebSocket server
    /// reference (This is not the owner)
    /// \param[in] connection_id my connection id (WebSocket server
    /// identifies this Html5_video_connection by this id).
    /// \param[in] socket A WebSocket connection. This object is the
    /// owner of this socket.
    Html5_video_connection(
        Nvindex_rendering_context*           irc_ref,
        nv::index_common::String_dict const* prj_ref,
        Html5_video_websocket_server*        html5_video_websocket_server_ref,
        mi::Sint32                           connection_id,        
        mi::http::IWeb_socket*               socket);

    /// Dtor
    ~Html5_video_connection();

    /// initializetion
    /// \return true when encoder creation succeeded.
    bool init();

    /// Data handler call back by websocket
    /// This callback function is executed when we receive data from the client browser  
    ///
    /// \param[in] websocket data handling websocket
    /// \param[in] buffer    data buffer
    /// \param[in] binary_frame flag indicating whether the data buffer is binary or text
    void data_handler(mi::http::IWeb_socket*  websocket, 
                      mi::neuraylib::IBuffer* buffer,
                      bool                    binary_frame);

    /// State handler call back by websocket
    /// This callback function is executed when an event occurs on the WebSocket occurs
    ///
    /// \param[in] websocket state changed websocket
    void state_handler(mi::http::IWeb_socket* websocket);

    /// get rendering context reference
    Nvindex_rendering_context* get_rendering_context_ref();

    /// get connection ID
    mi::Sint32 get_connection_id() const;

    /// get (retained) WebSocket
    mi::http::IWeb_socket* get_websocket();

    /// set WebSocket state to the encoder thread
    /// \param[in] websocket_state updated WebSocket state information
    void set_websocket_state(mi::http::IWeb_socket::State websocket_state);
    /// get WebSocket state from the encoder thread
    /// \return current WebSocket state
    mi::http::IWeb_socket::State get_websocket_state();

    /// get project
    /// \return project
    nv::index_common::String_dict const* get_project() const;

    /// shutdown the connection, stop the encode thread
    void shutdown_encodeing();    

private:
    /// prohibit default Ctor
    Html5_video_connection();
    /// reference to the nvindex_rendering_context (not owner)
    Nvindex_rendering_context*              m_irc_ref;
    /// project object reference (object const)
    nv::index_common::String_dict const*    m_project_ref;
    /// html5 video WebSocket server object reference.
    Html5_video_websocket_server*           m_html5_video_websocket_server_ref;
    /// connection ID. The WevSocket server recognize this object by this ID
    mi::Sint32                              m_connection_id;
    /// web socket for sending data 
    mi::base::Handle<mi::http::IWeb_socket> m_websocket;
    /// encoder thread object.
    Html5_video_encoder_thread*             m_encoder_thread;
    /// Lock for WebSocket access
    mi::base::Lock                          m_websocket_lock;
};

//======================================================================
/// Html5 video stream WebSocket server
/// 
/// When a client request to connect to web socket, this creates a
/// video connection
class Html5_video_websocket_server
{
public:
    /// ctor
    ///
    /// \param[in] irc_ref index rendering context reference
    /// \param[in] prj_ref project (option) reference
    Html5_video_websocket_server(
        Nvindex_rendering_context*           irc_ref,
        nv::index_common::String_dict const* prj_ref);

    /// dtor
    virtual ~Html5_video_websocket_server();

    /// connection request handle
    bool websocket_handler(mi::http::IWeb_socket* socket);

    /// close a connection by id
    void close_connection_by_id(mi::Sint32 connection_id);

    /// shutdown server, close all connections if they exist.
    void shutdown();

private:
    /// get the current ID and update the counter (need a lock)
    mi::Sint32 get_and_update_id();

private:
    /// reference to the nvindex_rendering_context (not owner)
    Nvindex_rendering_context*                    m_irc_ref;
    /// project object reference (object const)
    nv::index_common::String_dict const*          m_project_ref;
    /// connection ID
    mi::Sint32                                    m_connection_id;
    /// connection management
    std::map<mi::Sint32, Html5_video_connection*> m_connection_video_map;
    /// lock for connection id request
    mi::base::Lock                                m_connection_id_lock;
    /// lock for connection map access
    mi::base::Lock                                m_connection_map_lock;
};

//======================================================================
/// Html5 video stream connection server 
/// 
/// When a client request to connect to the server, this returns the 
/// specified html page.
class Html5_video_connection_http_server :
    public mi::base::Interface_implement<mi::http::IRequest_handler>
{
public:
    /// ctor
    Html5_video_connection_http_server();
    /// dtor
    virtual ~Html5_video_connection_http_server();
    /// connection request handle
    bool handle(mi::http::IConnection* connection);

private:
};

//======================================================================
/// Html5 video stream manager 
/// 
/// This manages html5 video stream related objects.
/// Currently the index rendering context owns this.
/// 
/// A simple class diagram
///
/// Html5_video_stream_manager
///   + has a --- Html5_video_connection_http_server
///   + has a --- Html5_video_websocket_server
///                 + has many --- Html5_video_connection
///
class Html5_video_stream_manager
{
public:
    /// Ctor
    Html5_video_stream_manager();
    /// Dtor
    virtual ~Html5_video_stream_manager();

    /// Initialize the video stream manager
    ///
    /// \param[in] irc_ref  index rendering context reference
    /// \param[in] prj_ref  project reference (options)
    /// \return true when succeeded.
    bool init(Nvindex_rendering_context*           irc_ref,    
              nv::index_common::String_dict const* prj_ref);

    /// Shutdown the html5 video stream
    void shutdown();

private:
    /// project object reference
    nv::index_common::String_dict const* m_project_ref;
    /// http server for html5 video stream
    mi::base::Handle<mi::http::IServer>  m_html5_http_server;
    /// html5 web socket server for multiple connections
    Html5_video_websocket_server*        m_html5_websocket_server;
};

//======================================================================
#endif  // NVIDIA_INDEX_HTML5_VIDEO_STREAM_H
