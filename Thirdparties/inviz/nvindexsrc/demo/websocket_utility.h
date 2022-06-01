/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief web socket handler utility

#ifndef NVIDIA_INDEX_WEBSOCKET_UTILITY_H
#define NVIDIA_INDEX_WEBSOCKET_UTILITY_H

#include <nv/index/iindex.h>

#include <string>

class Nvindex_rendering_context;

//======================================================================
/// Video data buffer for WebSocket
/// A wrapper of IBuffer for IVideo_data
class Websocket_video_data_buffer : public mi::base::Interface_implement<mi::neuraylib::IBuffer>
{
public:
    /// Ctor
    Websocket_video_data_buffer(mi::neuraylib::IVideo_data* data);
    /// Dtor
    virtual ~Websocket_video_data_buffer();
    
    const mi::Uint8* get_data() const;
    const mi::Uint8* get_data();
    mi::Size get_data_size() const;

private:
    mi::base::Handle<mi::neuraylib::IVideo_data> m_data;
};

//======================================================================
/// Web socket handler.
///
/// Using the owner class's handler.
/// \tparam T class of the object which handles the websocket callback
template <class T> 
class Websocket_handler :
    public mi::base::Interface_implement<mi::http::IWeb_socket_handler>
{
public:
    /// Ctor
    ///
    /// \param[in] url_path WebSocket URL path (path only, e.g., /myapp)
    /// \param[in] handler_obj_ref callback owner 
    /// \param[in] callback callback method
    Websocket_handler(
        const char* url_path, 
        T*          handler_obj_ref,
        bool (T::*callback)(mi::http::IWeb_socket* websocket))
        :
        m_url_path(url_path),
        m_handler_obj_ref(handler_obj_ref),
        m_callback(callback)
    {
        // empty
    }

    /// Dtor
    virtual ~Websocket_handler()
    {
        // empty
    }

    /// handler wrapper
    bool handle(mi::http::IWeb_socket* websocket)
    {
        if (websocket->get_url_path() != m_url_path)
        {
            return false;
        }
        return (*m_handler_obj_ref.*m_callback)(websocket);
    }

private:
    /// URL string 
    std::string m_url_path;
    /// Handler object reference (This obj should not delete this)
    T*          m_handler_obj_ref;
    /// Callback method
    bool (T::*m_callback)(mi::http::IWeb_socket* websocket);
};


//======================================================================
/// Web socket data handler.
///
/// Using the handler_obj_ref class's handler.
/// \tparam T class of the object which handles the websocket callback
template <class T>
class Websocket_data_handler
    : public mi::base::Interface_implement<mi::http::IWeb_socket_data_handler>
{
public:
    /// Ctor
    ///
    /// \param[in] handler_obj_ref    callback handler_obj_ref 
    /// \param[in] callback callback method
    Websocket_data_handler(
        T* handler_obj_ref,
        void (T::*callback)(mi::http::IWeb_socket* websocket, mi::neuraylib::IBuffer* buffer,
            bool binary_frame))
        : 
        m_handler_obj_ref(handler_obj_ref),
        m_callback(callback)
    {
        // empty
    }

    /// Dtor
    virtual ~Websocket_data_handler()
    {
        // empty
    }

    /// handler wrapper
    void handle(mi::http::IWeb_socket* websocket, mi::neuraylib::IBuffer* buffer, bool binary_frame)
    {
        (*m_handler_obj_ref.*m_callback)(websocket, buffer, binary_frame);
    }

private:
    /// Handler object reference (This obj should not delete this)
    T* m_handler_obj_ref;
    /// Callback method
    void (T::*m_callback)(mi::http::IWeb_socket* websocket, mi::neuraylib::IBuffer* buffer,
        bool binary_frame);
};

//======================================================================
/// Web socket state handler.
///
/// Using the handler_obj_ref class's handler.
/// \tparam T class of the object which handles the websocket callback
template <class T> class Websocket_state_handler
    : public mi::base::Interface_implement<mi::http::IWeb_socket_state_handler>
{
public:
    /// Ctor
    ///
    /// \param[in] handler_obj_ref    callback handler_obj_ref 
    /// \param[in] callback callback method
    Websocket_state_handler(
        T* handler_obj_ref,
        void (T::*callback)(mi::http::IWeb_socket* websocket))
        :
        m_handler_obj_ref(handler_obj_ref),
        m_callback(callback)
    {
        // empty
    }

    /// Dtor
    virtual ~Websocket_state_handler()
    {
        // empty
    }

    /// handler wrapper
    void handle(mi::http::IWeb_socket* websocket)
    {
        (*m_handler_obj_ref.*m_callback)(websocket);
    }

private:
    /// Handler object reference (This obj should not delete this)
    T* m_handler_obj_ref;
    /// Callback method
    void (T::*m_callback)(mi::http::IWeb_socket* websocket);
};

//======================================================================
#endif  // NVIDIA_INDEX_WEBSOCKET_UTILITY_H
