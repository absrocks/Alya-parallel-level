/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief RTMP server handler

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_RTMP_HANDLER_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_RTMP_HANDLER_H

#include "common/basic_rtmp_handler.h"

/// An RTMP play event handler that chooses the video codec and initializes it.
class Play_event_handler : public nv::index_common::Basic_RTMP_play_event_handler
{
public:
    Play_event_handler();

    virtual bool handle(
        bool                         is_start,
        mi::rtmp::IStream*           stream,
        mi::neuraylib::IVideo_data** out);

    // Returns the number of active video stream connections
    static mi::Sint32 get_nb_connections();

protected:
    virtual std::string get_video_codec_name();
    virtual std::string get_fallback_video_codec_name();

    virtual void init_video(
        mi::neuraylib::IVideo_encoder* codec,
        mi::rtmp::IStream*             stream);

    static mi::Sint32     s_nb_connections;
    static mi::base::Lock s_nb_connections_lock;
};

/// An RTMP frame event handler that encodes a frame and gives it to the RTMP server for
/// sending. Note that this event runs in another thread than the other event handlers, most
/// importantly the render handler, so care needs to be taken to avoid synchronization issues.
class Frame_event_handler : public nv::index_common::Basic_RTMP_frame_event_handler
{
public:
    Frame_event_handler(nv::index_common::Canvas_buffers* canvas_buffers);

    bool handle(
        mi::rtmp::IStream*           stream,
        mi::neuraylib::IVideo_data** out,
        bool                         send_queue_is_full);

public:
    /// Resize handling flag. 
    /// - This is class static for all connections.
    /// - is_resized() method is not static for the virtual method.
    /// We could make this implementation in the base class, but it is in common.
    static bool s_is_resized;

    /// set resize status 
    /// \param[in] is_resized true when resized
    static void set_resized(bool is_resized)
    {
        s_is_resized = is_resized;
    }

public:
    /// Resize handling.
    /// Note this is not static, but referring the class static for
    /// making this method virtual.  The base class needs this
    /// information for encoder stop/start.  
    /// \return true when resized.
    virtual bool is_resized() const 
    {
        return s_is_resized;
    }
};

/// An RTMP stream event handler that registers the play and frame event handlers above.
class Stream_event_handler : public nv::index_common::Basic_RTMP_stream_event_handler
{
public:
    Stream_event_handler(
        mi::base::Handle<mi::neuraylib::IScope> scope,
        nv::index_common::Canvas_buffers*       canvas_buffers);

protected:
    virtual mi::rtmp::IPlay_event_handler* create_play_event_handler();
    virtual mi::rtmp::IFrame_event_handler* create_frame_event_handler();
};

#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_RTMP_HANDLER_H
