/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Basic RTMP server implementation

#ifndef NVIDIA_INDEX_BIN_COMMON_BASIC_RTMP_HANDLER_H
#define NVIDIA_INDEX_BIN_COMMON_BASIC_RTMP_HANDLER_H

#include "canvas.h"

#include <mi/base/handle.h>
#include <mi/base/interface_implement.h>

#include <mi/dice.h>
#include <mi/math/vector.h>

#include <mi/neuraylib/rtmp.h>

#include "forwarding_logger.h"

#include <string>

namespace nv {
namespace index_common {

/// An RTMP play event handler that chooses the video codec and initializes it.
class Basic_RTMP_play_event_handler : public mi::base::Interface_implement<mi::rtmp::IPlay_event_handler>
{
public:
    Basic_RTMP_play_event_handler()
    {
        // empty
    }

    virtual bool handle(
        bool                         is_start,
        mi::rtmp::IStream*           stream,
        mi::neuraylib::IVideo_data** out)
    {
        std::string peer = stream->get_connection()->get_peer_address();

        if (is_start)
        {
            const std::string video_codec = get_video_codec_name();
            const std::string fallback_video_codec = get_fallback_video_codec_name();

            INFO_LOG << "Video started for " << peer << " with codec '" << video_codec << "'";

            // Initialize the video codec. This may log default settings (bitrate, framerate).
            if (!stream->use_codec(video_codec.c_str()))
            {
                if (fallback_video_codec.empty())
                {
                    ERROR_LOG << "Failed to load video codec '" << video_codec << "'.";
                    return false;
                }
                else
                {
                    // Try the fallback
                    INFO_LOG << "Could not load video codec '" << video_codec << "', "
                             << "trying fallback '" << fallback_video_codec << "'";
                    if (!stream->use_codec(fallback_video_codec.c_str()))
                    {
                        ERROR_LOG << "Failed to load video codec '" << fallback_video_codec << "'.";
                        return false;
                    }
                }
            }

            mi::base::Handle<mi::neuraylib::IVideo_encoder> codec(stream->get_video_codec());
            init_video(codec.get(), stream);
        }
        else
        {
            INFO_LOG << "Video stopped for " << peer;

            mi::base::Handle<mi::neuraylib::IVideo_encoder> codec(stream->get_video_codec());
            if (codec)
            {
                codec->close(out);
            }
        }

        return true;
    }

protected:
    /// Returns the codec name to use
    virtual std::string get_video_codec_name()
    {
        return "screen video";
    }

    /// Returns the codec name to use when initialization of the primary codec fails. When this is
    /// an empty string, the fallback will be disabled.
    virtual std::string get_fallback_video_codec_name()
    {
        return "";
    }

    virtual void init_video(
        mi::neuraylib::IVideo_encoder* codec,
        mi::rtmp::IStream*  stream)
    {
        // empty
    }
};

/// An RTMP frame event handler that encodes a frame and delivers it to the RTMP server for sending.
/// Note that this event runs in another thread than the other event handlers, most importantly the
/// render handler, so care needs to be taken to avoid synchronization issues.
class Basic_RTMP_frame_event_handler : public mi::base::Interface_implement<mi::rtmp::IFrame_event_handler>
{
public:
    Basic_RTMP_frame_event_handler(Canvas_buffers* canvas_buffers)
        : m_need_to_call_init(true),
          m_canvas_buffers(canvas_buffers)
    {
        // empty
    }

    virtual bool handle(
        mi::rtmp::IStream*           stream,
        mi::neuraylib::IVideo_data** out,
        bool                         send_queue_is_full)
    {
        if (send_queue_is_full) // we do not want to increase buffering
        {
            return true;
        }

        if (is_resized())
        {
            m_need_to_call_init = true; // re-initialization needed after resized
            return true;
        }

        mi::base::Handle<mi::neuraylib::IVideo_encoder> codec(stream->get_video_codec());
        if (m_need_to_call_init)
        {
            mi::math::Vector<mi::Uint32, 2> canvas_size = get_canvas_size();

            if (codec->init(canvas_size.x, canvas_size.y, out) != 0)
            {
                ERROR_LOG << "Failed to init codec.";
                return false;
            }

            m_need_to_call_init = false;
        }

        // Prepare the canvas for video encoding: The return value tells us whether a new frame is
        // available, but we ignore it here and just encode the old frame again if there is no new one.
        m_canvas_buffers->prepare_encoder_canvas();

        // The encoder canvas always contains a complete image and is only accessed by this thread.
        mi::base::Handle<nv::index_common::Canvas> canvas(m_canvas_buffers->get_encoder_canvas());
        if (!canvas.is_valid_interface())
        {
            return true;
        }

        // Run the video encoding, this can be an expensive operation
        const bool result = codec->encode_canvas(canvas.get(), out) == 0;
        return result;
    }

    /// Returns whether the canvas was resized.
    virtual bool is_resized() const 
    {
        return false; // No resize handling in the default implementation
    }

protected:
    virtual mi::math::Vector<mi::Uint32, 2> get_canvas_size()
    {
        Canvas* canvas = m_canvas_buffers->get_encoder_canvas();
        if (canvas)
        {
            return mi::math::Vector<mi::Uint32, 2>(canvas->get_resolution_x(), canvas->get_resolution_y());
        }
        else
        {
            return mi::math::Vector<mi::Uint32, 2>(0);
        }
    }

private:
    bool                 m_need_to_call_init;
    Canvas_buffers*      m_canvas_buffers;
};

/// An RTMP stream event handler that registers the play and frame event handlers.
class Basic_RTMP_stream_event_handler : public mi::base::Interface_implement<mi::rtmp::IStream_event_handler>
{
public:
    Basic_RTMP_stream_event_handler(
        mi::base::Handle<mi::neuraylib::IScope> scope,
        Canvas_buffers*                         canvas_buffers)
      : m_canvas_buffers(canvas_buffers)
    {
        // empty
    }

    virtual bool handle(
        bool               is_create,
        mi::rtmp::IStream* stream,
        const mi::IData*   command_arguments)
    {
        std::string peer = stream->get_connection()->get_peer_address();
        
        if (is_create)
        {
            INFO_LOG << "RTMP stream created for " << peer;
            
            mi::base::Handle<mi::rtmp::IPlay_event_handler> play_event_handler(
                create_play_event_handler());
            stream->register_play_event_handler(play_event_handler.get());
            
            mi::base::Handle<mi::rtmp::IFrame_event_handler> frame_event_handler(
                create_frame_event_handler());
            stream->register_frame_event_handler(frame_event_handler.get());
        }
        else
        {
            INFO_LOG << "RTMP stream removed for " << peer;
        }
        
        return true;
    }


protected:
    virtual mi::rtmp::IPlay_event_handler* create_play_event_handler()
    {
        return new Basic_RTMP_play_event_handler();
    }

    virtual mi::rtmp::IFrame_event_handler* create_frame_event_handler()
    {
        return new Basic_RTMP_frame_event_handler(get_canvas_buffers());
    }

    virtual Canvas_buffers* get_canvas_buffers()
    {
        return m_canvas_buffers;
    }

private:
    Canvas_buffers*                    m_canvas_buffers;
};

}} // namespace nv::index_common

#endif // NVIDIA_INDEX_BIN_COMMON_BASIC_RTMP_HANDLER_H
