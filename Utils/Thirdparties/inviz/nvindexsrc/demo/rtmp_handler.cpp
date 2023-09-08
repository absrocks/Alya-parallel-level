/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "rtmp_handler.h"

#include "common/common_utility.h"

#include "nvindex_appdata.h"

mi::Sint32 Play_event_handler::s_nb_connections = 0;
mi::base::Lock Play_event_handler::s_nb_connections_lock;

Play_event_handler::Play_event_handler()
  : Basic_RTMP_play_event_handler()
{
}

bool Play_event_handler::handle(
    bool                         is_start,
    mi::rtmp::IStream*           stream,
    mi::neuraylib::IVideo_data** out)
{
    bool result = Basic_RTMP_play_event_handler::handle(is_start, stream, out);

    if (result)
    {
        mi::base::Lock::Block block(&s_nb_connections_lock);

        if (is_start)
            s_nb_connections++;
        else
            s_nb_connections--;
    }

    return result;
}

mi::Sint32 Play_event_handler::get_nb_connections()
{
    mi::base::Lock::Block block(&s_nb_connections_lock);
    return s_nb_connections;
}

std::string Play_event_handler::get_video_codec_name()
{
    const nv::index_common::String_dict* dict = Nvindex_AppData::instance()->peek_app_proj();
    return dict->get("dice::rtmp_video_streaming::video_codec");
}

std::string Play_event_handler::get_fallback_video_codec_name()
{
    const nv::index_common::String_dict* dict = Nvindex_AppData::instance()->peek_app_proj();
    return dict->get("dice::rtmp_video_streaming::video_codec_fallback");
}

void Play_event_handler::init_video(
    mi::neuraylib::IVideo_encoder* codec,
    mi::rtmp::IStream*             stream)
{
    const nv::index_common::String_dict* dict = Nvindex_AppData::instance()->peek_app_proj();

    // Set video bitrate (if supported by the codec)
    const std::string video_bitrate = dict->get("dice::rtmp_video_streaming::video_bitrate");
    assert(video_bitrate != "");
    codec->set_parameter("bitrate", video_bitrate.c_str());

    // Set the frame rate for video encoding
    const std::string video_framerate = dict->get("dice::rtmp_video_streaming::video_framerate");
    assert(video_framerate != "");
    codec->set_parameter("framerate", video_framerate.c_str());
    
    // Set the preset for video encoding
    const std::string video_preset = dict->get("dice::rtmp_video_streaming::video_preset");
    if(video_preset != "")
        codec->set_parameter("preset", video_preset.c_str());

    if (get_video_codec_name() == "h264")
    {
        // Set the bitrate error for video encoding (x264 only)
        const std::string video_bitrate_error = dict->get("dice::rtmp_video_streaming::video_bitrate_error");
        assert(video_bitrate_error != "");
        codec->set_parameter("bitrate_error", video_bitrate_error.c_str());
    }

    // The "render_rate" indicates how often the RTMP server will try to call the registered
    // render event handler. This will be initialized with a default value when the video plugin
    // is started so we need to overwrite it again here.
    stream->set_property("render_rate", video_framerate.c_str());

    // Output the applied codec settings
    if (codec->get_parameter("bitrate"))
    {
        std::istringstream str(codec->get_parameter("bitrate"));
        mi::Sint32 bitrate = 0;
        str >> bitrate;
        INFO_LOG << "Video codec bitrate set to " << bitrate << " bits/s "
                 << "(" << (bitrate /  (8 * 1024)) << " kB/s)";
    }

    if (codec->get_parameter("bitrate_error") && !std::string(codec->get_parameter("bitrate_error")).empty())
        INFO_LOG << "Video codec bitrate error set to " << codec->get_parameter("bitrate_error");

    if (codec->get_parameter("framerate"))
        INFO_LOG << "Video codec framerate set to " << codec->get_parameter("framerate") << " fps";
        
    if (codec->get_parameter("preset"))
        INFO_LOG << "Video codec preset set to " << codec->get_parameter("preset");  
}

Frame_event_handler::Frame_event_handler(
    nv::index_common::Canvas_buffers* canvas_buffers)
  : Basic_RTMP_frame_event_handler(canvas_buffers)
{
    s_is_resized = false;
}

bool Frame_event_handler::handle(
    mi::rtmp::IStream*           stream,
    mi::neuraylib::IVideo_data** out,
    bool                         send_queue_is_full)
{
#ifdef DEBUG
    // Gather some statistics
    mi::Float64 encoding_time = nv::index_common::get_time();
#endif // DEBUG

    bool result = Basic_RTMP_frame_event_handler::handle(stream, out, send_queue_is_full);

#ifdef DEBUG
    encoding_time = nv::index_common::get_time() - encoding_time;

    mi::math::Vector<mi::Uint32, 2> size = get_canvas_size();
    DEBUG_LOG << "stat: Time for encoding a frame (" << size.x << "x" << size.y
              << "): " << encoding_time << " (ms)";
#endif // DEBUG

    return result;
}

/// definition of resizing status
bool Frame_event_handler::s_is_resized = false;

//----------------------------------------------------------------------
Stream_event_handler::Stream_event_handler(
    mi::base::Handle<mi::neuraylib::IScope> scope,
    nv::index_common::Canvas_buffers*       canvas_buffers)
  : nv::index_common::Basic_RTMP_stream_event_handler(scope, canvas_buffers)
{
}

mi::rtmp::IPlay_event_handler* Stream_event_handler::create_play_event_handler()
{
    return new Play_event_handler();
}

mi::rtmp::IFrame_event_handler* Stream_event_handler::create_frame_event_handler()
{
    return new Frame_event_handler(get_canvas_buffers());
}


