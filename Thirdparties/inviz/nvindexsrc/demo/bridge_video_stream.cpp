/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief NVIDIA IndeX DiCE Bridge Video Stream
#include "bridge_video_stream.h"

#include "nvindex_appdata.h"

//--------------------------------------------------------------------------------------------------
Bridge_video_stream::Data::Data(const void *ptr, mi::Size sz)
{
    if (sz != 0) {
        m_data = new mi::Uint8[sz];
        memcpy( m_data, ptr, sz );
        m_size = sz;
    }
}

//--------------------------------------------------------------------------------------------------
Bridge_video_stream::Data::~Data()
{
    if (m_data != 0)
        delete[] m_data;
}

//--------------------------------------------------------------------------------------------------
Bridge_video_stream::Bridge_video_stream(
    mi::bridge::IServer_session*        session,
    mi::Sint32                          video_context_id,
    Nvindex_rendering_context*          irc)
{
    m_cur_frame_id = 0;
    
    m_video_context = session->get_video_context(video_context_id);
    m_irc = irc;
    m_video_context->set_video_source(this);

    set_stream_format(std::string(Default_stream_format));
    set_stream_bitrate(Default_stream_bitrate);
    set_stream_fps(Default_stream_FPS);
    set_stream_gamma(Default_stream_gamma);
}

Bridge_video_stream::~Bridge_video_stream()
{
    if (m_video_context.is_valid_interface())
        m_video_context->close();
}


//--------------------------------------------------------------------------------------------------
void Bridge_video_stream::frame_ready(mi::Uint32 frame_id)
{
    m_cur_frame_id = frame_id;
    m_video_context->frame_ready();
}

//--------------------------------------------------------------------------------------------------
mi::Sint32 Bridge_video_stream::video_get_next_frame(
    mi::neuraylib::ICanvas** out_frame, 
    mi::neuraylib::IBuffer** out_data)
{
    // get canvas
    if(m_irc->m_canvas_buffers[0]->prepare_encoder_canvas())
    {
        *out_frame = m_irc->m_canvas_buffers[0]->get_encoder_canvas();
    }
    else
        *out_frame = 0;

    {
        mi::base::Lock::Block block(&Nvindex_AppData::instance()->m_performance_values_lock);
        if(Nvindex_AppData::instance()->m_performance_values)
        {
            nv::index::IPerformance_values* perf = Nvindex_AppData::instance()->m_performance_values.get();
            const mi::Float32 fps               = perf->get_time("frames_per_second");
            const mi::Float32 total_rendering   = perf->get_time("time_total_rendering");
            const mi::Float32 total_compositing = perf->get_time("time_total_compositing");
            const mi::Float32 rendering_only    = perf->get_time("time_rendering_only");
            const mi::Float32 gpu_upload        = perf->get_time("time_gpu_upload");
            
                // DICE suports only IString frame data, so encode it..
            std::ostringstream oss;
            oss << m_cur_frame_id 
                << fps
                << total_rendering
                << total_compositing
                << rendering_only
                << gpu_upload;
                // << std::endl; 

            // Copy into a data buffer and send.
            std::string s = oss.str();
            mi::base::Handle<mi::neuraylib::IBuffer> buffer( new Bridge_video_stream::Data( s.c_str(), s.length()+1 ) );
            buffer->retain();
            *out_data = buffer.get();
            
        }
    }
    
    return 0;
}

//--------------------------------------------------------------------------------------------------
void Bridge_video_stream::video_error(
    mi::Sint32 error_code, 
    const char* message)
{
    ERROR_LOG << "DiCE Bridge: Video_error " << error_code << ": " << message;
}
                
//--------------------------------------------------------------------------------------------------
void Bridge_video_stream::video_context_closed( mi::Sint32 reason )
{
    std::string errstr;
    errstr += "video_context closed: ";
    switch( reason ) 
    {
    case -1:
        errstr += "[-1] Network Error";
        ERROR_LOG << errstr;
        break;
    case 0:
        errstr += "[0] Closed by the server";
        WARN_LOG << errstr;
        break;
    case 1:
        errstr += "[1] Closed by the client";
        INFO_LOG << errstr;
        break;
    default:
        errstr += "Unknown";
        ERROR_LOG << errstr;
        break;
    }
    
}

//--------------------------------------------------------------------------------------------------
void Bridge_video_stream::set_stream_format( const std::string& format )
{
    m_stream_format = format;
    const std::string set_format = format == "auto" ? "h264" : format;
    m_video_context->set_video_format( set_format.c_str() );
    
    INFO_LOG << "DiCE Bridge: video stream format set to " << set_format;
}

//--------------------------------------------------------------------------------------------------
void Bridge_video_stream::set_stream_bitrate( mi::Uint32 bitrate )
{
    m_stream_bitrate = bitrate;
    m_video_context->set_bit_rate( bitrate );

    INFO_LOG << "DiCE Bridge: video stream bitrate set to " << bitrate;
}

//--------------------------------------------------------------------------------------------------
void Bridge_video_stream::set_stream_fps( mi::Uint32 fps )
{
    m_stream_fps = fps;
    m_video_context->set_max_frame_rate( fps );
    
    INFO_LOG << "DiCE Bridge: video stream fps set to " << fps;
}

void Bridge_video_stream::set_max_pending_frames( mi::Uint32 max_pending_frames)
{
    m_stream_max_pending_frames = max_pending_frames;
//    m_video_context->set_max_pending_frames( max_pending_frames );

    INFO_LOG << "DiCE Bridge: *Deprecated* video stream max pending frames set to " 
        << max_pending_frames;
}

//--------------------------------------------------------------------------------------------------
void Bridge_video_stream::set_stream_gamma( mi::Float32 gamma )
{
    m_stream_gamma = gamma;
    
    INFO_LOG << "DiCE Bridge: video stream gamma set to " << gamma;
    
}
