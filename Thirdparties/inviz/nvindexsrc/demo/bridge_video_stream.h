/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief NVIDIA IndeX DiCE Bridge Video Stream

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_BRIDGE_VIDEO_STREAM_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_BRIDGE_VIDEO_STREAM_H

static const int          Default_stream_FPS = 30;
static const int          Default_stream_bitrate = 5000000;
static const char* const  Default_stream_format = "h264";
static const float        Default_stream_gamma = 1.0f;

#include "nvindex_rendering_context.h"

class Bridge_video_stream
    : public mi::base::Interface_implement<mi::bridge::IVideo_source>
{
    public:

        Bridge_video_stream(
            mi::bridge::IServer_session*        session,
            mi::Sint32                          video_context_id,
            Nvindex_rendering_context*          irc);
            
        ~Bridge_video_stream();

        // Notify the video stream that a new frame should be read from the accum buffer
        void frame_ready(mi::Uint32 frame_id);

        // This is called by the video context when a frame_ready() is sent.
        virtual mi::Sint32 video_get_next_frame(
            mi::neuraylib::ICanvas** out_frame, 
            mi::neuraylib::IBuffer** out_data);

        virtual void video_error(
            mi::Sint32 error_code, 
            const char* message);
        
        virtual void video_context_closed( mi::Sint32 reason );

        void set_stream_format( const std::string& format );
        void set_stream_bitrate( mi::Uint32 bitrate );
        void set_stream_fps( mi::Uint32 fps );
        void set_max_pending_frames( mi::Uint32 max_pending_frames);
        void set_stream_gamma( mi::Float32 gamma );

    private:
        class Data : public
            mi::base::Interface_implement<mi::neuraylib::IBuffer>
        {
            public:
                Data(const void* ptr, mi::Size sz);
                ~Data();

                virtual const mi::Uint8* get_data() const
                    { return m_data; }
                virtual mi::Size get_data_size() const
                    { return m_size; }
            private:
                mi::Uint8   *m_data;
                mi::Size    m_size;
        };
        
        // // Stuff valid with creation
        mi::base::Handle<mi::bridge::IServer_video_context>     m_video_context;
        Nvindex_rendering_context*                              m_irc;

        // Stuff set through buffer attributes
        std::string     m_stream_format;
        mi::Uint32      m_stream_bitrate;
        mi::Uint32      m_stream_fps;
        mi::Uint32      m_stream_max_pending_frames;
        mi::Float32     m_stream_gamma;
        
        mi::Uint32      m_cur_frame_id;
};


#endif  //NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_BRIDGE_VIDEO_STREAM_H