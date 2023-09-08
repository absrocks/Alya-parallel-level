/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief NVIDIA IndeX rendering context

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_NVINDEX_RENDERING_CONTEXT_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_NVINDEX_RENDERING_CONTEXT_H

#include <mi/dice.h>

#include <nv/index/isession.h>

#include <vector>

#include "common/canvas.h"
#include "common/forwarding_logger.h"

#include "frame_info_callbacks.h"
#include "html5_video_stream.h"
#include "http_handler.h"
#include "span_renderer_if.h"


/// NVIDIA IndeX rendering context
class Nvindex_rendering_context
{
public:
    /// default constructor
    Nvindex_rendering_context()
        :
        m_current_scope(0),
        m_session_tag(mi::neuraylib::NULL_TAG),
        m_html5_video_stream_manager(0),
        m_frame_num(0),
        m_client_connection_issued(false),
        m_initialized(false),
        m_scene_setup_done(false),
        m_is_tsteps_playforward(false),
        m_is_tsteps_playbackward(false),
        m_is_pause_tsteps(false),
        m_is_tsteps_stepback(false),
        m_is_tsteps_stepforward(false)
    {
        // empty
    }
    /// destructor
    ~Nvindex_rendering_context()
    {
        // empty, though all the smart pointers will be destructed.
    }

    /// is valid this context?
    /// \return true when all the members are valid.
    bool is_valid() const
    {
        bool ret = true;
        ret = ret && m_database.is_valid_interface();
        ret = ret && m_global_scope.is_valid_interface();
        ret = ret && m_iindex_if.is_valid_interface();
        ret = ret && m_iindex_session.is_valid_interface();
        ret = ret && m_iindex_rendering.is_valid_interface();
        ret = ret && m_icluster_configuration.is_valid_interface();
        ret = ret && (m_session_tag.is_valid());
        ret = ret && m_progress_callback.is_valid_interface();
        // m_http_factory       may not exist (when rtmp video stream is off)
        // m_http_server        may not exist (when rtmp video stream is off)
        // http_request_handler may not exist (when rtmp video stream is off)
        // m_rtmp_factory       may not exist (when rtmp video stream is off)
        // m_rtmp_server        may not exist (when rtmp video stream is off)
        // m_connect_handler    may not exist (when rtmp video stream is off)
        ret = ret && m_iindex_scene_query.is_valid_interface();

        return ret;
    }

    /// call before DiCE shutdown
    void shutdown()
    {
        m_database                  = 0;
        m_global_scope              = 0;
        //   m_iindex_if   = 0;  // should not set to 0 here

        m_iindex_session            = 0;
        m_iindex_rendering          = 0;
        m_icluster_configuration    = 0;
        m_session_tag               = mi::neuraylib::NULL_TAG;
        m_progress_callback         = 0;
        m_frame_info_callbacks      = 0;

        m_http_factory              = 0;
        m_http_server               = 0;
        m_http_request_handler      = 0;
        m_rtmp_factory              = 0;
        m_rtmp_server               = 0;
        m_connect_handler           = 0;

        if (m_html5_video_stream_manager != 0)
        {
            m_html5_video_stream_manager->shutdown();
            delete m_html5_video_stream_manager;
            m_html5_video_stream_manager = 0;
        }

        // clear canvas buffers (handle)
        for (mi::Size i = 0; i < m_canvas_buffers.size(); ++i)
        {
            m_canvas_buffers[i] = 0;
        }
        m_canvas_buffers.clear();

        m_iindex_scene_query      = 0;

        m_cluster_change_callback = 0;
        m_span_buffer             = 0;
        m_frame_num               = 0; 
    }

    void set_current_scope(int nb)
    {
        if (nb >= 0 && nb < static_cast<int>(m_extra_scopes.size() + 1))
        {
            m_current_scope = nb;
        }
    }

    mi::base::Handle<mi::neuraylib::IScope> get_current_scope() const
    {
        const int current = m_current_scope;
        if (current == 0)
        {
            return m_global_scope;
        }
        else
        {
            return m_extra_scopes[current - 1];
        }
    }

public:
    /// database
    mi::base::Handle<mi::neuraylib::IDatabase> m_database;
    /// global scope
    mi::base::Handle<mi::neuraylib::IScope> m_global_scope;
    /// additional scopes (optional)
    std::vector<mi::base::Handle<mi::neuraylib::IScope> > m_extra_scopes;
    /// identifies the currently selected scope (0 means global scope)
    int m_current_scope;

    /// additional scopes for viewing scenarios (deprecated)
    std::vector<mi::base::Handle<mi::neuraylib::IScope> > m_viewing_scenario_scopes;
    /// nvindex library interface
    mi::base::Handle<nv::index::IIndex> m_iindex_if;
    /// iindex session
    mi::base::Handle<nv::index::IIndex_session> m_iindex_session;
    /// iindex rendering
    mi::base::Handle<nv::index::IIndex_rendering> m_iindex_rendering;
    /// rendering property query
    mi::base::Handle<nv::index::ICluster_configuration> m_icluster_configuration;
    /// session tag
    mi::neuraylib::Tag m_session_tag;
    /// progress callback
    mi::base::Handle<nv::index::IProgress_callback> m_progress_callback;
    /// frame info callback
    mi::base::Handle<Frame_info_callbacks> m_frame_info_callbacks;
    /// http factory
    mi::base::Handle<mi::http::IFactory> m_http_factory;
    /// http server
    mi::base::Handle<mi::http::IServer> m_http_server;
    /// request handler
    mi::base::Handle<HTTP_request_handler> m_http_request_handler;
    /// rtmp factory
    mi::base::Handle<mi::rtmp::IFactory> m_rtmp_factory;
    /// rtmp server
    mi::base::Handle<mi::rtmp::IServer> m_rtmp_server;
    /// connect handler
    mi::base::Handle<mi::rtmp::IConnect_event_handler> m_connect_handler;
    /// html5 video stream manager
    Html5_video_stream_manager*         m_html5_video_stream_manager;
    /// bridge server
    mi::base::Handle<mi::bridge::IBridge_server> m_bridge_server;
    /// bridge application
    mi::base::Handle<mi::bridge::IApplication> m_bridge_application;
    /// bridge application session handler
    mi::base::Handle<mi::bridge::IApplication_session_handler> m_bridge_application_session_handler;
    /// bridge video source
    mi::base::Handle<mi::bridge::IVideo_source> m_bridge_video_source;
    
    /// canvas buffers for rendering and video encoding
    std::vector<mi::base::Handle<nv::index_common::Canvas_buffers> >  m_canvas_buffers;
    /// scene query
    mi::base::Handle<nv::index::IIndex_scene_query> m_iindex_scene_query;
    /// cluster change callback 
    mi::base::Handle<nv::index::ICluster_change_callback> m_cluster_change_callback;
    /// span buffer renderer interface
    mi::base::Handle<Span_renderer_IF> m_span_buffer;
    /// frame number (set at the top of Index_application::render_frame())
    mi::Uint32 m_frame_num;

    /// Viewport list for multi-view mode
    mi::base::Handle<nv::index::IViewport_list> m_viewport_list;
    
    bool m_client_connection_issued;
    bool m_initialized;
    bool m_scene_setup_done;
    bool m_is_tsteps_playforward;
    bool m_is_tsteps_playbackward;
    bool m_is_pause_tsteps;
    bool m_is_tsteps_stepback;
    bool m_is_tsteps_stepforward;

private:
    /// copy constructor. prohibit until proved useful.
    Nvindex_rendering_context(Nvindex_rendering_context const &);
    /// operator=. prohibit until proved useful.
    Nvindex_rendering_context & operator=(Nvindex_rendering_context const &);
};

//----------------------------------------------------------------------
#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_NVINDEX_RENDERING_CONTEXT_H
