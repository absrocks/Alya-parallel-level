/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_RTMP_CONNECT_EVENT_HANDLER_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_RTMP_CONNECT_EVENT_HANDLER_H

#include <iterator>

#include <nv/index/iperformance_values.h>

#include "common/string_dict.h"
#include "nvindex_rendering_context.h"

/// An RTMP connect event handler that registers the stream and call event handlers above.
class Connect_event_handler : public mi::base::Interface_implement<mi::rtmp::IConnect_event_handler>
{
public:
    Connect_event_handler(
        Nvindex_rendering_context&           irc,
        const nv::index_common::String_dict& application_project,
        const std::string&                   session_cookie,
        const std::string&                   session_cookie_guest);

    bool handle(
        bool                    is_create,
        mi::rtmp::IConnection*  connection,
        const mi::IData*        command_arguments,
        const mi::IData*        user_arguments);

private:
    mi::base::Handle<mi::neuraylib::IScope>                     m_scope;
    Nvindex_rendering_context&                                  m_irc;
    /// copy of application_project (as at the construction time of this object)
    nv::index_common::String_dict                               m_application_project;
    std::string                                                 m_session_cookie;
    std::string                                                 m_session_cookie_guest;
};

#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_RTMP_CONNECT_EVENT_HANDLER_H
