/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief HTTP server handler

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_HTTP_HANDLER_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_HTTP_HANDLER_H

#include "common/basic_http_handler.h"

/// Our HTTP request handler
class HTTP_request_handler : public nv::index_common::Basic_HTTP_request_handler
{
public:
    HTTP_request_handler(
        const std::string&                              main_swf_file,
        std::vector<mi::math::Vector<mi::Uint32, 2> >&  resolutions,
        const std::string&                              rtmp_port = "");

    virtual bool handle(mi::http::IConnection* connection);

    void set_authentication(
        const std::string& auth_user,
        const std::string& auth_password);

    void set_guest_authentication(
        const std::string& auth_guest_user,
        const std::string& auth_guest_password);

    void set_resolutions(
        std::vector<mi::math::Vector<mi::Uint32, 2> >&  resolutions);

    std::string get_session_cookie() const;
    std::string get_session_cookie_guest() const;

private:
    virtual std::string get_page_title() const;

    virtual std::string get_favicon_url() const { return "/favicon.ico"; }

    virtual std::string get_extra_head_items() const {
        return "<link rel=\"shortcut icon\" href=\"" + get_favicon_url() + "\">\n";
    }

    virtual void get_flash_vars(
        mi::http::IRequest*                 request,
        std::map<std::string, std::string>& flash_vars);

    std::string m_auth_user;
    std::string m_auth_password;
    std::string m_session_cookie;

    std::string m_auth_guest_user;
    std::string m_auth_guest_password;
    std::string m_session_cookie_guest;

    mutable mi::base::Lock m_lock;
};

#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_HTTP_HANDLER_H
