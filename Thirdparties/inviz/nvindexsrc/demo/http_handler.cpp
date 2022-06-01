/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "http_handler.h"

#include <cstdio>
#include <cstdlib>

#include "common/common_utility.h"

#include "favicon.h"
#include "nvindex_appdata.h"

namespace {

static const char base64_chars[] = {
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz"
        "0123456789+/"};

std::string base64_encode(const unsigned char* bytes_to_encode, mi::Sint32 in_len) {
    std::string ret;
    mi::Sint32 i = 0;
    mi::Sint32 j = 0;
    unsigned char char_array_3[3];
    unsigned char char_array_4[4];
    while (in_len--)
    {
        char_array_3[i++] = *(bytes_to_encode++);

        if (i == 3)
        {
            char_array_4[0] =
                (char_array_3[0] & 0xfc) >> 2;
            char_array_4[1] =
                ((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
            char_array_4[2] =
                ((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);
            char_array_4[3] =
                char_array_3[2] & 0x3f;
            for (i = 0; (i < 4) ; i++)
                ret += base64_chars[char_array_4[i]];
            i = 0;
        }
    }

    if (i)
    {
        for (j = i; j < 3; j++)
            char_array_3[j] = '\0';
        char_array_4[0] =  (char_array_3[0] & 0xfc) >> 2;
        char_array_4[1] = ((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
        char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);
        char_array_4[3] =   char_array_3[2] & 0x3f;
        for (j = 0; (j < i + 1); j++)
            ret += base64_chars[char_array_4[j]];
        while ((i++ < 3))
            ret += '=';
    }
    return ret;
}

std::string generate_cookie()
{
    unsigned char buf[24];
#ifdef _WIN32
    srand(static_cast<int>(nv::index_common::get_time())); //FIXME: Doing this every time is not a good idea

    for (mi::Sint32 i = 0; i < sizeof(buf); ++i) {
        buf[i] = static_cast<mi::Uint8>((static_cast<mi::Float64>(rand()) / RAND_MAX) * 255.0);
    }
#else
    // Generate our random session cookie
    FILE* p_file = fopen("/dev/urandom", "r"); // that should be a good source of randomness
    if (p_file != 0)
    {
        fread(buf, sizeof(buf), 1, p_file);
        fclose(p_file);
    }
#endif // LINUX

    std::string cookie;
    const std::string s = base64_encode(buf, sizeof(buf));

    // Convert to modified Base64 for URL, so that it contains no characters that would need
    // encoding when used in URLs.
    for (size_t i=0; i < s.length(); ++i)
    {
        char c = s[i];
        if (c == '+')
            c = '-';
        else if (c == '/')
            c = '_';

        // Skip padding '=' altogether
        if (c != '=')
            cookie += c;
    }

    return cookie;
}

} // namespace

HTTP_request_handler::HTTP_request_handler(
    const std::string&                              main_swf_file,
    std::vector<mi::math::Vector<mi::Uint32, 2> >&  resolutions,
    const std::string&                              rtmp_port
) : Basic_HTTP_request_handler(main_swf_file, resolutions, rtmp_port)
{
    m_session_cookie       = generate_cookie();
    m_session_cookie_guest = generate_cookie();
}

std::string HTTP_request_handler::get_page_title() const
{
    const std::string default_title = "NVIDIA IndeX Viewer";
    const nv::index_common::String_dict* dict = Nvindex_AppData::instance()->peek_app_proj();
    return dict->get("app::page_title", default_title);
 }

void HTTP_request_handler::get_flash_vars(
    mi::http::IRequest*                 request,
    std::map<std::string, std::string>& flash_vars)
{
    Basic_HTTP_request_handler::get_flash_vars(request, flash_vars);

    // Session cookie
    const char* session_cookie = request->get_argument("session_cookie");
    if (session_cookie != 0 && !std::string(session_cookie).empty())
    {
        flash_vars["    session_cookie"] = "\"" + std::string(session_cookie) + "\"";
    }

    // Pass some settings from the project file to the flash application
    const nv::index_common::String_dict* dict = Nvindex_AppData::instance()->peek_app_proj();
    mi::Sint32 auto_reconnect = nv::index_common::get_sint32(dict->get("app::auto_reconnect", "0"));
    if (auto_reconnect > 0)
    {
        std::ostringstream s;
        s << auto_reconnect;
        flash_vars["    auto_reconnect"] = s.str();
    }
}

namespace {

bool check_credentials(
    const std::string& user,
    const std::string& password,
    const std::string& received_credentials)
{
    if (user.empty())
        return false;

    // Construct base64-encoded string from the correct credentials
    std::string good_credentials = user + ":" + password;
    good_credentials = base64_encode(
        reinterpret_cast<const unsigned char*>(good_credentials.c_str()),
        good_credentials.length());

    return (received_credentials == good_credentials);
}

} // namespace

bool HTTP_request_handler::handle(mi::http::IConnection* connection)
{
    mi::http::IRequest* request = connection->get_request();
    mi::http::IResponse* response = connection->get_response();
    const std::string url = request->get_url();

    const std::string auth_realm = "Viewer Authentication"; // Shown in the login dialog box in the browser

    // Do HTTP authentication, unless no username is set
    bool authenticated = false;

    const char* header_auth = request->get_header("Authorization");
    request->set_argument("session_cookie", 0);
    if (header_auth != 0 && !m_auth_user.empty())
    {
        // Header should contain "Basic <base64-encoded credentials>"
        const std::string auth = header_auth;
        const std::string auth_basic = "Basic ";
        if (auth.substr(0, auth_basic.length()) == auth_basic)
        {
            mi::base::Lock::Block block(&m_lock);
            const std::string received_credentials = auth.substr(auth_basic.length());

            // Authenticate if the credentials sent by the user match
            if (check_credentials(m_auth_user, m_auth_password, received_credentials))
            {
                authenticated = true;
                request->set_argument("session_cookie", m_session_cookie.c_str());
            }
            // Same for guest user
            else if (check_credentials(m_auth_guest_user, m_auth_guest_password, received_credentials))
            {
                authenticated = true;
                request->set_argument("session_cookie", m_session_cookie.c_str());
                request->set_argument("session_cookie", m_session_cookie_guest.c_str());
            }
            else
            {
                log_error(connection, "Authentication failed");
            }
        }
        else
        {
            // We only support "basic" authorization
            connection->print(html_status_page(response, 403, "Unsupported Authorization Scheme").c_str());
            log_error(connection, "Unsupported authorization scheme");
            return true;
        }
    }

    // Authentication is required unless no username is set
    if (!authenticated && !m_auth_user.empty())
    {
        response->set_header(
            "WWW-Authenticate",
            std::string("Basic realm=\"" + auth_realm + "\"").c_str());
        connection->print(html_status_page(response, 401, "Authorization Required").c_str());
        return true;
    }

    if (url == get_favicon_url())
    {
        // Icon that shows up in the web browser
        std::vector<mi::Uint8> data;
        data.assign(favicon, favicon + favicon_size);
        serve_buffer(connection, data, "image/x-icon");
    }
    else if (url == "/index_snapshot.png")
    {
        // Download the last snapshot image in PNG format (if available)
        serve_file(connection, "index_snapshot.png", "image/png");
    }
    else if (url == "/index_snapshot.ppm")
    {
        // Download the last snapshot image in PPM format (if available)
        serve_file(connection, "index_snapshot.ppm", "image/x-portable-pixmap");
    }
    else if (url == "/view" || url == "/index.html")
    {
        // Support deprecated URLs through a redirect
        response->set_header("Location", "/");
        connection->print(html_status_page(response, 307, "Temporary Redirect").c_str());
    }
    else if (url == "/crossdomain.xml")
    {
        // Cross-domain policy file to allow access to this server for Flash applications from any
        // domain. This is required so that a Flash viewer can query the status of other viewer
        // hosts.
        //
        // See also http://kb2.adobe.com/cps/142/tn_14213.html
        response->set_header("Content-Type", "text/xml");

        const std::string page =
            "<?xml version=\"1.0\"?>\n"
            "<!DOCTYPE cross-domain-policy SYSTEM "
            "\"http://www.macromedia.com/xml/dtds/cross-domain-policy.dtd\">\n"
            "<cross-domain-policy>\n"
            "<allow-access-from domain=\"*\" />\n"
            "</cross-domain-policy>\n";

        connection->print(page.c_str());
    }
    else
    {
        // Try default handler to handle the rest
        return Basic_HTTP_request_handler::handle(connection);
    }

    return true;
}

void HTTP_request_handler::set_authentication(
    const std::string& auth_user,
    const std::string& auth_password)
{
    mi::base::Lock::Block block(&m_lock);

    m_auth_user = auth_user;
    m_auth_password = auth_password;

    if (m_session_cookie.empty())
        ERROR_LOG << "Failed to create random session cookie, authentication for RTMP is disabled!";
}

void HTTP_request_handler::set_guest_authentication(
    const std::string& auth_guest_user,
    const std::string& auth_guest_password)
{
    mi::base::Lock::Block block(&m_lock);

    m_auth_guest_user = auth_guest_user;
    m_auth_guest_password = auth_guest_password;

    if (m_session_cookie_guest.empty())
        ERROR_LOG << "Failed to create random session cookie, authentication for RTMP is disabled!";
}

void HTTP_request_handler::set_resolutions(
    std::vector<mi::math::Vector<mi::Uint32, 2> >& resolutions)
{
    mi::base::Lock::Block block(&m_lock);
    m_resolutions = resolutions;
}

std::string HTTP_request_handler::get_session_cookie() const
{
    mi::base::Lock::Block block(&m_lock);
    return m_session_cookie;
}

std::string HTTP_request_handler::get_session_cookie_guest() const
{
    mi::base::Lock::Block block(&m_lock);
    return m_session_cookie_guest;
}
