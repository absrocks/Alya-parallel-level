/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Basic HTTP server implementation

#ifndef NVIDIA_INDEX_BIN_COMMON_BASIC_HTTP_HANDLER_H
#define NVIDIA_INDEX_BIN_COMMON_BASIC_HTTP_HANDLER_H

#include <mi/base/interface_implement.h>

#include <mi/neuraylib/http.h>
#include <mi/base/handle.h>
#include <mi/math/vector.h>

#include <map>
#include <string>
#include <vector>
#include <sstream>
#include <cstdio>

#include "forwarding_logger.h"
#include "swfobject_js.h"

namespace nv {
namespace index_common {

namespace {

/// A simple Uint8 buffer
class Simple_buffer : public mi::base::Interface_implement<mi::neuraylib::IBuffer>
{
  public:
    Simple_buffer(const std::vector<mi::Uint8>& content) { m_buffer = content; }

    const mi::Uint8* get_data() const { return &m_buffer[0]; }
    mi::Size get_data_size() const { return m_buffer.size(); }

  private:
    std::vector<mi::Uint8> m_buffer;
};

} // namespace

/// Basic HTTP server implementation
class Basic_HTTP_request_handler : public mi::base::Interface_implement<mi::http::IRequest_handler>
{
public:
    Basic_HTTP_request_handler(
        const std::string&                              main_viewer_file,
        std::vector<mi::math::Vector<mi::Uint32, 2> >&  resolutions,
        const std::string&                              rtmp_port = "")
        : m_main_viewer_file(main_viewer_file),
          m_resolutions(resolutions),
          m_rtmp_port(rtmp_port)
    {
        // empty
    }

    virtual bool handle(mi::http::IConnection* connection)
    {
        mi::http::IRequest* request = connection->get_request();
        mi::http::IResponse* response = connection->get_response();
        const std::string url = request->get_url();

        const std::string first_view = "/" + get_flash_viewer_url();
        const std::string second_view = "/1/" + get_flash_viewer_url();

        if (url == "/")
        {
            // Output the html page
            response->set_header("Content-Type", "text/html");
            response->set_header("Cache-Control", "no-cache");

            request->set_argument("viewing_id", "0"); // pass to get_flash_vars()
            const std::string page = html_page(request);
            connection->print(page.c_str());
        }
        else if (url == first_view)
        {
            // Output the flash file
            serve_file(connection, m_main_viewer_file, "application/x-shockwave-flash");
        }
        else if (url == second_view)
        {
            // Output the flash file
            serve_file(connection, m_main_viewer_file, "application/x-shockwave-flash");
        }
        else if (url == get_swfobject_url())
        {
            // Javascript for including flash
            std::vector<mi::Uint8> data;
            data.assign(swfobject_js, swfobject_js + sizeof(swfobject_js) - 1);
            serve_buffer(connection, data, "application/x-javascript");
        }
        else if (url == "/1/")
        {
            // Output the html page
            response->set_header("Content-Type", "text/html");
            response->set_header("Cache-Control", "no-cache");

            request->set_argument("viewing_id", "1"); // pass to get_flash_vars()
            const std::string page = html_page(request);
            connection->print(page.c_str());
        }
        else
        {
            // URL not found
            connection->print(html_status_page(response, 404, "Not Found").c_str());
            log_error(connection, "URL '" + url + "' not found");
        }

        return true;
    }

protected:
    /// Returns the string to use for the "title" tag.
    virtual std::string get_page_title() const { return "NVIDIA IndeX viewer"; }

    /// Returns the URL to the flash viewer (*.swf file).
    virtual std::string get_flash_viewer_url() const { return "viewer.swf"; }

    /// Returns the URL to the swfobject.js file.
    virtual std::string get_swfobject_url() const { return "/swfobject.js"; }

    /// Minimum required Flash Player version, will show a warning message box if not available
    virtual std::string get_required_flash_version() const { return "11.0"; }

    /// Returns a formatted list of the flash variables.
    virtual void get_flash_vars(
        mi::http::IRequest*                 request,
        std::map<std::string, std::string>& flash_vars)
    {
        mi::Uint32 viewing_id = 0;
        const char* viewing_id_str = request->get_argument("viewing_id");
        if (viewing_id_str != 0)
        {
            std::istringstream is(viewing_id_str);
            is >> viewing_id;
        }

        const mi::math::Vector<mi::Uint32, 2>& resolution = m_resolutions[viewing_id];

        // Pass video size to the flash player
        if (resolution.x > 0 && resolution.y > 0)
        {
            {
                std::ostringstream s;
                s << resolution.x;
                flash_vars["    video_width"] = s.str();
            }

            {
                std::ostringstream s;
                s << resolution.y;
                flash_vars["    video_height"] = s.str();
            }
        }

        // Pass non-default RTMP port to the flash player
        if (!m_rtmp_port.empty() && m_rtmp_port != "1935")
        {
            flash_vars["    rtmp_port"] = m_rtmp_port;
        }

        {
            std::ostringstream s;
            s << viewing_id;
            flash_vars["    viewing_session"] = s.str();
        }
    }

    /// Returns extra items to add in the "head".
    virtual std::string get_extra_head_items() const { return ""; }

    /// Returns the complete HTML page.
    virtual std::string html_page(mi::http::IRequest* request)
    {
        // Identifier, can be used by viewer application to verify that this is actually an NVIDIA
        // IndeX server.
        std::string identifier = "<meta name=\"generator\" content=\"NVIDIA IndeX Web Viewer\">";

        std::string html = std::string() +
            "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.01//EN\" \"http://www.w3.org/TR/html4/strict.dtd\">\n"
            "<html>\n"
            "<head>\n"
            "<meta http-equiv=\"Content-Type\" content=\"text/html;charset=UTF-8\">\n"
            "<title>" + get_page_title() + "</title>\n"
            + identifier + "\n"
            + get_extra_head_items() +

            // Set up CSS so that the flash plugin uses the entire page.
            // See http://dev.firefallpro.com/full-window-flash/ for details.
            "<style type=\"text/css\">\n"
            "<!--\n"
            "  html, body, div#content, div#content object { width: 100%; height: 100%; }\n"
            "  body { padding: 0; margin: 0; }\n"
            "  div#content, div#content object { overflow: hidden; }\n"
            "-->\n</style>\n"

            + html_include_flash(request) +

            "</head>\n"
            "<body>\n"
            "<div id=\"content\">\n"
            "  <div id=\"flash\">Flash Player " + get_required_flash_version() + " is required to view this website.</div>\n"
            "</div>\n"
            "</body>\n"
            "</html>\n";

        return html;
    }

    /// Returns the HTML code that includes the flash viewer.
    virtual std::string html_include_flash(mi::http::IRequest* request)
    {
        // Prepare variable to pass to the flash applet
        std::string flash_vars;
        std::map<std::string, std::string> flash_vars_map;
        get_flash_vars(request, flash_vars_map);
        std::map<std::string, std::string>::const_iterator it;
        for (it = flash_vars_map.begin(); it != flash_vars_map.end(); ++it)
        {
            flash_vars += it->first + ": " + it->second + ",\n";
        }

        std::string html = std::string() +
            // Include SWFObject JavaScript library for handling Flash content
            "<script type=\"text/javascript\" src=\"" + get_swfobject_url() + "\"></script>\n"

            "<script type=\"text/javascript\">\n"
            // Check Flash version
            "  var required = \"" + get_required_flash_version() + "\";\n"
            "  var available = swfobject.getFlashPlayerVersion();\n"
            "  if (!swfobject.hasFlashPlayerVersion(required))\n"
            "      alert(\"This application requires Flash Player \" + required + \n"
            "            \" but only version \" + available.major + \".\" + available.minor + \n"
            "            \" is available.\");\n"

            // Embed Flash Player using SWFObject
            "  var flashvars = {\n"
            + flash_vars +
            "    eof: true\n"
            "  };\n"
            // Only require Flash 10 here, as we want to try loading it anyway, after showing the
            // alert to the user
            "  swfobject.embedSWF(\"" + get_flash_viewer_url() + "\", \"flash\", \"100%\", \"100%\", "
            "\"10\", false, flashvars);\n"
            "</script>\n";

        return html;
    }

    /// Returns a HTML status page for showing errors.
    virtual std::string html_status_page(
        mi::http::IResponse* response,
        mi::Sint32           status_code,
        const std::string&   status_description,
        const std::string&   message = "")
    {
        response->set_header("Content-Type", "text/html");
        response->set_result_code(status_code, status_description.c_str());

        std::ostringstream s;
        s << "<html><head><title>" << status_code << " " << status_description << "</title></head>"
          << "<body>"
          << "<h1>" << status_code << " " << status_description << "</h1>"
          << message
          << "</body></html>";

        return s.str();
    }

    /// Logs an error (adding information about the client).
    virtual void log_error(
        mi::http::IConnection* connection,
        const std::string&     message)
    {
        ERROR_LOG << "HTTP request [" << connection->get_peer_address()  << "]: " << message;
    }

    /// Sends a data buffer to the client.
    virtual void serve_buffer(
        mi::http::IConnection*        connection,
        const std::vector<mi::Uint8>& data,
        const std::string&            content_type = "application/octet-stream")
    {
        mi::http::IResponse* response = connection->get_response();
        response->set_header("Content-Type", content_type.c_str());

        mi::base::Handle<mi::neuraylib::IBuffer> buffer(new Simple_buffer(data));
        connection->enqueue(buffer.get());
    }

    /// Sends a file to the client.
    virtual void serve_file(
        mi::http::IConnection* connection,
        const std::string&     file_name,
        const std::string&     content_type = "application/octet-stream")
    {
        mi::http::IResponse* response = connection->get_response();
        FILE * p_file = fopen(file_name.c_str(), "rb");
        if (p_file != 0){
            fseek(p_file, 0L, SEEK_END);
            mi::Sint64 const file_size = ftell(p_file);
            std::vector<mi::Uint8> data(file_size);
            rewind(p_file);
            fread((void*)&(data[0]), sizeof(mi::Uint8), file_size, p_file);
            fclose(p_file);

            serve_buffer(connection, data, content_type);
        }
        else
        {
            log_error(connection, "File '" + file_name + "' could not be found.");
            response->set_header("Content-Type", "text/html");
            response->set_result_code(500, "Internal Server Error");
            connection->print(
                html_status_page(
                    response,
                    500,
                    "Internal Server Error",
                    "System file could not be found.").c_str());
        }
    }

protected:
    std::string                                     m_main_viewer_file;
    std::vector<mi::math::Vector<mi::Uint32, 2> >   m_resolutions;
    std::string                                     m_rtmp_port;
};

}} // namespace nv::index_common

#endif // NVIDIA_INDEX_BIN_COMMON_BASIC_HTTP_HANDLER_H
