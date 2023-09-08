/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 ******************************************************************************/
/// \file

#include "xwin_fastxopendisplay.h"

#include <mi/base/lock.h>

#include "common/forwarding_logger.h"
#include "common/string_dict.h"

extern "C" {
#include "Xtrans.h"
#include <X11/Xlibint.h>
}

#include <string>
#include <map>
#include <cctype>
#include <cstdio>
#include <sstream>

using std::string;


// A map which is used as a cache to quickly check if a display is working or not. Along with a
// dedicated lock for the map.
static mi::base::Lock s_lock;
typedef std::map<std::string, bool> Display_map;
static Display_map s_displays;


// Tests if the specified connection is possible.
// Note that the screen has nothing to do with that and it is ignored
// (i.e. connection success does not mean that the screen exists).
bool test_connection(const Display_connection& connection)
{
    XtransConnInfo trans_conn = NULL;
    mi::Sint32 connect_stat = 0;

    std::ostringstream oss;
    if (!connection.protocol.empty()) {
        oss << connection.protocol << "/";
    }

    // the screen is ommited intentionally, Trans knows nothing about screens.
    oss << connection.hostname << ":" << connection.display;

    char address[256];
    snprintf(address, sizeof address, "%s", oss.str().c_str());

    if ((trans_conn = mi_XTransOpenCOTSClient(address)) == NULL) {
        return false;
    }

    if ((connect_stat = mi_XTransConnect(trans_conn, address)) < 0 ) {
        mi_XTransClose(trans_conn);
        trans_conn = NULL;

        // Although chances are slim, this is a point that our trick might fail.
        // TRANS_TRY_CONNECT_AGAIN means "try again after a delay" and since we
        // don't want to delay, we ignore it. If we ever run into the situation
        // that connecting to the X server is impossible while everything seems
        // properly set, but at the same time we get this debug message, we may
        // have to consider switching off the whole fast_XOpenDisplay mechanism.
        // Note that we WILL get TRANS_TRY_CONNECT_AGAIN if the protocol is TCP
        // and the connection fails, but TCP is a fallback for XOpenDisplay and
        // will only be used if UNIXCONN also failed, which could happen if for
        // example the application can't access "/tmp/.X11-unix/".
        if (connect_stat == TRANS_TRY_CONNECT_AGAIN) {
            INFO_LOG << "TRANS_TRY_CONNECT_AGAIN was ignored.";
        }
        return false;
    }
    mi_XTransDisconnect(trans_conn);
    mi_XTransClose(trans_conn);
    return true;
}


static bool parse_connection_string(
    const char *display_name, Display_connection *dc, bool *is_local=NULL)
{
    // based on code from libX11's ConnDis.c

    Display_connection ret;             // return value
    const char *lastp, *lastc, *p;      // char pointers 
    mi::Sint32 hostlen;                 // length tmp variable

    p = display_name;

    //
    // Step 0, find the protocol.  This is delimited by the optional
    // slash ('/').
    //
    for (lastp = p; *p && *p != ':' && *p != '/'; p++) ;
    if (!*p)
        return false;                     // must have a colon

    if (p != lastp && *p != ':') {        // protocol given?
        ret.protocol = string(lastp, p - lastp);
        p++;                              // skip the '/'
    }
    else {
        p = display_name;                 // reset the pointer in case no protocol was given
    }

    //
    // Step 1, find the hostname.  This is delimited by either one colon,
    // or two colons in the case of DECnet (DECnet Phase V allows a single
    // colon in the hostname).  (See note above regarding IPv6 numeric
    // addresses with triple colons or [] brackets.)
    //
    lastp = p;
    lastc = NULL;
    for (; *p; p++)
        if (*p == ':')
            lastc = p;

    if (!lastc)
        return false;      // must have a colon

    if ((lastp != lastc) && (*(lastc - 1) == ':')
#if defined(IPv6) && defined(AF_INET6)
      && ( ((lastc - 1) == lastp) || (*(lastc - 2) != ':'))
#endif  // defined(IPv6) && defined(AF_INET6)
       )
    {
        // DECnet display specified (we won't support that)
        return false;
    }
    else {
        hostlen = lastc - lastp;
    }

    if (hostlen > 0) {                  // hostname given?
        ret.hostname = string(lastp, hostlen);
    }

    p = lastc;

    // find if the connection is to the local host.
    if (is_local) {
        char localhost[256];
        *is_local = (_XGetHostname(localhost, sizeof localhost) > 0) && ret.hostname == localhost;
    }

    //
    // Step 2, find the display number.  This field is required and is
    // delimited either by a nul or a period, depending on whether or not
    // a screen number is present.
    //
    for (lastp = ++p; *p && isascii(*p) && isdigit(*p); p++) ;
    if ((p == lastp) ||                 // required field
        (*p != '\0' && *p != '.'))      // invalid non-digit terminator
    {
        return false;
    }
    ret.display = nv::index_common::get_sint32(string(lastp, p - lastp));

    //
    // Step 3, find the screen number.  This field is optional.  It is
    // present only if the display number was followed by a period (which
    // we've already verified is the only non-nul character).
    //
    if (*p) {
        for (lastp = ++p; *p && isascii(*p) && isdigit (*p); p++) ;
        if (p != lastp) {
            if (*p) {                   // non-digits
                return false;
            }
            ret.screen = nv::index_common::get_sint32(lastp);
        }
    }

    *dc = ret;
    return true;
}


// An implementation of XOpenDisplay which should not block for too long.
Display* fast_XOpenDisplay(const char* display_string)
{
    INFO_LOG << "Called fast_XOpenDisplay("<<display_string<<")";

    Display *disp = NULL;
    {
        bool first_time = true;
        bool previous_attempt_successful = false;
        bool connection_test_successful = false;

        mi::base::Lock::Block block(&s_lock);
        if (s_displays.find(display_string) != s_displays.end()) {
            // We know already whether this works or not from a previous attempt.
            first_time = false;
            previous_attempt_successful = s_displays[display_string];
        }
        else {
            // This is the first time we attempt to open the display.
            Display_connection connection;
            if (parse_connection_string(display_string, &connection)) {
                connection_test_successful = test_connection(connection);
                if (!connection_test_successful) {
                    // The connection failed, we assume that XOpenDisplay will also fail,
                    // we would rather not waste our time with it.
                    INFO_LOG << 
                        "test_connection("<<display_string<<
                        ") failed, it will be registered as non-working.";
                    s_displays[display_string] = false;
                }
            }
            else {
                INFO_LOG << 
                    "\""<<display_string<<"\" is not a valid string for an X Server connection.";
            }
        }
 
        if (!first_time) {
            if (previous_attempt_successful) {
                // An earlier attempt to open this display was successful.
                disp = XOpenDisplay(display_string);
                if (disp == NULL) {
                    // The X server probably crashed.
                    INFO_LOG << 
                        "XOpenDisplay("<<display_string<<") failed, but was successful previously. "
                        "Is the X server still alive?";
                }
            }
            else {
                INFO_LOG << 
                    "XOpenDisplay("<<display_string<<") will be skipped, as it failed previously.";
            }
        }
        else
        if (connection_test_successful) {
            // A connection is possible, but that's all we know so far.
            // Let's see if XOpenDisplay will also be successful.
            VERBOSE_LOG << "Attempting XOpenDisplay("<<display_string<<")";
            disp = XOpenDisplay(display_string);
            if (disp == NULL) {
                // test_connection() wasn't enough. This can happen if the screen does not exist,
                // but in this case XOpenDisplay returns immediately.
                INFO_LOG << 
                    "The connection test was successful, but XOpenDisplay failed.";
            }
            s_displays[display_string] = (disp != NULL);
        }

        if (disp != NULL) {
            if (disp->trans_conn != NULL)
                INFO_LOG << "XOpenDisplay("<<display_string<<") connected successfully using a "
                         << (mi_XTransIsLocal(disp->trans_conn)?"":"non-") << "local transport.";
            else
                INFO_LOG << "XOpenDisplay("<<display_string<<") connected, "
                         << "but transport information is not available.";
        }
    }

    return disp;
}

