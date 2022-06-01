/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 ******************************************************************************/
/// \file

#ifndef MI_XWIN_FASTXOPENDISPLAY
#define MI_XWIN_FASTXOPENDISPLAY

extern "C" {
#include <X11/Xlib.h>
}

#include <string>

/// X display connection information
struct Display_connection
{
    /// protocol string (empty for optimal local transport)
    std::string protocol;
    /// hostname string (empty for local host)
    std::string hostname;
    /// display number
    int display;
    /// screen number
    int screen;
};

/// Tests if the specified connection is possible.
/// Note that the screen has nothing to do with the connection and is ignored
/// (i.e. connection success does not mean that the screen exists).
bool test_connection(const Display_connection& connection);

/// An implementation of XOpenDisplay which should not block for too long.
Display *fast_XOpenDisplay(const char* connection_string);

#endif
