/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
///

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_COMMAND_HANDLER_COLORMAP_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_COMMAND_HANDLER_COLORMAP_H

#include "command_handler_base.h"

class Command_handler_colormap : public ICall_event_command_handler
{
public:
    Command_handler_colormap(
        Nvindex_rendering_context& irc,
        Nvindex_AppData*           appdata)
      : ICall_event_command_handler(irc, appdata)
    {}

    virtual bool handle(
        const std::string&    cmd,
        Call_event_arguments& args);
};

#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_COMMAND_HANDLER_COLORMAP_H
