/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
///

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_COMMAND_HANDLER_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_COMMAND_HANDLER_H

#include <map>
#include <string>

#include <mi/dice.h>

#include "nvindex_appdata.h"
#include "nvindex_rendering_context.h"

/// Wrapper around received arguments and outgoing responses of the RTMP Call_event_handler.
class Call_event_arguments
{
public:
    Call_event_arguments(
        const mi::base::Handle<const mi::IMap>& args);

    /// Returns the value of the argument with the given key, or an empty string if there is none.
    std::string get(const std::string& key);

    /// Returns whether the given key exists.
    bool has(const std::string& key);

    /// Adds the key-value pair to the output, which is returnd by generate_response().
    void set(
        const std::string& key,
        const std::string& value);

    /// Creates the response from the key-value pairs.
    mi::IMap* generate_response(
        const mi::base::Handle<nv::index::IIndex>& iindex_if) const;

    /// Removes all key-value pairs created by set().
    void clear_response();

    /// Copies the given key-value pair from input to output, if it exists in the input.
    void pass_through(const std::string& key);

private:
    const mi::base::Handle<const mi::IMap>     m_input;
    std::map<std::string, std::string>         m_input_cache;
    std::map<std::string, std::string>         m_output_cache;
};

/// Interface for handling commands received by the RTMP Call_event_handler.
class ICall_event_command_handler
{
public:
    ICall_event_command_handler(
        Nvindex_rendering_context& irc,
        Nvindex_AppData*           appdata);

    virtual ~ICall_event_command_handler();

    virtual bool handle(
        const std::string&    cmd,
        Call_event_arguments& args) = 0;

protected:
    Nvindex_rendering_context& m_irc;
    Nvindex_AppData*           m_appdata;
};

#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_COMMAND_HANDLER_H
