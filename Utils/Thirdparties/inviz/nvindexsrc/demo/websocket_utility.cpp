/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "websocket_utility.h"

#include "common/forwarding_logger.h"

#include "html5_video_stream.h"
#include "nvindex_rendering_context.h"

#include <string>
#include <cassert>

//======================================================================
Websocket_video_data_buffer::Websocket_video_data_buffer(mi::neuraylib::IVideo_data* data)
    :
    m_data(data)
{
    // empty
}

//----------------------------------------------------------------------
Websocket_video_data_buffer::~Websocket_video_data_buffer()
{
    m_data = 0;
}

//----------------------------------------------------------------------
const mi::Uint8* Websocket_video_data_buffer::get_data() const
{
    if(m_data.is_valid_interface())
    {
        return m_data->get_data();
    }
    else
    {
        return 0;
    }
}

//----------------------------------------------------------------------
const mi::Uint8* Websocket_video_data_buffer::get_data()
{
    if(m_data.is_valid_interface())
    {
        return m_data->get_data();
    }
    else
    {
        return 0;
    }
}

//----------------------------------------------------------------------
mi::Size Websocket_video_data_buffer::get_data_size() const
{
    if(m_data.is_valid_interface())
    {
        return m_data->get_data_size();
    }
    else
    {
        return 0;
    }
}

//======================================================================
