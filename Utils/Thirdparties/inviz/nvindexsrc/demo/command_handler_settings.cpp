/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/

#include "command_handler_settings.h"

#include "scene_utility.h"

using namespace nv::index_common;

bool Command_handler_settings::handle(
    const std::string&    cmd,
    Call_event_arguments& args)
{
    if (cmd == "get_video_codec")
    {
        args.set("cmd", cmd);

        String_dict* prj = Nvindex_AppData::instance()->peek_app_proj();
        args.set("codec", prj->get("dice::rtmp_video_streaming::video_codec"));
        args.set("bitrate", prj->get("dice::rtmp_video_streaming::video_bitrate"));

        return true;
    }
    else if (cmd == "set_video_codec")
    {
        args.set("cmd", cmd);

        String_dict* prj = Nvindex_AppData::instance()->peek_app_proj();

        // Video codec
        const std::string video_codec = args.get("codec");
        if (!video_codec.empty())
        {
            INFO_LOG << "Setting video codec to '" << video_codec << "'";
            prj->insert("dice::rtmp_video_streaming::video_codec", video_codec);
        }

        // Video bitrate
        const mi::Uint32 video_bitrate = get_uint32(args.get("bitrate"));
        if (video_bitrate > 0)
        {
            std::ostringstream os;
            os << (video_bitrate * 1000000); // Mbit -> bit
            INFO_LOG << "Setting video bitrate to " << os.str();
            prj->insert("dice::rtmp_video_streaming::video_bitrate", os.str());
        }

        return true;
    }

    return false;
}
