/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/

#include "command_handler_colormap.h"

#include "nv/index/icolormap.h"

#include "colormap_manager.h"

using namespace nv::index_common;

namespace {

std::string quote(const std::string& s)
{
    return '"' + escape_JSON(s) + '"';
}

} // namespace

bool Command_handler_colormap::handle(
    const std::string&    cmd,
    Call_event_arguments& args)
{
    if (cmd == "get_colormaps")
    {
        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            m_irc.get_current_scope()->create_transaction<mi::neuraylib::IDice_transaction>());

        mi::base::Handle<const nv::index::ISession> session(
            dice_transaction->access<nv::index::ISession>(m_irc.m_session_tag));

        const Colormap_manager* cmm = Nvindex_AppData::instance()->get_colormap_manager();
        const mi::Uint32 nb_colormaps = cmm->get_number_of_colormap();

        std::ostringstream os;
        os << "[";
        for (mi::Uint32 i=0; i < nb_colormaps; ++i)
        {
            if (i != 0)
                os << ",";

            mi::neuraylib::Tag tag = cmm->get_colormap_tag(i);

            std::ostringstream label;

            const char* tag_name = dice_transaction->tag_to_name(tag);
            if (tag_name != 0)
            {
                label << tag_name << " (" << tag << ")";
            }
            else
            {
                label << "colormap #" << i << " (" << tag << ")";
            }

            os << "{"
               << quote("label") << ":" << quote(label.str()) << ","
               << quote("tag") << ":" << tag
               << "}";
        }
        os << "]";

        args.set("cmd", cmd);
        args.set("colormaps", os.str().c_str());

        dice_transaction->commit();
        return true;
    }
    else if (cmd == "get_colormap")
    {
        mi::neuraylib::Tag tag;
        std::istringstream is(args.get("tag"));
        is >> tag.id;

        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            m_irc.get_current_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
        mi::base::Handle<const nv::index::ISession> session(
            dice_transaction->access<nv::index::ISession>(m_irc.m_session_tag));

        mi::base::Handle<const nv::index::IColormap> colormap(
            dice_transaction->access<nv::index::IColormap>(tag));
        if (colormap)
        {
            const mi::Uint32 nb_colormap_entries = static_cast<mi::Uint32>(colormap->get_number_of_entries());
            const mi::math::Color_struct* color_table = colormap->get_colormap();

            std::ostringstream os;
            for (mi::Uint32 i = 0; i < nb_colormap_entries; ++i)
            {
                if (i > 0)
                    os << " ";
                const mi::math::Color_struct& entry = color_table[i];
                os << entry.r << " " << entry.g << " " << entry.b << " " << entry.a;
            }
            args.set("values", os.str().c_str());
        }
        else
        {
            ERROR_LOG << "get_colormap: Could not access colormap with tag " << tag;
        }

        args.set("cmd", cmd);

        dice_transaction->commit();
        return true;
    }
    else if (cmd == "set_colormap")
    {
        mi::neuraylib::Tag tag;
        std::istringstream is(args.get("tag"));
        is >> tag.id;

        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            m_irc.get_current_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
        mi::base::Handle<const nv::index::ISession> session(
            dice_transaction->access<nv::index::ISession>(m_irc.m_session_tag));

        {
            mi::base::Handle<nv::index::IColormap> colormap(
                dice_transaction->edit<nv::index::IColormap>(tag));
            if (colormap)
            {
                const mi::Uint32 nb_colormap_entries = static_cast<mi::Uint32>(colormap->get_number_of_entries());
                std::vector<mi::math::Color_struct> color_table(nb_colormap_entries);

                const std::string values = args.get("values");
                std::istringstream is(values);

                for (mi::Uint32 i=0; i < nb_colormap_entries; ++i)
                {
                    mi::math::Color_struct& entry = color_table[i];
                    is >> entry.r >> entry.g >> entry.b >> entry.a;
                }
                colormap->set_colormap(&color_table[0], color_table.size());
            }
            else
            {
                ERROR_LOG << "set_colormap: Could not edit colormap with tag " << tag;
            }
        }
        args.set("cmd", cmd);

        dice_transaction->commit();
        return true;
    }

    return false;
}
