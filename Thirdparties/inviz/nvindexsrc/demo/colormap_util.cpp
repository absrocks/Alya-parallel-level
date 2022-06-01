/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "colormap_util.h"

#include <map>
#include <cassert>
#include <vector>
#include <algorithm>
#include <iostream>

#include <nv/index/icolormap.h>

#include "common/colormap_io.h"

#include "colormap_manager.h"
#include "nvindex_appdata.h"


//----------------------------------------------------------------------
void set_colormap_tag(mi::Uint32 colormap_idx, const mi::neuraylib::Tag& tag)
{
    Nvindex_AppData::instance()->get_colormap_manager()->set_tag(colormap_idx, tag);
}

//----------------------------------------------------------------------
mi::neuraylib::Tag get_colormap_tag(mi::Uint32 colormap_idx)
{
    Colormap_manager* cmm = Nvindex_AppData::instance()->get_colormap_manager();
    assert(cmm != 0);
    if (!cmm->is_valid_tag(colormap_idx))
    {
        ERROR_LOG << "get_colormap_tag: invalid colormap index (" << colormap_idx << ")";
        return mi::neuraylib::NULL_TAG;
    }

    const mi::neuraylib::Tag tag = cmm->get_colormap_tag(colormap_idx);

    return tag;
}

//----------------------------------------------------------------------
mi::Uint32 get_number_of_colormap()
{
    mi::Uint32 colormap_count = 
        Nvindex_AppData::instance()->get_colormap_manager()->get_number_of_colormap();
    return colormap_count;
}

//----------------------------------------------------------------------
bool load_all_colormap(
    const nv::index::IScene*             scene,
    const nv::index_common::String_dict& app_proj,
    mi::neuraylib::IDice_transaction*    dice_transaction)
{
    using nv::index::IColormap;

    const mi::Sint32 cmap_count = nv::index_common::get_sint32(app_proj.get("app::colormap::file_count", "0"));

    // find colormap filenames from projects
    std::vector<std::string>                       colmap_fname_vec;
    std::vector<mi::math::Vector<mi::Float32, 2> > colmap_domain_vec;
    std::vector<IColormap::Domain_boundary_mode>   colmap_domain_boundary_mode_vec;
    const std::string prefix = "app::colormap::file_";
    for (mi::Sint32 i = 0; i < cmap_count; ++i)
    {
        std::stringstream sstr;
        sstr << prefix << i;
        const std::string key   = sstr.str();
        const std::string fname = app_proj.get(key, "<notfound>");
        if (fname == "<notfound>")
        {
            ERROR_LOG << "missing colormap file entry in project file [" << key << "]";
            return false;
        }
        colmap_fname_vec.push_back(fname);

        const std::string domain = app_proj.get(key + "::domain", "0 1");
        colmap_domain_vec.push_back(nv::index_common::get_vec_float32_2(domain));

        const std::string mode_str = app_proj.get(key + "::domain_boundary_mode", "clamp_to_edge");
        IColormap::Domain_boundary_mode mode = IColormap::CLAMP_TO_EDGE;
        if (mode_str == "clamp_to_transparent")
        {
            mode = IColormap::CLAMP_TO_TRANSPARENT;
        }
        else if (mode_str != "clamp_to_edge")
        {
            ERROR_LOG << "Invalid domain_boundary_mode '" << mode_str << "' for colormap '" << key << "'";
        }
        colmap_domain_boundary_mode_vec.push_back(mode);
    }

    for (mi::Sint32 i = 0; i < cmap_count; ++i)
    {
        const std::string cmap_fname = colmap_fname_vec.at(i);

        std::vector<mi::math::Color_struct> colormap_entries;
        colormap_entries.reserve(256);

        std::string error_mes;
        if (!nv::index_common::load_colormap(cmap_fname, colormap_entries, error_mes))
        {
            ERROR_LOG << error_mes;
            ERROR_LOG << "Fail to load colormap from file [" << cmap_fname << "]. Returning invalid colormap.";
            // set black 256 entries.
            mi::math::Color default_col(0.0f, 0.0f, 0.0f, 1.0f);
            colormap_entries.resize(256, default_col);
        }

        // Store the colormap to the database
        mi::base::Handle<IColormap> colormap(scene->create_attribute<IColormap>());
        colormap->set_colormap(&(colormap_entries[0]), colormap_entries.size());

        const mi::math::Vector<mi::Float32, 2>& domain = colmap_domain_vec.at(i);
        colormap->set_domain(domain.x, domain.y);
        colormap->set_domain_boundary_mode(colmap_domain_boundary_mode_vec.at(i));

        std::ostringstream colormap_name;
        colormap_name << "colormap_" << i;

        mi::neuraylib::Tag tag = dice_transaction->name_to_tag(colormap_name.str().c_str());
        if (tag.is_valid())
        {
            // This will happen when reconnecting to a cluster in fail safety mode, where the
            // existing colormap data should not be overwritten
            INFO_LOG << "Colormap '" << colormap_name.str() << "' (tag " << tag << ") is already stored "
                     << "in database, using existing data";
        }
        else
        {
            tag = dice_transaction->store(
                colormap.get(), mi::neuraylib::Tag(), colormap_name.str().c_str());
            assert(tag.is_valid());
            if (!tag.is_valid())
            {
                ERROR_LOG << "load_all_colormap: Failed to store colormap job.";
                exit(1);
            }
        }

        Nvindex_AppData::instance()->get_colormap_manager()->set_tag(i, tag);
    }
    return true;
}

//----------------------------------------------------------------------
void dump_colormap(
    const nv::index::IColormap* cmap)
{
    mi::Uint32 const cmap_count = static_cast<mi::Uint32>(cmap->get_number_of_entries());
    std::string const vecname = "colormap_entry";
    std::cout << "std::vector< mi::math::Color > " << vecname << ";" << std::endl;

    for (mi::Uint32 i = 0; i < cmap_count; ++i)
    {
        mi::math::Color_struct const colst = cmap->get_color(i);
        std::cout << vecname << ".push_back(mi::math::Color(" << colst.r << ", " << colst.g << ", "
                  << colst.b << ", " << colst.a << "));" << std::endl;
    }
}

//----------------------------------------------------------------------
void dump_all_colormap(
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    Colormap_manager* cmm = Nvindex_AppData::instance()->get_colormap_manager();
    const mi::Uint32 cmap_count = cmm->get_number_of_colormap();
    for (mi::Uint32 i = 0; i < cmap_count; ++i)
    {
        const mi::neuraylib::Tag cmap_tag = cmm->get_colormap_tag(i);
        mi::base::Handle<const nv::index::IColormap> cmap(
            dice_transaction->access<nv::index::IColormap>(cmap_tag));
        if (!cmap.is_valid_interface())
        {
            WARN_LOG << "colormap id: " << i << " -> tag: " << cmap_tag.id << " is not valid.";
        }
        else
        {
            dump_colormap(cmap.get());
        }
    }
}

//----------------------------------------------------------------------
bool check_all_colormap_valid(
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    Colormap_manager* cmm = Nvindex_AppData::instance()->get_colormap_manager();
    mi::Uint32 const cmap_count = cmm->get_number_of_colormap();
    mi::Sint32 invalid_colormap_tag_count = 0;
    for (mi::Uint32 i = 0; i < cmap_count; ++i)
    {
        const mi::neuraylib::Tag cmap_tag = cmm->get_colormap_tag(i);
        mi::base::Handle<const nv::index::IColormap> cmap(
            dice_transaction->access<nv::index::IColormap>(cmap_tag));
        if (!cmap.is_valid_interface())
        {
            WARN_LOG << "colormap id: " << i << " -> tag: " << cmap_tag.id << " is not valid.";
            ++invalid_colormap_tag_count;
        }
    }

    if (invalid_colormap_tag_count > 0)
    {
        WARN_LOG << invalid_colormap_tag_count << " colormap tags are not valid";
        return false;
    }
    INFO_LOG << "Created " << cmap_count << " user-defined colormaps";

    // dump_all_colormap(dice_transaction);

    return true;
}

//----------------------------------------------------------------------
bool exec_userdef_colormap_operation(
    const mi::neuraylib::Tag          cmap_tag,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    {
        mi::base::Handle< nv::index::IColormap> cmap(
            dice_transaction->edit<nv::index::IColormap>(cmap_tag));
        assert(cmap.is_valid_interface());

        Colormap_manager* cmm = Nvindex_AppData::instance()->get_colormap_manager();
        const mi::Uint32 op_type = cmm->get_colormap_userdef_operation_type();
        if (op_type == 0)
        {
            // colormap reverse
            INFO_LOG << "exec_userdef_colormap_operation: colormap reverse.";
            mi::Sint32 const array_count = static_cast< mi::Sint32 >(cmap->get_number_of_entries());
            std::vector< mi::math::Color_struct > colmap_op_buf;
            for (mi::Sint32 i = array_count - 1; i >= 0; --i)
            {
                colmap_op_buf.push_back(cmap->get_color(i));
            }
            for (mi::Sint32 i = 0; i < array_count; ++i)
            {
                cmap->set_color(i, colmap_op_buf.at(i));
            }
        }
        else if (op_type == 1)
        {
            // colormap random shuffle
            INFO_LOG << "exec_userdef_colormap_operation: colormap random shuffle.";
            const mi::Sint32 array_count = static_cast< mi::Sint32 >(cmap->get_number_of_entries());
            std::vector< mi::math::Color_struct > colmap_op_buf;
            for (mi::Sint32 i = 0; i < array_count; ++i)
            {
                colmap_op_buf.push_back(cmap->get_color(i));
            }
            std::random_shuffle(colmap_op_buf.begin(), colmap_op_buf.end());
            for (mi::Sint32 i = 0; i < array_count; ++i)
            {
                cmap->set_color(i, colmap_op_buf.at(i));
            }
        }
        else
        {
            ERROR_LOG << "exec_userdef_colormap_operation: No such colormap operation id: "
                      << op_type << ", ignored.";
            return false;
        }
    }

    return true;
}

//----------------------------------------------------------------------
