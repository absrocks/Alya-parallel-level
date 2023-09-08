/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
///

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_SCENE_ELEMENT_IMPORTER_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_SCENE_ELEMENT_IMPORTER_H

#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <mi/dice.h>

#include <nv/index/iscene.h>

#include "common/forwarding_logger.h"
#include "common/string_dict.h"

class Scene_element_importer
{
public:
    Scene_element_importer(
        const mi::base::Handle<mi::neuraylib::IDice_transaction>& dice_transaction,
        const nv::index::IScene*                                  scene,
        mi::neuraylib::Tag                                        session_tag,
        std::map<std::string, mi::neuraylib::Tag>&                name_to_tag,
        std::vector<std::string>&                                 referencing_elements,
        const nv::index_common::String_dict&                      scene_dict);

    mi::neuraylib::Tag import(
        const std::string&             elem_name,
        const std::string&             elem_type,
        nv::index_common::String_dict& dict,
        bool&                          is_group,
        mi::neuraylib::Tag             old_tag = mi::neuraylib::Tag());

    static void get_supported_types(
        std::multimap<std::string, std::string>& types);

private:
    mi::base::Handle<nv::index::IScene_element> import_simple_shapes(
        const std::string&                   elem_name,
        const std::string&                   elem_type,
        const nv::index_common::String_dict& dict);

    mi::base::Handle<nv::index::IScene_element> import_raster_shapes(
        const std::string&                   elem_name,
        const std::string&                   elem_type,
        const nv::index_common::String_dict& dict);

    mi::base::Handle<nv::index::IScene_element> import_massive_shapes(
        const std::string&                                        elem_name,
        const std::string&                                        elem_type,
        nv::index_common::String_dict&                            dict, // intentionally non-const
        mi::neuraylib::Tag                                        old_tag,
        const mi::base::Handle<mi::neuraylib::IDice_transaction>& dice_transaction,
        bool&                                                     reuse_tag);

    mi::base::Handle<nv::index::IScene_element> import_attributes(
        const std::string&                   elem_name,
        const std::string&                   elem_type,
        const nv::index_common::String_dict& dict);

    mi::base::Handle<nv::index::IScene_element> import_raster_attributes(
        const std::string&                   elem_name,
        const std::string&                   elem_type,
        const nv::index_common::String_dict& dict);

    mi::base::Handle<nv::index::IScene_element> import_scene_groups(
        const std::string&                   elem_name,
        const std::string&                   elem_type,
        const nv::index_common::String_dict& dict);

    mi::neuraylib::Tag import_colormaps(
        const std::string&                   elem_name,
        const std::string&                   elem_type,
        const nv::index_common::String_dict& dict);

    inline mi::base::Handle<nv::index::IScene_element> null_result() const
    {
        return mi::base::Handle<nv::index::IScene_element>();
    }

    //
    // Shared helper functions
    //

    mi::neuraylib::Tag resolve_colormap(const std::string& id) const;

    mi::Float32 random(mi::Float32 min, mi::Float32 max) const;

    // Jet Color Map
    mi::math::Color_struct jetmap(
        mi::Float32 v,
        mi::Float32 vmin,
        mi::Float32 vmax) const;

    const mi::base::Handle<mi::neuraylib::IDice_transaction>& m_dice_transaction;
    const nv::index::IScene*                                  m_scene;
    mi::neuraylib::Tag                                        m_session_tag;
    std::map<std::string, mi::neuraylib::Tag>&                m_name_to_tag;
    std::vector<std::string>&                                 m_referencing_elements;
    const nv::index_common::String_dict&                      m_scene_dict;
    bool                                                      m_single_element;

    static const char*                                        s_supported_types[];
};

#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_SCENE_ELEMENT_IMPORTER_H
