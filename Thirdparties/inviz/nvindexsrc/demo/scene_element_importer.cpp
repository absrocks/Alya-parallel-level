/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/

#include "scene_element_importer.h"

#include <nv/index/icolormap.h>

#include "colormap_util.h"
#include "nvindex_appdata.h"

using namespace nv::index_common;

namespace {

const char* TYPE_ATTRIBUTE = "Attributes";
const char* TYPE_GROUPS    = "Scene Groups";
const char* TYPE_MASSIVE   = "Massive Data";
const char* TYPE_RASTER    = "Raster Shapes";
const char* TYPE_SHAPE     = "Shapes";

} // namespace

const char* Scene_element_importer::s_supported_types[] = {
    "circle",                               TYPE_RASTER,
    "clip_region",                          TYPE_ATTRIBUTE,
    "colormap",                             TYPE_ATTRIBUTE,
    "cone",                                 TYPE_SHAPE,
    "compute_technique",                    TYPE_ATTRIBUTE,
    "cylinder",                             TYPE_SHAPE,
    "depth_offset",                         TYPE_ATTRIBUTE,
    "directional_headlight",                TYPE_ATTRIBUTE,
    "directional_light",                    TYPE_ATTRIBUTE,
    "ellipse",                              TYPE_RASTER,
    "ellipsoid",                            TYPE_SHAPE,
    "heightfield",                          TYPE_MASSIVE,
    "heightfield_geometry",                 TYPE_ATTRIBUTE,
    "heightfield_wireframe_style",          TYPE_ATTRIBUTE,
    "horizon",                              TYPE_MASSIVE,   // Deprecated, use "heightfield" instead
    "horizon_geometry",                     TYPE_ATTRIBUTE, // Deprecated, use "heightfield_geometry" instead
    "icon_2d",                              TYPE_RASTER,
    "icon_3d",                              TYPE_RASTER,
    "intersection_highlighting",            TYPE_ATTRIBUTE,
    "irregular_volume",                     TYPE_MASSIVE,
    "label_2d",                             TYPE_RASTER,
    "label_3d",                             TYPE_RASTER,
    "label_font",                           TYPE_ATTRIBUTE,
    "label_layout",                         TYPE_ATTRIBUTE,
    "layered_graphics",                     TYPE_ATTRIBUTE,
    "line_path_2d",                         TYPE_RASTER,
    "line_path_3d",                         TYPE_RASTER,
    "line_set",                             TYPE_RASTER,
    "path_style",                           TYPE_ATTRIBUTE,
    "phong_gl",                             TYPE_ATTRIBUTE,
    "pipe_set",                             TYPE_MASSIVE,
    "plane",                                TYPE_SHAPE,
    "point_set",                            TYPE_RASTER,
    "polygon",                              TYPE_RASTER,
    "raster_benchmark_lines",               TYPE_RASTER,
    "raster_benchmark_points",              TYPE_RASTER,
    "rendering_kernel_program",             TYPE_ATTRIBUTE,
    "rendering_kernel_program_parameters",  TYPE_ATTRIBUTE,
    "reservoir_colormap",                   TYPE_ATTRIBUTE,
    "reservoir_grid",                       TYPE_MASSIVE,
    "shading_model",                        TYPE_ATTRIBUTE,
    "sparse_volume",                        TYPE_MASSIVE,
    "sparse_volume_rendering_properties",   TYPE_ATTRIBUTE,
    "sphere",                               TYPE_SHAPE,
    "static_scene_group",                   TYPE_GROUPS,
    "texture",                              TYPE_ATTRIBUTE,
    "texture_filter_mode",                  TYPE_ATTRIBUTE,
    "texturing_technique",                  TYPE_ATTRIBUTE,
    "transformed_scene_group",              TYPE_GROUPS,
    "time_step_assignment",                 TYPE_ATTRIBUTE,
    "triangle_mesh",                        TYPE_MASSIVE,
    "volume",                               TYPE_MASSIVE,
    "volume_rendering_properties",          TYPE_ATTRIBUTE,
    "volume_texture",                       TYPE_ATTRIBUTE,
    "wireframe_rendering_style",            TYPE_ATTRIBUTE,
    0 // This must be the last entry
};

Scene_element_importer::Scene_element_importer(
    const mi::base::Handle<mi::neuraylib::IDice_transaction>& dice_transaction,
    const nv::index::IScene*                                  scene,
    mi::neuraylib::Tag                                        session_tag,
    std::map<std::string, mi::neuraylib::Tag>&                name_to_tag,
    std::vector<std::string>&                                 referencing_elements,
    const nv::index_common::String_dict&                      scene_dict)
  : m_dice_transaction(dice_transaction),
    m_scene(scene),
    m_session_tag(session_tag),
    m_name_to_tag(name_to_tag),
    m_referencing_elements(referencing_elements),
    m_scene_dict(scene_dict),
    m_single_element(false)
{
}

mi::neuraylib::Tag Scene_element_importer::import(
    const std::string& elem_name,
    const std::string& elem_type,
    String_dict&       dict,
    bool&              is_group,
    mi::neuraylib::Tag old_tag)
{
    bool type_found = false;
    for (mi::Sint32 i = 0; s_supported_types[i] != 0; i += 2)
    {
        if (s_supported_types[i] == elem_type)
        {
            type_found = true;
        }
    }
    if (!type_found)
    {
        WARN_LOG << "The scene element type '" << elem_type << "' is not in the list of "
                 << "supported types.";
    }


    // Load additional scene element settings from a description file if specified. The current
    // settings from the dictionary take precedence over those from the description file.
    const std::string description_file = dict.get("description_file");
    if (!description_file.empty())
    {
        String_dict extra_dict;
        if (!load_application_project_file(description_file, extra_dict))
        {
            ERROR_LOG << "Could not load description file '" << description_file
                      << "' for element '" << elem_name << "'";
            return mi::neuraylib::Tag();
        }
        dict.insert_new(extra_dict);
    }

    mi::base::Handle<nv::index::IScene_element> element;

    is_group = false;
    m_single_element = old_tag.is_valid();
    bool reuse_tag = false;

    if (!element.is_valid_interface())
    {
        element = import_simple_shapes(elem_name, elem_type, dict);
    }

    if (!element.is_valid_interface())
    {
        element = import_raster_shapes(elem_name, elem_type, dict);
    }

    if (!element.is_valid_interface())
    {
        element = import_massive_shapes(elem_name, elem_type, dict, old_tag, m_dice_transaction, reuse_tag);
    }

    if (!element.is_valid_interface())
    {
        element = import_attributes(elem_name, elem_type, dict);
    }

    if (!element.is_valid_interface())
    {
        element = import_raster_attributes(elem_name, elem_type, dict);
    }

    if (!element.is_valid_interface())
    {
        element = import_scene_groups(elem_name, elem_type, dict);
        if (element.is_valid_interface())
        {
            is_group = true;
        }
    }

    mi::neuraylib::Tag new_tag;
    if (element.is_valid_interface())
    {
        // Common settings
        element->set_enabled(get_bool(dict.get("enabled", "true")));

        mi::base::Handle<nv::index::IShape> shape(
            element->get_interface<nv::index::IShape>());
        if (shape.is_valid_interface())
        {
            shape->set_pickable(get_bool(dict.get("pickable", "true")));
        }

        mi::base::Handle<nv::index::IDistributed_data> distributed_data(
            element->get_interface<nv::index::IDistributed_data>());
        if (distributed_data.is_valid_interface())
        {
            distributed_data->set_pickable(get_bool(dict.get("pickable", "true")));
        }

        // Disable volume rendering initially depending on global state
        //TODO: This is highly deprecated!
        if (elem_type == "volume" && !Nvindex_AppData::instance()->is_render_volume())
            element->set_enabled(false);

        if (reuse_tag)
        {
            // Reusing tag, edited existing element, no need to store again
            new_tag = old_tag;
        }
        else
        {
            mi::Sint32 privacy_level = mi::neuraylib::IDice_transaction::LOCAL_SCOPE;

            // Keep the privacy level of an existing scene element, it would be overwritten
            // otherwise. This effectively emulates the behavior of a normal edit operation, which
            // keeps the privacy level unchanged.
            if (old_tag.is_valid())
            {
                privacy_level = m_dice_transaction->get_privacy_level(old_tag);
            }

            // Store the newly created scene element
            new_tag = m_dice_transaction->store_for_reference_counting(
                element.get(), old_tag, elem_name.c_str(), privacy_level);
        }
    }
    else
    {
        // Colormaps have to be handled differently
        new_tag = import_colormaps(elem_name, elem_type, dict);
    }

    return new_tag;
}

void Scene_element_importer::get_supported_types(
    std::multimap<std::string, std::string>& types)
{
    for (mi::Sint32 i = 0; s_supported_types[i] != 0; i += 2)
    {
        const std::string name = s_supported_types[i];
        const std::string type = s_supported_types[i + 1];

        // Skip deprecated types
        if (name == "horizon" || name == "horizon_geometry")
        {
            continue;
        }

        types.insert(std::make_pair(type, name));
    }
}

mi::base::Handle<nv::index::IScene_element> Scene_element_importer::import_scene_groups(
    const std::string&                   elem_name,
    const std::string&                   elem_type,
    const nv::index_common::String_dict& dict)
{
    mi::base::Handle<nv::index::IScene_group> scene_group;

    if (elem_type == "static_scene_group")
    {
        mi::base::Handle<nv::index::IStatic_scene_group> group(
            m_scene->create_scene_group<nv::index::IStatic_scene_group>());
        scene_group = group;
    }
    else if (elem_type == "transformed_scene_group")
    {
        mi::base::Handle<nv::index::ITransformed_scene_group> group(
            m_scene->create_scene_group<nv::index::ITransformed_scene_group>());

        // Transformation
        mi::math::Matrix<mi::Float32, 4, 4> transform_mat(1.f);
        if (dict.is_defined("scale"))
        {
            mi::math::Vector<mi::Float32, 3> scaling = get_vec_float32_3(dict.get("scale"));
            mi::math::Matrix<mi::Float32, 4, 4> scale_matrix(1.f);
            scale_matrix.xx = scaling.x;
            scale_matrix.yy = scaling.y;
            scale_matrix.zz = scaling.z;
            transform_mat *= scale_matrix;
        }

        if (dict.is_defined("rotate"))
        {
            mi::math::Vector<mi::Float32, 3> rotation = get_vec_float32_3(dict.get("rotate"));
            transform_mat.rotate(rotation / 180.f * float(M_PI));
        }

        if (dict.is_defined("translate"))
            transform_mat.translate(get_vec_float32_3(dict.get("translate")));

        // Specifying the full transformation takes precedence
        if (dict.is_defined("transform"))
            transform_mat = get_mat_float32_4_4(dict.get("transform"));

        group->set_transform(transform_mat);
        scene_group = group;
    }

    if (!scene_group.is_valid_interface())
        return null_result();

    // Add the children tags directly if only the scene group is imported
    if (m_single_element)
    {
        std::istringstream children_list(dict.get("children"));
        while (!children_list.eof())
        {
            mi::neuraylib::Tag tag;
            children_list >> tag.id;
            if (tag.is_valid())
            {
                scene_group->append(tag, m_dice_transaction.get());
            }
        }
    }

    return scene_group;
}

mi::neuraylib::Tag Scene_element_importer::resolve_colormap(const std::string& id) const
{
    // When importing a single scene element then the id is actually the tag
    if (m_single_element)
    {
        mi::neuraylib::Tag tag;
        tag.id = get_uint32(id);
        return tag;
    }

    std::map<std::string, mi::neuraylib::Tag>::const_iterator it = m_name_to_tag.find(id);

    // Reuse existing colormap with the same id
    if (it != m_name_to_tag.end())
    {
        return it->second;
    }
    // Look for new inline colormap data with the given id
    else if (m_scene_dict.get(id + "::type", "") == "colormap")
    {
        String_dict dict;
        string_dict_key_prefix_filter(m_scene_dict, id + "::", dict, true);

        mi::Uint32 size = get_uint32(dict.get("size", "0"));

        // Read color values
        std::istringstream data(dict.get("data", ""));
        std::vector<mi::math::Color_struct> cm(size);
        mi::Uint32 i;
        for (i = 0; i < size; ++i)
        {
            if (data.eof() || !data)
                break;
            mi::math::Color_struct& c = cm[i];
            data >> c.r >> c.g >> c.b >> c.a;
        }

        if (i < size || !data)
        {
            ERROR_LOG << "Reading colormap data failed after " << i << " colors for '" << id << "'";
            return mi::neuraylib::Tag();
        }

        mi::base::Handle<nv::index::IColormap> colormap(
            m_scene->create_attribute<nv::index::IColormap>());
        colormap->set_colormap(cm.data(), cm.size());

        mi::math::Vector<mi::Float32, 2> domain = get_vec_float32_2(dict.get("domain", "0 1"));
        colormap->set_domain(domain.x, domain.y);

        const std::string mode_str = dict.get("domain_boundary_mode", "clamp_to_edge");
        nv::index::IColormap::Domain_boundary_mode mode = nv::index::IColormap::CLAMP_TO_EDGE;
        if (mode_str == "clamp_to_transparent")
        {
            mode = nv::index::IColormap::CLAMP_TO_TRANSPARENT;
        }
        else if (mode_str != "clamp_to_edge")
        {
            ERROR_LOG << "Invalid domain_boundary_mode '" << mode_str << "' for colormap '" << id << "'";
        }
        colormap->set_domain_boundary_mode(mode);

        const mi::Sint32 privacy_level = mi::neuraylib::IDice_transaction::LOCAL_SCOPE;
        mi::neuraylib::Tag tag = m_dice_transaction->store_for_reference_counting(
            colormap.get(), mi::neuraylib::NULL_TAG, id.c_str(), privacy_level);
        return tag;
    }

    // Nothing found, the id is probably an application-level colormap index
    return get_colormap_tag(get_uint32(id));
}

mi::Float32 Scene_element_importer::random(mi::Float32 min, mi::Float32 max) const
{
    mi::Float32 t = rand() / (1.f + RAND_MAX);

    return (min + (max - min) * t);
}

// Jet colormap
mi::math::Color_struct Scene_element_importer::jetmap(
    mi::Float32 v,
    mi::Float32 vmin,
    mi::Float32 vmax) const
{
    if (v < vmin)
        v = vmin;
    if (v > vmax)
        v = vmax;

    const mi::Float32 dv = vmax - vmin;
    mi::math::Color c(1.f);

    if (v < (vmin + 0.25f * dv))
    {
        c.r = 0.f;
        c.g = 4.f * (v - vmin) / dv;
    }
    else if (v < (vmin + 0.5f * dv))
    {
        c.r = 0.f;
        c.b = 1.f + 4.f * (vmin + 0.25f * dv - v) / dv;
    }
    else if (v < (vmin + 0.75f * dv))
    {
        c.r = 4.f * (v - vmin - 0.5f * dv) / dv;
        c.b = 0.f;
    }
    else
    {
        c.g = 1.f + 4.f * (vmin + 0.75f * dv - v) / dv;
        c.b = 0.f;
    }

    return c;
}
