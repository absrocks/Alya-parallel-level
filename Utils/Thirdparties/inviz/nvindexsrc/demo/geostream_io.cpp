/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "geostream_io.h"

#include <sstream>

#include <nv/index/iintersection_highlighting.h>
#include <nv/index/iregular_volume_texture.h>
#include <nv/index/iscene.h>

#include "common/scene_logger.h"
#include "common/common_utility.h"
#include "common/forwarding_logger.h"

#include "scene_element_importer.h"
#include "scene_utility.h"
#include "utilities.h"
#include "nvindex_appdata.h"

using namespace nv::index_common;

namespace {

void add_scene_graph_subtree(
    Scene_element_importer&                                   importer,
    const mi::base::Handle<mi::neuraylib::IDice_transaction>& dice_transaction,
    const nv::index::IScene*                                  scene,
    mi::neuraylib::Tag                                        scene_group_tag,
    const std::string&                                        scene_group_name,
    std::map<std::string, mi::neuraylib::Tag>&                name_to_tag,
    std::vector<std::string>&                                 referencing_elements,
    const String_dict&                                        scene_dict)
{
    std::istringstream children_list(scene_dict.get(scene_group_name + "::children"));
    while (!children_list.eof())
    {
        std::string elem_name;
        children_list >> elem_name;
        if (elem_name.empty())
            continue;

        // Reuse already created elements
        std::map<std::string, mi::neuraylib::Tag>::const_iterator it;
        it = name_to_tag.find(elem_name);
        if (it != name_to_tag.end())
        {
            // Just append the existing tag
            //TODO: Do not allow multiple for heightfield and volumes (check in library as well?)
            mi::base::Handle<nv::index::IScene_group> scene_group(
                dice_transaction->edit<nv::index::IScene_group>(scene_group_tag));

            scene_group->append(it->second, dice_transaction.get());
            continue;
        }

        String_dict dict;
        string_dict_key_prefix_filter(scene_dict, elem_name + "::", dict, true);

        // Determine the type of the element
        std::string elem_type = dict.get("type");
        if (elem_type.empty())
        {
            ERROR_LOG << "Element not found or type not specified for scene graph element '"
                      << elem_name << "'";
            continue;
        }

        // Try to import the element
        bool is_group;
        mi::neuraylib::Tag new_tag = importer.import(elem_name, elem_type, dict, is_group);

        if (new_tag.is_valid())
        {
            if (is_group)
            {
                // It was a group, traverse it
                add_scene_graph_subtree(
                    importer, dice_transaction, scene, new_tag, elem_name,
                    name_to_tag, referencing_elements, scene_dict);
            }

            // Add the new element to its parent scene group
            mi::base::Handle<nv::index::IScene_group> scene_group(
                dice_transaction->edit<nv::index::IScene_group>(scene_group_tag));
            scene_group->append(new_tag, dice_transaction.get());
            name_to_tag[elem_name] = new_tag;
        }
        else
        {
            ERROR_LOG << "Could not create scene graph element type '" << elem_type
                      << "' for element '" << elem_name << "'";
        }
    }

    // Resolve references at the end
    if (scene_group_name == "root")
    {
        for (size_t i=0; i < referencing_elements.size(); ++i)
        {
            const std::string elem_name = referencing_elements[i];

            String_dict dict;
            string_dict_key_prefix_filter(scene_dict, elem_name + "::", dict, true);
            const std::string elem_type = dict.get("type");
            const mi::neuraylib::Tag elem_tag = name_to_tag[elem_name];

            // Reference to intersecting shape in intersection highlighting attribute
            if (elem_type == "intersection_highlighting")
            {
                mi::base::Handle<nv::index::IIntersection_highlighting> attr(
                    dice_transaction->edit<nv::index::IIntersection_highlighting>(elem_tag));
                if (attr)
                {
                    std::string with = dict.get("intersect_with");
                    std::map<std::string, mi::neuraylib::Tag>::const_iterator it = name_to_tag.find(with);
                    if (it == name_to_tag.end())
                    {
                        ERROR_LOG << "Could not find scene graph element '" << with << "' which is "
                                  << "referenced by '" << elem_name << "'";
                    }
                    else
                    {
                        attr->set_intersection_shape(it->second);
                    }
                }
            }
            else if (elem_type == "volume_texture")
            {
                mi::base::Handle<nv::index::IRegular_volume_texture> attr(
                    dice_transaction->edit<nv::index::IRegular_volume_texture>(elem_tag));
                if (attr)
                {
                    std::string rvol_elem = dict.get("volume");
                    std::map<std::string, mi::neuraylib::Tag>::const_iterator it = name_to_tag.find(rvol_elem);
                    if (it == name_to_tag.end())
                    {
                        ERROR_LOG << "Could not find scene graph element '" << rvol_elem << "' which is "
                                  << "referenced by '" << elem_name << "'";
                    }
                    else
                    {
                        attr->set_volume_element(it->second);
                    }
                }
            }
        }
    }
}

} // namespace

mi::Sint32 add_scene_graph(
    const mi::base::Handle<mi::neuraylib::IDice_transaction>& dice_transaction,
    mi::neuraylib::Tag                                        session_tag,
    const String_dict&                                        dict)
{
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<nv::index::ISession>(session_tag));
    mi::base::Handle<const nv::index::IScene> scene(
        dice_transaction->access<nv::index::IScene>(session->get_scene()));

    String_dict scene_dict;
    string_dict_key_prefix_filter(dict, "app::scene::", scene_dict, true);

    std::map<std::string, mi::neuraylib::Tag> name_to_tag; // Maps scene graph element names to their tags
    std::vector<std::string> referencing_elements;

    Scene_element_importer importer(
        dice_transaction, scene.get(), session_tag, name_to_tag, referencing_elements, scene_dict);

    // Start traversal at the root node, i.e. the scene
    add_scene_graph_subtree(
        importer,
        dice_transaction,
        scene.get(),
        session->get_scene(),
        "root",
        name_to_tag,
        referencing_elements,
        scene_dict);

    // Just for verification:

    /*    const bool log_scene = get_bool(dict.get("app::scene::log", "no"));
    if (log_scene)
    {
        nv::index_common::Scene_logger* logger = new nv::index_common::Scene_logger;
        session->visit(static_cast<nv::index::IScene_visitor*>(logger), dice_transaction.get());
        delete logger;
	}*/
    const bool log_scene = get_bool(dict.get("app::scene::log", "yes"));

    if (log_scene)

      {

	nv::index_common::Scene_logger* logger = new nv::index_common::Scene_logger;

	session->visit(static_cast<nv::index::IScene_visitor*>(logger), dice_transaction.get());



	const mi::Uint32 nb_color_maps = Nvindex_AppData::instance()->get_colormap_manager()->get_number_of_colormap();

	const std::vector<mi::neuraylib::Tag>& colormaps = logger->get_color_maps();

	for (mi::Uint32 i = 0; i < colormaps.size(); ++i)

	  {

	    Nvindex_AppData::instance()->get_colormap_manager()->set_tag(nb_color_maps+i, colormaps[i]);

	  }



	delete logger;

      }
    return name_to_tag.size();
}
