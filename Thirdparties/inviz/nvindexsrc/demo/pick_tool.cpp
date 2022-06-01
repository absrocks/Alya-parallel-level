/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Picking tool implementation


#include "pick_tool.h"

#include "common/string_dict.h"

#include <nv/index/iheightfield_pick_result.h>
#include <nv/index/iintersection_highlight_pick_result.h>
#include <nv/index/ilight.h>
#include <nv/index/iline_set_pick_result.h>
#include <nv/index/imaterial.h>
#include <nv/index/ipath_query_results.h>
#include <nv/index/iplane.h>
#include <nv/index/iplane_pick_result.h>
#include <nv/index/ipoint_set.h>
#include <nv/index/ipoint_set_pick_result.h>
#include <nv/index/isparse_volume_scene_element.h>
#include <nv/index/itriangle_mesh_query_results.h>

#include "nvindex_appdata.h"
#include "nvindex_rendering_context.h"

#include <cassert>


//======================================================================
//----------------------------------------------------------------------
Pick_tool::Pick_tool()
{
    // empty
}
//----------------------------------------------------------------------
Pick_tool::~Pick_tool()
{
    // empty
}

//----------------------------------------------------------------------
mi::math::Vector<mi::Float32, 3> Pick_tool::get_pick_point_xyz(
    nv::index::IScene_pick_result* pick_result)
{
    assert(pick_result != 0);

    const mi::math::Vector<mi::Float32, 3> intersection_ijk = pick_result->get_intersection();
    const mi::math::Matrix<mi::Float32, 4, 4> transform     = pick_result->get_transform();
    const mi::math::Vector<mi::Float32, 3> intersection_xyz = 
        mi::math::transform_point(transform, intersection_ijk);

    return intersection_xyz;
}

//----------------------------------------------------------------------
std::string Pick_tool::get_pick_scene_element_tag_prefix()
{
    std::string pick_tag_prefix = "";
    const bool is_show_pick_points_in_scene_editor =
        nv::index_common::get_bool(Nvindex_AppData::instance()->peek_app_proj()->
                                   get("app::experimental::show_pick_point_in_scene_editor", "no"));
    if (is_show_pick_points_in_scene_editor)
    {
        pick_tag_prefix = "internal";
    }

    return pick_tag_prefix;
}

//----------------------------------------------------------------------
void Pick_tool::create_pick_point(
    nv::index::IScene_pick_result*    pick_result,
    const mi::neuraylib::Tag&         session_tag,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    assert(pick_result != 0);
    assert(dice_transaction != 0);

    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<nv::index::ISession>(session_tag));
    mi::base::Handle<nv::index::IScene> scene(
        dice_transaction->edit<nv::index::IScene>(session->get_scene()));

    // Create the containing scene group
    mi::base::Handle<nv::index::ITransformed_scene_group> scene_group(
        scene->create_scene_group<nv::index::ITransformed_scene_group>());

    // Create material attribute
    mi::base::Handle<nv::index::IPhong_gl> mat(
        scene->create_attribute<nv::index::IPhong_gl>());

    // Use white as base color, will be multiplied by the pick color
    mat->set_ambient(mi::math::Color(0.5f));
    mat->set_diffuse(mi::math::Color(0.4f));
    mat->set_specular(mi::math::Color(0.4f));
    mat->set_shininess(100.f);
    mat->set_opacity(1.f);

    // Create scene elements in global scope, so that no invalid tag access will happen when the
    // scene root is not localized
    const mi::neuraylib::Privacy_level GLOBAL_SCOPE = 0;

    const std::string pick_tag_prefix = get_pick_scene_element_tag_prefix();
    mi::neuraylib::Tag phong_tag = dice_transaction->store_for_reference_counting(
        mat.get(), mi::neuraylib::Tag(),
        (pick_tag_prefix + "__nv_viewer_pick_material").c_str(), GLOBAL_SCOPE);
    scene_group->append(phong_tag, dice_transaction);

    // Create light attribute
    mi::base::Handle<nv::index::IDirectional_headlight> light(
        scene->create_attribute<nv::index::IDirectional_headlight>());

    light->set_direction(mi::math::Vector<mi::Float32, 3>(1.f, -1.f, -1.f));
    light->set_intensity(mi::math::Color(1.f));

    mi::neuraylib::Tag light_tag = dice_transaction->store_for_reference_counting(
        light.get(), mi::neuraylib::Tag(),
        (pick_tag_prefix + "__nv_viewer_pick_light").c_str(), GLOBAL_SCOPE);
    scene_group->append(light_tag, dice_transaction);

    // Create the point set that represents the pick point
    std::vector<mi::math::Vector_struct<mi::Float32, 3> > vertices;
    const mi::math::Vector<mi::Float32, 3> intersection_xyz = get_pick_point_xyz(pick_result);
    vertices.push_back(intersection_xyz);
    vertices.push_back(intersection_xyz);

    INFO_LOG << "########## Adding pick point at " << intersection_xyz;

    std::vector<mi::math::Color_struct> per_vertex_color_values;
    mi::math::Color color = pick_result->get_color();
    per_vertex_color_values.push_back(color);
    color.a = 1.f;
    per_vertex_color_values.push_back(color); // second point is fully opaque

    std::vector<mi::Float32> per_vertex_radii;
    per_vertex_radii.push_back(30.f);
    per_vertex_radii.push_back(5.f);

    mi::base::Handle<nv::index::IPoint_set> pick_point(scene->create_shape<nv::index::IPoint_set>());
    // ... set the point styles
    pick_point->set_point_style(nv::index::IPoint_set::SHADED_CIRCLE);
    // ... and the vertices with colors and radii.
    pick_point->set_vertices(&vertices[0], vertices.size());
    pick_point->set_colors(&per_vertex_color_values[0], per_vertex_color_values.size());
    pick_point->set_radii(&per_vertex_radii[0], per_vertex_radii.size());

    // Do not allow picking the pick point
    pick_point->set_pickable(false);

    mi::neuraylib::Tag pick_point_tag = dice_transaction->store_for_reference_counting(
        pick_point.get(), mi::neuraylib::Tag(),
        (pick_tag_prefix + "__nv_viewer_pick_point").c_str(), GLOBAL_SCOPE);

    scene_group->append(pick_point_tag, dice_transaction);

    // Store the scene group and append it to the scene
    mi::neuraylib::Tag scene_group_tag = dice_transaction->store_for_reference_counting(
        scene_group.get(), mi::neuraylib::Tag(),
        (pick_tag_prefix + "__nv_viewer_pick_group").c_str(), GLOBAL_SCOPE);

    scene->append(scene_group_tag, dice_transaction);
}

//----------------------------------------------------------------------
void Pick_tool::remove_pick_point(
    const mi::neuraylib::Tag&         session_tag,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    assert(dice_transaction != 0);

    const std::string pick_tag_prefix = get_pick_scene_element_tag_prefix();

    const mi::neuraylib::Tag pick_point_group =
        dice_transaction->name_to_tag((pick_tag_prefix + "__nv_viewer_pick_group").c_str());
    if (pick_point_group.is_valid())
    {
        mi::base::Handle<const nv::index::ISession> session(
            dice_transaction->access<nv::index::ISession>(session_tag));
        mi::base::Handle<nv::index::IScene> scene(
            dice_transaction->edit<nv::index::IScene>(session->get_scene()));

        scene->remove(pick_point_group, dice_transaction);
    }
}

//----------------------------------------------------------------------
void Pick_tool::print_pick_result(
    Nvindex_rendering_context*        irc_ref,
    mi::Uint32                        intersect_num,
    nv::index::IScene_pick_result*    pick_result,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    assert(pick_result != 0);
    assert(dice_transaction != 0);

    std::stringstream log_out;
    log_out << "Intersection no. "    << intersect_num << "\n"
            << "\t Element (tag)    " << pick_result->get_scene_element().id << "\n";

    log_out << "\t Scene path tags: ";
    mi::base::Handle<const nv::index::IScene_path> scene_path(pick_result->get_scene_path());
    for (mi::Uint32 i = 0; i < scene_path->get_length(); ++i)
    {
        if (i > 0)
        {
            log_out << ", ";
        }
        log_out << scene_path->get_node(i);
    }
    log_out << "\n";

    mi::math::Vector<mi::Float32, 3> intersection_ijk = pick_result->get_intersection();
    mi::math::Vector<mi::Float32, 3> intersection_xyz = get_pick_point_xyz(pick_result);

    log_out << "\t Sub index:       " << pick_result->get_scene_element_sub_index() << "\n"
            << "\t Distance:        " << pick_result->get_distance() << "\n"
            << "\t IJK position:    " << intersection_ijk << "\n"
            << "\t XYZ position:    " << intersection_xyz << "\n"
            << "\t Color:           " << pick_result->get_color() << "\n";

    const mi::base::Handle<const nv::index::IData_sample> data_sample(pick_result->get_data_sample());
    if (data_sample != 0)
    {
        const mi::base::Handle<const nv::index::IData_sample_uint8>   ds_uint8  (data_sample->get_interface<const nv::index::IData_sample_uint8>());
        const mi::base::Handle<const nv::index::IData_sample_uint16>  ds_uint16 (data_sample->get_interface<const nv::index::IData_sample_uint16>());
        const mi::base::Handle<const nv::index::IData_sample_float32> ds_float32(data_sample->get_interface<const nv::index::IData_sample_float32>());
        const mi::base::Handle<const nv::index::IData_sample_rgba8>   ds_rgba8  (data_sample->get_interface<const nv::index::IData_sample_rgba8>());
        if (ds_uint8.is_valid_interface())
        {
            log_out << "\t Data sample:     " << mi::Uint32(ds_uint8->get_sample_value()) << " (uint8)\n";
        }
        else if (ds_uint16.is_valid_interface())
        {
            log_out << "\t Data sample:     " << ds_uint16->get_sample_value() << " (uint16)\n";
        }
        else if (ds_float32.is_valid_interface())
        {
            log_out << "\t Data sample:     " << ds_float32->get_sample_value() << " (float32)\n";
        }
        else if (ds_rgba8.is_valid_interface())
        {
            log_out << "\t Data sampler:     " << (mi::Uint32)ds_rgba8->get_sample_value().x << ", "
                    << (mi::Uint32)ds_rgba8->get_sample_value().y << ", "
                    << (mi::Uint32)ds_rgba8->get_sample_value().z << ", "
                    << (mi::Uint32)ds_rgba8->get_sample_value().w << " (rgba8)\n";
        }
    }
    INFO_LOG << log_out.str();

    if (pick_result->get_intersection_info_class() == nv::index::IPoint_set_pick_result::IID())
    {
        mi::base::Handle<const nv::index::IPoint_set_pick_result> point_pick_result(
            pick_result->get_interface<const nv::index::IPoint_set_pick_result>());
        if(point_pick_result.is_valid_interface())
        {
            INFO_LOG << "Point set specific pick results:";
            INFO_LOG << "\t Point index:     " << point_pick_result->get_point();
        }
    }

    if (pick_result->get_intersection_info_class() == nv::index::ILine_set_pick_result::IID())
    {
        mi::base::Handle<const nv::index::ILine_set_pick_result> line_pick_result(
            pick_result->get_interface<const nv::index::ILine_set_pick_result>());
        if (line_pick_result.is_valid_interface())
        {
            INFO_LOG << "Line set specific pick results:";
            INFO_LOG << "\t Segment index:   " << line_pick_result->get_segment();
        }
    }

    if (pick_result->get_intersection_info_class() == nv::index::IPath_pick_result::IID())
    {
        mi::base::Handle<const nv::index::IPath_pick_result> path_pick_result(
            pick_result->get_interface<const nv::index::IPath_pick_result>());
        if (path_pick_result.is_valid_interface())
        {
            INFO_LOG << "Path specific pick results:";
            INFO_LOG << "\t Segment index:   " << path_pick_result->get_segment();
        }
    }

    if (pick_result->get_intersection_info_class() == nv::index::ITriangle_mesh_pick_result::IID())
    {
        mi::Uint64 triangle_index = ~0ull;
        mi::base::Handle<const nv::index::ITriangle_mesh_pick_result> tri_mesh_pick_result(
            pick_result->get_interface<const nv::index::ITriangle_mesh_pick_result>());
        if (tri_mesh_pick_result.is_valid_interface())
        {
            triangle_index = tri_mesh_pick_result->get_triangle_index();
            INFO_LOG << "Triangle mesh specific pick results:" << "\n"
                     << "\t\t Triangle id:            " << tri_mesh_pick_result->get_triangle_index() << "\n"
                     << "\t\t Local triangle id:      " << tri_mesh_pick_result->get_local_triangle_index() << "\n"
                     << "\t\t Normal:                 " << tri_mesh_pick_result->get_normal() << "\n"
                     << "\t\t Color:                  " << tri_mesh_pick_result->get_color_value() << "\n"
                     << "\t\t Colormap color:         " << tri_mesh_pick_result->get_colormap_value() << "\n"
                     << "\t\t Barycentric coordinate: " << tri_mesh_pick_result->get_barycentric_coordinates();
        }
        // Valid triangle intersection, now query the triangle's details.
        if (triangle_index != ~0ull)
        {
            INFO_LOG << "Query the triangle details for triangle id: " << triangle_index << ".";
            // Query trinagle details
            mi::base::Handle<nv::index::IScene_lookup_result> entry_lookup_result(
                irc_ref->m_iindex_scene_query->entry_lookup(
                    triangle_index,                     // The global triangle index used to query the triangle's details.
                    pick_result->get_scene_element(),        // The tag that represents the triangle scene element.
                    irc_ref->m_session_tag, dice_transaction));

            if (entry_lookup_result.is_valid_interface())
            {
                mi::base::Handle<const nv::index::ITriangle_mesh_lookup_result> tri_mesh_lookup_result(
                    entry_lookup_result->get_interface<const nv::index::ITriangle_mesh_lookup_result>());
                if (tri_mesh_lookup_result.is_valid_interface())
                {
                    INFO_LOG << "Receiving the following per-vertex attributes:";
                    mi::math::Vector_struct<mi::Float32, 3> v0;
                    mi::math::Vector_struct<mi::Float32, 3> v1;
                    mi::math::Vector_struct<mi::Float32, 3> v2;
                    tri_mesh_lookup_result->get_vertices(v0, v1, v2);
                    INFO_LOG << "\tVertices:            (" << v0 << ","<< v1 << "," << v2 << ")";

                    tri_mesh_lookup_result->get_normals(v0, v1, v2);
                    INFO_LOG << "\tNormals:             (" << v0 << ","<< v1 << "," << v2 << ")";

                    mi::math::Vector_struct<mi::Float32, 2> st0;
                    mi::math::Vector_struct<mi::Float32, 2> st1;
                    mi::math::Vector_struct<mi::Float32, 2> st2;
                    tri_mesh_lookup_result->get_texture_coordinates(st0, st1, st2);
                    INFO_LOG << "\tTexture coordinates: (" << st0 << ","<< st1 << "," << st2 << ")";

                    mi::math::Color_struct c0;
                    mi::math::Color_struct c1;
                    mi::math::Color_struct c2;
                    tri_mesh_lookup_result->get_colors(c0, c1, c2);
                    INFO_LOG << "\tColors:              (" << c0 << ","<< c1 << "," << c2 << ")";

                    mi::Uint32 cm0, cm1, cm2;
                    tri_mesh_lookup_result->get_color_indices(cm0, cm1, cm2);
                    INFO_LOG << "\tColor index values:  (" << cm0 << ","<< cm1 << "," << cm2 << ")";
                }
            }
        }
    }

    if (pick_result->get_intersection_info_class() == nv::index::IHeightfield_pick_result::IID())
    {
        mi::base::Handle<const nv::index::IHeightfield_pick_result> compute_pick_result(
            pick_result->get_interface<const nv::index::IHeightfield_pick_result>());
        if (compute_pick_result.is_valid_interface() && compute_pick_result->is_computing_enabled())
        {
            INFO_LOG << "Compute heightfield specific pick results:";
            INFO_LOG << "\t Computed color values:   " << compute_pick_result->get_computed_color();
        }
    }

    if (pick_result->get_intersection_info_class() == nv::index::IPlane_pick_result::IID())
    {
        mi::base::Handle<const nv::index::IPlane_pick_result> compute_pick_result(
            pick_result->get_interface<const nv::index::IPlane_pick_result>());
        if (compute_pick_result.is_valid_interface())
        {
            INFO_LOG << "Plane specific pick results:";
            INFO_LOG << "\t Texture color values:    " << compute_pick_result->get_texture_color();
        }
    }

    if (pick_result->get_intersection_info_class() == nv::index::IIntersection_highlight_pick_result::IID())
    {
        mi::base::Handle<const nv::index::IIntersection_highlight_pick_result> highlight_pick_result(
            pick_result->get_interface<const nv::index::IIntersection_highlight_pick_result>());
        if (highlight_pick_result.is_valid_interface())
        {
            INFO_LOG << "Intersection highlight specific pick results:";
            INFO_LOG << "\t Intersection shape:      " << highlight_pick_result->get_intersection_shape();
        }
    }
}

//----------------------------------------------------------------------
mi::neuraylib::Tag Pick_tool::issue_pick_command(
    Nvindex_rendering_context*             irc_ref,
    bool                                   is_multi_view_mode,
    const mi::math::Vector<mi::Sint32, 2>& pick_location_on_canvas,
    Examiner_manipulator*                  examiner)
{
    assert(irc_ref  != 0);
    assert(examiner != 0);

    // Lock the scene while issuing the pick command
    mi::base::Lock::Block block(&Nvindex_AppData::instance()->m_scene_edit_lock);

    const nv::index::IIndex_canvas* pick_canvas = irc_ref->m_span_buffer.get();
    mi::base::Handle<nv::index::IScene_pick_results> scene_pick_results;

    if (is_multi_view_mode)
    {
        // multiview pick
        mi::base::Handle<nv::index::IScene_pick_results_list> scene_pick_results_list(
            irc_ref->m_iindex_scene_query->pick(
                pick_location_on_canvas,
                pick_canvas,
                irc_ref->m_viewport_list.get(),
                irc_ref->m_session_tag));

        for (mi::Size i = 0; i < scene_pick_results_list->size(); ++i)
        {
            mi::base::Handle<nv::index::IScene_pick_results> tmp_scene_pick_results(
                scene_pick_results_list->get(i));
            INFO_LOG << "Picked viewport " << i << ": " << tmp_scene_pick_results->get_viewport_index();
            if (i == 0)
            {
                scene_pick_results = scene_pick_results_list->get(i); // keep only the first one
            }
        }
    }
    else
    {
        // single view pick
        mi::base::Handle<mi::neuraylib::IScope> cur_scope(irc_ref->get_current_scope());
        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            cur_scope->create_transaction<mi::neuraylib::IDice_transaction>());
        assert(dice_transaction.is_valid_interface());
        {
            const mi::math::Vector<mi::Uint32, 2> pick_location_on_canvas_uint(
                static_cast<mi::Uint32>(pick_location_on_canvas.x),
                static_cast<mi::Uint32>(pick_location_on_canvas.y));

            scene_pick_results =
                irc_ref->m_iindex_scene_query->pick(
                    pick_location_on_canvas_uint, pick_canvas, irc_ref->m_session_tag, dice_transaction.get());
            // Remove old pick point scene group
            remove_pick_point(irc_ref->m_session_tag, dice_transaction.get());
        }
        dice_transaction->commit();
    }

    // Print the pick result
    if (!scene_pick_results.is_valid_interface() || scene_pick_results->get_nb_results() == 0)
    {
        INFO_LOG << "No intersection found.";
        examiner->set_pan_reference_point_valid(false);
        return mi::neuraylib::Tag();
    }

    mi::neuraylib::Tag selected_tag;
    const mi::Uint32 nb_results = scene_pick_results->get_nb_results();
    INFO_LOG << "Pick location on canvas: " << pick_location_on_canvas;
    INFO_LOG << "Number of intersections: " << nb_results;

    {
        mi::base::Handle<mi::neuraylib::IScope> cur_scope(irc_ref->get_current_scope());
        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            cur_scope->create_transaction<mi::neuraylib::IDice_transaction>());
        assert(dice_transaction.is_valid_interface());
        {
            for(mi::Uint32 i = 0; i < nb_results; ++i)
            {
                const mi::base::Handle<nv::index::IScene_pick_result> result(scene_pick_results->get_result(i));
                print_pick_result(irc_ref, i, result.get(), dice_transaction.get());
                if (i == 0)
                {
                    selected_tag = result->get_scene_element();
                }
            }

            // Render first pick point
            const bool show_pick = nv::index_common::get_bool(Nvindex_AppData::instance()->peek_app_proj()->get("app::show_pick_point", "true"));
            if (show_pick)
            {
                const mi::base::Handle<nv::index::IScene_pick_result> pick_result(
                    scene_pick_results->get_result(0));
                assert(pick_result.is_valid_interface());

                create_pick_point(pick_result.get(), irc_ref->m_session_tag, dice_transaction.get());
            }

            // Set rotation center to the examiner
            {
                // get the first pick point
                mi::base::Handle<nv::index::IScene_pick_result> pick_result(scene_pick_results->get_result(0));
                assert(pick_result.is_valid_interface());

                mi::math::Vector<mi::Float32, 3> pick_xyz_position = get_pick_point_xyz(pick_result.get());

                mi::base::Handle<const nv::index::ISession> session(
                    dice_transaction->access<const nv::index::ISession>(irc_ref->m_session_tag));
                assert(session.is_valid_interface());
                mi::base::Handle<const nv::index::IScene> scene(
                    dice_transaction->access<const nv::index::IScene>(session->get_scene()));
                assert(scene.is_valid_interface());
                const mi::math::Matrix<mi::Float32, 4, 4> scene_transform = scene->get_transform_matrix();
                const mi::math::Vector<mi::Float32, 3> pick_w_position =
                    mi::math::transform_point(scene_transform, pick_xyz_position);

                // for rotation
                examiner->set_examiner_rotation_center(pick_w_position);
                // for pan
                examiner->set_pan_reference_point_valid(true);
                examiner->set_pan_reference_point_object_space(pick_xyz_position);
            }
        }
        dice_transaction->commit();
    }

    return selected_tag;
}

//----------------------------------------------------------------------
