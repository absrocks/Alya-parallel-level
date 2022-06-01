/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/

#include "scene_element_importer.h"

#include <nv/index/icone.h>
#include <nv/index/icylinder.h>
#include <nv/index/iellipsoid.h>
#include <nv/index/iplane.h>
#include <nv/index/isphere.h>

using namespace nv::index_common;

mi::base::Handle<nv::index::IScene_element> Scene_element_importer::import_simple_shapes(
    const std::string& elem_name,
    const std::string& elem_type,
    const String_dict& dict)
{
    //
    // Sphere
    //
    if (elem_type == "sphere")
    {
        mi::base::Handle<nv::index::ISphere> sphere(
            m_scene->create_shape<nv::index::ISphere>());
        if (dict.is_defined("center"))
            sphere->set_center(get_vec_float32_3(dict.get("center")));
        if (dict.is_defined("radius"))
            sphere->set_radius(get_float32(dict.get("radius")));

        return sphere;
    }
    //
    // Ellipsoid
    //
    else if (elem_type == "ellipsoid")
    {
        mi::base::Handle<nv::index::IEllipsoid> ellipsoid(
            m_scene->create_shape<nv::index::IEllipsoid>());
        if (dict.is_defined("center"))
            ellipsoid->set_center(get_vec_float32_3(dict.get("center")));

        mi::math::Vector<mi::Float32, 3> a(1.f, 0.f, 0.f);
        mi::math::Vector<mi::Float32, 3> b(0.f, 1.f, 0.f);
        mi::Float32 c = 1.f;

        if (dict.is_defined("semi_axis_a"))
            a = get_vec_float32_3(dict.get("semi_axis_a"));
        if (dict.is_defined("semi_axis_b"))
            b = get_vec_float32_3(dict.get("semi_axis_b"));
        if (dict.is_defined("semi_axis_c"))
            c = get_float32(dict.get("semi_axis_c"));

        ellipsoid->set_semi_axes(a, b, c);

        return ellipsoid;
    }
    //
    // Cylinder
    //
    else if (elem_type == "cylinder")
    {
        mi::base::Handle<nv::index::ICylinder> cylinder(
            m_scene->create_shape<nv::index::ICylinder>());

        if (dict.is_defined("bottom"))
            cylinder->set_bottom(get_vec_float32_3(dict.get("bottom")));
        if (dict.is_defined("top"))
            cylinder->set_top(get_vec_float32_3(dict.get("top")));
        if (dict.is_defined("radius"))
            cylinder->set_radius(get_float32(dict.get("radius")));
        if (dict.is_defined("capped"))
            cylinder->set_capped(get_bool(dict.get("capped")));

        return cylinder;
    }
    //
    // Cone
    //
    else if (elem_type == "cone")
    {
        mi::base::Handle<nv::index::ICone> cone(
            m_scene->create_shape<nv::index::ICone>());

        if (dict.is_defined("center"))
            cone->set_center(get_vec_float32_3(dict.get("center")));
        if (dict.is_defined("tip"))
            cone->set_tip(get_vec_float32_3(dict.get("tip")));
        if (dict.is_defined("radius"))
            cone->set_radius(get_float32(dict.get("radius")));
        if (dict.is_defined("capped"))
            cone->set_capped(get_bool(dict.get("capped")));

        return cone;
    }
    //
    // Plane
    //
    else if (elem_type == "plane" || elem_type == "quad" || elem_type == "rectangle")
    {
        mi::base::Handle<nv::index::IPlane> plane(
            m_scene->create_shape<nv::index::IPlane>());

        plane->set_point(get_vec_float32_3(dict.get("position", "0 0 0")));
        plane->set_normal(get_vec_float32_3(dict.get("normal", "0 0 0")));
        plane->set_up(get_vec_float32_3(dict.get("up", "0 0 0")));

        mi::math::Vector<mi::Float32, 2> extent = get_vec_float32_2(dict.get("extent", "0 0"));
        plane->set_extent(extent);

        // Optionally move the plane along the normal (used in the scene editor)
        if (dict.is_defined("roam"))
        {
            const mi::Float32 amount = get_float32(dict.get("roam"));
            const mi::math::Vector<mi::Float32, 3> p = plane->get_point();
            const mi::math::Vector<mi::Float32, 3> n = plane->get_normal();
            plane->set_point(p + amount * n);
        }

        return plane;
    }

    return null_result();
}
