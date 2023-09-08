/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "shapes_utility.h"

#include <nv/index/isession.h>
#include <nv/index/iscene.h>
#include <nv/index/iline_set.h>
#include <nv/index/iregular_volume.h>
#include <nv/index/iline_set.h>

#include "common/common_utility.h"
#include "common/forwarding_logger.h"

#include <cassert>


static inline void get_line_segments(
    const mi::math::Bbox< mi::Float32, 3 >&                     bbox,
    std::vector< mi::math::Vector_struct< mi::Float32, 3> >&    line_segment_vertices)
{
    const mi::math::Vector< mi::Float32, 3 > min = bbox.min;
    const mi::math::Vector< mi::Float32, 3 > max = bbox.max;
    mi::math::Vector_struct< mi::Float32, 3> v;

    v.x = min.x; v.y = min.y; v.z = min.z;  // glVertex3f(min.x, min.y, min.z);
    line_segment_vertices.push_back(v);

                 v.y = max.y;               // glVertex3f(min.x, max.y, min.z);
    line_segment_vertices.push_back(v);
    
    v.x = max.x; v.y = min.y;               // glVertex3f(max.x, min.y, min.z);
    line_segment_vertices.push_back(v);
    
                 v.y = max.y;               // glVertex3f(max.x, max.y, min.z);
    line_segment_vertices.push_back(v);
    
                 v.y = min.y; v.z = max.z;  // glVertex3f(max.x, min.y, max.z);
    line_segment_vertices.push_back(v);
    
                 v.y = max.y;               // glVertex3f(max.x, max.y, max.z);
    line_segment_vertices.push_back(v);
    
    v.x = min.x; v.y = min.y;               // glVertex3f(min.x, min.y, max.z);
    line_segment_vertices.push_back(v);
    
                 v.y = max.y;               // glVertex3f(min.x, max.y, max.z);
    line_segment_vertices.push_back(v);

    // Line loop
    v.x = min.x; v.y = min.y; v.z = min.z;  // glVertex3f(min.x, min.y, min.z);
    line_segment_vertices.push_back(v);

    v.x = max.x;                            // glVertex3f(max.x, min.y, min.z);
    line_segment_vertices.push_back(v);
    line_segment_vertices.push_back(v);

    v.x = max.x;              v.z = max.z;  // glVertex3f(max.x, min.y, max.z);
    line_segment_vertices.push_back(v);
    line_segment_vertices.push_back(v);

    v.x = min.x;                            // glVertex3f(min.x, min.y, max.z)
    line_segment_vertices.push_back(v);
    line_segment_vertices.push_back(v);
    
                 v.y = min.y; v.z = min.z;  // glVertex3f(min.x, min.y, min.z); (1st vertex)
    line_segment_vertices.push_back(v);

    // Line loop
    v.x = min.x; v.y = max.y; v.z = min.z;  // glVertex3f(min.x, max.y, min.z);
    line_segment_vertices.push_back(v);

    v.x = max.x;                            // glVertex3f(max.x, max.y, min.z);
    line_segment_vertices.push_back(v);
    line_segment_vertices.push_back(v);

    v.x = max.x;              v.z = max.z;  // glVertex3f(max.x, max.y, max.z);
    line_segment_vertices.push_back(v);
    line_segment_vertices.push_back(v);

    v.x = min.x;                            // glVertex3f(min.x, max.y, max.z)
    line_segment_vertices.push_back(v);
    line_segment_vertices.push_back(v);
    
                              v.z = min.z;  // glVertex3f(min.x, max.y, min.z); (1st vertex)
    line_segment_vertices.push_back(v);
}

//----------------------------------------------------------------------
mi::neuraylib::Tag create_line_set(
    const mi::math::Bbox<mi::Float32, 3>&                       bbox,
    const mi::math::Color&                                      color,
    const mi::base::Handle<mi::neuraylib::IDice_transaction>&   dice_transaction,
    const mi::neuraylib::Tag&                                   session_tag,
    const std::string&                                          lines_name)
{
    if(!dice_transaction.is_valid_interface())
    {
        ERROR_LOG << "Invalid DiCE transaction (file: " << __FILE__ ", line: " << __LINE__ << ").";
        return mi::neuraylib::NULL_TAG;
    }

    if(!session_tag.is_valid())
    {
        ERROR_LOG << "Invalid session tag (file: " << __FILE__ ", line: " << __LINE__ << ").";
        return mi::neuraylib::NULL_TAG;
    }

    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<nv::index::ISession>(session_tag));
    assert(session.is_valid_interface());
    
    mi::base::Handle<const nv::index::IScene> scene(
        dice_transaction->access<const nv::index::IScene>(session->get_scene()));
    assert(scene.is_valid_interface());
    
    std::vector<mi::math::Vector_struct< mi::Float32, 3> > segment_vertices;
    std::vector<mi::math::Color_struct>                    color_per_segment;
    std::vector<mi::Float32>                               width_per_segment;
    get_line_segments(bbox, segment_vertices);
    color_per_segment.push_back(color);
    color_per_segment.push_back(color);
    color_per_segment.push_back(color);
    color_per_segment.push_back(color);
    color_per_segment.push_back(color);
    color_per_segment.push_back(color);
    color_per_segment.push_back(color);
    color_per_segment.push_back(color);
    color_per_segment.push_back(color);
    color_per_segment.push_back(color);
    color_per_segment.push_back(color);
    color_per_segment.push_back(color);
    width_per_segment.push_back(2.f);
    width_per_segment.push_back(2.f);
    width_per_segment.push_back(2.f);
    width_per_segment.push_back(2.f);
    width_per_segment.push_back(2.f);
    width_per_segment.push_back(2.f);
    width_per_segment.push_back(2.f);
    width_per_segment.push_back(2.f);
    width_per_segment.push_back(2.f);
    width_per_segment.push_back(2.f);
    width_per_segment.push_back(2.f);
    width_per_segment.push_back(2.f);

    // Create line set using the scene's factory 
    mi::base::Handle<nv::index::ILine_set> line_set(scene->create_shape<nv::index::ILine_set>());
    // ... set the line style such a dashed or dotted or solid linestyle (see ILine_set)
    line_set->set_line_style(nv::index::ILine_set::LINE_STYLE_DASHED);
    // ... set the line type. Currently only line segments are supported (see ILine_set)
    line_set->set_line_type(nv::index::ILine_set::LINE_TYPE_SEGMENTS);
    // ... and the lines with colors and widths all in line segment order.
    line_set->set_lines(&segment_vertices[0], segment_vertices.size());
    line_set->set_colors(&color_per_segment[0], color_per_segment.size());
    line_set->set_widths(&width_per_segment[0], width_per_segment.size());

    mi::neuraylib::Tag tag = dice_transaction->store_for_reference_counting(
            line_set.get(),
            mi::neuraylib::NULL_TAG,
            lines_name.c_str());
    
    return tag;
}

//----------------------------------------------------------------------
