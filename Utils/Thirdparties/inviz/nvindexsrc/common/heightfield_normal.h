/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief heightfield normal generator example implementation
#ifndef NVIDIA_INDEX_HEIGHTFIELD_NORMAL_H
#define NVIDIA_INDEX_HEIGHTFIELD_NORMAL_H

#include "common_utility.h"

#include <cassert>

#include <nv/index/iregular_heightfield_patch.h>

namespace nv {
namespace index_common {

namespace {

// outside patch is the same as the hole.
inline mi::Float32 get_height(
    const mi::Float32*      heightfield_data,
    mi::Sint64              xres,
    mi::Sint64              yres,
    mi::Sint64              x,
    mi::Sint64              y)
{
    if (x < 0 || x >= xres || y < 0 || y >= yres)
    {
        return nv::index::IRegular_heightfield_patch::get_hole_value();
    }
    else
    {
        return heightfield_data[xres * y + x];
    }
}

} // namespace

//----------------------------------------------------------------------
/// Generate heightfield vertex normal with a simple Laplacian.
///
/// The buffer must be pre-allocated.
///
/// \param[in] patch_raw_bbox raw patch bounding box
/// \param[in] height_data    height value data
/// \param[out] normal_data   computed normal vector data
inline void generate_vertex_normal(
    const mi::math::Bbox<mi::Sint64, 2>&           patch_raw_bbox,
    const mi::Float32*                             height_data,
    mi::math::Vector_struct<mi::Float32, 3>*       normal_data)
{
    assert(height_data != 0);
    assert(normal_data != 0);

    const mi::math::Vector<mi::Sint64, 2> patch_raw_dim = patch_raw_bbox.max - patch_raw_bbox.min;
    assert((patch_raw_dim.x > 0) && (patch_raw_dim.y > 0));
    const mi::Sint64 patch_raw_size = patch_raw_dim.x * patch_raw_dim.y;
    assert(patch_raw_size > 0);

    // Compute averaged normals.
    // Define the neighbors: center + 8 neighbors
    const mi::Sint32 vtx_count = 9;
    mi::math::Vector<mi::Sint32, 2> dxy[vtx_count];
    dxy[0] = mi::math::Vector<mi::Sint32, 2>( 0,  0);
    dxy[1] = mi::math::Vector<mi::Sint32, 2>( 1,  0);
    dxy[2] = mi::math::Vector<mi::Sint32, 2>( 1,  1);
    dxy[3] = mi::math::Vector<mi::Sint32, 2>( 0,  1);
    dxy[4] = mi::math::Vector<mi::Sint32, 2>(-1,  1);
    dxy[5] = mi::math::Vector<mi::Sint32, 2>(-1,  0);
    dxy[6] = mi::math::Vector<mi::Sint32, 2>(-1, -1);
    dxy[7] = mi::math::Vector<mi::Sint32, 2>( 0, -1);
    dxy[8] = mi::math::Vector<mi::Sint32, 2>( 1, -1);
    bool use_height[vtx_count];

    mi::Sint64 counter = 0;
    for(mi::Sint64 j = 0; j < patch_raw_dim.y; j++)
    {
        for(mi::Sint64 i = 0; i < patch_raw_dim.x; i++)
        {
            mi::Float32 sample_h[vtx_count];
            for (mi::Sint32 k = 0; k < vtx_count; ++k)
            {
                sample_h[k] = get_height(height_data, 
                                         patch_raw_dim.x, patch_raw_dim.y, 
                                         (i + dxy[k].x),  (j + dxy[k].y));

                use_height[k] = !nv::index::IRegular_heightfield_patch::is_hole(sample_h[k]);
                if (k == 0 && !use_height[k])
                {
                    break; // Early out
                }
            }

            if (!use_height[0])
            {
                // Holes do not get a normal
                normal_data[counter] = mi::math::Vector<mi::Float32, 3>(0.f);
                counter++;
                continue;
            }

            mi::math::Vector<mi::Float32, 3> edge_v[vtx_count];
            edge_v[0] = mi::math::Vector<mi::Float32, 3>(0.f); // entry 0 is not defined
            for (mi::Sint32 k = 1; k < vtx_count; ++k) // start with 1
            {
                // Center difference vector
                edge_v[k] = mi::math::Vector<mi::Float32, 3>(static_cast<mi::Float32>(dxy[k].x), 
                                                             static_cast<mi::Float32>(dxy[k].y), 
                                                             sample_h[k] - sample_h[0]);
                edge_v[k].normalize();
            }

            mi::math::Vector<mi::Float32, 3> nrm(0.f);
            for (mi::Sint32 k1 = 1; k1 < vtx_count; ++k1) // start with 1            
            {
                mi::Sint32 k2 = k1 + 1;
                if (k2 == vtx_count)
                {
                    k2 = 1;      // handle k1 == 8, k2 == 1 case
                }
                assert(k2 != 0);
                assert(k2 < vtx_count);

                if (use_height[k1] && use_height[k2]) // Ignore holes
                {
                    mi::math::Vector<mi::Float32, 3> face_n = mi::math::cross(edge_v[k1], edge_v[k2]);
                    face_n.normalize();
                    nrm += face_n;
                }
            }
            const bool is_nrm_finite = nrm.normalize();
            if (!is_nrm_finite)
            {
                // Set to a valid arbitrary normal instead of 0
                nrm = mi::math::Vector<mi::Float32, 3>(0.f, 0.f, 1.f);
            }
            
            assert(counter < patch_raw_size);
            no_unused_variable_warning_please(patch_raw_size);
            normal_data[counter] = nrm;
            counter++;
        }
    }
}

//----------------------------------------------------------------------
} // namespace index_common
} // namespace nv
#endif // NVIDIA_INDEX_HEIGHTFIELD_NORMAL_H
