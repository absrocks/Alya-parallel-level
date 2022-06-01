/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Bezier surface generation utility

#ifndef NVIDIA_INDEX_BIN_COMMON_BEZIER_SURFACE_UTILITY_H
#define NVIDIA_INDEX_BIN_COMMON_BEZIER_SURFACE_UTILITY_H

#include <mi/dice.h>
#include <cmath>
#include <map>
#include <vector>
#include "forwarding_logger.h"

namespace nv {
namespace index_common {

static const mi::math::Vector_struct<mi::Float32, 3> control_points[8][8] = { 
{
    {20,20,625},{500,20,750},{1000,20,20},{1500,20,20},{2000,20,20},{2500,20,20},{3000,20,1000},{3500,20,500}
},
{
    {20,500,625},{500,500,20},{1000,500,20},{1500,500,20},{2000,500,20},{2500,500,20},{ 3000,500,1000},{3500,500,500}
},
{
    {20,1000,625},{500,1000,20},{1000,1000,20}, {1500,1000,1500},{2000,1000,1500},{2500,1000,20},{3000,1000,1500},{3500,1000,500}
},
{
    {20,1500,625},{500,1500,20},{1000,1500,20}, {1500,1500,1500},{2000,1500,325},{2500,1500,20},{3000,1500,325},{3500,1500,20}
},
{
    {20,2000,625},{500,2000,20},{1000,2000,20},{1500,2000,1500},{2000,2000,1500},{2500,2000,20},{3000,2000,325},{3500,2000,500}
},
{
    {20,2500,625},{500,2500,20},{1000,2500,1500},{1500,2500,1500},{2000,2500,20},{2500,2500,20},{3000,2500,325},{3500,2500,500}
},
{
    {20,3000,625},{500,3000,20},{1000,3000,20},{1500,3000,20},{2000,3000,20},{2500,3000,20},{3000,3000,1000},{3500,3000,500}
},
{
    {20,3500,325},{500,3500,20},{1000,3500,20},{1500,3500,20},{2000,3500,20},{2500,3500,20},{3000,3500,250},{3500,3500,500}
}
};

/// intersection test between a cell and a patch
inline bool intersects_voxel(
    const mi::math::Vector_struct<mi::Float32, 3>&          cell_min,
    const mi::math::Vector_struct<mi::Float32, 3>&          cell_max,
    const mi::math::Bbox_struct<mi::Sint32, 3>&             patch_bounds)
{
    return (patch_bounds.min.x <= cell_max.x && patch_bounds.max.x >= cell_min.x &&
            patch_bounds.min.y <= cell_max.y && patch_bounds.max.y >= cell_min.y &&
            patch_bounds.min.z <= cell_max.z && patch_bounds.max.z >= cell_min.z);
        
} 

/// generate control points
inline void gen_control_points(
    std::map<mi::Uint32, std::vector<mi::math::Vector_struct<mi::Float32, 3> > >&           control_points_map)
{
    mi::Uint32 idx = 0;

    // For now simply have all the control points as they are
    // and later change it accommodate given bbox sizes
    for(mi::Uint32 i = 0; i < 8; ++i)
    {
        std::vector<mi::math::Vector_struct<mi::Float32, 3> > control_points_vec;
        for(mi::Uint32 j = 0; j < 8; ++j)
        {
            control_points_vec.push_back(control_points[i][j]);
            //INFO_LOG<<"point: "<<control_points[i][j];
        }
        if(!control_points_vec.empty())
            control_points_map[idx++] = control_points_vec;
    }
}

/// calculate a point position
inline mi::math::Vector_struct<mi::Float32, 3> calculate_point(
    mi::Float32                                                     t, 
    const std::vector<mi::math::Vector_struct<mi::Float32, 3> >&    pnts) 
{
    // Final v point
    mi::math::Vector_struct<mi::Float32, 3> p;

    // The t value inverted
    mi::Float32 it = 1.0f-t;

    // Calculate bernstein coefficients, currently only upto 8
    mi::Float32 b[8];
    b[7] = t*t*t*t*t*t*t;
    b[6] = 7*t*t*t*t*t*t*it;
    b[5] = 21*t*t*t*t*t*it*it;
    b[4] = 35*t*t*t*t*it*it*it;
    b[3] = 35*t*t*t*it*it*it*it;
    b[2] = 21*t*t*it*it*it*it*it;
    b[1] = 7*t*it*it*it*it*it*it;
    b[0] = it*it*it*it*it*it*it;
    
    // use the blending functions
    p.x = p.y = p.z = 0;
    for(mi::Uint32 i= 0; i < pnts.size(); ++i)
    {   
        p.x += b[i]*pnts[i].x;
        p.y += b[i]*pnts[i].y;
        p.z += b[i]*pnts[i].z;
    }
    return p;
}

/// estimate how many cells in the patch
inline mi::Uint32 estimate_patch_cells(
    const mi::Uint32&                                       layers,
    const mi::Uint32&                                       steps_u,
    const mi::Uint32&                                       steps_v,
    const mi::math::Bbox_struct<mi::Float32, 3>&            global_bbox,    
    const mi::math::Bbox_struct<mi::Sint32, 3>&             patch_bounds)
{
    //max columns and rows, add 1 to have right number of coordinates
    mi::Uint32 max_steps_u = steps_u+1;
    mi::Uint32 max_steps_v = steps_v+1;
    
    // Final row/column count
    mi::Uint32 row = 0;
    mi::Uint32 col = 0;
    
    mi::Float32 one_cell_distance_u = 1.0f/(static_cast<mi::Float32>(max_steps_u)-1.0f);
    mi::Float32 one_cell_distance_v = 1.0f/(static_cast<mi::Float32>(max_steps_v)-1.0f);
    mi::Float32 u_bound_min = (static_cast<mi::Float32>(patch_bounds.min.x) / global_bbox.max.x);
    mi::Float32 u_bound_max = (static_cast<mi::Float32>(patch_bounds.max.x) / global_bbox.max.x);
    mi::Float32 v_bound_min = (static_cast<mi::Float32>(patch_bounds.min.y) / global_bbox.max.y);
    mi::Float32 v_bound_max = (static_cast<mi::Float32>(patch_bounds.max.y) / global_bbox.max.y);
 
    for(mi::Uint32 i=0;i!=max_steps_v;++i)
    {
        // Calculate the parametric u value
        mi::Float32 v = static_cast<mi::Float32>(i)/(static_cast<mi::Float32>(max_steps_v)-1.0f);
        if(v >= v_bound_min - 2.0f * one_cell_distance_v && v <= v_bound_max + 2.0f * one_cell_distance_v)
        { 
            mi::Uint32 prev_col = 0;
            for(mi::Uint32 j=0;j!=max_steps_u;++j)
            {
                // Calculate the parametric v value
                mi::Float32 u = static_cast<mi::Float32>(j)/(static_cast<mi::Float32>(max_steps_u)-1.0f);
                if(u >= u_bound_min - 2.0f * one_cell_distance_u && u <= u_bound_max + 2.0f * one_cell_distance_u)
                {
                    prev_col++;
                }
            }
            if(prev_col > col)
                col = prev_col;
            row++;
        }
    }
    
    // Calculate indices for the each layer
    mi::Uint32 i = 0;
    mi::Uint32 total_cells = 0;
    for(mi::Uint32 layer = 1; layer <= layers; ++layer)
    {
        i = 0;
        while(true)
        {
            i++;
            if(i > (col - 1) * (row - 1)+(row - 2))
                break;
         
            if(i >= (col-1) && (i%col == 0))
                i++;
            total_cells++;
        }
    }

    return total_cells;
}

/// Generate a bezier_surface       
inline void generate_patch_cells(
    const mi::Uint32&                                       layers, 
    const mi::Uint32&                                       steps_u,
    const mi::Uint32&                                       steps_v,
    const mi::math::Bbox_struct<mi::Float32, 3>&            global_bbox,    
    const mi::math::Bbox_struct<mi::Sint32, 3>&             patch_bounds,
    std::vector<mi::math::Vector_struct<mi::Float32, 3> >&  cell_coordinates)
{

    // Generate control points based on the patch_bounds
    std::map<mi::Uint32, std::vector<mi::math::Vector_struct<mi::Float32, 3> > >  control_points_map;
    gen_control_points(control_points_map);
    
    //max columns and rows, add 1 to have right number of coordinates
    mi::Uint32 max_steps_u = steps_u+1;
    mi::Uint32 max_steps_v = steps_v+1;

    mi::Float32 one_cell_distance_u = 1.0f/(static_cast<mi::Float32>(max_steps_u)-1.0f);
    mi::Float32 one_cell_distance_v = 1.0f/(static_cast<mi::Float32>(max_steps_v)-1.0f);
    mi::Float32 u_bound_min = (static_cast<mi::Float32>(patch_bounds.min.x) / global_bbox.max.x);
    mi::Float32 u_bound_max = (static_cast<mi::Float32>(patch_bounds.max.x) / global_bbox.max.x);
    mi::Float32 v_bound_min = (static_cast<mi::Float32>(patch_bounds.min.y) / global_bbox.max.y);
    mi::Float32 v_bound_max = (static_cast<mi::Float32>(patch_bounds.max.y) / global_bbox.max.y);
 
    mi::Uint32 row = 0;
    mi::Uint32 col = 0;
    std::vector<mi::math::Vector_struct<mi::Float32, 3> > bezier_points_vec; 
    for(mi::Uint32 i=0;i!=max_steps_v;++i)
    {
        // Calculate the parametric u value
        mi::Float32 v = static_cast<mi::Float32>(i)/(static_cast<mi::Float32>(max_steps_v)-1.0f);
        if(v >= v_bound_min - 2.0f * one_cell_distance_v && v <= v_bound_max + 2.0f * one_cell_distance_v)
        { 
            mi::Uint32 prev_col = 0;
            for(mi::Uint32 j=0;j!=max_steps_u;++j)
            {
                // Calculate the parametric v value
                mi::Float32 u = static_cast<mi::Float32>(j)/(static_cast<mi::Float32>(max_steps_u)-1.0f);
                if(u >= u_bound_min - 2.0f * one_cell_distance_u && u <= u_bound_max + 2.0f * one_cell_distance_u)
                {
                    // Curve-sweep on u direction
                    std::vector<mi::math::Vector_struct<mi::Float32, 3> > sweep_points;
                    std::map<mi::Uint32, std::vector<mi::math::Vector_struct<mi::Float32, 3> > >::const_iterator itr = control_points_map.begin();
                    for(; itr!=control_points_map.end(); ++itr)
                    {
                        const std::vector<mi::math::Vector_struct<mi::Float32, 3> > control_points_vec = itr->second;
                        sweep_points.push_back(calculate_point(u, control_points_vec));
                    }                   
                    // Final point on the surface    
                    mi::math::Vector_struct<mi::Float32, 3> p = calculate_point(v, sweep_points);
                    bezier_points_vec.push_back(p);
                    prev_col++;
                }
            }
            if(prev_col > col)
                col = prev_col;
            row++;
        }
    }
    
    // Return if not enough coordinates are generated
    if(bezier_points_vec.empty() ||
        row < 2 || col < 2) return;
    
    // Resize to fit coordinates for all the layers
    mi::Uint32 grid_offset = col*row;
    bezier_points_vec.resize(grid_offset*(layers+1));
    
    // Stack bezier points to form cell layers
    mi::Uint32 zoffset_per_layer = 0;
    for(mi::Uint32 i = 0; i < grid_offset; ++i)
    {                
        for(mi::Uint32 layer = 1; layer <= layers; ++layer)
        {   
            zoffset_per_layer += 50;
            bezier_points_vec[i+grid_offset*layer] = bezier_points_vec[i];
            bezier_points_vec[i+grid_offset*layer].z = bezier_points_vec[i].z + static_cast<mi::Float32>(zoffset_per_layer);
        }
        zoffset_per_layer = 0;
    }
    
    // Calculate indices for the each layer
    mi::Uint32 i = 0;
    mi::Uint32 fixed_offset = 0;
    std::vector<mi::Uint32> grid_indices;
    for(mi::Uint32 layer = 1; layer <= layers; ++layer)
    {
        if(layer > 1)
            fixed_offset = grid_offset*(layer-1);
        
        i = 0;
        while(true)
        {
            i++;
            if(i > (col - 1) * (row - 1)+(row - 2))
                break;
         
            if(i >= (col-1) && (i%col == 0))
                i++;
            
            // Indices for each cell
            {
                grid_indices.push_back(i+fixed_offset);
                grid_indices.push_back(i+col+fixed_offset);
                grid_indices.push_back(i+col+grid_offset*layer);
                grid_indices.push_back(i+grid_offset*layer);
                grid_indices.push_back(i+1+fixed_offset);
                grid_indices.push_back(i+1+col+fixed_offset);
                grid_indices.push_back(i+1+col+grid_offset*layer);
                grid_indices.push_back(i+1+grid_offset*layer);                
            }
        }
    }

    for(mi::Uint32 i = 0; i < grid_indices.size(); ++i)
    {
        mi::Uint32  idx = grid_indices[i]-1;
        cell_coordinates.push_back(bezier_points_vec[idx]);
    }
}


//----------------------------------------------------------------------
}} // namespace nv::index_common
#endif // NVIDIA_INDEX_BIN_COMMON_BEZIER_SURFACE_UTILITY_H
