/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief external reservoir_grid data importer job


#include "synthetic_reservoir_patch_generator.h"

#include <nv/index/itriangle_mesh_subset.h>
#include <nv/index/itriangle_mesh_scene_element.h>
#include <zlib.h>
#include "bezier_surface_utility.h"
#include "forwarding_logger.h"
#include "common_utility.h"

namespace nv {
namespace index_common {

//----------------------------------------------------------------------
Synthetic_reservoir_patch_generator::Synthetic_reservoir_patch_generator(
    const mi::Uint32&                                   layers,
    const mi::Uint32&                                   cells_per_row,
    const mi::Uint32&                                   cells_per_column,
    const mi::math::Bbox_struct<mi::Float32, 3>&        global_bbox,
    const std::string &                                 cache_file,
    const std::string &                                 cache_type)
    :
    m_layers(layers),
    m_cells_per_row(cells_per_row),
    m_cells_per_column(cells_per_column),
    m_global_bbox(global_bbox),
    m_cache_file(cache_file),
    m_cache_type(cache_type)
{
    if(layers < 1)
    {
        WARN_LOG<<"Number of layers cannot be less than 1";
        m_layers = 1;
    }

    if(cells_per_row < 1 || cells_per_column < 1)
    {
        WARN_LOG<<"Number of cell count cannot be less than 1";
        m_cells_per_row = 1;
        m_cells_per_column = 1;
    }

    std::ostringstream s;
    s << "type=reservoir_grid\n"
      << "importer=artificial\n"
      << "layers=" << layers << "\n"
      << "cells_per_row=" << cells_per_row << "\n"
      << "cells_per_column=" << cells_per_column << "\n"
      << "cache_type=" << cache_type << "\n";
    m_configuration = s.str();
}

//----------------------------------------------------------------------
Synthetic_reservoir_patch_generator::Synthetic_reservoir_patch_generator(
): m_layers(0),
   m_cells_per_row(0),
   m_cells_per_column(0)

{
    // for serialization and copy only
}

//----------------------------------------------------------------------
Synthetic_reservoir_patch_generator::~Synthetic_reservoir_patch_generator()
{
    // empty
}

//----------------------------------------------------------------------
namespace
{
// Colormap generating utility function
static inline void generate_colormap(std::vector<mi::math::Color_struct>& color_entries)
{
   color_entries.resize(256);
   
   mi::Float32 vmin = 0.f; 
   mi::Float32 vmax = 0.9f;
   for(mi::Uint32 i=0; i<256; ++i)
   { 
       mi::Float32 v = (static_cast<mi::Float32>(i)+0.5f)/256.f;
       
       mi::math::Color_struct c = {1.f, 1.f, 1.f, 1.f};
       mi::Float32 dv;
        
       if (v < vmin)
          v = vmin;
       if (v > vmax)
          v = vmax;
       dv = vmax - vmin;

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
       color_entries[i] = c;
    }
}
} // anonymous namespace

//----------------------------------------------------------------------
mi::Size Synthetic_reservoir_patch_generator::estimate(
    const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{
    // Estimate number of cells without generating the actual points/cells
    std::vector<mi::math::Vector_struct<mi::Float32, 3> >  bezier_points_vec;
    mi::Uint32 total_cells = nv::index_common::estimate_patch_cells(
                                                m_layers,
                                                m_cells_per_row, 
                                                m_cells_per_column, 
                                                m_global_bbox, 
                                                bounding_box);
    // Total size in bytes
    mi::Size total_size = (total_cells * 
                            (8 * sizeof(mi::math::Vector_struct<mi::Float32, 3>)            // coordinates of the cell
                            + 6 * sizeof(mi::math::Vector<mi::Float32, 3>)                  // normals of the cell
                            + 36 * sizeof(mi::Uint32)                                       // triangle position indices of the cell
                            + 36 * sizeof(mi::Uint32)                                       // triangle normal indices of the cell
                            + 36 * sizeof(mi::Uint32)                                       // triangle colormap indices of the cell
                            + 12 * sizeof(nv::index::ITriangle_mesh_subset::Triflags))      // triangle wireframe flags indices of the cell    
                          );    
    //INFO_LOG <<"Memory consumption of "<< total_cells <<" cells in "<< bounding_box <<" is: "<< total_size <<" bytes";

    return total_size;

}

nv::index::IDistributed_data_subset* Synthetic_reservoir_patch_generator::create(
    const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
    nv::index::IData_subset_factory*                factory,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{
    std::ostringstream cache_file_name_str;
    cache_file_name_str << m_cache_file << "_"
                        << bounding_box.min.x << "." << bounding_box.max.x << "x"
                        << bounding_box.min.y << "." << bounding_box.max.y << "x"
                        << bounding_box.min.z << "." << bounding_box.max.z
                        << "_reservoir_synthetic_importer.tmp";

    std::string cache_file_name;
    cache_file_name = cache_file_name_str.str();

    std::ostringstream cache_file_name_str_gz;
    cache_file_name_str_gz << m_cache_file << "_"
                        << bounding_box.min.x << "." << bounding_box.max.x << "x"
                        << bounding_box.min.y << "." << bounding_box.max.y << "x"
                        << bounding_box.min.z << "." << bounding_box.max.z
                        << "_reservoir_synthetic_importer.tmp.gz";

    std::string cache_file_name_gz;
    cache_file_name_gz = cache_file_name_str_gz.str();

    FILE * cache_file = fopen(cache_file_name.c_str(), "rb");
    gzFile cache_file_gz = gzopen(cache_file_name_gz.c_str(), "rb");

    std::vector<mi::math::Vector_struct<mi::Float32, 2> > mesh_tex_coords;
    std::vector<mi::math::Color_struct > mesh_vertex_colors;
    std::vector<mi::Uint32> mesh_tex_coord_indices;
    std::vector<mi::Uint32> mesh_vertex_color_indices;
    std::vector<mi::Uint32> mesh_vertex_colormap_indices;

    // Materials
    std::vector<mi::Uint16> mesh_materials;
    mesh_materials.push_back(0);
    
    // triangle flags
    std::vector<nv::index::ITriangle_mesh_subset::Triflags> triangle_flags;
    
    // Get points for the reservoir based on bezier surface
    std::vector<mi::math::Vector_struct<mi::Float32, 3> > mesh_vectors;
    
    // Calculate Face normal per quad  
    std::vector<mi::math::Vector_struct<mi::Float32, 3> > mesh_normals;
    std::vector<mi::Uint32> mesh_normal_indices;
    
    // Triangulate each cell
    std::vector<mi::Uint32> mesh_pos_indices;
    
    //
    // Read data from the cache_file
    //
    if (cache_file && m_cache_type == "on")
    {
        INFO_LOG<<"Reading data from uncompressed cache file "<<cache_file_name;
        mi::Float64 start_time = nv::index_common::get_time();

        mi::Uint32 count = 0;
        if(fread(&count, sizeof(mi::Uint32), 1, cache_file))
            mesh_vectors.resize(count);
        if(count > 0 && fread(&mesh_vectors[0], sizeof(mi::Float32)*3, count, cache_file) != count)
            ERROR_LOG<<"Failed reading mesh_vectors from cache_file";

        count = 0;
        if(fread(&count, sizeof(mi::Uint32), 1, cache_file))
            mesh_normals.resize(count);
        if(count > 0 && fread(&mesh_normals[0], sizeof(mi::Float32), count*3, cache_file) != count*3)
            ERROR_LOG<<"Failed reading mesh_normals from cache_file";

        count = 0;
        if(fread(&count, sizeof(mi::Uint32), 1, cache_file))
            mesh_tex_coords.resize(count);
        if(count > 0 && fread(&mesh_tex_coords[0], sizeof(mi::Float32), count*2, cache_file) != count*2)
            ERROR_LOG<<"Failed reading mesh_tex_coords from cache_file";

        count = 0;
        if(fread(&count, sizeof(mi::Uint32), 1, cache_file))
            mesh_vertex_colors.resize(count);
        if(count > 0 && fread(&mesh_vertex_colors[0], sizeof(mi::Float32), count*4, cache_file) != count*4)
            ERROR_LOG<<"Failed reading mesh_vertex_colors from cache_file";

        count = 0;
        if(fread(&count, sizeof(mi::Uint32), 1, cache_file))
            mesh_pos_indices.resize(count);
        if(count > 0 && fread(&mesh_pos_indices[0], sizeof(mi::Uint32), count, cache_file) != count)
            ERROR_LOG<<"Failed reading mesh_pos_indices from cache_file";

        count = 0;
        if(fread(&count, sizeof(mi::Uint32), 1, cache_file))
            mesh_normal_indices.resize(count);
        if(count > 0 && fread(&mesh_normal_indices[0], sizeof(mi::Uint32), count, cache_file) != count)
            ERROR_LOG<<"Failed reading mesh_normal_indices from cache_file";

        count = 0;
        if(fread(&count, sizeof(mi::Uint32), 1, cache_file))
            mesh_tex_coord_indices.resize(count);
        if(count > 0 && fread(&mesh_tex_coord_indices[0], sizeof(mi::Uint32), count, cache_file) != count)
            ERROR_LOG<<"Failed reading mesh_tex_coord_indices from cache_file";

        count = 0;
        if(fread(&count, sizeof(mi::Uint32), 1, cache_file))
            mesh_vertex_color_indices.resize(count);
        if(count > 0 && fread(&mesh_vertex_color_indices[0], sizeof(mi::Uint32), count, cache_file) != count)
            ERROR_LOG<<"Failed reading mesh_vertex_color_indices from cache_file";

        count = 0;
        if(fread(&count, sizeof(mi::Uint32), 1, cache_file))
            mesh_vertex_colormap_indices.resize(count);
        if(count > 0 && fread(&mesh_vertex_colormap_indices[0], sizeof(mi::Uint32), count, cache_file) != count)
            ERROR_LOG<<"Failed reading mesh_vertex_colormap_indices from cache_file";

        count = 0;
        if(fread(&count, sizeof(mi::Uint32), 1, cache_file))
            mesh_materials.resize(count);
        if(count > 0 && fread(&mesh_materials[0], sizeof(mi::Uint16), count, cache_file) != count)
            ERROR_LOG<<"Failed reading mesh_materials from cache_file";

        count = 0;
        if(fread(&count, sizeof(mi::Uint32), 1, cache_file))
            triangle_flags.resize(count);
        if(count > 0 && fread(&triangle_flags[0], sizeof(nv::index::ITriangle_mesh_subset::Triflags), count, cache_file) != count)
            ERROR_LOG<<"Failed reading triangle_flags from cache_file";

        mi::math::Bbox_struct<mi::Float32, 3> bbox;
        fread(&bbox.min.x, sizeof(mi::Float32), 6, cache_file);

        mi::Float64 end_time_file = nv::index_common::get_time();

        fclose(cache_file);

        // ---------------------------------------------------------------------------------------------------
        // Create submesh (must be derived from ITriangle_mesh_subset)
        // and returns mesh to the caller. The caller takes ownership!
        mi::base::Handle<nv::index::ITriangle_mesh_subset> grid_mesh(
            factory->create<nv::index::ITriangle_mesh_subset>());
        if(grid_mesh.is_valid_interface() && !mesh_vectors.empty() )
        {
            // Generate global triangle IDs: 
            std::vector<mi::Uint64> global_triid;
            assert((mesh_pos_indices.size() % 3) == 0);
            const mi::Size nb_tris = mesh_pos_indices.size() / 3;
            global_triid.resize(nb_tris); 
            for(mi::Size i = 0; i < nb_tris; ++i){
                global_triid[i] = ~0UL; // set as invalid global triangle ID.
            }

            // check necessary data for creating a valid reservoir grid.
            if ((!mi::math::Bbox<mi::Float32, 3>(bbox).is_volume()) ||
                mesh_vectors.empty() || mesh_pos_indices.empty() || global_triid.empty()) {
                ERROR_LOG << "Cannot create a valid reservoir grid.";
                return 0;
            }

            const bool success = grid_mesh->initialize(
                bbox,
                &mesh_vectors[0],     mesh_vectors.size(),
                &mesh_pos_indices[0], mesh_pos_indices.size(),
                &global_triid[0],     global_triid.size(),

                mesh_normals.empty()                 ? 0 : &mesh_normals[0],                    mesh_normals.size(),
                mesh_tex_coords.empty()              ? 0 : &mesh_tex_coords[0],                 mesh_tex_coords.size(),
                mesh_vertex_colors.empty()           ? 0 : &mesh_vertex_colors[0],              mesh_vertex_colors.size(),

                mesh_normal_indices.empty()          ? 0 : &mesh_normal_indices[0],             
                mesh_tex_coord_indices.empty()       ? 0 : &mesh_tex_coord_indices[0],          
                mesh_vertex_color_indices.empty()    ? 0 : &mesh_vertex_color_indices[0],       
                mesh_vertex_colormap_indices.empty() ? 0 : &mesh_vertex_colormap_indices[0],

                mesh_materials.empty()               ? 0 : &mesh_materials[0],                  mesh_materials.size(),
                triangle_flags.empty()               ? 0 : &triangle_flags[0],                  triangle_flags.size());
            if(success)
            {
                INFO_LOG << "Took " << (end_time_file-start_time)
                         << " seconds to read " << mesh_pos_indices.size()/36
                         << " cells (" <<mesh_pos_indices.size()/3
                         << " triangles) " <<(nv::index_common::get_time()-end_time_file)
                         << " seconds to create the grid";
                grid_mesh->retain(); // since the handle will be out of scope.
                return grid_mesh.get();
            }
        }
        return 0;
    }
    else if(cache_file_gz && m_cache_type == "compressed")
    {
        INFO_LOG<<"Reading data from compressed cache file "<<cache_file_name_gz;
        mi::Float64 start_time = nv::index_common::get_time();

        mi::Uint32 count = 0;
        if(gzread(cache_file_gz, &count, sizeof(mi::Uint32)))
            mesh_vectors.resize(count);
        gzread(cache_file_gz, &mesh_vectors[0], count*sizeof(mi::Float32)*3);

        count = 0;
        if(gzread(cache_file_gz, &count, sizeof(mi::Uint32)))
            mesh_normals.resize(count);
        gzread(cache_file_gz, &mesh_normals[0], count*sizeof(mi::Float32)*3);

        count = 0;
        if(gzread(cache_file_gz, &count, sizeof(mi::Uint32)))
            mesh_tex_coords.resize(count);
        gzread(cache_file_gz, &mesh_tex_coords[0], count*sizeof(mi::Float32)*2);

        count = 0;
        if(gzread(cache_file_gz, &count, sizeof(mi::Uint32)))
            mesh_vertex_colors.resize(count);
        gzread(cache_file_gz, &mesh_vertex_colors[0], count*sizeof(mi::Float32)*4);

        count = 0;
        if(gzread(cache_file_gz, &count, sizeof(mi::Uint32)))
            mesh_pos_indices.resize(count);
        gzread(cache_file_gz, &mesh_pos_indices[0], count*sizeof(mi::Uint32));

        count = 0;
        if(gzread(cache_file_gz, &count, sizeof(mi::Uint32)))
            mesh_normal_indices.resize(count);
        gzread(cache_file_gz, &mesh_normal_indices[0], count*sizeof(mi::Uint32));

        count = 0;
        if(gzread(cache_file_gz, &count, sizeof(mi::Uint32)))
            mesh_tex_coord_indices.resize(count);
        gzread(cache_file_gz, &mesh_tex_coord_indices[0], count*sizeof(mi::Uint32));

        count = 0;
        if(gzread(cache_file_gz, &count, sizeof(mi::Uint32)))
            mesh_vertex_color_indices.resize(count);
        gzread(cache_file_gz, &mesh_vertex_color_indices[0], count*sizeof(mi::Uint32));

        count = 0;
        if(gzread(cache_file_gz, &count, sizeof(mi::Uint32)))
            mesh_vertex_colormap_indices.resize(count);
        gzread(cache_file_gz, &mesh_vertex_colormap_indices[0], count*sizeof(mi::Uint32));

        count = 0;
        if(gzread(cache_file_gz, &count, sizeof(mi::Uint32)))
            mesh_materials.resize(count);
        gzread(cache_file_gz, &mesh_materials[0], count*sizeof(mi::Uint16));

        count = 0;
        if(gzread(cache_file_gz, &count, sizeof(mi::Uint32)))
            triangle_flags.resize(count);
        gzread(cache_file_gz, &triangle_flags[0], count*sizeof(nv::index::ITriangle_mesh_subset::Triflags));

        mi::math::Bbox_struct<mi::Float32, 3> bbox;
        gzread(cache_file_gz, &bbox.min.x, sizeof(bbox));

        mi::Float64 end_time_file = nv::index_common::get_time();

        gzclose(cache_file_gz);

        // ---------------------------------------------------------------------------------------------------
        // Create submesh (must be derived from ITriangle_mesh_subset)
        // and returns mesh to the caller. The caller takes ownership!
        mi::base::Handle<nv::index::ITriangle_mesh_subset> grid_mesh(
            factory->create<nv::index::ITriangle_mesh_subset>());
        if(grid_mesh.is_valid_interface() && !mesh_vectors.empty())
        {
            // Generate global triangle IDs: 
            std::vector<mi::Uint64> global_triid;
            assert((mesh_pos_indices.size() % 3) == 0);
            const mi::Size nb_tris = mesh_pos_indices.size() / 3;
            global_triid.resize(nb_tris); 
            for(mi::Size i = 0; i < nb_tris; ++i){
                global_triid[i] = ~0UL; // set as invalid global triangle ID.
            }
            
            // check necessary data for creating a valid reservoir grid.
            if ((!mi::math::Bbox<mi::Float32, 3>(bbox).is_volume()) ||
                mesh_vectors.empty() || mesh_pos_indices.empty() || global_triid.empty()) {
                ERROR_LOG << "Cannot create a valid reservoir grid.";
                return 0;
            }

            const bool success = grid_mesh->initialize(
                bbox,
                &mesh_vectors[0],                    mesh_vectors.size(),
                &mesh_pos_indices[0],                mesh_pos_indices.size(),
                &global_triid[0],                    global_triid.size(),

                &mesh_normals[0],                    mesh_normals.size(),
                &mesh_tex_coords[0],                 mesh_tex_coords.size(),
                &mesh_vertex_colors[0],              mesh_vertex_colors.size(),

                &mesh_normal_indices[0],             
                &mesh_tex_coord_indices[0],          
                &mesh_vertex_color_indices[0],       
                &mesh_vertex_colormap_indices[0],
                
                &mesh_materials[0],                  mesh_materials.size(),
                &triangle_flags[0],                  triangle_flags.size());
            if(success)
            {
                INFO_LOG << "Took " << (end_time_file-start_time)
                         << " seconds to read " << mesh_pos_indices.size()/36
                         << " cells (" << mesh_pos_indices.size()/3
                         << " triangles) " <<(nv::index_common::get_time()-end_time_file)
                         << " seconds to create the grid";
                grid_mesh->retain(); // since the handle will be out of scope.
                return grid_mesh.get();
            }
        }

        return  NULL;
    }
            
    //const mi::Float64 start_time = current_system_time();
    nv::index_common::generate_patch_cells(m_layers, m_cells_per_row, m_cells_per_column, m_global_bbox, bounding_box, mesh_vectors);

    //mi::Float32 elapsed_time = current_system_time() - start_time;
     
    if (mesh_vectors.empty())
        WARN_LOG<< "No bezier coordinates generated!";
        
    //INFO_LOG<<"Took ("<<(elapsed_time)<<") seconds to generate the bezier coordinates";
    
    mi::Uint32 total_cells = mesh_vectors.size()/8;
    mi::Uint32 total_quads = total_cells*6;

    // For normals, since they are per face
    std::vector<mi::Uint32> quad_indices;
    
    // create visible edge mask....for this case add all edges except the edge(2,0)
    nv::index::ITriangle_mesh_subset::Triflags visible_edges = (nv::index::ITriangle_mesh_subset::Triflags)
                                    (nv::index::ITriangle_mesh_subset::TRIFLAGS_VISIBLE_EDGE_01 |
                                    nv::index::ITriangle_mesh_subset::TRIFLAGS_VISIBLE_EDGE_12);
                                    
    for(mi::Uint32 i=0; i < total_cells; ++i)
    {
        // Indices for the triangles
        mesh_pos_indices.push_back(0+8*i);
        mesh_pos_indices.push_back(1+8*i);
        mesh_pos_indices.push_back(5+8*i);

        mesh_pos_indices.push_back(5+8*i);
        mesh_pos_indices.push_back(4+8*i);
        mesh_pos_indices.push_back(0+8*i);
        
        mesh_pos_indices.push_back(6+8*i);
        mesh_pos_indices.push_back(2+8*i);
        mesh_pos_indices.push_back(3+8*i);
        
        mesh_pos_indices.push_back(3+8*i);
        mesh_pos_indices.push_back(7+8*i);
        mesh_pos_indices.push_back(6+8*i);
       
        mesh_pos_indices.push_back(3+8*i);
        mesh_pos_indices.push_back(0+8*i);
        mesh_pos_indices.push_back(4+8*i);
        
        mesh_pos_indices.push_back(4+8*i);
        mesh_pos_indices.push_back(7+8*i);
        mesh_pos_indices.push_back(3+8*i);

        mesh_pos_indices.push_back(1+8*i);
        mesh_pos_indices.push_back(2+8*i);
        mesh_pos_indices.push_back(6+8*i);
        
        mesh_pos_indices.push_back(6+8*i);
        mesh_pos_indices.push_back(5+8*i);
        mesh_pos_indices.push_back(1+8*i);
        
        mesh_pos_indices.push_back(2+8*i);
        mesh_pos_indices.push_back(1+8*i);
        mesh_pos_indices.push_back(0+8*i);
        
        mesh_pos_indices.push_back(0+8*i);
        mesh_pos_indices.push_back(3+8*i);
        mesh_pos_indices.push_back(2+8*i);
        
        mesh_pos_indices.push_back(4+8*i);
        mesh_pos_indices.push_back(5+8*i);
        mesh_pos_indices.push_back(6+8*i);
        
        mesh_pos_indices.push_back(6+8*i);
        mesh_pos_indices.push_back(7+8*i);
        mesh_pos_indices.push_back(4+8*i);
        
        // Also have indices for the quads
        quad_indices.push_back(0+8*i);
        quad_indices.push_back(1+8*i);
        quad_indices.push_back(5+8*i);
        quad_indices.push_back(4+8*i);
        
        quad_indices.push_back(7+8*i);
        quad_indices.push_back(6+8*i);
        quad_indices.push_back(2+8*i);
        quad_indices.push_back(3+8*i);
        
        quad_indices.push_back(3+8*i);
        quad_indices.push_back(0+8*i);
        quad_indices.push_back(4+8*i);
        quad_indices.push_back(7+8*i);
        
        quad_indices.push_back(1+8*i);
        quad_indices.push_back(2+8*i);
        quad_indices.push_back(6+8*i);
        quad_indices.push_back(5+8*i);
        
        quad_indices.push_back(3+8*i);
        quad_indices.push_back(2+8*i);
        quad_indices.push_back(1+8*i);
        quad_indices.push_back(0+8*i);
        
        quad_indices.push_back(4+8*i);
        quad_indices.push_back(5+8*i);
        quad_indices.push_back(6+8*i);
        quad_indices.push_back(7+8*i);
        
        // color_idx here mimicks property value of a cell in a reservoir
        mi::Uint32 color_idx = static_cast<mi::Uint32>(
            (255.0f*(mesh_vectors[0+8*i].z- m_global_bbox.min.z) /
             (m_global_bbox.max.z - m_global_bbox.min.z)));

        for(mi::Uint32 j = 0; j<36; ++j)
            mesh_vertex_colormap_indices.push_back(color_idx);
 
        // visible edges
        for(mi::Uint32 j = 0; j<12; ++j)
            triangle_flags.push_back(visible_edges);
    }

    for(mi::Uint32 i=0; i < total_quads; ++i)
    {
        mi::math::Vector<mi::Float32, 3> a(mesh_vectors[quad_indices[i*4+0]]);
        mi::math::Vector<mi::Float32, 3> b(mesh_vectors[quad_indices[i*4+1]]);
        mi::math::Vector<mi::Float32, 3> d(mesh_vectors[quad_indices[i*4+3]]);
        
        mi::math::Vector<mi::Float32, 3> v1 = b - a;
        mi::math::Vector<mi::Float32, 3> v2 = d - a;
        
        v1 = mi::math::cross(v1, v2);
        mesh_normals.push_back(v1);
        
        // Generate normal indices for the triangles
        for(mi::Uint32 j=0; j < 6; ++j)
            mesh_normal_indices.push_back(i);
    }
    
    // Create the final triangle mesh of the reservoir grid
    mi::math::Bbox_struct<mi::Float32, 3> bbox;
    bbox.min.x = static_cast<mi::Float32>(bounding_box.min.x);
    bbox.min.y = static_cast<mi::Float32>(bounding_box.min.y);
    bbox.min.z = static_cast<mi::Float32>(bounding_box.min.z);
    bbox.max.x = static_cast<mi::Float32>(bounding_box.max.x);
    bbox.max.y = static_cast<mi::Float32>(bounding_box.max.y);
    bbox.max.z = static_cast<mi::Float32>(bounding_box.max.z);

    INFO_LOG << "Loading (" << total_cells
             << ") out of (" << (m_cells_per_row)*(m_cells_per_column)*m_layers
             << ") total reservoir cells";
    
    if(m_cache_type == "on")
    {
        INFO_LOG<<"Cache file not found, writing cache file without compression "<<cache_file_name;

        mi::Float64 start_time_cache = nv::index_common::get_time();

        FILE * cache_file = fopen(cache_file_name.c_str(), "wb");

        mi::Uint32 count = mesh_vectors.size();
        fwrite(&count, sizeof(mi::Uint32), 1, cache_file);
        if (count > 0) fwrite(&mesh_vectors[0], sizeof(mi::Float32)*3, count, cache_file);

        count = mesh_normals.size();
        fwrite(&count, sizeof(mi::Uint32), 1, cache_file);
        if (count > 0) fwrite(&mesh_normals[0], sizeof(mi::Float32), count*3, cache_file);

        count = mesh_tex_coords.size();
        fwrite(&count, sizeof(mi::Uint32), 1, cache_file);
        if (count > 0) fwrite(&mesh_tex_coords[0], sizeof(mi::Float32), count*2, cache_file);

        count = mesh_vertex_colors.size();
        fwrite(&count, sizeof(mi::Uint32), 1, cache_file);
        if (count > 0) fwrite(&mesh_vertex_colors[0], sizeof(mi::Float32), count*4, cache_file);

        count = mesh_pos_indices.size();
        fwrite(&count, sizeof(mi::Uint32), 1, cache_file);
        if (count > 0) fwrite(&mesh_pos_indices[0], sizeof(mi::Uint32), count, cache_file);

        count = mesh_normal_indices.size();
        fwrite(&count, sizeof(mi::Uint32), 1, cache_file);
        if (count > 0) fwrite(&mesh_normal_indices[0], sizeof(mi::Uint32), count, cache_file);

        count = mesh_tex_coord_indices.size();
        fwrite(&count, sizeof(mi::Uint32), 1, cache_file);
        if (count > 0) fwrite(&mesh_tex_coord_indices[0], sizeof(mi::Uint32), count, cache_file);

        count = mesh_vertex_color_indices.size();
        fwrite(&count, sizeof(mi::Uint32), 1, cache_file);
        if (count > 0) fwrite(&mesh_vertex_color_indices[0], sizeof(mi::Uint32), count, cache_file);

        count = mesh_vertex_colormap_indices.size();
        fwrite(&count, sizeof(mi::Uint32), 1, cache_file);
        if (count > 0) fwrite(&mesh_vertex_colormap_indices[0], sizeof(mi::Uint32), count, cache_file);

        count = mesh_materials.size();
        fwrite(&count, sizeof(mi::Uint32), 1, cache_file);
        if (count > 0) fwrite(&mesh_materials[0], sizeof(mi::Uint16), count, cache_file);

        count = triangle_flags.size();
        fwrite(&count, sizeof(mi::Uint32), 1, cache_file);
        if (count > 0) fwrite(&triangle_flags[0], sizeof(nv::index::ITriangle_mesh_subset::Triflags), count, cache_file);

        fwrite(&bbox.min.x, sizeof(mi::Float32), 6, cache_file);

        fclose(cache_file);

        mi::Float64 end_time_cache = nv::index_common::get_time();
        INFO_LOG<<(end_time_cache-start_time_cache)<<" seconds to write the cache file";
    }
    else if(m_cache_type == "compressed")
    {
        INFO_LOG<<"Cache file not found, writing cache file with compression "<<cache_file_name_gz;

        mi::Float64 start_time_cache = nv::index_common::get_time();
        gzFile cache_file_gz = gzopen(cache_file_name_gz.c_str(), "wb6");

        mi::Uint32 count = mesh_vectors.size();
        gzwrite(cache_file_gz, &count, sizeof(mi::Uint32));
        gzwrite(cache_file_gz, &mesh_vectors[0], count*sizeof(mi::Float32)*3);

        count = mesh_normals.size();
        gzwrite(cache_file_gz, &count, sizeof(mi::Uint32));
        gzwrite(cache_file_gz, &mesh_normals[0], count*sizeof(mi::Float32)*3);

        count = mesh_tex_coords.size();
        gzwrite(cache_file_gz, &count, sizeof(mi::Uint32));
        gzwrite(cache_file_gz, &mesh_tex_coords[0], count*sizeof(mi::Float32)*2);

        count = mesh_vertex_colors.size();
        gzwrite(cache_file_gz, &count, sizeof(mi::Uint32));
        gzwrite(cache_file_gz, &mesh_vertex_colors[0], count*sizeof(mi::Float32)*4);

        count = mesh_pos_indices.size();
        gzwrite(cache_file_gz, &count, sizeof(mi::Uint32));
        gzwrite(cache_file_gz, &mesh_pos_indices[0], count*sizeof(mi::Uint32));

        count = mesh_normal_indices.size();
        gzwrite(cache_file_gz, &count, sizeof(mi::Uint32));
        gzwrite(cache_file_gz, &mesh_normal_indices[0], count*sizeof(mi::Uint32));

        count = mesh_tex_coord_indices.size();
        gzwrite(cache_file_gz, &count, sizeof(mi::Uint32));
        gzwrite(cache_file_gz, &mesh_tex_coord_indices[0], count*sizeof(mi::Uint32));

        count = mesh_vertex_color_indices.size();
        gzwrite(cache_file_gz, &count, sizeof(mi::Uint32));
        gzwrite(cache_file_gz, &mesh_vertex_color_indices[0], count*sizeof(mi::Uint32));

        count = mesh_vertex_colormap_indices.size();
        gzwrite(cache_file_gz, &count, sizeof(mi::Uint32));
        gzwrite(cache_file_gz, &mesh_vertex_colormap_indices[0], count*sizeof(mi::Uint32));

        count = mesh_materials.size();
        gzwrite(cache_file_gz, &count, sizeof(mi::Uint32));
        gzwrite(cache_file_gz, &mesh_materials[0], count*sizeof(mi::Uint16));

        count = triangle_flags.size();
        gzwrite(cache_file_gz, &count, sizeof(mi::Uint32));
        gzwrite(cache_file_gz, &triangle_flags[0], count*sizeof(nv::index::ITriangle_mesh_subset::Triflags));

        gzwrite(cache_file_gz, &bbox.min.x, sizeof(mi::Float32)*6);

        gzclose(cache_file_gz);

        mi::Float64 end_time_cache = nv::index_common::get_time();
        INFO_LOG<<(end_time_cache-start_time_cache)<<" seconds to write the cache file";
    }
    else
        INFO_LOG<<"Caching switched off, no cache files written";    
             
    // ---------------------------------------------------------------------------------------------------
    // Create submesh (must be derived from ITriangle_mesh_subset)
    // and returns mesh to the caller. The caller takes ownership!
    mi::base::Handle<nv::index::ITriangle_mesh_subset> grid_mesh(
        factory->create<nv::index::ITriangle_mesh_subset>());
    if(grid_mesh.is_valid_interface() && !mesh_vectors.empty())
    {
        // Generate global triangle IDs: 
        std::vector<mi::Uint64> global_triid;
        assert((mesh_pos_indices.size() % 3) == 0);
        const mi::Size nb_tris = mesh_pos_indices.size() / 3;
        global_triid.resize(nb_tris); 
        for(mi::Size i = 0; i < nb_tris; ++i){
            global_triid[i] = ~0UL; // set as invalid global triangle ID.
        }

        // check necessary data for creating a valid reservoir grid.
        if ((!mi::math::Bbox<mi::Float32, 3>(bbox).is_volume()) ||
            mesh_vectors.empty() || mesh_pos_indices.empty() || global_triid.empty()) {
            ERROR_LOG << "Cannot create a valid reservoir grid.";
            return 0;
        }

        const bool success = grid_mesh->initialize(
            bbox,
            &mesh_vectors[0],     mesh_vectors.size(),
            &mesh_pos_indices[0], mesh_pos_indices.size(),
            &global_triid[0],     global_triid.size(),

            mesh_normals.empty()                 ? 0 : &mesh_normals[0],                    mesh_normals.size(),
            mesh_tex_coords.empty()              ? 0 : &mesh_tex_coords[0],                 mesh_tex_coords.size(),
            mesh_vertex_colors.empty()           ? 0 : &mesh_vertex_colors[0],              mesh_vertex_colors.size(),

            mesh_normal_indices.empty()          ? 0 : &mesh_normal_indices[0],             
            mesh_tex_coord_indices.empty()       ? 0 : &mesh_tex_coord_indices[0],          
            mesh_vertex_color_indices.empty()    ? 0 : &mesh_vertex_color_indices[0],       
            mesh_vertex_colormap_indices.empty() ? 0 : &mesh_vertex_colormap_indices[0],
            
            mesh_materials.empty()               ? 0 : &mesh_materials[0],                  mesh_materials.size(),
            triangle_flags.empty()               ? 0 : &triangle_flags[0],                  triangle_flags.size());
        if(success)
        {
            grid_mesh->retain(); // since the handle will be out of scope.
            return grid_mesh.get();
        }
    }

    return NULL;
}

//----------------------------------------------------------------------
mi::base::Uuid Synthetic_reservoir_patch_generator::subset_id() const
{
    return nv::index::ITriangle_mesh_subset::IID();
}

//----------------------------------------------------------------------
const char* Synthetic_reservoir_patch_generator::get_configuration() const
{
    return m_configuration.c_str();
}

//----------------------------------------------------------------------
void Synthetic_reservoir_patch_generator::serialize(
    mi::neuraylib::ISerializer * serializer) const
{
    serializer->write(&m_layers, 1);
    serializer->write(&m_cells_per_row, 1);
    serializer->write(&m_cells_per_column, 1);
    serializer->write(&m_global_bbox.min.x, 6);
    
    mi::Uint32 nb_elements = mi::Uint32(m_cache_file.size());
    serializer->write(&nb_elements, 1);
    serializer->write(reinterpret_cast<const mi::Uint8*>(m_cache_file.c_str()), nb_elements);
    
    nb_elements = mi::Uint32(m_cache_type.size());
    serializer->write(&nb_elements, 1);
    serializer->write(reinterpret_cast<const mi::Uint8*>(m_cache_type.c_str()), nb_elements);
}

//----------------------------------------------------------------------
void Synthetic_reservoir_patch_generator::deserialize(
    mi::neuraylib::IDeserializer * deserializer)
{
    deserializer->read(&m_layers, 1);
    deserializer->read(&m_cells_per_row, 1);
    deserializer->read(&m_cells_per_column, 1);
    deserializer->read(&m_global_bbox.min.x, 6);
    
    mi::Uint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_cache_file.resize(nb_elements);
    deserializer->read(reinterpret_cast<mi::Uint8*>(&m_cache_file[0]), nb_elements);
    
    nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_cache_type.resize(nb_elements);
    deserializer->read(reinterpret_cast<mi::Uint8*>(&m_cache_type[0]), nb_elements);
}

//----------------------------------------------------------------------
}} // namespace nv::index_common
