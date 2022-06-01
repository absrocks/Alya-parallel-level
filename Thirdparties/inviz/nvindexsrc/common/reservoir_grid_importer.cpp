/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief external reservoir_grid data importer job


#include "reservoir_grid_importer.h"
#include <fstream>
#include <iterator>
#include <map>
#include <zlib.h>

#include <nv/index/itriangle_mesh_subset.h>
#include <nv/index/itriangle_mesh_scene_element.h>
#include <nv/index/iscene.h>
#include <nv/index/isession.h>
#include <nv/index/icolormap.h>

#include "common_utility.h"
#include "forwarding_logger.h"

namespace util {

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

/// Colormap generating utility function
void generate_colormap(std::vector<mi::math::Color_struct>& color_entries)
{
   color_entries.resize(256);

   mi::Float32 vmin = 0.f;
   mi::Float32 vmax = 0.9f;
   for(mi::Uint32 i=0; i<256; ++i)
   {
       mi::Float32 v = (i+0.5f)/256.f;

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

} // namespace util


namespace nv {
namespace index_common {

//----------------------------------------------------------------------
Reservoir_grid_importer::Reservoir_grid_importer(
        const std::string &                             geometry_file,
        const std::string &                             scalar_file,
        const mi::math::Bbox_struct<mi::Float32, 3>&    bbox,
        const std::string &                             cache_type)
    :
    m_geometry_file(geometry_file),
    m_scalar_file(scalar_file),
    m_global_bbox(bbox),
    m_cache_type(cache_type)
{
    std::ostringstream s;
    s << "type=reservoir_grid\n"
      << "importer=raw\n"
      << "geometry_file=" << geometry_file << "\n"
      << "scalar_file=" << scalar_file << "\n"
      << "cache_type=" << cache_type << "\n";
    m_configuration = s.str();
}

//----------------------------------------------------------------------
// constructor
Reservoir_grid_importer::Reservoir_grid_importer(
)
{
    // for serialization and copy only
}

//----------------------------------------------------------------------
// destructor
Reservoir_grid_importer::~Reservoir_grid_importer()
{
    // empty
}

//----------------------------------------------------------------------
const char* Reservoir_grid_importer::get_configuration() const
{
    return m_configuration.c_str();
}

namespace
{
//----------------------------------------------------------------------
// read geometry data from file
void read_reservoir_geometry(
    const std::string&                                       geometry_file,
    std::vector<mi::math::Vector_struct<mi::Float32, 3> >&   mesh_vectors,
    std::vector<mi::Uint64>&                                 cell_indices,
    std::map<mi::Uint32, mi::Uint64>&                        cells_per_part,
    const mi::math::Bbox_struct<mi::Sint32, 3>&              bounds)
{
    FILE * geometry_file_pointer = fopen(geometry_file.c_str(), "rb");

    if(geometry_file_pointer == NULL)
        ERROR_LOG<<"Reservoir geometry file not found";

    // C Binary
    char buffer[80];
    fread(buffer, sizeof(char), 80, geometry_file_pointer);

    // Description line 1
    fread(buffer, sizeof(char), 80, geometry_file_pointer);

    // Description line 2
    fread(buffer, sizeof(char), 80, geometry_file_pointer);

    // Node-id
    fread(buffer, sizeof(char), 80, geometry_file_pointer);

    // Element-id
    fread(buffer, sizeof(char), 80, geometry_file_pointer);

    mi::Uint64 offset = 0;
    while(!feof(geometry_file_pointer))
    {
        // Part
        fread(buffer, sizeof(char), 80, geometry_file_pointer);

        // Part number
        mi::Uint32 part_num = 0;
        fread(&part_num, sizeof(mi::Uint32), 1, geometry_file_pointer);

        if(part_num == 0)
            break;
        // Description line for the part
        fread(buffer, sizeof(char), 80, geometry_file_pointer);

        // Coordinates
        fread(buffer, sizeof(char), 80, geometry_file_pointer);

        // Total coordinates
        mi::Uint64 coordinate_size = 0;
        fread(&coordinate_size, sizeof(mi::Uint64), 1, geometry_file_pointer);

        if(coordinate_size == 0)
            WARN_LOG<<"Failed to read coordinate_size";

        std::vector<mi::Float32> xpos;
        std::vector<mi::Float32> ypos;
        std::vector<mi::Float32> zpos;

        xpos.resize(coordinate_size);
        fread(&xpos[0], sizeof(mi::Float32), coordinate_size, geometry_file_pointer);
        ypos.resize(coordinate_size);
        fread(&ypos[0], sizeof(mi::Float32), coordinate_size, geometry_file_pointer);
        zpos.resize(coordinate_size);
        fread(&zpos[0], sizeof(mi::Float32), coordinate_size, geometry_file_pointer);

        // Element type
        fread(buffer, sizeof(char), 80, geometry_file_pointer);

        // Total cells
        mi::Uint64 total_cells = 0;
        fread(&total_cells, sizeof(mi::Uint64), 1, geometry_file_pointer);

        if(total_cells == 0)
            WARN_LOG<<"Failed to read total_cells";

        cells_per_part[part_num] = total_cells;

        if(part_num > 1)
            offset += total_cells;

        // Read cell indices
        std::vector<mi::Uint64> file_cell_indices;
        file_cell_indices.resize(total_cells*8);
        fread(&file_cell_indices[0], sizeof(mi::Uint64), total_cells*8, geometry_file_pointer);

        for(mi::Uint64 i = 0; i < total_cells; ++i)
        {
            mi::math::Vector_struct<mi::Float32, 3> cell_pos_0;
            mi::Uint64 current_idx = file_cell_indices[i*8+0];
            cell_pos_0.x = xpos[current_idx];
            cell_pos_0.y = ypos[current_idx];
            cell_pos_0.z = zpos[current_idx];

            mi::math::Vector_struct<mi::Float32, 3> cell_pos_6;
            current_idx = file_cell_indices[i*8+6];
            cell_pos_6.x = xpos[current_idx];
            cell_pos_6.y = ypos[current_idx];
            cell_pos_6.z = zpos[current_idx];

            if(util::intersects_voxel(cell_pos_0, cell_pos_6, bounds))
            {
                mesh_vectors.push_back(cell_pos_0);
                current_idx = file_cell_indices[i*8+0];
                //cell_indices.push_back(offset+current_idx);

                mi::math::Vector_struct<mi::Float32, 3> cell_pos_1;
                current_idx = file_cell_indices[i*8+1];
                cell_pos_1.x = xpos[current_idx];
                cell_pos_1.y = ypos[current_idx];
                cell_pos_1.z = zpos[current_idx];
                mesh_vectors.push_back(cell_pos_1);
                //cell_indices.push_back(offset+current_idx);

                mi::math::Vector_struct<mi::Float32, 3> cell_pos_2;
                current_idx = file_cell_indices[i*8+2];
                cell_pos_2.x = xpos[current_idx];
                cell_pos_2.y = ypos[current_idx];
                cell_pos_2.z = zpos[current_idx];
                mesh_vectors.push_back(cell_pos_2);
                //cell_indices.push_back(offset+current_idx);

                mi::math::Vector_struct<mi::Float32, 3> cell_pos_3;
                current_idx = file_cell_indices[i*8+3];
                cell_pos_3.x = xpos[current_idx];
                cell_pos_3.y = ypos[current_idx];
                cell_pos_3.z = zpos[current_idx];
                mesh_vectors.push_back(cell_pos_3);
                //cell_indices.push_back(offset+current_idx);

                mi::math::Vector_struct<mi::Float32, 3> cell_pos_4;
                current_idx = file_cell_indices[i*8+4];
                cell_pos_4.x = xpos[current_idx];
                cell_pos_4.y = ypos[current_idx];
                cell_pos_4.z = zpos[current_idx];
                mesh_vectors.push_back(cell_pos_4);
                //cell_indices.push_back(offset+current_idx);

                mi::math::Vector_struct<mi::Float32, 3> cell_pos_5;
                current_idx = file_cell_indices[i*8+5];
                cell_pos_5.x = xpos[current_idx];
                cell_pos_5.y = ypos[current_idx];
                cell_pos_5.z = zpos[current_idx];
                mesh_vectors.push_back(cell_pos_5);
                //cell_indices.push_back(offset+current_idx);

                mesh_vectors.push_back(cell_pos_6);
                current_idx = file_cell_indices[i*8+6];
                //cell_indices.push_back(offset+current_idx);

                mi::math::Vector_struct<mi::Float32, 3> cell_pos_7;
                current_idx = file_cell_indices[i*8+7];
                cell_pos_7.x = xpos[current_idx];
                cell_pos_7.y = ypos[current_idx];
                cell_pos_7.z = zpos[current_idx];
                mesh_vectors.push_back(cell_pos_7);
                //cell_indices.push_back(offset+current_idx);

                cell_indices.push_back(offset+i);
            }
        }
    }
    fclose(geometry_file_pointer);
}

void read_reservoir_scalar(
    const std::string&                  scalar_file,
    std::map<mi::Uint32, mi::Uint64>&   cells_per_part,
    mi::Float32&                        scalar_min,
    mi::Float32&                        scalar_max,
    std::vector<mi::Float32>&           scalar_values)
{
    if(cells_per_part.empty())
        WARN_LOG<<"Cells per part information unavailable!";

    FILE * scalar_file_pointer = fopen(scalar_file.c_str(), "rb");

    if(scalar_file_pointer == NULL)
        ERROR_LOG<<"Scalar file not found";

    // Description line
    char buffer[80];
    fread(buffer, sizeof(char), 80, scalar_file_pointer);

    while(!feof(scalar_file_pointer))
    {
        // Part
        fread(buffer, sizeof(char), 80, scalar_file_pointer);

        // Part number
        mi::Uint32 part_num = 0;
        fread(&part_num, sizeof(mi::Uint32), 1, scalar_file_pointer);

        if(part_num == 0)
            break;

        //hexa8
        fread(buffer, sizeof(char), 80, scalar_file_pointer);

        mi::Uint64 part_cells = cells_per_part[part_num];
        std::vector<mi::Float32> part_scalar_values;

        part_scalar_values.resize(part_cells);

        if(fread(&part_scalar_values[0], sizeof(mi::Float32), part_cells, scalar_file_pointer) != part_cells)
            ERROR_LOG<<"Failed to read scalar values from the file";
        else
            std::copy(part_scalar_values.begin(), part_scalar_values.end(), std::back_inserter(scalar_values));

    }
    fclose(scalar_file_pointer);

    std::vector<mi::Float32>::iterator it = std::min_element(scalar_values.begin(), scalar_values.end());
    scalar_min = *it;
    it = std::max_element(scalar_values.begin(), scalar_values.end());
    scalar_max = *it;
}
} // namespace

//----------------------------------------------------------------------
mi::Size Reservoir_grid_importer::estimate(
    const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{
    // no estimate available since cells are generated on the fly per subcube
    return 0;
}

nv::index::IDistributed_data_subset* Reservoir_grid_importer::create(
    const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
    nv::index::IData_subset_factory*                factory,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{
    std::ostringstream cache_file_name_str;
    cache_file_name_str << m_geometry_file << "_"
                        << bounding_box.min.x << "." << bounding_box.max.x << "x"
                        << bounding_box.min.y << "." << bounding_box.max.y << "x"
                        << bounding_box.min.z << "." << bounding_box.max.z
                        << "_reservoir_raw_importer.tmp";

    std::string cache_file_name;
    cache_file_name = cache_file_name_str.str();

    std::ostringstream cache_file_name_str_gz;
    cache_file_name_str_gz << m_geometry_file << "_"
                        << bounding_box.min.x << "." << bounding_box.max.x << "x"
                        << bounding_box.min.y << "." << bounding_box.max.y << "x"
                        << bounding_box.min.z << "." << bounding_box.max.z
                        << "_reservoir_raw_importer.tmp.gz";

    std::string cache_file_name_gz;
    cache_file_name_gz = cache_file_name_str_gz.str();

    FILE * cache_file = fopen(cache_file_name.c_str(), "rb");
    gzFile cache_file_gz = gzopen(cache_file_name_gz.c_str(), "rb");

    std::vector<mi::math::Vector_struct<mi::Float32, 3> > mesh_vectors;
    std::vector<mi::math::Vector_struct<mi::Float32, 2> > mesh_tex_coords;
    std::vector<mi::math::Vector_struct<mi::Float32, 3> > mesh_normals;
    std::vector<mi::math::Color_struct > mesh_vertex_colors;
    std::vector<mi::Uint32> mesh_tex_coord_indices;
    std::vector<mi::Uint32> mesh_vertex_color_indices;
    std::vector<mi::Uint32> mesh_vertex_colormap_indices;
    std::vector<mi::Uint32> mesh_normal_indices;
    std::vector<mi::Uint32> mesh_pos_indices;

    // Materials
    std::vector<mi::Uint16> mesh_materials;
    mesh_materials.push_back(0);

    // triangle flags
    std::vector<nv::index::ITriangle_mesh_subset::Triflags> triangle_flags;

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
        mi::base::Handle<nv::index::ITriangle_mesh_subset>
            grid_mesh(factory->create<nv::index::ITriangle_mesh_subset>());
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
                grid_mesh->retain(); // since the handle will be out of scope next.
                return grid_mesh.get();
            }
        }

        return 0;
    }

    mi::Float64 start_time = nv::index_common::get_time();

    // Generate color values as a colormap
    // std::vector<mi::math::Color_struct> color_entries;
    // util::generate_colormap(color_entries);

    // if(color_entries.empty())
      // WARN_LOG<<"Colormap entries are empty!";

    // Get coordinates for the reservoir
    std::vector<mi::Uint64> cell_indices;
    std::map<mi::Uint32, mi::Uint64> cells_per_part;
    nv::index_common::read_reservoir_geometry(m_geometry_file, mesh_vectors, cell_indices, cells_per_part, bounding_box);

    if (mesh_vectors.empty())
        WARN_LOG<< "No coordinates indices were read!";

    // Get scalar values for the reservoir
    std::vector<mi::Float32> scalar_values;
    mi::Float32 scalar_min = 0;
    mi::Float32 scalar_max = 0;
    if(m_scalar_file != "")
        nv::index_common::read_reservoir_scalar(m_scalar_file, cells_per_part, scalar_max, scalar_min, scalar_values);

    if (scalar_values.empty())
        WARN_LOG<< "No scalar values read from the file!";

    // for(mi::Uint32 i = 0; i < scalar_values.size(); ++i)
        // INFO_LOG<<"scalar value: "<<scalar_values[i];

    mi::Uint32 total_cells = mesh_vectors.size()/8;
    mi::Uint32 total_quads = total_cells*6;

    std::vector<mi::Uint32> quad_indices;

    // create visible edge mask....for this case add all edges except the edge(2,0)
    nv::index::ITriangle_mesh_subset::Triflags visible_edges = (nv::index::ITriangle_mesh_subset::Triflags)
                                    (nv::index::ITriangle_mesh_subset::TRIFLAGS_VISIBLE_EDGE_01|
                                    nv::index::ITriangle_mesh_subset::TRIFLAGS_VISIBLE_EDGE_12);

    // Triangulate each cell
    for(mi::Uint32 i=0; i < total_cells; ++i)
    {
        // Indices for the triangles
        mesh_pos_indices.push_back(0+8*i); // 0
        mesh_pos_indices.push_back(1+8*i);
        mesh_pos_indices.push_back(5+8*i);

        mesh_pos_indices.push_back(5+8*i); // 1
        mesh_pos_indices.push_back(4+8*i);
        mesh_pos_indices.push_back(0+8*i);

        mesh_pos_indices.push_back(6+8*i); // 2
        mesh_pos_indices.push_back(2+8*i);
        mesh_pos_indices.push_back(3+8*i);

        mesh_pos_indices.push_back(3+8*i); // 3
        mesh_pos_indices.push_back(7+8*i);
        mesh_pos_indices.push_back(6+8*i);

        mesh_pos_indices.push_back(3+8*i); // 4
        mesh_pos_indices.push_back(0+8*i);
        mesh_pos_indices.push_back(4+8*i);

        mesh_pos_indices.push_back(4+8*i); // 5
        mesh_pos_indices.push_back(7+8*i);
        mesh_pos_indices.push_back(3+8*i);

        mesh_pos_indices.push_back(1+8*i); // 6
        mesh_pos_indices.push_back(2+8*i);
        mesh_pos_indices.push_back(6+8*i);

        mesh_pos_indices.push_back(6+8*i); // 7
        mesh_pos_indices.push_back(5+8*i);
        mesh_pos_indices.push_back(1+8*i);

        mesh_pos_indices.push_back(2+8*i); // 8
        mesh_pos_indices.push_back(1+8*i);
        mesh_pos_indices.push_back(0+8*i);

        mesh_pos_indices.push_back(0+8*i); // 9
        mesh_pos_indices.push_back(3+8*i);
        mesh_pos_indices.push_back(2+8*i);

        mesh_pos_indices.push_back(4+8*i); // 10
        mesh_pos_indices.push_back(5+8*i);
        mesh_pos_indices.push_back(6+8*i);

        mesh_pos_indices.push_back(6+8*i); // 11
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

        mi::Uint32 color_idx  = 0;
        if(!scalar_values.empty())
        {
            mi::Uint64 scalar_idx = cell_indices[i];
            mi::Float32 current_scalar = scalar_values[scalar_idx];
            color_idx = static_cast<mi::Uint32>(255.0f*(current_scalar - scalar_min) / (scalar_max - scalar_min));
            //INFO_LOG<<"scalar_idx is: "<<scalar_idx<<" current_scalar is: "<<current_scalar<<" color_idx: "<<color_idx;
        }
        else
        {
            color_idx = static_cast<mi::Uint32>(
                255.0f*(mesh_vectors[0+8*i].z- m_global_bbox.min.z) /
                (m_global_bbox.max.z - m_global_bbox.min.z));
        }

        // If no colormap, share same color across all the triangles of a cell
        //mesh_vertex_colors.push_back(color_entries[color_idx]);
        // for(mi::Uint32 j = 0; j<36; ++j)
            // mesh_vertex_color_indices.push_back(i);

        // Colormap indices for vertices
        for(mi::Uint32 j = 0; j<36; ++j)
            mesh_vertex_colormap_indices.push_back(color_idx);

        // Wireframe visible edges
        for(mi::Uint32 j = 0; j<12; ++j)
            triangle_flags.push_back(visible_edges);
    }

    // Calculate Face normal per quad
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

    //
    // Write cache file if not found
    //
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
    mi::base::Handle<nv::index::ITriangle_mesh_subset>
        grid_mesh(factory->create<nv::index::ITriangle_mesh_subset>());
    if(grid_mesh.is_valid_interface() && !mesh_vectors.empty())
    {
        // create global triangle ID based on cell_indices
        std::vector<mi::Uint64> global_triid;
        const mi::Size nb_cell = cell_indices.size();
        global_triid.resize(12L * nb_cell); // 12 triangle/cell
        for(mi::Size i = 0; i < nb_cell; ++i){
            for(mi::Size j = 0; j < 12; ++j){
                const mi::Size gid = i * 12 + j;
                global_triid[gid] = gid;
            }
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
            mi::Float64 end_time_file = nv::index_common::get_time();
            INFO_LOG << "Took " << (end_time_file-start_time)
                     << " seconds to read " << mesh_pos_indices.size()/36
                     << " cells (" << mesh_pos_indices.size()/3
                     << " triangles) " <<(nv::index_common::get_time()-end_time_file)
                     << " seconds to create the grid";
            grid_mesh->retain(); // since the handle will be out of scope.
            return grid_mesh.get();
        }
    }
    
    return 0;
}

//----------------------------------------------------------------------
mi::base::Uuid Reservoir_grid_importer::subset_id() const
{
    // reservoir grid importer generates ITriangle_mesh_subset
    return nv::index::ITriangle_mesh_subset::IID();
}

//----------------------------------------------------------------------
// serialize
void Reservoir_grid_importer::serialize(
    mi::neuraylib::ISerializer * serializer) const
{
    mi::Uint32 nb_elements = mi::Uint32(m_geometry_file.size());
    serializer->write(&nb_elements, 1);
    serializer->write(reinterpret_cast<const mi::Uint8*>(m_geometry_file.c_str()), nb_elements);

    nb_elements = mi::Uint32(m_scalar_file.size());
    serializer->write(&nb_elements, 1);
    serializer->write(reinterpret_cast<const mi::Uint8*>(m_scalar_file.c_str()), nb_elements);

    serializer->write(&m_global_bbox.min.x, 3);
    serializer->write(&m_global_bbox.max.x, 3);

    nb_elements = mi::Uint32(m_cache_type.size());
    serializer->write(&nb_elements, 1);
    serializer->write(reinterpret_cast<const mi::Uint8*>(m_cache_type.c_str()), nb_elements);
}

//----------------------------------------------------------------------
// deserialize
void Reservoir_grid_importer::deserialize(
    mi::neuraylib::IDeserializer * deserializer)
{
    mi::Uint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_geometry_file.resize(nb_elements);
    deserializer->read(reinterpret_cast<mi::Uint8*>(&m_geometry_file[0]), nb_elements);

    nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_scalar_file.resize(nb_elements);
    deserializer->read(reinterpret_cast<mi::Uint8*>(&m_scalar_file[0]), nb_elements);

    deserializer->read(&m_global_bbox.min.x, 3);
    deserializer->read(&m_global_bbox.max.x, 3);

    nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_cache_type.resize(nb_elements);
    deserializer->read(reinterpret_cast<mi::Uint8*>(&m_cache_type[0]), nb_elements);
}

//----------------------------------------------------------------------
}} // namespace nv::index_common
