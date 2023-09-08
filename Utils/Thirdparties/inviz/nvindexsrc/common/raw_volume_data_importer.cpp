/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Volume data importer reading a simple raw data format from file.

#include "forwarding_logger.h"

#include "raw_volume_data_importer.h"
#include "large_file_io.h"
#include "common_utility.h"

#include <nv/index/iregular_volume.h>
#include <nv/index/iregular_volume_brick.h>

#include <cstdio>

namespace nv {
namespace index_common {

mi::base::Lock Raw_volume_data_importer::s_file_access_lock;
mi::Uint64     Raw_volume_data_importer::s_total_read_bytes = 0;
mi::Float64    Raw_volume_data_importer::s_total_read_time  = 0.0;

//----------------------------------------------------------------------
Raw_volume_data_importer::Raw_volume_data_importer(
    const std::string&                            file_name,
    const mi::math::Vector_struct<mi::Uint32, 3>& size)
  : m_file_name(file_name),
    m_size(size),
    m_use_cache(true),
    m_cache_compression(0),
    m_serial_access(false),
    m_show_stats(false),
    m_demo_scale_hack_00(false)
{
    m_configuration = std::string() +
        "importer=raw\n" +
        "input_file=" + file_name + "\n";
}

Raw_volume_data_importer::Raw_volume_data_importer()
  : m_use_cache(true),
    m_cache_compression(0),
    m_serial_access(false),
    m_show_stats(false),
    m_demo_scale_hack_00(false)
{
}

//----------------------------------------------------------------------

namespace {

static inline void output_volume_import_stat(
    const std::string & header,
    mi::Uint64  volume_brick_size,
    mi::Float64 time_taken,
    mi::Uint64  total_read_bytes,
    mi::Float64 total_read_time)
{                        
    const mi::Float64 volume_brick_size_f64 = static_cast<mi::Float64>(volume_brick_size);
    const mi::Float64 total_read_bytes_f64  = static_cast<mi::Float64>(total_read_bytes);
    INFO_LOG << header     << (volume_brick_size_f64 / (1024.0 * 1024.0)) << " MB in " << time_taken << "s, "
             << "\ttotal " << (total_read_bytes_f64  / (1024.0 * 1024.0)) << " MB in " << total_read_time << "s, "
             << "\tavg "   << (total_read_bytes_f64  / total_read_time) / (1024.0 * 1024.0)  << " MB/s.";
}


static inline mi::Uint8 brick_value(
    const mi::math::Vector<mi::Uint32, 3>& brick_pos,
    const mi::Uint8*const                  brick_data,
    const mi::math::Vector<mi::Uint32, 3>& brick_dim)
{
    const mi::Size o =   static_cast<mi::Size>(brick_pos.z)
                       + static_cast<mi::Size>(brick_pos.y) * brick_dim.z
                       + static_cast<mi::Size>(brick_pos.x) * brick_dim.z * brick_dim.y;
    return brick_data[o];
}


} // namespace


mi::Size Raw_volume_data_importer::estimate(
    const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{
    // This is an 8-bit raw volume data importer:
    const mi::Size dx = bounding_box.max.x - bounding_box.min.x;
    const mi::Size dy = bounding_box.max.y - bounding_box.min.y;
    const mi::Size dz = bounding_box.max.z - bounding_box.min.z;
    const mi::Size volume_brick_size = dx * dy * dz;
    return volume_brick_size;
}

nv::index::IDistributed_data_subset* Raw_volume_data_importer::create(
    const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
    nv::index::IData_subset_factory*                factory,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{
    mi::base::Lock::Block block;
    const mi::math::Bbox<mi::Sint32, 3> bounds = bounding_box;

    if (m_serial_access)
    {
        block.set(&s_file_access_lock); // Acquire the lock
    }

    // The size for volumetric data brick that needs to be loaded
    mi::Sint64 dx = bounds.max.x - bounds.min.x;
    mi::Sint64 dy = bounds.max.y - bounds.min.y;
    mi::Sint64 dz = bounds.max.z - bounds.min.z;
    mi::Sint64 volume_brick_size = dx * dy * dz;

    INFO_LOG << "Raw volume data importer loads '" << m_file_name << "', bounds: " << bounds;

    mi::base::Handle<nv::index::IRegular_volume_brick_uint8> 
        volume_brick(factory->create<nv::index::IRegular_volume_brick_uint8>());
        
    if (!volume_brick.is_valid_interface())
    {
        ERROR_LOG << "Cannot create a volume brick.";
        return 0;
    }

    // allocate voxel data storage via volume_brick.
    mi::Uint8* voxel_data_storage = volume_brick->generate_voxel_storage(bounding_box);
    if (voxel_data_storage == 0)
    {
        ERROR_LOG << "Cannot generate voxel storage.";
        return 0;
    }

    // --- Demo hack parts >>>>>
    // We are assuming 508 subcube size with border 2
    // results in 1016 subcubes with border 4
    const mi::math::Vector<mi::Sint32, 3> demo_scale_factors   = mi::math::Vector<mi::Sint32, 3>(2);
    const mi::Size                        demo_scale_memory    =   demo_scale_factors.x
                                                                 * demo_scale_factors.y
                                                                 * demo_scale_factors.z;
    mi::math::Bbox<mi::Sint32, 3>         demo_scale_bounds_00 = bounds;

    if (m_demo_scale_hack_00)
    {
        demo_scale_bounds_00.min += mi::math::Vector<mi::Sint32, 3>(2) * demo_scale_factors;
        demo_scale_bounds_00.max -= mi::math::Vector<mi::Sint32, 3>(2) * demo_scale_factors;

        demo_scale_bounds_00.min /= demo_scale_factors;
        demo_scale_bounds_00.max /= demo_scale_factors;

        demo_scale_bounds_00.min -= mi::math::Vector<mi::Sint32, 3>(2);
        demo_scale_bounds_00.max += mi::math::Vector<mi::Sint32, 3>(2);

        INFO_LOG << "demo scale hack engaged... make it so number one!";
        INFO_LOG << "size:                 " << m_size;
        INFO_LOG << "bounds:               " << bounds;
        INFO_LOG << "demo_scale_bounds_00: " << demo_scale_bounds_00;
    }
    // <<<<< Demo hack parts ---

    // ---------------------------------------------------------------------------------------------------
    // Try importing the data from a temporary cache file possibly generated in a previous import process.
    std::string cache_file_name;
    io::File cache_file(m_cache_compression);
    bool deprecated_format = false;

    // First try the old naming scheme for compatibility
    if (m_use_cache)
    {
        // In the old naming scheme the bounding box did not contain the border
        mi::math::Bbox<mi::Sint32, 3> old_bounds = bounds;
        old_bounds.min += mi::math::Vector<mi::Sint32, 3>(1);
        old_bounds.max -= mi::math::Vector<mi::Sint32, 3>(1);

        std::ostringstream cache_file_name_str;
        cache_file_name_str << m_file_name << "_"
                            << old_bounds.min.x << "." << old_bounds.max.x << "x"
                            << old_bounds.min.y << "." << old_bounds.max.y << "x"
                            << old_bounds.min.z << "." << old_bounds.max.z
                            << "_external_xyz_importer.tmp";
       
        cache_file_name = cache_file_name_str.str();
        cache_file.open(cache_file_name.c_str(), std::ios_base::in);

        if (cache_file)
            deprecated_format = true;
    }

    if (m_use_cache && !cache_file)
    {    
        std::ostringstream cache_file_name_str;
        if (m_demo_scale_hack_00)
        {
            cache_file_name_str << m_file_name << "_"
                                << demo_scale_bounds_00.min.x << "." << demo_scale_bounds_00.max.x << "x"
                                << demo_scale_bounds_00.min.y << "." << demo_scale_bounds_00.max.y << "x"
                                << demo_scale_bounds_00.min.z << "." << demo_scale_bounds_00.max.z
                                << "_raw_importer.tmp";
        }
        else
        {
            cache_file_name_str << m_file_name << "_"
                                << bounds.min.x << "." << bounds.max.x << "x"
                                << bounds.min.y << "." << bounds.max.y << "x"
                                << bounds.min.z << "." << bounds.max.z
                                << "_raw_importer.tmp";
        }

        cache_file_name = cache_file_name_str.str();
        cache_file.open(cache_file_name.c_str(), std::ios_base::in);
    }

    const mi::Float64 time_start = nv::index_common::get_time();

    if (!cache_file)
    {
        if (m_use_cache)
            INFO_LOG << "Cache file '" << cache_file_name << "' not yet available.";
    }
    else
    {
        INFO_LOG << "Using cache file '" << cache_file_name << "'"
                 << (m_demo_scale_hack_00 ? " (scaled brick data)" : "")
                 << (deprecated_format ? " (deprecated naming scheme, consider deleting and regenerating the cache files)" : "");

        if (m_demo_scale_hack_00) {
            mi::Uint8* demo_scale_src_data
                = new mi::Uint8[volume_brick_size / demo_scale_memory]; // we only read an 8th of the size

            // Fetch data for the entire brick at once
            if (   cache_file.read(demo_scale_src_data, 0ll, volume_brick_size / demo_scale_memory)
                != volume_brick_size / demo_scale_memory)
            {
                ERROR_LOG << "Raw_volume_importer::create(): "
                          << "error reading from file '" << cache_file_name << "', "
                          << "unable to read " << volume_brick_size / demo_scale_memory << " bytes.";
            }
            cache_file.close();


            if (m_serial_access)
            {
                block.release(); // disk access is done
            }

            using namespace mi::math;

            using mi::math::clamp;
            using mi::math::floor;
            using mi::math::ceil;

            typedef mi::math::Vector<mi::Uint32,  3> Vec3ui;
            typedef mi::math::Vector<mi::Float32, 3> Vec3f;
            typedef mi::math::Vector<mi::Float64, 3> Vec3d;

            const Vec3ui src_size = Vec3ui(  demo_scale_bounds_00.max
                                           - demo_scale_bounds_00.min);
            const Vec3ui dst_size = Vec3ui(  bounds.max - bounds.min);

            for (mi::Uint32 x = 0; x < dst_size.x; ++x)
            {
                for (mi::Uint32 y = 0; y < dst_size.y; ++y)
                {
                    for (mi::Uint32 z = 0; z < dst_size.z; ++z)
                    {
                        const Vec3d dst_p = Vec3d(x, y, z);
                        const Vec3d src_p = dst_p / Vec3d(demo_scale_factors);

                        const Vec3ui f = Vec3ui(clamp(static_cast<mi::Uint32>(floor(src_p.x)), 0u, src_size.x - 1),
                                                clamp(static_cast<mi::Uint32>(floor(src_p.y)), 0u, src_size.y - 1),
                                                clamp(static_cast<mi::Uint32>(floor(src_p.z)), 0u, src_size.z - 1));
                        const Vec3ui c = Vec3ui(clamp(static_cast<mi::Uint32>(ceil(src_p.x)), 0u, src_size.x - 1),
                                                clamp(static_cast<mi::Uint32>(ceil(src_p.y)), 0u, src_size.y - 1),
                                                clamp(static_cast<mi::Uint32>(ceil(src_p.z)), 0u, src_size.z - 1));
                        const Vec3f  w = Vec3f(src_p - Vec3d(f));

                        const mi::Float32 v000 = brick_value(Vec3ui(f.x, f.y, f.z), demo_scale_src_data, src_size);
                        const mi::Float32 v001 = brick_value(Vec3ui(f.x, f.y, c.z), demo_scale_src_data, src_size);
                        const mi::Float32 v010 = brick_value(Vec3ui(f.x, c.y, f.z), demo_scale_src_data, src_size);
                        const mi::Float32 v011 = brick_value(Vec3ui(f.x, c.y, c.z), demo_scale_src_data, src_size);
                        const mi::Float32 v100 = brick_value(Vec3ui(c.x, f.y, f.z), demo_scale_src_data, src_size);
                        const mi::Float32 v101 = brick_value(Vec3ui(c.x, f.y, c.z), demo_scale_src_data, src_size);
                        const mi::Float32 v110 = brick_value(Vec3ui(c.x, c.y, f.z), demo_scale_src_data, src_size);
                        const mi::Float32 v111 = brick_value(Vec3ui(c.x, c.y, c.z), demo_scale_src_data, src_size);
                        
                        const mi::Float32 x00 = lerp(v000, v100, w.x);
                        const mi::Float32 x01 = lerp(v001, v101, w.x);
                        const mi::Float32 x10 = lerp(v010, v110, w.x);
                        const mi::Float32 x11 = lerp(v011, v111, w.x);

                        const mi::Float32 y0 = lerp(x00, x10, w.y);
                        const mi::Float32 y1 = lerp(x01, x11, w.y);

                        const mi::Float32 v_int = lerp(y0, y1, w.z);

                        const mi::Uint8 oval = static_cast<mi::Uint8>(clamp(v_int, 0.0f, 255.0f));
                        //const mi::Uint8 oval = static_cast<mi::Uint8>(clamp(v000, 0.0f, 255.0f));

                        const mi::Size dst_o =   static_cast<mi::Size>(z)
                                               + static_cast<mi::Size>(y) * dst_size.z
                                               + static_cast<mi::Size>(x) * dst_size.z * dst_size.y;

                        voxel_data_storage[dst_o] = oval;
                    }
                }
            }
            delete [] demo_scale_src_data;
        }
        else
        {
            // Fetch data for the entire brick at once
            if (cache_file.read(voxel_data_storage, 0ll, volume_brick_size) != static_cast<mi::Uint64>(volume_brick_size)) {
                ERROR_LOG << "Raw_volume_importer::create(): "
                          << "error reading from file '" << cache_file_name << "', "
                          << "unable to read " << volume_brick_size << " bytes.";
            }
            cache_file.close();
        }

        mi::Float64 time_end = nv::index_common::get_time();
        mi::Float64 time_taken = time_end - time_start;
        s_total_read_time += time_taken;
        s_total_read_bytes += volume_brick_size;
        if (m_show_stats)
        {
            output_volume_import_stat("Volume (cached): Read ", volume_brick_size, time_taken,
                                      s_total_read_bytes, s_total_read_time);
        }

        volume_brick->retain();
        return volume_brick.get();
    }

    // -------------------------------------------------------------------------------------------
    // Otherwise, fetch the data from the original dataset file.
    io::File in_file(m_file_name.c_str(), std::ios_base::in);
    if (!in_file)
    {
        // If the file is not available for reading then create a subset that
        // only contains invalid values mi::Uint8(0).
        ERROR_LOG << "Cannot open volume data file '"  << m_file_name << "'. "
                  << "Will create a database element with invalid value '0' (type: mi::Uint8).";
        return NULL;
    }

    mi::Sint32 x_dst = 0;
    for (mi::Sint32 x = bounds.min.x; x < bounds.max.x; ++x)
    {
        if (x < 0 || x >= (int)m_size.x)
        {
            // Skip the border voxels that are outside of the dataset, they will be filled with
            // duplicated boundary data later.
            x_dst++;
            continue;
        }

        mi::Sint32 y_dst = 0;
        for (mi::Sint32 y = bounds.min.y; y < bounds.max.y; ++y)
        {
            if (y < 0 || y >= (int)m_size.y)
            {
                // Skip the border voxels that are outside of the dataset, they will be filled with
                // duplicated boundary data later.
                y_dst++;
                continue;
            }

            // Data offsets, be sure to use 64-bit integers here.
            mi::Uint64 memory_offset
                = (mi::Uint64)dz * dy * x_dst + dz * y_dst;
            mi::Uint64 file_offset
                = (mi::Uint64)m_size.z * m_size.y * x + m_size.z * y + bounds.min.z;

            mi::Uint64 bytes_to_read = bounds.max.z - bounds.min.z;

            // Skip the minimum z border voxels that are outside of the dataset, they will be filled
            // with duplicated boundary data later.
            mi::Sint32 skip_min = -std::min(bounds.min.z, 0);
            memory_offset += skip_min;
            file_offset += skip_min;
            bytes_to_read -= skip_min;

            // Skip the maximum z border voxels that are outside of the dataset, they will be filled
            // with duplicated boundary data later.
            mi::Sint32 skip_max = std::max(bounds.max.z - (int)m_size.z, 0);
            bytes_to_read -= skip_max;

            if (bytes_to_read != in_file.read(&voxel_data_storage[memory_offset], file_offset, bytes_to_read))
            {
                ERROR_LOG << "Raw_volume_importer::create(): "
                          << "Error reading from file '" << m_file_name << "', "
                          << "unable to read " << volume_brick_size
                          << " bytes from offset " << file_offset << ".";
            }

            y_dst++;
        }
        x_dst++;
    }
    in_file.close();

    // --------------------------------------------------------------------------------------
    // Duplicate voxels at the dataset boundary (clamping) to fill up the border if necessary
    //

    // Calculate position of the voxels at the boundary of the current brick in memory. These are
    // the voxels that will be duplicated towards the border.
    const mi::math::Bbox<mi::Sint32, 3> boundary_dst(
        -bounds.min,
         mi::math::Vector<mi::Sint32, 3>(m_size) - bounds.min - mi::math::Vector<mi::Sint32, 3>(1));

    // Duplicate voxels at the minimum z boundary of the dataset (if that is part of the current
    // brick, i.e. the min component of the bounding box is zero)
    for (mi::Sint64 z = 0; z < -bounds.min.z; ++z)
    {
        for (mi::Sint64 x=0; x < dx; ++x)
        {
            for (mi::Sint64 y=0; y < dy; ++y)
            {
                voxel_data_storage[dz * dy * x + dz * y + z]
                    = voxel_data_storage[dz * dy * x + dz * y + boundary_dst.min.z];
            }
        }
    }

    // Duplicate voxels at the maximum z boundary of the dataset (if that is part of the current
    // brick, i.e. the max component of the bounding box is larger than the dataset)
    for (mi::Sint64 z = m_size.z; z < bounds.max.z; ++z)
    {
        for (mi::Sint64 x=0; x < dx; ++x)
        {
            for (mi::Sint64 y=0; y < dy; ++y)
            {
                voxel_data_storage[dz * dy * x + dz * y + (z - bounds.min.z)]
                    = voxel_data_storage[dz * dy * x + dz * y + boundary_dst.max.z];
            }
        }
    }

    // Duplicate min y
    for (mi::Sint64 y = 0; y < -bounds.min.y; ++y)
    {
        for (mi::Sint64 x=0; x < dx; ++x)
        {
            for (mi::Sint64 z=0; z < dz; ++z)
            {
                voxel_data_storage[dz * dy * x + (dz * y) + z]
                    = voxel_data_storage[dz * dy * x + (dz * boundary_dst.min.y) + z];
            }
        }
    }

    // Duplicate max y
    for (mi::Sint64 y = m_size.y; y < bounds.max.y; ++y)
    {
        for (mi::Sint64 x=0; x < dx; ++x)
        {
            for (mi::Sint64 z=0; z < dz; ++z)
            {
                voxel_data_storage[dz * dy * x + dz * (y - bounds.min.y) + z]
                    = voxel_data_storage[dz * dy * x + dz * boundary_dst.max.y + z];
            }
        }
    }

    // Duplicate min x
    for (mi::Sint64 x = 0; x < -bounds.min.x; ++x)
    {
        memcpy(&voxel_data_storage[dz * dy * x], &voxel_data_storage[dz * dy * boundary_dst.min.x], dz * dy);
    }

    // Duplicate max x
    for (mi::Sint64 x = m_size.x; x < bounds.max.x; ++x)
    {
        memcpy(&voxel_data_storage[dz * dy * (x - bounds.min.x)], &voxel_data_storage[dz * dy * boundary_dst.max.x], dz * dy);
    }

    mi::Float64 time_end = nv::index_common::get_time();
    mi::Float64 time_taken = time_end - time_start;
    s_total_read_time += time_taken;
    s_total_read_bytes += volume_brick_size;
    if (m_show_stats)
    {
        output_volume_import_stat("Volume: Read ", volume_brick_size, time_taken,
                                  s_total_read_bytes, s_total_read_time);
    }

    // Save data to cache file to accelerate the loading in the future.
    if (m_use_cache)
    {
        cache_file.open(cache_file_name, std::ios_base::out);
        if (!cache_file)
        {
            INFO_LOG << "Cache file '" << cache_file_name << "' cannot be opened for writing.";
        }
        else
        {
            INFO_LOG << "Writing volume data to cache file '" << cache_file_name << "'";
            if (cache_file.write(voxel_data_storage, 0ll, volume_brick_size) != static_cast<mi::Uint64>(volume_brick_size))
            {
                ERROR_LOG << "Raw_volume_importer::create(): "
                          << "Error writing to file '" << cache_file_name << "', "
                          << "unable to write " << volume_brick_size << " bytes.";
            }
            cache_file.close();
        }
    }

    volume_brick->retain();
    return volume_brick.get();
}

//----------------------------------------------------------------------
mi::base::Uuid Raw_volume_data_importer::subset_id() const
{
    // currently generate uint8 volume only
    return nv::index::IRegular_volume_brick_uint8::IID();
}

//----------------------------------------------------------------------
void Raw_volume_data_importer::serialize(
    mi::neuraylib::ISerializer* serializer) const
{
    mi::Uint32 nb_elements = mi::Uint32(m_file_name.size());
    serializer->write(&nb_elements, 1);
    serializer->write(reinterpret_cast<const mi::Uint8*>(m_file_name.c_str()), nb_elements);

    serializer->write(&m_size.x, 3);
    serializer->write(&m_use_cache, 1);
    serializer->write(&m_cache_compression, 1);
    serializer->write(&m_serial_access, 1);
    serializer->write(&m_show_stats, 1);
    serializer->write(&m_demo_scale_hack_00, 1);

    nb_elements = mi::Uint32(m_configuration.size());
    serializer->write(&nb_elements, 1);
    serializer->write(reinterpret_cast<const mi::Uint8*>(m_configuration.c_str()), nb_elements);
}

//----------------------------------------------------------------------
void Raw_volume_data_importer::deserialize(
    mi::neuraylib::IDeserializer* deserializer)
{
    mi::Uint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_file_name.resize(nb_elements);
    deserializer->read(reinterpret_cast<mi::Uint8*>(&m_file_name[0]), nb_elements);

    deserializer->read(&m_size.x, 3);
    deserializer->read(&m_use_cache, 1);
    deserializer->read(&m_cache_compression, 1);
    deserializer->read(&m_serial_access, 1);
    deserializer->read(&m_show_stats, 1);
    deserializer->read(&m_demo_scale_hack_00, 1);

    deserializer->read(&nb_elements, 1);
    m_configuration.resize(nb_elements);
    deserializer->read(reinterpret_cast<mi::Uint8*>(&m_configuration[0]), nb_elements);
}

//----------------------------------------------------------------------
const char* Raw_volume_data_importer::get_configuration() const
{
    return m_configuration.c_str();
}

}} // namespace nv::index_common
