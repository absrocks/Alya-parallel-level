/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Volume data importer reading a simple raw data format from file.

#include <cstdio>
#include <map>

#include <nv/index/iregular_volume.h>
#include <mi/base/lock.h>

#include "common_utility.h"
#include "forwarding_logger.h"
#include "large_file_io.h"

#include "short_value_raw_volume_data_importer.h"

namespace nv {
namespace index_common {

mi::base::Lock Short_value_raw_volume_data_importer::s_file_access_lock;
mi::Uint64     Short_value_raw_volume_data_importer::s_total_read_bytes = 0;
mi::Float64    Short_value_raw_volume_data_importer::s_total_read_time  = 0.0;

//----------------------------------------------------------------------
Short_value_raw_volume_data_importer::Short_value_raw_volume_data_importer(
    const std::string&                            file_name,
    const std::string&                            src_folder,
    const std::string&                            dst_folder,
    const mi::math::Vector_struct<mi::Uint32, 3>& size)
  : m_file_name(file_name),
    m_src_folder(src_folder),
    m_dst_folder(dst_folder),
    m_size(size),
    m_use_cache(false),
    m_cache_compression(0),
    m_serial_access(false),
    m_show_stats(false)
{
    m_configuration = std::string() +
        "importer=short\n" +
        "input_file=" + file_name + "\n";
}

Short_value_raw_volume_data_importer::Short_value_raw_volume_data_importer()
  : m_use_cache(false),
    m_cache_compression(0),
    m_serial_access(false),
    m_show_stats(false)
{
}

Short_value_raw_volume_data_importer::~Short_value_raw_volume_data_importer()
{
}
//----------------------------------------------------------------------

namespace {

static inline void output_volume_import_stat(
    const std::string&      header,
    mi::Uint64              volume_brick_size,
    mi::Float64             time_taken,
    mi::Uint64              total_read_bytes,
    mi::Float64             total_read_time)
{                        
    const mi::Float64 volume_brick_size_f64 = static_cast<mi::Float64>(volume_brick_size);
    const mi::Float64 total_read_bytes_f64  = static_cast<mi::Float64>(total_read_bytes);
    INFO_LOG << header     << (volume_brick_size_f64 / (1024.0 * 1024.0)) << " MB in " << time_taken << "s, "
             << "\ttotal " << (total_read_bytes_f64  / (1024.0 * 1024.0)) << " MB in " << total_read_time << "s, "
             << "\tavg "   << (total_read_bytes_f64  / total_read_time) / (1024.0 * 1024.0)  << " MB/s.";
}

} // namespace


void Short_value_raw_volume_data_importer::add_time_step(const std::string& file)
{
    m_time_steps.push_back(file);
}

mi::Size Short_value_raw_volume_data_importer::estimate(
    const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{
    // This is an 32-bit float volume data importer:
    const mi::Size dx = bounding_box.max.x - bounding_box.min.x;
    const mi::Size dy = bounding_box.max.y - bounding_box.min.y;
    const mi::Size dz = bounding_box.max.z - bounding_box.min.z;
    const mi::Size volume_brick_size = dx * dy * dz * sizeof(mi::Uint16);
    return volume_brick_size;
}

nv::index::IDistributed_data_subset* Short_value_raw_volume_data_importer::create(
    const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
    nv::index::IData_subset_factory*                factory,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{
    const std::string& file_name = m_file_name;

    return create(file_name, bounding_box, factory, dice_transaction);
}

nv::index::IDistributed_data_subset* Short_value_raw_volume_data_importer::create(
    const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
    mi::Uint32                                      time_step,
    nv::index::IData_subset_factory*                factory,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{
    const std::string& file_name = m_time_steps[time_step%m_time_steps.size()];

    return create(file_name, bounding_box, factory, dice_transaction);
}


nv::index::IDistributed_data_subset* Short_value_raw_volume_data_importer::create(
    const std::string&                              file_name,
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

    std::string src_filename = file_name;
    if (!m_src_folder.empty())
    {
        std::ostringstream file_name_str;
        file_name_str << m_src_folder << "/" << file_name;
        src_filename = file_name_str.str();
    }

    // The size for volumetric data brick that needs to be loaded
    mi::Sint64 dx = bounds.max.x - bounds.min.x;
    mi::Sint64 dy = bounds.max.y - bounds.min.y;
    mi::Sint64 dz = bounds.max.z - bounds.min.z;
    mi::Sint64 volume_brick_size = dx * dy * dz;

    INFO_LOG << "16-bit volume data importer loads '" << src_filename << "', bounds: " << bounds;

    mi::base::Handle<nv::index::IRegular_volume_brick_uint16> volume_brick(factory->create<nv::index::IRegular_volume_brick_uint16>());
        
    if (!volume_brick.is_valid_interface())
    {
        ERROR_LOG << "Cannot create a volume brick.";
        return 0;
    }

    // allocate voxel data storage via volume_brick.
    mi::Uint16* voxel_data_storage = volume_brick->generate_voxel_storage(bounding_box);
    if (voxel_data_storage == 0)
    {
        ERROR_LOG << "Cannot generate a short typed voxel storage.";
        return 0;
    }

    // ---------------------------------------------------------------------------------------------------
    // Try importing the data from a temporary cache file possibly generated in a previous import process.
    std::string cache_file_name;
    io::File cache_file(m_cache_compression);

    if (m_use_cache)
    {    
        std::ostringstream cache_file_name_str;
        {
            cache_file_name_str << m_dst_folder << "/" << file_name << "_"
                                << bounds.min.x << "." << bounds.max.x << "x"
                                << bounds.min.y << "." << bounds.max.y << "x"
                                << bounds.min.z << "." << bounds.max.z
                                << "_16-bit_cache" << ".tmp";
        }

        cache_file_name = cache_file_name_str.str();
        cache_file.open(cache_file_name.c_str(), std::ios_base::in);
    }

    const mi::Float64 time_start = nv::index_common::get_time();

    if (!cache_file)
    {
        if (m_use_cache)
        {
            INFO_LOG << "Cache file '" << cache_file_name << "' not yet available.";
        }
    }
    else
    {
        INFO_LOG << "Using cache file '" << cache_file_name << "'";
        {
            const mi::Uint64 bytes_to_read = static_cast<mi::Uint64>(volume_brick_size) * sizeof(mi::Uint16);
            const mi::Uint64 bytes_read = cache_file.read(voxel_data_storage, 0ll, bytes_to_read);
            // Fetch data for the entire brick at once
            if (bytes_read != bytes_to_read)
            {
                ERROR_LOG << "Short_value_raw_volume_data_importer::create(): "
                          << "error reading from file '" << cache_file_name << "', "
                          << "unable to read " << bytes_to_read << " bytes, got only " << bytes_read << ".";
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
    io::File in_file(src_filename.c_str(), std::ios_base::in);
    if (!in_file)
    {
        // If the file is not available for reading then create a subset that
        // only contains invalid values mi::Uint8(0).
        ERROR_LOG << "Cannot open volume data file '"  << src_filename << "'. "
                  << "Will create a database element with invalid value '0' (type: mi::Uint16).";
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

            mi::Uint64 bytes_to_read = (bounds.max.z - bounds.min.z);

            // Skip the minimum z border voxels that are outside of the dataset, they will be filled
            // with duplicated boundary data later.
            mi::Sint32 skip_min = (-std::min(bounds.min.z, 0));
            memory_offset += skip_min;
            file_offset += skip_min;
            bytes_to_read -= skip_min;

            // Skip the maximum z border voxels that are outside of the dataset, they will be filled
            // with duplicated boundary data later.
            mi::Sint32 skip_max = (std::max(bounds.max.z - (int)m_size.z, 0));
            bytes_to_read -= skip_max;

            if ((bytes_to_read*sizeof(mi::Uint16)) != in_file.read(&voxel_data_storage[memory_offset],
                                                    file_offset*sizeof(mi::Uint16),
                                                    bytes_to_read*sizeof(mi::Uint16)))
            {
                ERROR_LOG << "Short_value_raw_volume_data_importer::create(): "
                          << "Error reading from file '" << src_filename << "', "
                          << "unable to read " << (bytes_to_read*sizeof(mi::Uint16))
                          << " bytes from offset " << (file_offset*sizeof(mi::Uint16)) << ". Brick: " << bounds;
                return NULL;
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
            {
                INFO_LOG << "Writing volume values to cache file '" << cache_file_name << "'";
                const mi::Uint64 bytes_to_write = static_cast<mi::Uint64>(volume_brick_size) * sizeof(mi::Uint16);
                const mi::Uint64 bytes_written = cache_file.write(voxel_data_storage, 0ll, bytes_to_write);
                if (bytes_written != bytes_to_write)
                {
                    ERROR_LOG << "Short_value_raw_volume_data_importer::create(): "
                              << "Error writing to file '" << cache_file_name << "', "
                              << "unable to write " << bytes_to_write << " bytes, could only write " << bytes_written << ".";
                }
            }
            cache_file.close();
        }
    }

    volume_brick->retain();
    return volume_brick.get();
}

namespace {
void serialize_string(
    mi::neuraylib::ISerializer* serializer,
    const std::string&          the_string)
{
    mi::Uint32 nb_elements = mi::Uint32(the_string.size());
    serializer->write(&nb_elements, 1);
    serializer->write(reinterpret_cast<const mi::Uint8*>(the_string.c_str()), nb_elements);

}

void deserialize_string(
    mi::neuraylib::IDeserializer* deserializer,
    std::string&                  the_string)
{
    mi::Uint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    the_string.resize(nb_elements);
    deserializer->read(reinterpret_cast<mi::Uint8*>(&the_string[0]), nb_elements);
}
}

//----------------------------------------------------------------------
void Short_value_raw_volume_data_importer::serialize(
    mi::neuraylib::ISerializer* serializer) const
{
    // mi::Uint32 nb_elements = mi::Uint32(m_file_name.size());
    // serializer->write(&nb_elements, 1);
    // serializer->write(reinterpret_cast<const mi::Uint8*>(m_file_name.c_str()), nb_elements);
    serialize_string(serializer, m_file_name);
    serialize_string(serializer, m_src_folder);
    serialize_string(serializer, m_dst_folder);

    serializer->write(&m_size.x, 3);
    serializer->write(&m_range.x, 2);
    const mi::Uint32 nb_elements = m_time_steps.size();
    serializer->write(&nb_elements, 1);
    for(mi::Uint32 i=0; i<nb_elements; ++i)
    {
        serialize_string(serializer, m_time_steps[i]);
    }

    serializer->write(&m_use_cache, 1);
    serializer->write(&m_cache_compression, 1);
    serializer->write(&m_serial_access, 1);
    serializer->write(&m_show_stats, 1);

    //nb_elements = mi::Uint32(m_configuration.size());
    //serializer->write(&nb_elements, 1);
    //serializer->write(reinterpret_cast<const mi::Uint8*>(m_configuration.c_str()), nb_elements);
    serialize_string(serializer, m_configuration);
}

//----------------------------------------------------------------------
void Short_value_raw_volume_data_importer::deserialize(
    mi::neuraylib::IDeserializer* deserializer)
{
    // mi::Uint32 nb_elements = 0;
    // deserializer->read(&nb_elements, 1);
    // m_file_name.resize(nb_elements);
    // deserializer->read(reinterpret_cast<mi::Uint8*>(&m_file_name[0]), nb_elements);
    deserialize_string(deserializer, m_file_name);
    deserialize_string(deserializer, m_src_folder);
    deserialize_string(deserializer, m_dst_folder);

    deserializer->read(&m_size.x, 3);
    deserializer->read(&m_range.x, 2);
    mi::Uint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    for(mi::Uint32 i=0; i<nb_elements; ++i)
    {
        std::string brick_name;
        deserialize_string(deserializer, brick_name);
        m_time_steps.push_back(brick_name);
    }

    deserializer->read(&m_use_cache, 1);
    deserializer->read(&m_cache_compression, 1);
    deserializer->read(&m_serial_access, 1);
    deserializer->read(&m_show_stats, 1);

    // deserializer->read(&nb_elements, 1);
    // m_configuration.resize(nb_elements);
    // deserializer->read(reinterpret_cast<mi::Uint8*>(&m_configuration[0]), nb_elements);
    deserialize_string(deserializer, m_configuration);
}

//----------------------------------------------------------------------
const char* Short_value_raw_volume_data_importer::get_configuration() const
{
    return m_configuration.c_str();
}

}} // namespace nv::index_common
