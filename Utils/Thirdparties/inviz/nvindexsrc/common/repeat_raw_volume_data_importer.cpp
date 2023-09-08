/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Volume data importer reading a simple raw data format from file.

#include "forwarding_logger.h"

#include "repeat_raw_volume_data_importer.h"
#include "large_file_io.h"
#include "common_utility.h"

#include <nv/index/iregular_volume.h>
#include <nv/index/iregular_volume_brick.h>

#include <cstdio>

namespace nv {
namespace index_common {

//----------------------------------------------------------------------
// static member definition
mi::Uint8*          Repeat_raw_volume_data_importer::m_source = 0;
mi::base::Lock      Repeat_raw_volume_data_importer::m_source_lock;
mi::base::Atom32    Repeat_raw_volume_data_importer::m_cache_uses;

//----------------------------------------------------------------------
Repeat_raw_volume_data_importer::Repeat_raw_volume_data_importer(
    const std::string&                     file_name,
    const mi::math::Vector<mi::Uint32, 3>& source_size,
    const mi::math::Vector<mi::Uint32, 3>& repeat,
    bool                                   cache_source_data)
    : m_file_name(file_name),
      m_source_size(source_size),
      m_repeat(repeat),
      m_cache_source_data(cache_source_data)
{
    m_size = mi::math::Vector<mi::Uint32, 3>(m_source_size) * mi::math::Vector<mi::Uint32, 3>(m_repeat);

    if (m_cache_source_data && (m_source != 0))
    {
        ERROR_LOG << "Trying to create a new Repeat_raw_volume_data_importer with source caching "
                  << "while another one is still active. Disabling source caching.";
        m_cache_source_data = false;
    }

    std::ostringstream s;
    s << "importer=repeat\n"
      << "input_file=" << m_file_name << "\n"
      << "source size=" << m_source_size.x << " " << m_source_size.y << " " << m_source_size.z << "\n"
      << "repeat=" << m_repeat.x << " " << m_repeat.y << " " << m_repeat.z << "\n";
    m_configuration = s.str();

    DEBUG_LOG << m_configuration;
}

Repeat_raw_volume_data_importer::Repeat_raw_volume_data_importer()
    : m_file_name(),
      m_source_size(),
      m_repeat(),
      m_cache_source_data(false),
      m_size()
{}

//----------------------------------------------------------------------
Repeat_raw_volume_data_importer::~Repeat_raw_volume_data_importer()
{
    if (m_cache_source_data && (m_source != 0))
    {
        if(m_cache_uses==0)
        {
            // Make sure this only gets deleted once per host and only if not in use at all
            mi::base::Lock::Block block(&m_source_lock);
            delete[] m_source;
            m_source = 0;
        }
    }
}

//----------------------------------------------------------------------
mi::Uint8* Repeat_raw_volume_data_importer::load_source_data(
    const std::string& file_name,
    mi::Uint64         source_size) const
{
    INFO_LOG << "Repeat data importer loading volume data file '" << file_name << "'";

    io::File in_file(file_name, std::ios_base::in);
    if (!in_file)
    {
        ERROR_LOG << "Cannot open volume data file '"  << file_name << "'.";
        return 0;
    }

    mi::Uint8* source = new mi::Uint8[source_size];

    if (source_size != in_file.read(source, 0ull, sizeof(mi::Uint8) * source_size))
    {
        ERROR_LOG << "Repeat_raw_volume_data_importer::load_source_data(): fread() failed";
    }
    in_file.close();

    return source;
}

//----------------------------------------------------------------------
mi::Size Repeat_raw_volume_data_importer::estimate(
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

nv::index::IDistributed_data_subset* Repeat_raw_volume_data_importer::create(
    const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
    nv::index::IData_subset_factory*                factory,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{
    const mi::math::Bbox<mi::Sint32, 3> bounds = bounding_box;

    // The size for volumetric data brick
    mi::Uint64 dx = bounds.max.x - bounds.min.x;
    mi::Uint64 dy = bounds.max.y - bounds.min.y;
    mi::Uint64 dz = bounds.max.z - bounds.min.z;

    mi::Uint64 source_size = 
        static_cast<mi::Uint64>(m_source_size.x) *
        static_cast<mi::Uint64>(m_source_size.y) * 
        static_cast<mi::Uint64>(m_source_size.z);

    DEBUG_LOG << "Source size: " << m_source_size;
    DEBUG_LOG << "Size:        " << m_size;

    mi::Uint8* source = 0;

    // Count cache uses
    ++m_cache_uses;

    if (m_cache_source_data)
    {
        mi::base::Lock::Block block(&m_source_lock);
        if (m_source == 0)
        {
            m_source = load_source_data(m_file_name, source_size);
        }
        source = m_source;
    }
    else
    {
        DEBUG_LOG << "Loading non-cache.";
        source = load_source_data(m_file_name, source_size);
    }

    if (source == 0)
    {
        ERROR_LOG << "The file is not available for reading. Returning an invalid brick element.";
        --m_cache_uses;
        return NULL;
    }

    mi::base::Handle<nv::index::IRegular_volume_brick_uint8> volume_brick(
        factory->create<nv::index::IRegular_volume_brick_uint8>());
    
    if (!volume_brick.is_valid_interface()) {
        ERROR_LOG << "Cannot create a volume brick.";
        return 0;
    }

    mi::Uint8* voxel_data_storage = volume_brick->generate_voxel_storage(bounding_box);
    if (voxel_data_storage == 0) {
        ERROR_LOG << "Cannot generate a voxel storage.";
        return 0;
    }

    INFO_LOG << "Scaling brick loaded from file '" << m_file_name
             << "', with bounding box: " << bounds
             << " and repeat factor: " << m_repeat
             << ".";

    mi::Uint32 x_dst = 0;
    for (mi::Sint32 x = bounds.min.x; x < bounds.max.x; ++x)
    {
        if (x < 0 || x >= (mi::Sint32)m_size.x)
        {
            // Skip the border voxels that are outside of the dataset, they will be filled with
            // duplicated boundary data later.
            x_dst++;
            continue;
        }

        mi::Uint32 y_dst = 0;
        for (mi::Sint32 y = bounds.min.y; y < bounds.max.y; ++y)
        {
            if (y < 0 || y >= (mi::Sint32)m_size.y)
            {
                // Skip the border voxels that are outside of the dataset, they will be filled with
                // duplicated boundary data later.
                y_dst++;
                continue;
            }

            // Data offsets, be sure to use 64-bit integers here.
            mi::Uint64 memory_offset = (mi::Uint64)dz * dy * x_dst + dz * y_dst;

            mi::Uint64 source_offset
                = (mi::Uint64)m_source_size.z * m_source_size.y * (x / m_repeat.x) + m_source_size.z * (y / m_repeat.y);

            mi::Sint64 bytes_to_read = bounds.max.z - bounds.min.z;

            // Skip the minimum z border voxels that are outside of the dataset, they will be filled
            // with duplicated boundary data later.
            mi::Sint32 skip_min = -std::min(bounds.min.z, 0);

            // Skip the maximum z border voxels that are outside of the dataset, they will be filled
            // with duplicated boundary data later.
            mi::Sint32 skip_max = std::max(bounds.max.z - (mi::Sint32)m_size.z, 0);

            for (mi::Sint32 z=skip_min; z < bytes_to_read - skip_max; ++z)
            {
                const mi::Uint64 idx = source_offset + (bounds.min.z + z) / m_repeat.z;
                voxel_data_storage[memory_offset + z] = source[idx];
            }

            y_dst++;
        }
        x_dst++;
    }
    --m_cache_uses;

    //
    // Duplicate voxels at the dataset boundary (clamping) to fill up the border if necessary
    //

    // Calculate position of the voxels at the boundary of the current brick in memory. These are
    // the voxels that will be duplicated towards the border.
    const mi::math::Bbox<mi::Sint32, 3> boundary_dst(
        -bounds.min,
        mi::math::Vector<mi::Sint32, 3>(m_size) - bounds.min - mi::math::Vector<mi::Sint32, 3>(1));

    // Duplicate voxels at the minimum z boundary of the dataset (if that is part of the current
    // brick, i.e. the min component of the bounding box is zero)
    for (mi::Sint32 z = 0; z < -bounds.min.z; ++z)
    {
        for (mi::Uint32 x=0; x < dx; ++x)
        {
            for (mi::Uint32 y=0; y < dy; ++y)
            {
                voxel_data_storage[dz * dy * x + dz * y + z] = voxel_data_storage[dz * dy * x + dz * y + boundary_dst.min.z];
            }
        }
    }

    // Duplicate voxels at the maximum z boundary of the dataset (if that is part of the current
    // brick, i.e. the max component of the bounding box is larger than the dataset)
    for (mi::Sint32 z = m_size.z; z < bounds.max.z; ++z)
    {
        for (mi::Uint32 x=0; x < dx; ++x)
        {
            for (mi::Uint32 y=0; y < dy; ++y)
            {
                voxel_data_storage[dz * dy * x + dz * y + (z - bounds.min.z)] = voxel_data_storage[dz * dy * x + dz * y + boundary_dst.max.z];
            }
        }
    }

    // Duplicate min y
    for (mi::Sint32 y = 0; y < -bounds.min.y; ++y)
    {
        for (mi::Uint32 x=0; x < dx; ++x)
        {
            for (mi::Uint32 z=0; z < dz; ++z)
            {
                voxel_data_storage[dz * dy * x + (dz * y) + z] = voxel_data_storage[dz * dy * x + (dz * boundary_dst.min.y) + z];
            }
        }
    }

    // Duplicate max y
    for (mi::Sint32 y = m_size.y; y < bounds.max.y; ++y)
    {
        for (mi::Uint32 x=0; x < dx; ++x)
        {
            for (mi::Uint32 z=0; z < dz; ++z)
            {
                voxel_data_storage[dz * dy * x + dz * (y - bounds.min.y) + z] = voxel_data_storage[dz * dy * x + dz * boundary_dst.max.y + z];
            }
        }
    }

    // Duplicate min x
    for (mi::Sint32 x = 0; x < -bounds.min.x; ++x)
    {
        memcpy(&voxel_data_storage[dz * dy * x], &voxel_data_storage[dz * dy * boundary_dst.min.x], dz * dy);
    }

    // Duplicate max x
    for (mi::Sint32 x = m_size.x; x < bounds.max.x; ++x)
    {
        memcpy(&voxel_data_storage[dz * dy * (x - bounds.min.x)], &voxel_data_storage[dz * dy * boundary_dst.max.x], dz * dy);
    }

    if (!m_cache_source_data)
    {
        delete[] source;
        source = 0;
    }

    volume_brick->retain();     // since the handle will be out of scope next.
    return volume_brick.get();
}

//----------------------------------------------------------------------
void Repeat_raw_volume_data_importer::serialize(
    mi::neuraylib::ISerializer* serializer) const
{
    mi::Uint32 nb_elements = mi::Uint32(m_file_name.size());
    serializer->write(&nb_elements, 1);
    serializer->write(reinterpret_cast<const mi::Uint8*>(m_file_name.c_str()), nb_elements);

    serializer->write(&m_size.x, 3);
    serializer->write(&m_source_size.x, 3);
    serializer->write(&m_repeat.x, 3);
    serializer->write(&m_cache_source_data, 1);
}

//----------------------------------------------------------------------
void Repeat_raw_volume_data_importer::deserialize(
    mi::neuraylib::IDeserializer* deserializer)
{
    mi::Uint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_file_name.resize(nb_elements);
    deserializer->read(reinterpret_cast<mi::Uint8*>(&m_file_name[0]), nb_elements);

    deserializer->read(&m_size.x, 3);
    deserializer->read(&m_source_size.x, 3);
    deserializer->read(&m_repeat.x, 3);
    deserializer->read(&m_cache_source_data, 1);
}

//----------------------------------------------------------------------
mi::base::Uuid Repeat_raw_volume_data_importer::subset_id() const
{
    // currently generate uint8 volume only
    return nv::index::IRegular_volume_brick_uint8::IID();
}

//----------------------------------------------------------------------
const char* Repeat_raw_volume_data_importer::get_configuration() const
{
    return m_configuration.c_str();
}

//----------------------------------------------------------------------
}} // namespace nv::index_common
