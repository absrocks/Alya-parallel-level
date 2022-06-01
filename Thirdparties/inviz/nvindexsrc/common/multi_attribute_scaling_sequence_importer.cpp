/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Volume data importer reading a simple raw data format from file.

#include "forwarding_logger.h"

#include "multi_attribute_scaling_sequence_importer.h"
#include "large_file_io.h"
#include "common_utility.h"

#include <nv/index/iregular_volume.h>
#include <nv/index/iregular_volume_brick.h>

#include <limits>
#include <cstdio>

namespace nv {
namespace index_common {

mi::base::Lock Multi_attribute_scaling_sequence_importer::s_file_access_lock;
mi::Uint64     Multi_attribute_scaling_sequence_importer::s_total_read_bytes = 0;
mi::Float64    Multi_attribute_scaling_sequence_importer::s_total_read_time  = 0.0;

//----------------------------------------------------------------------
Multi_attribute_scaling_sequence_importer::Multi_attribute_scaling_sequence_importer(
    const std::string&                            attribute_file_1,
    const std::string&                            attribute_file_2,
    const std::string&                            output_file,
    const mi::math::Vector_struct<mi::Uint32, 3>& size)
  : m_attribute_file_1(attribute_file_1),
    m_attribute_file_2(attribute_file_2),
    m_output_file(output_file),
    m_size(size),
    m_time_interval(0),
    m_start_time(0),
    m_use_cache(true),
    m_cache_compression(0),
    m_serial_access(false),
    m_show_stats(false)
{
    m_configuration = std::string() +
        "importer=raw\n" +
        "input_file=" + attribute_file_1 + "\n";
}

Multi_attribute_scaling_sequence_importer::Multi_attribute_scaling_sequence_importer()
  : m_time_interval(0),
    m_start_time(0),
    m_use_cache(true),
    m_cache_compression(0),
    m_serial_access(false),
    m_show_stats(false)
{
    // empty
}

//----------------------------------------------------------------------

namespace {

static inline mi::Float32 round(mi::Float32 x)
{
    return std::floor(x + 0.5f);
}

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

mi::Size Multi_attribute_scaling_sequence_importer::estimate(
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

nv::index::IDistributed_data_subset* Multi_attribute_scaling_sequence_importer::create(
    const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
    mi::Uint32                                      time_step,
    nv::index::IData_subset_factory*                factory,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{
    std::string file_name_1;
    {
        mi::Uint32 file_sequence = 0;
        
        if (time_step == 0)
            file_sequence = m_start_time;
        else
            file_sequence = m_start_time + time_step*m_time_interval;
        
        std::ostringstream time_step_str;
            time_step_str << m_attribute_file_1 << file_sequence << ".raw";
        file_name_1 = time_step_str.str(); 
    }
    
    std::string file_name_2;
    {
        mi::Uint32 file_sequence = 0;
        
        if (time_step == 0)
            file_sequence = m_start_time;
        else
            file_sequence = m_start_time + time_step*m_time_interval;
        
        std::ostringstream time_step_str;
            time_step_str << m_attribute_file_2 << file_sequence << ".raw";
        file_name_2 = time_step_str.str(); 
    }

    std::string output_file_name;
    {
        mi::Uint32 file_sequence = 0;
        
        if (time_step == 0)
            file_sequence = m_start_time;
        else
            file_sequence = m_start_time + time_step*m_time_interval;
        
        std::ostringstream time_step_str;
            time_step_str << m_output_file << file_sequence << ".raw";
        output_file_name = time_step_str.str(); 
    }
    
    return create(file_name_1, file_name_2, output_file_name, bounding_box, factory, dice_transaction);
}

nv::index::IDistributed_data_subset* Multi_attribute_scaling_sequence_importer::create(
    const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
    nv::index::IData_subset_factory*                factory,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{
    std::string file_name_1;
    {
        std::ostringstream time_step_str;
            time_step_str << m_attribute_file_1 << m_start_time << ".raw";
        file_name_1 = time_step_str.str(); 
    }
    
    std::string file_name_2;
    {
        std::ostringstream time_step_str;
            time_step_str << m_attribute_file_2 << m_start_time << ".raw";
        file_name_2 = time_step_str.str(); 
    }
    
    std::string output_file_name;
    {
        std::ostringstream time_step_str;
            time_step_str << m_output_file << m_start_time << ".raw";
        output_file_name = time_step_str.str(); 
    }

    return create(file_name_1, file_name_2, output_file_name, bounding_box, factory, dice_transaction);
}

// TODO: Design this after GTC-2016 if required.
mi::Uint8 scale_attributes(
    const mi::Float32                                               attribute_value_1,
    const mi::Float32                                               attribute_value_2,
    const std::map<mi::Uint32, mi::math::Vector<mi::Float32, 2> >&  attribute_minmax_map)
{

    std::map<mi::Uint32, mi::math::Vector<mi::Float32, 2> >::const_iterator it_attr_1 = attribute_minmax_map.find(1);
    std::map<mi::Uint32, mi::math::Vector<mi::Float32, 2> >::const_iterator it_attr_2 = attribute_minmax_map.find(2);
    
    // // const mi::Uint8 max_uint = 255; //11111111
    // // const mi::Uint8 min_uint = 0; //00000000
    
    const mi::Uint8 max_uint_attribute_1 = 200; //11001000
    const mi::Uint8 min_uint_attribute_1 = 0; //00000000
    
    const mi::Uint8 max_uint_attribute_2 = 255; //11111111
    const mi::Uint8 min_uint_attribute_2 = 201; //11001001

    mi::Uint8 attribute_scaled_1 = ((attribute_value_1 - it_attr_1->second.x) * (max_uint_attribute_1 - min_uint_attribute_1) / (it_attr_1->second.y - it_attr_1->second.x) + min_uint_attribute_1);
    mi::Uint8 attribute_scaled_2 = (round(attribute_value_2 - it_attr_2->second.x) * (max_uint_attribute_2 - min_uint_attribute_2) / (it_attr_2->second.y - it_attr_2->second.x) + min_uint_attribute_2);

    return (attribute_scaled_2 - min_uint_attribute_2 + attribute_scaled_1);       
    
    //mi::Float32 attribute_scaled_1 = (attribute_value_1 - it_attr_1->second.x) / (it_attr_1->second.y - it_attr_1->second.x);
	//return static_cast<mi::Uint8>(round(attribute_scaled_1 * max_uint));
}

nv::index::IDistributed_data_subset* Multi_attribute_scaling_sequence_importer::create(
    const std::string&                              attribute_file_1,
    const std::string&                              attribute_file_2,
    const std::string&                              file_name,
    const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
    nv::index::IData_subset_factory*                factory,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{
    mi::base::Lock::Block block;
    if (m_serial_access)
    {
        block.set(&s_file_access_lock); // Acquire the lock
    }

    // The size for volumetric data brick that needs to be loaded
    const mi::math::Bbox<mi::Sint32, 3> bounds = bounding_box;
    mi::Sint64 dx = bounds.max.x - bounds.min.x;
    mi::Sint64 dy = bounds.max.y - bounds.min.y;
    mi::Sint64 dz = bounds.max.z - bounds.min.z;
    mi::Sint64 volume_brick_size = dx * dy * dz;

    INFO_LOG << "Multi_attribute_scaling_sequence_importer loads attribute-1: " 
        << attribute_file_1 << " with size: " << m_size << "', bounds: " << bounds;
    INFO_LOG << "Multi_attribute_scaling_sequence_importer loads attribute-2: " 
        << attribute_file_2 << " with size: " << m_size << "', bounds: " << bounds;

    // Create a volume brick
    mi::base::Handle<nv::index::IRegular_volume_brick_uint8> 
        volume_brick(factory->create<nv::index::IRegular_volume_brick_uint8>());
    if (!volume_brick.is_valid_interface()) {
        ERROR_LOG << "Cannot create a volume brick.";
        return 0;
    }

    // Allocate voxel data storage via volume_brick.
    mi::Uint8* voxel_data = volume_brick->generate_voxel_storage(bounding_box);
    if (voxel_data == 0) 
    {
        ERROR_LOG << "Cannot generate voxel storage.";
        return 0;
    }
    
    // ---------------------------------------------------------------------------------------------------
    // Try importing the data from a temporary cache file possibly generated in a previous import process.
    std::string cache_file_name;
    io::File cache_file(m_cache_compression);
    if (m_use_cache && !cache_file)
    {    
        std::ostringstream cache_file_name_str;
        {
            cache_file_name_str << file_name << "_"
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
        INFO_LOG << "Using cache file '" << cache_file_name;

        {
            // Fetch data for the entire brick at once
            if (cache_file.read(voxel_data, 0ll, volume_brick_size) != static_cast<mi::Uint64>(volume_brick_size))
            {
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

        volume_brick->retain(); // since the handle will be out of scope.
        return volume_brick.get();
    }

    // // -------------------------------------------------------------------------------------------
    // // Otherwise, fetch the data from the original dataset file.
    // io::File in_file_1(attribute_file_1.c_str(), std::ios_base::in);
    // if (!in_file_1)
    // {
        // // If the file is not available for reading then create a subset that
        // // only contains invalid values mi::Uint8(0).
        // ERROR_LOG << "Cannot open volume data file '"  << attribute_file_1 << "'. "
                  // << "Will create a database element with invalid value '0' (type: mi::Uint8).";
        // return NULL;
    // }
    
    // io::File in_file_2(attribute_file_2.c_str(), std::ios_base::in);
    // if (!in_file_2)
    // {
        // // If the file is not available for reading then create a subset that
        // // only contains invalid values mi::Uint8(0).
        // ERROR_LOG << "Cannot open volume data file '"  << attribute_file_2 << "'. "
                  // << "Will create a database element with invalid value '0' (type: mi::Uint8).";
        // return NULL;
    // }

    mi::Uint64 input_size = static_cast<mi::Uint64>(m_size.x)*m_size.y*m_size.z;
    mi::Float32* attribute_data_1 = new mi::Float32[input_size];
    FILE* attribute_fp_1 = fopen(attribute_file_1.c_str(), "rb");
    fread(attribute_data_1, sizeof(mi::Float32), input_size, attribute_fp_1);
    fclose(attribute_fp_1);
    
    mi::Float32* attribute_data_2 = new mi::Float32[input_size];
    FILE* attribute_fp_2 = fopen(attribute_file_2.c_str(), "rb");
    fread(attribute_data_2, sizeof(mi::Float32), input_size, attribute_fp_2);
    fclose(attribute_fp_2);
    
    mi::Uint64 x_dst = 0;
    for (mi::Sint32 x = bounds.min.x; x < bounds.max.x; ++x)
    {
        if (x < 0 || x >= (int)m_size.x)
        {
            // Skip the border voxels that are outside of the dataset, they will be filled with
            // duplicated boundary data later.
            x_dst++;
            continue;
        }

        mi::Uint64 y_dst = 0;
        for (mi::Sint32 y = bounds.min.y; y < bounds.max.y; ++y)
        {
            if (y < 0 || y >= (int)m_size.y)
            {
                // Skip the border voxels that are outside of the dataset, they will be filled with
                // duplicated boundary data later.
                y_dst++;
                continue;
            }

            mi::Uint64 z_dst = 0;
			for (mi::Sint64 z = bounds.min.z; z < bounds.max.z; ++z)
			{
                if (z < 0 || z >= (int)m_size.z)
				{
					z_dst++;
					continue;
				}
                mi::Uint64 memory_offset = static_cast<mi::Uint64>(dz) * dy * x_dst + dz * y_dst + z_dst;

                mi::Uint64 read_idx = static_cast<mi::Uint64>(x)
                    + (y)*m_size.x
                    + (z)*m_size.x*m_size.y;

                voxel_data[memory_offset] = scale_attributes(attribute_data_1[read_idx], attribute_data_2[read_idx], m_attribute_minmax_map);
                
                z_dst++;
            }
            y_dst++;
        }
        x_dst++;
    }
    delete[] attribute_data_1;
    delete[] attribute_data_2;


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
                voxel_data[dz * dy * x + dz * y + z]
                    = voxel_data[dz * dy * x + dz * y + boundary_dst.min.z];
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
                voxel_data[dz * dy * x + dz * y + (z - bounds.min.z)]
                    = voxel_data[dz * dy * x + dz * y + boundary_dst.max.z];
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
                voxel_data[dz * dy * x + (dz * y) + z]
                    = voxel_data[dz * dy * x + (dz * boundary_dst.min.y) + z];
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
                voxel_data[dz * dy * x + dz * (y - bounds.min.y) + z]
                    = voxel_data[dz * dy * x + dz * boundary_dst.max.y + z];
            }
        }
    }

    // Duplicate min x
    for (mi::Sint64 x = 0; x < -bounds.min.x; ++x)
    {
        memcpy(&voxel_data[dz * dy * x], &voxel_data[dz * dy * boundary_dst.min.x], dz * dy);
    }

    // Duplicate max x
    for (mi::Sint64 x = m_size.x; x < bounds.max.x; ++x)
    {
        memcpy(&voxel_data[dz * dy * (x - bounds.min.x)], &voxel_data[dz * dy * boundary_dst.max.x], dz * dy);
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
			if (cache_file.write(voxel_data, 0ll, volume_brick_size) != static_cast<mi::Uint64>(volume_brick_size))
			{
				ERROR_LOG << "Raw_volume_importer::create(): "
						  << "Error writing to file '" << cache_file_name << "', "
						  << "unable to write " << volume_brick_size << " bytes.";
			}
			cache_file.close();
		}
    }

    volume_brick->retain();     // since the handle will be out of scope.
    return volume_brick.get();
}

//----------------------------------------------------------------------
mi::base::Uuid Multi_attribute_scaling_sequence_importer::subset_id() const
{
    // currently generate uint8 volume only
    return nv::index::IRegular_volume_brick_uint8::IID();
}

//----------------------------------------------------------------------
const char* Multi_attribute_scaling_sequence_importer::get_configuration() const
{
    return m_configuration.c_str();
}

//----------------------------------------------------------------------
void Multi_attribute_scaling_sequence_importer::serialize(
    mi::neuraylib::ISerializer* serializer) const
{
    mi::Uint32 nb_elements = mi::Uint32(m_attribute_file_1.size());
    serializer->write(&nb_elements, 1);
    serializer->write(reinterpret_cast<const mi::Uint8*>(m_attribute_file_1.c_str()), nb_elements);
    
    nb_elements = mi::Uint32(m_attribute_file_2.size());
    serializer->write(&nb_elements, 1);
    serializer->write(reinterpret_cast<const mi::Uint8*>(m_attribute_file_2.c_str()), nb_elements);
    
    nb_elements = mi::Uint32(m_output_file.size());
    serializer->write(&nb_elements, 1);
    serializer->write(reinterpret_cast<const mi::Uint8*>(m_output_file.c_str()), nb_elements);
    
    {
		const mi::Size nb_elements = m_attribute_minmax_map.size();
		serializer->write(&nb_elements, 1);

		std::map<mi::Uint32, mi::math::Vector<mi::Float32, 2> >::const_iterator itr = m_attribute_minmax_map.begin();
		for (; itr != m_attribute_minmax_map.end(); ++itr)
		{
			serializer->write(&itr->first, 1);
			serializer->write(&itr->second.x, 2);
		}
	}

    serializer->write(&m_size.x, 3);
    serializer->write(&m_time_interval, 1);
    serializer->write(&m_start_time, 1);
    serializer->write(&m_use_cache, 1);
    serializer->write(&m_cache_compression, 1);
    serializer->write(&m_serial_access, 1);
    serializer->write(&m_show_stats, 1);

    nb_elements = mi::Uint32(m_configuration.size());
    nb_elements = mi::Uint32(m_configuration.size());
    serializer->write(&nb_elements, 1);
    serializer->write(reinterpret_cast<const mi::Uint8*>(m_configuration.c_str()), nb_elements);
}

//----------------------------------------------------------------------
void Multi_attribute_scaling_sequence_importer::deserialize(
    mi::neuraylib::IDeserializer* deserializer)
{
    mi::Uint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_attribute_file_1.resize(nb_elements);
    deserializer->read(reinterpret_cast<mi::Uint8*>(&m_attribute_file_1[0]), nb_elements);
    
    nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_attribute_file_2.resize(nb_elements);
    deserializer->read(reinterpret_cast<mi::Uint8*>(&m_attribute_file_2[0]), nb_elements);
    
    nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_output_file.resize(nb_elements);
    deserializer->read(reinterpret_cast<mi::Uint8*>(&m_output_file[0]), nb_elements);
    
    {
		mi::Size nb_elements = 0;
		deserializer->read(&nb_elements, 1);

		for (mi::Uint32 i = 0; i < nb_elements; ++i)
		{
			mi::Sint32 attribute_id = 0;
			deserializer->read(&attribute_id, 1);
			mi::math::Vector<mi::Float32, 2> minmax;
			deserializer->read(&minmax.x, 2);
			m_attribute_minmax_map[attribute_id] = minmax;
		}
	}

    deserializer->read(&m_size.x, 3);
    deserializer->read(&m_time_interval, 1);
    deserializer->read(&m_start_time, 1);
    deserializer->read(&m_use_cache, 1);
    deserializer->read(&m_cache_compression, 1);
    deserializer->read(&m_serial_access, 1);
    deserializer->read(&m_show_stats, 1);

    deserializer->read(&nb_elements, 1);
    m_configuration.resize(nb_elements);
    deserializer->read(reinterpret_cast<mi::Uint8*>(&m_configuration[0]), nb_elements);
}

//----------------------------------------------------------------------
}} // namespace nv::index_common
