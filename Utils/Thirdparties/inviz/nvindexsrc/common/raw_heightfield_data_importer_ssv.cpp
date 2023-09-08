/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Height field data importer reading raw binary formatted floats as elevation values. obsolete

#include "raw_heightfield_data_importer_ssv.h"

#include "common_utility.h"
#include "forwarding_logger.h"
#include "normal_encoder.h"
#include "ppm_io.h"
#include "raw_heightfield_data_type.h"

#include <cstdio>
#include <cassert>
#include <fstream>

#ifndef WIN32
#include <sys/time.h>
#endif

namespace nv {
namespace index_common {

mi::base::Lock Raw_heightfield_data_importer_ssv::s_file_access_lock;

mi::Uint64     Raw_heightfield_data_importer_ssv::s_total_read_bytes = 0;
mi::Float64    Raw_heightfield_data_importer_ssv::s_total_read_time  = 0.0;

//----------------------------------------------------------------------
static void output_heightfield_import_stat(
    const std::string & header,
    mi::Uint64  bytes_read,
    mi::Float64 time_taken,
    mi::Uint64  total_read_bytes,
    mi::Float64 total_read_time)
{
    const mi::Float64 bytes_read_f64       = static_cast<mi::Float64>(bytes_read);
    const mi::Float64 total_read_bytes_f64 = static_cast<mi::Float64>(total_read_bytes);
    INFO_LOG << header << (bytes_read_f64 / (1024.0 * 1024.0)) << " MB in " << time_taken << "s, "
             << "\ttotal " << (total_read_bytes_f64 / (1024 * 1024)) << " MB in " << total_read_time << "s, "
             << "\tavg " << (total_read_bytes_f64 / total_read_time) / (1024 * 1024) << " MB/s.";
}

//----------------------------------------------------------------------
Raw_heightfield_data_importer_ssv::Raw_heightfield_data_importer_ssv(
    const std::string&                            filename,
    const mi::math::Vector_struct<mi::Uint32, 2>& size)
{
    m_filename           = filename;
    m_heightfield_size   = size;
    m_use_cache          = false;
    m_serial_access      = false;
    m_show_stats         = false;
    m_image_verification = false;

    m_configuration = std::string() +
        "importer=raw\n" +
        "input_file=" + m_filename + "\n";
}

//----------------------------------------------------------------------
Raw_heightfield_data_importer_ssv::Raw_heightfield_data_importer_ssv()
{
    m_filename           = "";
    m_use_cache          = false;
    m_serial_access      = false;
    m_show_stats         = false;
    m_image_verification = false;
}

//----------------------------------------------------------------------
Raw_heightfield_data_importer_ssv::~Raw_heightfield_data_importer_ssv()
{
    // empty
}

//----------------------------------------------------------------------
const char* Raw_heightfield_data_importer_ssv::get_configuration() const
{
    return m_configuration.c_str();
}

//----------------------------------------------------------------------
mi::Size Raw_heightfield_data_importer_ssv::estimate(
    const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{
    // Size of the height field patch in memory
    const mi::Size range_i = bounding_box.max.x - bounding_box.min.x;
    const mi::Size range_j = bounding_box.max.y - bounding_box.min.y;
    const mi::Size patch_size = range_i*range_j;
    // ... elevation values + normal values per grid position:
    const mi::Size data_size = patch_size * (sizeof(mi::Float32) + sizeof(mi::Float32)*3);
    // note: color values not yet supported!
    return data_size;
}

//----------------------------------------------------------------------
nv::index::IDistributed_data_subset* Raw_heightfield_data_importer_ssv::create(
    const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
    nv::index::IData_subset_factory*                factory,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{
    mi::base::Lock::Block block;
    if (m_serial_access)
    {
        block.set(&s_file_access_lock); // Acquire the lock
    }

    // For more than 0.5G heightfield data entries, use 64bit.
    const mi::math::Vector<mi::Sint64, 2> heightfield_size(
        static_cast< mi::Sint64 >(m_heightfield_size.x),
        static_cast< mi::Sint64 >(m_heightfield_size.y));

    mi::math::Bbox<mi::Sint64, 2> patch_bbox(
        static_cast< mi::Sint64 >(bounding_box.min.x),
        static_cast< mi::Sint64 >(bounding_box.min.y),
        static_cast< mi::Sint64 >(bounding_box.max.x),
        static_cast< mi::Sint64 >(bounding_box.max.y));


    // Size of the height field patch in memory
    const mi::Sint64 range_i = patch_bbox.max.x - patch_bbox.min.x;
    const mi::Sint64 range_j = patch_bbox.max.y - patch_bbox.min.y;
    const mi::Sint64 grid_size = range_i * range_j;
    assert(range_i > 0);
    assert(range_j > 0);

    // Load a cache file if available
    std::string cache_file_name;
    if (m_use_cache)
    {
        nv::index::IRegular_heightfield_patch* heightfield =
            try_load_from_cache_file(bounding_box, grid_size, factory, cache_file_name);
        if (heightfield != 0)
        {
            // have read the heightfield from the cache file.
            return heightfield;
        }
    }

    // Open file
    FILE* p_file = fopen(m_filename.c_str(), "rb");
    if(p_file == 0)
    {
        ERROR_LOG << "Cannot open heightfield elevation file: " << m_filename;
        return 0;               // failed to load
    }

    // read the header
    if (!load_header(p_file))
    {
        fclose(p_file);
        return 0;
    }

    // read the data
    mi::Float32 *elevation_data = 0;
    mi::math::Vector_struct<mi::Float32, 3>* normal_vector_data = 0;
    mi::base::Handle<nv::index::IRegular_heightfield_patch> heightfield_patch(
        factory->create<nv::index::IRegular_heightfield_patch>());
    if (!heightfield_patch.is_valid_interface())
    {
        ERROR_LOG << "Cannot create IRegular_heightfield_patch.";
        return 0;
    }

    // Initialize patch and allocate the heightfield data storage
    {
        const mi::math::Bbox_struct<mi::Sint32, 2> bbox_patch_rect_st = {
            bounding_box.min.x, bounding_box.min.y,
            bounding_box.max.x, bounding_box.max.y,
        };
        if (!heightfield_patch->initialize(bbox_patch_rect_st))
        {
            ERROR_LOG << "Fail to initialize heightfield patch";
            return 0;
        }

        elevation_data = heightfield_patch->generate_elevation_storage();
        if (elevation_data == 0)
        {
            ERROR_LOG << "Cannot generate elevation storage.";
            return 0;
        }

        normal_vector_data = heightfield_patch->generate_normal_storage();
        if (normal_vector_data == 0)
        {
            elevation_data = 0; // invalidate
            ERROR_LOG << "Cannot generate normal storage.";
            return 0;
        }
    }

    // Read the data to fill the data storage
    load_data(p_file,
              cache_file_name,
              heightfield_size,
              bounding_box,
              patch_bbox,
              // file_type,
              elevation_data,
              normal_vector_data
        );

    heightfield_patch->retain(); // since the handle will be out of scope next.
    return heightfield_patch.get();
}

//----------------------------------------------------------------------
void Raw_heightfield_data_importer_ssv::serialize(
    mi::neuraylib::ISerializer * serializer) const
{
    mi::Uint32 nb_elements = mi::Uint32(m_filename.size());
    serializer->write(&nb_elements, 1);
    serializer->write(reinterpret_cast<const mi::Uint8*>(m_filename.c_str()), nb_elements);
    serializer->write(&m_heightfield_size.x, 2);
    serializer->write(&m_use_cache, 1);
    serializer->write(&m_serial_access, 1);
    serializer->write(&m_show_stats, 1);
    serializer->write(&m_image_verification, 1);
}

//----------------------------------------------------------------------
void Raw_heightfield_data_importer_ssv::deserialize(
    mi::neuraylib::IDeserializer * deserializer)
{
    mi::Uint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_filename.resize(nb_elements);
    deserializer->read(reinterpret_cast<mi::Uint8*>(&m_filename[0]), nb_elements);
    deserializer->read(&m_heightfield_size.x, 2);
    deserializer->read(&m_use_cache, 1);
    deserializer->read(&m_serial_access, 1);
    deserializer->read(&m_show_stats, 1);
    deserializer->read(&m_image_verification, 1);
}

//----------------------------------------------------------------------
nv::index::IRegular_heightfield_patch* Raw_heightfield_data_importer_ssv::try_load_from_cache_file(
    const mi::math::Bbox_struct<mi::Sint32, 3>& bounding_box,
    const mi::Sint64 grid_size,
    nv::index::IData_subset_factory* factory,
    std::string& cache_file_name
    ) const
{
    assert(m_use_cache);

    std::ostringstream cache_file_name_sstr;
    cache_file_name_sstr << m_filename << "_"
                         << bounding_box.min.x << "." << bounding_box.max.x << "x"
                         << bounding_box.min.y << "." << bounding_box.max.y
                         << "_raw_importer.tmp";

    INFO_LOG << "Trying to load cache file: " << cache_file_name_sstr.str();
    cache_file_name = cache_file_name_sstr.str();

    std::ifstream cache_file;
    cache_file.open(cache_file_name.c_str(), std::ios_base::in | std::ios_base::binary);

    if(!cache_file)
    {
        INFO_LOG << "Height field cache file '" << cache_file_name << "' is not yet available.";
    }
    else
    {
        // Create a patch and allocate data storage
        File_type file_type = FILE_TYPE_SSV;
        mi::Float32 *elevation_data = 0;
        mi::math::Vector_struct<mi::Float32, 3>* normal_vector_data = 0;
        mi::base::Handle<nv::index::IRegular_heightfield_patch> heightfield_patch(
            factory->create<nv::index::IRegular_heightfield_patch>());
        if (!heightfield_patch.is_valid_interface())
        {
            ERROR_LOG << "Cannot create IRegular_heightfield_patch.";
            return 0;
        }

        // Initialize patch and allocate the heightfield data storage
        {
            const mi::math::Bbox_struct<mi::Sint32, 2> bbox_patch_rect_st = {
                bounding_box.min.x, bounding_box.min.y,
                bounding_box.max.x, bounding_box.max.y,
            };
            if (!heightfield_patch->initialize(bbox_patch_rect_st))
            {
                ERROR_LOG << "Fail to initialize heightfield patch";
                return 0;
            }

            elevation_data = heightfield_patch->generate_elevation_storage();
            if (elevation_data == 0)
            {
                ERROR_LOG << "Cannot generate elevation storage.";
                return 0;
            }

            normal_vector_data = heightfield_patch->generate_normal_storage();
            if (normal_vector_data == 0)
            {
                elevation_data = 0; // invalidate
                ERROR_LOG << "Cannot generate normal storage.";
                return 0;
            }
        }

        {
            // tmp_encoded_normal_data scope
            mi::Float64 time_start = nv::index_common::get_time();
                
            // get file type
            cache_file.read(reinterpret_cast<char*>(&file_type), sizeof(File_type));

            // Fetch data for the entire patch at once
            cache_file.read(reinterpret_cast<char*>(elevation_data), grid_size * sizeof(mi::Float32));

            if (file_type != FILE_TYPE_SSV)
            {
                ERROR_LOG << "Cache read: not handled raw heightfield importer type: " << file_type;
            }

            mi::Uint32* tmp_encoded_normal_data = new mi::Uint32[grid_size];
            cache_file.read(reinterpret_cast<char*>(tmp_encoded_normal_data), grid_size * sizeof(mi::Uint32));
            for(mi::Uint32 i=0; i<grid_size; ++i)
            {
                normal_vector_data[i] = Normal_encoder::instance()->get_normal(tmp_encoded_normal_data[i]);
            }
            delete[] tmp_encoded_normal_data;
            tmp_encoded_normal_data = 0;

            mi::Uint64 bytes_read = grid_size * (sizeof(mi::Float32) + sizeof(mi::Uint32));

            mi::math::Vector<mi::Float32, 2> height_range;
            cache_file.read(reinterpret_cast<char*>(&height_range), sizeof(height_range));
            bytes_read += sizeof(height_range);

            if (!cache_file)
            {
                ERROR_LOG << "Height field importer: error reading from file '" << cache_file_name << "'";
            }
            cache_file.close();

            const mi::Float64 time_end = nv::index_common::get_time();
            const mi::Float64 time_taken = time_end - time_start;
            s_total_read_time += time_taken;
            s_total_read_bytes += bytes_read;

            if (m_show_stats)
            {
                output_heightfield_import_stat("Height field (cached): Read ", bytes_read, time_taken,
                                               s_total_read_bytes, s_total_read_time);
            }
        }
            
        heightfield_patch->retain(); // since the handle will be out of scope next.
        return heightfield_patch.get();
    }
    return 0;
}

//----------------------------------------------------------------------
bool Raw_heightfield_data_importer_ssv::load_header(FILE* p_file) const
{
    assert(p_file != 0);

    std::string magic_str;
    for(mi::Sint32 i = 0; (i < 3) && (!feof(p_file)); ++i)
    {
        const mi::Sint32 c = fgetc(p_file);
        if(c == -1)
        {
            ERROR_LOG << "Unexpected EOF when load a heightfield elevation file: " << m_filename;
            // fclose(p_file);
            return false;
        }
        magic_str += static_cast< char >(c);
    }

    if(magic_str != "ssv")
    {
        ERROR_LOG << m_filename << " is not a binary formatted heightfield elevation file.";
        return false;
    }

    // Read height field name
    char event_name[256];
    if(fgets(event_name, static_cast<int>((sizeof(event_name)/sizeof(char))), p_file) == 0)
    {
        ERROR_LOG << "Failed to get the heightfield name from the file [" << m_filename << "]";
        return false;
    }

    // Read IJ-dimensions of entire height field, and compare the value we
    // already have for sanity check.
    size_t const dim_count = 2;
    mi::Uint32 dim_ij[dim_count] = { 0, 0, };
    size_t const dim_read_count = fread((void*)&(dim_ij[0]), sizeof(mi::Uint32), dim_count, p_file);
    if (dim_count != dim_read_count)
    {
        // fclose(p_file);
        ERROR_LOG << "Failed to read the heightfield dimension [" << m_filename << "]";
        return false;        
    }
    
    if ((dim_ij[0] != m_heightfield_size.x) || (dim_ij[1] != m_heightfield_size.y))
    {
        ERROR_LOG << "The given size (" << m_heightfield_size.x << "," << m_heightfield_size.y
                  << ") doesn't match the heightfield dataset size (" << dim_ij[0] << "," << dim_ij[1] << ")";

        assert((dim_ij[0] == m_heightfield_size.x) && (dim_ij[1] == m_heightfield_size.y));
        // fclose(p_file);
        return false;
    }

    return true;
}

//----------------------------------------------------------------------
void Raw_heightfield_data_importer_ssv::load_data(
    FILE* p_file,
    std::string& cache_file_name,
    const mi::math::Vector<mi::Sint64, 2>& heightfield_size,
    const mi::math::Bbox_struct<mi::Sint32, 3>& bounding_box,
    const mi::math::Bbox<mi::Sint64, 2>& patch_bbox,
    mi::Float32* elevation_data,
    mi::math::Vector_struct<mi::Float32, 3>* normal_vector_data
    ) const
{
    assert(p_file != 0);

    const mi::Sint64 range_i = patch_bbox.max.x - patch_bbox.min.x;
    const mi::Sint64 range_j = patch_bbox.max.y - patch_bbox.min.y;
    const mi::Sint64 grid_size = range_i * range_j;

    mi::Uint32* tmp_encoded_normal_data = 0;
    {
        // tmp_encoded_normal_data scope
        const mi::Sint64 raw_data_begin = ftell(p_file);
        const mi::Sint64 normals_begin = raw_data_begin + sizeof(mi::Float32) * heightfield_size.x * heightfield_size.y;
        const mi::Float64 time_start = nv::index_common::get_time();

        // reading elevation values
        mi::Uint64 bytes_read = 0;
        for (mi::Sint64 j = 0; j < range_j; ++j)
        {
            const mi::Sint64 idx_file = ((patch_bbox.min.y + j) * heightfield_size.x) + patch_bbox.min.x;
            const mi::Sint64 idx_map  = j * range_i;

            fseek(p_file, raw_data_begin + idx_file * sizeof(mi::Float32), SEEK_SET);
            fread((void*)&(elevation_data[idx_map]), sizeof(mi::Float32), range_i, p_file);
            bytes_read += range_i * sizeof(mi::Float32);
        }
        
        // reading elevation normals
        tmp_encoded_normal_data = new mi::Uint32[grid_size];
                    
        for (mi::Sint64 j=0; j < range_j; ++j)
        {
            const mi::Sint64 idx_file = ((patch_bbox.min.y + j) * heightfield_size.x) + patch_bbox.min.x;
            const mi::Sint64 idx_map  = j * range_i;

            fseek(p_file, normals_begin + idx_file * sizeof(mi::Uint32), SEEK_SET);
            fread((void*)&(tmp_encoded_normal_data[idx_map]), sizeof(mi::Uint32), range_i, p_file);
            bytes_read += range_i * sizeof(mi::Uint32);
        }
                    
        for(mi::Uint32 i=0; i<grid_size; ++i)
        {
            normal_vector_data[i] = Normal_encoder::instance()->get_normal(tmp_encoded_normal_data[i]);
        }
                    
        fclose(p_file);
        
        const mi::Float64 time_end = nv::index_common::get_time();

        // Determine extent of bounding box in k-dimension from height values
        mi::Float32 min_k = mi::base::numeric_traits<mi::Float32>::max();
        mi::Float32 max_k = mi::base::numeric_traits<mi::Float32>::negative_max();
        for (mi::Sint64 i=0; i < grid_size; ++i)
        {
            mi::Float32& value = elevation_data[i];
            if (value > -1.f)
            {
                if (value < min_k)
                    min_k = value;
                if (value > max_k)
                    max_k = value;
            }
            else
            {
                // Convert to new format
                value = nv::index::IRegular_heightfield_patch::get_hole_value();
            }
        }

        // k-min and k-max should be different to avoid missing planar heightfields
        if (max_k == min_k)
        {
            max_k = min_k + 1.f;
        }

        const mi::math::Vector<mi::Float32, 2> height_range(min_k, max_k);

        INFO_LOG << "Loaded heightfield (raw_ssv) '" << m_filename << "' patch bounds: " << patch_bbox
                 << " height range: " << height_range << " bbox: " << bounding_box;

        mi::Float64 time_taken = time_end - time_start;
        s_total_read_time += time_taken;
        s_total_read_bytes += bytes_read;

        if (m_show_stats)
        {
            output_heightfield_import_stat("Heightfield: Read ", bytes_read, time_taken,
                                           s_total_read_bytes, s_total_read_time);
        }

        //
        // Save data to cache file to accelerate the loading in the future
        //
        if (m_use_cache)
        {
            std::ofstream cache_file;
            cache_file.open(cache_file_name.c_str(), std::ios_base::out | std::ios_base::binary);
            if (!cache_file)
            {
                INFO_LOG << "Height field cache file '" << cache_file_name << "' cannot be opened for writing.";
            }
            else
            {
                INFO_LOG << "Writing the height field patch data to a cache file '" << cache_file_name << "'";
                File_type file_type = FILE_TYPE_COUNT;
                cache_file.write(reinterpret_cast<char*>(&file_type), sizeof(File_type));
                if (file_type != FILE_TYPE_SSV)
                {
                    ERROR_LOG << "Wrong magic, is this ssv file? This may cause a problem.";
                }

                cache_file.write(reinterpret_cast<char*>(elevation_data), grid_size * sizeof(mi::Float32));
                cache_file.write(reinterpret_cast<char*>(tmp_encoded_normal_data), grid_size * sizeof(mi::Uint32));

                cache_file.write(reinterpret_cast<const char*>(&height_range), sizeof(height_range));
                cache_file.close();
            }
        }

        //
        // Verify the import by storing the elevation values in a ppm file.
        //
        if(m_image_verification)
        {
            std::vector<mi::Float32> elevation_values;
            for (mi::Uint32 i = 0; i < grid_size; ++i)
            {
                elevation_values.push_back(elevation_data[i]);
            }
            std::ostringstream ppm_file_sstr;
            ppm_file_sstr << "heightfield_elevation_values_"
                          << bounding_box.min.x << "." << bounding_box.max.x << "x"
                          << bounding_box.min.y << "." << bounding_box.max.y
                          << "_raw_importer.ppm";
            nv::index_common::write_ppm(
                ppm_file_sstr.str(),
                elevation_values,
                height_range,
                static_cast<mi::Sint32>(range_i),
                static_cast<mi::Sint32>(range_j));
        }
    }
    
    if (tmp_encoded_normal_data)
    {
        delete[] tmp_encoded_normal_data;
    }
}

//----------------------------------------------------------------------
}} // namespace nv::index_common
