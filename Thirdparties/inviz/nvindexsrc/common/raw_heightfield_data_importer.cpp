/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Heightfield data importer reading raw binary formatted elevation (float) and normal (float,3).

#include "raw_heightfield_data_importer.h"

#include "common_utility.h"
#include "forwarding_logger.h"
#include "ppm_io.h"

#include <cstdio>
#include <cassert>
#include <fstream>
#include <sstream>

#ifndef WIN32
#include <sys/time.h>
#endif

namespace nv {
namespace index_common {

mi::base::Lock Raw_heightfield_data_importer::s_file_access_lock;

mi::Uint64     Raw_heightfield_data_importer::s_total_read_bytes = 0;
mi::Float64    Raw_heightfield_data_importer::s_total_read_time  = 0.0;

namespace
{

//----------------------------------------------------------------------
void output_heightfield_import_stat(
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
/// check ssv file
/// \param[in] fname file name to check
/// \return true if the file has a ssv header
bool is_ssv_file(const std::string& fname)
{
    FILE * p_file = fopen(fname.c_str(), "rb");
    if(p_file == 0)
    {
        ERROR_LOG << "Cannot open raw heightfield file: '" << fname << "'";
        return false;
    }

    bool is_ssv_header = true;
    char const * const expected_header = "ssv";
    // get 'ssv' at top of the file 
    for (mi::Sint32 i = 0; i < 3; ++i)
    {
        const mi::Sint32 ch = getc(p_file);
        if (feof(p_file))
        {
            is_ssv_header = false; // too short
            break;
        }

        if (static_cast<char>(ch) != expected_header[i])
        {
            is_ssv_header = false; // not ssv
            break;
        }            
    }
    
    fclose(p_file);

    return is_ssv_header;
}

//----------------------------------------------------------------------
} // anonymous namespace

//----------------------------------------------------------------------
Raw_heightfield_data_importer::Raw_heightfield_data_importer(
    const std::string&                            filename,
    const mi::math::Vector_struct<mi::Uint32, 2>& size)
{
    m_filename           = filename;
    m_heightfield_size   = size;
    m_use_cache          = false;
    m_regenerate_normals = false;
    m_serial_access      = false;
    m_show_stats         = false;
    m_image_verification = false;

    m_configuration = std::string() +
        "importer=raw\n" +
        "input_file=" + m_filename + "\n";
}

//----------------------------------------------------------------------
Raw_heightfield_data_importer::Raw_heightfield_data_importer()
{
    m_filename           = "";
    m_use_cache          = false;
    m_regenerate_normals = false;
    m_serial_access      = false;
    m_show_stats         = false;
    m_image_verification = false;
}

//----------------------------------------------------------------------
Raw_heightfield_data_importer::~Raw_heightfield_data_importer()
{
    // empty
}

//----------------------------------------------------------------------
const char* Raw_heightfield_data_importer::get_configuration() const
{
    return m_configuration.c_str();
}

//----------------------------------------------------------------------
mi::Size Raw_heightfield_data_importer::estimate(
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
nv::index::IDistributed_data_subset* Raw_heightfield_data_importer::create(
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

    // Check the file magic
    Raw_heightfield_file_format_type file_type;
    if (!is_raw_heightfield_file_magic(m_filename, file_type))
    {
        ERROR_LOG << "Raw heightfield file has wrong magic: '" << m_filename << "'";
        return 0;               // wrong magic
    }

    // Open file
    FILE* p_file = fopen(m_filename.c_str(), "rb");
    if(p_file == 0)
    {
        ERROR_LOG << "Cannot open raw heightfield file: '" << m_filename << "'";
        return 0;               // failed to load
    }

    // read the header
    String_dict header_dict;
    if (!load_commented_header_part(p_file, header_dict))
    {
        fclose(p_file);
        return 0;
    }

    if (!check_header_consistency(header_dict, heightfield_size))
    {
        WARN_LOG << "Inconsistency detected between the given data and in the file header. "
                 << "This may cause a problem.";
    }

    // Create a heightfield patch
    mi::Float32*                             elevation_data     = 0;
    mi::math::Vector_struct<mi::Float32, 3>* normal_vector_data = 0;
    mi::base::Handle<nv::index::IRegular_heightfield_patch> heightfield_patch(
        factory->create<nv::index::IRegular_heightfield_patch>());
    if (!heightfield_patch.is_valid_interface())
    {
        ERROR_LOG << "Cannot create IRegular_heightfield_patch.";
        fclose(p_file);
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

        // Acquire storage for normals, unless they should be computed automatically by the library
        if (!m_regenerate_normals)
        {
            normal_vector_data = heightfield_patch->generate_normal_storage();

            if (normal_vector_data == 0)
            {
                elevation_data = 0; // invalidate
                ERROR_LOG << "Cannot generate normal storage.";
                return 0;
            }
        }
    }

    // Read the data to fill the data storage
    load_data(p_file,
              cache_file_name,
              file_type,
              heightfield_size,
              bounding_box,
              patch_bbox,
              elevation_data,
              normal_vector_data);

    heightfield_patch->retain(); // since the handle will be out of scope next.
    return heightfield_patch.get();
}

//----------------------------------------------------------------------
void Raw_heightfield_data_importer::serialize(
    mi::neuraylib::ISerializer * serializer) const
{
    mi::Uint32 nb_elements = mi::Uint32(m_filename.size());
    serializer->write(&nb_elements, 1);
    serializer->write(reinterpret_cast<const mi::Uint8*>(m_filename.c_str()), nb_elements);
    serializer->write(&m_heightfield_size.x, 2);
    serializer->write(&m_use_cache, 1);
    serializer->write(&m_regenerate_normals, 1);
    serializer->write(&m_serial_access, 1);
    serializer->write(&m_show_stats, 1);
    serializer->write(&m_image_verification, 1);
}

//----------------------------------------------------------------------
void Raw_heightfield_data_importer::deserialize(
    mi::neuraylib::IDeserializer * deserializer)
{
    mi::Uint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_filename.resize(nb_elements);
    deserializer->read(reinterpret_cast<mi::Uint8*>(&m_filename[0]), nb_elements);
    deserializer->read(&m_heightfield_size.x, 2);
    deserializer->read(&m_use_cache, 1);
    deserializer->read(&m_regenerate_normals, 1);
    deserializer->read(&m_serial_access, 1);
    deserializer->read(&m_show_stats, 1);
    deserializer->read(&m_image_verification, 1);
}

//----------------------------------------------------------------------
nv::index::IRegular_heightfield_patch* Raw_heightfield_data_importer::try_load_from_cache_file(
    const mi::math::Bbox_struct<mi::Sint32, 3>& bounding_box,
    const mi::Sint64                            grid_size,
    nv::index::IData_subset_factory*            factory,
    std::string&                                cache_file_name
    ) const
{
    assert(m_use_cache);

    std::ostringstream cache_file_name_str;
    cache_file_name_str << m_filename << "_"
                        << bounding_box.min.x << "." << bounding_box.max.x << "x"
                        << bounding_box.min.y << "." << bounding_box.max.y
                        << "_raw_importer.tmp";

	INFO_LOG << "Trying to load cache file: " << cache_file_name_str.str();
    cache_file_name = cache_file_name_str.str();

    std::ifstream cache_file;
    cache_file.open(cache_file_name.c_str(), std::ios_base::in | std::ios_base::binary);

    if(!cache_file)
    {
        INFO_LOG << "Heightfield cache file '" << cache_file_name << "' is not yet available.";
    }
    else
    {
        // Create a heightfield patch
        Raw_heightfield_file_format_type file_type = RHFFT_NONE;
        mi::Float32*                             elevation_data     = 0;
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

            elevation_data =
                heightfield_patch->generate_elevation_storage();
            if (elevation_data == 0)
            {
                ERROR_LOG << "Cannot generate elevation storage.";
                return 0;
            }

            // Acquire storage for normals, unless they should be computed automatically by the library
            if (!m_regenerate_normals)
            {
                normal_vector_data = heightfield_patch->generate_normal_storage();
                if (normal_vector_data == 0)
                {
                    ERROR_LOG << "Cannot generate normal storage.";
                    return 0;
                }
            }
        }

        {
            mi::Float64 time_start = nv::index_common::get_time();

            // get file type
            cache_file.read(reinterpret_cast<char*>(&file_type), sizeof(Raw_heightfield_file_format_type));

            if (file_type != RHFFT_INDEX_RAW_v1 && file_type != RHFFT_INDEX_RAW_v2)
            {
                ERROR_LOG << "Heightfield importer: unsupported file type version " << file_type
                          << " when reading from file '" << cache_file_name << "'";
                return 0;
            }

            // Fetch data for the entire patch at once
            cache_file.read(reinterpret_cast<char*>(elevation_data), grid_size * sizeof(mi::Float32));

            mi::Uint64 bytes_read = grid_size * sizeof(mi::Float32);

            if (normal_vector_data != 0)
            {
                cache_file.read(
                    reinterpret_cast<char*>(normal_vector_data),
                    grid_size * sizeof(mi::math::Vector_struct<mi::Float32, 3>));

                bytes_read += grid_size * sizeof(mi::math::Vector_struct<mi::Float32, 3>);
            }
            else
            {
                // Just seek over it
                cache_file.seekg(
                    grid_size * sizeof(mi::math::Vector_struct<mi::Float32, 3>),
                    std::ios_base::cur);
            }

            mi::math::Vector<mi::Float32, 2> height_range;
            cache_file.read(reinterpret_cast<char*>(&height_range), sizeof(height_range));
            bytes_read += sizeof(height_range);

            if (!cache_file)
            {
                ERROR_LOG << "Heightfield importer: error reading from file '" << cache_file_name << "'";
            }
            cache_file.close();

            const mi::Float64 time_end = nv::index_common::get_time();
            const mi::Float64 time_taken = time_end - time_start;
            s_total_read_time += time_taken;
            s_total_read_bytes += bytes_read;

            if (m_show_stats) 
            {
                output_heightfield_import_stat("Heightfield (cached): Read ", bytes_read, time_taken,
                                               s_total_read_bytes, s_total_read_time);
            }

            if (file_type == RHFFT_INDEX_RAW_v1)
            {
                // Holes are marked as -1 in the old format, need to convert data to new format
                static bool info_printed = false;
                if (!info_printed)
                {
                    WARN_LOG << "The cache file '" << cache_file_name << "' is using an old file format and "
                             << "needs to be converted during loading. For best performance, please regenerate "
                             << "(i.e. delete) the cache file. This message will only be printed once.";
                    info_printed = true;
                }

                for (mi::Sint64 i=0; i < grid_size; ++i)
                {
                    mi::Float32& value = elevation_data[i];
                    if (value == -1.f)
                    {
                        value = nv::index::IRegular_heightfield_patch::get_hole_value();
                    }
                }
            }
        }

        heightfield_patch->retain();
        return heightfield_patch.get();
    }
    return 0;
}

//----------------------------------------------------------------------
bool Raw_heightfield_data_importer::check_header_consistency(
    const String_dict& header_dict,
    const mi::math::Vector<mi::Sint64, 2> heightfield_size) const
{
    if (!header_dict.is_defined("size"))
    {
        ERROR_LOG << "Header has no 'size' item.";
        return false;
    }

    const mi::math::Vector<mi::Sint64, 2> size_in_header = get_vec_sint64_2(header_dict.get("size"));
    if (heightfield_size != size_in_header)
    {
        ERROR_LOG << "Inconsistent header (" << size_in_header << ") and given information ("
                  << heightfield_size << ").";
        return false;
    }

    return true;
}

//----------------------------------------------------------------------
void Raw_heightfield_data_importer::load_data(
    FILE*                                       p_file,
    const std::string&                          cache_file_name,
    Raw_heightfield_file_format_type            file_type,
    const mi::math::Vector<mi::Sint64, 2>&      heightfield_size,
    const mi::math::Bbox_struct<mi::Sint32, 3>& bounding_box,
    const mi::math::Bbox<mi::Sint64, 2>&        patch_bbox,
    mi::Float32*                                elevation_data,
    mi::math::Vector_struct<mi::Float32, 3>*    normal_vector_data) const
{
    assert(p_file != 0);

    const mi::Sint64 range_i = patch_bbox.max.x - patch_bbox.min.x;
    const mi::Sint64 range_j = patch_bbox.max.y - patch_bbox.min.y;
    const mi::Sint64 grid_size = range_i * range_j;

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

    if (normal_vector_data != 0)
    {
        for (mi::Sint64 j = 0; j < range_j; ++j)
        {
            const mi::Sint64 idx_file = ((patch_bbox.min.y + j) * heightfield_size.x) + patch_bbox.min.x;
            const mi::Sint64 idx_map  = j * range_i;

            fseek(p_file, normals_begin + idx_file * sizeof(mi::math::Vector_struct<mi::Float32, 3>), SEEK_SET);
            fread((void*)&(normal_vector_data[idx_map]), sizeof(mi::math::Vector_struct<mi::Float32, 3>), range_i, p_file);
            bytes_read += range_i * sizeof(mi::math::Vector_struct<mi::Float32, 3>);
        }
    }

    fclose(p_file);

    const mi::Float64 time_end = nv::index_common::get_time();

    // Determine extent of bounding box in k-dimension from height values
    mi::Float32 min_k = mi::base::numeric_traits<mi::Float32>::max();
    mi::Float32 max_k = mi::base::numeric_traits<mi::Float32>::negative_max();

    if (file_type == RHFFT_INDEX_RAW_v1)
    {
        for (mi::Sint64 i=0; i < grid_size; ++i)
        {
            mi::Float32& value = elevation_data[i];
            if (value > -1.0f)
            {
                if (value < min_k)
                {
                    min_k = value;
                }
                if (value > max_k)
                {
                    max_k = value;
                }
            }
            else
            {
                // Convert
                value = nv::index::IRegular_heightfield_patch::get_hole_value();
            }
        }
    }
    else if (file_type == RHFFT_INDEX_RAW_v2)
    {
        for (mi::Sint64 i=0; i < grid_size; ++i)
        {
            const mi::Float32 value = elevation_data[i];
            if (!nv::index::IRegular_heightfield_patch::is_hole(value))
            {
                if (value < min_k)
                {
                    min_k = value;
                }
                if (value > max_k)
                {
                    max_k = value;
                }
            }
        }
    }

    // k-min and k-max should be different to avoid missing planar heightfields
    if (max_k == min_k)
    {
        max_k = min_k + 1.0f;
    }

    const mi::math::Vector<mi::Float32, 2> height_range(min_k, max_k);

    INFO_LOG << "Loaded heightfield (raw) '" << m_filename << "' patch bounds: " << patch_bbox
             << " height range: " << height_range << " bbox: " << bounding_box
             << (file_type == RHFFT_INDEX_RAW_v1 ? " (old file format, not supporting negative height values)" : "");;

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
            INFO_LOG << "Heightfield cache file '" << cache_file_name << "' cannot be opened for writing.";
        }
        else
        {
            INFO_LOG << "Writing the heightfield patch data to a cache file '" << cache_file_name << "'";
            const Raw_heightfield_file_format_type file_type = RHFFT_INDEX_RAW_v2;
            cache_file.write(reinterpret_cast<const char*>(&file_type), sizeof(Raw_heightfield_file_format_type));

            cache_file.write(reinterpret_cast<char*>(elevation_data), grid_size * sizeof(mi::Float32));
            cache_file.write(reinterpret_cast<char*>(normal_vector_data), grid_size * sizeof(mi::math::Vector_struct<mi::Float32, 3>));
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
        std::ostringstream ppm_file;
        ppm_file << "heightfield_elevation_values_"
                 << bounding_box.min.x << "." << bounding_box.max.x << "x"
                 << bounding_box.min.y << "." << bounding_box.max.y
                 << "_raw_importer.ppm";
        nv::index_common::write_ppm(
            ppm_file.str(),
            elevation_values,
            height_range,
            static_cast<mi::Sint32>(range_i),
            static_cast<mi::Sint32>(range_j));
    }
}

//----------------------------------------------------------------------
//======================================================================
//----------------------------------------------------------------------
bool is_raw_heightfield_file_magic(
    const std::string&                infname,
    Raw_heightfield_file_format_type& file_type)
{
    const std::string raw_heightfield_magic_str = "index_raw_heightfield";
    mi::Sint32        file_type_int = -1;
    std::string       error_mes;
    const mi::Sint32  is_magic_ok = get_file_magic(
        infname, raw_heightfield_magic_str, file_type_int, &error_mes);
    if (is_magic_ok != 0)
    {
        if (is_ssv_file(infname))
        {
            // It seems old ssv file
            ERROR_LOG << "The file '" << infname << "' looks like  an old ssv raw heightfield file (before IndeX 1.3).\n"
                      << "This is deprecated and soon there will be no support anymore.\n"
                      << "You can try to use 'importer = raw_ssv' in your project file\n"
                      << "instead of 'importer = raw'.";
        }
        else
        {
            ERROR_LOG << "The file '" << infname << "' is not an NVIDIA IndeX raw heightfield file (wrong magic).\n"
                      << "The file header should start with '#! " << raw_heightfield_magic_str;
        }
        
        return false;
    }

    if (file_type_int == RHFFT_INDEX_RAW_v1 || file_type_int == RHFFT_INDEX_RAW_v2)
    {
        file_type = static_cast<Raw_heightfield_file_format_type>(file_type_int);
    }
    else
    {
        file_type = RHFFT_NONE;
    }

    if (file_type == RHFFT_NONE)
    {
        ERROR_LOG << "The file '" << infname << "' is not a raw heightfield file (wrong file type: '"
                  << file_type_int << "')\n";
        return false;
    }

    return true;
}

//----------------------------------------------------------------------
bool load_commented_header_part(FILE* p_file,
                                nv::index_common::String_dict & header_dict)
{
    if (p_file == 0)
    {
        ERROR_LOG << "invalid file pointer.";
        return false;
    }

    // max lines of the header part
    const mi::Sint32 max_line = 1024;
    const mi::Sint32 buf_size = 1024;
    char buf[buf_size];
    std::string header_str;

    // skip the first line (magic line)
    if (fgets(buf, buf_size, p_file) == 0)
    {
        ERROR_LOG << "cannot read the given file.";
        return false;
    }

    // read from the second line
    std::string end_marker = "#-end\n";
    mi::Sint32 line_count = 1;
    bool is_end = false;
    for (mi::Sint32 i = 0; i < max_line; ++i)
    {
        // can read?
        if (fgets(buf, buf_size, p_file) == 0)
        {
            ERROR_LOG << "cannot read the raw heightfield file at " << line_count << ". File truncated?";
            return false;
        }
        ++line_count;

        const std::string line(buf);

        // end?
        if(line == end_marker)
        {
            is_end = true;
            break;
        }

        // header?
        if(line[0] != '#')
        {
            ERROR_LOG << "Unexpected header line at " << line_count << ". '#' not found.";
            break;
        }
        assert(line.size() > 1);
        header_str += std::string(&(buf[1])); // skip the #
    }

    if (!is_end)
    {
        ERROR_LOG << "Cannot find the end marker (#-end) in the header section. Broken file?";
        return false;
    }

    header_dict.clear();
    std::stringstream sstr(header_str);
    header_dict.read(sstr);

    return true;
}

//----------------------------------------------------------------------
}} // namespace nv::index_common
