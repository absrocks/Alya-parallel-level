/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "ppm_heightfield_importer.h"

#include <nv/index/iregular_heightfield_patch.h>

#include "forwarding_logger.h"
#include "ppm_io.h"

#include <cstdio>

namespace nv {
namespace index_common {
//----------------------------------------------------------------------
PPM_heightfield_importer::PPM_heightfield_importer()
    :
    m_filename(),
    m_scale(1.0f),
    m_offset(0.0f),
    m_binary_mask_filename(),
    m_size()
{
    // empty
}

//----------------------------------------------------------------------
PPM_heightfield_importer::PPM_heightfield_importer(
    const std::string&                                          filename,
    mi::Float32                                                 scale,
    mi::Float32                                                 offset,
    const std::string&                                          binary_mask_filename,
    const mi::math::Vector_struct<mi::Uint32, 2>&               size)
  :
    m_filename(filename),
    m_scale(scale),
    m_offset(offset),
    m_binary_mask_filename(binary_mask_filename),
    m_size(size)
{
    std::ostringstream s;
    s << "importer=ppm\n"
      << "input_file=" << filename << "\n"
      << "importer_binary_mask=" << binary_mask_filename << "\n"
      << "importer_scale=" << scale << "\n"
      << "importer_offset=" << offset << "\n";
    m_configuration = s.str();
}

//----------------------------------------------------------------------
PPM_heightfield_importer::~PPM_heightfield_importer()
{
    // empty
}

//----------------------------------------------------------------------
namespace
{
static inline mi::Float32 get_height(
    const mi::Float32*      elevation_data,
    mi::Sint64              xres,
    mi::Sint64              yres,
    mi::Sint64              x,
    mi::Sint64              y)
{
    if (x < 0)
    {
        x = 0;
    }
    else if (x >= xres)
    {
        x = xres-1;
    }

    if (y < 0)
    {
        y = 0;
    }
    else if (y >= yres)
    {
        y = yres-1;
    }

    return elevation_data[xres * y + x];
}
}

//----------------------------------------------------------------------
mi::Size PPM_heightfield_importer::estimate(
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
nv::index::IDistributed_data_subset* PPM_heightfield_importer::create(
    const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
    nv::index::IData_subset_factory*                factory,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{
    DEBUG_LOG << "PPM file based heightfield importer[" 
              << m_filename 
              << "] for generating a heightfield patch in bounding box "
              << bounding_box;

    bool use_binary_mask = (m_binary_mask_filename != "");

    FILE* b_file = NULL;
    if(use_binary_mask)
    {
        // open the heightfield data file
        // binary mask is represent the hole in heightfield. (When the value is 0, it is hole.)
        b_file = fopen(m_binary_mask_filename.c_str(), "rb");
        if(b_file == 0)
        {
            WARN_LOG << "Cannot open binary mask file: " << m_binary_mask_filename;
            use_binary_mask = false;
        }
    }
    if(use_binary_mask)
    {
        // read the pgm file header
        std::string pgm_filemagic;
        mi::Sint32 img_width  = -1;
        mi::Sint32 img_height = -1;
        mi::Sint32 img_depth  = -1;
        std::string error_mes;
        if(!load_ppm_sub_header(b_file, m_binary_mask_filename, pgm_filemagic,
                                img_width, img_height, img_depth,
                                error_mes)){
            ERROR_LOG << "Can not recognize pgm file header: " << error_mes;
            fclose(b_file);            
            b_file = 0;
            return 0;
        }

        // only binary formatted portable pixmaps (PPMs) are supported
        if(pgm_filemagic != "P5")
        {
            fclose(b_file);
            b_file = 0;
            ERROR_LOG << m_binary_mask_filename << " is not a binary formatted portable graymap (PGM) heightfield file.";
            use_binary_mask = false;
        }
        //         DEBUG_LOG << "heightfield image importer [" << m_binary_mask_filename << "][" << pgm_filemagic
        //                   << "] resolution (" << img_width << " " << img_height << " " << img_depth << ")";
    }

    // open the heightfield data file
    FILE * p_file = fopen(m_filename.c_str(), "rb");
    if(p_file == 0)
    {
        ERROR_LOG << "Cannot open heightfield map file: " << m_filename;
        return 0;
    }

    std::string ppm_filemagic;
    mi::Sint32 img_width  = -1;
    mi::Sint32 img_height = -1;
    mi::Sint32 img_depth  = -1;
    std::string error_mes;
    if(!load_ppm_sub_header(p_file, m_filename, ppm_filemagic,
                            img_width, img_height, img_depth,
                            error_mes)){
        ERROR_LOG << "Can not recognize ppm file header: " << error_mes;
        fclose(p_file);
        p_file = 0;
        return 0;
    }

    // only binary formatted portable pixmaps (PPMs) are supported
    if(ppm_filemagic != "P6")
    {
        ERROR_LOG << m_filename << " is not a binary formatted portable pixmap (PPM) heightfield file. ["
                  << ppm_filemagic << "]";
        fclose(p_file);
        p_file = 0;
        return 0;
    }

    // for more than 0.5G heightfield data entries, use 64bit.
    mi::math::Bbox<mi::Sint64, 2> patch_bbox(
        static_cast< mi::Sint64 >(bounding_box.min.x),
        static_cast< mi::Sint64 >(bounding_box.min.y),
        static_cast< mi::Sint64 >(bounding_box.max.x),
        static_cast< mi::Sint64 >(bounding_box.max.y));

    // for more than 0.5G heightfield data entries, use 64bit.
    const mi::math::Vector<mi::Sint64, 2> heightfield_size(static_cast<mi::Sint64>(m_size.x),
                                                           static_cast<mi::Sint64>(m_size.y));

    // Size of the heightfield patch in memory
    const mi::Sint64 range_i = patch_bbox.max.x - patch_bbox.min.x; // bounding box
    const mi::Sint64 range_j = patch_bbox.max.y - patch_bbox.min.y;

    // Create a patch and allocate data storage
    mi::base::Handle<nv::index::IRegular_heightfield_patch> heightfield_patch(
        factory->create<nv::index::IRegular_heightfield_patch>());
    if (!heightfield_patch.is_valid_interface())
    {
        ERROR_LOG << "Cannot create IRegular_heightfield_patch.";
        return 0;
    }

    const mi::math::Bbox_struct<mi::Sint32, 2> bbox_patch_rect_st = {
        bounding_box.min.x, bounding_box.min.y,
        bounding_box.max.x, bounding_box.max.y,
    };
    if (!heightfield_patch->initialize(bbox_patch_rect_st))
    {
        ERROR_LOG << "Fail to initialize heightfield patch";
        return 0;
    }

    mi::Float32* elevation_data = heightfield_patch->generate_elevation_storage();
    if (elevation_data == 0)
    {
        ERROR_LOG << "Cannot generate elevation storage.";
        return 0;
    }

    mi::math::Vector_struct<mi::Float32, 3>* normal_data = heightfield_patch->generate_normal_storage();
    if (normal_data == 0)
    {
        ERROR_LOG << "Cannot generate normal storage.";
        return 0;
    }
    
    const mi::Sint64 grid_size = range_i * range_j;
    mi::Uint8* binary_mask = new mi::Uint8[grid_size];
    mi::Uint8* intermediate_heightfield_data = new mi::Uint8[grid_size*3];
    assert(binary_mask != 0);
    assert(intermediate_heightfield_data != 0);

    const mi::Sint64 raw_data_begin = ftell(p_file);
    //     const mi::Sint64 normals_begin  = raw_data_begin +
    //         sizeof(mi::Float32) * heightfield_size.x * heightfield_size.y;

    if(use_binary_mask)
    {
        const mi::Sint64 mask_data_begin = ftell(b_file);
        for(mi::Sint64 j=0; j < range_j; ++j)
        {
            const mi::Sint64 idx_file = ((patch_bbox.min.y + j) * heightfield_size.x) + patch_bbox.min.x;
            const mi::Sint64 idx_map  = j * range_i;

            fseek(b_file, mask_data_begin + idx_file * sizeof(mi::Uint8), SEEK_SET);
            fread((void*)&(binary_mask[idx_map]), sizeof(mi::Uint8), range_i, b_file);
        }
        fclose(b_file);
    }

    for(mi::Sint64 j=0; j < range_j; ++j)
    {
        const mi::Sint64 idx_file = ((patch_bbox.min.y + j) * heightfield_size.x) + patch_bbox.min.x;
        const mi::Sint64 idx_map  = j * range_i * 3;

        fseek(p_file, raw_data_begin + idx_file * sizeof(mi::Uint8)*3, SEEK_SET);
        fread((void*)&(intermediate_heightfield_data[idx_map]), sizeof(mi::Uint8)*3, range_i, p_file);
    }
    fclose(p_file);

    // set the elevation values first ...
    for(mi::Sint64 i=0; i<grid_size; ++i)
    {
        mi::Float32 elevation_value = static_cast<mi::Float32>(intermediate_heightfield_data[i*3]);
        elevation_data[i] = m_offset + (m_scale * elevation_value);
    }
    delete[] intermediate_heightfield_data;
    intermediate_heightfield_data = 0;

    // ... then compute averaged normals ...
    mi::Uint32 counter = 0;
    for(mi::Sint64 j=0; j < range_j; j++)
    {
        for(mi::Sint64 i=0; i < range_i; i++)
        {

            const mi::Float32 d0 = get_height(elevation_data, range_i, range_j, i,   j  );
            const mi::Float32 d1 = get_height(elevation_data, range_i, range_j, i+1, j  );
            const mi::Float32 d2 = get_height(elevation_data, range_i, range_j, i+1, j+1);
            const mi::Float32 d3 = get_height(elevation_data, range_i, range_j, i,   j+1);
            const mi::Float32 d4 = get_height(elevation_data, range_i, range_j, i-1, j+1);
            const mi::Float32 d5 = get_height(elevation_data, range_i, range_j, i-1, j  );
            const mi::Float32 d6 = get_height(elevation_data, range_i, range_j, i-1, j-1);
            const mi::Float32 d7 = get_height(elevation_data, range_i, range_j, i,   j-1);
            const mi::Float32 d8 = get_height(elevation_data, range_i, range_j, i+1, j-1);

            const mi::math::Vector<mi::Float32, 3> v1( 1, 0,  d1-d0);
            const mi::math::Vector<mi::Float32, 3> v2( 1, 1,  d2-d0);
            const mi::math::Vector<mi::Float32, 3> v3( 0, 1,  d3-d0);
            const mi::math::Vector<mi::Float32, 3> v4(-1, 1,  d4-d0);
            const mi::math::Vector<mi::Float32, 3> v5(-1, 0,  d5-d0);
            const mi::math::Vector<mi::Float32, 3> v6(-1,-1,  d6-d0);
            const mi::math::Vector<mi::Float32, 3> v7( 0,-1,  d7-d0);
            const mi::math::Vector<mi::Float32, 3> v8( 1,-1,  d8-d0);

            mi::math::Vector<mi::Float32, 3> n = mi::math::cross(v1, v2);
            n += mi::math::cross(v2, v3);
            n += mi::math::cross(v3, v4);
            n += mi::math::cross(v4, v5);
            n += mi::math::cross(v5, v6);
            n += mi::math::cross(v6, v7);
            n += mi::math::cross(v7, v8);
            n += mi::math::cross(v8, v1);
            n.normalize();

            normal_data[counter] = n;
            ++counter;
        }
    }

    // ... and finally do the masking after the normal calculations to not harm the normal averaging!
    for(mi::Sint64 i=0; i<grid_size; ++i)
    {
        mi::Float32 elevation_value = static_cast<mi::Float32>(elevation_data[i]);

        if(use_binary_mask)
        {
            elevation_value = binary_mask[i] == 0 ? nv::index::IRegular_heightfield_patch::get_hole_value() : elevation_value;
        }
        else
        {
            if(elevation_value == 0.f)
            {
                elevation_value = nv::index::IRegular_heightfield_patch::get_hole_value();
            }
        }
        elevation_data[i] = elevation_value;
    }

    delete[] binary_mask;
    binary_mask = 0;

    // determine extent of bounding box in k-dimension from height values
    mi::Float32 min_k = mi::base::numeric_traits<mi::Float32>::max();
    mi::Float32 max_k = mi::base::numeric_traits<mi::Float32>::negative_max();
    for(mi::Sint64 i=0; i < grid_size; ++i)
    {
        const mi::Float32 value = elevation_data[i];
        if (!nv::index::IRegular_heightfield_patch::is_hole(value))
        {
            if(value < min_k)
            {
                min_k = value;
            }
            if(value > max_k)
            {
                max_k = value;
            }
        }
    }

    // k-min and k-max should be different to avoid missing planar heightfield
    if (max_k == min_k)
    {
        max_k = min_k + 1.f;
    }

    // store bounding box results
    mi::math::Vector<mi::Float32, 2> height_range;
    height_range.x = min_k;
    height_range.y = max_k;

    heightfield_patch->retain(); // since the handle is out of scope next.
    return heightfield_patch.get();
}
//----------------------------------------------------------------------
mi::base::Uuid PPM_heightfield_importer::subset_id() const
{
    return nv::index::IRegular_heightfield_patch::IID();
}

//----------------------------------------------------------------------
const char* PPM_heightfield_importer::get_configuration() const
{
    return m_configuration.c_str();
}

//----------------------------------------------------------------------
void PPM_heightfield_importer::serialize(
    mi::neuraylib::ISerializer * serializer) const
{
    mi::Uint32 nb_elements = mi::Uint32(m_filename.size());
    serializer->write(&nb_elements, 1);
    serializer->write(reinterpret_cast<const mi::Uint8*>(m_filename.c_str()), nb_elements);

    serializer->write(&m_scale, 1);
    serializer->write(&m_offset, 1);

    nb_elements = mi::Uint32(m_binary_mask_filename.size());
    serializer->write(&nb_elements, 1);
    serializer->write(reinterpret_cast<const mi::Uint8*>(m_binary_mask_filename.c_str()), nb_elements);

    serializer->write(&m_size.x, 2);
}

//----------------------------------------------------------------------
void PPM_heightfield_importer::deserialize(
    mi::neuraylib::IDeserializer * deserializer)
{
    mi::Uint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_filename.resize(nb_elements);
    deserializer->read(reinterpret_cast<mi::Uint8*>(&m_filename[0]), nb_elements);

    deserializer->read(&m_scale, 1);
    deserializer->read(&m_offset, 1);

    nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_binary_mask_filename.resize(nb_elements);
    deserializer->read(reinterpret_cast<mi::Uint8*>(&m_binary_mask_filename[0]), nb_elements);

    deserializer->read(&m_size.x, 2);
}

//----------------------------------------------------------------------
}} // namespace nv::index_common
