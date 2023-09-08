/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "synthetic_heightfield_generator.h"

#include <cstdio>
#include <cassert>

#include <nv/index/iregular_heightfield_patch.h>

#include "forwarding_logger.h"
#include "encode_height.h"
#include "common_utility.h"
#include "heightfield_normal.h"



namespace nv {
namespace index_common {

//----------------------------------------------------------------------
Synthetic_heightfield_generator::Synthetic_heightfield_generator()
{
    // empty
}

//----------------------------------------------------------------------
Synthetic_heightfield_generator::Synthetic_heightfield_generator(
    const std::string&                              technique,
    const std::string&                              technique_param,
    const mi::math::Vector_struct<mi::Uint32, 2>&   size)
  :
    m_technique(technique),
    m_technique_param(technique_param),
    m_size(size)
{
    m_configuration = std::string() +
        "importer=synthetic\n" +
        "synthetic_type=" + m_technique + "\n" +
        "parameter=" + m_technique_param + "\n";
}

//----------------------------------------------------------------------
Synthetic_heightfield_generator::~Synthetic_heightfield_generator()
{
    // empty
}

//----------------------------------------------------------------------
mi::Size Synthetic_heightfield_generator::estimate(
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
nv::index::IDistributed_data_subset* Synthetic_heightfield_generator::create(
    const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
    nv::index::IData_subset_factory*                factory,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{
    DEBUG_LOG << "Synthetic heightfield generation using generation technique [" 
              << m_technique 
              << "] for generating data in bounding box "
              << bounding_box;

    // for more than 0.5G heightfield data entries, use 64bit.
    mi::math::Bbox<mi::Sint64, 2> patch_raw_bbox(
        static_cast< mi::Sint64 >(bounding_box.min.x),
        static_cast< mi::Sint64 >(bounding_box.min.y),
        static_cast< mi::Sint64 >(bounding_box.max.x),
        static_cast< mi::Sint64 >(bounding_box.max.y));

    // Size of the heightfield patch in memory
    const mi::math::Vector<mi::Sint64, 2> patch_raw_dim = patch_raw_bbox.max - patch_raw_bbox.min;
    assert((patch_raw_dim.x > 0) && (patch_raw_dim.y > 0));

    // Create a heightfield patch
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

    const mi::Sint64 patch_raw_size = patch_raw_dim.x * patch_raw_dim.y;
    DEBUG_LOG << "Size of the heightfield: " << m_size
              << ", patch_raw_bbox: " << patch_raw_bbox
              << ", patch_raw_size: " << patch_raw_size;

    mi::Sint32 syn_type = get_height_synthetic_type(m_technique);
    if(syn_type == HSM_Count)
    {
        ERROR_LOG << "No such technique identifier for a synthetic heightfield, use default type.";
        syn_type = HSM_Default;
    }

    if(syn_type == HSM_Perlin_noise)
    {
        this->generate_height_perlin_noise(patch_raw_bbox, elevation_data);
    }
    else if(syn_type == HSM_0)
    {
        this->generate_height_0(patch_raw_bbox, elevation_data); // with height parameter
    }
    else
    {
        this->generate_height_no_parameter(patch_raw_bbox, elevation_data, syn_type);
    }

    // generate the normal based on the height
    generate_vertex_normal(patch_raw_bbox, elevation_data, normal_data);

    // Determine extent of bounding box in k-dimension from height values
    mi::Float32 min_k = mi::base::numeric_traits<mi::Float32>::max();
    mi::Float32 max_k = mi::base::numeric_traits<mi::Float32>::negative_max();
    for(mi::Sint64 i = 0; i < patch_raw_size; ++i)
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

    // Store bounding box results
    mi::math::Vector<mi::Float32, 2> height_range;
    height_range.x = min_k;
    height_range.y = max_k;

    heightfield_patch->retain(); // since the handle will be out of scope.
    return heightfield_patch.get();
}

//----------------------------------------------------------------------
mi::base::Uuid Synthetic_heightfield_generator::subset_id() const
{
    return nv::index::IRegular_heightfield_patch::IID();
}

//----------------------------------------------------------------------
const char* Synthetic_heightfield_generator::get_configuration() const
{
    return m_configuration.c_str();
}

//----------------------------------------------------------------------
void Synthetic_heightfield_generator::serialize(
    mi::neuraylib::ISerializer * serializer) const
{
    mi::Uint32 nb_elements = mi::Uint32(m_technique.size());
    serializer->write(&nb_elements, 1);
    serializer->write(reinterpret_cast<const mi::Uint8*>(m_technique.c_str()), nb_elements);

    nb_elements = mi::Uint32(m_technique_param.size());
    serializer->write(&nb_elements, 1);
    serializer->write(reinterpret_cast<const mi::Uint8*>(m_technique_param.c_str()), nb_elements);

    serializer->write(&m_size.x, 2);
}

//----------------------------------------------------------------------
void Synthetic_heightfield_generator::deserialize(
    mi::neuraylib::IDeserializer * deserializer)
{
    mi::Uint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_technique.resize(nb_elements);
    deserializer->read(reinterpret_cast<mi::Uint8*>(&m_technique[0]), nb_elements);

    nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_technique_param.resize(nb_elements);
    deserializer->read(reinterpret_cast<mi::Uint8*>(&m_technique_param[0]), nb_elements);

    deserializer->read(&m_size.x, 2);
}

//----------------------------------------------------------------------
void Synthetic_heightfield_generator::generate_height_no_parameter(
    const mi::math::Bbox<mi::Sint64, 2>&    patch_raw_bbox,
    mi::Float32 * const                     p_elevation_data,
    const mi::Sint32                        syn_type) const
{
    // Size of the heightfield patch in memory
    const mi::math::Vector<mi::Sint64, 2> patch_raw_dim = patch_raw_bbox.max - patch_raw_bbox.min;
    no_unused_variable_warning_please(patch_raw_dim);
    assert((patch_raw_dim.x > 0) && (patch_raw_dim.y > 0));

    nv::index_common::String_dict key_value_map;
    std::string err_mes;
    if (!nv::index_common::get_string_dict_from_inline_parameter_string(
            m_technique_param, key_value_map, err_mes))
    {
        ERROR_LOG << "Failed to parse the parameter strings.\n" << err_mes;
    }

    // Scale all height values
    mi::Float32 scale = 1.0f;
    if (key_value_map.is_defined("scale"))
    {
        scale = nv::index_common::get_float32(key_value_map.get("scale"));
    }

    mi::math::Vector<mi::Sint64, 2> ij = patch_raw_bbox.min;
    for(ij.x = patch_raw_bbox.min.x; ij.x < patch_raw_bbox.max.x; ++ij.x)
    {
        for(ij.y = patch_raw_bbox.min.y; ij.y < patch_raw_bbox.max.y; ++ij.y)
        {
            const mi::Sint64 ij_idx = get_heightfield_index(ij, patch_raw_bbox);
            p_elevation_data[ij_idx] = get_height_synthetic_value(ij, syn_type) * scale;
        }
    }

    // Optionally add embedded geometry (isolated points and connecting lines)
    const bool geometry = nv::index_common::get_bool(key_value_map.get("geometry", "0"));

    // Optionally put a hole in every patch
    if (nv::index_common::get_bool(key_value_map.get("hole", "0")))
    {
        const mi::math::Vector<mi::Sint64, 2> border = patch_raw_dim / 4;
        const mi::math::Vector<mi::Sint64, 2> center = patch_raw_bbox.min + (patch_raw_dim / 2);
        for (ij.x = patch_raw_bbox.min.x + border.x; ij.x < patch_raw_bbox.max.x - border.x ; ++ij.x)
        {
            for (ij.y = patch_raw_bbox.min.y + border.y; ij.y < patch_raw_bbox.max.y - border.y; ++ij.y)
            {
                if (geometry)
                {
                    if (ij.x == center.x)
                    {
                        // Connecting line
                        continue;
                    }
                    else if (ij.y == center.y && (ij.x % 4 == 0))
                    {
                        // Isolated point
                        continue;
                    }
                }
                const mi::Sint64 ij_idx = get_heightfield_index(ij, patch_raw_bbox);
                p_elevation_data[ij_idx] = nv::index::IRegular_heightfield_patch::get_hole_value();
            }
        }
    }
}

//----------------------------------------------------------------------
void Synthetic_heightfield_generator::generate_height_0(
    const mi::math::Bbox<mi::Sint64, 2>&    patch_raw_bbox,
    mi::Float32 * const                     p_heightfield_data) const
{
    const mi::math::Vector<mi::Sint64, 2> patch_raw_dim = patch_raw_bbox.max - patch_raw_bbox.min;
    no_unused_variable_warning_please(patch_raw_dim);
    assert((patch_raw_dim.x > 0) && (patch_raw_dim.y > 0));

    nv::index_common::String_dict key_value_map;
    std::string err_mes;
    if (!nv::index_common::get_string_dict_from_inline_parameter_string(
            m_technique_param, key_value_map, err_mes))
    {
        ERROR_LOG << "Failed to parse the parameter strings.\n" << err_mes;
    }

    // get "height"
    mi::Float32 height = 5.0f;
    if(!key_value_map.is_defined("height"))
    {
        WARN_LOG << "User didn't specified the parameter 'height'. use default: " << height;
    }
    else
    {
        height = nv::index_common::get_float32(key_value_map.get("height"));
        INFO_LOG << "Heightfield generator HSM_0: height: " << height;
    }

    mi::math::Vector<mi::Sint64, 2> ij = patch_raw_bbox.min;
    for (ij.x = patch_raw_bbox.min.x; ij.x < patch_raw_bbox.max.x; ++ij.x)
    {
        for (ij.y = patch_raw_bbox.min.y; ij.y < patch_raw_bbox.max.y; ++ij.y)
        {
            const mi::Sint64 ij_idx = get_heightfield_index(ij, patch_raw_bbox);
            p_heightfield_data[ij_idx] = height;
        }
    }
}


//----------------------------------------------------------------------
void Synthetic_heightfield_generator::generate_height_perlin_noise(
    const mi::math::Bbox<mi::Sint64, 2>&    patch_raw_bbox,
    mi::Float32 * const                     p_elevation_data) const
{
    const mi::Float32 time = 0.0f;
    const mi::Sint64  patch_unit_i = 512;    
    const mi::Uint32 terms = 5;
    const mi::math::Vector<mi::Float32, 4> turbulence_weight(1.0f, 1.0f, 1.0f, 1.0f);
    const bool abs_noise = false;
    const bool ridged = true;
    const mi::Float32 height_scale  = 0.5f * 250.0f;
    const mi::Float32 height_offset = 10.0f;

    // Size of the heightfield patch in memory
    const mi::math::Vector<mi::Sint64, 2> patch_raw_dim = patch_raw_bbox.max - patch_raw_bbox.min;
    no_unused_variable_warning_please(patch_raw_dim);
    assert((patch_raw_dim.x > 0) && (patch_raw_dim.y > 0));

    mi::math::Vector<mi::Sint64, 2> ij = patch_raw_bbox.min;
    for(ij.x = patch_raw_bbox.min.x; ij.x < patch_raw_bbox.max.x; ++ij.x)
    {
        for(ij.y = patch_raw_bbox.min.y; ij.y < patch_raw_bbox.max.y; ++ij.y)
        {
            const mi::Sint64 ij_idx = get_heightfield_index(ij, patch_raw_bbox);
            const mi::Float32 val = encode_height_perlin_noise(ij, patch_unit_i, time, terms, 
                                                         turbulence_weight, abs_noise, ridged);
            const mi::Float32 height = (val * height_scale) + height_offset;
            p_elevation_data[ij_idx] = height;
        }
    }
}

//----------------------------------------------------------------------
}} // namespace nv::index_common
