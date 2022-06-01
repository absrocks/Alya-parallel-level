/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Volume data generator.

#include <cassert>

#include <mi/dice.h>

#include <nv/index/iregular_volume.h>
#include <nv/index/iregular_volume_brick.h>

#include "synthetic_volume_generator.h"
#include "forwarding_logger.h"
#include "common_utility.h"
#include "encode_voxel.h"
#include "perlin_noise_utility.h"

namespace nv {
namespace index_common {

//----------------------------------------------------------------------
namespace // anonymous namespace to hide
{
inline mi::Uint32 encode_voxel_transparent_test(const mi::math::Vector<mi::Sint64, 3>& ijk)
{
    // different value for the neighbors
    const mi::Uint32 ret = static_cast<mi::Uint32>(11 * ijk.x) + 16;
    return ret;
}

inline void encode_voxel_rgba_cube(
    mi::math::Vector_struct<mi::Uint8, 4>&  voxel,
    const mi::math::Bbox<mi::Uint32, 3>&    partial_bounds,
    const mi::math::Vector<mi::Sint64, 3>&  ijk)
{
    mi::Float32 rf = static_cast<mi::Float32>(ijk.x - partial_bounds.min.x) /
        static_cast<mi::Float32>(partial_bounds.max.x - partial_bounds.min.x);
    mi::Float32 gf = static_cast<mi::Float32>(ijk.y - partial_bounds.min.y) /
        static_cast<mi::Float32>(partial_bounds.max.y - partial_bounds.min.y);
    mi::Float32 bf = static_cast<mi::Float32>(ijk.z - partial_bounds.min.z) /
        static_cast<mi::Float32>(partial_bounds.max.z - partial_bounds.min.z);
                        
    voxel.x = static_cast<mi::Uint8>(rf * 255.0f + 0.5f);
    voxel.y = static_cast<mi::Uint8>(gf * 255.0f + 0.5f);
    voxel.z = static_cast<mi::Uint8>(bf * 255.0f + 0.5f);
    voxel.w = 255u;
}

// FIXME: Some code duplication with encode_voxel.h
inline mi::Uint8 get_voxel_as_uint8(
    const mi::math::Vector< mi::Sint64, 3 >&    ijk,
    mi::Uint32                                  synthesis_method)
{
    if(synthesis_method == Synthetic_volume_generator::DEFAULT) // default: something fancy
        return static_cast< mi::Uint8 >((2 * ijk.x + 4 * ijk.y + ijk.z) % 255);
    else if (synthesis_method == Synthetic_volume_generator::IJ)
        return static_cast< mi::Uint8 >(encode_voxel_ij_in_8bit(ijk));
    else if (synthesis_method == Synthetic_volume_generator::JK)
        return static_cast< mi::Uint8 >(encode_voxel_jk_in_8bit(ijk));
    else if (synthesis_method == Synthetic_volume_generator::KI)
        return static_cast< mi::Uint8 >(encode_voxel_ki_in_8bit(ijk));
    else if (synthesis_method == Synthetic_volume_generator::I)
        return static_cast< mi::Uint8 >(encode_voxel_i_in_8bit(ijk));
    else if (synthesis_method == Synthetic_volume_generator::J)
        return static_cast< mi::Uint8 >(encode_voxel_j_in_8bit(ijk));
    else if (synthesis_method == Synthetic_volume_generator::K)
        return static_cast< mi::Uint8 >(encode_voxel_k_in_8bit(ijk));
    else if (synthesis_method == Synthetic_volume_generator::SPHERE_0)
        return static_cast< mi::Uint8 >(encode_volume_sphere_value(ijk));
    else if (synthesis_method == Synthetic_volume_generator::TRANSPARENT_TEST)
        return static_cast< mi::Uint8 >(encode_voxel_transparent_test(ijk));
    else if (synthesis_method == Synthetic_volume_generator::PERLIN_NOISE)
    {
        ERROR_LOG << "Empty parameter string: the PERLIN_NOISE sythesis technique requires parameter.";
        assert(false);
        return mi::Uint8(0);
    }

    assert(false); // should not happened
    return mi::Uint8(0);
}

inline mi::math::Vector_struct<mi::Uint8, 4> get_voxel_as_rgba(
    const mi::math::Bbox<mi::Uint32, 3>&        partial_bounds,
    const mi::math::Vector< mi::Sint64, 3 >&    ijk,
    mi::Uint32                                  synthesis_method)
{
    mi::math::Vector_struct<mi::Uint8, 4> voxel;
    
    switch(synthesis_method)
    {
        default:
            encode_voxel_rgba_cube(voxel, partial_bounds, ijk);
            break;
    }
    
    return voxel;
        
}

inline mi::Float32 get_voxel_as_float32(
    const mi::math::Vector< mi::Sint64, 3 >&    ijk,
    mi::Uint32                                  synthesis_method,
    const mi::math::Vector<mi::Float32, 2>&     range,
    bool                                        normalize)
{
    float value;
    if (synthesis_method == Synthetic_volume_generator::DEFAULT ||
        synthesis_method == Synthetic_volume_generator::I)
    {
        value = static_cast<mi::Float32>(ijk.x);
    }
    else if (synthesis_method == Synthetic_volume_generator::J)
    {
        value = static_cast<mi::Float32>(ijk.y);
    }
    else if (synthesis_method == Synthetic_volume_generator::K)
    {
        value = static_cast<mi::Float32>(ijk.z);
    }
    else
    {
        ERROR_LOG << "Unsupported synthesis method " << synthesis_method << " for float32";
        return 0.f;
    }

    if (normalize)
    {
        // Translate, scale and then normalize to [0, 1] using fmod()
        return std::fmod(range.x + value / (range.y - range.x), 1.f);
    }
    else
    {
        // Just translate and scale
        return range.x + value / (range.y - range.x);
    }
}

} // namespace

//----------------------------------------------------------------------
Synthetic_volume_generator::Synthetic_volume_generator(
    const mi::math::Bbox_struct<mi::Uint32, 3>& bounds,
    Voxel_format                                voxel_format,    
    const std::string&                          synthesis_method,
    const std::string&                          parameters)
    :
    m_whole_bbox(bounds),
    m_synthetic_type(SPHERE_0),
    m_voxel_format(voxel_format),
    m_parameter_str(parameters)
{
    m_synthetic_type = nv::index_common::get_voxel_synthetic_type(synthesis_method);
    if (m_synthetic_type == COUNT)
    {
        ERROR_LOG << "No such synthesis method [" << synthesis_method << "], will use default type.";
        m_synthetic_type = DEFAULT;
    }

    std::string voxel_format_str;
    switch (voxel_format)
    {
      case VOXEL_FORMAT_UINT_8:   voxel_format_str = "uint8";   break;
      case VOXEL_FORMAT_RGBA_8:   voxel_format_str = "rgba8";   break;
      case VOXEL_FORMAT_UINT_16:  voxel_format_str = "uint16";  break;
      case VOXEL_FORMAT_FLOAT_32: voxel_format_str = "float32"; break;
      default:                    voxel_format_str = "unknown"; break;
    }
    
    std::ostringstream s;
    s << "importer=synthetic\n"
      << "voxel_format=" << voxel_format_str << "\n"
      << "synthetic_type=" << synthesis_method << "\n"
      << "parameter=" << parameters << "\n";
    m_configuration = s.str(); 
}

Synthetic_volume_generator::Synthetic_volume_generator()
    :
    m_whole_bbox(),
    m_synthetic_type(0),
    m_voxel_format(VOXEL_FORMAT_UINT_8),
    m_parameter_str()
{
    // empty
}


//----------------------------------------------------------------------
mi::Size Synthetic_volume_generator::estimate(
    const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{
    // This is an 8-bit raw volume data importer:
    const mi::Size dx = bounding_box.max.x - bounding_box.min.x;
    const mi::Size dy = bounding_box.max.y - bounding_box.min.y;
    const mi::Size dz = bounding_box.max.z - bounding_box.min.z;
    const mi::Size volume_brick_size = dx * dy * dz;
    
    switch (m_voxel_format)
    {
      case VOXEL_FORMAT_UINT_8:
          return volume_brick_size*sizeof(mi::Uint8);

      case VOXEL_FORMAT_RGBA_8:
          return volume_brick_size*sizeof(mi::Uint32);

      case VOXEL_FORMAT_UINT_16:
          return volume_brick_size*sizeof(mi::Uint16);

      case VOXEL_FORMAT_FLOAT_32:
          return volume_brick_size*sizeof(mi::Float32);

      default:
          ERROR_LOG << "Voxel format " << m_voxel_format << " not supported by estimate()";
          return volume_brick_size*sizeof(mi::Uint8);
    }
}



//----------------------------------------------------------------------
void Synthetic_volume_generator::generate_synthetic_data(
    mi::Uint8*                              brick_data,
    const mi::math::Bbox<mi::Sint64, 3>&    raw_bbox) const
{
    mi::math::Vector< mi::Sint64, 3 > ijk(0, 0, 0);
    for (ijk[0] = raw_bbox.min[0]; ijk[0] < raw_bbox.max[0]; ++ijk[0])
    {
        for (ijk[1] = raw_bbox.min[1]; ijk[1] < raw_bbox.max[1]; ++ijk[1])
        {
            for (ijk[2] = raw_bbox.min[2]; ijk[2] < raw_bbox.max[2]; ++ijk[2])
            {
                const mi::Sint64 vol_idx = get_volume_index(ijk, raw_bbox);
                assert(vol_idx >= 0);
                assert(vol_idx <  raw_bbox.volume());
                brick_data[vol_idx] = get_voxel_as_uint8(ijk, m_synthetic_type);
            }
        }
    }
}

void Synthetic_volume_generator::generate_synthetic_data(
    mi::math::Vector_struct<mi::Uint8, 4>*  brick_data,
    const mi::math::Bbox<mi::Sint64, 3>&    raw_bbox) const
{
    mi::math::Vector< mi::Sint64, 3 > ijk(0, 0, 0);
    for (ijk[0] = raw_bbox.min[0]; ijk[0] < raw_bbox.max[0]; ++ijk[0])
    {
        for (ijk[1] = raw_bbox.min[1]; ijk[1] < raw_bbox.max[1]; ++ijk[1])
        {
            for (ijk[2] = raw_bbox.min[2]; ijk[2] < raw_bbox.max[2]; ++ijk[2])
            {
                const mi::Sint64 vol_idx = get_volume_index(ijk, raw_bbox);
                assert(vol_idx >= 0);
                assert(vol_idx <  raw_bbox.volume());
                brick_data[vol_idx] = get_voxel_as_rgba(m_whole_bbox, ijk, m_synthetic_type);
            }
        }
    }
}

void Synthetic_volume_generator::generate_synthetic_data(
    mi::Float32*                            brick_data,
    const mi::math::Bbox<mi::Sint64, 3>&    raw_bbox,
    const mi::math::Vector<mi::Float32, 2>& range,
    bool                                    normalize) const
{
    mi::math::Vector< mi::Sint64, 3 > ijk(0, 0, 0);
    for (ijk[0] = raw_bbox.min[0]; ijk[0] < raw_bbox.max[0]; ++ijk[0])
    {
        for (ijk[1] = raw_bbox.min[1]; ijk[1] < raw_bbox.max[1]; ++ijk[1])
        {
            for (ijk[2] = raw_bbox.min[2]; ijk[2] < raw_bbox.max[2]; ++ijk[2])
            {
                const mi::Sint64 vol_idx = get_volume_index(ijk, raw_bbox);
                assert(vol_idx >= 0);
                assert(vol_idx <  raw_bbox.volume());
                brick_data[vol_idx] = get_voxel_as_float32(ijk, m_synthetic_type, range, normalize);
            }
        }
    }
}

//----------------------------------------------------------------------
nv::index::IDistributed_data_subset* Synthetic_volume_generator::create(
    const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
    nv::index::IData_subset_factory*                factory,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{
    // The size for volumetric data in this spatial area.
    const mi::math::Bbox< mi::Sint64, 3 > raw_bbox(
        bounding_box.min.x, bounding_box.min.y, bounding_box.min.z,
        bounding_box.max.x, bounding_box.max.y, bounding_box.max.z);
    const mi::Sint64 volume_brick_size = raw_bbox.volume();

    INFO_LOG << "Synthetic volume generator create type: " << m_synthetic_type
             << ", bounding box: " << bounding_box;
    
    switch (m_voxel_format)
    {
        default:
            return 0;
            
        case VOXEL_FORMAT_UINT_8:
        {
            // Create a volume brick and allocate the storage.
            mi::base::Handle<nv::index::IRegular_volume_brick_uint8> volume_brick(
                factory->create<nv::index::IRegular_volume_brick_uint8>());
            if (!volume_brick.is_valid_interface()) {
                ERROR_LOG << "Cannot allocate a volume brick.";
                return 0;
            }

            mi::Uint8* brick_data = volume_brick->generate_voxel_storage(bounding_box);
            if (brick_data == 0)
            {
                ERROR_LOG << "Cannot generate a voxel storage.";
                return 0;
            }
            memset(brick_data, 0, volume_brick_size * sizeof(mi::Uint8));

            // If synthetic type needs parameters, we need to handle them.
            if(m_synthetic_type == PERLIN_NOISE)
            {
                generate_perlin_noise_volume<mi::Uint8>(brick_data, raw_bbox, dice_transaction);
            }
            // Fill volume with zeroes
            else if(m_synthetic_type == ZERO)
            {
                memset(brick_data, 0, volume_brick_size * sizeof(mi::Uint8));
            }
            else
            {
                // Synthesis volume data
                generate_synthetic_data(brick_data, raw_bbox);
            }

            volume_brick->retain();     // since the handle will be out of scope.
            return volume_brick.get();
        }
        
        case VOXEL_FORMAT_RGBA_8:
        {
            // Create a volume brick and allocate the storage.
            mi::base::Handle<nv::index::IRegular_volume_brick_rgba8ui> volume_brick(
                factory->create<nv::index::IRegular_volume_brick_rgba8ui>());
            if (!volume_brick.is_valid_interface())
            {
                ERROR_LOG << "Cannot allocate a volume brick.";
                return 0;
            }

            mi::math::Vector_struct<mi::Uint8, 4>* brick_data = volume_brick->generate_voxel_storage(bounding_box);
            if (brick_data == 0)
            {
                ERROR_LOG << "Cannot generate a voxel storage.";
                return 0;
            }
            
            // If synthetic type needs parameters, we need to handle them.
            if(m_synthetic_type == PERLIN_NOISE)
            {
                generate_perlin_noise_volume<mi::math::Vector_struct<mi::Uint8, 4> >(brick_data, raw_bbox, dice_transaction);
            }
            else if(m_synthetic_type == ZERO)
            {
                // Fill volume with zeroes
                memset(brick_data, 0, volume_brick_size * sizeof(mi::Uint32));
            }
            else
            {
                generate_synthetic_data(brick_data, raw_bbox);
            }

            volume_brick->retain();     // since the handle will be out of scope.
            return volume_brick.get();
        }

        case VOXEL_FORMAT_UINT_16:
        {
            // Create a volume brick and allocate the storage.
            mi::base::Handle<nv::index::IRegular_volume_brick_uint16> volume_brick(
                factory->create<nv::index::IRegular_volume_brick_uint16>());
            if (!volume_brick.is_valid_interface())
            {
                ERROR_LOG << "Cannot allocate a volume brick.";
                return 0;
            }

            mi::Uint16* brick_data = volume_brick->generate_voxel_storage(bounding_box);
            if (brick_data == 0)
            {
                ERROR_LOG << "Cannot generate a voxel storage.";
                return 0;
            }
            
            // If synthetic type needs parameters, we need to handle them.
            if(m_synthetic_type == PERLIN_NOISE)
            {
                generate_perlin_noise_volume<mi::Uint16>(brick_data, raw_bbox, dice_transaction);
            }
            else if(m_synthetic_type == ZERO)
            {
                // Fill volume with zeroes
                memset(brick_data, 0, volume_brick_size * sizeof(mi::Uint16));
            }
            else
            {
                ERROR_LOG << "Unsupported synthetic type for voxel format uint16.";
                return 0;
            }

            volume_brick->retain();     // since the handle will be out of scope.
            return volume_brick.get();
        }

        case VOXEL_FORMAT_FLOAT_32:
        {
            // Create a volume brick and allocate the storage.
            mi::base::Handle<nv::index::IRegular_volume_brick_float32> volume_brick(
                factory->create<nv::index::IRegular_volume_brick_float32>());
            if (!volume_brick.is_valid_interface())
            {
                ERROR_LOG << "Cannot allocate a volume brick.";
                return 0;
            }

            mi::Float32* brick_data = volume_brick->generate_voxel_storage(bounding_box);
            if (brick_data == 0)
            {
                ERROR_LOG << "Cannot generate a voxel storage.";
                return 0;
            }
            
            // If synthetic type needs parameters, we need to handle them.
            if(m_synthetic_type == PERLIN_NOISE)
            {
                generate_perlin_noise_volume<mi::Float32>(brick_data, raw_bbox, dice_transaction);
            }
            else if(m_synthetic_type == ZERO)
            {
                // Fill volume with float zeroes (0.0f)
                mi::Float32* brick_data_iter = brick_data;
                for (mi::Sint64 i = 0; i < volume_brick_size; ++i)
                {
                    (*brick_data_iter) = 0.0f;
                    ++brick_data_iter;
                }
            }
            else
            {
                String_dict param_dict;
                std::string err_mes;
                get_string_dict_from_inline_parameter_string(
                    m_parameter_str, param_dict, err_mes);

                const mi::math::Vector<mi::Float32, 2> range = get_vec_float32_2(
                    param_dict.get("range", "0 1"));
                const bool normalize = get_bool(param_dict.get("normalize", "yes"));

                generate_synthetic_data(brick_data, raw_bbox, range, normalize);
            }

            volume_brick->retain();     // since the handle will be out of scope.
            return volume_brick.get();
        }
    }
    
}

//----------------------------------------------------------------------
mi::base::Uuid Synthetic_volume_generator::subset_id() const
{
    switch (m_voxel_format)
    {
        default:
            ERROR_LOG << "Unknown voxel format.";
            return nv::index::IDistributed_data_subset::IID();
        
        case VOXEL_FORMAT_UINT_8:
            return nv::index::IRegular_volume_brick_uint8::IID();        
            
        case VOXEL_FORMAT_RGBA_8:
            return nv::index::IRegular_volume_brick_rgba8ui::IID();

        case VOXEL_FORMAT_UINT_16:
            return nv::index::IRegular_volume_brick_uint16::IID();        

        case VOXEL_FORMAT_FLOAT_32:
            return nv::index::IRegular_volume_brick_float32::IID();        
    }
    
}

//----------------------------------------------------------------------
void Synthetic_volume_generator::serialize(mi::neuraylib::ISerializer* serializer) const
{
    serializer->write(&m_whole_bbox.min.x, 6);
    serializer->write(&m_synthetic_type, 1);

    mi::Uint32 nb_elements = mi::Uint32(m_parameter_str.size());
    serializer->write(&nb_elements, 1);
    if (nb_elements > 0)
    {
        serializer->write(reinterpret_cast<const mi::Uint8*>(m_parameter_str.c_str()), nb_elements);
    }

    nb_elements = mi::Uint32(m_configuration.size());
    serializer->write(&nb_elements, 1);
    if (nb_elements > 0)
    {
        serializer->write(reinterpret_cast<const mi::Uint8*>(m_configuration.c_str()), nb_elements);
    }

    const mi::Uint32 voxel_fmt = m_voxel_format;
    serializer->write(&voxel_fmt, 1);
}

//----------------------------------------------------------------------
void Synthetic_volume_generator::deserialize(mi::neuraylib::IDeserializer* deserializer)
{
    deserializer->read(&m_whole_bbox.min.x, 6);
    deserializer->read(&m_synthetic_type, 1);

    mi::Uint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_parameter_str.resize(nb_elements);
    if (nb_elements > 0)
    {
        deserializer->read(reinterpret_cast<mi::Uint8*>(&m_parameter_str[0]), nb_elements);
    }

    deserializer->read(&nb_elements, 1);
    m_configuration.resize(nb_elements);
    if (nb_elements > 0)
    {
        deserializer->read(reinterpret_cast<mi::Uint8*>(&m_configuration[0]), nb_elements);
    }

    mi::Uint32 voxel_fmt = VOXEL_FORMAT_UINT_8;
    deserializer->read(&voxel_fmt, 1);
    m_voxel_format = Voxel_format(voxel_fmt);
}

//----------------------------------------------------------------------
const char* Synthetic_volume_generator::get_configuration() const
{
    return m_configuration.c_str();
}


}} // namespace nv::index_common
