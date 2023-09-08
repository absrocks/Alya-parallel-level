/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief voxel value encoder for test

#ifndef NVIDIA_INDEX_BIN_COMMON_ENCODE_VOXEL_H
#define NVIDIA_INDEX_BIN_COMMON_ENCODE_VOXEL_H

#include <mi/dice.h>

#include <mi/base/types.h>
#include <mi/math/vector.h>
#include <mi/math/bbox.h>

#include "forwarding_logger.h"
#include "common_utility.h"
#include "perlin_noise_utility.h"

#include <string>
#include <cassert>
#include <iostream>

namespace nv {
namespace index_common {
//----------------------------------------------------------------------

/// Synthetic volume synthesis method
enum Volume_data_encode_e {
    /// default synthesis method
    VDE_Default,
    /// ijk encoding for i and j
    VDE_IJ,
    /// ijk encoding for j and k
    VDE_JK,
    /// ijk encoding for k and i
    VDE_KI,
    /// ijk encoding for i
    VDE_I,
    /// ijk encoding for j
    VDE_J,
    /// ijk encoding for k
    VDE_K,
    /// sphere type 0
    VDE_SPHERE_0,
    /// for a transparent test
    VDE_TRANSPARENT_TEST,
    /// Perlin noise (need parameters)
    VDE_PERLIN_NOISE,
    /// set all voxels to zero
    VDE_ZERO,
    /// sentinel
    VDE_Count
};

/// type of voxel (fixed, depends on library)
typedef mi::Uint8 VoxelType;

//----------------------------------------------------------------------
/// encode ij in 8 bit (the return value type is Sint32).
///
/// \param[in]  ijk ijk coordinates
/// \return     encoded ij in 8 bit
inline mi::Sint32 encode_voxel_ij_in_8bit(mi::math::Vector<mi::Sint64, 3> const & ijk)
{
    mi::Sint32 const ret = static_cast< mi::Sint32 >(((ijk.x & 15) << 4) + (ijk.y & 15));
    return ret;
}

//----------------------------------------------------------------------
/// encode jk in 8 bit (the return value type is Sint32).
///
/// \param[in]  ijk ijk coordinates
/// \return     encoded jk in 8 bit
inline mi::Sint32 encode_voxel_jk_in_8bit(mi::math::Vector<mi::Sint64, 3> const & ijk)
{
    mi::Sint32 const ret = static_cast< mi::Sint32 >(((ijk.y & 15) << 4) + (ijk.z & 15));
    return ret;
}

//----------------------------------------------------------------------
/// encode ki in 8 bit (the return value type is Sint32).
///
/// \param[in]  ijk ijk coordinates
/// \return     encoded ki in 8 bit
inline mi::Sint32 encode_voxel_ki_in_8bit(mi::math::Vector<mi::Sint64, 3> const & ijk)
{
    mi::Sint32 const ret = static_cast< mi::Sint32 >(((ijk.z & 15) << 4) + (ijk.x & 15));
    return ret;
}

//----------------------------------------------------------------------
/// encode i in 8 bit (the return value type is Sint32).
///
/// \param[in]  ijk ijk coordinates
/// \return     encoded i in 8 bit
inline mi::Sint32 encode_voxel_i_in_8bit(mi::math::Vector<mi::Sint64, 3> const & ijk)
{
    mi::Sint32 const ret = static_cast< mi::Sint32 >(ijk.x & 255);
    return ret;
}

//----------------------------------------------------------------------
/// encode j in 8 bit (the return value type is Sint32).
///
/// \param[in]  ijk ijk coordinates
/// \return     encoded j in 8 bit
inline mi::Sint32 encode_voxel_j_in_8bit(mi::math::Vector<mi::Sint64, 3> const & ijk)
{
    mi::Sint32 const ret = static_cast< mi::Sint32 >(ijk.y & 255);
    return ret;
}

//----------------------------------------------------------------------
/// encode k in 8 bit (the return value type is Sint32).
///
/// \param[in]  ijk ijk coordinates
/// \return     encoded k in 8 bit
inline mi::Sint32 encode_voxel_k_in_8bit(mi::math::Vector<mi::Sint64, 3> const & ijk)
{
    mi::Sint32 const ret = static_cast< mi::Sint32 >(ijk.z & 255);
    return ret;
}

//----------------------------------------------------------------------
/// encode spherical voxel value in 8 bit (the return value type is Sint32).
///
/// \param[in]  ijk ijk coordinates
/// \return     encoded sphere in 8 bit
inline mi::Sint32 encode_volume_sphere_value(const mi::math::Vector< mi::Sint64, 3 >& ijk)
{
    mi::Sint64  const radius   = 250;
    mi::Float64 const radiusf  = static_cast< mi::Float64 >(radius);
    mi::Sint64  const diameter = 2 * radius;
    mi::Float64 im = static_cast< mi::Float32 >(((ijk[0] % diameter) - radius)) / radiusf;
    mi::Float64 jm = static_cast< mi::Float32 >(((ijk[1] % diameter) - radius)) / radiusf;
    mi::Float64 km = static_cast< mi::Float32 >(((ijk[2] % diameter) - radius)) / radiusf;
    // 1.7320508 = \sqrt{3.0}
    mi::Float64 ampf = (200.0/1.7320508) * sqrt(im * im + jm * jm + km * km);

    return static_cast< VoxelType >(ampf);
}

//----------------------------------------------------------------------
/// encode voxel value for the transparent test
///
/// \param[in]  ijk ijk coordinates
/// \return     encoded voxel in 8 bit
inline mi::Sint32 encode_voxel_transparent_test(mi::math::Vector<mi::Sint64, 3> const & ijk)
{
    // different value for the neighbors
    mi::Sint32 const ret = static_cast< mi::Sint32 >(11 * ijk.x) + 16;
    return ret;
}

//----------------------------------------------------------------------
/// Encoding the Perlin noise values based on the given type to match various regular volume formats. 
template<typename T>
T encode_perlin_noise_float_value(mi::Float32 perlin_value)                  { T t; return t; }

/// Encoding the Perlin noise values for Uint8
template<>
inline mi::Uint8 encode_perlin_noise_float_value<mi::Uint8>(mi::Float32 perlin_value)
{
    return static_cast<mi::Uint8>((perlin_value+1.f) * 127.5f);
}

/// Encoding the Perlin noise values for Uint16
template<>
inline mi::Uint16 encode_perlin_noise_float_value<mi::Uint16>(mi::Float32 perlin_value)
{
    return static_cast<mi::Uint16>((perlin_value+1.f) * 0.5);
}

/// Encoding the Perlin noise values for Float32
template<>
inline mi::Float32 encode_perlin_noise_float_value<mi::Float32>(mi::Float32 perlin_value)
{
    return ((perlin_value+1.f) * 0.5f);
}

/// Encoding the Perlin noise values for Vector_struct<mi::Uint8, 4>
template<>
inline mi::math::Vector_struct<mi::Uint8, 4> encode_perlin_noise_float_value<mi::math::Vector_struct<mi::Uint8, 4> >(mi::Float32 perlin_value)
{
    perlin_value = ((perlin_value+1.f) * 127.5f);
    mi::Uint8 v = static_cast<mi::Uint8>(perlin_value);

    mi::math::Vector_struct<mi::Uint8, 4> value;
    value.x = v;
    value.y = v;
    value.z = v;
    value.w = v;

    return value;
}
//----------------------------------------------------------------------

//----------------------------------------------------------------------
/// Encode Perlin noise voxel value.
///
/// Inline parameter string
///  "cube_unit=512,time=0,terms=5,turbulence_weight=1 1 1 1,abs_noise=0,ridged=0"
///
/// \param[in] ijk ijk coordinates
/// \param[in] cube_unit_i repeating Perlin noise cube
/// unit. Recommended to set this to the volume size
/// \param[in] time Perlin noise generator time (0.0 for static)
/// \param[in] terms how many frequency terms are summed up
/// \param[in] turbulence_weight 
/// \param[in] abs_noise
/// \param[in] ridged
/// \return encoded Perlin noise in 8 bit
template<typename T>
inline T encode_perlin_noise_value(
    const mi::math::Vector< mi::Sint64, 3 >& ijk,
    mi::Sint64 cube_unit_i,
    mi::Float32 time,
    mi::Uint32 terms,
    const mi::math::Vector< mi::Float32, 4 >& turbulence_weight,
    bool abs_noise,
    bool ridged)
{
    const mi::Float64 cube_unit_f = static_cast< mi::Float64 >(cube_unit_i);

    const mi::math::Vector< mi::Float32, 3 > pos(
        static_cast< mi::Float32 >(static_cast< mi::Float64 >(ijk.x % cube_unit_i) / cube_unit_f),
        static_cast< mi::Float32 >(static_cast< mi::Float64 >(ijk.y % cube_unit_i) / cube_unit_f),
        static_cast< mi::Float32 >(static_cast< mi::Float64 >(ijk.z % cube_unit_i) / cube_unit_f));

    const mi::Float32 value = summed_perlin_noise(pos, time, terms, turbulence_weight, abs_noise, ridged);

    assert(-1.0f <= value);
    assert(value <= 1.0f);
    return encode_perlin_noise_float_value<T>(value);
}

/// Parallel perlin noise
template<typename T>
class Parallel_perlin_noise :
        public mi::neuraylib::Fragmented_job<0x1d29a35f,0x1786,0x4558,0x90,0x48,0x73,0x88,0x9d,0x8f,0x66,0x30>
{
public:
    // Need this typedef to work around an issue with commas in default arguments in gcc versions
    // earlier than 4.4, see https://gcc.gnu.org/bugzilla/show_bug.cgi?id=57
    typedef mi::math::Vector<mi::Float32, 4> Vector_Float32;

    Parallel_perlin_noise(
        T*                                      brick_data,
        const mi::math::Bbox< mi::Sint64, 3 >&  brick_bbox,
        mi::Sint64                              cube_unit         = 512,
        mi::Float32                             time              = 0.0f,
        mi::Uint32                              terms             = 5,
        mi::math::Vector< mi::Float32, 4>       turbulence_weight = Vector_Float32(1.0f, 1.0f, 1.0f, 1.0f),
        bool                                    abs_noise         = false,
        bool                                    ridged            = false
    ) 
    {
        m_brick_data        = brick_data;
        m_brick_bbox        = brick_bbox;
        m_cube_unit         = cube_unit;
        m_time              = time;
        m_terms             = terms;
        m_turbulence_weight = turbulence_weight;
        m_abs_noise         = abs_noise;
        m_ridged            = ridged;
    }

    /// Default or serialization only
    Parallel_perlin_noise() {}

    virtual void  execute_fragment(
        mi::neuraylib::IDice_transaction*               dice_transaction,
        mi::Size                                        index,
        mi::Size                                        count,
        const mi::neuraylib::IJob_execution_context*    context)
    {
        // INFO_LOG << "Launching fragment: " << index << "/" << count;
        mi::math::Vector< mi::Sint64, 3 > ijk(m_brick_bbox.min[0]+index, 0, 0);

        for (ijk[1] = m_brick_bbox.min[1]; ijk[1] < m_brick_bbox.max[1]; ++ijk[1])
        {
            for (ijk[2] = m_brick_bbox.min[2]; ijk[2] < m_brick_bbox.max[2]; ++ijk[2])
            {
                const mi::Sint64 vol_idx = get_volume_index(ijk, m_brick_bbox);
                assert(vol_idx >= 0);
                assert(vol_idx <  m_brick_bbox.volume());
                m_brick_data[vol_idx] = nv::index_common::encode_perlin_noise_value<T>(
                    ijk, m_cube_unit, m_time, m_terms, m_turbulence_weight, m_abs_noise, m_ridged);
            }
        }
    }

private:
    T*                                m_brick_data;
    mi::math::Bbox< mi::Sint64, 3>    m_brick_bbox;
    mi::Sint64                        m_cube_unit;
    mi::Float32                       m_time;
    mi::Uint32                        m_terms;
    mi::math::Vector<mi::Float32, 4>  m_turbulence_weight;
    bool                              m_abs_noise;
    bool                              m_ridged;
};

//----------------------------------------------------------------------
/// get synthesis method name
///
/// \param[in] synthesis_method  synthesis method identifier.
/// \return synthesis name
inline std::string get_voxel_synthetic_name(mi::Sint32 synthesis_method)
{
    char const * const p_name[11] = {
        "default",
        "ij",
        "jk",
        "ki",
        "i",
        "j",
        "k",
        "sphere_0",
        "transparent_test",
        "perlin_noise",
        NULL,
    };
    assert(synthesis_method >= 0);
    assert(synthesis_method <  9);

    return std::string(p_name[synthesis_method]);
}

//----------------------------------------------------------------------
/// get synthesis method ID
///
/// \param[in] synthesis_name   synthesis name
/// \return synthesis type enum. \see Synthetic_volume_generator.
inline mi::Sint32  get_voxel_synthetic_type(std::string const & synthesis_name)
{
    if (synthesis_name == "default")
    {
        return VDE_Default;
    }
    else if (synthesis_name == "ij")
    {
        return VDE_IJ;
    }
    else if (synthesis_name == "jk")
    {
        return VDE_JK;
    }
    else if (synthesis_name == "ki")
    {
        return VDE_KI;
    }
    else if (synthesis_name == "i"){
        return VDE_I;
    }
    else if (synthesis_name == "j")
    {
        return VDE_J;
    }
    else if (synthesis_name == "k")
    {
        return VDE_K;
    }
    else if (synthesis_name == "sphere_0")
    {
        return VDE_SPHERE_0;
    }
    else if (synthesis_name == "transparent_test")
    {
        return VDE_TRANSPARENT_TEST;
    }
    else if (synthesis_name == "perlin_noise")
    {
        return VDE_PERLIN_NOISE;
    }
    else if (synthesis_name == "zero")
    {
        return VDE_ZERO;
    }
    else{
        return VDE_Count;
    }
}

//----------------------------------------------------------------------
/// compute voxel value based on ijk coordinates.
inline VoxelType get_voxel_synthetic_value(
    mi::math::Vector< mi::Sint64, 3 > const & ijk,
    mi::Sint32 synthesis_method)
{
    // no switch here.
    if (synthesis_method == VDE_Default)
    {
        // default: something fancy
        return static_cast< VoxelType >((2 * ijk.x + 4 * ijk.y + ijk.z) % 255);
    }
    else if (synthesis_method == VDE_IJ)
    {
        return static_cast< VoxelType >(encode_voxel_ij_in_8bit(ijk));
    }
    else if (synthesis_method == VDE_JK)
    {
        return static_cast< VoxelType >(encode_voxel_jk_in_8bit(ijk));
    }
    else if (synthesis_method == VDE_KI)
    {
        return static_cast< VoxelType >(encode_voxel_ki_in_8bit(ijk));
    }
    else if (synthesis_method == VDE_I)
    {
        return static_cast< VoxelType >(encode_voxel_i_in_8bit(ijk));
    }
    else if (synthesis_method == VDE_J)
    {
        return static_cast< VoxelType >(encode_voxel_j_in_8bit(ijk));
    }
    else if (synthesis_method == VDE_K)
    {
        return static_cast< VoxelType >(encode_voxel_k_in_8bit(ijk));
    }
    else if (synthesis_method == VDE_SPHERE_0)
    {
        return static_cast< VoxelType >(encode_volume_sphere_value(ijk));
    }
    else if (synthesis_method == VDE_TRANSPARENT_TEST)
    {
        return static_cast< VoxelType >(encode_voxel_transparent_test(ijk));
    }
    else if (synthesis_method == VDE_PERLIN_NOISE)
    {
        std::cerr << "Internal Error! VDE_PERLIN_NOISE needs args, "
                  << "Should not be called by this." << std::endl;
        assert(false);
        return VoxelType(0);
    }
    else if (synthesis_method == VDE_ZERO)
    {
        std::cerr << "Internal Error! VDE_ZERO should not be called by this." << std::endl;
        assert(false);
        return VoxelType(0);
    }
    assert(false); // should not happened
    return VoxelType(0);
}

//----------------------------------------------------------------------
/// test voxel value generated by synthetic volume generator
///
/// \param[in] mes message for this test
/// \param[in] raw_bbox   raw bounding box (memory image bounding box)
/// \param[in] valid_bbox value valid bounding box
/// \param[in] ijk_encode_type ijk encode type. \see Synthetic_volume_generator
/// \param[in] p_data voxel data
/// \return true when all the voxel values are expected.
inline bool test_synthetic_voxel_value(
    std::string const & mes,
    mi::math::Bbox< mi::Sint64, 3 > const & raw_bbox,
    mi::math::Bbox< mi::Sint64, 3 > const & valid_bbox,
    mi::Sint32 ijk_encode_type,
    VoxelType const * const p_data
    )
{
    // INFO_LOG << "test_synthetic_voxel_value: " << mes;
    bool is_ok = true;
    mi::math::Vector< mi::Sint64, 3 > ijk = valid_bbox.min;
    for (ijk[0] = valid_bbox.min[0]; ijk[0] < valid_bbox.max[0]; ++ijk[0])
    {
        for (ijk[1] = valid_bbox.min[1]; ijk[1] < valid_bbox.max[1]; ++ijk[1])
        {
            for (ijk[2] = valid_bbox.min[2]; ijk[2] < valid_bbox.max[2]; ++ijk[2])
            {
                VoxelType expected_vval = get_voxel_synthetic_value(ijk, ijk_encode_type);
                const mi::Sint64 vol_idx = nv::index_common::get_volume_index(ijk, raw_bbox);
                const VoxelType  vval    = p_data[vol_idx];
                if (vval != expected_vval)
                {
                    ERROR_LOG << "ijk = " << ijk << ", idx = " << vol_idx
                              << ", access = " << static_cast< mi::Sint32 >(vval)
                              << " != "        << static_cast< mi::Sint32 >(expected_vval)
                              << " (expected).";
                    is_ok = false;
                }
            }
        }
    }
    // INFO_LOG << "test_synthetic_voxel_value: " << mes << "... check done";
    return is_ok;
}

//----------------------------------------------------------------------
/// partial volume copy. p_dst_data and p_src_data should be different
/// volumes.
///
/// \param[in,out] p_dst_data    destination volume data
/// \param[in]    dst_raw_bbox  destination data bounding box
/// \param[in]    p_src_data    source volume data
/// \param[in]    src_raw_bbox  source data buonding box
/// \param[in]    copy_src_bbox copy source volume region
/// \param[in]    dst_origin    destination min location position.
/// \return true when succeeded.
inline bool copy_partial_volume(
    VoxelType       *       p_dst_data,
    mi::math::Bbox< mi::Sint64, 3 >   const & dst_raw_bbox,
    VoxelType const * const p_src_data,
    mi::math::Bbox< mi::Sint64, 3 >   const & src_raw_bbox,
    mi::math::Bbox< mi::Sint64, 3 >   const & copy_src_bbox,
    mi::math::Vector< mi::Sint64, 3 > const & dst_origin
    )
{
    // check the parameter validity
    if (!dst_raw_bbox.is_volume())
    {
        ERROR_LOG << "illegal destination volume. dst_raw_bbox = " << dst_raw_bbox;
        return false;
    }
    if (!src_raw_bbox.is_volume())
    {
        ERROR_LOG << "illegal source volume. src_raw_bbox = " << src_raw_bbox;
        return false;
    }
    if (!copy_src_bbox.is_volume())
    {
        ERROR_LOG << "illegal copy region. copy_src_bbox = " << copy_src_bbox;
        return false;
    }
    if (!bbox_contains(src_raw_bbox, copy_src_bbox))
    {
        ERROR_LOG << "copy region is out of source range. src_raw_bbox = " << src_raw_bbox
                  << ", copy_src_bbox = " << copy_src_bbox;
        return false;
    }
    if (!dst_raw_bbox.contains(dst_origin))
    {
        ERROR_LOG << "destination origin is out of range. dst_raw_bbox = " << dst_raw_bbox
                  << ", dst_origin = " << dst_origin;
        return false;
    }
    mi::math::Vector< mi::Sint64, 3 > const dst_size_vec3  = dst_raw_bbox.max  - dst_raw_bbox.min;
    mi::math::Vector< mi::Sint64, 3 > const src_size_vec3  = src_raw_bbox.max  - src_raw_bbox.min;
    mi::math::Vector< mi::Sint64, 3 > const copy_size_vec3 = copy_src_bbox.max - copy_src_bbox.min;
    mi::math::Bbox< mi::Sint64, 3 >   const copy_dst_bbox(dst_origin, dst_origin + copy_size_vec3);
    if (!bbox_contains(dst_raw_bbox, copy_dst_bbox))
    {
        ERROR_LOG << "destination region is out of destination range. dst_raw_bbox = "
                  << dst_raw_bbox << ", copy_dst_bbox = " << copy_dst_bbox;
        return false;
    }
    mi::math::Vector< mi::Sint64, 3 > const offset_to_src_ijk =
        dst_origin - copy_src_bbox.min;

    mi::Sint64 const src_size  = src_size_vec3[0]  * src_size_vec3[1]  * src_size_vec3[2];
    mi::Sint64 const dst_size  = dst_size_vec3[0]  * dst_size_vec3[1]  * dst_size_vec3[2];
    mi::Sint64 const copy_size = copy_size_vec3[0] * copy_size_vec3[1] * copy_size_vec3[2];

    assert(src_size  == src_raw_bbox.volume());
    assert(dst_size  == dst_raw_bbox.volume());  no_unused_variable_warning_please(dst_size);
    assert(copy_size == copy_src_bbox.volume()); no_unused_variable_warning_please(copy_size);

    mi::math::Vector< mi::Sint64, 3 > ijk = copy_src_bbox.min; // for ijk[2] initialization
    for (ijk[0] = copy_src_bbox.min[0]; ijk[0] < copy_src_bbox.max[0]; ++(ijk[0]))
    {
        for (ijk[1] = copy_src_bbox.min[1]; ijk[1] < copy_src_bbox.max[1]; ++(ijk[1]))
        {
            mi::Sint64 const src_index = get_volume_index(ijk, src_raw_bbox);
            assert(src_index >= 0);
            if ((src_index + copy_size_vec3[2] - 1) >= src_size)
            {
                ERROR_LOG << "ijk = " << ijk
                          << ", src_index = "     << src_index
                          << ", src_raw_bbox = "  << src_raw_bbox
                          << ", copy_src_bbox = " << copy_src_bbox
                          << ", (src_index + copy_size_vec3[2] - 1) = "
                          << (src_index + copy_size_vec3[2] - 1)
                          << ", src_size = "      << src_size;
            }
            assert((src_index + copy_size_vec3[2] - 1) < src_size);

            mi::math::Vector< mi::Sint64, 3 > const dst_ijk = ijk + offset_to_src_ijk;

            mi::Sint64 const dst_index = get_volume_index(dst_ijk, dst_raw_bbox);
            assert(0 <= dst_index);
            assert((dst_index + copy_size_vec3[2] - 1) < dst_size);

            VoxelType const * const p_src = &(p_src_data[src_index]);
            VoxelType       *       p_dst = &(p_dst_data[dst_index]);
            memcpy(p_dst, p_src, sizeof(VoxelType) * copy_size_vec3[2]);
        }
    }

    return true;
}

//----------------------------------------------------------------------
/// expand border slices by copying valid voxels
///
/// \param[in,out] p_dst_data    destination volume data
/// \param[in] raw_bbox   raw (memory image) bounding box
/// \param[in] valid_bbox valid data bbox
/// \return    expanded result bbox
inline mi::math::Bbox< mi::Sint64, 3 > expand_border_by_copy(
    VoxelType * p_dst_data,
    mi::math::Bbox< mi::Sint64, 3 > const & raw_bbox,
    mi::math::Bbox< mi::Sint64, 3 > const & valid_bbox)
{
    // looking for the 6 boundaries
    // char const * const p_dir[] = {"i", "j", "k", 0, };

    // expand bounding box
    mi::math::Bbox< mi::Sint64, 3 > expand_bbox = valid_bbox;

    for (mi::Sint32 i = 0; i < 3; ++i)
    {
        if (raw_bbox.min[i] < expand_bbox.min[i])
        {
            // duplication source slice
            mi::math::Bbox< mi::Sint64, 3 > src_bbox = expand_bbox;
            src_bbox.max[i] = expand_bbox.min[i] + 1;
            mi::math::Vector< mi::Sint64, 3 > dst_origin = expand_bbox.min;
            dst_origin[i]  -= 1;

            // DEBUG_LOG << "copying boundary min." << p_dir[i]
            //           << ", raw = " << raw_bbox << ", expand_bbox = " << expand_bbox
            //           << "\n"
            //           << "src = "   << src_bbox
            //           << ", dst_origin = " << dst_origin;
            bool const success = copy_partial_volume(p_dst_data, raw_bbox,
                                                     p_dst_data, raw_bbox,
                                                     src_bbox,   dst_origin);
            if (!success) assert(0);

            // update expand bbox
            expand_bbox.min[i] -= 1;
        }

        if (raw_bbox.max[i] > expand_bbox.max[i])
        {
            // duplication source slice
            mi::math::Bbox< mi::Sint64, 3 > src_bbox = expand_bbox;
            src_bbox.min[i] = expand_bbox.max[i] - 1;

            mi::math::Vector< mi::Sint64, 3 > dst_origin = expand_bbox.min;
            dst_origin[i] = expand_bbox.max[i];

            // DEBUG_LOG << "copying boundary max." << p_dir[i]
            //           << ", raw = " << raw_bbox << ", expand_bbox = " << expand_bbox
            //           << "\n"
            //           << "src = "   << src_bbox
            //           << ", dst_origin = " << dst_origin;
            bool const success = copy_partial_volume(p_dst_data, raw_bbox,
                                                     p_dst_data, raw_bbox,
                                                     src_bbox,   dst_origin);
            if (!success) assert(0);

            // update expand bbox
            expand_bbox.max[i] += 1;
        }
    }

    return expand_bbox;
}

//----------------------------------------------------------------------
}} // namespace nv::index_common

#endif // NVIDIA_INDEX_BIN_COMMON_ENCODE_VOXEL_H
