/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief height value encoder for test

#ifndef NVIDIA_INDEX_BIN_COMMON_ENCODE_HEIGHT_H
#define NVIDIA_INDEX_BIN_COMMON_ENCODE_HEIGHT_H

#include <mi/base/types.h>
#include <mi/math/vector.h>
#include <mi/math/bbox.h>

#include "forwarding_logger.h"
#include "common_utility.h"
#include "perlin_noise_utility.h"

#include <cassert>
#include <string>

//----------------------------------------------------------------------
/// Synthetic height value synthesis method
enum Heightfield_synthesis_method {
    /// default synthesis method
    HSM_Default,
    /// ij encoding for i and j
    HSM_IJ,
    /// ij encoding for i
    HSM_I,
    /// ij encoding for j
    HSM_J,
    /// Perlin noise
    HSM_Perlin_noise,
    /// 0 (flat)
    HSM_0,
    /// sentinel
    HSM_Count
};

/// type of height value (fixed, depends on library)
typedef mi::Float32 HeightType;

//----------------------------------------------------------------------
/// encode height ij in [0,255]
///
/// \param[in]  ij ij coordinates
/// \return     encoded height at ij in [0, 255]
inline HeightType encode_height_ij(mi::math::Vector<mi::Sint64, 2> const & ij)
{
    HeightType ret = static_cast< HeightType >((ij.x + ij.y) & 255);
    return ret;
}

//----------------------------------------------------------------------
/// encode height i in [0,255]
///
/// \param[in]  ij ij coordinates
/// \return     encoded height at i in [0, 255]
inline HeightType encode_height_i(mi::math::Vector<mi::Sint64, 2> const & ij)
{
    HeightType const ret = static_cast< HeightType >(ij.x & 255);
    return ret;
}

//----------------------------------------------------------------------
/// encode height j in [0,255]
///
/// \param[in]  ij ij coordinates
/// \return     encoded height at j in [0, 255]
inline HeightType encode_height_j(mi::math::Vector<mi::Sint64, 2> const & ij)
{
    HeightType const ret = static_cast< HeightType >(ij.y & 255);
    return ret;
}

//----------------------------------------------------------------------
/// encode Perlin noise height value in HeightType
///
/// Inline parameter string
///  "patch_unit=512,time=0,terms=5,turbulence_weight=1 1 1 1,abs_noise=0,ridged=0"
///
/// \param[in] ij ij coordinates
/// \param[in] patch_unit_i repeating Perlin noise cube
/// unit. Recommended to set this to the volume size
/// \param[in] time Perlin noise generator time (0.0 for static)
/// \param[in] terms how many frequency terms are summed up
/// \param[in] turbulence_weight 
/// \param[in] abs_noise
/// \param[in] ridged
/// \return encoded Perlin noise in HeightType [0,1]
inline HeightType encode_height_perlin_noise(
    mi::math::Vector< mi::Sint64, 2 > const & ij,
    mi::Sint64 patch_unit_i,
    mi::Float32 time,
    mi::Uint32 terms,
    mi::math::Vector< mi::Float32, 4 > const & turbulence_weight,
    bool abs_noise,
    bool ridged)
{
    const mi::Float64 patch_unit_f = static_cast< mi::Float64 >(patch_unit_i);
    const mi::Float64 fixed_z = 0.1f;
    mi::math::Vector< mi::Float32, 3 > pos(
        static_cast< mi::Float32 >(static_cast< mi::Float64 >(ij.x % patch_unit_i) / patch_unit_f),
        static_cast< mi::Float32 >(static_cast< mi::Float64 >(ij.y % patch_unit_i) / patch_unit_f),
        static_cast<mi::Float32>(fixed_z));
    const mi::Float32 val = nv::index_common::summed_perlin_noise(pos, time, terms, 
                                                     turbulence_weight, 
                                                     abs_noise, ridged);
    assert(-1.0f <= val);
    assert(val   <= 1.0f);

    return 0.5f * (1.0f + val);
}

//----------------------------------------------------------------------
/// get synthesis method name
///
/// \param[in] synthesis_method  synthesis method identifier.
/// \return synthesis name
inline std::string get_height_synthetic_name(mi::Sint32 synthesis_method)
{
    char const * const p_name[7] = {
        "default",
        "ij",
        "i",
        "j",
        "perlin_noise",
        "0",
        0,
    };
    assert(synthesis_method >= 0);
    assert(synthesis_method <  7);

    return std::string(p_name[synthesis_method]);
}

//----------------------------------------------------------------------
/// get synthesis method ID
///
/// \param[in] synthesis_name   synthesis name
/// \return synthesis type enum. \see Synthetic_heightfield_generator
inline mi::Sint32 get_height_synthetic_type(std::string const & synthesis_name)
{
    if(synthesis_name == "default"){
        return HSM_Default;
    }else if(synthesis_name == "ij"){
        return HSM_IJ;
    }else if(synthesis_name == "i"){
        return HSM_I;
    }else if(synthesis_name == "j"){
        return HSM_J;
    }else if(synthesis_name == "perlin_noise"){
        return HSM_Perlin_noise;
    }else if(synthesis_name == "0"){
        return HSM_0;
    }else{
        return HSM_Count;
    }
}

//----------------------------------------------------------------------
/// compute voxel value based on ijk coordinates.
inline HeightType get_height_synthetic_value(
    mi::math::Vector< mi::Sint64, 2 > const & ij,
    mi::Sint32 synthesis_method)
{
    // no switch here.
    if(synthesis_method == HSM_Default){
        // default: ij
        return encode_height_ij(ij);
    }
    else if(synthesis_method == HSM_IJ){
        return encode_height_ij(ij);
    }
    else if(synthesis_method == HSM_I){
        return encode_height_i(ij);
    }
    else if(synthesis_method == HSM_J){
        return encode_height_j(ij);
    }
    assert(false); // should not happen (HSM_0, HSM_Perlin_noise should not be here.)
    return HeightType(0.0f);
}

//----------------------------------------------------------------------
/// test heightfield height value generated by synthetic heightfield generator
///
/// \param[in] mes message for this test
/// \param[in] heightfield_raw_bbox   raw bounding box (memory image bounding box)
/// \param[in] heightfield_valid_bbox value valid bounding box
/// \param[in] ij_encode_type ijk encode type.
/// \param[in] p_height_data height data
/// \return true when all the height values are expected.
inline bool test_synthetic_height_value(
    std::string const & mes,
    mi::math::Bbox< mi::Sint64, 2 > const & heightfield_raw_bbox,
    mi::math::Bbox< mi::Sint64, 2 > const & heightfield_valid_bbox,
    mi::Sint32 ij_encode_type,
    HeightType * p_height_data
    )
{
    // INFO_LOG << "test_synthetic_height_value: " << mes;
    bool is_ok = true;
    mi::math::Vector< mi::Sint64, 2 > ij(heightfield_valid_bbox.min.x, heightfield_valid_bbox.min.y);
    for(ij[0] = heightfield_valid_bbox.min[0]; ij[0] < heightfield_valid_bbox.max[0]; ++ij[0]){ // i
        for(ij[1] = heightfield_valid_bbox.min[1]; ij[1] < heightfield_valid_bbox.max[1]; ++ij[1]){ // j
            HeightType expected_hval = get_height_synthetic_value(ij, ij_encode_type);
            mi::Sint64 const hor_idx = nv::index_common::get_heightfield_index(ij, heightfield_raw_bbox);
            HeightType const hval = p_height_data[hor_idx];
            // ERROR_LOG << "test at " << hor_idx << ", hval: " << hval;
            if(hval != expected_hval){
                ERROR_LOG << "ij = " << ij << ", idx = " << hor_idx
                          << ", access = " << static_cast< mi::Sint32 >(hval)
                          << " != "        << static_cast< mi::Sint32 >(expected_hval)
                          << " (expected).";
                is_ok = false;
            }
        }
    }
    // INFO_LOG << "test_synthetic_height_value: " << mes << "... check done";
    return is_ok;
}

//----------------------------------------------------------------------

#endif // NVIDIA_INDEX_BIN_COMMON_ENCODE_HEIGHT_H
