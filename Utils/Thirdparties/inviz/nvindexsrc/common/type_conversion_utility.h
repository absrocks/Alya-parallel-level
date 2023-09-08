/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief type conversion utiliy

#ifndef COMMON_TYPE_CONVERSION_UTILITY_H
#define COMMON_TYPE_CONVERSION_UTILITY_H

#include <mi/dice.h>

#include <cassert>

namespace nv {
namespace index_common {
/// type of voxel (fixed, depends on library)
typedef mi::Uint8 VoxelType;

//----------------------------------------------------------------------
/// get mi::math::Color_struct<mi::Float32, 4> from mi::math::Color<mi::Float32, 4>
/// \param[in] col color (r, g, b, a)
/// \return    POD struct type color 3
inline mi::math::Color_struct get_color_struct(
    const mi::math::Color & col)
{
    mi::math::Color_struct col_st;

    col_st.r = col.r;
    col_st.g = col.g;
    col_st.b = col.b;
    col_st.a = col.a;

    return col_st;
}

//----------------------------------------------------------------------
/// get mi::math::Color_struct<mi::Float32, 4> from float rgba. Since
/// Color_struct is pure POD type and should not have a
/// constructor. This is a constructor function.
///
/// \param[in] r col color r
/// \param[in] g col color g
/// \param[in] b col color b
/// \param[in] a col color a
/// \return    POD struct type color
inline mi::math::Color_struct get_color_struct(
    mi::Float32 r, mi::Float32 g, mi::Float32 b, mi::Float32 a)
{
    const mi::math::Color_struct col_st = { r, g, b, a };
    return col_st;
}

//======================================================================
/// Bounding box type conversion
/// FromType, DIM bounding box to ToType, DIM bounding box conversion
/// \param[in] bbox_fromtype  <FromType, DIM> bounding box
/// \return    <ToType, DIM> bounding box
template < typename ToType, typename FromType,  mi::Size DIM>
inline mi::math::Bbox< ToType, DIM > convert_bbox_type(
    const mi::math::Bbox< FromType, DIM > & bbox_fromtype)
{
    assert(!bbox_fromtype.empty());
    mi::math::Bbox< ToType, DIM > bbox_totype;
    
    for(mi::Size i = 0; i < DIM; ++i){
        bbox_totype.min[i] = static_cast< ToType >(bbox_fromtype.min[i]);
        bbox_totype.max[i] = static_cast< ToType >(bbox_fromtype.max[i]);
    }
    
    return bbox_totype;
}

//======================================================================
// conversion: bbox_sint64_3 to bbox_*

//----------------------------------------------------------------------
/// Sint64, 3 bounding box to Uint32, 3 bounding box conversion
/// \param[in] bbox_sint64_3  Sint64, 3 bouding box
/// \return    Uint32, 3 bouding box
inline mi::math::Bbox< mi::Uint32, 3 > convert_bbox_sint64_3_to_bbox_uint32_3(
    const mi::math::Bbox< mi::Sint64, 3 > & bbox_sint64_3)
{
    assert(!bbox_sint64_3.empty());

    mi::math::Bbox< mi::Uint32, 3 > bbox_uint32_3;

    for(mi::Sint32 i = 0; i < 3; ++i){
        assert(bbox_sint64_3.min[i] >= 0);
        assert(bbox_sint64_3.min[i] <=
               static_cast< mi::Sint64 >(mi::base::numeric_traits<mi::Uint32>::max()));
        assert(bbox_sint64_3.max[i] >= 0);
        assert(bbox_sint64_3.max[i] <=
               static_cast< mi::Sint64 >(mi::base::numeric_traits<mi::Uint32>::max()));

        bbox_uint32_3.min[i] = static_cast< mi::Uint32 >(bbox_sint64_3.min[i]);
        bbox_uint32_3.max[i] = static_cast< mi::Uint32 >(bbox_sint64_3.max[i]);
    }

    return bbox_uint32_3;
}


//----------------------------------------------------------------------
}} // namespace nv::index_common
#endif // COMMON_TYPE_CONVERSION_UTILITY_H
