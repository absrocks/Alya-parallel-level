/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief interpolating 3D points

#ifndef NVIDIA_INDEX_GEOSPATIAL_BIN_COMMON_SPLINE_H
#define NVIDIA_INDEX_GEOSPATIAL_BIN_COMMON_SPLINE_H

#include <mi/math/vector.h>
#include <vector>

namespace nv
{
namespace index_common
{

class Cubic_hermite_spline
{
public:
    Cubic_hermite_spline()
        :
        m_tension(0.0f),
        m_bias(0.0f)
    {
        // empty
    }
    
    void build_curve(
        const mi::math::Vector<mi::Float32, 3> *points,
        mi::Size    num_points,
        mi::Float32 bias,
        mi::Float32 tension);
    
    mi::math::Vector<mi::Float32, 3> get_value(
        mi::Size idx,
        mi::Float32 t);
    
private:
    std::vector<mi::math::Vector<mi::Float32, 3> >   m_points;
    std::vector<mi::math::Vector<mi::Float32, 3> >   m_mk;
    mi::Float32                                      m_tension;
    mi::Float32                                      m_bias;
};

} // namespace index_common
} // namespace nv
#endif // NVIDIA_INDEX_GEOSPATIAL_BIN_COMMON_SPLINE_H
