/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief interpolating 3D points

#include "spline.h"

#include "forwarding_logger.h"

#include <cassert>

namespace nv
{
namespace index_common
{
//----------------------------------------------------------------------
void Cubic_hermite_spline::build_curve(
    const mi::math::Vector<mi::Float32, 3> *points,
    mi::Size num_points,
    mi::Float32 bias,
    mi::Float32 tension)
{
    assert(num_points >= 2);

    m_points.resize(num_points);
    
    // INFO_LOG << "SPLINE: num_points: " << num_points;
    for(mi::Size i = 0; i < num_points; ++i)
    {
        m_points[i] = points[i];
        // INFO_LOG << "SPLINE: p[" << i << "] = " << points[i];
        // INFO_LOG << "SPLINE: P[" << i << "] = " << points[i];
    }
    
    m_bias    = bias;
    m_tension = tension;
    mi::Float32 f = (1.f + bias)*(1.f - tension);
    
    m_mk.resize(num_points);
    
    // INFO_LOG << "SPLINE: bias: " << bias;
    // INFO_LOG << "SPLINE: tension: " << tension;

    m_mk[0] = (m_points[1] - m_points[0])*f;
    m_mk[num_points-1] = (m_points[num_points-1] - m_points[num_points-2])*f;
    
    for (mi::Size i = 1; i < num_points-1; ++i)
    {
        m_mk[i] = (m_points[i+1] - m_points[i-1])*(0.5*f);
    }
}

//----------------------------------------------------------------------
mi::math::Vector<mi::Float32, 3> Cubic_hermite_spline::get_value(mi::Size idx, mi::Float32 t)
{
    mi::math::Vector<mi::Float32, 3> v;

    if(idx > m_points.size() - 1)
    {
        WARN_LOG << "SPLINE: Out of range...";
    
        v.x = v.y = v.z = 0.f;
    }
    else
    {
        float tt = t*t;
        float ttt = tt*t;

        v = m_points[idx]*(2.f*ttt - 3.f*tt + 1.f) + m_mk[idx]*(ttt - 2.f*tt + t)
            + m_mk[idx+1]*(ttt - tt) + m_points[idx+1]*(-2.f*ttt + 3.f*tt);
            
        // INFO_LOG << "SPLINE: idx: " << idx << ", " << idx+1;
        // INFO_LOG << "SPLINE: points: " << m_points[idx] << ", " << m_points[idx+1];
        // INFO_LOG << "SPLINE: mk: " << m_mk[idx] << ", " << m_mk[idx+1];
        // INFO_LOG << "SPLINE: p(" << t << ")= " << v;
    }

    return v;
}
//----------------------------------------------------------------------

} // namespace index_common
} // namespace nv
