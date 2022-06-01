/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "color_space.h"

namespace nv {
namespace index_common {

//----------------------------------------------------------------------
mi::math::Color rgb_to_hsv(
    const mi::math::Color& rgb_color)
{
    const mi::Float32 r = rgb_color.r;
    const mi::Float32 g = rgb_color.g;
    const mi::Float32 b = rgb_color.b;

    const mi::Float32 max = std::max(std::max(r, g), b);
    const mi::Float32 min = std::min(std::min(r, g), b);
    
    mi::Float32 h = max - min;

    if (h > 0.0f) 
    {
        if (max == r)
        {
            h = (g - b) / h;
            if (h < 0.0f)
            {
                h += 6.0f;
            }
        }
        else if (max == g)
        {
            h = 2.0f + (b - r) / h;
        }
        else
        {
            h = 4.0f + (r - g) / h;
        }
    }
    h /= 6.0f;

    mi::Float32 s = (max - min);
    if (max != 0.0f)
    {
        s /= max;
    }
    mi::Float32 v = max;

    return mi::math::Color(h, s, v, rgb_color.a);
}

//----------------------------------------------------------------------
mi::math::Color hsv_to_rgb(
    const mi::math::Color& hsv_color)
{
    const mi::Float32 h = hsv_color[0];
    const mi::Float32 s = hsv_color[1];
    const mi::Float32 v = hsv_color[2];
    const mi::Float32 a = hsv_color[3];

    if (s == 0.0f) 
    {
        return mi::math::Color(v, v, v, a);
    }

    const mi::Float32 vh = h * 6.0f;
    const mi::Float32 vi = std::floor(vh);

    const mi::Float32 v1 = v * (1.0f - s);
    const mi::Float32 v2 = v * (1.0f - s * (vh - vi));
    const mi::Float32 v3 = v * (1.0f - s * (1.0f - (vh - vi)));

    switch (static_cast<mi::Sint32>(vi))
    {
    case 0:
        return mi::math::Color( v, v3, v1, a);
    case 1:
        return mi::math::Color(v2,  v, v1, a);
    case 2:
        return mi::math::Color(v1,  v, v3, a);
    case 3:
        return mi::math::Color(v1, v2,  v, a);
    case 4:
        return mi::math::Color(v3, v1,  v, a);
    case 5:
        return mi::math::Color( v, v1, v2, a);
    default:
        return mi::math::Color( v, v3, v1, a);
    }

    return mi::math::Color(0.0f, 0.0f, 0.0f, a); // not really needed
}

//----------------------------------------------------------------------
} // namespace index_common
} // namespace nv
