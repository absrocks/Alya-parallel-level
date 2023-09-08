/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "colormap_ppm.h"

#include "forwarding_logger.h"

namespace nv {
namespace index_common {
//----------------------------------------------------------------------
bool colormap_to_ppm(
    const std::vector<mi::math::Color>& colormap_entries,
    PPM_image& ppm_image)
{
    if (colormap_entries.empty())
    {
        ERROR_LOG << "colormap_to_ppm: empty colormap.";
        return false;
    }

    const mi::Sint32 height = 256;    // alpha value range [0,255]
    mi::Sint32 width = static_cast<mi::Sint32>(colormap_entries.size());

    ppm_image.resize_buffer(width, height);
    mi::math::Color_struct fillcol = { 0.0f, 0.0f, 0.0f, 1.0f, };
    ppm_image.clear_buffer(fillcol);

    for (mi::Sint32 x = 0; x < width; ++x)
    {
        const mi::math::Color col = colormap_entries[x];
        const mi::Sint32 alpha_height = static_cast<mi::Sint32>(255.0f * col[3]); // expand [0,1] -> [0,255]
        // ERROR_LOG << "alpha_height: " << alpha_height << ", " << col[3];
        assert(alpha_height >= 0);
        assert(alpha_height <  height);

        for (mi::Sint32 y = 0; y < alpha_height; ++y) 
        {
            const mi::Sint32 pixel_y = height - y - 1;
            assert(pixel_y >= 0);
            ppm_image.put_color(x, pixel_y, col);
        }
    }

    return true;
}

//----------------------------------------------------------------------
} // namespace index_common
} // namespace nv
