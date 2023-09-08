/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief color space conversion

#ifndef NVIDIA_INDEX_BIN_COMMON_COLOR_SPACE_H
#define NVIDIA_INDEX_BIN_COMMON_COLOR_SPACE_H

#include <mi/math/color.h>

namespace nv {
namespace index_common {
//----------------------------------------------------------------------
/// Get HSV color from RGB. Alpha will be copied.
///
/// \param[in] rgb_color  RGB color
/// \return converted HSV color (hue, saturation, value, alpha)
mi::math::Color rgb_to_hsv(
    const mi::math::Color& rgb_color);

//----------------------------------------------------------------------
/// Get HSV color from RGB. Alpha will be copied.
///
/// \param[in] hsv_color  HSV color (hue, saturation, value, alpha)
/// \return converted RGB color 
mi::math::Color hsv_to_rgb(
    const mi::math::Color& hsv_color);

//----------------------------------------------------------------------

} // namespace index_common
} // namespace nv
#endif // NVIDIA_INDEX_BIN_COMMON_COLOR_SPACE_H
