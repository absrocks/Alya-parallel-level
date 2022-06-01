/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief colormap to ppm convert

#ifndef NVIDIA_INDEX_BIN_COMMON_COLORMAP_PPM_H
#define NVIDIA_INDEX_BIN_COMMON_COLORMAP_PPM_H

#include <vector>

#include <mi/math/color.h>

#include "ppm_io.h"

namespace nv {
namespace index_common {
//----------------------------------------------------------------------
/// Colormap to ppm format convert
///
/// Note: only handle discretized colormap (color array) 
///
/// \param[in]  colormap_entries colormap: array of colors
/// \param[out] ppm_image (output) ppm image object
/// \return true when succeeded
bool colormap_to_ppm(
    const std::vector<mi::math::Color>& colormap_entries,
    PPM_image& ppm_image);

//----------------------------------------------------------------------
} // namespace index_common
} // namespace nv
#endif // NVIDIA_INDEX_BIN_COMMON_COLORMAP_PPM_H
