/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief colormap IO routine

#ifndef NVIDIA_INDEX_BIN_COMMON_COLORMAP_IO_H
#define NVIDIA_INDEX_BIN_COMMON_COLORMAP_IO_H

#include <vector>
#include <string>
#include <cassert>

#include <mi/math/color.h>

namespace nv {
namespace index_common {
//----------------------------------------------------------------------
/// Load colormap from a file
///
/// \param[in]  colormap_fname colormap filename
/// \param[out] colormap_entries (output) colormap: array of colors
/// \param[out] error_mes        (output) error message when return false.
/// \return true when load succeeded.
bool load_colormap(
    const std::string &                     colormap_fname,
    std::vector< mi::math::Color_struct > & colormap_entries,
    std::string & error_mes);

//----------------------------------------------------------------------
/// return extension from the fname.
/// \param[in] fname file name
/// \return extension part of the fname including the '.'.
std::string get_filename_extension(const std::string & fname);

//----------------------------------------------------------------------
/// check this is 16bit colormap file
/// \param[in] colormap_fname colormap filename
/// \return true if the colormap filename has 16 bit type extension (".cmap")
bool is_16bit_colormap_file(const std::string & colormap_fname);

//----------------------------------------------------------------------
/// load 16 bit colormap file
/// \param[in] colormap_fname    colormap file name for read
/// \param[out] colormap_entries (output) colormap entries
/// \param[out] error_mes        (output) error message when return false
/// \return true when load succeeded
bool load_16bit_colormap(
    const std::string &                     colormap_fname,
    std::vector< mi::math::Color_struct > & colormap_entries,
    std::string & error_mes);

//----------------------------------------------------------------------
/// write 16 bit colormap file
/// \param[in] colormap_fname    colormap file name for write
/// \param[in] colormap_entries  colormap entries
/// \param[out] error_mes        (output) error message when return false
/// \return true when load succeeded
bool write_16bit_colormap(
    std::string const &                           colormap_fname,
    std::vector< mi::math::Color_struct > const & colormap_entries,
    std::string & error_mes);

//----------------------------------------------------------------------

}} // namespace nv::index_common
#endif // NVIDIA_INDEX_BIN_COMMON_COLORMAP_IO_H
