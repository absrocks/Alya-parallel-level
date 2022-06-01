/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief colormap utilities

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_COLORMAP_UTIL_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_COLORMAP_UTIL_H

#include <nv/index/isession.h>
#include <nv/index/iscene.h>

#include <map>
#include <cassert>
#include <vector>

#include "common/forwarding_logger.h"
#include "common/string_dict.h"

//----------------------------------------------------------------------

/// get colormap tag from an index
///
/// \param[in] colormap_idx colormap index
/// \return    colormap tag corresponds colormap_idx
mi::neuraylib::Tag get_colormap_tag(mi::Uint32 colormap_idx);

//----------------------------------------------------------------------
/// set colormap tag from an index
///
/// \param[in] colormap_idx colormap index
/// \param[in] tag          colormap tag corresponds colormap_idx
void set_colormap_tag(mi::Uint32                colormap_idx, 
                      const mi::neuraylib::Tag& tag);

//----------------------------------------------------------------------
/// get number of colormaps
///
/// \return    number of colormaps (not each number of colormap's entries)
mi::Uint32 get_number_of_colormap();

//----------------------------------------------------------------------
/// load all colormap
///
/// \param[in] scene             scene object
/// \param[in] app_config        application configuration (colormap project)
/// \param[in] dice_transaction  dice transaction to access dice database
/// \return true when succeeded
bool load_all_colormap(
    const nv::index::IScene*             scene,
    const nv::index_common::String_dict& app_config,
    mi::neuraylib::IDice_transaction*    dice_transaction);

//----------------------------------------------------------------------
/// sanity check all the colormaps are valid.
///
/// Note: Only works when not defined USE_SYSTEM_COLORMAP
///
/// \param[in] dice_transaction  dice transaction to access dice database
bool check_all_colormap_valid(
    mi::neuraylib::IDice_transaction* dice_transaction);

//----------------------------------------------------------------------
/// operate user defined colormap operation
///
/// \param[in] cmap_tag         colormap tag to operate
/// \param[in] dice_transaction dice transaction to access dice database
/// \return true when succeeded.
bool exec_userdef_colormap_operation(
    const mi::neuraylib::Tag          cmap_tag,
    mi::neuraylib::IDice_transaction* dice_transaction);

//----------------------------------------------------------------------

#endif // #ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_COLORMAP_UTIL_H
