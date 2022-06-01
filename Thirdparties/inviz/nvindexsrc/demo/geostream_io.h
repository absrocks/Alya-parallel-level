/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief geostream viewer data IO (input/output)

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_GEOSTREAM_IO_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_GEOSTREAM_IO_H

#include <mi/dice.h>

#include "common/string_dict.h"

#include <string>
//----------------------------------------------------------------------
/// get volume "filter" generation configuration
/// \param[in] dice_transaction dice transaction
/// \param[in] session_tag current session
/// \param[in] volume_dict  copy generation parameters in a dict
/// \return import configuration tag
mi::neuraylib::Tag get_filter_volume_config_tag(
    const mi::base::Handle<mi::neuraylib::IDice_transaction>& dice_transaction,
    const mi::neuraylib::Tag& session_tag,
    nv::index_common::String_dict & volume_dict);

//----------------------------------------------------------------------
/// Add the scene graph specified in the project files
///
/// \param[in] dice_transaction dice db transaction
/// \param[in] session_tag session tag
/// \param[in] dict  scene graph dict data
/// \return number of add scene graph nodes.
mi::Sint32 add_scene_graph(
    const mi::base::Handle<mi::neuraylib::IDice_transaction>& dice_transaction,
    mi::neuraylib::Tag                                        session_tag,
    const nv::index_common::String_dict&                      dict);

//----------------------------------------------------------------------

#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_GEOSTREAM_IO_H
