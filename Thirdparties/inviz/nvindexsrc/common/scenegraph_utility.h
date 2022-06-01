/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief volume related utilities

#ifndef NVIDIA_INDEX_UTILITY_SCENEGRAPH_H
#define NVIDIA_INDEX_UTILITY_SCENEGRAPH_H

#include <mi/dice.h>

#include <string>

namespace nv {
namespace index_common {

//----------------------------------------------------------------------
/// ...
/// ...
///
/// \param[in] session_tag              Session operating IndeX.
///
/// \param[in] group_name               The name of the scenegraph group stored in the database.
///
/// \param[in] dice_transaction         The transaction used for the operation.
///
/// \return tag of the newly created volume
///
mi::neuraylib::Tag get_scene_group(
    const mi::neuraylib::Tag&         session_tag,
    const std::string&                group_name,
    mi::neuraylib::IDice_transaction* dice_transaction);


/// ...
/// ...
///
/// \param[in] session_tag              Session operating IndeX.
///
/// \param[in] group_name               The name of the scenegraph group.
///
/// \param[in] dice_transaction         The transaction used for the operation.
///
/// \return tag of the newly created volume
///
mi::neuraylib::Tag create_scene_group(
    const mi::neuraylib::Tag&         session_tag,
    const std::string&                group_name,
    mi::neuraylib::IDice_transaction* dice_transaction);

//----------------------------------------------------------------------

}} // namespace nv::index_common

#endif // NVIDIA_INDEX_UTILITY_SCENEGRAPH_H
