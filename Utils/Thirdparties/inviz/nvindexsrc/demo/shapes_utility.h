/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief volume related utilities

#ifndef NVIDIA_INDEX_UTILITY_SHAPES_H
#define NVIDIA_INDEX_UTILITY_SHAPES_H

#include <mi/dice.h>

#include <string>


//----------------------------------------------------------------------

/// ...
/// ...
///
/// \param[in] dice_transaction         The transaction used for the operation.
///
/// \param[in] session_tag              Session operating IndeX.
///
/// \param[in] group_name               The name of the scenegraph group.
///
/// \return tag                         The tag identifies the newly created line set
///
mi::neuraylib::Tag create_line_set(
    const mi::math::Bbox<mi::Float32, 3>&                       bbox,
    const mi::math::Color&                                      color,
    const mi::base::Handle<mi::neuraylib::IDice_transaction>&   dice_transaction,
    const mi::neuraylib::Tag&                                   session_tag,
    const std::string&                                          group_name);

//----------------------------------------------------------------------


#endif // NVIDIA_INDEX_UTILITY_SHAPES_H
