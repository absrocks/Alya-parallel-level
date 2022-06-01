/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief heightfield seed manipuration

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_HEIGHTFIELD_SEED_MANIPURATION_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_HEIGHTFIELD_SEED_MANIPURATION_H

#include <mi/dice.h>
#include <string>

//----------------------------------------------------------------------
/// Loading a set of seed lines given in the heightfield's local ij space.
///
/// \param[in] scope the scope
/// \param[in] session_tag              session tag
/// \param[in] heightfield_scene_element_id heightfield scene element id
/// \param[in] seed_lines_file_name seed lines file name
void import_seed_lines(
    mi::base::Handle<mi::neuraylib::IScope>&    scope,
    const mi::neuraylib::Tag&                   session_tag,
    mi::Uint32                                  heightfield_scene_element_id,
    const std::string&                          seed_lines_file_name);

//----------------------------------------------------------------------
/// Remove an arbitrary seed line from the heightfield dataset (for testing purposes).
///
/// \param[in] scope the scope
/// \param[in] session_tag              session tag
/// \param[in] heightfield_scene_element_id heightfield scene element id
void remove_arbitrary_seed_line(
    mi::base::Handle<mi::neuraylib::IScope>&    scope,
    const mi::neuraylib::Tag&                   session_tag,
    mi::Uint32                                  heightfield_scene_element_id);

//----------------------------------------------------------------------
/// Remove an arbitrary seed line from the heightfield dataset (for testing purposes).
///
/// \param[in] scope the scope
/// \param[in] session_tag              session tag
/// \param[in] heightfield_scene_element_id heightfield scene element id
void remove_arbitrary_seed_point(
    mi::base::Handle<mi::neuraylib::IScope>&    scope,
    const mi::neuraylib::Tag&                   session_tag,
    mi::Uint32                                  heightfield_scene_element_id);

#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_HEIGHTFIELD_SEED_MANIPURATION_H
