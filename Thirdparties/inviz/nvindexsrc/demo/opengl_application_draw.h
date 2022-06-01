/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief opengl application draw

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_OPENGL_APPLICATION_DRAW_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_OPENGL_APPLICATION_DRAW_H

#include <mi/dice.h>
#include <mi/math/bbox.h>


// forward declaration
namespace nv {
namespace index_common {
class String_dict;
}} // namespace
class Opengl_application_buffer;

//----------------------------------------------------------------------
/// application OpenGL initialization.
///
/// Application side OpenGL related set up. Mainly, static primitives
/// setup.
///
/// \param[in] p_app_buf     opengl application buffer
/// \param[in] window_resolution  window (buffer) resolution
void gl_application_initialize(Opengl_application_buffer * p_app_buf,
                               const mi::math::Vector<mi::Sint32,2>& window_resolution);

//----------------------------------------------------------------------
/// application OpenGL draw call.
///
/// application can draw any opaque object that override the z-buffer
/// here.  When this is called, model view and camera matrices have
/// already set up. It is your responsibility if the matrices changed.
///
/// \param[in,out] p_app_buf OpenGL application buffer. Should be
/// initialized, the RGBA and Z buffer is written in this function.
void gl_application_draw_object(
    Opengl_application_buffer* p_app_buf);


//----------------------------------------------------------------------
/// application OpenGL drawing example: create wells.
///
/// \param[in] roi_xyz  region of interest in xyz coordinates
/// \param[in] well_opt well creation option
void gl_application_example_well_creation(
    const mi::math::Bbox< mi::Float32, 3 >& roi_xyz,
    const nv::index_common::String_dict&    well_opt);

//----------------------------------------------------------------------
/// application OpenGL shutdown.
///
/// Application side OpenGL related shutdown. Mainly, static primitives
/// deletion.
///
void gl_application_shutdown();

#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_OPENGL_APPLICATION_DRAW_H
