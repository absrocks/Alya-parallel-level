/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief opengl application draw for no opengl (stubs)

#include "opengl_application_draw.h"

//----------------------------------------------------------------------
void gl_application_initialize(Opengl_application_buffer*,
                               const mi::math::Vector<mi::Sint32,2>& main_window_size)
{
    // empty
}

//----------------------------------------------------------------------
void gl_application_draw_object(Opengl_application_buffer*)
{
    // empty
}

//----------------------------------------------------------------------
void gl_application_example_well_creation(
    const mi::math::Bbox< mi::Float32, 3 >& roi_xyz,
    const nv::index_common::String_dict& well_opt)
{
    // empty
}

//----------------------------------------------------------------------
void gl_application_shutdown()
{
    // empty
}

//----------------------------------------------------------------------

