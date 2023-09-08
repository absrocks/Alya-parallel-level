/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "opengl_appdata.h"
#include "opengl_application_buffer.h"
#include "opengl_drawing_utilities.h"

#include "common/forwarding_logger.h"

#ifdef USE_OPENGL
#include "opengl_offscreen_context.h"
#endif

//----------------------------------------------------------------------
OpenGL_AppData::~OpenGL_AppData()
{
    gl_delete_offscreen_context(m_p_context);
    m_p_context = 0;

    if(m_p_opengl_application_buffer != 0)
    {
        delete m_p_opengl_application_buffer;
        m_p_opengl_application_buffer = 0;
    }
}

//----------------------------------------------------------------------
void OpenGL_AppData::set_offscreen_context(Offscreen_context * p_context)
{
    m_p_context = p_context;
}

//----------------------------------------------------------------------
Opengl_application_buffer * OpenGL_AppData::get_opengl_application_buffer_ptr()
{
    if (m_p_opengl_application_buffer == 0)
    {
        m_p_opengl_application_buffer = new Opengl_application_buffer();
    }
    return m_p_opengl_application_buffer;
}

//----------------------------------------------------------------------
OpenGL_AppData::OpenGL_AppData()
    :
    m_is_show_color_table(false),
    m_color_table_position(mi::math::Vector< mi::Sint32, 2 >(60, 20)),
    m_colormap_scale(1.0f),
    m_axis_display_list(0),
    m_quad_display_list(0),
    m_cube_display_list(0),
    m_p_context(0),
    m_p_opengl_application_buffer(0)
{
    // empty
}

//----------------------------------------------------------------------
#ifdef USE_OPENGL
//----------------------------------------------------------------------
void OpenGL_AppData::resize_offscreen_context(mi::Sint32 width, mi::Sint32 height)
{
    if (m_p_context == 0)
    {
        ERROR_LOG << "OpenGL_AppData::resize_offscreen_context: no offscreen context ready.";
        return;
    }
    m_p_context->resize(width, height);
}
//----------------------------------------------------------------------
#else
//----------------------------------------------------------------------
void OpenGL_AppData::resize_offscreen_context(mi::Sint32 width, mi::Sint32 height)
{
    // empty
}

//----------------------------------------------------------------------
#endif

