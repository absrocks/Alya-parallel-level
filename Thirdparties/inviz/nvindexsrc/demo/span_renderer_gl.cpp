/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "span_renderer_gl.h"

#include <GL/glew.h>
#include <cassert>

#include "common/forwarding_logger.h"

#include "opengl_drawing_utilities.h"

//----------------------------------------------------------------------
Span_renderer_gl::Span_renderer_gl()
  : 
    Span_renderer_IF(),
    m_pbo_enabled(false),
    m_pbo(0)
{
    m_buffer_resolution.x = -1;
    m_buffer_resolution.y = -1;

    m_pbo_enabled = glewGetExtension("GL_ARB_pixel_buffer_object");
    if (m_pbo_enabled)
    {
        glGenBuffersARB(1, &m_pbo);
    }
    else{
        ERROR_LOG << "PBO not supported!";
    }

}

//----------------------------------------------------------------------
Span_renderer_gl::~Span_renderer_gl()
{
    // empty
}

//----------------------------------------------------------------------
void Span_renderer_gl::prepare()
{
    m_screen_space_subdivision.clear();

    // Set blending mode for pre-multiplied alpha (which means that the RGB values were already
    // multiplied by the alpha). 
    glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);

    // Disable depth buffer writes and depth test, as we are doing 2D operations without
    // valid depth information
    glDepthMask(GL_FALSE); 
    glDisable(GL_DEPTH_TEST);
}

//----------------------------------------------------------------------
void Span_renderer_gl::receive_tile(
    mi::Uint8*                                  buffer,
    mi::Uint32                                  buffer_size,
    const mi::math::Bbox_struct<mi::Uint32, 2>& area)
{
    const mi::Uint32 x_range = area.max.x - area.min.x;
    const mi::Uint32 y_range = area.max.y - area.min.y;
    if(x_range==0 || y_range==0)
    {
        DEBUG_LOG << "No image space area to render!";
        return;
    }

    if(!buffer)
    {
        DEBUG_LOG << "No image buffer contents to be renderered - possibly due to internal optimizations";
        return;
    }

    // Acquire the lock
    {
        mi::base::Lock::Block block(&m_receive_tile_lock);
        m_screen_space_subdivision.push_back(area);
    }

    // Always use blending to support correct display of wireframe geometry to visualize the
    // bounding boxes or region of interest.
    glEnable(GL_BLEND);
    
    if(m_pbo_enabled)
    {
        // Bind the PBO
        glBindBufferARB(GL_PIXEL_UNPACK_BUFFER, m_pbo);

        // Create and initialize the PBO data store with the received tile buffer
        glBufferDataARB(GL_PIXEL_UNPACK_BUFFER, buffer_size*sizeof(GLubyte), (void**)buffer, GL_STATIC_READ);

        // Draw the contents of the PBO data store
        glRasterPos2i(area.min.x, area.min.y);
        glDrawPixels(x_range, y_range, GL_RGBA, GL_UNSIGNED_BYTE, 0);

        // Unbind the PBO
        glBindBufferARB(GL_PIXEL_UNPACK_BUFFER_ARB, 0);
    }
    else
    {
        glRasterPos2i(area.min.x, area.min.y);
        glDrawPixels(x_range, y_range, GL_RGBA, GL_UNSIGNED_BYTE, buffer);
    }

    glDisable(GL_BLEND);
}

//----------------------------------------------------------------------
void Span_renderer_gl::receive_tile_blend(
    mi::Uint8*                                      buffer,
    mi::Uint32                                      buffer_size,
    const mi::math::Bbox_struct<mi::Uint32, 2>&     area)
{
    // Blending is used by default
    receive_tile(buffer, buffer_size, area);
}

//----------------------------------------------------------------------
void Span_renderer_gl::finish()
{
    // Restore default settings
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);
}

//----------------------------------------------------------------------
void Span_renderer_gl::set_buffer_resolution(
    const mi::math::Vector_struct<mi::Sint32, 2> & buffer_resolution)
{
    assert(buffer_resolution.x > 0);
    assert(buffer_resolution.y > 0);

    m_buffer_resolution = buffer_resolution;
    gl_initialize_viewport(m_buffer_resolution.x, m_buffer_resolution.y, this->get_background_color());
}

//----------------------------------------------------------------------
mi::math::Vector_struct<mi::Sint32, 2> Span_renderer_gl::get_buffer_resolution() const
{
    return m_buffer_resolution;
}

//----------------------------------------------------------------------
void Span_renderer_gl::get_screen_space_subdivision(
    std::vector<mi::math::Bbox_struct<mi::Uint32, 2> >& screen_space_subdivision)
{
    mi::base::Lock::Block block(&m_receive_tile_lock);
    screen_space_subdivision = m_screen_space_subdivision;
}

//----------------------------------------------------------------------
mi::Size Span_renderer_gl::get_nb_of_screen_space_subdivision()
{
    mi::base::Lock::Block block(&m_receive_tile_lock);
    return m_screen_space_subdivision.size();
}

//----------------------------------------------------------------------
void Span_renderer_gl::copy_pixel(mi::Sint32 width, mi::Sint32 height,
                                  void * data_ptr)
{
    glReadPixels(0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, (GLvoid *)data_ptr);
}

//----------------------------------------------------------------------
