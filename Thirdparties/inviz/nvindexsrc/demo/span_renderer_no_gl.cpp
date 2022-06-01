/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "span_renderer_no_gl.h"

#include <iostream>

#include "common/forwarding_logger.h"

//----------------------------------------------------------------------
Span_renderer_no_gl::Span_renderer_no_gl()
    :
    Span_renderer_IF(),
    m_p_framebuffer(0),
    m_screen_space_subdivision()
{
    m_buffer_resolution.x = -1;
    m_buffer_resolution.y = -1;
}

//----------------------------------------------------------------------
Span_renderer_no_gl::~Span_renderer_no_gl()
{
    DEBUG_LOG << "Span_renderer_no_gl::~Span_renderer_no_gl()";
    if(m_p_framebuffer != 0){
        delete [] m_p_framebuffer;
        m_p_framebuffer = 0;
    }
    m_buffer_resolution.x = -1;
    m_buffer_resolution.y = -1;
    m_screen_space_subdivision.clear();
}

//----------------------------------------------------------------------
void Span_renderer_no_gl::prepare()
{
    m_screen_space_subdivision.resize(0);
}

//----------------------------------------------------------------------
void Span_renderer_no_gl::receive_tile(
    mi::Uint8* span_buffer,
    mi::Uint32 buffer_size,
    const mi::math::Bbox_struct<mi::Uint32, 2>& covered_area)
{
    const mi::Uint32 x_width  = covered_area.max.x - covered_area.min.x;
    const mi::Uint32 y_height = covered_area.max.y - covered_area.min.y;
    if(x_width==0 || y_height==0){
        DEBUG_LOG << "No image space area to render!";
        return;
    }

    {
        mi::base::Lock::Block block(&m_receive_tile_lock);
        m_screen_space_subdivision.push_back(covered_area);
    }

    if(!span_buffer)
    {
        DEBUG_LOG << "No image buffer contents to be renderered - possibly due to internal optimizations";
        return;
    }

    if(m_p_framebuffer == 0){
        DEBUG_LOG << "Span_renderer_no_gl: buffer has no size.";
        return;
    }

    // No blending necessary: Since the incoming image data uses pre-multiplied alpha (which means
    // that the RGB values were already multiplied by the alpha), a simple copy is sufficient.
    const mi::Uint8* p_span_buffer_iter = span_buffer;
    const mi::Uint32 copy_byte = x_width * 4;

    for (mi::Uint32 y = covered_area.min.y; y < covered_area.max.y; ++y)
    {
        // Copy scanline with memcpy
        assert(covered_area.min.x < covered_area.max.x);
        mi::Uint8* p_dest = m_p_framebuffer + get_index(covered_area.min.x, y, 0);
        std::memcpy(p_dest, p_span_buffer_iter, copy_byte);
        p_span_buffer_iter += copy_byte;
    }
}

//----------------------------------------------------------------------
void Span_renderer_no_gl::receive_tile_blend(
    mi::Uint8* span_buffer,
    mi::Uint32 buffer_size,
    const mi::math::Bbox_struct<mi::Uint32, 2>& covered_area)
{
    const mi::Uint32 x_width  = covered_area.max.x - covered_area.min.x;
    const mi::Uint32 y_height = covered_area.max.y - covered_area.min.y;
    if(x_width==0 || y_height==0)
    {
        DEBUG_LOG << "No image space area to render!";
        return;
    }

    {
        mi::base::Lock::Block block(&m_receive_tile_lock);
        m_screen_space_subdivision.push_back(covered_area);
    }

    if(!span_buffer)
    {
        DEBUG_LOG << "No image buffer contents to be renderered - possibly due to internal optimizations";
        return;
    }

    if(m_p_framebuffer == 0){
        DEBUG_LOG << "Span_renderer_no_gl: buffer has no size.";
        return;
    }

    // Blend received data (source) over existing image (destination)
    const mi::Uint8* p_src = span_buffer;
    for (mi::Uint32 y = covered_area.min.y; y < covered_area.max.y; ++y)
    {
        mi::Uint8* p_dst = m_p_framebuffer + get_index(covered_area.min.x, y, 0);
        for (mi::Uint32 x = covered_area.min.x; x < covered_area.max.x; ++x)
        {
            const mi::Uint8 src_alpha = p_src[3];
            if (src_alpha == 0)
            {
                // Source is fully transparent, destination unchanged
                p_src += 4;
                p_dst += 4;                
            }
            else if (src_alpha == 255)
            {
                // Source is fully opaque, just overwrite
                *p_dst++ = *p_src++;
                *p_dst++ = *p_src++;
                *p_dst++ = *p_src++;
                *p_dst++ = *p_src++;
            }
            else
            {
                // Blend just as glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA) would do, to blend the
                // image with pre-multiplied alpha.
                mi::Float32 r_dst = p_dst[0] / 255.f;
                mi::Float32 g_dst = p_dst[1] / 255.f;
                mi::Float32 b_dst = p_dst[2] / 255.f;
                mi::Float32 a_dst = p_dst[3] / 255.f;

                mi::Float32 r_src = p_src[0] / 255.f;
                mi::Float32 g_src = p_src[1] / 255.f;
                mi::Float32 b_src = p_src[2] / 255.f;
                mi::Float32 a_src = p_src[3] / 255.f;
                p_src += 4;

                const mi::Float32 a_one_minus_src = 1.f - a_src;
                *p_dst++ = static_cast<mi::Uint8>(255 * mi::math::saturate(r_src + r_dst * a_one_minus_src));
                *p_dst++ = static_cast<mi::Uint8>(255 * mi::math::saturate(g_src + g_dst * a_one_minus_src));
                *p_dst++ = static_cast<mi::Uint8>(255 * mi::math::saturate(b_src + b_dst * a_one_minus_src));
                *p_dst++ = static_cast<mi::Uint8>(255 * mi::math::saturate(a_src + a_dst * a_one_minus_src));
            }
        }
    }
}

//----------------------------------------------------------------------
void Span_renderer_no_gl::finish()
{
    // empty
}

//----------------------------------------------------------------------
void Span_renderer_no_gl::set_buffer_resolution(
    const mi::math::Vector_struct<mi::Sint32, 2> & screen_resolution_size)
{
    assert(screen_resolution_size.x > 0);
    assert(screen_resolution_size.y > 0);

    this->resize_buffer(screen_resolution_size.x, screen_resolution_size.y);

    // initialize the framebuffer
    const mi::math::Color_struct bg = this->get_background_color();
    if(bg.r == bg.g && bg.g == bg.b)
    {
        mi::Uint8 value = static_cast<mi::Uint8>(255 * bg.r);
        memset(m_p_framebuffer, value, this->get_framebuffer_size());
    }
    else
    {
        mi::Uint8* p_dst = &m_p_framebuffer[0];
        const mi::Uint32 canvas_size = screen_resolution_size.x * screen_resolution_size.y;
        for(mi::Uint32 i=0; i<canvas_size; ++i)
        {
            *p_dst++ = static_cast<mi::Uint8>(255 * bg.r);
            *p_dst++ = static_cast<mi::Uint8>(255 * bg.g);
            *p_dst++ = static_cast<mi::Uint8>(255 * bg.b);
            *p_dst++ = static_cast<mi::Uint8>(0);
        }
    }
}

//----------------------------------------------------------------------
mi::math::Vector_struct<mi::Sint32, 2> Span_renderer_no_gl::get_buffer_resolution() const
{
    return m_buffer_resolution;
}

//----------------------------------------------------------------------
void Span_renderer_no_gl::get_screen_space_subdivision(
    std::vector<mi::math::Bbox_struct<mi::Uint32, 2> >& screen_space_subdivision)
{
    mi::base::Lock::Block block(&m_receive_tile_lock);
    screen_space_subdivision = m_screen_space_subdivision;
}

//----------------------------------------------------------------------
mi::Size Span_renderer_no_gl::get_nb_of_screen_space_subdivision()
{
    mi::base::Lock::Block block(&m_receive_tile_lock);
    return m_screen_space_subdivision.size();
}

//----------------------------------------------------------------------
void Span_renderer_no_gl::copy_pixel(mi::Sint32 width, mi::Sint32 height,
                                     void * data_ptr)
{
    memcpy(data_ptr, m_p_framebuffer, width * height * 4);
}

//----------------------------------------------------------------------
void Span_renderer_no_gl::resize_buffer(mi::Sint32 width, mi::Sint32 height)
{
    assert(width  >= 0);
    assert(height >= 0);

    if((width == m_buffer_resolution.x) && (height == m_buffer_resolution.y)){
        // no need to resize
        return;
    }

    // first need to set the width and height for get_framebuffer_size()
    m_buffer_resolution.x = width;
    m_buffer_resolution.y = height;

    delete[] m_p_framebuffer;
    m_p_framebuffer = new mi::Uint8[this->get_framebuffer_size()];
}

//----------------------------------------------------------------------
