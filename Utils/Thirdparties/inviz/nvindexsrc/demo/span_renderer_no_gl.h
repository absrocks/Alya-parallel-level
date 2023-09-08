/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief span renderer with no-OpenGL functions

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_SPAN_RENDERER_NO_GL_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_SPAN_RENDERER_NO_GL_H

#include <nv/index/iindex.h>
#include <mi/base/interface_implement.h>
#include <mi/base/lock.h>

#include <vector>
#include <cassert>

#include "span_renderer_if.h"


/// This class implement the span buffer rendering which enables
/// user-defined rendering of span data and shall be passed to the
/// rendering call ('iindex_rendering->render').
///
/// This particular blends all spans into an user framebuffer without
/// any OpenGL functionality.
class Span_renderer_no_gl : public Span_renderer_IF
{
public:
    /// constructor
    Span_renderer_no_gl();
    /// destructor
    virtual ~Span_renderer_no_gl();

    /// prepare the receive tile
    virtual void prepare();

    /// receive a tile
    ///
    /// One frame is constructed as the following pseudo code by the
    /// IndeX library
    ///
    /// \code
    /// span_buffer->prepare();
    /// while(sending_tiles){
    ///     span_buffer->receive_tile(...);
    /// }
    /// span_buffer->finish();
    /// \endcode
    ///
    /// \param[in] span_buffer  span buffer
    /// \param[in] buffer_size  span buffer sieze
    /// \param[in] covered_area area covered by the span buffer
    virtual void receive_tile(
        mi::Uint8* span_buffer,
        mi::Uint32 buffer_size,
        const mi::math::Bbox_struct<mi::Uint32, 2>& covered_area);

    virtual void receive_tile_blend(
        mi::Uint8* span_buffer,
        mi::Uint32 buffer_size,
        const mi::math::Bbox_struct<mi::Uint32, 2>& covered_area);

    virtual void finish();

    virtual bool is_multi_thread_capable() const
    {
        // Both serial and parallel rendering is supported, so let the user decide
        return Span_renderer_IF::get_use_parallel_rendering();
    }

    virtual void set_buffer_resolution(
        const mi::math::Vector_struct<mi::Sint32, 2> & buffer_resolution);

    virtual mi::math::Vector_struct<mi::Sint32, 2> get_buffer_resolution() const;

public:
    /// get screen space subdivision information for visualization
    ///
    /// \param[out] screen_space_subdivision (output)
    /// screen_space_subdivision scnreen space subdicition bounding
    /// box list
    void get_screen_space_subdivision(
        std::vector<mi::math::Bbox_struct<mi::Uint32, 2> >& screen_space_subdivision);

    /// get number of screen space subdivision information for
    /// performance measurement
    ///
    /// \return number of spans
    virtual mi::Size get_nb_of_screen_space_subdivision();

    /// copy the internal buffer to data_ptr
    ///
    /// \param[in]  width  viewport width
    /// \param[in]  height viewport height
    /// \param[out] data_ptr canvas data head (\see mi::neuraylib::ITile::get_data())
    virtual void copy_pixel(mi::Sint32 width, mi::Sint32 height, void * data_ptr);

    /// get the class name
    /// \return class name
    virtual std::string get_class_name() const {
        return std::string("Span_renderer_no_gl");
    }

private:
    /// resize the internal buffer and initialize it
    ///
    /// \param[in] width  framebuffer width
    /// \param[in] height framebuffer height
    void resize_buffer(mi::Sint32 width, mi::Sint32 height);

    /// get the index of the position x,y
    ///
    /// \param[in] x pixel position of x
    /// \param[in] y pixel position of y
    /// \param[in] z pixel position of z, RGBA [0,3]
    /// \return index of the framebuffer
    inline mi::Sint32 get_index(mi::Sint32 x, mi::Sint32 y, mi::Sint32 z){
        assert((m_buffer_resolution.x > 0) && (m_buffer_resolution.y > 0));
        assert((0 <= x) && (x < m_buffer_resolution.x));
        assert((0 <= y) && (y < m_buffer_resolution.y));
        assert((0 <= z) && (z <  4)); // RGBA

        mi::Sint32 const idx = (4 * ((m_buffer_resolution.x * y) + x)) + z;
        assert((0 <= idx) && (idx < this->get_framebuffer_size()));

        return idx;
    }

    /// get the framebuffer size
    ///
    /// \return size of the framebuffer
    mi::Sint32 get_framebuffer_size() const
    {
        return m_buffer_resolution.x * m_buffer_resolution.y * 4;
    }

private:
    /// complete framebuffer
    mi::Uint8 * m_p_framebuffer;
    /// frame buffer size information
    mi::math::Vector_struct<mi::Sint32, 2> m_buffer_resolution;

    std::vector<mi::math::Bbox_struct<mi::Uint32, 2> > m_screen_space_subdivision;

    /// receive tile lock
    mi::base::Lock m_receive_tile_lock;
};

//----------------------------------------------------------------------
#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_SPAN_RENDERER_NO_GL_H
