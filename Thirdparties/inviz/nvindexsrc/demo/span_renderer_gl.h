/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief span renderer with OpenGL implementation

#ifndef MI_SUBSURFACE_VIEWER_OPENGL_BASED_SPAN_RENDERER_H
#define MI_SUBSURFACE_VIEWER_OPENGL_BASED_SPAN_RENDERER_H

#include <nv/index/iindex.h>
#include <mi/base/interface_implement.h>
#include <mi/base/lock.h>

#include <GL/glew.h>
#include <vector>

#include "span_renderer_if.h"



/// This class implement the span buffer rendering which enables user-defined rendering
/// of span data and shall be passed to the rendering call ('iindex_rendering->render').
/// This particular blends all spans into an OpenGL framebuffer.
class Span_renderer_gl : public Span_renderer_IF
{
public:
    /// constructor
    Span_renderer_gl();
    /// destructor
    virtual ~Span_renderer_gl();

public:
    //----------------------------------------------------------------------
    // Implemented of IIndex_canvas
    //----------------------------------------------------------------------

    virtual void prepare();

    virtual void receive_tile(
        mi::Uint8*                              span_buffer,
        mi::Uint32                              buffer_size,
        const mi::math::Bbox_struct<mi::Uint32, 2>& covered_area);

    virtual void receive_tile_blend(
        mi::Uint8*                              span_buffer,
        mi::Uint32                              buffer_size,
        const mi::math::Bbox_struct<mi::Uint32, 2>& covered_area);

    virtual void finish();

    virtual bool is_multi_thread_capable() const { return false; }

    virtual void set_buffer_resolution(
        const mi::math::Vector_struct<mi::Sint32, 2> & screen_resolution_size);

    virtual mi::math::Vector_struct<mi::Sint32, 2> get_buffer_resolution() const;

public:
    //----------------------------------------------------------------------
    // implementation of Span_renderer_IF
    //----------------------------------------------------------------------

    /// get screen space subdivision information for visualization
    ///
    /// \param[out] screen_space_subdivision (output)
    /// screen_space_subdivision scnreen space subdicition bounding
    /// box list
    virtual void get_screen_space_subdivision(
        std::vector<mi::math::Bbox_struct<mi::Uint32, 2> >& screen_space_subdivision);

    /// get number of screen space subdivision information for
    /// performance measurement
    ///
    /// \return number of spans
    virtual mi::Size get_nb_of_screen_space_subdivision();

    /// copy the OpenGL buffer to data_ptr
    ///
    /// \param[in]  width  viewport width
    /// \param[in]  height viewport height
    /// \param[out] data_ptr canvas data head (\see mi::neuraylib::ITile::get_data())
    virtual void copy_pixel(mi::Sint32 width, mi::Sint32 height, void * data_ptr);

    /// get the class name
    /// \return class name
    virtual std::string get_class_name() const {
        return std::string("Span_renderer_gl");
    }

private:
    GLboolean       m_pbo_enabled;
    GLuint          m_pbo;

    std::vector<mi::math::Bbox_struct<mi::Uint32, 2> > m_screen_space_subdivision;

    /// receive tile lock
    mi::base::Lock m_receive_tile_lock;

    /// screen resolution size
    mi::math::Vector_struct<mi::Sint32, 2> m_buffer_resolution;
};

#endif // MI_SUBSURFACE_VIEWER_OPENGL_BASED_SPAN_RENDERER_H
