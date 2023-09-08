/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief span renderer interface for OpenGL and non-OpenGL implementation

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_SPAN_RENDERER_IF_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_SPAN_RENDERER_IF_H

#include <nv/index/iindex.h>
#include <mi/base/interface_implement.h>

#include <vector>
#include <string>

namespace nv {
namespace index_common {
class String_dict;
}} // namespace

/// Span renderer application side interface.
///
/// We have now two kind of span renderer. Both has the same
/// interface, but the different implementation. This is the interface
/// for them.
class Span_renderer_IF :
        public mi::base::Interface_implement<nv::index::IIndex_canvas>
{
public:
    /// constructor
    Span_renderer_IF()
      : m_use_parallel_rendering(true)
    {
        m_background_color.r = 0.0f;
        m_background_color.g = 0.0f;
        m_background_color.b = 0.0f;
        m_background_color.a = 1.0f;
    }
    /// destructor
    virtual ~Span_renderer_IF()
    {
        // empty
    }

public:
    /// get screen space subdivision information for visualization
    ///
    /// \param[out] screen_space_subdivision (output)
    /// screen_space_subdivision scnreen space subdicition bounding
    /// box list
    virtual void get_screen_space_subdivision(
        std::vector<mi::math::Bbox_struct<mi::Uint32, 2> >& screen_space_subdivision) = 0;

    /// get number of screen space subdivision information for
    /// performance measurement
    ///
    /// \return number of spans
    virtual mi::Size get_nb_of_screen_space_subdivision() = 0;

    /// set current background color
    /// \param[in] background_color buffer background color 
    virtual void set_background_color(const mi::math::Color_struct & background_color)
    {
        m_background_color = background_color;
    }

    /// get current background color
    /// \return current background color 
    virtual mi::math::Color_struct get_background_color() const
    {
        return m_background_color;
    }

    /// copy the internal buffer to data_ptr
    ///
    /// \param[in]  width  viewport width
    /// \param[in]  height viewport height
    /// \param[out] data_ptr canvas data head (\see mi::neuraylib::ITile::get_data())
    virtual void copy_pixel(mi::Sint32 width, mi::Sint32 height, void * data_ptr) = 0;

    /// get the class name
    /// \return class name
    virtual std::string get_class_name() const = 0;

    /// Controls whether to use parallel rendering, if supported by the span renderer. If it is not
    /// supported, this setting will have no effect.
    virtual void set_use_parallel_rendering(bool enable)
    {
        m_use_parallel_rendering = enable;
    }

    /// use parallel rendering?
    /// \return true when run in parallel.
    virtual bool get_use_parallel_rendering() const
    {
        return m_use_parallel_rendering;
    }

private:
    bool m_use_parallel_rendering;
    mi::math::Color_struct m_background_color;
};

//----------------------------------------------------------------------
#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_SPAN_RENDERER_IF_H
