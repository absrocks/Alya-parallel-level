/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief span renderer with no-OpenGL functions

#ifndef NVIDIA_INDEX_BIN_COMMON_IMAGE_FILE_CANVAS_H
#define NVIDIA_INDEX_BIN_COMMON_IMAGE_FILE_CANVAS_H

#include <mi/base/interface_implement.h>
#include <mi/base/lock.h>

#include <nv/index/iindex.h>
#include <nv/index/iviewport.h>

#include <vector>
#include <string>
#include <cassert>

#include "forwarding_logger.h"
#include "ppm_io.h"

namespace nv {
namespace index_common {

/// This class implement the span buffer rendering which enables
/// user-defined rendering of span data and shall be passed to the
/// rendering call ('iindex_rendering->render').
///
/// This particular blends all spans into an user framebuffer without
/// any OpenGL functionality.
class Image_file_canvas : public mi::base::Interface_implement<nv::index::IIndex_canvas>
{
public:
    /// constructor
    Image_file_canvas()
        :
        m_p_framebuffer(0),
        m_width(-1),
        m_height(-1),
        m_screen_space_subdivision(),
        m_save_fname("imagefile.ppm")
    {
        // empty
    }

    /// constructor
    ///
    /// \param[in] width  viewport width
    /// \param[in] height viewport height
    /// \param[in] fname  save file name
    Image_file_canvas(mi::Sint32 width, mi::Sint32 height, std::string const & fname)
        :
        m_p_framebuffer(0),
        m_width(-1),
        m_height(-1),
        m_screen_space_subdivision(),
        m_save_fname("imagefile.ppm")
    {
        this->set_save_filename(fname);
        this->set_buffer_resolution(mi::math::Vector<mi::Sint32, 2>(width, height)); // FIXME
    }

    /// destructor
    virtual ~Image_file_canvas()
    {
        DEBUG_LOG << "Image_file_canvas::~Image_file_canvas()";
    }

    /// shutdown this
    void shutdown()
    {
        if(m_p_framebuffer != 0)
        {
            delete [] m_p_framebuffer;
            m_p_framebuffer = 0;
        }
        m_width  = -1;
        m_height = -1;
        m_screen_space_subdivision.clear();
    }


    /// set save filename
    /// \param[in] fname filename to save
    void set_save_filename(std::string const & fname)
    {
        m_save_fname = fname;
    }

    /// get save filename
    /// \return current save filename
    std::string get_save_filename() const
    {
        return m_save_fname;
    }

    /// get width
    /// \return width size
    mi::Sint32 get_width() const
    {
        return m_width;
    }

    /// get height
    /// \return height size
    mi::Sint32 get_height() const
    {
        return m_height;
    }

    /// get ppm color buffer.
    /// Note: expensive
    ///
    /// \param[out] ppm_color_buf (output) ppm color buffer
    /// \return true when succeeded.
    bool get_ppm_color_buffer(
        std::vector< mi::math::Color_struct > & ppm_color_buf) const
    {
        mi::Sint32 const width  = this->get_width();
        mi::Sint32 const height = this->get_height();
        if((width <= 0) || (height <= 0) || (m_p_framebuffer == 0))
        {
            return false;
        }

        bool const is_success =
            get_ppm_color_buffer_from_pixel_buffer(width, height, m_p_framebuffer,
                                                   ppm_color_buf);
        return is_success;
    }

public:
    //------------------------------------------------------------
    // IIndex_canvas implementation
    //------------------------------------------------------------

    virtual void prepare()
    {
        m_screen_space_subdivision.resize(0);
        // initialize the framebuffer
        memset(m_p_framebuffer, 0, this->get_framebuffer_size());
    }

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
        
        if(m_p_framebuffer == 0)
        {
            ERROR_LOG << "Image_file_canvas: buffer has no size. Have you called set_buffer_resolution()?";
            return;
        }
        
        // No blending necessary: Since the incoming image data uses pre-multiplied alpha (which means
        // that the RGB values were already multiplied by the alpha), a simple copy is sufficient.
        mi::Uint8 * p_span_buffer_iter = span_buffer;
        const size_t copy_byte = x_width * 4;
        
        for (mi::Uint32 y = covered_area.min.y; y < covered_area.max.y; ++y)
        {
            // Copy scanline with memcpy
            assert(covered_area.min.x < covered_area.max.x);
            mi::Uint8 * p_dest = m_p_framebuffer + get_index(covered_area.min.x, y, 0);
            memcpy(p_dest, p_span_buffer_iter, copy_byte);
            p_span_buffer_iter += copy_byte;
        }
    }
    
    virtual void receive_tile_blend(
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
        
        if(m_p_framebuffer == 0)
        {
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

    virtual void finish()
    {
        assert(m_p_framebuffer != 0);
        if (this->get_save_filename().empty())
        {
            DEBUG_LOG << "Image_file_canvas::finish: empty output file name... no image saved.";
            return;
        }
        
        bool const ret = save_pixel_buffer_to_ppm_file(m_width, m_height,
                                                       reinterpret_cast< mi::Uint8 *>(m_p_framebuffer),
                                                       this->get_save_filename());
        if (ret)
        {
            DEBUG_LOG << "Image_file_canvas::finish: save image ["
                      << this->get_save_filename() + "] ... done.";
        }
        else
        {
            ERROR_LOG << "Image_file_canvas::finish: save image ["
                      << this->get_save_filename() << "] ... failed.";
        }
    }
    
    virtual bool is_multi_thread_capable() const
    {
        // Both serial and parallel rendering is supported, so let the user decide
        return true;
    }

    virtual void set_buffer_resolution(
        const mi::math::Vector_struct<mi::Sint32, 2> & screen_resolution_size)
    {
        assert(screen_resolution_size.x > 0);
        assert(screen_resolution_size.y > 0);

        this->resize_buffer(screen_resolution_size.x, screen_resolution_size.y);
     
        // initialize the framebuffer
        memset(m_p_framebuffer, 0, this->get_framebuffer_size());
    }

    virtual mi::math::Vector_struct<mi::Sint32, 2> get_buffer_resolution() const
    {
        // TODO: FIXME
        const mi::math::Vector_struct<mi::Sint32, 2> ret = {m_width, m_height,};
        return ret;
    }

public:
    /// get screen space subdivision information for visualization
    ///
    /// \param[out] screen_space_subdivision (output)
    /// screen_space_subdivision scnreen space subdicition bounding
    /// box list
    void get_screen_space_subdivision(
        std::vector<mi::math::Bbox_struct<mi::Uint32, 2> >& screen_space_subdivision) const
    {
        screen_space_subdivision = m_screen_space_subdivision;
    }
    
    /// get number of screen space subdivision information for
    /// performance measurement
    ///
    /// \return number of spans
    virtual mi::Size get_nb_of_screen_space_subdivision() const
    {
        return m_screen_space_subdivision.size();
    }
    
    /// copy the internal buffer to data_ptr
    ///
    /// \param[in]  width  viewport width
    /// \param[in]  height viewport height
    /// \param[out] data_ptr canvas data head (\see mi::neuraylib::ITile::get_data())
    virtual void copy_pixel(mi::Sint32 width, mi::Sint32 height,
                            void * data_ptr)
    {
        mi::Uint8 * uchar_ptr = reinterpret_cast< mi::Uint8 * >(data_ptr);

        // We can optimize the x, y == 0 case with memcopy like this.
        // mi::Sint32 const xmin = 0;
        // mi::Sint32 const xmax = width;
        // mi::Sint32 const ymin = 0;
        // mi::Sint32 const ymax = height;

        // copy whole buffer
        memcpy(uchar_ptr, m_p_framebuffer, this->get_framebuffer_size());
    }
    
    /// get the class name
    /// \return class name
    virtual std::string get_class_name() const
    {
        return std::string("Image_file_canvas");
    }
    
private:
    /// resize the internal buffer and initialize it
    ///
    /// \param[in] width  framebuffer width
    /// \param[in] height framebuffer height
    void resize_buffer(mi::Sint32 width, mi::Sint32 height)
    {
        assert(width  > 0);
        assert(height > 0);
        
        if ((width == m_width) && (height == m_height))
        {
            // no need to resize
            return;
        }
        
        if (m_p_framebuffer != 0)
        {
            delete [] m_p_framebuffer;
            m_p_framebuffer = 0;
        }
        
        // first need to set the width and height for get_framebuffer_size()
        m_width  = width;
        m_height = height;
        
        m_p_framebuffer = new mi::Uint8[this->get_framebuffer_size()];
    }
    
    /// get the index of the position x,y
    ///
    /// \param[in] x pixel position of x
    /// \param[in] y pixel position of y
    /// \param[in] z pixel position of z, RGBA [0,3]
    /// \return index of the framebuffer
    inline mi::Sint32 get_index(mi::Sint32 x, mi::Sint32 y, mi::Sint32 z)
    {
        assert((m_width > 0) && (m_height > 0));
#ifdef DEBUG
        if(!(((0 <= x) && (x < m_width))  && 
             ((0 <= y) && (y < m_height)) && 
             ((0 <= z) && (z <  4))))
        {
            ERROR_LOG << "get_index: accessing outside of the boundary. ("
                      << x << ", " << y << ", " << z << ") of [" << m_width << ", " << m_height << "]";
        }
#endif
        assert((0 <= x) && (x < m_width));
        assert((0 <= y) && (y < m_height));
        assert((0 <= z) && (z <  4)); // RGBA
        
        mi::Sint32 const idx = (4 * ((m_width * y) + x)) + z;
        assert((0 <= idx) && (idx < this->get_framebuffer_size()));
        
        return idx;
    }
    
    /// get the framebuffer size
    ///
    /// \return size of the framebuffer
    mi::Sint32 get_framebuffer_size() const
    {
        return m_width * m_height * 4;
    }
    
private:
    /// complete framebuffer
    mi::Uint8 * m_p_framebuffer;
    /// frame buffer size information: width  FIXME: use m_buffer_resolution
    mi::Sint32  m_width;
    /// frame buffer size information: height
    mi::Sint32  m_height;
    
    std::vector<mi::math::Bbox_struct<mi::Uint32, 2> > m_screen_space_subdivision;
    
    /// save file name of this canvas
    std::string m_save_fname;
    
    /// receive tile lock
    mi::base::Lock m_receive_tile_lock;
};

//----------------------------------------------------------------------
}} // namespace nv::index_common
#endif // NVIDIA_INDEX_BIN_COMMON_IMAGE_FILE_CANVAS_H
