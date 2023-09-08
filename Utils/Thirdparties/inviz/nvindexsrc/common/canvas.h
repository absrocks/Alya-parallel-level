/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Example implementation of Time and Canvas

#ifndef NVIDIA_INDEX_BIN_COMMON_CANVAS_H
#define NVIDIA_INDEX_BIN_COMMON_CANVAS_H

#include <mi/base/types.h>
#include <mi/base/handle.h>
#include <mi/base/interface_implement.h>
#include <mi/base/lock.h>

#include <mi/math/vector.h>

#include <mi/neuraylib/icanvas.h>
#include <mi/neuraylib/itile.h>

#include <cassert>

namespace nv {
namespace index_common {

/// Example Tile implementation used by Canvas
class Tile : public mi::base::Interface_implement<mi::neuraylib::ITile>
{
  public:
    /// Constructor. Creates a tile of the given width and height.
    Tile(
        mi::Uint32 width,
        mi::Uint32 height)
      : m_width(width),
        m_height(height),
        m_data(0)
    {
        m_data = new mi::Uint8[m_width * m_height * 4];
    }

    /// Destructor
    virtual ~Tile()
    {
        delete[] m_data;
        m_data = 0;
    }

public:
    // Implement the interface of mi::neuraylib::ITile
    virtual void get_pixel(
        mi::Uint32          x_offset,
        mi::Uint32          y_offset,
        mi::Float32*        floats) const
    {
        mi::Uint8* position = &m_data[(x_offset + y_offset * m_width) * 4];
        floats[0] = static_cast<mi::Float32>(position[0]) * mi::Float32(1.0/255.0);
        floats[1] = static_cast<mi::Float32>(position[1]) * mi::Float32(1.0/255.0);
        floats[2] = static_cast<mi::Float32>(position[2]) * mi::Float32(1.0/255.0);
        floats[3] = static_cast<mi::Float32>(position[3]) * mi::Float32(1.0/255.0);
    }
    virtual void set_pixel(
        mi::Uint32          x_offset,
        mi::Uint32          y_offset,
        const mi::Float32*  floats)
    {
        mi::Uint8* position = &m_data[(x_offset + y_offset * m_width) * 4];
        position[0] = static_cast<mi::Uint8>( mi::math::clamp( floats[0], 0.0f, 1.0f) * 255.f);
        position[1] = static_cast<mi::Uint8>( mi::math::clamp( floats[1], 0.0f, 1.0f) * 255.f);
        position[2] = static_cast<mi::Uint8>( mi::math::clamp( floats[2], 0.0f, 1.0f) * 255.f);
        position[3] = static_cast<mi::Uint8>( mi::math::clamp( floats[3], 0.0f, 1.0f) * 255.f);
    }

    virtual const char* get_type()          const   { return "Rgba"; }
    virtual mi::Uint32 get_resolution_x()   const   { return m_width; }
    virtual mi::Uint32 get_resolution_y()   const   { return m_height; }
    virtual const void* get_data()          const   { return m_data; }
    virtual void* get_data()                        { return m_data; }

private:
    // unused default constructor
    Tile();
private:
    mi::Uint32 m_width;
    mi::Uint32 m_height;
    mi::Uint8* m_data;
};

/// Example canvas class
class Canvas : public mi::base::Interface_implement<mi::neuraylib::ICanvas>
{
public:
    /// Constructor. Creates a Canvas of the given width and height.
    Canvas(
        mi::Uint32 width,
        mi::Uint32 height)
        : 
        m_width(0),
        m_height(0),
        m_gamma(2.2f),
        m_tile()
    {
        m_tile = 0;
        resize_buffer(width, height);
    }

    // Destructor
    virtual ~Canvas()
    {
        m_tile = 0;
    }

    // Implement the interface of mi::neuraylib::ICanvas
    mi::Uint32 get_resolution_x()       const { return m_width;     }
    mi::Uint32 get_resolution_y()       const { return m_height;    }
    mi::Uint32 get_tile_resolution_x()  const { return m_width;     }
    mi::Uint32 get_tile_resolution_y()  const { return m_height;    }
    mi::Uint32 get_tiles_size_x()       const { return 1;           }
    mi::Uint32 get_tiles_size_y()       const { return 1;           }
    mi::Uint32 get_layers_size()        const { return 1;           }

    const mi::neuraylib::ITile* get_tile(
        mi::Uint32 pixel_x,
        mi::Uint32 pixel_y,
        mi::Uint32 layer = 0) const
    {
        //TODO: LOCK?
        m_tile->retain();
        return m_tile.get();
    }

    mi::neuraylib::ITile* get_tile(
        mi::Uint32 pixel_x,
        mi::Uint32 pixel_y,
        mi::Uint32 layer = 0)
    {
        //TODO: LOCK?
        m_tile->retain();
        return m_tile.get();
    }

    const char* get_type() const { return "Rgba"; }
    mi::Float32 get_gamma() const { return m_gamma; }
    void set_gamma(mi::Float32 gamma) { m_gamma = gamma; }

    /// resize the buffer
    /// \param[in] width  new buffer width, must be > 0
    /// \param[in] height new buffer height, must be > 0
    void resize_buffer(mi::Uint32 width, mi::Uint32 height)
    {
        assert(width  > 0);
        assert(height > 0);

        if ((m_width == width) && (m_height == height))
        {
            return;             // same size to the current size, no resize.
        }

        // clean up and re-allocation
        m_tile   = 0;
        m_width  = width;        
        m_height = height;
        m_tile = new Tile(m_width, m_height);
    }

private:
    /// unused default constructor
    Canvas();

private:
    mi::Uint32 m_width;
    mi::Uint32 m_height;
    mi::Float32 m_gamma;

    /// The tiles of this canvas
    mi::base::Handle<Tile>  m_tile;
};

/// Encapsulate triple buffering of Canvas buffers.
///
/// By making sure a buffer only ever gets accessed by one thread, locking is only necessary when
/// swapping buffers, which should be very fast.
class Canvas_buffers : public mi::base::Interface_implement<mi::base::IInterface>
{
  public:
    /// Constructs the internal canvas buffers with the given size.
    Canvas_buffers(
        mi::Uint32 width,
        mi::Uint32 height)
        :
        m_render_frame_number(-1)
    {
        update_buffer_size(width, height);
    }

    /// destructor
    virtual ~Canvas_buffers()
    {
        m_render_canvas   = 0;
        m_finished_canvas = 0;
        m_encoder_canvas  = 0;
    }

    /// Changes the size of the internal canvas buffers.
    void update_buffer_size(
        mi::Uint32 width,
        mi::Uint32 height)
    {
        m_render_canvas   = new Canvas(width, height);
        m_finished_canvas = new Canvas(width, height);
        m_encoder_canvas  = new Canvas(width, height);
        
        m_new_finished_canvas_available = false;
    }

    /// Should be called by the rendering thread after it has finished rendering an image.
    void rendering_finished()
    {
        mi::base::Lock::Block block(&m_lock_canvas);

        // Swap buffers and tell the video encoder that there is a new image available
        std::swap(m_render_canvas, m_finished_canvas);
        m_new_finished_canvas_available = true;
        ++m_render_frame_number;
    }

    /// Should be called by the video encoder thread before it accesses and starts encoding an
    /// image.
    bool prepare_encoder_canvas()
    {
        mi::base::Lock::Block block(&m_lock_canvas);

        if (m_new_finished_canvas_available)
        {
            // The video encoder is now responsible for encoding the new image
            std::swap(m_encoder_canvas, m_finished_canvas);
            m_new_finished_canvas_available = false;
            return true;
        }
        else
        {
            // The finished canvas was not changed
            return false;
        }
    }

    /// Returns the canvas the rendering thread should write to.
    /// The rendering thread should call rendering_finished() when the image is ready to be
    /// delivered to the video encoder.
    Canvas* get_render_canvas()
    {
        if(m_render_canvas)
            m_render_canvas->retain();
        
        return m_render_canvas.get();
    }

    /// Returns the canvas the video encoder should read from.
    /// The video encoder thread should call prepare_encoder_canvas() before calling this method, to
    /// make sure it gets the most recent image.
    Canvas* get_encoder_canvas()
    {
        if(m_encoder_canvas)
            m_encoder_canvas->retain();
        
        return m_encoder_canvas.get();
    }

    /// Get the render canvas resolution.
    /// 
    /// \return render canvas resolution
    mi::math::Vector<mi::Uint32, 2> get_render_canvas_resolution() const
    {
        return mi::math::Vector<mi::Uint32, 2>(m_render_canvas->get_resolution_x(),
                                               m_render_canvas->get_resolution_y());
    }

    /// Get the rendered frame number
    ///
    /// \return the current render frame number
    mi::Sint64 get_render_frame_number() const
    {
        return m_render_frame_number;
    }

private:
    /// unused default constructor
    Canvas_buffers();

private:
    /// The render canvas is only written to by the render thread.
    mi::base::Handle<Canvas>                  m_render_canvas;

    /// The finished canvas always contains a complete rendered image, no read or write operations
    /// happen to it directly.
    mi::base::Handle<Canvas>                  m_finished_canvas;

    /// The encoder canvas is only read by the video encoder thread.
    mi::base::Handle<Canvas>                  m_encoder_canvas;

    /// The lock prevents access to the canvases while they are swapped.
    mi::base::Lock                            m_lock_canvas;

    /// Flag tells that the renderer thread has finished an image and that the video encoder has not
    /// yet processed it.
    bool                                      m_new_finished_canvas_available;
    
    /// rendered frame number
    mi::Sint64                                m_render_frame_number;
};

}} // namespace nv::index_common

#endif // NVIDIA_INDEX_BIN_COMMON_CANVAS_H
