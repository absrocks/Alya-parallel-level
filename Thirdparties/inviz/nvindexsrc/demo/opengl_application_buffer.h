/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief opengl application buffer

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_OPENGL_APPLICATION_BUFFER_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_OPENGL_APPLICATION_BUFFER_H

#include <nv/index/iopengl_application_buffer.h>

#include <cassert>
#include <string>

//----------------------------------------------------------------------
/// opengl application buffer
///
/// This implementation has no dependency to OpenGL.
class Opengl_application_buffer :
        public mi::base::Interface_implement<nv::index::IOpengl_application_buffer>
{
public:
    /// default constructor
    Opengl_application_buffer();
    /// destructor
    virtual ~Opengl_application_buffer();

public:
    /// \addgroup iopegl_application_buffer_implementation iopegl_application_buffer_implementation
    /// \brief implementation of IOpengl_application_buffer
    /// @{

    /// get the buffer resolution.
    ///
    /// \return resolution (width,height) in pixels. may (-1, -1) if not initialized.
    virtual mi::math::Vector_struct< mi::Sint32, 2 > get_resolution() const;

    /// get the pointer of Z-buffer. The same structure to OpenGL
    /// Z-buffer. May 0 if not initialized.
    /// This returns writable raw pointer. Use with care.
    /// Depth value is [0, maxuint].
    /// \return pointer to the top address of mi::Uint32[pixel_count]
    virtual mi::Uint32* get_z_buffer_ptr();

    /// @}

    /// Set Z-buffer precision (bits)
    ///
    /// \param[in] precision of z-buffer in bits
    void set_z_buffer_precision(mi::Uint32 precision);

    /// Get Z-buffer precision (bits)
    ///
    /// \return precision of z-buffer in bits
    mi::Uint32 get_z_buffer_precision() const;

    /// resize the buffer. usually causing memory allocation.
    ///
    /// \param[in] new_resolution new resolution of this buffer.
    void resize_buffer(mi::math::Vector_struct< mi::Sint32, 2 > const & new_resolution);

    /// clear buffer
    void clear_buffer();

    /// is buffer allocated
    /// \return true when the buffer is allocated.
    bool is_buffer_allocated() const ;

public:
    /// get pixel index. Assume buffer size is less than 2G.
    static inline mi::Sint32 get_pixel_index(mi::Sint32 width,
                                             mi::Sint32 height,
                                             mi::Sint32 channel,
                                             mi::Sint32 x,
                                             mi::Sint32 y,
                                             mi::Sint32 rgba)
    {
        assert((0 <= x)    && (x    < width));
        assert((0 <= y)    && (y    < height));
        assert((0 <= rgba) && (rgba < channel));
        return ((channel * (x + (y * width))) + rgba);
    }

public:
    /// debug: dump RGBA and Z to ppm image (no alpha)
    ///
    /// \param[in] fname output file basename. {fname}.rgba.ppm, and
    /// {fname}.z.ppm file will be created.
    /// \return true when succeeded.
    bool debug_write_buffer_to_file(std::string const & fname);

private:
    /// deallocate memory
    void delete_memory();

    /// allocate memory accuoding to the new_resolution.
    /// Deallocation/initialization must be done before this call.
    /// \param[in] resolution new resolution of this buffer.
    void allocate_memory(mi::math::Vector< mi::Sint32, 2 > const & new_resolution);

private:
    /// buffer resolution
    mi::math::Vector< mi::Sint32, 2 > m_resolution;
    /// Z-buffer
    mi::Uint32* m_z_buffer;
    /// Z-buffer precision (bits)
    mi::Uint32 m_z_buffer_precision;

private:
    /// copy constructor. prohibit until proved useful.
    Opengl_application_buffer(Opengl_application_buffer const &);
    /// operator=. prohibit until proved useful.
    Opengl_application_buffer const & operator=(Opengl_application_buffer const &);
};

//----------------------------------------------------------------------


#endif // #ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_OPENGL_APPLICATION_BUFFER_H
