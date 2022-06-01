/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "opengl_application_buffer.h"

#include "common/colormap_io.h"
#include "common/forwarding_logger.h"
#include "common/ppm_io.h"

#include "utilities.h"

#include <vector>

//----------------------------------------------------------------------
Opengl_application_buffer::Opengl_application_buffer()
  : m_resolution(-1, -1),
    m_z_buffer(0),
    m_z_buffer_precision(24)
{
    // empty
}

//----------------------------------------------------------------------
Opengl_application_buffer::~Opengl_application_buffer()
{
    this->delete_memory();
    DEBUG_LOG << "Opengl_application_buffer::~Opengl_application_buffer() called.";
}

//----------------------------------------------------------------------
mi::math::Vector_struct< mi::Sint32, 2 > Opengl_application_buffer::get_resolution() const
{
    return m_resolution;
}

//----------------------------------------------------------------------
mi::Uint32* Opengl_application_buffer::get_z_buffer_ptr()
{
    return m_z_buffer;
}

//----------------------------------------------------------------------
void Opengl_application_buffer::set_z_buffer_precision(mi::Uint32 precision)
{
    if (precision != 24 && precision != 32)
    {
        ERROR_LOG << "Trying to set unsupported z-buffer precision " << precision
                  << ". Please check X server status when the precision is reported as 0."; 
        return;
    }

    m_z_buffer_precision = precision;
}

//----------------------------------------------------------------------
mi::Uint32 Opengl_application_buffer::get_z_buffer_precision() const
{
    return m_z_buffer_precision;
}

//----------------------------------------------------------------------
void Opengl_application_buffer::resize_buffer(
    mi::math::Vector_struct< mi::Sint32, 2 > const & new_resolution)
{
    mi::math::Vector< mi::Sint32, 2 > new_resolution_vec(new_resolution);
    if(m_resolution != new_resolution_vec){
        this->delete_memory();
        this->allocate_memory(new_resolution_vec);
        this->clear_buffer();
    }
}

//----------------------------------------------------------------------
void Opengl_application_buffer::clear_buffer()
{
    if(m_z_buffer == 0){
        return;
    }

    assert(this->is_buffer_allocated());

    mi::Sint64 const pixel_count = m_resolution.x * m_resolution.y;
    assert(pixel_count > 0);

    memset(m_z_buffer, 0, (pixel_count * sizeof(mi::Uint32)));
}

//----------------------------------------------------------------------
bool Opengl_application_buffer::is_buffer_allocated() const
{
    if((m_resolution.x > 0) && (m_resolution.y > 0) &&
       (m_z_buffer    != 0))
    {
        return true;
    }
    return false;
}

//----------------------------------------------------------------------
bool Opengl_application_buffer::debug_write_buffer_to_file(
    std::string const & fname)
{
    if(fname.size() == 0){
        ERROR_LOG << "Opengl_application_buffer::debug_write_buffer_to_file: empty filename.";
        return false;
    }
    if(!this->is_buffer_allocated()){
        ERROR_LOG << "Opengl_application_buffer::debug_write_buffer_to_file: no buffer. "
                  << "(no resized.)";
        return false;
    }

    mi::math::Vector< mi::Sint32, 2 > const res = this->get_resolution();
    assert((res.x > 0) && (res.y > 0));
    assert(static_cast< mi::Sint64 >(res.x) * static_cast< mi::Sint64 >(res.y) <
           mi::base::numeric_traits< mi::Sint32 >::max());
    mi::Sint32 const pixel_count = res.x * res.y;
    assert(pixel_count > 0);

/*    
    // write RGBA buffer
    mi::Float64 const color_scale = 1.0/255.0;
    std::vector< mi::math::Color_struct > ppm_rgba_buf(pixel_count);
    for(mi::Sint32 y = 0; y < res.y; ++y){
        for(mi::Sint32 x = 0; x < res.x; ++x){
            mi::Sint32 const gl_rgba_idx =
                Opengl_application_buffer::get_pixel_index(res.x, res.y, 4, x, y, 0);
            mi::Sint32 const ppm_idx = (res.y * y) + x;
            mi::math::Vector< mi::Float64, 3 > dcol(
                color_scale * static_cast< mi::Float64 >(m_rgba_buffer[gl_rgba_idx + 0]),
                color_scale * static_cast< mi::Float64 >(m_rgba_buffer[gl_rgba_idx + 1]),
                color_scale * static_cast< mi::Float64 >(m_rgba_buffer[gl_rgba_idx + 2]));
            ppm_rgba_buf.at(ppm_idx).r = static_cast< mi::Float32 >(dcol.x);
            ppm_rgba_buf.at(ppm_idx).g = static_cast< mi::Float32 >(dcol.y);
            ppm_rgba_buf.at(ppm_idx).b = static_cast< mi::Float32 >(dcol.z);
        }
    }

    std::string const rgba_fname = fname + ".rgba.ppm";
    bool ret = write_ppm(rgba_fname, ppm_rgba_buf, res.y, res.x);
    if(!ret){
        ERROR_LOG << "fail to save [" << rgba_fname << "]";
        return false;
    }
    
    INFO_LOG << "Writing rgba buffer to file: " << rgba_fname;
*/

    // write Z buffer
    // OpenGL depth buffer read scales to 32bit uint value
    mi::Float64 const depth_scale = 1.0 / (m_z_buffer_precision == 24 ? 16777215.0 : 4294967295.0);
    std::vector< mi::math::Color_struct > ppm_z_buf(pixel_count);
    
    for(mi::Sint32 y = 0; y < res.y; ++y)
    {
        for(mi::Sint32 x = 0; x < res.x; ++x)
        {
            mi::Sint32 const gl_z_idx = Opengl_application_buffer::get_pixel_index(res.x, res.y, 1, x, y, 0);
            mi::Sint32 const ppm_idx = (res.y * y) + x;
            mi::Uint32 depth = static_cast< mi::Uint32 >(static_cast< mi::Float64 >(m_z_buffer[gl_z_idx]) * depth_scale);
            
            assert((0.0 <= depth) && (depth <= 1.0));
            
            ppm_z_buf.at(ppm_idx).r = static_cast<mi::Float32>(depth);
            ppm_z_buf.at(ppm_idx).g = static_cast<mi::Float32>(depth);
            ppm_z_buf.at(ppm_idx).b = static_cast<mi::Float32>(depth);
        }
    }

    std::string const z_fname = fname + ".z.ppm";
    bool ret = nv::index_common::write_ppm(z_fname, ppm_z_buf, res.y, res.x);
    INFO_LOG << "Writing z buffer to file: " << z_fname;
    if(!ret){
        ERROR_LOG << "failed to save [" << z_fname << "]";
        return false;
    }

    return true;
}

//----------------------------------------------------------------------
void Opengl_application_buffer::delete_memory()
{
    m_resolution = mi::math::Vector< mi::Sint32, 2 >(-1, -1);
    if(m_z_buffer != 0){
        delete [] m_z_buffer;
        m_z_buffer = 0;
    }
}

//----------------------------------------------------------------------
void Opengl_application_buffer::allocate_memory(
    mi::math::Vector< mi::Sint32, 2 > const & new_resolution)
{
    if((new_resolution.x <= 0) || (new_resolution.y <= 0)){
        ERROR_LOG << "Opengl_application_buffer::allocate_memory: "
                  << "cannot allocate negative memory size.";
        return;
    }

    assert(m_z_buffer    == 0);

    m_resolution = new_resolution;
    mi::Sint64 const pixel_count = m_resolution.x * m_resolution.y;
    assert(pixel_count > 0);

    m_z_buffer = new mi::Uint32[pixel_count];
}

//----------------------------------------------------------------------

