/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "image_tile.h"

#include "common/colormap_io.h"
#include "common/ppm_io.h"

//----------------------------------------------------------------------
Image_tile::Image_tile()
    :
    m_ij_bbox(),
    m_color_buffer()
{
    // empty
}

//----------------------------------------------------------------------
Image_tile::~Image_tile()
{
    m_ij_bbox.     clear();
    m_color_buffer.clear();
}

//----------------------------------------------------------------------
Image_tile::Image_tile(Image_tile const & rhs)
{
    m_ij_bbox      = rhs.m_ij_bbox;
    m_color_buffer = rhs.m_color_buffer;
}

//----------------------------------------------------------------------
Image_tile& Image_tile::operator=(Image_tile const & rhs)
{
    if(this != & rhs){
        m_ij_bbox      = rhs.m_ij_bbox;
        m_color_buffer = rhs.m_color_buffer;
    }
    return *this;
}

//----------------------------------------------------------------------
bool Image_tile::resize_buffer(mi::math::Bbox<mi::Sint64, 2 > const & ij_bbox)
{
    if(!ij_bbox.is_plane()){
        DEBUG_LOG << "Image_tile::resize_buffer: invalid ij_bbox: " << ij_bbox;
        return false;
    }

    m_ij_bbox = ij_bbox;

    mi::math::Vector< mi::Sint64, 2 > buf_size = this->get_buffer_size();
    // DEBUG_LOG << "resize to " << buf_size;
    m_color_buffer.resize(buf_size.x * buf_size.y);

    return true;
}

//----------------------------------------------------------------------
bool Image_tile::is_valid_buffer() const
{
    if(m_ij_bbox.is_plane()){
        return true;
    }
    return false;
}

//----------------------------------------------------------------------
void Image_tile::clear_buffer(mi::math::Color const & col)
{
    assert(this->is_valid_buffer());
    size_t buf_size = m_color_buffer.size();

    for(size_t i = 0; i < buf_size; ++i){
        m_color_buffer.at(i) = col;
    }
}

//----------------------------------------------------------------------
mi::math::Vector< mi::Sint64, 2 > Image_tile::get_buffer_size() const
{
    mi::math::Vector< mi::Sint64, 2 > buf_size = this->m_ij_bbox.max - m_ij_bbox.min;
    return buf_size;
}

//----------------------------------------------------------------------
mi::math::Bbox<mi::Sint64, 2 > Image_tile::get_bounding_box() const
{
    return m_ij_bbox;
}

//----------------------------------------------------------------------
void Image_tile::put_color(mi::math::Vector< mi::Sint64, 2 > const & ij, mi::math::Color const & col)
{
    m_color_buffer.at(this->get_idx(ij)) = col;
}

//----------------------------------------------------------------------
mi::math::Color Image_tile::get_color(mi::math::Vector< mi::Sint64, 2 > const & ij) const
{
    return m_color_buffer.at(this->get_idx(ij));
}

//----------------------------------------------------------------------
bool Image_tile::copy_partial_image_tile(
    Image_tile const & src_image_tile,
    mi::math::Bbox<   mi::Sint64, 2 > const & copy_src_bbox,
    mi::math::Vector< mi::Sint64, 2 > const & dst_origin
    )
{
    // check the parameter validity
    if(!src_image_tile.is_valid_buffer()){
        ERROR_LOG << "illegal source tile image.";
        return false;
    }
    if(!this->is_valid_buffer()){
        ERROR_LOG << "illegal this tile image.";
        return false;
    }
    if(!nv::index_common::bbox_contains(src_image_tile.get_bounding_box(), copy_src_bbox)){
        ERROR_LOG << "copy region is out of source image tile range. src_bbox: "
                  << src_image_tile.get_bounding_box()
                  << ", copy_src_bbox: " << copy_src_bbox;
        return false;
    }
    if(!this->is_inside(dst_origin)){
        ERROR_LOG << "destination origin is out of range. dst_raw_bbox: "
                  << this->get_bounding_box()
                  << ", dst_origin: " << dst_origin;
        return false;
    }

    mi::math::Bbox< mi::Sint64, 2 >   const dst_raw_bbox   = this->get_bounding_box();
    mi::math::Bbox< mi::Sint64, 2 >   const src_raw_bbox   = src_image_tile.get_bounding_box();

    mi::math::Vector< mi::Sint64, 2 > const dst_size_vec2  = dst_raw_bbox.max  - dst_raw_bbox.min;
    mi::math::Vector< mi::Sint64, 2 > const src_size_vec2  = src_raw_bbox.max  - src_raw_bbox.min;
    mi::math::Vector< mi::Sint64, 2 > const copy_size_vec2 = copy_src_bbox.max - copy_src_bbox.min;

    mi::math::Bbox< mi::Sint64, 2 >   const copy_dst_bbox(dst_origin, dst_origin + copy_size_vec2);
    if(!dst_raw_bbox.contains(copy_dst_bbox.max)){
        ERROR_LOG << "destination region is out of destination range. dst_raw_bbox = "
                  << dst_raw_bbox << ", copy_dst_bbox = " << copy_dst_bbox;
        return false;
    }

    mi::math::Vector< mi::Sint64, 2 > const offset_to_src_ij =
        dst_origin - copy_src_bbox.min;

    mi::Sint64 const dst_size  = dst_size_vec2[0]  * dst_size_vec2[1];
    mi::Sint64 const src_size  = src_size_vec2[0]  * src_size_vec2[1];
    mi::Sint64 const copy_size = copy_size_vec2[0] * copy_size_vec2[1];
    assert(dst_size  == dst_raw_bbox.volume());  nv::index_common::no_unused_variable_warning_please(dst_size);
    assert(src_size  == src_raw_bbox.volume());  nv::index_common::no_unused_variable_warning_please(src_size);
    assert(copy_size == copy_src_bbox.volume()); nv::index_common::no_unused_variable_warning_please(copy_size);

    // Note: no assumption on memory layout. Some optimization may be possible.
    mi::math::Vector< mi::Sint64, 2 > ij = copy_src_bbox.min; // for ij[1] initialization
    for(ij[0] = copy_src_bbox.min[0]; ij[0] < copy_src_bbox.max[0]; ++(ij[0]))
    {
        for(ij[1] = copy_src_bbox.min[1]; ij[1] < copy_src_bbox.max[1]; ++(ij[1]))
        {
            mi::math::Color const col = src_image_tile.get_color(ij);
            mi::math::Vector< mi::Sint64, 2 > const dst_ij = ij + offset_to_src_ij;
            this->put_color(dst_ij, col);
        }
    }

    return true;
}

//----------------------------------------------------------------------
bool Image_tile::save_buffer(std::string const & fname)
{
    assert(this->is_valid_buffer());
    mi::math::Vector< mi::Sint64, 2 > const buf_size = this->get_buffer_size();

    // convert buffer type
    std::vector< mi::math::Color_struct > tmp_color_st_buffer;
    size_t const buf_len = m_color_buffer.size();
    tmp_color_st_buffer.resize(buf_len);
    for(size_t i = 0; i < buf_len; ++i){
        mi::math::Color const & col = m_color_buffer.at(i);
        mi::math::Color_struct col_st = { col.r, col.g, col.b, col.a };
        tmp_color_st_buffer.at(i) = col_st;
    }

    bool const is_ok = nv::index_common::write_ppm(
        fname, tmp_color_st_buffer, static_cast<int>(buf_size[0]), static_cast<int>(buf_size[1]));
    return is_ok;
}

//----------------------------------------------------------------------
bool Image_tile::load_buffer(std::string const & fname)
{
    std::vector< mi::math::Color_struct > ppm_color_st_buf;
    mi::Sint32 img_width  = -1;
    mi::Sint32 img_height = -1;
    std::string error_mes;

    if(!nv::index_common::load_ppm(fname, ppm_color_st_buf, img_width, img_height, error_mes)){
        ERROR_LOG << "Image_tile::load_buffer: failed to load[" << fname << "].\n" << error_mes;
        return false;
    }
    mi::math::Bbox<mi::Sint64, 2 > const ij_bbox(0, 0, img_width, img_height);
    if(!this->resize_buffer(ij_bbox)){
        ERROR_LOG << "Image_tile::load_buffer: failed to resize to " << ij_bbox;
        return false;
    }

    size_t const buf_len = ppm_color_st_buf.size();
    for(size_t i = 0; i < buf_len; ++i){
        m_color_buffer.at(i) = ppm_color_st_buf.at(i);
    }

    return true;
}

//----------------------------------------------------------------------
void Image_tile::serialize(mi::neuraylib::ISerializer* serializer) const
{
    assert(this->is_valid_buffer());

    // serialize bounding box
    serializer->write(&(m_ij_bbox.min.x), 6);

    // serialize color buffer
    mi::Size const color_buffer_size = m_color_buffer.size();
    serializer->write(&color_buffer_size, 1);

    for(mi::Size i = 0; i < color_buffer_size; ++i){
        mi::math::Color const col = m_color_buffer.at(i);
        serializer->write(&(col.r), 4);
    }
}

//----------------------------------------------------------------------
void Image_tile::deserialize(mi::neuraylib::IDeserializer* deserializer)
{
    // deserialize bounding box
    mi::math::Bbox<mi::Sint64, 2 > ij_bbox;
    deserializer->read(&(ij_bbox.min.x), 6);
    this->resize_buffer(ij_bbox);

    // deserialize color buffer
    mi::Size color_buffer_size = 0;
    deserializer->read(&color_buffer_size, 1);
    mi::math::Vector< mi::Sint64, 2 > buf_size = this->get_buffer_size();
    assert(color_buffer_size == static_cast< mi::Size >(buf_size[0] * buf_size[1]));
    assert(color_buffer_size == m_color_buffer.size());
    nv::index_common::no_unused_variable_warning_please(buf_size);

    for(mi::Size i = 0; i < color_buffer_size; ++i){
        mi::math::Color col;
        deserializer->read(&(col.r), 4);
        m_color_buffer.at(i) = col;
    }
}

//----------------------------------------------------------------------

