/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief ppm file IO routine

#ifndef NVIDIA_INDEX_BIN_COMMON_PPM_IO_H
#define NVIDIA_INDEX_BIN_COMMON_PPM_IO_H

#include <vector>
#include <string>
#include <cassert>

#include <mi/math/color.h>

namespace nv {
namespace index_common {
//----------------------------------------------------------------------
/// Load ppm file (P6 file)
///
/// \param[in]  ppm_fname   ppm filename
/// \param[out] ppm_color   (output) color contents
/// \param[out] img_width   (output) image width  size
/// \param[out] img_height  (output) image height size
/// \param[out] error_mes   (output) error message when returns false
/// \return true when succeeded
bool load_ppm(
    const std::string &                     ppm_fname,
    std::vector< mi::math::Color_struct > & ppm_color,
    mi::Sint32 & img_width,
    mi::Sint32 & img_height,
    std::string & error_mes);


//----------------------------------------------------------------------
/// Load ppm file sub: loading header of ppm file
///
/// \param[in]  pfp read opened file stream
/// \param[in]  fname    a reading filename for message out
/// \param[out] ppm_filemagic (output) file magic of ppm
/// \param[out] img_width     (output) image width
/// \param[out] img_height    (output) image height
/// \param[out] img_depth     (output) image depth
/// \param[out] error_mes     (output) error message when returns false
/// \return true when success
bool load_ppm_sub_header(FILE * pfp,
                         std::string const & fname,
                         std::string & ppm_filemagic,
                         mi::Sint32 & img_width,
                         mi::Sint32 & img_height,
                         mi::Sint32 & img_depth,
                         std::string & error_mes);

//----------------------------------------------------------------------
/// get a ppm color buffer from a raw pixel buffer
///
/// \param[in] width  canvas width
/// \param[in] height canvas height
/// \param[in] uchar_pixel_ptr pixel buffer's first pointer (assumed RGBA)
/// \param[out] ppm_color_buf  (output) ppm color buffer
/// \return truen when succeeded
bool get_ppm_color_buffer_from_pixel_buffer(
    mi::Uint32 width,
    mi::Uint32 height,
    mi::Uint8 const * const uchar_pixel_ptr,
    std::vector< mi::math::Color_struct > & ppm_color_buf);

//----------------------------------------------------------------------
/// write a snapshot image to a ppm file.
///
/// \param[in] width  canvas width
/// \param[in] height canvas height
/// \param[in] uchar_pixel_ptr pixel buffer's first pointer (assumed RGBA)
/// \param[in] ofname output filename
/// \return truen when succeeded
bool save_pixel_buffer_to_ppm_file(
    mi::Uint32 width,
    mi::Uint32 height,
    mi::Uint8 const *  const uchar_pixel_ptr,
    std::string const & ofname);

//----------------------------------------------------------------------
/// write ppm file (P6 file)
///
/// \param[in] ppm_fname   ppm filename
/// \param[in] ppm_color   color contents. The range of Color_struct is [0,1]
/// \param[in] img_width   image width  size
/// \param[in] img_height  image height size
/// \return true when succeeded
bool write_ppm(
    std::string const &                           ppm_fname,
    std::vector< mi::math::Color_struct > const & ppm_color,
    mi::Sint32 const img_width,
    mi::Sint32 const img_height);

//----------------------------------------------------------------------
/// write ppm file (P6 file)
///
/// \param[in] ppm_fname        ppm filename
/// \param[in] heights          height values typically used in an elevation map/height field
/// \param[in] height_range     the range of height values [min/max]
/// \param[in] img_width        image width  size
/// \param[in] img_height       image height size
/// \return true when succeeded
bool write_ppm(
    const std::string&                      ppm_fname,
    const std::vector<mi::Float32>&         heights,
    const mi::math::Vector<mi::Float32, 2>& height_range,
    const mi::Sint32                        img_width,
    const mi::Sint32                        img_height);

//----------------------------------------------------------------------
/// color image helper class
class PPM_image
{
public:
    /// default constructor
    PPM_image();
    /// destructor
    virtual ~PPM_image();

    /// resize buffer
    void resize_buffer(mi::Sint32 width, mi::Sint32 height);
    /// get width
    /// \return width of the buffer
    mi::Sint32 get_width() const;
    /// get height
    /// \return height of the buffer
    mi::Sint32 get_height() const;
    /// is valid
    bool is_valid_buffer() const;
    /// clear buffer
    /// \param[in] col color to fill
    void clear_buffer(mi::math::Color_struct const & col);
    /// put color
    void put_color(mi::Sint32 x, mi::Sint32 y, mi::math::Color_struct const & col);
    /// get color
    mi::math::Color_struct get_color(mi::Sint32 x, mi::Sint32 y);
    /// save buffer
    /// \param[in] fname filename to save
    bool save_buffer(std::string const & fname);
    /// load buffer
    /// \param[in] fname filename to load
    bool load_buffer(std::string const & fname);
    /// set buffer by copy
    /// \param[in] width  image width
    /// \param[in] height image height
    /// \param[in] color_buffer color buffer
    /// \return true when succeeded
    bool set_buffer(mi::Sint32 width, mi::Sint32 height,
                     std::vector< mi::math::Color_struct > const & color_buffer);
    /// is equal contents?
    /// \param[in] other_image other PPM_image.
    /// \param[in] threshold   consider not equal if more than this percentag of pixels differ
    /// \return true when the buffer contents are the same.
    bool is_equal(PPM_image const & other_image, mi::Float32 threshold);

private:
    /// get index
    inline mi::Sint32 get_idx(mi::Sint32 x, mi::Sint32 y){
        assert(this->is_valid_buffer());
        assert((0 <= x) && (x < m_width));
        assert((0 <= y) && (y < m_height));
        return (m_width * y) + x;
    }

private:
    /// pixels
    std::vector< mi::math::Color_struct > m_color_buffer;
    /// image width
    mi::Sint32 m_width;
    /// image height
    mi::Sint32 m_height;

private:
    /// copy constructor. prohibit until proved useful.
    PPM_image(PPM_image const &);
    /// operator=. prohibit until proved useful.
    PPM_image const & operator=(PPM_image const &);
};

//----------------------------------------------------------------------
}} // namespace nv::index_common
#endif // NVIDIA_INDEX_BIN_COMMON_PPM_IO_H
