/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "ppm_io.h"

#include <cstdio>
#include <cassert>
#include <sstream>

#include "forwarding_logger.h"

namespace nv {
namespace index_common {
//----------------------------------------------------------------------
/// get ppm color from char buffer
///
/// \param[in] p_color_buf color buffer pointer. The size must >= 3.
/// \return a color made by the first three component in the buffer.
static mi::math::Color_struct get_ppm_color_from_buffer(unsigned char * p_buf)
{
    assert(p_buf != 0);

    mi::Float32 col[3] = {0.0f, 0.0f, 0.0f};
    for(mi::Sint32 i = 0; i < 3; ++i){
        col[i] = static_cast< mi::Float32 >(p_buf[i]);
    }

    const mi::Float32 RECIPR = 1/255.0f; // ppm color values are in [0, 255]
    const mi::Float32 ALPHA  = 1.0f;
    mi::math::Color_struct nrcol = { col[0] * RECIPR, col[1] * RECIPR, col[2] * RECIPR, ALPHA };

    return nrcol;
}

//----------------------------------------------------------------------
bool load_ppm_sub_header(FILE * pfp,
                         std::string const & fname,
                         std::string & ppm_filemagic,
                         mi::Sint32 & img_width,
                         mi::Sint32 & img_height,
                         mi::Sint32 & img_depth,
                         std::string & error_mes)
{
    const mi::Sint32 BUFSIZE = 1024;
    char  line[BUFSIZE];

    // get ppm magic
    if(fgets(line, BUFSIZE, pfp) == 0){
        error_mes = "Error! Can not read a ppm file magic [" + fname + "] (r)";
        return false;
    }

    std::string filemagic(line);
    if(filemagic.size() != 3){
        error_mes = "Error! broken file magic infile [" + fname + "] (r)";
        return false;
    }
    filemagic = filemagic.substr(0, 2);

    // skip comment and empty lines
    while(fgets(line, BUFSIZE, pfp) != 0){
        if(line[0] != '#' && line[0] != '\n'){
            break;
        }
    }

    // get the size of image
    mi::Sint32 width  = -1;
    mi::Sint32 height = -1;
    {
        std::stringstream instr(line);
        instr >> width >> height;
        if(instr.fail()){
            error_mes = "Error! unrecognized ppm header. [" + std::string(line) + "]";
            return false;
        }
    }

    // get the depth
    if(fgets(line, BUFSIZE, pfp) == 0){
        error_mes = "Error! Can not read a depth of ppm file [" + fname + "] (r)";
        return false;
    }

    mi::Sint32 depth = -1;
    {
        std::stringstream instr(line);
        instr >> depth;
        if(instr.fail()){
            error_mes = "Error! unrecognized ppm header. [" + std::string(line) + "]";
            return false;
        }
    }

    if(!((width > 0) && (height > 0) && (depth > 0))){
        std::stringstream sstr;
        sstr << "Error! illegal image parameter. width, height, depth = [" << width
             << ", " << height << ", " << depth << "]";
        error_mes = sstr.str();
        return false;
    }

    // now we know the correct parameters, set them.
    ppm_filemagic = filemagic;
    img_width     = width;
    img_height    = height;
    img_depth     = depth;

    return true;
}

//----------------------------------------------------------------------
/// Load ppm file "P6" sub
///
/// \param[in]  pfp read opened file stream
/// \param[in]  fname    a reading filename for message out
/// \param[out] color_storage (output) colormap strage
/// \param[out] img_width     (output) image width
/// \param[out] img_height    (output) image height
/// \param[out] img_depth     (output) image depth
/// \param[out] error_mes     (output) error message when returns false
/// \return true when success
static bool load_ppm_p6_sub(FILE * pfp,
                            std::string const & fname,
                            std::vector< mi::math::Color_struct > & color_storage,
                            mi::Sint32 img_width,
                            mi::Sint32 img_height,
                            mi::Sint32 img_depth,
                            std::string & error_mes)
{
    assert(img_width  > 0);
    assert(img_height > 0);

    if(img_depth != 255){
        error_mes = "Error! Unexpected depth. The ppm file (P6) depth should be 255. ["
            + fname + "]";
        return false;
    }

    bool is_ok = true;
    size_t const BUFSIZE = img_width * img_height * 3; // RGB
    unsigned char * p_imgbuf = new unsigned char[BUFSIZE];
    {
        size_t const read_count = fread(p_imgbuf, 1, BUFSIZE, pfp);
        if(read_count == BUFSIZE){
            for(size_t i = 0; i < BUFSIZE; i = i + 3){
                mi::math::Color_struct const col = get_ppm_color_from_buffer(&(p_imgbuf[i]));
                color_storage.push_back(col);
            }
        }
        else{
            std::stringstream os;
            os << "Error! read failed, file truncated "
               << "(bytes read: " << read_count
               << ", expected bytes read: " << BUFSIZE
               << ")? [" << fname << "]";
            error_mes = os.str();
            is_ok = false;
        }
    }
    delete [] p_imgbuf;
    p_imgbuf = 0;

    return is_ok;
}

//----------------------------------------------------------------------
/// Load ppm file sub
///
/// \param[in]  pfp read opened file stream
/// \param[in]  fname    a reading filename for message out
/// \param[out] color_storage (output) colormap strage
/// \param[out] img_width     (output) image width
/// \param[out] img_height    (output) image height
/// \param[out] img_depth     (output) image depth
/// \param[out] error_mes     (output) error message when returns false
/// \return true when success
static bool load_ppm_sub(FILE * pfp,
                         std::string const & fname,
                         std::vector< mi::math::Color_struct > & color_storage,
                         mi::Sint32 & img_width,
                         mi::Sint32 & img_height,
                         mi::Sint32 & img_depth,
                         std::string & error_mes)
{
    assert(pfp != 0);

    std::string ppm_filemagic("");
    bool const is_header_ok = load_ppm_sub_header(pfp, fname,
                                                  ppm_filemagic,
                                                  img_width, img_height, img_depth,
                                                  error_mes);

    if(!is_header_ok){
        return false;
    }

    bool is_ppm_ok = false;
    if(ppm_filemagic == "P6"){
        is_ppm_ok = load_ppm_p6_sub(pfp, fname,
                                    color_storage, img_width, img_height,img_depth,
                                    error_mes);
    }
    else{
        error_mes = "Error! unsupported ppm file format [" + ppm_filemagic + "], P6 only.";
    }

    return is_ppm_ok;
}

//----------------------------------------------------------------------
bool load_ppm(
    const std::string &                     ppm_fname,
    std::vector< mi::math::Color_struct > & ppm_color,
    mi::Sint32 & img_width,
    mi::Sint32 & img_height,
    std::string & error_mes)
{
    //-- Open
    FILE * pfp = fopen(ppm_fname.c_str(), "rb");
    if(pfp == 0){
        error_mes = "Error! Can not open a ppm file [" + ppm_fname + "] (r)";
        return false;
    }

    std::vector< mi::math::Color_struct > color_storage;
    mi::Sint32 width  = -1;
    mi::Sint32 height = -1;
    mi::Sint32 depth  = -1;
    bool const ret = load_ppm_sub(pfp, ppm_fname, color_storage, width, height, depth,
                                  error_mes);

    //-- Close
    fclose(pfp);

    if(ret){
        // the values are right, return them, otherwise, keep as is (don't mess up).
        ppm_color  = color_storage;
        img_width  = width;
        img_height = height;
    }

    return ret;
}

//----------------------------------------------------------------------
bool get_ppm_color_buffer_from_pixel_buffer(
    mi::Uint32 width,
    mi::Uint32 height,
    mi::Uint8 const * const uchar_pixel_ptr,
    std::vector< mi::math::Color_struct > & ppm_color_buf)
{
    if((height <= 0) || (width <= 0) || (uchar_pixel_ptr == 0)){
        return false;
    }

    mi::Sint32 const swidth  = static_cast< mi::Sint32 >(width);
    mi::Sint32 const sheight = static_cast< mi::Sint32 >(height);

    for(mi::Sint32 h = sheight - 1; h >= 0; --h){
        for(mi::Sint32 w = 0; w < swidth; ++w){
            mi::Sint32 const bidx = (h * width + w) * 4;
            mi::Float32 const r = static_cast< mi::Float32 >(uchar_pixel_ptr[bidx    ]) / 255.0f;
            mi::Float32 const g = static_cast< mi::Float32 >(uchar_pixel_ptr[bidx + 1]) / 255.0f;
            mi::Float32 const b = static_cast< mi::Float32 >(uchar_pixel_ptr[bidx + 2]) / 255.0f;
            mi::Float32 const a = static_cast< mi::Float32 >(uchar_pixel_ptr[bidx + 3]) / 255.0f;

            mi::math::Color_struct col = {
                r, g, b, a,
            };
            ppm_color_buf.push_back(col);
        }
    }

    return true;
}

//----------------------------------------------------------------------
bool save_pixel_buffer_to_ppm_file(mi::Uint32 width,
                                   mi::Uint32 height,
                                   mi::Uint8 const * const uchar_pixel_ptr,
                                   std::string const & ofname)
{
    std::vector< mi::math::Color_struct > ppm_color_buf;
    bool ret = get_ppm_color_buffer_from_pixel_buffer(width, height, uchar_pixel_ptr, ppm_color_buf);
    if(!ret){
        return false;
    }

    ret = write_ppm(ofname, ppm_color_buf, width, height);

    return ret;
}

//----------------------------------------------------------------------
bool write_ppm(
    std::string const &                           ppm_fname,
    std::vector< mi::math::Color_struct > const & ppm_color,
    mi::Sint32 const img_width,
    mi::Sint32 const img_height)
{
    if((img_width <= 0)  ||
       (img_height <= 0) ||
       (ppm_color.size() != static_cast< size_t >(img_width * img_height)))
    {
        ERROR_LOG << "write_ppm(): illegal parameters for " << ppm_fname;
        return false;
    }

    //-- Open
    FILE * pfp = fopen(ppm_fname.c_str(), "wb");
    if(pfp == 0)
    {
        ERROR_LOG << "write_ppm(): Cannot open file " << ppm_fname;
        return false;
    }
    clearerr(pfp);

    // write P6 ppm

    // header
    fprintf(pfp, "P6\n");
    fprintf(pfp, "# CREATOR: nv::index::colormap_io 1.0\n");
    fprintf(pfp, "%d %d\n", img_width, img_height);
    fprintf(pfp, "255\n");

    mi::Float32 COEF = 255.0f;
    for(mi::Sint32 y = 0; y < img_height; ++y)
    {
        for(mi::Sint32 x = 0; x < img_width; ++x)
        {
            mi::Sint32 idx = (img_width * y) + x;
            mi::math::Color_struct const col = ppm_color.at(idx);
            unsigned char ucol_r = static_cast< unsigned char >(COEF * col.r);
            unsigned char ucol_g = static_cast< unsigned char >(COEF * col.g);
            unsigned char ucol_b = static_cast< unsigned char >(COEF * col.b);
            fprintf(pfp, "%c%c%c", ucol_r, ucol_g, ucol_b);
        }
    }

    mi::Sint32 const err = ferror(pfp);
    //-- Close
    fclose(pfp);

    if(err != 0)
    {
        ERROR_LOG << "write_ppm(): Write failed for " << ppm_fname;
        return false;
    }

    return true;
}

bool write_ppm(
    const std::string&                      ppm_fname,
    const std::vector<mi::Float32>&         heights,
    const mi::math::Vector<mi::Float32, 2>& height_range,
    const mi::Sint32                        img_width,
    const mi::Sint32                        img_height)
{
    if((img_width <= 0)  ||
       (img_height <= 0) ||
       (heights.size() != static_cast< size_t >(img_width * img_height)))
    {
        ERROR_LOG << "write_ppm(): illegal parameters for " << ppm_fname;
        return false;
    }

    //-- Open
    FILE * pfp = fopen(ppm_fname.c_str(), "wb");
    if(pfp == 0)
    {
        ERROR_LOG << "write_ppm(): Cannot open file " << ppm_fname;
        return false;
    }
    clearerr(pfp);

    // write P6 ppm

    // header
    fprintf(pfp, "P6\n");
    fprintf(pfp, "# CREATOR: nv::index::colormap_io 1.0\n");
    fprintf(pfp, "%d %d\n", img_width, img_height);
    fprintf(pfp, "255\n");

    mi::Float32 COEF = 255.0f;

    mi::Float32 denomiator = 1.f / (height_range.y - height_range.x);
    for(mi::Sint32 y = 0; y < img_height; ++y)
    {
        for(mi::Sint32 x = 0; x < img_width; ++x)
        {
            mi::Sint32 idx = (img_width * y) + x;
            mi::Float32 col = (heights.at(idx)-height_range.x)*denomiator;
            unsigned char ucol_r = static_cast< unsigned char >(COEF * col);
            unsigned char ucol_g = static_cast< unsigned char >(COEF * col);
            unsigned char ucol_b = static_cast< unsigned char >(COEF * col);
            fprintf(pfp, "%c%c%c", ucol_r, ucol_g, ucol_b);
        }
    }

    mi::Sint32 const err = ferror(pfp);
    //-- Close
    fclose(pfp);

    if(err != 0)
    {
        ERROR_LOG << "write_ppm(): Write failed for " << ppm_fname;
        return false;
    }

    return true;
}

//----------------------------------------------------------------------
//======================================================================
//----------------------------------------------------------------------
PPM_image::PPM_image()
    :
    m_width(-1),
    m_height(-1)
{
    // empty
}

//----------------------------------------------------------------------
PPM_image::~PPM_image()
{
    m_color_buffer.clear();
    m_width  = -1;
    m_height = -1;
}

//----------------------------------------------------------------------
void PPM_image::resize_buffer(mi::Sint32 width, mi::Sint32 height)
{
    assert(width  > 0);
    assert(height > 0);
    m_width  = width;
    m_height = height;

    m_color_buffer.resize(m_width * m_height);
}

//----------------------------------------------------------------------
mi::Sint32 PPM_image::get_width() const
{
    return m_width;
}

//----------------------------------------------------------------------
mi::Sint32 PPM_image::get_height() const
{
    return m_height;
}

//----------------------------------------------------------------------
bool PPM_image::is_valid_buffer() const
{
    if((m_width > 0) && (m_height > 0)){
        assert((m_width * m_height) == static_cast< mi::Sint32 >(m_color_buffer.size()));
        return true;
    }
    return false;
}

//----------------------------------------------------------------------
void PPM_image::clear_buffer(mi::math::Color_struct const & col)
{
    assert(this->is_valid_buffer());
    size_t buf_size = m_color_buffer.size();
    for(size_t i = 0; i < buf_size; ++i){
        m_color_buffer.at(i) = col;
    }
}

//----------------------------------------------------------------------
void PPM_image::put_color(mi::Sint32 x, mi::Sint32 y, mi::math::Color_struct const & col)
{
    m_color_buffer.at(this->get_idx(x, y)) = col;
}

//----------------------------------------------------------------------
mi::math::Color_struct PPM_image::get_color(mi::Sint32 x, mi::Sint32 y)
{
    return m_color_buffer.at(this->get_idx(x, y));
}

//----------------------------------------------------------------------
bool PPM_image::save_buffer(std::string const & fname)
{
    assert(this->is_valid_buffer());
    bool const is_ok = write_ppm(fname, m_color_buffer, m_width, m_height);
    return is_ok;
}

//----------------------------------------------------------------------
bool PPM_image::load_buffer(std::string const & fname)
{
    std::vector< mi::math::Color_struct > ppm_color;
    mi::Sint32 img_width  = -1;
    mi::Sint32 img_height = -1;
    std::string error_mes;
    bool const is_success = load_ppm(fname, ppm_color, img_width, img_height,
                                     error_mes);
    if(!is_success){
        ERROR_LOG << "Fail to PPM_image:::load_buffer [" << fname << "]\n" << error_mes;
        return false;
    }

    assert(img_width  > 0);
    assert(img_height > 0);

    this->resize_buffer(img_width, img_height);
    assert(ppm_color.size() == m_color_buffer.size());

    m_color_buffer = ppm_color;

    return true;
}

//----------------------------------------------------------------------
bool PPM_image::set_buffer(mi::Sint32 width, mi::Sint32 height,
                            std::vector< mi::math::Color_struct > const & color_buffer)
{
    if((width <= 0) || (height <= 0)){
        return false;
    }
    this->resize_buffer(width, height);

    if((width * height) != static_cast< mi::Sint32 >(color_buffer.size())){
        return false;
    }

    m_color_buffer = color_buffer;
    return true;
}

//----------------------------------------------------------------------
/// L1 difference of RGB
///
/// \param[in] c0 color0 to compare
/// \param[in] c1 color1 to compare
/// \return L1 difference
inline mi::Float32 diff_rgb_L1(mi::math::Color const & c0,
                               mi::math::Color const & c1)
{
    mi::Float32 const rgb_L1 = fabsf(c0[0] - c1[0]) + fabsf(c0[1] - c1[1]) + fabsf(c0[2] - c1[2]);
    return rgb_L1;
}

//----------------------------------------------------------------------
bool PPM_image::is_equal(PPM_image const & other_image, mi::Float32 threshold)
{
    if(this->get_width() != other_image.get_width()){
        DEBUG_LOG << "width differs.";
        return false;
    }
    if(this->get_height() != other_image.get_height()){
        DEBUG_LOG << "height differs.";
        return false;
    }
    if(!(this->is_valid_buffer())){
        DEBUG_LOG << "this is not a valid PPM_image.";
        return false;
    }
    if(!(other_image.is_valid_buffer())){
        DEBUG_LOG << "other is not a valid PPM_image.";
        return false;
    }

    mi::Sint32 differpix_count = 0;
    mi::Sint32 abs_differpix_count = 0;
    size_t const bufsz = m_color_buffer.size();
    for(size_t i = 0; i < bufsz; ++i){
        mi::math::Color c0 = this->m_color_buffer[i];
        mi::math::Color c1 = other_image.m_color_buffer[i];
        // Don't compare alpha since ppm doesn't have it.
        // truncate them since ppm has only integer value.
        // clamp since ppm file doesn't handle float
        c0 *= 255.0f;
        c0.r = static_cast< mi::Float32 >(static_cast< unsigned char >(c0.r));
        c0.g = static_cast< mi::Float32 >(static_cast< unsigned char >(c0.g));
        c0.b = static_cast< mi::Float32 >(static_cast< unsigned char >(c0.b));

        c1 *= 255.0f;
        c1.r = static_cast< mi::Float32 >(static_cast< unsigned char >(c1.r));
        c1.g = static_cast< mi::Float32 >(static_cast< unsigned char >(c1.g));
        c1.b = static_cast< mi::Float32 >(static_cast< unsigned char >(c1.b));

        mi::Float32 const diff_L1 = diff_rgb_L1(c0, c1);
        if(diff_L1 != 0.0f){
            ++abs_differpix_count;
        }
        if(diff_L1 > 3.0f){
            ++differpix_count;
        }
        // mi::Sint32 const x = i / this->get_width();
        // mi::Sint32 const y = i % this->get_height();
        //  std::cout << "at [" << x << ", " << y << ", " << j << "]: " << c0_j
        //            << " != " << c1_j << std::endl;
    }

    if(abs_differpix_count > 0){
        std::stringstream sstr;
        sstr << "Absolute different pixel count: " << abs_differpix_count << "/" << bufsz;
        INFO_LOG << sstr.str();
    }

    if(differpix_count > 0){
        std::stringstream sstr;
        sstr << "Different pixels: " << differpix_count << "/" << bufsz;
        mi::Float32 percentage
            = (static_cast<mi::Float32>(differpix_count) / static_cast<mi::Float32>(bufsz)) * 100.f;
        if(percentage > threshold){
            sstr << ", " << percentage << "%, more than threshold " << threshold << "%.";
            ERROR_LOG << sstr.str();
            return false;
        }
        INFO_LOG << sstr.str();
    }
    return true;
}

//----------------------------------------------------------------------
}} // namespace nv::index_common svv
