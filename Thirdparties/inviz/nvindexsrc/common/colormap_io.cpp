/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "colormap_io.h"

#include <cstdio>
#include <cassert>
#include <deque>
#include <sstream>

#include "ppm_io.h"
#include "type_conversion_utility.h"

namespace nv {
namespace index_common {
//----------------------------------------------------------------------
std::string get_filename_extension(const std::string & fname)
{
    // get the extension of the filename
    size_t foundpos = fname.find_last_of(".", fname.size());
    if(foundpos == std::string::npos){
        return std::string(""); // extension not found
    }

    return fname.substr(foundpos, fname.size() - foundpos);
}

//----------------------------------------------------------------------
bool is_16bit_colormap_file(const std::string & colormap_fname)
{
    // get the extension of the filename
    const std::string ext = get_filename_extension(colormap_fname);
    if(ext.empty()){
        return false;           // no extension, not the 16bit colormap file
    }

    if(ext == ".cmap"){
        return true;
    }

    return false;
}

//----------------------------------------------------------------------
bool load_16bit_colormap(
    const std::string &                     colormap_fname,
    std::vector< mi::math::Color_struct > & colormap_entries,
    std::string & error_mes)
{
    FILE * pfp = fopen(colormap_fname.c_str(), "rb");
    if(pfp == 0)
    {
        error_mes = "Can not open the following 16-bit colormap file: \"" + colormap_fname + "\".";
        return false;
    }

    const mi::Float32 RECIPR = 1.0f/65535.0f; // 16bit color values are given as [0, 65535]
    const mi::Sint32 BUFSIZE = 1024;
    char line[BUFSIZE];
    mi::Sint32 line_num = 0;
    std::vector< mi::math::Color_struct > colormap_storage;
    bool is_ok = true;

    while(fgets(line, BUFSIZE, pfp) != 0){
        ++line_num;
        std::string line_str(line);
        if(line_str.size() == 0){
            // skip an empty line: just for sure for the last line.
            continue;
        }
        mi::Float32 col[4] = { 0.0f, 0.0f, 0.0f, 1.0f, };
        // some static code analyzer don't aware that we limit the
        // buffer size so this sscanf has no huge input data.
        const mi::Sint32 ret_elem = sscanf(line, "%20f %20f %20f %20f",
                                           &(col[0]), &(col[1]), &(col[2]), &(col[3]));
        if(ret_elem != 4){
            std::stringstream sstr;
            sstr << "Illegal line [" << colormap_fname << "]:" << line_num
                 << ", file collapsed? abort.";
            error_mes = sstr.str();
            is_ok = false;
            break;
        }
        const mi::math::Color_struct col_vec =
            get_color_struct(col[0] * RECIPR,
                             col[1] * RECIPR,
                             col[2] * RECIPR,
                             col[3] * RECIPR);
        colormap_storage.push_back(mi::math::clamp(col_vec, 0.0f, 1.0f));
    }

    if(is_ok){
        colormap_entries = colormap_storage;
    }
    fclose(pfp);

    return is_ok;
}

//----------------------------------------------------------------------
bool write_16bit_colormap(
    const std::string &                           colormap_fname,
    std::vector< mi::math::Color_struct > const & colormap_entries,
    std::string & error_mes)
{
    const size_t colormap_count = colormap_entries.size();
    if(colormap_count == 0){
        error_mes = "Empty colormap entry.";
        return false;
    }

    const std::string ext = get_filename_extension(colormap_fname);
    if(ext != ".cmap"){
        error_mes = "The file extension is not for 16bit colormap. "
            "[" + ext + "] instead of [.cmap]. But continue.";

    }

    FILE * pfp = fopen(colormap_fname.c_str(), "w");
    if(pfp == 0){
        error_mes = "Can not open a 16bit colormap [" + colormap_fname + "] (w)";
        return false;
    }

    const mi::Float32 FACTER = 65535.0f; // 16bit color values are given as [0, 65535]
    bool is_ok = true;

    for(size_t i = 0; i < colormap_count; ++i){
        mi::math::Color_struct col = colormap_entries.at(i);

        mi::Sint32 ret = fprintf(pfp, "%f %f %f %f\n",
                                 FACTER * col.r, FACTER * col.g, FACTER * col.b, FACTER * col.a);
        if(ret < 7){
            error_mes = "Can not write a file.";
            is_ok = false;
            break;
        }
    }

    fclose(pfp);

    return is_ok;
}

//----------------------------------------------------------------------
bool load_colormap(
    const std::string &                   colormap_fname,
    std::vector<mi::math::Color_struct> & colormap_entries,
    std::string & error_mes)
{
    const bool is_16bit_colormap = is_16bit_colormap_file(colormap_fname);

    bool ret = false;
    if(is_16bit_colormap){
        ret = load_16bit_colormap(colormap_fname, colormap_entries, error_mes);
    }
    else{
        mi::Sint32 img_w = 0;
        mi::Sint32 img_h = 0;
        ret = load_ppm(colormap_fname, colormap_entries, img_w, img_h, error_mes);
    }

    return ret;
}

//----------------------------------------------------------------------
}} // namespace nv::index_common svv
