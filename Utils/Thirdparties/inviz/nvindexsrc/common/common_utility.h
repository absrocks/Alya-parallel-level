/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief common utility

#ifndef NVIDIA_INDEX_BIN_COMMON_COMMON_UTILITY_H
#define NVIDIA_INDEX_BIN_COMMON_COMMON_UTILITY_H

#include <mi/dice.h>

#include <cassert>
#include <fstream>
#include <sstream>
#include <string>

#include "forwarding_logger.h"
#include "string_dict.h"
#include "tokenizer.h"

#ifdef LINUX
#include <unistd.h>
#endif

namespace nv {
namespace index_common {
//----------------------------------------------------------------------
/// High resolution sleep
/// \param[in] seconds sleep seconds (e.g., 0.01 sec)
inline void sleep(mi::Float32 seconds)
{
#ifdef _WIN32
    // The windows version expects milliseconds here.
    mi::Uint32 millis = (mi::Uint32)(seconds * 1000);
    if (millis == 0)
        millis = 1;
    ::Sleep(millis);
#else
    useconds_t micros = (useconds_t)(seconds * 1000000);
    ::usleep((useconds_t) micros);
#endif  // _WIN32
}

//----------------------------------------------------------------------
/// get current time
/// \return get Float64 epoch time in second
inline mi::Float64 get_time()
{
#ifdef _WIN32
    static bool init = false;
    static mi::Float64 frequency;
    if(!init)
    {
        //_tzset();
        LARGE_INTEGER freq;
        QueryPerformanceFrequency(&freq);
        frequency = (mi::Float64)freq.QuadPart;
        init      = true;
    }
    LARGE_INTEGER counter;
    QueryPerformanceCounter(&counter);

    return (mi::Float64)counter.QuadPart / frequency;
#else
    timeval tv;
    gettimeofday(&tv, NULL);
    return static_cast<mi::Float64>(tv.tv_sec) + (static_cast<mi::Float64>(tv.tv_usec) * 1.0e-6);
#endif
}

//----------------------------------------------------------------------
/// retrieves the host name of this machine
/// \return hostname of this machine.
inline std::string get_host_name()
{
    std::string host_name = "unknown";
#ifdef LINUX
    char buf[256];
    if (gethostname(buf, sizeof(buf)) == 0)
        host_name = buf;
#else
    char* host_name_env = getenv("HOSTNAME");
    if (host_name_env == NULL)
        host_name_env = getenv("HOST");

    if (host_name_env != NULL)
        host_name = host_name_env;
    //else
    //    INFO_LOG << "Environment variable 'HOSTNAME' or 'HOST' not set on host.";
#endif // LINUX

    return host_name;
}

//----------------------------------------------------------------------
/// get current working directory as a string
/// \return current working directory
inline std::string get_current_directory_name_string()
{
    std::string cwd_name = "no current directory name available.";
#ifdef LINUX
    char *newed_buffer = get_current_dir_name();
    assert(newed_buffer != 0);
    cwd_name = newed_buffer;

#elif _WIN32
    DWORD length = GetCurrentDirectory(0, 0);
    cwd_name = "";
    if (length > 0)
    {
        mi::Size length_size = static_cast<mi::Size>(length);
        char *newed_buffer = new char[length_size];
        assert(newed_buffer != 0);
        for (mi::Size i = 0; i < length_size; ++i)
        {
            newed_buffer[i] = '\0';
        }
        if (GetCurrentDirectory(length, newed_buffer) > 0)
        {
            cwd_name = newed_buffer;
            assert(cwd_name.length() == length_size);
        }
    }

#else
    // Nothing to do for other platform
#endif

    return cwd_name;
}

//----------------------------------------------------------------------
/// get volume index
/// \param[in] ijk             ijk coordinate of volume data access
/// \param[in] volume_raw_bbox volume data bounding box.
/// \return The index of volume data value
inline mi::Sint64 get_volume_index(
    const mi::math::Vector< mi::Sint64, 3 >& ijk,
    const mi::math::Bbox<   mi::Sint64, 3 >& volume_raw_bbox)
{
#ifdef DEBUG
    for(mi::Sint32 i = 0; i < 3; ++i)
    {
        if((!(volume_raw_bbox.min[i] <= ijk[i])) ||
           (!(volume_raw_bbox.max[i] >  ijk[i])))
        {
            ERROR_LOG << "get_volume_index: volume_raw_bbox = " << volume_raw_bbox
                      << ", ijk = " << ijk;
        }
        assert(volume_raw_bbox.min[i] <= ijk[i]);
        assert(volume_raw_bbox.max[i] >  ijk[i]);
    }
#endif  // DEBUG
    const mi::math::Vector< mi::Sint64, 3 > vol_size_vec3 = volume_raw_bbox.max - volume_raw_bbox.min;
    const mi::math::Vector< mi::Sint64, 3 > offset_idx = ijk - volume_raw_bbox.min;

    // Horner
    const mi::Sint64 index =
        vol_size_vec3[2] * ((offset_idx[0] * vol_size_vec3[1]) + offset_idx[1]) + offset_idx[2];

    assert(index >= 0);
    assert(index < (vol_size_vec3[0] * vol_size_vec3[1] * vol_size_vec3[2]));

    return index;
}

//----------------------------------------------------------------------
/// get heightfield index
/// \param[in] ij                   ij coordinate of heightfield patch
/// \param[in] heightfield_raw_bbox heightfield patch bounding box.
/// \return height value index of the heightfield
inline mi::Sint64 get_heightfield_index(
    const mi::math::Vector< mi::Sint64, 2 >& ij,
    const mi::math::Bbox<   mi::Sint64, 2 >& heightfield_raw_bbox)
{
#ifdef DEBUG
    for(mi::Sint32 i = 0; i < 2; ++i)
    {
        if((!(heightfield_raw_bbox.min[i] <= ij[i])) ||
           (!(heightfield_raw_bbox.max[i] >  ij[i])))
        {
            ERROR_LOG << "get_heightfield_index: heightfield_raw_bbox = " << heightfield_raw_bbox
                      << ", ij = " << ij;
        }
        assert(heightfield_raw_bbox.min[i] <= ij[i]);
        assert(heightfield_raw_bbox.max[i] >  ij[i]);
    }
#endif  // DEBUG
    const mi::math::Vector< mi::Sint64, 2 > heightfield_size_vec2 = heightfield_raw_bbox.max - heightfield_raw_bbox.min;
    const mi::math::Vector< mi::Sint64, 2 > offset_idx = ij - heightfield_raw_bbox.min;

    // Horner
    // mi::Sint64 const index = (offset_idx[0] * heightfield_size_vec2[1]) + offset_idx[1];  // j direction continuous.
    const mi::Sint64 index = (offset_idx[1] * heightfield_size_vec2[0]) + offset_idx[0]; // i direction continuous.

    assert(index >= 0);
    assert(index < (heightfield_size_vec2[0] * heightfield_size_vec2[1]));

    return index;
}

//----------------------------------------------------------------------
/// volume containing test
///
/// \param[in] bbox0 volume bounding box 0
/// \param[in] bbox1 volume bounding box 1
/// \return true when volume 0 contains volume 1.
inline bool bbox_contains(mi::math::Bbox< mi::Sint64, 3 > const & bbox0,
                          mi::math::Bbox< mi::Sint64, 3 > const & bbox1)
{
    if(bbox0.contains(bbox1.min) && bbox0.contains(bbox1.max)){
        return true;
    }
    return false;
}

//----------------------------------------------------------------------
/// patch containing test
///
/// \param[in] patchbbox0 patch bounding box 0
/// \param[in] patchbbox1 patch bounding box 1
/// \return true when patch 0 contains patch 1.
inline bool bbox_contains(mi::math::Bbox< mi::Sint64, 2 > const & patchbbox0,
                          mi::math::Bbox< mi::Sint64, 2 > const & patchbbox1)
{
    if(patchbbox0.contains(patchbbox1.min) && patchbbox0.contains(patchbbox1.max)){
        return true;
    }
    return false;
}

//----------------------------------------------------------------------
/// Access a given variable in order to avoid "unused variable" warnings.
template <class T>
inline void no_unused_variable_warning_please(T & /* ignored */)
{
    // empty
}

//----------------------------------------------------------------------
/// get string mi::Sint32 list to vector< mi::Sint32 >.
inline std::vector< mi::Sint32 > get_sint32_vector(std::string const & int_str_list)
{
    std::vector< std::string > tokens;
    nv::index_common::Tokenizer::parse(int_str_list, " ", tokens);

    std::vector< mi::Sint32 > int_vec;
    for(std::vector< std::string >::const_iterator ti = tokens.begin(); ti != tokens.end(); ++ti){
        if((*ti) == ""){
            // skip empty string (== str.strip())
            continue;
        }
        int_vec.push_back(get_sint32(*ti));
    }
    return int_vec;
}


//----------------------------------------------------------------------
/// get an bbox-and result as an bbox, Sint64
inline mi::math::Bbox< mi::Sint64, 3 > bbox_and(mi::math::Bbox< mi::Sint64, 3 > const & bbox0,
                                                mi::math::Bbox< mi::Sint64, 3 > const & bbox1)
{
    mi::math::Bbox< mi::Sint64, 3 > res_bbox;

    if(!bbox0.intersects(bbox1)){
        // no and region
        // DEBUG_LOG << "bbox_and: no and region.";
        return mi::math::Bbox< mi::Sint64, 3 >();
    }

    for(mi::Sint32 i = 0; i < 3; ++i){
        res_bbox.min[i] = std::max(bbox0.min[i], bbox1.min[i]);
        res_bbox.max[i] = std::min(bbox0.max[i], bbox1.max[i]);
    }

    return res_bbox;
}

//----------------------------------------------------------------------
/// get an bbox-and result as an bbox, Float32
inline mi::math::Bbox< mi::Float32, 3 > bbox_and(mi::math::Bbox< mi::Float32, 3 > const & bbox0,
                                                 mi::math::Bbox< mi::Float32, 3 > const & bbox1)
{
    mi::math::Bbox< mi::Float32, 3 > res_bbox;

    if(!bbox0.intersects(bbox1)){
        // no and region
        // DEBUG_LOG << "bbox_and: no and region.";
        return mi::math::Bbox< mi::Float32, 3 >();
    }

    for(mi::Sint32 i = 0; i < 3; ++i){
        res_bbox.min[i] = std::max(bbox0.min[i], bbox1.min[i]);
        res_bbox.max[i] = std::min(bbox0.max[i], bbox1.max[i]);
    }

    return res_bbox;
}


//----------------------------------------------------------------------
/// get a map from a inline parameter string
///
/// Inline parameter string is one line string with comma separated key=value.
/// E.g., "key=value,key=value,..." 
///
/// \param[in]  inline_pstring inline parameter string
/// \param[out] out_map output std::map
/// \param[out] err_mes error message if returns false
/// \return true when file parsing succeeded
inline bool get_map_from_inline_parameter_string(
    std::string const & inline_pstring,
    std::map< std::string, std::string > & out_map,
    std::string & err_mes
    )
{
    std::string fn = "get_map_from_inline_parameter_string: ";

    std::vector< std::string > key_value_vec;
    std::string const separator = ",";

    nv::index_common::Tokenizer::parse(inline_pstring, separator, key_value_vec);
    std::map< std::string, std::string > out_map_tmp;

    // get each key_value
    size_t len = key_value_vec.size();
    for(size_t i = 0; i < len; ++i){
        std::vector< std::string > key_value;
        nv::index_common::Tokenizer::parse(key_value_vec[i], "=", key_value);
        if(key_value.size() != 2){
            // no key=value (or missing one of them)
            std::stringstream sstr;
            sstr << "Error! " << fn << "[" << i << "]-th element [" 
                 << key_value_vec[i] << "] is not a 'key=value'.";
            err_mes = sstr.str();
            return false;
        }
        if(out_map_tmp.find(key_value[0]) != out_map_tmp.end()){
            // duplicated key has no sense for this. Cast an error.
            std::stringstream sstr;
            sstr << "Error! " << fn << "[" << i << "]-th element [" 
                 << key_value_vec[i] << "]'s key is duplicated.";
            err_mes = sstr.str();
            return false;
        }
        out_map_tmp[key_value[0]] = key_value[1];
    }

    out_map = out_map_tmp;

    return true;
}

inline
bool read_file(const std::string& file_path, std::string& out_file_string)
{
    std::ifstream   src_file;

    src_file.open(file_path.c_str(), std::ios::in | std::ios::binary | std::ios::ate);

    if (!src_file.is_open()) {
        return false;
    }
    const std::ifstream::pos_type src_file_size = src_file.tellg();
    src_file.seekg(std::ios::beg);

    out_file_string.resize(src_file_size, '\0');
    src_file.read(&out_file_string[0], src_file_size);

    if (!src_file) {
        out_file_string.clear();
        return false;
    }

    src_file.close();

    return true;
}

//----------------------------------------------------------------------
}} // namespace nv::index_common
#endif // NVIDIA_INDEX_BIN_COMMON_COMMON_UTILITY_H
