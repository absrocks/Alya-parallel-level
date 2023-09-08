/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "utilities.h"
#include <limits>
#include <mi/math/vector.h>

#include "common/string_dict.h"
#include "common/common_utility.h"

#include "ring_buffer.h"

#include <iomanip>
#include <istream>
#include <ostream>
#include <sstream>

#ifdef _WIN32
#    ifndef WIN32_LEAN_AND_MEAN
#      define WIN32_LEAN_AND_MEAN 1
#    endif
#  include <windows.h>
#endif

#ifdef LINUX
#include <unistd.h>
#endif

//----------------------------------------------------------------------
std::string concat_string_with_separator(
    std::vector< std::string > const & str_vec,
    std::string const & sep_str)
{
    size_t const nb_heightfield = str_vec.size();
    std::ostringstream os;
    for(size_t i = 0; i < nb_heightfield; ++i){
        os << str_vec.at(i);
        if(i < (nb_heightfield - 1)){
            os << sep_str;
        }
    }
    return os.str();
}

//----------------------------------------------------------------------
void get_region_vector(const mi::math::Vector<mi::Sint32, 2> & tatal_range_bbox,
                       mi::Sint32 interval_size,
                       std::vector< mi::math::Vector<mi::Sint32, 2> > & region_vec)
{
    assert(interval_size > 0);

    region_vec.clear();

    // get n-slices from a patial volume
    const mi::Sint32 start_k = tatal_range_bbox.x;
    const mi::Sint32 max_k   = tatal_range_bbox.y;
    const mi::Sint32 total_slices     = max_k - start_k;
    const mi::Sint32 n_slice_reminder = total_slices % interval_size;
    const mi::Sint32 k_n_slices       = total_slices / interval_size;
    assert((k_n_slices * interval_size + n_slice_reminder) == total_slices);

    for(mi::Sint32 i = 0; i < k_n_slices; ++i){
        mi::math::Vector<mi::Sint32, 2> range_bbox;
        range_bbox.x = start_k +  i      * interval_size;
        range_bbox.y = start_k + (i + 1) * interval_size;
        region_vec.push_back(range_bbox);
    }

    if(n_slice_reminder > 0){
        const mi::Sint32 reminder_start = start_k + (k_n_slices * interval_size);

        mi::math::Vector<mi::Sint32, 2> range_bbox;
        range_bbox.x = reminder_start;
        range_bbox.y = max_k;
        region_vec.push_back(range_bbox);
    }
}

//----------------------------------------------------------------------
mi::Sint32 meshgrid4_get_idx(mi::Sint32 ix, mi::Sint32 iy, mi::Sint32 iz)
{
    return ix + ((iy + (iz * 4)) * 4);
}

//----------------------------------------------------------------------
void meshgrid4(mi::math::Vector<mi::Sint64, 4> const & xrange,
               mi::math::Vector<mi::Sint64, 4> const & yrange,
               mi::math::Vector<mi::Sint64, 4> const & zrange,
               mi::math::Vector<mi::Sint64, 3> * meshgrid4_result)
{
    // check the inputs. (only for computing grid, not general case.)
    for(mi::Sint32 i = 0; i < 3; ++i){
        if((!(xrange[i] < xrange[i + 1])) ||
           (!(yrange[i] < yrange[i + 1])) ||
           (!(zrange[i] < zrange[i + 1]))){
            ERROR_LOG << "The range is not ascendant.";
            assert(false);
        }
    }
    for(mi::Sint32 ix = 0; ix < 4; ++ix){
        for(mi::Sint32 iy = 0; iy < 4; ++iy){
            for(mi::Sint32 iz = 0; iz < 4; ++iz){
                mi::Sint32 const idx = meshgrid4_get_idx(ix, iy, iz);
                meshgrid4_result[idx].x = xrange[ix];
                meshgrid4_result[idx].y = yrange[iy];
                meshgrid4_result[idx].z = zrange[iz];
            }
        }
    }
}

//----------------------------------------------------------------------
bool check_necessary_key(nv::index_common::String_dict const & dict,
                         char const * const p_key[])
{
    std::vector< std::string > undef_list;
    bool const is_all_key = nv::index_common::is_all_keys_defined(dict, p_key, &undef_list);
    if(!is_all_key){
        std::stringstream sstr;
        sstr << "following undefined keys in the project file.";
        for(std::vector< std::string >::const_iterator si = undef_list.begin();
            si != undef_list.end(); ++si){
            sstr << "\n[" << *si << "]";
        }
        ERROR_LOG << sstr.str();
        return false;
    }
    return true;
}

//----------------------------------------------------------------------
mi::math::Matrix_struct<mi::Float32, 4, 4>
get_translation_matrix(mi::math::Vector_struct<mi::Float32, 3> const & translation_vec)
{
    mi::math::Matrix<mi::Float32, 4, 4> translation_mat(
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        translation_vec.x, translation_vec.y, translation_vec.z, 1);
    return translation_mat;

}

//----------------------------------------------------------------------
mi::math::Matrix_struct<mi::Float32, 4, 4>
get_scaling_matrix(mi::math::Vector_struct<mi::Float32, 3> const & scaling_vec)
{
    mi::math::Matrix<mi::Float32, 4, 4> scaling_mat(
        scaling_vec.x, 0, 0, 0,
        0, scaling_vec.y, 0, 0,
        0, 0, scaling_vec.z, 0,
        0, 0, 0, 1);
    return scaling_mat;
}

//----------------------------------------------------------------------
mi::math::Matrix_struct<mi::Float32, 4, 4>
get_rotation_matrix_about_x(mi::Float32 angle)
{
    mi::Float32 sinp, cosp;
    mi::math::sincos(angle, sinp, cosp);

    mi::math::Matrix<mi::Float32, 4, 4> rotation_mat(
        1, 0, 0, 0,
        0, cosp, sinp, 0,
        0, -sinp, cosp, 0,
        0, 0, 0, 1);
    return rotation_mat;
}

//----------------------------------------------------------------------
mi::math::Matrix_struct<mi::Float32, 4, 4>
get_rotation_matrix_about_y(mi::Float32 angle)
{
    mi::Float32 sinp, cosp;
    mi::math::sincos(angle, sinp, cosp);

    mi::math::Matrix<mi::Float32, 4, 4> rotation_mat(
        cosp, 0, -sinp, 0,
        0, 1, 0, 0,
        sinp, 0, cosp, 0,
        0, 0, 0, 1);
    return rotation_mat;
}

//----------------------------------------------------------------------
mi::math::Matrix_struct<mi::Float32, 4, 4>
get_rotation_matrix_about_z(mi::Float32 angle)
{
    mi::Float32 sinp, cosp;
    mi::math::sincos(angle, sinp, cosp);

    mi::math::Matrix<mi::Float32, 4, 4> rotation_mat(
        cosp, sinp, 0, 0,
        -sinp, cosp, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1);
    return rotation_mat;

}

//----------------------------------------------------------------------
std::string current_system_iso_date_str()
{
#ifdef WIN32
    SYSTEMTIME st;
    ::GetSystemTime(&st);
    /*
      WORD wYear;
      WORD wMonth;
      WORD wDayOfWeek;
      WORD wDay;
      WORD wHour;
      WORD wMinute;
      WORD wSecond;
      WORD wMilliseconds;
    */
    std::ostringstream os;
    // "1993-6-30"
    os << st.wYear << "-" << st.wMonth << "-" << st.wDay;

    return os.str();
#else
    // Linux
    time_t timep;
    time(&timep);
    struct tm time_struct_data;
    localtime_r(&timep, &time_struct_data);

    std::ostringstream os;
    os << (time_struct_data.tm_year + 1900) << "-" 
       << (time_struct_data.tm_mon + 1)     << "-" 
       << time_struct_data.tm_mday;
    return os.str();
#endif  // WIN32
}

//----------------------------------------------------------------------
std::string current_system_iso_time_str()
{
#ifdef WIN32
    SYSTEMTIME st;
    ::GetSystemTime(&st);
    /*
      WORD wYear;
      WORD wMonth;
      WORD wDayOfWeek;
      WORD wDay;
      WORD wHour;
      WORD wMinute;
      WORD wSecond;
      WORD wMilliseconds;
    */
    std::ostringstream os;
    // "1993-6-30"
    os << st.wHour << ":" << st.wMinute << ":" << st.wSecond;

    return os.str();
#else
    // Linux
    time_t timep;
    time(&timep);
    struct tm time_struct_data;
    localtime_r(&timep, &time_struct_data);

    // mi::Sint32 tm_sec;			/* Seconds.	[0-60] (1 leap second) */
    // mi::Sint32 tm_min;			/* Minutes.	[0-59] */
    // mi::Sint32 tm_hour;			/* Hours.	[0-23] */

    const mi::Size BUFSZ = 256;
    char buf[BUFSZ];
    snprintf(buf, BUFSZ, "%02d:%02d:%02d", 
             time_struct_data.tm_hour,
             time_struct_data.tm_min,
             time_struct_data.tm_sec);
    const std::string ret_str = buf;

    return ret_str;
#endif  // WIN32
}


//----------------------------------------------------------------------
std::string current_system_calender_str()
{
#ifdef WIN32
//#  error "Not implemented this on WIN32"
    SYSTEMTIME st;
    ::GetSystemTime(&st);
    /*
  WORD wYear;
  WORD wMonth;
  WORD wDayOfWeek;
  WORD wDay;
  WORD wHour;
  WORD wMinute;
  WORD wSecond;
  WORD wMilliseconds;
  */
    std::ostringstream os;

    // "Wed Jun 30 21:49:08 1993\n"
    switch (st.wDay) {
        case 0: os << "Sun"; break;
        case 1: os << "Mon"; break;
        case 2: os << "Tue"; break;
        case 3: os << "Wed"; break;
        case 4: os << "Thu"; break;
        case 5: os << "Fru"; break;
        case 6: os << "Sat"; break;
        default: os << "???";
    }
    os << " ";
    switch (st.wMonth) {
        case  1: os << "Jan"; break;
        case  2: os << "Feb"; break;
        case  3: os << "Mar"; break;
        case  4: os << "Apr"; break;
        case  5: os << "May"; break;
        case  6: os << "Jun"; break;
        case  7: os << "Jul"; break;
        case  8: os << "Aug"; break;
        case  9: os << "Sep"; break;
        case 10: os << "Oct"; break;
        case 12: os << "Nov"; break;
        case 13: os << "Dec"; break;
    }
    os << " " << st.wHour
       << ":" << st.wMinute
       << ":" << st.wSecond
       << " " << st.wYear;

    return os.str();
#else
    // Linux
    mi::Sint32 const BUF_LEN = 64;     // see man ctime()
    char buf[BUF_LEN];
    time_t time_st;
    time(&time_st);
    ctime_r(&time_st, buf);
    std::string chopbuf(buf);
    const size_t fpos = chopbuf.find("\n");
    std::string retbuf = chopbuf;
    if(fpos != std::string::npos){
        retbuf = chopbuf.substr(0, chopbuf.size() - 1);
    }
    return retbuf;
#endif  // WIN32
}

//----------------------------------------------------------------------
bool bbox_inclusive_contain(
    mi::math::Bbox< mi::Sint64, 3 > const & bbox0,
    mi::math::Bbox< mi::Sint64, 3 > const & bbox1)
{
    if(bbox0.empty()){
        DEBUG_LOG << "bbox_inclusive_contain: bbox0 is empty.";
        return false;
    }
    if(bbox1.empty()){
        DEBUG_LOG << "bbox_inclusive_contain: bbox1 is empty.";
        return false;
    }

    mi::math::Bbox<mi::Sint64, 3> const and_bbox = nv::index_common::bbox_and(bbox0, bbox1);
    if(bbox1 == and_bbox){
        return true;
    }
    return false;
}

//----------------------------------------------------------------------
bool patch_bbox_inclusive_contain(
    mi::math::Bbox< mi::Sint64, 3 > const & patch_bbox0,
    mi::math::Bbox< mi::Sint64, 3 > const & patch_bbox1)
{
    if(patch_bbox0.empty()){
        DEBUG_LOG << "patch_bbox_inclusive_contain: patch_bbox0 is empty.";
        return false;
    }
    if(patch_bbox1.empty()){
        DEBUG_LOG << "patch_bbox_inclusive_contain: patch_bbox1 is empty.";
        return false;
    }

    mi::math::Bbox<mi::Sint64, 3> const and_bbox = nv::index_common::bbox_and(patch_bbox0, patch_bbox1);
    bool inc_contain = true;
    for(mi::Sint32 i = 0; i < 2; ++i){
        if(patch_bbox1.min[i] != and_bbox.min[i]){
            inc_contain = false;
            break;
        }
        if(patch_bbox1.max[i] != and_bbox.max[i]){
            inc_contain = false;
            break;
        }
    }
    return inc_contain;
}

//----------------------------------------------------------------------
mi::math::Bbox<mi::Float32, 3> get_transformed_bbox(
    const mi::math::Bbox<mi::Float32, 3>&      bbox,
    const mi::math::Matrix<mi::Float32, 4, 4>& transform_mat)
{
    const mi::math::Vector<mi::Float32, 3> new_min(transform_point(transform_mat, bbox.min));
    const mi::math::Vector<mi::Float32, 3> new_max(transform_point(transform_mat, bbox.max));
    const mi::math::Bbox< mi::Float32, 3 > res_bbox(new_min, new_max);
    return res_bbox;
}

//----------------------------------------------------------------------
bool is_bool_state_change_to(bool b_old, bool b_new, bool b_expected)
{
    if((b_old != b_new) && (b_new == b_expected)){ // changed && expected
        return true;
    }
    return false;
}

//----------------------------------------------------------------------
std::string get_app_build_info()
{
    std::string app_build_info_str;
    // Output build environment (when available). Please do not change
    // this code (also the macro in the Makefile), this is for
    // detecting build environment dependent problem.
#ifdef APP_BUILD_ENVIRONMENT
    app_build_info_str = std::string(MI_BASE_STRINGIZE(APP_BUILD_ENVIRONMENT)) + " ";
#else
    app_build_info_str = "user ";
#endif  // APP_BUILD_ENVIRONMENT

    // It seems intel compiler also defines __GNUC__.
#ifdef __GNUC__
#  ifdef __INTEL_COMPILER
    app_build_info_str += "icc ";
#  else
    app_build_info_str += "gcc ";
#  endif
#elif _MSC_VER
    app_build_info_str += "vc ";
#else
    app_build_info_str += "(unknown compiler) ";
#endif

    return app_build_info_str;
}

//----------------------------------------------------------------------
nv::index_common::String_dict get_prefix_entry_key(
    std::string const & prefix,
    std::vector< nv::index_common::String_dict > const & option_vec)
{
    nv::index_common::String_dict extracted_opt;
    for(std::vector< nv::index_common::String_dict >::const_iterator si = option_vec.begin();
        si != option_vec.end();
        ++si)
    {
        nv::index_common::String_dict one_opt;
        mi::Sint32 const item_count =
            string_dict_key_prefix_filter(*si, prefix, one_opt);
        if(item_count == 0){
            continue;           // no prefix entry
        }
        extracted_opt.insert_all(one_opt);
    }

    return extracted_opt;
}

namespace {

std::string expand_environment_variables(
    const std::string&     s,
    bool&                  found,
    std::string::size_type pos = 0)
{
    const std::string::size_type pos_start = s.find("${", pos);
    if (pos_start == std::string::npos)
    {
        return s;
    }

    const std::string pre  = s.substr(0, pos_start);
    std::string post = s.substr(pos_start + 2);

    const std::string::size_type pos_end = post.find("}");
    if (pos_end == std::string::npos)
    {
        return s;
    }

    const std::string name = post.substr(0, pos_end);

    // Only allow access to environment variable with this prefix as a basic safety measure
    const std::string prefix = "NVINDEX_";
    if (name.substr(0, prefix.size()) == prefix)
    {
        post = post.substr(pos_end + 1);

        std::string val = "";
        const char* env = getenv(name.c_str());
        if (env)
        {
            val = std::string(env);
        }
        found = true;

        return expand_environment_variables(pre + val + post, found, pre.size() + val.size());
    }
    else
    {
        ERROR_LOG << "Variable '" << name << "' ignored because it does not start "
                  << "with '" << prefix << "', no further variable expansion will be applied on this line.";
        return s;
    }
}

} // namespace

//----------------------------------------------------------------------
bool load_application_project_file(std::string const & app_prj_fname,
                                   nv::index_common::String_dict  & app_prj_dict)
{
    mi::Sint32 file_version = -1;
    std::string err_mes;
    bool const ret = get_string_dict_with_magic(app_prj_fname,
                                                "index_app_project",
                                                app_prj_dict, file_version, &err_mes);
    if(!ret){
        ERROR_LOG << "Cannot load application project file [" << app_prj_fname << "].";
        ERROR_LOG << err_mes;
        return false;
    }
    mi::Sint32 const app_version = 0;
    if(file_version != app_version){
        ERROR_LOG << "Wrong project file version. version = " << file_version
                  << ", instead of " << app_version;
        return false;
    }

    // Expand environment variables such as ${NVINDEX_FOO}, but only if explicitly enabled
    if (nv::index_common::get_bool(app_prj_dict.get("app::variable_expansion", "no")))
    {
        using nv::index_common::String_dict;
        for (String_dict::const_iterator it = app_prj_dict.begin(); it != app_prj_dict.end(); ++it)
        {
            bool found = false;
            const std::string s = expand_environment_variables(it->second, found);
            if (found)
            {
                app_prj_dict.set(it->first, s);
            }
        }
    }

    return true;
}

//----------------------------------------------------------------------
nv::index_common::String_dict get_extra_project_option(
    nv::index_common::String_dict const & project_file_opt)
{
    std::string const prefix = "app::project_file::";
    nv::index_common::String_dict file_opt;
    mi::Sint32 const item_count =
        string_dict_key_prefix_filter(project_file_opt, prefix, file_opt);
    if(item_count == 0){
        nv::index_common::String_dict empty_string_opt;
        return empty_string_opt; // no extra project option
    }

    std::vector< std::string > keyvec;
    for(nv::index_common::String_dict::const_iterator si = file_opt.begin(); si != file_opt.end(); ++si){
        keyvec.push_back(si->first);
    }
    // get the lex-ordered vector
    std::sort(keyvec.begin(), keyvec.end());

    nv::index_common::String_dict ext_opt;
    for(std::vector< std::string >::const_iterator ki = keyvec.begin(); ki != keyvec.end(); ++ki){
        std::string const prj_fname = file_opt.get(*ki);
        INFO_LOG << "Loading extra project file [" << (*ki) << "] = [" << prj_fname <<"]";

        nv::index_common::String_dict opt_buf;
        if(!load_application_project_file(prj_fname, opt_buf)){
            ERROR_LOG << "Failed to load extra project file [" << prj_fname << "]. skip.";
        }
        ext_opt.insert_all(opt_buf);
    }

    //     if(get_bool(com_line_opt.get("dice::verbose", "0"))){
    //         ext_opt.write(INFO_LOG, "[verbose] extprj: ");
    //     }

    return ext_opt;
}

//----------------------------------------------------------------------
void check_project_compatibility(nv::index_common::String_dict & app_proj)
{
    std::vector<std::string> key_vec;
    std::vector<std::string> mes_vec;
    std::vector<std::string> rep_vec;

    key_vec.push_back("dice::network::per_ip_config::number_of_remote_host");
    mes_vec.push_back("This key is no longer supported and ignored.");
    rep_vec.push_back("");      // no replacement

    key_vec.push_back("app::examiner::view_all_view_type");
    mes_vec.push_back("This key has been deprecated."); // 2015-11
    rep_vec.push_back("app::examiner::predefined_view function.");

    key_vec.push_back("app::examiner::view_all_slack_factor");
    mes_vec.push_back("This key is no longer supported and ignored.");
    rep_vec.push_back("");      // no replacement 2015-11

    key_vec.push_back("app::examiner::view_all::slack_factor");
    mes_vec.push_back("This key is no longer supported and ignored.");
    rep_vec.push_back("");      // no replacement 2015-11

    key_vec.push_back("app::examiner::view_all_scene_bbox");
    mes_vec.push_back("This key has been deprecated. Please use app::examiner::predefined_view function.");
    rep_vec.push_back("");      // 2015-11

    key_vec.push_back("app::examiner::translation_auto_center_mode");
    mes_vec.push_back("This key is no longer supported and ignored.");
    rep_vec.push_back("");      // no replacement 2015-11

    key_vec.push_back("index::window_resolution");
    mes_vec.push_back("This key is deprecated, please use index::canvas_resolution.");
    rep_vec.push_back("index::canvas_resolution");

    key_vec.push_back("index::camera::window_resolution");
    mes_vec.push_back("This key is deprecated, please use index::canvas_resolution.");
    rep_vec.push_back("index::canvas_resolution");

    key_vec.push_back("index::image_file_window_resolution");
    mes_vec.push_back("This key is deprecated, please use index::image_file_canvas_resolution.");
    rep_vec.push_back("index::image_file_canvas_resolution");
     
    key_vec.push_back("index::CUDA_volume_memory");
    mes_vec.push_back("This key is no longer supported and ignored.");
    rep_vec.push_back("");      // no replacement

    key_vec.push_back("index::camera::to");
    mes_vec.push_back("This key is deprecated. Please update as index::camera::dir value.");
    rep_vec.push_back(""); // 2015-11

    key_vec.push_back("app::examiner::front_view::view");
    mes_vec.push_back("This key is deprecated. Please update with app::examiner::predefined_view method.");
    rep_vec.push_back(""); // 2015-11

    key_vec.push_back("app::examiner::front_view::up");
    mes_vec.push_back("This key is deprecated. Please update with app::examiner::predefined_view method.");
    rep_vec.push_back(""); // 2015-11

    key_vec.push_back("app::examiner::side_view::view");
    mes_vec.push_back("This key is deprecated. Please update with app::examiner::predefined_view method.");
    rep_vec.push_back(""); // 2015-11

    key_vec.push_back("app::examiner::side_view::up");
    mes_vec.push_back("This key is deprecated. Please update with app::examiner::predefined_view method.");
    rep_vec.push_back(""); // 2015-11

    key_vec.push_back("app::examiner::top_view::view");
    mes_vec.push_back("This key is deprecated. Please update with app::examiner::predefined_view method.");
    rep_vec.push_back(""); // 2015-11

    key_vec.push_back("app::examiner::top_view::up");
    mes_vec.push_back("This key is deprecated. Please update with app::examiner::predefined_view method.");
    rep_vec.push_back(""); // 2015-11

    key_vec.push_back("app::examiner::quarter_view::view");
    mes_vec.push_back("This key is deprecated. Please update with app::examiner::predefined_view method.");
    rep_vec.push_back(""); // 2015-11

    key_vec.push_back("app::examiner::quarter_view::up");
    mes_vec.push_back("This key is deprecated. Please update with app::examiner::predefined_view method.");
    rep_vec.push_back(""); // 2015-11

    key_vec.push_back("app::examiner::is_navigation_dive_in_when_zoom");
    mes_vec.push_back("This key is no longer supported and ignored.");
    rep_vec.push_back(""); // 2015-11

    key_vec.push_back("app::examiner::is_navigation_pick_to_lookat");
    mes_vec.push_back("This key is no longer supported and ignored.");
    rep_vec.push_back(""); // 2015-11

    key_vec.push_back("app::examiner::startup_view_all");
    mes_vec.push_back("This key is deprecated. Please use app::examiner::predefined_view::startup_update.");
    rep_vec.push_back("app::examiner::predefined_view::startup_update"); // 2015-12

    key_vec.push_back("app::examiner::initial_rotation_center");
    mes_vec.push_back("This key is no longer supported and ignored. Please use app::examiner::initial_rotation_center::type.");
    rep_vec.push_back(""); // 2015-12 Can be removed milestone 3 (never exposed)

    key_vec.push_back("view_all::predefined_view_all::enable");
    mes_vec.push_back("This key is no longer supported and ignored.");
    rep_vec.push_back(""); // 2015-12 Can be removed milestone 3 (never exposed)

    key_vec.push_back("app::examiner::rotation_center_any");
    mes_vec.push_back("This key is no longer supported and ignored.");
    rep_vec.push_back(""); // 2015-12 Can be removed milestone 3 (never exposed)

    // RTMP video stream key update due to other video stream (html5, bridge) support
    key_vec.push_back("dice::video_streaming");
    mes_vec.push_back("This key has been deprecated, please use dice::rtmp_video_streaming::enabled due to multiple kinds of video stream support.");
    rep_vec.push_back("dice::rtmp_video_streaming::enabled"); // 2016-08 

    key_vec.push_back("dice::video_codec");
    mes_vec.push_back("This key has been deprecated, please use dice::rtmp_video_streaming::video_codec due to multiple kinds of video stream support.");
    rep_vec.push_back("dice::rtmp_video_streaming::video_codec"); // 2016-08 

    key_vec.push_back("dice::video_codec_fallback");
    mes_vec.push_back("This key has been deprecated, please use dice::rtmp_video_streaming::video_codec_fallback due to multiple kinds of video stream support.");
    rep_vec.push_back("dice::rtmp_video_streaming::video_codec_fallback"); // 2016-08 

    key_vec.push_back("dice::video_bitrate");
    mes_vec.push_back("This key has been deprecated, please use dice::rtmp_video_streaming::video_codec_bitrate due to multiple kinds of video stream support.");
    rep_vec.push_back("dice::rtmp_video_streaming::video_bitrate"); // 2016-08 

    key_vec.push_back("dice::video_framerate");
    mes_vec.push_back("This key has been deprecated, please use dice::rtmp_video_streaming::video_framerate due to multiple kinds of video stream support.");
    rep_vec.push_back("dice::rtmp_video_streaming::video_framerate"); // 2016-08 
    
    key_vec.push_back("dice::video_bitrate_error");
    mes_vec.push_back("This key has been deprecated, please use dice::rtmp_video_streaming::video_bitrate_error due to multiple kinds of video stream support.");
    rep_vec.push_back("dice::rtmp_video_streaming::video_bitrate_error"); // 2016-08 

    key_vec.push_back("dice::video_preset");
    mes_vec.push_back("This key has been deprecated, please use dice::rtmp_video_streaming::video_preset due to multiple kinds of video stream support.");
    rep_vec.push_back("dice::rtmp_video_streaming::video_preset"); // 2016-08 

    key_vec.push_back("dice::rtmp_port");
    mes_vec.push_back("This key has been deprecated, please use dice::rtmp_video_streaming::port due to multiple kinds of video stream support.");
    rep_vec.push_back("dice::rtmp_video_streaming::port"); // 2016-08 

    key_vec.push_back("dice::rtmp_listen");
    mes_vec.push_back("This key has been deprecated, please use dice::rtmp_video_streaming::listen due to multiple kinds of video stream support.");
    rep_vec.push_back("dice::rtmp_video_streaming::listen"); // 2016-08 

    key_vec.push_back("dice::flash_client");
    mes_vec.push_back("This key has been deprecated, please use dice::rtmp_video_streaming::flash_client due to multiple kinds of video stream support.");
    rep_vec.push_back("dice::rtmp_video_streaming::flash_client"); // 2016-08 

    assert(key_vec.size() == mes_vec.size());

    const mi::Uint32 key_count = key_vec.size();
    for (mi::Uint32 i = 0; i < key_count; ++i)
    {
        if (app_proj.is_defined(key_vec.at(i)))
        {
            WARN_LOG << "Project key [" << key_vec.at(i) << "]: " << mes_vec.at(i);
            if (!(rep_vec.at(i).empty()))
            {
                const std::string old_var = app_proj.get(key_vec.at(i));
                app_proj.insert(rep_vec.at(i), old_var);
                INFO_LOG << "Transfer the value from old key ["
                         << key_vec.at(i) << "] to the new key [" << rep_vec.at(i) 
                         << "] with value [" << old_var << "].";
            }
        }
    }

    // check the indeX key still there
    for (nv::index_common::String_dict::const_iterator ai = app_proj.begin(); ai != app_proj.end(); ++ai)
    {
        if (ai->first.find("indeX") == 0)
        {
            ERROR_LOG << "Found [" << ai->first << "] key that has 'indeX::'! This is obsoleted key. "
                      << "Use index (all lower case key) instead.";
        }
    }
}


//======================================================================
// Ring_buffer_stat
//----------------------------------------------------------------------
void Ring_buffer_stat::initialize()
{
    m_max = 0.0f;
    m_ave = 0.0f;
    m_var = 0.0f;
    m_sample_count = 0;
}

//----------------------------------------------------------------------
std::string Ring_buffer_stat::to_string() const
{
    std::stringstream sstr;
    sstr << "(max, ave, var, samples) = ("
         <<  m_max << ", " << m_ave << ", " << m_var << ", "
         << m_sample_count << ")";
    return sstr.str();
}

//----------------------------------------------------------------------
bool get_ring_buffer_stat(Ring_buffer * p_rbuf, Ring_buffer_stat & rbs)
{
    rbs.initialize();

    mi::Float32 maxfps = 0.0f;
    mi::Float32 avefps = 0.0f;
    // mi::Float32 varfps = 0.0f;

    if(p_rbuf->size() == 0){
        return false;
    }

    for(Ring_buffer::iterator ri = p_rbuf->begin(); ri != p_rbuf->end(); ++ri){
        mi::Float32 const curval = static_cast<mi::Float32>(*ri);

        // max
        if(maxfps < curval){
            maxfps = curval;
        }
        // ave
        avefps += curval;
        // var
        // NIN;
    }
    avefps /= static_cast< mi::Float32 >(p_rbuf->size());
    rbs.m_max = maxfps;
    rbs.m_ave = avefps;
    rbs.m_sample_count = p_rbuf->size();

    return true;
}

//----------------------------------------------------------------------
Cumulative_stat::Cumulative_stat()
    :
    m_max(0.0f),
    m_sum(0.0f),
    m_sample_count(0)
{
    this->clear();
}

//----------------------------------------------------------------------
Cumulative_stat::~Cumulative_stat()
{
    // empty
}

//----------------------------------------------------------------------
void Cumulative_stat::clear()
{
    m_max = -std::numeric_limits< mi::Float32 >::max();
    m_sum = 0.0;
    m_sample_count = 0;
}

//----------------------------------------------------------------------
bool Cumulative_stat::is_valid_stat() const
{
    if(m_sample_count == 0){
        return false;
    }
    return true;
}

//----------------------------------------------------------------------
void Cumulative_stat::add_sample(mi::Float64 dat)
{
    if(m_max < dat){
        m_max = dat;
    }
    m_sum += dat;
    ++m_sample_count;
    // std::cout << "DEBUG: sum = " << m_sum << ", samples = " << m_sample_count << std::endl;
}

//----------------------------------------------------------------------
mi::Float64 Cumulative_stat::get_max() const
{
    return m_max;
}

//----------------------------------------------------------------------
mi::Float64 Cumulative_stat::get_ave() const
{
    if(m_sample_count == 0){
        ERROR_LOG << "no samples, no average.";
        return 0.0;
    }
    mi::Float64 const ave = m_sum / static_cast< mi::Float64 >(m_sample_count);
    return ave;
}

//----------------------------------------------------------------------
mi::Sint32  Cumulative_stat::get_sample_count() const
{
    return m_sample_count;
}

//----------------------------------------------------------------------
std::string Cumulative_stat::to_string() const
{
    if(!this->is_valid_stat()){
        return std::string("Error! no samples.");
    }

    std::stringstream sstr;
    sstr << "(max, ave, samples) = ("
         <<  this->get_max() << ", " << this->get_ave() << ", "
         << this->get_sample_count() << ")";
    return sstr.str();
}

//----------------------------------------------------------------------

std::string escape_JSON(const std::string& unescaped)
{
    std::string esc;
    esc.reserve(unescaped.length()); // Escaped string will be at least this long
    for (size_t i=0; i < unescaped.size(); ++i)
    {
        const unsigned char c = unescaped[i];
        switch (c)
        {
          // Use two-character sequence escape representation, but do not perform the optional
          // escaping of the slash (solidus) character
          case '"':  esc += "\\\""; break; // quotation mark
          case '\\': esc += "\\\\"; break; // backslash (reverse solidus)
          case '\b': esc += "\\b";  break; // backspace
          case '\f': esc += "\\f";  break; // form feed
          case '\n': esc += "\\n";  break; // newline (line feed)
          case '\r': esc += "\\r";  break; // carriage return
          case '\t': esc += "\\t";  break; // tab
          default:
              if (c > 0x1f)
              {
                  // Character doesn't need to be escaped
                  esc += c;
              }
              else
              {
                  // Escape control characters in hex, as in "\u001C"
                  std::ostringstream os;
                  os << "\\u"
                     << std::hex << std::uppercase << std::setw(4) << std::setfill('0')
                     << static_cast<int>(c);
                  esc += os.str();
              }
        }
    }

    return esc;
}

//----------------------------------------------------------------------
