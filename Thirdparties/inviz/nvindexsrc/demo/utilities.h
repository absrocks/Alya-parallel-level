/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief utility functions

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_UTILITIES_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_UTILITIES_H

#include <mi/base/types.h>
#include <mi/math/vector.h>
#include <mi/math/matrix.h>
#include <mi/math/color.h>
#include <mi/math/bbox.h>

#include <string>
#include <vector>
#include <cassert>

#include "common/forwarding_logger.h"

// forward declaration
namespace nv {
namespace index_common {
class String_dict;
}} // namespace

//======================================================================
// conversion: bbox_sint64_3 to bbox_*

//----------------------------------------------------------------------
/// convert bounding box (Sint64_3) to bbox (Sint64_2)
/// \param[in] ijk_bbox_64_3 volume bounding box
/// \return patch bounding box (only x, y)
inline mi::math::Bbox<mi::Sint64, 2> convert_bbox_sint64_3_to_sint64_2(
    mi::math::Bbox<mi::Sint64, 3> const & ijk_bbox_64_3)
{
    mi::math::Bbox<mi::Sint64, 2> ij_bbox_64_2(
        ijk_bbox_64_3.min.x, ijk_bbox_64_3.min.y,
        ijk_bbox_64_3.max.x, ijk_bbox_64_3.max.y);

    return ij_bbox_64_2;
}

//----------------------------------------------------------------------
/// Concatenate names with a separator. The separator is only inserted
/// inbetween strings.
///
/// ex.
///   concat_string_with_separator({"foo", "bar", "baz"}, "|") -> "foo|bar|baz".
///
/// \param[in] str_vec  string vector to be concatenated.
/// \param[in] sep_str  separator string
/// \return concatenated names
std::string concat_string_with_separator(
    std::vector< std::string > const & str_vec,
    std::string const & sep_str);

//----------------------------------------------------------------------
/// get region with interval size.
///
/// Input:
///        [a,b) and an interval size s, where s > 0, b-a > 0
///
/// Output: {r}
///
///    (1) r in [a+ks,a+s*(k+1)) for k = {0, 1, .. (b-a)/s}
///        if (b-a) mod s = 0
///
///    (2) r in [a+ks,a+s*(k+1)) for k = {0, 1, .. floor((b-a)/s)}
///              with [a+floor((b-a)/s)*s,b) otherwise.
///
/// Example: Convert [region_start, region_max) to
///  [region_start,                  region_start +     interval_size),
///  [region_start + interval_size , region_start + 2 * interval_size),
///  ...
///  [region_start + k   * interval_size , region_start + (k+1) * interval_size),
///  ...
///  [region_start + n-1 * interval_size , region_max)   // reminder
///
/// \param[in]  total_range_bbox    [region_start, region_max)
/// \param[in]  interval_size interval size. should be > 0.
/// \param[out] region_vec    (output) result of regions.
void get_region_vector(mi::math::Vector<mi::Sint32, 2> const & total_range_bbox,
                       mi::Sint32 interval_size,
                       std::vector< mi::math::Vector<mi::Sint32, 2> > & region_vec);


//----------------------------------------------------------------------
/// meshgrid function for vector 4
///
/// \param[in]  xrange x range
/// \param[in]  yrange y range
/// \param[in]  zrange z range
/// \param[out] meshgrid4_result meshgrid result.
void meshgrid4(mi::math::Vector<mi::Sint64, 4> const & xrange,
               mi::math::Vector<mi::Sint64, 4> const & yrange,
               mi::math::Vector<mi::Sint64, 4> const & zrange,
               mi::math::Vector<mi::Sint64, 3> * meshgrid4_result);

//----------------------------------------------------------------------
/// meshgrid support function for vector 4
///
/// \param[in]  ix x range
/// \param[in]  iy y range
/// \param[in]  iz z range
/// \return     index of array 64
mi::Sint32 meshgrid4_get_idx(mi::Sint32 ix, mi::Sint32 iy, mi::Sint32 iz);

//----------------------------------------------------------------------
/// check necessary keys in the String_dict. If something is missing,
/// prints out an error and returns false.
bool check_necessary_key(nv::index_common::String_dict const & dict,
                         char const * const p_key[]);

//----------------------------------------------------------------------

// ---------------------------------------------------------------------
inline bool intersection(
    mi::Float32                                     Y,
    const mi::math::Vector_struct<mi::Float32, 2>&  from,
    const mi::math::Vector_struct<mi::Float32, 2>&  to,
    mi::Float32&                                    x_value)
{
    if(from.y<Y && to.y<Y)
        return false;

    if(from.y>Y && to.y>Y)
        return false;

    x_value = to.x + (((Y - to.y)/(to.y-from.y) ) * (to.x-from.x));

    return true;
}

// ---------------------------------------------------------------------
inline void compute_intersection(
    mi::Float32                                                     Y,
    const std::vector<mi::math::Vector_struct<mi::Float32, 2> >&    polygon,
    std::vector<mi::Float32>&                                       intersecton_spans)
{
    const mi::Size nb_vertices = polygon.size();
    mi::Float32 x_value = 0.0f;
    assert(nb_vertices > 0);    // before does -1
    for(mi::Size i=0; i<nb_vertices-1; ++i)
    {
        const bool intersected = intersection(Y, polygon[i], polygon[i+1], x_value);
        if(intersected)
        {
            intersecton_spans.push_back(x_value);
        }
    }

    const bool intersected = intersection(Y, polygon[nb_vertices-1], polygon[0], x_value);
    if(intersected)
    {
        intersecton_spans.push_back(x_value);
    }

    std::sort(
        intersecton_spans.begin(),
        intersecton_spans.end());
}

//----------------------------------------------------------------------
/// get translation matrix from a translation vector
/// \param[in] translation_vec tranlation vector
/// \return translation matrix
mi::math::Matrix_struct<mi::Float32, 4, 4>
get_translation_matrix(mi::math::Vector_struct<mi::Float32, 3> const & translation_vec);

//----------------------------------------------------------------------
/// get scaling matrix from a scaling vector
/// \param[in] scaling_vec scaling vector
/// \return scaling matrix
mi::math::Matrix_struct<mi::Float32, 4, 4>
get_scaling_matrix(mi::math::Vector_struct<mi::Float32, 3> const & scaling_vec);

//----------------------------------------------------------------------
/// get rotation matrix about X axis for a given angle
/// \param[in] angle
/// \return rotation matrix about X axis
mi::math::Matrix_struct<mi::Float32, 4, 4>
get_rotation_matrix_about_x(mi::Float32 angle);

//----------------------------------------------------------------------
/// get rotation matrix about Y axis for a given angle
/// \param[in] angle
/// \return rotation matrix about Y axis
mi::math::Matrix_struct<mi::Float32, 4, 4>
get_rotation_matrix_about_y(mi::Float32 angle);

//----------------------------------------------------------------------
/// get rotation matrix about Z axis for a given angle
/// \param[in] angle
/// \return rotation matrix about Z axis
mi::math::Matrix_struct<mi::Float32, 4, 4>
get_rotation_matrix_about_z(mi::Float32 angle);

//----------------------------------------------------------------------
/// get current system date as a ISO date (e.g, 2014-6-11)
/// \return current system time as ISO date
std::string current_system_iso_date_str();

//----------------------------------------------------------------------
/// get current system time as a ISO date (e.g, 15:20:09)
/// \return current system time as ISO time
std::string current_system_iso_time_str();

//----------------------------------------------------------------------
/// get current system time as a calender string.
/// \return current system time as a calender string..
std::string current_system_calender_str();

//----------------------------------------------------------------------
/// bounding box 'inclusive_contain' operation
///
/// bbox0 inclusive_contain bbox1 means, bbox1's min, max are inside
/// or also equal on the border of the bbox0.
///
/// \param[in] bbox0 bounding box 0
/// \param[in] bbox1 bounding box 1
/// \return true if bbox0 inclusive contain bbox1
/// if one of the volume is empty, return false and cast message when
/// debug build. (The message is only visible in debug build.)
bool bbox_inclusive_contain(
    mi::math::Bbox< mi::Sint64, 3 > const & bbox0,
    mi::math::Bbox< mi::Sint64, 3 > const & bbox1);

//----------------------------------------------------------------------
/// patch bounding box 'inclusive_contain' operation (as a patch, so
/// only x and y are considered.)
///
/// patch_bbox0 inclusive_contain patch_bbox1 means, x and y of
/// patch_bbox1's min, max are inside or also equal on the border of
/// the bbox0.
///
/// \param[in] bbox0 bounding box 0
/// \param[in] bbox1 bounding box 1
/// \return true if bbox0 inclusive contain bbox1
/// if one of the volume is empty, return false and cast message when
/// debug build. (The message is only visible in debug build.)
bool patch_bbox_inclusive_contain(
    mi::math::Bbox< mi::Sint64, 3 > const & patch_bbox0,
    mi::math::Bbox< mi::Sint64, 3 > const & patch_bbox1);

//----------------------------------------------------------------------
/// get a transformed bbox
///
/// \param[in]  bbox          a bounding box
/// \param[in]  transform_mat a transformation matrix
mi::math::Bbox<mi::Float32, 3> get_transformed_bbox(
    const mi::math::Bbox<mi::Float32, 3>&      bbox,
    const mi::math::Matrix<mi::Float32, 4, 4>& transform_mat);

//----------------------------------------------------------------------
/// is bool state change to a specific state?
///
///
/// \param[in] b_old old state
/// \param[in] b_new new state
/// \param[in] b_expected expected state
/// \return When a bool state is changed b_old to b_new and b_new is
/// b_expected, true.
bool is_bool_state_change_to(bool b_old, bool b_new, bool b_expected);

//----------------------------------------------------------------------
/// get the build environment information
///
/// \return build environment information string
std::string get_app_build_info();

//----------------------------------------------------------------------
/// get the prefix entry keys from multiple String_dicts
///
/// Get the entries that their key has the prefix are extracted with a
/// priority. The latter entry has higher priority. The use case is
/// that getting different priority options, e.g., project file and
/// command line option.
///
/// \param[in] prefix prefix for the extra project file entry. E.g., "project_file::".
/// \param[in] option_vec   option vector
/// \return prefixed entries
nv::index_common::String_dict get_prefix_entry_key(
    std::string const & prefix,
    std::vector< nv::index_common::String_dict > const & option_vec);

//----------------------------------------------------------------------
/// get the extra project file entry
///
/// Look up the "app::project_file::*" keys and load these files and
/// accumulate the option entries.
///
/// \param[in] project_file_opt project file options
/// \return accumulated items
nv::index_common::String_dict get_extra_project_option(
    nv::index_common::String_dict const & project_file_opt);

//----------------------------------------------------------------------
/// load application project file. Check the header and version.
///
/// \param[in]  app_prj_fname application project file
/// \param[out] app_prj_dict  (output) application project contents dictionary
/// \return true when load succeeded.
bool load_application_project_file(std::string const & app_prj_fname,
                                   nv::index_common::String_dict  & app_prj_dict);

//----------------------------------------------------------------------
/// check old project entry and only warn it. No change the app_proj.
///
/// \param[in,out] app_prj application project string dictionary
void check_project_compatibility(nv::index_common::String_dict & app_proj);

//----------------------------------------------------------------------
/// Escape special characters to create a valid JSON string representation.
/// See https://tools.ietf.org/html/rfc7158#section-7
///
/// \param[in] unescaped String that should be escaped
/// \return Escaped string
std::string escape_JSON(const std::string& unescaped);

//----------------------------------------------------------------------

/// ring buffer stat computation
class Ring_buffer_stat
{
public:
    /// default constructor
    Ring_buffer_stat()
        :
        m_max(0.0),
        m_ave(0.0),
        m_var(-1.0),            // invalid
        m_sample_count(-1)
    {
        // empty;
    }

    /// initialize to all zero
    void initialize();

    /// destructor
    ~Ring_buffer_stat()
    {
        // empty
    }

    /// get string representation
    std::string to_string() const;

public:
    /// maximal value
    mi::Float32 m_max;
    /// average value
    mi::Float32 m_ave;
    /// variance value
    mi::Float32 m_var;
    /// sample count
    mi::Sint32  m_sample_count;

private:
    /// copy constructor. prohibit until proved useful.
    Ring_buffer_stat(Ring_buffer_stat const &);
    /// operator=. prohibit until proved useful.
    Ring_buffer_stat const & operator=(Ring_buffer_stat const &);
};

class Ring_buffer;

/// get ring buffer statistics
/// \param[in]  p_rbuf ring buffer
/// \param[out] rbs ring buffer statistics result.
/// \return false if no data in the ring buffer
bool get_ring_buffer_stat(Ring_buffer * p_rbuf, Ring_buffer_stat & rbs);

//----------------------------------------------------------------------

/// cumulative statistics
class Cumulative_stat
{
public:
    /// default constructor
    Cumulative_stat();

    /// destructor
    ~Cumulative_stat();

    /// clear the state (The object becomes invalid since no samples.)
    void clear();

    /// is statistics valid
    bool is_valid_stat() const;

    /// add one sample
    /// \param[in] dat one sample data
    void add_sample(mi::Float64 dat);

    /// get max
    /// \return current max value. -1 if not valid (e.g., no samples)
    mi::Float64 get_max() const;

    /// get average
    /// \return current average value. -1 if not valid (e.g., no samples)
    mi::Float64 get_ave() const;

    /// get number of samples
    /// \return current sample number. 0 if no samples.
    mi::Sint32  get_sample_count() const;

    /// get string representation
    std::string to_string() const;

public:
    /// maximal value
    mi::Float64 m_max;
    /// sample sum value
    mi::Float64 m_sum;
    /// sample count
    mi::Sint32  m_sample_count;

private:
    /// copy constructor. prohibit until proved useful.
    Cumulative_stat(Cumulative_stat const &);
    /// operator=. prohibit until proved useful.
    Cumulative_stat const & operator=(Cumulative_stat const &);
};

//----------------------------------------------------------------------

#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_UTILITIES_H
