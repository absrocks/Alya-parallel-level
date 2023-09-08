/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief camera utlities

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_CAMERA_UTIL_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_CAMERA_UTIL_H

#include "stereo_camera.h"

#include <nv/index/icamera.h>

#include "common/string_dict.h"

#include <vector>

class Nvindex_rendering_context;

//----------------------------------------------------------------------
//======================================================================
//----------------------------------------------------------------------
/// Camera parameter container: FIXME to class, put in common
class Camera_parameter
{
public:
    /// default constructor of camera parameter container
    Camera_parameter();

    /// set camera parameters by string dict
    ///
    /// \param[in] cam_param          a string dict that contains camera parameters
    /// \param[in] string_dict_prefix camera parameter's string dict prefix (ex. index::camera::)
    /// \return number of parameters found in cam_param
    mi::Sint32 set_parameter_by_string_dict(
        const nv::index_common::String_dict& cam_param, 
        const std::string& string_dict_prefix);

    /// set camera parameters into a string dict
    ///
    /// \param[in] string_dict_prefix camera parameter's string dict prefix (ex. index::camera::)
    /// \param[out] cam_param          a string dict that contains camera parameters
    /// FIXME: can be a const method
    void get_parameter_to_string_dict(
        const std::string& string_dict_prefix,
        nv::index_common::String_dict& cam_param);

    /// copy this parameters to a camera
    void copy_this_parameter_to_camera(nv::index::ICamera *cam) const;

    /// copy a camera parameter to this
    void copy_camera_parameter_to_this(const nv::index::ICamera *cam);

    /// get a vector parameters.
    /// 
    /// \param[in]  key keyword to get a Vector
    /// \param[out] retval return value
    /// \return true when valid key
    bool get_vector_value(const std::string& key, mi::math::Vector<mi::Float32, 3>& retval) const;

    /// get a scalar parameters
    /// 
    /// \param[in]  key keyword to get a Float64
    /// \param[out] retval return value
    /// \return true when valid key
    bool get_scalar_value(const std::string& key, mi::Float64& retval) const;

    /// get a bool parameters
    /// 
    /// \param[in]  key keyword to get a bool
    /// \param[out] retval return value
    /// \return true when valid key
    bool get_bool_value(const std::string& key, bool& retval) const;

    /// get a string parameters
    /// 
    /// \param[in]  key keyword to get a string
    /// \param[out] retval return value
    /// \return true when valid key
    bool get_string_value(const std::string& key, std::string& retval) const;

private:
    /// vector parameter keys
    std::vector<std::string> m_vector_key_vec;
    /// scalar parameter keys
    std::vector<std::string> m_scalar_key_vec;
    /// bool parameter keys
    std::vector<std::string> m_bool_key_vec;
    /// string parameter keys
    std::vector<std::string> m_string_key_vec;

    /// from, dir, up vectors
    std::map<std::string, mi::math::Vector<mi::Float32, 3> > m_vector_param_map;
    /// aspect, aperture, focal, clip_min, clip_max
    std::map<std::string, mi::Float64>                       m_scalar_param_map;
    /// orthographic, view_all
    std::map<std::string, mi::Sint32>                        m_bool_param_map;
    /// name
    std::map<std::string, std::string>                       m_string_param_map;

};

/// Camera toolkit methods
class Camera_tool
{
public:
    /// Ctor
    Camera_tool()
    {
        // empty
    }
    /// Dtor
    ~Camera_tool()
    {
        // empty
    }

    /// print camera parameter from camara object
    ///
    /// \param[in] cam camera object
    static void print_camera_param(
        const nv::index::ICamera* cam);

    /// print camera parameter from tag and transaction
    ///
    /// \param[in] camera_tag        camera to be used
    /// \param[in] dice_transaction  dice transaction for edit the camera
    static void print_camera_param(
        const mi::neuraylib::Tag& camera_tag,
        mi::neuraylib::IDice_transaction* dice_transaction);

    /// get camera parameter from string
    ///
    /// \param[in]  cam_param_str camera parameter string
    /// \param[out] eye    (output) eye position
    /// \param[out] dir    (output) direction vector
    /// \param[out] up     (output) up vector
    /// \param[out] error_mes (output) error message when returns false
    /// \return true when cam_param_str is valid.
    static bool get_camera_param_from_string(
        const std::string& cam_param_str,
        mi::math::Vector< mi::Float32, 3 >& eye,
        mi::math::Vector< mi::Float32, 3 >& dir,
        mi::math::Vector< mi::Float32, 3 >& up,
        std::string & error_mes);
 
    /// compute the camera basis vectors
    ///
    /// \param[in]  view_dir camera viewing direction
    /// \param[in]  up_dir   camera up direction
    /// \param[out] xdir     (output) basis x
    /// \param[out] ydir     (output) basis y
    /// \param[out] zdir     (output) basis z
    ///
    static void compute_camera_basis(
        const mi::math::Vector<mi::Float32, 3>& view_dir,
        const mi::math::Vector<mi::Float32, 3>& up_dir,
        mi::math::Vector<mi::Float32, 3>&       xdir,
        mi::math::Vector<mi::Float32, 3>&       ydir,
        mi::math::Vector<mi::Float32, 3>&       zdir);
 
    /// Get gluLookAt vectors
    ///
    /// \param[in]  cam    ICamera object
    /// \param[out] eye    (output) eye position vector
    /// \param[out] center (output) lookat center position vector
    /// \param[out] up     (output) up vector
    ///
    static void get_glu_lookat_vector(
        const nv::index::ICamera*         cam,
        mi::math::Vector<mi::Float32, 3>& eye,
        mi::math::Vector<mi::Float32, 3>& center,
        mi::math::Vector<mi::Float32, 3>& up);
 
    /// get the model view matrix
    ///
    /// \param[in] cam                 a camera
    /// \param[in] scene_transform     scene transformation
    /// \param[in] viewport_resolution viewport resolution
    /// \return MV matrix
    static mi::math::Matrix<mi::Float64, 4, 4> get_mv_matrix(
        nv::index::ICamera*                        cam,
        const mi::math::Matrix<mi::Float32, 4, 4>& scene_transform,
        const mi::math::Vector<mi::Sint32, 2>&     viewport_resolution);

};

//----------------------------------------------------------------------
#endif // #ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_CAMERA_UTIL_H
 
