/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "camera_utility.h"

#include <cassert>
#include <iostream>
#include <iomanip>

#include <mi/base/handle.h>

#include <nv/index/icamera.h>

#include "common/common_utility.h"
#include "common/forwarding_logger.h"

#include "examiner_manipulator.h"

//----------------------------------------------------------------------
// Camera_parameter
//----------------------------------------------------------------------
Camera_parameter::Camera_parameter()
{
    // camera vector parameter keys
    m_vector_key_vec.push_back("from");
    m_vector_key_vec.push_back("dir");
    m_vector_key_vec.push_back("up");

    // camera scalar parameter keys
    m_scalar_key_vec.push_back("aspect");
    m_scalar_key_vec.push_back("aperture");
    m_scalar_key_vec.push_back("focal");
    m_scalar_key_vec.push_back("clip_max");
    m_scalar_key_vec.push_back("clip_min");

    // camera bool parameter keys
    m_bool_key_vec.push_back("orthographic");
    m_bool_key_vec.push_back("view_all");

    // camera string parameter keys
    m_string_key_vec.push_back("name");               // for GUI (predefined_view)
    m_string_key_vec.push_back("view_all_bbox_type"); // for predefined_view

    // set default values
    m_vector_param_map["from"] = mi::math::Vector<mi::Float32, 3>(0.0f, 0.0f, -5.0f);
    m_vector_param_map["dir"]  = mi::math::Vector<mi::Float32, 3>(0.0f, 0.0f,  0.0f);
    m_vector_param_map["up"]   = mi::math::Vector<mi::Float32, 3>(0.0f, 1.0f,  0.0f);

    m_scalar_param_map["aspect"]   = 1.0;
    m_scalar_param_map["aperture"] = 0.033;
    m_scalar_param_map["focal"]    = 0.03;
    m_scalar_param_map["clip_min"] = 0.01;
    m_scalar_param_map["clip_max"] = 1000.0;

    m_bool_param_map["orthographic"] = 0; // bool is represented by Sint32
    m_bool_param_map["view_all"]     = 0; // for predefined_view

    m_string_param_map["name"]       = "Predef_view";        // for predefined_view
    m_string_param_map["view_all_bbox_type"]  = "scene_roi"; // for predefined_view
}

//----------------------------------------------------------------------
mi::Sint32 Camera_parameter::set_parameter_by_string_dict(
    const nv::index_common::String_dict& cam_param,
    const std::string& string_dict_prefix)
{
    mi::Sint32 found_param_count = 0;
    const mi::Size nb_vector_param = m_vector_key_vec.size();
    for (mi::Size i = 0; i < nb_vector_param; ++i)
    {
        const std::string cur_key = m_vector_key_vec[i];
        if (cam_param.is_defined(string_dict_prefix + cur_key))
        {
            m_vector_param_map[cur_key] =
                nv::index_common::get_vec_float32_3(cam_param.get(string_dict_prefix + cur_key));
            ++found_param_count;
        }
    }

    const mi::Size nb_scalar_param = m_scalar_key_vec.size();
    for (mi::Size i = 0; i < nb_scalar_param; ++i)
    {
        const std::string cur_key = m_scalar_key_vec[i];
        if (cam_param.is_defined(string_dict_prefix + cur_key))
        {
            m_scalar_param_map[cur_key] =
                nv::index_common::get_float64(cam_param.get(string_dict_prefix + cur_key));
            ++found_param_count;
        }
    }

    const mi::Size nb_bool_param = m_bool_key_vec.size();
    for (mi::Size i = 0; i < nb_bool_param; ++i)
    {
        const std::string cur_key = m_bool_key_vec[i];
        if (cam_param.is_defined(string_dict_prefix + cur_key))
        {
            m_bool_param_map[cur_key] =
                nv::index_common::get_bool(cam_param.get(string_dict_prefix + cur_key));
            ++found_param_count;
        }
    }

    const mi::Size nb_string_param = m_string_key_vec.size();
    for (mi::Size i = 0; i < nb_string_param; ++i)
    {
        const std::string cur_key = m_string_key_vec[i];
        if (cam_param.is_defined(string_dict_prefix + cur_key))
        {
            m_string_param_map[cur_key] = cam_param.get(string_dict_prefix + cur_key);
            ++found_param_count;
        }
    }

    // sanity check
    const mi::Float64 clip_min = m_scalar_param_map["clip_min"];
    const mi::Float64 clip_max = m_scalar_param_map["clip_max"];
    if ((clip_min > 0.0) && (clip_max > 0.0))
    {
        if(clip_max <= clip_min)
        {
            ERROR_LOG << "set_parameter_by_string_dict: illegal clipping plane, clip_min must be < clip_max. "
                      << "Setting to 0.01, 10000.0";
            m_scalar_param_map["clip_min"] = 0.01;
            m_scalar_param_map["clip_max"] = 10000.0;
        }
    }
    // else negative, do nothing

    return found_param_count;
}

//----------------------------------------------------------------------
void Camera_parameter::get_parameter_to_string_dict(
    const std::string& string_dict_prefix,
    nv::index_common::String_dict& cam_param)
{
    const mi::Size nb_vector_param = m_vector_key_vec.size();
    for (mi::Size i = 0; i < nb_vector_param; ++i)
    {
        const std::string cur_key = m_vector_key_vec[i];
        cam_param.insert(string_dict_prefix + cur_key,
                         nv::index_common::to_string(m_vector_param_map[cur_key]));
    }

    const mi::Size nb_scalar_param = m_scalar_key_vec.size();
    for (mi::Size i = 0; i < nb_scalar_param; ++i)
    {
        const std::string cur_key = m_scalar_key_vec[i];
        cam_param.insert(string_dict_prefix + cur_key,
                         nv::index_common::to_string(m_scalar_param_map[cur_key]));
    }

    const mi::Size nb_bool_param = m_bool_key_vec.size();
    for (mi::Size i = 0; i < nb_bool_param; ++i)
    {
        const std::string cur_key = m_bool_key_vec[i];
        cam_param.insert(string_dict_prefix + cur_key,
                         nv::index_common::to_string(m_bool_param_map[cur_key]));
    }

    const mi::Size nb_string_param = m_string_key_vec.size();
    for (mi::Size i = 0; i < nb_string_param; ++i)
    {
        const std::string cur_key = m_string_key_vec[i];
        cam_param.insert(string_dict_prefix + cur_key, m_string_param_map[cur_key]);
    }
}

//----------------------------------------------------------------------
static void debug_message_if_lt_0(mi::Float64 val,
                            const std::string& param_name, 
                            const std::string& message)
{
    if (val <= 0.0)
    {
        DEBUG_LOG << param_name << ": " << val << " is <= 0.0. " << message;
    }
}

void Camera_parameter::copy_this_parameter_to_camera(nv::index::ICamera *cam) const
{
    const std::string mn = "Camera_parameter::copy_this_parameter_to_camera: ";
    if (cam == 0)
    {
        ERROR_LOG << mn << "invalid camera.";
        return;
    }

    // vectors
    assert(m_vector_param_map.find("from") != m_vector_param_map.end());
    assert(m_vector_param_map.find("dir")  != m_vector_param_map.end());
    assert(m_vector_param_map.find("up")   != m_vector_param_map.end());

    const mi::math::Vector<mi::Float32, 3> from = m_vector_param_map.find("from")->second;
    const mi::math::Vector<mi::Float32, 3> dir  = m_vector_param_map.find("dir") ->second;
    const mi::math::Vector<mi::Float32, 3> up   = m_vector_param_map.find("up")  ->second;

    // scalars
    assert(m_scalar_param_map.find("aspect")   != m_scalar_param_map.end());
    assert(m_scalar_param_map.find("aperture") != m_scalar_param_map.end());
    assert(m_scalar_param_map.find("focal")    != m_scalar_param_map.end());
    assert(m_scalar_param_map.find("clip_min") != m_scalar_param_map.end());
    assert(m_scalar_param_map.find("clip_max") != m_scalar_param_map.end());

    const mi::Float64 aspect   = m_scalar_param_map.find("aspect")  ->second;
    const mi::Float64 aperture = m_scalar_param_map.find("aperture")->second;
    mi::Float64       clip_min = m_scalar_param_map.find("clip_min")->second;
    mi::Float64       clip_max = m_scalar_param_map.find("clip_max")->second;

    if ((clip_min > 0.0) && (clip_max > 0.0))
    {
        if (clip_max <= clip_min)
        {
            ERROR_LOG << "copy_this_parameter_to_camera: illegal clipping plane, clip_min (" << clip_min
                      << ") must be < clip_max (" << clip_max << "). Setting to 0.01, 100.0";
            clip_min = 0.01;
            clip_max = 100.0;
        }
    }
    else
    {
        // if one of clip_min, clip_max < 0.0, set both -1.0
        clip_min = -1.0;
        clip_max = -1.0;
        DEBUG_LOG << "One of the clip_min, clip_max is < 0.0, set both -1.0.";
    }

    // always set the following three parameters
    cam->set_eye_point(from);
    cam->set_view_direction(dir);
    cam->set_up_direction(up);

    mi::base::Handle<nv::index::IPerspective_camera>  perspective_camera(
        cam->get_interface<nv::index::IPerspective_camera>());
    mi::base::Handle<nv::index::IOrthographic_camera> ortho_camera(
        cam->get_interface<nv::index::IOrthographic_camera>());
    if (perspective_camera.is_valid_interface())
    {
        if (aspect > 0.0)
        {
            perspective_camera->set_aspect(aspect);
        }

        if (aperture > 0.0)
        {
            perspective_camera->set_aperture(aperture);
        }

        const mi::Float64 focal = m_scalar_param_map.find("focal")->second;
        if (focal > 0.0)
        {
            perspective_camera->set_focal(focal);
        }

        if (clip_min > 0.0)
        {
            perspective_camera->set_clip_min(clip_min);
        }

        if (clip_max > 0.0)
        {
            perspective_camera->set_clip_max(clip_max);
        }
        
        debug_message_if_lt_0(aspect,   "aspect",   "ignored.");
        debug_message_if_lt_0(aperture, "aperture", "ignored.");
        debug_message_if_lt_0(focal,    "focal",    "ignored.");
        debug_message_if_lt_0(clip_min, "clip_min", "ignored.");
        debug_message_if_lt_0(clip_max, "clip_max", "ignored.");
    }
    else if (ortho_camera.is_valid_interface())
    {
        if (aspect > 0.0)
        {
            ortho_camera->set_aspect(aspect);
        }

        if (aperture > 0.0)
        {
            ortho_camera->set_aperture(aperture);
        }

        if (clip_min > 0.0)
        {
            ortho_camera->set_clip_min(clip_min);
        }

        if (clip_max > 0.0)
        {
            ortho_camera->set_clip_max(clip_max);
        }

        debug_message_if_lt_0(aspect,   "aspect",   "ignored.");
        debug_message_if_lt_0(aperture, "aperture", "ignored.");
        debug_message_if_lt_0(clip_min, "clip_min", "ignored.");
        debug_message_if_lt_0(clip_max, "clip_max", "ignored.");
    }
    else
    {
        ERROR_LOG << "Unknown camera type.";
    }
}

//----------------------------------------------------------------------
void Camera_parameter::copy_camera_parameter_to_this(const nv::index::ICamera *cam)
{
    const std::string mn = "Camera_parameter::copy_camera_parameter_to_this: ";
    if (cam == 0)
    {
        ERROR_LOG << mn << "invalid camera.";
        return;
    }

    const mi::math::Vector<mi::Float32, 3> from = cam->get_eye_point();
    const mi::math::Vector<mi::Float32, 3> dir  = cam->get_view_direction();
    const mi::math::Vector<mi::Float32, 3> up   = cam->get_up_direction();

    // vectors
    m_vector_param_map["from"] = from;
    m_vector_param_map["dir"]  = dir;
    m_vector_param_map["up"]   = up;

    mi::base::Handle<const nv::index::IPerspective_camera>  perspective_camera(
        cam->get_interface<const nv::index::IPerspective_camera>());
    mi::base::Handle<const nv::index::IOrthographic_camera> ortho_camera(
        cam->get_interface<const nv::index::IOrthographic_camera>());

    // scalars
    if (perspective_camera.is_valid_interface())
    {
        m_scalar_param_map["aspect"]   = perspective_camera->get_aspect();
        m_scalar_param_map["aperture"] = perspective_camera->get_aperture();
        m_scalar_param_map["focal"]    = perspective_camera->get_focal();
        m_scalar_param_map["clip_min"] = perspective_camera->get_clip_min();
        m_scalar_param_map["clip_max"] = perspective_camera->get_clip_max();
    }
    else if (ortho_camera.is_valid_interface())
    {
        m_scalar_param_map["aspect"]   = ortho_camera->get_aspect();
        m_scalar_param_map["aperture"] = ortho_camera->get_aperture();
        m_scalar_param_map["clip_min"] = ortho_camera->get_clip_min();
        m_scalar_param_map["clip_max"] = ortho_camera->get_clip_max();
    }
    else
    {
        ERROR_LOG << "Unknown camera type.";
    }
}

//----------------------------------------------------------------------
bool Camera_parameter::get_vector_value(const std::string& key,
                                        mi::math::Vector<mi::Float32, 3>& retval) const
{
    const std::string mn = "Camera_parameter::get_vector_value: ";
    std::map<std::string, mi::math::Vector<mi::Float32, 3> >::const_iterator vi =
        m_vector_param_map.find(key);
    if (vi == m_vector_param_map.end())
    {
        ERROR_LOG << mn << "No such key: " << key;
        return false;
    }
    
    retval = vi->second;
    return true;
}

//----------------------------------------------------------------------
bool Camera_parameter::get_scalar_value(const std::string& key,
                                        mi::Float64& retval) const
{
    const std::string mn = "Camera_parameter::get_scalar_value: ";
    std::map<std::string, mi::Float64>::const_iterator si =
        m_scalar_param_map.find(key);
    if (si == m_scalar_param_map.end())
    {
        ERROR_LOG << mn << "No such key: " << key;
        return false;
    }

    retval = si->second;    
    return true;
}

//----------------------------------------------------------------------
bool Camera_parameter::get_bool_value(const std::string& key, 
                                      bool& retval) const
{
    const std::string mn = "Camera_parameter::get_bool_value: ";
    std::map<std::string, mi::Sint32>::const_iterator bi =
        m_bool_param_map.find(key);
    if (bi == m_bool_param_map.end())
    {
        ERROR_LOG << mn << "No such key: " << key;
        return false;
    }
    retval = (bi->second != 0); // Sint32 to bool conversion (We don't use an std bit vector)
    
    return true;
}

//----------------------------------------------------------------------
bool Camera_parameter::get_string_value(const std::string& key,
                                        std::string& retval) const
{
    const std::string mn = "Camera_parameter::get_string_value: ";
    std::map<std::string, std::string>::const_iterator si =
        m_string_param_map.find(key);
    if (si == m_string_param_map.end())
    {
        ERROR_LOG << mn << "No such key: " << key;
        return false;
    }

    retval = si->second;
    return true;
}

//----------------------------------------------------------------------
//======================================================================
//----------------------------------------------------------------------
void Camera_tool::print_camera_param(const nv::index::ICamera* cam)
{
    assert(cam != 0);

    // Simple camera parameters
    mi::math::Vector< mi::Float32, 3 > eye = cam->get_eye_point();
    mi::math::Vector< mi::Float32, 3 > dir = cam->get_view_direction();
    mi::math::Vector< mi::Float32, 3 > up  = cam->get_up_direction();

    mi::base::Handle<const nv::index::IPerspective_camera>  perspective_camera(cam->get_interface<const nv::index::IPerspective_camera>());
    mi::base::Handle<const nv::index::IOrthographic_camera> ortho_camera(cam->get_interface<const nv::index::IOrthographic_camera>());

    // Print directly to stdout, for easy cut-and-paste
    std::cout << "Camera parameters:\n"
              << std::setprecision (15) // Need full precision to be able to correctly restore the camera
              << "index::camera::from = "<< eye.x << " " << eye.y << " " << eye.z    << "\n"
              << "index::camera::dir = " << dir.x << " " << dir.y << " " << dir.z << "\n"
              << "index::camera::up = "  << up.x  << " " << up.y  << " " << up.z     << "\n";

    if (perspective_camera)
    {
        std::cout << std::setprecision(15)
                  << "index::camera::aspect = "   << perspective_camera->get_aspect()   << "\n"
                  << "index::camera::aperture = " << perspective_camera->get_aperture() << "\n"
                  << "index::camera::focal = "    << perspective_camera->get_focal()    << "\n"
                  << "index::camera::clip_min = " << perspective_camera->get_clip_min() << "\n"
                  << "index::camera::clip_max = " << perspective_camera->get_clip_max() << "\n"
                  << "index::camera::orthographic = " << "no" << "\n";
    }
    if (ortho_camera)
    {
        std::cout << std::setprecision(15)
                  << "index::camera::aspect = "   << ortho_camera->get_aspect()   << "\n"
                  << "index::camera::aperture = " << ortho_camera->get_aperture() << "\n"
                  << "index::camera::clip_min = " << ortho_camera->get_clip_min() << "\n"
                  << "index::camera::clip_max = " << ortho_camera->get_clip_max() << "\n"
                  << "index::camera::orthographic = " << "yes" << "\n";
    }
    std::cout << "app::examiner::predefined_view::startup_update = no\n";
    std::cout << "\n";

    // Print again in regression test format
    std::cout << std::setprecision (15)
              << "ICamera::Set:\n"
              << "    from: "<< eye.x << " " << eye.y << " " << eye.z << "\n"
              << "    dir: " << dir.x << " " << dir.y << " " << dir.z << "\n"
              << "    up: "  << up.x  << " " << up.y  << " " << up.z  << "\n";

    if (perspective_camera)
    {
        std::cout << std::setprecision(15)
              << "    aspect: "   << perspective_camera->get_aspect()   << "\n"
              << "    aperture: " << perspective_camera->get_aperture() << "\n"
              << "    focal: "    << perspective_camera->get_focal()    << "\n"
              << "    clip_min: " << perspective_camera->get_clip_min() << "\n"
              << "    clip_max: " << perspective_camera->get_clip_max() << "\n"
              << "    orthographic: " << "no" << "\n";
    }
    if (ortho_camera)
    {
        std::cout << std::setprecision(15)
              << "    aspect: "   << ortho_camera->get_aspect()   << "\n"
              << "    aperture: " << ortho_camera->get_aperture() << "\n"
              << "    clip_min: " << ortho_camera->get_clip_min() << "\n"
              << "    clip_max: " << ortho_camera->get_clip_max() << "\n"
              << "    orthographic: " << "yes" << "\n";
    }

    std::cout << std::endl;
}

//----------------------------------------------------------------------
void Camera_tool::print_camera_param(
    const mi::neuraylib::Tag&         camera_tag,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    const mi::base::Handle< const nv::index::ICamera > cam(
        dice_transaction->access< nv::index::ICamera >(camera_tag));
    assert(cam.is_valid_interface());

    print_camera_param(cam.get());
}

//----------------------------------------------------------------------
bool Camera_tool::get_camera_param_from_string(
    const std::string&                  cam_param_str,
    mi::math::Vector< mi::Float32, 3 >& eye,
    mi::math::Vector< mi::Float32, 3 >& dir,
    mi::math::Vector< mi::Float32, 3 >& up,
    std::string &                       error_mes)
{
    std::istringstream isstr(cam_param_str);
    isstr >> eye.x >> eye.y >> eye.z
          >> dir.x >> dir.y >> dir.z
          >> up.x  >> up.y  >> up.z;

    //     INFO_LOG << "DEBUG: cam_pram_str[" << cam_param_str << "]\n"
    //              << "eye: " << eye.x << " " << eye.y << " " << eye.z << "\n"
    //              << "dir: " << dir.x << " " << dir.y << " " << dir.z << "\n"
    //              << "up:  " << up.x  << " " << up.y  << " " << up.z;
    if (!isstr)
    {
        error_mes = "get_camera_param_from_string: invalid cam_param_str. ["
            + cam_param_str + "]";
        return false;
    }

    return true;
}

//----------------------------------------------------------------------
void Camera_tool::compute_camera_basis(
    const mi::math::Vector<mi::Float32, 3>& view_dir,
    const mi::math::Vector<mi::Float32, 3>& up_dir,
    mi::math::Vector<mi::Float32, 3>&       xdir,
    mi::math::Vector<mi::Float32, 3>&       ydir,
    mi::math::Vector<mi::Float32, 3>&       zdir)
{
    // Camera orthonormal basis, left handed
    bool is_normalized = false;

    zdir = view_dir;
    is_normalized = zdir.normalize();
    assert(is_normalized);

    xdir = mi::math::cross(zdir, up_dir);
    is_normalized = xdir.normalize();
    assert(is_normalized);

    ydir = mi::math::cross(xdir, zdir);
    is_normalized = ydir.normalize();
    assert(is_normalized);

    nv::index_common::no_unused_variable_warning_please(is_normalized);
}

//----------------------------------------------------------------------
void Camera_tool::get_glu_lookat_vector(
    const nv::index::ICamera*         cam,
    mi::math::Vector<mi::Float32, 3>& eye,
    mi::math::Vector<mi::Float32, 3>& center,
    mi::math::Vector<mi::Float32, 3>& up)
{
    assert(cam != 0);

    const mi::math::Vector<mi::Float32, 3> view_dir = cam->get_view_direction();
    // This should be some vector with length.
    // The ideal test is assert(length(view_dir) == 1.0f);
    // but it is considered numerical problem.
    assert(length(view_dir) > 0.9f);

    eye    = cam->get_eye_point();
    center = eye + view_dir;
    up     = cam->get_up_direction();
}

//----------------------------------------------------------------------
static mi::math::Matrix<mi::Float64, 4, 4> get_orthographic_matrix(
    nv::index::IOrthographic_camera* ocam)
{
    assert(ocam != 0);

    const mi::Float64 film_width  = ocam->get_aperture();
    const mi::Float64 film_height = film_width / ocam->get_aspect();

    const mi::Float64 left    = -film_width / 2.0;
    const mi::Float64 right   =  film_width / 2.0;
    const mi::Float64 bottom  = -film_height / 2.0;
    const mi::Float64 top     =  film_height / 2.0;
    const mi::Float64 nearVal = ocam->get_clip_min();
    const mi::Float64 farVal  = ocam->get_clip_max();

    const mi::math::Vector<mi::Float64, 3> t(
        -(right + left) / (right - left),
        -(top + bottom) / (top - bottom),
        -(farVal + nearVal) / (farVal - nearVal));

    const mi::math::Matrix<mi::Float64, 4, 4> ortho_mat(
        2.0 / (right - left), 0.0,                  0.0,                       0.0,
        0.0,                  2.0 / (top - bottom), 0.0,                       0.0,
        0.0,                  0.0,                  -2.0 / (farVal - nearVal), 0.0,
        t.x,                  t.y,                  t.z,                       1.0);
    return ortho_mat;
}

//----------------------------------------------------------------------
static mi::math::Matrix<mi::Float64, 4, 4> get_perspective_matrix(
    nv::index::IPerspective_camera* pcam)
{
    assert(pcam != 0);

    const mi::Float64 fovy_rad = pcam->get_fov_y_rad();
    const mi::Float64 aspect   = pcam->get_aspect();
    const mi::Float64 znear    = pcam->get_clip_min();
    const mi::Float64 zfar     = pcam->get_clip_max();
    const mi::Float64 fr       = 1.0 / tan(fovy_rad / 2.0);

    const mi::Float64 clip_a = (zfar + znear) / (znear - zfar);
    const mi::Float64 clip_b = (2.0 * zfar * znear) / (znear - zfar);

    const mi::math::Matrix<mi::Float64, 4, 4> perspective_mat(
        (fr / aspect),     0.0,       0.0,        0.0,
        0.0,               fr,        0.0,        0.0,
        0.0,               0.0,       clip_a,    -1.0,
        0.0,               0.0,       clip_b,     0.0);
    return perspective_mat;
}

//----------------------------------------------------------------------
static mi::math::Matrix<mi::Float64, 4, 4> get_camera_matrix(
    nv::index::ICamera* cam)
{
    assert(cam != 0);

    mi::math::Vector<mi::Float32, 3> e = cam->get_eye_point();
    mi::math::Vector<mi::Float64, 3> eyepos(e.x, e.y, e.z);
    mi::math::Vector<mi::Float64, 3> up_direction(mi::math::Vector<mi::Float32, 3>(cam->get_up_direction()));
    
    // Camera orthonormal basis. gluLookat manner (right hand)
    mi::math::Vector<mi::Float64, 3> zdir = -mi::math::Vector<mi::Float64, 3>(mi::math::Vector<mi::Float32, 3>(cam->get_view_direction()));
    zdir.normalize();
    mi::math::Vector<mi::Float64, 3> xdir = mi::math::cross(up_direction, zdir);
    xdir.normalize();
    mi::math::Vector<mi::Float64, 3> ydir = mi::math::cross(zdir, xdir);
    ydir.normalize();

    mi::math::Matrix<mi::Float64, 4,4 > cam_mat(
        xdir.x, ydir.x, zdir.x, 0.0,
        xdir.y, ydir.y, zdir.y, 0.0,
        xdir.z, ydir.z, zdir.z, 0.0,
        -mi::math::dot(xdir, eyepos),
        -mi::math::dot(ydir, eyepos),
        -mi::math::dot(zdir, eyepos),
        1.0);
    return cam_mat;
}

//----------------------------------------------------------------------
static mi::math::Matrix<mi::Float64, 4, 4> get_projection_matrix(
    nv::index::ICamera* cam)
{
    assert(cam != 0);

    mi::base::Handle<nv::index::IPerspective_camera>  perspective_camera(
        cam->get_interface<nv::index::IPerspective_camera>());
    if (perspective_camera.is_valid_interface())
    {
        return get_perspective_matrix(perspective_camera.get());
    }

    mi::base::Handle<nv::index::IOrthographic_camera> ortho_camera(
        cam->get_interface<nv::index::IOrthographic_camera>());
    if (ortho_camera.is_valid_interface())
    {
        return get_orthographic_matrix(ortho_camera.get());
    }

    ERROR_LOG << "Unknown type of camera, return I matrix.";
    return mi::math::Matrix<mi::Float64, 4, 4>(1.0);
}

//----------------------------------------------------------------------
mi::math::Matrix<mi::Float64, 4, 4> Camera_tool::get_mv_matrix(
    nv::index::ICamera*                        cam,
    const mi::math::Matrix<mi::Float32, 4, 4>& scene_transform,
    const mi::math::Vector<mi::Sint32, 2>&     viewport_resolution)
{
    const mi::math::Matrix<mi::Float64, 4, 4> scene_mat(scene_transform);
    const mi::math::Matrix<mi::Float64, 4, 4> camera_mat = get_camera_matrix(cam);

    const mi::math::Matrix<mi::Float64, 4, 4> mv_mat = scene_mat * camera_mat;

    const mi::math::Matrix<mi::Float64, 4, 4> proj_mat = get_projection_matrix(cam);
    const mi::math::Matrix<mi::Float64, 4, 4> mvp_mat = mv_mat * proj_mat;

    static const mi::math::Matrix<mi::Float64, 4, 4> transform_to_raster_mat =
        mi::math::Matrix<mi::Float64, 4, 4> (
            1.0, 0.0, 0.0, 0.0,
            0.0, 1.0, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0,
            1.0, 1.0, 0.0, 1.0)  *
        mi::math::Matrix<mi::Float64, 4, 4>(
            0.5, 0.0, 0.0, 0.0,
            0.0, 0.5, 0.0, 0.0,
            0.0, 0.0, 1.0, 0.0,
            0.0, 0.0, 0.0, 1.0);

    assert(viewport_resolution.x > 0);
    assert(viewport_resolution.y > 0);

    const mi::math::Matrix<mi::Float64, 4, 4> scale_to_window_mat(
        static_cast<mi::Float64>(viewport_resolution.x), 0.0, 0.0, 0.0,
        0.0, static_cast<mi::Float64>(viewport_resolution.y), 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0);

    const mi::math::Matrix<mi::Float64, 4, 4> object_to_screen
        = mvp_mat * transform_to_raster_mat * scale_to_window_mat;

    return object_to_screen;
}

//----------------------------------------------------------------------
