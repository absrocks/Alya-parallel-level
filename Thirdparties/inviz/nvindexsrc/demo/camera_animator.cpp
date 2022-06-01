/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "camera_animator.h"

#include <nv/index/icamera.h>

#include <cassert>
#include <iterator>

#include "common/forwarding_logger.h"
#include "common/tokenizer.h"
#include "common/string_dict.h"

#include "nvindex_appdata.h"

//======================================================================
// Camera_animator_data
//======================================================================
//----------------------------------------------------------------------
Camera_animator_data::Camera_animator_data()
    :
    m_eye   (0.0f, 0.0f, -1.0f),
    m_dir   (0.0f, 0.0f,  1.0f),
    m_up    (0.0f, 1.0f,  0.0f),
    m_aspect(1.0f),
    m_aperture(0.033f),
    m_focal(0.03),
    m_clip_min(0.01f),
    m_clip_max(10000.0f),
    m_nb_transit_frame(10)
{
    // empty
}

//----------------------------------------------------------------------
Camera_animator_data::Camera_animator_data(
    const mi::math::Vector< mi::Float32, 3 >& eye,
    const mi::math::Vector< mi::Float32, 3 >& dir,
    const mi::math::Vector< mi::Float32, 3 >& up,
    mi::Float64 aspect,
    mi::Float64 aperture,
    mi::Float64 focal,
    mi::Float64 clip_min,
    mi::Float64 clip_max,
    mi::Sint32  nb_transit_frame)
    :
    m_eye   (eye),
    m_dir   (dir),
    m_up    (up),
    m_aspect(aspect),
    m_aperture(aperture),
    m_focal(focal),
    m_clip_min(clip_min),
    m_clip_max(clip_max),
    m_nb_transit_frame(nb_transit_frame)
{
    // empty
}

//----------------------------------------------------------------------
std::string Camera_animator_data::to_string(mi::Size frame_number) const
{
    const std::string prefix_str = "app::camera_animator::key::";
    std::stringstream sstr;

    sstr << prefix_str << frame_number << "::from = "           << m_eye      << "\n"
         << prefix_str << frame_number << "::dir = "            << m_dir      << "\n"
         << prefix_str << frame_number << "::up = "             << m_up       << "\n"
         << prefix_str << frame_number << "::aspect = "         << m_aspect   << "\n"
         << prefix_str << frame_number << "::aperture = "       << m_aperture << "\n"
         << prefix_str << frame_number << "::focal = "          << m_focal    << "\n"
         << prefix_str << frame_number << "::clip_min = "       << m_clip_min << "\n"
         << prefix_str << frame_number << "::clip_max = "       << m_clip_max << "\n"
         << prefix_str << frame_number << "::transit_frames = " << m_nb_transit_frame;

    return sstr.str();
}

//----------------------------------------------------------------------
//======================================================================
// ICamera_animator
//======================================================================
//----------------------------------------------------------------------
Camera_animator::Camera_animator()
    :
    m_cur_data_idx(0),
    m_interp_step(50),      // each view to view has m_interp_step + 1 steps so far.
    m_cur_interp_step_idx(0)
{
    // empty
}

//----------------------------------------------------------------------
void Camera_animator::clear()
{
    m_cam_dat_vec.clear();
}

//----------------------------------------------------------------------
void Camera_animator::append(const Camera_animator_data& adat)
{
    m_cam_dat_vec.push_back(adat);
}

//----------------------------------------------------------------------
void Camera_animator::reset(bool is_reset)
{
    m_cur_interp_step_idx = 0;
    if(is_reset){
        m_cur_data_idx = 0;
    }
    else{
        m_cur_data_idx = static_cast< mi::Sint32 >(m_cam_dat_vec.size() + 1);
    }
}

//----------------------------------------------------------------------
bool Camera_animator::is_valid()
{
    if(m_cam_dat_vec.size() < 2){
        // INFO_LOG << "Camera_animator: m_cam_dat_vec(" << m_cam_dat_vec.size() << ") < 2";
        return false;
    }
    if(m_cur_data_idx + 1 >= static_cast< mi::Sint32 >(m_cam_dat_vec.size())){
        // INFO_LOG << "Camera_animator: cur index " << m_cur_data_idx
        //          << " + 1 >= m_cam_dat_vec(" << m_cam_dat_vec.size() << ")";
        return false;       // can not interpolate, no more data.
    }

    return true;
}

//----------------------------------------------------------------------
bool Camera_animator::next()
{
    if(!this->is_valid()){
        return false;
    }

    //  INFO_LOG << "DEBUG: curdata = " << m_cur_data_idx << "/" << m_cam_dat_vec.size()
    //           << ", step = "<< m_cur_interp_step_idx << "/" << m_interp_step;

    if(m_cur_interp_step_idx + 1 <= m_interp_step){
        // can step interval
        ++m_cur_interp_step_idx;
        return true;
    }
    else{
        // next data
        ++m_cur_data_idx;
        m_cur_interp_step_idx = 0;

        if(m_cur_data_idx + 1 >= static_cast< mi::Sint32 >(m_cam_dat_vec.size())){
            return false;       // can not interpolate, no more data.
        }
    }

    return true;
}

//----------------------------------------------------------------------
mi::Size Camera_animator::size() const
{
    return m_cam_dat_vec.size();
}

//----------------------------------------------------------------------
Camera_animator_data Camera_animator::get_ith_data(mi::Size idx) const
{
    if(idx >= this->size()){
        ERROR_LOG << "Camera_animator::get_ith_data: idx: " << idx << " is out of range.";
        return Camera_animator_data();
    }

    return m_cam_dat_vec.at(idx);
}

//----------------------------------------------------------------------
void Camera_animator::dump() const
{
    for(mi::Size i = 0; i < this->size(); ++i){
        INFO_LOG << this->get_ith_data(i).to_string(i);
    }
}

//----------------------------------------------------------------------
Camera_animator_data Camera_animator::get_current_interp_camera_parameter() const
{
    assert(m_cur_data_idx + 1 < static_cast< mi::Sint32 >(m_cam_dat_vec.size()));
    Camera_animator_data src = m_cam_dat_vec.at(m_cur_data_idx);
    Camera_animator_data dst = m_cam_dat_vec.at(m_cur_data_idx + 1);

    Camera_animator_data interp_res =
        this->interp(src, dst, m_interp_step, m_cur_interp_step_idx);
    return interp_res;
}

//----------------------------------------------------------------------
Camera_animator_data Camera_animator::interp(const Camera_animator_data & src,
                                             const Camera_animator_data & dst,
                                             mi::Sint32 max_idx,
                                             mi::Sint32 cur_idx) const
{
    assert(max_idx >  0);
    assert(cur_idx >= 0);
    assert(cur_idx <= max_idx);
    Camera_animator_data res;

    mi::Float64 alpha = (static_cast<mi::Float64>(cur_idx))/(static_cast<mi::Float64>(max_idx));

    res.m_eye = ((1.0 - alpha) * src.m_eye) + (alpha * dst.m_eye);
    res.m_dir = ((1.0 - alpha) * src.m_dir) + (alpha * dst.m_dir);
    if (!res.m_dir.normalize())
    {
        ERROR_LOG << "Linear interpolation is break down for camera direction. Set to source.";
        res.m_dir = src.m_dir;
    }

    mi::math::Vector< mi::Float32, 3 > up = ((1.0 - alpha) * src.m_up) + (alpha * dst.m_up);
    if (length(up) < 1.0e-6)
    {
        WARN_LOG << "Linear interpolation is break down... set up to Y up.";
        up.x = 0.0; up.y = 1.0; up.z = 0.0;
    }
    up.normalize();
    res.m_up = up;
    // INFO_LOG << "alpha = " << alpha;

    // FIXME: no other value interpolation: aperture, clip, ... (may have no meaning)

    return res;
}

//======================================================================

// Forward declaration
static bool get_camera_animator_from_project(
    const nv::index_common::String_dict & animate_prj_dict,
    Camera_animator & camera_anim);

Camera_animator get_default_camera_animator(
    const nv::index_common::String_dict* app_prj)
{
    Camera_animator anim;

    if (app_prj->is_defined("app::camera_animator::key::number"))
        get_camera_animator_from_project(*app_prj, anim);

    return anim;
}

//----------------------------------------------------------------------
/// update up camera by camera animator
///
/// \param[in]  camera_tag camera tag to update the camera parameter
/// \param[in]  cam_anim   reference to the camera animator 
/// \param[in]  dice_transaction dice transaction for the camera edit
void update_camera_by_animator(
    const mi::neuraylib::Tag & camera_tag,
    Camera_animator * cam_anim,
    mi::neuraylib::IDice_transaction * dice_transaction)
{
    assert(camera_tag.is_valid());
    assert(cam_anim != 0);
    assert(dice_transaction != 0);

    // when camera animator is enabled, start the animation
    if(cam_anim->is_valid()){
        Camera_animator_data cad = cam_anim->get_current_interp_camera_parameter();
        cam_anim->next();
        // if(!cam_anim->is_valid()){
        //     cam_anim->reset(); // loop the animation
        // }
        {
            { // handle perspective camera
                mi::base::Handle< nv::index::IPerspective_camera > cam(
                    dice_transaction->edit< nv::index::IPerspective_camera >(camera_tag));

                if (cam.is_valid_interface())
                {
                    cam->set_eye_point(     cad.m_eye);
                    cam->set_view_direction(cad.m_dir);
                    cam->set_up_direction(  cad.m_up);
                    cam->set_aspect(        cad.m_aspect);
                    cam->set_aperture(      cad.m_aperture);
                    cam->set_focal(         cad.m_focal);
                    cam->set_clip_min(      cad.m_clip_min);
                    cam->set_clip_max(      cad.m_clip_max);
                }
            }
            { // handle orthographic camera
                mi::base::Handle< nv::index::IOrthographic_camera > cam(
                    dice_transaction->edit< nv::index::IOrthographic_camera >(camera_tag));

                if (cam.is_valid_interface())
                {
                    cam->set_eye_point(     cad.m_eye);
                    cam->set_view_direction(cad.m_dir);
                    cam->set_up_direction(  cad.m_up);
                    cam->set_aspect(        cad.m_aspect);
                    cam->set_aperture(      cad.m_aperture);
                    cam->set_clip_min(      cad.m_clip_min);
                    cam->set_clip_max(      cad.m_clip_max);
                }
            }
        }
    }
}

//----------------------------------------------------------------------
/// load an animation camera project file and returns it as a project object.
///
/// \param[in[  animate_prj_fname animation camera project file
/// \param[out] animate_prj_dict  (output) string dict contains the project
/// \return true when load succeeded.
static bool load_animation_camera_project_file(const std::string & animate_prj_fname,
                                               nv::index_common::String_dict & animate_prj_dict)
{
    mi::Sint32 file_version = -1;
    std::string err_mes;
    bool const ret = get_string_dict_with_magic(animate_prj_fname,
                                                "animate_camera_project",
                                                animate_prj_dict, file_version, &err_mes);
    if(!ret){
        ERROR_LOG << "Cannot load animate camera project file [" << animate_prj_fname << "].";
        ERROR_LOG << err_mes;
        return false;
    }
    mi::Sint32 const app_version = 0;
    if(file_version != app_version){
        ERROR_LOG << "Wrong animate project file version. version = " << file_version
                  << ", instead of " << app_version;
        return false;
    }
    return true;
}

//----------------------------------------------------------------------
/// project to camera animator conversion function
///
/// \param[in]  animate_prj_dict project contains animation information
/// \param[out] camera_animator  camera animator
/// \return true when succeeded
static bool get_camera_animator_from_project(
    const nv::index_common::String_dict & animate_prj_dict,
    Camera_animator & camera_anim)
{
    const std::string fn = "get_camera_animator_from_project: ";
    const std::string prefix_str = "app::camera_animator::";

    if(!(animate_prj_dict.is_defined(prefix_str + "key::number"))){
        ERROR_LOG << fn << "no " + prefix_str + "key::number. Ignored.";
        return false;
    }

    const mi::Sint32 key_count = nv::index_common::get_sint32(animate_prj_dict.get(prefix_str + "key::number"));
    if (key_count <= 0)
    {
        ERROR_LOG << fn << "invalid key frame count. Ignored.";
        return false;
    }

    Camera_animator tmp_ca;
    for(mi::Sint32 i = 0; i < key_count; ++i){
        const std::string cur_prefix = prefix_str + "key::" + nv::index_common::to_string(i) + "::";
        std::vector< std::string > key_list;
        key_list.push_back(cur_prefix + "from");
        // key_list.push_back(cur_prefix + "to");
        key_list.push_back(cur_prefix + "up");
        key_list.push_back(cur_prefix + "aspect");
        key_list.push_back(cur_prefix + "aperture");
        key_list.push_back(cur_prefix + "focal");
        key_list.push_back(cur_prefix + "clip_min");
        key_list.push_back(cur_prefix + "clip_max");
        key_list.push_back(cur_prefix + "transit_frames");

        std::vector< std::string > undef_list;
        bool const is_keys_ok = is_all_keys_defined(animate_prj_dict, key_list, &undef_list);
        std::stringstream sstr;
        if(!is_keys_ok){
            std::copy(undef_list.begin(), undef_list.end(), std::ostream_iterator< std::string >(sstr, " "));
            ERROR_LOG << fn << "Undefined camera key entries:\n" << sstr.str();
            return false;
        }

        mi::math::Vector<mi::Float32, 3> from = nv::index_common::get_vec_float32_3(animate_prj_dict.get(cur_prefix + "from"));
        mi::math::Vector<mi::Float32, 3> dir  = nv::index_common::get_vec_float32_3(animate_prj_dict.get(cur_prefix + "dir", "0 0 1"));

        // Transition code: "to" is deprecated, use "dir": TODO
        if (animate_prj_dict.is_defined(cur_prefix + "to"))
        {
            WARN_LOG << "camera parameter: " << (cur_prefix + "to") << " is obsolete, please use " 
                     << (cur_prefix + "dir") << " instead.";
            const mi::math::Vector<mi::Float32, 3> to = nv::index_common::get_vec_float32_3(animate_prj_dict.get(cur_prefix + "to"));
            dir = from - to;
            INFO_LOG << "camera parameter: " << (cur_prefix + "dir") << " is calculated from 'from - to'  : " << dir;
        }

        mi::math::Vector<mi::Float32, 3> up   = nv::index_common::get_vec_float32_3(animate_prj_dict.get(cur_prefix + "up"));
        mi::Float64 aspect         = nv::index_common::get_float64(animate_prj_dict.get(cur_prefix + "aspect"));
        mi::Float64 aperture       = nv::index_common::get_float64(animate_prj_dict.get(cur_prefix + "aperture"));
        mi::Float64 focal          = nv::index_common::get_float64(animate_prj_dict.get(cur_prefix + "focal"));
        mi::Float64 clip_min       = nv::index_common::get_float64(animate_prj_dict.get(cur_prefix + "clip_min"));
        mi::Float64 clip_max       = nv::index_common::get_float64(animate_prj_dict.get(cur_prefix + "clip_max"));
        mi::Sint32  transit_frames = nv::index_common::get_sint32 (animate_prj_dict.get(cur_prefix + "transit_frames"));
        if(clip_min >= clip_max){
            ERROR_LOG << "inconsistent clip_min and clip_max [" << clip_min << clip_max << "]";
            return false;
        }
        // app::camera_animator::key::0::from = 0.703629493713379 -0.0909369438886642 3.95295476913452
        // app::camera_animator::key::0::to = 0.0514815151691437 0.011218941770494 0.812628209590912
        // app::camera_animator::key::0::up = -0.0160848833620548 0.999264001846313 0.034823652356863
        // app::camera_animator::key::0::aspect = 1
        // app::camera_animator::key::0::aperture = 0.0329999998211861
        // app::camera_animator::key::0::focal = 0.0299999993294477
        // app::camera_animator::key::0::clip_min = 0.321217000484467
        // app::camera_animator::key::0::clip_max = 32.1217002868652
        // app::camera_animator::key::0::transit_frames = 100
        
        tmp_ca.append(Camera_animator_data(from, dir, up, 
                                           aspect, aperture, focal,
                                           clip_min, clip_max, transit_frames));
    }

    camera_anim = tmp_ca;
    
    return true;
}

//----------------------------------------------------------------------
bool animate_camera_by_command_str(const std::string & command, std::string & ret_mes)
{
    // parse the command, get the filename
    std::string separators = " ";
    std::vector< std::string > tokens;
    nv::index_common::Tokenizer::parse(command, separators, tokens);
    
    const mi::Size token_size = tokens.size();
    if((token_size <= 1) || (tokens[0] != "animate_camera")){
        ret_mes = "Invalid command: animate_camera [animate.prj|?]";
        return false;
    }

    // help
    if(tokens[1] == "?"){
        INFO_LOG << "animate_camera [animate.prj|?]\n"
                 << "  animate.prj ... project file for animatiing camera.\n" 
                 << "  ?           ... This help message." ;
        ret_mes = "Usage: animate_camera [animate.prj|?]";
        return false;
    }

    // load the file
    const std::string animate_file_str = tokens[1];
    nv::index_common::String_dict animate_prj_dict;
    bool is_success = load_animation_camera_project_file(animate_file_str, 
                                                         animate_prj_dict);
    if(!is_success){
        ERROR_LOG << "failed to load camera animation file [" << animate_file_str << "]";
        return false;
    }
    animate_prj_dict.write(std::cout, ""); // debug


    // convert a string_dict to camera_animator.
    Camera_animator ca;
    is_success = get_camera_animator_from_project(animate_prj_dict, ca);
    if(!is_success){
        ERROR_LOG << "failed to construct a camera animator object from [" << animate_file_str << "]";
        return false;
    }

    INFO_LOG << "Load camera animation project [" << animate_file_str 
             << "], size: " << ca.size() << ", " 
             << (ca.is_valid() ? "valid" : "invalid");
    ca.dump();                  // DEBUG

    Nvindex_AppData::instance()->set_camera_animator_data(ca);
    Nvindex_AppData::instance()->peek_camera_animator_data()->reset(true);

    return true;
}

//----------------------------------------------------------------------
