/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "stereo_camera.h"

#include <cassert>
#include <sstream>

#include <nv/index/icamera.h>

#include "examiner_manipulator.h"
#include "nvindex_appdata.h"
#include "camera_utility.h"

//----------------------------------------------------------------------
Stereo_camera::Stereo_camera()
    :
    // m_cam_tag[]
    m_current_cam_idx(SC_Left),
    m_eye_separation_dist(1.0)
{
    m_cam_tag[SC_Left]  = mi::neuraylib::NULL_TAG;
    m_cam_tag[SC_Right] = mi::neuraylib::NULL_TAG;
}

//----------------------------------------------------------------------
Stereo_camera::Stereo_camera(mi::neuraylib::Tag left_cam_tag,
                             mi::neuraylib::Tag right_cam_tag)
    :
    // m_cam_tag[]
    m_current_cam_idx(SC_Left),
    m_eye_separation_dist(1.0)
{
    this->set_camera_tag(left_cam_tag, right_cam_tag);
}

//----------------------------------------------------------------------
void Stereo_camera::set_camera_tag(mi::neuraylib::Tag left_cam_tag,
                                   mi::neuraylib::Tag right_cam_tag)
{
    assert(left_cam_tag .is_valid());
    assert(right_cam_tag.is_valid());

    m_cam_tag[SC_Left]  = left_cam_tag;
    m_cam_tag[SC_Right] = right_cam_tag;
}

//----------------------------------------------------------------------
mi::neuraylib::Tag Stereo_camera::get_camera_tag(mi::Uint32 cidx) const
{
    assert(this->is_valid_index(cidx));

    return m_cam_tag[cidx];
}

//----------------------------------------------------------------------
void Stereo_camera::set_current_camera(mi::Uint32 cidx)
{
    assert(this->is_valid_index(cidx));
    m_current_cam_idx = cidx;
}

//----------------------------------------------------------------------
mi::neuraylib::Tag Stereo_camera::get_current_camera() const
{
    return m_cam_tag[m_current_cam_idx];
}

//----------------------------------------------------------------------
void Stereo_camera::swap_current_camera()
{
    if (m_current_cam_idx == SC_Left)
    {
        m_current_cam_idx = SC_Right;
    }
    else
    {
        m_current_cam_idx = SC_Left;
    }
}

//----------------------------------------------------------------------
void Stereo_camera::update_right_camera_and_swap_current_camera(
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    // always update even only left to see, since there is a direct
    // access method.
    this->update_right_camera(dice_transaction);
    this->swap_current_camera();
}

//----------------------------------------------------------------------
void Stereo_camera::set_eye_separation(mi::Float64 separation_dist)
{
    m_eye_separation_dist = separation_dist;
}


//----------------------------------------------------------------------
mi::Float64 Stereo_camera::get_eye_separation() const
{
    return m_eye_separation_dist;
}

//----------------------------------------------------------------------
void Stereo_camera::update_right_camera(
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    assert(dice_transaction != 0);
    assert(this->is_valid());

    mi::base::Handle< const nv::index::IPerspective_camera > left_cam(
        dice_transaction->access< nv::index::IPerspective_camera >(m_cam_tag[SC_Left]));
    mi::base::Handle< nv::index::IPerspective_camera > right_cam(
        dice_transaction->edit< nv::index::IPerspective_camera >(m_cam_tag[SC_Right]));
    assert(left_cam.is_valid_interface());
    assert(right_cam.is_valid_interface());

    // synchronize the camera
    right_cam->assign(left_cam.get());

    // get camera position parameter and compute camera coords basis
    // currently two cameras share the position
    mi::math::Vector<mi::Float32, 3> leye (left_cam->get_eye_point());
    mi::math::Vector<mi::Float32, 3> llook = 
        Nvindex_AppData::instance()->get_user_interaction(0)->
        get_examiner()->get_examiner_rotation_center();
    mi::math::Vector<mi::Float32, 3> lup  (left_cam->get_up_direction());
    mi::math::Vector<mi::Float32, 3> lview(left_cam->get_view_direction());
    mi::math::Vector<mi::Float32, 3> l_cam_xdir = mi::math::cross(lview, lup);
    l_cam_xdir.normalize();

    mi::math::Vector<mi::Float32, 3> reye = leye + (this->get_eye_separation() * l_cam_xdir);
    right_cam->set_eye_point(reye);

    mi::math::Vector<mi::Float32, 3> rview = llook - reye;
    rview.normalize();
    right_cam->set_view_direction(rview);    
    right_cam->set_up_direction(lup); // the same as the left up.
}

//----------------------------------------------------------------------
bool Stereo_camera::is_valid() const
{
    if ((m_cam_tag[SC_Left] .is_valid()) &&
        (m_cam_tag[SC_Right].is_valid()))
    {
        return true;
    }
    return false;
}

//----------------------------------------------------------------------
bool Stereo_camera::is_valid_index(mi::Uint32 cidx) const
{
    if ((cidx == SC_Left) || (cidx == SC_Right))
    {
        return true;
    }
    return false;
}

//----------------------------------------------------------------------
std::string Stereo_camera::to_string() const
{
    std::stringstream sstr;
    sstr << "camtag [L]: " << m_cam_tag[Stereo_camera::SC_Left].id
         << ", [R]: "      << m_cam_tag[Stereo_camera::SC_Right].id
         << ", current cam: "  << m_current_cam_idx
         << ", eye distance: " << m_eye_separation_dist;

    return sstr.str();
}

//----------------------------------------------------------------------
