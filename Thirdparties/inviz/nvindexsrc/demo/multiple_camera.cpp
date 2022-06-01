/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "multiple_camera.h"

#include <cassert>

#include "common/forwarding_logger.h"

//----------------------------------------------------------------------
Multiple_stereo_camera::Multiple_stereo_camera()
    :
    m_main_cam(),
    m_use_otho_camera(false),
    m_ortho_camera_tag(mi::neuraylib::NULL_TAG)
{
    // empty
}

//----------------------------------------------------------------------
void Multiple_stereo_camera::set_main_camera(const Stereo_camera & cam)
{
    assert(cam.is_valid());
    m_main_cam = cam;
}

//----------------------------------------------------------------------
Stereo_camera Multiple_stereo_camera::get_main_camera() const
{
    return m_main_cam;
}

//----------------------------------------------------------------------
Stereo_camera * Multiple_stereo_camera::peek_main_camera()
{
    return &m_main_cam;
}

//----------------------------------------------------------------------
mi::neuraylib::Tag Multiple_stereo_camera::get_current_main_camera_tag() const
{
    const mi::neuraylib::Tag current_cam_tag = m_main_cam.get_current_camera();

    return current_cam_tag;
}

//----------------------------------------------------------------------
mi::neuraylib::Tag Multiple_stereo_camera::get_main_base_camera_tag() const
{
    const mi::neuraylib::Tag base_cam_tag = m_main_cam.get_camera_tag(Stereo_camera::SC_Left);

    return base_cam_tag;
}

//----------------------------------------------------------------------
void Multiple_stereo_camera::update_swap_main_camera_view(
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    m_main_cam.update_right_camera_and_swap_current_camera(dice_transaction);
}

//----------------------------------------------------------------------
void Multiple_stereo_camera::append_camera(const Stereo_camera & cam)
{
    assert(cam.is_valid());
    m_stereo_cam_vec.push_back(cam);
}

//----------------------------------------------------------------------
Stereo_camera Multiple_stereo_camera::get_storage_camera(mi::Sint32 cam_idx) const
{
    if(!this->is_valid_index(cam_idx)){
        return Stereo_camera();
    }
    Stereo_camera ret_cam = m_stereo_cam_vec.at(cam_idx);
    assert(ret_cam.is_valid());

    return ret_cam;
}

//----------------------------------------------------------------------
Stereo_camera * Multiple_stereo_camera::peek_storage_camera(mi::Sint32 cam_idx)
{
    if(!this->is_valid_index(cam_idx)){
        return 0;
    }
    Stereo_camera * p_ret_cam = &(m_stereo_cam_vec.at(cam_idx));
    assert(p_ret_cam->is_valid());

    return p_ret_cam;
}

//----------------------------------------------------------------------
void Multiple_stereo_camera::initialize_storage_camera_by_main_camera(
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    if(!m_main_cam.is_valid()){
        ERROR_LOG << "Multiple_stereo_camera::initialize_storage_camera_by_main_camera: "
                  << "invalid main camera. cannot init.";
        return;
    }

    Stereo_camera * const p_maincam = this->peek_main_camera();
    for (size_t i = 0; i < m_stereo_cam_vec.size(); ++i)
    {
        Stereo_camera * p_dstcam  = this->peek_storage_camera(i);
        this->copy_stereo_camera(p_dstcam, p_maincam, dice_transaction);
    }
}


//----------------------------------------------------------------------
bool Multiple_stereo_camera::save_camera_state(
    mi::Sint32 dest_cam_idx,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    if(!(this->get_main_camera().is_valid())){
        ERROR_LOG << "Multiple_stereo_camera::save_camera_state: no main camera.";
        return false;
    }
    if(!this->is_valid_index(dest_cam_idx)){
        return false;
    }

    Stereo_camera * p_maincam = this->peek_main_camera();
    Stereo_camera * p_dstcam  = this->peek_storage_camera(dest_cam_idx);
    const bool ret = this->copy_stereo_camera(p_dstcam, p_maincam, dice_transaction);

    return ret;
}

//----------------------------------------------------------------------
bool Multiple_stereo_camera::load_camera_state(
    mi::Sint32         src_cam_idx,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    if(!(this->get_main_camera().is_valid())){
        ERROR_LOG << "Multiple_stereo_camera::laod_camera_state: no main camera.";
        return false;
    }
    if(!this->is_valid_index(src_cam_idx)){
        return false;
    }

    Stereo_camera * p_srccam  = this->peek_storage_camera(src_cam_idx);
    Stereo_camera * p_maincam = this->peek_main_camera();
    const bool ret = this->copy_stereo_camera(p_maincam, p_srccam, dice_transaction);

    return ret;
}

//----------------------------------------------------------------------
mi::Sint32 Multiple_stereo_camera::get_camera_count() const
{
    return static_cast< mi::Sint32 >(m_stereo_cam_vec.size());
}

//----------------------------------------------------------------------
bool Multiple_stereo_camera::is_valid_index(mi::Sint32 idx) const
{
    if((idx < 0) || (idx >= this->get_camera_count())){
        ERROR_LOG << "Multiple_stereo_camera::is_valid_index: idx (" << idx << ") is not in range. [0, "
                  << this->get_camera_count() << ")";
        return false;
    }

    return true;
}

//----------------------------------------------------------------------
bool Multiple_stereo_camera::copy_stereo_camera(
    Stereo_camera * p_dst, Stereo_camera * p_src,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    if((p_src == 0) || (!p_src->is_valid())){
        ERROR_LOG << "Multiple_stereo_camera::copy_stereo_camera: invalid copy source camera. ignored.";
        return false;
    }
    if((p_dst == 0) || (!p_dst->is_valid())){
        ERROR_LOG << "Multiple_stereo_camera::copy_stereo_camera: invalid copy destination camera. ignored.";
        return false;
    }
    if(p_dst == p_src){
        WARN_LOG << "Multiple_stereo_camera::copy_stereo_camera: identical copy source and destination. ignored.";
        return false;
    }

    mi::base::Handle< const nv::index::IPerspective_camera > s_l_cam(
        dice_transaction->access< nv::index::IPerspective_camera >(
            p_src->get_camera_tag(Stereo_camera::SC_Left)));
    mi::base::Handle< const nv::index::IPerspective_camera > s_r_cam(
        dice_transaction->access< nv::index::IPerspective_camera >(
            p_src->get_camera_tag(Stereo_camera::SC_Right)));
    assert(s_l_cam.is_valid_interface());
    assert(s_r_cam.is_valid_interface());

    mi::base::Handle< nv::index::IPerspective_camera > d_l_cam(
        dice_transaction->edit< nv::index::IPerspective_camera >(
            p_dst->get_camera_tag(Stereo_camera::SC_Left)));
    mi::base::Handle< nv::index::IPerspective_camera > d_r_cam(
        dice_transaction->edit< nv::index::IPerspective_camera >(
            p_dst->get_camera_tag(Stereo_camera::SC_Right)));
    assert(d_l_cam.is_valid_interface());
    assert(d_r_cam.is_valid_interface());

    d_l_cam->assign(s_l_cam.get());
    d_r_cam->assign(s_r_cam.get());

    return true;
}

mi::neuraylib::Tag Multiple_stereo_camera::get_ortho_camera_tag() const
{
    return m_ortho_camera_tag;
}

void Multiple_stereo_camera::set_ortho_camera_tag(const mi::neuraylib::Tag& t)
{
    assert(t.is_valid());
    m_ortho_camera_tag = t;
}

bool Multiple_stereo_camera::get_use_ortho_camera() const
{
    return m_use_otho_camera;
}

void Multiple_stereo_camera::set_use_ortho_camera(bool u)
{
    m_use_otho_camera = u;
}

mi::neuraylib::Tag Multiple_stereo_camera::get_current_camera_tag_considering_ortho() const
{
    if (get_use_ortho_camera())
    {
        return get_ortho_camera_tag();
    }
    
    return get_current_main_camera_tag();
}

//----------------------------------------------------------------------
