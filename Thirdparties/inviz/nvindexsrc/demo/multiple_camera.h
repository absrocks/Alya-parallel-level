/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief handling multiple (stereo) cameras

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_MULTIPLE_STEREO_CAMERA_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_MULTIPLE_STEREO_CAMERA_H

#include "stereo_camera.h"

#include <nv/index/icamera.h>

#include <vector>

/// A multiple stereo camera
///
/// This camera has:
///   - main camera. This is used for navigation.
///   - storage camera. save and load the view to these cameras.
///
/// Code Usage
/// \code
/// Multiple_stereo_camera multi_cam;
/// // set main camera
/// multi_cam.set_main_camera_tag(
///     session.get()->create_camera(dice_transaction));
/// // add three storage cameras
/// for(mi::Sint32 i = 0; i < 3; ++i)
/// {
///     multi_cam.append_camera_tag(
///               session.get()->create_camera(dice_transaction));
/// }
///
/// // save main camera to camera 1
/// current_camera_tag = multi_cam.save_camera_state(1, dice_transaction);
///
/// // restore camera 1 to main camera
/// multi_cam.load_camera_state(1, dice_transaction.get());
/// \endcode
class Multiple_stereo_camera
{
public:
    /// default constructor
    Multiple_stereo_camera();

    /// set main camera
    ///
    /// \param[in] cam main camera, must be valid camera
    void set_main_camera(const Stereo_camera & cam);

    /// get main camera (copy)
    ///
    /// \return main camera, if no camera, invalid camera returns
    Stereo_camera get_main_camera() const;

    /// peek main camera (reference)
    ///
    /// \return reference to the main camera, if no camera, reference
    /// to an invalid camera.
    Stereo_camera * peek_main_camera();

    /// get main camera's current camera tag. for convenient.
    ///
    /// current camera depends on stereo mode. It is Left camera or
    /// Right camera depends on the main camera's state.
    ///
    /// \return main camera's current camera tag.
    mi::neuraylib::Tag get_current_main_camera_tag() const;

    /// get main camera's base camera tag. for convenient.
    ///
    /// base camera is always the main left camera. Because all the
    /// manipuration and stereo camera computation is based on the left camera.
    ///
    /// \return main camera's base camera tag.
    mi::neuraylib::Tag get_main_base_camera_tag() const;

    /// update camera state and swap the view (left and right) of
    /// current main camera.  for stereo view.
    ///
    /// \param[in] dice_transaction dice transaction to update the ICamera instance
    void update_swap_main_camera_view(
        mi::neuraylib::IDice_transaction* dice_transaction);

    /// append storage camera
    ///
    /// \param[in] cam camera, must be valid camera
    void append_camera(const Stereo_camera & cam);

    /// get storage camera (copy)
    ///
    /// \param[in] cam_idx camera index
    /// \return    storage camera
    Stereo_camera get_storage_camera(mi::Sint32 cam_idx) const;

    /// peek storage camera (reference)
    ///
    /// \param[in] cam_idx camera index
    /// \return    reference to the storage camera
    Stereo_camera * peek_storage_camera(mi::Sint32 cam_idx);

    /// initialize all storage cameras by the main camera
    void initialize_storage_camera_by_main_camera(
        mi::neuraylib::IDice_transaction* dice_transaction);

    /// save main camera state to camera cam_idx
    ///
    /// \param[in] dest_cam_idx copy destination camera index
    /// \param[in] dice_transaction dice transaction to update the ICamera instance
    /// \return true when succeeded
    bool save_camera_state(
        mi::Sint32         dest_cam_idx,
        mi::neuraylib::IDice_transaction* dice_transaction);

    /// load main camera state from camera cam_idx
    ///
    /// \param[in] src_cam_idx copy source camera tag index
    /// \param[in] dice_trans dice transaction to update the ICamera instance
    /// \return true when succeeded
    bool load_camera_state(
        mi::Sint32         src_cam_idx,
        mi::neuraylib::IDice_transaction* dice_transaction);

    /// get number of cameras
    ///
    /// \return number of cameras stored
    mi::Sint32 get_camera_count() const;

    /// get the ortho camera's tag.
    /// \return the ortho camera's tag.
    mi::neuraylib::Tag get_ortho_camera_tag() const;
    /// get the ortho camera's tag.
    /// \return the ortho camera's tag.
    void set_ortho_camera_tag(const mi::neuraylib::Tag& t);

    /// query if the ortho camera is used.
    bool get_use_ortho_camera() const;
    /// set if ortho camera is used.
    void set_use_ortho_camera(bool u);

    /// get current camera tag considering ortho camera
    mi::neuraylib::Tag get_current_camera_tag_considering_ortho() const;

private:
    /// check the index validity
    /// \param[in] idx index of the camera tag stored
    /// \return true when valid
    bool is_valid_index(mi::Sint32 idx) const;

    /// copy stereo camera
    ///
    /// \param[in] p_dst_cam pointer to the destination stereo camera
    /// \param[in] p_src_cam pointer to the source stereo camera
    /// \param[in] dice_transaction  the transaction
    /// \return true when succeeded
    bool copy_stereo_camera(
        Stereo_camera * p_dst, Stereo_camera * p_src,
        mi::neuraylib::IDice_transaction* dice_transaction);

private:
    /// main camera tag
    Stereo_camera                m_main_cam;
    /// stereo camera vector
    std::vector< Stereo_camera > m_stereo_cam_vec;

    bool                         m_use_otho_camera;
    mi::neuraylib::Tag           m_ortho_camera_tag;
};

//----------------------------------------------------------------------
#endif // #ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_MULTIPLE_STEREO_CAMERA_H
