/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief a stereo camera example (two cameras)

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_STEREO_CAMERA_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_STEREO_CAMERA_H

#include <mi/dice.h>

#include <vector>
#include <string>


/// a simple stereo camera
///
/// This camera has two ICameras, SC_Left and SC_Right. In this
/// implementation, the left camera is the main camera, this means you
/// move the left camera only, and the right camera parameters are
/// computed. (Another implementation is of course possible.)
///
class Stereo_camera
{
public:
    /// left or right camera
    enum Stereo_camera_LR_e {
        /// left camera
        SC_Left,
        /// right camera
        SC_Right,
        /// number of cameras in the stereo camera
        SC_Cemara_count
    };

public:
    /// default constructor. create an invelid stereo camera
    Stereo_camera();

    /// constructor with left and right camera
    /// \param[in] left_cam_tag  left camera tag
    /// \param[in] right_cam_tag right camera tag
    Stereo_camera(mi::neuraylib::Tag left_cam_tag,
                  mi::neuraylib::Tag right_cam_tag);

    /// set main camera tag
    ///
    /// \param[in] left_cam_tag  left  camera tag, must be valid camera tags
    /// \param[in] right_cam_tag right camera tag, must be valid camera tags
    void set_camera_tag(mi::neuraylib::Tag left_cam_tag,
                        mi::neuraylib::Tag right_cam_tag);

    /// get camera tag
    ///
    /// \param[in] cidx camera index. {SC_Left|SC_Right}
    /// \return main camera index, if no camera, NULL_TAG
    mi::neuraylib::Tag get_camera_tag(mi::Uint32 cidx = SC_Left) const;

    /// set current camara
    ///
    /// \param[in] cidx camera index. {SC_Left|SC_Right}
    void set_current_camera(mi::Uint32 cidx);

    /// get current camara tag
    ///
    /// \return current camera tag
    mi::neuraylib::Tag get_current_camera() const;

    /// swap current camera: left <-> right.
    void swap_current_camera();

    /// swap current camera: left <-> right.
    ///   - synchronize left and right camera parameters
    ///   - update right position
    void update_right_camera_and_swap_current_camera(
        mi::neuraylib::IDice_transaction* dice_transaction);

    /// set eye separation distance
    ///
    /// \param[in] separation_dist eye separation distance
    void set_eye_separation(mi::Float64 separation_dist);

    /// get eye separation distance
    ///
    /// \return get current eye separation distance
    mi::Float64 get_eye_separation() const;

    /// update the right camera state based on left camera state.
    ///
    /// \param[in] dice_transaction dice transaction to update the camera
    /// state
    void update_right_camera(
        mi::neuraylib::IDice_transaction* dice_transaction);

    /// the camera is valid?
    /// \return true when camera tags are not NULL_TAG.
    bool is_valid() const;

    /// get string representation of this object
    std::string to_string() const;

private:
    /// check the camera index validity
    /// \param[in] cidx left or right camera
    /// \return true when valid
    bool is_valid_index(mi::Uint32 cidx) const;

private:
    /// camera tags
    mi::neuraylib::Tag m_cam_tag[Stereo_camera::SC_Cemara_count];
    /// current camera index
    mi::Uint32         m_current_cam_idx;
    /// eye separation distance. No unit so far. This depends on the scene.
    mi::Float64        m_eye_separation_dist;
};

#endif // #ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_STEREO_CAMERA_H
