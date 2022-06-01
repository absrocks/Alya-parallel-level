/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief camera animator for a dynamic performance test

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_DEBUG_VIEWER_CAMERA_ANIMATOR_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_DEBUG_VIEWER_CAMERA_ANIMATOR_H

#include <mi/dice.h>

#include "common/string_dict.h"

#include <string>
#include <vector>

/// camera animator's data.
/// All camera data except the resolution. (They will be deprecated anyway.)
class Camera_animator_data
{
public:
    /// default constructor
    Camera_animator_data();

    /// constructor
    ///
    /// \param[in] eye              camera eye position
    /// \param[in] dir              camera direction vector
    /// \param[in] up               camera up vector
    /// \param[in] aspect           aspect ratio  
    /// \param[in] aperture         aperture
    /// \param[in] focal            focal length
    /// \param[in] clip_min clip    plane min distance
    /// \param[in] clip_max clip    plane max distance
    /// \param[in] nb_transit_frame number of transit frames from the current to the next
    Camera_animator_data(
        const mi::math::Vector< mi::Float32, 3 >& eye,
        const mi::math::Vector< mi::Float32, 3 >& dir,
        const mi::math::Vector< mi::Float32, 3 >& up,
        mi::Float64 aspect,
        mi::Float64 aperture,
        mi::Float64 focal,
        mi::Float64 clip_min,
        mi::Float64 clip_max,
        mi::Sint32  nb_transit_frame);

    /// get the object string representation
    /// \param[in] frame_number frame number 
    std::string to_string(mi::Size frame_number) const;

public:
    /// camera eye position
    mi::math::Vector< mi::Float32, 3 > m_eye;
    /// camera direction vector
    mi::math::Vector< mi::Float32, 3 > m_dir;
    /// camera up vector
    mi::math::Vector< mi::Float32, 3 > m_up;
    /// aspect ratio
    mi::Float64                        m_aspect;
    /// aperture
    mi::Float64                        m_aperture;
    /// focal length
    mi::Float64                        m_focal;
    /// clip plane min distance
    mi::Float64                        m_clip_min;
    /// clip plane max distance
    mi::Float64                        m_clip_max;
    /// number of transit frames
    mi::Sint32                         m_nb_transit_frame;
};

/// camera animator
///
/// This has a list of camera parameter and returns the interpolated
/// parameters each get_current_interp_camera_parameter() call.
class Camera_animator
{
public:
    /// constructor
    Camera_animator();

    /// clear the data
    void clear();

    /// append the Camera_animator_data
    /// \param adat a camera animation data
    void append(const Camera_animator_data & adat);

    /// reset the interpolation status
    /// \param[in] is_reset when true, reset the state, false, make
    /// the state invalid.
    void reset(bool is_reset);

    /// current interpolation is valid?
    /// \return if no more next, return false.
    bool is_valid();

    /// step forward the interpolation
    /// \return if no more next, return false.
    bool next();

    /// get current interpolated camera parameter
    Camera_animator_data get_current_interp_camera_parameter() const;

    /// get the size of the animator data
    /// \return camera animator vector size.
    mi::Size size() const;

    /// i-th Camera_animator_data
    /// \param[in] idx data index
    /// \return i-th Camera_animator_data. Default object when idx is out of range.
    Camera_animator_data get_ith_data(mi::Size idx) const;

    /// dump the contents to INFO_LOG
    void dump() const;

private:
    /// just a simple linear interpolation
    Camera_animator_data interp(const Camera_animator_data & src,
                                const Camera_animator_data & dst,
                                mi::Sint32 max_idx,
                                mi::Sint32 cur_idx) const;

private:
    /// data vector
    std::vector< Camera_animator_data > m_cam_dat_vec;
    /// current data index
    mi::Sint32 m_cur_data_idx;
    /// interpolation step
    mi::Sint32 m_interp_step;
    /// interpolation step index
    mi::Sint32 m_cur_interp_step_idx;
};

//----------------------------------------------------------------------
/// get Camera_animator
/// \return a default Camera_animator instance.
Camera_animator get_default_camera_animator(
    const nv::index_common::String_dict* app_prj);

//----------------------------------------------------------------------
/// Update a camera by a camera animator
///
/// \param[in]  camera_tag camera tag to update the camera parameter
/// \param[in]  cam_anim   reference to the camera animator 
/// \param[in]  dice_transaction dice transaction for the camera edit
void update_camera_by_animator(
    const mi::neuraylib::Tag& camera_tag,
    Camera_animator* cam_anim,
    mi::neuraylib::IDice_transaction* dice_transaction);

//----------------------------------------------------------------------
/// Update a camera by a camera animator
///
/// command: animate_camera animation_project.prj
///
/// \param[in]   command    command string
/// \param[out]  ret_str    return result string 
/// \return true when success
bool animate_camera_by_command_str(const std::string & command, std::string & ret_str);

//----------------------------------------------------------------------

#endif // #ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_DEBUG_VIEWER_CAMERA_ANIMATOR_H
