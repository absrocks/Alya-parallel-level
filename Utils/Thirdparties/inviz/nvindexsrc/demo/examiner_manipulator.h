/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief examiner manipulator (camera movement by mouse/keyboard input)

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_EXAMINER_MANIPULATOR_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_EXAMINER_MANIPULATOR_H

#include <mi/math/vector.h>
#include <mi/math/matrix.h>
#include <mi/base/types.h>

#include <cassert>
#include <string>
#include <vector>
#include <nv/index/icamera.h>
#include <nv/index/isession.h>
#include <nv/index/iscene.h>

#include "camera_utility.h"
#include "multiple_camera.h"

//----------------------------------------------------------------------

/// IndeX scene examiner manipulator.
///
/// This is so called examiner type scene manipulator. There are three
/// modes of this manipulator.
///
/// - trackball rotation mode: You can fly around the rotation center point.
/// - pan mode: camera pan movement
/// - zoom mode: camera zoom movement
///
class Examiner_manipulator
{
public:
    /// manipulation modes
    enum Manipulation_mode {
        /// trackball rotation (examine) mode
        EM_TrackballRotation,
        /// pan mode
        EM_Pan,
        /// zoom mode
        EM_Zoom,
        /// sentinel
        EM_Manipulation_mode_count
    };

    /// Mouse button
    enum Mouse_button {
        /// Left mouse button
        Left_button   = 0,
        /// Middle mouse button
        Middle_button = 1,
        /// Right mouse button
        Right_button  = 2,
        /// sentinel
        Mouse_button_count
    };

    /// Mouse action
    enum Mouse_action {
        /// Mouse button action down
        Button_down,
        /// Mouse button action up
        Button_up,
        /// sentinel
        Mouse_action_count
    };

public:
    /// Constructor
    Examiner_manipulator();

    /// initialize by dict (convenient method)
    ///
    /// \param[in] dict initi this with app::examiner::* key and its values
    void initialize_by_dict(const nv::index_common::String_dict& dict);

    /// set manipulation mode
    /// \param[in] manip_mode
    void set_manipulation_mode(Examiner_manipulator::Manipulation_mode manip_mode);
    /// get manipulation mode
    /// \return current manipulation mode
    Examiner_manipulator::Manipulation_mode get_manipulation_mode() const;

    /// Set camera pan mode screen match
    /// \param[in] is_match_screen true when screen match mode
    void set_pan_screen_match_mode(bool is_match_screen);
    /// Get camera pan mode screen match
    /// \return true when screen match mode
    bool is_pan_screen_match_mode() const;

    /// set pick navigation behavior
    /// \param[in] is_pick_to_lookat when true, camera will dive in
    /// the scene when zoom. Otherwise the loookat point is fixed.
    void set_navigation_pick_to_lookat(bool is_pick_to_lookat);
    /// get pick navigation behavior
    /// \return true when in the 'pick to lookat' mode
    bool is_navigation_pick_to_lookat() const;

    /// set main window_resolution
    /// \param[in] win_resolution main window resolution in pixel
    void set_main_window_resolution(const mi::math::Vector<mi::Sint32, 2>& win_resolution);
    /// get main window resolution
    /// \return current examiner window resolution
    mi::math::Vector<mi::Sint32, 2> get_main_window_resolution() const;
    /// get window_resolution uint32
    /// \return current examiner window resolution
    mi::math::Vector<mi::Uint32, 2> get_main_window_resolution_uint32() const;

    /// set scene bounding box (for navigation: pan and zoom distance computation)
    ///
    /// \param[in] scene_bbox scene bounding box
    void set_scene_bbox(const mi::math::Bbox<mi::Float32, 3>& scene_bbox);
    /// get scene bounding box
    ///
    /// \return current scene bounding box
    mi::math::Bbox<mi::Float32, 3> get_scene_bbox() const;

    /// Set examiner rotation center.
    ///
    /// This is the center of the camera rotation in world coordinate
    /// when the examiner is in the trackball mode.
    ///
    /// \param[in] rot_center examiner rotation center
    void set_examiner_rotation_center(
        const mi::math::Vector<mi::Float32, 3>& rot_center);
    /// Get examiner rotation center.
    mi::math::Vector<mi::Float32, 3> get_examiner_rotation_center() const;

    /// Set pan reference point
    ///
    /// This is the pan movement reference point.
    ///
    /// \param[in] pan_ref  pan reference point in object space
    void set_pan_reference_point_object_space(
        const mi::math::Vector<mi::Float32, 3>& pan_ref);
    /// Get pan reference point in object space
    mi::math::Vector<mi::Float32, 3> get_pan_reference_point_object_space() const;

    /// Set pan reference point validity
    ///
    /// \param[in] is_valid true when valid
    void set_pan_reference_point_valid(bool is_valid);
    /// Get pan reference point
    bool is_pan_reference_point_valid() const;

    /// Set reference point zoom factor
    /// \param[in] zoom_factor zoom factor, must be > 1.0
    void set_reference_point_zoom_factor(mi::Float32 zoom_factor);
    /// Get reference point zoom factor
    mi::Float32 get_reference_point_zoom_factor() const;

public:
    /// mouse press event handling
    ///
    /// \param[in] button        Which button was pressed
    /// \param[in] mouse_action  mouse action, {Button_press, Button_release}
    /// \param[in] mouse_x mouse position screen coordinate x
    /// \param[in] mouse_y mouse position screen coordinate y
    /// \param[in] camera_tag camera tag
    /// \param[in] dice_transaction dice transaction to change the ICamera instance
    void mouse_press_event(mi::Sint32                        button,
                           Mouse_action                      mouse_action,
                           mi::Sint32                        mouse_x,
                           mi::Sint32                        mouse_y,
                           const mi::neuraylib::Tag&         camera_tag,
                           mi::neuraylib::IDice_transaction* dice_transaction);

    /// mouse motion event handling
    ///
    /// \param[in] button        Which button was pressed
    /// \param[in] mouse_x mouse position screen coordinate x
    /// \param[in] mouse_y mouse position screen coordinate y
    /// \param[in] camera_tag camera tag
    /// \param[in] scene_tag scene for global transform access
    /// \param[in] dice_transaction dice transaction to change the ICamera instance
    void mouse_motion_event(mi::Sint32                        button,
                            mi::Sint32                        mouse_x,
                            mi::Sint32                        mouse_y,
                            const mi::neuraylib::Tag&         camera_tag,
                            const mi::neuraylib::Tag&         scene_tag,
                            mi::neuraylib::IDice_transaction* dice_transaction);

    /// mouse wheel event handling
    ///
    /// \param[in] wheel_delta mouse wheel delta (+1 or -1)
    /// \param[in] camera_tag camera tag
    /// \param[in] scene_tag scene for global transform access
    /// \param[in] dice_transaction dice transaction to change the ICamera instance
    void mouse_wheel_event(mi::Sint32                        wheel_delta,
                           const mi::neuraylib::Tag&         camera_tag,
                           const mi::neuraylib::Tag&         scene_tag,
                           mi::neuraylib::IDice_transaction* dice_transaction);



    /// set last mouse 2d position
    /// \param[in] mouse_x last x mouse position
    /// \param[in] mouse_y last y mouse position
    void set_last_mouse2d(mi::Sint32 mouse_x, mi::Sint32 mouse_y);

    /// set last 3d position on a sphere (trackball movement)
    /// \param[in] point_3d last position point.
    void set_last_mouse3d(const mi::math::Vector<mi::Float32, 3>& point_3d);

    /// Direct camera zoom
    void camera_zoom(
        mi::Float32                       zoom_factor,
        const mi::neuraylib::Tag&         camera_tag,
        const mi::neuraylib::Tag&         scene_tag,
        mi::neuraylib::IDice_transaction* dice_transaction);

    /// camera zoom to the reference point (to the pan point)
    ///
    /// \param[in] is_zoom_in zoom in when true, zoom out otherwise
    /// \param[in] scene      scene object
    /// \param[in] cam        camera object
    void camera_zoom_to_pan_reference_point(
        bool                              is_zoom_in,
        const nv::index::IScene*          scene,
        nv::index::ICamera*               cam);

public:
    // predefined camera parameter interface

    /// clear predefined camera parameters
    void clear_predefined_camera();
    /// append a predefined camera parameter
    void append_predefined_camera(Camera_parameter& cam_param);
    /// get current number of camera parameters
    mi::Size nb_predefined_camera_parameter() const;
    /// get camera parameter by index
    Camera_parameter get_predefined_camera_parameter(mi::Size idx);
    /// get predefined view name string vector
    std::vector<std::string> get_predefined_view_name_vec() const;

    /// set current predefined view index
    /// \param[in] idx index of current predefined view
    void set_predefined_view_index(mi::Size idx);
    /// get current predefined view index
    /// \return index of current predefined view
    mi::Size get_predefined_view_index() const;

    /// Set the camera parameter to the current predefined view
    /// index's camera parameter
    ///
    /// \param[in]  session          current session (for view_all)
    /// \param[in]  scene            current scene   (for view_all)
    /// \param[in]  ms_camera        set parameters to this multiple stereo camera
    /// \param[in]  dice_transaction dice db transaction
    void set_predefined_view_parameter_to_camera(
        const nv::index::ISession*        session,
        const nv::index::IScene*          scene,
        Multiple_stereo_camera*           ms_camera,
        mi::neuraylib::IDice_transaction* dice_transaction);

    /// set a pan movement heuristics by pick. use last valid pick point
    /// \param[in] is_enabled true when this heuristics enabled.
    void set_pan_reference_point_use_pick_point(bool is_enabled);
    /// get a pan movement heuristics by pick.
    /// \return true when this heuristics enabled.
    bool is_pan_reference_point_use_pick_point() const;

    /// set a pan movement heuristics by bisect point.
    /// \param[in] is_enabled true when this heuristics enabled.
    void set_pan_reference_point_roi_bisect_and_solid(bool is_enabled);
    /// get a pan movement heuristics by bisect point.
    /// \return true when this heuristics enabled.
    bool is_pan_reference_point_roi_bisect_and_solid() const;

private:
    /// Initialize predefined views.
    /// A subroutine of initialize_by_dict
    ///
    /// \param[in] dict option that has predefined view parameters
    void initialize_predefined_view_by_dict(const nv::index_common::String_dict& dict);


    /// map the mouse position on an unit sphere.
    ///
    /// \param[in]  mouse_pos mouse position in screen in pixel
    /// \param[out] p_3d_pos mapped position on a unit sphere
    /// \return truen when the p_3d_pos is valid
    bool map_to_sphere(const mi::math::Vector<mi::Sint32, 2>& mouse_pos,
                       mi::math::Vector<mi::Float32, 3>*      p_3d_pos);

    /// move camera as trackball rotation mode
    ///
    /// \param[in]  new_point_3d updated point in 3d
    /// \param[in]  camera_tag   camera tag
    /// \param[in]  dice_transaction   dice transaction to update the ICamera instance
    void move_camera_trackballrotate_mode(
        const mi::math::Vector<mi::Float32, 3>& new_point_3d,
        const mi::neuraylib::Tag&               camera_tag,
        mi::neuraylib::IDice_transaction*       dice_transaction);

    /// camera pan
    ///
    /// \param[in]  new_point_3d updated point in 3d
    /// \param[in]  camera_tag   camera tag
    /// \param[in]  scene_tag    scene tag for global transform access
    /// \param[in]  dice_transaction   dice transaction to update the ICamera instance
    void move_camera_pan_mode(
        const mi::math::Vector<mi::Sint32, 2>& new_point_2d,
        const mi::neuraylib::Tag&              camera_tag,
        const mi::neuraylib::Tag&              scene_tag,
        mi::neuraylib::IDice_transaction*      dice_transaction);

    /// camera zoom
    ///
    /// \param[in] new_point_3d updated point in 3d
    /// \param[in] camera_tag   camera tag
    /// \param[in] scene_tag    scene tag for global transform access
    /// \param[in] dice_transaction   dice transaction to update the ICamera instance
    void move_camera_zoom_mode(
        const mi::math::Vector<mi::Sint32, 2>& new_point_2d,
        const mi::neuraylib::Tag&              camera_tag,
        const mi::neuraylib::Tag&              scene_tag,
        mi::neuraylib::IDice_transaction*      dice_transaction);

    /// clear last mouse and coordinate information. Invalidation.
    void clear_last_mouse_info();

    /// get normalized screen coordinates
    ///
    /// normalized screen coordinate is norm_x = [-0.5, 0.5], norm_y = [-0.5, 0.5],
    /// (0, 0) is the center, Y is up direction.
    ///
    ///   +-----+-----+ (0.5, 0.5)
    ///   +     |     +
    ///   +     |     +
    ///   +-----+-----+
    ///   +     |(0,0)+
    ///   +     |     +
    ///   +-----+-----+
    ///(-0.5,-0.5)
    ///
    ///
    /// \param[in]  scr_x_i  screen coordinate x
    /// \param[in]  scr_y_i  screen coordinate y
    /// \param[out] norm_x   normalized screen coordinate x
    /// \param[out] norm_y   normalized screen coordinate y
    /// \return true when the coordinate is valid.
    bool get_normalized_screen_coords(mi::Sint32   scr_x_i,
                                      mi::Sint32   scr_y_i,
                                      mi::Float64& norm_x,
                                      mi::Float64& norm_y);

    /// View the whole boundaing box. Entry point for both perspective
    /// camera and orthographic camera.
    ///
    /// \param[in] world_bbox bounding box in the world space to see
    ///            (transformed)
    /// \param[in] viewdir   camera viewing direction to be set
    /// \param[in] up        camera up vector to be set
    /// \param[in] cam camera object to be updated
    void view_bbox_all(const mi::math::Bbox<mi::Float32, 3>&   world_bbox,
                       const mi::math::Vector<mi::Float32, 3>& viewdir,
                       const mi::math::Vector<mi::Float32, 3>& up,
                       nv::index::ICamera*                     cam);

    /// View the whole boundaing box for perspective camera
    ///
    /// \param[in]  world_bbox bounding box in the world space to see
    ///             (transformed)
    /// \param[in]  viewdir   camera viewing direction to be set
    /// \param[in]  up        camera up vector to be set
    /// \param[out] pcam      perspective camera object to be updated
    void view_bbox_all_perspective(
        const mi::math::Bbox<mi::Float32, 3>&   world_bbox,
        const mi::math::Vector<mi::Float32, 3>& viewdir,
        const mi::math::Vector<mi::Float32, 3>& up,
        nv::index::IPerspective_camera*         pcam);

    /// View the whole boundaing box for orthographic camera
    ///
    /// \param[in]  world_bbox bounding box in the world space to see
    ///             (transformed)
    /// \param[in]  viewdir   camera viewing direction to be set
    /// \param[in]  up        camera up vector to be set
    /// \param[out] ocam      orthographic camera object to be updated
    void view_bbox_all_orthographic(
        const mi::math::Bbox<mi::Float32, 3>&   world_bbox,
        const mi::math::Vector<mi::Float32, 3>& viewdir,
        const mi::math::Vector<mi::Float32, 3>& up,
        nv::index::IOrthographic_camera*        ocam);

    /// set camera parameter to the current predefined view index's
    /// camera parameter. Subroutine of
    /// set_predefined_view_parameter_to_camera.
    ///
    /// \param[in]  world_bbox bounding box in the world space
    /// \param[in]  cam_param  camera parameters
    /// \param[out] cam        the camera object to be updated
    void set_predefined_view_parameter_to_camera_sub(
        const mi::math::Bbox<mi::Float32, 3>& world_bbox,
        const Camera_parameter&               cam_param,
        nv::index::ICamera*                   cam);


    /// Pan the camera (screen match)
    ///
    /// \param[in] scene scene object
    /// \param[in] last_nx last mouse x position in normalized coordinate
    /// \param[in] last_ny last mouse y position in normalized coordinate
    /// \param[in] cur_nx  current mouse x position in normalized coordinate
    /// \param[in] cur_ny  current mouse y position in normalized coordinate
    /// \param[in,out] cam camera to be updated
    void pan_screen_match_camera(
        const nv::index::IScene* scene,
        mi::Float64              last_nx,
        mi::Float64              last_ny,
        mi::Float64              cur_nx,
        mi::Float64              cur_ny,
        nv::index::ICamera*      cam);

private:
    /// main window resolution
    mi::math::Vector<mi::Sint32, 2>  m_main_window_resolution;
    /// examiner mode
    Manipulation_mode m_manip_mode;
    /// Examiner's rotation center xyz position
    mi::math::Vector<mi::Float32, 3> m_examiner_rotation_center;
    /// last mouse position
    mi::math::Vector<mi::Sint32, 2>  m_last_point_2d;
    /// last point 3d on a sphere for trackball movement
    mi::math::Vector<mi::Float32, 3> m_last_point_3d;
    /// scene bounding box
    mi::math::Bbox< mi::Float32, 3 > m_scene_bbox;

    /// predefined camera parameters
    std::vector<Camera_parameter> m_predef_cam_param;
    /// current predefined view index
    mi::Size m_current_predefined_view_index;

    /// pan reference point heuristics mode: pick
    bool m_pan_reference_point_use_pick_point;
    /// pan reference point
    mi::math::Vector<mi::Float32, 3> m_pan_reference_point_object_space;
    /// pan reference point validity
    bool m_pan_reference_point_valid;

    /// pan reference point heuristics mode: bisect
    bool m_pan_reference_point_roi_bisect_and_solid;

    /// reference point zoom factor
    mi::Float32 m_reference_point_zoom_factor;
};

//----------------------------------------------------------------------

#endif /// NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_EXAMINER_MANIPULATOR_H
