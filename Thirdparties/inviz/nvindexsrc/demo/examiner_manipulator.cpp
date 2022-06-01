/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "examiner_manipulator.h"

#include <nv/index/iscene.h>

#include <iostream>

#include "multiple_camera.h"    // Multiple stereo camera test
#include "nvindex_appdata.h"
#include "scene_utility.h"
#include "utilities.h"

#include "common/forwarding_logger.h"
#include "common/common_utility.h"
#include "common/string_dict.h"

//======================================================================
//----------------------------------------------------------------------
/// adjustment of the reference point ROI bisect and solid
/// (clip_min) by a heuristics
///
/// \param[in] ref_center  reference point center
/// \param[in] bbox_radius the radius of the bounding box
/// \param[in] obj_2_world object to world space matrix
/// \param[in] ocam        camera object
/// \return updated reference center position
static mi::math::Vector<mi::Float64, 3> adjust_reference_point_by_bisect_for_pan_move(
    const mi::math::Vector<mi::Float64, 3>&    ref_center_obj_d,
    const mi::Float64                          bbox_radius_obj_d,
    const mi::math::Matrix<mi::Float64, 4, 4>& obj_2_world,
    nv::index::ICamera*                        cam)
{
    mi::math::Matrix<mi::Float64, 4, 4> world_2_obj = obj_2_world;
    if (!world_2_obj.invert())
    {
        ERROR_LOG << "cannot get world to object space matrix.";
        return ref_center_obj_d;
    }

    const mi::math::Vector<mi::Float32, 3> eye_point_ws_f = cam->get_eye_point();
    const mi::math::Vector<mi::Float64, 3> eye_point_ws_d(eye_point_ws_f);
    const mi::math::Vector<mi::Float64, 3> eye_point_obj_d =
        mi::math::transform_point(world_2_obj, eye_point_ws_d);
    const mi::math::Vector<mi::Float32, 3> viewdir_ws_f   = cam->get_view_direction();
    const mi::math::Vector<mi::Float64, 3> viewdir_ws_d(viewdir_ws_f);
    mi::math::Vector<mi::Float64, 3> viewdir_obj_d =
        mi::math::transform_vector(world_2_obj, viewdir_ws_d);
    viewdir_obj_d.normalize();

    mi::Float64 clip_min = 0.0001; // heuristic epsilon for the orthographic camera
    mi::base::Handle<nv::index::IPerspective_camera> perspective_camera(
        cam->get_interface<nv::index::IPerspective_camera>());
    if (perspective_camera.is_valid_interface())
    {
        clip_min = perspective_camera->get_clip_min();
    }

    // edge point
    const mi::math::Vector<mi::Float64, 3> edge_point_obj_d =
        ref_center_obj_d + bbox_radius_obj_d * viewdir_obj_d;

    // Is edge point visible?
    if (mi::math::dot(edge_point_obj_d - eye_point_obj_d, viewdir_obj_d) <= 0.0)
    {
        // non visible. put the reference point with minimal distance
        // in front of the camera.
        const mi::math::Vector<mi::Float64, 3> new_ref_center_obj_d =
            eye_point_obj_d + clip_min * viewdir_obj_d;
        return new_ref_center_obj_d;
    }

    // visible and enough far
    if (length(ref_center_obj_d - eye_point_obj_d) >=  bbox_radius_obj_d)
    {
        return ref_center_obj_d;
    }

    // try bisect position
    const mi::math::Vector<mi::Float64, 3> bisect_point_obj_d =
        0.5 * (eye_point_obj_d + edge_point_obj_d);
    if (length(eye_point_obj_d - bisect_point_obj_d) < clip_min)
    {
        // too close
        const mi::math::Vector<mi::Float64, 3> new_ref_center_obj_d =
            eye_point_obj_d + clip_min * viewdir_obj_d;
        return new_ref_center_obj_d;
    }

    // use bisect point
    return bisect_point_obj_d;
}

//----------------------------------------------------------------------
//======================================================================
/// Rotates the camera around a point.
///
/// \param[in] angle_rad Rotation angle in radians
/// \param[in] rot_axis  Vector defining the rotation axis
/// \param[in,out] cam camera object
void rotate_camera(
    const mi::math::Vector<mi::Float32, 3>& rotate_center,
    mi::Float32                             angle_rad,
    const mi::math::Vector<mi::Float32, 3>& rot_axis_in,
    nv::index::ICamera*                     cam)
{
    assert(cam != 0);

    const mi::math::Vector<mi::Float32, 3> up_dir        = cam->get_up_direction();
    const mi::math::Vector<mi::Float32, 3> eye_point     = cam->get_eye_point();
    const mi::math::Vector<mi::Float32, 3> view_dir      = cam->get_view_direction();

    mi::math::Vector<mi::Float32, 3> ev[3];
    Camera_tool::compute_camera_basis(view_dir, up_dir, ev[0], ev[1], ev[2]);

    mi::math::Vector<mi::Float32, 3> rot_axis(rot_axis_in);
    rot_axis.normalize();       // for Matrix::set_rotation

    // OpenGL
    rot_axis = rot_axis.x * ev[0] + rot_axis.y * ev[1] - rot_axis.z * ev[2];

    // rotate camera
    mi::math::Matrix<mi::Float32, 4, 4> rotmat(1.0f); // = I
    rotmat.set_rotation(rot_axis, (mi::Float64)angle_rad);

    // Move camera position to origin (offset)
    const mi::math::Vector<mi::Float32, 3> epos = eye_point - rotate_center;

    // Rotate the camera position and back to the offset
    mi::math::Vector<mi::Float32, 3> new_eye_pos    = transform_point(rotmat, epos) + rotate_center;
    mi::math::Vector<mi::Float32, 3> new_view_dir   = transform_vector(rotmat, view_dir);
    mi::math::Vector<mi::Float32, 3> new_up_vec     = transform_vector(rotmat, up_dir);

    // Update the camera parameters
    cam->set_eye_point(new_eye_pos);
    cam->set_view_direction(new_view_dir);
    cam->set_up_direction(new_up_vec);
}

//----------------------------------------------------------------------
//======================================================================
//----------------------------------------------------------------------
Examiner_manipulator::Examiner_manipulator()
    :
    m_main_window_resolution(1024, 1024),
    m_manip_mode(EM_TrackballRotation),
    m_examiner_rotation_center(0.0f, 0.0f, 0.0f),
    m_last_point_2d(-1, -1),
    m_last_point_3d(0.0, 0.0, 0.0),
    m_scene_bbox(-1.0f, -1.0f, -1.0f, 1.0f, 1.0f, 1.0f),
    m_current_predefined_view_index(0),
    m_pan_reference_point_use_pick_point(false),
    m_pan_reference_point_object_space(0.0f, 0.0f, 0.0f),
    m_pan_reference_point_valid(false),
    m_pan_reference_point_roi_bisect_and_solid(false),
    m_reference_point_zoom_factor(1.5)
{
    // empty
}

//----------------------------------------------------------------------
void Examiner_manipulator::initialize_by_dict(const nv::index_common::String_dict& dict)
{
    const std::string prefix = "app::examiner::";

    // rotate mode and center
    {
        const std::string center_type = dict.get(prefix + "initial_rotation_center::type", "roi_center");

        if (center_type == "world_coordinates")
        {
            const mi::math::Vector<mi::Float32, 3> initial_rotation_center =
                nv::index_common::get_vec_float32_3(
                    dict.get(prefix + "initial_rotation_center::world_coordinates", "0 0 0"));
            set_examiner_rotation_center(initial_rotation_center);
            INFO_LOG << "Examiner: initial_rotation_center: " << initial_rotation_center;
        }
        else if (center_type == "roi_center")
        {
            // This initialization will be delayed after scene
            // setup. (This needs the global transformation to
            // calculate the center of the roi in the world space.)
        }
        else
        {
            WARN_LOG << "No initial rotation center define method initial_rotation_center::type = "
                     << center_type << ", set it to origin.";
            const mi::math::Vector<mi::Float32, 3> default_initial_rotation_center(0.0f, 0.0f, 0.0f);
            set_examiner_rotation_center(default_initial_rotation_center);
        }
    }        

    // predefined view
    if (nv::index_common::get_bool(dict.get(prefix + "predefined_view::enable", "no")))
    {
        initialize_predefined_view_by_dict(dict);
    }
    else
    {
        // default top view setting
        nv::index_common::String_dict default_predef;
        default_predef.insert("app::examiner::predefined_view::view_count",  "1");
        default_predef.insert("app::examiner::predefined_view::0::name",     "Default Top view");
        default_predef.insert("app::examiner::predefined_view::0::view_all", "yes");
        default_predef.insert("app::examiner::predefined_view::0::from",     "0 0 0");
        default_predef.insert("app::examiner::predefined_view::0::dir",      "0 0 -1");
        default_predef.insert("app::examiner::predefined_view::0::up",       "0 1 0");
        default_predef.insert("app::examiner::predefined_view::0::aspect",   "1");
        default_predef.insert("app::examiner::predefined_view::0::aperture", "0.033");
        default_predef.insert("app::examiner::predefined_view::0::focal",    "0.03");
        default_predef.insert("app::examiner::predefined_view::0::clip_min", "0.001");
        default_predef.insert("app::examiner::predefined_view::0::clip_max", "10000");
        default_predef.insert("app::examiner::predefined_view::0::orthographic", "no");
        default_predef.insert("app::examiner::predefined_view::0::view_all_bbox_type", "scene_roi");

        initialize_predefined_view_by_dict(default_predef);
    }

    // experimental pan heuristic
    {
        const bool is_heuristics_pick = nv::index_common::get_bool(dict.get(prefix + "pan_reference_point_use_pick_point", "no"));
        set_pan_reference_point_use_pick_point(is_heuristics_pick);

        const bool is_heuristics_bisect = nv::index_common::get_bool(dict.get(prefix + "pan_reference_point_roi_bisect_and_solid", "no"));
        set_pan_reference_point_roi_bisect_and_solid(is_heuristics_bisect);
    }

    // experimental key zoom
    {
        const mi::Float32 zoom_factor =
            nv::index_common::get_float32(dict.get(prefix + "reference_point_zoom_factor", "2.0"));
        set_reference_point_zoom_factor(zoom_factor);
    }
}

//----------------------------------------------------------------------
void Examiner_manipulator::set_manipulation_mode(
    Examiner_manipulator::Manipulation_mode manip_mode)
{
    m_manip_mode = manip_mode;
}

//----------------------------------------------------------------------
Examiner_manipulator::Manipulation_mode Examiner_manipulator::get_manipulation_mode() const
{
    return m_manip_mode;
}

//----------------------------------------------------------------------
void Examiner_manipulator::set_main_window_resolution(const mi::math::Vector<mi::Sint32, 2> & win_resolution)
{
    assert(win_resolution.x > 0);
    assert(win_resolution.y > 0);
    m_main_window_resolution = win_resolution;
}

//----------------------------------------------------------------------
mi::math::Vector<mi::Sint32, 2> Examiner_manipulator::get_main_window_resolution() const
{
    return m_main_window_resolution;
}

//----------------------------------------------------------------------
mi::math::Vector<mi::Uint32, 2> Examiner_manipulator::get_main_window_resolution_uint32() const
{
    return mi::math::Vector<mi::Uint32, 2>(static_cast<mi::Uint32>(m_main_window_resolution.x),
                                           static_cast<mi::Uint32>(m_main_window_resolution.y));
}

//----------------------------------------------------------------------
void Examiner_manipulator::set_scene_bbox(const mi::math::Bbox<mi::Float32, 3>& scene_bbox)
{
    m_scene_bbox = scene_bbox;
}

//----------------------------------------------------------------------
mi::math::Bbox<mi::Float32, 3> Examiner_manipulator::get_scene_bbox() const
{
    return m_scene_bbox;
}

//----------------------------------------------------------------------
void Examiner_manipulator::set_examiner_rotation_center(
    const mi::math::Vector<mi::Float32, 3>& rot_center)
{
    m_examiner_rotation_center = rot_center;
    // INFO_LOG << "Set rotation center: " << m_examiner_rotation_center;
}

//----------------------------------------------------------------------
mi::math::Vector<mi::Float32, 3> Examiner_manipulator::get_examiner_rotation_center() const
{
    return m_examiner_rotation_center;
}

//----------------------------------------------------------------------
void Examiner_manipulator::set_pan_reference_point_object_space(
    const mi::math::Vector<mi::Float32, 3>& pan_ref)
{
    m_pan_reference_point_object_space = pan_ref;
}

//----------------------------------------------------------------------
mi::math::Vector<mi::Float32, 3> Examiner_manipulator::get_pan_reference_point_object_space() const
{
    return m_pan_reference_point_object_space;
}

//----------------------------------------------------------------------
void Examiner_manipulator::set_pan_reference_point_valid(bool is_valid)
{
    m_pan_reference_point_valid = is_valid;
}

//----------------------------------------------------------------------
bool Examiner_manipulator::is_pan_reference_point_valid() const
{
    return m_pan_reference_point_valid;
}

//----------------------------------------------------------------------
void Examiner_manipulator::set_reference_point_zoom_factor(mi::Float32 zoom_factor)
{
    if (zoom_factor <= 1.0)
    {
        ERROR_LOG << "Invalid zoom factor: " << zoom_factor << ", should be > 1.0.";
        return;
    }
    m_reference_point_zoom_factor = zoom_factor;
}

//----------------------------------------------------------------------
mi::Float32 Examiner_manipulator::get_reference_point_zoom_factor() const
{
    return m_reference_point_zoom_factor;
}

//----------------------------------------------------------------------
void Examiner_manipulator::mouse_press_event(
    mi::Sint32                        button,
    Mouse_action                      mouse_action,
    mi::Sint32                        mouse_x,
    mi::Sint32                        mouse_y,
    const mi::neuraylib::Tag&         camera_tag,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    if (button == Left_button)
    {
        if (mouse_action == Button_down)
        {
            // init/invalidate the mouse info
            this->clear_last_mouse_info();

            if (this->get_manipulation_mode() == EM_TrackballRotation)
            {
                this->set_last_mouse2d(mouse_x, mouse_y);
                // FIXME: m_win_width and m_win_height must be update when
                // window size is changed.
                mi::math::Vector<mi::Float32, 3> point_3d(0.0, 0.0, 0.0);
                const bool r = map_to_sphere(m_last_point_2d, &point_3d);
                if (!r)
                {
                    WARN_LOG << "Fail to map mouse press position on the sphere. Ignored.";
                    return;
                }
                this->set_last_mouse3d(point_3d);
            }
            else if (this->get_manipulation_mode() == EM_Pan)
            {
                this->set_last_mouse2d(mouse_x, mouse_y); // record the start position
            }
            else if (this->get_manipulation_mode() == EM_Zoom)
            {
                this->set_last_mouse2d(mouse_x, mouse_y); // record the start position
            }
            else
            {
                assert(false);  // never happen
            }
        }
    }
    else if (button == Middle_button)
    {
        if (mouse_action == Button_down)
        {
            // init/invalidate the mouse info
            this->clear_last_mouse_info();

            this->set_last_mouse2d(mouse_x, mouse_y); // record the start position
        }
    }
}

//----------------------------------------------------------------------
void Examiner_manipulator::mouse_motion_event(
    mi::Sint32                        button,
    mi::Sint32                        mouse_x,
    mi::Sint32                        mouse_y,
    const mi::neuraylib::Tag&         camera_tag,
    const mi::neuraylib::Tag&         scene_tag,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    mi::math::Vector<mi::Sint32, 2> new_point_2d(mouse_x, mouse_y);
    bool const is_inside =
        (new_point_2d.x >= 0) && (new_point_2d.x < m_main_window_resolution.x) &&
        (new_point_2d.y >= 0) && (new_point_2d.y < m_main_window_resolution.y);
    if (!is_inside)
    {
        return;             // mouse pointer is out of range, do nothing.
    }

    if (button == Left_button)
    {
        if (this->get_manipulation_mode() == EM_TrackballRotation)
        {
            mi::math::Vector<mi::Float32, 3> new_point_3d(0.0, 0.0, 0.0);
            const bool r = map_to_sphere(new_point_2d, &new_point_3d);
            if (!r)
            {
                WARN_LOG << "Fail to map mouse motion position on the sphere. Ignored.";
                return;
            }

            this->move_camera_trackballrotate_mode(new_point_3d, camera_tag, dice_transaction);

            this->set_last_mouse2d(mouse_x, mouse_y);
            this->set_last_mouse3d(new_point_3d);
        }
        else if (this->get_manipulation_mode() == EM_Pan)
        {
            this->move_camera_pan_mode(new_point_2d, camera_tag, scene_tag, dice_transaction);
            this->set_last_mouse2d(mouse_x, mouse_y);
        }
        else if (this->get_manipulation_mode() == EM_Zoom)
        {
            this->move_camera_zoom_mode(new_point_2d, camera_tag, scene_tag, dice_transaction);
            this->set_last_mouse2d(mouse_x, mouse_y);
        }
        else
        {
            assert(false);  // never happen
        }
    }
    else if(button == Middle_button)
    {
        // Pan with the middle button
        this->move_camera_pan_mode(new_point_2d, camera_tag, scene_tag, dice_transaction);
        this->set_last_mouse2d(mouse_x, mouse_y);
    }
}

//----------------------------------------------------------------------
void Examiner_manipulator::mouse_wheel_event(
    mi::Sint32                        wheel_delta,
    const mi::neuraylib::Tag&         camera_tag,
    const mi::neuraylib::Tag&         scene_tag,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    mi::base::Handle<const nv::index::IScene> scene(
        dice_transaction->access<nv::index::IScene>(scene_tag));
    assert(scene.is_valid_interface());

    // zoom factor is relative to the bbox
    const mi::math::Bbox< mi::Float32, 3 > transformed_bbox =
        get_transformed_bbox(m_scene_bbox, scene->get_transform_matrix());
    const mi::Float64 scr_radius  = transformed_bbox.diagonal_length() * 0.5;
    const mi::Float64 z_dist      = wheel_delta * scr_radius;
    const mi::Float64 zoom_factor = 0.4; // Heuristics
    const mi::Float64 delta_z     = zoom_factor * z_dist;
    if (fabs(delta_z) < 1.0e-8)
    {
        DEBUG_LOG << "|zoom distance| < 1.0e-8.";
    }

    { // this scope is important to have the edit-handle run out-of-scope before the next edit
        mi::base::Handle< nv::index::IOrthographic_camera > ortho_cam(
            dice_transaction->edit< nv::index::IOrthographic_camera >(camera_tag));
        if (ortho_cam)
        {
            mi::math::Vector<mi::Float32, 3> eye_poi   (ortho_cam->get_eye_point());
            mi::math::Vector<mi::Float32, 3> view_dir  (ortho_cam->get_view_direction());

            mi::Float64 ap = std::max(0.00001, ortho_cam->get_aperture() - delta_z);
            ortho_cam->set_aperture(ap);
        }
    }

    {
        mi::base::Handle< nv::index::IPerspective_camera > perspective_cam(
            dice_transaction->edit< nv::index::IPerspective_camera >(camera_tag));
        if (perspective_cam)
        {
            mi::math::Vector<mi::Float32, 3> eye_poi (perspective_cam->get_eye_point());
            mi::math::Vector<mi::Float32, 3> view_dir(perspective_cam->get_view_direction());
            eye_poi = eye_poi + delta_z * view_dir;
            perspective_cam->set_eye_point(eye_poi);
        }
    }
}

//----------------------------------------------------------------------
void Examiner_manipulator::set_last_mouse2d(mi::Sint32 mouse_x, mi::Sint32 mouse_y)
{
    m_last_point_2d = mi::math::Vector<mi::Sint32, 2>(mouse_x, mouse_y);
}

//----------------------------------------------------------------------
void Examiner_manipulator::set_last_mouse3d(const mi::math::Vector<mi::Float32, 3> & point_3d)
{
    m_last_point_3d = point_3d;
}

//----------------------------------------------------------------------
void Examiner_manipulator::clear_predefined_camera()
{
    m_predef_cam_param.clear();
}

//----------------------------------------------------------------------
void Examiner_manipulator::append_predefined_camera(Camera_parameter& cam_param)
{
    m_predef_cam_param.push_back(cam_param);
}

//----------------------------------------------------------------------
mi::Size Examiner_manipulator::nb_predefined_camera_parameter() const
{
    return m_predef_cam_param.size();
}

//----------------------------------------------------------------------
Camera_parameter Examiner_manipulator::get_predefined_camera_parameter(mi::Size idx)
{
    const std::string mn = "Examiner_manipulator::get_camera_parameter: ";
    if (idx >= m_predef_cam_param.size())
    {
        ERROR_LOG << mn << "index is out of range";
        return Camera_parameter();
    }

    return m_predef_cam_param[idx];
}

//----------------------------------------------------------------------
std::vector<std::string> Examiner_manipulator::get_predefined_view_name_vec() const
{
    std::vector<std::string> predef_name_vec;
    for (std::vector<Camera_parameter>::const_iterator ci = m_predef_cam_param.begin();
         ci != m_predef_cam_param.end();
         ++ci)
    {
        std::string view_name = "<no_name>";
        const bool is_valid_key = ci->get_string_value("name", view_name);
        if (!is_valid_key)
        {
            ERROR_LOG << "Cannot find the predefined view";
        }

        // GUI limitation
        if (view_name.find(',') != std::string::npos)
        {
            const std::string orig_name = view_name;
            const mi::Size nb_view_name = view_name.length();
            for (mi::Size i = 0; i < nb_view_name; ++i)
            {
                if (view_name[i] == ',')
                {
                    view_name[i] = ' ';
                }
            }
            ERROR_LOG << "Cannot use character ',' in the name. Replaced ',' with ' ': "
                      << orig_name << " -> " << view_name;
        }
        predef_name_vec.push_back(view_name);
    }

    return predef_name_vec;
}

//----------------------------------------------------------------------
void Examiner_manipulator::set_predefined_view_index(mi::Size idx)
{
    const std::string mn = "Examiner_manipulator::set_predefined_view_index: ";
    if (idx >= m_predef_cam_param.size())
    {
        ERROR_LOG << mn << "index out of range: " << idx << ", ignored.";
        return;
    }

    m_current_predefined_view_index = idx;
}

//----------------------------------------------------------------------
mi::Size Examiner_manipulator::get_predefined_view_index() const
{
    return m_current_predefined_view_index;
}

//----------------------------------------------------------------------
static mi::math::Bbox<mi::Float32, 3> get_bbox_by_type(
    const std::string&                bbox_type,
    const nv::index::ISession*        session,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    assert(dice_transaction != 0);

    mi::math::Bbox<mi::Float32, 3> bbox;
    if (bbox_type == "volume")
    {
        // volume bbox
        const std::vector<mi::neuraylib::Tag> volume_tag_vec = Nvindex_AppData::instance()->get_volume_tag_vec();

        if (!(volume_tag_vec.empty()))
        {
            // Note: currently we use the first volume only
            mi::base::Handle<const nv::index::IRegular_volume> volume(
                dice_transaction->access<const nv::index::IRegular_volume>(volume_tag_vec.at(0)));
            assert(volume.is_valid_interface());
            bbox = volume->get_XYZ_clipped_bounding_box();
        }
        else
        {
            WARN_LOG << "Cannot retrieve the first volume bounding box.";
        }
    }
    else if (bbox_type == "scene_roi")
    {
        assert(session != 0);
        // xyz region of interest bbox
        bbox = get_XYZ_global_region_of_interest_bbox(session, dice_transaction);
    }
    else
    {
        ERROR_LOG << "Unknown bbox_type: " << bbox_type << ", use scene_roi.";

        // fallback
        assert(session != 0);
        bbox = get_XYZ_global_region_of_interest_bbox(session, dice_transaction);
    }

    return bbox;
}


//----------------------------------------------------------------------
void Examiner_manipulator::set_predefined_view_parameter_to_camera(
    const nv::index::ISession*        session,
    const nv::index::IScene*          scene,
    Multiple_stereo_camera*           ms_camera,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    assert(session   != 0);
    assert(scene     != 0);
    assert(ms_camera != 0);
    assert(dice_transaction != 0);

    const std::string mn = "Examiner_manipulator::set_predefined_view_parameter_to_camera: ";
    const mi::Size predef_idx = get_predefined_view_index();
    if (predef_idx >= nb_predefined_camera_parameter())
    {
        ERROR_LOG << mn  << "predefined view index is out of range: " << predef_idx
                  << "/" << nb_predefined_camera_parameter();
        return;
    }
    Camera_parameter cam_param = get_predefined_camera_parameter(predef_idx);

    std::string bbox_type = "scene";
    {
        const bool is_valid_key = cam_param.get_string_value("view_all_bbox_type", bbox_type);
        assert(is_valid_key);
        nv::index_common::no_unused_variable_warning_please(is_valid_key);
    }

    bool is_orthographic_camera = false;
    {
        const bool is_valid_key = cam_param.get_bool_value("orthographic", is_orthographic_camera);
        assert(is_valid_key);
        nv::index_common::no_unused_variable_warning_please(is_valid_key);
    }
    ms_camera->set_use_ortho_camera(is_orthographic_camera);

    const mi::math::Bbox<mi::Float32, 3> current_bbox =
        get_bbox_by_type(bbox_type, session, dice_transaction);
    const mi::math::Bbox< mi::Float32, 3 > translated_bbox =
        get_transformed_bbox(current_bbox, scene->get_transform_matrix());

    // Set the parameters to both camera to have an illusion of the
    // same camera has a different mode.
    {
        mi::base::Handle<nv::index::ICamera> ortho_cam(
            dice_transaction->edit<nv::index::ICamera>(
                ms_camera->get_ortho_camera_tag()));

        set_predefined_view_parameter_to_camera_sub(translated_bbox,
                                                    cam_param,
                                                    ortho_cam.get());
    }
    {
        mi::base::Handle<nv::index::ICamera> main_cam(
            dice_transaction->edit<nv::index::ICamera>(
                ms_camera->get_main_base_camera_tag()));

        set_predefined_view_parameter_to_camera_sub(translated_bbox,
                                                    cam_param,
                                                    main_cam.get());
    }

}

//----------------------------------------------------------------------
void Examiner_manipulator::set_pan_reference_point_use_pick_point(bool is_enabled)
{
    m_pan_reference_point_use_pick_point = is_enabled;
}

//----------------------------------------------------------------------
bool Examiner_manipulator::is_pan_reference_point_use_pick_point() const
{
    return m_pan_reference_point_use_pick_point;
}

//----------------------------------------------------------------------
void Examiner_manipulator::set_pan_reference_point_roi_bisect_and_solid(bool is_enabled)
{
    m_pan_reference_point_roi_bisect_and_solid = is_enabled;
}

//----------------------------------------------------------------------
bool Examiner_manipulator::is_pan_reference_point_roi_bisect_and_solid() const
{
    return m_pan_reference_point_roi_bisect_and_solid;
}

//----------------------------------------------------------------------
void Examiner_manipulator::initialize_predefined_view_by_dict(
    const nv::index_common::String_dict& dict)
{
    const std::string prefix = "app::examiner::";

    const mi::Size predefined_view_count =
        nv::index_common::get_sint32(dict.get(prefix + "predefined_view::view_count", "0"));
    if (predefined_view_count <= 0)
    {
        ERROR_LOG << "illegal predefined view_count: " << predefined_view_count
                  << "No predefined view is set.";

        return;
    }

    for (mi::Size i = 0; i < predefined_view_count; ++i)
    {
        std::stringstream sstr;
        sstr << prefix << "predefined_view::" << i << "::";
        Camera_parameter cam_param;
        const mi::Sint32 param_count =
            cam_param.set_parameter_by_string_dict(dict, sstr.str());
        if (param_count <= 0)
        {
            ERROR_LOG << "no predefined camera parameters set for predefined camera: " << i;
        }

        {
            // debug out
            nv::index_common::String_dict cam_out;
            cam_param.get_parameter_to_string_dict(sstr.str(), cam_out);
            std::stringstream osstr;
            cam_out.write(osstr, "");
            DEBUG_LOG << "predefined_view[" << i << "]\n" << osstr.str();
        }
        append_predefined_camera(cam_param);
    }

    // current index can be set only after the parameters has set
    mi::Sint32 predefined_view_index_s32 =
        nv::index_common::get_sint32(dict.get(prefix + "predefined_view::view_index", "0"));
    if (predefined_view_index_s32 < 0)
    {
        ERROR_LOG << "illegal predefined_view::view_index: " << predefined_view_index_s32
                  << ", set index 0.";
        predefined_view_index_s32 = 0;
    }
    set_predefined_view_index(static_cast<mi::Size>(predefined_view_index_s32));
}

//----------------------------------------------------------------------
bool Examiner_manipulator::map_to_sphere(
    const mi::math::Vector<mi::Sint32, 2>& mouse_pos,
    mi::math::Vector<mi::Float32, 3>*      p_3d_pos)
{
    mi::Float64 norm_x = 0.0;
    mi::Float64 norm_y = 0.0;
    if (!this->get_normalized_screen_coords(mouse_pos.x, mouse_pos.y, norm_x, norm_y))
    {
        (*p_3d_pos)[0] = 0.0;
        (*p_3d_pos)[1] = 0.0;
        (*p_3d_pos)[2] = 0.0;
        return false;
    }

    // normalization into the rotational angle from GMU.
    mi::Float64 sx = sin(0.5 * M_PI * norm_x);
    mi::Float64 sy = sin(0.5 * M_PI * norm_y);
    mi::Float64 sx2sy2 = (sx * sx) + (sy * sy);

    mi::Float64 sz = 0.0;
    if (sx2sy2 < 1.0)
    {
        sz = sqrt(1.0 - sx2sy2);
    }

    (*p_3d_pos)[0] = static_cast<mi::Float32>(sx);
    (*p_3d_pos)[1] = static_cast<mi::Float32>(sy);
    (*p_3d_pos)[2] = static_cast<mi::Float32>(sz);

    //         fprintf(stderr, "map_to_sphere: mouse %g %g (%g %g) [%g %g] -> 3D %f %f %f\n",
    //                 wx, wy, norm_x, norm_y, win_width, win_height,
    //                 (*p_3d_pos)[0], (*p_3d_pos)[1], (*p_3d_pos)[2]);

    return true;
}

//----------------------------------------------------------------------
void Examiner_manipulator::move_camera_trackballrotate_mode(
    const mi::math::Vector<mi::Float32, 3>& new_point_3d,
    const mi::neuraylib::Tag&               camera_tag,
    mi::neuraylib::IDice_transaction*       dice_transaction)
{
    if (m_last_point_3d.z == 0.0)
    {
        return;             // the last state is illegal, ignore.
    }

    if (new_point_3d.z == 0.0)
    {
        return;             // the new state is illegal, ignore.
    }

    mi::Float64 rotate_angle_d = 0.0f;
    mi::math::Vector<mi::Float32, 3> rot_axis  = cross(new_point_3d, m_last_point_3d);
    const mi::Float64                cos_angle = dot(m_last_point_3d, new_point_3d);
    if (fabs(cos_angle) < 1.0)
    {
        rotate_angle_d = 2.0 * acos(cos_angle); // radian, spinor rotation
    }
    const mi::Float32 rotate_angle_f = static_cast<mi::Float32>(rotate_angle_d);

    //  DEBUG_LOG << "last: "    << m_last_point_3d << ", new: " << new_point_3d
    //            << ", axis: "  << rot_axis << ", cos: " << cos_angle;

    // update camera
    {
        mi::base::Handle< nv::index::ICamera > cam(
            dice_transaction->edit< nv::index::ICamera >(camera_tag));
        assert(cam.is_valid_interface());

        const mi::math::Vector<mi::Float32, 3> rotation_center =
            get_examiner_rotation_center();
        rotate_camera(rotation_center, rotate_angle_f, rot_axis, cam.get());

        // nv::index_common::print_camera_param(cam);
    }
}

//----------------------------------------------------------------------
void Examiner_manipulator::move_camera_pan_mode(
    const mi::math::Vector<mi::Sint32, 2>& new_point_2d,
    const mi::neuraylib::Tag&              camera_tag,
    const mi::neuraylib::Tag&              scene_tag,
    mi::neuraylib::IDice_transaction*      dice_transaction)
{
    mi::Float64 cur_nx  = 0.0, cur_ny  = 0.0, last_nx = 0.0, last_ny = 0.0;
    if (!get_normalized_screen_coords(new_point_2d.x, new_point_2d.y, cur_nx, cur_ny))
    {
        return;       // invalid mouse position or screen, ignore.
    }

    if (!this->get_normalized_screen_coords(m_last_point_2d.x, m_last_point_2d.y, last_nx, last_ny))
    {
        return;       // invalid mouse position or screen, ignore.
    }

    // update the camera
    {
        mi::base::Handle< nv::index::ICamera > cam(
            dice_transaction->edit< nv::index::ICamera >(camera_tag));
        assert(cam.is_valid_interface());

        mi::base::Handle<const nv::index::IScene> scene(
            dice_transaction->access<nv::index::IScene>(scene_tag));
        assert(scene.is_valid_interface());

        // screen and mouse pointer movement will match in the pan mode.
        this->pan_screen_match_camera(scene.get(), last_nx, last_ny, cur_nx, cur_ny, cam.get());
    }
}

//----------------------------------------------------------------------
void Examiner_manipulator::move_camera_zoom_mode(
    const mi::math::Vector<mi::Sint32, 2>& new_point_2d,
    const mi::neuraylib::Tag&              camera_tag,
    const mi::neuraylib::Tag&              scene_tag,
    mi::neuraylib::IDice_transaction*      dice_transaction)
{
    mi::Float64 cur_nx = 0.0;
    mi::Float64 cur_ny = 0.0;
    if (!get_normalized_screen_coords(new_point_2d.x, new_point_2d.y, cur_nx, cur_ny))
    {
        // invalid mouse position or screen, ignore.
        return;
    }

    mi::Float64 last_nx = 0.0;
    mi::Float64 last_ny = 0.0;
    if (!this->get_normalized_screen_coords(m_last_point_2d.x, m_last_point_2d.y, last_nx, last_ny))
    {
        // invalid mouse position or screen, ignore.
        return;
    }

    // update camera
    {
        mi::base::Handle<const nv::index::IScene> scene(
            dice_transaction->access<nv::index::IScene>(scene_tag));
        assert(scene.is_valid_interface());

        const mi::math::Bbox< mi::Float32, 3 > transformed_bbox =
            get_transformed_bbox(m_scene_bbox, scene->get_transform_matrix());
        const mi::Float64 scr_radius = transformed_bbox.diagonal_length() * 0.5;
        const mi::Float64 z_dist      = 2.0 * (last_ny - cur_ny) * scr_radius;
        const mi::Float64 zoom_factor = 2.0;
        const mi::Float64 delta_z     = zoom_factor * z_dist;
        if (fabs(delta_z) < 1.0e-8)
        {
            DEBUG_LOG << "|zoom distance| < 1.0e-8.";
        }

        { // this scope is important to have the edit-handle run out-of-scope before the next edit
            mi::base::Handle< nv::index::IOrthographic_camera > ortho_cam(
                dice_transaction->edit< nv::index::IOrthographic_camera >(camera_tag));
            if (ortho_cam)
            {
                mi::math::Vector<mi::Float32, 3> eye_poi   (ortho_cam->get_eye_point());
                mi::math::Vector<mi::Float32, 3> view_dir  (ortho_cam->get_view_direction());

                mi::Float64 ap = std::max(0.00001, ortho_cam->get_aperture() - delta_z);
                ortho_cam->set_aperture(ap);
            }
        }

        {
            mi::base::Handle< nv::index::IPerspective_camera > perspective_cam(
                dice_transaction->edit< nv::index::IPerspective_camera >(camera_tag));
            if (perspective_cam)
            {
                mi::math::Vector<mi::Float32, 3> eye_poi (perspective_cam->get_eye_point());
                mi::math::Vector<mi::Float32, 3> view_dir(perspective_cam->get_view_direction());
                eye_poi = eye_poi + delta_z * view_dir;
                perspective_cam->set_eye_point(eye_poi);
            }
        }
    }
}

//----------------------------------------------------------------------
void Examiner_manipulator::camera_zoom(
    mi::Float32                       zoom_factor,
    const mi::neuraylib::Tag&         camera_tag,
    const mi::neuraylib::Tag&         scene_tag,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    // update camera
    {
        mi::base::Handle< nv::index::ICamera > cam(
            dice_transaction->edit< nv::index::ICamera >(camera_tag));
        assert(cam.is_valid_interface());
        mi::base::Handle<const nv::index::IScene> scene(
            dice_transaction->access<nv::index::IScene>(scene_tag));
        assert(scene.is_valid_interface());

        mi::math::Vector<mi::Float32, 3>       eye_poi (cam->get_eye_point());
        const mi::math::Vector<mi::Float32, 3> view_dir(cam->get_view_direction());

        eye_poi = eye_poi + zoom_factor * view_dir;
        cam->set_eye_point(eye_poi);
    }
}

//----------------------------------------------------------------------
void Examiner_manipulator::camera_zoom_to_pan_reference_point(
    bool                     is_zoom_in,
    const nv::index::IScene* scene,
    nv::index::ICamera*      cam)
{
    assert(scene != 0);
    assert(cam != 0);

    if (!is_pan_reference_point_valid())
    {
        INFO_LOG << "No valid zoom in/out reference point. Pick a point.";
        return;
    }

    const mi::math::Matrix<mi::Float32, 4, 4> obj_2_world_f = scene->get_transform_matrix();
    const mi::math::Matrix<mi::Float64, 4, 4> obj_2_world_d(obj_2_world_f);
    const mi::math::Vector<mi::Float64, 3>    ref_point_obj_d(get_pan_reference_point_object_space());
    const mi::math::Vector<mi::Float64, 3>    ref_point_ws_d = mi::math::transform_point(obj_2_world_d, ref_point_obj_d);
    const mi::math::Vector<mi::Float64, 3>    cam_ws_d(mi::math::Vector<mi::Float32, 3>(cam->get_eye_point()));
    const mi::math::Vector<mi::Float64, 3>    vdir_d(mi::math::Vector<mi::Float32, 3>(cam->get_view_direction()));
    const mi::math::Vector<mi::Float64, 3>    rdir_d(ref_point_ws_d - cam_ws_d);
    mi::Float64 clip_min = 0.0001; // heuristic epsilon for the orthographic camera
    mi::base::Handle<nv::index::IPerspective_camera> perspective_camera(
        cam->get_interface<nv::index::IPerspective_camera>());
    if (perspective_camera.is_valid_interface())
    {
        clip_min = perspective_camera->get_clip_min();
    }

    if (length(rdir_d) < clip_min)
    {
        WARN_LOG << "too close to the reference point for computing zoom in/out.";
        return;
    }

    mi::math::Vector<mi::Float64, 3>  n_vdir_d = vdir_d;
    mi::math::Vector<mi::Float64, 3>  n_rdir_d = rdir_d;
    if (!n_vdir_d.normalize())
    {
        ERROR_LOG << "Invalid view vector.";
        return;
    }

    if (!n_rdir_d.normalize())
    {
        ERROR_LOG << "Invalid reference-point-to vector.";
        return;
    }

    const mi::Float64 cos_vr = mi::math::dot(n_vdir_d, n_rdir_d);
    if (mi::math::abs(cos_vr) < 0.3)
    {
        WARN_LOG << "camera does not looking at the reference direction.";
        return;
    }

    const mi::Float64 r_dist = mi::math::length(rdir_d);
    const mi::Float64 along_view_dist = r_dist / cos_vr;

    const mi::Float64 zoom_factor = get_reference_point_zoom_factor();
    if (zoom_factor <= 1.0)
    {
        ERROR_LOG << "Zoom factor must be > 1.0: " << zoom_factor;
        return;
    }

    const mi::Float64 zoom_in_factor  = (zoom_factor - 1.0) / zoom_factor;
    const mi::Float64 zoom_out_factor = zoom_factor - 1.0;

    if (is_zoom_in)
    {
        // zoom in
        const mi::math::Vector<mi::Float64, 3> zoom_in_eye_pos_d =
            cam_ws_d + zoom_in_factor * along_view_dist * vdir_d;
        const mi::math::Vector<mi::Float32, 3> zoom_in_eye_pos_f(zoom_in_eye_pos_d);
        cam->set_eye_point(zoom_in_eye_pos_f);
    }
    else
    {
        // zoom out
        const mi::math::Vector<mi::Float64, 3> zoom_out_eye_pos_d =
            cam_ws_d - (zoom_out_factor * along_view_dist * vdir_d);
        const mi::math::Vector<mi::Float32, 3> zoom_out_eye_pos_f(zoom_out_eye_pos_d);
        cam->set_eye_point(zoom_out_eye_pos_f);
    }
}

//----------------------------------------------------------------------
void Examiner_manipulator::clear_last_mouse_info()
{
    m_last_point_2d = mi::math::Vector<mi::Sint32, 2>(-1, -1);
    m_last_point_3d = mi::math::Vector<mi::Float32, 3>(0.0, 0.0, 0.0);
}

//----------------------------------------------------------------------
bool Examiner_manipulator::get_normalized_screen_coords(mi::Sint32 scr_x_i,
                                                        mi::Sint32 scr_y_i,
                                                        mi::Float64 & norm_x,
                                                        mi::Float64 & norm_y)
{
    // screen should be exist
    mi::math::Vector<mi::Sint32,2> window_resolution = this->get_main_window_resolution();
    assert((window_resolution.x > 0) && (window_resolution.y > 0));

    const mi::Float64 win_width  = (mi::Float64)(window_resolution.x);
    const mi::Float64 win_height = (mi::Float64)(window_resolution.y);
    const mi::Float64 wx         = (mi::Float64)(scr_x_i);
    const mi::Float64 wy         = (mi::Float64)(scr_y_i);
    norm_x = (wx - (win_width / 2.0))  / win_width;
    norm_y = ((win_height / 2.0) - wy) / win_height;

    if ((wx < 0.0) || (wx >= win_width) || (wy < 0.0) || (wy >= win_height))
    {
        WARN_LOG << "Out of range coordinates. [" << wx << ", " << wy << "] should be in ["
                 << win_width << " " << win_height << "]";
        return false;
    }

    return true;
}

//----------------------------------------------------------------------
void Examiner_manipulator::view_bbox_all(
    const mi::math::Bbox<mi::Float32, 3>&   world_bbox,
    const mi::math::Vector<mi::Float32, 3>& viewdir,
    const mi::math::Vector<mi::Float32, 3>& up,
    nv::index::ICamera*                     cam)
{
    if (world_bbox.empty())
    {
        WARN_LOG << "view_bbox_all: uninitialized bbox, ignored.";
        return;
    }

    // perspective camera view bbox all
    mi::base::Handle<nv::index::IPerspective_camera> perspective_camera(
        cam->get_interface<nv::index::IPerspective_camera>());
    if (perspective_camera.is_valid_interface())
    {
        view_bbox_all_perspective(world_bbox, viewdir, up, perspective_camera.get());
        return;
    }

    // orthographic camera view bbox all
    mi::base::Handle<nv::index::IOrthographic_camera> ortho_cam(
        cam->get_interface<nv::index::IOrthographic_camera>());
    if (ortho_cam.is_valid_interface())
    {
        view_bbox_all_orthographic(world_bbox, viewdir, up, ortho_cam.get());
        return;
    }
}

//----------------------------------------------------------------------
void Examiner_manipulator::view_bbox_all_perspective(
    const mi::math::Bbox<mi::Float32, 3>&   world_bbox,
    const mi::math::Vector<mi::Float32, 3>& viewdir,
    const mi::math::Vector<mi::Float32, 3>& up,
    nv::index::IPerspective_camera*         pcam)
{
    const std::string mn = "Examiner_manipulator::view_bbox_all_perspective: ";

    const mi::math::Vector<mi::Float32, 3> min    = world_bbox.min;
    const mi::math::Vector<mi::Float32, 3> max    = world_bbox.max;
    const mi::math::Vector<mi::Float32, 3> center = world_bbox.center();
    const mi::Float32 rad = mi::math::euclidean_distance(max, min) / 2.0f;

    mi::Float32 fov_rad_2 = static_cast<mi::Float32>(pcam->get_fov_y_rad() / 2.0);
    assert((0.0 < fov_rad_2) && (fov_rad_2 < (M_PI/2.0)));
    mi::Float32 dist = rad / tan(fov_rad_2);
    // DEBUG_LOG << "rad = " << rad << ", dist = " << dist << ", tan(" << fov_rad_2
    //           << ") = " << tan(fov_rad_2);

    mi::math::Vector<mi::Float32, 3> eyepos = -(dist * viewdir) + center;

    const mi::Float32 near_dist_from_eye = dist - rad;
    const mi::Float32 cam_clip_min = static_cast<mi::Float32>(pcam->get_clip_min());
    if ((near_dist_from_eye > 0.0f) && // not dive in the volume
        (near_dist_from_eye < cam_clip_min))
    {
        INFO_LOG << mn << "clip_min: " << cam_clip_min << " is currently clipped out the front. "
                 << "Too see all the bounding box, you may need to adjust the clip_min value.";
    }

    const mi::Float32 cam_clip_max = static_cast<mi::Float32>(pcam->get_clip_max());
    const mi::Float32 far_point_from_eye = dist + rad;
    if (cam_clip_max < far_point_from_eye)
    {
        INFO_LOG << mn << "clip_max: "<< cam_clip_max << " is currently near the far point. "
                 << "Too see all the bounding box, you may need to adjust the clip_max value.";
    }

    pcam->set_eye_point(eyepos);
    pcam->set_view_direction(viewdir);
    pcam->set_up_direction(up);

    // INFO_LOG << "perspective cam: world_bbox = " << world_bbox << ", center bbox = " << center
    //          << ", min: " << min << ", max: " << max
    //          << ", radius: "    << rad
    //          << ", tan("        << fov_rad_2
    //          << ") = "          << tan(fov_rad_2)
    //          << ", distance: "  << dist
    //     ;
}


//----------------------------------------------------------------------
void Examiner_manipulator::view_bbox_all_orthographic(
    const mi::math::Bbox< mi::Float32, 3 >& world_bbox,
    const mi::math::Vector<mi::Float32, 3>& viewdir,
    const mi::math::Vector<mi::Float32, 3>& up,
    nv::index::IOrthographic_camera*        ocam)
{
    const std::string mn = "Examiner_manipulator::view_bbox_all_orthographic: ";

    const mi::math::Vector<mi::Float32, 3> min    = world_bbox.min;
    const mi::math::Vector<mi::Float32, 3> max    = world_bbox.max;
    const mi::math::Vector<mi::Float32, 3> center = world_bbox.center();
    const mi::Float32 rad = mi::math::euclidean_distance(max, min) / 2.0f;
    const mi::Float32 dist = rad;
    const mi::math::Vector<mi::Float32, 3> eyepos = -(dist * viewdir) + center;

    const mi::Float32 near_dist_from_eye = dist - rad;
    const mi::Float32 cam_clip_min = static_cast<mi::Float32>(ocam->get_clip_min());
    if ((near_dist_from_eye > 0.0f) && // not dive in the volume
        (near_dist_from_eye < cam_clip_min))
    {
        INFO_LOG << mn << "clip_min: " << cam_clip_min << " is currently clipped out the front. "
                 << "Too see all the bounding box, you may need to adjust the clip_min value.";
    }

    const mi::Float32 cam_clip_max = static_cast<mi::Float32>(ocam->get_clip_max());
    const mi::Float32 far_point_from_eye = (dist + rad);
    if (cam_clip_max < far_point_from_eye)
    {
        INFO_LOG << mn << "clip_max: "<< cam_clip_max << " is currently near the far point. "
                 << "Too see all the bounding box, you may need to adjust the clip_max value.";
    }

    ocam->set_eye_point(eyepos);
    ocam->set_view_direction(viewdir);
    ocam->set_up_direction(up);
    const mi::Float64 aperture = 2.0f * rad; // diameter of the bbox
    ocam->set_aperture(aperture);
}

//----------------------------------------------------------------------
void Examiner_manipulator::set_predefined_view_parameter_to_camera_sub(
    const mi::math::Bbox<mi::Float32, 3>& world_bbox,
    const Camera_parameter&               cam_param,
    nv::index::ICamera*                   cam)
{
    assert(cam != 0);
    cam_param.copy_this_parameter_to_camera(cam);

    // view all?
    bool is_view_all  = false;
    bool is_valid_key = cam_param.get_bool_value("view_all", is_view_all);
    assert(is_valid_key);
    nv::index_common::no_unused_variable_warning_please(is_valid_key);
    if (is_view_all)
    {
        mi::math::Vector<mi::Float32, 3> predef_viewdir;
        is_valid_key = cam_param.get_vector_value("dir", predef_viewdir);
        assert(is_valid_key);

        mi::math::Vector<mi::Float32, 3> predef_up;
        is_valid_key = cam_param.get_vector_value("up", predef_up);
        assert(is_valid_key);

        view_bbox_all(world_bbox, predef_viewdir, predef_up, cam);
    }
}

//----------------------------------------------------------------------
void Examiner_manipulator::pan_screen_match_camera(
    const nv::index::IScene* scene,
    mi::Float64              last_nx,
    mi::Float64              last_ny,
    mi::Float64              cur_nx,
    mi::Float64              cur_ny,
    nv::index::ICamera*      cam)
{
    assert(cam      != 0);
    assert(scene    != 0);

    // screen pan mode: note not for the multiview since the canvas size is used.
    // FIXME: accessing app_prj in a rendering loop should be avoided.
    nv::index_common::String_dict* app_prj = Nvindex_AppData::instance()->peek_app_proj();
    const mi::math::Vector<mi::Sint32, 2> canvas_size =
        nv::index_common::get_vec_sint32_2(app_prj->get("index::canvas_resolution", "1024 1024"));
    const mi::math::Matrix<mi::Float32, 4, 4> scene_transform  = scene->get_transform_matrix();
    const mi::math::Matrix<mi::Float64, 4, 4> scene_transform_d(scene_transform);
    const mi::math::Matrix<mi::Float64, 4, 4> object_to_screen = 
        Camera_tool::get_mv_matrix(cam, scene_transform, canvas_size);
    // INFO_LOG << "screen pan mode: transform:\n" << scene_transform << "\nobject to screen:\n" << object_to_screen;

    mi::math::Matrix<mi::Float64, 4, 4> screen_to_object(object_to_screen);
    if (!screen_to_object.invert())
    {
        WARN_LOG << "Cannot invert object to screen matrix";
    }

    // Use ROI for the bbox
    const mi::math::Bbox<mi::Float32, 3> roi_bbox = scene->get_clipped_bounding_box();
    const mi::math::Vector<mi::Float32, 3> scene_center_f = roi_bbox.center();
    // INFO_LOG << "scene bbox: " << m_scene_bbox << ", center: " << scene_center_f;

    const mi::math::Vector<mi::Float32, 3> up_dir   = cam->get_up_direction();
    const mi::math::Vector<mi::Float32, 3> view_dir = cam->get_view_direction();
    mi::math::Vector<mi::Float32, 3> ev[3];
    Camera_tool::compute_camera_basis(view_dir, up_dir, ev[0], ev[1], ev[2]);

    // reference point
    mi::math::Vector<mi::Float64, 3> scene_center_d(static_cast<mi::Float64>(scene_center_f.x),
                                                    static_cast<mi::Float64>(scene_center_f.y),
                                                    static_cast<mi::Float64>(scene_center_f.z));

    bool is_adjusted = false;
    if (is_pan_reference_point_use_pick_point())
    {
        if (is_pan_reference_point_valid())
        {
            scene_center_d = mi::math::Vector<mi::Float64, 3>(get_pan_reference_point_object_space());
            is_adjusted = true;
        }
    }

    if ((!is_adjusted) && is_pan_reference_point_roi_bisect_and_solid())
    {
        const mi::Float64 bbox_rad_d = 0.5 * static_cast<mi::Float64>(roi_bbox.diagonal_length());
        const mi::math::Vector<mi::Float64, 3> new_scene_center_d =
            adjust_reference_point_by_bisect_for_pan_move(scene_center_d, bbox_rad_d, scene_transform_d, cam);
        if (scene_center_d != new_scene_center_d)
        {
            scene_center_d = new_scene_center_d;
            // WARN_LOG << "Updated: " << scene_center_d << " -> " << new_scene_center_d;
        }
        is_adjusted = true;
    }

    const mi::math::Vector<mi::Float64, 3> center_screen = mi::math::transform_point(object_to_screen, scene_center_d);
    const mi::math::Vector<mi::Float64, 3> center_world  = mi::math::transform_point(scene_transform,  scene_center_d);

    // on screen, get +x and +y
    const mi::math::Vector<mi::Float64, 3> center_px_screen = center_screen + mi::math::Vector<mi::Float64, 3>(1.0, 0.0, 0.0);
    const mi::math::Vector<mi::Float64, 3> center_py_screen = center_screen + mi::math::Vector<mi::Float64, 3>(0.0, 1.0, 0.0);

    // screen to object to world
    const mi::math::Matrix<mi::Float64, 4, 4> screen_to_world = screen_to_object * scene_transform_d;
    const mi::math::Vector<mi::Float64, 3> center_px_world = mi::math::transform_point(screen_to_world, center_px_screen);
    const mi::math::Vector<mi::Float64, 3> center_py_world = mi::math::transform_point(screen_to_world, center_py_screen);

    // get dx and dy in the world space
    const mi::math::Vector<mi::Float64, 3> dx_world = center_px_world - center_world;
    const mi::math::Vector<mi::Float64, 3> dy_world = center_py_world - center_world;

    // pan move by pixels
    const mi::Float64 x_dist_pixel = (last_nx - cur_nx) * static_cast<mi::Float64>(canvas_size.x);
    const mi::Float64 y_dist_pixel = (last_ny - cur_ny) * static_cast<mi::Float64>(canvas_size.y);

    // move the camera in the world space
    const mi::math::Vector<mi::Float64, 3> delta_world_d = x_dist_pixel * dx_world + y_dist_pixel * dy_world;
    const mi::math::Vector<mi::Float32, 3> delta_world_f(static_cast<mi::Float32>(delta_world_d.x),
                                                         static_cast<mi::Float32>(delta_world_d.y),
                                                         static_cast<mi::Float32>(delta_world_d.z));

    mi::math::Vector<mi::Float32, 3> eye_point = cam->get_eye_point();
    eye_point += delta_world_f;
    cam->set_eye_point(eye_point);
}

//----------------------------------------------------------------------
