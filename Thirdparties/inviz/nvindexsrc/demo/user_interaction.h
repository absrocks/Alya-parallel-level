/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief common user interaction code

#ifndef NVIDIA_INDEX_BIN_USER_INTERACTION_H
#define NVIDIA_INDEX_BIN_USER_INTERACTION_H

#include <mi/dice.h>

#include "examiner_manipulator.h"
#include "multiple_camera.h"

class Nvindex_rendering_context;

/// User interaction handling class
///
/// RTMP, dice bridge, HTML5 have own command stream, but at the end,
/// mouse motion and so on are the same inetraction. This class handle
/// such part.
///
class User_interaction
{
public:
    /// default constructor
    User_interaction();
    /// destructor
    ~User_interaction();

    /// initialize by dict (convenient method)
    ///
    /// \param[in] dict initialize this object (including the examiner)
    void initialize_by_dict(const nv::index_common::String_dict& dict);

    /// Get this user interaction's examiner
    ///
    /// \return an examiner 
    Examiner_manipulator* get_examiner();

    /// Get this user interaction's multiple stereo camera
    ///
    /// \return a multiple stereo camera
    Multiple_stereo_camera* get_multiple_stereo_camera();

    /// Handle the mouse motion
    ///
    /// \param[in] button           mouse button (left, middle, right)
    /// \param[in] mouse_x          mouse pixel position X
    /// \param[in] mouse_y          mouse pixel position Y
    /// \param[in] scene_tag        scene tag
    /// \param[in] dice_transaction dice transaction
    void mouse_motion(
        mi::Sint32                        button,
        mi::Sint32                        mouse_x,
        mi::Sint32                        mouse_y,
        const mi::neuraylib::Tag&         scene_tag,
        mi::neuraylib::IDice_transaction* dice_transaction);

    /// Handle the mouse button action (press, release)
    ///
    /// \param[in] button       mouse button (left, middle, right)
    /// \param[in] mouse_action mouse action (press, release)
    /// \param[in] mouse_x      mouse pixel position X
    /// \param[in] mouse_y      mouse pixel position Y
    /// \param[in] dice_transaction  dice transaction
    /// \return the pick position on canvas, no pick if negative pick
    /// position coordinate
    mi::math::Vector<mi::Sint32, 2> mouse_button_action(
        mi::Sint32                         button,
        Examiner_manipulator::Mouse_action mouse_action,
        mi::Sint32                         mouse_x,
        mi::Sint32                         mouse_y,
        mi::neuraylib::IDice_transaction*  dice_transaction);

    /// Handle pick command issue 
    ///
    /// \param[in] irc_ref            index rendering context reference
    /// \param[in] is_multi_view_mode true when multi-view pick
    /// \param[in] pick_position      pick position on the canvas
    /// \return Tag of the first picked scene element, if any.
    mi::neuraylib::Tag issue_pick_command(
        Nvindex_rendering_context*             irc_ref,
        bool                                   is_multi_view_mode,
        const mi::math::Vector<mi::Sint32, 2>& pick_position);

    /// Handle mouse wheel 
    ///
    /// \param[in] irc_ref   index rendering context reference
    /// \param[in] scene_tag scene tag
    /// \param[in] mouse_x   current mouse position x
    /// \param[in] mouse_y   current mouse position y
    /// \param[in] delta     wheel delta (+1 or -1)
    /// \param[in] dice_transaction  dice transaction
    void mouse_wheel_action(
        Nvindex_rendering_context*        irc_ref,
        const mi::neuraylib::Tag&         scene_tag,
        mi::Sint32                        mouse_x,
        mi::Sint32                        mouse_y,
        mi::Sint32                        delta,
        mi::neuraylib::IDice_transaction* dice_transaction);
    
    /// Handle key down action
    ///
    /// This is called when pushed:
    ///  - a key with a special key
    ///  - only a special key
    ///
    /// \param[in] irc_ref      index rendering context reference
    /// \param[in] keycode      keycode
    /// \param[in] is_shift_on  true when shift is on
    /// \param[in] is_ctrl_on   true when ctrl is on
    /// \param[in] is_alt_on    true when alt is on
    /// \param[in] is_meta_on   true when meta is on
    ///
    void key_down_action(
        Nvindex_rendering_context* irc_ref,
        const mi::Sint32           keycode, 
        const bool                 is_shift_on,
        const bool                 is_ctrl_on,
        const bool                 is_alt_on,
        const bool                 is_meta_on);

    /// Handle key press action (an ascii key is pushed)
    ///
    /// This is called when pushed:
    ///  - A key without a special key
    ///
    /// \param[in] irc_ref      index rendering context reference
    /// \param[in] keycode      keycode
    /// \param[in] is_shift_on  true when shift is on
    /// \param[in] is_ctrl_on   true when ctrl is on
    /// \param[in] is_alt_on    true when alt is on
    /// \param[in] is_meta_on   true when meta is on
    ///
    /// \return true when the key was suported, false not supported
    /// yet. This is for rtmp_call_event_handler handling. This
    /// key_down_action doesn't support a specific key set (e.g.,
    /// viewing_scenario), thus this returns false and need to be
    /// handled in the rtmp_call_event_handler.
    ///
    bool key_press_action(
        Nvindex_rendering_context* irc_ref,
        const mi::Sint32           keycode, 
        const bool                 is_shift_on,
        const bool                 is_ctrl_on,
        const bool                 is_alt_on,
        const bool                 is_meta_on);

    /// Handle key up action (any key)
    ///
    /// \param[in] irc_ref      index rendering context reference
    /// \param[in] keycode      keycode
    /// \param[in] is_shift_on  true when shift is on
    /// \param[in] is_ctrl_on   true when ctrl is on
    /// \param[in] is_alt_on    true when alt is on
    /// \param[in] is_meta_on   true when meta is on
    ///
    void key_up_action(
        Nvindex_rendering_context* irc_ref,
        const mi::Sint32           keycode, 
        const bool                 is_shift_on,
        const bool                 is_ctrl_on,
        const bool                 is_alt_on,
        const bool                 is_meta_on);

    /// Handle video resize action
    ///
    /// \param[in] irc_ref      index rendering context reference
    /// \param[in] video_width  new video width
    /// \param[in] video_height new video height
    ///
    void video_resize_action(
        Nvindex_rendering_context* irc_ref,
        const mi::Sint32           video_width, 
        const mi::Sint32           video_height);

    /// Handle shutdown action
    ///
    /// \param[in] irc_ref      index rendering context reference
    ///
    void shutdown_server_action(
        Nvindex_rendering_context* irc_ref);


private:
    /// Set the rotation center to the volume (currently, first volume
    /// in the scene) in ROI.
    ///
    /// \param[in] irc_ref index rendering context reference
    void set_rotation_center_to_volume_roi_center(
        Nvindex_rendering_context* irc_ref);

    /// camera zoom in/out to the reference point
    /// 
    /// \param[in] irc_ref    index rendering context reference
    /// \param[in] is_zoom_in true when zoom in, false when zoom out.
    void camera_zoom_in_reference_point(
        Nvindex_rendering_context* irc_ref,
        bool                       is_zoom_in);

private:
    /// last mouse position X
    mi::Sint32             m_last_press_mouse_x;
    /// last mouse position Y
    mi::Sint32             m_last_press_mouse_y;
    /// This is the user interaction's examiner
    Examiner_manipulator   m_examiner;
    /// Multiple stereo cameras for this user interaction
    Multiple_stereo_camera m_multiple_camera;
};

#endif // NVIDIA_INDEX_BIN_USER_INTERACTION_H
