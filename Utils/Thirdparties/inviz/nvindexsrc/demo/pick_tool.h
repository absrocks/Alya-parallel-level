/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Picking tool

#ifndef NVIDIA_INDEX_BIN_PICK_TOOL_H
#define NVIDIA_INDEX_BIN_PICK_TOOL_H

#include <mi/dice.h>

#include <nv/index/iscene_query_results.h>

#include <string>

class Nvindex_rendering_context;
class Examiner_manipulator;

/// Pick tool
///
/// When a user pick on a canvas, this tool shows a sphere on the
/// picking point. Also this tool show the picking information.
///
class Pick_tool
{
public:
    /// default constructor
    Pick_tool();
    /// destructor
    ~Pick_tool();

    /// Issue pick command.
    ///
    ///  - Put a sphere at the picked point
    ///  - Print picked scene element tag and pick information
    ///
    /// \param[in] irc_ref                 index rendering context       
    /// \param[in] is_multi_view_mode      true when multiview mode
    /// \param[in] pick_location_on_canvas pick location on canvas
    /// \param[in] examiner                The examiner operated by this picking
    ///
    /// \return Tag of the first picked scene element, if any.
    ///
    static mi::neuraylib::Tag issue_pick_command(
        Nvindex_rendering_context*             irc_ref,
        bool                                   is_multi_view_mode,
        const mi::math::Vector<mi::Sint32, 2>& pick_location_on_canvas,
        Examiner_manipulator*                  examiner);

private:
    /// Get the pick point in 3D space
    static mi::math::Vector<mi::Float32, 3> get_pick_point_xyz(
        nv::index::IScene_pick_result* pick_result);

    /// Get the pick scene element tag prefix for scene editor
    /// visibility control
    static std::string get_pick_scene_element_tag_prefix();

    /// Create a pick point as a sphere
    static void create_pick_point(
        nv::index::IScene_pick_result*    pick_result,
        const mi::neuraylib::Tag&         session_tag,
        mi::neuraylib::IDice_transaction* dice_transaction);

    /// Remove a pick point sphere
    static void remove_pick_point(
        const mi::neuraylib::Tag&         session_tag,
        mi::neuraylib::IDice_transaction* dice_transaction);

    /// Print pick result to console
    static void print_pick_result(
        Nvindex_rendering_context*        irc_ref,
        mi::Uint32                        intersect_num,
        nv::index::IScene_pick_result*    pick_result,
        mi::neuraylib::IDice_transaction* dice_transaction);

};

#endif // NVIDIA_INDEX_BIN_PICK_TOOL_H
