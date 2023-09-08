/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief heightfield workflow example implementation

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_HEIGHTFIELD_WORKFLOW_FUNCTIONALITY_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_HEIGHTFIELD_WORKFLOW_FUNCTIONALITY_H

#include <mi/math/vector.h>
#include <mi/math/matrix.h>
#include <mi/base/types.h>

#include <nv/index/iindex.h>

#include <cassert>

#include "common/forwarding_logger.h"

#include "utilities.h"
#include "snapping_tool.h"

/// IndeX scene workflow functionality.
///
/// Singleton
///
class Heightfield_workflow_functionality
{
public:
    /// Workflows
    enum Workflow {
        NONE                = 0,    ///
        PICK_OPERATION      = 1,    ///
        POLYGON_DELETION    = 2,    ///
        ELEVATION_CHANGE    = 3,    ///
        GRIDDING            = 4,    ///
        AUTOTRACKING        = 5     ///
    };

    /// clear the workflow status
    void clear();

    /// set heightfield workflow type
    void set_workflow(Workflow workflow_type);
    /// get heightfield workflow type
    Workflow get_workflow() const { return m_current_workflow; }

    /// set selected heightfield
    void set_heightfield(const mi::neuraylib::Tag& heightfield_tag)
    {
        m_heightfield_tag = heightfield_tag;
        m_bounding_polygon.clear();
    }
    /// get selected heightfield
    /// \return selected heightfield tag.
    const mi::neuraylib::Tag& get_heightfield() const { return m_heightfield_tag; }

    /// set selected volume
    void set_volume(const mi::neuraylib::Tag& volume_tag)
    {
        m_volume_tag = volume_tag;
    }
    /// get selected volume
    const mi::neuraylib::Tag& get_volume() const { return m_volume_tag; }

    /// Receive the next steps explanations for ui display
    std::string get_next_step() const;

public:
    // Manual Picking Process -----------------------------------------------
    /// Process, i.e., apply functionality
    ///
    /// \param[in] pixel_position   picking pixel position
    /// \param[in] pick_canvas      canvas to pick
    /// \param[in] snapping_tool    pick tool (a reference, this instance is not an owner of this object)
    /// \param[in] dice_transaction db transaction
    void process_pick(
        const mi::math::Vector_struct<mi::Uint32, 2>&   pixel_position,
        const nv::index::IIndex_canvas*                 pick_canvas,        
        Snapping_tool*                                  snapping_tool,
        mi::neuraylib::IDice_transaction*               dice_transaction) const;

    /// is heightfield manual pick ready?
    /// \return true if heightfield manual pick workflow state.
    bool is_heightfield_workflow_manual_pick() const;

public:
    // Polygon Delete Operation ---------------------------------------------

    /// in IJK space of the heightfield.
    ///
    /// There is a simple degenerated case filter. But it is only see
    /// the cases locally. No global degenerated case considered.
    ///
    /// \return true the vertex is added. false the adding vertex is a
    /// degenerated case
    bool add_vertex_to_bounding_polygon(
        const mi::math::Vector<mi::Float32, 3>& ijk_coord);

    /// in IJK space of the heightfield
    const std::vector<mi::math::Vector<mi::Float32, 3> >& get_bounding_polygon() const;
    void get_bounding_polygon_convex_hull(std::vector<mi::math::Vector<mi::Float32, 3> >& convex_hull) const;

    void delete_polygon(
        bool                              use_convex_hull,
        mi::neuraylib::IDice_transaction* dice_transaction) const;

    /// is heightfield delete operation ready?
    /// \return true if heightfield delete operation workflow state.
    bool is_heightfield_workflow_delete_operation() const;

public:
    // Polygon Elevation Change Operation ---------------------------------------------

    /// update heightfield elevation values
    ///
    /// \param[in] new_elevation_value new elevation value (or scaling
    ///            factor based on is_scale_op)
    /// \param[in] is_scale_op         scale the current value
    ///            instead of set the new height value when true.
    /// \param[in] use_convex_hull use convex hull of the polygon when true.
    /// \param[in] dice_transaction db transaction
    void update_heightfield_elevation_values(
        mi::Float32                       new_elevation_value,
        bool                              is_scale_op,
        bool                              use_convex_hull,
        mi::neuraylib::IDice_transaction* dice_transaction) const;

    /// is heightfield elevation change operation ready?
    /// \return true if any heightfield elevation change operation workflow state.
    bool is_heightfield_workflow_elevation_change_operation() const;

public:
    // Gridding Operation ---------------------------------------------
    void gridding(
        bool                              use_convex_hull,
        mi::neuraylib::IDice_transaction* dice_transaction) const;

    /// is heightfield gridding change operation ready?
    /// \return true if any heightfield gridding change operation workflow state.
    bool is_heightfield_workflow_gridding_operation() const;
    // Gridding Operation ---------------------------------------------

public:
    // Autotacking Operation ------------------------------------------
    void autotracking(
        mi::base::Handle<mi::neuraylib::IScope>&    scope) const;

    /// is heightfield autotracking operation ready?
    /// \return true if any heightfield autotracking operation workflow state.
    bool is_autotracking_workflow_operation() const;

    /// Set the associated volume tag.
    void set_associated_scene_element(
        const mi::neuraylib::Tag& associated_volume_tag)
    {
        m_associated_volume_tag = associated_volume_tag;
        m_bounding_polygon.clear();
    }
    /// get associated volume scene element
    /// \return selected volume tag.
    const mi::neuraylib::Tag& get_associated_scene_element() const { return m_associated_volume_tag; }

public:
    //-----------------------------------------------------------------
    // for testing (unit tests)
    //-----------------------------------------------------------------

    /// set the xyz world pick position for test.
    ///
    /// \param[in] xyz_world xyz world coordinate position
    void set_last_xyz_world_pick_position(const mi::math::Vector<mi::Float32, 3>& xyz_world);

    /// get the last input pick position in the workflow for test.
    ///
    /// \return input xyz world pick position
    const mi::math::Vector<mi::Float32, 3>& get_last_xyz_world_pick_position() const;

    /// set the xyz volume pick position.
    ///
    /// \param[in] xyz_volume xyz volume coordinate position
    void set_last_xyz_volume_pick_position(const mi::math::Vector<mi::Float32, 3>& xyz_volume);

    /// get the last result pick position on the volume. If
    /// the pick position can not ideitify the volume data, it
    /// returns (-1,-1,-1).
    ///
    /// \return result xyz pick position of the volume. (-1,-1,-1)
    /// when no volume at the pick psition.
    const mi::math::Vector<mi::Float32, 3>& get_last_xyz_volume_pick_position() const;

    /// clear the pick position
    void clear_pick_position();

public:
    /// Access workflow instance
    static Heightfield_workflow_functionality* instance()
    {
        if (g_heightfield_workflow_functionality == 0)
        {
            ERROR_LOG << "Heightfield workflow has not been init.";
        }
        return g_heightfield_workflow_functionality;
    }

    /// initialize the heightfield workflow
    static void init(
        mi::base::Handle<nv::index::IIndex_session>& iindex_session,
        const mi::neuraylib::Tag&                    session_tag)
    {
        g_heightfield_workflow_functionality =
            new Heightfield_workflow_functionality(iindex_session, session_tag);
    }

    /// delete the heightfield workflow instance
    static void delete_instance()
    {
        if (g_heightfield_workflow_functionality != 0)
        {
            delete g_heightfield_workflow_functionality;
            g_heightfield_workflow_functionality = 0;
        }
    }

private:
    /// Singleton instance
    static Heightfield_workflow_functionality* g_heightfield_workflow_functionality;

    /// Constructor
    Heightfield_workflow_functionality()
    {
        /* - not allowed - */
        abort();
    }

    /// constructor
    Heightfield_workflow_functionality(
        mi::base::Handle<nv::index::IIndex_session>& iindex_session,
        const mi::neuraylib::Tag&                    session_tag)
        :
        m_iindex_session(),
        m_session_tag(session_tag),
        m_current_workflow(NONE),
        m_heightfield_tag(mi::neuraylib::NULL_TAG),
        m_volume_tag     (mi::neuraylib::NULL_TAG),
        m_snapping_tool(0),
        m_bounding_polygon(),
        m_associated_volume_tag(mi::neuraylib::NULL_TAG)
    {
        m_iindex_session = iindex_session; // ref count up
    }

    /// destructor
    ~Heightfield_workflow_functionality()
    {
        m_iindex_session        = 0; // ref count down
        m_session_tag           = mi::neuraylib::NULL_TAG;
        m_current_workflow      = NONE;
        m_heightfield_tag       = mi::neuraylib::NULL_TAG;
        m_volume_tag            = mi::neuraylib::NULL_TAG;
        m_snapping_tool         = 0; // This is not a owner, just dereference
        // m_bounding_polygon
        m_associated_volume_tag = mi::neuraylib::NULL_TAG;
    }

private:
    mi::base::Handle<nv::index::IIndex_session>             m_iindex_session;
    mi::neuraylib::Tag                                      m_session_tag;
    Workflow                                                m_current_workflow;

    // Workflow details
    mi::neuraylib::Tag                                      m_heightfield_tag;  // heightfield to operate on
    mi::neuraylib::Tag                                      m_volume_tag;       // volume to operate on

    // Pick Operatioon
    Snapping_tool*                                          m_snapping_tool;

    // Polygon for deleting elevation values or changing their heights
    std::vector<mi::math::Vector<mi::Float32, 3> >          m_bounding_polygon;  // should be non-self-intersecting - will not be tested!

    // Assigned associated volume for autotracking
    mi::neuraylib::Tag                                      m_associated_volume_tag;


    //----------------------------------------------------------------------
    // for test
    /// the last pick position on the world coordinate.
    mi::math::Vector<mi::Float32, 3> m_last_xyz_world_pick_position;
    /// the last result pick position on the volume.
    mi::math::Vector<mi::Float32, 3> m_last_xyz_volume_pick_position;

};

#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_HEIGHTFIELD_WORKFLOW_FUNCTIONALITY_H
