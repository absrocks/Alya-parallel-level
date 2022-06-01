/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "heightfield_workflow_functionality.h"

#include <mi/math/vector.h>
#include <mi/math/matrix.h>
#include <mi/base/types.h>

#include <nv/index/iheightfield_interaction.h>
#include <nv/index/isession.h>
#include <nv/index/iregular_heightfield.h>
#include <nv/index/idistributed_compute_algorithm.h>

#include <cassert>

#include "common/common_utility.h"
#include "common/distributed_heightfield_elevation_change.h"
#include "common/forwarding_logger.h"

#include "autotracking_workflow_operation.h"
#include "distributed_heightfield_elevation_delete.h"
#include "simple_gridding_operation.h"
#include "utilities.h"

namespace
{
//----------------------------------------------------------------------
static inline mi::Float32 ccw_triangle_size(
    const mi::math::Vector<mi::Float32, 3>& p0,
    const mi::math::Vector<mi::Float32, 3>& p1,
    const mi::math::Vector<mi::Float32, 3>& p2)
{
    return (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p2[0] - p0[0]) * (p1[1] - p0[1]);
}

//----------------------------------------------------------------------
class Angular_sort_criteria
{
public:
    Angular_sort_criteria(
        const mi::math::Vector<mi::Float32, 3>& p) : m_p(p)
    {
        // empty
    }

    mi::Sint32 operator()(
        const mi::math::Vector<mi::Float32, 3>& p0,
        const mi::math::Vector<mi::Float32, 3>& p1)
    {
        const mi::Float32 is_left_value = ccw_triangle_size(m_p, p0, p1);
        return static_cast<mi::Sint32>(is_left_value>0);
    }

private:
    mi::math::Vector<mi::Float32, 3> m_p;
};

//----------------------------------------------------------------------
/// compute the 2D convex hull of the vertices using the 'Graham Scan' algorithm
void compute_convex_hull(
    const std::vector<mi::math::Vector<mi::Float32, 3> >&   polygon_in,
    std::vector<mi::math::Vector<mi::Float32, 3> >&         convex_hull_verts)
{
    std::vector<mi::math::Vector<mi::Float32, 3> > polygon = polygon_in;
    const mi::Uint32 nb_vertices = mi::Uint32(polygon.size());

    /// Compute initial vertex of the convex hull which is the rightmost vertex
    mi::math::Vector<mi::Float32, 3> initial_hull_value = polygon[0];

    mi::Uint32 initial_idx = 0;
    for (mi::Uint32 i = 1; i < nb_vertices; ++i)
    {
        if(polygon[i][0] > initial_hull_value[0])   /// rightmost vertex check
        {
            initial_hull_value = mi::math::Vector<mi::Float32, 3>(
                polygon[i][0], polygon[i][1], polygon[i][2]);
            initial_idx = i;
        }
    }
    polygon[initial_idx] = polygon[0];
    polygon[0] = initial_hull_value;

    /// Initialize stack
    convex_hull_verts.push_back(initial_hull_value);
    std::sort(
        polygon.begin()+1,
        polygon.end(),
        Angular_sort_criteria(initial_hull_value));
    convex_hull_verts.push_back(polygon[1]);

    /// Compute remaining convex hull values
    for (mi::Uint32 i = 2; i < nb_vertices; )
    {
        const mi::math::Vector<mi::Float32, 3> pi = polygon[i];

        const mi::Uint32 stack_size = mi::Uint32(convex_hull_verts.size());
        assert(stack_size >= 2);
        const mi::math::Vector<mi::Float32, 3> p1 = convex_hull_verts[stack_size-2];
        const mi::math::Vector<mi::Float32, 3> p2 = convex_hull_verts[stack_size-1];

        const mi::Float32 left = ccw_triangle_size(p1, p2, pi);
        if (left > 0.0f)
        {
            convex_hull_verts.push_back(pi);    /// push on stack
            ++i;
        }
        else if (left == 0.0f)
        {
            // degenerated, ignored.
            WARN_LOG << "Degenerated case detected. The result may be wrong.";
            ++i;
        }
        else
        {
            convex_hull_verts.pop_back();  /// pop last element from stack
        }
    }
}
//----------------------------------------------------------------------
}   // anonymous namespace

Heightfield_workflow_functionality* Heightfield_workflow_functionality::g_heightfield_workflow_functionality = 0;

//----------------------------------------------------------------------
void Heightfield_workflow_functionality::clear()
{
    m_current_workflow = NONE;
    m_heightfield_tag  = mi::neuraylib::NULL_TAG;
    m_volume_tag       = mi::neuraylib::NULL_TAG;
    m_bounding_polygon.clear();
    m_associated_volume_tag = mi::neuraylib::NULL_TAG;
}

//----------------------------------------------------------------------
void Heightfield_workflow_functionality::set_workflow(Workflow workflow_type)
{
    m_current_workflow = workflow_type;
    m_heightfield_tag = mi::neuraylib::NULL_TAG;
    m_bounding_polygon.clear();
}

//----------------------------------------------------------------------
std::string Heightfield_workflow_functionality::get_next_step() const
{
    if (m_current_workflow == NONE)
    {
        return "Choose a workflow.";
    }
    else if (!m_heightfield_tag.is_valid())
    {
        return "Choose a heightfield to operate on.";
    }

    if (m_current_workflow == PICK_OPERATION)
    {
        return "Manually pick a slice.\nIf the first intersection represents a slice then "
               "the ray/slice intersection will be used to compute a pick location based "
               "on the build-in user-defined snapping tool.";
    }
    else if(    m_current_workflow == POLYGON_DELETION
             || m_current_workflow == ELEVATION_CHANGE
             || m_current_workflow == GRIDDING)
    {
        return "'Pick' or 'mouse over' on a heightfield.\nIf the first intersection represents a heightfield then "
               "the ray/slice intersection will be used to compute a intersection point "
               "and appends the intersection point to the polygon's vertex list.";
    }
    else if (Heightfield_workflow_functionality::AUTOTRACKING)
    {
        if (!m_associated_volume_tag.is_valid())
        {
            return "Choose a volume scene element.";
        }
    }

    return "Done";
}

//----------------------------------------------------------------------
void Heightfield_workflow_functionality::process_pick(
    const mi::math::Vector_struct<mi::Uint32, 2>&   pixel_position,
    const nv::index::IIndex_canvas*                 pick_canvas,
    Snapping_tool*                                  snapping_tool,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{
    assert(this->is_heightfield_workflow_manual_pick()); // check the workflow is in the manual pick state.

    mi::base::Handle<nv::index::IHeightfield_interaction> heightfield_interaction(
        m_iindex_session->create_heightfield_interaction_interface(
            m_heightfield_tag, m_session_tag, dice_transaction));
    assert(heightfield_interaction != 0);

    heightfield_interaction->pick(pixel_position, pick_canvas, m_volume_tag, snapping_tool, dice_transaction);
}

bool Heightfield_workflow_functionality::is_heightfield_workflow_manual_pick() const
{
    if ((this->get_heightfield().is_valid()) &&
        (this->get_workflow() == Heightfield_workflow_functionality::PICK_OPERATION))
    {
        return true;
    }

    return false;
}

//----------------------------------------------------------------------
bool Heightfield_workflow_functionality::add_vertex_to_bounding_polygon(
    const mi::math::Vector<mi::Float32, 3>& ijk_coord)
{
    if (m_bounding_polygon.empty())
    {
        m_bounding_polygon.push_back(ijk_coord);
        return true;
    }

    if (m_bounding_polygon.back() == ijk_coord)
    {
        // WARN_LOG << "Rejected to add the same point.";
        return false;
    }

    m_bounding_polygon.push_back(ijk_coord);
    return true;
}

const std::vector<mi::math::Vector<mi::Float32, 3> >& Heightfield_workflow_functionality::get_bounding_polygon() const
{
    return m_bounding_polygon;
}

void Heightfield_workflow_functionality::get_bounding_polygon_convex_hull(
    std::vector<mi::math::Vector<mi::Float32, 3> >& convex_hull) const
{
    assert(     Heightfield_workflow_functionality::instance()->is_heightfield_workflow_delete_operation()
             || Heightfield_workflow_functionality::instance()->is_heightfield_workflow_elevation_change_operation()
             || Heightfield_workflow_functionality::instance()->is_heightfield_workflow_gridding_operation() );

    if (m_bounding_polygon.size() > 2)
    {
        compute_convex_hull(m_bounding_polygon, convex_hull);
    }
}

void Heightfield_workflow_functionality::delete_polygon(
    bool                              use_convex_hull,
    mi::neuraylib::IDice_transaction* dice_transaction) const
{
    assert(this->is_heightfield_workflow_delete_operation()); // check the workflow is in the polygon deletion state.
    if (m_bounding_polygon.size() <= 2)
    {
        ERROR_LOG << "Not enough points for the polygon.";
        return;
    }

    mi::base::Handle<nv::index::IHeightfield_interaction> heightfield_interaction(
        m_iindex_session->create_heightfield_interaction_interface(
            m_heightfield_tag, m_session_tag, dice_transaction));
    assert(heightfield_interaction.is_valid_interface());

    std::vector<mi::math::Vector<mi::Float32, 3> > polygon;
    if (use_convex_hull)
    {
        compute_convex_hull(m_bounding_polygon, polygon);
    }
    else
    {
        const mi::Uint32 nb_polygon_vertices = m_bounding_polygon.size();
        for (mi::Uint32 i = 0; i < nb_polygon_vertices; ++i)
        {
            polygon.push_back(m_bounding_polygon[i]);
        }
    }
    const mi::Uint32 nb_vertices = polygon.size();

    std::vector<mi::math::Vector_struct<mi::Float32, 2> > polygon_2d;
    for (mi::Uint32 i = 0; i < nb_vertices; ++i)
    {
        const mi::math::Vector<mi::Float32, 3>& vertex = polygon[i];
        mi::math::Vector_struct<mi::Float32, 2> vertex_2d;
        vertex_2d.x = vertex.x;
        vertex_2d.y = vertex.y;
        polygon_2d.push_back(vertex_2d);
    }

    // Access session
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<const nv::index::ISession>(m_session_tag));

    // Query the heightfield extent/bounds
    mi::base::Handle<const nv::index::IRegular_heightfield> heightfield(
        dice_transaction->access<const nv::index::IRegular_heightfield>(m_heightfield_tag));
    const mi::math::Bbox_struct<mi::Float32, 3> heightfield_bounds = heightfield->get_IJK_bounding_box();
    mi::math::Bbox_struct<mi::Uint32, 2> heightfield_area;
    heightfield_area.min.x = static_cast< mi::Uint32 >(heightfield_bounds.min.x);
    heightfield_area.min.y = static_cast< mi::Uint32 >(heightfield_bounds.min.y);
    heightfield_area.max.x = static_cast< mi::Uint32 >(heightfield_bounds.max.x);
    heightfield_area.max.y = static_cast< mi::Uint32 >(heightfield_bounds.max.y);

    // Access the distribution scheme
    const mi::neuraylib::Tag& dist_layout_tag = session->get_distribution_layout();
    mi::base::Handle<const nv::index::IData_distribution> distribution_layout(
        dice_transaction->access<const nv::index::IData_distribution>(dist_layout_tag));
    assert(distribution_layout.is_valid_interface());

    // heightfield distribution layout
    nv::index::IRegular_heightfield_data_locality* data_locality =
        distribution_layout->retrieve_data_locality_for_editing(
            m_heightfield_tag,
            heightfield_area,
            dice_transaction);

    std::vector<mi::Uint32> cluster_host_ids;
    const mi::Uint32 nb_cluster_hosts = data_locality->get_nb_cluster_nodes();
    assert(nb_cluster_hosts > 0); // at least one host know the heightfield.
    
    for (mi::Uint32 i = 0; i < nb_cluster_hosts; ++i)
    {
        cluster_host_ids.push_back(data_locality->get_cluster_node(i));
    }

    // Choose a distributed computing algorithm
    mi::base::Handle< nv::index::IDistributed_compute_algorithm >
        distributed_compute_algorithm(new Distributed_heightfield_elevation_delete(
                                          m_session_tag,
                                          m_heightfield_tag,
                                          polygon_2d,
                                          cluster_host_ids));
    heightfield_interaction->invoke_computing(distributed_compute_algorithm.get(), dice_transaction);
}

//----------------------------------------------------------------------
bool Heightfield_workflow_functionality::is_heightfield_workflow_delete_operation() const
{
    if ((this->get_heightfield().is_valid()) &&
        (this->get_workflow() == Heightfield_workflow_functionality::POLYGON_DELETION))
    {
        return true;
    }
    return false;
}

//----------------------------------------------------------------------
void Heightfield_workflow_functionality::update_heightfield_elevation_values(
    mi::Float32                       new_elevation_value,
    bool                              is_scaling_op,
    bool                              use_convex_hull,
    mi::neuraylib::IDice_transaction* dice_transaction) const
{
    assert(this->is_heightfield_workflow_elevation_change_operation()); // check the workflow is in the elevation change state.

    // heightfield_interaction is new-ed
    mi::base::Handle<nv::index::IHeightfield_interaction>
        heightfield_interaction(m_iindex_session->create_heightfield_interaction_interface(
                                m_heightfield_tag, m_session_tag, dice_transaction));
    assert(heightfield_interaction.is_valid_interface());

    std::vector<mi::math::Vector<mi::Float32, 3> > polygon;
    if (use_convex_hull)
    {
        compute_convex_hull(m_bounding_polygon, polygon);
    }
    else
    {
        const mi::Uint32 nb_polygon_vertices = m_bounding_polygon.size();
        for (mi::Uint32 i = 0; i < nb_polygon_vertices; ++i)
        {
            polygon.push_back(m_bounding_polygon[i]);
        }
    }
    const mi::Uint32 nb_vertices = polygon.size();

    std::vector<mi::math::Vector_struct<mi::Float32, 2> > polygon_2d;
    for (mi::Uint32 i = 0; i < nb_vertices; ++i)
    {
        const mi::math::Vector<mi::Float32, 3>& vertex = polygon[i];
        mi::math::Vector_struct<mi::Float32, 2> vertex_2d;
        vertex_2d.x = vertex.x;
        vertex_2d.y = vertex.y;
        polygon_2d.push_back(vertex_2d);
    }

    // Access session
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<const nv::index::ISession>(m_session_tag));

    // Query the heightfield extent/bounds
    mi::base::Handle<const nv::index::IRegular_heightfield> heightfield(
        dice_transaction->access<const nv::index::IRegular_heightfield>(m_heightfield_tag));
    const mi::math::Bbox_struct<mi::Float32, 3> heightfield_bounds = heightfield->get_IJK_bounding_box();
    mi::math::Bbox_struct<mi::Uint32, 2> heightfield_area;
    heightfield_area.min.x = static_cast< mi::Uint32 >(heightfield_bounds.min.x);
    heightfield_area.min.y = static_cast< mi::Uint32 >(heightfield_bounds.min.y);
    heightfield_area.max.x = static_cast< mi::Uint32 >(heightfield_bounds.max.x);
    heightfield_area.max.y = static_cast< mi::Uint32 >(heightfield_bounds.max.y);

    // Access the distribution scheme
    const mi::neuraylib::Tag& dist_layout_tag = session->get_distribution_layout();
    mi::base::Handle<const nv::index::IData_distribution> distribution_layout(
        dice_transaction->access<const nv::index::IData_distribution>(dist_layout_tag));

    // Heightfield distribution layout
    nv::index::IRegular_heightfield_data_locality* data_locality = 
        distribution_layout->retrieve_data_locality_for_editing(
            m_heightfield_tag,
            heightfield_area,
            dice_transaction);

    std::vector<mi::Uint32> cluster_host_ids;
    const mi::Uint32 nb_cluster_hosts = data_locality->get_nb_cluster_nodes();
    for (mi::Uint32 i = 0; i < nb_cluster_hosts; ++i)
    {
        cluster_host_ids.push_back(data_locality->get_cluster_node(i));
    }

    // Choose a distributed computing algorithm
    mi::base::Handle<nv::index::IDistributed_compute_algorithm>
        distributed_compute_algorithm(new nv::index_common::Distributed_heightfield_elevation_change(
                                          m_session_tag,
                                          m_heightfield_tag,
                                          is_scaling_op,
                                          new_elevation_value,
                                          polygon_2d,
                                          cluster_host_ids));
    heightfield_interaction->invoke_computing(distributed_compute_algorithm.get(), dice_transaction);
}

bool Heightfield_workflow_functionality::is_heightfield_workflow_elevation_change_operation() const
{
    if( (this->get_heightfield().is_valid()) &&
        (this->get_workflow() == Heightfield_workflow_functionality::ELEVATION_CHANGE) )
    {
        return true;
    }
    return false;
}

//----------------------------------------------------------------------
void Heightfield_workflow_functionality::gridding(
    bool                              use_convex_hull,
    mi::neuraylib::IDice_transaction* dice_transaction) const
{
    assert(this->is_heightfield_workflow_gridding_operation()); // check the workflow is in the elevation change state.

    mi::base::Handle<nv::index::IHeightfield_interaction> heightfield_interaction(
        m_iindex_session->create_heightfield_interaction_interface(
            m_heightfield_tag, m_session_tag, dice_transaction));
    assert(heightfield_interaction.is_valid_interface());

    std::vector<mi::math::Vector<mi::Float32, 3> > polygon;
    if(use_convex_hull)
    {
        compute_convex_hull(m_bounding_polygon, polygon);
    }
    else
    {
        const mi::Uint32 nb_polygon_vertices = m_bounding_polygon.size();
        for(mi::Uint32 i = 0; i<nb_polygon_vertices; ++i)
        {
            polygon.push_back(m_bounding_polygon[i]);
        }
    }
    const mi::Uint32 nb_vertices = polygon.size();

    std::vector<mi::math::Vector_struct<mi::Float32, 2> > polygon_2d;
    for (mi::Uint32 i = 0; i < nb_vertices; ++i)
    {
        const mi::math::Vector<mi::Float32, 3>& vertex = polygon[i];
        mi::math::Vector_struct<mi::Float32, 2> vertex_2d;
        vertex_2d.x = vertex.x;
        vertex_2d.y = vertex.y;
        polygon_2d.push_back(vertex_2d);
    }

    // Access session
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<const nv::index::ISession>(m_session_tag));

    // Query the heightfield extent/bounds
    mi::base::Handle<const nv::index::IRegular_heightfield> heightfield(
        dice_transaction->access<const nv::index::IRegular_heightfield>(m_heightfield_tag));
    const mi::math::Bbox_struct<mi::Float32, 3> heightfield_bounds = heightfield->get_IJK_bounding_box();
    mi::math::Bbox_struct<mi::Uint32, 2> heightfield_area;
    heightfield_area.min.x = static_cast< mi::Uint32 >(heightfield_bounds.min.x);
    heightfield_area.min.y = static_cast< mi::Uint32 >(heightfield_bounds.min.y);
    heightfield_area.max.x = static_cast< mi::Uint32 >(heightfield_bounds.max.x);
    heightfield_area.max.y = static_cast< mi::Uint32 >(heightfield_bounds.max.y);

    // Access the distribution scheme
    const mi::neuraylib::Tag& dist_layout_tag = session->get_distribution_layout();
    mi::base::Handle<const nv::index::IData_distribution> distribution_layout(
        dice_transaction->access<const nv::index::IData_distribution>(dist_layout_tag));

    // Heightfield distribution layout
    nv::index::IRegular_heightfield_data_locality* data_locality = 
        distribution_layout->retrieve_data_locality_for_editing(
            m_heightfield_tag,
            heightfield_area,
            dice_transaction);

    std::vector<mi::Uint32> cluster_host_ids;
    const mi::Uint32 nb_cluster_hosts = data_locality->get_nb_cluster_nodes();
    for(mi::Uint32 i = 0; i < nb_cluster_hosts; ++i)
    {
        cluster_host_ids.push_back(data_locality->get_cluster_node(i));
    }

    // Choose a distributed computing algorithm
    mi::base::Handle<Simple_gridding_operation> distributed_compute_algorithm(
        new Simple_gridding_operation(
            m_session_tag, m_heightfield_tag, polygon_2d));

    heightfield_interaction->invoke_computing(distributed_compute_algorithm.get(), dice_transaction);
}

bool Heightfield_workflow_functionality::is_heightfield_workflow_gridding_operation() const
{
    if( (this->get_heightfield().is_valid()) &&
        (this->get_workflow() == Heightfield_workflow_functionality::GRIDDING) )
    {
        return true;
    }
    return false;
}

//----------------------------------------------------------------------
void Heightfield_workflow_functionality::autotracking(
    mi::base::Handle<mi::neuraylib::IScope>&    scope) const
{
    // get dice transaction
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        scope->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());

    mi::base::Handle<Autotracking_workflow_operation> autotracking;
    mi::base::Handle<nv::index::IHeightfield_interaction> heightfield_interaction;
    {
        mi::base::Handle<nv::index::IHeightfield_interaction> hi(
            m_iindex_session->create_heightfield_interaction_interface(
                m_heightfield_tag, m_session_tag, dice_transaction.get()));
        heightfield_interaction.swap(hi);

        // Access session
        mi::base::Handle<const nv::index::ISession> session(
            dice_transaction->access<const nv::index::ISession>(m_session_tag));

        // Query the heightfield extent/bounds
        mi::base::Handle<const nv::index::IRegular_heightfield> heightfield(
            dice_transaction->access<const nv::index::IRegular_heightfield>(m_heightfield_tag));
        const mi::math::Bbox_struct<mi::Float32, 3> heightfield_bounds = heightfield->get_IJK_bounding_box();
        mi::math::Bbox_struct<mi::Uint32, 2> heightfield_area;
        heightfield_area.min.x = static_cast< mi::Uint32 >(heightfield_bounds.min.x);
        heightfield_area.min.y = static_cast< mi::Uint32 >(heightfield_bounds.min.y);
        heightfield_area.max.x = static_cast< mi::Uint32 >(heightfield_bounds.max.x);
        heightfield_area.max.y = static_cast< mi::Uint32 >(heightfield_bounds.max.y);

        // Access the distribution scheme
        const mi::neuraylib::Tag& dist_layout_tag = session->get_distribution_layout();
        mi::base::Handle<const nv::index::IData_distribution> distribution_layout(
            dice_transaction->access<const nv::index::IData_distribution>(dist_layout_tag));

        // Heightfield distribution layout
        nv::index::IRegular_heightfield_data_locality* data_locality =
            distribution_layout->retrieve_data_locality_for_editing(
                m_heightfield_tag,
                heightfield_area,
                dice_transaction.get());

        std::vector<mi::Uint32> cluster_host_ids;
        const mi::Uint32 nb_cluster_hosts = data_locality->get_nb_cluster_nodes();
        for(mi::Uint32 i = 0; i < nb_cluster_hosts; ++i)
        {
            cluster_host_ids.push_back(data_locality->get_cluster_node(i));
        }

        // Choose a distributed computing algorithm
        mi::base::Handle<Autotracking_workflow_operation> at(
            new Autotracking_workflow_operation(
                m_session_tag, m_heightfield_tag, m_associated_volume_tag,
                true, /* 4-way (false) or 8-way (true) seeding*/
                5000    /* max number of iterations */ ));
        autotracking.swap(at);
    }
    // commit changes
    dice_transaction->commit();

    assert(autotracking.is_valid_interface());
    assert(heightfield_interaction.is_valid_interface());

    mi::Uint32 nb_iterations = 1;
    while (autotracking->needs_iteration())
    {
        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            scope->create_transaction<mi::neuraylib::IDice_transaction>());
        assert(dice_transaction.is_valid_interface());
        {
            if ((nb_iterations % 100) == 0)
            {
                autotracking->enable_elevation_updates();
            }
            heightfield_interaction->invoke_computing(autotracking.get(), dice_transaction.get());
        }

        // commit changes
        dice_transaction->commit();

        if ((nb_iterations % 100) == 0)
        {
            autotracking->disable_elevation_updates();
            nv::index_common::sleep(1.0f);   // take a break and let the renderer and compositor do their job .... TODO: remove sleep!
        }
        nb_iterations++;
    }
}

//----------------------------------------------------------------------
bool Heightfield_workflow_functionality::is_autotracking_workflow_operation() const
{
    if( (this->get_heightfield().is_valid()) &&
        (this->m_associated_volume_tag.is_valid()) &&
        (this->get_workflow() == Heightfield_workflow_functionality::AUTOTRACKING) )
    {
        return true;
    }
    return false;
}

//----------------------------------------------------------------------
void Heightfield_workflow_functionality::set_last_xyz_world_pick_position(
    const mi::math::Vector<mi::Float32, 3>& xyz_world)
{
    m_last_xyz_world_pick_position = xyz_world;
}

//----------------------------------------------------------------------
const mi::math::Vector<mi::Float32, 3>& Heightfield_workflow_functionality::get_last_xyz_world_pick_position() const
{
    return m_last_xyz_world_pick_position;
}

//----------------------------------------------------------------------
void Heightfield_workflow_functionality::set_last_xyz_volume_pick_position(
    const mi::math::Vector<mi::Float32, 3>& xyz_volume)
{
    m_last_xyz_volume_pick_position = xyz_volume;
}

//----------------------------------------------------------------------
const mi::math::Vector<mi::Float32, 3>& Heightfield_workflow_functionality::get_last_xyz_volume_pick_position() const
{
    return m_last_xyz_volume_pick_position;
}

//----------------------------------------------------------------------
void Heightfield_workflow_functionality::clear_pick_position()
{
    m_last_xyz_world_pick_position  = mi::math::Vector<mi::Float32, 3>(-1.0f, -1.0f, -1.0f);
    m_last_xyz_volume_pick_position = mi::math::Vector<mi::Float32, 3>(-1.0f, -1.0f, -1.0f);
}

//----------------------------------------------------------------------

