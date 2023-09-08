/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "autotracking_workflow_operation.h"

#include <cassert>

#include <nv/index/isession.h>
#include <nv/index/iregular_heightfield.h>
#include <nv/index/iregular_heightfield_compute_task.h>
#include <nv/index/idistributed_data_access.h>

#include "distributed_heightfield_elevation_delete.h"
#include "utilities.h"

#include "common/forwarding_logger.h"

namespace
{

static inline void print_seeds(
    const std::vector<mi::math::Vector_struct<mi::Float32, 3> >& seed_points)
{
    const mi::Uint32 nb_seed_points = seed_points.size();
    for(mi::Uint32 i=0; i<nb_seed_points; ++i)
    {
        INFO_LOG << "Seed point (" << i << ")" << " = " << seed_points[i];
    }
}

// --------------------------------------------------------------------------------------------
static inline void determine_seed_points(
    const mi::neuraylib::Tag&                               session_tag,
    const mi::neuraylib::Tag&                               heightfield_tag,
    mi::neuraylib::IDice_transaction*                       dice_transaction,
    std::vector<mi::math::Vector_struct<mi::Float32, 3> >&  seed_points)
{
    // Access session
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<const nv::index::ISession>(session_tag));
    // const mi::neuraylib::Tag& dist_layout_tag = session->get_distribution_layout();

    // Query the heightfield's seed points
    mi::base::Handle<const nv::index::IRegular_heightfield> heightfield(
        dice_transaction->access<const nv::index::IRegular_heightfield>(heightfield_tag));
    const mi::Uint32 nb_seed_points = heightfield->get_nb_seed_points();
    seed_points.resize(nb_seed_points);
    for(mi::Uint32 i=0; i<nb_seed_points; ++i)
    {
        heightfield->get_seed_point(i, seed_points[i]);
        DEBUG_LOG << "Seed point (" << i << ")" << " = " << seed_points[i];
    }
}

// --------------------------------------------------------------------------------------------
static inline nv::index::IRegular_heightfield_data_locality* compute_distribution_layout(
    const mi::neuraylib::Tag&          session_tag,
    const mi::neuraylib::Tag&          heightfield_tag,
    mi::neuraylib::IDice_transaction*  dice_transaction,
    bool                               editing)
{
    // Access session
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<const nv::index::ISession>(session_tag));
    const mi::neuraylib::Tag& dist_layout_tag = session->get_distribution_layout();

    // Query the heightfield extent/bounds
    mi::base::Handle<const nv::index::IRegular_heightfield> heightfield(
        dice_transaction->access<const nv::index::IRegular_heightfield>(heightfield_tag));
    const mi::math::Bbox_struct<mi::Float32, 3> heightfield_bounds = heightfield->get_IJK_bounding_box();
    mi::math::Bbox_struct<mi::Uint32, 2> heightfield_area;
    heightfield_area.min.x = static_cast< mi::Uint32 >(heightfield_bounds.min.x);
    heightfield_area.min.y = static_cast< mi::Uint32 >(heightfield_bounds.min.y);
    heightfield_area.max.x = static_cast< mi::Uint32 >(heightfield_bounds.max.x);
    heightfield_area.max.y = static_cast< mi::Uint32 >(heightfield_bounds.max.y);

    // Access the distribution scheme
    mi::base::Handle<const nv::index::IData_distribution> distribution_layout(
        dice_transaction->access<const nv::index::IData_distribution>(dist_layout_tag));

    // Heightfield distribution layout
    if(editing)
    {
        nv::index::IRegular_heightfield_data_locality* data_locality =
            distribution_layout->retrieve_data_locality_for_editing(
                heightfield_tag, heightfield_area, dice_transaction);
        return data_locality;
    }
    else
    {
        nv::index::IRegular_heightfield_data_locality* data_locality =
            distribution_layout->retrieve_data_locality(
                heightfield_tag, heightfield_area, dice_transaction);
        return data_locality;
    }
}
// --------------------------------------------------------------------------------------------
} // anonymous namespace

// --------------------------------------------------------------------------------------------
Autotracking_workflow_operation::Autotracking_workflow_operation(
    const mi::neuraylib::Tag&                       session_tag,
    const mi::neuraylib::Tag&                       heightfield_tag,
    const mi::neuraylib::Tag&                       volume_tag,
    bool                                            eight_way_seeding,
    mi::Uint32                                      max_nb_interations)
 :  m_session_tag(session_tag),
    m_heightfield_tag(heightfield_tag),
    m_volume_tag(volume_tag),
    m_eight_way_seeding(eight_way_seeding ? 1 : 0),
    m_max_nb_interations(max_nb_interations),
    m_nb_iterations(0),
    m_seeding_and_evaluating(NULL),
    m_heightfield_seeding_updates(NULL),
    m_elevation_updates(false)
{
}

Autotracking_workflow_operation::~Autotracking_workflow_operation()
{
    if(m_seeding_and_evaluating)
        delete m_seeding_and_evaluating;
        
    if(m_heightfield_seeding_updates)
        delete m_heightfield_seeding_updates;
}

bool Autotracking_workflow_operation::needs_iteration() const
{
    if(!m_seeding_and_evaluating || ! m_heightfield_seeding_updates)
        return true;

    const mi::Uint32 nb_seed_points = m_seeding_and_evaluating->get_nb_seed_points();
    return (nb_seed_points != 0 && m_nb_iterations<m_max_nb_interations);
}

/// FIXME: need test: value range to bbox
void Autotracking_workflow_operation::get_updated_bounding_box(
    mi::math::Bbox_struct<mi::Float32, 3>& bbox) const
{
    if(m_heightfield_seeding_updates && m_elevation_updates)
    {
        mi::Float32 min_height = -1.f;
        mi::Float32 max_height = -1.f;
        m_heightfield_seeding_updates->get_min_max_heights(min_height, max_height);

        if(bbox.min.z > min_height){
            bbox.min.z = min_height;
        }
        if(bbox.max.z < max_height){
            bbox.max.z = max_height;
        }
        bbox.min.z = (bbox.min.z < 0.0f) ? 0.0f : bbox.min.z;
            
        INFO_LOG << "Updated heightfield heights: " << bbox.min.z << ", " << bbox.max.z;
    }
}

// --------------------------------------------------------------------------------------------
// The fragment execute locally manages the workflow operations.
void Autotracking_workflow_operation::execute_fragment(
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Size                                        index,
    mi::Size                                        count,
    const mi::neuraylib::IJob_execution_context*    context)
{
    if(m_nb_iterations == 0)
    {
        // Determine the seed points first.
        std::vector<mi::math::Vector_struct<mi::Float32, 3> > seed_points;
        determine_seed_points(m_session_tag, m_heightfield_tag, dice_transaction, seed_points);

        if(seed_points.empty())
        {
            INFO_LOG << "No seed point given. Terminate autotracking workflow early.";
            return;
        }

        // Query the heightfield extent/bounds
        mi::base::Handle<const nv::index::IRegular_heightfield> heightfield(
            dice_transaction->access<const nv::index::IRegular_heightfield>(m_heightfield_tag));
        const mi::math::Bbox_struct<mi::Float32, 3> heightfield_bounds = heightfield->get_IJK_bounding_box();
        mi::math::Vector_struct<mi::Uint32, 2> heightfield_size;
        heightfield_size.x = static_cast< mi::Uint32 >(heightfield_bounds.max.x - heightfield_bounds.min.x);
        heightfield_size.y = static_cast< mi::Uint32 >(heightfield_bounds.max.y - heightfield_bounds.min.y);

        // Create data locality ....
        mi::base::Handle<nv::index::IRegular_heightfield_data_locality> data_locality_for_editing(
            compute_distribution_layout(m_session_tag, m_heightfield_tag, dice_transaction, true));
        mi::base::Handle<nv::index::IRegular_heightfield_data_locality> data_locality_for_seeding(
            compute_distribution_layout(m_session_tag, m_heightfield_tag, dice_transaction, false));

        // Cluster id referencing hosts that host heightfield data
        std::vector<mi::Uint32> cluster_host_ids_for_editing;
        const mi::Uint32 nb_cluster_hosts_for_editing = data_locality_for_editing->get_nb_cluster_nodes();
        for(mi::Uint32 i=0; i<nb_cluster_hosts_for_editing; ++i)
        {
            cluster_host_ids_for_editing.push_back(data_locality_for_editing->get_cluster_node(i));
            INFO_LOG << "Data locality for editing: " << i << " on host: " << data_locality_for_editing->get_cluster_node(i);
        }

        // Cluster id referencing hosts that host heightfield data
        std::vector<mi::Uint32> cluster_host_ids_for_seeding;
        const mi::Uint32 nb_cluster_hosts_for_seeding = data_locality_for_seeding->get_nb_cluster_nodes();
        for(mi::Uint32 i=0; i<nb_cluster_hosts_for_seeding; ++i)
        {
            cluster_host_ids_for_seeding.push_back(data_locality_for_seeding->get_cluster_node(i));
            INFO_LOG << "Data locality for seeding: " << i << " on host: " << data_locality_for_editing->get_cluster_node(i);
        }

        // 1.Step: Delete all heightfield elevation values within the heightfield's *and* the global's region of interests (set to '-1').
        Distributed_heightfield_elevation_delete delete_operation(
            m_session_tag, m_heightfield_tag, cluster_host_ids_for_editing);
        dice_transaction->execute_fragmented(&delete_operation, cluster_host_ids_for_editing.size());

        // 2.Step: Set all heightfield elevation values marked as hole within the polygon area.
        m_seeding_and_evaluating = new Distributed_seed_point_evaluation(
            m_session_tag,
            m_heightfield_tag,
            m_volume_tag,
            seed_points,
            cluster_host_ids_for_seeding,
            heightfield_size,
            m_eight_way_seeding != 0);
        print_seeds(seed_points);
        
        m_heightfield_seeding_updates = new Distributed_heightfield_seedings_updates(
            m_session_tag, m_heightfield_tag, cluster_host_ids_for_editing);
    }

    //
    const bool needs_execution = this->needs_iteration();
    if(needs_execution)
    {
        INFO_LOG << "SEEDING .......... (iteration = " << m_nb_iterations << ")";

        if(m_elevation_updates)
            m_seeding_and_evaluating->enable_elevation_updates();

        dice_transaction->execute_fragmented(m_seeding_and_evaluating, m_seeding_and_evaluating->get_nb_of_fragments());
        m_seeding_and_evaluating->update_seed_points();
        
        if(m_elevation_updates)
        {
            m_seeding_and_evaluating->disable_elevation_updates();
            const std::vector<std::vector<std::pair<mi::math::Bbox_struct<mi::Uint32, 2>, mi::neuraylib::Tag> > >& patch_tag_lists = 
                m_seeding_and_evaluating->get_patch_tag_lists();
            m_heightfield_seeding_updates->set_patch_tag_lists(patch_tag_lists);

            dice_transaction->execute_fragmented(m_heightfield_seeding_updates, m_heightfield_seeding_updates->get_nb_of_fragments());
        }
    }
    m_nb_iterations++;

    // Final call ...........
    if(!this->needs_iteration())
        finish(dice_transaction);
}

void Autotracking_workflow_operation::finish(
    mi::neuraylib::IDice_transaction*   dice_transaction)
{
    // Finally update the heightfield elevation values
    m_elevation_updates = true;
    m_seeding_and_evaluating->enable_elevation_updates();
    dice_transaction->execute_fragmented(m_seeding_and_evaluating, m_seeding_and_evaluating->get_nb_of_fragments());
    m_seeding_and_evaluating->disable_elevation_updates();
    
    dice_transaction->execute_fragmented(m_heightfield_seeding_updates, m_heightfield_seeding_updates->get_nb_of_fragments());
}


// --------------------------------------------------------------------------------------------

