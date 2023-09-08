/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief auto tracking workflow operation
#ifndef NVIDIA_INDEX_AUTOTRACKING_WORKFLOW_OPERATION_H
#define NVIDIA_INDEX_AUTOTRACKING_WORKFLOW_OPERATION_H

#include <mi/base/types.h>
#include <mi/dice.h>
#include <mi/base/interface_declare.h>

#include <nv/index/idistributed_compute_algorithm.h>
#include <nv/index/idistributed_data_locality.h>

#include "distributed_seed_point_evaluation.h"
#include "distributed_heightfield_seeding_updates.h"

#include <mi/base/lock.h>

#include <vector>

/// Proof of concept autotracking compute algorithm that ...
class Autotracking_workflow_operation :
    public nv::index::Distributed_compute_algorithm<0xce2ddd71,0xe6b9,0x44fd,0x83,0x39,0xa0,0x6b,0x8e,0x94,0xd7,0x53>
{
public:
    Autotracking_workflow_operation(
        const mi::neuraylib::Tag& session,
        const mi::neuraylib::Tag& heightfield_tag,
        const mi::neuraylib::Tag& volume_tag,
        bool                      eight_way_seeding,
        mi::Uint32                max_nb_interations);

    Autotracking_workflow_operation()
        :
        m_session_tag(mi::neuraylib::NULL_TAG),
        m_heightfield_tag(mi::neuraylib::NULL_TAG),
        m_volume_tag(mi::neuraylib::NULL_TAG),
        m_eight_way_seeding(0),
        m_max_nb_interations(0),
        m_nb_iterations(0),
        m_seeding_and_evaluating(0),
        m_heightfield_seeding_updates(0),
        m_elevation_updates(false)
    {
        // for serialization only
    }

    virtual ~Autotracking_workflow_operation();
    
    /// 
    bool needs_iteration() const;
    
    void enable_elevation_updates()
    {
        m_elevation_updates = true;
    }
    void disable_elevation_updates()
    {
        m_elevation_updates = false;
    }
    
    /// Number of fragments to start the compute algorithm
    virtual mi::Size get_nb_of_fragments() const { return 1; /* only one fragment, which starts sub jobs! */ }
    
    /// updated bounding box after compute
    virtual void get_updated_bounding_box(
        mi::math::Bbox_struct<mi::Float32, 3>& bbox) const;


    // -----------------------------------------------------------------------------------------
    /// Parallel compute on the local cluster host.
    /// \param dice_transaction         The transaction the job is executing in.
    /// \param index                    The index identifies for which fragment the execute function is called.
    /// \param count                    The number of fragments in which the job is split.
    virtual void execute_fragment(
        mi::neuraylib::IDice_transaction*               dice_transaction,
        mi::Size                                        index,
        mi::Size                                        count,
        const mi::neuraylib::IJob_execution_context*    context);
    
private:
    void finish(mi::neuraylib::IDice_transaction*   dice_transaction);

private:
    mi::neuraylib::Tag                      m_session_tag;
    mi::neuraylib::Tag                      m_heightfield_tag;
    mi::neuraylib::Tag                      m_volume_tag;
    mi::Uint32                              m_eight_way_seeding;
    mi::Uint32                              m_max_nb_interations;
    
    mi::Uint32                              m_nb_iterations;
    Distributed_seed_point_evaluation*      m_seeding_and_evaluating;
    Distributed_heightfield_seedings_updates* m_heightfield_seeding_updates;
    bool                                    m_elevation_updates;
};

// ----------------------------------------------------------------------------------------------


#endif // NVIDIA_INDEX_AUTOTRACKING_WORKFLOW_OPERATION_H
