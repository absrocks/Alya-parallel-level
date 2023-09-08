/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "distributed_volume_bordermap_generation.h"

#include <mi/dice.h>
#include <nv/index/isession.h>
#include <nv/index/iconfig_settings.h>

#include <cassert>

#include "scene_utility.h"
#include "utilities.h"
#include "volume_bordermap_database_element.h"
#include "volume_filter_compute_task.h"

#include "common/forwarding_logger.h"
#include "common/type_conversion_utility.h"

//----------------------------------------------------------------------
Distributed_volume_bordermap_generation::Distributed_volume_bordermap_generation(
    const mi::neuraylib::Tag& session_tag,
    const mi::neuraylib::Tag& volume_tag,
    mi::Sint32 filter_type,
    mi::math::Bbox< mi::Uint32, 3 > const & active_ijk_bbox,
    std::vector< mi::Uint32 > const &       cluster_hosts)
    :
    m_session_tag(session_tag),
    m_volume_tag(volume_tag),
    m_bordermap_element_tag(mi::neuraylib::NULL_TAG),
    m_filter_type(filter_type),
    m_active_ijk_bbox(active_ijk_bbox),
    m_cluster_hosts(cluster_hosts)
{
    m_result_job_bordermap_tag_vec.resize(0);

    INFO_LOG << "Distributed_volume_bordermap_generation: creating border set "
             << active_ijk_bbox <<  ".";
}

//----------------------------------------------------------------------
Distributed_volume_bordermap_generation::~Distributed_volume_bordermap_generation()
{
    DEBUG_LOG << "~Distributed_volume_bordermap_generation";
}

//----------------------------------------------------------------------
std::vector< mi::neuraylib::Tag > const &
Distributed_volume_bordermap_generation::get_job_result_bordermap_tag_vec() const
{
    return m_result_job_bordermap_tag_vec;
}

//----------------------------------------------------------------------
void Distributed_volume_bordermap_generation::assign_fragments_to_hosts(
    mi::Uint32* slots,
    mi::Size    nr_slots)
{
    if (nr_slots != m_cluster_hosts.size()){
        ERROR_LOG << "Distributed_volume_bordermap_generation: The number of host ids present ("
                  << m_cluster_hosts.size()
                  << ") and requested (" << nr_slots << ") doesn't match the requested count.";
    }
    else{
        for (mi::Size i=0; i < nr_slots; ++i){
            slots[i] = m_cluster_hosts[i];
        }
    }
}

//----------------------------------------------------------------------
void Distributed_volume_bordermap_generation::serialize(
    mi::neuraylib::ISerializer* serializer) const
{
    serializer->write(&m_session_tag.id, 1);
    serializer->write(&m_volume_tag.id, 1);
    serializer->write(&m_bordermap_element_tag.id, 1);
    serializer->write(&m_filter_type, 1);
    serializer->write(&m_active_ijk_bbox.min.x, 1);
    serializer->write(&m_active_ijk_bbox.min.y, 1);
    serializer->write(&m_active_ijk_bbox.min.z, 1);
    serializer->write(&m_active_ijk_bbox.max.x, 1);
    serializer->write(&m_active_ijk_bbox.max.y, 1);
    serializer->write(&m_active_ijk_bbox.max.z, 1);
    mi::Uint32 nb_elements = mi::Uint32(m_cluster_hosts.size());
    serializer->write(&nb_elements, 1);
    for (mi::Uint32 i=0; i < nb_elements; ++i){
        serializer->write(&m_cluster_hosts[i], 1);
    }
}

//----------------------------------------------------------------------
void Distributed_volume_bordermap_generation::deserialize(
    mi::neuraylib::IDeserializer* deserializer)
{
    deserializer->read(&m_session_tag.id, 1);
    deserializer->read(&m_volume_tag.id, 1);
    deserializer->read(&m_bordermap_element_tag.id, 1);
    deserializer->read(&m_filter_type, 1);
    deserializer->read(&m_active_ijk_bbox.min.x, 1);
    deserializer->read(&m_active_ijk_bbox.min.y, 1);
    deserializer->read(&m_active_ijk_bbox.min.z, 1);
    deserializer->read(&m_active_ijk_bbox.max.x, 1);
    deserializer->read(&m_active_ijk_bbox.max.y, 1);
    deserializer->read(&m_active_ijk_bbox.max.z, 1);
    mi::Uint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_cluster_hosts.resize(nb_elements);
    for (mi::Uint32 i=0; i < nb_elements; ++i){
        deserializer->read(&m_cluster_hosts[i], 1);
    }
}

//----------------------------------------------------------------------
void Distributed_volume_bordermap_generation::create_border(
    mi::neuraylib::IDice_transaction*         dice_transaction,
    nv::index::IRegular_volume_data_locality* locality,
    mi::Uint32                                host_id,
    mi::neuraylib::Tag                        job_local_bordermap_tag)
{
    mi::Uint32 const nb_brick_edits = static_cast<mi::Uint32>(locality->get_nb_bounding_box(host_id));
    assert(job_local_bordermap_tag.is_valid());

    for(mi::Uint32 i=0; i < nb_brick_edits; ++i)
    {
        mi::base::Handle<nv::index::IRegular_volume_data_edit> volume_data_edit(
            locality->create_data_edit(dice_transaction, host_id, i));
        assert(volume_data_edit.is_valid_interface());

        mi::base::Handle<Volume_bordermap_generation> gen_opration(
            new Volume_bordermap_generation(
                nv::index_common::convert_bbox_type<mi::Sint64, mi::Uint32, 3>(m_active_ijk_bbox),
                m_session_tag,
                m_volume_tag,
                job_local_bordermap_tag
                ));
        {
            volume_data_edit->execute_compute_task(gen_opration.get(), dice_transaction);
        }
    }
}

//----------------------------------------------------------------------
void Distributed_volume_bordermap_generation::execute_fragment(
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Size                                        index,
    mi::Size                                        count,
    const mi::neuraylib::IJob_execution_context*    context)
{
    mi::Uint32 const cluster_host = m_cluster_hosts[index];
    assert(count > 0);

    // Note: actually, we need a lock this. But most probably we would
    // not see any problem.
    m_result_job_bordermap_tag_vec.resize(count, mi::neuraylib::NULL_TAG);

    mi::base::Handle< nv::index::IRegular_volume_data_locality > locality(
        new_volume_distribution_layout(m_session_tag,
                                       m_volume_tag,
                                       m_active_ijk_bbox,
                                       dice_transaction));

    // create the local border set and store to the DB
    mi::neuraylib::Tag job_bmap_tag = mi::neuraylib::NULL_TAG;
    {
        mi::base::Handle< Volume_bordermap_database_element > job_bmap(
            new Volume_bordermap_database_element());
        job_bmap_tag = dice_transaction->store(job_bmap.get());
        assert(job_bmap_tag.is_valid());
    }

    INFO_LOG << "Distributed_volume_bordermap_generation: create border on local host: "
             << cluster_host;
    this->create_border(dice_transaction, locality.get(), cluster_host,
                        job_bmap_tag);

    assert(index < m_result_job_bordermap_tag_vec.size());
    m_result_job_bordermap_tag_vec[index] = job_bmap_tag;
}

//----------------------------------------------------------------------
void Distributed_volume_bordermap_generation::execute_fragment_remote(
    mi::neuraylib::ISerializer*                     serializer,
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Size                                        index,
    mi::Size                                        count,
    const mi::neuraylib::IJob_execution_context*    context)
{
    mi::Uint32 const cluster_host = m_cluster_hosts[index];
    assert(m_result_job_bordermap_tag_vec.empty());

    mi::base::Handle< nv::index::IRegular_volume_data_locality > locality(
        new_volume_distribution_layout(m_session_tag,
                                       m_volume_tag,
                                       m_active_ijk_bbox,
                                       dice_transaction));

    // create the local border set and store to the DB
    mi::neuraylib::Tag job_bmap_tag = mi::neuraylib::NULL_TAG;
    {
        mi::base::Handle< Volume_bordermap_database_element > job_bmap(
            new Volume_bordermap_database_element());
        job_bmap_tag = dice_transaction->store(job_bmap.get());
        assert(job_bmap_tag.is_valid());
    }

    DEBUG_LOG << "Distributed_volume_bordermap_generation: create border on remote host: "
              << cluster_host;
    this->create_border(dice_transaction, locality.get(), cluster_host,
                        job_bmap_tag);
    // send to the host
    serializer->write(&(job_bmap_tag.id), 1);
}

//----------------------------------------------------------------------
void Distributed_volume_bordermap_generation::receive_remote_result(
    mi::neuraylib::IDeserializer*       deserializer,
    mi::neuraylib::IDice_transaction*   transaction,
    mi::Size                            index,
    mi::Size                            count)
{
    // If this failed, execute_fragment really needs a lock there.
    assert(m_result_job_bordermap_tag_vec.size() == count);

    DEBUG_LOG << "Distributed_volume_bordermap_generation::Receiving result of fragment "
              << index << " of " << count << " from a remote." << std::endl;
    mi::neuraylib::Tag result_bordermap_tag = mi::neuraylib::NULL_TAG;

    deserializer->read(&(result_bordermap_tag.id), 1);

    assert(index < m_result_job_bordermap_tag_vec.size());
    m_result_job_bordermap_tag_vec[index] = result_bordermap_tag;
}

//----------------------------------------------------------------------
