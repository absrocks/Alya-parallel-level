/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "distributed_volume_filter.h"

#include <mi/dice.h>
#include <nv/index/isession.h>
#include <nv/index/iconfig_settings.h>

#include <cassert>
#include <iterator>

#include "compute_utility.h"
#include "scene_utility.h"
#include "utilities.h"
#include "volume_filter_compute_task.h"

#include "common/forwarding_logger.h"
#include "common/type_conversion_utility.h"

//----------------------------------------------------------------------
Distributed_volume_filter::Distributed_volume_filter(
    const mi::neuraylib::Tag& session_tag,
    const mi::neuraylib::Tag& volume_tag,
    const mi::neuraylib::Tag& bordermap_element_tag,
    mi::Sint32 filter_type,
    mi::math::Bbox< mi::Uint32, 3 > const & active_ijk_bbox,
    std::vector< mi::Uint32 > const & cluster_hosts,
    mi::Uint32 parallel_count)
    :
    m_session_tag(session_tag),
    m_volume_tag(volume_tag),
    m_bordermap_element_tag(bordermap_element_tag),
    m_filter_type(filter_type),
    m_active_ijk_bbox(active_ijk_bbox),
    m_cluster_hosts(cluster_hosts),
    m_parallel_count(parallel_count),
    m_fragment_to_host_id(),
    m_fragment_to_host_local_thread_id()
{
    assert(m_bordermap_element_tag.is_valid());
    assert(m_parallel_count > 0);
    INFO_LOG << "Distributed_volume_filter: creating border set "
             << active_ijk_bbox <<  ".";

    setup_parallel_fragment(m_cluster_hosts,
                            m_parallel_count,
                            m_fragment_to_host_id,
                            m_fragment_to_host_local_thread_id);
}

//----------------------------------------------------------------------
Distributed_volume_filter::~Distributed_volume_filter()
{
    DEBUG_LOG << "~Distributed_volume_filter";
}

//----------------------------------------------------------------------
mi::Size Distributed_volume_filter::get_nb_of_fragments() const
{
    assert(m_fragment_to_host_id.size() == m_fragment_to_host_local_thread_id.size());

    mi::Size const frag_counts = m_fragment_to_host_id.size();
    INFO_LOG << "number of fragments: " << frag_counts;
    return frag_counts;
}


//----------------------------------------------------------------------
mi::neuraylib::IFragmented_job::Scheduling_mode
Distributed_volume_filter::get_scheduling_mode() const
{
    return mi::neuraylib::IFragmented_job::USER_DEFINED;
}

//----------------------------------------------------------------------
void Distributed_volume_filter::assign_fragments_to_hosts(
    mi::Uint32* slots,
    mi::Size    nr_slots)
{
    INFO_LOG << "assign_fragments_to_hosts: nr_slots = " << nr_slots;

    if (nr_slots != m_fragment_to_host_id.size()){
        ERROR_LOG << "Distributed_volume_filter: The number of host ids present ("
                  << m_cluster_hosts.size()
                  << ") and requested (" << nr_slots << ") doesn't match the requested count.";
    }
    else{
        for(mi::Size i=0; i < nr_slots; ++i){
            slots[i] = m_fragment_to_host_id[i];
        }
    }
}

//----------------------------------------------------------------------
void Distributed_volume_filter::serialize(
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
    {
        mi::Uint32 nb_elements = mi::Uint32(m_cluster_hosts.size());
        serializer->write(&nb_elements, 1);
        for (mi::Uint32 i=0; i < nb_elements; ++i){
            serializer->write(&m_cluster_hosts[i], 1);
        }
    }
    serializer->write(&m_parallel_count, 1);
    {
        mi::Uint32 nb_elements = mi::Uint32(m_fragment_to_host_id.size());
        serializer->write(&nb_elements, 1);
        for (mi::Uint32 i=0; i < nb_elements; ++i){
            serializer->write(&m_fragment_to_host_id[i], 1);
        }
    }
    {
        mi::Uint32 nb_elements = mi::Uint32(m_fragment_to_host_local_thread_id.size());
        serializer->write(&nb_elements, 1);
        for (mi::Uint32 i=0; i < nb_elements; ++i){
            serializer->write(&m_fragment_to_host_local_thread_id[i], 1);
        }
    }
}

//----------------------------------------------------------------------
void Distributed_volume_filter::deserialize(
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
    {
        mi::Uint32 nb_elements = 0;
        deserializer->read(&nb_elements, 1);
        m_cluster_hosts.resize(nb_elements);
        for (mi::Uint32 i=0; i < nb_elements; ++i){
            deserializer->read(&m_cluster_hosts[i], 1);
        }
    }
    deserializer->read(&m_parallel_count, 1);
    {
        mi::Uint32 nb_elements = 0;
        deserializer->read(&nb_elements, 1);
        m_fragment_to_host_id.resize(nb_elements);
        for (mi::Uint32 i=0; i < nb_elements; ++i){
            deserializer->read(&m_fragment_to_host_id[i], 1);
        }
    }
    {
        mi::Uint32 nb_elements = 0;
        deserializer->read(&nb_elements, 1);
        m_fragment_to_host_local_thread_id.resize(nb_elements);
        for (mi::Uint32 i=0; i < nb_elements; ++i){
            deserializer->read(&m_fragment_to_host_local_thread_id[i], 1);
        }
    }
}

//----------------------------------------------------------------------
void Distributed_volume_filter::apply_filter(
    mi::neuraylib::IDice_transaction*         dice_transaction,
    nv::index::IRegular_volume_data_locality* locality,
    mi::Uint32                                host_id,
    mi::Uint32                                fragment_index)
{
    assert(m_parallel_count > 0);
    mi::Uint32 const nb_brick_edits = static_cast<mi::Uint32>(locality->get_nb_bounding_box(host_id));

    assert(fragment_index < m_fragment_to_host_local_thread_id.size());
    mi::Sint32 const tid = m_fragment_to_host_local_thread_id[fragment_index];

    for(mi::Sint32 i = tid; i < static_cast< mi::Sint32 >(nb_brick_edits); i += m_parallel_count)
    {
        DEBUG_LOG << "apply_filter: fragment_index = " << fragment_index
                  << ", tid = " << tid << ", i = " << i;
        mi::base::Handle<nv::index::IRegular_volume_data_edit> volume_data_edit(
            locality->create_data_edit(dice_transaction, host_id, i));
        assert(volume_data_edit.is_valid_interface());

        mi::base::Handle<Volume_filter_apply> filter(
            new Volume_filter_apply(
                nv::index_common::convert_bbox_type<mi::Sint64, mi::Uint32, 3>(m_active_ijk_bbox),
                m_session_tag,
                m_volume_tag,
                m_bordermap_element_tag,
                m_filter_type,
                fragment_index));
        {
            volume_data_edit->execute_compute_task(filter.get(), dice_transaction);
        }
    }
}


//----------------------------------------------------------------------
void Distributed_volume_filter::execute_fragment(
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Size                                        index,
    mi::Size                                        count,
    const mi::neuraylib::IJob_execution_context*    context)
{
    mi::Uint32 const host_idx = static_cast<mi::Uint32>(index / m_parallel_count);
    assert(host_idx < m_cluster_hosts.size());
    mi::Uint32 const cluster_host = m_cluster_hosts[host_idx];

    mi::base::Handle< nv::index::IRegular_volume_data_locality > locality(
        new_volume_distribution_layout(m_session_tag,
                                       m_volume_tag,
                                       m_active_ijk_bbox,
                                       dice_transaction));

    INFO_LOG << "Distributed_volume_filter: apply filter on local host: "
             << cluster_host << ", " << index << "/" << count;
    this->apply_filter(dice_transaction, locality.get(), cluster_host, static_cast<mi::Uint32>(index));
}

//----------------------------------------------------------------------
void Distributed_volume_filter::execute_fragment_remote(
    mi::neuraylib::ISerializer*                     serializer,
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Size                                        index,
    mi::Size                                        count,
    const mi::neuraylib::IJob_execution_context*    context)
{
    mi::Uint32 const host_idx = static_cast<mi::Uint32>(index / m_parallel_count);
    assert(host_idx < m_cluster_hosts.size());
    mi::Uint32 const cluster_host = m_cluster_hosts[host_idx];

    mi::base::Handle< nv::index::IRegular_volume_data_locality > locality(
        new_volume_distribution_layout(m_session_tag,
                                       m_volume_tag,
                                       m_active_ijk_bbox,
                                       dice_transaction));

    INFO_LOG << "Distributed_volume_filter: apply filter on remote host: "
             << cluster_host << ", " << index << "/" << count;
    this->apply_filter(dice_transaction, locality.get(), cluster_host, static_cast<mi::Uint32>(index));
}

//----------------------------------------------------------------------
