/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/

#include "irregular_volume_data_editing.h"

#include <cassert>
#include <map>
#include <cstdlib>
#include <iostream>
#include <nv/index/iirregular_volume_compute_task.h>
#include <nv/index/iirregular_volume_scene_element.h>
#include <nv/index/iirregular_volume_subset.h>
#include <nv/index/isession.h>

#include "common/forwarding_logger.h"
#include "common/irregular_volume_importer.h"

#include "../alya/alyaglobal.h"

//----------------------------------------------------------------------
/// Edit irregular volume compute task
class Scala_editing_operation : public nv::index::Irregular_volume_compute_task
{
public:
    /// \param ijk_bbox         Specifies the part of the volume where the operation should be applied.
    /// \param volume_size      Size of the entire volume.
    /// \param border           volume border size.
    /// \param bricks_locality  A vector with all bricks locality information.
    Scala_editing_operation(
        mi::Uint32                              time_step,
        const mi::math::Bbox<mi::Float32, 3>&   ijk_bbox)
      : m_time_step(time_step), m_active_ijk_bbox(ijk_bbox)
    {
        // Not used
    }

  virtual bool edit(
        const mi::math::Bbox_struct<mi::Float32, 3>&    bbox,
        nv::index::IIrregular_volume_subset*            subset_data,
        mi::neuraylib::IDice_transaction*               dice_transaction) const
    {
      // Not used
      //INFO_LOG << "Editing call for subset contained in " << bbox;
      
        const mi::Uint32 nb_attributes = subset_data->get_nb_attributes();
        if (nb_attributes > 0)
        {
            nv::index::IIrregular_volume_subset::Attribute_storage      attribute_store;
            nv::index::IIrregular_volume_subset::Attribute_parameters   attribute_parameter;

            subset_data->get_attribute_parameters(0, attribute_parameter);
            subset_data->get_attribute(0, attribute_store);
            const mi::Uint32 nb_attrib_values = attribute_parameter.nb_attrib_values;
            if (attribute_parameter.type == nv::index::IIrregular_volume_subset::ATTRIB_TYPE_UINT8)
            {
                DEBUG_LOG << "type: Uint8, " << nb_attrib_values << " elements - subset/bbox: " << bbox;
                mi::Uint8* values = static_cast<mi::Uint8*>(attribute_store.attrib_values);
                for (mi::Uint32 i = 0; i < nb_attrib_values; ++i)
                {
                    mi::Uint8 value = values[i];
                    if (value == 255) value = 0;
                    else value += 1;
                    values[i] = value;
                }
                return true;
            }
            if (attribute_parameter.type == nv::index::IIrregular_volume_subset::ATTRIB_TYPE_UINT16)
            {
                DEBUG_LOG << "type: Uint16, " << nb_attrib_values << " elements - subset/bbox: " << bbox;
                mi::Uint16* values = static_cast<mi::Uint16*>(attribute_store.attrib_values);
                for (mi::Uint32 i = 0; i < nb_attrib_values; ++i)
                {
                    mi::Uint16 value = values[i];
                    if (value > mi::base::numeric_traits<mi::Uint16>::max()-100) value = 0;
                    else value += 100;
                    values[i] = value;
                }
                return true;
            }
            if (attribute_parameter.type == nv::index::IIrregular_volume_subset::ATTRIB_TYPE_FLOAT32)
            {
                DEBUG_LOG << "type: Float32, " << nb_attrib_values << " elements - subset/bbox: " << bbox;
                mi::Float32* values = static_cast<mi::Float32*>(attribute_store.attrib_values);
		for (mi::Uint32 i = 0; i < nb_attrib_values; ++i)
		  {
		    mi::Float32 value = values[i];
		    if (value > mi::base::numeric_traits<mi::Float32>::max()+1.f) value = 0;
		    else value += 1.f;
		    if (value > 99.f) value = -50.f;
		    else value += 1.f;
		  
		    //mi::Uint32 idx = global_to_local_vtx_idx_map.at( global_numbering[i] - 1 );
		    //values[idx] = pt_to_scalardata[i];
		  }
                return true;
            }
        }
        return false;
    }

    virtual nv::index::IIrregular_volume_subset* edit(
        const mi::math::Bbox_struct<mi::Float32, 3>& bbox,
        nv::index::IData_subset_factory*             factory,
        mi::neuraylib::IDice_transaction*            dice_transaction) const
    {
        std::ostringstream folder_and_basic_name;
        // \\NETAPP-DEB03\bigdata\medical\bsc_heart\300kheart150_ts\heart-intra-000142.ts
        folder_and_basic_name << "\\\\NETAPP-DEB03\\bigdata\\medical\\bsc_heart\\300kheart150_ts\\heart-intra-";
        if(m_time_step<10)
            folder_and_basic_name << "00000" << m_time_step << ".ts";
        else if (m_time_step<100)
            folder_and_basic_name << "0000" << m_time_step << ".ts";
        else if (m_time_step<1000)
            folder_and_basic_name << "000" << m_time_step << ".ts";

        const std::string filename = folder_and_basic_name.str();
        
        nv::index_common::Irregular_volume_importer importer(filename, false);

        nv::index::IIrregular_volume_subset* subset = static_cast<nv::index::IIrregular_volume_subset*>(
            importer.create(bbox, factory, dice_transaction));

        return subset;
    }

    virtual Operation_mode get_operation_mode() const
    {
        return OPERATION_MODE_SCALA_VALUE_EDITING;
        //return OPERATION_MODE_TOPOLOGY_EDITING;
    }

private:
    mi::Uint32                          m_time_step;
    mi::math::Bbox<mi::Float32, 3>      m_active_ijk_bbox;
};

//----------------------------------------------------------------------
/// Fragmented job which starts one fragment per GPU and executes the compute tasks for all bricks
/// that are stored on that GPU.
class Local_compute_scheduler :
    public mi::neuraylib::Fragmented_job<0x65f080d,0xe51d,0x437c,0x8f,0x6a,0x37,0xeb,0xb7,0x61,0xfa,0x59>
{
public:
    Local_compute_scheduler(
        mi::Uint32                                      host_id,
        nv::index::IIrregular_volume_data_locality*     data_locality,
        nv::index::Irregular_volume_compute_task*       compute_task)
        : m_host_id(host_id),
          m_data_locality(data_locality),
          m_compute_task(compute_task)
    {
    }

    /// Number of fragments to start the job
    mi::Size get_nb_of_fragments() const
    {
        const mi::Uint32 nb_subsets = m_data_locality->get_nb_bounding_box(m_host_id);
        return nb_subsets;
    }

    /// Process the bricks which are all located on a single device
    virtual void execute_fragment(
        mi::neuraylib::IDice_transaction*               dice_transaction,
        mi::Size                                        index,
        mi::Size                                        count,
        const mi::neuraylib::IJob_execution_context*    context)
    {
        mi::base::Handle<nv::index::IIrregular_volume_data_edit> data_editing(
            m_data_locality->create_data_edit(dice_transaction, m_host_id, index));
        data_editing->execute_compute_task(m_compute_task, dice_transaction);
    }

    /// Returns the CPU load per fragment of the fragmented job.
    virtual mi::Float32 get_cpu_load() const { return 1.f; }

    /// Returns the GPU load per fragment of the fragmented job.
    virtual mi::Float32 get_gpu_load() const { return 0.f; }

private:
    mi::Uint32                                              m_host_id;
    nv::index::IIrregular_volume_data_locality*             m_data_locality;
    nv::index::Irregular_volume_compute_task*               m_compute_task;
};

//----------------------------------------------------------------------
Irregular_volume_data_editing::Irregular_volume_data_editing(
    mi::Uint32                                                      time_step,
    const mi::neuraylib::Tag&                                       session,
    const mi::neuraylib::Tag&                                       volume_scene_element,
    const mi::math::Bbox<mi::Float32, 3>&                           active_ijk_bbox,
    nv::index::IIrregular_volume_data_locality*                     data_locality)
  : m_time_step(time_step),
    m_session_tag(session),
    m_scene_element_tag(volume_scene_element),
    m_active_ijk_bbox(active_ijk_bbox),
    m_data_locality(data_locality),
    m_nb_bricks(0)
{
    const mi::Size nb_hosts = m_data_locality->get_nb_cluster_nodes();
    for(mi::Size j=0; j<nb_hosts; ++j)
    {
        const mi::Uint32 host = m_data_locality->get_cluster_node(j);
        const mi::Uint32 nb_bounding_boxes = m_data_locality->get_nb_bounding_box(host);
        m_nb_bricks += nb_bounding_boxes;
    }
}

//----------------------------------------------------------------------
Irregular_volume_data_editing::Irregular_volume_data_editing()
  : m_time_step(0),
    m_session_tag(),
    m_scene_element_tag(),
    m_active_ijk_bbox(),
    m_data_locality(0),
    m_nb_bricks(0)
{
    // for serialization only
}

// -----------------------------------------------------------------------------------------
mi::Size Irregular_volume_data_editing::get_nb_of_fragments() const
{
    // One fragment per host
    return m_data_locality->get_nb_cluster_nodes();
}

// -----------------------------------------------------------------------------------------
void Irregular_volume_data_editing::create_remote_data_locality(
    mi::neuraylib::IDice_transaction*   dice_transaction)
{
    mi::base::Lock::Block block(&m_lock);

    if (m_data_locality == 0)
    {
        // Access the session as the main resource hook in NVIDIA IndeX:
        mi::base::Handle<const nv::index::ISession> session(
            dice_transaction->access<nv::index::ISession>(m_session_tag));

        // Access the sessions's distribution scheme
        const mi::neuraylib::Tag data_distribution_tag = session->get_distribution_layout();
        mi::base::Handle<const nv::index::IData_distribution> distribution(
            dice_transaction->access<nv::index::IData_distribution>(data_distribution_tag));

        // Retrieve the data locality information for the volume dataset.
        // If set the bounding box restricts the data locality information to a sub volume.
        m_data_locality = distribution->retrieve_data_locality(
            m_scene_element_tag, m_active_ijk_bbox, dice_transaction);
    }
}

// -----------------------------------------------------------------------------------------
void Irregular_volume_data_editing::assign_fragments_to_hosts(
    mi::Uint32* slots,
    mi::Size    nr_slots)
{
    // Assign one fragment to each host
    const mi::Size nb_hosts = m_data_locality->get_nb_cluster_nodes();
    if (nr_slots != nb_hosts)
    {
        ERROR_LOG << "The number of hosts (" << nb_hosts << ") and the number of the requested slots "
                  << "(" << nr_slots << ") don't match when assigning "
                  << "fragments (job tasks) to the cluster machines.";
    }

    m_host_ids.resize(nb_hosts);
    
    for (mi::Size i=0; i < nb_hosts; ++i)
    {
        const mi::Uint32 host_id = m_data_locality->get_cluster_node(i);
        // Schedule the job execution of fragment i on this host
        slots[i] = host_id;
        // Also store the host id locally
        m_host_ids[i] = host_id;
    }

}

// --------------------------------------------------------------------------------------------
void Irregular_volume_data_editing::invoke_compute_tasks_per_host(
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Uint32                                      host_id)
{
    // Edit the irregular volume subset data on the given host by running a fragmented job
    // that makes sure the compute task for each local subset is executed.
    Scala_editing_operation op(m_time_step, m_active_ijk_bbox);
    Local_compute_scheduler scheduler(host_id, m_data_locality, &op);
    dice_transaction->execute_fragmented(&scheduler, scheduler.get_nb_of_fragments());
}

//----------------------------------------------------------------------
void Irregular_volume_data_editing::invoke_compute_tasks_per_brick(
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Uint32                                      host_id,
    mi::Uint32                                      subset_id)
{
    // Edit the brick by applying a compute task
    mi::base::Handle<nv::index::IIrregular_volume_data_edit> brick_data_editing(
        m_data_locality->create_data_edit(dice_transaction, host_id, subset_id));

    Scala_editing_operation op(m_time_step, m_active_ijk_bbox);
    brick_data_editing->execute_compute_task(&op, dice_transaction);
}

// --------------------------------------------------------------------------------------------
void Irregular_volume_data_editing::execute_fragment(
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Size                                        index,
    mi::Size                                        count,
    const mi::neuraylib::IJob_execution_context*    context)
{
    const mi::Uint32 host_id = m_host_ids[index];

    // INFO_LOG << "GPU job with index " << index << " (of " << count << ") "
             // << "runs on host " << host_id;
             
    invoke_compute_tasks_per_host(dice_transaction, host_id);
}

//----------------------------------------------------------------------
void Irregular_volume_data_editing::execute_fragment_remote(
    mi::neuraylib::ISerializer*                     serializer,
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Size                                        index,
    mi::Size                                        count,
    const mi::neuraylib::IJob_execution_context*    context)
{
    const mi::Uint32 host_id = m_host_ids[index];

    // Create the data locality on the remote machine
    create_remote_data_locality(dice_transaction);

    // INFO_LOG << "GPU job with index " << index << " (of " << count << ") "
             // << "runs on remote host " << host_id;
             
    invoke_compute_tasks_per_host(dice_transaction, host_id);
}

// --------------------------------------------------------------------------------------------
void Irregular_volume_data_editing::serialize(mi::neuraylib::ISerializer* serializer) const
{
    serializer->write(&m_time_step, 1);
    serializer->write(&m_session_tag.id, 1);
    serializer->write(&m_scene_element_tag.id, 1);
    serializer->write(&m_active_ijk_bbox.min.x, 6);

    serializer->write(&m_nb_bricks, 1);

    mi::Uint32 nb_elements = static_cast<mi::Uint32>(m_host_ids.size());
    serializer->write(&nb_elements, 1);
    for (mi::Uint32 i = 0; i < nb_elements; ++i)
    {
        serializer->write(&m_host_ids[i], 1);
    }

    nb_elements = static_cast<mi::Uint32>(m_brick_ids.size());
    serializer->write(&nb_elements, 1);
    for (mi::Uint32 i = 0; i < nb_elements; ++i)
    {
        serializer->write(&m_brick_ids[i], 1);
    }
}

//----------------------------------------------------------------------
void Irregular_volume_data_editing::deserialize(mi::neuraylib::IDeserializer* deserializer)
{
    deserializer->read(&m_time_step, 1);
    deserializer->read(&m_session_tag.id, 1);
    deserializer->read(&m_scene_element_tag.id, 1);
    deserializer->read(&m_active_ijk_bbox.min.x, 6);

    deserializer->read(&m_nb_bricks, 1);

    mi::Uint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_host_ids.resize(nb_elements);
    for (mi::Uint32 i = 0; i < nb_elements; ++i)
    {
        deserializer->read(&m_host_ids[i], 1);
    }

    nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_brick_ids.resize(nb_elements);
    for (mi::Uint32 i = 0; i < nb_elements; ++i)
    {
        deserializer->read(&m_brick_ids[i], 1);
    }
}
