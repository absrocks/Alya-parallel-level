/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/

#include "cuda_ipc_volume_editing.h"

#include <cassert>
#include <map>

#include <nv/index/iconfig_settings.h>
#include <nv/index/iregular_volume.h>
#include <nv/index/iregular_volume_compute_task.h>
#include <nv/index/isession.h>

#include "common/forwarding_logger.h"

#include <cuda_runtime.h>
#include "cuda_ipc_volume_editing_kernel.cuh"

//----------------------------------------------------------------------
/// Edit volume compute task
class Voxel_value_set_operation :
    public nv::index::Regular_volume_compute_task
{
public:
    /// \param ijk_bbox         Specifies the part of the volume where the operation should be applied.
    /// \param volume_size      Size of the entire volume.
    /// \param border           volume border size.
    /// \param bricks_locality  A vector with all bricks locality information.
    Voxel_value_set_operation(
        const mi::math::Bbox<mi::Uint32, 3>&   ijk_bbox,
        const mi::math::Vector<mi::Uint32, 3>& volume_size,
        mi::Uint32                             border,
        const std::vector<Brick_locality>&     bricks_locality)
      : m_active_ijk_bbox(ijk_bbox),
        m_volume_size(volume_size),
        m_border(border),
        m_bricks_locality(bricks_locality)
    {}

    /// Perform computation on the given volume data.
    /// \param brick_ijk_bbox               Local space 3D bounding box of the current volume brick.
    /// \param voxel_values                 The voxel data of the given volume brick.
    /// \param dice_transaction             Current transaction.
    virtual bool compute(
        const mi::math::Bbox_struct<mi::Uint32, 3>& brick_ijk_bbox_struct,
        nv::index::IRegular_volume_data*            voxel_data,
        mi::neuraylib::IDice_transaction*           dice_transaction) const
    {
        // Not used
        return true;
    }

    /// Perform computation on the given volume data using GPUs.
    /// \param brick_ijk_bbox               Local space 3D bounding box of the current volume brick.
    /// \param voxel_values                 The voxel data of the given volume brick.
    /// \param dice_transaction             Current transaction.
    virtual bool compute_gpu(
        const mi::math::Bbox_struct<mi::Uint32, 3>&  brick_ijk_bbox_struct,
        cudaArray*                                   voxel_values,
        mi::neuraylib::IDice_transaction*            dice_transaction) const
    {
        // lookup brick
        const mi::math::Bbox<mi::Sint32, 3> brick_ijk_bbox = mi::math::Bbox<mi::Sint32, 3>(
            mi::math::Bbox<mi::Uint32, 3>(brick_ijk_bbox_struct)); // convert to signed
        
        mi::Sint32 idx = -1;
        for(mi::Uint32 i=0; i<m_bricks_locality.size(); ++i)
        {
            if(m_bricks_locality[i].bbox.min == brick_ijk_bbox.min)
            {
                idx = i;
                break;
            }
        }
            
        const mi::math::Bbox<mi::Sint32, 3> work_bbox = get_work_bbox(brick_ijk_bbox);
        if (!work_bbox.is_volume())
        {
            return false;
        }
        
        // DEVICE_MEMORY_OFFSET in number of elements of the voxel type.
        const size_t DEVICE_MEMORY_OFFSET = 0; 
        // SRC_BRICK_HALO for when source volume bricks contain app compute border voxels that are not wished for rendering
        const mi::math::Vector<mi::Sint32, 3> SRC_BRICK_HALO(0, 0, 0); 
        
        // Size of the work bounding box
        const mi::math::Vector<mi::Sint32, 3> work_size = work_bbox.max - work_bbox.min;
        const mi::math::Vector<mi::Sint32, 3> src_work_size = work_size + 2*SRC_BRICK_HALO;

        cuda::set_voxel_values(
            voxel_values, 
            work_size.begin(), 
            m_bricks_locality[idx].ipc_handle, 
            m_bricks_locality[idx].dev_ptr, 
            src_work_size.begin(),
            DEVICE_MEMORY_OFFSET,
            m_bricks_locality[idx].device_id);

        return true;
    }

    virtual Compute_mode get_compute_mode() const
    {
        return COMPUTE_MODE_GPU;
    }

private:
    mi::math::Bbox<mi::Sint32, 3> get_work_bbox(const mi::math::Bbox<mi::Sint32, 3>& brick_ijk_bbox) const
    {
        if(m_border == 0)
            return brick_ijk_bbox;
        
        // Work bbox should include boundary data
        mi::math::Bbox<mi::Sint32, 3> work_bbox(
            brick_ijk_bbox.min - mi::math::Vector<mi::Sint32, 3>(m_border),
            brick_ijk_bbox.max + mi::math::Vector<mi::Sint32, 3>(m_border));

        // Intersect with the active IJK bbox, which specifies where the operation should be applied
        work_bbox.min = mi::math::elementwise_max(work_bbox.min, mi::math::Vector<mi::Sint32, 3>(m_active_ijk_bbox.min));
        work_bbox.max = mi::math::elementwise_min(work_bbox.max, mi::math::Vector<mi::Sint32, 3>(m_active_ijk_bbox.max));

        // Boundary data on the outer border of the volume is duplicated, so make sure it is
        // updated as well (even when it is outside of the active IJK bbox).
        for (mi::Uint32 i=0; i < 3; ++i)
        {
            if (work_bbox.min[i] == 0)
            {
                work_bbox.min[i]--;
            }

            if (work_bbox.max[i] == static_cast<mi::Sint32>(m_volume_size[i]))
            {
                work_bbox.max[i]++;
            }
        }

        return work_bbox;
    }

    mi::math::Bbox<mi::Uint32, 3>        m_active_ijk_bbox;
    mi::math::Vector<mi::Uint32, 3>      m_volume_size;
    mi::Uint32                           m_border;
    std::vector<Brick_locality>          m_bricks_locality;
};

//----------------------------------------------------------------------
/// Fragmented job which starts one fragment per GPU and executes the compute tasks for all bricks
/// that are stored on that GPU.
class GPU_volume_compute_scheduler :
    public mi::neuraylib::Fragmented_job<0x65f080d,0xe68d,0x437c,0x8f,0x6a,0x37,0xeb,0xb7,0x61,0xfa,0x59>
{
public:
    /// \param host_id           Specifies the current host.
    /// \param data_locality     Contains information on which GPU the data for a certain brick is stored.
    /// \param compute_task      Compute task that should be applied to the bricks
    /// \param dice_transaction  Current transaction.
    GPU_volume_compute_scheduler(
        mi::Uint32                                      host_id,
        nv::index::IRegular_volume_data_locality*       data_locality,
        nv::index::Regular_volume_compute_task*         compute_task,
        mi::neuraylib::IDice_transaction*               dice_transaction)
      : m_host_id(host_id),
        m_data_locality(data_locality),
        m_compute_task(compute_task)
    {
        const mi::Uint32 nb_bounding_boxes = m_data_locality->get_nb_bounding_box(host_id);

        // Create a map that holds which bricks are stored on which device
        std::map<mi::Sint32, std::vector<mi::Uint32> > bricks_per_device_map;
        for (mi::Uint32 i=0; i < nb_bounding_boxes; ++i)
        {
            mi::base::Handle<nv::index::IRegular_volume_data_edit> brick_data_editing(
                m_data_locality->create_data_edit(dice_transaction, host_id, i));

            mi::Sint32 device = brick_data_editing->get_assigned_device(dice_transaction);
            if (device >= 0)
            {
                bricks_per_device_map[device].push_back(i);
            }
        }

        // Convert the map into a vector of vectors
        std::map<mi::Sint32, std::vector<mi::Uint32> >::const_iterator it;
        for (it = bricks_per_device_map.begin(); it != bricks_per_device_map.end(); ++it)
        {
            m_bricks_per_device.push_back(it->second);
        }
    }

    /// Number of fragments to start the job
    mi::Size get_nb_of_fragments() const
    {
        return m_bricks_per_device.size(); // One fragment per device
    }

    /// Process the bricks which are all located on a single device
    virtual void execute_fragment(
        mi::neuraylib::IDice_transaction*               dice_transaction,
        mi::Size                                        index,
        mi::Size                                        count,
        const mi::neuraylib::IJob_execution_context*    context)
    {
        assert(index < m_bricks_per_device.size());

        // Iterate over bricks stored on the current device
        const std::vector<mi::Uint32>& bricks = m_bricks_per_device[index];
        for (size_t i=0; i < bricks.size(); ++i)
        {
            mi::base::Handle<nv::index::IRegular_volume_data_edit> brick_data_editing(
                m_data_locality->create_data_edit(dice_transaction, m_host_id, bricks[i]));

            brick_data_editing->execute_compute_task(m_compute_task, dice_transaction);
        }
    }

    /// Returns the CPU load per fragment of the fragmented job.
    virtual mi::Float32 get_cpu_load() const { return 0.f; }

    /// Returns the GPU load per fragment of the fragmented job.
    virtual mi::Float32 get_gpu_load() const { return 1.f; }

private:
    mi::Uint32                                              m_host_id;
    nv::index::IRegular_volume_data_locality*               m_data_locality;
    nv::index::Regular_volume_compute_task*                 m_compute_task;
    std::vector<std::vector<mi::Uint32> >                   m_bricks_per_device;
};

//----------------------------------------------------------------------
CUDA_IPC_volume_editing::CUDA_IPC_volume_editing(
    const mi::neuraylib::Tag&                                       session,
    const mi::neuraylib::Tag&                                       volume_scene_element,
    const mi::math::Bbox<mi::Uint32, 3>&                            active_ijk_bbox,
    const mi::math::Vector<mi::Uint32, 3>&                          volume_size,
    mi::Uint32                                                      border,
    nv::index::IRegular_volume_data_locality*                       data_locality,
    const std::vector<Brick_locality>&                              bricks_locality)
  : m_session_tag(session),
    m_volume_scene_element_tag(volume_scene_element),
    m_active_ijk_bbox(active_ijk_bbox),
    m_volume_size(volume_size),
    m_border(border),
    m_data_locality(data_locality),
    m_nb_bricks(0),
    m_bricks_locality(bricks_locality)
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
CUDA_IPC_volume_editing::CUDA_IPC_volume_editing()
  : m_session_tag(),
    m_volume_scene_element_tag(),
    m_active_ijk_bbox(),
    m_volume_size(),
    m_border(0),
    m_data_locality(0),
    m_nb_bricks(0)
{
    // for serialization only
}

// -----------------------------------------------------------------------------------------
mi::Size CUDA_IPC_volume_editing::get_nb_of_fragments() const
{
    // One fragment per host
    return m_data_locality->get_nb_cluster_nodes();
}

// -----------------------------------------------------------------------------------------
void CUDA_IPC_volume_editing::create_remote_data_locality(
    mi::neuraylib::IDice_transaction*   dice_transaction)
{
    // Unfortunately, we cannot serialize the data locality information (abstract interface)
    // and thus have to re-create it.
    //
    // TODO: Design a factory scheme that can be used anywhere
    //       without the need to have a dice transaction at hand.
    //
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
            m_volume_scene_element_tag, m_active_ijk_bbox, dice_transaction);
    }
}

// -----------------------------------------------------------------------------------------
void CUDA_IPC_volume_editing::assign_fragments_to_hosts(
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
void CUDA_IPC_volume_editing::invoke_compute_tasks_per_host(
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Uint32                                      host_id)
{
    Voxel_value_set_operation op(m_active_ijk_bbox, m_volume_size, m_border, m_bricks_locality);

    // Edit the bricks on the given host by running a fragmented job that makes sure the compute
    // task for each brick is executed on the GPU where the brick data is stored
    GPU_volume_compute_scheduler scheduler(host_id, m_data_locality, &op, dice_transaction);
    dice_transaction->execute_fragmented(&scheduler, scheduler.get_nb_of_fragments());
}

//----------------------------------------------------------------------
void CUDA_IPC_volume_editing::invoke_compute_tasks_per_brick(
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Uint32                                      host_id,
    mi::Uint32                                      brick_id)
{
    // Edit the brick by applying a compute task
    mi::base::Handle<nv::index::IRegular_volume_data_edit> brick_data_editing(
        m_data_locality->create_data_edit(dice_transaction, host_id, brick_id));

    Voxel_value_set_operation op(m_active_ijk_bbox, m_volume_size, m_border, m_bricks_locality);
    brick_data_editing->execute_compute_task(&op, dice_transaction);
}

// --------------------------------------------------------------------------------------------
void CUDA_IPC_volume_editing::execute_fragment(
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
void CUDA_IPC_volume_editing::execute_fragment_remote(
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
void CUDA_IPC_volume_editing::serialize(mi::neuraylib::ISerializer* serializer) const
{
    serializer->write(&m_session_tag.id, 1);
    serializer->write(&m_volume_scene_element_tag.id, 1);
    serializer->write(&m_active_ijk_bbox.min.x, 6);
    serializer->write(&m_volume_size.x, 3);

    serializer->write(&m_nb_bricks, 1);
    
    mi::Uint32 nb_elements = static_cast<mi::Uint32>(m_host_ids.size());
    serializer->write(&nb_elements, 1);
    for (mi::Uint32 i=0; i < nb_elements; ++i)
    {
        serializer->write(&m_host_ids[i], 1);
    }
        
    nb_elements = static_cast<mi::Uint32>(m_brick_ids.size());
    serializer->write(&nb_elements, 1);
    for (mi::Uint32 i=0; i < nb_elements; ++i)
    {
        serializer->write(&m_brick_ids[i], 1);
    }
    
    nb_elements = static_cast<mi::Uint32>(m_bricks_locality.size());
    serializer->write(&nb_elements, 1);
    const mi::Uint8* bricks_info_bytes = reinterpret_cast<const mi::Uint8*>(&m_bricks_locality[0]);
    serializer->write(bricks_info_bytes, sizeof(Brick_locality)*nb_elements);
    
}

//----------------------------------------------------------------------
void CUDA_IPC_volume_editing::deserialize(mi::neuraylib::IDeserializer* deserializer)
{
    deserializer->read(&m_session_tag.id, 1);
    deserializer->read(&m_volume_scene_element_tag.id, 1);
    deserializer->read(&m_active_ijk_bbox.min.x, 6);
    deserializer->read(&m_volume_size.x, 3);

    deserializer->read(&m_nb_bricks, 1);
    
    mi::Uint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_host_ids.resize(nb_elements);
    for (mi::Uint32 i=0; i < nb_elements; ++i)
    {
        deserializer->read(&m_host_ids[i], 1);
    }

    nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_brick_ids.resize(nb_elements);
    for (mi::Uint32 i=0; i < nb_elements; ++i)
    {
        deserializer->read(&m_brick_ids[i], 1);
    }
    
    nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_bricks_locality.resize(nb_elements);
    mi::Uint8* bricks_info_bytes = reinterpret_cast<mi::Uint8*>(&m_bricks_locality[0]);
    deserializer->read(bricks_info_bytes, sizeof(Brick_locality)*nb_elements);
}

void CUDA_IPC_volume_editing::init(const nv::index_common::String_dict* prj)
{

}

void CUDA_IPC_volume_editing::step(
    mi::neuraylib::Tag                                           volume_tag,
    mi::neuraylib::Tag                                           session_tag,
    mi::neuraylib::IDice_transaction*                            dice_transaction,
    const std::vector<Brick_locality>&                           bricks_locality)
{
    ERROR_LOG << "Initializing CUDA IPC volume editing....................................";
    
    // Access the session
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<nv::index::ISession>(session_tag));

    // Access the volume scene element
    mi::base::Handle<const nv::index::IRegular_volume> volume(
        dice_transaction->access<nv::index::IRegular_volume>(volume_tag));

    // The entire region of interest is used for computing the volume data
    const mi::math::Bbox<mi::Uint32, 3> query_bbox = volume->get_IJK_bounding_box();

    // Access the distribution scheme
    const mi::neuraylib::Tag dist_layout_tag = session->get_distribution_layout();
    mi::base::Handle<const nv::index::IData_distribution> distribution_layout(
        dice_transaction->access<nv::index::IData_distribution>(dist_layout_tag));

    // Retrieve configuration to get the border size
    mi::base::Handle<const nv::index::IConfig_settings> config_settings(
        dice_transaction->access<nv::index::IConfig_settings>(session->get_config()));

    if (config_settings->get_subcube_border_size() != 0)
    {
        ERROR_LOG << "CUDA IPC volume editing requires border size 0 but got "
                  << config_settings->get_subcube_border_size();
        return;
    }

    // Distribution layout
    mi::base::Handle<nv::index::IRegular_volume_data_locality> data_locality(
        distribution_layout->retrieve_data_locality(volume_tag, query_bbox, dice_transaction));

    // Set up distributed editing process.
    mi::base::Handle<nv::index::IDistributed_compute_algorithm> task(
        new CUDA_IPC_volume_editing(
            session_tag,
            volume_tag,
            query_bbox,
            volume->get_volume_size(),
            0,
            data_locality.get(),
            bricks_locality));

    // Start the fragmented job
    dice_transaction->execute_fragmented(task.get(), task->get_nb_of_fragments());
}

mi::Uint64 CUDA_IPC_volume_editing::get_allocated_device_memory()
{
    return 0;
}

