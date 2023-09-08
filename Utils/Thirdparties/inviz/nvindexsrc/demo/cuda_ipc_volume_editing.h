/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Volume editing via CUDA IPC 

#ifndef NVIDIA_INDEX_CUDA_IPC_VOLUME_EDITING_H
#define NVIDIA_INDEX_CUDA_IPC_VOLUME_EDITING_H

#include <mi/base/types.h>
#include <mi/dice.h>
#include <mi/base/interface_declare.h>

#include <cuda_runtime.h>

#include <nv/index/idistributed_compute_algorithm.h>
#include <nv/index/idistributed_data_locality.h>

#include <vector>

#include "common/string_dict.h"


/// Brick locality information
struct Brick_locality
{
    mi::math::Bbox<mi::Sint32, 3>   bbox;
    cudaIpcMemHandle_t              ipc_handle;
    mi::Uint8*                      dev_ptr;                                   
    mi::Uint32                      brick_id;
    mi::Uint32                      host_id;
    mi::Uint32                      device_id;
};

/// Brick locality information
struct Brick_affinity
{
    mi::math::Bbox<mi::Float32, 3>  bbox;
    mi::Uint32                      host_id;
    mi::Uint32                      device_id;
};

/// Simple distributed/parallel volume editing operation that is used to transfer
/// volume bricks allocated by the application to IndeX internal volume brick storage
/// without using host memory. All transfers are done device to device using CUDA IPC.

class CUDA_IPC_volume_editing :
    public nv::index::Distributed_compute_algorithm<0xc2743f80,0x2a72,0x45f8,0x9a,0x86,0xf9,0x92,0x23,0x18,0x5f,0xbe>
{
public:
    /// \param session                Tag of the session
    /// \param volume_tag             Tag of the volume scene element
    /// \param active_ijk_bbox        Specifies the part of the volume where the operation should be applied
    /// \param volume_size            Size of the entire volume
    /// \param border                 Size of boundary around each subcube
    /// \param data_locality          Locality information of the volume data
    CUDA_IPC_volume_editing(
        const mi::neuraylib::Tag&                   session,
        const mi::neuraylib::Tag&                   volume_tag,
        const mi::math::Bbox<mi::Uint32, 3>&        active_ijk_bbox,
        const mi::math::Vector<mi::Uint32, 3>&      volume_size,
        mi::Uint32                                  border,
        nv::index::IRegular_volume_data_locality*   data_locality,
        const std::vector<Brick_locality>&          bricks_locality);

    // Default constructor for serialization only
    CUDA_IPC_volume_editing();

public:
    virtual mi::Size get_nb_of_fragments() const;

    virtual mi::neuraylib::IFragmented_job::Scheduling_mode get_scheduling_mode() const
    {
        return mi::neuraylib::IFragmented_job::USER_DEFINED;
    }

    virtual void assign_fragments_to_hosts(
        mi::Uint32* slots,
        mi::Size    nr_slots);

    virtual void execute_fragment(
        mi::neuraylib::IDice_transaction*               dice_transaction,
        mi::Size                                        index,
        mi::Size                                        count,
        const mi::neuraylib::IJob_execution_context*    context);

    virtual void execute_fragment_remote(
        mi::neuraylib::ISerializer*                     serializer,
        mi::neuraylib::IDice_transaction*               dice_transaction,
        mi::Size                                        index,
        mi::Size                                        count,
        const mi::neuraylib::IJob_execution_context*    context);

    virtual void serialize(mi::neuraylib::ISerializer *serializer) const;
    virtual void deserialize(mi::neuraylib::IDeserializer* deserializer);
    virtual void get_updated_bounding_box(mi::math::Bbox_struct<mi::Float32, 3>& bbox) const
    {
        // empty
    }
    
public:
    static void init(const nv::index_common::String_dict* prj);

    static void step(
        mi::neuraylib::Tag                                           volume_tag,
        mi::neuraylib::Tag                                           session_tag,
        mi::neuraylib::IDice_transaction*                            dice_transaction,
        const std::vector<Brick_locality>&                           bricks_locality);

    static mi::Uint64 get_allocated_device_memory();
    

private:
    void create_remote_data_locality(
        mi::neuraylib::IDice_transaction*                   dice_transaction);

    void invoke_compute_tasks_per_host(
        mi::neuraylib::IDice_transaction*                   dice_transaction,
        mi::Uint32                                          host_id);

    void invoke_compute_tasks_per_brick(
        mi::neuraylib::IDice_transaction*                   dice_transaction,
        mi::Uint32                                          host_id,
        mi::Uint32                                          brick_id);

private:
    mi::neuraylib::Tag                                      m_session_tag;
    mi::neuraylib::Tag                                      m_volume_scene_element_tag;
    mi::Uint8                                               m_voxel_value;
    mi::math::Bbox<mi::Uint32, 3>                           m_active_ijk_bbox;
    mi::math::Vector<mi::Uint32, 3>                         m_volume_size;
    mi::Uint32                                              m_border;

    // User-defined scheduling.
    nv::index::IRegular_volume_data_locality*               m_data_locality;
    mi::Uint32                                              m_nb_bricks;

    std::vector<mi::Uint32>                                 m_host_ids;
    std::vector<mi::Uint32>                                 m_brick_ids;

    // Need a lock when creating a data locality.
    mi::base::Lock                                          m_lock;
    
    // brick information
    std::vector<Brick_locality>                             m_bricks_locality;
    
};

#endif // NVIDIA_INDEX_CUDA_IPC_VOLUME_EDITING_H
