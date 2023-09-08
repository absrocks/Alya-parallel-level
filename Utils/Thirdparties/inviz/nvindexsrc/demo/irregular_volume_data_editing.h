/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/

#ifndef NVIDIA_INDEX_IRREGULAR_VOLUME_DATA_EDITING_H
#define NVIDIA_INDEX_IRREGULAR_VOLUME_DATA_EDITING_H

#include <mi/base/types.h>
#include <mi/dice.h>
#include <mi/base/interface_declare.h>

#include <nv/index/idistributed_compute_algorithm.h>
#include <nv/index/idistributed_data_locality.h>

#include <vector>

class Irregular_volume_data_editing :
    public nv::index::Distributed_compute_algorithm<0xa2115911,0x73ef,0x48cb,0x97,0x0f,0x5a,0x48,0xfc,0x1e,0x72,0x29>
{
public:
    Irregular_volume_data_editing(
        mi::Uint32                                  time_step,
        const mi::neuraylib::Tag&                   session,
        const mi::neuraylib::Tag&                   volume_tag,
        const mi::math::Bbox<mi::Float32, 3>&       active_ijk_bbox,
        nv::index::IIrregular_volume_data_locality* data_locality);

    // Default constructor for serialization only
    Irregular_volume_data_editing();

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

    mi::Uint32                                              m_time_step;
    mi::neuraylib::Tag                                      m_session_tag;
    mi::neuraylib::Tag                                      m_scene_element_tag;
    mi::math::Bbox<mi::Float32, 3>                          m_active_ijk_bbox;

    // User-defined scheduling.
    nv::index::IIrregular_volume_data_locality*             m_data_locality;
    mi::Uint32                                              m_nb_bricks;

    std::vector<mi::Uint32>                                 m_host_ids;
    std::vector<mi::Uint32>                                 m_brick_ids;

    // Need a lock when creating a data locality.
    mi::base::Lock                                          m_lock;
};

#endif // NVIDIA_INDEX_IRREGULAR_VOLUME_DATA_EDITING_H
