/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief voxel value change computation on a volume data 

#ifndef NVIDIA_INDEX_DISTRIBUTED_VOXEL_VALUE_CHANGE_H
#define NVIDIA_INDEX_DISTRIBUTED_VOXEL_VALUE_CHANGE_H

#include <vector>

#include <mi/base/types.h>
#include <mi/dice.h>
#include <mi/base/interface_declare.h>

#include <nv/index/idistributed_compute_algorithm.h>
#include <nv/index/idistributed_data_locality.h>

namespace nv {
namespace index_common {

/// Simple computing algorithm that sets all voxel values in a given
/// region to a constant value.
class Distributed_voxel_value_change :
    public mi::neuraylib::Fragmented_job<0x92442ea6,0x9cfe,0x4c87,0x8c,0xeb,0xc3,0x71,0x60,0xd2,0x81,0x57,
                                         nv::index::IDistributed_compute_algorithm>
{
public:
    /// Constructor
    /// \param session_tag         Tag of the session
    /// \param volume_tag          Tag of the volume scene element
    /// \param is_host_assign_mode Host assignment mode when true, and ignore the voxel value
    /// \param voxel_value     Value to set when is_host_assign_mode is false. 
    /// \param active_ijk_bbox     Specifies the part of the volume where the operation should be applied
    /// \param volume_size         Size of the entire volume
    /// \param cluster_hosts
    Distributed_voxel_value_change(
        const mi::neuraylib::Tag&              session_tag,
        const mi::neuraylib::Tag&              volume_tag,
        bool                                   is_host_assign_mode,
        mi::Uint8                              voxel_value,
        const mi::math::Bbox<mi::Uint32, 3>&   active_ijk_bbox,
        const mi::math::Vector<mi::Uint32, 3>& volume_size,
        const std::vector<mi::Uint32>&         cluster_hosts);

    Distributed_voxel_value_change()
        :
        m_session_tag(mi::neuraylib::NULL_TAG),
        m_volume_tag(mi::neuraylib::NULL_TAG),
        m_is_host_assign_mode(false),
        m_voxel_value(0),
        m_active_ijk_bbox(),
        m_volume_size(),
        m_cluster_hosts()
    {
        // for serialization only
    }

public:
    // Implemented Fragmented_job

    virtual mi::Size get_nb_of_fragments() const { return m_cluster_hosts.size(); }

    virtual void get_updated_bounding_box(mi::math::Bbox_struct<mi::Float32, 3>& bbox) const { /* empty */ }

    virtual mi::neuraylib::IFragmented_job::Scheduling_mode get_scheduling_mode() const { return mi::neuraylib::IFragmented_job::USER_DEFINED;}

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

    virtual mi::base::Uuid get_class_id() const { return IID(); }

    virtual void serialize(mi::neuraylib::ISerializer *serializer) const;

    virtual void deserialize(mi::neuraylib::IDeserializer* deserializer);

private:
    void edit_regular_volume_data(
        mi::neuraylib::IDice_transaction*               dice_transaction,
        nv::index::IRegular_volume_data_locality*       locality,
        mi::Uint32                                      host_id);

private:
    mi::neuraylib::Tag              m_session_tag;
    mi::neuraylib::Tag              m_volume_tag;
    bool                            m_is_host_assign_mode;
    mi::Uint8                       m_voxel_value;
    mi::math::Bbox<mi::Uint32, 3>   m_active_ijk_bbox;
    mi::math::Vector<mi::Uint32, 3> m_volume_size;

    // User-defined scheduling.
    std::vector<mi::Uint32>         m_cluster_hosts;
};

}} // namespace nv::index_common

#endif // NVIDIA_INDEX_DISTRIBUTED_VOXEL_VALUE_CHANGE_H
