/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include <cassert>

#include <nv/index/isession.h>
#include <nv/index/iregular_volume.h>
#include <nv/index/iregular_volume_compute_task.h>

#include "forwarding_logger.h"
#include "distributed_voxel_value_change.h"

namespace nv {
namespace index_common {

/// Volume data value set compute task
class Volume_voxel_value_set_operation : public nv::index::Regular_volume_compute_task
{
public:
    /// Constructor
    ///
    /// \param[in] voxel_value  Voxel value to set
    /// \param[in] ijk_bbox     Specifies the part of the volume where the operation should be applied.
    /// \param[in] volume_size  Size of the entire volume.
    Volume_voxel_value_set_operation(
        mi::Sint8                              voxel_value,
        const mi::math::Bbox<mi::Uint32, 3>&   ijk_bbox,
        const mi::math::Vector<mi::Uint32, 3>& volume_size)
        :
        m_voxel_value(voxel_value),
        m_active_ijk_bbox(ijk_bbox),
        m_volume_size(volume_size)
    {
        // empty
    }

    /// Perform computation on the given volume data.
    ///
    /// \param[in] brick_ijk_bbox_struct IJK bounding box of the current volume brick.
    /// \param[in] voxel_data            Voxel data of the volume brick.
    /// \param[in] dice_transaction      Current dice db transaction.
    virtual bool compute(
        const mi::math::Bbox_struct<mi::Uint32, 3>& brick_ijk_bbox_struct,
        nv::index::IRegular_volume_data*            voxel_data,
        mi::neuraylib::IDice_transaction*           dice_transaction) const
    {
        assert(voxel_data != 0);
        assert(dice_transaction != 0);
        
        mi::base::Handle<nv::index::IRegular_volume_data_uint8> volume_data_uint8(voxel_data->get_interface<nv::index::IRegular_volume_data_uint8>());
        mi::base::Handle<nv::index::IRegular_volume_data_rgba8> volume_data_rgba8(voxel_data->get_interface<nv::index::IRegular_volume_data_rgba8>());
        mi::base::Handle<nv::index::IRegular_volume_data_uint16> volume_data_uint16(voxel_data->get_interface<nv::index::IRegular_volume_data_uint16>());
        mi::base::Handle<nv::index::IRegular_volume_data_float32> volume_data_float32(voxel_data->get_interface<nv::index::IRegular_volume_data_float32>());
    
        const mi::math::Bbox<mi::Sint32, 3> brick_ijk_bbox(
            brick_ijk_bbox_struct.min.x, brick_ijk_bbox_struct.min.y, brick_ijk_bbox_struct.min.z,
            brick_ijk_bbox_struct.max.x, brick_ijk_bbox_struct.max.y, brick_ijk_bbox_struct.max.z);

        // Work bbox should include the boundary data
        mi::math::Bbox<mi::Sint32, 3> work_bbox(
            brick_ijk_bbox.min - mi::math::Vector<mi::Sint32, 3>(1),        // FIXME: The boundary is not necessary just 1.
            brick_ijk_bbox.max + mi::math::Vector<mi::Sint32, 3>(1));       // FIXME: The boundary is not necessary just 1.

        // Intersect with the active IJK bbox, which specifies where the operation should be applied
        work_bbox.min = mi::math::elementwise_max(work_bbox.min, mi::math::Vector<mi::Sint32, 3>(m_active_ijk_bbox.min));
        work_bbox.max = mi::math::elementwise_min(work_bbox.max, mi::math::Vector<mi::Sint32, 3>(m_active_ijk_bbox.max));

        // The boundary data at the outer border of the regular volume dataset is duplicated, so make sure it is
        // updated as well even if it is outside of the active IJK bbox.
        for (mi::Sint32 i=0; i < 3; ++i)
        {
            if (work_bbox.min[i] == 0)
            {
                work_bbox.min[i]--;
            }

            if (work_bbox.max[i] == (mi::Sint32)m_volume_size[i])
            {
                work_bbox.max[i]++;
            }
        }

        if (!work_bbox.is_volume())
        {
            return false;
        }

        // Size of the work bounding box
        const mi::math::Vector<mi::Sint32, 3> work = work_bbox.max - work_bbox.min;

        // Offset into the brick data
        const mi::math::Vector<mi::Sint32, 3> offset = work_bbox.min - brick_ijk_bbox.min + mi::math::Vector<mi::Sint32, 3>(1); // FIXME: The boundary is not necessary just 1.

        // Size of the brick (for addressing)
        const mi::math::Vector<mi::Sint32, 3> d = brick_ijk_bbox.max - brick_ijk_bbox.min + mi::math::Vector<mi::Sint32, 3>(2); // FIXME: The boundary is not necessary just 1.
        
        if(volume_data_uint8)
        {
            mi::Uint8* voxel_values = volume_data_uint8->get_voxel_data_mutable();
            
            // Now assign the voxel to all selected voxels
            for (mi::Sint32 x=0; x < work.x; ++x)
            {
                for (mi::Sint32 y=0; y < work.y; ++y)
                {
                    for (mi::Sint32 z=0; z < work.z; ++z)
                    {
                        mi::Uint32 p = ((x + offset.x) * d.y * d.z) + ((y + offset.y) * d.z) + z + offset.z;
                        voxel_values[p] = m_voxel_value;
                    }
                }
            }

            return true;
        }
        else if(volume_data_rgba8)
        {
            mi::math::Vector_struct<mi::Uint8, 4>* voxel_values = volume_data_rgba8->get_voxel_data_mutable();
            
            // Now assign the voxel value to all selected voxels
            for (mi::Sint32 x=0; x < work.x; ++x)
            {
                for (mi::Sint32 y=0; y < work.y; ++y)
                {
                    for (mi::Sint32 z=0; z < work.z; ++z)
                    {
                        mi::Uint32 p = ((x + offset.x) * d.y * d.z) + ((y + offset.y) * d.z) + z + offset.z;
                        voxel_values[p].x = m_voxel_value;
                        voxel_values[p].y = m_voxel_value;
                        voxel_values[p].z = m_voxel_value;
                        voxel_values[p].w = 255u;
                    }
                }
            }

            return true;
        }
        else if(volume_data_uint16)
        {
            mi::Uint16* voxel_values = volume_data_uint16->get_voxel_data_mutable();
            
            // Now assign the voxel to all selected voxels
            for (mi::Sint32 x=0; x < work.x; ++x)
            {
                for (mi::Sint32 y=0; y < work.y; ++y)
                {
                    for (mi::Sint32 z=0; z < work.z; ++z)
                    {
                        mi::Uint32 p = ((x + offset.x) * d.y * d.z) + ((y + offset.y) * d.z) + z + offset.z;
                        voxel_values[p] = m_voxel_value;
                    }
                }
            }

            return true;
        }
        else if(volume_data_float32)
        {
            mi::Float32* voxel_values = volume_data_float32->get_voxel_data_mutable();
            
            // Now assign the voxel to all selected voxels
            for (mi::Sint32 x=0; x < work.x; ++x)
            {
                for (mi::Sint32 y=0; y < work.y; ++y)
                {
                    for (mi::Sint32 z=0; z < work.z; ++z)
                    {
                        mi::Uint32 p = ((x + offset.x) * d.y * d.z) + ((y + offset.y) * d.z) + z + offset.z;
                        voxel_values[p] = static_cast<mi::Float32>(m_voxel_value);
                    }
                }
            }

            return true;
        }
        else
        {
            ERROR_LOG << "Regular volume editing operation failed because the voxel datatype is neither uint8 nor rgba8.";
            return false;
        }
        
    }

private:
    mi::Sint8                            m_voxel_value;
    mi::math::Bbox<mi::Uint32, 3>        m_active_ijk_bbox;
    mi::math::Vector<mi::Uint32, 3>      m_volume_size;
};

//----------------------------------------------------------------------
Distributed_voxel_value_change::Distributed_voxel_value_change(
    const mi::neuraylib::Tag&              session_tag,
    const mi::neuraylib::Tag&              volume_tag,
    bool                                   is_host_assign_mode,
    mi::Uint8                              voxel_value,
    const mi::math::Bbox<mi::Uint32, 3>&   active_ijk_bbox,
    const mi::math::Vector<mi::Uint32, 3>& volume_size,
    const std::vector<mi::Uint32>&         cluster_hosts)
  : m_session_tag(session_tag),
    m_volume_tag(volume_tag),
    m_is_host_assign_mode(is_host_assign_mode),
    m_voxel_value(voxel_value),
    m_active_ijk_bbox(active_ijk_bbox),
    m_volume_size(volume_size),
    m_cluster_hosts(cluster_hosts)
{
    INFO_LOG << "The operation applied the distributed regular volume's voxel data will set all voxel values contained in "
             << active_ijk_bbox << " to the following value: " << ((mi::Uint32)voxel_value) << ".";
}

// -----------------------------------------------------------------------------------------

void Distributed_voxel_value_change::assign_fragments_to_hosts(
    mi::Uint32* slots,
    mi::Size    nr_slots)
{
    if (nr_slots != m_cluster_hosts.size())
    {
        ERROR_LOG << "The distributed regular volume's voxel data cannot be edited. "
                  << "The number of hosts in the cluster (" << m_cluster_hosts.size()
                  << ") and the number of hosts that were requested to perform the operation (" << nr_slots << ") doesn't match.";
    }
    else
    {
        for (mi::Size i=0; i < nr_slots; ++i)
        {
            slots[i] = m_cluster_hosts[i];
        }
    }
}

// --------------------------------------------------------------------------------------------

void Distributed_voxel_value_change::serialize(
    mi::neuraylib::ISerializer* serializer) const
{
    serializer->write(&m_session_tag.id, 1);
    serializer->write(&m_volume_tag.id, 1);
    serializer->write(&m_is_host_assign_mode, 1);
    serializer->write(&m_voxel_value, 1);
    serializer->write(&m_active_ijk_bbox.min.x, 6);
    serializer->write(&m_volume_size.x, 3);

    const mi::Uint32 nb_elements = mi::Uint32(m_cluster_hosts.size());
    serializer->write(&nb_elements, 1);
    for (mi::Uint32 i=0; i < nb_elements; ++i)
    {
        serializer->write(&m_cluster_hosts[i], 1);
    }
}

void Distributed_voxel_value_change::deserialize(
    mi::neuraylib::IDeserializer* deserializer)
{
    deserializer->read(&m_session_tag.id, 1);
    deserializer->read(&m_volume_tag.id, 1);
    deserializer->read(&m_is_host_assign_mode, 1);
    deserializer->read(&m_voxel_value, 1);
    deserializer->read(&m_active_ijk_bbox.min.x, 6);
    deserializer->read(&m_volume_size.x, 3);

    mi::Uint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_cluster_hosts.resize(nb_elements);
    for (mi::Uint32 i=0; i < nb_elements; ++i)
    {
        deserializer->read(&m_cluster_hosts[i], 1);
    }
}

// --------------------------------------------------------------------------------------------

void Distributed_voxel_value_change::edit_regular_volume_data(
    mi::neuraylib::IDice_transaction*               dice_transaction,
    nv::index::IRegular_volume_data_locality*       locality,
    mi::Uint32                                      host_id)
{
    const mi::Uint32 nb_brick_edits = static_cast<mi::Uint32>(locality->get_nb_bounding_box(host_id));

    for (mi::Uint32 i=0; i < nb_brick_edits; ++i)
    {
        // const mi::math::Bbox_struct<mi::Uint32, 3>& roi = locality->get_bounding_box(host_id, i);

        mi::base::Handle<nv::index::IRegular_volume_data_edit> volume_data_edit(
            locality->create_data_edit(dice_transaction, host_id, i));

        mi::Sint8 setting_voxel_value = 0;
        if(m_is_host_assign_mode)
        {
            // use host_id related value for m_is_host_assign_mode == true
            setting_voxel_value = static_cast<mi::Sint8>(host_id % 256);
        }
        else
        {
            setting_voxel_value = m_voxel_value;
        }

        mi::base::Handle<Volume_voxel_value_set_operation> op(
            new Volume_voxel_value_set_operation(setting_voxel_value, m_active_ijk_bbox, m_volume_size));
        {
            volume_data_edit->execute_compute_task(op.get(), dice_transaction);
        }
    }
}

// --------------------------------------------------------------------------------------------

namespace {

nv::index::IRegular_volume_data_locality* get_regular_volume_distribution_layout(
    const mi::neuraylib::Tag&                   session_tag,
    const mi::neuraylib::Tag&                   volume_tag,
    const mi::math::Bbox_struct<mi::Uint32, 3>& ijk_query_bounds,
    mi::neuraylib::IDice_transaction *          dice_transaction)
{
    assert(session_tag.is_valid());
    assert(volume_tag.is_valid());
    assert(dice_transaction != 0);

    mi::base::Handle<const nv::index::ISession> session(dice_transaction->access<nv::index::ISession>(session_tag));

    // Access the distribution scheme
    const mi::neuraylib::Tag dist_layout_tag = session->get_distribution_layout();
    mi::base::Handle<const nv::index::IData_distribution> distribution_layout(dice_transaction->access<nv::index::IData_distribution>(dist_layout_tag));

    // Regular volume data distribution layout.
    nv::index::IRegular_volume_data_locality* data_locality = distribution_layout->retrieve_data_locality(
        volume_tag, ijk_query_bounds, dice_transaction);

    return data_locality;
}

} // namespace

//----------------------------------------------------------------------

void Distributed_voxel_value_change::execute_fragment(
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Size                                        index,
    mi::Size                                        count,
    const mi::neuraylib::IJob_execution_context*    context)
{
    mi::base::Handle<nv::index::IRegular_volume_data_locality> locality(
        get_regular_volume_distribution_layout(m_session_tag, m_volume_tag, m_active_ijk_bbox, dice_transaction));

    const mi::Uint32 cluster_host = m_cluster_hosts[index];
    DEBUG_LOG << "The distributed data editing of the regular volume started on local host (id) " << cluster_host << ".";
    edit_regular_volume_data(dice_transaction, locality.get(), cluster_host);
}

void Distributed_voxel_value_change::execute_fragment_remote(
    mi::neuraylib::ISerializer*                     serializer,
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Size                                        index,
    mi::Size                                        count,
    const mi::neuraylib::IJob_execution_context*    context)
{
    mi::base::Handle<nv::index::IRegular_volume_data_locality> locality(
        get_regular_volume_distribution_layout(m_session_tag, m_volume_tag, m_active_ijk_bbox, dice_transaction));

    const mi::Uint32 cluster_host = m_cluster_hosts[index];
    DEBUG_LOG << "The distributed data editing of the regular volume started on remote host (id) " << cluster_host << ".";
    edit_regular_volume_data(dice_transaction, locality.get(), cluster_host);
}

}} // namespace nv::index_common
