/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief heightfiled elevation change example

#ifndef NVIDIA_INDEX_DISTRIBUTED_HEIGHTFIELD_ELEVATION_CHANGE_H
#define NVIDIA_INDEX_DISTRIBUTED_HEIGHTFIELD_ELEVATION_CHANGE_H

#include <mi/base/types.h>
#include <mi/dice.h>
#include <mi/base/interface_declare.h>

#include <nv/index/idistributed_compute_algorithm.h>
#include <nv/index/idistributed_data_locality.h>

#include <vector>

namespace nv {
namespace index_common {

/// Proof of concept computing algorithm that simply changes
/// the elevation values of the given heightfield.
class Distributed_heightfield_elevation_change :
    public nv::index::Distributed_compute_algorithm<0x0d03bab9,0x4d21,0x4e59,0x80,0x95,0x43,0x4b,0x4b,0x7e,0xe9,0xff>
{
public:
    Distributed_heightfield_elevation_change(
        const mi::neuraylib::Tag&                                       session_tag,
        const mi::neuraylib::Tag&                                       heightfield_tag,
        bool                                                            scale_operation,
        mi::Float32                                                     elevation_manipulation_value,
        const std::vector<mi::math::Vector_struct<mi::Float32, 2> >&    polygon,
        const std::vector<mi::Uint32>&                                  cluster_hosts);

    Distributed_heightfield_elevation_change()
        :
        m_session_tag(mi::neuraylib::NULL_TAG),
        m_heightfield_tag(mi::neuraylib::NULL_TAG),
        m_scale_operation(false),
        m_elevation_manipulation_value(0.0f),
        m_polygon_of_interest(),
        m_cluster_hosts(),
        m_fragment_bounding_boxes(),
        m_data_lock()
    {
        // for serialization only
    }

    virtual ~Distributed_heightfield_elevation_change();

    /// Number of fragments to start the compute algorithm
    virtual mi::Size get_nb_of_fragments() const { return m_cluster_hosts.size(); }

    virtual void get_updated_bounding_box(mi::math::Bbox_struct<mi::Float32, 3>& bbox) const;

public:
    // Implemented Fragmented_job

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

    virtual void receive_remote_result(
        mi::neuraylib::IDeserializer*       deserializer,
        mi::neuraylib::IDice_transaction*   dice_transaction,
        mi::Size                            index,
        mi::Size                            count);

    virtual void serialize(mi::neuraylib::ISerializer *serializer) const;

    virtual void deserialize(mi::neuraylib::IDeserializer* deserializer);

private:
    mi::math::Bbox<mi::Float32, 3> edit_heightfield_data(
        mi::neuraylib::IDice_transaction*               dice_transaction,
        nv::index::IRegular_heightfield_data_locality*  locality,
        mi::Uint32                                      host_id);

private:
    mi::neuraylib::Tag                                      m_session_tag;
    mi::neuraylib::Tag                                      m_heightfield_tag;
    bool                                                    m_scale_operation;
    mi::Float32                                             m_elevation_manipulation_value;
    std::vector<mi::math::Vector_struct<mi::Float32, 2> >   m_polygon_of_interest;

    // First implemented the user-defined scheduling.
    // Later used to access resp. edit the heightfield data on local or remote host.
    std::vector<mi::Uint32>                                 m_cluster_hosts;

    // Do not serialize/deserialize below
    std::vector<mi::math::Bbox<mi::Float32, 3> >            m_fragment_bounding_boxes;
    mi::base::Lock                                          m_data_lock;
};

}} // namespace nv::index_common

#endif // NVIDIA_INDEX_DISTRIBUTED_HEIGHTFIELD_ELEVATION_CHANGE_H
