/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief heightfield elevation value delete

#ifndef NVIDIA_INDEX_DISTRIBUTED_HORIZON_SEEDING_UPDATES_H
#define NVIDIA_INDEX_DISTRIBUTED_HORIZON_SEEDING_UPDATES_H

#include <mi/base/types.h>
#include <mi/dice.h>
#include <mi/base/interface_declare.h>

#include <nv/index/idistributed_compute_algorithm.h>
#include <nv/index/idistributed_data_locality.h>

#include <mi/base/lock.h>

#include <vector>


/// Proof of concept computing algorithm that simply changes
/// the elevation values of the given heightfield.
class Distributed_heightfield_seedings_updates :
    public nv::index::Distributed_compute_algorithm<0x68c13b6a,0x9f5d,0x40a1,0x86,0x8b,0x2b,0x1b,0xd9,0x96,0xb7,0xe0>
{
public:
    Distributed_heightfield_seedings_updates(
        const mi::neuraylib::Tag&                                       session_tag,
        const mi::neuraylib::Tag&                                       heightfield_tag,
        const std::vector<mi::Uint32>&                                  cluster_hosts);
    Distributed_heightfield_seedings_updates()
        :
        m_session_tag(mi::neuraylib::NULL_TAG),
        m_heightfield_tag(mi::neuraylib::NULL_TAG),
        m_cluster_hosts(),
        m_patch_tag_pair_vector(),
        m_compute_lock(),
        m_max_value(0.0f),
        m_min_value(0.0f)
    {
        // for serialization only
    }

    virtual ~Distributed_heightfield_seedings_updates();

    // Number of fragments to start the compute algorithm
    virtual mi::Size get_nb_of_fragments() const { return m_cluster_hosts.size(); }

    /// Updates the bounding box associated with this compute
    /// algorithm according to the results of the computation.
    /// \deprecated
    /// \param[in] bbox Current bounding box, will be modified if necessary.
    virtual void get_updated_bounding_box(mi::math::Bbox_struct<mi::Float32, 3>& bbox) const
    {
        // empty
    }

    ///
    void get_min_max_heights(mi::Float32& min_height, mi::Float32& max_height) const
    {
        min_height = m_min_value;
        max_height = m_max_value;
    }

    void set_patch_tag_lists(
        const std::vector<std::vector<std::pair<mi::math::Bbox_struct<mi::Uint32, 2>, mi::neuraylib::Tag> > >& patch_tag_lists);

    // -----------------------------------------------------------------------------------------
    /// Job scheduling is user-defined
    virtual mi::neuraylib::IFragmented_job::Scheduling_mode get_scheduling_mode() const
    {
        return mi::neuraylib::IFragmented_job::USER_DEFINED;
    }

    /// Override the normal scheduling of fragments to host by providing an array of host ids
    /// corresponding to the index for the fragment to assign to that host. For example:
    /// 1 2 2 2 3 would mean fragment 0 is assigned to host 1, fragment 1 to host 2 etc.
    ///
    /// \param slots        Pointer to an array of host ids in the order of which
    ///                     fragment should be assigned to which host.
    /// \param nr_slots     The number of host ids in the array that the pointer \p slots
    ///                     points to.
    virtual void assign_fragments_to_hosts(
        mi::Uint32* slots,
        mi::Size    nr_slots);

    /// Parallel compute on the local cluster host.
    /// \param dice_transaction         The transaction the job is executing in.
    /// \param index                    The index identifies for which fragment the execute function is called.
    /// \param count                    The number of fragments in which the job is split.
    virtual void execute_fragment(
        mi::neuraylib::IDice_transaction*               dice_transaction,
        mi::Size                                        index,
        mi::Size                                        count,
        const mi::neuraylib::IJob_execution_context*    context);

    /// Execute one fragment of the job on a remote host. This is executed on a different host than
    /// calling host and given to the receive_remote_result function.
    /// \param dice_transaction         The transaction the job is executing in.
    /// \param index                    The index identifies for which fragment the execute function is called.
    /// \param count                    The number of fragments in which the job is split.
    virtual void execute_fragment_remote(
		mi::neuraylib::ISerializer*                     serializer,
        mi::neuraylib::IDice_transaction*               dice_transaction,
        mi::Size                                        index,
        mi::Size                                        count,
        const mi::neuraylib::IJob_execution_context*    context);

    /// Receive a job result from a remote host and integrate it in the end result. This function
    /// the same effect as the call of execute_fragment on the original host.
    /// \param dice_transaction         The transaction the job is executing in.
    /// \param index                    The index identifies for which fragment the execute function is called.
    /// \param count                    The number of fragments in which the job is split.
    virtual void receive_remote_result(
        mi::neuraylib::IDeserializer*       deserializer,
        mi::neuraylib::IDice_transaction*   dice_transaction,
        mi::Size                            index,
        mi::Size                            count);

    // --------------------------------------------------------------------------------------------
    /// Serialize the class to the given serializer.
    /// \param serializer Write to this serializer.
    virtual void serialize(mi::neuraylib::ISerializer *serializer) const;

    /// Deserialize the class from the given deserializer.
    /// \param serializer Read from this deserializer.
    virtual void deserialize(mi::neuraylib::IDeserializer* deserializer);

private:
    /// Compute the heightfield data layout ....
    nv::index::IRegular_heightfield_data_locality* get_distribution_layout(
        mi::neuraylib::IDice_transaction*  dice_transaction);

    void update_heightfield_surface(
        mi::neuraylib::IDice_transaction*                                                         dice_transaction,
        mi::Uint32                                                                                host_id,
        const std::vector<std::pair<mi::math::Bbox_struct<mi::Uint32, 2>, mi::neuraylib::Tag> >&  patch_tag_pair_vector,
        mi::math::Vector_struct<mi::Float32, 2>&                                                  min_max_heights);

    void set_min_max_heights(
        mi::Float32 min_height,
        mi::Float32 max_height);

private:
    mi::neuraylib::Tag                                      m_session_tag;
    mi::neuraylib::Tag                                      m_heightfield_tag;

    // First implemented the user-defined scheduling.
    // Later used to access resp. edit the heightfield data on local or remote host.
    std::vector<mi::Uint32>                                                             m_cluster_hosts;
    std::vector<std::pair<mi::math::Bbox_struct<mi::Uint32, 2>, mi::neuraylib::Tag> >   m_patch_tag_pair_vector;

    // Locking when computing the returned height values
    mi::base::Lock                                                      m_compute_lock;
    mi::Float32                                                         m_max_value;
    mi::Float32                                                         m_min_value;
};


#endif // NVIDIA_INDEX_DISTRIBUTED_HORIZON_SEEDING_UPDATES_H
