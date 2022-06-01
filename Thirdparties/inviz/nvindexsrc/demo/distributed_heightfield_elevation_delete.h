/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief heightfield elevation value delete

#ifndef NVIDIA_INDEX_DISTRIBUTED_HEIGHTFIELD_ELEVATION_DELETE_H
#define NVIDIA_INDEX_DISTRIBUTED_HEIGHTFIELD_ELEVATION_DELETE_H

#include <mi/base/types.h>
#include <mi/dice.h>
#include <mi/base/interface_declare.h>

#include <nv/index/idistributed_compute_algorithm.h>
#include <nv/index/idistributed_data_locality.h>

#include <vector>


/// Proof of concept computing algorithm that simply changes
/// the elevation values of the given heightfield.
class Distributed_heightfield_elevation_delete :
    public nv::index::Distributed_compute_algorithm<0x7b864170,0x2db9,0x4b1a,0xbd,0xbb,0x23,0x13,0x5a,0xc9,0x35,0x39>
{
public:
    Distributed_heightfield_elevation_delete(
        const mi::neuraylib::Tag&                                       session_tag,
        const mi::neuraylib::Tag&                                       heightfield_tag,
        const std::vector<mi::Uint32>&                                  cluster_hosts);
    Distributed_heightfield_elevation_delete(
        const mi::neuraylib::Tag&                                       session_tag,
        const mi::neuraylib::Tag&                                       heightfield_tag,
        const std::vector<mi::math::Vector_struct<mi::Float32, 2> >&    polygon,
        const std::vector<mi::Uint32>&                                  cluster_hosts);

    Distributed_heightfield_elevation_delete()
        :
        m_session_tag(mi::neuraylib::NULL_TAG),
        m_heightfield_tag(mi::neuraylib::NULL_TAG),
        m_cluster_hosts(),
        m_polygon_of_interest()
    {
        // for serialization only
    }

    virtual ~Distributed_heightfield_elevation_delete();

    // Number of fragments to start the compute algorithm
    virtual mi::Size get_nb_of_fragments() const { return m_cluster_hosts.size(); }

    /// Updates the bounding box associated with this compute
    /// algorithm according to the results of the computation.
    /// \deprecated
    /// \param[in,out] bbox Current bbox, will be modified if necessary.
    virtual void get_updated_bounding_box(mi::math::Bbox_struct<mi::Float32, 3>& bbox) const
    {
        // empty
    }


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
		mi::neuraylib::ISerializer*  serializer,
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
    void edit_heightfield_data(
        mi::neuraylib::IDice_transaction*               dice_transaction,
        nv::index::IRegular_heightfield_data_locality*  locality,
        mi::Uint32                                      host_id);

private:
    mi::neuraylib::Tag                                      m_session_tag;
    mi::neuraylib::Tag                                      m_heightfield_tag;

    // First implemented the user-defined scheduling.
    // Later used to access resp. edit the heightfield data on local or remote host.
    std::vector<mi::Uint32>                                 m_cluster_hosts;
    std::vector<mi::math::Vector_struct<mi::Float32, 2> >   m_polygon_of_interest;
};


#endif // NVIDIA_INDEX_DISTRIBUTED_HEIGHTFIELD_ELEVATION_DELETE_H
