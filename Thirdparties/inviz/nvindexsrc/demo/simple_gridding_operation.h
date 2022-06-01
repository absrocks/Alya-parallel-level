/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief gridding operation example

#ifndef NVIDIA_INDEX_SIMPLE_GRIDDING_OPERATION_H
#define NVIDIA_INDEX_SIMPLE_GRIDDING_OPERATION_H

#include <mi/base/types.h>
#include <mi/dice.h>
#include <mi/base/interface_declare.h>

#include <nv/index/idistributed_compute_algorithm.h>
#include <nv/index/idistributed_data_locality.h>

#include <mi/base/lock.h>

#include <vector>

//----------------------------------------------------------------------
/// Proof of concept computing algorithm that simply changes
/// the elevation values of the given heightfield.
class Simple_gridding_operation :
    public nv::index::Distributed_compute_algorithm<0xc0074e9a,0x0b0d,0x479b,0xb1,0x92,0x3d,0xf9,0x58,0xa4,0x58,0x2e>
{
public:
    /// constructor
    Simple_gridding_operation(
        const mi::neuraylib::Tag&                                    session_tag,
        const mi::neuraylib::Tag&                                    heightfield_tag,
        const std::vector<mi::math::Vector_struct<mi::Float32, 2> >& polygon);

    /// default constructor
    Simple_gridding_operation()
        :
        m_session_tag(mi::neuraylib::NULL_TAG),
        m_heightfield_tag(mi::neuraylib::NULL_TAG),
        m_polygon_of_interest(),
        m_mean_value(0.0f)
    {
        // for serialization only
    }

    /// destructor
    virtual ~Simple_gridding_operation();

public:
    // Implemented Distributed_compute_algorithm

    virtual mi::Size get_nb_of_fragments() const
    {
        // only one fragment, which starts sub jobs!
        return 1;
    }

    virtual void get_updated_bounding_box(mi::math::Bbox_struct<mi::Float32, 3>& bbox) const
    {
        if (bbox.min.z > m_mean_value)
        {
            bbox.min.z = m_mean_value;
        }
        if (bbox.max.z < m_mean_value)
        {
            bbox.max.z = m_mean_value;
        }
    }

    virtual void execute_fragment(
        mi::neuraylib::IDice_transaction*               dice_transaction,
        mi::Size                                        index,
        mi::Size                                        count,
        const mi::neuraylib::IJob_execution_context*    context);

private:
    mi::neuraylib::Tag                                      m_session_tag;
    mi::neuraylib::Tag                                      m_heightfield_tag;
    std::vector<mi::math::Vector_struct<mi::Float32, 2> >   m_polygon_of_interest;
    mi::Float32                                             m_mean_value;
};

// ----------------------------------------------------------------------------------------------
/// compute mean height value
class Gridder_elevation_mean_compute :
    public mi::neuraylib::Fragmented_job<0x556aa9e0,0x2094,0x4dc0,0x94,0xdf,0x7b,0xf8,0xc3,0x6c,0x6f,0x58>
{
public:
    /// constructor
    Gridder_elevation_mean_compute(
        const mi::neuraylib::Tag&                                    session_tag,
        const mi::neuraylib::Tag&                                    heightfield_tag,
        const std::vector<mi::math::Vector_struct<mi::Float32, 2> >& polygon,
        const std::vector<mi::Uint32>&                               cluster_hosts);

    /// default constructor
    Gridder_elevation_mean_compute()
        :
        m_session_tag(mi::neuraylib::NULL_TAG),
        m_heightfield_tag(mi::neuraylib::NULL_TAG),
        m_polygon_of_interest(),
        m_cluster_hosts(),
        m_result_nb_elevation_values(0),
        m_result_accumulated_elevation_values(),
        m_compute_lock()
    {
        // for serialization only
    };

    /// destructor
    virtual ~Gridder_elevation_mean_compute();

    /// Mean elevation value
    mi::Float32 get_mean_value() const;

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

    /// Serialize the class to the given serializer.
    /// \param serializer Write to this serializer.
    virtual void serialize(mi::neuraylib::ISerializer *serializer) const;

    /// Deserialize the class from the given deserializer.
    /// \param serializer Read from this deserializer.
    virtual void deserialize(mi::neuraylib::IDeserializer* deserializer);

private:
    /// Compute the mean based on cluster results.
    void compute_local_mean(
        mi::neuraylib::IDice_transaction*              dice_transaction,
        nv::index::IRegular_heightfield_data_locality* locality,
        mi::Uint32                                     host_id,
        mi::Uint32                                     patch_id,
        mi::Float32&                                   elevation,
        mi::Uint32&                                    nb_elevations);

    /// Compute the mean based on cluster results.
    void compute_mean(mi::Float32 elevation, mi::Uint32 nb_elevations);

private:
    mi::neuraylib::Tag                                      m_session_tag;
    mi::neuraylib::Tag                                      m_heightfield_tag;
    std::vector<mi::math::Vector_struct<mi::Float32, 2> >   m_polygon_of_interest;

    // To implemented the user-defined scheduling.
    // Later used to access resp. edit the heightfield data on local or remote host.
    std::vector<mi::Uint32>                                 m_cluster_hosts;

    mi::Uint32                                              m_result_nb_elevation_values;
    std::vector<mi::Float32>                                m_result_accumulated_elevation_values;

    // Locking when computing the returned results
    mi::base::Lock                                          m_compute_lock;
};

//----------------------------------------------------------------------

#endif // NVIDIA_INDEX_SIMPLE_GRIDDING_OPERATION_H
