/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief distributed volume bordermap generation job

#ifndef NVIDIA_INDEX_DISTRIBUTED_VOLUME_BORDERMAP_GENERATION_H
#define NVIDIA_INDEX_DISTRIBUTED_VOLUME_BORDERMAP_GENERATION_H

#include <mi/dice.h>
#include <mi/base/types.h>

#include <nv/index/idistributed_compute_algorithm.h>
#include <nv/index/idistributed_data_locality.h>

#include <vector>

/// Example of inset amplitude value filter job: step 1 create volume bordermap.
class Distributed_volume_bordermap_generation :
    public nv::index::Distributed_compute_algorithm<0xb623a43a,0xba6f,0x4c6f,0x8d,0x89,0x03,0x22,0x53,0xe1,0xdc,0xe8>
{
public:
    /// constructor
    ///
    /// \param[in] session_tag     session tag
    /// \param[in] volume_tag      volume scene element tag
    /// \param[in] filter_type     type name of the filter
    /// \param[in] active_ijk_bbox Specifies the part of the volume
    /// where the operation should be applied
    /// \param[in] cluster_hosts
    Distributed_volume_bordermap_generation(
        const mi::neuraylib::Tag& session_tag,
        const mi::neuraylib::Tag& volume_tag,
        mi::Sint32 filter_type,
        mi::math::Bbox< mi::Uint32, 3 > const & active_ijk_bbox,
        std::vector< mi::Uint32 > const &       cluster_hosts);

    /// constructor
    Distributed_volume_bordermap_generation()
        :
        m_session_tag(mi::neuraylib::NULL_TAG),
        m_volume_tag(mi::neuraylib::NULL_TAG),
        m_bordermap_element_tag(mi::neuraylib::NULL_TAG),
        m_filter_type(0),
        m_active_ijk_bbox(),
        m_cluster_hosts(),
        m_result_job_bordermap_tag_vec()
    {
        // empty. for serialization only
    };

    /// destructor
    virtual ~Distributed_volume_bordermap_generation();

    /// get the job result
    /// \return created bordermap tag vector
    std::vector< mi::neuraylib::Tag > const & get_job_result_bordermap_tag_vec() const;

public:
    // implementation of IDistributed_compute_algorithm

    virtual mi::Size get_nb_of_fragments() const
    {
        return m_cluster_hosts.size();
    }

    virtual void get_updated_bounding_box(mi::math::Bbox_struct<mi::Float32, 3>& bbox) const
    {
        // empty
    }

public:

    /// Job scheduling is user-defined
    virtual mi::neuraylib::IFragmented_job::Scheduling_mode get_scheduling_mode() const
    {
        return mi::neuraylib::IFragmented_job::USER_DEFINED;
    }

    /// Override the normal scheduling of fragments to host by providing an array of host ids
    /// corresponding to the index for the fragment to assign to that host. For example:
    /// 1 2 2 2 3 would mean fragment 0 is assigned to host 1, fragment 1 to host 2 etc.
    ///
    /// \param[in] slots        Pointer to an array of host ids in the order of which
    ///                     fragment should be assigned to which host.
    /// \param[in] nr_slots     The number of host ids in the array that the pointer \p slots
    ///                     points to.
    virtual void assign_fragments_to_hosts(
        mi::Uint32* slots,
        mi::Size    nr_slots);

    /// Parallel compute on the local cluster host.
    /// \param[in] dice_transaction  The transaction the job is executing in.
    /// \param[in] index             The index identifies for which fragment the
    /// execute function is called.
    /// \param[in] count             The number of fragments in which the job is split.
    /// \param[in] context
    virtual void execute_fragment(
        mi::neuraylib::IDice_transaction*               dice_transaction,
        mi::Size                                        index,
        mi::Size                                        count,
        const mi::neuraylib::IJob_execution_context*    context);

    /// Execute one fragment of the job on a remote host. This is executed on a different host than
    /// calling host and given to the receive_remote_result function.
    ///
    /// \param[in] serializer               serializer
    /// \param[in] dice_transaction         The transaction the job is executing in.
    /// \param[in] index The index identifies for which fragment the
    /// execute function is called.
    /// \param[in] count                    The number of fragments in which the job is split.
    /// \param[in] context
    virtual void execute_fragment_remote(
        mi::neuraylib::ISerializer*                     serializer,
        mi::neuraylib::IDice_transaction*               dice_transaction,
        mi::Size                                        index,
        mi::Size                                        count,
        const mi::neuraylib::IJob_execution_context*    context);

    /// get the result from remote. On the local host, the result from
    /// the fragment is read from the deserializer and stored in the
    /// result array.
    virtual void receive_remote_result(
        mi::neuraylib::IDeserializer*       deserializer,
        mi::neuraylib::IDice_transaction*   transaction,
        mi::Size                            index,
        mi::Size                            count);


    /// Return the unique class id of this class
    virtual mi::base::Uuid get_class_id() const
    {
        return IID();
    }

    /// Serialize the class to the given serializer.
    virtual void serialize(mi::neuraylib::ISerializer *serializer) const;

    /// Deserialize the class from the given deserializer.
    virtual void deserialize(mi::neuraylib::IDeserializer* deserializer);

private:
    /// create border
    /// \param[in] dice_transaction dice database transaction
    /// \param[in] locality data locality
    /// \param[in] host_id  host id
    /// \param[in] job_local_bordermap_tag job local bordermap tag
    void create_border(
        mi::neuraylib::IDice_transaction*         dice_transaction,
        nv::index::IRegular_volume_data_locality* locality,
        mi::Uint32                                host_id,
        mi::neuraylib::Tag                        job_local_bordermap_tag);

private:
    /// session
    mi::neuraylib::Tag m_session_tag;
    /// scene element (volume data tag to be processed)
    mi::neuraylib::Tag m_volume_tag;
    /// bordermap database element tag
    mi::neuraylib::Tag m_bordermap_element_tag;
    /// Volume_data_filter's filter type
    mi::Sint32 m_filter_type;
    /// computing IJK region of interesr bounding box
    mi::math::Bbox<mi::Uint32, 3> m_active_ijk_bbox;
    /// User-defined scheduling.
    std::vector<mi::Uint32> m_cluster_hosts;

    /// results of job border map vector (not serialized)
    std::vector< mi::neuraylib::Tag > m_result_job_bordermap_tag_vec;
};


#endif // NVIDIA_INDEX_DISTRIBUTED_VOLUME_BORDERMAP_GENERATION_H
