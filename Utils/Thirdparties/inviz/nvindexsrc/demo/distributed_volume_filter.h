/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief distributed volume filter job: step 2 filter the data

#ifndef NVIDIA_INDEX_DISTRIBUTED_VOLUME_FILTER_H
#define NVIDIA_INDEX_DISTRIBUTED_VOLUME_FILTER_H

#include <mi/dice.h>
#include <mi/base/types.h>

#include <nv/index/idistributed_compute_algorithm.h>
#include <nv/index/idistributed_data_locality.h>

#include <vector>

/// Example of inset amplitude value filter job: step 2 filter the data.
class Distributed_volume_filter :
        public nv::index::Distributed_compute_algorithm<0x36b808ca,0x1d1b,0x4003,0x84,0x57,0xe7,0xe7,0xa5,0x6b,0x4e,0x3d>
{
public:
    /// constructor
    ///
    /// \param[in] session_tag   Tag of the session
    /// \param[in] volume_tag    Tag of the volume scene element
    /// \param[in] bordermap_element_tag Tag of the bordermap element
    /// \param[in] filter_type type of the filter
    /// \param[in] active_ijk_bbox Specifies the part of the volume
    /// where the operation should be applied
    /// \param[in] cluster_hosts host id list of the current cluster
    /// \param[in] parallel_count how many threads runs for processing.
    Distributed_volume_filter(
        const mi::neuraylib::Tag& session_tag,
        const mi::neuraylib::Tag& volume_tag,
        const mi::neuraylib::Tag& bordermap_element_tag,
        mi::Sint32 filter_type,
        mi::math::Bbox< mi::Uint32, 3 > const & active_ijk_bbox,
        std::vector< mi::Uint32 > const & cluster_hosts,
        mi::Uint32 parallel_count);

    /// constructor
    Distributed_volume_filter()
        :
        m_session_tag(mi::neuraylib::NULL_TAG),
        m_volume_tag(mi::neuraylib::NULL_TAG),
        m_bordermap_element_tag(mi::neuraylib::NULL_TAG),
        m_filter_type(0),
        m_active_ijk_bbox(),
        m_cluster_hosts(),
        m_parallel_count(0),
        m_fragment_to_host_id(),
        m_fragment_to_host_local_thread_id()
    {
        // empty. for serialization only
    };

    /// destructor
    virtual ~Distributed_volume_filter();

public:
    // implementation of IDistributed_compute_algorithm

    virtual mi::Size get_nb_of_fragments() const;

    virtual void get_updated_bounding_box(mi::math::Bbox_struct<mi::Float32, 3>& bbox) const
    {
        // empty
    }

public:
    virtual mi::neuraylib::IFragmented_job::Scheduling_mode get_scheduling_mode() const;

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

private:
    /// apply filter
    /// \param[in] dice_transaction dice database transaction
    /// \param[in] locality data locality
    /// \param[in] host_id  host id
    /// \param[in] fragment_index fragment index
    void apply_filter(
        mi::neuraylib::IDice_transaction*         dice_transaction,
        nv::index::IRegular_volume_data_locality* locality,
        mi::Uint32                                host_id,
        mi::Uint32                                fragment_index);

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
    /// User-defined scheduling. available cluster hostid list
    std::vector< mi::Uint32 > m_cluster_hosts;
    /// using number of m_parallel_count threads for filtering.
    mi::Uint32 m_parallel_count;
    /// fragment host map.   map: fragment_index -> host id
    std::vector< mi::Sint32 > m_fragment_to_host_id;
    /// fragment thread map: map: fragment_index -> host local thread id
    std::vector< mi::Sint32 > m_fragment_to_host_local_thread_id;

};


#endif // #ifndef NVIDIA_INDEX_DISTRIBUTED_VOLUME_FILTER_H
