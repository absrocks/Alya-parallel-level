/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief distributed amplitude value map generation job at the heightfield intersection.

#ifndef NVIDIA_INDEX_TEMPLATE_DISTRIBUTED_AMPLITUDE_MAP_JOB_H
#define NVIDIA_INDEX_TEMPLATE_DISTRIBUTED_AMPLITUDE_MAP_JOB_H

#include <mi/base/types.h>
#include <mi/dice.h>
#include <mi/base/interface_declare.h>

#include <nv/index/idistributed_compute_algorithm.h>
#include <nv/index/idistributed_data_locality.h>

#include <vector>

#include "heightfield_data_retrieval_appjob.h"
#include "image_tile.h"

#include "common/encode_voxel.h"
#include "common/encode_height.h"

/// Distributed amplitude value map at the heightfield intersection points.
///
/// Create an image of volume amplitude value that is on the heightfield data.
class Distributed_amplitude_map_job :
    public nv::index::Distributed_compute_algorithm<0xec206d63,0x7fe6, 0x47ae,0x9a,0x72,0xc9,0x05,0x78,0x36,0xc6,0x95>
{
public:
    /// constructor
    ///
    /// Pass the parameters for this job.
    ///
    /// \param[in] session_tag     session tag
    /// \param[in] heightfield_tag heightfield tag
    /// \param[in] volume_tag      volume tag
    /// \param[in] cluster_hosts   host assignment
    Distributed_amplitude_map_job(
        const mi::neuraylib::Tag&      session_tag,
        const mi::neuraylib::Tag&      heightfield_tag,
        const mi::neuraylib::Tag&      volume_tag,
        std::vector<mi::Uint32> const & cluster_hosts);

    /// default constructor. Only for serialization.
    Distributed_amplitude_map_job()
        :
        m_session_tag(mi::neuraylib::NULL_TAG),
        m_heightfield_tag(mi::neuraylib::NULL_TAG),
        m_volume_tag(mi::neuraylib::NULL_TAG),
        m_cluster_hosts(),
        m_whole_amplitude_map(),
        m_partial_amplitude_map_vec(),
        m_data_lock()
    {
        // empty
    }

    /// destructor
    virtual ~Distributed_amplitude_map_job();

    /// Number of fragments to start the compute algorithm
    virtual mi::Size get_nb_of_fragments() const { return m_cluster_hosts.size(); }

    /// Updates the heightfield bounding box according to
    /// the results of the computation.
    /// \deprecated
    /// \param[in,out] bbox Current heightfield bbox, will be modified if necessary.
    virtual void get_updated_bounding_box(mi::math::Bbox_struct<mi::Float32, 3>& bbox) const;

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
        mi::neuraylib::IJob_execution_context const *   context);

    /// Execute one fragment of the job on a remote host. This is executed on a different host than
    /// calling host and given to the receive_remote_result function.
    /// \param dice_transaction         The transaction the job is executing in.
    /// \param index                    The index identifies for which fragment the execute function is called.
    /// \param count                    The number of fragments in which the job is split.
    virtual void execute_fragment_remote(
        mi::neuraylib::ISerializer*        serializer,
        mi::neuraylib::IDice_transaction*             dice_transaction,
        mi::Size                                      index,
        mi::Size                                      count,
        mi::neuraylib::IJob_execution_context const * context);

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
    /// \param[in] deserializer Read from this deserializer.
    virtual void deserialize(mi::neuraylib::IDeserializer* deserializer);

public:
    /// save the result amplitude map to file
    ///
    /// \param[in] fname amplitude map filename to be saved
    /// \return true when succeeded
    bool save_amplitude_map(std::string const & fname);

private:
    /// get global xyz ROI
    /// \param[in]  dice_transaction dice transaction
    /// \return global xyz ROI
    mi::math::Bbox< mi::Float32, 3 > get_global_xyz_roi(
        mi::base::Handle< mi::neuraylib::IDice_transaction > & dice_transaction);

    /// get heightfield patches in the global ROI
    ///
    /// \param[in]  host_id host id of this node
    /// \param[out] patch_ij_bbox_vec  (output) patch ij bounding box vector
    /// \param[in]  dice_transaction dice transaction
    /// \return true when success
    bool get_heightfield_bbox_in_ROI(
        mi::Uint32 host_id,
        std::vector< mi::math::Bbox<mi::Sint32, 2> > & patch_ij_bbox_vec,
        mi::base::Handle< mi::neuraylib::IDice_transaction > & dice_transaction);

    /// get heightfield ijk bbox
    ///
    /// \param[in] heightfield_data_access heightfield data access
    /// \return heightfield ijk bounding box
    mi::math::Bbox<mi::Uint32, 3> get_heightfield_ijk_bbox(
        mi::base::Handle<nv::index::IRegular_heightfield_data_access>& heightfield_data_access);

    /// get heightfield patch height range as integer
    ///
    /// \param[in] patch_bbox_s64_2 patch bbox
    /// \param[in] p_height_values  height value array
    /// \return height range (exclude -1 values).
    mi::math::Vector< mi::Sint64, 2 > get_heightfield_patch_height_range(
        mi::math::Bbox< mi::Sint64, 2 > const & patch_bbox_s64_2,
        mi::Float32* p_height_values);

    /// allocate the result amplitude map
    /// \param[in] dice_transaction dice transaction
    void allocate_result_map(
        mi::base::Handle< mi::neuraylib::IDice_transaction > & dice_transaction);


    /// generate amplitude map
    /// \param[in] dice_transaction dice transaction
    /// \param[in] host_id          host id
    /// \param[in] is_result_host   true when this is the result collecting host
    bool generate_amplitude_map(
        mi::base::Handle< mi::neuraylib::IDice_transaction > & dice_transaction,
        mi::Uint32                                             host_id,
        bool                                                   is_result_host);

private:
    /// session tag
    mi::neuraylib::Tag                           m_session_tag;
    /// heightfield tag to manipulate
    mi::neuraylib::Tag                           m_heightfield_tag;
    /// volume tag to manipulate
    mi::neuraylib::Tag                           m_volume_tag;

    // First implemented the user-defined scheduling.
    // Later used to access resp. edit the heightfield data on local or remote host.
    std::vector<mi::Uint32>                      m_cluster_hosts;

private:
    /// the whole amplitude map (only the host allocates the buffer)
    Image_tile m_whole_amplitude_map;

    /// the local computed partial amplitude map vector
    std::vector< Image_tile > m_partial_amplitude_map_vec;

private:
    // non serialized objects

    /// a lock
    mi::base::Lock                               m_data_lock;
};


#endif // NVIDIA_INDEX_TEMPLATE_DISTRIBUTED_AMPLITUDE_MAP_JOB_H
