/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "distributed_amplitude_map_job.h"

#include <cassert>

#include <nv/index/isession.h>
#include <nv/index/iregular_heightfield.h>
#include <nv/index/iregular_heightfield_compute_task.h>

#include "scene_utility.h"
#include "utilities.h"
#include "volume_data_retrieval_appjob.h"

#include "common/forwarding_logger.h"
#include "common/type_conversion_utility.h"

//----------------------------------------------------------------------
// FIXME: put this in scene_utility.
nv::index::IRegular_heightfield_data_locality * new_heightfield_bbox_distribution_layout(
    const mi::neuraylib::Tag&            session_tag,
    const mi::neuraylib::Tag&            heightfield_tag,
    mi::math::Bbox_struct<mi::Uint32, 2> const & ij_query_bounds,
    bool                                  is_edit,
    mi::neuraylib::IDice_transaction *    dice_transaction)
{
    assert(session_tag.is_valid());
    assert(heightfield_tag.is_valid());
    assert(dice_transaction != 0);

    // Access session
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<const nv::index::ISession>(session_tag));
    assert(session.is_valid_interface());
    const mi::neuraylib::Tag& dist_layout_tag = session->get_distribution_layout();
    assert(dist_layout_tag.is_valid());

    // Access the distribution scheme
    mi::base::Handle<const nv::index::IData_distribution> distribution_layout(
        dice_transaction->access<const nv::index::IData_distribution>(dist_layout_tag));
    assert(distribution_layout.is_valid_interface());

    // Heightfield distribution layout: this is new-ed. Use handle at the
    // caller side. (We can not use handle here since it is
    // immidiately out of scope.)
    nv::index::IRegular_heightfield_data_locality* p_data_locality = 0;
    if(is_edit){
        p_data_locality = distribution_layout->retrieve_data_locality_for_editing(
            heightfield_tag, ij_query_bounds, dice_transaction);
    }
    else{
        p_data_locality = distribution_layout->retrieve_data_locality(
            heightfield_tag, ij_query_bounds, dice_transaction);
    }

    return p_data_locality;
}


//----------------------------------------------------------------------
/// extend heightfield bbox in Z direction
/// \param[in,out] heightfield_ijk_bbox extend this bbox
static void extend_heightfield_bbox_z(mi::math::Bbox<mi::Uint32, 3> & heightfield_ijk_bbox)
{
    if(heightfield_ijk_bbox.min.z > 0){
        // if possible try min to zero first
        --heightfield_ijk_bbox.min.z;
    }
    // also add two more to max
    ++heightfield_ijk_bbox.max.z;
}

//======================================================================
Distributed_amplitude_map_job::Distributed_amplitude_map_job(
    const mi::neuraylib::Tag&      session_tag,
    const mi::neuraylib::Tag&      heightfield_tag,
    const mi::neuraylib::Tag&      volume_tag,
    std::vector<mi::Uint32> const & cluster_hosts)
    :
    m_session_tag(session_tag),
    m_heightfield_tag(heightfield_tag),
    m_volume_tag(volume_tag),
    m_cluster_hosts(cluster_hosts)
{
    // empty
}

//----------------------------------------------------------------------
Distributed_amplitude_map_job::~Distributed_amplitude_map_job()
{
    // empty
}

//----------------------------------------------------------------------
void Distributed_amplitude_map_job::get_updated_bounding_box(
    mi::math::Bbox_struct<mi::Float32, 3>& bbox) const
{
    // empty
}

//----------------------------------------------------------------------
void Distributed_amplitude_map_job::assign_fragments_to_hosts(
    mi::Uint32* slots,
    mi::Size    nr_slots)
{
    if(nr_slots != m_cluster_hosts.size())
    {
        ERROR_LOG << "The number of host ids present ("
                  << m_cluster_hosts.size()
                  << ") and requested ("
                  << nr_slots
                  << ") doesn't match the requested count.";
        assert(m_cluster_hosts.size()==nr_slots);
    }
    else
    {
        for(mi::Size i=0; i<nr_slots; ++i)
        {
            slots[i] = m_cluster_hosts[i];
        }
    }
}

//----------------------------------------------------------------------
void Distributed_amplitude_map_job::execute_fragment(
    mi::neuraylib::IDice_transaction*              dice_transaction,
    mi::Size                                       index,
    mi::Size                                       count,
    mi::neuraylib::IJob_execution_context const *  context)
{
    assert(index < m_cluster_hosts.size());
    mi::Uint32 const cluster_host = m_cluster_hosts[index];

    mi::base::Handle<mi::neuraylib::IDice_transaction>
        dice_trans_hnd(dice_transaction, mi::base::DUP_INTERFACE);

    bool const is_result_host = true;
    this->generate_amplitude_map(dice_trans_hnd, cluster_host, is_result_host);

    {
        // scope for the lock
        mi::base::Lock::Block block(&m_data_lock);
        assert(m_whole_amplitude_map.is_valid_buffer());

        for(std::vector< Image_tile >::const_iterator mi = m_partial_amplitude_map_vec.begin();
            mi != m_partial_amplitude_map_vec.end();
            ++mi)
        {
            m_whole_amplitude_map.copy_partial_image_tile(*mi,
                                                          mi->get_bounding_box(),
                                                          mi->get_bounding_box().min);
            DEBUG_LOG << "execute_fragment: partial map: " << mi->get_bounding_box();
        }
    }
}

//----------------------------------------------------------------------
void Distributed_amplitude_map_job::execute_fragment_remote(
    mi::neuraylib::ISerializer*                   serializer,
    mi::neuraylib::IDice_transaction*             dice_transaction,
    mi::Size                                      index,
    mi::Size                                      count,
    mi::neuraylib::IJob_execution_context const * context)
{
    assert(index < m_cluster_hosts.size());
    mi::Uint32 const cluster_host = m_cluster_hosts[index];

    mi::base::Handle<mi::neuraylib::IDice_transaction>
        dice_trans_hnd(dice_transaction, mi::base::DUP_INTERFACE);

    bool const is_result_host = false;
    this->generate_amplitude_map(dice_trans_hnd, cluster_host, is_result_host);


    // serialize the result.
    mi::Size const map_count = m_partial_amplitude_map_vec.size();
    serializer->write(&map_count, 1);
    for(std::vector< Image_tile >::const_iterator mi = m_partial_amplitude_map_vec.begin();
        mi != m_partial_amplitude_map_vec.end();
        ++mi)
    {
        mi->serialize(serializer);
        DEBUG_LOG << "execute_fragment_remote[" << index << "]: " << mi->get_bounding_box();
    }
}

//----------------------------------------------------------------------
void Distributed_amplitude_map_job::receive_remote_result(
    mi::neuraylib::IDeserializer*                   deserializer,
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Size                                        index,
    mi::Size                                        count)
{

    std::vector< Image_tile > receive_tile_vec;
    mi::Size map_count = 0;
    deserializer->read(&map_count, 1);

    for(mi::Size i = 0; i < map_count; ++i){
        Image_tile receive_tile;
        receive_tile.deserialize(deserializer);
        receive_tile_vec.push_back(receive_tile);
    }
    {
        // scope for the lock
        mi::base::Lock::Block block(&m_data_lock);

        for(std::vector< Image_tile >::const_iterator mi = receive_tile_vec.begin();
            mi != receive_tile_vec.end();
            ++mi)
        {
            assert(mi->is_valid_buffer());
            m_whole_amplitude_map.copy_partial_image_tile(*mi,
                                                          mi->get_bounding_box(),
                                                          mi->get_bounding_box().min);
            DEBUG_LOG << "receive_remote_result copy[" << index << "]: " << mi->get_bounding_box();
        }
    }
}

//----------------------------------------------------------------------
void Distributed_amplitude_map_job::serialize(
    mi::neuraylib::ISerializer* serializer) const
{
    serializer->write(&m_session_tag.id,   1);
    serializer->write(&m_heightfield_tag.id, 1);
    serializer->write(&m_volume_tag.id, 1);

    mi::Uint32 nb_elements = mi::Uint32(m_cluster_hosts.size());
    serializer->write(&nb_elements, 1);
    for(mi::Uint32 i=0; i<nb_elements; ++i){
        serializer->write(&m_cluster_hosts[i], 1);
    }
}

//----------------------------------------------------------------------
void Distributed_amplitude_map_job::deserialize(
    mi::neuraylib::IDeserializer* deserializer)
{
    deserializer->read(&m_session_tag.id,   1);
    deserializer->read(&m_heightfield_tag.id, 1);
    deserializer->read(&m_volume_tag.id, 1);

    mi::Uint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_cluster_hosts.resize(nb_elements);
    for(mi::Uint32 i=0; i<nb_elements; ++i){
        deserializer->read(&m_cluster_hosts[i], 1);
    }
}

//----------------------------------------------------------------------
void Distributed_amplitude_map_job::allocate_result_map(
    mi::base::Handle< mi::neuraylib::IDice_transaction > & dice_transaction)
{
    if(m_whole_amplitude_map.is_valid_buffer()){
        ERROR_LOG << "duplicate allocation? bbox: " << m_whole_amplitude_map.get_bounding_box();
        return;
    }

    mi::math::Bbox< mi::Sint64, 3 > const heightfield_ijk_bbox = get_heightfield_local_IJK_ROI_bbox(
        dice_transaction.get(), m_heightfield_tag);
    mi::math::Bbox< mi::Sint64, 2 > const patch_ij_bbox =
        convert_bbox_sint64_3_to_sint64_2(heightfield_ijk_bbox);
    if(!patch_ij_bbox.is_plane()){
        ERROR_LOG << "Distributed_amplitude_map_job::allocate_result_map: illegal heightfield bbox.";
        return;
    }
    m_whole_amplitude_map.resize_buffer(patch_ij_bbox);
    INFO_LOG << "amplitude map size: " << patch_ij_bbox;
}

//----------------------------------------------------------------------
bool Distributed_amplitude_map_job::generate_amplitude_map(
    mi::base::Handle< mi::neuraylib::IDice_transaction > & dice_transaction,
    mi::Uint32                                             host_id,
    bool                                                   is_result_host)
{
    if(is_result_host){
        // If this is the result collecting host, allocate the final result map.
        this->allocate_result_map(dice_transaction);
    }

    // get local heightfield patch bbox
    std::vector< mi::math::Bbox<mi::Sint32, 2> > patch_ij_bbox_vec;
    if(!this->get_heightfield_bbox_in_ROI(host_id, patch_ij_bbox_vec, dice_transaction)){
        ERROR_LOG << "can not get patch value ranges.";
        return false;
    }

    // retrieve heightfield data
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<const nv::index::ISession>(
            m_session_tag));
    assert(session.is_valid_interface());
    mi::base::Handle<const nv::index::IDistributed_data_access_factory>
        subsurface_data_access_factory(
            dice_transaction->access<const nv::index::IDistributed_data_access_factory>(
                session->get_data_access_factory()));
    assert(subsurface_data_access_factory.is_valid_interface());

    Heightfield_data_retrieval_appjob heightfield_retrieval_job(
        m_heightfield_tag,
        patch_ij_bbox_vec,
        subsurface_data_access_factory);
    dice_transaction->execute_fragmented(&heightfield_retrieval_job, 1);


    // process patch pidx
    size_t const patch_count = patch_ij_bbox_vec.size();
    // DEBUG_LOG << "retrieved heightfield data. patch_count = " << patch_count;

    for(size_t pidx = 0; pidx < patch_count; ++pidx)
    {
        // get heightfield bbox
        mi::base::Handle< nv::index::IRegular_heightfield_data_access > heightfield_data_access =
            heightfield_retrieval_job.get_heightfield_data(pidx);

        mi::math::Bbox<mi::Uint32, 3> heightfield_ijk_bbox =
            this->get_heightfield_ijk_bbox(heightfield_data_access);
        // DEBUG_LOG << "heightfield_ijk_bbox: " << heightfield_ijk_bbox;
        extend_heightfield_bbox_z(heightfield_ijk_bbox);
        // DEBUG_LOG << "extended heightfield_ijk_bbox: " << heightfield_ijk_bbox;
        assert(heightfield_ijk_bbox.is_volume());

        // get the volume data.
        //
        // Here, we assume we have enough memory for processing a patch.
        // A patch may need one column of subregion at most.
        //
        // Note: For simplicity of example, this implementation
        // processes patch by patch. But this may trigger volume data
        // transfer. In real implementation, patch should be
        // transferred for efficiency since patch transfer is less
        // expensive in general compare to volume transfer.
        mi::base::Handle<nv::index::IRegular_volume_data_access> volume_data_access(
            subsurface_data_access_factory->create_regular_volume_data_access(m_volume_tag));
        assert(volume_data_access.is_valid_interface());
        volume_data_access->access(heightfield_ijk_bbox, dice_transaction.get());

        const mi::base::Handle<const nv::index::IRegular_volume_data> volume_data(
            volume_data_access->get_volume_data());
        const mi::base::Handle<const nv::index::IRegular_volume_data_uint8> volume_data_uint8(
            volume_data->get_interface<const nv::index::IRegular_volume_data_uint8>());

        if (!volume_data_uint8) {
            ERROR_LOG << "Amplitude-map generation failed: data access on non-uint8 volume type.";
            return false;
        }

        mi::math::Bbox<mi::Sint64, 3> const heightfield_ijk_bbox_64_3(
            nv::index_common::convert_bbox_type<mi::Sint64, mi::Uint32, 3>(heightfield_ijk_bbox));
        mi::math::Bbox<mi::Uint32, 3> const effective_bbox = volume_data_access->get_bounding_box();
        mi::math::Bbox<mi::Sint64, 3> const effective_bbox_64_3(
            nv::index_common::convert_bbox_type<mi::Sint64, mi::Uint32, 3>(effective_bbox));
        // DEBUG_LOG << "effective volume bbox: " << effective_bbox;
        if(effective_bbox.empty()){
            ERROR_LOG << "effective_bbox has no volume.";
        }

        // create a sub image
        nv::index_common::VoxelType  const * const p_volume_dat  = volume_data_uint8->get_voxel_data();
        HeightType const * const p_height_dat = heightfield_data_access->get_elevation_values();
        mi::math::Bbox< mi::Sint64, 2 > const heightfield_raw_bbox_64_2 =
            convert_bbox_sint64_3_to_sint64_2(heightfield_ijk_bbox_64_3);


        mi::math::Color const hole_col(1.0, 0.0, 0.0, 1.0);
        mi::math::Color amp_col (0.0, 0.0, 0.0, 1.0);

        Image_tile amp_map_tile;
        amp_map_tile.resize_buffer(heightfield_raw_bbox_64_2);

        mi::math::Vector< mi::Sint64, 2 > ij(heightfield_ijk_bbox_64_3.min.x, heightfield_ijk_bbox_64_3.min.y);
        for(ij[1] = heightfield_ijk_bbox_64_3.min[1]; ij[1] < heightfield_ijk_bbox_64_3.max[1]; ++ij[1])
        {
            for(ij[0] = heightfield_ijk_bbox_64_3.min[0]; ij[0] < heightfield_ijk_bbox_64_3.max[0]; ++ij[0])
            {
                mi::Sint64 const hor_idx = nv::index_common::get_heightfield_index(ij, heightfield_raw_bbox_64_2);
                mi::Sint64 const hval = static_cast< mi::Sint64 >(p_height_dat[hor_idx]);
                mi::math::Vector< mi::Sint64, 3 > ijk(ij[0], ij[1], hval);
                if(hval < 0)
                {
                    amp_map_tile.put_color(ij, hole_col);
                }
                else
                {
                    mi::Sint64 const vol_idx = nv::index_common::get_volume_index(ijk, effective_bbox_64_3);
                    nv::index_common::VoxelType const amplitude_val = p_volume_dat[vol_idx];
                    mi::Float32 const mono_col = static_cast< mi::Float32 >(amplitude_val) / 255.0f;
                    amp_col.r = amp_col.g = amp_col.b = mono_col;
                    amp_map_tile.put_color(ij, amp_col);
                }
            }
        }
        m_partial_amplitude_map_vec.push_back(amp_map_tile);
        INFO_LOG << "processed: patch " << pidx << "/" << patch_count;
    }

    return true;
}

//----------------------------------------------------------------------
bool Distributed_amplitude_map_job::save_amplitude_map(std::string const & fname)
{
    if(!m_whole_amplitude_map.is_valid_buffer()){
        ERROR_LOG << "invalid result amplitude map buffer.";
        return false;
    }
    bool const is_success = m_whole_amplitude_map.save_buffer(fname);

    return is_success;
}

//----------------------------------------------------------------------
mi::math::Bbox< mi::Float32, 3 > Distributed_amplitude_map_job::get_global_xyz_roi(
        mi::base::Handle< mi::neuraylib::IDice_transaction > & dice_transaction)
{
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<const nv::index::ISession>(
            m_session_tag));
    assert(session.is_valid_interface());

    mi::math::Bbox<mi::Float32, 3> const global_xyz_roi =
        get_XYZ_global_region_of_interest_bbox(session.get(), dice_transaction.get());

    return global_xyz_roi;
}

//----------------------------------------------------------------------
bool Distributed_amplitude_map_job::get_heightfield_bbox_in_ROI(
    mi::Uint32 host_id,
    std::vector< mi::math::Bbox<mi::Sint32, 2> > & patch_ij_bbox_vec,
    mi::base::Handle< mi::neuraylib::IDice_transaction > & dice_transaction)
{
    // get heightfield data locality
    const mi::math::Bbox<mi::Float32, 3> global_xyz_roi = this->get_global_xyz_roi(dice_transaction);
    // DEBUG_LOG << "global_xyz_roi = " << global_xyz_roi;

    for (mi::Sint32 i = 0; i < 2; ++i)
    {
        if (floor(global_xyz_roi.min[i]) < 0.0)
        {
            ERROR_LOG << "This example did not handle minus min ROI.";
            return false;
        }

        if (ceil(global_xyz_roi.max[i]) < 0.0)
        {
            ERROR_LOG << "This example did not handle minus max ROI.";
            return false;
        }
    }

    mi::math::Bbox_struct<mi::Uint32, 2> ij_heightfield_roi;
    ij_heightfield_roi.min.x = static_cast< mi::Uint32 >(floor(global_xyz_roi.min.x));
    ij_heightfield_roi.min.y = static_cast< mi::Uint32 >(floor(global_xyz_roi.min.y));
    ij_heightfield_roi.max.x = static_cast< mi::Uint32 >(ceil (global_xyz_roi.max.x));
    ij_heightfield_roi.max.y = static_cast< mi::Uint32 >(ceil (global_xyz_roi.max.y));

    bool const is_edit = false;
    mi::base::Handle<nv::index::IRegular_heightfield_data_locality> heightfield_data_locality(
        new_heightfield_bbox_distribution_layout(m_session_tag,
                                                 m_heightfield_tag,
                                                 ij_heightfield_roi,
                                                 is_edit,
                                                 dice_transaction.get()));

    // get this node local heightfield bboxes
    // mi::Uint32 const h_node_count     = heightfield_data_locality->get_nb_cluster_nodes();
    const mi::Uint32 patch_bbox_count = 
        static_cast<mi::Uint32>(heightfield_data_locality->get_nb_bounding_box(host_id));
    // DEBUG_LOG << "hostid: " << host_id << "/" << h_node_count << ", bbox_count = "
    //           << patch_bbox_count;

    patch_ij_bbox_vec.clear();
    for(mi::Uint32 i = 0; i < patch_bbox_count; ++i)
    {
        const mi::math::Bbox<mi::Sint32, 3> patch_3D_bbox = heightfield_data_locality->get_bounding_box(host_id, i);
        mi::math::Bbox<mi::Sint32, 2> patch_bbox;
        patch_bbox.min.x = patch_3D_bbox.min.x;
        patch_bbox.min.y = patch_3D_bbox.min.y;
        patch_bbox.max.x = patch_3D_bbox.max.x;
        patch_bbox.max.y = patch_3D_bbox.max.y;
        patch_ij_bbox_vec.push_back(patch_bbox);
    }

    return true;
}

//----------------------------------------------------------------------
mi::math::Bbox<mi::Uint32, 3> Distributed_amplitude_map_job::get_heightfield_ijk_bbox(
    mi::base::Handle< nv::index::IRegular_heightfield_data_access >& heightfield_data_access)
{
    assert(heightfield_data_access.is_valid_interface());

    const mi::math::Bbox_struct<mi::Uint32, 2>& patch_bbox = heightfield_data_access->get_patch_bounding_box();

    // Retrieved height values
    const mi::math::Bbox< mi::Sint64, 2> patch_bbox_s64_2(
        nv::index_common::convert_bbox_type<mi::Sint64, mi::Uint32, 2>(
            mi::math::Bbox< mi::Uint32, 2>(patch_bbox)));

    // first get patch min max value
    const mi::math::Vector< mi::Sint64, 2> height_range =
        this->get_heightfield_patch_height_range(patch_bbox_s64_2,
                                                 heightfield_data_access->get_elevation_values());

    const mi::math::Bbox<mi::Uint32, 3> heightfield_bbox(patch_bbox.min.x,
                                                   patch_bbox.min.y,
                                                   static_cast< mi::Uint32 >(height_range[0]),
                                                   patch_bbox.max.x,
                                                   patch_bbox.max.y,
                                                   static_cast< mi::Uint32 >(height_range[1]));
    return heightfield_bbox;
}


//----------------------------------------------------------------------
mi::math::Vector< mi::Sint64, 2 > Distributed_amplitude_map_job::get_heightfield_patch_height_range(
    mi::math::Bbox< mi::Sint64, 2 > const & patch_bbox_s64_2,
    mi::Float32* p_height_values)
{
    mi::Float32 height_min = 1.0e23f;
    mi::Float32 height_max = 0.0f;
    for(mi::Sint64 i = patch_bbox_s64_2.min[0]; i < patch_bbox_s64_2.max[0]; ++i){
        for(mi::Sint64 j = patch_bbox_s64_2.min[1]; j < patch_bbox_s64_2.max[1]; ++j){
            mi::Sint64 const idx =
                nv::index_common::get_heightfield_index(mi::math::Vector< mi::Sint64, 2 >(i, j), patch_bbox_s64_2);
            if(p_height_values[idx] < 0.0f){
                continue;
            }
            if(p_height_values[idx] > height_max){
                height_max = p_height_values[idx];
            }
            if(p_height_values[idx] < height_min){
                height_min = p_height_values[idx];
            }
        }
    }
    // DEBUG_LOG << "height_min = " << height_min << ", height_max = " << height_max;
    return mi::math::Vector< mi::Sint64, 2 >(static_cast< mi::Sint64 >(floor(height_min)),
                                             static_cast< mi::Sint64 >(ceil( height_max)));
}

//----------------------------------------------------------------------
