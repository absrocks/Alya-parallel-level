/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Volume data importer reading a simple raw data format from file.

#include "volume_replica_generation.h"

#include "forwarding_logger.h"
#include "common_utility.h"
#include "type_conversion_utility.h"
#include "encode_voxel.h"

#include <nv/index/iregular_volume_brick.h>

#include <nv/index/idistributed_data_locality.h>
#include <nv/index/isession.h>

#include <cstdio>

namespace nv {
namespace index_common {

//----------------------------------------------------------------------
Volume_replica_generation::Volume_replica_generation(
    const mi::neuraylib::Tag&                   session_tag,
    const mi::neuraylib::Tag&                   input_volume_tag,
    const mi::math::Bbox_struct<mi::Uint32, 3>& ijk_partial_bbox) 
    :
    m_session_tag(session_tag),
    m_input_volume_tag(input_volume_tag),
    m_ijk_partial_bbox(ijk_partial_bbox)
{
    // TODO: m_configuration ;
}

//----------------------------------------------------------------------
Volume_replica_generation::Volume_replica_generation()
    : 
    m_session_tag(mi::neuraylib::NULL_TAG),
    m_input_volume_tag(mi::neuraylib::NULL_TAG)
{
    // empty
}

//----------------------------------------------------------------------
mi::Size Volume_replica_generation::estimate(
    const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{
    // This is an 8-bit raw volume data importer:
    const mi::Size dx = bounding_box.max.x - bounding_box.min.x;
    const mi::Size dy = bounding_box.max.y - bounding_box.min.y;
    const mi::Size dz = bounding_box.max.z - bounding_box.min.z;
    const mi::Size volume_brick_size = dx * dy * dz;
    return volume_brick_size;
}

//----------------------------------------------------------------------
nv::index::IDistributed_data_subset* Volume_replica_generation::create(
    const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
    nv::index::IData_subset_factory*                factory,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{
    const mi::math::Bbox<mi::Sint32, 3> bounds = bounding_box;
    
    // Verify if the the scene element tag is valid.
    if(!m_input_volume_tag.is_valid())
    {
        ERROR_LOG << "Invalid scene element tag. No brick generated for " << bounds << ".";
        return 0;
    }

    // The size for volume data per subregion.
    mi::math::Bbox< mi::Sint64, 3 > raw_dst_bbox(
        nv::index_common::convert_bbox_type<mi::Sint64, mi::Sint32, 3>(bounds));

    INFO_LOG << "Volume_copy_generation copy. Output bbox: " << bounds
             << ", raw_dst_bbox = " << raw_dst_bbox;

    // Access the session
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<const nv::index::ISession>(m_session_tag));
    assert(session.is_valid_interface());

    // Access the data access factory
    const mi::neuraylib::Tag& data_access_tag = session->get_data_access_factory();
    mi::base::Handle<const nv::index::IDistributed_data_access_factory> access_factory(
        dice_transaction->access<const nv::index::IDistributed_data_access_factory>(data_access_tag));
    assert(access_factory.is_valid_interface());

    // Query the volume data
    mi::base::Handle<nv::index::IRegular_volume_data_access> volume_data_access(
        access_factory->create_regular_volume_data_access(m_input_volume_tag));
    assert(volume_data_access.is_valid_interface());

    // Create a volume brick
    mi::base::Handle<nv::index::IRegular_volume_brick_uint8> volume_brick(
        factory->create<nv::index::IRegular_volume_brick_uint8>());
    if (!volume_brick.is_valid_interface())
    {
        ERROR_LOG << "Cannot create a volume brick in Volume_replica_generation.";
        return 0;
    }

    // Allocate storage and fill it with the source data
    if (!this->retrieve_source_data(bounds, volume_data_access.get(), volume_brick.get(), dice_transaction))
    {
        ERROR_LOG << "Failed to retrieve the voxel data.";
        return 0;
    }

    // INFO_LOG << "Data copied.";

    volume_brick->retain();     // Since the handle will be out of scope. 
    return volume_brick.get();
}

//----------------------------------------------------------------------
mi::base::Uuid Volume_replica_generation::subset_id() const
{
    return nv::index::IRegular_volume_brick_uint8::IID();
}

//----------------------------------------------------------------------
const char* Volume_replica_generation::get_configuration() const
{
    return m_configuration.c_str();
}

//----------------------------------------------------------------------
void Volume_replica_generation::serialize(
    mi::neuraylib::ISerializer* serializer) const
{
    serializer->write(&m_session_tag.id, 1);
    serializer->write(&m_input_volume_tag.id, 1);
    serializer->write(&m_ijk_partial_bbox.min.x, 6);
}

//----------------------------------------------------------------------
void Volume_replica_generation::deserialize(
    mi::neuraylib::IDeserializer* deserializer)
{
    deserializer->read(&m_session_tag.id, 1);
    deserializer->read(&m_input_volume_tag.id, 1);
    deserializer->read(&m_ijk_partial_bbox.min.x, 6);
}

//----------------------------------------------------------------------
bool Volume_replica_generation::retrieve_source_data(
    const mi::math::Bbox<mi::Sint32, 3>&    bounds,
    nv::index::IRegular_volume_data_access* volume_data_access,
    nv::index::IRegular_volume_brick_uint8* volume_brick,
    mi::neuraylib::IDice_transaction*       dice_transaction
    ) const
{
    assert(volume_data_access != 0);
    assert(volume_brick       != 0);
    assert(dice_transaction   != 0);

    DEBUG_LOG << "Source volume ijk ROI bbox = " << m_ijk_partial_bbox
              << ", brick bbox = " << bounds;

    //----------------------------------------------------------------------
    // Allocate the voxel storage
    mi::Uint8* voxel_data = volume_brick->generate_voxel_storage(bounds);
    if (voxel_data == 0)
    {
        ERROR_LOG << "Cannot generate voxel storage in Volume_replica_generation.";
        return false;
    }

    //----------------------------------------------------------------------
    // compute source volume bbox
    mi::math::Bbox< mi::Sint64, 3 > output_bbox(
        nv::index_common::convert_bbox_type<mi::Sint64, mi::Sint32, 3>(bounds));
    const mi::math::Bbox< mi::Sint64, 3 > whole_src_roi_bbox(
        nv::index_common::convert_bbox_type<mi::Sint64, mi::Uint32, 3>(m_ijk_partial_bbox));
    const mi::math::Vector< mi::Sint64, 3 > whole_src_size_vec3 =
        whole_src_roi_bbox.max - whole_src_roi_bbox.min;
    const mi::math::Bbox< mi::Sint64, 3 > whole_src_size_bbox(
        mi::math::Vector< mi::Sint64, 3 >(0, 0, 0), whole_src_size_vec3);

    mi::math::Bbox< mi::Sint64, 3 > copy_src_roi_bbox(
        whole_src_roi_bbox.min + output_bbox.min,
        whole_src_roi_bbox.min + output_bbox.max);

    // clip to the valid range
    copy_src_roi_bbox = bbox_and(copy_src_roi_bbox, whole_src_roi_bbox);

    //----------------------------------------------------------------------
    // compute destination (brick) volume bbox
    mi::math::Bbox< mi::Sint64, 3 > raw_dst_bbox = output_bbox;
    const mi::math::Bbox< mi::Sint64, 3 >  dst_valid_bbox = bbox_and(raw_dst_bbox, whole_src_size_bbox);

    // acces to the source volume data
    const mi::math::Bbox_struct< mi::Uint32, 3 > copy_src_roi_bbox_st =
        convert_bbox_type<mi::Uint32, mi::Sint64, 3>(copy_src_roi_bbox);
    INFO_LOG << "Volume_copy_generation::access_source_data: copy source bounding box = "
             << copy_src_roi_bbox_st;
    if (volume_data_access->access(copy_src_roi_bbox_st, dice_transaction) < 0)
    {
        ERROR_LOG << "Failed to volume data access. bbox: " << copy_src_roi_bbox_st;
        return false;
    }

    const mi::base::Handle<const nv::index::IRegular_volume_data> volume_data(
        volume_data_access->get_volume_data());
    const mi::base::Handle<const nv::index::IRegular_volume_data_uint8> volume_data_uint8(
        volume_data->get_interface<const nv::index::IRegular_volume_data_uint8>());

    if (!volume_data_uint8) {
        ERROR_LOG << "Failed to access volume data: data access on non-uint8 volume type.";
        return false;
    }

    // Copy the voxel data to brick's data storage
    const mi::Uint8* const src_data = volume_data_uint8->get_voxel_data();
    const bool success = 
        nv::index_common::copy_partial_volume(voxel_data,        raw_dst_bbox,
                                              src_data,          copy_src_roi_bbox,
                                              copy_src_roi_bbox, dst_valid_bbox.min);
    if (!success)
    {
        ERROR_LOG << "Failed to copy the voxel data.";
        return false;
    }

    // FIXME: fill the brick border
    // mi::math::Bbox< mi::Sint64, 3 > const expanded_bbox =
    //     nv::index_common::expand_border_by_copy(brick_data, raw_dst_bbox, dst_valid_bbox);
    // assert(expanded_bbox == raw_dst_bbox);
    // nv::index_common::no_unused_variable_warning_please(expanded_bbox);

    // #ifdef DEBUG
    //         // This test is correct iff the source volume is synthetic and
    //         // its type is SDGSART_IJ type.
    //         bool const r_src =
    //             test_syhthetic_voxel_value("volume_data_gen_copy_job: src data",
    //                                         src_bbox,   src_bbox,
    //                                         SDVE_IJ,    src_data);
    //         // did not construct the border. only check the raw_query_bbox range.
    //         bool const r_dst =
    //             test_syhthetic_voxel_value("volume_data_gen_copy_job: dst data",
    //                                         raw_dst_bbox, src_bbox,
    //                                         SDVE_IJ,      brick_data);
    //         assert(r_src);
    //         assert(r_dst);
    // #endif // DEBUG

    return true;
}

//----------------------------------------------------------------------
}} // namespace nv::index_common
