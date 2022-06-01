/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Public volume data importer based on a simple raw data format.

#ifndef NVIDIA_INDEX_BIN_COMMON_VOLUME_REPLICA_GENERATION_H
#define NVIDIA_INDEX_BIN_COMMON_VOLUME_REPLICA_GENERATION_H

#include <mi/dice.h>
#include <mi/base/interface_declare.h>
#include <mi/base/interface_implement.h>

#include <nv/index/idistributed_data_access.h>
#include <nv/index/idistributed_data_import_callback.h>

#include <string>

namespace nv {
namespace index {
class IRegular_volume_brick_uint8;
}
}

namespace nv {
namespace index_common {

/// Volume replica generation.
/// This copies a volume using importer facility.
class Volume_replica_generation :
    public nv::index::Distributed_discrete_data_import_callback<0xa2c6327d,0xa271,0x41f4,0xb9,0x75,0x05,0x28,0x34,0x91,0x15,0xa9>
{
public:
    /// Constructor
    Volume_replica_generation(
        const mi::neuraylib::Tag&                   session_tag,
        const mi::neuraylib::Tag&                   input_volume_tag,
        const mi::math::Bbox_struct<mi::Uint32, 3>& ijk_partial_bounds);

    /// The default constructor is required for serialization only.
    Volume_replica_generation();

    // -------------------------------------------------------------------------
    // Implements the methods of the importer callback interface \c IDistributed_data_import_callback to
    // estimate memory consumption and to create a volume brick.
    virtual mi::Size estimate(
        const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
        mi::neuraylib::IDice_transaction*               dice_transaction) const;
    
    virtual nv::index::IDistributed_data_subset* create(
        const mi::math::Bbox_struct<mi::Sint32, 3>&     bbox,
        nv::index::IData_subset_factory*                factory,
        mi::neuraylib::IDice_transaction*               dice_transaction) const;

    virtual mi::base::Uuid subset_id() const;

    // -------------------------------------------------------------------------------------------
    /// Returns configuration settings that may be used by the library for the session export
    /// mechanism provided by ISession::export_session().
    virtual const char* get_configuration() const;

    // -------------------------------------------------------------------------------------------
    /// Serialize the class to the given serializer.
    /// \param[in] serializer Write to this serializer.
    virtual void serialize(mi::neuraylib::ISerializer* serializer) const;

    /// Deserialize the class from the given deserializer.
    /// \param[in] deserializer Read from this deserializer.
    virtual void deserialize(mi::neuraylib::IDeserializer* deserializer);

private:
    /// Retrieve the source data to the voxel_data.
    ///
    /// \param[in] bounds             volume data bounding box
    /// \param[in] volume_data_access data volume data access object
    /// \param[in,out] volume_brick   (in,out) volume data to be filled
    /// \param[in] dice_transaction   the db transaction
    /// \return true when succeeded
    bool retrieve_source_data(
        const mi::math::Bbox<mi::Sint32, 3>&    bounds,
        nv::index::IRegular_volume_data_access* volume_data_access,
        nv::index::IRegular_volume_brick_uint8* volume_brick,
        mi::neuraylib::IDice_transaction*       dice_transaction
        ) const;

private:
    // Session that allows retrieving the input volume.
    mi::neuraylib::Tag                      m_session_tag;
    
    // Input volume that shall be replicated.
    mi::neuraylib::Tag                      m_input_volume_tag;

    // ..
    mi::math::Bbox_struct<mi::Uint32, 3>    m_ijk_partial_bbox;

    /// configuration information for the session exporter
    std::string                             m_configuration;
};

}} // namespace nv::index_common

#endif // NVIDIA_INDEX_BIN_COMMON_VOLUME_REPLICA_GENERATION_H
