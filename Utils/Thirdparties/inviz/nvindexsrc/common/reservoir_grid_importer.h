/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief External reservoir grid data importer example implementation

#ifndef NVIDIA_INDEX_SDK_RESERVOIR_GRID_IMPORTER_H
#define NVIDIA_INDEX_SDK_RESERVOIR_GRID_IMPORTER_H

#include <mi/dice.h>

#include <nv/index/idistributed_data_import_callback.h>

#include <string>

namespace nv {
namespace index_common {

/// External reservoir grid data importer example implementation
class Reservoir_grid_importer :
    public nv::index::Distributed_discrete_data_import_callback<0x71f9f1a3,0x014e,0x4791,0xb8,0x12,0x17,0x6a,0x56,0x88,0x1c,0x7d>
{
public:
    /// constructor
    Reservoir_grid_importer(
        const std::string &                             geometry_file,
        const std::string &                             scalar_file,
        const mi::math::Bbox_struct<mi::Float32, 3>&    bbox,
        const std::string &                             cache_type);

    /// default constructor
    Reservoir_grid_importer();

    /// destructor
    virtual ~Reservoir_grid_importer();

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

    /// Returns configuration settings that may be used by the library for the session export
    /// mechanism provided by ISession::export_session().
    virtual const char* get_configuration() const;

    // ----------------------------------------------------------------------------
    /// Serialize the class to the given serializer.
    virtual void serialize(mi::neuraylib::ISerializer * serializer) const;

    /// Deserialize the class from the given deserializer.
    virtual void deserialize(mi::neuraylib::IDeserializer * deserializer);

private:
    std::string                             m_geometry_file;
    std::string                             m_scalar_file;
    mi::math::Bbox_struct<mi::Float32, 3>   m_global_bbox;
    std::string                             m_cache_type;

    /// configuration information for the session exporter
    std::string                             m_configuration;
};

}} // namespace nv::index_common

#endif // NVIDIA_INDEX_SDK_RESERVOIR_GRID_IMPORTER_H
