/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief external reservoir grid data generator job example implementation

#ifndef NVIDIA_INDEX_SDK_SYNTHETIC_RESERVOIRGRID_GENERATOR_H
#define NVIDIA_INDEX_SDK_SYNTHETIC_RESERVOIRGRID_GENERATOR_H

#include <mi/dice.h>

#include <nv/index/idistributed_data_import_callback.h>

#include <string>

namespace nv {
namespace index_common {

/// Synthetic reservoir patch generator.
/// An example implementation of massive data importer.
class Synthetic_reservoir_patch_generator :
    public nv::index::Distributed_discrete_data_import_callback<0x12511f11,0x84e4,0x4247,0xad,0xb7,0x3e,0x2a,0x9b,0x33,0xd6,0xc1>
{
public:
    /// constructor
    Synthetic_reservoir_patch_generator(
        const mi::Uint32& layers,
        const mi::Uint32& m_cells_per_row,
        const mi::Uint32& m_cells_per_column,
        const mi::math::Bbox_struct<mi::Float32, 3>& global_bbox,
        const std::string &                          cache_file,
        const std::string &                          cache_type);

    /// default constructor
    Synthetic_reservoir_patch_generator();

    /// destructor
    virtual ~Synthetic_reservoir_patch_generator();

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
    mi::Uint32                              m_layers;
    mi::Uint32                              m_cells_per_row;
    mi::Uint32                              m_cells_per_column;
    mi::math::Bbox_struct<mi::Float32, 3>   m_global_bbox;
    std::string                             m_cache_file;
    std::string                             m_cache_type;

    /// configuration information for the session exporter
    std::string                             m_configuration;
};

}} // namespace nv::index_common

#endif // NVIDIA_INDEX_SDK_SYNTHETIC_RESERVOIRGRID_GENERATOR_H
