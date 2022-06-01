/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief External triangle-mesh data import job.

#ifndef NVIDIA_INDEX_BIN_TRIANGLE_MESH_IMPORTER_H
#define NVIDIA_INDEX_BIN_TRIANGLE_MESH_IMPORTER_H

#include <mi/dice.h>

#include <nv/index/idistributed_data_import_callback.h>

#include <string>

namespace nv {
namespace index_common {

/// Triangle mesh importer implementation.
class Triangle_mesh_importer :
    public nv::index::Distributed_discrete_data_import_callback<0xd6f2c6e5,0x6c42,0x42a7,0xbf,0x0c,0x0f,0xb0,0x2a,0xdf,0x78,0xa6>
{
public:
    /// constructor
    Triangle_mesh_importer(const std::string& filename);

    /// default constructor
    Triangle_mesh_importer();

    /// destructor
    virtual ~Triangle_mesh_importer();

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

    // -------------------------------------------------------------------------
    /// Returns configuration settings that may be used by the library for the session export
    /// mechanism provided by ISession::export_session().
    virtual const char* get_configuration() const;

    // --------------------------------------------------------------------
    /// Serialize the class to the given serializer.
    /// \param serializer Write to this serializer.
    virtual void serialize(mi::neuraylib::ISerializer* serializer) const;

    /// Deserialize the class from the given deserializer.
    /// \param deserializer Read from this deserializer.
    virtual void deserialize(mi::neuraylib::IDeserializer* deserializer);

private:
    /// Triangle mesh's data filename
    std::string                             m_filename;

    /// configuration information for the session exporter
    std::string                             m_configuration;
};

}} // namespace nv::index_common

#endif // NVIDIA_INDEX_BIN_TRIANGLE_MESH_IMPORTER_H
