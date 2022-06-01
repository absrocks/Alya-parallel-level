/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief external heightfield ppm data importer

#ifndef NVIDIA_INDEX_PPM_HEIGHTFIELD_IMPORTER_H
#define NVIDIA_INDEX_PPM_HEIGHTFIELD_IMPORTER_H

#include <mi/dice.h>

#include <nv/index/idistributed_data_import_callback.h>

#include <string>

namespace nv {
namespace index_common {
/// Heightfield ppm file importer.
///
class PPM_heightfield_importer :
        public nv::index::Distributed_discrete_data_import_callback<0x5762c097,0x713d,0x47e0,0xaf,0x32,0x37,0x1f,0xcd,0x74,0x3b,0xc5>
{
public:
    /// default constructor
    PPM_heightfield_importer();

    /// constructor
    ///
    /// elevation value = scale * (file value) + offset
    ///
    /// \param[in] filename     heightfield ppm file name
    /// \param[in] scale        elevation value scaling factor
    /// \param[in] offset       elevation value offset factor
    /// \param[in] binary_mask_filename binary mask file name for holes
    /// \param[in] size         heightfield patch whole size
    ///
    PPM_heightfield_importer(
        const std::string&                                          filename,
        mi::Float32                                                 scale,
        mi::Float32                                                 offset,
        const std::string&                                          binary_mask_filename,
        const mi::math::Vector_struct<mi::Uint32, 2>&               size);

    /// destructor
    virtual ~PPM_heightfield_importer();

    // -------------------------------------------------------------------------
    // Implements the methods of the importer callback interface \c IDistributed_data_import_callback to
    // estimate memory consuptions and to create a volume brick.
    virtual mi::Size estimate(
        const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
        mi::neuraylib::IDice_transaction*               dice_transaction) const;
    
    virtual nv::index::IDistributed_data_subset* create(
        const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
        nv::index::IData_subset_factory*                factory,
        mi::neuraylib::IDice_transaction*               dice_transaction) const;

    virtual mi::base::Uuid subset_id() const;

    /// -------------------------------------------------------------------------------------------
    /// Returns configuration settings that may be used by the library for the session export
    /// mechanism provided by ISession::export_session().
    virtual const char* get_configuration() const;

    /// -------------------------------------------------------------------------------------------
    /// Serialize the class to the given serializer.
    virtual void serialize(mi::neuraylib::ISerializer * serializer) const;

    /// Deserialize the class from the given deserializer.
    virtual void deserialize(mi::neuraylib::IDeserializer * deserializer);

private:
    /// Heightfield ppm filename
    std::string                                                     m_filename;
    /// Elevation value scaling factor
    mi::Float32                                                     m_scale;
    /// Elevation value offset
    mi::Float32                                                     m_offset;
    /// Binary mask file name for holes
    std::string                                                     m_binary_mask_filename;
    /// Heightfield size
    mi::math::Vector<mi::Uint32, 2>                                 m_size;

    // Configuration information for the session export
    std::string                                                     m_configuration;
};

//----------------------------------------------------------------------
}} // namespace nv::index_common
#endif // NVIDIA_INDEX_PPM_HEIGHTFIELD_IMPORTER_H
