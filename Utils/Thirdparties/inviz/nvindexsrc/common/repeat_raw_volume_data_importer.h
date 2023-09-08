/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Public volume data importer based on a simple raw data format.

#ifndef NVIDIA_INDEX_BIN_COMMON_REPEAT_RAW_VOLUME_DATA_IMPORTER_H
#define NVIDIA_INDEX_BIN_COMMON_REPEAT_RAW_VOLUME_DATA_IMPORTER_H

#include <mi/dice.h>
#include <mi/base/interface_declare.h>
#include <mi/base/interface_implement.h>

#include <nv/index/idistributed_data_import_callback.h>

#include <string>

namespace nv {
namespace index_common {

/// The raw volume data importer reads volume brick data from file.
/// The class implement the \c IDistributed_data_import_callback interfaces class.
/// The implementation represents a very simplistic example of NVIDIA IndeX distributed
/// data import scheme and shall (1) demonstrate how to leverage the import callback to
/// generate a subset of a large-scale dataset in a parallel and distributed way and
/// (2) enabled application-writers and NVIDIA IndeX integrators to get started instantly.
///
/// The importer allows for cached brick data generated and stored in a previous sessions
/// in order to speed up a parallel import. Furthermore, the imported allows scaling up
/// a given dataset to "any" size to create larger datasets, e.g., for benchmarking.
///
class Repeat_raw_volume_data_importer :
    public nv::index::Distributed_discrete_data_import_callback<0x85500ff0,0x17fb,0x484b,0xb2,0x92,0x43,0x94,0xf9,0x36,0x95,0x69>
{
public:
    /// The raw volume data importer reads 8-bit volume brick data from file.
    /// The importer requires the file name/location as well as the size/extent of the
    /// dataset. 
    ///
    /// \param[in] file_name            The dataset's file name.
    /// \param[in] source_size          The size of the original volume dataset.
    /// \param[in] repeat               The scaling factor to scale the original dataset in each direction (ijk).
    /// \param[in] cache_source_data    If set to true, the entire source dataset will only be loaded once
    ///                                 per host and kept in memory. This can significantly speed up the
    ///                                 import process. However, the source dataset should not be too large
    ///                                 and this can only be used for a single dataset at the same time.
    ///
    Repeat_raw_volume_data_importer(
        const std::string&                          file_name,
        const mi::math::Vector<mi::Uint32, 3>&      source_size,
        const mi::math::Vector<mi::Uint32, 3>&      repeat,
        bool                                        cache_source_data);

    /// The default constructor is required for serialization only.
    Repeat_raw_volume_data_importer();
    ~Repeat_raw_volume_data_importer();

    /// The repeat importer modifies the size of the volume based on the size of the source dataset
    /// and the repeat factor. Since NVIDIA IndeX requires the size of the dataset that will be
    /// uploaded to the cluster the following method provides the modified size.
    ///
    /// \return         Returns the calculated size of the volume scene element.
    ///
    const mi::math::Vector<mi::Uint32, 3>& get_volume_size() const { return m_size; }
    
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
    /// \param serializer Write to this serializer.
    virtual void serialize(mi::neuraylib::ISerializer* serializer) const;

    /// Deserialize the class from the given deserializer.
    /// \param deserializer Read from this deserializer.
    virtual void deserialize(mi::neuraylib::IDeserializer* deserializer);

private:
    mi::Uint8* load_source_data(
        const std::string& file_name,
        mi::Uint64         source_size) const;

    /// Source data file name
    std::string                             m_file_name;
    
    /// Source volume size
    mi::math::Vector<mi::Uint32, 3>         m_source_size;
    
    /// Repeat parameter. Magnify factor for each dimension.
    mi::math::Vector<mi::Uint32, 3>         m_repeat;
    
    /// If 'true' then the data is cached in local memory
    bool                                    m_cache_source_data;
    
    /// Volume size after applying the repeat operation to the source dataset
    mi::math::Vector<mi::Uint32, 3>         m_size;

    /// Configuration information for the session exporter
    std::string                             m_configuration;

private:
    // Used for caching the source data and speed up the repeat-import.
    static mi::Uint8*                       m_source;
    static mi::base::Lock                   m_source_lock;
    static mi::base::Atom32                 m_cache_uses;
};

}} // namespace nv::index_common

#endif // NVIDIA_INDEX_BIN_COMMON_REPEAT_RAW_VOLUME_DATA_IMPORTER_H
