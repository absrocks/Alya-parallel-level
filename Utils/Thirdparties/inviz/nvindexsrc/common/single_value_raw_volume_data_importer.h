/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Public volume data importer based on a simple scalar data format.

#ifndef NVIDIA_INDEX_BIN_COMMON_SINGLE_VALUE_RAW_VOLUME_DATA_IMPORTER_H
#define NVIDIA_INDEX_BIN_COMMON_SINGLE_VALUE_RAW_VOLUME_DATA_IMPORTER_H

#include <string>

#include <mi/dice.h>
#include <mi/base/interface_declare.h>
#include <mi/base/interface_implement.h>

#include <nv/index/idistributed_data_import_callback.h>
#include <nv/index/iregular_volume_brick.h>

namespace nv {
namespace index_common {

/// The scalar volume data importer reads volume brick data from file.
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
class Single_value_raw_volume_data_importer :
    public nv::index::Distributed_discrete_data_import_callback<0xd2fe0d1b,0xc9dd,0x4e40,0xbc,0x51,0x0b,0x23,0x0f,0x81,0xc4,0x82>
{
public:
    /// The raw volume data importer reads 8-bit volume brick data from file.
    /// The importer requires the file name/location as well as the size/extent of the
    /// dataset. 
    ///
    /// \param[in] file_name    The file name/location of the volume dataset file.
    /// \param[in] src_folder   source dataset file directory
    /// \param[in] dst_folder   destination dataset file directory
    /// \param[in] size         The size/extent of the entire volume dataset.
    ///
    Single_value_raw_volume_data_importer(
        const std::string&                            file_name,
        const std::string&                            src_folder,
        const std::string&                            dst_folder,
        const mi::math::Vector_struct<mi::Uint32, 3>& size);

    /// The default constructor is required for serialization only.
    Single_value_raw_volume_data_importer();
    virtual ~Single_value_raw_volume_data_importer();
    
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

    virtual nv::index::IDistributed_data_subset* create(
        const mi::math::Bbox_struct<mi::Sint32, 3>&     bbox,
        mi::Uint32                                      time_step,
        nv::index::IData_subset_factory*                factory,
        mi::neuraylib::IDice_transaction*               dice_transaction) const;

    virtual mi::base::Uuid subset_id() const { return nv::index::IRegular_volume_brick_float32::IID(); }
    // -------------------------------------------------------------------------

    /// Normalize the voxel values to the range [0,1]
    virtual void set_range(const mi::math::Vector<mi::Float32, 2>& range) { m_range = range; }

    // Adding additional time steps for a time lapsed animation.
    virtual void add_time_step(const std::string& file);

    // 
    virtual void set_b_frame_offset(mi::Uint32 b_frame_offset) { m_b_frame_offset=b_frame_offset; }
    virtual void enable_difference_encoding(bool use)          { m_use_difference_encoding = use; }

    /// get stored file name
    /// \return volume data file name
    virtual const char* get_file_name() const { return m_file_name.c_str(); }

    /// get volume data size
    /// \return volume data size
    virtual const mi::math::Vector_struct<mi::Uint32, 3>& get_volume_size() const { return m_size; }

    /// use cache files to improve read performance
    virtual void set_use_cache(bool use_cache) { m_use_cache = use_cache; }
    virtual bool get_use_cache() const { return m_use_cache; }

    /// use gzip to compress cache files, 0 means uncompresed
    virtual void set_cache_compression(mi::Uint32 compression_level) { m_cache_compression = compression_level; }
    virtual bool get_cache_compression() const { return m_cache_compression != 0u; }

    /// only allow one thread at a time to read a file, to improve performance with spinning disks
    virtual void set_serial_access(bool serial_access) { m_serial_access = serial_access; }
    virtual bool get_serial_access() const { return m_serial_access; }

    /// show some statistics when loading
    virtual void set_show_stats(bool show_stats) { m_show_stats = show_stats; }
    virtual bool get_show_stats() const { return m_show_stats; }

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
    nv::index::IDistributed_data_subset* create(
        const std::string&                              file_name,
        bool                                            is_b_frame,
        const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
        nv::index::IData_subset_factory*                factory,
        mi::neuraylib::IDice_transaction*               dice_transaction) const;

    std::string                             m_file_name;
    std::string                             m_src_folder;
    std::string                             m_dst_folder;
    mi::math::Vector<mi::Uint32, 3>         m_size;

    mi::math::Vector<mi::Float32, 2>        m_range;
    std::vector<std::string>                m_time_steps;
    mi::Uint32                              m_b_frame_offset;
    bool                                    m_use_difference_encoding;

    /// use cache files to improve read performance
    bool                                    m_use_cache;

    /// cache gzip compression level, 0 means uncompresed
    mi::Uint32                              m_cache_compression;

    /// only allow one thread at a time to read a file, to improve performance with spinning disks
    bool                                    m_serial_access;

    /// show some statistics when loading
    bool                                    m_show_stats;

    /// configuration information for the session exporter
    std::string                             m_configuration;

    // static variables.
    static mi::base::Lock                   s_file_access_lock;
    static mi::Uint64                       s_total_read_bytes;
    static mi::Float64                      s_total_read_time;
};

}} // namespace nv::index_common

#endif // NVIDIA_INDEX_BIN_COMMON_SINGLE_VALUE_RAW_VOLUME_DATA_IMPORTER_H
