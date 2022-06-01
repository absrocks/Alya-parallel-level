/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Public volume data importer based on a simple raw data format.

#ifndef NVIDIA_INDEX_BIN_COMMON_MULTI_ATTRIBUTE_SCALING_SEQUENCE_IMPORTER_H
#define NVIDIA_INDEX_BIN_COMMON_MULTI_ATTRIBUTE_SCALING_SEQUENCE_IMPORTER_H

#include <mi/dice.h>
#include <mi/base/interface_declare.h>
#include <mi/base/interface_implement.h>

#include <nv/index/idistributed_data_import_callback.h>

#include <map>
#include <string>

namespace nv {
namespace index_common {

/// The Multi_attribute_scaling_sequence_importer reads different attributes and combines
/// them into one attribute.
/// The class implement the \c IDistributed_data_import_callback interfaces class.

///
class Multi_attribute_scaling_sequence_importer :
    public nv::index::Distributed_discrete_data_import_callback<0xa8c7b57a,0x1c3d,0x4d68,0xb1,0x05,0xf3,0x09,0x6e,0xca,0x06,0x1a>
{
public:
    Multi_attribute_scaling_sequence_importer(
        const std::string&                            attribute_file_1,
        const std::string&                            attribute_file_2,
        const std::string&                            output_file,
        const mi::math::Vector_struct<mi::Uint32, 3>& size);

    /// The default constructor is required for serialization only.
    Multi_attribute_scaling_sequence_importer();
    
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

    virtual mi::base::Uuid subset_id() const;

    /// Returns configuration settings that may be used by the library for the session export
    /// mechanism provided by ISession::export_session().
    virtual const char* get_configuration() const;

    /// get stored file name
    /// \return volume data file name
    virtual const char* get_file_name() const { return m_attribute_file_1.c_str(); }

    /// get volume data size
    /// \return volume data size
    virtual const mi::math::Vector_struct<mi::Uint32, 3>& get_volume_size() const { return m_size; }
    
    /// set min max values of an attribute
    virtual void set_attribute_minmax(mi::Uint32 attribute_id, mi::math::Vector<mi::Float32, 2> min_max) 
    { 
        m_attribute_minmax_map[attribute_id] = min_max; 
    }
    
    /// set/get time interval between timesteps
    virtual void set_time_interval(mi::Uint32 time_interval) { m_time_interval = time_interval; }
    virtual mi::Uint32 get_time_interval() const { return m_time_interval; }
    
    /// set/get start time of the timesteps
    virtual void set_start_time(mi::Uint32 start_time) { m_start_time = start_time; }
    virtual mi::Uint32 get_start_time() const { return m_start_time; }
    
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

    // -------------------------------------------------------------------------------------------
    /// Serialize the class to the given serializer.
    /// \param serializer Write to this serializer.
    virtual void serialize(mi::neuraylib::ISerializer* serializer) const;

    /// Deserialize the class from the given deserializer.
    /// \param deserializer Read from this deserializer.
    virtual void deserialize(mi::neuraylib::IDeserializer* deserializer);

private:
    nv::index::IDistributed_data_subset* create(
        const std::string&                              attribute_file_1,
        const std::string&                              attribute_file_2,
        const std::string&                              output_file,
        const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
        nv::index::IData_subset_factory*                factory,
        mi::neuraylib::IDice_transaction*               dice_transaction) const;

    std::string                             m_attribute_file_1;
    std::string                             m_attribute_file_2;
    std::string                             m_output_file;
    mi::math::Vector<mi::Uint32, 3>         m_size;
    
    /// Attribute min/max values
    std::map<mi::Uint32, mi::math::Vector<mi::Float32, 2> > m_attribute_minmax_map;
    
    /// Time interval between each timestep
    mi::Uint32                              m_time_interval;
    
    /// Starting time of the timesteps
    mi::Uint32                              m_start_time;
    
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
#endif // NVIDIA_INDEX_BIN_COMMON_MULTI_ATTRIBUTE_SCALING_SEQUENCE_IMPORTER_H
