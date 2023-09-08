/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief external sparse volume data import job. It loads a subset of an sparse volume 
/// contained in an given sub-region.

#ifndef NVIDIA_INDEX_BIN_SPARSE_VOLUME_IMPORTER_H
#define NVIDIA_INDEX_BIN_SPARSE_VOLUME_IMPORTER_H

#include <string>

#include <mi/dice.h>

#include <nv/index/idistributed_data_import_callback.h>
#include <nv/index/isparse_volume_subset.h>

// disable OpenVDB warnings
#ifdef WIN_NT
#pragma warning(push)
#pragma warning(disable: 4800 4244 4146)
#endif

#ifdef LINUX
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#endif

#include <openvdb/openvdb.h>

#ifdef WIN_NT
#pragma warning(pop)
#endif

#ifdef LINUX
#pragma GCC diagnostic pop
#endif


namespace nv {
namespace index_common {

/// External sparse volume importer implementation
class Sparse_volume_importer :
    public nv::index::Distributed_discrete_data_import_callback<0xea29f10f,0x6914,0x4a32,0x98,0x0,0x51,0x89,0x14,0x52,0x44,0x5a>
{
public:
    struct OpenVDB_cache_file_info
    {
        std::string             m_file_name;
        openvdb::io::File*      m_file;
        openvdb::GridBase::Ptr  m_grid_attrib_density;

        void clear();

        OpenVDB_cache_file_info();
        ~OpenVDB_cache_file_info();
    };

public:
    Sparse_volume_importer();
    Sparse_volume_importer(
        const std::string& input_dir,
        const std::string& cache_out_dir,
        const std::string& file_base_name,
        const std::string& file_ext,
        const std::string& field_name,
        const mi::Uint32   ts_enum_len,
        const mi::Uint32   ts_enum_stride,
        const mi::Uint32   ts_enum_off,
        bool               use_cache = false);
    virtual ~Sparse_volume_importer();
    
    // Implements the methods of the importer callback interface \c IDistributed_data_import_callback to
    // estimate memory consumption and to create a volume brick.
    virtual mi::Size estimate(
        const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
        mi::neuraylib::IDice_transaction*               dice_transaction) const;
    
    virtual nv::index::IDistributed_data_subset* create(
        const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
        nv::index::IData_subset_factory*                factory,
        mi::neuraylib::IDice_transaction*               dice_transaction) const;

    virtual nv::index::IDistributed_data_subset* create(
        const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
        mi::Uint32                                      time_step,
        nv::index::IData_subset_factory*                factory,
        mi::neuraylib::IDice_transaction*               dice_transaction) const;

    virtual mi::base::Uuid subset_id() const
    {
        return nv::index::ISparse_volume_subset::IID();
    }

    /// Returns configuration settings that may be used by the library for the session export
    /// mechanism provided by ISession::export_session().
    virtual const char* get_configuration() const;

    /// Serialize the class to the given serializer.
    /// \param serializer Write to this serializer.
    virtual void serialize(mi::neuraylib::ISerializer* serializer) const;

    /// De-serialize the class from the given deserializer.
    /// \param deserializer Read from this deserializer.
    virtual void deserialize(mi::neuraylib::IDeserializer* deserializer);

protected:
    // Implements the core file loading. Will be used by both create() variants internally.
    nv::index::IDistributed_data_subset* create(
        const std::string&                              filename,
        const std::string&                              filename_cache,
        const std::string&                              field_name,
        const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
        nv::index::IData_subset_factory*                factory,
        mi::neuraylib::IDice_transaction*               dice_transaction) const;

    bool load_from_vdb(
        const std::string&                          filename,
        const mi::math::Bbox_struct<mi::Sint32, 3>& bounding_box,
        const std::string&                          data_field_name,
        nv::index::ISparse_volume_subset*           svol_subset) const;

protected:
    std::string                             m_input_directory;          ///< Sparse volume data input directory filename
    std::string                             m_cache_output_directory;   ///< Directory where the generated cache files are directed to
    std::string                             m_file_base_name;           ///< 
    std::string                             m_file_extension;           ///< 

    std::string                             m_data_field_name;          ///< the data field to load from vdb files

    mi::Uint32                              m_timestep_enumeration_length;
    mi::Uint32                              m_timestep_enumeration_stride;
    mi::Uint32                              m_timestep_enumeration_offset;

    std::string                             m_configuration;            ///< Configuration information for the session exporter

    bool                                    m_use_cache;        ///< Use cache file

    static mi::base::Lock                   m_cached_input_lock;
    static OpenVDB_cache_file_info          m_cached_input_file;
};

} // namespace index_common
} // namespace nv

#endif // NVIDIA_INDEX_BIN_SPARSE_VOLUME_IMPORTER_H
