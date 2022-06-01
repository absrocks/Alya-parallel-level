/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Heightfield data importer reading raw binary formatted elevation (float) and normal (float,3).

#ifndef NVIDIA_INDEX_BIN_COMMON_RAW_HEIGHTFIELD_DATA_IMPORTER_H
#define NVIDIA_INDEX_BIN_COMMON_RAW_HEIGHTFIELD_DATA_IMPORTER_H

#include <mi/dice.h>
#include <mi/base/interface_declare.h>
#include <mi/base/interface_implement.h>

#include <nv/index/idistributed_data_import_callback.h>
#include <nv/index/iregular_heightfield_patch.h>

#include <string>

#include "raw_heightfield_data_type.h"
#include "string_dict.h"

namespace nv {
namespace index_common {

/// Heightfield raw data importer.
///
/// - height value (Float32)
/// - normal       (Float32,3)
///
/// File structure.
///---
/// #! index_raw_heightfield 1
/// # eventname = heightfield_test
/// # size = 1024 1024
/// #... (key = value, project format)
/// #-end
/// height value binary data...
/// normal value binary data...
///---
class Raw_heightfield_data_importer :
        public nv::index::Distributed_discrete_data_import_callback<0x79c7e58e,0x59f0,0x45d5,0x83,0x4c,0x22,0x6a,0x2f,0xb9,0x8f,0xc3>
{
public:
    /// constructor
    ///
    /// \param[in] filename raw heightfield file name
    /// \param[in] size     heightfield patch size
    Raw_heightfield_data_importer(
        const std::string&                            filename,
        const mi::math::Vector_struct<mi::Uint32, 2>& size);

    /// default constructor
    Raw_heightfield_data_importer();

    /// destructor
    virtual ~Raw_heightfield_data_importer();

    // -------------------------------------------------------------------------
    // Implements the methods of the importer callback interface \c IDistributed_data_import_callback to
    // estimate memory consumption and to create a volume brick.
    virtual mi::Size estimate(
        const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
        mi::neuraylib::IDice_transaction*               dice_transaction) const;
    
    virtual nv::index::IDistributed_data_subset* create(
        const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
        nv::index::IData_subset_factory*                factory,
        mi::neuraylib::IDice_transaction*               dice_transaction) const;

    virtual mi::base::Uuid subset_id() const
    {
        return nv::index::IRegular_heightfield_patch::IID();
    }

    // -------------------------------------------------------------------------
    /// Returns configuration settings that may be used by the library for the session export
    /// mechanism provided by ISession::export_session().
    virtual const char* get_configuration() const;

    //----------------------------------------------------------------------
    virtual const char* get_file_name() const { return m_filename.c_str(); }
    virtual const mi::math::Vector_struct<mi::Uint32, 2> & get_heightfield_size() const { return m_heightfield_size; }

    /// Use cache files to improve read performance
    virtual void set_use_cache(bool use_cache) { m_use_cache = use_cache; }
    virtual bool get_use_cache() const { return m_use_cache; }

    /// Recompute normals instead of loading them from file
    virtual void set_regenerate_normals(bool regenerate) { m_regenerate_normals = regenerate; }
    virtual bool get_regenerate_normals() const { return m_regenerate_normals; }

    /// Only allow one thread at a time to read a file, to improve performance with spinning disks
    virtual void set_serial_access(bool serial_access) { m_serial_access = serial_access; }
    virtual bool get_serial_access() const { return m_serial_access; }

    /// Print statistics when importing the height field data.
    virtual void set_show_stats(bool show_stats) { m_show_stats = show_stats; }
    virtual bool get_show_stats() const { return m_show_stats; }
    //----------------------------------------------------------------------

    /// Serialize the class to the given serializer.
    virtual void serialize(mi::neuraylib::ISerializer * serializer) const;

    /// Deserialize the class from the given deserializer.
    virtual void deserialize(mi::neuraylib::IDeserializer * deserializer);

private:
    /// Try to load a cache file.
    ///
    /// Try importing the data from a tempory cache file possibly
    /// generated beforehand if possible. 0, we have not yet had a
    /// cache file.
    ///
    /// \param[in]  bounding_box    bounding box of the patch
    /// \param[in]  factory         data (heightfield) factory
    /// \param[out] cache_file_name raw heightfield cache file name 
    /// 
    /// \return heightfield patch, 0 when failed to load the
    /// heightfield from the cache file.
    nv::index::IRegular_heightfield_patch* try_load_from_cache_file(
        const mi::math::Bbox_struct<mi::Sint32, 3>& bounding_box,
        const mi::Sint64                            grid_size,
        nv::index::IData_subset_factory*            factory,
        std::string&                                cache_file_name) const;
    
    /// Check the consistency between header file and the given information
    ///
    /// \param[in] header_dict  file header information in dict format
    /// \param[in] bounding_box given bounding box
    /// \return true when the header and the given info are consistent.
    bool check_header_consistency(
        const String_dict& header_dict, 
        const mi::math::Vector<mi::Sint64, 2> heightfield_size) const;

    /// Load the data
    void load_data(
        FILE*                                       p_file,
        const std::string&                          cache_file_name,
        Raw_heightfield_file_format_type            file_type,
        const mi::math::Vector<mi::Sint64, 2>&      heightfield_size,
        const mi::math::Bbox_struct<mi::Sint32, 3>& bounding_box,
        const mi::math::Bbox<mi::Sint64, 2>&        patch_bbox,
        mi::Float32*                                elevation_data,
        mi::math::Vector_struct<mi::Float32, 3>*    normal_vector_data) const;

    /// heightfield raw patch filename
    std::string                             m_filename;

    /// total heightfield size
    mi::math::Vector<mi::Uint32, 2>         m_heightfield_size;

    /// configuration information for the session exporter
    std::string                             m_configuration;

    /// use cache files to improve read performance
    bool                                    m_use_cache;
    /// recompute normals instead of loading them from file
    bool                                    m_regenerate_normals;
    /// only allow one thread at a time to read a file, to improve performance with spinning disks
    bool                                    m_serial_access;
    /// show some statistics when loading
    bool                                    m_show_stats;
    /// Verify imported elevation values using ppm file export
    bool                                    m_image_verification;

    static mi::Uint64                       s_total_read_bytes;
    static mi::Float64                      s_total_read_time;

    static mi::base::Lock                   s_file_access_lock;
};

//----------------------------------------------------------------------
/// Check raw heightfield file magic
///
/// \param[in] infname    input file name
/// \param[out] file_type file type in the header
/// \return    true when valid heightfield file magic
bool is_raw_heightfield_file_magic(
    const std::string&                infname,
    Raw_heightfield_file_format_type& file_type);

//----------------------------------------------------------------------
/// Load raw heightfield header.
/// Note: The file magic should have been checked.
///
///-----
/// #! index_raw_heightfield 1
/// #eventname = heightfield_test
/// #size = 1024 1024
/// #... (key = value, project format)
/// #-end
/// heightfield binary data Float32           ...
/// normal binary data      Vector<Float32,3> ...
///-----
/// Currently raw data file_type is RHFFT_INDEX_RAW_v2 (== 2) at the first magic line
/// #! index_raw_heightfield 2
///
bool load_commented_header_part(FILE* p_file,
                                nv::index_common::String_dict & header_dict);

//----------------------------------------------------------------------
} // namespace index_common
} // namespace nv
#endif // NVIDIA_INDEX_BIN_COMMON_RAW_HEIGHTFIELD_DATA_IMPORTER_H
