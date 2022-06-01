/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief external irregular volume data import job for VTK files. 
/// It loads a subset of an irregular volume contained in an given voxel

#ifndef NVIDIA_INDEX_BIN_IRREGULAR_VOLUME_IMPORTER_VTK_H
#define NVIDIA_INDEX_BIN_IRREGULAR_VOLUME_IMPORTER_VTK_H

#include <string>
#include <vector>

#include <mi/dice.h>

#include <nv/index/idistributed_data_import_callback.h>
#include <nv/index/iirregular_volume_subset.h>

namespace nv {
namespace index_common {

/// External irregular volume importer implementation
class Irregular_volume_importer_vtk :
    public nv::index::Distributed_continuous_data_import_callback<0x48e5c2ed,0x2be7,0x4f7e,0xba,0x7b,0xe1,0xf8,0xf1,0x3a,0x32,0x68>
{
public:
    /// constructor
    Irregular_volume_importer_vtk(
        const std::string&  filename, 
        bool                use_cache = false);

    /// default constructor
    Irregular_volume_importer_vtk();

    /// destructor
    virtual ~Irregular_volume_importer_vtk();
    
    // -------------------------------------------------------------------------
    // Implements the methods of the importer callback interface \c IDistributed_data_import_callback to
    // estimate memory consumption and to create a volume brick.
    virtual mi::Size estimate(
        const mi::math::Bbox_struct<mi::Float32, 3>&    bounding_box,
        mi::neuraylib::IDice_transaction*               dice_transaction) const;
    
    virtual nv::index::IDistributed_data_subset* create(
        const mi::math::Bbox_struct<mi::Float32, 3>&    bbox,
        nv::index::IData_subset_factory*                factory,
        mi::neuraylib::IDice_transaction*               dice_transaction) const;

    virtual mi::base::Uuid subset_id() const
    {
        return nv::index::IIrregular_volume_subset::IID();
    }

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
    virtual bool load_irregular_volume_from_cache(
        std::ifstream                                   &ifile,
        nv::index::IIrregular_volume_subset*            irregular_volume_subset) const;

    virtual bool save_irregular_volume_to_cache(
        std::ofstream                                           &ofile,
        const nv::index::IIrregular_volume_subset::Mesh_storage &mesh_storage,
        const nv::index::IIrregular_volume_subset*              irregular_volume_subset) const;

    /// Irregular volume's data filename
    std::string                             m_filename;
    
    /// configuration information for the session exporter
    std::string                             m_configuration;
    
    /// use cache file
    bool                                    m_use_cache;

    /// local caching of input file
    struct IVOL_ts_mesh
    {
        std::vector<mi::math::Vector<mi::Float32, 3> >      vertices;
        mi::math::Bbox<mi::Float32, 3>                      bbox;

        std::vector<mi::math::Vector<mi::Uint32, 4> >       tetrahedron_indices;

        std::vector<mi::Float32>                            attributes;
        mi::math::Vector<mi::Float32, 2>                    attributes_value_range;

        IVOL_ts_mesh();
    };

    static mi::base::Lock                   m_cached_input_lock;
    static std::string                      m_cached_input_file;
    static IVOL_ts_mesh                     m_cached_input_ivol_mesh;

};

}} // namespace nv::index_common

#endif // NVIDIA_INDEX_BIN_IRREGULAR_VOLUME_IMPORTER_VTK_H

