/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief external irregular volume data import job. It loads a subset of an irregular volume 
//  contained in an given voxel (returned in memory allocated in execute method).

#ifndef NVIDIA_INDEX_BIN_IRREGULAR_VOLUME_IMPORTER_H
#define NVIDIA_INDEX_BIN_IRREGULAR_VOLUME_IMPORTER_H

#include <string>
#include <vector>

#include <mi/dice.h>

#include <nv/index/idistributed_data_import_callback.h>
#include <nv/index/iirregular_volume_subset.h>

namespace nv {
namespace index_common {

/// EnSight Gold file utilities
namespace esg
{

/// Geometry element type enum
enum Type_element
{
    TYPE_ELEMENT_INVALID = 0,
    TYPE_ELEMENT_TETRAHEDRON,
    TYPE_ELEMENT_HEXAHEDRON,
};

/// Geometry information of irregular volume
struct Geometry
{
    Geometry();
    
    void rebuild_subset_indices();
    
    Type_element    type;
    mi::Uint32      part_id;
    mi::Uint32      nb_vertices;
    mi::Uint32      nb_cells;
    
    std::vector<mi::math::Vector_struct<mi::Float32, 3> >   vertices;
    std::vector<mi::math::Vector_struct<mi::Float32, 3> >   displ;
    std::vector<mi::Uint32>                                 cell_vtx_ids;
    std::vector<mi::Float32>                                scalar;
};

/// list of Geometry type
typedef std::vector<Geometry> Geometry_list;

/// file type detection and read. It automatically detects if it's binary or ASCII
void read_case_from_file(
    const mi::math::Bbox_struct<mi::Float32, 3>& bbox,
    const std::string&	                         geometry_file,
    const std::string&	                         scalar_file,
    Geometry_list&                               geometry_list);
    
/// file type detection and read. It automatically detects if it's binary or ASCII
void convert_case_to_ts(
    const mi::math::Bbox_struct<mi::Float32, 3>& bbox,
    const std::string&	                         geometry_file,
    const std::string&	                         scalar_file,
    const std::string&	                         displ_file,
    const std::string&	                         ts_file,
    mi::Uint32      	                         nb_frames,
    Geometry_list&                               geometry_list);

/// read geometry: ascii
void read_geometry_from_ascii_file(
    const mi::math::Bbox_struct<mi::Float32, 3>& bbox,
    const std::string&	                         geometry_file,
    Geometry_list&                               geometry_list);

/// read geometry: binary
void read_geometry_from_binary_file(
    const mi::math::Bbox_struct<mi::Float32, 3>& bbox,
    const std::string&	                         geometry_file,
    Geometry_list&                               geometry_list);

/// read scalar: ascii
void read_scalar_from_ascii_file(
    const std::string&	scalar_file,
    Geometry_list &geometry_list);

/// read scalar: binary
void read_scalar_from_binary_file(
    const std::string&	scalar_file,
    Geometry_list &geometry_list);
    
/// read displacement data: ascii
void read_displacement_from_ascii_file(
    const std::string&	displ_file,
    Geometry_list &geometry_list);
    
}

/// External irregular volume importer implementation
class Irregular_volume_importer :
    public nv::index::Distributed_continuous_data_import_callback<0xa41773dd,0x4c9b,0x467d,0xa2,0x73,0x0,0xea,0xb0,0xff,0x98,0x9b>
{
public:
    /// constructor
    Irregular_volume_importer(const std::string& filename, bool use_cache = false);

    /// Temporary constructor for EnSight geometry and scalar file.
    Irregular_volume_importer(
        const std::string&                      geometry_file, 
        const std::string&                      scalar_file,
        bool                                    use_cache = false);
        
    /// Temporary constructor for EnSight to export to ts.
    Irregular_volume_importer(
        const std::string&                      geometry_file, 
        const std::string&                      scalar_file,
        const std::string&                      displ_file,
        const std::string&                      ts_file,
        mi::Uint32                              nb_frames,
        bool                                    use_cache = false);

    /// default constructor
    Irregular_volume_importer();

    /// destructor
    virtual ~Irregular_volume_importer();
    
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

    virtual bool load_irregular_volume_from_ts(
        const mi::math::Bbox_struct<mi::Float32, 3>&    bounding_box,
        nv::index::IIrregular_volume_subset*            irregular_volume_subset) const;

    virtual bool load_irregular_volume_from_ensight_gold(
        const mi::math::Bbox_struct<mi::Float32, 3>&    bounding_box,
        nv::index::IIrregular_volume_subset*            irregular_volume_subset) const;
        
    virtual bool load_irregular_volume_from_cache(
        std::ifstream                                   &ifile,
        nv::index::IIrregular_volume_subset*            irregular_volume_subset) const;

    virtual bool save_irregular_volume_to_cache(
        std::ofstream                                           &ofile,
        const nv::index::IIrregular_volume_subset::Mesh_storage &mesh_storage,
        const nv::index::IIrregular_volume_subset*              irregular_volume_subset) const;

    /// Irregular volume's data filename
    std::string                             m_filename;

    /// Temporary filename for Ensight scalar files
    std::string                             m_scalar_filename;

    /// Temporary filename for Ensight displacement files
    std::string                             m_displ_filename;

    /// Temporary filename for Ensight displacement files
    std::string                             m_ts_filename;
    
    /// Temporary frame number for Ensight dataset with time varying meshes
    mi::Uint32                             m_nb_frames;
    
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

#endif // NVIDIA_INDEX_BIN_IRREGULAR_VOLUME_IMPORTER_H

