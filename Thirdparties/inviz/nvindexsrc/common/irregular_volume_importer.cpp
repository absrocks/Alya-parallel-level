/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief external irregular volume data importer job

#include "irregular_volume_importer.h"
#include <iostream>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <map>
#include <typeinfo>
#include <vector>
#include <nv/index/iirregular_volume_subset.h>
#include "forwarding_logger.h"
#include "../alya/alyaglobal.h"

namespace nv {
namespace index_common {

//----------------------------------------------------------------------
Irregular_volume_importer::IVOL_ts_mesh::IVOL_ts_mesh()
  : bbox(mi::math::Vector<mi::Float32, 3>( (std::numeric_limits<float>::max)()),
         mi::math::Vector<mi::Float32, 3>(-(std::numeric_limits<float>::max)()))
  , attributes_value_range( (std::numeric_limits<float>::max)(),
                           -(std::numeric_limits<float>::max)())
{
}

mi::base::Lock                          Irregular_volume_importer::m_cached_input_lock;
std::string                             Irregular_volume_importer::m_cached_input_file;
Irregular_volume_importer::IVOL_ts_mesh Irregular_volume_importer::m_cached_input_ivol_mesh;

//----------------------------------------------------------------------
Irregular_volume_importer::Irregular_volume_importer(
    const std::string&  filename,
    bool                use_cache)
{
    m_use_cache = use_cache;

    m_filename = filename;

    m_nb_frames = 1;
    m_configuration = std::string() +
        "importer=raw\n" +
        "input_file=" + m_filename + "\n";
}

Irregular_volume_importer::Irregular_volume_importer(
    const std::string&                          geometry_file, 
    const std::string&                          scalar_file,
    bool                                        use_cache)
{
    m_use_cache = use_cache;
    
    m_filename = geometry_file;
    m_scalar_filename = scalar_file;

    m_nb_frames = 1;
    
    m_configuration = std::string() +
        "importer=raw\n" +
        "input_file=" + geometry_file + "\n" +
        "input_scalar_file=" + scalar_file + "\n";
}

Irregular_volume_importer::Irregular_volume_importer(
        const std::string&                      geometry_file, 
        const std::string&                      scalar_file,
        const std::string&                      displ_file,
        const std::string&                      ts_file,
        mi::Uint32                              nb_frames,
        bool                                    use_cache)
{
    m_use_cache = use_cache;
    
    m_filename = geometry_file;
    m_scalar_filename = scalar_file;
    m_displ_filename = displ_file;
    m_ts_filename = ts_file;

    m_nb_frames = nb_frames;
    
    m_configuration = std::string() +
        "importer=raw\n" +
        "input_file=" + geometry_file + "\n" +
        "input_scalar_file=" + scalar_file + "\n" + 
        "input_displ_file=" + displ_file + "\n";
}


//----------------------------------------------------------------------
// constructor
Irregular_volume_importer::Irregular_volume_importer()
{
    m_use_cache = false;
}

//----------------------------------------------------------------------
// destructor
Irregular_volume_importer::~Irregular_volume_importer()
{
}

//----------------------------------------------------------------------
mi::Size Irregular_volume_importer::estimate(
    const mi::math::Bbox_struct<mi::Float32, 3>&    bounding_box,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{
    return 0;
}

//----------------------------------------------------------------------
nv::index::IDistributed_data_subset* Irregular_volume_importer::create(
    const mi::math::Bbox_struct<mi::Float32, 3>&    bounding_box,
    nv::index::IData_subset_factory*                factory,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{
    INFO_LOG << "Irregular volume importer loads '" << m_filename << "', bounds: " << bounding_box;

    const std::string file_extension =
          m_filename.find_last_of(".") != std::string::npos
        ? m_filename.substr(m_filename.find_last_of(".") + 1)
        : "";
      
    // TODO: Temporary using 'geom' and 'scalar' instead of 'case' for EnSight files 
    if(     file_extension == ""
       || !(file_extension == "ts" || file_extension == "geom"))
    {
        ERROR_LOG << "Unknown file type: '" << m_filename << "'.";
        return NULL;
    }
    
    mi::base::Handle<nv::index::IIrregular_volume_subset>
        irregular_volume_subset(factory->create<nv::index::IIrregular_volume_subset>());

    if (!irregular_volume_subset.is_valid_interface()) {
        ERROR_LOG << "Cannot create an irregular volume subset.";
        return NULL;
    }
                    
    if(file_extension == "ts") {
        if (!load_irregular_volume_from_ts(bounding_box, irregular_volume_subset.get())) {
            return NULL;
        }
    }
    // TODO: Temporary using 'geom' and 'scalar' instead of 'case' for EnSight files 
    else if(file_extension == "geom") {
        if (!load_irregular_volume_from_ensight_gold(bounding_box, irregular_volume_subset.get())) {
            return NULL;
        }
    }
                    
    irregular_volume_subset->retain();
    return irregular_volume_subset.get();
}

//----------------------------------------------------------------------
const char* Irregular_volume_importer::get_configuration() const
{
    return m_configuration.c_str();
}

//----------------------------------------------------------------------
void Irregular_volume_importer::serialize(
    mi::neuraylib::ISerializer * serializer) const
{
    serializer->write(&m_use_cache, 1);
    
    mi::Uint32 nb_elements = mi::Uint32(m_filename.size());
    serializer->write(&nb_elements, 1);
    serializer->write(reinterpret_cast<const mi::Uint8*>(m_filename.c_str()), nb_elements);
    
    nb_elements = mi::Uint32(m_scalar_filename.size());
    serializer->write(&nb_elements, 1);
    serializer->write(reinterpret_cast<const mi::Uint8*>(m_scalar_filename.c_str()), nb_elements);
    
}

//----------------------------------------------------------------------
void Irregular_volume_importer::deserialize(
    mi::neuraylib::IDeserializer * deserializer)
{
    deserializer->read(&m_use_cache, 1);
    
    mi::Uint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_filename.resize(nb_elements);
    deserializer->read(reinterpret_cast<mi::Uint8*>(&m_filename[0]), nb_elements);
    
    nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_scalar_filename.resize(nb_elements);
    deserializer->read(reinterpret_cast<mi::Uint8*>(&m_scalar_filename[0]), nb_elements);
}

//----------------------------------------------------------------------
bool Irregular_volume_importer::load_irregular_volume_from_ts(
    const mi::math::Bbox_struct<mi::Float32, 3>&    bounding_box,
    nv::index::IIrregular_volume_subset*            irregular_volume_subset) const
{
    // select the attribute type
    typedef mi::Float32 IVOL_attrib_type;
    using mi::math::Vector;
    using mi::math::Vector_struct;

    // typedef Vector<mi::Uint32, 3>           Vec3ui;
    // typedef Vector_struct<mi::Uint32, 3>    Vec3ui_struct;
    typedef Vector<mi::Uint32, 4>           Vec4ui;
    // typedef Vector_struct<mi::Uint32, 4>    Vec4ui_struct;

    typedef Vector<mi::Float32, 3>          Vec3f;
    // typedef Vector_struct<mi::Float32, 3>   Vec3f_struct;
    // typedef Vector<mi::Float32, 4>          Vec4f;
    // typedef Vector_struct<mi::Float32, 4>   Vec4f_struct;

    std::ostringstream cache_file_name_str;
    std::ifstream      ifile_cache;

    if (m_use_cache) {
        cache_file_name_str << m_filename << "_"
                            << bounding_box.min.x << "." << bounding_box.max.x << "x"
                            << bounding_box.min.y << "." << bounding_box.max.y << "x"
                            << bounding_box.min.z << "." << bounding_box.max.z
                            << ".ivc";
                            
        ifile_cache.open(cache_file_name_str.str().c_str(), std::ios_base::in | std::ios_base::binary);
    }

    // read irregular volume subset from cache file
    if (m_use_cache && ifile_cache) {
        return load_irregular_volume_from_cache(ifile_cache, irregular_volume_subset);
    }
    else {
        const IVOL_ts_mesh* ivol_mesh = 0;

        {
            mi::base::Lock::Block block(&m_cached_input_lock);

            if (m_cached_input_file != m_filename) {
                std::ifstream ifile(m_filename.c_str(), std::ios_base::in | std::ios_base::binary);

                if (!ifile.is_open()) {
                    ERROR_LOG << "Irregular volume importer (.ts file): error opening input file '" << m_filename << "'.";
                    return false;
                }
        
                mi::Uint32  nb_global_vertices     = 0u;
                mi::Uint32  nb_global_tetrahedrons = 0u;

                ifile >> nb_global_vertices;
                ifile >> nb_global_tetrahedrons;
        
                if (nb_global_vertices == 0 || nb_global_tetrahedrons == 0) {
                    ERROR_LOG << "Irregular volume importer (.ts file): "
                              << "invalid mesh data read (invalid number of vertices or tetrahedrons) "
                              << "file: '" << m_filename << "'";
                    return false;
                }

                // reset the cached mesh
                m_cached_input_ivol_mesh = IVOL_ts_mesh();

                // read tetrahedron's mesh vertices and attributes
                m_cached_input_ivol_mesh.vertices.resize(nb_global_vertices);
                m_cached_input_ivol_mesh.attributes.resize(nb_global_vertices);
                m_cached_input_ivol_mesh.tetrahedron_indices.resize(nb_global_tetrahedrons);

                // read vertices
                mi::math::Bbox<mi::Float32, 3> dataset_extent;

                for (mi::Uint32 v = 0u; v < nb_global_vertices; ++v) {
                    Vec3f&       cur_vert   = m_cached_input_ivol_mesh.vertices[v];
                    mi::Float32& cur_attrib = m_cached_input_ivol_mesh.attributes[v];

                    ifile >> cur_vert.x;
                    ifile >> cur_vert.y;
                    ifile >> cur_vert.z;
                    ifile >> cur_attrib;
		    cur_attrib = v ;
                    dataset_extent.min.x = mi::math::min(dataset_extent.min.x, cur_vert.x);
                    dataset_extent.max.x = mi::math::max(dataset_extent.max.x, cur_vert.x);
                    dataset_extent.min.y = mi::math::min(dataset_extent.min.y, cur_vert.y);
                    dataset_extent.max.y = mi::math::max(dataset_extent.max.y, cur_vert.y);
                    dataset_extent.min.z = mi::math::min(dataset_extent.min.z, cur_vert.z);
                    dataset_extent.max.z = mi::math::max(dataset_extent.max.z, cur_vert.z);

                    m_cached_input_ivol_mesh.bbox.min = mi::math::elementwise_min(m_cached_input_ivol_mesh.bbox.min, cur_vert);
                    m_cached_input_ivol_mesh.bbox.max = mi::math::elementwise_max(m_cached_input_ivol_mesh.bbox.max, cur_vert);

                    m_cached_input_ivol_mesh.attributes_value_range.x =
                        mi::math::min(m_cached_input_ivol_mesh.attributes_value_range.x, cur_attrib);
                    m_cached_input_ivol_mesh.attributes_value_range.y =
                        mi::math::max(m_cached_input_ivol_mesh.attributes_value_range.y, cur_attrib);
                }

                // ERROR_LOG << dataset_extent;

                // read tetrahedron indices
                for (mi::Uint32 t = 0u; t < nb_global_tetrahedrons; ++t) {
                    ifile >> m_cached_input_ivol_mesh.tetrahedron_indices[t].x;
                    ifile >> m_cached_input_ivol_mesh.tetrahedron_indices[t].y;
                    ifile >> m_cached_input_ivol_mesh.tetrahedron_indices[t].z;
                    ifile >> m_cached_input_ivol_mesh.tetrahedron_indices[t].w;
                }

                ifile.close();
            }

            m_cached_input_file = m_filename;
            ivol_mesh           = &m_cached_input_ivol_mesh; 
        }

        if (false)
        {
            mi::Float32 range_min = (std::numeric_limits<float>::max)();
            mi::Float32 range_max = -(std::numeric_limits<float>::max)();
            {
                std::ifstream ifile("value_range.txt", std::ios_base::in);
                ifile >> range_min;
                ifile >> range_max;
                WARN_LOG << "Read value range: " << range_min << ", " << range_max;
            }

            {
                WARN_LOG << "Current value range: " << m_cached_input_ivol_mesh.attributes_value_range.x
                    << ", " << m_cached_input_ivol_mesh.attributes_value_range.y;

                range_min = mi::math::min(m_cached_input_ivol_mesh.attributes_value_range.x, range_min);
                range_max = mi::math::max(m_cached_input_ivol_mesh.attributes_value_range.y, range_max);

                std::ofstream ofile("value_range.txt", std::ios_base::out | std::ios_base::trunc);
                ofile << range_min << " " << range_max;
                WARN_LOG << "Written value range: " << range_min << ", " << range_max;
            }
        }
        
        // read tetrahedrons and collect only those that intersects the subset bounding box
        mi::math::Bbox<mi::Float32, 3>      subset_bbox(static_cast<mi::Float32>(bounding_box.min.x),
                                                        static_cast<mi::Float32>(bounding_box.min.y),
                                                        static_cast<mi::Float32>(bounding_box.min.z),
                                                        static_cast<mi::Float32>(bounding_box.max.x),
                                                        static_cast<mi::Float32>(bounding_box.max.y),
                                                        static_cast<mi::Float32>(bounding_box.max.z));

        mi::Uint32                          nb_subset_vertices     = 0u;
        mi::Uint32                          nb_subset_tetrahedrons = 0u;

        std::vector<Vec4ui>                 subset_tetrahedrons;
        std::vector<Vec3f>                  subset_vertices;
        std::vector<IVOL_attrib_type>       subset_attributes;
        
        // generate the irregular volume mesh subset (all tetrahedra contained in the subset-bbox)
        // * also generate the global max edge length required by the renderer

        const std::vector<Vec4ui>&      global_tet_indices = ivol_mesh->tetrahedron_indices;
        const std::vector<Vec3f>&       global_vertices    = ivol_mesh->vertices;
        const std::vector<mi::Float32>& global_attribs     = ivol_mesh->attributes;

        mi::Float32 max_edge_len_sqr = -(std::numeric_limits<mi::Float32>::max)();
        for (mi::Size t = 0u; t < ivol_mesh->tetrahedron_indices.size(); ++t) {
            const Vec4ui&                   tet_vtx_indices = global_tet_indices[t];
            mi::math::Bbox<mi::Float32, 3>  tet_bbox;

            const Vec3f a = Vec3f(global_vertices[tet_vtx_indices.x]);
            const Vec3f b = Vec3f(global_vertices[tet_vtx_indices.y]);
            const Vec3f c = Vec3f(global_vertices[tet_vtx_indices.z]);
            const Vec3f d = Vec3f(global_vertices[tet_vtx_indices.w]);

            // TODO: use triangle-box intersection for all four faces, much more precise
            tet_bbox.clear();
            tet_bbox.insert(a);
            tet_bbox.insert(b);
            tet_bbox.insert(c);
            tet_bbox.insert(d);
            
            if (tet_bbox.intersects(subset_bbox)) {
                subset_tetrahedrons.push_back(tet_vtx_indices);
            }

            // calculate max edge length
            const Vec3f e0 = b - a;
            const Vec3f e1 = c - a;
            const Vec3f e2 = d - a;
            const Vec3f e3 = c - b;
            const Vec3f e4 = d - b;
            const Vec3f e5 = d - c;

            max_edge_len_sqr = mi::math::max(max_edge_len_sqr, mi::math::square_length(e0));
            max_edge_len_sqr = mi::math::max(max_edge_len_sqr, mi::math::square_length(e1));
            max_edge_len_sqr = mi::math::max(max_edge_len_sqr, mi::math::square_length(e2));
            max_edge_len_sqr = mi::math::max(max_edge_len_sqr, mi::math::square_length(e3));
            max_edge_len_sqr = mi::math::max(max_edge_len_sqr, mi::math::square_length(e4));
            max_edge_len_sqr = mi::math::max(max_edge_len_sqr, mi::math::square_length(e5));
        }
        
        // build subset vertex list and remap indices
        nb_subset_tetrahedrons = subset_tetrahedrons.size();

	//SImulation in-situ registry
	int simboxid = registry_sim.bbox_exist(bounding_box.min.x,bounding_box.min.y,bounding_box.min.z,bounding_box.max.x,bounding_box.max.y,bounding_box.max.z);
	
	if (simboxid == -1)
	  {
	    simboxid = registry_sim.add_bbox(bounding_box.min.x,bounding_box.min.y,bounding_box.min.z,bounding_box.max.x,bounding_box.max.y,bounding_box.max.z); 
	  }
	
        for (mi::Uint32 t = 0u; t < nb_subset_tetrahedrons; ++t) {
            Vec4ui& tet_vtx_indices = subset_tetrahedrons[t];
            for(mi::Uint32 j = 0; j < 4u; ++j) {
                mi::Uint32& vtx_index = tet_vtx_indices[j];
		
                std::map<unsigned int,unsigned int>::const_iterator kt = registry_sim.localfind(simboxid,vtx_index);
		if(kt != registry_sim.localend(simboxid) ) {
		  vtx_index = kt->second;
                }
                else {
		  const mi::Uint32 new_vtx_idx = subset_vertices.size();
		  registry_sim.localset(simboxid,vtx_index,new_vtx_idx);
		  registry_sim.globalset(simboxid,new_vtx_idx,vtx_index);
		  subset_vertices.push_back(global_vertices[vtx_index]);
		  subset_attributes.push_back(global_attribs[vtx_index]);
		  vtx_index = new_vtx_idx;
		}
	    }
	}
	
        nb_subset_vertices                 = subset_vertices.size();
        mi::Uint32 nb_cells                = nb_subset_tetrahedrons;
        mi::Uint32 nb_cell_face_indices    = nb_cells * 4u;

        mi::Uint32 nb_faces                = nb_cell_face_indices;
        mi::Uint32 nb_face_vtx_indices     = nb_faces * 3u;

        nv::index::IIrregular_volume_subset::Mesh_parameters mesh_params;

        // general mesh geometry and topology info
        mesh_params.nb_vertices            = nb_subset_vertices;
        mesh_params.nb_face_vtx_indices    = nb_face_vtx_indices;
        mesh_params.nb_faces               = nb_faces;
        mesh_params.nb_cell_face_indices   = nb_cell_face_indices;
        mesh_params.nb_cells               = nb_cells;

        // required mesh information for the renderer
        mesh_params.global_max_edge_length = mi::math::sqrt(max_edge_len_sqr);

        nv::index::IIrregular_volume_subset::Mesh_storage mesh_storage;
        if (!irregular_volume_subset->generate_mesh_storage(mesh_params, mesh_storage)) {
            ERROR_LOG << "Irregular volume importer (.ts file): "
                      << "unable to generate irregular volume mesh storage, "
                      << "file: '" << m_filename << "'";
            return false;
        }

        if (mesh_params.nb_vertices > 0) {
            nv::index::IIrregular_volume_subset::Attribute_parameters attrib_params;
            attrib_params.affiliation        = nv::index::IIrregular_volume_subset::ATTRIB_AFFIL_PER_VERTEX;
            attrib_params.type               = nv::index::IIrregular_volume_subset::ATTRIB_TYPE_FLOAT32;
            attrib_params.nb_attrib_values   = nb_subset_vertices;

            nv::index::IIrregular_volume_subset::Attribute_storage attribute_storage;
            if (!irregular_volume_subset->generate_attribute_storage(0u, attrib_params, attribute_storage)) {
                ERROR_LOG << "Irregular volume importer (.ts file): "
                          << "unable to generate irregular volume attribute storage, "
                          << "file: '" << m_filename << "'";
                return false;
            }

            IVOL_attrib_type* subset_attrib_values = reinterpret_cast<IVOL_attrib_type*>(attribute_storage.attrib_values);
        
            // copy vertices and attributes
            for (mi::Uint32 v = 0u; v < nb_subset_vertices; ++v) {
                mesh_storage.vertices[v] = subset_vertices[v];
                subset_attrib_values[v]  = subset_attributes[v];
            }
        
            mi::Uint32 next_vidx = 0u;
            mi::Uint32 next_fidx = 0u;

            // generate cells, cell's face indices, faces, face's vertex indices 
            for (mi::Uint32 t = 0u; t < nb_subset_tetrahedrons; ++t) {
                const mi::Uint32 a   = subset_tetrahedrons[t].x;
                const mi::Uint32 b   = subset_tetrahedrons[t].y;
                const mi::Uint32 c   = subset_tetrahedrons[t].z;
                const mi::Uint32 d   = subset_tetrahedrons[t].w;

                const Vec3f& av      = subset_vertices[a];
                const Vec3f& bv      = subset_vertices[b];
                const Vec3f& cv      = subset_vertices[c];
                const Vec3f& dv      = subset_vertices[d];

                const Vec3f centroid = (av + bv + cv + dv) * 0.25f;
            
                // TODO: when available, use C++11 lambda
                // * adds a face to the ivol mesh storage
                // * tries to orient faces to have correct CCW vertex ordering
#define IVOL_ADD_TET_FACE(i0, i1, i2, p0, p1, p2)                                   \
                do {                                                                \
                    const Vec3f e1 = ((p1) - (p0));                                 \
                    const Vec3f e2 = ((p2) - (p0));                                 \
                                                                                    \
                    /* face plane */                                                \
                    Vec3f       n = cross(e1, e2);                                  \
                    n.normalize();                                                  \
                    const mi::Float32 dst = -(dot(n, (p0)));                        \
                                                                                    \
                    /* tetrahedron centroid distance to face plane */               \
                    const mi::Float32 cd = dot(n, centroid) + dst;                  \
                                                                                    \
                    if (cd < 0.0f) {                                                \
                        /* correct ordering */                                      \
                        mesh_storage.face_vtx_indices[next_vidx + 0u] = (i0);       \
                        mesh_storage.face_vtx_indices[next_vidx + 1u] = (i1);       \
                        mesh_storage.face_vtx_indices[next_vidx + 2u] = (i2);       \
                    }                                                               \
                    else {                                                          \
                        /* invert ordering */                                       \
                        mesh_storage.face_vtx_indices[next_vidx + 0u] = (i0);       \
                        mesh_storage.face_vtx_indices[next_vidx + 1u] = (i2);       \
                        mesh_storage.face_vtx_indices[next_vidx + 2u] = (i1);       \
                    }                                                               \
                                                                                    \
                    mesh_storage.faces[next_fidx].nb_vertices        = 3u;          \
                    mesh_storage.faces[next_fidx].start_vertex_index = next_vidx;   \
                                                                                    \
                    next_fidx += 1;                                                 \
                    next_vidx += 3;                                                 \
                } while (0)

                // cell's face indices
                // * this might seem redundant here, but consider datasets with shared faces among cells
                mesh_storage.cell_face_indices[next_fidx + 0] = next_fidx + 0;
                mesh_storage.cell_face_indices[next_fidx + 1] = next_fidx + 1;
                mesh_storage.cell_face_indices[next_fidx + 2] = next_fidx + 2;
                mesh_storage.cell_face_indices[next_fidx + 3] = next_fidx + 3;

                // fill cell data
                mesh_storage.cells[t].nb_faces         = 4u;
                mesh_storage.cells[t].start_face_index = next_fidx;

                // needs to happen down here as it increases next_fidx and next_vidx
                IVOL_ADD_TET_FACE(a, c, d, av, cv, dv);
                IVOL_ADD_TET_FACE(a, b, c, av, bv, cv);
                IVOL_ADD_TET_FACE(a, d, b, av, dv, bv);
                IVOL_ADD_TET_FACE(b, d, c, bv, dv, cv);

#undef IVOL_ADD_TET_FACE

                // TODO: determine max edge length
                // * needs to be a global measure
                //{
                //    ivol_mesh->m_max_edge_length = max(length(bv - av), ivol_mesh->m_max_edge_length);
                //    ivol_mesh->m_max_edge_length = max(length(cv - av), ivol_mesh->m_max_edge_length);
                //    ivol_mesh->m_max_edge_length = max(length(dv - av), ivol_mesh->m_max_edge_length);

                //    ivol_mesh->m_max_edge_length = max(length(cv - bv), ivol_mesh->m_max_edge_length);
                //    ivol_mesh->m_max_edge_length = max(length(dv - bv), ivol_mesh->m_max_edge_length);

                //    ivol_mesh->m_max_edge_length = max(length(dv - cv), ivol_mesh->m_max_edge_length);
                //}
            }
        }

        // save irregular volume subset to cache file
        if (m_use_cache) {
            std::ofstream ofile_cache(cache_file_name_str.str().c_str(), std::ios_base::out | std::ios_base::binary);
            if (ofile_cache) {
                save_irregular_volume_to_cache(ofile_cache, mesh_storage, irregular_volume_subset);
            }
        }
    }
    
    return true;
}

//----------------------------------------------------------------------
bool Irregular_volume_importer::load_irregular_volume_from_ensight_gold(
    const mi::math::Bbox_struct<mi::Float32, 3>&    bounding_box,
    nv::index::IIrregular_volume_subset*            irregular_volume_subset) const
{
    if(m_scalar_filename == "")
    {
        ERROR_LOG << "Scalar file name is empty.";
        return false;
    }

    std::ostringstream cache_file_name_str;
    std::ifstream ifile_cache;
    if(m_use_cache)
    {
        cache_file_name_str << m_filename << "_"
                            << bounding_box.min.x << "." << bounding_box.max.x << "x"
                            << bounding_box.min.y << "." << bounding_box.max.y << "x"
                            << bounding_box.min.z << "." << bounding_box.max.z
                            << ".ivc";
                            
        ifile_cache.open(cache_file_name_str.str().c_str(), std::ios_base::in | std::ios_base::binary);
    }

    // read irregular volume subset from cache file
    if(m_use_cache && ifile_cache)
    {
        return load_irregular_volume_from_cache(ifile_cache, irregular_volume_subset);
    }
    else
    {
        /// read geometry from Ensight file
        esg::Geometry_list geometry_list;
        
        if(m_nb_frames > 1)
            esg::convert_case_to_ts(bounding_box, m_filename, m_scalar_filename, m_displ_filename, m_ts_filename, m_nb_frames, geometry_list);
        else
            esg::read_case_from_file(bounding_box, m_filename, m_scalar_filename, geometry_list);
        
        /// rebuild indices and collapse all geometry to only one.
        mi::Uint32 nb_geometries = geometry_list.size();
        mi::Uint32 nb_vertices = 0u;
        mi::Uint32 nb_faces = 0u;
        mi::Uint32 nb_face_vtx_indices = 0u;
        mi::Uint32 nb_cells = 0u;
        mi::Uint32 nb_cell_face_indices = 0u;
        
        std::vector<mi::math::Vector_struct<mi::Float32, 3> >   vertices;
        std::vector<mi::Uint32>                                 cell_vtx_ids;
        std::vector<mi::Float32>                                scalar;
        
        esg::Type_element type_element = esg::TYPE_ELEMENT_INVALID;
        
        for(mi::Uint32 i=0; i<nb_geometries; ++i)
        {
            // skip invalid geometry
            if(geometry_list[i].type == esg::TYPE_ELEMENT_INVALID)
                continue;
                
            type_element = geometry_list[i].type;

            // rebuild cell vertex indices
            geometry_list[i].rebuild_subset_indices();

            // remap cell vertex indices to global list
            if(nb_vertices > 0)
            {
                std::vector<mi::Uint32>::iterator jt = geometry_list[i].cell_vtx_ids.begin();
                for(; jt!= geometry_list[i].cell_vtx_ids.end(); ++jt)
                    *jt += nb_vertices;
            }
            
            // collapse data
            vertices.insert(vertices.end(), geometry_list[i].vertices.begin(), geometry_list[i].vertices.end());
            geometry_list[i].vertices.clear();
            
            cell_vtx_ids.insert(cell_vtx_ids.end(), geometry_list[i].cell_vtx_ids.begin(), geometry_list[i].cell_vtx_ids.end());
            geometry_list[i].cell_vtx_ids.clear();
            
            scalar.insert(scalar.end(), geometry_list[i].scalar.begin(), geometry_list[i].scalar.end());
            geometry_list[i].scalar.clear();
            
            nb_vertices += geometry_list[i].nb_vertices;
            nb_cells += geometry_list[i].nb_cells;
            
            switch(geometry_list[i].type)
            {
                case esg::TYPE_ELEMENT_TETRAHEDRON:
                    nb_faces += geometry_list[i].nb_cells*4u;
                    nb_face_vtx_indices += geometry_list[i].nb_cells*12u;
                    break;

                case esg::TYPE_ELEMENT_HEXAHEDRON:
                    nb_faces += geometry_list[i].nb_cells*6u;
                    nb_face_vtx_indices += geometry_list[i].nb_cells*24u;
                    break;
                    
                default:
                    break;
            }
                
            
        }
        
        nb_cell_face_indices = nb_faces;
        
        nv::index::IIrregular_volume_subset::Mesh_parameters mesh_params;
        mesh_params.nb_vertices          = nb_vertices;
        mesh_params.nb_face_vtx_indices  = nb_face_vtx_indices;
        mesh_params.nb_faces             = nb_faces;
        mesh_params.nb_cell_face_indices = nb_cell_face_indices;
        mesh_params.nb_cells             = nb_cells;

        nv::index::IIrregular_volume_subset::Mesh_storage mesh_storage;
        if(!irregular_volume_subset->generate_mesh_storage(mesh_params, mesh_storage))
            return false;
            
        nv::index::IIrregular_volume_subset::Attribute_parameters attrib_params;
        attrib_params.affiliation      = nv::index::IIrregular_volume_subset::ATTRIB_AFFIL_PER_VERTEX;
        attrib_params.type             = nv::index::IIrregular_volume_subset::ATTRIB_TYPE_UINT8;
        attrib_params.nb_attrib_values = nb_vertices;


        nv::index::IIrregular_volume_subset::Attribute_storage attribute_storage;
        if(!irregular_volume_subset->generate_attribute_storage(0u, attrib_params, attribute_storage))
            return false;

        memcpy(mesh_storage.vertices, &vertices[0], sizeof(mi::Float32)*3*nb_vertices);
        memcpy(attribute_storage.attrib_values, &scalar[0], sizeof(mi::Float32)*nb_vertices);
                
        // read cells, cell's face indices, faces, face's vertex indices 
        mi::Uint32 faces_per_cell = 0u;
        mi::Uint32 vertices_per_face = 0u;
        mi::Uint32 vertices_per_cell = 0u;
        
        switch(type_element)
        {
            case esg::TYPE_ELEMENT_TETRAHEDRON:
                faces_per_cell = 4u;
                vertices_per_face = 3u;
                vertices_per_cell = 4u;
                break;
                
            case esg::TYPE_ELEMENT_HEXAHEDRON:
                faces_per_cell = 6u;
                vertices_per_face = 4u;
                vertices_per_cell = 8u;
                break;
                
            default:
                break;
        }
        
        
        for(mi::Uint32 t=0u; t < nb_cells; ++t) 
        {
            // fill cell data
            mesh_storage.cells[t].nb_faces = faces_per_cell;
            mesh_storage.cells[t].start_face_index = t*faces_per_cell;
            
            // fill cell's face data
            for(mi::Uint32 i=0u; i<faces_per_cell; i++)
            {
                const mi::Uint32 cur_face_idx = t*faces_per_cell + i;
                
                // cell's face indices
                mesh_storage.cell_face_indices[cur_face_idx] = cur_face_idx;
                
                // faces
                mesh_storage.faces[cur_face_idx].nb_vertices = vertices_per_face;
                mesh_storage.faces[cur_face_idx].start_vertex_index = cur_face_idx*vertices_per_face;
            }
            
            // face's vertex indices
            mi::Uint32 *vp = &cell_vtx_ids[t*vertices_per_cell];
            mi::Uint32* fvi = &mesh_storage.face_vtx_indices[t*faces_per_cell*vertices_per_face];
            
            switch(type_element)
            {
                case esg::TYPE_ELEMENT_TETRAHEDRON:
                    fvi[ 0] = vp[0];    fvi[ 1] = vp[1];    fvi[ 2] = vp[3];
                    fvi[ 3] = vp[1];    fvi[ 4] = vp[2];    fvi[ 5] = vp[3];
                    fvi[ 6] = vp[2];    fvi[ 7] = vp[0];    fvi[ 8] = vp[3];
                    fvi[ 9] = vp[0];    fvi[10] = vp[2];    fvi[11] = vp[1];
                    break;

                case esg::TYPE_ELEMENT_HEXAHEDRON:
                    fvi[ 0] = vp[0];    fvi[ 1] = vp[1];    fvi[ 2] = vp[5];    fvi[ 3] = vp[4];
                    fvi[ 4] = vp[1];    fvi[ 5] = vp[2];    fvi[ 6] = vp[6];    fvi[ 7] = vp[5];
                    fvi[ 8] = vp[2];    fvi[ 9] = vp[3];    fvi[10] = vp[7];    fvi[11] = vp[6];
                    fvi[12] = vp[3];    fvi[13] = vp[0];    fvi[14] = vp[4];    fvi[15] = vp[7];
                    fvi[16] = vp[4];    fvi[17] = vp[5];    fvi[18] = vp[6];    fvi[19] = vp[7];
                    fvi[20] = vp[3];    fvi[21] = vp[2];    fvi[22] = vp[1];    fvi[23] = vp[0];
                    break;
                    
                default:
                    break;
                
                    
            }
            
        }
        
        // save irregular volume subset to cache file
        if(m_use_cache)
        {
            std::ofstream ofile_cache(cache_file_name_str.str().c_str(), std::ios_base::out | std::ios_base::binary);
            if(ofile_cache)
                save_irregular_volume_to_cache(ofile_cache, mesh_storage, irregular_volume_subset);
        }
    }
    return true;
}

//----------------------------------------------------------------------
bool Irregular_volume_importer::load_irregular_volume_from_cache(
    std::ifstream                                   &ifile,
    nv::index::IIrregular_volume_subset*            irregular_volume_subset) const
{
    // read mesh storage params
    nv::index::IIrregular_volume_subset::Mesh_parameters mesh_params;
    ifile.read(reinterpret_cast<char*>(&mesh_params), 
        sizeof(nv::index::IIrregular_volume_subset::Mesh_parameters));
            
    // create and mesh storage arrays
    nv::index::IIrregular_volume_subset::Mesh_storage mesh_storage;
    if(!irregular_volume_subset->generate_mesh_storage(mesh_params, mesh_storage))
        return false;
        
    ifile.read(reinterpret_cast<char*>(mesh_storage.vertices),
        sizeof(mi::math::Vector_struct<mi::Float32, 3>)*mesh_params.nb_vertices);
        
    ifile.read(reinterpret_cast<char*>(mesh_storage.face_vtx_indices),
        sizeof(mi::Uint32)*mesh_params.nb_face_vtx_indices);
        
    ifile.read(reinterpret_cast<char*>(mesh_storage.faces),
        sizeof(nv::index::IIrregular_volume_subset::Face)*mesh_params.nb_faces);
        
    ifile.read(reinterpret_cast<char*>(mesh_storage.cell_face_indices),
        sizeof(mi::Uint32)*mesh_params.nb_cell_face_indices);
        
    ifile.read(reinterpret_cast<char*>(mesh_storage.cells),
        sizeof(nv::index::IIrregular_volume_subset::Cell)*mesh_params.nb_cells);
    
    // read attribute parameters
    mi::Uint32 nb_attributes = 0u;
    ifile.read(reinterpret_cast<char*>(&nb_attributes), sizeof(mi::Uint32));

    for(mi::Uint32 i=0; i<nb_attributes; ++i)
    {
        mi::Uint32 attribute_index = 0u;
        ifile.read(reinterpret_cast<char*>(&attribute_index), sizeof(mi::Uint32));
        
        nv::index::IIrregular_volume_subset::Attribute_parameters attrib_params;
        ifile.read(reinterpret_cast<char*>(&attrib_params), 
            sizeof(nv::index::IIrregular_volume_subset::Attribute_parameters));

        // create and read attribute storage
        nv::index::IIrregular_volume_subset::Attribute_storage attribute_storage;
        if(!irregular_volume_subset->generate_attribute_storage(attribute_index, attrib_params, attribute_storage))
            return false;
            
        switch(attrib_params.type)
        {
            case nv::index::IIrregular_volume_subset::ATTRIB_TYPE_UINT8:
            
                ifile.read(reinterpret_cast<char*>(attribute_storage.attrib_values),
                    sizeof(mi::Uint8)*attrib_params.nb_attrib_values);
                break;
                
            case nv::index::IIrregular_volume_subset::ATTRIB_TYPE_UINT16:
            
                ifile.read(reinterpret_cast<char*>(attribute_storage.attrib_values),
                    sizeof(mi::Uint16)*attrib_params.nb_attrib_values);
                break;

            case nv::index::IIrregular_volume_subset::ATTRIB_TYPE_FLOAT32:
            
                ifile.read(reinterpret_cast<char*>(attribute_storage.attrib_values),
                    sizeof(mi::Float32)*attrib_params.nb_attrib_values);
                break;
        }
    }
    
    return true;
}

//----------------------------------------------------------------------
bool Irregular_volume_importer::save_irregular_volume_to_cache(
    std::ofstream                                           &ofile,
    const nv::index::IIrregular_volume_subset::Mesh_storage &mesh_storage,
    const nv::index::IIrregular_volume_subset*              irregular_volume_subset) const
{
    nv::index::IIrregular_volume_subset::Mesh_parameters mesh_params = 
        irregular_volume_subset->get_mesh_parameters();

    // write mesh storage parameters
    ofile.write(reinterpret_cast<char*>(&mesh_params), 
        sizeof(nv::index::IIrregular_volume_subset::Mesh_parameters));
        
    // write mesh storage
    ofile.write(reinterpret_cast<char*>(mesh_storage.vertices),
        sizeof(mi::math::Vector_struct<mi::Float32, 3>)*mesh_params.nb_vertices);
        
    ofile.write(reinterpret_cast<char*>(mesh_storage.face_vtx_indices),
        sizeof(mi::Uint32)*mesh_params.nb_face_vtx_indices);
        
    ofile.write(reinterpret_cast<char*>(mesh_storage.faces),
        sizeof(nv::index::IIrregular_volume_subset::Face)*mesh_params.nb_faces);
        
    ofile.write(reinterpret_cast<char*>(mesh_storage.cell_face_indices),
        sizeof(mi::Uint32)*mesh_params.nb_cell_face_indices);
        
    ofile.write(reinterpret_cast<char*>(mesh_storage.cells),
        sizeof(nv::index::IIrregular_volume_subset::Cell)*mesh_params.nb_cells);
        
    // write attribute storage
    mi::Uint32 nb_attributes = irregular_volume_subset->get_nb_attributes();
    for(mi::Uint32 i=0; i<nb_attributes; ++i)
    {
        // attribute index
        ofile.write(reinterpret_cast<char*>(&i), sizeof(mi::Uint32));
        
        //attribute params
        nv::index::IIrregular_volume_subset::Attribute_parameters attrib_params;
        irregular_volume_subset->get_attribute_parameters(i, attrib_params);
        ofile.write(reinterpret_cast<char*>(&attrib_params), 
            sizeof(nv::index::IIrregular_volume_subset::Attribute_parameters));
            
        nv::index::IIrregular_volume_subset::Attribute_storage attribute_storage;
        irregular_volume_subset->get_attribute(i, attribute_storage);
        
        switch(attrib_params.type)
        {
            case nv::index::IIrregular_volume_subset::ATTRIB_TYPE_UINT8:
            
                ofile.write(reinterpret_cast<char*>(attribute_storage.attrib_values),
                    sizeof(mi::Uint8)*attrib_params.nb_attrib_values);
                break;
                
            case nv::index::IIrregular_volume_subset::ATTRIB_TYPE_UINT16:
            
                ofile.write(reinterpret_cast<char*>(attribute_storage.attrib_values),
                    sizeof(mi::Uint16)*attrib_params.nb_attrib_values);
                break;

            case nv::index::IIrregular_volume_subset::ATTRIB_TYPE_FLOAT32:
            
                ofile.write(reinterpret_cast<char*>(attribute_storage.attrib_values),
                    sizeof(mi::Float32)*attrib_params.nb_attrib_values);
                break;
        }
    }
    return true;
}


namespace esg
{

//----------------------------------------------------------------------
/// read a line
inline mi::Uint32 read_line(
    std::ifstream&  ifile, 
    char            result[256])
{
    ifile.getline(result, 256);
    if (ifile.fail())
    {
        ifile.clear();
        return 0;
    }

    return 1;
}

//----------------------------------------------------------------------
/// read a next line
inline mi::Uint32 read_next_data_line(
    std::ifstream&  ifile, 
    char            result[256])
{
    bool is_comment = true;
    mi::Uint32 value = 1;

    while( is_comment && value )
    {

        value = read_line(ifile, result);
        if( *result && result[0] != '#' )
        {
            size_t len = strlen( result );
            unsigned int i = 0;
            while( i < len && (static_cast<unsigned int>(result[i]) <= 255) && isspace(result[i]) )
            {
                ++i;
            }
            // If there was only space characters this is a comment, thus skip it
            if( i != len )
            {
                // The line was not empty, not beginning by '#' and not composed
                // of only white space, this is not a comment
                is_comment = false;
            }
        }
    }

    return value;
}

//----------------------------------------------------------------------
Geometry::Geometry()
{
    type = TYPE_ELEMENT_INVALID;
    part_id = ~0u;
    nb_vertices = 0u;
    nb_cells = 0u;
}

void Geometry::rebuild_subset_indices()
{
    // build vertex list subset and remap indices
    std::vector<mi::math::Vector_struct<mi::Float32, 3> > subset_vertices;
    std::vector<mi::Float32> subset_scalar;
    std::map<mi::Uint32, mi::Uint32> vidx_to_vjdx;

    std::vector<mi::Uint32>::iterator it;
    
    for(it=cell_vtx_ids.begin(); it!=cell_vtx_ids.end(); ++it)
    {
        std::map<mi::Uint32, mi::Uint32>::const_iterator kt = vidx_to_vjdx.find(*it);
        if(kt != vidx_to_vjdx.end())
        {
            *it = kt->second;
        }
        else
        {
            mi::Uint32 vjdx = subset_vertices.size();
            vidx_to_vjdx[*it] = vjdx;
            subset_vertices.push_back(vertices[*it]);
            subset_scalar.push_back(scalar[*it]);
            *it = vjdx;
        }
    }
    
    vertices = subset_vertices;
    scalar = subset_scalar;
}

//----------------------------------------------------------------------
void read_case_from_file(
    const mi::math::Bbox_struct<mi::Float32, 3>& bbox,
    const std::string&                           geometry_file,
    const std::string&	                         scalar_file,
    Geometry_list&                               geometry_list)
{
    std::string file_extension = geometry_file.find_last_of(".") != std::string::npos ? 
                    geometry_file.substr(geometry_file.find_last_of(".") + 1) : "";
                    
    if(file_extension == "" || file_extension != "geom")
    {
        ERROR_LOG << "Error opening geometry file: '" << geometry_file << "'.";
        ERROR_LOG << "Unknown file extension.";
        return;
    }

    // Opening the text file as binary. If not, the reader fails to read
    // files with Unix line endings on Windows machines.
    std::ifstream ifile(geometry_file.c_str(), std::ios_base::in|std::ios_base::binary);
    if (!ifile) 
    {
        ERROR_LOG << "Error opening geometry file: '" << geometry_file << "'.";
        return;
    }
    
    char line[256];
    char sub_line[256];

    // find out weather is binary or text ensight gold file format.
    // Description line 1 or 'C binary'
    read_next_data_line(ifile, line);
    sscanf(line, " %*s %s", sub_line);
    ifile.close();

    // is a binary file? Then use the binary parser
    if (strncmp(sub_line, "Binary", 6) == 0)
    {
        read_geometry_from_binary_file(bbox, geometry_file, geometry_list);
        read_scalar_from_binary_file(scalar_file, geometry_list);
    }
    else
    {
        read_geometry_from_ascii_file(bbox, geometry_file, geometry_list);
        read_scalar_from_ascii_file(scalar_file, geometry_list);
    }
}

//----------------------------------------------------------------------
void convert_case_to_ts(
    const mi::math::Bbox_struct<mi::Float32, 3>& bbox,
    const std::string&	                         geometry_file,
    const std::string&	                         scalar_file,
    const std::string&	                         displ_file,
    const std::string&	                         ts_file,
    mi::Uint32      	                         nb_frames,
    Geometry_list&                               geometry_list)
{
    ERROR_LOG << "Scalar: " << scalar_file;
    ERROR_LOG << "Displacement: " << displ_file;
    
    std::string file_extension = geometry_file.find_last_of(".") != std::string::npos ? 
                    geometry_file.substr(geometry_file.find_last_of(".") + 1) : "";
                    
    if(file_extension == "" || file_extension != "geom")
    {
        ERROR_LOG << "Error opening geometry file: '" << geometry_file << "'.";
        ERROR_LOG << "Unknown file extension.";
        return;
    }

    // Opening the text file as binary. If not, the reader fails to read
    // files with Unix line endings on Windows machines.
    std::ifstream ifile(geometry_file.c_str(), std::ios_base::in|std::ios_base::binary);
    if (!ifile) 
    {
        ERROR_LOG << "Error opening geometry file: '" << geometry_file << "'.";
        return;
    }
    
    char line[256];
    char sub_line[256];

    // find out weather is binary or text ensight gold file format.
    // Description line 1 or 'C binary'
    read_next_data_line(ifile, line);
    sscanf(line, " %*s %s", sub_line);
    ifile.close();

    // is a binary file? Then use the binary parser
    if (strncmp(sub_line, "Binary", 6) == 0)
    {
        read_geometry_from_binary_file(bbox, geometry_file, geometry_list);
        read_scalar_from_binary_file(scalar_file, geometry_list);
    }
    else
    {
        const mi::Float32 mult = 1.0f;
        // std::string scalar_base_name ("/h/bigdata/medical/bsc_heart/ensimerge1M/heart.ensi.INTRA-");
        // std::string displ_base_name ("/h/bigdata/medical/bsc_heart/ensimerge1M/heart.ensi.DISPL-");
        char buff[256];
        
        mi::math::Bbox_struct<mi::Float32, 3> bbox2;
        bbox2.min.x = bbox2.min.y = bbox2.min.z =  -10000.0f;
        bbox2.max.x = bbox2.max.y = bbox2.max.z =  10000.0f;

        read_geometry_from_ascii_file(bbox2, geometry_file, geometry_list);
        
        
        for(mi::Uint32 frame=1; frame <= nb_frames; frame++)
        {
            
            sprintf(buff, "%06d", frame);            
            // std::string ofname = std::string("heart-intra-") + std::string(buff) + std::string(".ts"); 
            // std::string ofname = std::string("coarse_heart-") + std::string(buff) + std::string(".ts"); 
            std::string ofname = ts_file + "-" + std::string(buff) + std::string(".ts"); 
            std::ofstream ofile(ofname.c_str(), std::ios_base::out | std::ios_base::binary);
            
            if(ofile)
            {
                geometry_list[0].scalar.clear();
                
                std::string scalar_file_enum = scalar_file + "-" + std::string(buff);                
                ERROR_LOG << "Scalar file: " << scalar_file_enum;
                read_scalar_from_ascii_file(scalar_file_enum, geometry_list);
                
                std::string displ_file_enum = displ_file + "-" + std::string(buff);                
                ERROR_LOG << "Displacement file: " << displ_file_enum;
                read_displacement_from_ascii_file(displ_file_enum, geometry_list);
                
                ofile << geometry_list[0].nb_vertices << " " << geometry_list[0].nb_cells << std::endl;
                
                for(mi::Uint32 i=0; i < geometry_list[0].nb_vertices; i++)
                {
                    ofile << 
                        geometry_list[0].vertices[i].x + geometry_list[0].displ[i].x*mult << " " <<
                        geometry_list[0].vertices[i].y + geometry_list[0].displ[i].y*mult << " " <<
                        geometry_list[0].vertices[i].z + geometry_list[0].displ[i].z*mult << " " <<
                        geometry_list[0].scalar[i] << std::endl;
                }
                
                for(mi::Uint32 i=0; i < geometry_list[0].nb_cells; i++)
                {
                    mi::Uint32 idx = 4*i;
                    ofile << 
                        geometry_list[0].cell_vtx_ids[idx] << " " <<
                        geometry_list[0].cell_vtx_ids[idx+1] << " " <<
                        geometry_list[0].cell_vtx_ids[idx+2] << " " <<
                        geometry_list[0].cell_vtx_ids[idx+3] << std::endl;
                }
            }
            
          
            ofile.close();
        }
      
        ERROR_LOG << "Conversion from Ensight to TS done!!!!!!!!!!!!";
        // read_geometry_from_ascii_file(bbox, geometry_file, geometry_list);
        // read_scalar_from_ascii_file(scalar_file, geometry_list);
    }
}


//----------------------------------------------------------------------
void read_geometry_from_ascii_file(
    const mi::math::Bbox_struct<mi::Float32, 3>& bbox,
    const std::string&                           geometry_file,
    Geometry_list&                               geometry_list)
{
    std::string file_extension = geometry_file.find_last_of(".") != std::string::npos ? 
                    geometry_file.substr(geometry_file.find_last_of(".") + 1) : "";
                    
    if(file_extension == "" || file_extension != "geom")
    {
        ERROR_LOG << "Error opening geometry file: '" << geometry_file << "'.";
        ERROR_LOG << "Unknown file extension.";
        return;
    }
    
    const mi::math::Bbox<mi::Float32, 3> subset_bbox(bbox);

    // Opening the text file as binary. If not, the reader fails to read
    // files with Unix line endings on Windows machines.
    std::ifstream ifile(geometry_file.c_str(), std::ios_base::in|std::ios_base::binary);
    if (!ifile) 
    {
        ERROR_LOG << "Error opening geometry file: '" << geometry_file << "'.";
        return;
    }
    
    mi::Uint32 line_count;
    char line[256];
    char sub_line[256];

    // Description line 1.
    read_next_data_line(ifile, line);
    
    // Description line 2. it could be blank
    read_line(ifile, line);
    
    // Node-id
    bool node_id_listed = false;
    read_next_data_line(ifile, line);
    sscanf(line, " %*s %*s %s", sub_line);
    if (strncmp(sub_line, "given", 5) == 0 || strncmp(sub_line, "ignore", 6) == 0)
        node_id_listed = true;

    // Element-id
    bool element_id_listed = false;
    read_next_data_line(ifile, line);
    sscanf(line, " %*s %*s %s", sub_line);
    if (strncmp(sub_line, "given", 5) == 0 || strncmp(sub_line, "ignore", 6) == 0)
        element_id_listed = true;
    
    line_count = read_next_data_line(ifile, line); // "extents" or "part"
    if (strncmp(line, "extents", 7) == 0)
    {
        // Skipping the extent lines for now.
        read_next_data_line(ifile, line);
        read_next_data_line(ifile, line);
        read_next_data_line(ifile, line);

        line_count = read_next_data_line(ifile, line); // "part"
    }
    
    while (line_count && strncmp(line, "part", 4) == 0)
    {
        // Part-id
        read_next_data_line(ifile, line);
        mi::Uint32 part_id = static_cast<mi::Uint32>(atoi(line) - 1); // EnSight starts #ing at 1.

        // part description line or interface?
        read_next_data_line(ifile, line); 
        if (strncmp(line, "interface", 9) == 0)
        {
            ERROR_LOG << "Error reading geometry file: '" << geometry_file << "'.";
            ERROR_LOG << "Interface sections are not supported.";
            return;
        }

        // Coordinate or block?
        read_next_data_line(ifile, line);
        if (strncmp(line, "coordinates", 11) != 0)
        {
            ERROR_LOG << "Error reading geometry file: '" << geometry_file << "'.";
            ERROR_LOG << "Expecting coordinates.";
            return;
        }
        
        geometry_list.push_back(Geometry());
        Geometry &geometry = geometry_list.back();
        geometry.part_id  = part_id;

        // Number of coordinates
        read_next_data_line(ifile, line);
        geometry.nb_vertices = static_cast<mi::Uint32>(atoi(line));

        // coordinates ids are provided when node-id are listed. Not used for now
        // std::vector<mi::Uint32> coord_ids;
        if(node_id_listed)
        {
            // coord_ids.resize(geometry.nb_vertices);
            for(mi::Uint32 i=0; i<geometry.nb_vertices; ++i)
            {
                read_next_data_line(ifile, line);
                // coord_ids[i] = static_cast<mi::Uint32>(atoi(line));
            }
        }

        // Vertex coordinates
        mi::math::Bbox<mi::Float32, 3> dataset_extent;

        geometry.vertices.resize(geometry.nb_vertices);

        for(mi::Uint32 i=0; i<geometry.nb_vertices; ++i)
        {
            read_next_data_line(ifile, line);
            geometry.vertices[i].x = static_cast<mi::Float32>(atof(line));
            dataset_extent.min.x = mi::math::min(dataset_extent.min.x, geometry.vertices[i].x);
            dataset_extent.max.x = mi::math::max(dataset_extent.max.x, geometry.vertices[i].x);
        }

        for(mi::Uint32 i=0; i<geometry.nb_vertices; i++)
        {
            read_next_data_line(ifile, line);
            geometry.vertices[i].y = static_cast<mi::Float32>(atof(line));
            dataset_extent.min.y = mi::math::min(dataset_extent.min.y, geometry.vertices[i].y);
            dataset_extent.max.y = mi::math::max(dataset_extent.max.y, geometry.vertices[i].y);
        }

        for(mi::Uint32 i=0; i<geometry.nb_vertices; ++i)
        {
            read_next_data_line(ifile, line);
            geometry.vertices[i].z = static_cast<mi::Float32>(atof(line));
            dataset_extent.min.z = mi::math::min(dataset_extent.min.z, geometry.vertices[i].z);
            dataset_extent.max.z = mi::math::max(dataset_extent.max.z, geometry.vertices[i].z);
        }

        // INFO_LOG << dataset_extent;

        // Element type
        read_next_data_line(ifile, line);
        
        // Hexahedron 8 nodes
        if (strncmp(line, "hexa8", 5) == 0)
        {
            geometry.type = TYPE_ELEMENT_HEXAHEDRON;
            geometry.nb_cells = 0u;
            
            // Number of elements (cells)
            read_next_data_line(ifile, line);
            mi::Uint32 nb_cells = static_cast<mi::Uint32>(atoi(line));
                
            // Element-ids. Provided when Element-ids are listed. Not used for now
            // std::vector<mi::Uint32> element_ids;
            if(element_id_listed)
            {
                // element_ids.resize(nb_cells);
                for(mi::Uint32 i=0u; i<nb_cells; ++i)
                {
                    read_next_data_line(ifile, line);
                    // element_ids[i] = static_cast<mi::Uint32>(atoi(line));
                }
            }

            // face vertex indices
            mi::math::Bbox<mi::Float32, 3> cur_bbox;
            std::vector<mi::Uint32> ids(8);
            for(mi::Uint32 i=0u; i<nb_cells; ++i)
            {
                read_next_data_line(ifile, line);
                sscanf(line, " %d %d %d %d %d %d %d %d", &ids[0], &ids[1], &ids[2], &ids[3], 
                    &ids[4], &ids[5], &ids[6], &ids[7]);
                
                // Ensight start #ing at 1, move it to zeroing
                for(mi::Uint32 j=0u; j<8u; ++j)
                    ids[j] -= 1u;
                    
                // bounding box test
                cur_bbox.clear();
                for(mi::Uint32 j=0u; j<8u; ++j)
                    cur_bbox.insert(geometry.vertices[ids[j]]);
                
                if (cur_bbox.intersects(subset_bbox))
                {
                    geometry.cell_vtx_ids.insert(geometry.cell_vtx_ids.end(), ids.begin(), ids.end());
                    geometry.nb_cells++;
                }
            }
            
        }
        // Tetrahedron 4 nodes
        else if (strncmp(line, "tetra4", 6) == 0)
        {
            geometry.type = TYPE_ELEMENT_TETRAHEDRON;
            geometry.nb_cells = 0u;
        
            // Number of elements (faces)
            read_next_data_line(ifile, line);
            mi::Uint32 nb_cells = static_cast<mi::Uint32>(atoi(line));
                
            // Element-ids. Provided when Element-ids are listed. Not used for now
            // std::vector<mi::Uint32> element_ids;
            if(element_id_listed)
            {
                // element_ids.resize(nb_cells);
                for(mi::Uint32 i=0u; i<nb_cells; ++i)
                {
                    read_next_data_line(ifile, line);
                    // element_ids[i] = static_cast<mi::Uint32>(atoi(line));
                }
            }

            // face vertex indices
            mi::math::Bbox<mi::Float32, 3> cur_bbox;
            std::vector<mi::Uint32> ids(4);
            for(mi::Uint32 i=0u; i<nb_cells; ++i)
            {
                read_next_data_line(ifile, line);
                sscanf(line, " %d %d %d %d", &ids[0], &ids[1], &ids[2], &ids[3]);                
                
                // Ensight start #ing at 1, move it to zeroing
                for(mi::Uint32 j=0u; j<4u; ++j)
                    ids[j] -= 1u;
                    
                // bounding box test
                cur_bbox.clear();
                for(mi::Uint32 j=0u; j<4u; ++j)
                    cur_bbox.insert(geometry.vertices[ids[j]]);
                
                if (cur_bbox.intersects(subset_bbox))
                {
                    geometry.cell_vtx_ids.insert(geometry.cell_vtx_ids.end(), ids.begin(), ids.end());
                    geometry.nb_cells++;
                }
                
            }
        }
        // Quadrangle 4 nodes. Not used yet
        else if (strncmp(line, "quad4", 5) == 0)
        {
            // free vertices, we are not supporting quad objects.
            // but keep number of vertices that is needed to read the scalar values
            geometry.vertices.clear();
            geometry.nb_cells = 0u;
           
            // Number of elements (faces)
            read_next_data_line(ifile, line);
            mi::Uint32 nb_cells = static_cast<mi::Uint32>(atoi(line));
                
            // Element-ids. Provided when Element-ids are listed. Not used for now
            // std::vector<mi::Uint32> element_ids;
            if(element_id_listed)
            {
                // element_ids.resize(nb_cells);
                for(mi::Uint32 i=0; i<nb_cells; ++i)
                {
                    read_next_data_line(ifile, line);
                    // element_ids[i] = static_cast<mi::Uint32>(atoi(line));
                }
            }

            // face vertex indices
            // std::vector<mi::Uint32> cell_vtx_ids;
            // cell_vtx_ids.resize(nb_cells*4);

            mi::Uint32 dummy_ids[4];
            for(mi::Uint32 i=0; i<nb_cells; ++i)
            {
                read_next_data_line(ifile, line);
                sscanf(line, " %d %d %d %d", dummy_ids, dummy_ids+1, dummy_ids+2, dummy_ids+3);                
            }
        }
        else
        {
            ERROR_LOG << "Error reading geometry file: '" << geometry_file << "'.";
            ERROR_LOG << "Element type not supported.";
            return;
        }

        line_count = read_next_data_line(ifile, line); // "part or eof"
    }
    
    ifile.close();
}

//----------------------------------------------------------------------
void read_geometry_from_binary_file(
    const mi::math::Bbox_struct<mi::Float32, 3>& bbox,
    const std::string&                           geometry_file,
    Geometry_list&                               geometry_list)
{
}

void read_scalar_from_ascii_file(
    const std::string&	scalar_file,
    Geometry_list &geometry_list)
{
    // Opening the text file as binary. If not, the reader fails to read
    // files with Unix line endings on Windows machines.
    std::ifstream ifile(scalar_file.c_str(), std::ios_base::in|std::ios_base::binary);
    if (!ifile) 
    {
        ERROR_LOG << "Error opening scalar file: '" << scalar_file << "'.";
        return;
    }
    
    char line[256];
    
    // Description line
    read_next_data_line(ifile, line);

    // for each part
    while(read_next_data_line(ifile, line) && strncmp(line, "part", 4) == 0)
    {
        // Part-id
        read_next_data_line(ifile, line);
        mi::Uint32 part_id = static_cast<mi::Uint32>(atoi(line) - 1); // EnSight starts #ing at 1.
        
        //find associated geometry
        mi::Uint32 geom_idx = ~0u;
        mi::Uint32 nb_geometries = geometry_list.size();
        
        for(mi::Uint32 i=0u; i<nb_geometries; ++i)
        {
            if(geometry_list[i].part_id == part_id)
            {
                geom_idx = i;
                break;
            }
        }
        
        if(geom_idx == ~0u)
        {
            ERROR_LOG << "Error reading scalar file: '" << scalar_file << "'.";
            ERROR_LOG << "No part-id match between geometry and scalar files.";
            return;
        }
        
        Geometry &geometry = geometry_list[geom_idx];
        
        // // Element type. IGNORE
        // read_next_data_line(ifile, line);
        
        // "coordinates" or "block"
        read_next_data_line(ifile, line);
        
        if(strncmp(line, "coordinates", 11) != 0)
        {
            ERROR_LOG << "Error reading scalar file: '" << scalar_file << "'.";
            ERROR_LOG << "Expecting 'coordinates' keyword.";
            return;
        }
        
        // read scalars
        geometry.scalar.resize(geometry.nb_vertices);
        for(mi::Uint32 i=0; i<geometry.nb_vertices; ++i)
        {
            read_next_data_line(ifile, line);            
            geometry.scalar[i] = static_cast<mi::Float32>(atof(line));
        }
    }
    
    ifile.close();

}

void read_scalar_from_binary_file(
    const std::string&	scalar_file,
    Geometry_list &geometry_list)
{
}

/// read displacement data: ascii
void read_displacement_from_ascii_file(
    const std::string&	displ_file,
    Geometry_list &geometry_list)
{
    // Opening the text file as binary. If not, the reader fails to read
    // files with Unix line endings on Windows machines.
    std::ifstream ifile(displ_file.c_str(), std::ios_base::in|std::ios_base::binary);
    if (!ifile) 
    {
        ERROR_LOG << "Error opening displacement file: '" << displ_file << "'.";
        return;
    }
    
    char line[256];
    
    // Description line
    read_next_data_line(ifile, line);

    // for each part
    while(read_next_data_line(ifile, line) && strncmp(line, "part", 4) == 0)
    {
        // Part-id
        read_next_data_line(ifile, line);
        mi::Uint32 part_id = static_cast<mi::Uint32>(atoi(line) - 1); // EnSight starts #ing at 1.
        
        //find associated geometry
        mi::Uint32 geom_idx = ~0u;
        mi::Uint32 nb_geometries = geometry_list.size();
        
        for(mi::Uint32 i=0u; i<nb_geometries; ++i)
        {
            if(geometry_list[i].part_id == part_id)
            {
                geom_idx = i;
                break;
            }
        }
        
        if(geom_idx == ~0u)
        {
            ERROR_LOG << "Error reading scalar file: '" << displ_file << "'.";
            ERROR_LOG << "No part-id match between geometry and displacement files.";
            return;
        }
        
        Geometry &geometry = geometry_list[geom_idx];
        
        // // Element type. IGNORE
        // read_next_data_line(ifile, line);
        
        // "coordinates" or "block"
        read_next_data_line(ifile, line);
        
        if(strncmp(line, "coordinates", 11) != 0)
        {
            ERROR_LOG << "Error reading displacement file: '" << displ_file << "'.";
            ERROR_LOG << "Expecting 'coordinates' keyword.";
            return;
        }
        
        // read scalars
        geometry.displ.resize(geometry.nb_vertices);
        for(mi::Uint32 i=0; i<geometry.nb_vertices; ++i)
        {
            read_next_data_line(ifile, line);            
            geometry.displ[i].x = static_cast<mi::Float32>(atof(line));
        }
        
        for(mi::Uint32 i=0; i<geometry.nb_vertices; ++i)
        {
            read_next_data_line(ifile, line);            
            geometry.displ[i].y = static_cast<mi::Float32>(atof(line));
        }
        
        for(mi::Uint32 i=0; i<geometry.nb_vertices; ++i)
        {
            read_next_data_line(ifile, line);            
            geometry.displ[i].z = static_cast<mi::Float32>(atof(line));
        }
        
    }
    
    ifile.close();
    
}



} // namespace nv::index_common::esg
//----------------------------------------------------------------------
}} // namespace nv::index_common
 
