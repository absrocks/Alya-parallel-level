/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief external sparse volume data importer job

#include "sparse_volume_importer.h"

#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <vector>

#include <nv/index/isparse_volume_subset.h>

#include "forwarding_logger.h"

// disable OpenVDB warnings
#if WIN_NT
#pragma warning(push)
#pragma warning(disable: 4800 4244 4146)
#endif

#include <openvdb/openvdb.h>

typedef openvdb::tree::Tree<
            openvdb::tree::RootNode<
                openvdb::tree::InternalNode<
                    openvdb::tree::InternalNode<
                        openvdb::tree::InternalNode<
                            openvdb::tree::LeafNode<openvdb::Vec3f,4>,
                        3>,
                    3>,
                3>
            >
        > VDB_3334_tree_vec3f; 
typedef openvdb::Grid<VDB_3334_tree_vec3f>   VDB_3334_grid_vec3f;

typedef openvdb::tree::Tree<
            openvdb::tree::RootNode<
                openvdb::tree::InternalNode<
                    openvdb::tree::InternalNode<
                        openvdb::tree::InternalNode<
                            openvdb::tree::LeafNode<float,4>,
                        3>,
                    3>,
                3>
            >
        > VDB_3334_tree_fp32; 

typedef openvdb::Grid<VDB_3334_tree_fp32>   VDB_3334_grid_fp32;

// Todo chris:
// * inside_bbox: is it an intersection when b.min == a.max? they then share no voxels...
// * data->setBackground? not defined anymore in OpenVDB 3.1.0?
// * vector<>.data() not available as we are external code and not using C++11ff yet

namespace nv {
namespace index_common {

// some typedefs
typedef mi::math::Vector<mi::Sint32, 3>         Vec3i;
typedef mi::math::Vector<mi::Uint32, 3>         Vec3u;

typedef mi::math::Vector_struct<mi::Sint32, 3>  Vec3i_struct;
typedef mi::math::Vector_struct<mi::Uint32, 3>  Vec3u_struct;

typedef mi::math::Bbox<mi::Sint32, 3>           Bbox3i;
typedef mi::math::Bbox<mi::Uint32, 3>           Bbox3u;

typedef mi::math::Bbox_struct <mi::Sint32, 3>   Bbox3i_struct;
typedef mi::math::Bbox_struct <mi::Uint32, 3>   Bbox3u_struct;

namespace util {

struct VDBSetup {
    VDBSetup() {
        openvdb::initialize();
        VDB_3334_grid_vec3f::registerGrid();
        VDB_3334_grid_fp32::registerGrid();
    }
    ~VDBSetup() {
        openvdb::uninitialize();
    }
};

/// \returns    true if the given brick falls within the bounding box.
static bool inside_bbox(
    const openvdb::Coord&   corner,
    const Vec3i&            extend,
    const Bbox3i_struct&    bbox)
{
    // The brick is a range, of course.
    const Vec3i  in_min(corner.x(), corner.y(), corner.z());
    const Vec3i  in_max = in_min + extend;

    // The brick might intersect the bbox but not be contained wholly within
    // it.  We're conservative for now and say an intersection is good enough.
    return    (bbox.min.x < in_max.x && in_min.x < bbox.max.x)
           && (bbox.min.y < in_max.y && in_min.y < bbox.max.y)
           && (bbox.min.z < in_max.z && in_min.z < bbox.max.z);
}

static size_t nbricks(const std::vector<float>& bkdata, const size_t bsize[3])
{
    const size_t nvoxels = bsize[0]*bsize[1]*bsize[2];
    assert(bkdata.size() % nvoxels == 0); // incorrect brick size?
    return bkdata.size() / nvoxels;
}

static bool is_vdb(const std::string& filename)
{
    std::ifstream ifs(filename.c_str());
    
    if(!ifs)
    {
        ERROR_LOG << "Sparse volume importer: could not access '" << filename << "'.";
        return false;
    }

    mi::Uint32 file_magic = 0u;
    ifs.read(reinterpret_cast<char*>(&file_magic), sizeof(file_magic));

    if (ifs.gcount() != sizeof(file_magic))
    {
        ERROR_LOG << "Sparse volume importer: could not read magic from '" << filename << "'.";
        ifs.close();
        return false;
    }
    ifs.close();

    return file_magic == openvdb::OPENVDB_MAGIC;
}

template<typename GRID_TYPE>
struct VDB_grid_types
{
    // intentionally empty
};

template<>
struct VDB_grid_types<openvdb::FloatGrid>
{
    typedef openvdb::FloatGrid              grid_type;
    typedef grid_type::TreeType             tree_type;
    typedef tree_type::LeafNodeType         leaf_node_type;

    static const mi::Sint32                 grid_leaf_extent = 8;
};

template<>
struct VDB_grid_types<VDB_3334_grid_fp32>
{
    typedef VDB_3334_grid_fp32              grid_type;
    typedef VDB_3334_grid_fp32::TreeType    tree_type;
    typedef tree_type::LeafNodeType         leaf_node_type;

    static const mi::Sint32                 grid_leaf_extent = 16;
};

template<typename GRID_TYPE>
void load_vdb_grid_leaves(
    openvdb::GridBase::Ptr&     vdata,
    std::vector<float>&         bkdata,
    std::vector<Vec3i_struct>&  bkpos,
    const Bbox3i_struct&        bounding_box)
{
    using namespace openvdb;

    const Vec3i leaf_extent(VDB_grid_types<GRID_TYPE>::grid_leaf_extent,
                            VDB_grid_types<GRID_TYPE>::grid_leaf_extent,
                            VDB_grid_types<GRID_TYPE>::grid_leaf_extent);

    typedef typename VDB_grid_types<GRID_TYPE>::grid_type vdb_grid_type;
    typename vdb_grid_type::Ptr data = GridBase::grid<vdb_grid_type>(vdata);
    //data->setBackground(0.f);

    Bbox3i data_bbox;

    typedef typename VDB_grid_types<GRID_TYPE>::tree_type::LeafCIter LeafIter;
    for(LeafIter iter = data->tree().cbeginLeaf(); iter; ++iter)
    {
        const typename VDB_grid_types<GRID_TYPE>::leaf_node_type& leaf = *iter;
        // what brick is this?
        Coord origin;
        leaf.getOrigin(origin);

        { // test the data bbox
            const Vec3i n_min(origin.x(), origin.y(), origin.z());
            const Vec3i n_max = n_min + leaf_extent;
            data_bbox.insert(n_min);
            data_bbox.insert(n_max);
        }

        if (!util::inside_bbox(origin, leaf_extent, bounding_box))
        {
            continue;
        }

        const size_t dims[3] = { static_cast<size_t>(leaf_extent.x),
                                 static_cast<size_t>(leaf_extent.y),
                                 static_cast<size_t>(leaf_extent.z) };
        const size_t bsize   = dims[0]*dims[1]*dims[2];
        const size_t nb      = nbricks(bkdata, dims);
        // Reserve some space for the brick's (meta)data.
        bkdata.resize(bsize*(nb+1));
        bkpos.resize(nb+1);

        // Copy the data into our buffer.
        const size_t offset = bsize*nb;
        const typename VDB_grid_types<GRID_TYPE>::leaf_node_type::Buffer& b = leaf.buffer();

        const mi::Size brick_size = bsize * sizeof(mi::Float32);
        mi::Float32* dst_data = bkdata.data() + offset;
        memcpy(dst_data, b.data(), brick_size);

        // Copy over metadata of where the brick lies.
        Vec3i_struct& pos = bkpos[nb];
        pos.x = origin[0];
        pos.y = origin[1];
        pos.z = origin[2];
    }

    INFO_LOG << "bbox: " << data_bbox << ", bkgrnd: " << data->tree().getBackgroundValue()->str();
}

template<typename GRID_TYPE>
mi::Sint32 generate_attrib_brick_pool(
    openvdb::GridBase::Ptr&             vdata,
    std::vector<float>&                 bkdata,
    std::vector<Vec3i_struct>&          bkpos,
    const Bbox3i_struct&                bounding_box,
    nv::index::ISparse_volume_subset*   svol_subset)
{
    // No ghost data in current files (... I think, needs testing).
    const mi::Uint32 nghost = 0u;
    // All data we have are float.
    const nv::index::Sparse_volume_voxel_format fmt = nv::index::SPARSE_VOLUME_VOXEL_FORMAT_FLOAT32;

    const Vec3u bdims = Vec3u(VDB_grid_types<GRID_TYPE>::grid_leaf_extent,
                              VDB_grid_types<GRID_TYPE>::grid_leaf_extent,
                              VDB_grid_types<GRID_TYPE>::grid_leaf_extent);

    return svol_subset->generate_attribute_brick_pool(bdims, nghost, fmt, bkpos.size());
}

// insert the brick into the pool:
// foreach(brick in m_bkdata):
//   const bsize1 = bs[0]*bs[1]*bs[2];
//   pool->insert_brick(m_bkpos[brick], m_bkdata.data()+brick*bsize1)

static void insert_bricks(
    const std::vector<float>&                       bkdata,
    const size_t                                    bsize[3],
    const std::vector<Vec3i_struct>&                bkpos,
    nv::index::ISparse_volume_subset_brick_pool*    pool)
{
    const size_t nb = bkpos.size();
    assert(nb == nbricks(bkdata, bsize));

    for (size_t i = 0; i < nb; ++i)
    {
        const size_t offset = i * bsize[0] * bsize[1] * bsize[2];
        if (pool->insert_brick(bkpos[i], bkdata.data() + offset) == -1)
        {
            ERROR_LOG << "Sparse volume importer: error inserting brick " << i << " into cache!";
            return;
        }
    }
}

} // namespace util

mi::base::Lock                                      Sparse_volume_importer::m_cached_input_lock;
Sparse_volume_importer::OpenVDB_cache_file_info     Sparse_volume_importer::m_cached_input_file;

Sparse_volume_importer::OpenVDB_cache_file_info::OpenVDB_cache_file_info()
  : m_file(NULL)
{
}

Sparse_volume_importer::OpenVDB_cache_file_info::~OpenVDB_cache_file_info()
{
    clear();
}

void Sparse_volume_importer::OpenVDB_cache_file_info::clear()
{
    m_grid_attrib_density.reset();
    if (m_file != NULL) {
        m_file->close();
        delete m_file;
        m_file = NULL;
    }
    m_file_name.clear();
}

Sparse_volume_importer::Sparse_volume_importer()
  : m_use_cache(false)
{
}

Sparse_volume_importer::Sparse_volume_importer(
        const std::string& input_dir,
        const std::string& cache_out_dir,
        const std::string& file_base_name,
        const std::string& file_ext,
        const std::string& field_name,
        const mi::Uint32   ts_enum_len,
        const mi::Uint32   ts_enum_stride,
        const mi::Uint32   ts_enum_off,
        bool               use_cache)
  : m_input_directory(input_dir)
  , m_cache_output_directory(cache_out_dir)
  , m_file_base_name(file_base_name)
  , m_file_extension(file_ext)
  , m_data_field_name(field_name)
  , m_timestep_enumeration_length(ts_enum_len)
  , m_timestep_enumeration_stride(ts_enum_stride)
  , m_timestep_enumeration_offset(ts_enum_off)
  , m_use_cache(use_cache)
{
    std::stringstream cfg;

    cfg << "importer=vdb\n"
        << "input_directory=" << m_input_directory << "\n"
        << "input_file_base_name=" << m_file_base_name << "\n"
        << "input_file_extension=" << m_file_extension << "\n"
        << "input_file_field=" << m_data_field_name << "\n"
        << "output_cache_directory=" << m_cache_output_directory << "\n"
        << "input_time_step_enum_length=" << m_timestep_enumeration_length << "\n"
        << "input_time_step_enum_stride=" << m_timestep_enumeration_stride << "\n"
        << "input_time_step_enum_offset=" << m_timestep_enumeration_offset << "\n"
        << "cache=" << (m_use_cache ? "true" : "false");

    m_configuration = cfg.str();
}

Sparse_volume_importer::~Sparse_volume_importer()
{
}

mi::Size Sparse_volume_importer::estimate(
    const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{
    return 0ull;
}

nv::index::IDistributed_data_subset* Sparse_volume_importer::create(
    const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
    nv::index::IData_subset_factory*                factory,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{
    return create(bounding_box, 0u, factory, dice_transaction);
}

nv::index::IDistributed_data_subset* Sparse_volume_importer::create(
    const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
    mi::Uint32                                      time_step,
    nv::index::IData_subset_factory*                factory,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{

    std::stringstream file_name;
    file_name << m_file_base_name;

    if (m_timestep_enumeration_length > 0) {
        const mi::Uint32 ts = m_timestep_enumeration_offset + time_step * m_timestep_enumeration_stride;
        file_name << std::setfill('0') << std::setw(m_timestep_enumeration_length) << ts;
    }
    
    file_name << m_file_extension;

    std::stringstream in_file_name;
    in_file_name << m_input_directory << '/' << file_name.str();

    std::ostringstream out_file_name;
    if (m_cache_output_directory.empty())
    {
        out_file_name << m_input_directory;
    }
    else
    {
        out_file_name << m_cache_output_directory;
    }
    
    out_file_name << '/' << file_name.str()
                  << "_field_"  << m_data_field_name
                  << "_bounds_" << bounding_box.min.x << "." << bounding_box.max.x << "x"
                                << bounding_box.min.y << "." << bounding_box.max.y << "x"
                                << bounding_box.min.z << "." << bounding_box.max.z
                  << ".idx_svol_part";


    //TODO generate time-step file names
    return create(in_file_name.str(), out_file_name.str(), m_data_field_name, bounding_box, factory, dice_transaction);
}

nv::index::IDistributed_data_subset* Sparse_volume_importer::create(
    const std::string&                              filename,
    const std::string&                              filename_cache,
    const std::string&                              field_name,
    const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
    nv::index::IData_subset_factory*                factory,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{
    INFO_LOG << "Sparse volume importer: loads '" << filename << "', bounds: " << bounding_box;

    mi::base::Handle<nv::index::ISparse_volume_subset>
        sparse_volume_subset(factory->create<nv::index::ISparse_volume_subset>());

    if (!sparse_volume_subset.is_valid_interface())
    {
        ERROR_LOG << "Sparse volume importer: unable to create a sparse volume subset.";
        return NULL;
    }

    bool cache_file_read = false;
    
    if (m_use_cache) {
        bool cache_file_exists = false;
        {
            std::ifstream in_file(filename_cache.c_str(), std::ios_base::in | std::ios_base::binary);
            cache_file_exists = in_file.is_open();
        }

        if (cache_file_exists) {
            cache_file_read = sparse_volume_subset->load_internal_data_representation(filename_cache.c_str());
        }

        if (cache_file_read) {
            INFO_LOG << "Successfully read cache file: " << filename_cache;
        }
    }

    if (!cache_file_read) {
        const std::string file_extension =
              filename.find_last_of(".") != std::string::npos
            ? filename.substr(filename.find_last_of(".") + 1)
            : "";

        const bool is_vdb_file = (file_extension == "vdb") && util::is_vdb(filename);

        // check for supported file type
        if (!is_vdb_file)
        {
            ERROR_LOG << "Sparse volume importer: unknown or unsupported file type: '" << filename << "'.";
            return NULL;
        }

        if (is_vdb_file)
        {
            if (!load_from_vdb(filename, bounding_box, field_name, sparse_volume_subset.get()))
            {
                return NULL;
            }
        }

        if (m_use_cache) {
            sparse_volume_subset->finalize();
            if (!sparse_volume_subset->store_internal_data_representation(filename_cache.c_str())) {
                ERROR_LOG << "Sparse volume importer: failed to write subset cache file '" << filename_cache << "'.";
            }
        }
    }

    sparse_volume_subset->retain();
    return sparse_volume_subset.get();
}

const char* Sparse_volume_importer::get_configuration() const
{
    return m_configuration.c_str();
}

void serialize_string(
    const std::string&          s,
    mi::neuraylib::ISerializer* serializer)
{
    mi::Uint32 nb_elements = mi::Uint32(s.size());
    serializer->write(&nb_elements, 1);
    serializer->write(reinterpret_cast<const mi::Uint8*>(s.c_str()), nb_elements);
}

void deserialize_string(
    std::string&                  s,
    mi::neuraylib::IDeserializer* deserializer)
{
    mi::Uint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    s.resize(nb_elements);
    deserializer->read(reinterpret_cast<mi::Uint8*>(&s[0]), nb_elements);
}

void Sparse_volume_importer::serialize(
    mi::neuraylib::ISerializer * serializer) const
{
    serializer->write(&m_use_cache, 1);
    
    serialize_string(m_input_directory,         serializer);
    serialize_string(m_cache_output_directory,  serializer);
    serialize_string(m_file_base_name,          serializer);
    serialize_string(m_file_extension,          serializer);

    serialize_string(m_data_field_name,         serializer);

    serializer->write(&m_timestep_enumeration_length, 1);
    serializer->write(&m_timestep_enumeration_stride, 1);
    serializer->write(&m_timestep_enumeration_offset, 1);
}

void Sparse_volume_importer::deserialize(
    mi::neuraylib::IDeserializer * deserializer)
{
    deserializer->read(&m_use_cache, 1);
    
    deserialize_string(m_input_directory,        deserializer);
    deserialize_string(m_cache_output_directory, deserializer);
    deserialize_string(m_file_base_name,         deserializer);
    deserialize_string(m_file_extension,         deserializer);

    deserialize_string(m_data_field_name,        deserializer);

    deserializer->read(&m_timestep_enumeration_length, 1);
    deserializer->read(&m_timestep_enumeration_stride, 1);
    deserializer->read(&m_timestep_enumeration_offset, 1);
}

bool Sparse_volume_importer::load_from_vdb(
    const std::string&                          filename,
    const mi::math::Bbox_struct<mi::Sint32, 3>& bounding_box,
    const std::string&                          data_field_name,
    nv::index::ISparse_volume_subset*           svol_subset) const
{
    using namespace openvdb;
    const util::VDBSetup    vdb_setup_guard;

    {
        mi::base::Lock::Block file_cache_lock(&m_cached_input_lock);

        if (m_cached_input_file.m_file_name != filename) {
            m_cached_input_file.clear();

            m_cached_input_file.m_file = new io::File(filename.c_str());
            if (!m_cached_input_file.m_file->open())
            {
                ERROR_LOG << "Sparse volume importer: could not open VDB '" << filename << "'.";
                m_cached_input_file.clear();
                return false;
            }

            { // read and output meta data for all grids contained in the file
                // JESUS... I want auto here
                const GridPtrVecPtr grids_meta = m_cached_input_file.m_file->readAllGridMetadata();
                
                INFO_LOG << "File '" << filename << "' contains " << grids_meta->size() << "grids: \n";
                for (mi::Size i = 0; i < grids_meta->size(); ++i) {
                    const GridBase::Ptr grid_meta = grids_meta->at(i);
                    std::stringstream os;
                    grid_meta->print(os);
                    INFO_LOG << " *** grid " << i << ": \n" << os.str() << "\n";
                }
            }

            // This class should really have an m_field instead of hard-coding density.
            //const BBoxd rb = BBoxd(Vec3d(0, 0, 0), Vec3d(500, 1630, 1630));
            m_cached_input_file.m_grid_attrib_density = m_cached_input_file.m_file->readGrid(data_field_name/*, rb*/);
            if (false) {
                const GridBase::Ptr grid_meta = m_cached_input_file.m_grid_attrib_density;
                std::stringstream os;
                grid_meta->print(os);
                INFO_LOG << " * grid name:     " << grid_meta->getName() << "\n"
                         << "   - info:        " << os.str() << "\n"
                         << "   - mem_usage:   " << std::fixed << std::setprecision(3)
                                                 << double(grid_meta->memUsage()) / (1024.0 * 1024.0) << "MiB\n"
                         << "   - voxels size: " << grid_meta->voxelSize().x() << ", "
                                                 << grid_meta->voxelSize().y() << ", "
                                                 << grid_meta->voxelSize().z() << "\n"
                         << "   - eval bbox:   " << grid_meta->evalActiveVoxelBoundingBox().min().x() << ", "
                                                 << grid_meta->evalActiveVoxelBoundingBox().min().y() << ", "
                                                 << grid_meta->evalActiveVoxelBoundingBox().min().z() << " x "
                                                 << grid_meta->evalActiveVoxelBoundingBox().max().x() << ", "
                                                 << grid_meta->evalActiveVoxelBoundingBox().max().y() << ", "
                                                 << grid_meta->evalActiveVoxelBoundingBox().max().z() << "\n"
                         << "   - eval dim :   " << grid_meta->evalActiveVoxelDim().x() << ", "
                                                 << grid_meta->evalActiveVoxelDim().y() << ", "
                                                 << grid_meta->evalActiveVoxelDim().z() << "\n";
            }

            if (m_cached_input_file.m_grid_attrib_density->valueType() != "float")
            {
                ERROR_LOG << "Sparse volume importer: we presently only support scalar FP fields.";
                m_cached_input_file.clear();
                return false;
            }
            m_cached_input_file.m_file_name = filename;
        }
    }

    std::vector<mi::Float32>    density_brick_data;
    std::vector<Vec3i_struct>   density_brick_positions;

    size_t bsize[3] = {8,8,8};
    mi::Sint32 poolID = -1;
    if(m_cached_input_file.m_grid_attrib_density->isType<FloatGrid>()) {
        util::load_vdb_grid_leaves<FloatGrid>(
                            m_cached_input_file.m_grid_attrib_density,
                            density_brick_data,
                            density_brick_positions,
                            bounding_box);
        if (density_brick_positions.empty()) {
            // no brick data found intersecting the sub-region
            return true;
        }

        poolID = util::generate_attrib_brick_pool<FloatGrid>(
                            m_cached_input_file.m_grid_attrib_density,
                            density_brick_data,
                            density_brick_positions,
                            bounding_box,
                            svol_subset);
        bsize[0] = bsize[1] = bsize[2] = 8;
    }
    else if(m_cached_input_file.m_grid_attrib_density->isType<VDB_3334_grid_fp32>()) {
        util::load_vdb_grid_leaves<VDB_3334_grid_fp32>(
                            m_cached_input_file.m_grid_attrib_density,
                            density_brick_data,
                            density_brick_positions,
                            bounding_box);
        if (density_brick_positions.empty()) {
            // no brick data found intersecting the sub-region
            return true;
        }

        poolID = util::generate_attrib_brick_pool<VDB_3334_grid_fp32>(
                            m_cached_input_file.m_grid_attrib_density,
                            density_brick_data,
                            density_brick_positions,
                            bounding_box,
                            svol_subset);
        bsize[0] = bsize[1] = bsize[2] = 16;
    } else {
        ERROR_LOG << "Sparse volume importer: unknown VDB tree type, cannot load.";
        return false;
    }

    nv::index::ISparse_volume_subset_brick_pool* pool = svol_subset->get_attribute_brick_pool(poolID);
    if (pool == NULL)
    {
        ERROR_LOG << "Sparse volume importer: pool generation failed?";
        return false;
    }

    util::insert_bricks(density_brick_data, bsize, density_brick_positions, pool);

    return true;
}

} // namespace index_common
} // namespace nv

#if WIN_NT
#pragma warning(pop)
#endif
