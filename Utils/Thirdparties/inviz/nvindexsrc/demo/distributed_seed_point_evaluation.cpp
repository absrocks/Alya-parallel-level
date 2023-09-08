/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "distributed_seed_point_evaluation.h"

#include <cassert>
#include <map>

#include <nv/index/isession.h>
#include <nv/index/iregular_heightfield.h>
#include <nv/index/iregular_heightfield_compute_task.h>
#include <nv/index/iregular_volume_brick.h>

#include "common/forwarding_logger.h"
#include "volume_brick_element.h"

#include "utilities.h"

namespace // anonymous
{

//
mi::Uint8 access_amplitude_value(
    mi::neuraylib::IDice_transaction*               dice_transaction,
    const mi::math::Vector_struct<mi::Float32, 4>&  seed,
    const nv::index::IRegular_volume*               volume,
    const nv::index::IRegular_volume_data_locality* volume_data_locality,
    mi::Uint32                                      cluster_host_id)
{
    // Transform manual pick into the volume local space (ijk space).
    mi::math::Vector<mi::Float32, 3> xyz_pick_position(seed.x, seed.y, seed.z);
    mi::math::Vector<mi::Float32, 3> transformed_pick = xyz_pick_position; // transform_point(m_xyz_to_ijk, xyz_pick_position);
    // transformed_pick.x = mi::math::round(transformed_pick.x);
    // transformed_pick.y = mi::math::round(transformed_pick.y);
    // transformed_pick.z = mi::math::round(transformed_pick.z);

    const mi::Uint32 nb_local_accesses = static_cast<mi::Uint32>(volume_data_locality->get_nb_bounding_box(cluster_host_id));
    mi::math::Bbox_struct<mi::Uint32, 3> brick_ijk_box;
    mi::Uint8 amplitude_value = 128;

    for(mi::Uint32 i=0; i<nb_local_accesses; ++i)
    {
        const mi::Uint8* volume_amplitude_data = volume_data_locality->access_local_data(
            dice_transaction, cluster_host_id, i, brick_ijk_box);

        if(volume_amplitude_data)
        {
            const mi::Uint32 range_y = brick_ijk_box.max.y - brick_ijk_box.min.y;
            const mi::Uint32 range_z = brick_ijk_box.max.z - brick_ijk_box.min.z;
            const mi::Uint32 range_zy = range_z*range_y;

            if(    brick_ijk_box.min.x <= transformed_pick.x && transformed_pick.x <= brick_ijk_box.max.x
                && brick_ijk_box.min.y <= transformed_pick.y && transformed_pick.y <= brick_ijk_box.max.y)
            {
                if(mi::Sint32(brick_ijk_box.min.z) <= transformed_pick.z
                             && transformed_pick.z <= mi::Sint32(brick_ijk_box.max.z) )
                {
                    const mi::Sint32 index =
                        static_cast< mi::Sint32 >(
                            (transformed_pick.x-brick_ijk_box.min.x) * range_zy
                            + (transformed_pick.y-brick_ijk_box.min.y) * range_z
                            + (transformed_pick.z-brick_ijk_box.min.z));
                    assert(index >= 0);
                    amplitude_value = volume_amplitude_data[index];

                    break;
                }
            }
        }
    }
    return amplitude_value;
}

// class Lex_compare :
//     public std::binary_function<
//         const mi::math::Vector_struct<mi::Float32, 3>&,
//         const mi::math::Vector_struct<mi::Float32, 3>,
//         bool>
// {
// public:
//     Lex_compare(mi::Uint32 x_size)
//         :  m_x_size(x_size){}

//     bool operator()(
//         const mi::math::Vector_struct<mi::Float32, 3>& vec_a,
//         const mi::math::Vector_struct<mi::Float32, 3>& vec_b)
//     {
//         mi::Uint32 val_a = static_cast< mi::Uint32 >(vec_a.y * m_x_size + vec_a.x);
//         mi::Uint32 val_b = static_cast< mi::Uint32 >(vec_b.y * m_x_size + vec_b.x);

//         return val_a < val_b;
//     }

// private:
//     mi::Uint32 m_x_size;
// };

static inline void print_samples(
    const std::vector<mi::math::Vector_struct<mi::Float32, 4> >& seed_points)
{
    const mi::Uint32 nb_seed_points = seed_points.size();
    for(mi::Uint32 i=0; i<nb_seed_points; ++i)
    {
        INFO_LOG << "sample no. " << i << ":\t" << seed_points[i];
    }
}

// --------------------------------------------------------------------------------------------
static inline nv::index::IRegular_volume_data_locality* compute_volume_distribution_layout(
    const mi::neuraylib::Tag&               session_tag,
    const mi::neuraylib::Tag&               volume_tag,
    mi::neuraylib::IDice_transaction*       dice_transaction,
    mi::math::Bbox_struct<mi::Uint32, 2>&   heightfield_patch_bbox)
{
    // Access session
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<const nv::index::ISession>(session_tag));
    const mi::neuraylib::Tag& dist_layout_tag = session->get_distribution_layout();

    // Query the heightfield extent/bounds
    mi::base::Handle<const nv::index::IRegular_volume> volume(
        dice_transaction->access<const nv::index::IRegular_volume>(volume_tag));
    const mi::math::Bbox_struct<mi::Float32, 3> volume_bbox = volume->get_IJK_region_of_interest();

    mi::math::Bbox_struct<mi::Uint32, 3> query_bounds;
    query_bounds.min.x = heightfield_patch_bbox.min.x;
    query_bounds.min.y = heightfield_patch_bbox.min.y;
    query_bounds.min.z = static_cast< mi::Uint32 >(volume_bbox.min.z);
    query_bounds.max.x = heightfield_patch_bbox.max.x;
    query_bounds.max.y = heightfield_patch_bbox.max.y;
    query_bounds.max.z = static_cast< mi::Uint32 >(volume_bbox.max.z);

    // Access the distribution scheme
    mi::base::Handle<const nv::index::IData_distribution> distribution_layout(
        dice_transaction->access<const nv::index::IData_distribution>(dist_layout_tag));

    // Volume distribution layout
    nv::index::IRegular_volume_data_locality* data_locality = distribution_layout->retrieve_data_locality(
        volume_tag, query_bounds, dice_transaction);

    return data_locality;
}

// --------------------------------------------------------------------------------------------
static inline nv::index::IRegular_heightfield_data_locality* compute_heightfield_distribution_layout(
    const mi::neuraylib::Tag&               session_tag,
    const mi::neuraylib::Tag&               heightfield_tag,
    mi::neuraylib::IDice_transaction*       dice_transaction,
    mi::math::Bbox_struct<mi::Uint32, 2>&   heightfield_area)
{
    // Access session
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<const nv::index::ISession>(session_tag));
    const mi::neuraylib::Tag& dist_layout_tag = session->get_distribution_layout();

    // Query the heightfield extent/bounds
    mi::base::Handle<const nv::index::IRegular_heightfield> heightfield(
        dice_transaction->access<const nv::index::IRegular_heightfield>(heightfield_tag));
    const mi::math::Bbox_struct<mi::Float32, 3> heightfield_bounds = heightfield->get_IJK_bounding_box();
    heightfield_area.min.x = static_cast< mi::Uint32 >(heightfield_bounds.min.x);
    heightfield_area.min.y = static_cast< mi::Uint32 >(heightfield_bounds.min.y);
    heightfield_area.max.x = static_cast< mi::Uint32 >(heightfield_bounds.max.x);
    heightfield_area.max.y = static_cast< mi::Uint32 >(heightfield_bounds.max.y);

    // Access the distribution scheme
    mi::base::Handle<const nv::index::IData_distribution> distribution_layout(
        dice_transaction->access<const nv::index::IData_distribution>(dist_layout_tag));

    // Heightfield distribution layout
    nv::index::IRegular_heightfield_data_locality* data_locality =
        distribution_layout->retrieve_data_locality(
            heightfield_tag, heightfield_area, dice_transaction);
    return data_locality;
}

} // anonymous namespace
// --------------------------------------------------------------------------------------------

Flow_grid::Flow_grid(
    const mi::math::Bbox_struct<mi::Uint32, 2>&   range)
  : m_range(range)
{
    m_max_value = -1.f;
    m_min_value = -1.f;

    m_x_range = m_range.max.x-m_range.min.x + 2 + 1;
    m_y_range = m_range.max.y-m_range.min.y + 2 + 1;

    const mi::Uint32 size = m_x_range * m_y_range;
    m_flow_grid.resize(size, -2.f); // set the size of the grid and init each grid location with '-2', i.e., not yet processed.

    // const mi::Uint32 last_row_first_index = m_x_range*(m_y_range-1);
}

Flow_grid::~Flow_grid()
{
    // empty
}

void Flow_grid::set_border_sample(
    const mi::math::Vector_struct<mi::Float32, 4>& border) const
{
    mi::math::Bbox_struct<mi::Uint32, 2> range;
    range.min.x = (m_range.min.x==0) ? m_range.min.x : m_range.min.x-1;
    range.min.y = (m_range.min.y==0) ? m_range.min.y : m_range.min.y-1;
    range.max.x = m_range.max.x+1;
    range.max.y = m_range.max.y+1;

    if(    border.x >= range.min.x
        && border.x <= range.max.x
        && border.y >= range.min.y
        && border.y <= range.max.y )
    {
        const mi::Sint32 local_i = static_cast< mi::Sint32 >(border.x - range.min.x);
        const mi::Sint32 local_j = static_cast< mi::Sint32 >(border.y - range.min.y);

        // INFO_LOG << "BORDER TO SET: " << border << " in " << m_range << " using (" << local_i << "," << local_j << ")";

        mi::Sint32 offset_i = 0;
        if(m_range.min.x == 0)
        {
            offset_i = 1;
        }
        mi::Sint32 offset_j = 0;
        if(m_range.min.y == 0)
        {
            offset_j = 1;
        }

        const mi::Sint32 sample_position = (local_j+offset_j) * (m_x_range) + (local_i+offset_i);
        m_flow_grid[sample_position] = border.z;
    }
}

mi::Float32 Flow_grid::get_sample(
    mi::Uint32 i_value,
    mi::Uint32 j_value) const
{
    const mi::Sint32 local_i = i_value;
    const mi::Sint32 local_j = j_value;

    mi::Sint32 offset_i = 0;
    if(m_range.min.x==0)
    {
        offset_i = 1;
    }
    mi::Sint32 offset_j = 0;
    if(m_range.min.y == 0)
    {
        offset_j = 1;
    }

    mi::Sint32 sample_position = (local_j+offset_j) * (m_x_range) + (local_i+offset_i);

    if(sample_position >= (mi::Sint32)m_flow_grid.size())
    {
        DEBUG_LOG << "LOCAL: " << local_i << "," << local_j << " with offset (" << offset_i << "," << offset_i << ") at patch: " << m_range;
        return -1.f;
    }
    else
    {
        mi::Float32 sample_value = m_flow_grid[sample_position];
        if(sample_value == -2)
            return -1;
        else
            return sample_value;
    }
}

mi::Float32 Flow_grid::get_sample_value(
    mi::Uint32 i_value,
    mi::Uint32 j_value) const
{
    const mi::Uint32 local_i = i_value - m_range.min.x;
    const mi::Uint32 local_j = j_value - m_range.min.y;

    const mi::Uint32 sample_position = (local_j+1) * m_x_range + (local_i+1);

    if(sample_position >= m_flow_grid.size())
    {
        DEBUG_LOG << "get: IJ value: " << i_value << "," << j_value  << " in " << m_range;
        DEBUG_LOG << "get: IJ location: " << local_i << "," << local_j
                 << " HAS SAMPLE VALUE: " << m_flow_grid[sample_position];

        assert(false);
        return -1.f;
    }

    return m_flow_grid[sample_position];
}

/// EDITING
void Flow_grid::set_sample_value(
    mi::Uint32 i_value,
    mi::Uint32 j_value,
    mi::Float32 k_value) const
{
    const mi::Uint32 local_i = i_value - m_range.min.x;
    const mi::Uint32 local_j = j_value - m_range.min.y;

    const mi::Uint32 sample_position = (local_j+1) * m_x_range + (local_i+1);

    if(sample_position >= m_flow_grid.size())
    {
        DEBUG_LOG << "set: IJ value: " << i_value << "," << j_value  << " in " << m_range;
        DEBUG_LOG << "set: IJ location: " << local_i << "," << local_j << " TO: " << k_value;

        assert(false);
        return;
    }

    m_flow_grid[sample_position] = k_value;

    if(k_value>=0.f)
    {
        if(m_max_value==-1 || m_max_value<k_value)
            m_max_value = k_value;
        if(m_min_value==-1 || m_min_value>k_value)
            m_min_value = k_value;
    }
}

void Flow_grid::serialize(mi::neuraylib::ISerializer* serializer) const
{
    serializer->write(&m_range.min.x, 4);

    serializer->write(&m_x_range, 1);
    serializer->write(&m_y_range, 1);

    const mi::Uint32 nb_values = mi::Uint32(m_flow_grid.size());
    serializer->write(&m_flow_grid[0], nb_values);

    serializer->write(&m_max_value, 1);
    serializer->write(&m_min_value, 1);
}

void Flow_grid::deserialize(mi::neuraylib::IDeserializer* deserializer)
{
    deserializer->read(&m_range.min.x, 4);

    deserializer->read(&m_x_range, 1);
    deserializer->read(&m_y_range, 1);

    const mi::Uint32 size = m_x_range*m_y_range;
    m_flow_grid.resize(size);
    deserializer->read(&m_flow_grid[0], size);

    deserializer->read(&m_max_value, 1);
    deserializer->read(&m_min_value, 1);
}

mi::neuraylib::IElement* Flow_grid::copy() const
{
    Flow_grid* grid = new Flow_grid();
    grid->m_range     = this->m_range;
    grid->m_x_range   = this->m_x_range;
    grid->m_y_range   = this->m_y_range;

    const mi::Uint32 size = m_x_range * m_y_range;
    grid->m_flow_grid.resize(size, -2.f); // set the size of the grid and init each grid location with '-2', i.e., not yet processed.
    for(mi::Uint32 i=0; i<size; ++i)
        grid->m_flow_grid[i] = this->m_flow_grid[i];

    grid->m_max_value = this->m_max_value;
    grid->m_min_value = this->m_min_value;

    return grid;
}

// ------------------------------------------------------------------------------------------

Distributed_seed_point_evaluation::Distributed_seed_point_evaluation(
    const mi::neuraylib::Tag&                                       session_tag,
    const mi::neuraylib::Tag&                                       heightfield_tag,
    const mi::neuraylib::Tag&                                       volume_tag,
    const std::vector<mi::math::Vector_struct<mi::Float32, 3> >&    seed_points,
    const std::vector<mi::Uint32>&                                  cluster_hosts,
    const mi::math::Vector_struct<mi::Uint32, 2>&                   heightfield_size,
    bool                                                            eight_way_seedings)
 :  m_session_tag(session_tag),
    m_heightfield_tag(heightfield_tag),
    m_volume_tag(volume_tag),
    m_seed_points(0),
    m_cluster_hosts(cluster_hosts),
    m_update_heightfield(0),
    m_dangling_points(0)
{
    const mi::Uint32 nb_seeds = seed_points.size();
    for(mi::Uint32 i=0; i<nb_seeds; ++i)
    {
        mi::math::Vector_struct<mi::Float32, 4> seed;
        seed.x = seed_points[i].x;
        seed.y = seed_points[i].y;
        seed.z = seed_points[i].z;
        seed.w = -1.f;
        m_seed_points.push_back(seed);
    }
    m_samples.resize(m_cluster_hosts.size());
    m_dangling_samples.resize(m_cluster_hosts.size());
    m_border_samples.resize(m_cluster_hosts.size());
    m_flow_grid_tags.resize(m_cluster_hosts.size());
    m_volume_brick_tags.resize(m_cluster_hosts.size());
    m_patch_tag_vectors.resize(m_cluster_hosts.size());
    m_heightfield_size = heightfield_size;
    m_eight_way_seedings = eight_way_seedings ? 1 : 0;
}

Distributed_seed_point_evaluation::~Distributed_seed_point_evaluation()
{
}


// -----------------------------------------------------------------------------------------
void Distributed_seed_point_evaluation::assign_fragments_to_hosts(
    mi::Uint32* slots,
    mi::Size    nr_slots)
{
    if(nr_slots != m_cluster_hosts.size())
    {
        ERROR_LOG << "The number of host ids present ("
                  << m_cluster_hosts.size()
                  << ") and requested ("
                  << nr_slots
                  << ") doesn't match the requested count.";
        assert(m_cluster_hosts.size()==nr_slots);
    }
    else
    {
        for(mi::Size i=0; i<nr_slots; ++i)
        {
            slots[i] = m_cluster_hosts[i];
        }
    }
}

// --------------------------------------------------------------------------------------------
void Distributed_seed_point_evaluation::serialize(
    mi::neuraylib::ISerializer* serializer) const
{
    serializer->write(&m_session_tag.id, 1);
    serializer->write(&m_heightfield_tag.id, 1);
    serializer->write(&m_volume_tag.id, 1);

    const mi::Uint32 nb_seeds = mi::Uint32(m_seed_points.size());
    serializer->write(&nb_seeds, 1);
    for(mi::Uint32 i=0; i<nb_seeds; ++i)
        serializer->write(&m_seed_points[i].x, 4);

    const mi::Uint32 nb_danglings = mi::Uint32(m_dangling_points.size());
    serializer->write(&nb_danglings, 1);
    for(mi::Uint32 i=0; i<nb_danglings; ++i)
        serializer->write(&m_dangling_points[i].x, 4);

    const mi::Uint32 nb_borders = mi::Uint32(m_border_points.size());
    serializer->write(&nb_borders, 1);
    for(mi::Uint32 i=0; i<nb_borders; ++i)
        serializer->write(&m_border_points[i].x, 4);

    serializer->write(&m_update_heightfield, 1);

    const mi::Uint32 nb_tag_lists = mi::Uint32(m_flow_grid_tags.size());
    serializer->write(&nb_tag_lists, 1);
    for(mi::Uint32 i=0; i<nb_tag_lists; ++i)
    {
        const std::vector<mi::neuraylib::Tag>& flow_grid_tags = m_flow_grid_tags[i];
        const mi::Uint32 nb_tags = flow_grid_tags.size();
        serializer->write(&nb_tags, 1);
        for(mi::Uint32 j=0; j<nb_tags; ++j)
            serializer->write(&flow_grid_tags[j].id, 1);
    }

    const mi::Uint32 nb_elements = mi::Uint32(m_cluster_hosts.size());
    serializer->write(&nb_elements, 1);
    for(mi::Uint32 i=0; i<nb_elements; ++i)
        serializer->write(&m_cluster_hosts[i], 1);


    // ... and then tag lists
    const mi::Uint32 nb_lists = m_volume_brick_tags.size();
    serializer->write(&nb_lists, 1);
    for(mi::Uint32 k=0; k<nb_lists; ++k)
    {
        const std::vector<std::vector<mi::neuraylib::Tag> >& volume_brick_tags = m_volume_brick_tags[k];
        const mi::Uint32 nb_tag_lists = volume_brick_tags.size();
        serializer->write(&nb_tag_lists, 1);

        for(mi::Uint32 i=0; i<nb_tag_lists; ++i)
        {
            const std::vector<mi::neuraylib::Tag>& tags = volume_brick_tags[i];
            const mi::Uint32 nb_tags = tags.size();
            serializer->write(&nb_tags, 1);
            for(mi::Uint32 j=0; j<nb_tags; ++j)
            {
                serializer->write(&tags[j].id, 1);
            }
        }
    }

    serializer->write(&m_eight_way_seedings, 1);
}

//----------------------------------------------------------------------
void Distributed_seed_point_evaluation::deserialize(
    mi::neuraylib::IDeserializer* deserializer)
{
    deserializer->read(&m_session_tag.id, 1);
    deserializer->read(&m_heightfield_tag.id, 1);
    deserializer->read(&m_volume_tag.id, 1);

    mi::Uint32 nb_seeds = 0;
    deserializer->read(&nb_seeds, 1);
    m_seed_points.resize(nb_seeds);
    for(mi::Uint32 i=0; i<nb_seeds; ++i)
        deserializer->read(&m_seed_points[i].x, 4);

    mi::Uint32 nb_danglings = 0;
    deserializer->read(&nb_danglings, 1);
    m_dangling_points.resize(nb_danglings);
    for(mi::Uint32 i=0; i<nb_danglings; ++i)
        deserializer->read(&m_dangling_points[i].x, 4);

    mi::Uint32 nb_borders = 0;
    deserializer->read(&nb_borders, 1);
    m_border_points.resize(nb_borders);
    for(mi::Uint32 i=0; i<nb_borders; ++i)
        deserializer->read(&m_border_points[i].x, 4);

    deserializer->read(&m_update_heightfield, 1);

    mi::Uint32 nb_tag_lists = 0;
    deserializer->read(&nb_tag_lists, 1);
    m_flow_grid_tags.resize(nb_tag_lists);
    for(mi::Uint32 i=0; i<nb_tag_lists; ++i)
    {
        std::vector<mi::neuraylib::Tag>& flow_grid_tags = m_flow_grid_tags[i];
        mi::Uint32 nb_tags = 0;
        deserializer->read(&nb_tags, 1);
        flow_grid_tags.resize(nb_tags);
        for(mi::Uint32 j=0; j<nb_tags; ++j)
            deserializer->read(&flow_grid_tags[j].id, 1);
    }

    mi::Uint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_cluster_hosts.resize(nb_elements);
    for(mi::Uint32 i=0; i<nb_elements; ++i)
        deserializer->read(&m_cluster_hosts[i], 1);

    // Tag lists ...
    mi::Uint32 nb_lists = 0;
    deserializer->read(&nb_lists, 1);
    m_volume_brick_tags.resize(nb_lists);
    for(mi::Uint32 k=0; k<nb_lists; ++k)
    {
        std::vector<std::vector<mi::neuraylib::Tag> >& volume_brick_tags = m_volume_brick_tags[k];
        mi::Uint32 nb_tag_lists = 0;
        deserializer->read(&nb_tag_lists, 1);
        volume_brick_tags.resize(nb_tag_lists);
        for(mi::Uint32 i=0; i<nb_tag_lists; ++i)
        {
            std::vector<mi::neuraylib::Tag>& tags = volume_brick_tags[i];
            mi::Uint32 nb_tags = 0;
            deserializer->read(&nb_tags, 1);
            tags.resize(nb_tags);
            for(mi::Uint32 j=0; j<nb_tags; ++j)
            {
                deserializer->read(&tags[j].id, 1);
            }
        }
    }

    // Eight way ...
    deserializer->read(&m_eight_way_seedings, 1);
}

// --------------------------------------------------------------------------------------------
namespace { // anonymous

class Gate_element
{
public:
    Gate_element()
    {
        m_amplitude_value = 0;
        m_k = 0.f;
    }
    Gate_element(
        mi::Uint8   value,
        mi::Float32 k)
    {
        m_amplitude_value = value;
        m_k = k;
    }

    bool operator<(const Gate_element& rhs) const
    {
        if(this->m_amplitude_value < rhs.m_amplitude_value)
            return true;
        else if(this->m_amplitude_value > rhs.m_amplitude_value)
            return false;
        else
        {
            return ( (this->m_k*this->m_k) <= (rhs.m_k*rhs.m_k) );
        }
    }

    bool operator>(const Gate_element& rhs) const
    {
        if(this->m_amplitude_value > rhs.m_amplitude_value)
            return true;
        else if(this->m_amplitude_value < rhs.m_amplitude_value)
            return false;
        else
        {
            return ( (this->m_k*this->m_k) <= (rhs.m_k*rhs.m_k) );
        }
    }

    mi::Uint8   m_amplitude_value;
    mi::Float32 m_k;
};

// Computing height values based on the given volume dataset.
class Compute_height_value : public ICompute_height_value
{
public:
    Compute_height_value(
        mi::Uint32                                         gate_size,
        const mi::math::Bbox_struct<mi::Uint32, 2>&        patch_area,
        const mi::neuraylib::Tag&                          volume_tag,
        const nv::index::IRegular_volume*    volume,
        const nv::index::IRegular_volume_data_locality*    volume_data_locality,
        mi::Uint32                                         cluster_host_id,
        std::vector<mi::neuraylib::Tag>&                   volume_brick_tags,
        const nv::index::IDistributed_data_access_factory* volume_data_access_factory)
      : m_gate_size(gate_size),
        m_patch_area(patch_area),
        m_volume_tag(volume_tag),
        m_volume(volume),
        m_cluster_host_id(cluster_host_id),
        m_volume_data_locality(volume_data_locality),
        m_volume_brick_tags(volume_brick_tags),
        m_volume_data_access_factory(volume_data_access_factory)
    {
        m_xyz_to_ijk = m_volume->get_transform();
        m_ijk_to_xyz = m_xyz_to_ijk;
        m_ijk_to_xyz.invert();
    }

    void process_volume_brick(
        const mi::Uint8*                                volume_amplitude_data,
        const mi::math::Bbox_struct<mi::Sint32, 3>&     brick_ijk_box,
        const mi::math::Vector<mi::Float32, 3>&         transformed_pick,
        mi::Uint32                                      nb_gate_elements,
        std::vector<Gate_element>&                      gate,
        std::map<mi::Sint32, bool>&                     gate_values_found,
        mi::Sint32                                      gate_from,
        mi::Sint32                                      gate_to)
    {
        if(volume_amplitude_data)
        {
            // const mi::Uint32 range_x = brick_ijk_box.max.x - brick_ijk_box.min.x;
            const mi::Uint32 range_y = brick_ijk_box.max.y - brick_ijk_box.min.y;
            const mi::Uint32 range_z = brick_ijk_box.max.z - brick_ijk_box.min.z;
            const mi::Uint32 range_zy = range_z*range_y;

            if(    brick_ijk_box.min.x <= transformed_pick.x && transformed_pick.x <= brick_ijk_box.max.x
                && brick_ijk_box.min.y <= transformed_pick.y && transformed_pick.y <= brick_ijk_box.max.y)
            {
                const mi::Sint32 index =
                    static_cast< mi::Sint32 >((transformed_pick.x-brick_ijk_box.min.x) * range_zy
                                              + (transformed_pick.y-brick_ijk_box.min.y) * range_z);

                for(mi::Sint32 j=gate_from; j<gate_to; ++j)
                {
                    if(gate_values_found[j])
                        continue;

                    const mi::Sint32 K = static_cast< mi::Sint32 >(transformed_pick.z) + j;

                    if(         mi::Sint32(brick_ijk_box.min.z) <= K
                        && K <= mi::Sint32(brick_ijk_box.max.z) )
                    {
                        mi::Sint32 value_access = index + (K-brick_ijk_box.min.z);
                        Gate_element gate_element(volume_amplitude_data[value_access], static_cast<mi::Float32>(j));

                        // INFO_LOG << "gate element: " << gate_element.m_k << " K-brick_ijk_box.min.z: " << (K-brick_ijk_box.min.z) << " -> " << mi::Uint32(gate_element.m_amplitude_value);
                        gate.push_back(gate_element);
                        gate_values_found[j] = true;
                    }

                    if(gate.size() == nb_gate_elements)
                        return;
                }
            }
        }
    }

    //
    virtual void compute_height_value(
        mi::neuraylib::IDice_transaction*               dice_transaction,
        const mi::math::Vector_struct<mi::Float32, 4>&  input_seed,
        mi::math::Vector_struct<mi::Float32, 4>&        output_seed)
    {
        if(m_gate_size>0 && input_seed.z>=0)
        {
            mi::math::Bbox_struct<mi::Uint32, 3> volume_bbox = m_volume->get_IJK_bounding_box();

            // Transform manual pick into the volume's object space (ijk space).
            mi::math::Vector<mi::Float32, 3> xyz_pick_position(input_seed.x, input_seed.y, input_seed.z);
            mi::math::Vector<mi::Float32, 3> transformed_pick = xyz_pick_position;
            transform_point(m_xyz_to_ijk, xyz_pick_position);
            transformed_pick.x = mi::math::round(transformed_pick.x);
            transformed_pick.y = mi::math::round(transformed_pick.y);
            transformed_pick.z = mi::math::round(transformed_pick.z);

            if(    volume_bbox.min.x > transformed_pick.x
                || volume_bbox.max.x < transformed_pick.x
                || volume_bbox.min.y > transformed_pick.y
                || volume_bbox.max.y < transformed_pick.y
                || volume_bbox.min.z > transformed_pick.z
                || volume_bbox.max.z < transformed_pick.z )
            {
                INFO_LOG << "No data acquired for pick adjustment. Returning original pick position.";
                output_seed.z = -1.f;
                output_seed.w = -1.f;
                return;
            }

            mi::Sint32 gate_from = -mi::Sint32(m_gate_size);
            if(static_cast< mi::Sint32 >(transformed_pick.z)+gate_from <=
               static_cast< mi::Sint32 >(volume_bbox.min.z))
            {
                gate_from = volume_bbox.min.z - static_cast< mi::Sint32 >(transformed_pick.z);
            }

            mi::Sint32 gate_to =  mi::Sint32(m_gate_size+1);
            if(static_cast< mi::Sint32 >(transformed_pick.z)+gate_to >=
               static_cast< mi::Sint32 >(volume_bbox.max.z))
            {
                gate_to = volume_bbox.max.z - static_cast< mi::Sint32 >(transformed_pick.z);
            }

            // Number of gate elements to be gathered.
            const mi::Uint32 nb_gate_elements = mi::math::abs(gate_from) + gate_to;

            // INFO_LOG << "GATE: " << gate_from << " to " << gate_to << " = " << nb_gate_elements << " elements";

            std::map<mi::Sint32, bool> gate_values_found;
            for(mi::Sint32 j=gate_from; j<gate_to; ++j)
            {
                gate_values_found[j] = false;
            }

            std::vector<Gate_element> gate;
            mi::Uint8 k_amplitude_value = mi::Uint8(input_seed.w);

            // First accessing data stored locally on the present host avoiding any data copies and
            // volume data momory transfers
            const mi::Uint32 nb_local_accesses = static_cast<mi::Uint32>(m_volume_data_locality->get_nb_bounding_box(m_cluster_host_id));
            for(mi::Uint32 i=0; i<nb_local_accesses; ++i)
            {
                mi::math::Bbox_struct<mi::Uint32, 3> brick_ijk_box_uint;
                const mi::Uint8* volume_amplitude_data = m_volume_data_locality->access_local_data(
                    dice_transaction, m_cluster_host_id, i, brick_ijk_box_uint);
                mi::math::Bbox<mi::Sint32, 3> brick_ijk_box(
                    brick_ijk_box_uint.min.x, brick_ijk_box_uint.min.y, brick_ijk_box_uint.min.z,
                    brick_ijk_box_uint.max.x, brick_ijk_box_uint.max.y, brick_ijk_box_uint.max.z);

                process_volume_brick(
                    volume_amplitude_data, brick_ijk_box, transformed_pick,
                    nb_gate_elements, gate, gate_values_found, gate_from, gate_to);
            }

            // Second accessing data fetched through the network previously and now stored locally on the present host
            // avoiding additional volume data momory copies and data transfers
            if(gate.size()<nb_gate_elements)
            {
                const mi::Uint32 nb_tags = m_volume_brick_tags.size();
                for(mi::Uint32 i=0; i<nb_tags; ++i)
                {
                    // INFO_LOG << "Accessing volume data: tag: " << m_volume_brick_tags[i].id;
                    mi::base::Handle<const Volume_brick_element> volume_brick(
                        dice_transaction->access<const Volume_brick_element>(m_volume_brick_tags[i]));

                    const mi::math::Bbox<mi::Sint32, 3> brick_ijk_box = volume_brick->get_bounding_box();
                    const mi::Uint8* volume_amplitude_data = volume_brick->get_voxel_data();
                    // INFO_LOG << "Accessed data tag: " << m_volume_brick_tags[i].id << " for box: " << brick_ijk_box;

                    process_volume_brick(
                        volume_amplitude_data, brick_ijk_box, transformed_pick,
                        nb_gate_elements, gate, gate_values_found, gate_from, gate_to);
                }
            }

            // Third, if not all gate values have been gathered yet (i.e., some are missing) then fetch volume data
            // through the network and store them in the distributed database. That is, this operation causes both
            // network traffic and memory copies.
            //
            // Later data accesses to the fetched volume data (see 'second' above) avoid network traffic and
            // memory copies.
            if(gate.size()<nb_gate_elements)
            {
                // INFO_LOG << "gate size: " << gate.size() << " ----> elements required: " << nb_gate_elements
                //          << " pick: " << transformed_pick << " volume bbox: " << volume_bbox;
                mi::Sint32 min_missing = 0;
                mi::Sint32 max_missing = 0;
                bool missing_values = false;
                std::map<mi::Sint32, bool>::iterator map_itr = gate_values_found.begin();
                for(; map_itr!=gate_values_found.end(); ++map_itr)
                {
                    if(map_itr->second == false)
                    {
                        if(missing_values == false)
                        {
                            missing_values = true;
                            min_missing = map_itr->first;
                            max_missing = map_itr->first;
                        }
                        else
                        {
                            if(min_missing > map_itr->first)
                                min_missing = map_itr->first;

                            if(max_missing < map_itr->first)
                                max_missing = map_itr->first;
                        }
                    }
                }

                if(missing_values)
                {
                    mi::Sint32 brick_from = min_missing - mi::Sint32(m_gate_size*3); // make large enough to avoid unneeded network traffic later
                    if(static_cast< mi::Sint32 >(transformed_pick.z)+brick_from <=
                       static_cast< mi::Sint32 >(volume_bbox.min.z))
                    {
                        brick_from = volume_bbox.min.z - static_cast< mi::Sint32 >(transformed_pick.z);
                    }

                    mi::Sint32 brick_to = max_missing + mi::Sint32(m_gate_size*3); // make large enough to avoid unneeded network traffic later
                    if(static_cast< mi::Sint32 >(transformed_pick.z)+brick_to >=
                       static_cast< mi::Sint32 >(volume_bbox.max.z))
                    {
                        brick_to = volume_bbox.max.z - static_cast< mi::Sint32 >(transformed_pick.z);
                    }

                    mi::math::Bbox_struct<mi::Uint32, 3> query_bounds;
                    query_bounds.min.x = (m_patch_area.min.x>0) ? m_patch_area.min.x-1 : m_patch_area.min.x;
                    query_bounds.min.y = (m_patch_area.min.y>0) ? m_patch_area.min.y-1 : m_patch_area.min.y;
                    query_bounds.min.z = static_cast< mi::Uint32 >(transformed_pick.z) + brick_from;
                    query_bounds.max.x = m_patch_area.max.x+1;
                    query_bounds.max.y = m_patch_area.max.y+1;
                    query_bounds.max.z = static_cast< mi::Uint32 >(transformed_pick.z) + brick_to;

                    mi::base::Handle<nv::index::IRegular_volume_data_access> volume_data_access(
                        m_volume_data_access_factory->create_regular_volume_data_access(m_volume_tag));
                    volume_data_access->access(query_bounds, dice_transaction);

                    const mi::base::Handle<const nv::index::IRegular_volume_data> volume_data(
                        volume_data_access->get_volume_data());
                    const mi::base::Handle<const nv::index::IRegular_volume_data_uint8> volume_data_uint8(
                        volume_data->get_interface<const nv::index::IRegular_volume_data_uint8>());

                    if (!volume_data_uint8) {
                        ERROR_LOG << "Failed to access volume data: data access on non-uint8 volume type.";
                        return;
                    }

                    const mi::math::Bbox_struct<mi::Uint32, 3>& brick_ijk_box_uint = volume_data_access->get_bounding_box();
                    mi::math::Bbox<mi::Sint32, 3> brick_ijk_box(
                        brick_ijk_box_uint.min.x, brick_ijk_box_uint.min.y, brick_ijk_box_uint.min.z,
                        brick_ijk_box_uint.max.x, brick_ijk_box_uint.max.y, brick_ijk_box_uint.max.z);
                    
                    const mi::Uint8* volume_amplitude_data = volume_data_uint8->get_voxel_data();

                    mi::base::Handle< Volume_brick_element > volume_brick
                        (new Volume_brick_element(brick_ijk_box, volume_amplitude_data));
                    // ... gate value computation again ....
                    process_volume_brick(
                        volume_brick->get_voxel_data(), volume_brick->get_bounding_box(), transformed_pick,
                        nb_gate_elements, gate, gate_values_found, gate_from, gate_to); // with original gate of course

                    // Storing volume brick data in distributed database
                    const mi::neuraylib::Tag volume_brick_tag = dice_transaction->store(volume_brick.get());
                    // INFO_LOG << "CREATED and STORED data tag: " << volume_brick_tag.id << " for box: " << brick_ijk_box;
                    m_volume_brick_tags.push_back(volume_brick_tag); // .. also available for trace gate in the neighorhood of the seed
                }
            }

            // Early out
            if(gate.empty())
            {
                DEBUG_LOG << "No data acquired for pick adjustment. Returning original pick position.";
                output_seed.z = -1.f;
                output_seed.w = -1.f;
                return;
            }

            Gate_element current_reference = gate[0];
            if(mi::Uint32(k_amplitude_value)>127)
            {
                for(mi::Uint32 i=0; i<gate.size(); ++i)
                {
                    if(current_reference<gate[i])
                        current_reference = gate[i];
                }
            }
            else
            {
                for(mi::Uint32 i=0; i<gate.size(); ++i)
                {
                    if(current_reference>gate[i])
                        current_reference = gate[i];
                }
            }

            // Transform new pick location back to xyz space and retrun
            transformed_pick.z += current_reference.m_k;
            xyz_pick_position = transform_point(m_ijk_to_xyz, transformed_pick);

            output_seed.z = xyz_pick_position.z;
            output_seed.w = current_reference.m_amplitude_value;
        }
        else
        {
            DEBUG_LOG << "TERMINATED EARLY. Input seed has the following amplitude value:" << input_seed.w;
            output_seed.z = -1.f;   // hole
            output_seed.w = -1.f;   // nothing valid
        }
    }

    const std::vector<mi::neuraylib::Tag>& get_tags() const
    {
        return m_volume_brick_tags;
    }

private:
    mi::Uint32                                              m_gate_size;
    mi::math::Bbox_struct<mi::Uint32, 2>                    m_patch_area;
    mi::neuraylib::Tag                                      m_volume_tag;
    const nv::index::IRegular_volume*                       m_volume;
    mi::Uint32                                              m_cluster_host_id;
    const nv::index::IRegular_volume_data_locality*         m_volume_data_locality;
    std::vector<mi::neuraylib::Tag>                         m_volume_brick_tags;
    const nv::index::IDistributed_data_access_factory*      m_volume_data_access_factory;

    mi::math::Matrix<mi::Float32, 4, 4>                     m_xyz_to_ijk;
    mi::math::Matrix<mi::Float32, 4, 4>                     m_ijk_to_xyz;
};

// SIMPLE TEST CALL TO VERIFY THE SEEDING PROCESS!! (used in evaluate_single_seed_point(..) )
#if 0
mi::Float32 compute_height_value(mi::Float32  height)
{
    return height + 0.2f;
}
#endif

//
void evaluate_single_dangling_point(
    mi::neuraylib::IDice_transaction*                       dice_transaction,
    const mi::math::Vector_struct<mi::Float32, 4>&          input_dangling_point,
    const Flow_grid*                                        flow_grid,
    ICompute_height_value*                                  height_compute,
    mi::math::Vector_struct<mi::Float32, 4>&                output_point,
    std::vector<mi::math::Vector_struct<mi::Float32, 4> >&  result_samples)
{
    const mi::Float32 sample_value =
        flow_grid->get_sample_value(static_cast< mi::Uint32 >(input_dangling_point.x),
                                    static_cast< mi::Uint32 >(input_dangling_point.y));
    if(sample_value == -2)
    {
        mi::math::Vector_struct<mi::Float32, 4> output_sample = input_dangling_point;
        // output_sample.z = compute_height_value(input_seed.z);
        height_compute->compute_height_value(dice_transaction, input_dangling_point, output_sample);
        flow_grid->set_sample_value(static_cast< mi::Uint32 >(output_sample.x),
                                    static_cast< mi::Uint32 >(output_sample.y),
                                    output_sample.z);
        result_samples.push_back(output_sample);

        output_point = output_sample;
    }
}

//
void evaluate_single_seed_point(
    mi::neuraylib::IDice_transaction*                       dice_transaction,
    const mi::math::Vector_struct<mi::Float32, 4>&          input_seed,
    const Flow_grid*                                        flow_grid,
    ICompute_height_value*                                  height_compute,
    std::vector<mi::math::Vector_struct<mi::Float32, 4> >&  result_samples)
{
    const mi::Float32 sample_value =
        flow_grid->get_sample_value(static_cast< mi::Uint32 >(input_seed.x),
                                    static_cast< mi::Uint32 >(input_seed.y));
    if(sample_value == -2)
    {
        mi::math::Vector_struct<mi::Float32, 4> output_sample = input_seed;
        // output_sample.z = compute_height_value(input_seed.z);
        height_compute->compute_height_value(dice_transaction, input_seed, output_sample);

        flow_grid->set_sample_value(static_cast< mi::Uint32 >(output_sample.x),
                                    static_cast< mi::Uint32 >(output_sample.y),
                                    output_sample.z);
        result_samples.push_back(output_sample);
    }
}


void compute_single_seed_flow(
    mi::neuraylib::IDice_transaction*                               dice_transaction,
    bool                                                            eight_way_seedings,
    const mi::math::Vector_struct<mi::Float32, 4>&                  seed,
    const mi::math::Bbox_struct<mi::Uint32, 2>&                     patch,
    const Flow_grid*                                                flow_grid,
    ICompute_height_value*                                          height_compute,
    std::vector<mi::math::Vector_struct<mi::Float32, 4> >&          result_samples,
    std::vector<mi::math::Vector_struct<mi::Float32, 4> >&          result_dangling_samples)
{
    std::vector<mi::math::Vector_struct<mi::Float32, 4> > potential_seeds;
    // 4-way seeding
    mi::math::Vector_struct<mi::Float32, 4> top   = seed; top.y   += 1.f; potential_seeds.push_back(top);
    mi::math::Vector_struct<mi::Float32, 4> left  = seed; left.x  -= 1.f; potential_seeds.push_back(left);
    mi::math::Vector_struct<mi::Float32, 4> down  = seed; down.y  -= 1.f; potential_seeds.push_back(down);
    mi::math::Vector_struct<mi::Float32, 4> right = seed; right.x += 1.f; potential_seeds.push_back(right);

    if(eight_way_seedings==1) // 8-way seeding
    {
        mi::math::Vector_struct<mi::Float32, 4> top_left   = seed; top_left.x   -= 1.f; top_left.y   += 1.f; potential_seeds.push_back(top_left);
        mi::math::Vector_struct<mi::Float32, 4> down_left  = seed; down_left.x  -= 1.f; down_left.y  -= 1.f; potential_seeds.push_back(down_left);
        mi::math::Vector_struct<mi::Float32, 4> down_right = seed; down_right.x += 1.f; down_right.y -= 1.f; potential_seeds.push_back(down_right);
        mi::math::Vector_struct<mi::Float32, 4> top_right  = seed; top_right.x  += 1.f; top_right.y  += 1.f; potential_seeds.push_back(top_right);
    }

    mi::Uint32 nb_potentials = potential_seeds.size();
    // mi::Uint32 nb_resulting_seeds = 0;
    for(mi::Uint32 i=0; i<nb_potentials; ++i)
    {
        const mi::math::Vector_struct<mi::Float32, 4>& potential = potential_seeds[i];

        // INFO_LOG << "Seed " << seed << " potential -> " << potential;

        if(    patch.min.x <= potential.x && potential.x <= patch.max.x
            && patch.min.y <= potential.y && potential.y <= patch.max.y)
        {
            evaluate_single_seed_point(
                dice_transaction,
                potential,
                flow_grid,
                height_compute,
                result_samples);
        }
        else
        {
            result_dangling_samples.push_back(potential);
        }
    }
}

} // anonymous namespace

// --------------------------------------------------------------------------------------------
void Distributed_seed_point_evaluation::evaluate_seed_flow(
    mi::neuraylib::IDice_transaction*                       dice_transaction,
    const mi::math::Bbox_struct<mi::Uint32, 2>&             patch,
    const Flow_grid*                                        flow_grid,
    ICompute_height_value*                                  height_compute,
    std::vector<mi::math::Vector_struct<mi::Float32, 4> >&  result_samples,
    std::vector<mi::math::Vector_struct<mi::Float32, 4> >&  result_dangling_samples) const
{
    const mi::Uint32 nb_seeds = m_seed_points.size();
    for(mi::Uint32 i=0; i<nb_seeds; ++i)
    {
        const mi::math::Vector_struct<mi::Float32, 4>& seed = m_seed_points[i];
        if(    patch.min.x <= seed.x && seed.x <= patch.max.x
            && patch.min.y <= seed.y && seed.y <= patch.max.y)
        {
            compute_single_seed_flow(
                dice_transaction,
                m_eight_way_seedings != 0,
                seed,                   // seed value
                patch,
                flow_grid,
                height_compute,
                result_samples,
                result_dangling_samples);
        }
    }
}

// --------------------------------------------------------------------------------------------
void Distributed_seed_point_evaluation::evaluate_seed_flow(
    mi::neuraylib::IDice_transaction*                      dice_transaction,
    nv::index::IRegular_heightfield_data_locality*         heightfield_locality,
    nv::index::IRegular_volume_data_locality*              volume_data_locality,
    mi::Uint32                                             host_id,
    std::vector<mi::neuraylib::Tag>&                       flow_grid_tags,
    std::vector<std::vector<mi::neuraylib::Tag> >&         volume_brick_tags,
    std::vector<mi::math::Vector_struct<mi::Float32, 4> >& result_samples,
    std::vector<mi::math::Vector_struct<mi::Float32, 4> >& result_dangling_samples,
    std::vector<mi::math::Vector_struct<mi::Float32, 4> >& result_border_samples)
{
    // Access session and regular volume data access factory
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<const nv::index::ISession>(m_session_tag));

    const mi::neuraylib::Tag& access_factory_tag = session->get_data_access_factory();

    mi::base::Handle<const nv::index::IDistributed_data_access_factory> access_factory(
        dice_transaction->access<const nv::index::IDistributed_data_access_factory>(access_factory_tag));

    mi::base::Handle<const nv::index::IRegular_volume> volume(
        dice_transaction->access<const nv::index::IRegular_volume>(m_volume_tag));

    // Initialize the seeding process per host and for all assigned patches
    // Done only once!
    const mi::Uint32 nb_patch_edits = static_cast<mi::Uint32>(heightfield_locality->get_nb_bounding_box(host_id));
    if(flow_grid_tags.empty())
    {
        // initialize the tag list per patch
        volume_brick_tags.resize(nb_patch_edits);

        // initialize the flow grids per patch
        for(mi::Uint32 i=0; i<nb_patch_edits; ++i)
        {
            const mi::math::Bbox<mi::Sint32, 3> patch_3D_bbox = heightfield_locality->get_bounding_box(host_id, i);
            mi::math::Bbox<mi::Uint32, 2> patch; // x and y are always non-negative
            patch.min.x = static_cast<mi::Uint32>(patch_3D_bbox.min.x);
            patch.min.y = static_cast<mi::Uint32>(patch_3D_bbox.min.y);
            patch.max.x = static_cast<mi::Uint32>(patch_3D_bbox.max.x);
            patch.max.y = static_cast<mi::Uint32>(patch_3D_bbox.max.y);

            mi::base::Handle< Flow_grid >flow_grid(new Flow_grid(patch));

            // Initialize the flow grid with the initial seed points
            const mi::Uint32 nb_seeds = m_seed_points.size();
            for(mi::Uint32 j=0; j<nb_seeds; ++j)
            {
                mi::math::Vector_struct<mi::Float32, 4>& seed = m_seed_points[j];
                if(    patch.min.x <= seed.x && seed.x <= patch.max.x
                    && patch.min.y <= seed.y && seed.y <= patch.max.y)
                {
                    const mi::Uint8 amplitude_value = access_amplitude_value(
                        dice_transaction,
                        seed,
                        volume.get(),
                        volume_data_locality,
                        host_id);
                    flow_grid->set_sample_value(static_cast< mi::Uint32 >(seed.x),
                                                static_cast< mi::Uint32 >(seed.y),
                                                seed.z);
                    seed.w = amplitude_value;
                }
            }

            const mi::neuraylib::Tag flow_grid_tag = dice_transaction->store(flow_grid.get());
            flow_grid_tags.push_back(flow_grid_tag);
        }
    }

    // Invoke the seeding process ........
    for(mi::Uint32 i=0; i<nb_patch_edits; ++i)
    {
        // Flow grid storing the elevation values at each ij position but also steers the seeding (flow) process
        mi::base::Handle<const Flow_grid> flow_grid(dice_transaction->access<const Flow_grid>(flow_grid_tags[i]));
        const mi::math::Bbox<mi::Sint32, 3> patch_3D_bbox = heightfield_locality->get_bounding_box(host_id, i);
        mi::math::Bbox<mi::Uint32, 2> patch; // x and y are always non-negative
        patch.min.x = static_cast<mi::Uint32>(patch_3D_bbox.min.x);
        patch.min.y = static_cast<mi::Uint32>(patch_3D_bbox.min.y);
        patch.max.x = static_cast<mi::Uint32>(patch_3D_bbox.max.x);
        patch.max.y = static_cast<mi::Uint32>(patch_3D_bbox.max.y);

        // Technique that computes the k position for each neighboring trace
        Compute_height_value height_compute(
            10,
            patch,
            m_volume_tag,
            volume.get(),
            volume_data_locality,
            host_id,
            volume_brick_tags[i],
            access_factory.get());

        // Handle dangling samples first.
        const mi::Uint32 nb_dangling_samples = m_dangling_points.size();
        for(mi::Uint32 j=0; j<nb_dangling_samples; ++j)
        {
            // A dangling sample can only be on the border
            const mi::math::Vector_struct<mi::Float32, 4>& dangling = m_dangling_points[j];
            if(    patch.min.x <= dangling.x && dangling.x <= patch.max.x
                && patch.min.y <= dangling.y && dangling.y <= patch.max.y)
            {
                mi::math::Vector_struct<mi::Float32, 4> border_point = dangling;
                evaluate_single_dangling_point(
                    dice_transaction,
                    dangling,
                    flow_grid.get(),
                    &height_compute,
                    border_point,
                    result_samples);

                result_border_samples.push_back(border_point);
            }
        }

        const mi::Uint32 nb_border_samples = m_border_points.size();
        for(mi::Uint32 j=0; j<nb_border_samples; ++j)
        {
            // A dangling sample can only be on the border
            const mi::math::Vector_struct<mi::Float32, 4>& border = m_border_points[j];
            flow_grid->set_border_sample(border);
        }

        // Handle all seed points.
        evaluate_seed_flow(
            dice_transaction,
            patch,
            flow_grid.get(),
            &height_compute,
            result_samples,
            result_dangling_samples);

        const std::vector<mi::neuraylib::Tag>& tag_list = height_compute.get_tags();
        const mi::Uint32 nb_tags = tag_list.size();
        if(nb_tags>volume_brick_tags[i].size() )
        {
            volume_brick_tags[i].clear();
            for(mi::Uint32 k=0; k<nb_tags; ++k)
            {
                volume_brick_tags[i].push_back(tag_list[k]);
            }
        }
    }
}

// --------------------------------------------------------------------------------------------
// Parallel compute on the local cluster host.
void Distributed_seed_point_evaluation::update_seed_points()
{
    m_seed_points.clear();
    const mi::Uint32 nb_seed_sets = m_samples.size();
    for(mi::Uint32 i=0; i<nb_seed_sets; ++i)
    {
        m_seed_points.insert(m_seed_points.end(), m_samples[i].begin(), m_samples[i].end());
        m_samples[i].clear();
    }
    // std::sort(m_seed_points.begin(), m_seed_points.end(), Lex_compare(m_heightfield_size.x));
    // INFO_LOG << "UPDATED SEED POINTS";
    // INFO_LOG << "===================";
    // print_samples(m_seed_points);

    m_dangling_points.clear();
    const mi::Uint32 nb_dangling_sets = m_dangling_samples.size();
    for(mi::Uint32 i=0; i<nb_dangling_sets; ++i)
    {
        m_dangling_points.insert(m_dangling_points.end(), m_dangling_samples[i].begin(), m_dangling_samples[i].end());
        m_dangling_samples[i].clear();
    }
    // std::sort(m_dangling_points.begin(), m_dangling_points.end(), Lex_compare(m_heightfield_size.x));
    // INFO_LOG << "UPDATED DANGLING POINTS";
    // INFO_LOG << "=======================";
    // print_samples(m_dangling_points);

    m_border_points.clear();
    const mi::Uint32 nb_border_sets = m_border_samples.size();
    for(mi::Uint32 i=0; i<nb_border_sets; ++i)
    {
        m_border_points.insert(m_border_points.end(), m_border_samples[i].begin(), m_border_samples[i].end());
        m_border_samples[i].clear();
    }
    // std::sort(m_border_points.begin(), m_border_points.end(), Lex_compare(m_heightfield_size.x));
    // INFO_LOG << "UPDATED BORDER POINTS";
    // INFO_LOG << "=====================";
    // print_samples(m_border_points);
}

mi::Uint32 Distributed_seed_point_evaluation::get_nb_seed_points() const
{
    return m_seed_points.size()+m_dangling_points.size()+m_border_points.size();
}

// --------------------------------------------------------------------------------------------
// Parallel compute on the local cluster host.
void Distributed_seed_point_evaluation::execute_fragment(
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Size                                        index,
    mi::Size                                        count,
    const mi::neuraylib::IJob_execution_context*    context)
{
    const mi::Uint32 cluster_host = m_cluster_hosts[index];
    mi::math::Bbox_struct<mi::Uint32, 2> heightfield_area;
    mi::base::Handle<nv::index::IRegular_heightfield_data_locality> heightfield_data_locality(
        compute_heightfield_distribution_layout(m_session_tag, m_heightfield_tag, dice_transaction, heightfield_area));

    mi::base::Handle<nv::index::IRegular_volume_data_locality> volume_data_locality(
        compute_volume_distribution_layout(m_session_tag, m_volume_tag, dice_transaction, heightfield_area));

    evaluate_seed_flow(
        dice_transaction,
        heightfield_data_locality.get(),
        volume_data_locality.get(),
        cluster_host,
        m_flow_grid_tags[index],
        m_volume_brick_tags[index],
        m_samples[index],
        m_dangling_samples[index],
        m_border_samples[index]);

    //
    if(m_update_heightfield==1) // minor optimization
    {
        const std::vector<mi::neuraylib::Tag>& flow_grid_tags = m_flow_grid_tags[index];

        const mi::Uint32 nb_patch_edits = static_cast<mi::Uint32>(heightfield_data_locality->get_nb_bounding_box(cluster_host));
        std::vector<std::pair<mi::math::Bbox_struct<mi::Uint32, 2>, mi::neuraylib::Tag> >& patch_tag_vector = m_patch_tag_vectors[index];
        patch_tag_vector.clear();

        for(mi::Uint32 i=0; i<nb_patch_edits; ++i)
        {
            // Flow grid edit needs to be triggered ;(
            mi::base::Handle<const Flow_grid> flow_grid(dice_transaction->edit<const Flow_grid>(flow_grid_tags[i]));

            const mi::math::Bbox<mi::Sint32, 3> patch_3D_bbox = heightfield_data_locality->get_bounding_box(cluster_host, i);
            mi::math::Bbox_struct<mi::Uint32, 2> patch; // x and y are always non-negative
            patch.min.x = static_cast<mi::Uint32>(patch_3D_bbox.min.x);
            patch.min.y = static_cast<mi::Uint32>(patch_3D_bbox.min.y);
            patch.max.x = static_cast<mi::Uint32>(patch_3D_bbox.max.x);
            patch.max.y = static_cast<mi::Uint32>(patch_3D_bbox.max.y);

            std::pair<mi::math::Bbox_struct<mi::Uint32, 2>, mi::neuraylib::Tag> pair(patch, flow_grid_tags[i]);
            patch_tag_vector.push_back(pair);
        }
    }
}

// Execute one fragment of the job on a remote host. This is executed on a different host than
// calling host and given to the receive_remote_result function.
void Distributed_seed_point_evaluation::execute_fragment_remote(
    mi::neuraylib::ISerializer*                     serializer,
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Size                                        index,
    mi::Size                                        count,
    const mi::neuraylib::IJob_execution_context*    context)
{
    const mi::Uint32 cluster_host = m_cluster_hosts[index];
    mi::math::Bbox_struct<mi::Uint32, 2> heightfield_area;
    mi::base::Handle<nv::index::IRegular_heightfield_data_locality> heightfield_data_locality(
        compute_heightfield_distribution_layout(m_session_tag, m_heightfield_tag, dice_transaction, heightfield_area));

    mi::base::Handle<nv::index::IRegular_volume_data_locality> volume_data_locality(
        compute_volume_distribution_layout(m_session_tag, m_volume_tag, dice_transaction, heightfield_area));

    // Tag referencing the flow grids
    std::vector<mi::neuraylib::Tag>& flow_grid_tags = m_flow_grid_tags[index];
    // Tags referening the volume subcubes
    std::vector<std::vector<mi::neuraylib::Tag> >& volume_brick_tags = m_volume_brick_tags[index];

    std::vector<mi::math::Vector_struct<mi::Float32, 4> > samples;
    std::vector<mi::math::Vector_struct<mi::Float32, 4> > dangling_samples;
    std::vector<mi::math::Vector_struct<mi::Float32, 4> > border_samples;

    evaluate_seed_flow(
        dice_transaction,
        heightfield_data_locality.get(),
        volume_data_locality.get(),
        cluster_host,
        flow_grid_tags,
        volume_brick_tags,
        samples,
        dangling_samples,
        border_samples);

    // Serialize results
    // .... patch-tag list first is required
    if(m_update_heightfield==1)  // minor optimization
    {
        const mi::Uint32 nb_patch_edits = static_cast<mi::Uint32>(heightfield_data_locality->get_nb_bounding_box(cluster_host));

        serializer->write(&nb_patch_edits, 1);
        for(mi::Uint32 i=0; i<nb_patch_edits; ++i)
        {
            // Flow grid edit needs to be triggered ;(
            mi::base::Handle<const Flow_grid> flow_grid(dice_transaction->edit<const Flow_grid>(flow_grid_tags[i]));

            const mi::math::Bbox<mi::Sint32, 3> patch_3D_bbox = heightfield_data_locality->get_bounding_box(cluster_host, i);
            mi::math::Bbox<mi::Uint32, 2> patch; // x and y are always non-negative
            patch.min.x = static_cast<mi::Uint32>(patch_3D_bbox.min.x);
            patch.min.y = static_cast<mi::Uint32>(patch_3D_bbox.min.y);
            patch.max.x = static_cast<mi::Uint32>(patch_3D_bbox.max.x);
            patch.max.y = static_cast<mi::Uint32>(patch_3D_bbox.max.y);

            serializer->write(&patch.min.x, 2);
            serializer->write(&patch.max.x, 2);
            serializer->write(&flow_grid_tags[i].id, 1);
        }
    }

    // ... then the actual samples
    const mi::Uint32 nb_samples = mi::Uint32(samples.size());
    serializer->write(&nb_samples, 1);
    for(mi::Uint32 i=0; i<nb_samples; ++i)
        serializer->write(&samples[i].x, 4);

    // ... then dangling samples
    const mi::Uint32 nb_dangling_samples = mi::Uint32(dangling_samples.size());
    serializer->write(&nb_dangling_samples, 1);
    for(mi::Uint32 i=0; i<nb_dangling_samples; ++i)
        serializer->write(&dangling_samples[i].x, 4);

    // ... then border samples
    const mi::Uint32 nb_border_samples = mi::Uint32(border_samples.size());
    serializer->write(&nb_border_samples, 1);
    for(mi::Uint32 i=0; i<nb_border_samples; ++i)
        serializer->write(&border_samples[i].x, 4);

    // ... and then flow grid tags
    const mi::Uint32 nb_flow_grid_tags = mi::Uint32(flow_grid_tags.size());
    serializer->write(&nb_flow_grid_tags, 1);
    for(mi::Uint32 i=0; i<nb_flow_grid_tags; ++i)
        serializer->write(&flow_grid_tags[i].id, 1);

    // ... and then tag lists
    const mi::Uint32 nb_tag_lists = volume_brick_tags.size();
    serializer->write(&nb_tag_lists, 1);
    for(mi::Uint32 i=0; i<nb_tag_lists; ++i)
    {
        const std::vector<mi::neuraylib::Tag>& tags = volume_brick_tags[i];
        const mi::Uint32 nb_tags = tags.size();
        serializer->write(&nb_tags, 1);
        for(mi::Uint32 j=0; j<nb_tags; ++j)
        {
            serializer->write(&tags[j].id, 1);
        }
    }
}

// Receive a job result from a remote host and integrate it in the end result. This function
// the same effect as the call of execute_fragment on the original host.
void Distributed_seed_point_evaluation::receive_remote_result(
    mi::neuraylib::IDeserializer*                   deserializer,
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Size                                        index,
    mi::Size                                        count)
{
    // Deserialize results
    // .... patch-tag list first is required
    if(m_update_heightfield==1)
    {
        std::vector<std::pair<mi::math::Bbox_struct<mi::Uint32, 2>, mi::neuraylib::Tag> >& patch_tag_vector = m_patch_tag_vectors[index];
        patch_tag_vector.clear();

        mi::Uint32 nb_patch_edits = 0;
        deserializer->read(&nb_patch_edits, 1);

        for(mi::Uint32 i=0; i<nb_patch_edits; ++i)
        {
            std::pair<mi::math::Bbox_struct<mi::Uint32, 2>, mi::neuraylib::Tag> pair;
            deserializer->read(&pair.first.min.x, 2);
            deserializer->read(&pair.first.max.x, 2);
            deserializer->read(&pair.second.id, 1);

            patch_tag_vector.push_back(pair);
        }
    }


    // ... then the actual samples
    std::vector<mi::math::Vector_struct<mi::Float32, 4> >& samples = m_samples[index];
    samples.clear();
    mi::Uint32 nb_samples = 0;
    deserializer->read(&nb_samples, 1);
    samples.resize(nb_samples);
    for(mi::Uint32 i=0; i<nb_samples; ++i)
        deserializer->read(&samples[i].x, 4);

    // ... then dangling samples
    std::vector<mi::math::Vector_struct<mi::Float32, 4> >& dangling_samples = m_dangling_samples[index];
    mi::Uint32 nb_dangling_samples = 0;
    deserializer->read(&nb_dangling_samples, 1);
    dangling_samples.resize(nb_dangling_samples);
    for(mi::Uint32 i=0; i<nb_dangling_samples; ++i)
        deserializer->read(&dangling_samples[i].x, 4);

    // ... then border samples
    std::vector<mi::math::Vector_struct<mi::Float32, 4> >& border_samples = m_border_samples[index];
    mi::Uint32 nb_border_samples = 0;
    deserializer->read(&nb_border_samples, 1);
    border_samples.resize(nb_border_samples);
    for(mi::Uint32 i=0; i<nb_border_samples; ++i)
        deserializer->read(&border_samples[i].x, 4);

    // ... and then flow grid tags
    std::vector<mi::neuraylib::Tag>& flow_grid_tags = m_flow_grid_tags[index];
    flow_grid_tags.clear();
    mi::Uint32 nb_flow_grid_tags = 0;
    deserializer->read(&nb_flow_grid_tags, 1);
    flow_grid_tags.resize(nb_flow_grid_tags);
    for(mi::Uint32 i=0; i<nb_flow_grid_tags; ++i)
        deserializer->read(&flow_grid_tags[i].id, 1);

    // ... and then the tag lists
    std::vector<std::vector<mi::neuraylib::Tag> >& volume_brick_tags = m_volume_brick_tags[index];
    volume_brick_tags.clear();
    mi::Uint32 nb_tag_lists = 0;
    deserializer->read(&nb_tag_lists, 1);
    for(mi::Uint32 i=0; i<nb_tag_lists; ++i)
    {
        std::vector<mi::neuraylib::Tag> tags;
        mi::Uint32 nb_tags = 0;
        deserializer->read(&nb_tags, 1);
        tags.resize(nb_tags);
        for(mi::Uint32 j=0; j<nb_tags; ++j)
        {
            deserializer->read(&tags[j].id, 1);
        }
        volume_brick_tags.push_back(tags);
    }
}

void Distributed_seed_point_evaluation::clean_up(
    mi::neuraylib::IDice_transaction*   dice_transaction)
{
    // Clean up the flow grid tags first ...
    const mi::Uint32 nb_tag_lists = m_flow_grid_tags.size();
    for(mi::Uint32 i=0; i<nb_tag_lists; ++i)
    {
        std::vector<mi::neuraylib::Tag>& flow_grid_tags = m_flow_grid_tags[i];
        const mi::Uint32 nb_tags = flow_grid_tags.size();

        for(mi::Uint32 j=0; j<nb_tags; ++j)
        {
            dice_transaction->remove(flow_grid_tags[j]);
        }

        flow_grid_tags.clear();
    }

    // ... and then volume tags.
    const mi::Uint32 nb_lists = m_volume_brick_tags.size();
    for(mi::Uint32 k=0; k<nb_lists; ++k)
    {
        std::vector<std::vector<mi::neuraylib::Tag> >& volume_brick_tags = m_volume_brick_tags[k];
        const mi::Uint32 nb_tag_lists = volume_brick_tags.size();

        for(mi::Uint32 i=0; i<nb_tag_lists; ++i)
        {
            std::vector<mi::neuraylib::Tag>& tags = volume_brick_tags[i];
            const mi::Uint32 nb_tags = tags.size();

            for(mi::Uint32 j=0; j<nb_tags; ++j)
            {
                dice_transaction->remove(tags[j]);
            }
            tags.clear();
        }
        volume_brick_tags.clear();
    }
    m_volume_brick_tags.clear();
}


