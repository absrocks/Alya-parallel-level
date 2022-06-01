/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "distributed_heightfield_seeding_updates.h"

#include <cassert>

#include <nv/index/isession.h>
#include <nv/index/iregular_heightfield.h>
#include <nv/index/iregular_heightfield_compute_task.h>

#include "distributed_seed_point_evaluation.h"
#include "utilities.h"

#include "common/forwarding_logger.h"
#include "common/type_conversion_utility.h"


//----------------------------------------------------------------------
namespace // anonymous namespace
{
class Heightfield_elevation_update_operation :
    public mi::base::Interface_implement<nv::index::IRegular_heightfield_compute_task>
{
public:
    Heightfield_elevation_update_operation(
        const mi::math::Bbox_struct<mi::Sint32, 2>&     roi,
        mi::neuraylib::Tag                              flow_grid_tag)
      : m_roi(roi),
        m_flow_grid_tag(flow_grid_tag)
    {
        m_min_max_heights.x = -1;
        m_min_max_heights.y = -1;
    }

    void get_min_max_heights(
        mi::math::Vector_struct<mi::Float32, 2>& min_max_heights) const
    {
        if(min_max_heights.x > m_min_max_heights.x)
            min_max_heights.x = m_min_max_heights.x;
        if(min_max_heights.y < m_min_max_heights.y)
            min_max_heights.y = m_min_max_heights.y;
    }

public:
    // Implemented IRegular_heightfield_compute_task

    virtual void get_region_of_interest_for_compute(
        mi::math::Bbox_struct<mi::Sint32, 2>& roi) const
    {
        roi = m_roi;
    }

    // FIXME: Updated to bounding box access, but not tested
    virtual bool compute(
        const mi::math::Bbox_struct<mi::Sint32, 2>& ij_patch_bbox,
        mi::Float32*                                elevation_values,
        mi::neuraylib::IDice_transaction*           dice_transaction) const
    {
        if (m_flow_grid_tag != mi::neuraylib::NULL_TAG)
        {
            mi::base::Handle<const Flow_grid> flow_grid(dice_transaction->access<const Flow_grid>(m_flow_grid_tag));

            mi::Float32* dst = &elevation_values[0];
            for (mi::Sint32 j = ij_patch_bbox.min.y; j < ij_patch_bbox.max.y; ++j)
            {
                for (mi::Sint32 i = ij_patch_bbox.min.x; i < ij_patch_bbox.max.x; ++i)
                {
                    const mi::Float32 value = flow_grid->get_sample(i-ij_patch_bbox.min.x, j-ij_patch_bbox.min.y);
                    *dst++ = value;

                    if (value < m_min_max_heights.x)
                    {
                        m_min_max_heights.x = value;
                    }

                    if (value > m_min_max_heights.y)
                    {
                        m_min_max_heights.y = value;
                    }
                }
            }
            return true;
        }
        else
        {
            return false;
        }
    }

    virtual bool compute(
        const mi::math::Bbox_struct<mi::Sint32, 2>&     ij_patch_bbox,
        mi::Float32*                                    elevation_values,
        mi::math::Vector_struct<mi::Float32, 3>*        normal_values,
        mi::neuraylib::IDice_transaction*               dice_transaction) const
    {
        return false;
    }

    // This override is not necessary, but without this override, ICC
    // compiler warns.
    virtual bool user_defined_normal_computation() const { return false; }

private:
    mi::math::Bbox<mi::Sint32, 2> m_roi;
    mi::neuraylib::Tag            m_flow_grid_tag;

    mutable mi::math::Vector_struct<mi::Float32, 2>  m_min_max_heights;
};

} // anonymous namespace
//----------------------------------------------------------------------

Distributed_heightfield_seedings_updates::Distributed_heightfield_seedings_updates(
    const mi::neuraylib::Tag&                                       session_tag,
    const mi::neuraylib::Tag&                                       heightfield_tag,
    const std::vector<mi::Uint32>&                                  cluster_hosts)
 :  m_session_tag(session_tag),
    m_heightfield_tag(heightfield_tag),
    m_cluster_hosts(cluster_hosts),
    m_max_value(0.0f),
    m_min_value(0.0f)
{
}

//----------------------------------------------------------------------
Distributed_heightfield_seedings_updates::~Distributed_heightfield_seedings_updates()
{
    // empty
}


//-----------------------------------------------------------------------------------------
void Distributed_heightfield_seedings_updates::assign_fragments_to_hosts(
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

//----------------------------------------------------------------------
void Distributed_heightfield_seedings_updates::serialize(
    mi::neuraylib::ISerializer* serializer) const
{
    serializer->write(&m_session_tag.id, 1);
    serializer->write(&m_heightfield_tag.id, 1);

    const mi::Uint32 nb_elements = mi::Uint32(m_cluster_hosts.size());
    serializer->write(&nb_elements, 1);
    for(mi::Uint32 i=0; i<nb_elements; ++i)
        serializer->write(&m_cluster_hosts[i], 1);

    const mi::Uint32 nb_patch_tag_pairs = m_patch_tag_pair_vector.size();
    serializer->write(&nb_patch_tag_pairs, 1);
    for(mi::Uint32 j=0; j<nb_patch_tag_pairs; ++j)
    {
        serializer->write(&m_patch_tag_pair_vector[j].first.min.x, 2);
        serializer->write(&m_patch_tag_pair_vector[j].first.max.x, 2);
        serializer->write(&m_patch_tag_pair_vector[j].second.id, 1);
    }
}

//----------------------------------------------------------------------
void Distributed_heightfield_seedings_updates::deserialize(
    mi::neuraylib::IDeserializer* deserializer)
{
    deserializer->read(&m_session_tag.id, 1);
    deserializer->read(&m_heightfield_tag.id, 1);

    mi::Uint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_cluster_hosts.resize(nb_elements);
    for(mi::Uint32 i=0; i<nb_elements; ++i)
        deserializer->read(&m_cluster_hosts[i], 1);

    mi::Uint32 nb_patch_tag_pairs = 0;
    deserializer->read(&nb_patch_tag_pairs, 1);
    m_patch_tag_pair_vector.resize(nb_patch_tag_pairs);
    for(mi::Uint32 j=0; j<nb_patch_tag_pairs; ++j)
    {
        deserializer->read(&m_patch_tag_pair_vector[j].first.min.x, 2);
        deserializer->read(&m_patch_tag_pair_vector[j].first.max.x, 2);
        deserializer->read(&m_patch_tag_pair_vector[j].second.id, 1);
    }
}

nv::index::IRegular_heightfield_data_locality* Distributed_heightfield_seedings_updates::get_distribution_layout(
    mi::neuraylib::IDice_transaction*  dice_transaction)
{
    // Access session
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<const nv::index::ISession>(m_session_tag));
    const mi::neuraylib::Tag& dist_layout_tag = session->get_distribution_layout();

    // Query the heightfield extent/bounds
    mi::base::Handle<const nv::index::IRegular_heightfield> heightfield(
        dice_transaction->access<const nv::index::IRegular_heightfield>(m_heightfield_tag));
    const mi::math::Bbox_struct<mi::Float32, 3> heightfield_bounds = heightfield->get_IJK_bounding_box();
    mi::math::Bbox_struct<mi::Uint32, 2> heightfield_area;
    heightfield_area.min.x = static_cast< mi::Uint32 >(heightfield_bounds.min.x);
    heightfield_area.min.y = static_cast< mi::Uint32 >(heightfield_bounds.min.y);
    heightfield_area.max.x = static_cast< mi::Uint32 >(heightfield_bounds.max.x);
    heightfield_area.max.y = static_cast< mi::Uint32 >(heightfield_bounds.max.y);

    // Access the distribution scheme
    mi::base::Handle<const nv::index::IData_distribution> distribution_layout(
        dice_transaction->access<const nv::index::IData_distribution>(dist_layout_tag));

    // Heightfield distribution layout
    nv::index::IRegular_heightfield_data_locality* data_locality = distribution_layout->retrieve_data_locality_for_editing(
        m_heightfield_tag,
        heightfield_area,
        dice_transaction);

    return data_locality;
}

//----------------------------------------------------------------------
void Distributed_heightfield_seedings_updates::set_patch_tag_lists(
    const std::vector<std::vector<std::pair<mi::math::Bbox_struct<mi::Uint32, 2>, mi::neuraylib::Tag> > >& patch_tag_lists)
{
    const mi::Uint32 nb_lists = patch_tag_lists.size();
    m_patch_tag_pair_vector.clear();

    for(mi::Uint32 i=0; i<nb_lists; ++i)
    {
        const std::vector<std::pair<mi::math::Bbox_struct<mi::Uint32, 2>, mi::neuraylib::Tag> >& input_list = patch_tag_lists[i];

        const mi::Uint32 nb_elements = input_list.size();
        for(mi::Uint32 j=0; j<nb_elements; ++j)
        {
            m_patch_tag_pair_vector.push_back(input_list[j]);
        }
    }
}

//----------------------------------------------------------------------
void Distributed_heightfield_seedings_updates::update_heightfield_surface(
    mi::neuraylib::IDice_transaction*                                                         dice_transaction,
    mi::Uint32                                                                                host_id,
    const std::vector<std::pair<mi::math::Bbox_struct<mi::Uint32, 2>, mi::neuraylib::Tag> >&  patch_tag_pair_vector,
    mi::math::Vector_struct<mi::Float32, 2>&                                                  min_max_heights)
{
    mi::base::Handle<nv::index::IRegular_heightfield_data_locality> locality(
        this->get_distribution_layout(dice_transaction));

    const mi::Uint32 nb_patches = patch_tag_pair_vector.size();
    // INFO_LOG << "NUMBER OF PATCHES OPERATED ON " << nb_patches;

    const mi::Uint32 nb_editing_patches = static_cast<mi::Uint32>(locality->get_nb_bounding_box(host_id));
    // INFO_LOG << "NUMBER OF PATCHES " << nb_editing_patches << " ON HOST: " << host_id;

    for(mi::Uint32 i=0; i<nb_patches; ++i)
    {
        const mi::math::Bbox_struct<mi::Uint32, 2>& patch = patch_tag_pair_vector[i].first;

        for(mi::Uint32 j=0; j<nb_editing_patches; ++j)
        {
            // x and y are always non-negative
            const mi::math::Bbox<mi::Sint32, 3> patch_3D_bbox = locality->get_bounding_box(host_id, j);

            if (   static_cast<mi::Uint32>(patch_3D_bbox.min.x) == patch.min.x
                && static_cast<mi::Uint32>(patch_3D_bbox.max.x) == patch.max.x
                && static_cast<mi::Uint32>(patch_3D_bbox.min.y) == patch.min.y
                && static_cast<mi::Uint32>(patch_3D_bbox.max.y) == patch.max.y)
            {
                mi::base::Handle<nv::index::IRegular_heightfield_data_edit> heightfield_data_edit(
                    locality->create_data_edit(dice_transaction, host_id, j));

                Heightfield_elevation_update_operation update_operation(
                    nv::index_common::convert_bbox_type<mi::Sint32, mi::Uint32, 2>(patch), // FIXME: remove the conversion by using Sint32 bbox
                    patch_tag_pair_vector[i].second);
                heightfield_data_edit->edit(&update_operation, dice_transaction);

                // Get min/max height values within the grid
                update_operation.get_min_max_heights(min_max_heights);
                // INFO_LOG << "RETURNED HEIGHTS: " << min_max_heights;
            }
        }
    }
}


//----------------------------------------------------------------------
void Distributed_heightfield_seedings_updates::execute_fragment(
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Size                                        index,
    mi::Size                                        count,
    const mi::neuraylib::IJob_execution_context*    context)
{
    const mi::Uint32 cluster_host = m_cluster_hosts[index];
    mi::math::Vector_struct<mi::Float32, 2> min_max_heights;

    update_heightfield_surface(dice_transaction, cluster_host, m_patch_tag_pair_vector, min_max_heights);

    set_min_max_heights(min_max_heights.x, min_max_heights.y);
}

//----------------------------------------------------------------------
void Distributed_heightfield_seedings_updates::execute_fragment_remote(
    mi::neuraylib::ISerializer*                     serializer,
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Size                                        index,
    mi::Size                                        count,
    const mi::neuraylib::IJob_execution_context*    context)
{
    const mi::Uint32 cluster_host = m_cluster_hosts[index];
    mi::math::Vector_struct<mi::Float32, 2> min_max_heights;

    update_heightfield_surface(dice_transaction, cluster_host, m_patch_tag_pair_vector, min_max_heights);

    // Serialize the maximum heights ...
    serializer->write(&min_max_heights.x, 2);
}

//----------------------------------------------------------------------
void Distributed_heightfield_seedings_updates::receive_remote_result(
    mi::neuraylib::IDeserializer*                   deserializer,
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Size                                        index,
    mi::Size                                        count)
{
    // Deserialize the height values ...
    mi::math::Vector_struct<mi::Float32, 2> min_max_heights;
    deserializer->read(&min_max_heights.x, 2);
    set_min_max_heights(min_max_heights.x, min_max_heights.y);
}

//----------------------------------------------------------------------
void Distributed_heightfield_seedings_updates::set_min_max_heights(
    mi::Float32 min_height,
    mi::Float32 max_height)
{
    mi::base::Lock::Block block(&m_compute_lock);
    if(m_max_value<max_height)
        m_max_value = max_height;
    if(m_min_value>min_height)
        m_min_value = min_height;
}

//----------------------------------------------------------------------
