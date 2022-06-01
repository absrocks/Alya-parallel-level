/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "simple_gridding_operation.h"

#include <cassert>

#include <nv/index/isession.h>
#include <nv/index/iregular_heightfield.h>
#include <nv/index/iregular_heightfield_compute_task.h>
#include <nv/index/iregular_heightfield_patch.h>
#include <nv/index/idistributed_data_access.h>

#include "common/forwarding_logger.h"
#include "common/distributed_heightfield_elevation_change.h"

#include "utilities.h"

namespace // anonymous namespace
{
//----------------------------------------------------------------------
static nv::index::IRegular_heightfield_data_locality* compute_distribution_layout(
    const mi::neuraylib::Tag&         session_tag,
    const mi::neuraylib::Tag&         heightfield_tag,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    assert(session_tag.    is_valid());
    assert(heightfield_tag.is_valid());
    assert(dice_transaction != 0);

    // Access session
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<const nv::index::ISession>(session_tag));
    assert(session.is_valid_interface());

    const mi::neuraylib::Tag& dist_layout_tag = session->get_distribution_layout();
    assert(dist_layout_tag.is_valid());

    // Query the heightfield extent/bounds
    mi::base::Handle<const nv::index::IRegular_heightfield> heightfield(
        dice_transaction->access<const nv::index::IRegular_heightfield>(heightfield_tag));
    assert(heightfield.is_valid_interface());

    const mi::math::Bbox_struct<mi::Float32, 3> heightfield_bounds = heightfield->get_IJK_bounding_box();
    mi::math::Bbox_struct<mi::Uint32, 2> heightfield_area;
    heightfield_area.min.x = static_cast<mi::Uint32>(heightfield_bounds.min.x);
    heightfield_area.min.y = static_cast<mi::Uint32>(heightfield_bounds.min.y);
    heightfield_area.max.x = static_cast<mi::Uint32>(heightfield_bounds.max.x);
    heightfield_area.max.y = static_cast<mi::Uint32>(heightfield_bounds.max.y);

    // Access the distribution scheme
    mi::base::Handle<const nv::index::IData_distribution> distribution_layout(
        dice_transaction->access<const nv::index::IData_distribution>(dist_layout_tag));
    assert(distribution_layout.is_valid_interface());

    // Heightfield distribution layout
    nv::index::IRegular_heightfield_data_locality* data_locality =
        distribution_layout->retrieve_data_locality_for_editing(
            heightfield_tag, heightfield_area, dice_transaction);
    assert(data_locality);

    return data_locality;
}
//----------------------------------------------------------------------
} // anonymous namespace


//----------------------------------------------------------------------
Simple_gridding_operation::Simple_gridding_operation(
    const mi::neuraylib::Tag&                                    session_tag,
    const mi::neuraylib::Tag&                                    heightfield_tag,
    const std::vector<mi::math::Vector_struct<mi::Float32, 2> >& polygon)
    :
    m_session_tag(session_tag),
    m_heightfield_tag(heightfield_tag),
    m_polygon_of_interest(polygon)
{
    m_mean_value = -1.0f;
}

//----------------------------------------------------------------------
Simple_gridding_operation::~Simple_gridding_operation()
{
    // empty
}

//----------------------------------------------------------------------
void Simple_gridding_operation::execute_fragment(
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Size                                        index,
    mi::Size                                        count,
    const mi::neuraylib::IJob_execution_context*    context)
{
    assert(dice_transaction != 0);

    // Create data locality ....
    mi::base::Handle<nv::index::IRegular_heightfield_data_locality> data_locality(
        compute_distribution_layout(
            m_session_tag, m_heightfield_tag, dice_transaction));
    assert(data_locality.is_valid_interface());

    // Cluster id referencing hosts that host heightfield data
    std::vector<mi::Uint32> cluster_host_ids;
    const mi::Uint32 nb_cluster_hosts = data_locality->get_nb_cluster_nodes();
    for (mi::Uint32 i = 0; i < nb_cluster_hosts; ++i)
    {
        cluster_host_ids.push_back(data_locality->get_cluster_node(i));
        INFO_LOG << "Heightfield data on host: " << data_locality->get_cluster_node(i) << " (index: " << i << ")";
    }

    // 1.Step: Query and compute mean heightfield elevation value (within the polygon area).
    Gridder_elevation_mean_compute mean_compute(
        m_session_tag,
        m_heightfield_tag,
        m_polygon_of_interest,
        cluster_host_ids);

    dice_transaction->execute_fragmented(&mean_compute, cluster_host_ids.size());
    m_mean_value = mean_compute.get_mean_value();

    // 2.Step: Set all heightfield elevation values marked as hole within the polygon area.
    const bool is_scaling_op = false;
    nv::index_common::Distributed_heightfield_elevation_change
        gridding_fill(m_session_tag,
                      m_heightfield_tag,
                      is_scaling_op,
                      1.0f, //m_mean_value,
                      m_polygon_of_interest,
                      cluster_host_ids);
    dice_transaction->execute_fragmented(&gridding_fill, cluster_host_ids.size());
}

//----------------------------------------------------------------------
Gridder_elevation_mean_compute::Gridder_elevation_mean_compute(
    const mi::neuraylib::Tag&                                       session_tag,
    const mi::neuraylib::Tag&                                       heightfield_tag,
    const std::vector<mi::math::Vector_struct<mi::Float32, 2> >&    polygon,
    const std::vector<mi::Uint32>&                                  cluster_hosts)
    :
    m_session_tag(session_tag),
    m_heightfield_tag(heightfield_tag),
    m_polygon_of_interest(polygon),
    m_cluster_hosts(cluster_hosts),
    m_result_nb_elevation_values(0)
{
    // empty
}

//----------------------------------------------------------------------
Gridder_elevation_mean_compute::~Gridder_elevation_mean_compute()
{
    // empty
}

//----------------------------------------------------------------------
void Gridder_elevation_mean_compute::assign_fragments_to_hosts(
    mi::Uint32* slots,
    mi::Size    nr_slots)
{
    if (nr_slots != m_cluster_hosts.size())
    {
        ERROR_LOG << "The number of host ids present ("
                  << m_cluster_hosts.size()
                  << ") and requested ("
                  << nr_slots
                  << ") doesn't match the requested count.";
        assert(m_cluster_hosts.size() == nr_slots);
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
void Gridder_elevation_mean_compute::serialize(
    mi::neuraylib::ISerializer* serializer) const
{
    serializer->write(&m_session_tag.id, 1);
    serializer->write(&m_heightfield_tag.id, 1);

    mi::Uint32 nb_elements = mi::Uint32(m_polygon_of_interest.size());
    serializer->write(&nb_elements, 1);
    for (mi::Uint32 i = 0; i < nb_elements; ++i)
    {
        serializer->write(&m_polygon_of_interest[i].x, 2);
    }

    nb_elements = static_cast<mi::Uint32>(m_cluster_hosts.size());
    serializer->write(&nb_elements, 1);
    for (mi::Uint32 i = 0; i < nb_elements; ++i)
    {
        serializer->write(&m_cluster_hosts[i], 1);
    }
}

//----------------------------------------------------------------------
void Gridder_elevation_mean_compute::deserialize(
    mi::neuraylib::IDeserializer* deserializer)
{
    deserializer->read(&m_session_tag.id, 1);
    deserializer->read(&m_heightfield_tag.id, 1);

    mi::Uint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_polygon_of_interest.resize(nb_elements);
    for (mi::Uint32 i = 0; i < nb_elements; ++i)
    {
        deserializer->read(&m_polygon_of_interest[i].x, 2);
    }

    nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_cluster_hosts.resize(nb_elements);
    for (mi::Uint32 i = 0; i < nb_elements; ++i)
    {
        deserializer->read(&m_cluster_hosts[i], 1);
    }
}

//----------------------------------------------------------------------
// FIXME: need tests due to value range -> bbox change
void Gridder_elevation_mean_compute::compute_local_mean(
    mi::neuraylib::IDice_transaction*              dice_transaction,
    nv::index::IRegular_heightfield_data_locality* locality,
    mi::Uint32                                     host_id,
    mi::Uint32                                     patch_id,
    mi::Float32&                                   elevation,
    mi::Uint32&                                    nb_elevations)
{
    assert(dice_transaction != 0);
    assert(locality != 0);
    assert(m_session_tag.is_valid());
    assert(m_heightfield_tag.is_valid());

    elevation = 0.0f;
    nb_elevations = 0;


    // const mi::Uint32 nb_patches = locality->get_nb_bounding_box(host_id);
    // Access session
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<const nv::index::ISession>(m_session_tag));
    assert(session.is_valid_interface());

    const mi::neuraylib::Tag& access_factory_tag = session->get_data_access_factory();
    assert(access_factory_tag.is_valid());

    mi::base::Handle<const nv::index::IDistributed_data_access_factory> access_factory(
        dice_transaction->access<const nv::index::IDistributed_data_access_factory>(access_factory_tag));
    assert(access_factory.is_valid_interface());

    const mi::math::Bbox<mi::Sint32, 3> patch_3D_bbox = locality->get_bounding_box(host_id, patch_id);
    mi::math::Bbox<mi::Uint32, 2> roi; // x and y are always non-negative
    roi.min.x = static_cast<mi::Uint32>(patch_3D_bbox.min.x);
    roi.min.y = static_cast<mi::Uint32>(patch_3D_bbox.min.y);
    roi.max.x = static_cast<mi::Uint32>(patch_3D_bbox.max.x);
    roi.max.y = static_cast<mi::Uint32>(patch_3D_bbox.max.y);

    mi::base::Handle<nv::index::IRegular_heightfield_data_access> heightfield_data_access(
        access_factory->create_regular_heightfield_data_access(m_heightfield_tag));
    assert(heightfield_data_access.is_valid_interface());

    heightfield_data_access->access(roi, dice_transaction);
    const mi::math::Bbox_struct<mi::Uint32, 2>& patch_bbox =
        heightfield_data_access->get_patch_bounding_box();
    mi::Float32* elevation_values = heightfield_data_access->get_elevation_values();
    assert(elevation_values != 0);

    mi::Uint32 nb_elevation_values = (patch_bbox.max.x - patch_bbox.min.x) * (patch_bbox.max.y - patch_bbox.min.y);
    if(m_polygon_of_interest.empty())
    {
        for (mi::Uint32 i = 0; i < nb_elevation_values; ++i)
        {
            if (!nv::index::IRegular_heightfield_patch::is_hole(*elevation_values))
            {
                elevation += (*elevation_values);
                nb_elevations++;
            }
            elevation_values++;
        }
    }
    else
    {
        mi::math::Vector<mi::Float32, 2> min_value =
            mi::math::Vector<mi::Float32, 2>(static_cast<mi::Float32>(patch_bbox.max.x),
                                             static_cast<mi::Float32>(patch_bbox.max.y));
        mi::math::Vector<mi::Float32, 2> max_value =
            mi::math::Vector<mi::Float32, 2>(static_cast<mi::Float32>(patch_bbox.min.x),
                                             static_cast<mi::Float32>(patch_bbox.min.y));

        const mi::Uint32 nb_polygon_verts = m_polygon_of_interest.size();

        for (mi::Uint32 i = 0; i < nb_polygon_verts; ++i)
        {
            const mi::math::Vector<mi::Float32, 2>& polygon_vert = m_polygon_of_interest[i];
            if (polygon_vert.x < min_value.x)
            {
                min_value.x = polygon_vert.x;
            }
            if (polygon_vert.y < min_value.y)
            {
                min_value.y = polygon_vert.y;
            }
            if (polygon_vert.x > max_value.x)
            {
                max_value.x = polygon_vert.x;
            }
            if (polygon_vert.y > max_value.y)
            {
                max_value.y = polygon_vert.y;
            }
        }

        // Clamp to bounding box
        if (min_value.x < patch_bbox.min.x)
        {
            min_value.x = static_cast<mi::Float32>(patch_bbox.min.x);
        }
        if (min_value.y < patch_bbox.min.y)
        {
            min_value.y = static_cast<mi::Float32>(patch_bbox.min.y);
        }
        if (max_value.x > patch_bbox.max.x)
        {
            max_value.x = static_cast<mi::Float32>(patch_bbox.max.x);
        }
        if (max_value.y > patch_bbox.max.y)
        {
            max_value.y = static_cast<mi::Float32>(patch_bbox.max.y);
        }

        const mi::Uint32 min_value_y =
            static_cast< mi::Uint32 >(mi::math::ceil(min_value.y) - patch_bbox.min.y);
        const mi::Uint32 max_value_y =
            static_cast< mi::Uint32 >(mi::math::floor(max_value.y) - patch_bbox.min.y);
        std::vector<mi::Float32> intersecton_spans;

        for (mi::Uint32 y = min_value_y; y <= max_value_y; ++y)
        {
            compute_intersection(static_cast<mi::Float32>(y+patch_bbox.min.y), m_polygon_of_interest, intersecton_spans);
            mi::Float32* vertices = &elevation_values[(patch_bbox.max.x-patch_bbox.min.x+1) * y];

            const mi::Uint32 nb_intersections = intersecton_spans.size();
            if ((nb_intersections % 2) == 0)
            {
                for(mi::Uint32 i = 0; i < nb_intersections; i += 2)
                {
                    mi::Sint32 from = static_cast< mi::Sint32 >(intersecton_spans[i]);
                    mi::Sint32 to   = static_cast< mi::Sint32 >(intersecton_spans[i+1]);

                    mi::Float32* elevation_values = vertices;

                    for (mi::Sint32 k = patch_bbox.min.x;
                         k <= static_cast< mi::Sint32 >(patch_bbox.max.x); ++k)
                    {
                        if ((k >= from) && (k <= to))
                        {
                            if (!nv::index::IRegular_heightfield_patch::is_hole(*elevation_values))
                            {
                                elevation += (*elevation_values);
                                nb_elevations++;
                            }
                        }
                        elevation_values++;
                    }
                }
            }
            intersecton_spans.clear();
        }
    }
}

//----------------------------------------------------------------------
void Gridder_elevation_mean_compute::execute_fragment(
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Size                                        index,
    mi::Size                                        count,
    const mi::neuraylib::IJob_execution_context*    context)
{
    const mi::Uint32 cluster_host = m_cluster_hosts[index];
    mi::base::Handle<nv::index::IRegular_heightfield_data_locality> locality(
        compute_distribution_layout(m_session_tag, m_heightfield_tag, dice_transaction));
    assert(locality.is_valid_interface());

    const mi::Uint32 nb_patches = static_cast<mi::Uint32>(locality->get_nb_bounding_box(cluster_host));

    for (mi::Uint32 i = 0; i < nb_patches; ++i)
    {
        mi::Uint32 nb_elevation_values = 0;
        mi::Float32 accumulated_elevation = 0.0f;
        compute_local_mean(dice_transaction, locality.get(), cluster_host, i, accumulated_elevation, nb_elevation_values);

        this->compute_mean(accumulated_elevation, nb_elevation_values);
    }
}

//----------------------------------------------------------------------
void Gridder_elevation_mean_compute::execute_fragment_remote(
    mi::neuraylib::ISerializer*                     serializer,
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Size                                        index,
    mi::Size                                        count,
    const mi::neuraylib::IJob_execution_context*    context)
{
    const mi::Uint32 cluster_host = m_cluster_hosts[index];
    mi::base::Handle<nv::index::IRegular_heightfield_data_locality> locality(
        compute_distribution_layout(
            m_session_tag, m_heightfield_tag, dice_transaction));
    assert(locality.is_valid_interface());

    const mi::Uint32 nb_patches = static_cast<mi::Uint32>(locality->get_nb_bounding_box(cluster_host));

    serializer->write(&nb_patches, 1);
    for(mi::Uint32 i = 0; i < nb_patches; ++i)
    {
        mi::Uint32 nb_elevation_values = 0;
        mi::Float32 accumulated_elevation = 0.0f;
        compute_local_mean(dice_transaction, locality.get(), cluster_host, i, accumulated_elevation, nb_elevation_values);

        // Writing the results to the stream
        serializer->write(&nb_elevation_values, 1);
        serializer->write(&accumulated_elevation, 1);
    }
}

//----------------------------------------------------------------------
void Gridder_elevation_mean_compute::receive_remote_result(
    mi::neuraylib::IDeserializer*                   deserializer,
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Size                                        index,
    mi::Size                                        count)
{
    // Reading the results from the stream
    mi::Uint32 nb_patches = 0;
    deserializer->read(&nb_patches, 1);
    for(mi::Uint32 i = 0; i < nb_patches; ++i)
    {
        mi::Uint32 nb_elevations = 0;
        deserializer->read(&nb_elevations, 1);
        mi::Float32 elevation = 0.f;
        deserializer->read(&elevation, 1);

        this->compute_mean(elevation, nb_elevations);
    }
}

//----------------------------------------------------------------------
void Gridder_elevation_mean_compute::compute_mean(
    mi::Float32 elevation,
    mi::Uint32  nb_elevations)
{
    // Lock mean compute setup.
    //
    // Could also be removed if array size is set in advance and atomic
    // counter is used for the number of accumulated heightfield elevation values.
    //
    mi::base::Lock::Block block(&m_compute_lock);

    // The number of accumulated heightfield elevation values.
    m_result_nb_elevation_values += nb_elevations;
    // Accumulated heightfield values per patch.
    m_result_accumulated_elevation_values.push_back(elevation);
}

//----------------------------------------------------------------------
mi::Float32 Gridder_elevation_mean_compute::get_mean_value() const
{
    // Actually compute the mean value on the fly.
    // Accumulate the weighted per patch elevation values.
    mi::Float32 elevation_mean = 0.f;
    const mi::Uint32 nb_elevations = m_result_accumulated_elevation_values.size();
    for (mi::Uint32 i = 0; i < nb_elevations; ++i)
    {
        assert(m_result_nb_elevation_values != 0);
        elevation_mean += (m_result_accumulated_elevation_values[i] / mi::Float32(m_result_nb_elevation_values));
    }
    return elevation_mean;
}

//----------------------------------------------------------------------
