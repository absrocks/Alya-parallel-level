/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "distributed_heightfield_elevation_change.h"

#include <cassert>

#include <nv/index/isession.h>
#include <nv/index/iregular_heightfield.h>
#include <nv/index/iregular_heightfield_compute_task.h>
#include <nv/index/iregular_heightfield_patch.h>

#include "forwarding_logger.h"

namespace nv {
namespace index_common {

namespace {

bool intersection(
    mi::Float32                                     Y,
    const mi::math::Vector_struct<mi::Float32, 2>&  from,
    const mi::math::Vector_struct<mi::Float32, 2>&  to,
    mi::Float32&                                    x_value)
{
    if(from.y<Y && to.y<Y)
        return false;

    if(from.y>Y && to.y>Y)
        return false;

    x_value = to.x + (((Y - to.y)/(to.y-from.y) ) * (to.x-from.x));

    return true;
}

void compute_intersection(
    mi::Float32                                                     Y,
    const std::vector<mi::math::Vector_struct<mi::Float32, 2> >&    polygon,
    std::vector<mi::Float32>&                                       intersection_spans)
{
    const mi::Size nb_vertices = polygon.size();
    mi::Float32 x_value = 0.0f;
    assert(nb_vertices > 0);
    for(mi::Size i=0; i<nb_vertices-1; ++i)
    {
        const bool intersected = intersection(Y, polygon[i], polygon[i+1], x_value);
        if(intersected)
        {
            intersection_spans.push_back(x_value);
        }
    }

    const bool intersected = intersection(Y, polygon[nb_vertices-1], polygon[0], x_value);
    if(intersected)
    {
        intersection_spans.push_back(x_value);
    }

    std::sort(
        intersection_spans.begin(),
        intersection_spans.end());
}

} // namespace

/// Heightfield elevation value change compute task.
class Heightfield_elevation_scale_operation :
        public mi::base::Interface_implement<nv::index::IRegular_heightfield_compute_task>
{
public:
    /// constructor
    Heightfield_elevation_scale_operation(
        bool                                                            scale_operation,
        mi::Float32                                                     elevation_manipulation_value,
        const std::vector<mi::math::Vector_struct<mi::Float32, 2> >&    polygon)
        : m_scale_operation(scale_operation),
          m_elevation_manipulation_value(elevation_manipulation_value),
          m_polygon_of_interest(polygon)
    {
        // empty
    }

    /// 
    virtual void get_region_of_interest_for_compute(
        mi::math::Bbox_struct<mi::Sint32, 2>& roi) const
    {
        const mi::Size nb_polygon_verts = m_polygon_of_interest.size();
        if (nb_polygon_verts > 0)
        {
            const mi::math::Vector<mi::Float32, 2>& polygon_vert = m_polygon_of_interest[0];
            roi.min.x = static_cast< mi::Sint32 >(mi::math::floor(polygon_vert.x));
            roi.max.x = static_cast< mi::Sint32 >(mi::math::ceil (polygon_vert.x));
            roi.min.y = static_cast< mi::Sint32 >(mi::math::floor(polygon_vert.y));
            roi.max.y = static_cast< mi::Sint32 >(mi::math::ceil (polygon_vert.y));

            for (mi::Uint32 i = 1; i < nb_polygon_verts; ++i)
            {
                const mi::math::Vector<mi::Float32, 2>& polygon_vert = m_polygon_of_interest[i];
                if (polygon_vert.x < roi.min.x)
                {
                    roi.min.x = static_cast<mi::Sint32>(mi::math::floor(polygon_vert.x));
                }
                if (polygon_vert.x > roi.max.x)
                {
                    roi.max.x = static_cast<mi::Sint32>(mi::math::floor(polygon_vert.x));
                }
                if (polygon_vert.y < roi.min.y)
                {
                    roi.min.y = static_cast<mi::Sint32>(mi::math::floor(polygon_vert.y));
                }
                if (polygon_vert.y > roi.max.y)
                {
                    roi.max.y = static_cast<mi::Sint32>(mi::math::floor(polygon_vert.y));
                }
            }
        }
    }

    /// compute: with normal version    
    virtual bool compute(
        const mi::math::Bbox_struct<mi::Sint32, 2>&     ijk_patch_bbox,
        mi::Float32*                                    elevation_values,
        mi::math::Vector_struct<mi::Float32, 3>*        normal_values,
        mi::neuraylib::IDice_transaction*               dice_transaction) const
    {
        mi::math::Vector<mi::Float32, 2> min_value =
            mi::math::Vector<mi::Float32, 2>(static_cast<mi::Float32>(ijk_patch_bbox.max.x),
                                             static_cast<mi::Float32>(ijk_patch_bbox.max.y));
        mi::math::Vector<mi::Float32, 2> max_value =
            mi::math::Vector<mi::Float32, 2>(static_cast<mi::Float32>(ijk_patch_bbox.min.x),
                                             static_cast<mi::Float32>(ijk_patch_bbox.min.y));

        const mi::Sint32 nb_elevation_values = 
            (ijk_patch_bbox.max.x - ijk_patch_bbox.min.x) * (ijk_patch_bbox.max.y - ijk_patch_bbox.min.y);
        assert(nb_elevation_values >= 0);

        // Predefined normal
        mi::math::Vector_struct<mi::Float32, 3> up = { 0.0f, 0.0f, 1.0f, };

        if (m_polygon_of_interest.empty())
        {
            if (m_scale_operation)
            {
                for (mi::Sint32 i = 0; i < nb_elevation_values; ++i)
                {
                    if (!nv::index::IRegular_heightfield_patch::is_hole(elevation_values[i]))
                    {
                        elevation_values[i] *= m_elevation_manipulation_value;
                        normal_values[i] = up;
                    }
                }
            }
            else
            {
                for (mi::Sint32 i=0; i < nb_elevation_values; ++i)
                {
                    if (!nv::index::IRegular_heightfield_patch::is_hole(elevation_values[i]))
                    {
                        elevation_values[i] = m_elevation_manipulation_value;
                        normal_values[i] = up;
                    }
                }
            }
        }
        else
        {
            const mi::Size nb_polygon_verts = m_polygon_of_interest.size();
            
            for (mi::Size i = 0; i < nb_polygon_verts; ++i)
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
            if (min_value.x < ijk_patch_bbox.min.x)
            {
                min_value.x = static_cast<mi::Float32>(ijk_patch_bbox.min.x);
            }
            if (min_value.y < ijk_patch_bbox.min.y)
            {
                min_value.y = static_cast<mi::Float32>(ijk_patch_bbox.min.y);
            }
            if (max_value.x > ijk_patch_bbox.max.x)
            {
                max_value.x = static_cast<mi::Float32>(ijk_patch_bbox.max.x);
            }
            if (max_value.y > ijk_patch_bbox.max.y)
            {
                max_value.y = static_cast<mi::Float32>(ijk_patch_bbox.max.y);
            }

            const mi::Sint32 min_value_y =
                static_cast< mi::Sint32 >(mi::math::ceil(min_value.y) - 
                                          static_cast<mi::Float32>(ijk_patch_bbox.min.y));
            const mi::Sint32 max_value_y =
                static_cast< mi::Sint32 >(mi::math::floor(max_value.y) - 
                                          static_cast<mi::Float32>(ijk_patch_bbox.min.y));
            std::vector<mi::Float32> intersecton_spans;
            
            for(mi::Sint32 y = min_value_y; y < max_value_y; ++y)
            {
                compute_intersection(static_cast<mi::Float32>(y+ijk_patch_bbox.min.y), m_polygon_of_interest, intersecton_spans);
                mi::Float32* vertices = &elevation_values[(ijk_patch_bbox.max.x-ijk_patch_bbox.min.x) * y];
                mi::math::Vector_struct<mi::Float32, 3>* normals = &normal_values[(ijk_patch_bbox.max.x-ijk_patch_bbox.min.x) * y];

                const mi::Size nb_intersections = intersecton_spans.size();
                for(mi::Size i = 0; (i + 1) < nb_intersections; i += 2)
                {
                    mi::Sint32 from = static_cast<mi::Sint32>(intersecton_spans[i]);
                    mi::Sint32 to   = static_cast<mi::Sint32>(intersecton_spans[i+1]);

                    if(!m_scale_operation)
                    {
                        mi::Float32* VERT = vertices;
                        mi::math::Vector_struct<mi::Float32, 3>* NORM = normals;
                        
                        for (mi::Sint32 span_x = ijk_patch_bbox.min.x; span_x < static_cast<mi::Sint32>(ijk_patch_bbox.max.x); ++span_x)
                        {
                            if((span_x >= from) && (span_x <= to))
                            {
                                *VERT = m_elevation_manipulation_value;
                                *NORM = up;
                            }
                            VERT++;
                            NORM++;
                        }
                    }
                    else
                    {
                        mi::Float32* VERT = vertices;
                        mi::math::Vector_struct<mi::Float32, 3>* NORM = normals;

                        for (mi::Sint32 span_x = ijk_patch_bbox.min.x; span_x < static_cast<mi::Sint32>(ijk_patch_bbox.max.x); ++span_x)
                        {
                            if((span_x >= from) && (span_x <= to))
                            {
                                *VERT *= m_elevation_manipulation_value;
                                *NORM = up;
                            }
                            VERT++;
                            NORM++;
                        }
                    }
                }
                intersecton_spans.clear();
            }
        }

        return true;
    }

    /// compute: without normal version
    virtual bool compute(
        const mi::math::Bbox_struct<mi::Sint32, 2>&     ijk_patch_bbox,
        mi::Float32*                                    elevation_values,
        mi::neuraylib::IDice_transaction*               dice_transaction) const
    {
        mi::math::Vector<mi::Float32, 2> min_value = 
            mi::math::Vector<mi::Float32, 2>(static_cast<mi::Float32>(ijk_patch_bbox.max.x),
                                             static_cast<mi::Float32>(ijk_patch_bbox.max.y));
        mi::math::Vector<mi::Float32, 2> max_value =
            mi::math::Vector<mi::Float32, 2>(static_cast<mi::Float32>(ijk_patch_bbox.min.x),
                                             static_cast<mi::Float32>(ijk_patch_bbox.min.y));

        const mi::Sint32 nb_elevation_values =
            (ijk_patch_bbox.max.x - ijk_patch_bbox.min.x) * (ijk_patch_bbox.max.y - ijk_patch_bbox.min.y);
        assert(nb_elevation_values >= 0);

        if (m_polygon_of_interest.empty())
        {
            if (m_scale_operation)
            {
                for(mi::Sint32 i = 0; i < nb_elevation_values; ++i)
                {
                    if (!nv::index::IRegular_heightfield_patch::is_hole(elevation_values[i]))
                    {
                        elevation_values[i] *= m_elevation_manipulation_value;
                    }
                }
            }
            else
            {
                for (mi::Sint32 i = 0; i < nb_elevation_values; ++i)
                {
                    if (!nv::index::IRegular_heightfield_patch::is_hole(elevation_values[i]))
                    {
                        elevation_values[i] = m_elevation_manipulation_value;
                    }
                }
            }
        }
        else
        {
            const mi::Size nb_polygon_verts = m_polygon_of_interest.size();

            for(mi::Size i = 0; i < nb_polygon_verts; ++i)
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
            if (min_value.x < ijk_patch_bbox.min.x)
            {
                min_value.x = static_cast<mi::Float32>(ijk_patch_bbox.min.x);
            }
            if (min_value.y < ijk_patch_bbox.min.y)
            {
                min_value.y = static_cast<mi::Float32>(ijk_patch_bbox.min.y);
            }
            if (max_value.x > ijk_patch_bbox.max.x)
            {
                max_value.x = static_cast<mi::Float32>(ijk_patch_bbox.max.x);
            }
            if (max_value.y > ijk_patch_bbox.max.y)
            {
                max_value.y = static_cast<mi::Float32>(ijk_patch_bbox.max.y);
            }

            const mi::Sint32 min_value_y =
                static_cast<mi::Sint32>(mi::math::ceil(min_value.y) - 
                                        static_cast<mi::Float32>(ijk_patch_bbox.min.y));
            const mi::Sint32 max_value_y =
                static_cast<mi::Sint32>(mi::math::floor(max_value.y) - 
                                        static_cast<mi::Float32>(ijk_patch_bbox.min.y));
            std::vector<mi::Float32> intersecton_spans;
            
            for(mi::Sint32 y = min_value_y; y < max_value_y; ++y)
            {
                compute_intersection(static_cast<mi::Float32>(y + ijk_patch_bbox.min.y), m_polygon_of_interest, intersecton_spans);
                mi::Float32* vertices = &elevation_values[(ijk_patch_bbox.max.x - ijk_patch_bbox.min.x) * y];

                const mi::Size nb_intersections = intersecton_spans.size();
                for(mi::Size i = 0; (i + 1) < nb_intersections; i += 2)
                {
                    mi::Sint32 from = static_cast<mi::Sint32>(intersecton_spans[i]);
                    mi::Sint32 to   = static_cast<mi::Sint32>(intersecton_spans[i+1]);

                    if(!m_scale_operation)
                    {
                        mi::Float32* VERT = vertices;
                        for (mi::Sint32 span_x = ijk_patch_bbox.min.x;
                             span_x < static_cast<mi::Sint32>(ijk_patch_bbox.max.x); ++span_x)
                        {
                            if ((span_x >= from) && (span_x <= to))
                            {
                                *VERT = m_elevation_manipulation_value;
                            }
                            VERT++;
                        }
                    }
                    else
                    {
                        mi::Float32* VERT = vertices;
                        for(mi::Sint32 span_x = ijk_patch_bbox.min.x;
                            span_x < static_cast<mi::Sint32>(ijk_patch_bbox.max.x); ++span_x)
                        {
                            if ((span_x >= from) && (span_x <= to))
                            {
                                *VERT *= m_elevation_manipulation_value;
                            }
                            VERT++;
                        }
                    }
                }
                intersecton_spans.clear();
            }
        }

        return true;
    }

    virtual bool user_defined_normal_computation() const { return true; }

private:
    bool                                                    m_scale_operation;
    mi::Float32                                             m_elevation_manipulation_value;
    std::vector<mi::math::Vector_struct<mi::Float32, 2> >   m_polygon_of_interest;
};

//-----------------------------------------------------------------------------------------
Distributed_heightfield_elevation_change::Distributed_heightfield_elevation_change(
    const mi::neuraylib::Tag&                                       session_tag,
    const mi::neuraylib::Tag&                                       heightfield_tag,
    bool                                                            scale_operation,
    mi::Float32                                                     elevation_manipulation_value,
    const std::vector<mi::math::Vector_struct<mi::Float32, 2> >&    polygon,
    const std::vector<mi::Uint32>&                                  cluster_hosts)
    :  m_session_tag(session_tag),
       m_heightfield_tag(heightfield_tag),
       m_scale_operation(scale_operation),
       m_elevation_manipulation_value(elevation_manipulation_value),
       m_polygon_of_interest(polygon),
       m_cluster_hosts(cluster_hosts)
{
    // empty
}

//-----------------------------------------------------------------------------------------
Distributed_heightfield_elevation_change::~Distributed_heightfield_elevation_change()
{
    // empty
}

//-----------------------------------------------------------------------------------------
void Distributed_heightfield_elevation_change::assign_fragments_to_hosts(
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
        assert(m_cluster_hosts.size()==nr_slots);
    }
    else
    {
        for (mi::Size i = 0; i < nr_slots; ++i)
        {
            slots[i] = m_cluster_hosts[i];
        }
    }
}

//----------------------------------------------------------------------
void Distributed_heightfield_elevation_change::serialize(
    mi::neuraylib::ISerializer* serializer) const
{
    serializer->write(&m_session_tag.id, 1);
    serializer->write(&m_heightfield_tag.id, 1);
    serializer->write(&m_scale_operation, 1);
    serializer->write(&m_elevation_manipulation_value, 1);

    mi::Uint32 nb_elements = mi::Uint32(m_polygon_of_interest.size());
    serializer->write(&nb_elements, 1);
    for (mi::Uint32 i = 0; i < nb_elements; ++i)
    {
        serializer->write(&m_polygon_of_interest[i].x, 2);
    }

    nb_elements = mi::Uint32(m_cluster_hosts.size());
    serializer->write(&nb_elements, 1);
    for (mi::Uint32 i = 0; i < nb_elements; ++i)
    {
        serializer->write(&m_cluster_hosts[i], 1);
    }
}

//----------------------------------------------------------------------
void Distributed_heightfield_elevation_change::deserialize(
    mi::neuraylib::IDeserializer* deserializer)
{
    deserializer->read(&m_session_tag.id, 1);
    deserializer->read(&m_heightfield_tag.id, 1);
    deserializer->read(&m_scale_operation, 1);
    deserializer->read(&m_elevation_manipulation_value, 1);

    mi::Uint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_polygon_of_interest.resize(nb_elements);
    for (mi::Uint32 i=0; i<nb_elements; ++i)
    {
        deserializer->read(&m_polygon_of_interest[i].x, 2);
    }

    nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_cluster_hosts.resize(nb_elements);
    for (mi::Uint32 i=0; i<nb_elements; ++i)
    {
        deserializer->read(&m_cluster_hosts[i], 1);
    }
}

//----------------------------------------------------------------------
mi::math::Bbox<mi::Float32, 3> Distributed_heightfield_elevation_change::edit_heightfield_data(
    mi::neuraylib::IDice_transaction*               dice_transaction,
    nv::index::IRegular_heightfield_data_locality*  locality,
    mi::Uint32                                      host_id)
{
    // Overall bounding box of all changed patches
    mi::math::Bbox<mi::Float32, 3> bbox;

    const mi::Uint32 nb_bbox = static_cast<mi::Uint32>(locality->get_nb_bounding_box(host_id));
    for (mi::Uint32 i = 0; i < nb_bbox; ++i)
    {
        mi::base::Handle<nv::index::IRegular_heightfield_data_edit> heightfield_data_edit(
            locality->create_data_edit(dice_transaction, host_id, i));
        assert(heightfield_data_edit.is_valid_interface());
        
        Heightfield_elevation_scale_operation scale_operation(
            m_scale_operation, m_elevation_manipulation_value, m_polygon_of_interest);

        // Execute the scale operation
        heightfield_data_edit->edit(&scale_operation, dice_transaction);

        // Retrieve the updated bounding box and add it to the overall bbox
        mi::math::Bbox_struct<mi::Float32, 3> updated_bbox;
        heightfield_data_edit->get_updated_bounding_box(updated_bbox);
        bbox.push_back(updated_bbox);
    }

    return bbox;
}

// --------------------------------------------------------------------------------------------

namespace {

nv::index::IRegular_heightfield_data_locality * new_heightfield_distribution_layout(
    const mi::neuraylib::Tag&                   session_tag,
    const mi::neuraylib::Tag&                   heightfield_tag,
    const mi::math::Bbox_struct<mi::Uint32, 2>& ij_query_bounds,
    bool                                        is_edit,
    mi::neuraylib::IDice_transaction*           dice_transaction)
{
    assert(session_tag  .is_valid());
    assert(heightfield_tag.is_valid());
    assert(dice_transaction != 0);

    // Access session
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<const nv::index::ISession>(session_tag));
    assert(session.is_valid_interface());
    const mi::neuraylib::Tag& dist_layout_tag = session->get_distribution_layout();
    assert(dist_layout_tag.is_valid());

    // Access the distribution scheme
    mi::base::Handle<const nv::index::IData_distribution> distribution_layout(
        dice_transaction->access<const nv::index::IData_distribution>(dist_layout_tag));
    assert(distribution_layout.is_valid_interface());

    // Heightfield distribution layout: this is new-ed. Use handle at the
    // caller side. (We can not use handle here since it is
    // immidiately out of scope.)
    nv::index::IRegular_heightfield_data_locality* p_data_locality = 0;
    if (is_edit)
    {
        p_data_locality = distribution_layout->retrieve_data_locality_for_editing(
            heightfield_tag, ij_query_bounds, dice_transaction);
    }
    else
    {
        p_data_locality = distribution_layout->retrieve_data_locality(
            heightfield_tag, ij_query_bounds, dice_transaction);
    }

    return p_data_locality;
}

nv::index::IRegular_heightfield_data_locality* new_whole_heightfield_distribution_layout(
    const mi::neuraylib::Tag&             session_tag,
    const mi::neuraylib::Tag&             heightfield_tag,
    bool                                  is_edit,
    mi::neuraylib::IDice_transaction*     dice_transaction,
    mi::math::Bbox_struct<mi::Uint32, 3>& heightfield_bbox)
{
    assert(session_tag  .is_valid());
    assert(heightfield_tag.is_valid());
    assert(dice_transaction != 0);

    // Query the heightfield bbox
    mi::base::Handle<const nv::index::IRegular_heightfield> heightfield(
        dice_transaction->access<const nv::index::IRegular_heightfield>(heightfield_tag));
    assert(heightfield.is_valid_interface());

    // return the ijk heightfield bbox
    mi::math::Bbox<mi::Float32, 3> const & ijk_float_bbox = heightfield->get_IJK_bounding_box();

    // convert to heightfield bbox
    for (mi::Sint32 i = 0; i < 3; ++i)
    {
        // check only
        assert(floor(ijk_float_bbox.min[i]) >= 0.0);
        assert(ceil( ijk_float_bbox.max[i]) >= 0.0);
    }
    heightfield_bbox.min.x = static_cast<mi::Uint32>(floor(ijk_float_bbox.min.x));
    heightfield_bbox.max.x = static_cast<mi::Uint32>(ceil( ijk_float_bbox.max.x));
    heightfield_bbox.min.y = static_cast<mi::Uint32>(floor(ijk_float_bbox.min.y));
    heightfield_bbox.max.y = static_cast<mi::Uint32>(ceil( ijk_float_bbox.max.y));
    heightfield_bbox.min.z = static_cast<mi::Uint32>(floor(ijk_float_bbox.min.z));
    heightfield_bbox.max.z = static_cast<mi::Uint32>(ceil( ijk_float_bbox.max.z));

    // convert to heightfield area (IJ) for query
    mi::math::Bbox_struct<mi::Uint32, 2> heightfield_ij_bbox;
    heightfield_ij_bbox.min.x = heightfield_bbox.min.x;
    heightfield_ij_bbox.max.x = heightfield_bbox.max.x;
    heightfield_ij_bbox.min.y = heightfield_bbox.min.y;
    heightfield_ij_bbox.max.y = heightfield_bbox.max.y;

    nv::index::IRegular_heightfield_data_locality * p_locality = new_heightfield_distribution_layout(
        session_tag,
        heightfield_tag,
        heightfield_ij_bbox,
        is_edit,
        dice_transaction);

    return p_locality;
}

}

//----------------------------------------------------------------------
void Distributed_heightfield_elevation_change::execute_fragment(
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Size                                        index,
    mi::Size                                        count,
    const mi::neuraylib::IJob_execution_context*    context)
{
    const mi::Uint32 cluster_host = m_cluster_hosts[index];

    mi::math::Bbox_struct<mi::Uint32, 3> heightfield_bbox_result;
    bool const is_edit = true;
    mi::base::Handle<nv::index::IRegular_heightfield_data_locality> heightfield_data_locality(
        new_whole_heightfield_distribution_layout(m_session_tag,
                                                  m_heightfield_tag,
                                                  is_edit,
                                                  dice_transaction,
                                                  heightfield_bbox_result));
    assert(heightfield_data_locality.is_valid_interface());

    mi::math::Bbox<mi::Float32, 3> bbox = edit_heightfield_data(dice_transaction,
                                                                heightfield_data_locality.get(), cluster_host);

    mi::base::Lock::Block block(&m_data_lock);
    m_fragment_bounding_boxes.push_back(bbox);
}

//----------------------------------------------------------------------
void Distributed_heightfield_elevation_change::execute_fragment_remote(
    mi::neuraylib::ISerializer*                     serializer,
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Size                                        index,
    mi::Size                                        count,
    const mi::neuraylib::IJob_execution_context*    context)
{
    const mi::Uint32 cluster_host = m_cluster_hosts[index];

    mi::math::Bbox_struct<mi::Uint32, 3> heightfield_bbox_result;
    bool const is_edit = true;
    mi::base::Handle<nv::index::IRegular_heightfield_data_locality> heightfield_data_locality(
        new_whole_heightfield_distribution_layout(m_session_tag,
                                                  m_heightfield_tag,
                                                  is_edit,
                                                  dice_transaction,
                                                  heightfield_bbox_result));
    assert(heightfield_data_locality.is_valid_interface());

    mi::math::Bbox<mi::Float32, 3> bbox = edit_heightfield_data(dice_transaction,
                                                                heightfield_data_locality.get(), cluster_host);

    // Serialize bounding box result
    serializer->write(&bbox.min.x, 6);
}

//----------------------------------------------------------------------
void Distributed_heightfield_elevation_change::receive_remote_result(
    mi::neuraylib::IDeserializer*                   deserializer,
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Size                                        index,
    mi::Size                                        count)
{
    // Receive bounding box result from remote host
    mi::math::Bbox<mi::Float32, 3> bbox;
    deserializer->read(&bbox.min.x, 6);

    mi::base::Lock::Block block(&m_data_lock);
    m_fragment_bounding_boxes.push_back(bbox);
}

void Distributed_heightfield_elevation_change::get_updated_bounding_box(
    mi::math::Bbox_struct<mi::Float32, 3>& bbox_struct) const
{
    mi::math::Bbox<mi::Float32, 3> bbox(bbox_struct);

    // This will enlarge the bounding box so that the resulting heightfield fits in
    for (mi::Size i = 0; i < m_fragment_bounding_boxes.size(); ++i)
    {
        bbox.insert(m_fragment_bounding_boxes[i]);
    }
    bbox_struct = bbox;
}

//----------------------------------------------------------------------
}} // namespace nv::index_common
