/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "distributed_heightfield_elevation_delete.h"

#include <cassert>

#include <nv/index/isession.h>
#include <nv/index/iregular_heightfield.h>
#include <nv/index/iregular_heightfield_compute_task.h>
#include <nv/index/iregular_heightfield_patch.h>

#include "common/forwarding_logger.h"
#include "common/type_conversion_utility.h"

#include "scene_utility.h"
#include "utilities.h"

class Heightfield_polygon_based_delete_operation :
    public mi::base::Interface_implement<nv::index::IRegular_heightfield_compute_task>
{
public:
    Heightfield_polygon_based_delete_operation(
        const std::vector<mi::math::Vector_struct<mi::Float32, 2> >& polygon)
      : m_polygon_of_interest(polygon)
    {
        // empty
    }

public:
    // Implemented IRegular_heightfield_compute_task

    virtual void get_region_of_interest_for_compute(
        mi::math::Bbox_struct<mi::Sint32, 2>& roi) const
    {
        const mi::Size nb_polygon_verts = m_polygon_of_interest.size();
        if (nb_polygon_verts > 0)
        {
            const mi::math::Vector<mi::Float32, 2>& polygon_vert = m_polygon_of_interest[0];
            roi.min.x = static_cast<mi::Sint32>(mi::math::floor(polygon_vert.x));
            roi.max.x = static_cast<mi::Sint32>(mi::math::ceil(polygon_vert.x));
            roi.min.y = static_cast<mi::Sint32>(mi::math::floor(polygon_vert.y));
            roi.max.y = static_cast<mi::Sint32>(mi::math::ceil(polygon_vert.y));

            for (mi::Size i = 1; i < nb_polygon_verts; ++i)
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

    // FIXME: Updated to bounding box access, but not tested
    virtual bool compute(
        const mi::math::Bbox_struct<mi::Sint32, 2>& ij_patch_bbox,
        mi::Float32*                                elevation_values,
        mi::neuraylib::IDice_transaction*           dice_transaction) const
    {
        mi::math::Vector<mi::Float32, 2> min_value =
            mi::math::Vector<mi::Float32, 2>(static_cast<mi::Float32>(ij_patch_bbox.max.x),
                                             static_cast<mi::Float32>(ij_patch_bbox.max.y));
        mi::math::Vector<mi::Float32, 2> max_value =
            mi::math::Vector<mi::Float32, 2>(static_cast<mi::Float32>(ij_patch_bbox.min.x),
                                             static_cast<mi::Float32>(ij_patch_bbox.min.y));

        const mi::Size nb_polygon_verts = m_polygon_of_interest.size();
        
        for (mi::Size i = 0; i < nb_polygon_verts; ++i)
        {
            const mi::math::Vector<mi::Float32, 2>& polygon_vert = m_polygon_of_interest[i];
            if (polygon_vert.x < min_value.x) min_value.x = polygon_vert.x;
            if (polygon_vert.y < min_value.y) min_value.y = polygon_vert.y;
            if (polygon_vert.x > max_value.x) max_value.x = polygon_vert.x;
            if (polygon_vert.y > max_value.y) max_value.y = polygon_vert.y;
        }

        // Clamp to bounding box
        if (min_value.x < ij_patch_bbox.min.x) min_value.x = static_cast<mi::Float32>(ij_patch_bbox.min.x);
        if (min_value.y < ij_patch_bbox.min.y) min_value.y = static_cast<mi::Float32>(ij_patch_bbox.min.y);
        if (max_value.x > ij_patch_bbox.max.x) max_value.x = static_cast<mi::Float32>(ij_patch_bbox.max.x);
        if (max_value.y > ij_patch_bbox.max.y) max_value.y = static_cast<mi::Float32>(ij_patch_bbox.max.y);

        const mi::Sint32 min_value_y =
            static_cast< mi::Uint32 >(mi::math::ceil(min_value.y)  - ij_patch_bbox.min.y);
        const mi::Sint32 max_value_y =
            static_cast< mi::Uint32 >(mi::math::floor(max_value.y) - ij_patch_bbox.min.y);
        std::vector<mi::Float32> intersecton_spans;

        for(mi::Sint32 y=min_value_y; y<max_value_y; ++y)
        {
            compute_intersection(static_cast<mi::Float32>(y+ij_patch_bbox.min.y), m_polygon_of_interest, intersecton_spans);
            mi::Float32* vertices = &elevation_values[(ij_patch_bbox.max.x-ij_patch_bbox.min.x) * y];

            const mi::Size nb_intersections = intersecton_spans.size();
            for(mi::Size i=0; i+1<nb_intersections; i+=2)
            {
                mi::Sint32 from = static_cast< mi::Sint32 >(intersecton_spans[i]);
                mi::Sint32 to   = static_cast< mi::Sint32 >(intersecton_spans[i+1]);

                mi::Float32* VERT = vertices;
                // for(mi::Sint32 scan_x = ij_patch_bbox.min.x; 
                //     scan_x <= static_cast< mi::Sint32 >(ij_patch_bbox.max.x); ++k)
                // {
                //     if(k>=from && k<=to)
                //         *VERT = -1;
                //     VERT++;
                // }
                // DELETEME value range -> bbox
                // TODO: Need a small test (an exact polygon should specified)
                for(mi::Sint32 scan_x = ij_patch_bbox.min.x; 
                    scan_x < static_cast< mi::Sint32 >(ij_patch_bbox.max.x); ++scan_x)
                {
                    if((scan_x>=from) && (scan_x<=to)){
                        *VERT = nv::index::IRegular_heightfield_patch::get_hole_value();
                    }
                    VERT++;
                }
            }
            intersecton_spans.clear();
        }

        return true;
    }

    virtual bool compute(
        const mi::math::Bbox_struct<mi::Uint32, 2>&     ij_patch_bbox,
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
    std::vector<mi::math::Vector_struct<mi::Float32, 2> >   m_polygon_of_interest;
};

class Heightfield_patch_delete_operation :
    public mi::base::Interface_implement<nv::index::IRegular_heightfield_compute_task>
{
public:
    Heightfield_patch_delete_operation(
        const mi::math::Bbox_struct<mi::Sint32, 2>& roi) : m_roi(roi)
    {
        // empty
    }

public:
    // Implemented IRegular_heightfield_compute_task

    virtual void get_region_of_interest_for_compute(
        mi::math::Bbox_struct<mi::Sint32, 2>& roi) const
    {
        roi = m_roi;
    }

    // FIXME: Updated to bounding box access, but not tested
    // TODO: need a simple exact test
    virtual bool compute(
        const mi::math::Bbox_struct<mi::Sint32, 2>&     ijk_patch_bbox,
        mi::Float32*                                    elevation_values,
        mi::neuraylib::IDice_transaction*               dice_transaction) const
    {
        const mi::Sint32 nb_elevation_values =
            (ijk_patch_bbox.max.x - ijk_patch_bbox.min.x) * (ijk_patch_bbox.max.y - ijk_patch_bbox.min.y);
        assert(nb_elevation_values >= 0);
        
        for (mi::Sint32 i = 0; i < nb_elevation_values; ++i)
        {
            elevation_values[i] = nv::index::IRegular_heightfield_patch::get_hole_value();
        }

        return true;
    }

    // This override is not necessary, but without this override, ICC
    // compiler warns.
    virtual bool compute(
        const mi::math::Bbox_struct<mi::Sint32, 2>& ij_patch_bbox,
        mi::Float32*                                elevation_values,
        mi::math::Vector_struct<mi::Float32, 3>*    normal_values,
        mi::neuraylib::IDice_transaction*           dice_transaction) const
    {
        return false;
    }

    // This override is not necessary, but without this override, ICC
    // compiler warns.
    virtual bool user_defined_normal_computation() const { return false; }

private:
    mi::math::Bbox_struct<mi::Sint32, 2>    m_roi;
};

//----------------------------------------------------------------------
Distributed_heightfield_elevation_delete::Distributed_heightfield_elevation_delete(
    const mi::neuraylib::Tag&      session_tag,
    const mi::neuraylib::Tag&      heightfield_tag,
    const std::vector<mi::Uint32>& cluster_hosts)
 :  m_session_tag(session_tag),
    m_heightfield_tag(heightfield_tag),
    m_cluster_hosts(cluster_hosts)
{
    m_polygon_of_interest.resize(0);
}

//----------------------------------------------------------------------
Distributed_heightfield_elevation_delete::Distributed_heightfield_elevation_delete(
    const mi::neuraylib::Tag&                                       session_tag,
    const mi::neuraylib::Tag&                                       heightfield_tag,
    const std::vector<mi::math::Vector_struct<mi::Float32, 2> >&    polygon,
    const std::vector<mi::Uint32>&                                  cluster_hosts)
 :  m_session_tag(session_tag),
    m_heightfield_tag(heightfield_tag),
    m_cluster_hosts(cluster_hosts),
    m_polygon_of_interest(polygon)
{
    // empty
}

//----------------------------------------------------------------------
Distributed_heightfield_elevation_delete::~Distributed_heightfield_elevation_delete()
{
    // empty
}

//-----------------------------------------------------------------------------------------
void Distributed_heightfield_elevation_delete::assign_fragments_to_hosts(
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
void Distributed_heightfield_elevation_delete::serialize(
    mi::neuraylib::ISerializer* serializer) const
{
    serializer->write(&m_session_tag.id, 1);
    serializer->write(&m_heightfield_tag.id, 1);

    mi::Uint32 nb_elements = mi::Uint32(m_polygon_of_interest.size());
    serializer->write(&nb_elements, 1);
    for(mi::Uint32 i=0; i<nb_elements; ++i)
        serializer->write(&m_polygon_of_interest[i].x, 2);

    nb_elements = mi::Uint32(m_cluster_hosts.size());
    serializer->write(&nb_elements, 1);
    for(mi::Uint32 i=0; i<nb_elements; ++i)
        serializer->write(&m_cluster_hosts[i], 1);
}

//----------------------------------------------------------------------
void Distributed_heightfield_elevation_delete::deserialize(
    mi::neuraylib::IDeserializer* deserializer)
{
    deserializer->read(&m_session_tag.id, 1);
    deserializer->read(&m_heightfield_tag.id, 1);

    mi::Uint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_polygon_of_interest.resize(nb_elements);
    for(mi::Uint32 i=0; i<nb_elements; ++i)
        deserializer->read(&m_polygon_of_interest[i].x, 2);

    nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_cluster_hosts.resize(nb_elements);
    for(mi::Uint32 i=0; i<nb_elements; ++i)
        deserializer->read(&m_cluster_hosts[i], 1);
}

//----------------------------------------------------------------------
// FIXME: need tests due to value range -> bounding box change
void Distributed_heightfield_elevation_delete::edit_heightfield_data(
    mi::neuraylib::IDice_transaction*              dice_transaction,
    nv::index::IRegular_heightfield_data_locality* locality,
    mi::Uint32                                     host_id)
{
    const mi::Uint32 nb_patch_edits = static_cast<mi::Uint32>(locality->get_nb_bounding_box(host_id));

    for(mi::Uint32 i=0; i<nb_patch_edits; ++i)
    {
        mi::math::Bbox<mi::Sint32, 3> roi_3D = locality->get_bounding_box(host_id, i);
        mi::math::Bbox<mi::Sint32, 2> roi;
        roi.min.x = roi_3D.min.x;
        roi.min.y = roi_3D.min.y;
        roi.max.x = roi_3D.max.x;
        roi.max.y = roi_3D.max.y;

        mi::base::Handle<nv::index::IRegular_heightfield_data_edit> heightfield_data_edit(
            locality->create_data_edit(dice_transaction, host_id, i));

        DEBUG_LOG << "Deleting heightfield elevations on host: " << host_id;
        if(!m_polygon_of_interest.empty())
        {
            Heightfield_polygon_based_delete_operation delete_operation(m_polygon_of_interest);
            heightfield_data_edit->edit(&delete_operation, dice_transaction);
        }
        else
        {
            Heightfield_patch_delete_operation delete_operation(roi);
            heightfield_data_edit->edit(&delete_operation, dice_transaction);
        }
    }
}


//----------------------------------------------------------------------
void Distributed_heightfield_elevation_delete::execute_fragment(
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

    edit_heightfield_data(dice_transaction, heightfield_data_locality.get(), cluster_host);
}

//----------------------------------------------------------------------
void Distributed_heightfield_elevation_delete::execute_fragment_remote(
    mi::neuraylib::ISerializer*          serializer,
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

    edit_heightfield_data(dice_transaction, heightfield_data_locality.get(), cluster_host);
}

//----------------------------------------------------------------------
void Distributed_heightfield_elevation_delete::receive_remote_result(
    mi::neuraylib::IDeserializer*        deserializer,
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Size                                        index,
    mi::Size                                        count)
{
    // nothing to do .....
}

//----------------------------------------------------------------------
