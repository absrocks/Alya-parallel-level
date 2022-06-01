/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "snapping_tool.h"
#include "heightfield_workflow_functionality.h"

//----------------------------------------------------------------------
bool Snapping_tool::adjust_pick(
    const mi::math::Vector_struct<mi::Float32, 3>& pick_position,
    nv::index::IRegular_volume_data_access*        subsurface_data_access,
    mi::neuraylib::IDice_transaction*              dice_transaction,
    mi::math::Vector_struct<mi::Float32, 3>&       result_pick) const
{
    // record the pick coordinate for test
    Heightfield_workflow_functionality::instance()->set_last_xyz_world_pick_position(pick_position);
    Heightfield_workflow_functionality::instance()->set_last_xyz_volume_pick_position(
        mi::math::Vector<mi::Float32, 3>(-1.0f, -1.0f, -1.0f));

    if(subsurface_data_access)
    {
        // Access the volume scene element
        const mi::neuraylib::Tag& scene_element_tag = subsurface_data_access->get_scene_element();
        mi::base::Handle<const nv::index::IRegular_volume> volume_tag(
            dice_transaction->access<const nv::index::IRegular_volume>(scene_element_tag));

        // Transform manual pick into the volume's object space (ijk space).
        mi::math::Vector<mi::Float32, 3> ijk_pick_position = pick_position;
        ijk_pick_position.x = mi::math::round(ijk_pick_position.x);
        ijk_pick_position.y = mi::math::round(ijk_pick_position.y);
        ijk_pick_position.z = mi::math::round(ijk_pick_position.z);

        // Compute the query bounds, i.e., the column below and above the intersection with the slice.
        const mi::math::Bbox_struct<mi::Float32, 3>& volume_bbox = volume_tag->get_IJK_region_of_interest();
        mi::math::Bbox_struct<mi::Uint32, 3> query_bounds;
        query_bounds.min.x = static_cast< mi::Uint32 >(ijk_pick_position.x);
        query_bounds.min.y = static_cast< mi::Uint32 >(ijk_pick_position.y);
        query_bounds.min.z = 0;
        query_bounds.max.x = static_cast< mi::Uint32 >(ijk_pick_position.x + 1);
        query_bounds.max.y = static_cast< mi::Uint32 >(ijk_pick_position.y + 1);
        query_bounds.max.z = static_cast< mi::Uint32 >(volume_bbox.max.z);

        // Bounding boxes must have a volume
        assert(query_bounds.min.x < query_bounds.max.x);
        assert(query_bounds.min.y < query_bounds.max.y);
        assert(query_bounds.min.z < query_bounds.max.z);
        
        subsurface_data_access->access(query_bounds, dice_transaction);
        const mi::math::Bbox_struct<mi::Uint32, 3>& effective_bounds = subsurface_data_access->get_bounding_box();

        DEBUG_LOG << "Query bounds:     " << query_bounds;
        DEBUG_LOG << "Effective bounds: " << effective_bounds;

        const mi::base::Handle<const nv::index::IRegular_volume_data> volume_data(
            subsurface_data_access->get_volume_data());
        const mi::base::Handle<const nv::index::IRegular_volume_data_uint8> volume_data_uint8(
            volume_data->get_interface<const nv::index::IRegular_volume_data_uint8>());

        if (!volume_data_uint8) {
            ERROR_LOG << "Snapping tool failed: data access on non-uint8 volume type.";
            return false;
        }


        if (volume_data_uint8->get_voxel_data())
        {
            const mi::Sint32 k_snap_min =
                static_cast< mi::Sint32 >(((ijk_pick_position.z - m_k_range_min)<=0)
                                          ? 0 : (ijk_pick_position.z - m_k_range_min));
            const mi::Sint32 k_snap_max =
                static_cast< mi::Sint32 >(((ijk_pick_position.z + m_k_range_max)>=volume_bbox.max.z)
                                          ? volume_bbox.max.z : (ijk_pick_position.z + m_k_range_max));

            const mi::Uint8* ptr = volume_data_uint8->get_voxel_data();
            mi::Uint8 amplitude_value = ptr[mi::Uint32(ijk_pick_position.z)];
            mi::Uint32 k_position = static_cast< mi::Sint32 >(ijk_pick_position.z);

            for (mi::Sint32 k = k_snap_min; k <= k_snap_max; ++k)
            {
                const mi::Uint8 value = ptr[k];
                DEBUG_LOG << "K: " << k << " value: " << mi::Uint32(value);
                if(value<amplitude_value)
                {
                    amplitude_value = value;
                    k_position = k;
                }
            }
            ijk_pick_position.z = static_cast<mi::Float32>(k_position);
        }

        // Transform new pick location back to xyz space and retrun
        INFO_LOG << "Snapping tool returns a computed position: " << ijk_pick_position;

        Heightfield_workflow_functionality::instance()->set_last_xyz_volume_pick_position(
            ijk_pick_position); // record the computed result for test

        result_pick.x = ijk_pick_position.x;
        result_pick.y = ijk_pick_position.y;
        result_pick.z = ijk_pick_position.z;
        return true;
    }
    else
    {
        // Simply snap the intersection to the nearest integer valued position
        result_pick.x = mi::math::round(pick_position.x);
        result_pick.y = mi::math::round(pick_position.y);
        result_pick.z = mi::math::round(pick_position.z);
        INFO_LOG << "Snapping tool returns a snapped position: " << result_pick;
        return false;
    }
}

//----------------------------------------------------------------------

