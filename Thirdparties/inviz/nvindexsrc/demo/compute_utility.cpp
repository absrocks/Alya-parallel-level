/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "compute_utility.h"

#include <nv/index/iscene.h>
#include <nv/index/iregular_heightfield.h>
#include <nv/index/idistributed_data_locality.h>

#include <cassert>
#include <sstream>
#include <iterator>
#include <iostream>
#include <cstdio>

#include "common/common_utility.h"
#include "common/distributed_voxel_value_change.h"
#include "common/forwarding_logger.h"
#include "common/string_dict.h"
#include "common/type_conversion_utility.h"

#include "volume_bordermap_database_element.h"
#include "distributed_amplitude_map_job.h"
#include "distributed_volume_bordermap_generation.h"
#include "distributed_volume_filter.h"
#include "heightfield_workflow_functionality.h"
#include "nvindex_appdata.h"
#include "scene_utility.h"
#include "snapping_tool.h"
#include "volume_data_filter.h"

#include "common/volume_replica_generation.h"

//----------------------------------------------------------------------
std::string get_unique_dice_name(const std::string & base_dice_name, 
                                 mi::neuraylib::IDice_transaction * dice_transaction)
{
    assert(!base_dice_name.empty());
    assert(dice_transaction != 0);    

    const mi::Sint32 max_i = 1000000;
    for(mi::Sint32 i = 0; i < max_i; ++i){
        std::stringstream sstr;
        sstr << base_dice_name << "_" << i;
        mi::neuraylib::Tag check_tag = dice_transaction->name_to_tag(sstr.str().c_str());
        if(!check_tag.is_valid()){
            return sstr.str();
        }
    }
    ERROR_LOG << "Can't find unique name of [" << base_dice_name << "_0-" << max_i << "]";

    return base_dice_name;
}

//----------------------------------------------------------------------
void setup_parallel_fragment(std::vector< mi::Uint32 > const & cluster_hosts,
                             mi::Uint32 const parallel_count,
                             std::vector< mi::Sint32 > & fragment_to_host_id,
                             std::vector< mi::Sint32 > & fragment_to_host_local_thread_id)

{
    std::stringstream sstr;

    assert(!cluster_hosts.empty());
    assert(parallel_count       > 0);

    // set up M: fragment index ->host id
    fragment_to_host_id.clear();
    for(std::vector< mi::Uint32 >::const_iterator hi = cluster_hosts.begin();
        hi != cluster_hosts.end(); ++hi)
    {
        assert((*hi) != 0);
        for(mi::Uint32 i = 0; i < parallel_count; ++i){
            fragment_to_host_id.push_back(*hi);
        }
    }
    sstr << "M: fi -> hi = ";
    std::copy(fragment_to_host_id.begin(), fragment_to_host_id.end(),
              std::ostream_iterator< mi::Sint32 >(sstr, " "));
    sstr << "\n";

    // set up M: fragment index -> host local thread id
    fragment_to_host_local_thread_id.clear();
    size_t const host_count = cluster_hosts.size();
    for(size_t i = 0; i < host_count; ++i){
        for(mi::Uint32 i = 0; i < parallel_count; ++i){
            fragment_to_host_local_thread_id.push_back(i);
        }
    }
    sstr << "debug: M: fi -> ti = ";
    std::copy(fragment_to_host_local_thread_id.begin(),
              fragment_to_host_local_thread_id.end(),
              std::ostream_iterator< mi::Sint32 >(sstr, " "));

    DEBUG_LOG << sstr.str();
}

//----------------------------------------------------------------------
void heightfield_workflow_manual_picking(
    mi::base::Handle< mi::neuraylib::IScope > & scope,
    mi::math::Vector_struct< mi::Uint32, 2 > const & pick_location,
    const nv::index::IIndex_canvas*                  pick_canvas
    )
{
    assert(scope.is_valid_interface());
    // check the workflow is in the manual pick state.
    assert(Heightfield_workflow_functionality::instance()->is_heightfield_workflow_manual_pick());

    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        scope->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    {
        // Apply heightfield manual pick operation using a user-defined
        // snapping algorithm (a very simple one).
        Snapping_tool snapping_tool;
        Heightfield_workflow_functionality::instance()->
            process_pick(pick_location, pick_canvas, &snapping_tool, dice_transaction.get());
    }
    dice_transaction->commit();
}


//----------------------------------------------------------------------
void heightfield_workflow_append_vertex_to_bounding_polygon(
    mi::base::Handle< mi::neuraylib::IScope >&        scope,
    mi::base::Handle<nv::index::IIndex_scene_query>&  iindex_scene_query,
    const mi::neuraylib::Tag&                         session_tag,
    const mi::math::Vector_struct< mi::Uint32, 2 >&   pick_location,
    const nv::index::IIndex_canvas*                   pick_canvas
)
{
    assert(scope.is_valid_interface());
    // check the workflow is in the  polygon delete state.
    assert(     Heightfield_workflow_functionality::instance()->is_heightfield_workflow_delete_operation()
             || Heightfield_workflow_functionality::instance()->is_heightfield_workflow_elevation_change_operation()
             || Heightfield_workflow_functionality::instance()->is_heightfield_workflow_gridding_operation() );

    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        scope->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    mi::base::Handle<nv::index::IScene_pick_results> scene_pick_results(
        iindex_scene_query->pick(
            pick_location,
            pick_canvas,
            session_tag,
            dice_transaction.get()));
    dice_transaction->commit();

    const mi::Uint32 nb_results = scene_pick_results->get_nb_results();
    if(nb_results>0)
    {
        // INFO_LOG << "Number of intersections: " << nb_results;
        for(mi::Uint32 i=0; i<nb_results; ++i)
        {
            mi::base::Handle<nv::index::IScene_pick_result> result(scene_pick_results->get_result(i));

            if(Heightfield_workflow_functionality::instance()->get_heightfield() == result->get_scene_element())
            {
                const mi::math::Vector_struct<mi::Float32, 3>& intersection = result->get_intersection();

                if(Heightfield_workflow_functionality::instance()->add_vertex_to_bounding_polygon(intersection)){
                    INFO_LOG << "Adding a vertex to bounding polygon: \n"
                             << "\t Heightfield (tag)    " << result->get_scene_element().id << "\n"
                             << "\t IJK position:    " << intersection;
                }
                else{
                    WARN_LOG << "Rejected to add the same vertex to the last one.";
                }

                return;
            }
        }
    }
    INFO_LOG << "No valid heightfield intersection, i.e., no vertex added to heightfield delete polygon.";
}

//----------------------------------------------------------------------
void heightfield_workflow_delete_polygon(
    mi::base::Handle<mi::neuraylib::IScope>&    scope,
    bool                                        use_convex_hull)
{
    assert(scope.is_valid_interface());

    // check the workflow is in the polygon delete state.
    assert(Heightfield_workflow_functionality::instance()->is_heightfield_workflow_delete_operation());

    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        scope->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    {
        Heightfield_workflow_functionality::instance()->delete_polygon(use_convex_hull, dice_transaction.get());
    }
    dice_transaction->commit();
}

//----------------------------------------------------------------------
void heightfield_workflow_change_elevation_values(
    mi::base::Handle<mi::neuraylib::IScope>& scope,
    mi::Float32                              new_elevation_value,
    bool                                     is_scale_op,
    bool                                     use_convex_hull)
{
    assert(scope.is_valid_interface());

    // check the workflow is in the elevation change state.
    assert(Heightfield_workflow_functionality::instance()->is_heightfield_workflow_elevation_change_operation());

    // get dice transaction
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        scope->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    {
        // update heightfield data
        Heightfield_workflow_functionality::instance()->update_heightfield_elevation_values(
            new_elevation_value,
            is_scale_op,
            use_convex_hull,
            dice_transaction.get());
    }
    dice_transaction->commit();
}

//----------------------------------------------------------------------
void heightfield_workflow_gridding(
    mi::base::Handle<mi::neuraylib::IScope>&    scope,
    bool                                        use_convex_hull)
{
    assert(scope.is_valid_interface());

    // check the workflow is in the polygon delete state.
    assert(Heightfield_workflow_functionality::instance()->is_heightfield_workflow_gridding_operation());

    // get dice transaction
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        scope->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    {
        // delete heightfield data
        Heightfield_workflow_functionality::instance()->gridding(
            use_convex_hull,
            dice_transaction.get());
    }
    // commit changes
    dice_transaction->commit();
}

//----------------------------------------------------------------------
void heightfield_workflow_autotracking(
    mi::base::Handle<mi::neuraylib::IScope>& scope)
{
    // autotracking start
    Heightfield_workflow_functionality::instance()->autotracking(scope);
}

//----------------------------------------------------------------------
bool generate_heightfield_amplitude_map(
    const std::string&                           com_str,
    mi::base::Handle<mi::neuraylib::IScope>&     scope,
    mi::base::Handle<nv::index::IIndex_session>& iindex_session,
    const mi::neuraylib::Tag&                    session_tag)
{
    std::string const func_name = "generate_heightfield_amplitude_map:";
    std::string const com_name = "ha_map ";

    std::string const param = (com_str.substr(com_name.size()));
    if(param.size() == 0){
        ERROR_LOG << ("empty parameter: " + com_str);
        return false;
    }
    mi::Sint32 heightfield_idx = -1;
    mi::Sint32 volume_idx  = -1;
    {
        std::stringstream instr(param);
        instr >> heightfield_idx >> volume_idx;
        if(instr.fail()){
            ERROR_LOG << ("invalid parameters: " + com_str);
            return false;
        }
    }

    DEBUG_LOG << func_name << com_name
              << ", heightfield idx = " << heightfield_idx << ", volume_idx = " << volume_idx;

    assert(scope.is_valid_interface());
    assert(session_tag.is_valid());

    mi::neuraylib::Tag const heightfield_tag = Nvindex_AppData::instance()->get_heightfield_tag_from_idx(heightfield_idx);
    if(!heightfield_tag.is_valid()){
        ERROR_LOG << func_name << "no such heightfield_idx [" << heightfield_idx << "]";
        return false;
    }

    mi::neuraylib::Tag const volume_tag = Nvindex_AppData::instance()->get_volume_tag_from_idx(volume_idx);
    if(!volume_tag.is_valid()){
        ERROR_LOG << func_name << "no such volume_idx [" << volume_idx << "]";
        return false;
    }
    DEBUG_LOG << func_name
              << ", heightfield tag = " << heightfield_tag.id << ", volume_tag = " << volume_tag.id;

    // get dice transaction
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        scope->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());

    {
        // heightfield_interaction is new-ed.
        mi::base::Handle< nv::index::IHeightfield_interaction > heightfield_interaction(
            iindex_session->create_heightfield_interaction_interface(
                heightfield_tag, session_tag, dice_transaction.get()));
        assert(heightfield_interaction.is_valid_interface() != 0);

        mi::math::Bbox_struct<mi::Uint32, 3> heightfield_bbox_result;
        bool const is_edit = true;
        mi::base::Handle<nv::index::IRegular_heightfield_data_locality> heightfield_data_locality(
            new_whole_heightfield_distribution_layout(session_tag,
                                                      heightfield_tag,
                                                      is_edit,
                                                      dice_transaction.get(),
                                                      heightfield_bbox_result));
        assert(heightfield_data_locality.is_valid_interface());

        std::vector<mi::Uint32> cluster_host_ids;
        mi::Uint32 const nb_cluster_hosts = heightfield_data_locality->get_nb_cluster_nodes();
        DEBUG_LOG << "nb_cluster_hosts = " << nb_cluster_hosts;
        for(mi::Uint32 i=0; i<nb_cluster_hosts; ++i)
        {
            DEBUG_LOG << "host id [" << i << "] = " << heightfield_data_locality->get_cluster_node(i);
            cluster_host_ids.push_back(heightfield_data_locality->get_cluster_node(i));
        }

        // Choose a distributed computing algorithm
        Distributed_amplitude_map_job distributed_amplitude_map_job(
            session_tag,
            heightfield_tag,
            volume_tag,
            cluster_host_ids);
        heightfield_interaction->invoke_computing(&distributed_amplitude_map_job, dice_transaction.get());
        distributed_amplitude_map_job.save_amplitude_map("amplitude_map.ppm");
    }
    dice_transaction->commit();

    return true;
}

//----------------------------------------------------------------------
bool volume_set_amplitude_value(
    bool                      is_host_assign_mode,
    mi::Uint8                 amplitude_value,
    const mi::neuraylib::Tag& session_tag,
    const mi::neuraylib::Tag& volume_tag,
    mi::neuraylib::IDice_transaction * dice_transaction)
{
    assert(dice_transaction != 0);
    assert(session_tag.is_valid());    
    assert(volume_tag.is_valid());    

    // Work on the first volume only
    mi::math::Bbox<mi::Float32, 3> roi;
    mi::math::Vector_struct<mi::Uint32, 3> volume_size;
    {
        // get_cluster_host_containing_volume() access to the volume,
        // therefore if there is no scope and edit<> causes an error.
        // (Here is an access, so you can actually remove this scope,
        // but it is an error prone style.)
        mi::base::Handle<const nv::index::IRegular_volume> volume(
            dice_transaction->access<const nv::index::IRegular_volume>(volume_tag));
        if(!volume.is_valid_interface()){
            ERROR_LOG << "volume_set_amplitude_value: invalid volume_tag: " << volume_tag;
            return false;
        }
    
        // Use current per-volume region of interest as bounds for this operation
        roi = volume->get_IJK_region_of_interest();
        volume_size = volume->get_volume_size();
    }
    const mi::math::Bbox<mi::Uint32, 3> bounds(roi);

    // get all host ids in the cluster
    std::vector<mi::Uint32> cluster_host_ids;
    get_cluster_host_containing_volume(session_tag, volume_tag, bounds, cluster_host_ids,
                                       dice_transaction);
    
    // Set up distributed algorithm
    mi::base::Handle< nv::index::IDistributed_compute_algorithm >
        algo(new nv::index_common::Distributed_voxel_value_change(
                 session_tag, volume_tag, is_host_assign_mode,  amplitude_value,
                 bounds, volume_size, cluster_host_ids));

    // Start the fragmented job
    dice_transaction->execute_fragmented(algo.get(), algo->get_nb_of_fragments());

    return true;
}

//----------------------------------------------------------------------
bool volume_inset_filter_amplitude_values(
    const mi::neuraylib::Tag & session_tag,
    const mi::neuraylib::Tag & volume_tag,
    const std::string & filter_type_str,
    const mi::Sint32 parallel_count, // FIXME: may not needed. Use slice number (or is_parallel_execution boolean)
    mi::neuraylib::IDice_transaction * dice_transaction
    )
{
    mi::neuraylib::Tag bordermap_tag = mi::neuraylib::NULL_TAG;
    mi::math::Bbox<mi::Uint32, 3> clipped_bbox;
    std::vector<mi::Uint32> cluster_host_ids;

    assert(dice_transaction != 0);

    // necessary info to process volume
    if(!volume_tag.is_valid()){
        ERROR_LOG << "volume_inset_filter_amplitude_values: invalid volume tag ["
                  << volume_tag << "].";
        return false;
    }
    mi::math::Bbox< mi::Sint64, 3 > const clipped_src_ijk_roi_bbox =
        get_volume_local_IJK_ROI_bbox(volume_tag, dice_transaction);

    if(clipped_src_ijk_roi_bbox.empty()){
        ERROR_LOG << "[volume_data_filter]: can not find source volume ROI.";
        return false;
    }

    clipped_bbox = nv::index_common::convert_bbox_sint64_3_to_bbox_uint32_3(clipped_src_ijk_roi_bbox);

    // get all host ids in the cluster
    get_cluster_host_containing_volume(session_tag, volume_tag, clipped_bbox, cluster_host_ids,
                                       dice_transaction);
    if(cluster_host_ids.empty()){
        ERROR_LOG << "No host found. either volume_tag = " << volume_tag
                  << ", bbox = " << clipped_bbox
                  << " is not valid, or data is not loaded.";
        return false;
    }

    mi::Sint32 const ftype = Volume_data_filter::get_filter_type(filter_type_str);
    if(ftype == Volume_data_filter::VDFT_Count){
        ERROR_LOG << "volume_inset_filter_amplitude_values: unknown filter ["
                  << filter_type_str << "].";
        return false;
    }

    // distributed algorithm: step 1. bordermap creation job
    mi::base::Handle< nv::index::IDistributed_compute_algorithm >
        gen_border(new Distributed_volume_bordermap_generation(
                       session_tag,
                       volume_tag,
                       ftype,
                       clipped_bbox,
                       cluster_host_ids));
    assert(gen_border.is_valid_interface());

    // Start the fragmented job
    dice_transaction->execute_fragmented(gen_border.get(), gen_border->get_nb_of_fragments());

    INFO_LOG << "Distributed_volume_bordermap_generation: create bordermap done";

    // create the border set and store to the DB
    {
        mi::base::Handle< Volume_bordermap_database_element > bordermap(
            new Volume_bordermap_database_element());
        bordermap_tag = dice_transaction->store(bordermap.get(), mi::neuraylib::NULL_TAG,
                                                "bordermap");
        assert(bordermap_tag.is_valid());
    }
    // INFO_LOG << "global bordermap_tag = " << bordermap_tag.id;

    // get the created border result from the job
    {
        mi::base::Handle< Volume_bordermap_database_element > bordermap(
            dice_transaction->edit< Volume_bordermap_database_element >(bordermap_tag));
        assert(bordermap.is_valid_interface());

        std::vector< mi::neuraylib::Tag > result_tag_vec =
            static_cast< Distributed_volume_bordermap_generation * >(
                gen_border.get())->get_job_result_bordermap_tag_vec();
        for(std::vector< mi::neuraylib::Tag >::const_iterator ti = result_tag_vec.begin();
            ti != result_tag_vec.end(); ++ti)
        {
            mi::base::Handle< Volume_bordermap_database_element const > job_bmap(
                dice_transaction->access< Volume_bordermap_database_element const >(*ti));
            assert(job_bmap.is_valid_interface());
            bordermap->insert_bordermap(job_bmap.get());
        }
        // DEBUG_LOG << "dump the global bordermap: " << bordermap->to_string();
    }

    mi::Float64 const start_time = nv::index_common::get_time();
    {
        mi::Sint32 const ftype = Volume_data_filter::get_filter_type(filter_type_str);
        assert(parallel_count > 0);

        DEBUG_LOG << "parallel_count = " << parallel_count << ", cluster_host_ids.size() = "
                  << cluster_host_ids.size();

        // distributed algorithm: step 2. apply inset filter with bordermap
        mi::base::Handle< nv::index::IDistributed_compute_algorithm >
            vol_filter(new Distributed_volume_filter(
                           session_tag,
                           volume_tag,
                           bordermap_tag,
                           ftype,
                           clipped_bbox,
                           cluster_host_ids,
                           parallel_count));
        assert(vol_filter.is_valid_interface());

        // Start the fragmented job
        dice_transaction->execute_fragmented(vol_filter.get(), vol_filter->get_nb_of_fragments());

        INFO_LOG << "Distributed_volume_bordermap_generation: apply filter done";
    }
    mi::Float64 const elapsed_time = nv::index_common::get_time() - start_time;
    INFO_LOG << "filter elapsed time = " << elapsed_time;

    // clean up the border elements
    assert(bordermap_tag.is_valid());
    {
        mi::base::Handle< Volume_bordermap_database_element > bordermap(
            dice_transaction->edit< Volume_bordermap_database_element >(
                bordermap_tag));
        assert(bordermap.is_valid_interface());
        bordermap->clear_bordermap(dice_transaction);
    }

    return true;
}

//----------------------------------------------------------------------
mi::neuraylib::Tag generate_volume_by_copy(
    const mi::neuraylib::Tag & session_tag,
    const mi::neuraylib::Tag & src_volume_tag,
    const std::string & dst_volume_name,
    const mi::math::Vector<mi::Float32, 3> & tr_vec,
    mi::neuraylib::IDice_transaction * dice_transaction)
{
    assert(session_tag.is_valid());
    assert(src_volume_tag.is_valid());
    assert(!dst_volume_name.empty());
    assert(dice_transaction != 0);
    
    // get source volume size
    mi::base::Handle<const nv::index::IRegular_volume> src_volume(
        dice_transaction->access<const nv::index::IRegular_volume>(src_volume_tag));
    if(!src_volume.is_valid_interface()){
        ERROR_LOG << "No such volume in the scene, tag: " << src_volume_tag;
        return mi::neuraylib::NULL_TAG;
    }
    const mi::math::Vector_struct<mi::Uint32, 3> src_volume_size = src_volume->get_volume_size();
    const mi::math::Vector_struct<mi::Uint32, 3> zero = { 0, 0, 0, };
    mi::math::Bbox_struct<mi::Uint32, 3>   src_ijk_bbox;
    src_ijk_bbox.min = zero;
    src_ijk_bbox.max = src_volume_size;

    // create copy importer
    nv::index_common::Volume_replica_generation* volume_generator_cb = 
        new nv::index_common::Volume_replica_generation(session_tag, src_volume_tag, src_ijk_bbox);

    // get the scene for edit
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<const nv::index::ISession>(session_tag));
    assert(session.is_valid_interface());
    mi::base::Handle<nv::index::IScene> scene(
        dice_transaction->edit<nv::index::IScene>(session->get_scene()));
    assert(scene.is_valid_interface());

    // Create a volume scene element that represents a regular volume.
    const mi::math::Vector_struct<mi::Float32, 3> scale = { 1.0f, 1.0f, 1.0f, };
    const mi::Float32 rotate_k = 0.0f;
    mi::base::Handle<nv::index::IRegular_volume> dst_volume(
        scene->create_volume(scale, rotate_k, tr_vec, src_volume_size, volume_generator_cb, dice_transaction));
    assert(dst_volume.is_valid_interface());

    // Assign the same colormap to the source volume
    dst_volume->assign_colormap(src_volume->assigned_colormap());
    dst_volume->set_name(dst_volume_name.c_str());
    dst_volume->set_enabled(true);

    const std::string dst_volume_dice_name = get_unique_dice_name(dst_volume_name, dice_transaction);
                
    const mi::neuraylib::Tag dst_volume_tag = 
        dice_transaction->store_for_reference_counting(dst_volume.get(), 
                                                       mi::neuraylib::NULL_TAG, 
                                                       dst_volume_dice_name.c_str());
    assert(dst_volume_tag.is_valid());

    return dst_volume_tag;
}

//----------------------------------------------------------------------
mi::neuraylib::Tag add_volume_to_the_scene_root(
    const mi::neuraylib::Tag & session_tag,
    const mi::neuraylib::Tag & adding_volume_tag,
    mi::neuraylib::IDice_transaction * dice_transaction)
{
    assert(session_tag.is_valid());
    assert(adding_volume_tag.is_valid());
    assert(dice_transaction != 0);

    // get the scene for edit
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<const nv::index::ISession>(session_tag));
    assert(session.is_valid_interface());
    mi::base::Handle<nv::index::IScene> scene_edit(
        dice_transaction->edit<nv::index::IScene>(session->get_scene()));
    assert(scene_edit.is_valid_interface());

    mi::base::Handle<nv::index::IStatic_scene_group> static_group_node(
        scene_edit->create_scene_group<nv::index::IStatic_scene_group>());
    assert(static_group_node.is_valid_interface());

    static_group_node->append(adding_volume_tag, dice_transaction);

    const std::string group_dice_name = get_unique_dice_name("IStatic_scene_group", dice_transaction);

    const mi::neuraylib::Tag group_tag =
        dice_transaction->store_for_reference_counting(static_group_node.get(),
                                                       mi::neuraylib::NULL_TAG,
                                                       group_dice_name.c_str());

    // append the group to the scene root
    scene_edit->append(group_tag, dice_transaction);

    return group_tag;
}

//----------------------------------------------------------------------
bool test_verify_synthetic_volume_contents(
    std::string const & param,
    mi::base::Handle<mi::neuraylib::IScope> & scope,
    const mi::neuraylib::Tag&                session_tag)
{
    // start transaction
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        scope->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    {
        // Access the session
        mi::base::Handle<const nv::index::ISession> session(
            dice_transaction->
            access<const nv::index::ISession>(session_tag));
        assert(session.is_valid_interface());
        // Access the data access factory
        const mi::neuraylib::Tag& data_access_tag = session->get_data_access_factory();
        mi::base::Handle<const nv::index::IDistributed_data_access_factory> access_factory(
            dice_transaction->
            access<const nv::index::IDistributed_data_access_factory>(data_access_tag));
        assert(access_factory.is_valid_interface());

        mi::Sint32 const volume_idx = 0; // always access to the index 0 volume
        mi::neuraylib::Tag const volume_tag = Nvindex_AppData::instance()->get_volume_tag_from_idx(volume_idx);
        const mi::math::Vector<mi::Uint32, 3> vol_size = get_volume_size(volume_tag, dice_transaction.get());

        mi::base::Handle<const nv::index::IScene> scene(
            dice_transaction->access<const nv::index::IScene>(session->get_scene()));
        assert(scene.is_valid_interface());

        const std::vector<mi::neuraylib::Tag> volume_tag_vec = Nvindex_AppData::instance()->get_volume_tag_vec();        
        assert(!volume_tag_vec.empty());

        mi::base::Handle<const nv::index::IRegular_volume> volume_element(
            dice_transaction->access<const nv::index::IRegular_volume>(volume_tag_vec.at(0))); // fix to volume 0
        assert(volume_element.is_valid_interface());

        DEBUG_LOG << "global: get_clipped_bounding_box:   "   << scene->get_clipped_bounding_box();
        DEBUG_LOG << "local:  get_XYZ_bounding_box: "         << volume_element->get_XYZ_bounding_box();
        DEBUG_LOG << "local:  get_XYZ_clipped_bounding_box: " << volume_element->get_XYZ_clipped_bounding_box();

        mi::math::Bbox<mi::Uint32, 3> vol_bbox(mi::math::Vector<mi::Uint32, 3>(0, 0,  0), vol_size);

        mi::base::Handle<nv::index::IRegular_volume_data_access> volume_data_access(
            access_factory->create_regular_volume_data_access(volume_tag));
        assert(volume_data_access.is_valid_interface());
        volume_data_access->access(vol_bbox, dice_transaction.get());

        mi::math::Bbox<mi::Uint32, 3> const effective_bbox = volume_data_access->get_bounding_box();
        INFO_LOG << "query_bbox: " << vol_bbox << ", effective bbox: " << effective_bbox;
        assert(vol_bbox == effective_bbox);

        // test the contents
//         bool const is_verified =
//             test_synthetic_voxel_value("Test synthetic volume contents.",
//                                        convert_bbox_uint32_3_to_bbox_sint64_3(effective_bbox),
//                                        convert_bbox_uint32_3_to_bbox_sint64_3(effective_bbox),
//                                        SDVE_IJ,
//                                        volume_data_access->get_amplitude_values());
//         INFO_LOG << "verified contents: " << (is_verified ? "OK." : "failed.");
//         assert(is_verified);
    }
    dice_transaction->commit();

    return true;
}

//----------------------------------------------------------------------

