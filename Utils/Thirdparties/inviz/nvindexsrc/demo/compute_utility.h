/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief computation task related utility functions
#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_COMPUTE_UTILITY_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_COMPUTE_UTILITY_H

#include <mi/dice.h>
#include <nv/index/isession.h>
#include <nv/index/iindex.h>

#include <vector>
#include <string>

namespace nv {
namespace index_common {
class String_dict;
}} // namespace

//----------------------------------------------------------------------
/// get unique dice name from the base_dice_name.
///
/// Note: Using linear search. Fine a unique base_name + number.
///
/// \param[in] base_dice_name base dice name 
/// \param[in] dice_transaction dice transaction
/// \return return the unique dice name
std::string get_unique_dice_name(const std::string & base_dice_name, 
                                 mi::neuraylib::IDice_transaction * dice_transaction);

//----------------------------------------------------------------------

/// setup parallel fragment jobs
///
/// \param[in] cluster_hosts cluster host list. e.g., (1 2 3 4)
/// \param[in] parallel_count number of parallel fragment execution on
/// a host.
/// \param[out] fragment_to_host_id (output) fragment to host id map
/// (1 1 2 2 3 3 4 4) for two concurrence
/// \param[out] fragment_to_host_local_thread_id fragment fragment to
/// thread id map. ((0 1) (0 1) (0 1) (0 1)) means
/// ((host_1's thread 0 host_1's thread 1) (host_2's thread 0 host_2's thread 1)
///  (host_3's thread 0 host_3's thread 1) (host_4's thread 0 host_4's thread 1))
void setup_parallel_fragment(std::vector< mi::Uint32 > const & cluster_hosts,
                             mi::Uint32 const parallel_count,
                             std::vector< mi::Sint32 > & fragment_to_host_id,
                             std::vector< mi::Sint32 > & fragment_to_host_local_thread_id);

//----------------------------------------------------------------------
/// heightfield workflow manual picking
///
/// \param[in] scope the scope
/// \param[in] pick_location pick location on the screen
/// \param[in] pick_canvas   pick canvas
void heightfield_workflow_manual_picking(
    mi::base::Handle< mi::neuraylib::IScope > & scope,
    mi::math::Vector_struct< mi::Uint32, 2 > const & pick_location,
    const nv::index::IIndex_canvas*                  pick_canvas
    );

//----------------------------------------------------------------------
/// heightfield workflow appending a vertex to bounding polygon.
///
/// \param[in] scope the scope
/// \param[in] iindex_scene_query scene query
/// \param[in] session_tag session tag
/// \param[in] pick_location pick location on the screen
/// \param[in] pick_canvas   pick canvas
void heightfield_workflow_append_vertex_to_bounding_polygon(
    mi::base::Handle<mi::neuraylib::IScope>&          scope,
    mi::base::Handle<nv::index::IIndex_scene_query>&  iindex_scene_query,
    const mi::neuraylib::Tag&                         session_tag,
    const mi::math::Vector_struct< mi::Uint32, 2 >&   pick_location,
    const nv::index::IIndex_canvas*                   pick_canvas
    );

//----------------------------------------------------------------------
/// hightfield workflow polygon deletion
///
/// \param[in] scope the scope
/// \param[in] use_convex_hull if the convex hull or the polygon loop shall be used.
///
void heightfield_workflow_delete_polygon(
    mi::base::Handle< mi::neuraylib::IScope >&  scope,
    bool                                        use_convex_hull);

//----------------------------------------------------------------------
/// Heightfield workflow elevation change
///
/// \param[in] scope               the scope
/// \param[in] new_elevation_value new elevation value (or scaling
///            factor depends on is_scaling_op)
/// \param[in] is_scaling_op scaling operation when true instead of
///            set elevation value
/// \param[in] use_convex_hull if the convex hull or the polygon loop shall be used.
void heightfield_workflow_change_elevation_values(
    mi::base::Handle<mi::neuraylib::IScope>& scope,
    mi::Float32                              new_elevation_value,
    bool                                     is_scaling_op,
    bool                                     use_convex_hull);

//----------------------------------------------------------------------
/// Heightfield workflow gridding
///
/// \param[in] scope the scope
/// \param[in] use_convex_hull if the convex hull or the polygon loop shall be used.
///
void heightfield_workflow_gridding(
    mi::base::Handle<mi::neuraylib::IScope>&    scope,
    bool                                        use_convex_hull);

//----------------------------------------------------------------------
/// Heightfield autotracking workflow
///
/// \param[in] scope the scope
///
void heightfield_workflow_autotracking(
    mi::base::Handle<mi::neuraylib::IScope>& scope);

//----------------------------------------------------------------------
/// heightfield amplitude value map generation example
///
/// FIXME: Need test due to valuerange/bbox IF update
///
/// \param[in] com_str command string. expected "test_template_heightfield param".
/// \param[in] scope the db scope
/// \param[in] iindex_session IndeX session
/// \param[in] session_tag session tag
/// \return true when test succeeded
bool generate_heightfield_amplitude_map(
    const std::string&                                     com_str,
    mi::base::Handle<mi::neuraylib::IScope> &              scope,
    mi::base::Handle<nv::index::IIndex_session>&           iindex_session,
    const mi::neuraylib::Tag&                              session_tag);

//----------------------------------------------------------------------

/// Assign the given amplitude value to all voxels in the
/// region-of-interest of the given volume scene element.
///
/// \param[in] is_host_assign_mode show host assignment instead of amplitude value
/// \param[in] amplitude_value  voxel amplitude value to be set when is_host_assign_mode false
/// \param[in] session_tag      session tag to access to the data locality
/// \param[in] volume_tag       volume tag of the volume to be filled
/// \param[in] dice_transaction dice db transaction
/// \return true when success
bool volume_set_amplitude_value(
    bool                      is_host_assign_mode,
    mi::Uint8                 amplitude_value,
    const mi::neuraylib::Tag& session_tag,
    const mi::neuraylib::Tag& volume_tag,
    mi::neuraylib::IDice_transaction * dice_transaction);

//----------------------------------------------------------------------
/// apply volume inset filter of amplitude values
///
/// \param[in] session_tag      session tag
/// \param[in] volume_tag       volume tag to be filtered
/// \param[in] filter_type_str  filter type name
/// \param[in] parallel_count   
/// \param[in] dice_transaction dice db transaction
/// \return true when success
bool volume_inset_filter_amplitude_values(
    const mi::neuraylib::Tag & session_tag,
    const mi::neuraylib::Tag & volume_tag,
    const std::string & filter_type_str,
    const mi::Sint32 parallel_count, // FIXME: may not needed. Use slice number (or is_parallel_execution boolean)
    mi::neuraylib::IDice_transaction * dice_transaction
    );

//----------------------------------------------------------------------
/// generate volume data by copy.
///
/// \param[in] session_tag      session tag
/// \param[in] src_volume_tag   source volume tag
/// \param[in] dst_volume_name  generated destination volume name. If
///                             the name is not unique in dice DB, the
///                             generated volume name may be altered.
/// \param[in] tr_vec           translation vector for the generated volume
/// \param[in] dice_transaction dice translation
/// \return the tag of the generated volume
mi::neuraylib::Tag generate_volume_by_copy(
    const mi::neuraylib::Tag & session_tag,
    const mi::neuraylib::Tag & src_volume_tag,
    const std::string & dst_volume_name,
    const mi::math::Vector<mi::Float32, 3> & tr_vec,
    mi::neuraylib::IDice_transaction * dice_transaction
    );

//----------------------------------------------------------------------
/// create default static group for a volume and append it to the scene root
///
/// \param[in] session_tag       session tag
/// \param[in] adding_volume_tag a volume tag adding to the static group node
/// \param[in] dice_transaction dice translation
/// \return the tag of the generated static group node
mi::neuraylib::Tag add_volume_to_the_scene_root(
    const mi::neuraylib::Tag & session_tag,
    const mi::neuraylib::Tag & adding_volume_tag,
    mi::neuraylib::IDice_transaction * dice_transaction);

//----------------------------------------------------------------------
/// verify the synthetic volume contents
///
/// \param[in] param test parameter
/// \param[in] scope the db scope
/// \param[in] session_tag session tag
/// \return true when test succeeded
bool test_verify_synthetic_volume_contents(
    std::string const & param,
    mi::base::Handle<mi::neuraylib::IScope> & scope,
    const mi::neuraylib::Tag&                session_tag);

//----------------------------------------------------------------------
#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_COMPUTE_UTILITY_H
