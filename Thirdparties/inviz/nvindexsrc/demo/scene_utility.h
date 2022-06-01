/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief scene element access utility functions
///
/// For the test purpose, these are put in this file. Unit tests
/// access these functions also.
///
#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_SCENE_UTILITY_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_SCENE_UTILITY_H

#include <mi/dice.h>
#include <cassert>
#include <vector>
#include <string>

#include <nv/index/isession.h>
#include <nv/index/iscene_element.h>
#include <nv/index/iscene_group.h>
#include <nv/index/iindex.h>

#include "common/forwarding_logger.h"

#include "nvindex_rendering_context.h"

namespace nv {
namespace index_common {
class String_dict;
}} // namespace

class Nvindex_AppData;

//----------------------------------------------------------------------
/// Interface of scene traverse strategy.
class Scene_traverse_strategy_if
{
public:
    /// Default constructor
    Scene_traverse_strategy_if()
    {}

    /// Destructor
    virtual ~Scene_traverse_strategy_if()
    {}

    /// Controls whether disabled scene element nodes are visited.
    /// \return true if this strategy should only visit enable nodes,
    /// when false all nodes are visited.
    virtual bool visit_only_enabled_scene_elements() const
    {
        // Visit also disabled scene elements by default
        return false;
    }

    /// Controls whether disabled group nodes are visited.
    /// \return true if this strategy only visit enable group nodes,
    /// when false all groups are visited.
    virtual bool visit_only_enabled_groups() const
    {
        // Visit also disabled groups by default
        return false;
    }

    /// Called before doing recursion.
    /// \param[in] scene_element_tag current scene element tag
    /// \param[in] dice_transaction  db transaction
    virtual void run_before_recursion(
        mi::neuraylib::Tag                scene_element_tag,
        mi::neuraylib::IDice_transaction* dice_transaction) = 0;
    
    /// Called after doing recursion.
    /// \param[in] scene_element_tag current scene element tag
    /// \param[in] dice_transaction  db transaction
    virtual void run_after_recursion(
        mi::neuraylib::Tag                scene_element_tag,
        mi::neuraylib::IDice_transaction* dice_transaction)
    {
        // empty by default
    }

private:
    /// Prohibit copy constructor.
    Scene_traverse_strategy_if(const Scene_traverse_strategy_if& rhs);
    /// Prohibit operator=.
    const Scene_traverse_strategy_if& operator=(const Scene_traverse_strategy_if& rhs);
};

//----------------------------------------------------------------------
/// Example traverse strategy, print all scene element name.
class Scene_traverse_print_element : public Scene_traverse_strategy_if
{
public:
    /// default constructor
    Scene_traverse_print_element()
        :
        m_counter(0)
    {
        // empty
    }
    /// destructor
    virtual ~Scene_traverse_print_element()
    {
        // empty
    }

public:
    // Implemented Scene_traverse_strategy_if
    
    virtual void run_before_recursion(
        mi::neuraylib::Tag scene_element_tag,
        mi::neuraylib::IDice_transaction* dice_transaction)
    {
        assert(scene_element_tag.is_valid());
        assert(dice_transaction != 0);
        mi::base::Handle<const nv::index::IScene_element> scene_element(
            dice_transaction->access<const nv::index::IScene_element>(scene_element_tag));
        assert(scene_element.is_valid_interface());

        // print out this node.
        INFO_LOG << "[" << m_counter
                 << "], class name = " << scene_element->get_class_name()
                 << ", "               << (scene_element->get_enabled() ? "enabled" : "disabled");
        ++m_counter;
    }

public:
    /// get node count. To get a result after all nodes are traversed.
    mi::Sint32 get_node_count() const
    {
        return m_counter;
    }

private:
    /// a node counter
    mi::Sint32 m_counter;

private:
    /// copy constructor. prohibit until proved useful.
    Scene_traverse_print_element(const Scene_traverse_print_element & rhs);
    /// operator=. prohibit until proved useful.
    const Scene_traverse_print_element & operator=(const Scene_traverse_print_element & rhs);
};

//----------------------------------------------------------------------
/// Get the parent tag traverse
class Scene_traverse_get_parent : public Scene_traverse_strategy_if
{
public:
    /// default constructor
    /// \patam[in] child_tag a child tag who is looking for the parent.
    Scene_traverse_get_parent(const mi::neuraylib::Tag & child_tag)
        :
        m_child_tag(child_tag),
        m_parent_tag(mi::neuraylib::NULL_TAG)
    {
        // empty
        ERROR_LOG << m_child_tag << " is looking for a parent.";
    }
    /// destructor
    virtual ~Scene_traverse_get_parent()
    {
        // empty
    }

public:
    // Implemented Scene_traverse_strategy_if

    virtual void run_before_recursion(
        mi::neuraylib::Tag                scene_element_tag,
        mi::neuraylib::IDice_transaction* dice_transaction)
    {
        assert(scene_element_tag.is_valid());
        assert(dice_transaction != 0);

        mi::base::Handle<const nv::index::IScene_group> group(
            dice_transaction->access<const nv::index::IScene_group>(scene_element_tag));
        if(group.is_valid_interface()){
            // this is group, will check it out
            mi::Uint32 nb_childern = group->nb_elements();
            for(mi::Uint32 i = 0; i < nb_childern; ++i){
                const mi::neuraylib::Tag ch_tag = group->get_scene_element(i);
                if(m_child_tag == ch_tag){
                    // scene_element_tag is m_child_tag's parent
                    if(m_parent_tag.is_valid()){
                        WARN_LOG << "the child (" << m_child_tag << ") already found parent ("
                                 << m_parent_tag << "), overwritten by " << scene_element_tag;
                    }
                    m_parent_tag = scene_element_tag;
                    INFO_LOG << "Found parent (" << m_parent_tag << ") of ("
                             << m_child_tag << ")." ;
                    break;
                }
            }
        }
    }

public:
    /// get child tag
    mi::neuraylib::Tag get_child_tag() const
    {
        return m_child_tag;
    }
    /// get parent tag
    mi::neuraylib::Tag get_parent_tag() const
    {
        return m_parent_tag;
    }

private:
    /// child tag who is looking for the parent.
    mi::neuraylib::Tag m_child_tag;
    /// resulting parent
    mi::neuraylib::Tag m_parent_tag;

private:
    /// copy constructor. prohibit until proved useful.
    Scene_traverse_get_parent(const Scene_traverse_get_parent & rhs);
    /// operator=. prohibit until proved useful.
    const Scene_traverse_get_parent & operator=(const Scene_traverse_get_parent & rhs);
};


//----------------------------------------------------------------------
/// volume and heightfield data collection strategy
class Scene_traverse_volume_heightfield_element : public Scene_traverse_strategy_if
{
public:
    /// default constructor
    Scene_traverse_volume_heightfield_element();

    /// destructor
    virtual ~Scene_traverse_volume_heightfield_element();

    /// clear
    void initialize();

    /// get enabled number of volume scene element
    mi::Sint32 get_nb_volumes() const;
    /// get enabled number of sparse volume scene element
    mi::Sint32 get_nb_sparse_volumes() const;
    /// get enabled number of heightfield scene element
    mi::Sint32 get_nb_heightfield() const;
    /// get enabled number of rendering kernel parameter buffer scene element
    mi::Sint32 get_nb_rtc_parameter_buffer() const;

    /// get volume scene element tag vector
    std::vector<mi::neuraylib::Tag> get_volume_tag_vec() const;
    /// get sparse volume scene element tag vector
    std::vector<mi::neuraylib::Tag> get_sparse_volume_tag_vec() const;
    /// get heightfield scene element tag vector
    std::vector<mi::neuraylib::Tag> get_heightfield_tag_vec() const;
    /// get rendering kernel parameter buffer scene element tag vector
    std::vector<mi::neuraylib::Tag> get_rtc_param_buffer_tag_vec() const;

public:
    // Implemented Scene_traverse_strategy_if

    virtual bool visit_only_enabled_groups() const
    {
        // disabled groups are not visited
        return true;
    }

    virtual void run_before_recursion(
        mi::neuraylib::Tag                scene_element_tag,
        mi::neuraylib::IDice_transaction* dice_transaction);

private:
    /// volume scene element tag vector
    std::vector<mi::neuraylib::Tag> m_volume_tag_vec;
    /// sparse volume scene element tag vector
    std::vector<mi::neuraylib::Tag> m_sparse_volume_tag_vec;
    /// heightfield scene element tag vector
    std::vector<mi::neuraylib::Tag> m_heightfield_tag_vec;
    /// rendering kernel parameter buffer scene element tag vector
    std::vector<mi::neuraylib::Tag> m_rtc_param_buffer_tag_vec;

private:
    /// copy constructor. prohibit until proved useful.
    Scene_traverse_volume_heightfield_element(const Scene_traverse_volume_heightfield_element & rhs);
    /// operator=. prohibit until proved useful.
    const Scene_traverse_volume_heightfield_element & operator=(const Scene_traverse_volume_heightfield_element & rhs);
};

//----------------------------------------------------------------------
/// tag data collection strategy
class Scene_traverse_collect_tag : public Scene_traverse_strategy_if
{
public:
    /// default constructor
    Scene_traverse_collect_tag();

    /// destructor
    virtual ~Scene_traverse_collect_tag();

    /// clear
    void initialize();

    /// get number of tags in the scene
    mi::Size get_nb_tags() const;
    /// get tag vector
    std::vector<mi::neuraylib::Tag> get_tag_vec() const;

public:
    // Implemented Scene_traverse_strategy_if

    virtual bool visit_only_enabled_groups() const
    {
        // disabled groups are not visited
        return true;
    }

    virtual void run_before_recursion(
        mi::neuraylib::Tag                scene_element_tag,
        mi::neuraylib::IDice_transaction* dice_transaction);

private:
    /// all tags in the scene vector
    std::vector<mi::neuraylib::Tag> m_tag_vec;

private:
    /// copy constructor. prohibit until proved useful.
    Scene_traverse_collect_tag(const Scene_traverse_collect_tag & rhs);
    /// operator=. prohibit until proved useful.
    const Scene_traverse_collect_tag & operator=(const Scene_traverse_collect_tag & rhs);
};

//----------------------------------------------------------------------
/// Get the bounding box by traverse
class Scene_traverse_scene_bounding_box : public Scene_traverse_strategy_if
{
public:
    /// default constructor
    Scene_traverse_scene_bounding_box();

    /// destructor
    virtual ~Scene_traverse_scene_bounding_box();

    /// clear
    void initialize();

    /// Get the bounding box
    ///
    /// \return bounding box of the traverse result
    mi::math::Bbox<mi::Float32, 3> get_bounding_box() const;

public:
    // Implemented Scene_traverse_strategy_if

    virtual bool visit_only_enabled_groups() const
    {
        // disabled groups are not visited
        return true;
    }

    virtual void run_before_recursion(
        mi::neuraylib::Tag                scene_element_tag,
        mi::neuraylib::IDice_transaction* dice_transaction);

    virtual void run_after_recursion(
        mi::neuraylib::Tag                scene_element_tag,
        mi::neuraylib::IDice_transaction* dice_transaction);

private:
    /// Scene bounding box
    mi::math::Bbox<mi::Float32, 3>  m_scene_bbox;

    /// transformation matrix stack (multiplied leaves)
    std::vector< mi::math::Matrix<mi::Float32, 4, 4> > m_transform_mat_stack;

private:
    /// copy constructor. prohibit until proved useful.
    Scene_traverse_scene_bounding_box(const Scene_traverse_scene_bounding_box & rhs);
    /// operator=. prohibit until proved useful.
    const Scene_traverse_scene_bounding_box & operator=(const Scene_traverse_scene_bounding_box & rhs);
};

//----------------------------------------------------------------------
/// Get volume size. A convenient function.
///
/// \param[in]  volume_tag  a volume tag. This should be a valid volume tag.
/// \param[in]  dice_transaction dice transaction
/// \return size of the volume
mi::math::Vector<mi::Uint32, 3> get_volume_size(
    const mi::neuraylib::Tag &        volume_tag,
    mi::neuraylib::IDice_transaction* dice_transaction);

//----------------------------------------------------------------------
/// Get the current volume region of interest bounding box as Uint32
/// bbox. A convenient function.
///
/// \param[in]  volume_tag  a volume tag. This should be a valid volume tag.
/// \param[in]  dice_transaction dice transaction
/// \return region of interest bounding box in Uint32 bbox.
mi::math::Bbox<mi::Sint32, 3> get_volume_roi_bbox(
    const mi::neuraylib::Tag&         volume_tag,
    mi::neuraylib::IDice_transaction* dice_transaction);

//----------------------------------------------------------------------
/// Get the current heightfield region of interest bounding box as Sint32
/// bbox. A convenient function.
///
/// \param[in]  heightfield_tag  a heightfield tag. This should be a heightfield tag.
/// \param[in]  dice_transaction dice transaction
/// \return region of interest bounding box in Sint32 bbox.
mi::math::Bbox<mi::Sint32, 3> get_heightfield_roi_bbox(
    const mi::neuraylib::Tag&         heightfield_tag,
    mi::neuraylib::IDice_transaction* dice_transaction);


//----------------------------------------------------------------------
/// Get the whole scene bounding box by traverse the scene.
mi::math::Bbox<mi::Float32, 3> get_scene_bounding_box(
    mi::neuraylib::Tag                session_tag,
    mi::neuraylib::IDice_transaction* dice_transaction);

//----------------------------------------------------------------------
/// set examiner window resolution. Only change the examiner's
/// resolution, no camera resolution change.
/// \deprecated
/// \param[in] win_res window resolution (width, height)
void set_examiner_window_resolution(mi::Sint32_2 const & win_res);

//----------------------------------------------------------------------
/// request canvas resolution. You can call this from any threads.
///
/// \param[in] new_canvas_res requesting new canvas resolution
void request_canvas_resolution(
    mi::Sint32_2 const & new_canvas_res);

//----------------------------------------------------------------------
/// update window resolution to the requested resolution.
/// Resolutions for both the examiner and the main camera.
/// The resolution must be requested before.
///
/// \param[in] irc               IndeX rendering context
void update_canvas_resolution(
    Nvindex_rendering_context& irc);

//----------------------------------------------------------------------
/// change image canvas resolution. You can call this from any threads.
///
/// \param[in] new_image_res new image resolution
void set_image_canvas_resolution(
    const mi::Sint32_2& new_image_res);

//----------------------------------------------------------------------
/// get current image canvas resolution. You can call this from any threads.
///
/// \return get current image resolution
mi::Sint32_2 get_image_canvas_resolution();

//----------------------------------------------------------------------
/// get snapshot filename by mode (sequence index or frame number)
///
/// \param[in] appdata nvindex application data
/// \param[in] frame_num       frame number (given by the application)
/// \return snapshot filename
std::string get_snapshot_filename_by_mode(
    Nvindex_AppData* appdata,
    mi::Uint32       frame_num);

//----------------------------------------------------------------------
/// copy span buffer to canvas (e.g., encoding cavnas)
/// If necessary canvas will be resized.
///
/// \param[in] span_buffer span buffer pointer
/// \param[in] canvas      a canvas (e.g., encoding canvas)
void copy_result_pixel_to_canvas(Span_renderer_IF*         span_buffer,
                                 nv::index_common::Canvas* canvas);

//----------------------------------------------------------------------
/// get cluster host id they contain the volume with query bbox
///
/// \param[in] session_tag session tag
/// \param[in] volume_tag  accessing volume tag
/// \param[in] query_bbox  query bounding box of volume in IJK space
/// \param[out] cluster_host_id_vec (output) hosts who have volume in bounds
/// \param[in] dice_transaction dice db transaction
void get_cluster_host_containing_volume(
    const mi::neuraylib::Tag&            session_tag,
    const mi::neuraylib::Tag&            volume_tag,
    const mi::math::Bbox<mi::Uint32, 3>& query_bbox,
    std::vector<mi::Uint32>&             cluster_host_id_vec,
    mi::neuraylib::IDice_transaction*    dice_transaction);

//----------------------------------------------------------------------
/// get the global ROI and volume local ROI
///
/// \param[in] session session object
/// \param[in] dice_transaction transaction of DB
/// \return global region of interest bounding box
mi::math::Bbox<mi::Float32, 3> get_XYZ_global_region_of_interest_bbox(
    const nv::index::ISession*        session,
    mi::neuraylib::IDice_transaction* dice_transaction);

//----------------------------------------------------------------------
/// set XYZ global region of interest
///
/// \param[in] xyz_bbox xyz global region of interest bounding box
/// \param[in] session session object
/// \param[in] dice_transaction transaction of DB
void set_XYZ_global_region_of_interest_bbox(
    const mi::math::Bbox<mi::Float32, 3>& xyz_bbox,
    const nv::index::ISession*            session,
    mi::neuraylib::IDice_transaction*     dice_transaction);

//----------------------------------------------------------------------
/// get volume IJK ROI bounding box
///
/// \param[in] dice_transaction dice transaction
/// \param[in] volume_tag
/// \return volume local IJK bounding box
mi::math::Bbox< mi::Sint64, 3 > get_volume_local_IJK_ROI_bbox(
    mi::neuraylib::Tag const&         volume_tag,
    mi::neuraylib::IDice_transaction* dice_transaction);

//----------------------------------------------------------------------
/// set volume IJK ROI bounding box
///
/// \param[in] volume_tag volume tag
/// \param[in] volume_ijk_roi volume's IJK ROI
/// \param[in] dice_transaction dice transaction
void set_volume_local_IJK_ROI_bbox(
    const mi::neuraylib::Tag&              volume_tag,
    const mi::math::Bbox< mi::Sint64, 3 >& volume_ijk_roi,
    mi::neuraylib::IDice_transaction*      dice_transaction);

//----------------------------------------------------------------------
/// get heightfield IJK ROI bounding box
///
/// \param[in] dice_transaction dice transaction
/// \param[in] heightfield_tag heightfield tag
/// \return heightfield local IJK bounding box
mi::math::Bbox< mi::Sint64, 3 > get_heightfield_local_IJK_ROI_bbox(
    mi::neuraylib::IDice_transaction* dice_transaction,
    const mi::neuraylib::Tag&         heightfield_tag);

//----------------------------------------------------------------------
/// set heightfield IJK ROI bounding box
///
/// \param[in] heightfield_tag heightfield tag
/// \param[in] heightfield_ijk_roi heightfield's IJK ROI
/// \param[in] dice_transaction dice transaction
void set_heightfield_local_IJK_ROI_bbox(
    const mi::neuraylib::Tag& heightfield_tag,
    mi::math::Bbox< mi::Sint64, 3 > const & heightfield_ijk_roi,
    mi::neuraylib::IDice_transaction * dice_transaction);

//----------------------------------------------------------------------

/// Scene element type enum
enum Scene_element_type_e {
    /// volume
    SE_volume,
    /// heightfield
    SE_heightfield,
    /// point set
    SE_point,
    /// line set
    SE_line,
    /// plane 
    SE_plane,
    /// None
    SE_none
};

/// get scene element type (limited)
/// 
/// \param[in] scene_element_tag scene element tag.
/// \param[out] result_index when found the element, returns the
/// index. -1 when not found.
/// \param[in] dice_transaction dice transaction
/// \return type 
Scene_element_type_e get_scene_element_type(
    const mi::neuraylib::Tag&         scene_element_tag,
    mi::Sint32&                       result_index,
    mi::neuraylib::IDice_transaction* dice_transaction);

//----------------------------------------------------------------------
/// set slice change
///
/// \param[in] scope the scope
/// \param[in] current_volume_index current volume index
/// \param[in] slice_id slice id {INLINE_SECTION, CROSS_LINE_SECTION,
/// HORIZONTAL_SECTION, VERTICAL_PROFILE}
/// \param[in] slice_position slice position (what coordinate?)
/// \param[in] slice_colormap_tag slice colormap
/// \param[in] is_slice_enabled visualize slice when true
void set_slice_change(
    mi::neuraylib::IScope*    scope,
    mi::Uint32                current_volume_index,
    mi::Uint32                slice_id,
    mi::Float32               slice_position,
    const mi::neuraylib::Tag& slice_colormap_tag,
    bool                      is_slice_enabled);

//----------------------------------------------------------------------
/// get new-ed heightfield patch distribution layout for whole heightfield (data_locality)
///
/// FIXME: need test for valuerange/bbox update
///
/// \param[in]  session_tag               session tag
/// \param[in]  heightfield_tag           heightfield scene element tag
/// \param[in]  is_edit                   true for editing computation, false otherwise
/// \param[in]  dice_transaction          dice transaction
/// \param[out] heightfield_bbox          (output) heightfield bounding box
/// \return data locality
nv::index::IRegular_heightfield_data_locality* new_whole_heightfield_distribution_layout(
    const mi::neuraylib::Tag&             session_tag,
    const mi::neuraylib::Tag&             heightfield_tag,
    bool                                  is_edit,
    mi::neuraylib::IDice_transaction *    dice_transaction,
    mi::math::Bbox_struct<mi::Uint32, 3>& heightfield_bbox);

//----------------------------------------------------------------------
/// get new-ed heightfield patch distribution layout of ij_query_bbox (data_locality)
///
/// \param[in]  session_tag               session tag
/// \param[in]  heightfield_tag           heightfield_scene element tag
/// \param[in]  ij_query_bbox          ij query bounding box of heightfield data
/// \param[in]  is_edit                   true for editing computation, false otherwise
/// \param[in]  dice_transaction          dice transaction
/// \return data locality
nv::index::IRegular_heightfield_data_locality* new_heightfield_distribution_layout(
    const mi::neuraylib::Tag&            session_tag,
    const mi::neuraylib::Tag&            heightfield_tag,
    mi::math::Bbox_struct<mi::Uint32, 2> const & ij_query_bbox,
    bool                                  is_edit,
    mi::neuraylib::IDice_transaction *    dice_transaction);

//----------------------------------------------------------------------
/// get new-ed volume distribution layout (data_locality)
///
/// \param[in]  session_tag      session tag
/// \param[in]  volume_tag       volume tag
/// \param[in]  ijk_query_bbox   ijk query bounding box of volume data
/// \param[in]  dice_transaction dice transaction
/// \return new-ed data locality. recommended to receive it with Handle.
nv::index::IRegular_volume_data_locality* new_volume_distribution_layout(
    const mi::neuraylib::Tag & session_tag,
    const mi::neuraylib::Tag & volume_tag,
    const mi::math::Bbox_struct<mi::Uint32, 3> & ijk_query_bbox,
    mi::neuraylib::IDice_transaction * dice_transaction);

//----------------------------------------------------------------------
/// scene hierarchy traverse
///
/// \param[in] scene_group_tag traverse root group tag
/// \param[in] p_strategy      scene traverse strategy
/// \param[in] dice_transaction database transaction
/// \return true when success
bool traverse_scene_hierarchy_by_tag(
    const mi::neuraylib::Tag           scene_group_tag,
    Scene_traverse_strategy_if * const p_strategy,
    mi::neuraylib::IDice_transaction*  dice_transaction);

//----------------------------------------------------------------------
/// Wrapper to export the session to a file or stdout, adding some application-specifiy settings
void export_session(
    const Nvindex_rendering_context&  irc,
    mi::neuraylib::IDice_transaction* dice_transaction,
    const std::string&                output_filename = "");

//----------------------------------------------------------------------
/// Get scene status by traversing the scene. Traverse implementation.
///
/// \param[in] p_scene_traverse_strategy traverse strategy
/// \param[in] session_tag               session tag
/// \param[in] dice_transaction          dice db transaction
/// \return true when success
bool get_scene_status(
    Scene_traverse_strategy_if* p_scene_traverse_strategy,
    mi::neuraylib::Tag          session_tag,
    mi::neuraylib::IDice_transaction* dice_transaction);

//----------------------------------------------------------------------
/// Get center of volume inside of region of interest
///
/// \param[in]  irc_ref    index rendering context reference
/// \param[in]  volume_tag volume tag to considering
/// \param[out] volume_center volume center in world space
/// \return true when succeeded
bool get_center_of_volume_in_roi(
    Nvindex_rendering_context*        irc_ref,
    const mi::neuraylib::Tag&         volume_tag,
    mi::math::Vector<mi::Float32, 3>& volume_center);
//----------------------------------------------------------------------

#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_SCENE_UTILITY_H
