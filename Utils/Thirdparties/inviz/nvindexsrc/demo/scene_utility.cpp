/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief scene element access utility functions
///
/// For the test purpose, these are put in this file. Unit tests
/// access these functions also.
///

#include "scene_utility.h"

#include <nv/index/idistributed_data_locality.h>
#include <nv/index/iplane.h>
#include <nv/index/ipoint_set.h>
#include <nv/index/irendering_kernel_programs.h>
#include <nv/index/iregular_heightfield.h>
#include <nv/index/iregular_volume.h>
#include <nv/index/iscene.h>
#include <nv/index/isparse_volume_scene_element.h>
#include <nv/index/iviewport.h>

#include "common/colormap_io.h"
#include "common/common_utility.h"
#include "common/string_dict.h"
#include "common/type_conversion_utility.h"

#include "examiner_manipulator.h"
#include "nvindex_appdata.h"
#include "opengl_appdata.h"
#include "utilities.h"

//----------------------------------------------------------------------
Scene_traverse_volume_heightfield_element::Scene_traverse_volume_heightfield_element()
    :
    m_volume_tag_vec(),
    m_sparse_volume_tag_vec(),
    m_heightfield_tag_vec()
{
    // empty
}

//----------------------------------------------------------------------
Scene_traverse_volume_heightfield_element::~Scene_traverse_volume_heightfield_element()
{
    // empty
}

//----------------------------------------------------------------------
void Scene_traverse_volume_heightfield_element::initialize()
{
    m_volume_tag_vec.clear();
    m_heightfield_tag_vec.clear();
}

//----------------------------------------------------------------------
std::vector<mi::neuraylib::Tag> Scene_traverse_volume_heightfield_element::get_volume_tag_vec() const
{
    return m_volume_tag_vec;
}

//----------------------------------------------------------------------
std::vector<mi::neuraylib::Tag> Scene_traverse_volume_heightfield_element::get_sparse_volume_tag_vec() const
{
    return m_sparse_volume_tag_vec;
}

//----------------------------------------------------------------------
std::vector<mi::neuraylib::Tag> Scene_traverse_volume_heightfield_element::get_heightfield_tag_vec() const
{
    return m_heightfield_tag_vec;
}

std::vector<mi::neuraylib::Tag> Scene_traverse_volume_heightfield_element::get_rtc_param_buffer_tag_vec() const
{
    return m_rtc_param_buffer_tag_vec;
}

//----------------------------------------------------------------------
mi::Sint32 Scene_traverse_volume_heightfield_element::get_nb_volumes() const
{
    return static_cast<mi::Sint32>(m_volume_tag_vec.size());
}

//----------------------------------------------------------------------
mi::Sint32 Scene_traverse_volume_heightfield_element::get_nb_sparse_volumes() const
{
    return static_cast<mi::Sint32>(m_sparse_volume_tag_vec.size());
}

//----------------------------------------------------------------------
mi::Sint32 Scene_traverse_volume_heightfield_element::get_nb_heightfield() const
{
    return static_cast<mi::Sint32>(m_heightfield_tag_vec.size());
}

//----------------------------------------------------------------------
mi::Sint32 Scene_traverse_volume_heightfield_element::get_nb_rtc_parameter_buffer() const
{
    return static_cast<mi::Sint32>(m_rtc_param_buffer_tag_vec.size());
}

//----------------------------------------------------------------------
void Scene_traverse_volume_heightfield_element::run_before_recursion(
    mi::neuraylib::Tag                scene_element_tag,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    assert(scene_element_tag.is_valid());
    assert(dice_transaction != 0);

    mi::base::Handle<const nv::index::IScene_element> scene_element(
        dice_transaction->access<const nv::index::IScene_element>(scene_element_tag));
    assert(scene_element.is_valid_interface());

    mi::base::Handle<const nv::index::IRegular_volume> volume_element(
        scene_element->get_interface<const nv::index::IRegular_volume>());
    mi::base::Handle<const nv::index::IRegular_heightfield> heightfield_element(
        scene_element->get_interface<const nv::index::IRegular_heightfield>());
    mi::base::Handle<const nv::index::ISparse_volume_scene_element> svol_element(
        scene_element->get_interface<const nv::index::ISparse_volume_scene_element>());
    mi::base::Handle<const nv::index::IRendering_kernel_program_parameters> rtc_param_element(
        scene_element->get_interface<const nv::index::IRendering_kernel_program_parameters>());
    if (volume_element.is_valid_interface())
    {
        m_volume_tag_vec.push_back(scene_element_tag); // found volume
    }
    else if (svol_element.is_valid_interface())
    {
        m_sparse_volume_tag_vec.push_back(scene_element_tag); // found heightfield
    }
    else if (heightfield_element.is_valid_interface())
    {
        m_heightfield_tag_vec.push_back(scene_element_tag); // found heightfield
    }
    else if (rtc_param_element.is_valid_interface())
    {
        m_rtc_param_buffer_tag_vec.push_back(scene_element_tag); // found heightfield
    }
}

//----------------------------------------------------------------------
//======================================================================
//----------------------------------------------------------------------
Scene_traverse_collect_tag::Scene_traverse_collect_tag()
{
    // empty
}

//----------------------------------------------------------------------
Scene_traverse_collect_tag::~Scene_traverse_collect_tag()
{
    // empty
}

//----------------------------------------------------------------------
void Scene_traverse_collect_tag::initialize()
{
    m_tag_vec.clear();
}

//----------------------------------------------------------------------
mi::Size Scene_traverse_collect_tag::get_nb_tags() const
{
    return m_tag_vec.size();
}

//----------------------------------------------------------------------
std::vector<mi::neuraylib::Tag> Scene_traverse_collect_tag::get_tag_vec() const
{
    return m_tag_vec;
}

//----------------------------------------------------------------------
void Scene_traverse_collect_tag::run_before_recursion(
    mi::neuraylib::Tag                scene_element_tag,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    assert(scene_element_tag.is_valid());
    assert(dice_transaction != 0);

    m_tag_vec.push_back(scene_element_tag);
}

//----------------------------------------------------------------------
//======================================================================
//----------------------------------------------------------------------
Scene_traverse_scene_bounding_box::Scene_traverse_scene_bounding_box()
{
    initialize();
}

//----------------------------------------------------------------------
Scene_traverse_scene_bounding_box::~Scene_traverse_scene_bounding_box()
{
    // empty
}

//----------------------------------------------------------------------
void Scene_traverse_scene_bounding_box::initialize()
{
    m_scene_bbox.clear();

    m_transform_mat_stack.clear();
    mi::math::Matrix<mi::Float32, 4, 4> eye(1.0f);
    m_transform_mat_stack.push_back(eye);
}

//----------------------------------------------------------------------
mi::math::Bbox<mi::Float32, 3> Scene_traverse_scene_bounding_box::get_bounding_box() const
{
    return m_scene_bbox;
}

//----------------------------------------------------------------------
void Scene_traverse_scene_bounding_box::run_before_recursion(
    mi::neuraylib::Tag                scene_element_tag,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    assert(scene_element_tag.is_valid());
    assert(dice_transaction != 0);

    mi::base::Handle<const nv::index::IScene_element> scene_element(
        dice_transaction->access<nv::index::IScene_element>(scene_element_tag));
    if (!(scene_element_tag.is_valid()))
    {
        ERROR_LOG << "Failed to access to the scene element.";
        return;
    }

    mi::base::Handle<const nv::index::ITransformed_scene_group> transformed_scene_group(
        scene_element->get_interface<const nv::index::ITransformed_scene_group>());
    mi::base::Handle<const nv::index::IDistributed_data> distributed_data(
        scene_element->get_interface<const nv::index::IDistributed_data>());
    mi::base::Handle<const nv::index::IObject_space_shape> object_space_shape(
        scene_element->get_interface<const nv::index::IObject_space_shape>());

    assert(!m_transform_mat_stack.empty());
    if (transformed_scene_group.is_valid_interface())
    {
        mi::math::Matrix<mi::Float32, 4, 4> cur_mat = transformed_scene_group->get_transform();
        mi::math::Matrix<mi::Float32, 4, 4> stacked_mat = m_transform_mat_stack.back() * cur_mat;
        m_transform_mat_stack.push_back(stacked_mat);        
        // INFO_LOG << "scene_bbox: depth: " << m_transform_mat_stack.size() << ", cur_mat: " << cur_mat;
    }
    else if (distributed_data.is_valid_interface())
    {
        const mi::math::Bbox<mi::Float32, 3> cur_bbox = distributed_data->get_bounding_box();
        const mi::math::Bbox<mi::Float32, 3> tr_bbox = mi::math::transform_point(m_transform_mat_stack.back(), cur_bbox);
        // const mi::math::Bbox<mi::Float32, 3> debug_last = m_scene_bbox;
        m_scene_bbox.insert(tr_bbox);
        // INFO_LOG << "scene_bbox: distributed bbox: " << debug_last << " -> " << m_scene_bbox;
    }
    else if (object_space_shape.is_valid_interface())
    {
        const mi::math::Bbox<mi::Float32, 3> cur_bbox = object_space_shape->get_bounding_box();
        const mi::math::Bbox<mi::Float32, 3> tr_bbox = mi::math::transform_point(m_transform_mat_stack.back(), cur_bbox);
        // const mi::math::Bbox<mi::Float32, 3> debug_last = m_scene_bbox;
        m_scene_bbox.insert(tr_bbox);
        // INFO_LOG << "scene_bbox: shape bbox: " << debug_last << " -> " << m_scene_bbox;
    }
    else
    {
        // no transformation and not a 3D object
    }
}

//----------------------------------------------------------------------
void Scene_traverse_scene_bounding_box::run_after_recursion(
    mi::neuraylib::Tag                scene_element_tag,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    assert(scene_element_tag.is_valid());
    assert(dice_transaction != 0);

    mi::base::Handle<const nv::index::IScene_element> scene_element(
        dice_transaction->access<nv::index::IScene_element>(scene_element_tag));
    if (!(scene_element_tag.is_valid()))
    {
        ERROR_LOG << "Failed to access to the scene element.";
        return;
    }

    mi::base::Handle<const nv::index::ITransformed_scene_group> transformed_scene_group(
        scene_element->get_interface<const nv::index::ITransformed_scene_group>());

    assert(!m_transform_mat_stack.empty());
    if (transformed_scene_group.is_valid_interface())
    {
        m_transform_mat_stack.pop_back();        
    }

    assert(!m_transform_mat_stack.empty());
}

//----------------------------------------------------------------------
//======================================================================
//----------------------------------------------------------------------
mi::math::Vector<mi::Uint32, 3> get_volume_size(
    const mi::neuraylib::Tag&         volume_tag,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    assert(dice_transaction != 0);
    assert(volume_tag.is_valid());

    mi::base::Handle<const nv::index::IRegular_volume> volume_element(
        dice_transaction->access<const nv::index::IRegular_volume>(volume_tag));
    assert(volume_element.is_valid_interface());

    const mi::math::Vector<mi::Uint32, 3> ret_size = volume_element->get_volume_size();

    return ret_size;
}

//----------------------------------------------------------------------
mi::math::Bbox<mi::Sint32, 3> get_volume_roi_bbox(
    const mi::neuraylib::Tag&         volume_tag,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    assert(volume_tag.is_valid());
    assert(dice_transaction != 0);

    mi::base::Handle<const nv::index::IRegular_volume> volume_element(
        dice_transaction->access<const nv::index::IRegular_volume>(volume_tag));
    assert(volume_element.is_valid_interface());

    const mi::math::Bbox_struct<mi::Float32, 3> volume_ijk_roi_bbox_float =
        volume_element->get_IJK_region_of_interest();

    const mi::math::Bbox<mi::Sint32, 3> volume_ijk_roi_bbox_uint(
        nv::index_common::convert_bbox_type<mi::Sint32, mi::Float32, 3>(volume_ijk_roi_bbox_float));

    return volume_ijk_roi_bbox_uint;
}

//----------------------------------------------------------------------
mi::math::Bbox<mi::Sint32, 3> get_heightfield_roi_bbox(
    const mi::neuraylib::Tag&         heightfield_tag,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    assert(heightfield_tag.is_valid());
    assert(dice_transaction != 0);

    mi::base::Handle<const nv::index::IRegular_heightfield> heightfield_element(
        dice_transaction->access<const nv::index::IRegular_heightfield>(heightfield_tag));
    assert(heightfield_element.is_valid_interface());

    const mi::math::Bbox_struct<mi::Float32, 3> heightfield_ijk_roi_bbox_float =
        heightfield_element->get_IJK_region_of_interest();

    const mi::math::Bbox<mi::Sint32, 3> heightfield_ijk_roi_bbox_sint(
        nv::index_common::convert_bbox_type<mi::Sint32, mi::Float32, 3>(heightfield_ijk_roi_bbox_float));

    return heightfield_ijk_roi_bbox_sint;
}

//----------------------------------------------------------------------
mi::math::Bbox<mi::Float32, 3> get_scene_bounding_box(
    mi::neuraylib::Tag                session_tag,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    Scene_traverse_scene_bounding_box scene_bbox_strategy;

    const bool is_succeed = get_scene_status(&scene_bbox_strategy, session_tag, dice_transaction);
    if (is_succeed)
    {
        return scene_bbox_strategy.get_bounding_box();
    }

    return mi::math::Bbox<mi::Float32, 3>();
}

//----------------------------------------------------------------------
void set_examiner_window_resolution(const mi::Sint32_2& win_res)
{
    if (win_res.x > 0)
    {
        assert(win_res.y > 0);
        // To tell the window resolution to nvindex library, use
        // ICamera. In this example application, window resolution is
        // kept in the examiner manipulator.
        Nvindex_AppData::instance()->get_user_interaction(0)->
            get_examiner()->set_main_window_resolution(win_res);

        // Flash player 9 has bitmap size limitation 2880x2880
        // See http://kb2.adobe.com/cps/496/cpsid_49662.html
        //         if ((win_res.x > 2880) || (win_res.y > 2880)){
        //             WARN_LOG << "Too large window resolution for Flash Player 9. "
        //                      << "Adobe Flash Player 9 can show width or height up to 2880 pixels. "
        //                      << "See http://kb2.adobe.com/cps/496/cpsid_49662.html ";
        //         }
    }
    else
    {
        ERROR_LOG << "Illegal window resolution " << win_res << ". ignored.";
    }
}

//----------------------------------------------------------------------
void request_canvas_resolution(
    mi::Sint32_2 const & new_canvas_res)
{
    if ((new_canvas_res.x <= 0) || (new_canvas_res.y <= 0))
    {
        ERROR_LOG << "Cannot request negative resolution: " << new_canvas_res << ", ignored.";
        return;
    }

    if ((new_canvas_res.x > 4096) || (new_canvas_res.y > 4096))
    {
        WARN_LOG << "Requested canvas size seems too large: " << new_canvas_res
                 << ", may cause out of memory error.";
    }

    if ((new_canvas_res.x >= 10000) || (new_canvas_res.y >= 10000))
    {
        ERROR_LOG << "Requested canvas size seems too large: " << new_canvas_res
                  << ", may cause out of memory error, ignored";
        return;
    }

    // When there are multiple controls, please lock this. Currently
    // we don't have that in this application.
    Nvindex_AppData::instance()->set_request_canvas_resolution(new_canvas_res);
}

//----------------------------------------------------------------------
void update_canvas_resolution(Nvindex_rendering_context& irc)
{
    assert(irc.m_span_buffer.is_valid_interface());
    nv::index::IIndex_canvas* const index_canvas = irc.m_span_buffer.get();

    const mi::Sint32_2 requested_res = Nvindex_AppData::instance()->get_request_canvas_resolution();
    assert(requested_res.x > 0);
    assert(requested_res.y > 0);

    // We have 5 canvases/buffers to update.

    // 1. examiner resolution update
    set_examiner_window_resolution(requested_res);

    // 2. index canvas (span buffer) resolution update
    index_canvas->set_buffer_resolution(requested_res);

    // (3. opengl offscreen context if USE_OPENGL implementation)
    Nvindex_AppData::instance()->get_opengl_appdata()->resize_offscreen_context(requested_res.x, requested_res.y);

    // 4. http canvas update
    if (irc.m_http_request_handler.is_valid_interface())
    {
        std::vector<mi::math::Vector<mi::Uint32, 2> > res_vec;
        res_vec.push_back(mi::math::Vector<mi::Uint32, 2>(static_cast<mi::Uint32>(requested_res.x),
                                                          static_cast<mi::Uint32>(requested_res.y)));
        res_vec.push_back(mi::math::Vector<mi::Uint32, 2>(1024, 1024)); // second view, fixed resolution
        irc.m_http_request_handler->set_resolutions(res_vec);
    }

    // 5. encoder canvas will be updated when it copies the data from
    // index_canvas So is will not be updated at here.
}

//----------------------------------------------------------------------
void set_image_canvas_resolution(const mi::Sint32_2& new_image_res)
{
    assert(new_image_res.x > 0);
    assert(new_image_res.y > 0);

    Nvindex_AppData::instance()->peek_app_proj()->
        insert("index::image_file_canvas_resolution",
               nv::index_common::to_string(new_image_res));
}

//----------------------------------------------------------------------
mi::Sint32_2 get_image_canvas_resolution()
{
    const mi::Sint32_2 cur_img_res =
        nv::index_common::get_vec_sint32_2(
            Nvindex_AppData::instance()->peek_app_proj()->
            get("index::image_file_canvas_resolution", "1024 1024"));

    return cur_img_res;
}

//----------------------------------------------------------------------
std::string get_snapshot_filename_by_mode(
    Nvindex_AppData* appdata,
    mi::Uint32       frame_num)
{
    assert(appdata != 0);

    mi::Sint32 snapshot_idx = -1;

    if (appdata->is_snapshot_sequence_mode())
    {
        // use snapshot sequence for the snapshot index
        snapshot_idx = appdata->get_snapshot_index();
        appdata->inc_snapshot_index();
    }
    else
    {
        // use rendering frame number for the snapshot index
        snapshot_idx = static_cast<mi::Sint32>(frame_num);
    }
    assert(snapshot_idx >= 0);

    const std::string fname = appdata->get_snapshot_filename(snapshot_idx);

    return fname;
}

//----------------------------------------------------------------------
void copy_result_pixel_to_canvas(Span_renderer_IF*         span_buffer,
                                 nv::index_common::Canvas* canvas)
{
    assert(span_buffer != 0);
    assert(canvas      != 0);

    const mi::math::Vector<mi::Sint32, 2> span_buffer_size = span_buffer->get_buffer_resolution();

    if ((span_buffer_size.x != static_cast<mi::Sint32>(canvas->get_resolution_x())) ||
        (span_buffer_size.y != static_cast<mi::Sint32>(canvas->get_resolution_y())))
    {
        canvas->resize_buffer(span_buffer_size.x, span_buffer_size.y);
    }

    const mi::base::Handle<mi::neuraylib::ITile> tile(canvas->get_tile(0, 0, 0));
    span_buffer->copy_pixel(canvas->get_resolution_x(),
                            canvas->get_resolution_y(),
                            tile->get_data());
}

//----------------------------------------------------------------------
void get_cluster_host_containing_volume(
    const mi::neuraylib::Tag&              session_tag,
    const mi::neuraylib::Tag&              volume_tag,
    const mi::math::Bbox< mi::Uint32, 3 >& query_bounds,
    std::vector<mi::Uint32>&               cluster_host_id_vec,
    mi::neuraylib::IDice_transaction*      dice_transaction)
{
    assert(dice_transaction != 0);

    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<const nv::index::ISession>(
            session_tag));
    assert(session.is_valid_interface());

    // Access the distribution scheme
    const mi::neuraylib::Tag dist_layout_tag = session->get_distribution_layout();
    assert(dist_layout_tag.is_valid());
    mi::base::Handle<const nv::index::IData_distribution>
        distribution_layout(
            dice_transaction->access<const nv::index::IData_distribution>(
                dist_layout_tag));
    assert(distribution_layout.is_valid_interface());

    // Distribution layout
    mi::base::Handle<nv::index::IRegular_volume_data_locality> data_locality(
        distribution_layout->retrieve_data_locality(
            volume_tag, query_bounds, dice_transaction));
    assert(data_locality.is_valid_interface());

    // Determine all host ids in the cluster
    for (mi::Uint32 i = 0; i < data_locality->get_nb_cluster_nodes(); ++i)
    {
        // Host id 0 specifies that this subregion has not yet been loaded.
        // This should not happen.
        assert(data_locality->get_cluster_node(i) != 0);
        cluster_host_id_vec.push_back(data_locality->get_cluster_node(i));
    }
}

//----------------------------------------------------------------------
mi::math::Bbox<mi::Float32, 3> get_XYZ_global_region_of_interest_bbox(
    const nv::index::ISession*        session,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    assert(dice_transaction != 0);
    assert(session != 0);
    mi::base::Handle<const nv::index::IScene> scene(dice_transaction->access<const nv::index::IScene>(session->get_scene()));
    assert(scene.is_valid_interface());
    const mi::math::Bbox<mi::Float32, 3> global_roi_bbox = scene->get_clipped_bounding_box();

    return global_roi_bbox;
}

//----------------------------------------------------------------------
void set_XYZ_global_region_of_interest_bbox(
    const mi::math::Bbox<mi::Float32, 3>& xyz_bbox,
    const nv::index::ISession*            session,
    mi::neuraylib::IDice_transaction*     dice_transaction)
{
    assert(session != 0);
    assert(dice_transaction != 0);

    mi::base::Handle<nv::index::IScene> scene(
        dice_transaction->edit<nv::index::IScene>(
            session->get_scene()));
    assert(scene.is_valid_interface());

    scene->set_clipped_bounding_box(xyz_bbox);
}

//----------------------------------------------------------------------
mi::math::Bbox< mi::Sint64, 3 > get_volume_local_IJK_ROI_bbox(
    const mi::neuraylib::Tag&         volume_tag,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    assert(dice_transaction != 0);
    assert(volume_tag.is_valid());

    mi::math::Bbox< mi::Sint64, 3 > ijk_volume_roi_bbox;
    {
        mi::base::Handle<const nv::index::IRegular_volume> volume_element(
            dice_transaction->access<const nv::index::IRegular_volume>(volume_tag));
        assert(volume_element.is_valid_interface());

        ijk_volume_roi_bbox = nv::index_common::convert_bbox_type<mi::Sint64, mi::Float32, 3>(volume_element->get_IJK_region_of_interest());

        // sanity check. ROI in volume_size (need this constructor
        // since return value is a reference not to make a reference
        // of a reference.)
        assert(bbox_inclusive_contain(
                   nv::index_common::convert_bbox_type<mi::Sint64, mi::Uint32, 3>(
                       mi::math::Bbox< mi::Uint32, 3 >(volume_element->get_IJK_bounding_box())),
                   ijk_volume_roi_bbox));
    }

    return ijk_volume_roi_bbox;
}

//----------------------------------------------------------------------
void set_volume_local_IJK_ROI_bbox(
    const mi::neuraylib::Tag&            volume_tag,
    const mi::math::Bbox<mi::Sint64, 3>& volume_ijk_roi,
    mi::neuraylib::IDice_transaction*    dice_transaction)
{
    assert(dice_transaction != 0);
    assert(volume_tag.is_valid());

    if (volume_ijk_roi.empty())
    {
        WARN_LOG << "can not set empty volume_ijk_roi, ignored.";
        return;
    }

    mi::base::Handle< nv::index::IRegular_volume> volume_element(
        dice_transaction->edit< nv::index::IRegular_volume>(volume_tag));
    assert(volume_element.is_valid_interface());

    volume_element->set_IJK_region_of_interest(
        nv::index_common::convert_bbox_type<mi::Float32, mi::Sint64, 3>(volume_ijk_roi));
}

//----------------------------------------------------------------------
mi::math::Bbox< mi::Sint64, 3 > get_heightfield_local_IJK_ROI_bbox(
    mi::neuraylib::IDice_transaction* dice_transaction,
    const mi::neuraylib::Tag&         heightfield_tag)
{
    assert(dice_transaction != 0);
    assert(heightfield_tag.is_valid());

    mi::math::Bbox< mi::Sint64, 3 > ijk_heightfield_roi_bbox;
    {
        mi::base::Handle<const nv::index::IRegular_heightfield> heightfield(
            dice_transaction->access<const nv::index::IRegular_heightfield>(heightfield_tag));
        assert(heightfield.is_valid_interface());

        ijk_heightfield_roi_bbox =
            nv::index_common::convert_bbox_type<mi::Sint64, mi::Float32, 3>(
                heightfield->get_IJK_region_of_interest());

        // sanity check. ROI in volume_size (need this constructor
        // since return value is a reference not to make a reference
        // of a reference.)

        // bbox_inclusive_contain
        mi::math::Bbox< mi::Sint64, 3 > const hbbox =
            nv::index_common::convert_bbox_type<mi::Sint64, mi::Float32, 3>(
                mi::math::Bbox< mi::Float32, 3 >(heightfield->get_IJK_bounding_box()));
        assert(patch_bbox_inclusive_contain(hbbox, ijk_heightfield_roi_bbox));
        nv::index_common::no_unused_variable_warning_please(hbbox);
    }

    return ijk_heightfield_roi_bbox;
}

//----------------------------------------------------------------------
void set_heightfield_local_IJK_ROI_bbox(
    const mi::neuraylib::Tag&            heightfield_tag,
    const mi::math::Bbox<mi::Sint64, 3>& heightfield_ijk_roi,
    mi::neuraylib::IDice_transaction*    dice_transaction)
{
    assert(heightfield_tag.is_valid());
    assert(dice_transaction != 0);

    if (heightfield_ijk_roi.empty())
    {
        WARN_LOG << "can not set empty heightfield_ijk_roi, ignored.";
        return;
    }

    mi::base::Handle< nv::index::IRegular_heightfield> heightfield(
        dice_transaction->edit< nv::index::IRegular_heightfield>(heightfield_tag));
    assert(heightfield.is_valid_interface());

    heightfield->set_IJK_region_of_interest(
        nv::index_common::convert_bbox_type<mi::Float32, mi::Sint64, 3>(heightfield_ijk_roi));
}


//----------------------------------------------------------------------
Scene_element_type_e get_scene_element_type(
    const mi::neuraylib::Tag&         scene_element_tag,
    mi::Sint32&                       return_idx,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    assert(dice_transaction != 0);
    assert(scene_element_tag.is_valid());
    return_idx = -1;            // not found yet

    // FIXME: get_interface can do these.
    {
        const std::vector<mi::neuraylib::Tag> volume_tag_vec = Nvindex_AppData::instance()->get_volume_tag_vec();
        mi::Uint32 const                      volume_count   = volume_tag_vec.size();
        for (mi::Uint32 i = 0; i < volume_count; ++i)
        {
            mi::neuraylib::Tag const          vol_tag        = volume_tag_vec.at(i);
            if (vol_tag == scene_element_tag)
            {
                // found
                return_idx                                   = static_cast< mi::Sint32 >(i);
                // DEBUG_LOG << "tag: " << scene_element_tag.id << " is volume: " << i;
                return SE_volume;
            }
        }

        const std::vector<mi::neuraylib::Tag> hf_tag_vec     = Nvindex_AppData::instance()->get_heightfield_tag_vec();
        mi::Uint32 const                      nb_heightfield = hf_tag_vec.size();
        for (mi::Uint32 i = 0; i < nb_heightfield; ++i)
        {
            mi::neuraylib::Tag const          hf_tag         = hf_tag_vec.at(i);
            if (hf_tag == scene_element_tag)
            {
                // found
                return_idx                                   = static_cast< mi::Sint32 >(i);
                // DEBUG_LOG << "tag: " << scene_element_tag.id << " is heightfield: " << i;
                return SE_heightfield;
            }
        }

        {
            mi::base::Handle<const nv::index::IPoint_set> point(
                dice_transaction->access<const nv::index::IPoint_set>(
                    scene_element_tag));

            if (point)
            {
                return_idx = -1;
                return SE_point;
            }
        }

        {
            mi::base::Handle<const nv::index::IPlane> plane(
                dice_transaction->access<const nv::index::IPlane>(
                    scene_element_tag));

            if (plane)
            {
                return_idx = -1;
                return SE_plane;
            }
        }

        // TODO: to support slices. add slice check code here.
    }
    DEBUG_LOG << "tag: " << scene_element_tag.id << " is neither volume, heightfield, point, plane. "
              << "(slice is not yet supported in this function, but possible to implement. "
              << "See the comment in the code.)";

    return SE_none;
}

//----------------------------------------------------------------------
void set_slice_change(
    mi::neuraylib::IScope* scope,
    mi::Uint32                                current_seismic_volume_index,
    mi::Uint32                                slice_id,
    mi::Float32                               slice_position,
    const mi::neuraylib::Tag&                 slice_colormap_tag,
    bool is_slice_enabled
    )
{
    assert(scope != 0);

    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        scope->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    {
        bool                                  done           = false;
        const std::vector<mi::neuraylib::Tag> volume_tag_vec = Nvindex_AppData::instance()->get_volume_tag_vec();
        if (current_seismic_volume_index < volume_tag_vec.size())
        {
            mi::base::Handle<const nv::index::IRegular_volume> volume(
                dice_transaction->access<const nv::index::IRegular_volume>(
                    volume_tag_vec.at(current_seismic_volume_index)));
            assert(volume.is_valid_interface());

            if (slice_id < volume->get_nb_slices())
            {
                done                               = true;
                mi::neuraylib::Tag const slice_tag = volume->get_slice(slice_id);
                assert(slice_tag.is_valid());

                mi::base::Handle<nv::index::ISlice_scene_element> slice(
                    dice_transaction->edit<nv::index::ISlice_scene_element>(slice_tag));
                assert(slice.is_valid_interface());

                mi::base::Handle<nv::index::ISection_scene_element> section_slice(
                    slice->get_interface<nv::index::ISection_scene_element>());
                if (section_slice.get())
                {
                    section_slice->set_position(slice_position);
                }

                if (slice_colormap_tag.is_valid())
                {
                    slice->assign_colormap(slice_colormap_tag);
                }
                slice->set_enabled(is_slice_enabled);
            }
        }

        // Handle interaction group
        const nv::index_common::String_dict* dict       = Nvindex_AppData::instance()->peek_app_proj();
        const std::string                    group_name = dict->get("app::interaction_group::group");
        if (!done && !group_name.empty())
        {
            mi::neuraylib::Tag group_tag = dice_transaction->name_to_tag(group_name.c_str());
            if (group_tag)
            {
                mi::base::Handle<nv::index::ITransformed_scene_group> group(
                    dice_transaction->edit<nv::index::ITransformed_scene_group>(group_tag));
                if (group)
                {
                    const mi::math::Vector<mi::Float32, 3> dir
                        = nv::index_common::get_vec_float32_3(dict->get("app::interaction_group::translate", "1 0 0"));
                    const mi::Float32 amount
                        = nv::index_common::get_float32(dict->get("app::interaction_group::amount", "1"));

                    mi::math::Matrix<mi::Float32, 4, 4> transform(1.f);
                    transform.translate(dir * amount * slice_position);
                    group->set_transform(transform);
                }
                else
                {
                    ERROR_LOG << "Could not access an ITransformed_scene_group using the name '"
                              << group_name << "' as specified in app::interaction_group, tag "
                              << group_tag << ".";
                }
            }
            else
            {
                ERROR_LOG << "Could not find a database element with the name '" << group_name << "' "
                          << "as specified in app::interaction_group.";
            }
        }
    }
    dice_transaction->commit();
}

//----------------------------------------------------------------------
nv::index::IRegular_heightfield_data_locality * new_whole_heightfield_distribution_layout(
    const mi::neuraylib::Tag&             session_tag,
    const mi::neuraylib::Tag&             heightfield_tag,
    bool                                  is_edit,
    mi::neuraylib::IDice_transaction *    dice_transaction,
    mi::math::Bbox_struct<mi::Uint32, 3>& heightfield_bbox)
{
    assert(session_tag.is_valid());
    assert(heightfield_tag.is_valid());
    assert(dice_transaction != 0);

    // Query the heightfield extent/bounds
    mi::base::Handle<const nv::index::IRegular_heightfield> heightfield(
        dice_transaction->access<const nv::index::IRegular_heightfield>(heightfield_tag));
    assert(heightfield.is_valid_interface());

    // return the heightfield bbox
    mi::math::Bbox<mi::Float32, 3> const & float_bbox = heightfield->get_IJK_bounding_box();

    // convert to heightfield bbox
    for (mi::Sint32 i = 0; i < 3; ++i)
    {
        // check only
        assert(floor(float_bbox.min[i]) >= 0.0);
        assert(ceil( float_bbox.max[i]) >= 0.0);
    }
    heightfield_bbox.min.x = static_cast< mi::Uint32 >(floor(float_bbox.min.x));
    heightfield_bbox.max.x = static_cast< mi::Uint32 >(ceil( float_bbox.max.x));
    heightfield_bbox.min.y = static_cast< mi::Uint32 >(floor(float_bbox.min.y));
    heightfield_bbox.max.y = static_cast< mi::Uint32 >(ceil( float_bbox.max.y));
    heightfield_bbox.min.z = static_cast< mi::Uint32 >(floor(float_bbox.min.z));
    heightfield_bbox.max.z = static_cast< mi::Uint32 >(ceil( float_bbox.max.z));

    // convert to heightfield area for query
    mi::math::Bbox_struct<mi::Uint32, 2> heightfield_area;
    heightfield_area.min.x = heightfield_bbox.min.x;
    heightfield_area.max.x = heightfield_bbox.max.x;
    heightfield_area.min.y = heightfield_bbox.min.y;
    heightfield_area.max.y = heightfield_bbox.max.y;

    nv::index::IRegular_heightfield_data_locality * p_locality = new_heightfield_distribution_layout(
        session_tag,
        heightfield_tag,
        heightfield_area,
        is_edit,
        dice_transaction);

    return p_locality;
}


//----------------------------------------------------------------------
nv::index::IRegular_heightfield_data_locality * new_heightfield_distribution_layout(
    const mi::neuraylib::Tag&                   session_tag,
    const mi::neuraylib::Tag&                   heightfield_tag,
    const mi::math::Bbox_struct<mi::Uint32, 2>& ij_query_bounds,
    bool                                        is_edit,
    mi::neuraylib::IDice_transaction*           dice_transaction)
{
    assert(session_tag.is_valid());
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

    // Heightfiled distribution layout: this is new-ed. Use handle at the
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

//----------------------------------------------------------------------
nv::index::IRegular_volume_data_locality* new_volume_distribution_layout(
    const mi::neuraylib::Tag&                   session_tag,
    const mi::neuraylib::Tag&                   volume_tag,
    const mi::math::Bbox_struct<mi::Uint32, 3>& ijk_query_bounds,
    mi::neuraylib::IDice_transaction*           dice_transaction)
{
    assert(session_tag.is_valid());
    assert(volume_tag.is_valid());
    assert(dice_transaction != 0);

    // Access session
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<const nv::index::ISession>(session_tag));
    assert(session.is_valid_interface());

    // Access to the distribution scheme
    const mi::neuraylib::Tag & dist_layout_tag = session->get_distribution_layout();
    assert(dist_layout_tag.is_valid());
    mi::base::Handle<const nv::index::IData_distribution> distribution_layout(
        dice_transaction->access<const nv::index::IData_distribution>(dist_layout_tag));
    assert(distribution_layout.is_valid_interface());

    // Seismic distribution layout
    nv::index::IRegular_volume_data_locality* data_locality = distribution_layout->retrieve_data_locality(
        volume_tag,
        ijk_query_bounds,
        dice_transaction);

    return data_locality;
}

//----------------------------------------------------------------------
static void traverse_scene_hierarchy(
    const mi::neuraylib::Tag&          scene_element_tag,
    Scene_traverse_strategy_if * const p_strategy,
    mi::neuraylib::IDice_transaction*  dice_transaction)
{
    assert(scene_element_tag.is_valid());
    assert(p_strategy != 0);
    assert(dice_transaction != 0);

    mi::base::Handle< const nv::index::IScene_element > scene_element(
        dice_transaction->access<const nv::index::IScene_element>(scene_element_tag));
    assert(scene_element.is_valid_interface());

    mi::base::Handle< const nv::index::IScene_group > scene_group(
        scene_element->get_interface<const nv::index::IScene_group>());

    if (scene_group.is_valid_interface())
    {
        // group
        if (p_strategy->visit_only_enabled_groups())
        {
            if (!(scene_group->get_enabled()))
            {
                return;
            }
        }
    }
    else
    {
        // no group
        if (p_strategy->visit_only_enabled_scene_elements())
        {
            if (!(scene_element->get_enabled()))
            {
                return;
            }
        }
    }

    p_strategy->run_before_recursion(scene_element_tag, dice_transaction);

    if (scene_group.is_valid_interface())
    {
        const mi::Sint32 nb_element = scene_group->nb_elements();
        for (mi::Sint32 i = 0; i < nb_element; ++i)
        {
            const mi::neuraylib::Tag child_scene_element_tag = scene_group->get_scene_element(i);
            {
                mi::base::Handle<const nv::index::IScene_element> child_scene_element(
                    dice_transaction->access<const nv::index::IScene_element>(child_scene_element_tag));
                if (!child_scene_element.is_valid_interface())
                {
                    WARN_LOG << "traverse_scene_hierarchy: child scene element ("
                             << child_scene_element_tag.id << ") is not a type of IScene_element.";
                    return;
                }
            }
            traverse_scene_hierarchy(child_scene_element_tag, p_strategy, dice_transaction);
        }
    }

    p_strategy->run_after_recursion(scene_element_tag, dice_transaction);
}

//----------------------------------------------------------------------
bool traverse_scene_hierarchy_by_tag(
    const mi::neuraylib::Tag           scene_group_tag,
    Scene_traverse_strategy_if * const p_strategy,
    mi::neuraylib::IDice_transaction*  dice_transaction)
{
    assert(scene_group_tag.is_valid());
    assert(p_strategy != 0);
    assert(dice_transaction != 0);

    mi::base::Handle<const nv::index::IScene_group> scene_group(
        dice_transaction->access<const nv::index::IScene_group>(scene_group_tag));
    if (!scene_group.is_valid_interface())
    {
        WARN_LOG << "traverse_scene_hierarchy_by_tag: scene_group_tag ("
                 << scene_group_tag.id << ") is not a type of IScene_group. Can not traverse.";
        return false;;
    }

    // start traverse
    traverse_scene_hierarchy(scene_group_tag,
                             p_strategy,
                             dice_transaction);
    return true;
}

//----------------------------------------------------------------------
void export_session(
    const Nvindex_rendering_context&  irc,
    mi::neuraylib::IDice_transaction* dice_transaction,
    const std::string&                output_filename)
{
    assert(dice_transaction != 0);

    //TODO: The session exporter should be modified so that the scene update lock is not needed
    mi::base::Lock::Block block_update(&Nvindex_AppData::instance()->m_scene_update_lock);

    mi::base::Handle<mi::neuraylib::IFactory> factory(
        irc.m_iindex_if->get_api_component<mi::neuraylib::IFactory>());
    mi::base::Handle<mi::IString> str(factory->create<mi::IString>("String"));

    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<nv::index::ISession>(irc.m_session_tag));

    const nv::index::IViewport_list* viewport_list = 0;
    if (irc.m_viewport_list.is_valid_interface() && irc.m_viewport_list->size() > 0)
    {
        viewport_list = irc.m_viewport_list.get();
    }

    session->export_session(
        nv::index::ISession::EXPORT_DEFAULT | nv::index::ISession::EXPORT_USER_DATA,
        str.get(), dice_transaction, viewport_list);

    std::ostringstream s;
    s << str->get_c_str();

    const nv::index_common::String_dict* p_app_proj = Nvindex_AppData::instance()->peek_app_proj();

    //
    // Special handling for application-side setttings
    //

    // Interaction group
    const std::string interact = p_app_proj->get("app::interaction_group::group");
    if (!interact.empty())
    {
        s << "\n# Interaction group (application side)\n"
          << "app::interaction_group::group = "
          << p_app_proj->get("app::interaction_group::group") << "\n"
          << "app::interaction_group::translate = "
          << p_app_proj->get("app::interaction_group::translate", "0 0 0") << "\n"
          << "app::interaction_group::amount = "
          << p_app_proj->get("app::interaction_group::amount", "0") << "\n"
          << "app::interaction_group::min = "
          << p_app_proj->get("app::interaction_group::min", "-1000") << "\n"
          << "app::interaction_group::max = "
          << p_app_proj->get("app::interaction_group::max", "-1000") << "\n";
    }

    if (output_filename.empty() || output_filename  == "stdout")
    {
        std::cout << "\n"
                  << "----------------8<-------------[ cut here ]------------------\n"
                  << s.str()
                  << "----------------8<-------------[ cut here ]------------------\n" << std::endl;
    }
    else
    {
        INFO_LOG << "Writing export to file '" << output_filename << "'";
        std::ofstream f(output_filename.c_str());
        f << s.str();
        f.close();
    }
}

//----------------------------------------------------------------------
bool get_scene_status(
    Scene_traverse_strategy_if* p_scene_traverse_strategy,
    mi::neuraylib::Tag          session_tag,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    assert(p_scene_traverse_strategy != 0);
    assert(session_tag.is_valid());
    assert(dice_transaction != 0);

    bool is_success = false;
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<nv::index::ISession>(session_tag));
    assert(session.is_valid_interface());

    const mi::neuraylib::Tag scene_tag = session->get_scene();
    assert(scene_tag.is_valid());
    is_success = traverse_scene_hierarchy_by_tag(scene_tag,
                                                 p_scene_traverse_strategy,
                                                 dice_transaction);
    assert(is_success);

    return is_success;
}

//----------------------------------------------------------------------
bool get_center_of_volume_in_roi(
    Nvindex_rendering_context*        irc_ref,
    const mi::neuraylib::Tag&         volume_tag,
    mi::math::Vector<mi::Float32, 3>& volume_center)
{
    const std::string fn = "get_center_of_volume_in_roi:";
    assert(irc_ref != 0);

    if (!volume_tag.is_valid())
    {
        ERROR_LOG << "Given volume tag is invalid.";
        return false;
    }

    bool is_succeded = true;

    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        irc_ref->get_current_scope()->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    {
        mi::base::Handle<const nv::index::IRegular_volume> volume(
            dice_transaction->access<const nv::index::IRegular_volume>(volume_tag));
        assert(volume.is_valid_interface());

        // This is a volume local clipped bounding box (no relation with global ROI)
        const mi::math::Bbox<mi::Float32, 3> bbox = volume->get_XYZ_clipped_bounding_box();

        // get the scene transformation (FIXME: only get the
        // global transformation, didn't see the static_group_node's transformation)
        assert(irc_ref->m_session_tag.is_valid());
        mi::base::Handle<const nv::index::ISession> session(
            dice_transaction->access<const nv::index::ISession>(irc_ref->m_session_tag));
        assert(session.is_valid_interface());

        const mi::math::Bbox<mi::Float32, 3> global_roi_bbox =
            get_XYZ_global_region_of_interest_bbox(session.get(),
                                                   dice_transaction.get());
        const bool is_intersect = global_roi_bbox.intersects(bbox);
        const mi::math::Bbox<mi::Float32, 3> int_bbox =
            nv::index_common::bbox_and(global_roi_bbox, bbox);

        if (bbox.empty())
        {
            ERROR_LOG << fn << "Cannot have a volume center. The first volume has an empty bbox.";
            is_succeded = false;
        }
        else if (global_roi_bbox.empty())
        {
            ERROR_LOG << fn << "Cannot have a bbox center. Global ROI has empty bbox.";
            is_succeded = false;
        }
        else if (!is_intersect)
        {
            ERROR_LOG << fn << "Cannot have a bbox center. ROI and volume_bbox has no intersect.";
            is_succeded = false;
        }
        else if (!int_bbox.is_volume())
        {
            ERROR_LOG << fn << "Cannot have a bbox center. The intersection has no volume.";
            is_succeded = false;
        }
        else
        {
            // success
            mi::base::Handle<const nv::index::IScene> scene(
                dice_transaction->access<const nv::index::IScene>(session->get_scene()));
            assert(scene.is_valid_interface());
            const mi::math::Bbox<mi::Float32, 3> transformed_bbox =
                get_transformed_bbox(int_bbox, scene->get_transform_matrix());
            volume_center = transformed_bbox.center();;

        }
    }
    dice_transaction->commit();

    return is_succeded;
}

//----------------------------------------------------------------------
