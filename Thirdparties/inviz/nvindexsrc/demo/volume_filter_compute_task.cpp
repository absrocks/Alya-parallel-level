/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "volume_filter_compute_task.h"

#include <nv/index/isession.h>
#include <nv/index/iconfig_settings.h>

#include "common/common_utility.h"
#include "common/forwarding_logger.h"
#include "common/type_conversion_utility.h"

#include "utilities.h"
#include "volume_bordermap_database_element.h"
#include "volume_brick_element.h"
#include "volume_data_filter.h"

#include <cassert>
#include <vector>

// self test mode
// #define IS_SELF_TEST 0

//----------------------------------------------------------------------
/// get border bboxes
/// \param[in]  src_bbox   compute source bounding box
/// \param[in]  brick_bbox brick bounding box. compute destination bounding box
/// \param[out] result_border_vec border bounding box list.
static void get_border_bbox(
    const mi::math::Bbox< mi::Sint64, 3 >&          src_bbox,
    const mi::math::Bbox< mi::Sint64, 3 >&          brick_bbox,
    std::vector< mi::math::Bbox< mi::Sint64, 3 > >& result_border_vec)
{
    // looking for the 6 boundaries
    // char const * const p_dir[] = {"i", "j", "k", 0, };

    // DEBUG_LOG << "src: " << src_bbox << ", brick: " << brick_bbox;
    for (mi::Sint32 i = 0; i < 3; ++i)
    {
        if (src_bbox.min[i] < brick_bbox.min[i])
        {
            mi::math::Bbox< mi::Sint64, 3 > border_bbox = src_bbox;
            border_bbox.max[i] = src_bbox.min[i] + 2; // need two wide boundary
            result_border_vec.push_back(border_bbox);
            // DEBUG_LOG << "border " << p_dir[i] << "-:" << border_bbox;
        }
        if (src_bbox.max[i] > brick_bbox.max[i])
        {
            mi::math::Bbox< mi::Sint64, 3 > border_bbox = src_bbox;
            border_bbox.min[i] = src_bbox.max[i] - 2;
            result_border_vec.push_back(border_bbox);
            // DEBUG_LOG << "border " << p_dir[i] << "+:" << border_bbox;
        }
    }
    // DEBUG_LOG << "done";
}

//----------------------------------------------------------------------
/// test the border elements
void verify_border_elements(
    mi::neuraylib::IDice_transaction* dice_transaction,
    const mi::neuraylib::Tag          bordermap_element_tag,
    const std::string&                mes)
{
    // access to the bordermap
    assert(bordermap_element_tag.is_valid());
    mi::base::Handle< Volume_bordermap_database_element > bordermap(
        dice_transaction->edit< Volume_bordermap_database_element >(bordermap_element_tag));
    assert(bordermap.is_valid_interface());

    const Volume_bordermap_database_element::Border_map* p_map = bordermap->get_bordermap_ref();
    for (Volume_bordermap_database_element::Border_map::const_iterator ti = p_map->begin();
        ti != p_map->end(); ++ti)
    {
        const mi::neuraylib::Tag btag = ti->second;
        assert(btag.is_valid());
        mi::base::Handle< const Volume_brick_element > border_elem(
            dice_transaction->access< Volume_brick_element >(btag));
        const std::string tag_to_key =
            nv::index_common::to_string(
                nv::index_common::convert_bbox_type<mi::Sint64, mi::Sint32, 3>(border_elem->get_bounding_box()));
        if (tag_to_key == ti->first)
        {
            INFO_LOG << "verified key: " << tag_to_key;
        }
        else
        {
            ERROR_LOG << "Stored key: " << ti->first
                      << " and key from tag->element: " << tag_to_key
                      << "doesn't match.";
        }
        const mi::math::Bbox< mi::Sint64, 3 > bbox(
            nv::index_common::convert_bbox_type<mi::Sint64, mi::Sint32, 3>(border_elem->get_bounding_box()));
        const bool is_verified = test_synthetic_voxel_value(
            mes, bbox, bbox, nv::index_common::VDE_IJ, border_elem->get_voxel_data());
        if (!is_verified)
        {
            ERROR_LOG << "Fail to access the border element at " << mes << ".";
        }
    }
}


//----------------------------------------------------------------------
Volume_bordermap_generation::Volume_bordermap_generation(
    const mi::math::Bbox< mi::Sint64, 3 >& whole_ijk_roi_bbox,
    const mi::neuraylib::Tag&              session_tag,
    const mi::neuraylib::Tag&              input_scene_element_tag,
    const mi::neuraylib::Tag&              job_local_bordermap_element_tag)
    :
    m_whole_ijk_roi_bbox(whole_ijk_roi_bbox),
    m_session_tag(session_tag),
    m_input_volume_tag(input_scene_element_tag),
    m_bordermap_element_tag(job_local_bordermap_element_tag)
{
    // empty
}

//----------------------------------------------------------------------
bool Volume_bordermap_generation::compute(
    const mi::math::Bbox_struct<mi::Uint32, 3>& brick_ijk_bbox_struct,
    nv::index::IRegular_volume_data*            voxel_data,
    mi::neuraylib::IDice_transaction*           dice_transaction) const
{
    //----------------------------------------------------------------------
    // get the volume_data_access

    // Access the session
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->
        access<const nv::index::ISession>(m_session_tag));
    assert(session.is_valid_interface());
    // Access the data access factory
    const mi::neuraylib::Tag& data_access_tag = session->get_data_access_factory();
    mi::base::Handle<const nv::index::IDistributed_data_access_factory> access_factory(
        dice_transaction->
        access<const nv::index::IDistributed_data_access_factory>(data_access_tag));
    assert(access_factory.is_valid_interface());

    // check the subcube size
    mi::base::Handle<const nv::index::IConfig_settings> config_settings(
        dice_transaction->access<nv::index::IConfig_settings>(session->get_config()));
    assert(config_settings.is_valid_interface());
    mi::math::Vector_struct<mi::Uint32, 3> subcube_size;
    config_settings->get_subcube_size(subcube_size);
    if (!Volume_data_filter::is_valid_subcube_size(subcube_size))
    {
        const mi::Sint32 need_w = Volume_data_filter::VDF_extra_comp_width;
        ERROR_LOG << "Can not get source volume. required subcube_size >= ("
                  << need_w << ", " << need_w << ", " << need_w << "), "
                  << subcube_size;        // This is an assumption.
        return false;
    }

    //----------------------------------------------------------------------
    // compute source volume bbox for a 3x3x3 filter kernel.
    const mi::math::Vector< mi::Sint64, 3 > OneVec3(1, 1, 1);

    const mi::math::Bbox< mi::Sint64, 3 > output_bbox(
        nv::index_common::convert_bbox_type<mi::Sint64, mi::Uint32, 3>(brick_ijk_bbox_struct));
    const mi::math::Bbox< mi::Sint64, 3 > valid_src_bbox =
        Volume_data_filter::get_valid_compute_source_bbox(output_bbox, m_whole_ijk_roi_bbox);

    if (valid_src_bbox.empty())
    {
        DEBUG_LOG << "no compute task for this node (out of volume IJK ROI)";
        return true;
    }

    // check the local duplication
    std::vector< mi::math::Bbox< mi::Sint64, 3 > > border_bbox_vec;
    {
        // get the border bbox candidate.
        std::vector< mi::math::Bbox< mi::Sint64, 3 > > border_cand_bbox_vec;
        get_border_bbox(valid_src_bbox, output_bbox, border_cand_bbox_vec);

        assert(m_bordermap_element_tag.is_valid());
        mi::base::Handle< const Volume_bordermap_database_element > bordermap(
            dice_transaction->access< const Volume_bordermap_database_element >(
                m_bordermap_element_tag));
        assert(bordermap.is_valid_interface());

        for (size_t i = 0; i < border_cand_bbox_vec.size(); ++i)
        {
            const mi::math::Bbox< mi::Sint64, 3 > bbox = border_cand_bbox_vec[i];
            const mi::neuraylib::Tag btag = bordermap->get_element_tag(bbox);
            if (btag.is_valid())
            {
                // found duplication
                DEBUG_LOG << "found local duplication: " << bbox;
                continue;
            }
            border_bbox_vec.push_back(bbox);
        }
    }

    // access and store the border
    std::vector< mi::neuraylib::Tag > added_border_tag_vec;
    for (size_t i = 0; i < border_bbox_vec.size(); ++i)
    {
        mi::math::Bbox< mi::Sint64, 3 > const border_bbox = border_bbox_vec[i];

        // get the volume data
        mi::base::Handle<nv::index::IRegular_volume_data_access> volume_data_access(
            access_factory->create_regular_volume_data_access(m_input_volume_tag));
        assert(volume_data_access.is_valid_interface());

        mi::math::Bbox_struct< mi::Uint32, 3 > border_bbox_st(
            nv::index_common::convert_bbox_type<mi::Uint32, mi::Sint64, 3>(border_bbox));

        volume_data_access->access(border_bbox_st, dice_transaction);
        assert((border_bbox == nv::index_common::convert_bbox_type<mi::Sint64, mi::Uint32, 3>(
                    volume_data_access->get_bounding_box())));

        const mi::base::Handle<const nv::index::IRegular_volume_data> volume_data(
            volume_data_access->get_volume_data());
        const mi::base::Handle<const nv::index::IRegular_volume_data_uint8> volume_data_uint8(
            volume_data->get_interface<const nv::index::IRegular_volume_data_uint8>());

        if (!volume_data_uint8)
        {
            ERROR_LOG << "Volume border-map generation failed: data access on non-uint8 volume type.";
            return false;
        }

        mi::math::Bbox<mi::Uint32, 3> data_bbox_u = volume_data_access->get_bounding_box();
        mi::math::Bbox<mi::Sint32, 3> data_bbox(
            data_bbox_u.min.x, data_bbox_u.min.y, data_bbox_u.min.z,
            data_bbox_u.max.x, data_bbox_u.max.y, data_bbox_u.max.z);
        mi::base::Handle< Volume_brick_element > border_elem
            (new Volume_brick_element(data_bbox, volume_data_uint8->get_voxel_data()));

#ifdef IS_SELF_TEST
        {
            const mi::math::Bbox< mi::Sint64, 3 > bbox =
                nv::index_common::convert_bbox_type<mi::Sint64, mi::Uint32, 3>(
                    volume_data_access->get_bounding_box());
            std::string const mes = "creation_border";
            const bool is_verified = test_synthetic_voxel_value(
                mes, bbox, bbox, nv::index_common::VDE_IJ, volume_data_access->get_volume_data());
            if (!is_verified)
            {
                ERROR_LOG << "Fail to access the border element at " << mes << ".";
            }
        }
#endif  // IS_SELF_TEST

        const mi::neuraylib::Tag border_tag = dice_transaction->store(border_elem.get());
        assert(border_tag.is_valid());
        added_border_tag_vec.push_back(border_tag);
    }

    {
        // store the bordermap
        assert(m_bordermap_element_tag.is_valid());
        mi::base::Handle< Volume_bordermap_database_element > job_local_bmap(
            dice_transaction->edit< Volume_bordermap_database_element >(m_bordermap_element_tag));
        assert(job_local_bmap.is_valid_interface());
        assert(border_bbox_vec.size() == added_border_tag_vec.size());
        for (size_t i = 0; i < added_border_tag_vec.size(); ++i)
        {
            job_local_bmap->insert_bbox_volume_tag(border_bbox_vec[i], added_border_tag_vec[i]);
        }
    }

#ifdef IS_SELF_TEST
    verify_border_elements(dice_transaction, m_bordermap_element_tag, "store and access");
#endif  // IS_SELF_TEST

    return true;
}

//----------------------------------------------------------------------
//======================================================================
//----------------------------------------------------------------------
Volume_filter_apply::Volume_filter_apply(
    const mi::math::Bbox< mi::Sint64, 3 >& whole_ijk_roi_bbox,
    const mi::neuraylib::Tag&              session_tag,
    const mi::neuraylib::Tag&              input_scene_element_tag,
    const mi::neuraylib::Tag&              bordermap_element_tag,
    mi::Sint32                             filter_type,
    mi::Uint32                             fragment_index)
    :
    m_whole_ijk_roi_bbox(whole_ijk_roi_bbox),
    m_session_tag(session_tag),
    m_input_volume_tag(input_scene_element_tag),
    m_job_local_bordermap_tag(bordermap_element_tag),
    m_filter_type(filter_type),
    m_fragment_index(fragment_index)
{
    assert(m_job_local_bordermap_tag.is_valid());
}

//----------------------------------------------------------------------
bool Volume_filter_apply::compute(
    const mi::math::Bbox_struct<mi::Uint32, 3>&  brick_ijk_bbox_struct,
    nv::index::IRegular_volume_data*             voxel_data,
    mi::neuraylib::IDice_transaction*            dice_transaction) const
{
    //----------------------------------------------------------------------
    // compute source volume bbox for a 3x3x3 filter kernel.
    mi::math::Vector< mi::Sint64, 3 > const OneVec3(1, 1, 1);

    const mi::math::Bbox< mi::Sint64, 3 > output_bbox(
        nv::index_common::convert_bbox_type<mi::Sint64, mi::Uint32, 3>(brick_ijk_bbox_struct));

    const mi::math::Bbox< mi::Sint64, 3 > valid_src_bbox =
        Volume_data_filter::get_valid_compute_source_bbox(
            output_bbox, m_whole_ijk_roi_bbox);
    if (valid_src_bbox.empty())
    {
        DEBUG_LOG << "no computing region (may out of IJK ROI).";
        return true;
    }
    
    mi::base::Handle<nv::index::IRegular_volume_data_uint8> amplitude_data(
        voxel_data->get_interface<nv::index::IRegular_volume_data_uint8>());
    if (!amplitude_data)
    {
        ERROR_LOG << "Volume filter operation failed: input data not uint8 data.";
        return false;
    }

    mi::Uint8* dest_amplitude_values = amplitude_data->get_voxel_data_mutable();

    //----------------------------------------------------------------------
    // get the volume_data_access

    // Access the session
    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->
        access<const nv::index::ISession>(m_session_tag));
    assert(session.is_valid_interface());
    // Access the data access factory
    const mi::neuraylib::Tag& data_access_tag = session->get_data_access_factory();
    mi::base::Handle<const nv::index::IDistributed_data_access_factory> access_factory(
        dice_transaction->
        access<const nv::index::IDistributed_data_access_factory>(data_access_tag));
    assert(access_factory.is_valid_interface());

    //----------------------------------------------------------------------
    // perform filter
    // source data storage
    const mi::Sint64 vsize = valid_src_bbox.volume();
    nv::index_common::VoxelType * p_volume = new nv::index_common::VoxelType[vsize];
    {
        DEBUG_LOG << "allocated compute source volume: " << vsize
                  << ", volume IJK ROI = " << m_whole_ijk_roi_bbox;

        // collect the data needed for the computation
        this->access_volume_data(access_factory.get(), valid_src_bbox, dice_transaction, p_volume);
        this->fill_volume_border_data(output_bbox, valid_src_bbox, dice_transaction, p_volume);

        const mi::math::Bbox< mi::Sint64, 3 > raw_dst_bbox(
            output_bbox.min - OneVec3, output_bbox.max + OneVec3);

        // filter computation
        Volume_data_filter sfilter(m_filter_type);
        sfilter.apply_filter(dest_amplitude_values, raw_dst_bbox,
                             p_volume, valid_src_bbox,
                             m_whole_ijk_roi_bbox);
    }
    delete [] p_volume;
    p_volume = 0;

    return true;
}

//----------------------------------------------------------------------
void Volume_filter_apply::access_volume_data(
    const nv::index::IDistributed_data_access_factory* access_factory,
    const mi::math::Bbox<mi::Sint64, 3>&               valid_src_bbox,
    mi::neuraylib::IDice_transaction*                  dice_transaction,
    nv::index_common::VoxelType*                       p_volume) const
{
    assert(!valid_src_bbox.empty());
    const mi::math::Bbox_struct< mi::Uint32, 3 > valid_src_bbox_st(
        nv::index_common::convert_bbox_type<mi::Uint32, mi::Sint64, 3>(valid_src_bbox));

    mi::base::Handle<nv::index::IRegular_volume_data_access> volume_data_access(
        access_factory->create_regular_volume_data_access(m_input_volume_tag));
    assert(volume_data_access.is_valid_interface());
    volume_data_access->access(valid_src_bbox_st, dice_transaction);
    assert(valid_src_bbox == (nv::index_common::convert_bbox_type<mi::Sint64, mi::Uint32, 3>(
                                  volume_data_access->get_bounding_box())));

    const mi::base::Handle<const nv::index::IRegular_volume_data> volume_data(
        volume_data_access->get_volume_data());
    const mi::base::Handle<const nv::index::IRegular_volume_data_uint8> volume_data_uint8(
        volume_data->get_interface<const nv::index::IRegular_volume_data_uint8>());

    if (!volume_data_uint8)
    {
        ERROR_LOG << "Volume-filter application failed: data access on non-uint8 volume type.";
        return;
    }

    const mi::math::Vector< mi::Sint64, 3 > dst_origin = valid_src_bbox.min;
    const bool is_success =
        nv::index_common::copy_partial_volume(p_volume, valid_src_bbox,
                                              volume_data_uint8->get_voxel_data(), valid_src_bbox,
                                              valid_src_bbox, dst_origin);
    assert(is_success); nv::index_common::no_unused_variable_warning_please(is_success);

#ifdef IS_SELF_TEST
    const std::string mes = "volume valid data";
    // only output_bbox_bbox part is guaranteed to be correct in this state.
    const bool is_verified = test_synthetic_voxel_value(
        mes, valid_src_bbox, output_bbox, nv::index_common::VDE_IJ, p_volume);
    if (!is_verified)
    {
        ERROR_LOG << "Fail to access the border element at " << mes << ".";
    }
#endif  // IS_SELF_TEST

#ifdef IS_SELF_TEST
    verify_border_elements(dice_transaction, m_job_local_bordermap_tag, "prepare computing volume");
#endif  // IS_SELF_TEST
}

//----------------------------------------------------------------------
void Volume_filter_apply::fill_volume_border_data(
    const mi::math::Bbox< mi::Sint64, 3 >& output_bbox,
    const mi::math::Bbox< mi::Sint64, 3 >& valid_src_bbox,
    mi::neuraylib::IDice_transaction*      dice_transaction,
    nv::index_common::VoxelType*           p_volume) const
{
    std::vector< mi::math::Bbox< mi::Sint64, 3 > > border_bbox_vec;
    get_border_bbox(valid_src_bbox, output_bbox, border_bbox_vec);

    assert(m_job_local_bordermap_tag.is_valid());
    mi::base::Handle< const Volume_bordermap_database_element > bordermap(
        dice_transaction->access< const Volume_bordermap_database_element >(
            m_job_local_bordermap_tag));
    assert(bordermap.is_valid_interface());

    for (size_t i = 0; i < border_bbox_vec.size(); ++i)
    {
        const mi::math::Bbox< mi::Sint64, 3 > bbox = border_bbox_vec[i];
        const mi::neuraylib::Tag btag = bordermap->get_element_tag(bbox);
        // border map should be exist in create bordermap job
        if (!btag.is_valid())
        {
            DEBUG_LOG << "query bbox: " << bbox
                      << "dump the bordermap: "    << bordermap->to_string() << "\n"
                      << "global bordermap_tag = " << m_job_local_bordermap_tag.id;
        }
        assert(btag.is_valid());

        // access the border data
        mi::base::Handle< const Volume_brick_element > border_volume_data(
            dice_transaction->access< const Volume_brick_element >(btag));
        assert(border_volume_data.is_valid_interface());
        assert(bbox == (nv::index_common::convert_bbox_type<mi::Sint64, mi::Sint32, 3>(
                            border_volume_data->get_bounding_box())));

        const mi::math::Vector< mi::Sint64, 3 > dst_origin = bbox.min;
        const bool is_success =
            nv::index_common::copy_partial_volume(p_volume, valid_src_bbox,
                                                  border_volume_data->get_voxel_data(), bbox,
                                                  bbox, dst_origin);
        assert(is_success);
        nv::index_common::no_unused_variable_warning_please(is_success);
    }

#ifdef IS_SELF_TEST
    const std::string mes = "computation source";
    const bool is_verified = test_synthetic_voxel_value(
        mes, valid_src_bbox, valid_src_bbox, nv::index_common::VDE_IJ, p_volume);
    if (!is_verified)
    {
        ERROR_LOG << "Fail to access the border element at " << mes << ".";
    }
#endif  // IS_SELF_TEST
}

//----------------------------------------------------------------------
