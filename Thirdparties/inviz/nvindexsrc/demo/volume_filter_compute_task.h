/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief volume filter computation task: generating bordermap, apply filter

#ifndef NVIDIA_INDEX_VOLUME_FILTER_COMPUTE_TASK_H
#define NVIDIA_INDEX_VOLUME_FILTER_COMPUTE_TASK_H

#include <mi/dice.h>
#include <nv/index/iregular_volume_compute_task.h>
#include <nv/index/idistributed_data_access.h>

#include "common/encode_voxel.h"

/// volume amplitude data filter preprocess: create bordermap
class Volume_bordermap_generation :
        public nv::index::Regular_volume_compute_task
{
public:
    /// constructor
    ///
    /// \param[in] whole_ijk_roi_bbox specifies the part of the
    /// volume where the operation should be applied.
    /// \param[in] session_tag session tag
    /// \param[in] volume_tag input volume data scene element tag
    /// \param[in] job_local_bordermap_element_tag job local bordermap tag
    Volume_bordermap_generation(
        const mi::math::Bbox< mi::Sint64, 3 >& whole_ijk_roi_bbox,
        const mi::neuraylib::Tag&              session_tag,
        const mi::neuraylib::Tag&              volume_tag,
        const mi::neuraylib::Tag&              job_local_bordermap_element_tag);

    /// Perform computation on the given volume data.
    ///
    /// \param[in] brick_ijk_bbox_struct IJK bounding box of the current
    /// volume brick.
    /// \param[in] dest_amplitude_values Amplitude data of the volume
    /// brick. This is the result destination.
    /// \param[in] dice_transaction      Current transaction.
    virtual bool compute(
        const mi::math::Bbox_struct<mi::Uint32, 3>& brick_ijk_bbox_struct,
        nv::index::IRegular_volume_data*            voxel_data,
        mi::neuraylib::IDice_transaction*           dice_transaction) const;

private:
    /// whole computation IJK region of interest bounding box
    mi::math::Bbox< mi::Sint64, 3 > m_whole_ijk_roi_bbox;
    /// session tag
    mi::neuraylib::Tag m_session_tag;
    /// input volume tag
    mi::neuraylib::Tag m_input_volume_tag;
    /// bordermap element tag
    mi::neuraylib::Tag m_bordermap_element_tag;
};


/// volume amplitude data filter: apply a volume filter
class Volume_filter_apply :
    public nv::index::Regular_volume_compute_task
{
public:
    /// constructor
    ///
    /// \param[in] whole_ijk_roi_bbox specifies the part of the
    /// volume where the operation should be applied.
    /// \param[in] session_tag session tag
    /// \param[in] volume_tag input volume tag
    /// \param[in] bordermap_element_tag bordermap db element tag
    /// \param[in] filter_type filter type identifier
    /// \param[in] fragment_index fragment job index (thread ID)
    Volume_filter_apply(
        const mi::math::Bbox< mi::Sint64, 3 >& whole_ijk_roi_bbox,
        const mi::neuraylib::Tag&              session_tag,
        const mi::neuraylib::Tag&              volume_tag,
        const mi::neuraylib::Tag&              bordermap_element_tag,
        mi::Sint32                             filter_type,
        mi::Uint32                             fragment_index);

    /// Perform computation on the given volume data.
    ///
    /// \param[in] brick_ijk_bbox_struct IJK bounding box of the current
    /// volume brick.
    /// \param[in] dest_amplitude_values Amplitude data of the volume
    /// brick. This is the result destination.
    /// \param[in] dice_transaction      Current transaction.
    virtual bool compute(
        const mi::math::Bbox_struct<mi::Uint32, 3>&  brick_ijk_bbox_struct,
        nv::index::IRegular_volume_data*             voxel_data,
        mi::neuraylib::IDice_transaction*            dice_transaction) const;

private:
    /// access to the volume data, non border part
    ///
    /// \param[in] access_factory     access factory for the volume data
    /// \param[in] valid_src_bbox     valid computing source bounding box
    /// \param[in] dice_transaction   dice database transaction
    /// \param[out] (output) p_volume voxel data to be fill in (the central part)
    void access_volume_data(
        const nv::index::IDistributed_data_access_factory* access_factory,
        const mi::math::Bbox<mi::Sint64, 3>&               valid_src_bbox,
        mi::neuraylib::IDice_transaction*                  dice_transaction,
        nv::index_common::VoxelType*                       p_volume) const;

    /// fill in border data in dice database
    /// \param[in] output_bbox        output bounding box. The brick of this job
    /// \param[in] valid_src_bbox     valid computing source bounding box
    /// \param[in] dice_transaction   dice database transaction
    /// \param[out] (output) p_volume voxel data to be fill in (the border part)
    void fill_volume_border_data(
        const mi::math::Bbox<mi::Sint64, 3>& output_bbox,
        const mi::math::Bbox<mi::Sint64, 3>& valid_src_bbox,
        mi::neuraylib::IDice_transaction*    dice_transaction,
        nv::index_common::VoxelType*         p_volume) const;

private:
    /// whole computation IJK region of interest bounding box
    mi::math::Bbox< mi::Sint64, 3 > m_whole_ijk_roi_bbox;
    /// session tag
    mi::neuraylib::Tag m_session_tag;
    /// input volume tag
    mi::neuraylib::Tag m_input_volume_tag;
    /// bordermap element tag
    mi::neuraylib::Tag m_job_local_bordermap_tag;
    /// filter type
    mi::Sint32 m_filter_type;
    /// fragment index
    mi::Uint32 m_fragment_index;
};


#endif // #ifndef NVIDIA_INDEX_VOLUME_FILTER_COMPUTE_TASK_H
