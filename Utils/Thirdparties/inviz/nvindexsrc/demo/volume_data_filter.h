/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief simple example of volume data filter

#ifndef GEOSPATIAL_BIN_GEOSPATIAL_STREAM_VIEWER_VOLUME_DATA_FILTER_H
#define GEOSPATIAL_BIN_GEOSPATIAL_STREAM_VIEWER_VOLUME_DATA_FILTER_H

#include "common/encode_voxel.h"

#include <cassert>

/// simple example of volume data filter
class Volume_data_filter
{
public:
    /// kernel size constants
    enum VDF_kernel_size_e {
        /// 2d kernel (3x3)
        VDF_kernel_2d_size = 9,
        /// 3d kernel (3x3x3)
        VDF_kernel_3d_size = 27,
        /// extra computation width for the computation kernel. 2 for 3x3x3.
        VDF_extra_comp_width = 2
    };

    /// Filter method. Synchronize get_filter_name(),
    /// get_filter_type().
    enum Volume_data_filter_type_e {
        /// Gaussian IJ
        VDFT_Gaussian_IJ,
        /// Gaussian JK
        VDFT_Gaussian_JK,
        /// Gaussian KI
        VDFT_Gaussian_KI,
        /// Sobel IJ
        VDFT_Sobel_IJ,
        /// Sobel JK
        VDFT_Sobel_JK,
        /// Sobel KI
        VDFT_Sobel_KI,
        /// DIP median filter
        VDFT_DIP_MEDIAN,
        /// sentinel
        VDFT_Count
    };


public:
    /// get filter name
    /// \param[in] filter_type filter type id
    /// \return    filter name string
    static std::string get_filter_name(mi::Sint32 filter_type);

    /// get filter ID
    /// \param[in] filter_name filter name
    /// \return filter type ID
    static mi::Sint32  get_filter_type(std::string const & filter_name);

    /// is valid filter type?
    /// \param[in] filter_type filter typre ID
    /// \return true when filter type ID is valid
    static bool  is_valid_filter_type(mi::Sint32 filter_type);

    /// required octree subcube size
    /// \param[in] subcube_size octree subcube size
    /// \return true when subcube size satisfied the requirement
    static bool is_valid_subcube_size(
        mi::math::Vector_struct<mi::Uint32, 3> const & subcube_size);

    /// get valid compute source bounding box for 3d filtering
    ///
    /// \param[in] output_bbox output volume bounding box
    /// \param[in] whole_src_roi_bbox whole source region of interest bounding box
    /// \return valid compute source bounding box
    static mi::math::Bbox< mi::Sint64, 3 > get_valid_compute_source_bbox(
        mi::math::Bbox< mi::Sint64, 3 > const & output_bbox,
        mi::math::Bbox< mi::Sint64, 3 > const & whole_src_roi_bbox
        );

public:
    /// constructor with fiter type
    ///
    /// \param[in] filter_type filter type
    Volume_data_filter(mi::Sint32 filter_type);

    /// destructor
    virtual ~Volume_data_filter();

    /// apply the filter.
    /// The result is in the p_res_data.
    ///
    /// \param[in] p_dst_data   destination result data voxel array
    /// \param[in] raw_dst_bbox destination raw bounding box
    /// \param[in] p_src_data   computation source data voxel array
    /// \param[in] valid_src_bbox source data bounding box. This might
    /// be a partial bbox for local computing.
    /// \param[in] whole_src_roi_bbox whole source region of interest
    /// bounding box
    void apply_filter(nv::index_common::VoxelType * p_dst_data,
                      mi::math::Bbox< mi::Sint64, 3 > const & raw_dst_bbox,
                      nv::index_common::VoxelType * p_src_data,
                      mi::math::Bbox< mi::Sint64, 3 > const & valid_src_bbox,
                      mi::math::Bbox< mi::Sint64, 3 > const & whole_src_roi_bbox) const;

private:
    /// initiaize filter (kernal) according to the filter_type.
    /// \param[in] filter_type filter type
    void initialize_filter(mi::Sint32 filter_type);

    /// filter one voxel by a 2d kernel
    ///
    /// \param[in] p_src_data input source data
    /// \param[in] valid_src_bbox source bounding box
    /// \param[in] center_ijk  ijk coordinates of computing voxel
    /// \return filter result
    mi::Float32 filter_2d_one(
        nv::index_common::VoxelType * p_src_data,
        mi::math::Bbox< mi::Sint64, 3 >   const & valid_src_bbox,
        mi::math::Vector< mi::Sint64, 3 > const & center_ijk) const;

    /// filter one voxel by a 3d kernel
    ///
    /// \param[in] p_src_data input source data
    /// \param[in] valid_src_bbox source bounding box
    /// \param[in] center_ijk  ijk coordinates of computing voxel
    /// \return filter result
    mi::Float32 filter_3d_one(
        nv::index_common::VoxelType * p_src_data,
        mi::math::Bbox< mi::Sint64, 3 >   const & valid_src_bbox,
        mi::math::Vector< mi::Sint64, 3 > const & center_ijk) const;

    /// DIP median filter
    ///
    /// \param[in] p_src_data input source data
    /// \param[in] valid_src_bbox source bounding box
    /// \param[in] center_ijk ijk coordinates of computing voxel
    /// \return DIP median filter result
    mi::Float32 dip_median_3d_one(
        nv::index_common::VoxelType * p_src_data,
        mi::math::Bbox< mi::Sint64, 3 >   const & valid_src_bbox,
        mi::math::Vector< mi::Sint64, 3 > const & center_ijk) const;


private:
    /// get 2d 3x3 kernel index
    /// \param[in] x
    /// \param[in] y
    /// \return index of 3x3 kernel
    static mi::Sint32 get_2d_33_kernel_index(mi::Sint32 x, mi::Sint32 y);

    /// get 3d 3x3x3 kernel index
    /// \param[in] x
    /// \param[in] y
    /// \param[in] z
    /// \return index of 3x3x3 kernel
    static mi::Sint32 get_3d_33_kernel_index(mi::Sint32 x, mi::Sint32 y, mi::Sint32 z);

private:
    /// default constructor. prohivit this.
    Volume_data_filter();

private:
    // filter type
    mi::Sint32                           m_filter_type;
    // 2d kernel
    mi::Float32                          m_kernel_2d[VDF_kernel_2d_size];
    // 3d kernel
    mi::Float32                          m_kernel_3d[VDF_kernel_3d_size];
};


#endif // #ifndef GEOSPATIAL_BIN_GEOSPATIAL_STREAM_VIEWER_VOLUME_DATA_FILTER_H
