/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "volume_data_filter.h"

#include "common/common_utility.h"
#include "common/encode_voxel.h"
#include "common/forwarding_logger.h"

#include "utilities.h"


//----------------------------------------------------------------------
std::string Volume_data_filter::get_filter_name(mi::Sint32 filter_type)
{
    char const * const p_name[VDFT_Count + 1] = {
        "gaussian_ij",
        "gaussian_jk",
        "gaussian_ki",
        "sobel_ij",
        "sobel_jk",
        "sobel_ki",
        "dip_median",
        NULL,
    };
    assert(VDFT_Count  == 7);
    assert(is_valid_filter_type(filter_type));

    return std::string(p_name[filter_type]);
}

//----------------------------------------------------------------------
mi::Sint32  Volume_data_filter::get_filter_type(std::string const & filter_name)
{
    if(filter_name == "gaussian_ij"){
        return VDFT_Gaussian_IJ;
    }else if(filter_name == "gaussian_jk"){
        return VDFT_Gaussian_JK;
    }else if(filter_name == "gaussian_ki"){
        return VDFT_Gaussian_KI;
    }else if(filter_name == "sobel_ij"){
        return VDFT_Sobel_IJ;
    }else if(filter_name == "sobel_jk"){
        return VDFT_Sobel_JK;
    }else if(filter_name == "sobel_ki"){
        return VDFT_Sobel_KI;
    }else if(filter_name == "dip_median"){
        return VDFT_DIP_MEDIAN;
    }else{
        return VDFT_Count;
    }
}

//----------------------------------------------------------------------
bool Volume_data_filter::is_valid_filter_type(mi::Sint32 filter_type)
{
    if((filter_type < 0) || (filter_type >= VDFT_Count)){
        return false;
    }
    return true;
}

//----------------------------------------------------------------------
bool Volume_data_filter::is_valid_subcube_size(
    mi::math::Vector_struct<mi::Uint32, 3> const & subcube_size)
{
    if((subcube_size.x >= VDF_extra_comp_width) &&
       (subcube_size.y >= VDF_extra_comp_width) &&
       (subcube_size.z >= VDF_extra_comp_width))
    {
        return true;
    }
    return false;
}

//----------------------------------------------------------------------
mi::math::Bbox< mi::Sint64, 3 >
Volume_data_filter::get_valid_compute_source_bbox(
    mi::math::Bbox< mi::Sint64, 3 > const & output_bbox,
    mi::math::Bbox< mi::Sint64, 3 > const & whole_src_roi_bbox
    )
{
    mi::math::Vector< mi::Sint64, 3 > const whole_src_size_vec3 =
        whole_src_roi_bbox.max - whole_src_roi_bbox.min;
    mi::math::Bbox< mi::Sint64, 3 > const whole_src_size_bbox(
        mi::math::Vector< mi::Sint64, 3 >(0, 0, 0), whole_src_size_vec3);

    mi::math::Bbox< mi::Sint64, 3 > valid_src_roi_bbox(
        whole_src_roi_bbox.min + output_bbox.min,
        whole_src_roi_bbox.min + output_bbox.max);

    // add extra slices for computation. Depends on the filter.
    // As an example, we use always 3x3x3 kernel.
    mi::math::Vector< mi::Sint64, 3 > ONEVec3(1, 1, 1);
    valid_src_roi_bbox.min -= (VDF_extra_comp_width * ONEVec3);
    valid_src_roi_bbox.max += (VDF_extra_comp_width * ONEVec3);

    // clip to the valid range
    valid_src_roi_bbox = nv::index_common::bbox_and(valid_src_roi_bbox, whole_src_roi_bbox);

    return valid_src_roi_bbox;
}

//----------------------------------------------------------------------
Volume_data_filter::Volume_data_filter()
    :
    m_filter_type(VDFT_Gaussian_IJ)
    // m_kernel_2d
    // m_kernel_3d
{
#ifdef DEBUG
    this->initialize_filter(VDFT_Gaussian_IJ);
#endif
}

//----------------------------------------------------------------------
Volume_data_filter::Volume_data_filter(
    mi::Sint32 filter_type)
    :
    m_filter_type(filter_type)
{
    this->initialize_filter(m_filter_type);

    assert((sizeof(m_kernel_2d)/sizeof(mi::Float32)) == VDF_kernel_2d_size);
    assert((sizeof(m_kernel_3d)/sizeof(mi::Float32)) == VDF_kernel_3d_size);
}

//----------------------------------------------------------------------
Volume_data_filter::~Volume_data_filter()
{
    // empty
}

//----------------------------------------------------------------------
void Volume_data_filter::initialize_filter(mi::Sint32 filter_type)
{
    for(mi::Sint32 i = 0; i < VDF_kernel_2d_size; ++i){
        m_kernel_2d[i] = 0.0f;
    }
    for(mi::Sint32 i = 0; i < VDF_kernel_3d_size; ++i){
        m_kernel_3d[i] = 0.0f;
    }

    if((filter_type == VDFT_Gaussian_IJ) ||
       (filter_type == VDFT_Gaussian_JK) ||
       (filter_type == VDFT_Gaussian_KI)){
        mi::Float32 const gauss[VDF_kernel_2d_size] = {  0.0113f, 0.0838f, 0.0113f,
                                                          0.0838f, 0.6193f, 0.0838f,
                                                          0.0113f, 0.0838f, 0.0113f, };
        for(mi::Sint32 i = 0; i < VDF_kernel_2d_size; ++i){
            m_kernel_2d[i] = gauss[i];
        }
    }
    else if((filter_type == VDFT_Sobel_IJ) ||
            (filter_type == VDFT_Sobel_JK) ||
            (filter_type == VDFT_Sobel_KI)){
        mi::Float32 const sobel[VDF_kernel_2d_size] = { -1.0f, 0.0f, 1.0f,
                                                         -2.0f, 0.0f, 2.0f,
                                                         -1.0f, 0.0f, 1.0f, };
        for(mi::Sint32 i = 0; i < VDF_kernel_2d_size; ++i){
            m_kernel_2d[i] = sobel[i];
        }
    }
    else if(filter_type == VDFT_DIP_MEDIAN){
        // no preparation
    }
    else{
        ERROR_LOG << "Volume_data_filter::initialize_filter: "
                  << "Unknown filter type.";
    }
}

//----------------------------------------------------------------------
static mi::math::Vector< mi::Sint64, 3 > const IJK_OFFSET_TAB[4] = {
    mi::math::Vector< mi::Sint64, 3 >(1, 1, 0), // IJ
    mi::math::Vector< mi::Sint64, 3 >(0, 1, 1), // JK
    mi::math::Vector< mi::Sint64, 3 >(1, 0, 1), // KI
    mi::math::Vector< mi::Sint64, 3 >(1, 1, 1), // IJK
};

// IJK offset index table: filter type -> index
static mi::Sint32 const OFFSET_IDX[Volume_data_filter::VDFT_Count] = {
    0, 1, 2, 0, 1, 2, 3,
};

/// show filter operation progress
static void show_progress_debug(mi::math::Bbox< mi::Sint64, 3 >   const & bbox,
                                mi::math::Vector< mi::Sint64, 3 > const & ijk,
                                std::string const & mes)
{
    mi::Sint32 const dim_idx  = 0;
    mi::Sint32 const rep_freq = 10;
    mi::Sint64  const isize = bbox.max[dim_idx] - bbox.min[dim_idx];
    mi::Float32 const fsize = static_cast< mi::Float32 >(isize);
    if((isize <= rep_freq) ||
       ((ijk[dim_idx] - bbox.min[dim_idx]) % (isize / rep_freq)) == 0)
    {
        DEBUG_LOG << mes << floor((static_cast< mi::Float32 >(ijk[dim_idx] - bbox.min[dim_idx])
                                   / fsize) * 100.0f)
                  << "% done";
    }
}

//----------------------------------------------------------------------
void Volume_data_filter::apply_filter(
    nv::index_common::VoxelType * p_dst_data,
    mi::math::Bbox< mi::Sint64, 3 > const & raw_dst_bbox,
    nv::index_common::VoxelType * p_src_data,
    mi::math::Bbox< mi::Sint64, 3 > const & valid_src_bbox,
    mi::math::Bbox< mi::Sint64, 3 > const & whole_src_roi_bbox
    ) const
{
    assert(is_valid_filter_type(m_filter_type));
    mi::math::Vector< mi::Sint64, 3 > const ijk_offset =
        IJK_OFFSET_TAB[OFFSET_IDX[m_filter_type]];

    // get valid filter-able bbox from the source
    mi::math::Bbox< mi::Sint64, 3 > conv_bbox = valid_src_bbox;
    conv_bbox.min += ijk_offset;
    conv_bbox.max -= ijk_offset;

    // result bbox also should be considered.
    conv_bbox = nv::index_common::bbox_and(conv_bbox, raw_dst_bbox);
    //     DEBUG_LOG << "filter = "  << Volume_data_filter::get_filter_name(m_filter_type)
    //               << ", rowdst:"  << raw_dst_bbox
    //               << ", update: " << conv_bbox
    //               << ", src: "    << valid_src_bbox;

    mi::math::Vector< mi::Sint64, 3 > ijk = valid_src_bbox.min;

    if(m_filter_type == VDFT_DIP_MEDIAN){
        // 3d filter
        for(ijk[0] = conv_bbox.min[0]; ijk[0] < conv_bbox.max[0]; ++(ijk[0])){
            show_progress_debug(conv_bbox, ijk, "filter3d ");
            for(ijk[1] = conv_bbox.min[1]; ijk[1] < conv_bbox.max[1]; ++(ijk[1])){
                for(ijk[2] = conv_bbox.min[2]; ijk[2] < conv_bbox.max[2]; ++(ijk[2])){
                    mi::Float32 const res =
                        this->dip_median_3d_one(p_src_data, valid_src_bbox, ijk);
                    mi::Sint64 const dst_idx = nv::index_common::get_volume_index(ijk, raw_dst_bbox);
                    p_dst_data[dst_idx] = static_cast< nv::index_common::VoxelType >(res);
                }
            }
        }
    }
    else{
        // 2d filter
        for(ijk[0] = conv_bbox.min[0]; ijk[0] < conv_bbox.max[0]; ++(ijk[0])){
            show_progress_debug(conv_bbox, ijk, "filter2d ");
            for(ijk[1] = conv_bbox.min[1]; ijk[1] < conv_bbox.max[1]; ++(ijk[1])){
                for(ijk[2] = conv_bbox.min[2]; ijk[2] < conv_bbox.max[2]; ++(ijk[2])){
                    mi::Float32 const res =
                        this->filter_2d_one(p_src_data, valid_src_bbox, ijk);
                    mi::Sint64 const dst_idx = nv::index_common::get_volume_index(ijk, raw_dst_bbox);
                    p_dst_data[dst_idx] = static_cast< nv::index_common::VoxelType >(res);
                }
            }
        }
    }

    mi::math::Bbox< mi::Sint64, 3 > const valid_bbox =
        nv::index_common::bbox_and(whole_src_roi_bbox, conv_bbox);

    // fill the border. expanding border until max(valid_bbox +-1, raw_dst_bbox)
    mi::math::Vector< mi::Sint64, 3 > const ONEVec3(1, 1, 1);
    mi::math::Bbox< mi::Sint64, 3 > const expand_valid_bbox(valid_bbox.min - ONEVec3,
                                                            valid_bbox.max + ONEVec3);
    mi::math::Bbox< mi::Sint64, 3 > const expand_dst_bbox =
        nv::index_common::bbox_and(expand_valid_bbox, raw_dst_bbox);
    mi::math::Bbox< mi::Sint64, 3 > expanding_bbox = conv_bbox;
    assert(bbox_inclusive_contain(raw_dst_bbox, expand_dst_bbox));

    while(expanding_bbox != expand_dst_bbox){
        std::stringstream sstr;
        sstr << "expanding border: " << expanding_bbox << ", dst = " << expand_dst_bbox;
        expanding_bbox = nv::index_common::expand_border_by_copy(p_dst_data, raw_dst_bbox, expanding_bbox);
        sstr << " -> " << expanding_bbox;
        // DEBUG_LOG << sstr.str();
    }
}

//----------------------------------------------------------------------
mi::Float32 Volume_data_filter::filter_2d_one(
    nv::index_common::VoxelType * p_src_data,
    mi::math::Bbox< mi::Sint64, 3 >   const & valid_src_bbox,
    mi::math::Vector< mi::Sint64, 3 > const & center_ijk) const
{
    mi::Float32 res = 0.0f;

    mi::math::Vector< mi::Sint64, 3 > uvw(0, 0, 0);
    for(mi::Sint64 u = -1; u < 2; ++u){
        for(mi::Sint64 v = -1; v < 2; ++v){
            uvw[ OFFSET_IDX[m_filter_type]         ] = u;
            uvw[(OFFSET_IDX[m_filter_type] + 1) % 3] = v;
            mi::math::Vector< mi::Sint64, 3 > const v_ijk = center_ijk + uvw;
            mi::Sint64 const src_idx =  nv::index_common::get_volume_index(v_ijk, valid_src_bbox);
            mi::Float32 const sdat = static_cast< mi::Float32 >(p_src_data[src_idx]);
            mi::Float32 const conv = sdat * m_kernel_2d[get_2d_33_kernel_index(static_cast<mi::Sint32>(u),
                                                                               static_cast<mi::Sint32>(v))];
            res += conv;
        }
    }
    return fabs(res);
}

//----------------------------------------------------------------------
mi::Float32 Volume_data_filter::filter_3d_one(
    nv::index_common::VoxelType * p_src_data,
    mi::math::Bbox< mi::Sint64, 3 >   const & valid_src_bbox,
    mi::math::Vector< mi::Sint64, 3 > const & center_ijk) const
{
    mi::Float32 res = 0.0f;

    mi::math::Vector< mi::Sint64, 3 > uvw(0, 0, 0);
    for(mi::Sint64 u = -1; u < 2; ++u){
        for(mi::Sint64 v = -1; v < 2; ++v){
            for(mi::Sint64 w = -1; w < 2; ++w){
                // Not implemented any 3D filter
                // uvw = ;
                // mi::math::Vector< mi::Sint64, 3 > const v_ijk = center_ijk + uvw;
                // mi::Sint64 const src_idx =  get_volume_index(v_ijk, valid_src_bbox);
                // mi::Float32 const sdat = static_cast< mi::Float32 >(p_src_data[src_idx]);
                // mi::Float32 const conv = sdat * m_kernel_3d[get_3d_33_kernel_index(u, v, w)];
                // res += conv;
            }
        }
    }
    return fabs(res);
}

//----------------------------------------------------------------------
mi::Float32 Volume_data_filter::dip_median_3d_one(
    nv::index_common::VoxelType * p_src_data,
    mi::math::Bbox< mi::Sint64, 3 >   const & valid_src_bbox,
    mi::math::Vector< mi::Sint64, 3 > const & center_ijk) const
{
    const size_t bufsize    = 27;
    const size_t median_idx = 13;
    std::vector< mi::Float32 > float_buffer(bufsize, 0.0f);

    mi::math::Vector< mi::Sint64, 3 > uvw(0, 0, 0);
    size_t fb_idx = 0;
    for(uvw[0] = -1; uvw[0] < 2; ++uvw[0]){
        for(uvw[1] = -1; uvw[1] < 2; ++uvw[1]){
            for(uvw[2] = -1; uvw[2] < 2; ++uvw[2]){
                mi::math::Vector< mi::Sint64, 3 > const v_ijk = center_ijk + uvw;
                mi::Sint64 const src_idx =  nv::index_common::get_volume_index(v_ijk, valid_src_bbox);
                assert(fb_idx < bufsize);
                float_buffer[fb_idx] = static_cast< mi::Float32 >(p_src_data[src_idx]);
                ++fb_idx;
            }
        }
    }
    std::nth_element(float_buffer.begin(),
                     float_buffer.begin() + median_idx,
                     float_buffer.end());
    mi::Float32 const res = float_buffer[median_idx];
    return res;
}

//----------------------------------------------------------------------
mi::Sint32 Volume_data_filter::get_2d_33_kernel_index(mi::Sint32 x, mi::Sint32 y)
{
    assert(-1 <= x);
    assert(x  <= 1);
    assert(-1 <= y);
    assert(y  <= 1);

    mi::Sint32 const idx = (x + 1) + (y + 1) * 3;
    assert(0  <= idx);
    assert(idx<  9);

    return idx;
}

//----------------------------------------------------------------------
mi::Sint32 Volume_data_filter::get_3d_33_kernel_index(
    mi::Sint32 x, mi::Sint32 y, mi::Sint32 z)
{
    assert(-1 <= x);
    assert(x  <= 1);
    assert(-1 <= y);
    assert(y  <= 1);
    assert(-1 <= z);
    assert(z  <= 1);

    mi::Sint32 const idx = (x + 1) + ((y + 1) + ((z + 1) * 3) * 3);
    assert(0   <= idx);
    assert(idx <  27);

    return idx;
}

//----------------------------------------------------------------------
