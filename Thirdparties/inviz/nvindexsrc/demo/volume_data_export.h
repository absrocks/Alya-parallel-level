/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief external volume data exporter example

#ifndef GEOSPATIAL_BIN_GEOSPATIAL_STREAM_VIEWER_VOLUME_DATA_EXPORT_H
#define GEOSPATIAL_BIN_GEOSPATIAL_STREAM_VIEWER_VOLUME_DATA_EXPORT_H

#include <mi/base/types.h>
#include <mi/dice.h>
#include <mi/base/handle.h>
#include <nv/index/idistributed_data_locality.h>
#include <nv/index/idistributed_data_access.h>
#include <nv/index/iindex.h>

#include <cstdio>
#include <string>
#include <vector>

#include "common/encode_voxel.h"


/// volume data exporter
class Volume_data_exporter
{
public:
    /// constructor
    ///
    /// \param[in] volume_tag volume scene element tag
    /// \param[in] ijk_bbox exporting bounding box (bbox.max == volume size value)
    /// \param[in] data_distribution subsurface data distribution
    /// \param[in] data_access_factory data access factory
    Volume_data_exporter(
        const mi::neuraylib::Tag&                                            volume_tag,
        const mi::math::Bbox<mi::Sint32, 3>&                                 ijk_bbox,
        mi::base::Handle<const nv::index::IData_distribution>&               data_distribution,
        mi::base::Handle<const nv::index::IDistributed_data_access_factory>& data_access_factory);

    /// Export the volume data to a local file, whereas local means
    /// a file on a local cluster machine.
    /// \param[in] fpath The name of the local file.
    /// \param[in] dice_transaction transaction
    /// \return true when succeeded
    bool export_to_local_file(
        const std::string & fpath,
        mi::neuraylib::IDice_transaction*  dice_transaction) const;

private:
    /// export to local file with n slices at once.
    /// Assume data can be fit in the main memory
    ///
    /// \param[in] dice_transaction dice transaction
    /// \param[in] filename         exported volume data filename
    /// \param[in] slices_transfer_once_count number of slices to
    /// transfer at once from remote to host.
    /// \return true when succeeded
    bool export_to_local_file_with_n_slices(
        const std::string & fpath,
        mi::Sint32 slices_transfer_once_count,
        mi::neuraylib::IDice_transaction * dice_transaction) const;

    /// check the exporting region of interest and query bounding box
    ///
    bool is_valid_export_region(
        mi::math::Bbox<mi::Sint32, 3> const & volume_bbox,
        std::vector< mi::math::Bbox<mi::Sint32, 3> > const & export_all_bbox_vec)
        const;


    /// export partial volume.
    ///
    /// To avoid all the data copy to local, we peocess the volume partially.
    ///
    /// \param[in] p_export_fp      exporting FILE*.
    /// \param[in] whole_volume_bbox total volume bounding box
    /// \param[in] bbox_partial_vec partial exporting volume bounding box.
    /// \param[in] number_of_slices_transfer number of slices to be
    /// transfered at once. dice's current limitation is receiving
    /// data size. Keep the receiving size less than some amount. (1MB
    /// shuold be good... If the octree leaf size if 509^3, 4 is a good number.)
    /// \param[in] dice_transaction dice transaction
    /// \return true when succeeded
    bool export_partial_volume(
        FILE * p_export_fp,
        mi::math::Bbox<mi::Sint32, 3> const & whole_volume_bbox,
        std::vector< mi::math::Bbox<mi::Sint32, 3> > const & bbox_partial_vec,
        mi::Sint32 number_of_slices_transfer,
        mi::neuraylib::IDice_transaction * dice_transaction) const;

    /// export slices.
    ///
    /// \param[in] dice_transaction dice transaction
    /// \param[in] p_export_fp      exporting FILE*.
    /// \param[in] whole_volume_bbox total volume bounding box
    /// \param[in] slice_bbox_vec   slice bounding box vector
    /// \return true when succeeded
    bool export_slices(
        mi::neuraylib::IDice_transaction * dice_transaction,
        FILE * p_export_fp,
        mi::math::Bbox<mi::Sint32, 3> const & whole_volume_bbox,
        std::vector< mi::math::Bbox<mi::Sint32, 3> > const & slice_bbox_vec) const;

    /// get exporting bounding boxes (bounding box ... [min, max])
    ///
    /// \param[in]  dice_transaction dice transaction
    /// \param[out] export_bounding_box_vec exporting bounding box vector
    void get_exporting_bbox(
        mi::neuraylib::IDice_transaction*  dice_transaction,
        std::vector< mi::math::Bbox<mi::Sint32, 3> > & export_bounding_box_vec) const;

    /// write one bounding box (of subvolume) to the file
    ///
    /// \param[in] p_ofp                output file pointer
    /// \param[in] subvolume_bbox ijk   bounding box of the subvolume
    /// \param[in] whole_volume_bbox    whole volume bbox (this is
    /// actually not needed, but useful for check the correctness.)
    /// \param[in] p_data amplitude raw data to be written (the size
    /// is the size of subvolume)
    /// \return true when succeeded
    bool write_one_bbox_to_file(
        FILE * p_outfile,
        mi::math::Bbox<mi::Sint32, 3> const & subvolume_bbox,
        mi::math::Bbox<mi::Sint32, 3> const & whole_volume_bbox,
        const nv::index_common::VoxelType*    p_data) const;

    /// map subvolume bounding box ijk coordinate to data memory index
    ///
    /// subvolume bbox ijk coordinate shares the region of interest
    /// bounding box coordinate, but the value of subvolume bbox is
    /// equal or smaller than roi bbox.
    ///
    /// \param[in] i i of subvolume ijk
    /// \param[in] j j of subvolume ijk
    /// \param[in] k k of subvolume ijk
    /// \param[in] subvol_bbox exporting subvolume bounding box [min, max]
    /// \return linear index of the subvolume ijk.
    mi::Sint64 subvol_bbox_ijk_to_memory_index(
        mi::Sint32 i,
        mi::Sint32 j,
        mi::Sint32 k,
        mi::math::Bbox<mi::Sint32, 3> const & subvol_bbox) const;

    /// region of interest bounding box ijk coordinate to export file
    /// position index.
    ///
    /// ROI coordinates shareds subvolume coordinates. But, ROI
    /// coordinates values should includes subvolume coordinates
    /// values.
    ///
    /// \param[in] i i of ijk (ROI coordinate)
    /// \param[in] j j of ijk (ROI coordinate)
    /// \param[in] k k of ijk (ROI coordinate)
    /// \param[in] roi_bbox exporting region of interest bbox [min,max]
    /// \param[in] subvol_bbox volume bonding box [min, max)
    /// \return file position
    mi::Sint64 ROI_bbox_ijk_to_file_position_index(
        mi::Sint32 i,
        mi::Sint32 j,
        mi::Sint32 k,
        mi::math::Bbox<mi::Sint32, 3> const & roi_bbox,
        mi::math::Bbox<mi::Sint32, 3> const & subvol_bbox) const;


private:
    mi::neuraylib::Tag                                                   m_volume_tag;
    mi::math::Bbox<mi::Sint32, 3>                                        m_ijk_bbox;
    mi::base::Handle<const nv::index::IData_distribution>                m_data_distribution;
    mi::base::Handle<const nv::index::IDistributed_data_access_factory>  m_data_access_factory;
};

//----------------------------------------------------------------------
/// export a volume data.
///
/// This needs a transaction, please care the transaction nesting (not
/// possible in dice) when you use this.
///
/// \param[in] volume_tag           volume tag to be exported
/// \param[in] export_bbox          exporting volume ijk bounding box
/// \param[in] export_fname         export file name
/// \param[in] session_tag          session tag
/// \param[in] dice_transaction     transaction database scope
/// \return true when success
bool export_volume_data(
    const mi::neuraylib::Tag &            volume_tag,
    const mi::math::Bbox<mi::Sint32, 3> & ijk_bbox,
    const std::string&                    export_fname,
    const mi::neuraylib::Tag&             session_tag,
    mi::neuraylib::IDice_transaction *    dice_transaction);

//----------------------------------------------------------------------
#endif // GEOSPATIAL_BIN_GEOSPATIAL_STREAM_VIEWER_VOLUME_DATA_EXPORT_H
