/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief heightfield data exporter example

#ifndef BIN_GEOSPATIAL_STREAM_VIEWER_HEIGHTFIELD_DATA_EXPORT_H
#define BIN_GEOSPATIAL_STREAM_VIEWER_HEIGHTFIELD_DATA_EXPORT_H

#include <mi/dice.h>

#include <nv/index/idistributed_data_locality.h>
#include <nv/index/idistributed_data_access.h>
#include <nv/index/iindex.h>

#include <cstdio>
#include <string>
#include <vector>

//----------------------------------------------------------------------

/// Heightfield data exporter
class Heightfield_data_exporter
{
public:
    /// constructor
    ///
    /// \param[in] heightfield_tag heightfield scene element tag
    /// \param[in] ij_bbox exporting value range (value range = [min, max))
    /// \param[in] data_distribution
    /// \param[in] data_access_factory
    Heightfield_data_exporter(
        const mi::neuraylib::Tag&                                            heightfield_tag,
        const mi::math::Bbox<mi::Sint32, 2>&                                 ij_bbox,
        mi::base::Handle<const nv::index::IData_distribution>&               data_distribution,
        mi::base::Handle<const nv::index::IDistributed_data_access_factory>& data_access_factory);

    /// Export the heightfield data to a local file, whereas local means
    /// a file on a local cluster machine.
    /// \param[in] filename The name of the local file.
    /// \param[in] dice_transaction
    /// \return true when succeeded
    bool export_to_local_file(
        const std::string & filename,        
        mi::neuraylib::IDice_transaction * dice_transaction) const;

private:
    /// get exporting subregion bounding box
    void get_exporting_bbox(
        std::vector< mi::math::Bbox<mi::Sint32, 2> >& export_bbox_vec,
        mi::neuraylib::IDice_transaction *                   dice_transaction) const;
    /// check the exporting bbox is in (inclusive) the heightfield bbox
    bool is_exporting_bbox_in_heightfield(
        mi::neuraylib::IDice_transaction * dice_transaction) const;
    /// serialize heightfield data to a file
    ///
    /// \param[in] outfile             output file pointer.
    /// \param[in] header_size         size of header (bytes)
    /// \param[in] end_of_header_fpos  end of header file position.
    /// \param[in] heightfield_data_access data access for the exporting heightfield
    /// \return true when succeeded
    bool serialized_write_to_file(
        FILE* outfile,
        const mi::Size   header_size,
        const mi::Sint64 end_of_header_fpos,
        mi::base::Handle<nv::index::IRegular_heightfield_data_access>& heightfield_data_access) const;

private:
    mi::neuraylib::Tag                                                   m_heightfield_tag;
    mi::math::Bbox<mi::Sint32, 2>                                        m_ij_bbox;
    mi::base::Handle<const nv::index::IData_distribution>                m_data_distribution;
    mi::base::Handle<const nv::index::IDistributed_data_access_factory>  m_data_access_factory;
};

//----------------------------------------------------------------------

/// export a heightfield data.
///
/// \param[in] heightfield_tag heightfield tag to be exported
/// \param[in] ij_bbox         exporting IJ bounding box
/// \param[in] export_fname    export file name
/// \param[in] session_tag     session tag
/// \param[in] dice_transaction     transaction database scope
/// \return true when success
bool export_heightfield_data(
    const mi::neuraylib::Tag &            heightfield_tag,
    const mi::math::Bbox<mi::Sint32, 2> & ij_bbox,
    const std::string&                    export_fname,
    const mi::neuraylib::Tag&             session_tag,
    mi::neuraylib::IDice_transaction *    dice_transaction);


//----------------------------------------------------------------------
#endif // BIN_GEOSPATIAL_STREAM_VIEWER_HEIGHTFIELD_DATA_EXPORT_H
