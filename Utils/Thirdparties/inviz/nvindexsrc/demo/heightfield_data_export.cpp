/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "heightfield_data_export.h"

#include <nv/index/iregular_heightfield.h>
#include <nv/index/isession.h>
#include <nv/index/iscene.h>

#include <cassert>

#include "heightfield_data_retrieval_appjob.h"
#include "nvindex_appdata.h"
#include "scene_utility.h"

#include "common/common_utility.h"
#include "common/forwarding_logger.h"
#include "common/normal_encoder.h"
#include "common/type_conversion_utility.h"

//----------------------------------------------------------------------
Heightfield_data_exporter::Heightfield_data_exporter(
    const mi::neuraylib::Tag&                                            heightfield_tag,
    const mi::math::Bbox<mi::Sint32, 2>&                                 ij_bbox,
    mi::base::Handle<const nv::index::IData_distribution>&               data_distribution,
    mi::base::Handle<const nv::index::IDistributed_data_access_factory>& data_access_factory)
    :
    m_heightfield_tag(heightfield_tag),
    m_ij_bbox(ij_bbox),
    m_data_distribution(data_distribution),
    m_data_access_factory(data_access_factory)
{
    // empty
}

//----------------------------------------------------------------------
void Heightfield_data_exporter::get_exporting_bbox(
    std::vector< mi::math::Bbox<mi::Sint32, 2> >& export_bbox_vec,
    mi::neuraylib::IDice_transaction*             dice_transaction) const
{
    assert(dice_transaction != 0);
    export_bbox_vec.clear();

    // Figure out data locality.
    INFO_LOG << "Data locality call for bbox: " << m_ij_bbox;
    mi::base::Handle<nv::index::IRegular_heightfield_data_locality> data_locality(
        m_data_distribution->retrieve_data_locality(
            m_heightfield_tag,
            nv::index_common::convert_bbox_type<mi::Uint32, mi::Sint32, 2>(m_ij_bbox),
            dice_transaction));

    const mi::math::Bbox<mi::Sint32, 2> ij_intersect_bbox = m_ij_bbox;
    const mi::Uint32 nb_cluster_nodes = data_locality->get_nb_cluster_nodes();

    for (mi::Uint32 i = 0; i < nb_cluster_nodes; ++i)
    {
        const mi::Uint32 cluster_node_id = data_locality->get_cluster_node(i);
        INFO_LOG << "Cluster node " << i << " has id:" << cluster_node_id;
        const mi::Uint32 nb_bbox = static_cast<mi::Uint32>(data_locality->get_nb_bounding_box(cluster_node_id));

        for (mi::Uint32 bidx = 0; bidx < nb_bbox; ++bidx)
        {
            // EXAMPLE: simply pick the ones the you feel make sense for a clever export of heightfield.
            // 'Clever' should take into account:
            //  - the main memory restriction
            //  - the J-first and I-last export scheme
            //  - the number of data access calls (possibly)

            // The following is *not* very clever ;) - especially for large dataset sizes .....
            const mi::math::Bbox<mi::Sint32, 3> patch_3D_bbox(data_locality->get_bounding_box(cluster_node_id, bidx));
            mi::math::Bbox<mi::Sint32, 2> bbox(patch_3D_bbox.min.x, patch_3D_bbox.min.y, 
                                               patch_3D_bbox.max.x, patch_3D_bbox.max.y);

            // bbox intersection
            if (bbox.intersects(ij_intersect_bbox))
            {
                INFO_LOG << "get_exporting_bbox: bounding box = " << bbox;
                export_bbox_vec.push_back(bbox);
            }
        }
    }
}

//----------------------------------------------------------------------
bool Heightfield_data_exporter::is_exporting_bbox_in_heightfield(
    mi::neuraylib::IDice_transaction * dice_transaction) const
{
    assert(dice_transaction != 0);
    assert(m_heightfield_tag.is_valid());

    mi::base::Handle<const nv::index::IRegular_heightfield> heightfield(
        dice_transaction->access<const nv::index::IRegular_heightfield>(m_heightfield_tag));
    assert(heightfield.is_valid_interface());

    mi::math::Bbox<mi::Float32, 3> hf_bbox_f3 = heightfield->get_IJK_bounding_box();

    mi::math::Bbox<mi::Sint32, 2> hf_bbox_u2(static_cast<mi::Sint32>(floor(hf_bbox_f3.min.x)),
                                             static_cast<mi::Sint32>(floor(hf_bbox_f3.min.y)),
                                             static_cast<mi::Sint32>(ceil( hf_bbox_f3.max.x)),
                                             static_cast<mi::Sint32>(ceil( hf_bbox_f3.max.y)));
    if (!hf_bbox_u2.contains(m_ij_bbox.min))
    {
        ERROR_LOG << "bbox: " << m_ij_bbox << " min is not inclusive inside the heightfield bbox: "
                  << hf_bbox_u2;
        return false;
    }

    if (!hf_bbox_u2.contains(m_ij_bbox.max))
    {
        ERROR_LOG << "bbox: " << m_ij_bbox << " max is not inclusive inside the heightfield bbox: "
                  << hf_bbox_u2;
        return false;
    }

    return true;
}

//----------------------------------------------------------------------
bool Heightfield_data_exporter::serialized_write_to_file(
    FILE* outfile,
    const mi::Size   header_size,
    const mi::Sint64 end_of_header_fpos,
    mi::base::Handle<nv::index::IRegular_heightfield_data_access>& heightfield_data_access) const
{
    assert(outfile != 0);
    assert(heightfield_data_access.is_valid_interface());

    const mi::math::Bbox<mi::Sint32, 2> patch_bbox = 
        nv::index_common::convert_bbox_type<mi::Sint32, mi::Uint32, 2>(
            heightfield_data_access->get_patch_bounding_box());

    INFO_LOG << "Writing patch bbox: " << patch_bbox << " of " << m_ij_bbox;

    const mi::math::Vector<mi::Sint32, 2> whole_ij_bbox_dim = m_ij_bbox.max - m_ij_bbox.min;
    assert(whole_ij_bbox_dim.x > 0);
    assert(whole_ij_bbox_dim.y > 0);
    INFO_LOG << "Exporting whole dimension: " << whole_ij_bbox_dim << ", size: " << (whole_ij_bbox_dim.x * whole_ij_bbox_dim.y);    

    const mi::math::Vector<mi::Sint32,2> ij_patch_dim = patch_bbox.max - patch_bbox.min;
    const mi::Sint32 ij_patch_size = ij_patch_dim.x * ij_patch_dim.y;
    INFO_LOG << "Exporting patch: " << patch_bbox << ", size: " << ij_patch_size;   

    // Retrieved height values
    mi::Float32* height_values = heightfield_data_access->get_elevation_values();
    mi::math::Vector_struct<mi::Float32, 3> * normal_vector_values = heightfield_data_access->get_normal_values();

    const mi::Uint32 x_offset = 0; // patch_bbox.min.x; // 0; // FIXME for range export

    const mi::Uint32 patch_orig_x_in_file = patch_bbox.min.x - m_ij_bbox.min.x;
    const mi::Uint32 patch_orig_y_in_file = patch_bbox.min.y - m_ij_bbox.min.y;

    INFO_LOG << "patch orig: " << patch_orig_x_in_file << "," << patch_orig_y_in_file;

    // keep the current file position
    const mi::Sint64 height_begin_fpos = end_of_header_fpos;
    const mi::Sint64 normal_begin_fpos = height_begin_fpos + sizeof(mi::Float32) * (whole_ij_bbox_dim.x * whole_ij_bbox_dim.y);

    // Number of elements to write in each loop
    // const mi::Uint32 nb_elements = patch_range_i;
    mi::Size write_byte_count = 0;

    const mi::Sint32 nb_elements = ij_patch_dim.x;
    for (mi::Sint32 j = 0; j < ij_patch_dim.y; ++j)
    {
        const mi::Uint64 patch_idx = j * ij_patch_dim.x + x_offset;

        mi::Float32* heights = &height_values[patch_idx];
        mi::math::Vector_struct<mi::Float32, 3>* normals = &(normal_vector_values[patch_idx]);

        const mi::Sint64 global_height_idx = (patch_orig_y_in_file + j) * whole_ij_bbox_dim.x + patch_orig_x_in_file;
        const mi::Sint64 global_normal_idx = global_height_idx; // the index are the same between height and normal

        const mi::Sint64 height_cur_fpos = height_begin_fpos + global_height_idx * sizeof(mi::Float32);
        const mi::Sint64 normal_cur_fpos = normal_begin_fpos + global_normal_idx * sizeof(mi::math::Vector_struct<mi::Float32, 3>);

        // Checking the range
        assert(height_cur_fpos >= height_begin_fpos);
        assert(normal_cur_fpos >= normal_begin_fpos);
        assert(height_cur_fpos <  height_begin_fpos + static_cast<mi::Sint64>((whole_ij_bbox_dim.x * whole_ij_bbox_dim.y) * sizeof(mi::Float32)));
        assert(normal_cur_fpos <  normal_begin_fpos + static_cast<mi::Sint64>((whole_ij_bbox_dim.x * whole_ij_bbox_dim.y) * sizeof(mi::math::Vector_struct<mi::Float32, 3 >)));

        // Writing the height values first
        mi::Sint64 fseek_success = fseek(outfile, height_cur_fpos, SEEK_SET);
        if (fseek_success != 0)
        {
            ERROR_LOG << "fseek for writing height values failed.";
            return false;
        }

        // Assuming continuous in the memory
        size_t write_success = fwrite(heights, sizeof(mi::Float32), nb_elements, outfile);
        if (write_success != static_cast<size_t>(nb_elements))
        {
            ERROR_LOG << "Writing height values failed.";
            return false;
        }

        // Then writing the height values
        fseek_success = fseek(outfile, normal_cur_fpos, SEEK_SET);
        if (fseek_success != 0)
        {
            ERROR_LOG << "fseek for writing height values failed.";
            return false;
        }
        // Assuming continuous in the memory
        write_success = fwrite(normals, sizeof(mi::math::Vector_struct<mi::Float32, 3 >), nb_elements, outfile);
        if (write_success != static_cast<size_t>(nb_elements))
        {
            ERROR_LOG << "Writing normal values failed.";
            return false;
        }

        write_byte_count += nb_elements * (sizeof(mi::Float32) + sizeof(mi::math::Vector_struct<mi::Float32, 3>));
    }

// #define REPORT_WRITE_BYTE_COUNT
#ifdef REPORT_WRITE_BYTE_COUNT
    mi::Sint64 estimated_item_count = whole_ij_bbox_dim.x * whole_ij_bbox_dim.y;
    mi::Size   file_byte_size = header_size +
        (sizeof(mi::Float32) + sizeof(mi::math::Vector_struct<mi::Float32, 3>)) * static_cast<mi::Size>(estimated_item_count);
    mi::Size   written_byte_size = header_size + write_byte_count;
    INFO_LOG << "estimated_file_size: " << file_byte_size 
             << " (header: " << header_size << "), written: " << written_byte_size;
#else
    nv::index_common::no_unused_variable_warning_please(write_byte_count);
#endif    

    return true;
}

//----------------------------------------------------------------------
bool Heightfield_data_exporter::export_to_local_file(
    const std::string & filename,
    mi::neuraylib::IDice_transaction*  dice_transaction) const
{
    assert(dice_transaction != 0);
    if (filename.empty())
    {
        ERROR_LOG << "empty output filename. no export.";
        return false;
    }

    if (m_ij_bbox.empty())
    {
        ERROR_LOG << "empty exporting bounding box. nothing will be exported.";
        return false;
    }
    INFO_LOG << "ij bounding box [min,max) = " << m_ij_bbox;

    if (!this->is_exporting_bbox_in_heightfield(dice_transaction))
    {
        return false;
    }
    
    // get exporting bounding boxes
    std::vector< mi::math::Bbox<mi::Sint32, 2> > export_bbox_vec;
    this->get_exporting_bbox(export_bbox_vec, dice_transaction);

    // DEBUG
    for (mi::Size i = 0; i < export_bbox_vec.size(); ++i)
    {
        // Found BUG here: twice
        // WARN_LOG << "DEBUG: export_bbox_vec[" << i << "]: " << export_bbox_vec.at(i);
    }


    // Start the fragmented job that runs only on the local machine.
    Heightfield_data_retrieval_appjob data_retrieval_job(
        m_heightfield_tag,
        export_bbox_vec,
        m_data_access_factory);
    const mi::Size nb_query_bbox = export_bbox_vec.size();
    if (nb_query_bbox == 0)
    {
        ERROR_LOG << "exporting bounding box list is empty. The exporting bbox may be not valid. no export.";        
        return false;
    }

    dice_transaction->execute_fragmented(&data_retrieval_job, 1);

    // Open a file for exporting the heightfield data
    assert(!filename.empty());
    FILE* export_file = fopen(filename.c_str(), "wb+");
    if (export_file == 0)
    {
        ERROR_LOG << "Failed to open exporting file [" << filename << "]. no export.";
        return false;
    }
    INFO_LOG << "Exporting heightfield to file = [" << filename << "]";

    // Writing the header
    const std::string heightfield_name = "exported_heightfield";
    const mi::math::Vector<mi::Sint32,2> whole_ij_bbox_dim = m_ij_bbox.max - m_ij_bbox.min;
    INFO_LOG << "Heightfield size: (" << whole_ij_bbox_dim;

    std::stringstream sstr;
    sstr << "#! index_raw_heightfield 2\n"
         << "#name = " << heightfield_name << "\n"
         << "#size = " << whole_ij_bbox_dim.x << " " << whole_ij_bbox_dim.y << "\n"
         << "#-end\n";

    fprintf(export_file, "%s", sstr.str().c_str());
    fflush(export_file);
    const mi::Size header_size = sstr.str().length();

    // Keep end of header position now 
    const mi::Sint64 end_of_header_fpos = ftell(export_file);
    
    // Now, the returned heightfield data can be stored as they are local now ....
    for (mi::Size i = 0; i < nb_query_bbox; ++i)
    {
        mi::base::Handle<nv::index::IRegular_heightfield_data_access> heightfield_data_access =
            data_retrieval_job.get_heightfield_data(static_cast<mi::Uint32>(i));

        INFO_LOG << "Received heightfield data patch bbox: " << 
            heightfield_data_access->get_patch_bounding_box();
        bool success = serialized_write_to_file(export_file, header_size, end_of_header_fpos, heightfield_data_access);
        if (!success)
        {
            fclose(export_file);
            return false;
        }
    }

    const mi::Sint32 success_fclose = fclose(export_file);
    if (success_fclose != 0)
    {
        // fclose failed
        ERROR_LOG << "Cannot close the file [" << filename << "]";
        return false;
    }

    return true;
}

//----------------------------------------------------------------------
bool export_heightfield_data(
    const mi::neuraylib::Tag &            heightfield_tag,
    const mi::math::Bbox<mi::Sint32, 2> & ij_bbox,
    const std::string&                    export_fname,
    const mi::neuraylib::Tag&             session_tag,
    mi::neuraylib::IDice_transaction *    dice_transaction)
{
    assert(heightfield_tag.is_valid());
    assert(session_tag.is_valid());
    assert(dice_transaction != 0);

    if (export_fname.empty())
    {
        ERROR_LOG << "empty output heightfield filename. no export.";
        return false;
    }

    if (ij_bbox.empty())
    {
        ERROR_LOG << "ij bbox is empty. no export.";        
        return false;
    }
    
    if (!ij_bbox.is_plane())
    {
        ERROR_LOG << "ij bbox is not a plane. no export.";        
        return false;
    }


    bool is_success = false;
    {
        mi::base::Handle<const nv::index::ISession> session(
            dice_transaction->access<const nv::index::ISession>(session_tag));
        assert(session.is_valid_interface());
        
        mi::base::Handle<const nv::index::IData_distribution> subsurface_data_distribution(
            dice_transaction->access<const nv::index::IData_distribution>(
                session->get_distribution_layout()));
        assert(subsurface_data_distribution.is_valid_interface());

        mi::base::Handle<const nv::index::IDistributed_data_access_factory>
            subsurface_data_access_factory(
                dice_transaction->access<const nv::index::IDistributed_data_access_factory>(
                    session->get_data_access_factory()));
        assert(subsurface_data_access_factory.is_valid_interface());

        {
            mi::base::Handle<const nv::index::IRegular_heightfield> heightfield(
                dice_transaction->access<const nv::index::IRegular_heightfield>(heightfield_tag));
            if (!heightfield_tag.is_valid())
            {
                ERROR_LOG << "invalid heightfield_tag: " << heightfield_tag;
                return false;
            }
            assert(heightfield.is_valid_interface());

            // check the ROI validity
            mi::math::Bbox<mi::Float32, 3> xyz_roi_float =
                get_XYZ_global_region_of_interest_bbox(session.get(), dice_transaction);
            mi::math::Bbox<mi::Float32, 3> ijk_bbox_float = heightfield->get_IJK_bounding_box();
            if (!(xyz_roi_float.contains(ijk_bbox_float.min) && xyz_roi_float.contains(ijk_bbox_float.max)))
            {
                // ijk_bbox is not inside (or on) the xyz_roi
                WARN_LOG << "The exporting heightfield ijk_bbox is outside of the global xyz ROI. "
                         << "This will give a corrupted result. global ROI:"
                         << xyz_roi_float << ", ijk_bbox: " << ijk_bbox_float
                         << ". All data should be loaded.";
            }

            Heightfield_data_exporter exporter(
                heightfield_tag,
                ij_bbox,
                subsurface_data_distribution,
                subsurface_data_access_factory);
            
            INFO_LOG << "Exporting heightfield data in bbox " << ij_bbox 
                     << " to [" << export_fname << "].";

            is_success = exporter.export_to_local_file(export_fname, dice_transaction);
            INFO_LOG << "Finished exporting heightfield data."; 
        }
    }

    return is_success;
}

//----------------------------------------------------------------------
