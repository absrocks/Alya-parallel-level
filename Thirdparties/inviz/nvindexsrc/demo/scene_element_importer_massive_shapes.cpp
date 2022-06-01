/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/

#include "scene_element_importer.h"

#include "common/triangle_mesh_importer.h"
#include "common/synthetic_reservoir_patch_generator.h"
#include "common/reservoir_grid_importer.h"

#include "common/ppm_heightfield_importer.h"
#include "common/raw_heightfield_data_importer.h"
#include "common/raw_heightfield_data_importer_ssv.h"
#include "common/synthetic_heightfield_generator.h"

#include "common/multi_attribute_scaling_sequence_importer.h"
#include "common/raw_volume_data_importer.h"
#include "common/raw_volume_data_sequence_importer.h"
#include "common/repeat_raw_volume_data_importer.h"
#include "common/rgba_raw_volume_data_importer.h"
#include "common/single_value_raw_volume_data_importer.h"
#include "common/short_value_raw_volume_data_importer.h"
#include "common/synthetic_volume_generator.h"

#include "common/irregular_volume_importer.h"
#include "common/irregular_volume_importer_vtk.h"

#ifdef NV_IDX_USE_OPENVDB_IMPORTER
#include "common/sparse_volume_importer.h"
#endif

#include "colormap_util.h"
#include "nvindex_appdata.h"
#include "scene_utility.h"
#include "utilities.h"
#include "volume_data_filter.h"

using namespace nv::index_common;

namespace {

//----------------------------------------------------------------------
// importer/generator configuration
//----------------------------------------------------------------------
/// get "raw" importer configuration
/// \param[in] volume_fname volume file name
/// \param[in] volume_size  volume size
/// \return Returns the raw volume data importer.
Raw_volume_data_importer* get_raw_volume_data_importer(const nv::index_common::String_dict& volume_dict)
{
    char const * const p_key[] = { "input_file", "size", 0, };
    if (!check_necessary_key(volume_dict, p_key))
    {
        ERROR_LOG << "get_raw_volume_data_importer: return invalid importer.";
        return NULL;
    }

    const std::string volume_fname = volume_dict.get("input_file");
    const mi::math::Vector<mi::Uint32, 3> volume_size = get_vec_uint32_3(volume_dict.get("size"));

    Raw_volume_data_importer* importer = new Raw_volume_data_importer(volume_fname.c_str(), volume_size);

    importer->set_use_cache(get_bool(volume_dict.get("use_cache", "true")));
    importer->set_cache_compression(get_uint32(volume_dict.get("cache_compression", "0")));
    importer->set_serial_access(get_bool(volume_dict.get("serial_access", "false")));
    importer->set_show_stats(get_bool(volume_dict.get("stats", "false")));
    importer->set_demo_scale_hack_00(get_bool(volume_dict.get("demo_brick_scale_2x2x2", "false")));

    return importer;
}

//----------------------------------------------------------------------
// importer/generator configuration
//----------------------------------------------------------------------
/// get "raw" importer configuration
/// \param[in] volume_fname volume file name
/// \param[in] volume_size  volume size
// \return Returns the raw volume data importer.
Rgba_raw_volume_data_importer* get_raw_rgba_volume_data_importer(const nv::index_common::String_dict& volume_dict)
{
    char const * const p_key[] = { "input_file", "size", 0, };
    if (!check_necessary_key(volume_dict, p_key))
    {
        ERROR_LOG << "get_raw_rgba_volume_data_importer: return invalid importer.";
        return NULL;
    }

    const std::string volume_fname = volume_dict.get("input_file");
    const std::string src_folder = volume_dict.get("src_folder", "");
    const std::string dst_folder = volume_dict.get("dst_folder", "");
    const mi::math::Vector<mi::Uint32, 3> volume_size = get_vec_uint32_3(volume_dict.get("size"));

    Rgba_raw_volume_data_importer* importer = new Rgba_raw_volume_data_importer(
        volume_fname.c_str(),
        src_folder,
        dst_folder,
        volume_size);

    importer->set_use_cache(get_bool(volume_dict.get("use_cache", "false")));
    importer->set_cache_compression(get_uint32(volume_dict.get("cache_compression", "0")));
    const mi::Uint32 b_frame_offset = get_uint32(volume_dict.get("b_frame_offset", "0"));
    importer->set_use_difference_encoding(get_bool(volume_dict.get("use_difference_encoding", "false")));
    importer->set_b_frame_offset(b_frame_offset);
    const mi::Uint32 nb_time_steps = get_uint32(volume_dict.get("nb_time_steps", "0"));
    for(mi::Uint32 time_step=0; time_step<nb_time_steps; time_step++)
    {
        std::ostringstream time_step_str;
        time_step_str << "time_step_" << time_step;
        const std::string time_step_file = volume_dict.get(time_step_str.str());
        importer->add_time_step(time_step_file);
    }
    importer->set_serial_access(get_bool(volume_dict.get("serial_access", "false")));
    importer->set_show_stats(get_bool(volume_dict.get("stats", "false")));

    return importer;
}


//----------------------------------------------------------------------
// importer/generator configuration
//----------------------------------------------------------------------
/// get "raw" importer configuration
/// \param[in] volume_fname volume file name
/// \param[in] volume_size  volume size
/// \return Returns the raw volume data importer.
Single_value_raw_volume_data_importer* get_raw_float_volume_data_importer(const nv::index_common::String_dict& volume_dict)
{
    char const * const p_key[] = { "input_file", "size", 0, };
    if (!check_necessary_key(volume_dict, p_key))
    {
        ERROR_LOG << "get_raw_float_volume_data_importer: return invalid importer.";
        return NULL;
    }

    const std::string volume_fname = volume_dict.get("input_file");
    const std::string src_folder = volume_dict.get("src_folder", "");
    const std::string dst_folder = volume_dict.get("dst_folder", "");
    const mi::math::Vector<mi::Uint32, 3> volume_size = get_vec_uint32_3(volume_dict.get("size"));
    const mi::math::Vector<mi::Float32, 2> range = get_vec_float32_2(volume_dict.get("range", "0 0"));

    Single_value_raw_volume_data_importer* importer = new Single_value_raw_volume_data_importer(
        volume_fname.c_str(),
        src_folder,
        dst_folder,
        volume_size);

    importer->set_range(range);
    importer->set_use_cache(get_bool(volume_dict.get("use_cache", "false")));
    importer->set_cache_compression(get_uint32(volume_dict.get("cache_compression", "0")));
    const mi::Uint32 b_frame_offset = get_uint32(volume_dict.get("b_frame_offset", "0"));
    importer->enable_difference_encoding(get_bool(volume_dict.get("use_difference_encoding", "false")));
    importer->set_b_frame_offset(b_frame_offset);
    const mi::Uint32 nb_time_steps = get_uint32(volume_dict.get("nb_time_steps", "0"));
    for(mi::Uint32 time_step=0; time_step<nb_time_steps; time_step++)
    {
        std::ostringstream time_step_str;
        time_step_str << "time_step_" << time_step;
        const std::string time_step_file = volume_dict.get(time_step_str.str());
        importer->add_time_step(time_step_file);
    }
    importer->set_serial_access(get_bool(volume_dict.get("serial_access", "false")));
    importer->set_show_stats(get_bool(volume_dict.get("stats", "false")));

    return importer;
}

//----------------------------------------------------------------------
// importer/generator configuration
//----------------------------------------------------------------------
/// get "raw" importer configuration
/// \param[in] volume_fname volume file name
/// \param[in] volume_size  volume size
/// \return Returns the raw volume data importer.
Short_value_raw_volume_data_importer* get_raw_short_volume_data_importer(const nv::index_common::String_dict& volume_dict)
{
    char const * const p_key[] = { "input_file", "size", 0, };
    if (!check_necessary_key(volume_dict, p_key))
    {
        ERROR_LOG << "get_raw_short_volume_data_importer: return invalid importer.";
        return NULL;
    }

    const std::string volume_fname = volume_dict.get("input_file");
    const std::string src_folder = volume_dict.get("src_folder", "");
    const std::string dst_folder = volume_dict.get("dst_folder", "");
    const mi::math::Vector<mi::Uint32, 3> volume_size = get_vec_uint32_3(volume_dict.get("size"));

    Short_value_raw_volume_data_importer* importer = new Short_value_raw_volume_data_importer(
        volume_fname.c_str(),
        src_folder,
        dst_folder,
        volume_size);

    importer->set_use_cache(get_bool(volume_dict.get("use_cache", "false")));
    importer->set_cache_compression(get_uint32(volume_dict.get("cache_compression", "0")));
    const mi::Uint32 nb_time_steps = get_uint32(volume_dict.get("nb_time_steps", "0"));
    for(mi::Uint32 time_step=0; time_step<nb_time_steps; time_step++)
    {
        std::ostringstream time_step_str;
        time_step_str << "time_step_" << time_step;
        const std::string time_step_file = volume_dict.get(time_step_str.str());
        importer->add_time_step(time_step_file);
    }
    importer->set_serial_access(get_bool(volume_dict.get("serial_access", "false")));
    importer->set_show_stats(get_bool(volume_dict.get("stats", "false")));

    return importer;
}

//----------------------------------------------------------------------
// importer/generator configuration
//----------------------------------------------------------------------
/// get "raw" importer configuration
/// \param[in] volume_fname volume file name
/// \param[in] volume_size  volume size
/// \return Returns the raw volume data importer.
Raw_volume_data_sequence_importer* get_raw_volume_data_sequence_importer(
    const nv::index_common::String_dict& volume_dict)
{
    char const * const p_key[] = { "input_file", "size", 0, };
    if (!check_necessary_key(volume_dict, p_key))
    {
        ERROR_LOG << "get_raw_volume_data_importer: return invalid importer.";
        return NULL;
    }

    const std::string volume_fname = volume_dict.get("input_file");
    const mi::math::Vector<mi::Uint32, 3> volume_size = get_vec_uint32_3(volume_dict.get("size"));

    Raw_volume_data_sequence_importer* importer = new Raw_volume_data_sequence_importer(volume_fname.c_str(), volume_size);

	importer->set_time_interval(get_uint32(volume_dict.get("time_interval", "0")));
	importer->set_start_time(get_uint32(volume_dict.get("start_time", "0")));
	importer->set_use_cache(get_bool(volume_dict.get("use_cache", "true")));
    importer->set_cache_compression(get_uint32(volume_dict.get("cache_compression", "0")));
    importer->set_serial_access(get_bool(volume_dict.get("serial_access", "false")));
    importer->set_show_stats(get_bool(volume_dict.get("stats", "false")));
    importer->set_demo_scale_hack_00(get_bool(volume_dict.get("demo_brick_scale_2x2x2", "false")));

    return importer;
}

//----------------------------------------------------------------------
/// get "repeat" importer configuration
/// \param[in,out] volume_dict (input/output) repeat parameter in a
/// dict, size will be updated.
/// \return Returns the repeat importer.
Multi_attribute_scaling_sequence_importer* get_multi_attribute_scaling_sequence_importer(
    const nv::index_common::String_dict& volume_dict)
{
    char const * const p_key[] = { "attribute_file_1", "attribute_file_2", "output_file", "size", 0, };
    if (!check_necessary_key(volume_dict, p_key))
    {
        ERROR_LOG << "get_multi_attribute_scaling_sequence_importer: returned invalid importer.";
        return NULL;
    }

    const std::string volume_fname_1 = volume_dict.get("attribute_file_1");
    const std::string volume_fname_2 = volume_dict.get("attribute_file_2");
    const std::string output_fname = volume_dict.get("output_file");
    const mi::math::Vector<mi::Uint32, 3> volume_size = get_vec_uint32_3(volume_dict.get("size"));

    Multi_attribute_scaling_sequence_importer* importer = new Multi_attribute_scaling_sequence_importer(
        volume_fname_1.c_str(), 
        volume_fname_2.c_str(), 
        output_fname,
        volume_size);

	importer->set_attribute_minmax(1, get_vec_float32_2(volume_dict.get("attribute_min_max_1")));
	importer->set_attribute_minmax(2, get_vec_float32_2(volume_dict.get("attribute_min_max_2")));
	importer->set_time_interval(get_uint32(volume_dict.get("time_interval", "0")));
	importer->set_start_time(get_uint32(volume_dict.get("start_time", "0")));
	importer->set_use_cache(get_bool(volume_dict.get("use_cache", "true")));
    importer->set_cache_compression(get_uint32(volume_dict.get("cache_compression", "0")));
    importer->set_serial_access(get_bool(volume_dict.get("serial_access", "false")));
    importer->set_show_stats(get_bool(volume_dict.get("stats", "false")));

    return importer;
}

//----------------------------------------------------------------------
/// get "repeat" importer configuration
/// \param[in,out] volume_dict (input/output) repeat parameter in a
/// dict, size will be updated.
/// \return Returns the repeat importer.
Repeat_raw_volume_data_importer* get_repeat_volume_importer(nv::index_common::String_dict& volume_dict)
{
    const char* const p_key[] = { "input_file", "size", "repeat", 0, };
    if (!check_necessary_key(volume_dict, p_key))
    {
        ERROR_LOG << "get_repeat_volume_importer: return invalid importer.";
        return NULL;
    }

    const std::string volume_fname = volume_dict.get("input_file");
    mi::math::Vector<mi::Uint32, 3> volume_size = get_vec_uint32_3(volume_dict.get("size"));
    const mi::math::Vector<mi::Uint32, 3> repeat =
        get_vec_uint32_3(volume_dict.get("repeat", "1 1 1"));
    bool const cache_source_data = get_bool(volume_dict.get("cache_source", "true"));

    Repeat_raw_volume_data_importer* importer = new Repeat_raw_volume_data_importer(
        volume_fname.c_str(), volume_size, repeat, cache_source_data);

    // This importer modifies the volume size, so update the size here.
    volume_size = importer->get_volume_size();
    std::stringstream sstr;
    sstr << volume_size.x << " " << volume_size.y << " " << volume_size.z;
    volume_dict.insert("size", sstr.str());

    return importer;
}

//----------------------------------------------------------------------
/// get "synthetic" volume generation configuration
/// \param[in] volume_dict  The dictionary provides key-value pairs passed with the scene files.
/// \return                 An instance of the synthetic volume generator class.
Synthetic_volume_generator* get_synthetic_volume_generator(const nv::index_common::String_dict& volume_dict)
{
    const char* const p_key[] = { "size", 0, };
    if (!check_necessary_key(volume_dict, p_key))
    {
        ERROR_LOG << "get_synthetic_volume_generator: return invalid.";
        return NULL;
    }

    const mi::math::Vector<mi::Uint32, 3> volume_size = get_vec_uint32_3(volume_dict.get("size"));
    const mi::math::Bbox<mi::Uint32, 3> ijk_bbox(mi::math::Vector<mi::Uint32, 3>(0), volume_size);

    const std::string synthetic_type_str = volume_dict.get("synthetic_type", "default");
    const std::string parameter_str      = volume_dict.get("parameter", "");
    const std::string voxel_format_str   = volume_dict.get("voxel_format", "uint8");
    
    INFO_LOG << "[generator]: synthetic: synthetic_type [" << synthetic_type_str
             << "], parameter [" << parameter_str << "], voxel_format [" << voxel_format_str << "]";

    const mi::math::Vector<mi::Float32, 3> trans = get_vec_float32_3(volume_dict.get("translate", "0 0 0"));
    if (trans != mi::math::Vector<mi::Float32, 3>(0.f))
    {
        INFO_LOG << "----------------------------------------------------------------------";
        INFO_LOG << "The synthetic volume generator does not support volume translation, "
                 << "results may be incorrect.";
        INFO_LOG << "----------------------------------------------------------------------";
    }

    Synthetic_volume_generator::Voxel_format voxel_format = Synthetic_volume_generator::VOXEL_FORMAT_UINT_8;
    if (voxel_format_str == "uint8")
    {
        voxel_format = Synthetic_volume_generator::VOXEL_FORMAT_UINT_8;
    }
    else if (voxel_format_str == "rgba8")
    {
        voxel_format = Synthetic_volume_generator::VOXEL_FORMAT_RGBA_8;
    }
    else if (voxel_format_str == "uint16")
    {
        voxel_format = Synthetic_volume_generator::VOXEL_FORMAT_UINT_16;
    }
    else if (voxel_format_str == "float32")
    {
        voxel_format = Synthetic_volume_generator::VOXEL_FORMAT_FLOAT_32;
    }
    else
    {
        ERROR_LOG << "Unknown voxel format '" << voxel_format_str << "' for Synthetic_volume_generator";
        return NULL;
    }

    Synthetic_volume_generator* generator = new Synthetic_volume_generator(
            ijk_bbox,
            voxel_format,
            synthetic_type_str,
            parameter_str);
    assert(generator);

    return generator;
}

bool check_scope(
    const mi::base::Handle<mi::neuraylib::IScope>& scope,
    const nv::index_common::String_dict&           dict,
    mi::neuraylib::Tag                             old_tag)
{
    if (scope->get_privacy_level() == 0)
    {
        // In global scope, all fine
        return true;
    }

    if (!old_tag.is_valid())
    {
        ERROR_LOG << "Massive scene elements may only be created in global scope";
        return false;
    }
    else if (get_bool(dict.get("allow_reloading", "no")))
    {
        ERROR_LOG << "Massive scene element importers may only be modified in global scope";
        return false;
    }

    return true;
}

} // namespace

mi::base::Handle<nv::index::IScene_element> Scene_element_importer::import_massive_shapes(
    const std::string&                                        elem_name,
    const std::string&                                        elem_type,
    nv::index_common::String_dict&                            dict,
    mi::neuraylib::Tag                                        old_tag,
    const mi::base::Handle<mi::neuraylib::IDice_transaction>& dice_transaction,
    bool&                                                     reuse_tag)
{
    mi::base::Handle<mi::neuraylib::IScope> scope(dice_transaction->get_scope());

    //
    // Volume
    //
    if (elem_type == "volume")
    {
        // Volume size in three dimensions
        mi::math::Vector<mi::Uint32, 3> volume_size = get_vec_uint32_3(dict.get("size", "0 0 0"));

        if (volume_size.x == 0 && volume_size.y == 0 && volume_size.z == 0)
        {
            ERROR_LOG << "No size given for volume '" << elem_name << "'";
            return null_result();
        }

        //
        // Prepare the volume importer
        //
        mi::base::Handle<nv::index::IRegular_volume> volume;

        if (!check_scope(scope, dict, old_tag))
        {
            ERROR_LOG << "Invalid scope for operation on volume " << old_tag;
        }
        else if (old_tag.is_valid() && !get_bool(dict.get("allow_reloading", "no")))
        {
            // Just edit the existing volume instead of loading it again
            volume = dice_transaction->edit<nv::index::IRegular_volume>(old_tag);
            if (volume)
            {
                reuse_tag = true;
            }
        }
        else
        {
            // Create new volume
            const std::string importer = dict.get("importer", "raw");
            nv::index::IDistributed_discrete_data_import_callback* importer_callback = NULL;
            if (importer == "raw")
            {
                importer_callback = get_raw_volume_data_importer(dict);
            }
            else if (importer == "raw_rgba")
            {
                importer_callback = get_raw_rgba_volume_data_importer(dict);
            }
            else if (importer == "raw_float")
            {
                importer_callback = get_raw_float_volume_data_importer(dict);
            }
            else if(importer == "raw_short")
            {
                importer_callback = get_raw_short_volume_data_importer(dict);
            }
            else if (importer == "raw_sequence")
            {
                importer_callback = get_raw_volume_data_sequence_importer(dict);
            }
            else if (importer == "multi_attribute_scaling_sequence")
            {
                importer_callback = get_multi_attribute_scaling_sequence_importer(dict);
            }
            else if (importer == "repeat")
            {
                importer_callback = get_repeat_volume_importer(dict);
            }
            else if (importer == "artificial")
            {
                WARN_LOG << "importer type 'artificial' is obsolete. Please use 'synthetic'";
                importer_callback = get_synthetic_volume_generator(dict);
            }
            else if (importer == "synthetic")
            {
                importer_callback = get_synthetic_volume_generator(dict);
            }
            else
            {
                ERROR_LOG << "No such registered importer for volume: " << importer;
                return null_result();
            }

            // The importer might have modified the volume size (e.g. repeat), so update here
            volume_size = get_vec_uint32_3(dict.get("size"));

            //
            // Define the volume
            //
            mi::math::Vector<mi::Float32, 3> scale = get_vec_float32_3(dict.get("scale", "1 1 1"));

            mi::Float32 rotate_k = get_float32(dict.get("rotate", "0")) / 180.f * static_cast<mi::Float32>(MI_PI);
            // Value in radians takes precedence (to prevent numerical inaccuracies during conversion)
            if (dict.is_defined("rotate_rad"))
            {
                rotate_k = get_float32(dict.get("rotate_rad"));
            }

            mi::math::Vector<mi::Float32, 3> translate = get_vec_float32_3(dict.get("translate", "0 0 0"));

            // Create the volume scene element
            volume = m_scene->create_volume(
                scale, rotate_k, translate, volume_size,
                importer_callback, m_dice_transaction.get());
        }

        if (volume == 0)
        {
            ERROR_LOG << "Failed to create volume '" << elem_name << "'";
            return null_result();
        }

        // Set properties
        volume->set_name(dict.get("name", "unnamed volume").c_str());

        // Assign colormap
        mi::neuraylib::Tag colormap_tag;
        if (dict.is_defined("colormap"))
        {
            colormap_tag = resolve_colormap(dict.get("colormap"));
        }
        else
        {
            colormap_tag = get_colormap_tag(Nvindex_AppData::instance()->get_current_colormap_index());
        }
        volume->assign_colormap(colormap_tag);

        // Clipping
        if (dict.is_defined("clip"))
        {
            volume->set_IJK_region_of_interest(get_bbox_float32_3(dict.get("clip")));
        }

        //
        // Create slices through the volume
        //

        if (scope->get_privacy_level() > 0)
        {
            //TODO: Implement support for this
            INFO_LOG << "Note: Editing volume " << old_tag << " while not in global scope, "
                     << "changes to section slices will be ignored";
        }
        else
        {
            // First clear any existing slices
            while (volume->get_nb_slices() > 0)
            {
                volume->remove_slice(volume->get_slice(0));
            }

            if (dict.is_defined("nb_slices"))
            {
                const mi::Uint32 nb_slices = get_uint32(dict.get("nb_slices"));
                for (mi::Uint32 i = 0; i < nb_slices; ++i)
                {
                    std::ostringstream s;
                    s << "slice" << i << "::";
                    const std::string header = s.str();

                    mi::neuraylib::Tag slice_tag;
                    const bool enabled = get_bool(dict.get(header + "enabled", "yes"));
                    const mi::neuraylib::Tag colormap_tag = resolve_colormap(
                        dict.get(header + "colormap"));

                    if (dict.get(header + "type", "") == "vertical_profile")
                    {
                        std::vector<mi::math::Vector<mi::Float32, 2> > vertices;
                        const mi::Uint32 nb_vertices = get_uint32(dict.get(header + "nb_vertices", "0"));
                        for (mi::Uint32 j = 0; j < nb_vertices; ++j)
                        {
                            std::ostringstream name;
                            name << header << "v" << j;
                            if (dict.is_defined(name.str()))
                            {
                                vertices.push_back(get_vec_float32_2(dict.get(name.str())));
                            }
                        }

                        slice_tag = volume->add_vertical_profile(
                            m_dice_transaction.get(),
                            &vertices[0],
                            vertices.size(),
                            colormap_tag,
                            enabled);
                    }
                    else
                    {
                        nv::index::ISection_scene_element::Slice_orientation orientation;
                        const std::string s = dict.get(header + "orientation");
                        if (s == "inline")
                            orientation = nv::index::ISection_scene_element::INLINE_SECTION;
                        else if (s == "cross")
                            orientation = nv::index::ISection_scene_element::CROSS_LINE_SECTION;
                        else
                            orientation = nv::index::ISection_scene_element::HORIZONTAL_SECTION;

                        slice_tag = volume->add_section(
                            m_dice_transaction.get(),
                            orientation,
                            get_float32(dict.get(header + "position", "0")),
                            colormap_tag,
                            enabled);
                    }

                    mi::base::Handle<nv::index::IShape> slice(
                        dice_transaction->edit<nv::index::IShape>(slice_tag));
                    if (slice.is_valid_interface())
                    {
                        slice->set_pickable(get_bool(dict.get(header + "pickable", "yes")));
                    }
                }
            }
            else
            {
                //
                // Default: Three section slices
                //
                const mi::Uint32 colormap_count = get_number_of_colormap();
                const mi::neuraylib::Tag colormap_tag_0 = get_colormap_tag(std::min(3u, colormap_count - 1));

                const mi::neuraylib::Tag inline_section_tag = volume->add_section(
                    m_dice_transaction.get(), nv::index::ISection_scene_element::INLINE_SECTION,
                    mi::Float32(volume_size.x / 4.f), colormap_tag_0, false);
                assert(inline_section_tag.is_valid());
                nv::index_common::no_unused_variable_warning_please(inline_section_tag);

                const mi::neuraylib::Tag colormap_tag_1 = get_colormap_tag(std::min(5u, colormap_count - 1));
                const mi::neuraylib::Tag crossline_section_tag = volume->add_section(
                    m_dice_transaction.get(), nv::index::ISection_scene_element::CROSS_LINE_SECTION,
                    mi::Float32(volume_size.y / 4.f), colormap_tag_1, false);
                assert(colormap_tag_1       .is_valid());
                assert(crossline_section_tag.is_valid());
                nv::index_common::no_unused_variable_warning_please(crossline_section_tag);

                const mi::neuraylib::Tag colormap_tag_2 = get_colormap_tag(std::min(7u, colormap_count - 1));
                const mi::neuraylib::Tag horizontal_section_tag = volume->add_section(
                    m_dice_transaction.get(), nv::index::ISection_scene_element::HORIZONTAL_SECTION,
                    mi::Float32(volume_size.z / 4.f), colormap_tag_2, false);
                assert(horizontal_section_tag.is_valid());
                nv::index_common::no_unused_variable_warning_please(horizontal_section_tag);

                // Vertical Profile
                const mi::neuraylib::Tag colormap_tag_3 = get_colormap_tag(std::min(9u, colormap_count - 1));

                std::vector<mi::math::Vector<mi::Float32, 2> > v(5);
                v[0].x = volume_size.x * 0.7f; v[0].y = volume_size.y * 0.1f;
                v[1].x = volume_size.x * 0.2f; v[1].y = volume_size.y * 0.3f;
                v[2].x = volume_size.x * 0.1f; v[2].y = volume_size.y * 0.8f;
                v[3].x = volume_size.x * 0.6f; v[3].y = volume_size.y * 0.6f;
                v[4].x = volume_size.x * 0.3f; v[4].y = volume_size.y * 1.0f;

                const mi::neuraylib::Tag vertical_profile_tag = volume->add_vertical_profile(
                    m_dice_transaction.get(), &v[0], v.size(), colormap_tag_3, false);
                assert(vertical_profile_tag.is_valid());
                nv::index_common::no_unused_variable_warning_please(vertical_profile_tag);
            }
        }

        return volume;
    }
    //
    // Heightfield
    //
    else if (elem_type == "horizon" || elem_type == "heightfield")
    {
        if (elem_type == "horizon")
        {
            WARN_LOG << "'type = horizon' is obsolete, please use 'type = heightfield' instead." ;
        }

        // Heightfield dataset file
        const std::string heightfield_file = dict.get("input_file");

        // Heightfield size in two dimensions
        const mi::math::Vector<mi::Uint32, 2> heightfield_size(get_vec_uint32_2(dict.get("size", "0 0")));

        if (heightfield_size.x == 0 && heightfield_size.y == 0)
        {
            ERROR_LOG << "No size given for heightfield '" << elem_name << "'";
            return null_result();
        }

        //
        // Prepare the heightfield importer
        //
        mi::base::Handle<nv::index::IRegular_heightfield> heightfield;

        if (!check_scope(scope, dict, old_tag))
        {
            ERROR_LOG << "Invalid scope for operation on heightfield " << old_tag;
        }
        else if (old_tag.is_valid() && !get_bool(dict.get("allow_reloading", "no")))
        {
            // Just edit the existing heightfield instead of loading it again
            heightfield = dice_transaction->edit<nv::index::IRegular_heightfield>(old_tag);
            if (heightfield)
                reuse_tag = true;
        }
        else
        {
            // Create new heightfield
            const std::string heightfield_importer = dict.get("importer", "raw");
            nv::index::IDistributed_discrete_data_import_callback* importer_callback = NULL;

            if (heightfield_importer == "raw")
            {
                nv::index_common::Raw_heightfield_data_importer* importer =
                    new nv::index_common::Raw_heightfield_data_importer(heightfield_file, heightfield_size);

                importer->set_use_cache(get_bool(dict.get("cache", "false")));
                importer->set_regenerate_normals(get_bool(dict.get("regenerate_normals", "false")));
                importer->set_serial_access(get_bool(dict.get("serial_access", "false")));
                importer->set_show_stats(get_bool(dict.get("stats", "false")));

                importer_callback = importer;

                if (heightfield_file.empty())
                {
                    if (!dict.is_defined("importer"))
                    {
                        INFO_LOG << "Heightfield section has no [importer], use raw importer.";
                    }

                    ERROR_LOG << "Heightfield raw importer misses [input_file] item.";
                }
            }
            else if (heightfield_importer == "raw_ssv") // FIXME: obsolete
            {
                nv::index_common::Raw_heightfield_data_importer_ssv* importer =
                    new nv::index_common::Raw_heightfield_data_importer_ssv(heightfield_file, heightfield_size);

                importer->set_use_cache(get_bool(dict.get("cache", "false")));
                importer->set_serial_access(get_bool(dict.get("serial_access", "false")));
                importer->set_show_stats(get_bool(dict.get("stats", "false")));

                importer_callback = importer;

                if (heightfield_file.empty())
                {
                    ERROR_LOG << "Heightfield raw_ssv importer misses [input_file] item.";
                }
            }
            else if (heightfield_importer == "ppm")
            {
                const mi::Float32 scale  = get_float32(dict.get("importer_scale",  "1"));
                const mi::Float32 offset = get_float32(dict.get("importer_offset", "0"));

                // Binary mask - must be PBM (portable graymap) formatted form and must
                // have the same dimension as the heightfield dataset file!
                const std::string binary_mask_file = dict.get("importer_binary_mask", "");

                nv::index_common::PPM_heightfield_importer* importer =
                    new nv::index_common::PPM_heightfield_importer(heightfield_file,
                                                                   scale,
                                                                   offset,
                                                                   binary_mask_file,   // PBM (portable bitmap) that defines holes in the surface
                                                                   heightfield_size);
            
                importer_callback = importer;

                if (heightfield_file.empty())
                {
                    ERROR_LOG << "Heightfield ppm importer misses [input_file] item.";
                }
            }
            else if (heightfield_importer == "artificial" || heightfield_importer == "synthetic")
            {
                if(heightfield_importer == "artificial")
                {
                    WARN_LOG << "heightfield importer type 'artificial' is obsolete. Please use 'synthetic' instead.";
                }

                const std::string tech_name  = dict.get("synthetic_type", "default");
                const std::string tech_param = dict.get("parameter",      "");
                nv::index_common::Synthetic_heightfield_generator* importer =
                    new nv::index_common::Synthetic_heightfield_generator(tech_name, tech_param, heightfield_size);

                importer_callback = importer;
            }
            else
            {
                ERROR_LOG << "The heightfield importer '"
                          << heightfield_importer
                          << "' is not supported through the scene file.";
                return null_result();
            }

            //
            // Define the heightfield
            //
            mi::math::Vector<mi::Float32, 3> scale = get_vec_float32_3(dict.get("scale", "1 1 1"));

            mi::Float32 rotate_k = get_float32(dict.get("rotate", "0")) / 180.f * static_cast<mi::Float32>(MI_PI);
            // Value in radians takes precedence (to prevent numerical inaccuracies during conversion)
            if (dict.is_defined("rotate_rad"))
            {
                rotate_k = get_float32(dict.get("rotate_rad"));
            }

            mi::math::Vector<mi::Float32, 3> translate = get_vec_float32_3(dict.get("translate", "0 0 0"));

            // The range in which the k values are in (required by the distribution scheme)
            if (!dict.is_defined("range"))
            {
                WARN_LOG << "****** No value range given for heightfield '" << elem_name << "', "
                         << "using defaults [0 2000].";
            }
            mi::math::Vector<mi::Float32, 2> elevation_range(get_vec_float32_2(dict.get("range", "0 2000")));

            // Create the heightfield scene element
            heightfield = m_scene->create_regular_heightfield(
                scale, rotate_k, translate,
                heightfield_size,
                elevation_range,
                importer_callback,
                dice_transaction.get());
            assert(heightfield != 0);
        }

        if (heightfield == 0)
        {
            ERROR_LOG << "Failed to create heightfield '" << elem_name << "'";
            return null_result();
        }

        // Set properties
        heightfield->set_name(dict.get("name", "unnamed heightfield").c_str());

        if (dict.is_defined("color"))
        {
            const mi::math::Color color = get_color(dict.get("color", "1 1 1 1"));
            WARN_LOG << "Color setting '" << color << "' is deprecated for heightfield '"
                     << elem_name << "' and will be ignored. Please use a scene graph material instead.";
        }

        // Assign colormap/volume texturing
        if (dict.is_defined("colormap"))
        {
            mi::neuraylib::Tag colormap_tag = resolve_colormap(dict.get("colormap"));
            heightfield->assign_colormap(colormap_tag);
            heightfield->set_colormap_mapping(true);
        }

        // Define a clip area (region of interest)
        if (dict.is_defined("clip"))
        {
            heightfield->set_IJK_region_of_interest(get_bbox_float32_3(dict.get("clip")));
        }

        return heightfield;
    }
    //
    // Pipe set
    //
    else if (elem_type == "pipe_set")
    {
#ifdef PIPE_SET_SUPPORT
        // Pipe dataset file
        const std::string pipe_file = dict.get("input_file");
        const mi::math::Bbox<mi::Float32, 3> pipe_bbox = get_bbox_float32_3(dict.get("box"));

        //
        // Prepare the pipe importer
        //
        mi::neuraylib::Tag import_strategy_tag;

        // Reuse existing importer unless reloading is requested
        if (!check_scope(scope, dict, old_tag))
        {
            ERROR_LOG << "Invalid scope for operation on pipe set " << old_tag;
        }
        else if (old_tag.is_valid() && !get_bool(dict.get("allow_reloading", "no")))
        {
            mi::base::Handle<const nv::index::IPipe_set_scene_element> pipe(
                dice_transaction->access<nv::index::IPipe_set_scene_element>(old_tag));
            if (pipe)
                import_strategy_tag = pipe->get_import_strategy();
        }

        if (!import_strategy_tag.is_valid())
        {
            const std::string pipe_importer = dict.get("importer", "raw");
            if (pipe_importer == "raw")
            {
                // Deprecated: This import scheme doesn't fit for the load balancing.

                // mi::base::Handle<nv::index_common::Pipe_import_config> import_strategy(
                //     new nv::index_common::Pipe_import_config(pipe_file, pipe_bbox));
                // import_strategy_tag = m_dice_transaction->store_for_reference_counting(
                //     import_strategy.get());
            }
            else
            {
                ERROR_LOG << "No such registered importer for pipe: " << pipe_importer;
            }
            assert(import_strategy_tag.is_valid());
        }

        // Create the pipe scene element
        mi::base::Handle<nv::index::IPipe_set_scene_element> pset(
            m_scene->create_pipe_set(import_strategy_tag, pipe_bbox));
        assert(pset != 0);

        return NULL;
#else
        WARN_LOG << "Current version has no IPipe_set_scene_element import support.";
        return null_result();
#endif
    }
    //
    // Triangle mesh
    //
    else if (elem_type == "triangle_mesh")
    {
        mi::base::Handle<nv::index::ITriangle_mesh_scene_element> mesh;

        if (!check_scope(scope, dict, old_tag))
        {
            ERROR_LOG << "Invalid scope for operation on triangle mesh " << old_tag;
        }
        else if (old_tag.is_valid() && !get_bool(dict.get("allow_reloading", "no")))
        {
            // Just edit the existing triangle mesh instead of loading it again
            mesh = dice_transaction->edit<nv::index::ITriangle_mesh_scene_element>(old_tag);
            if (mesh)
                reuse_tag = true;
        }
        else
        {
            // Triangle mesh dataset file
            const std::string mesh_file = dict.get("input_file");
            if(!dict.is_defined("box"))
            {
                ERROR_LOG << "The triangle mesh has no bounding box information (box). The triangle mesh data loading might fail.";
            }

            const mi::math::Bbox<mi::Float32, 3> mesh_bbox = get_bbox_float32_3(dict.get("box"));

            //
            // Prepare the triangle mesh importer
            //
            nv::index::IDistributed_discrete_data_import_callback* importer_callback = NULL;

            const std::string mesh_importer = dict.get("importer", "raw");
            if (mesh_importer == "raw")
            {
                nv::index_common::Triangle_mesh_importer* importer =
                    new nv::index_common::Triangle_mesh_importer(mesh_file);
                importer_callback = importer;
            }
            else
            {
                ERROR_LOG << "No such registered importer for triangle meshes: " << mesh_importer;
            }
            assert(importer_callback!=NULL);

            // Create the triangle mesh scene element
            mesh = m_scene->create_triangle_mesh(mesh_bbox, importer_callback, dice_transaction.get());
            assert(mesh != NULL);
        }

        return mesh;
    }
    //
    // Reservoir grids
    //
    else if (elem_type == "reservoir_grid")
    {
        const mi::math::Bbox<mi::Float32, 3> reservoir_grid_bbox = get_bbox_float32_3(dict.get("box"));

        nv::index::IDistributed_discrete_data_import_callback* importer_callback = NULL;

        const std::string reservoir_grid_importer = dict.get("importer", "artificial");
        if (reservoir_grid_importer == "artificial")
        {
            // Reservoir grid dataset file
            const mi::Uint32 layers = get_uint32(dict.get("layers", "1"));
            const mi::Uint32 cells_per_row = get_uint32(dict.get("cells_per_row", "1"));
            const mi::Uint32 cells_per_column = get_uint32(dict.get("cells_per_column", "1"));
            const std::string cache_file = dict.get("cache_file");
            const std::string cache_type   = dict.get("cache_type");

            nv::index_common::Synthetic_reservoir_patch_generator* generator =
                new nv::index_common::Synthetic_reservoir_patch_generator(
                    layers, cells_per_row, cells_per_column, reservoir_grid_bbox, cache_file,cache_type);
            importer_callback = generator;
        }
        else if (reservoir_grid_importer == "raw")
        {
            const std::string reservoir_file = dict.get("geometry_file");
            const std::string scalar_file = dict.get("scalar_file");
            const std::string cache_type   = dict.get("cache_type");

            if (reservoir_file == "")
                WARN_LOG<<"Missing input file for reservoir grid importer";

            if (scalar_file == "")
                WARN_LOG << "Missing scalar file for reservoir grid importer, "
                         << "cells are colored based on the depth value";

            nv::index_common::Reservoir_grid_importer* importer =
                new nv::index_common::Reservoir_grid_importer(
                    reservoir_file, scalar_file, reservoir_grid_bbox, cache_type);
            importer_callback = importer;
        }
        else
        {
            ERROR_LOG << "No such registered importer for triangle mesh: " << reservoir_grid_importer;
        }
        assert(importer_callback!=NULL);

        // Create the triangle mesh scene element for the reservoir_grid
        mi::base::Handle<nv::index::ITriangle_mesh_scene_element> mesh_reservoir(
            m_scene->create_triangle_mesh(reservoir_grid_bbox, importer_callback, dice_transaction.get()));
        assert(mesh_reservoir != 0);

        return mesh_reservoir;
    }

    //
    // Irregular volumes
    //
    else if (elem_type == "irregular_volume")
    {
        // Unused warning
        // const mi::math::Bbox<mi::Float32, 3> ivol_mesh_bbox = get_bbox_float32_3(dict.get("box"));

        mi::base::Handle<nv::index::IIrregular_volume_scene_element> ivol_mesh;

        if (!check_scope(scope, dict, old_tag))
        {
            ERROR_LOG << "Invalid scope for operation on irregular volume " << old_tag;
        }
        else if (old_tag.is_valid() && !get_bool(dict.get("allow_reloading", "no")))
        {
            // Just edit the existing triangle mesh instead of loading it again
            ivol_mesh = dice_transaction->edit<nv::index::IIrregular_volume_scene_element>(old_tag);
            if (ivol_mesh)
            {
                reuse_tag = true;
            }
        }
        else
        {
            const std::string ivol_mesh_file = dict.get("input_file", "");
            const std::string ivol_scalar_file = dict.get("scalar_file", "");
            const std::string ivol_displ_file = dict.get("displ_file", "");
            const std::string ivol_ts_file = dict.get("ts_file", "");
            const mi::Uint32 ivol_nb_frames = get_uint32(dict.get("gtc2007_converter_nb_frames", "1"));

            if(!dict.is_defined("box"))
            {
                ERROR_LOG << "The irregular volume dataset has no bounding box information (box)."
                          << "The loading process might fail.";
            }

            const mi::math::Bbox<mi::Float32, 3> ivol_mesh_bbox = get_bbox_float32_3(dict.get("box"));

            // Prepare importer
            nv::index::IDistributed_continuous_data_import_callback* importer_callback = NULL;

            const std::string mesh_importer = dict.get("importer", "raw");
            if (mesh_importer == "raw")
            {
                nv::index_common::Irregular_volume_importer* importer = NULL;
                
                // EnSight dataset?
                if(ivol_scalar_file != "")
                {
                    // in converting to ts mode?
                    if(ivol_nb_frames > 1)
                    {
                        importer =
                            new nv::index_common::Irregular_volume_importer(
                                ivol_mesh_file, ivol_scalar_file, ivol_displ_file, ivol_ts_file, ivol_nb_frames);
                    }
                    else
                    {
                        importer =
                            new nv::index_common::Irregular_volume_importer(ivol_mesh_file, ivol_scalar_file);
                    }
                }
                else
                {
                    importer =
                        new nv::index_common::Irregular_volume_importer(ivol_mesh_file);
                }
                
                importer_callback = importer;
            }
            else if (mesh_importer == "vtk")
            {
                importer_callback = new nv::index_common::Irregular_volume_importer_vtk(ivol_mesh_file);
            }
            else
            {
                ERROR_LOG << "No such registered importer for irregular volume datasets: " << mesh_importer;
            }
            assert(importer_callback != NULL);

            // Create the irregular volume scene element
            ivol_mesh = m_scene->create_irregular_volume(ivol_mesh_bbox, importer_callback, dice_transaction.get());
        }

        return ivol_mesh;
    }

    //
    // Sparse volumes
    //
    else if (elem_type == "sparse_volume")
    {
        // Unused warning
        // const mi::math::Bbox<mi::Float32, 3> ivol_mesh_bbox = get_bbox_float32_3(dict.get("box"));

        mi::base::Handle<nv::index::ISparse_volume_scene_element> svol_elem;

        if (old_tag.is_valid() && !get_bool(dict.get("allow_reloading", "no")))
        {
            // Just edit the existing triangle mesh instead of loading it again
            svol_elem = dice_transaction->edit<nv::index::ISparse_volume_scene_element>(old_tag);
            if (svol_elem)
            {
                reuse_tag = true;
            }
        }
        else
        {
            if (!dict.is_defined("box"))
            {
                ERROR_LOG << "The sparse volume dataset has no bounding box information (box)."
                          << "The loading process might fail.";
            }

            const mi::math::Bbox<mi::Float32, 3> svol_bbox = get_bbox_float32_3(dict.get("box"));

            // Prepare importer
            nv::index::IDistributed_discrete_data_import_callback* importer_callback = NULL;

            const std::string svol_importer = dict.get("importer", "vdb");
#ifdef NV_IDX_USE_OPENVDB_IMPORTER
            const bool        svol_use_cache      = get_bool(dict.get("cache", "false"));

            const std::string svol_in_dir         = dict.get("input_directory",      "");
            const std::string svol_in_file_base   = dict.get("input_file_base_name", "");
            const std::string svol_in_file_ext    = dict.get("input_file_extension", "");
            const std::string svol_in_field       = dict.get("input_file_field",     "density");
                                                  
            const std::string svol_out_dir        = dict.get("output_cache_directory", "");

            const mi::Uint32  svol_ts_enum_len    = get_uint32(dict.get("input_time_step_enum_length", "0"));
            const mi::Uint32  svol_ts_enum_stride = get_uint32(dict.get("input_time_step_enum_stride", "1"));
            const mi::Uint32  svol_ts_enum_offset = get_uint32(dict.get("input_time_step_enum_offset", "0"));

            if (svol_importer == "vdb")
            {
                nv::index_common::Sparse_volume_importer* importer = NULL;
                importer =
                    new nv::index_common::Sparse_volume_importer(
                            svol_in_dir,
                            svol_out_dir,
                            svol_in_file_base,
                            svol_in_file_ext,
                            svol_in_field,
                            svol_ts_enum_len,
                            svol_ts_enum_stride,
                            svol_ts_enum_offset,
                            svol_use_cache);

                importer_callback = importer;
            }
            else
#endif
            {
                ERROR_LOG << "No such registered importer for sparse volume datasets: " << svol_importer;
            }
            assert(importer_callback != NULL);

            // Create the sparse volume scene element
            svol_elem = m_scene->create_sparse_volume(svol_bbox, importer_callback, dice_transaction.get());
        }

        return svol_elem;
    }

    return null_result();
}
