/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/

#include "scene_element_importer.h"

#include <nv/index/iclip_region.h>
#include <nv/index/icolormap.h>
#include <nv/index/idepth_offset.h>
#include <nv/index/ifont.h>
#include <nv/index/iintersection_highlighting.h>
#include <nv/index/ilight.h>
#include <nv/index/ilabel.h>
#include <nv/index/imaterial.h>
#include <nv/index/ipath.h>
#include <nv/index/ipipe.h>
#include <nv/index/iregular_volume_rendering_properties.h>
#include <nv/index/iregular_volume_texture.h>
#include <nv/index/irendering_kernel_programs.h>
#include <nv/index/iscene_group.h>
#include <nv/index/ishading_model.h>
#include <nv/index/isparse_volume_rendering_properties.h>
#include <nv/index/itexture.h>
#include <nv/index/itexture_filter_mode.h>
#include <nv/index/itriangle_mesh_scene_element.h>
#include <nv/index/iwireframe_rendering_style.h>

#include "common/common_utility.h"
#include "common/distributed_compute_techniques.h"
#include "common/ppm_io.h"

#include "colormap_util.h"
#include "nvindex_appdata.h"
#include "rtc_parameter_buffer_manip.h"

using namespace nv::index_common;

namespace {

// Colormap generating utility function
void generate_colormap(std::vector<mi::math::Color_struct>& color_entries)
{
   color_entries.resize(256);

   mi::Float32 vmin = 0.0f;
   mi::Float32 vmax = 0.9f;
   for (mi::Uint32 i=0; i < 256; ++i)
   {
       mi::Float32 v = (i + 0.5f) / 256.f;

       mi::math::Color c(1.f, 1.f, 1.f, 1.f);
       mi::Float32 dv;

       if (v < vmin)
          v = vmin;
       if (v > vmax)
          v = vmax;
       dv = vmax - vmin;

       if (v < (vmin + 0.25f * dv))
       {
          c.r = 0.f;
          c.g = 4.f * (v - vmin) / dv;
       }
       else if (v < (vmin + 0.5f * dv))
       {
          c.r = 0.f;
          c.b = 1.f + 4.f * (vmin + 0.25f * dv - v) / dv;
       }
       else if (v < (vmin + 0.75f * dv))
       {
          c.r = 4.f * (v - vmin - 0.5f * dv) / dv;
          c.b = 0.f;
       }
       else
       {
          c.g = 1.f + 4.f * (vmin + 0.75f * dv - v) / dv;
          c.b = 0.f;
       }
       color_entries[i] = c;
    }
}

} //namespace

mi::base::Handle<nv::index::IScene_element> Scene_element_importer::import_attributes(
    const std::string& elem_name,
    const std::string& elem_type,
    const String_dict& dict)
{
    //
    // Material
    //
    if (elem_type == "phong_gl")
    {
        mi::base::Handle<nv::index::IPhong_gl> mat(
            m_scene->create_attribute<nv::index::IPhong_gl>());

        if (dict.is_defined("ambient"))
            mat->set_ambient(get_color(dict.get("ambient")));
        if (dict.is_defined("diffuse"))
            mat->set_diffuse(get_color(dict.get("diffuse")));
        if (dict.is_defined("specular"))
            mat->set_specular(get_color(dict.get("specular")));
        if (dict.is_defined("shininess"))
            mat->set_shininess(get_float32(dict.get("shininess")));
        if (dict.is_defined("opacity"))
            mat->set_opacity(get_float32(dict.get("opacity")));

        return mat;
    }

    //
    // Time step Assignment
    //
    if (elem_type == "time_step_assignment")
    {
        mi::base::Handle<nv::index::ITime_step_assignment> assignment(
            m_scene->create_attribute<nv::index::ITime_step_assignment>());

        if (dict.is_defined("nb_time_steps"))
        {
            const mi::Uint64 nb_time_steps = get_uint64(dict.get("nb_time_steps"));
            assignment->set_nb_time_steps(nb_time_steps);
            ERROR_LOG << "Number of time steps: " << nb_time_steps;
        }
        if (dict.is_defined("interval"))
        {
            const mi::math::Vector<mi::Float32, 2> interval = get_vec_float32_2(dict.get("interval"));
            assignment->set_time_interval(interval.x, interval.y);
            ERROR_LOG << "Time Interval: " << interval;
        }

        return assignment;
    }
    //
    // Lights
    //
    else if (elem_type == "directional_light")
    {
        mi::base::Handle<nv::index::IDirectional_light> light(
            m_scene->create_attribute<nv::index::IDirectional_light>());

        if (dict.is_defined("direction"))
            light->set_direction(get_vec_float32_3(dict.get("direction")));
        if (dict.is_defined("intensity"))
            light->set_intensity(get_color(dict.get("intensity")));

        return light;
    }
    else if (elem_type == "directional_headlight")
    {
        mi::base::Handle<nv::index::IDirectional_headlight> light(
            m_scene->create_attribute<nv::index::IDirectional_headlight>());

        if (dict.is_defined("direction"))
            light->set_direction(get_vec_float32_3(dict.get("direction")));
        if (dict.is_defined("intensity"))
            light->set_intensity(get_color(dict.get("intensity")));

        return light;
    }
    //
    // Per-heightfield geometry style
    //
    else if (elem_type == "heightfield_geometry" || elem_type == "horizon_geometry")
    {
        mi::base::Handle<nv::index::IHeightfield_geometry_settings> cfg(
            m_scene->create_attribute<nv::index::IHeightfield_geometry_settings>());

        if (dict.is_defined("visible"))
            cfg->set_visible(get_bool(dict.get("visible")));

        if (dict.is_defined("mask"))
        {
            mi::Uint32 types = 0;
            std::istringstream types_list(dict.get("mask"));
            while (!types_list.eof())
            {
                std::string type;
                types_list >> type;
                if (type == "isolated_points")
                    types |= nv::index::IHeightfield_geometry_settings::TYPE_ISOLATED_POINTS;
                else if (type == "connecting_lines")
                    types |= nv::index::IHeightfield_geometry_settings::TYPE_CONNECTING_LINES;
                else if (type == "seed_points")
                    types |= nv::index::IHeightfield_geometry_settings::TYPE_SEED_POINTS;
                else if (type == "seed_lines")
                    types |= nv::index::IHeightfield_geometry_settings::TYPE_SEED_LINES;
                else if (type == "all")
                    types |= nv::index::IHeightfield_geometry_settings::TYPE_ALL;
                else if (!type.empty())
                    ERROR_LOG << "Invalid mask type '" << type << "' for element '" << elem_name << "'";
            }
            cfg->set_type_mask(types);
        }

        if (dict.is_defined("color"))
            cfg->set_color(get_color(dict.get("color")));

        if (dict.is_defined("mode"))
        {
            std::string mode = dict.get("mode");
            if (mode == "fixed")
                cfg->set_color_mode(nv::index::IHeightfield_geometry_settings::MODE_FIXED);
            else if (mode == "material")
                cfg->set_color_mode(nv::index::IHeightfield_geometry_settings::MODE_MATERIAL);
            else if (mode == "volume")
                cfg->set_color_mode(nv::index::IHeightfield_geometry_settings::MODE_VOLUME_TEXTURE);
            else if (mode == "computed")
                cfg->set_color_mode(nv::index::IHeightfield_geometry_settings::MODE_COMPUTED_TEXTURE);
            else
                ERROR_LOG << "Invalid mode '" << mode << "' for element '" << elem_name << "'";
        }

        if (dict.is_defined("render"))
        {
            std::string mode = dict.get("render");
            if (mode == "default")
                cfg->set_render_mode(nv::index::IHeightfield_geometry_settings::RENDER_DEFAULT);
            else if (mode == "z-axis" || mode == "flat")
                cfg->set_render_mode(nv::index::IHeightfield_geometry_settings::RENDER_Z_AXIS_ALIGNED);
            else if (mode == "screen")
                cfg->set_render_mode(nv::index::IHeightfield_geometry_settings::RENDER_SCREEN_ALIGNED);
            else if (mode == "rasterized")
                cfg->set_render_mode(nv::index::IHeightfield_geometry_settings::RENDER_RASTERIZED);
            else if (mode == "none")
                cfg->set_render_mode(nv::index::IHeightfield_geometry_settings::RENDER_NONE);
            else
                ERROR_LOG << "Invalid render type '" << mode << "' for element '" << elem_name << "'";
        }

        if (dict.is_defined("offset"))
            cfg->set_render_offset(get_float32(dict.get("offset")));

        if (dict.is_defined("size"))
        {
            mi::Float32 size = get_float32(dict.get("size"));
            cfg->set_geometry_size(size);
            if (cfg->get_geometry_size() != size)
                WARN_LOG <<  "Invalid geometry size " << size << " for element '" << elem_name << "', "
                         << "using size " << cfg->get_geometry_size() << " instead";
        }

        if (dict.is_defined("opacity_ramp_distances"))
        {
            mi::math::Vector<mi::Float32, 2> dists = get_vec_float32_2(dict.get("opacity_ramp_distances"));
            cfg->set_opacity_ramp_distances(dists);
        }

        return cfg;
    }
    //
    // Computed texture for planes and heightfields
    //
    else if (elem_type == "texturing_technique")
    {
        ERROR_LOG << "The 'texturing_technique' compute generator has been removed in favor of the improved "
                  << "'compute_technique' distributed compute generators.";
    }
    //
    // Distributed compute techniques
    //
    else if (elem_type == "compute_technique")
    {
        const mi::math::Vector<mi::Float32, 2> extent = get_vec_float32_2(dict.get("resolution", "100 100"));
        const bool enabled = get_bool(dict.get("enabled", "true"));

        const std::string generator = dict.get("generator", "");
        if (generator == "checkerboard")
        {
            nv::index::IDistributed_compute_destination_buffer_2d_texture::Buffer_format buffer_format
                = nv::index::IDistributed_compute_destination_buffer_2d_texture::RGBA_UINT8;

            const std::string format = dict.get("format", "rgba_float32");
            if (format == "intensity_uint8")
                buffer_format = nv::index::IDistributed_compute_destination_buffer_2d_texture::INTENSITY_UINT8;
            else if (format == "rgba_float32")
                buffer_format = nv::index::IDistributed_compute_destination_buffer_2d_texture::RGBA_FLOAT32;
            else if (format == "rgba_uint8")
                buffer_format = nv::index::IDistributed_compute_destination_buffer_2d_texture::RGBA_UINT8;
            else
                ERROR_LOG << "Invalid format '" << format << "' for element '" << elem_name << "'";

            mi::base::Handle<Distributed_compute_checkerboard_2d> tex(new Distributed_compute_checkerboard_2d(extent, buffer_format));

            if (dict.is_defined("color0"))
                tex->set_color_0(get_color(dict.get("color0")));
            if (dict.is_defined("color1"))
                tex->set_color_1(get_color(dict.get("color1")));

            tex->set_enabled(enabled);

            return tex;
        }
        if (generator == "checkerboard_volume")
        {
            const mi::math::Vector<mi::Uint32, 3> checker_size = get_vec_uint32_3(dict.get("checker_size", "10 10 10"));

            mi::base::Handle<Distributed_compute_checkerboard_3d> tex(new Distributed_compute_checkerboard_3d(checker_size));

            tex->set_enabled(enabled);

            return tex;
        }
        else if (generator == "mandelbrot")
        {
            mi::base::Handle<Distributed_compute_mandelbrot_2d> tex(new Distributed_compute_mandelbrot_2d(extent));

            tex->set_enabled(enabled);

            return tex;
        }
        else if (generator == "bitmap")
        {
            std::string fname = dict.get("filename", "").c_str();
            bool use_cache = get_bool(dict.get("use_cache", "no"));
            
            INFO_LOG << "GEO|Creating bitmap";

            mi::base::Handle<Distributed_compute_bitmap_mapping_2d> tex(new Distributed_compute_bitmap_mapping_2d(fname, use_cache, extent));

            tex->set_enabled(enabled);

            return tex;
        }
        else if (generator == "height_color")
        {
            const mi::math::Vector<mi::Float32, 2> height_range = get_vec_float32_2(
                dict.get("height_range", "0 256"));

            mi::base::Handle<Distributed_compute_height_color_mapping_2d> tex(
                new Distributed_compute_height_color_mapping_2d(height_range));

            tex->set_colors(
                get_color(dict.get("color0", "0 1 0")),
                get_color(dict.get("color1", "1 0 0")));

            tex->set_enabled(enabled);

            return tex;
        }
        else
        {
            ERROR_LOG << "Invalid texture generator '" << generator << "' for element '" << elem_name << "'";
        }
    }
    //
    // Path style
    //
    else if (elem_type == "path_style")
    {
        mi::base::Handle<nv::index::IPath_style>
            path_style(m_scene->create_attribute<nv::index::IPath_style>());

        std::string cap_str = dict.get("cap_style", "flat");
        if(cap_str == "none")
            path_style->set_cap_style(nv::index::IPath_style::CAP_STYLE_NONE);
        else if(cap_str == "flat")
            path_style->set_cap_style(nv::index::IPath_style::CAP_STYLE_FLAT);
        else if(cap_str == "round")
            path_style->set_cap_style(nv::index::IPath_style::CAP_STYLE_ROUND);

        std::string interpolation_str = dict.get("interpolation", "segment");

        if(interpolation_str == "segment")
            path_style->set_interpolation(nv::index::IPath_style::INTERPOLATION_SEGMENT);
        else if(interpolation_str == "nearest")
            path_style->set_interpolation(nv::index::IPath_style::INTERPOLATION_NEAREST);
        else if(interpolation_str == "linear")
            path_style->set_interpolation(nv::index::IPath_style::INTERPOLATION_LINEAR);


        std::string source_str = dict.get("color_source", "both");

        if(source_str == "none")
            path_style->set_color_source(nv::index::IPath_style::COLOR_SOURCE_NONE);
        else if(source_str == "colormap_only")
            path_style->set_color_source(nv::index::IPath_style::COLOR_SOURCE_COLORMAP_ONLY);
        else if(source_str == "rgba_only")
            path_style->set_color_source(nv::index::IPath_style::COLOR_SOURCE_RGBA_ONLY);
        else if(source_str == "both")
            path_style->set_color_source(nv::index::IPath_style::COLOR_SOURCE_BOTH);

        if (dict.is_defined("upsampling"))
        {
            const mi::Uint32 upsampling = get_uint32(dict.get("upsampling", "2"));
            mi::Float32 tension = get_float32(dict.get("tension", "0.0"));

            path_style->set_upsampling(true, upsampling, tension);
        }

        return path_style;
    }
    //
    // Clipping region
    //
    else if (elem_type == "clip_region")
    {
        mi::base::Handle<nv::index::IClip_region> attr(
            m_scene->create_attribute<nv::index::IClip_region>());

        if (dict.is_defined("bbox"))
        {
            const mi::math::Bbox<mi::Float32, 3> bbox = get_bbox_float32_3(dict.get("bbox"));
            attr->set_clip_bounding_box(bbox);
        }

        return attr;
    }
    //
    // Depth offset
    //
    else if (elem_type == "depth_offset")
    {
        mi::base::Handle<nv::index::IDepth_offset> depth_offset_attr(
            m_scene->create_attribute<nv::index::IDepth_offset>());

        mi::Float32 depth_offset = get_float32(dict.get("offset", "0"));
        depth_offset_attr->set_depth_offset(depth_offset);

        return depth_offset_attr;
    }
    //
    // Intersection highlighting
    //
    else if (elem_type == "intersection_highlighting")
    {
        mi::base::Handle<nv::index::IIntersection_highlighting> attr(
            m_scene->create_attribute<nv::index::IIntersection_highlighting>());

        if (dict.is_defined("intersect_with"))
        {
            std::istringstream is(dict.get("intersect_with"));
            mi::neuraylib::Tag intersect_tag;
            if (is >> intersect_tag.id)
            {
                // Interpret numerical value as tag
                attr->set_intersection_shape(intersect_tag);
            }
            else
            {
                // It must be a scene element id, the reference will be resolved later
                m_referencing_elements.push_back(elem_name);
            }
        }

        if (dict.is_defined("color"))
            attr->set_color(get_color(dict.get("color")));

        if (dict.is_defined("width"))
            attr->set_width(get_float32(dict.get("width")));

        if (dict.is_defined("smooth"))
            attr->set_smoothness(get_float32(dict.get("smooth")));

        if (dict.is_defined("discontinuity_limit"))
            attr->set_discontinuity_limit(get_float32(dict.get("discontinuity_limit")));

        return attr;
    }
    //
    // Regular volume texturing
    //
    else if (elem_type == "volume_texture")
    {
        mi::base::Handle<nv::index::IRegular_volume_texture> attr(
            m_scene->create_attribute<nv::index::IRegular_volume_texture>());

        if (dict.is_defined("volume"))
        {
            std::istringstream is(dict.get("volume"));
            mi::neuraylib::Tag rvol_tag;
            if (is >> rvol_tag.id)
            {
                // Interpret numerical value as tag
                attr->set_volume_element(rvol_tag);
            }
            else
            {
                // It must be a scene element id, the reference will be resolved later
                m_referencing_elements.push_back(elem_name);
            }
        }

        const std::string mode = dict.get("boundary_mode", "clamp_to_transparent");
        if (mode == "clamp_to_transparent")
        {
            attr->set_volume_boundary_mode(nv::index::IRegular_volume_texture::CLAMP_TO_TRANSPARENT);
        }
        else if (mode == "clamp_to_opaque")
        {
            attr->set_volume_boundary_mode(nv::index::IRegular_volume_texture::CLAMP_TO_OPAQUE);
        }
        else
        {
            ERROR_LOG << "Invalid boundary_mode '" << mode << "' for element " << elem_name;
        }

        return attr;
    }
    //
    // Wireframe style
    //
    else if (elem_type == "wireframe_rendering_style")
    {
        std::string wire_mode_str = dict.get("style", "outline");
        mi::Float32 wire_width = get_float32(dict.get("width", "1.0"));

        mi::base::Handle<nv::index::IWireframe_rendering_style> wireframe(
            m_scene->create_attribute<nv::index::IWireframe_rendering_style>());

        wireframe->set_wireframe_width(wire_width);
        wireframe->set_wireframe_color(get_color(dict.get("color", "1 0 0 1")));

        if(wire_mode_str == "outline")
        {
            wireframe->set_wireframe_style(nv::index::IWireframe_rendering_style::WIREFRAME_STYLE_OUTLINE);
            return wireframe;
        }
        else if(wire_mode_str == "wireframe")
        {
            wireframe->set_wireframe_style(nv::index::IWireframe_rendering_style::WIREFRAME_STYLE_WIREFRAME);
            return wireframe;
        }
        else
        {
            ERROR_LOG << "Invalid wireframe style: " << wire_mode_str;
        }
    }
    //
    // Heightfield wireframe style
    //
    else if (elem_type == "heightfield_wireframe_style")
    {
        bool error = false;

        std::string wire_mode_str = dict.get("style", "outline");
        std::string wire_topo_str = dict.get("topology", "quads");
        mi::Float32 wire_width = get_float32(dict.get("width", "1.0"));
        mi::Uint32 wire_resolution = get_uint32(dict.get("resolution", "1"));

        mi::base::Handle<nv::index::IHeightfield_wireframe_style> wireframe(
            m_scene->create_attribute<nv::index::IHeightfield_wireframe_style>());

        wireframe->set_wireframe_resolution(wire_resolution);
        wireframe->set_wireframe_width(wire_width);
        wireframe->set_wireframe_color(get_color(dict.get("color", "1 0 0 1")));

        if(wire_mode_str == "outline")
        {
            wireframe->set_wireframe_style(nv::index::IWireframe_rendering_style::WIREFRAME_STYLE_OUTLINE);
        }
        else if(wire_mode_str == "wireframe")
        {
            wireframe->set_wireframe_style(nv::index::IWireframe_rendering_style::WIREFRAME_STYLE_WIREFRAME);
        }
        else
        {
            ERROR_LOG << "Invalid wireframe style: " << wire_mode_str;
            error = true;
        }

        if(wire_topo_str == "triangles")
        {
            wireframe->set_wireframe_topology(
                nv::index::IHeightfield_wireframe_style::WIREFRAME_TOPOLOGY_TRIANGLES);
        }
        else if(wire_topo_str == "quads")
        {
            wireframe->set_wireframe_topology(
                nv::index::IHeightfield_wireframe_style::WIREFRAME_TOPOLOGY_QUADS);
        }
        else
        {
            ERROR_LOG << "Invalid wireframe topology: " << wire_topo_str;
            error = true;
        }

        if(!error)
        {
            return wireframe;
        }
    }
    //
    // Shading mode
    //
    else if (elem_type == "shading_model")
    {
        const bool two_sided_lighting = get_bool(dict.get("two_sided_lighting", "off"));
        const bool backface_culling   = get_bool(dict.get("backface_culling", "off"));

        nv::index::IShading_model::Front_face_vertex_order ff_order = nv::index::IShading_model::FRONT_FACE_CCW;
        const std::string ff_order_entry = dict.get("front_face", "ccw");
        if (ff_order_entry == "ccw") {
            // redundant
        }
        else if (ff_order_entry == "cw") {
            ff_order = nv::index::IShading_model::FRONT_FACE_CW;
        }
        else {
            ERROR_LOG << "No such face vertex-ordering mode for shading_model: " << ff_order_entry;
        }

        const std::string mode = dict.get("mode", "phong");
        if (mode == "phong")
        {
            mi::base::Handle<nv::index::IPhong_shading> shading_mode(
                m_scene->create_attribute<nv::index::IPhong_shading>());

            shading_mode->set_double_sided_lighting(two_sided_lighting);
            shading_mode->set_backface_culling(backface_culling);
            shading_mode->set_front_face_vertex_order(ff_order);

            return shading_mode;
        }
        else if (mode == "flat")
        {
            mi::base::Handle<nv::index::IFlat_shading> shading_mode(
                m_scene->create_attribute<nv::index::IFlat_shading>());

            shading_mode->set_double_sided_lighting(two_sided_lighting);
            shading_mode->set_backface_culling(backface_culling);
            shading_mode->set_front_face_vertex_order(ff_order);

            return shading_mode;
        }
        else
        {
            ERROR_LOG << "No such registered mode for shading_model: " << mode;
        }
    }
    // Texture filter mode
    else if (elem_type == "texture_filter_mode")
    {
        const std::string   tf_mode = dict.get("mode", "nearest");
        if (tf_mode == "nearest")
        {
            mi::base::Handle<nv::index::ITexture_filter_mode_nearest_neighbor> filter_mode(
                m_scene->create_attribute<nv::index::ITexture_filter_mode_nearest_neighbor>());

            return filter_mode;
        }
        else if (tf_mode == "linear")
        {
            mi::base::Handle<nv::index::ITexture_filter_mode_linear> filter_mode(
                m_scene->create_attribute<nv::index::ITexture_filter_mode_linear>());

            return filter_mode;
        }
    }
    else if (elem_type == "volume_rendering_properties")
    {
        const std::string   shade_mode = dict.get("mode", "no_lighting");

        mi::base::Handle<nv::index::IRegular_volume_rendering_properties> render_prop(
            m_scene->create_attribute<nv::index::IRegular_volume_rendering_properties>());

        if (shade_mode == "no_lighting")
        {
            render_prop->set_shading_mode(nv::index::IRegular_volume_rendering_properties::NO_LIGHTING);
        }
        else if (shade_mode == "phong_lighting")
        {
            render_prop->set_shading_mode(nv::index::IRegular_volume_rendering_properties::PHONG_MATERIAL_LIGHTING);
        }

        const mi::Float32 shading_grad_threshold = get_float32(dict.get("gradient_threshold", "0.0"));
        render_prop->set_shading_gradient_threshold(shading_grad_threshold);

        const mi::Float32 reference_step_size = get_float32(dict.get("reference_step_size", "1.0"));
        render_prop->set_reference_step_size(reference_step_size);

        return render_prop;
    }
    else if (elem_type == "sparse_volume_rendering_properties")
    {
        const std::string   filter_mode = dict.get("filter_mode", "nearest");

        mi::base::Handle<nv::index::ISparse_volume_rendering_properties> render_prop(
            m_scene->create_attribute<nv::index::ISparse_volume_rendering_properties>());

        if (filter_mode == "nearest")
        {
            render_prop->set_filter_mode(nv::index::SPARSE_VOLUME_FILTER_NEAREST);
        }
        else if (filter_mode == "trilinear")
        {
            render_prop->set_filter_mode(nv::index::SPARSE_VOLUME_FILTER_TRILINEAR_POST);
        }

        render_prop->set_sampling_distance(             get_float32(dict.get("sampling_distance",               "1.0")));
        render_prop->set_preintegrated_volume_rendering(get_bool(   dict.get("preintegrated_volume_rendering",  "false")));
        render_prop->set_debug_visualization_option(    get_uint32( dict.get("debug_vis_option",                "0")));

        return render_prop;
    }
    else if (elem_type == "rendering_kernel_program")
    {
        mi::base::Handle<nv::index::IRendering_kernel_program> rtc_program;

        const std::string rtc_target = dict.get("target", "volume_sample_program");
        if (rtc_target == "volume_sample_program")
        {
            rtc_program = m_scene->create_attribute<nv::index::IVolume_sample_program>();
        }
        else if (rtc_target == "surface_sample_program")
        {
            rtc_program = m_scene->create_attribute<nv::index::ISurface_sample_program>();
        }
        else
        {
            ERROR_LOG << "Invalid target '" << rtc_target << "' for element " << elem_name;
            return null_result();
        }

        std::string rtc_program_source;
        const std::string rtc_program_source_file = dict.get("source_file");
        if (!rtc_program_source_file.empty())
        {
            // load from file
            INFO_LOG << "Reading rendering kernel program source for element '" << elem_name << "' from file '"
                     << rtc_program_source_file << "'...";
            if (!nv::index_common::read_file(rtc_program_source_file, rtc_program_source))
            {
                ERROR_LOG << "Error reading file '" << rtc_program_source_file << "'.";
                return null_result();
            }
        }
        else
        {
            rtc_program_source = dict.get("source_string");
        }

        rtc_program->set_program_source(rtc_program_source.c_str());
        return rtc_program;
    }
    else if (elem_type == "rendering_kernel_program_parameters")
    {
        mi::base::Handle<nv::index::IRendering_kernel_program_parameters> rtc_params(
            m_scene->create_attribute<nv::index::IRendering_kernel_program_parameters>());

        RTC_parameter_buffer_manip::instance()->add_parameter_buffer_data(rtc_params);

        return rtc_params;
    }

    return null_result();
}

mi::neuraylib::Tag Scene_element_importer::import_colormaps(
    const std::string& elem_name,
    const std::string& elem_type,
    const String_dict& dict)
{
    mi::neuraylib::Tag new_tag;

    if (elem_type == "colormap")
    {
        const std::string type_str = dict.get("map_type");
        if (type_str == "procedural")
        {
            mi::base::Handle<nv::index::IColormap> colormap(
                m_scene->create_attribute<nv::index::IColormap>());

            const std::string map_name_str = dict.get("map_name", "jet");
            const mi::Uint32 map_size = get_uint32(dict.get("map_size", "128"));

            std::vector<mi::math::Color_struct> map;
            for (mi::Uint32 i=0; i<map_size; i++)
            {
                mi::Float32 t = (i + 0.5f)/map_size;

                if(map_name_str == "jet")
                    map.push_back(jetmap(t, 0.f, 1.f));
            }

            colormap->set_colormap(&map[0], map_size);

            new_tag = m_dice_transaction->store_for_reference_counting(
                colormap.get(), mi::neuraylib::NULL_TAG, elem_name.c_str());
        }
        else if (type_str == "lookup_table")
        {
            // assign existing colormap
            mi::neuraylib::Tag colormap_tag;
            if (dict.is_defined("map_index"))
            {
                colormap_tag = get_colormap_tag(get_uint32(dict.get("map_index")));
            }
            else
            {
                colormap_tag = get_colormap_tag(Nvindex_AppData::instance()->get_current_colormap_index());
            }
            new_tag = colormap_tag;
        }
        else if (type_str == "data")
        {
            new_tag = resolve_colormap(elem_name);
        }
        else if (type_str == "reference")
        {
            // This type is only used when exporting a single scene element, i.e. for the scene
            // editor. Just reuse the existing colormap with the given tag.
            new_tag = m_dice_transaction->name_to_tag(dict.get("colormap").c_str());
        }
        else
        {
            ERROR_LOG << "Invalid colormap map_type: '" << type_str << "' for element " << elem_name;
        }

        if (new_tag)
        {
            mi::base::Handle<nv::index::IColormap> colormap(
                m_dice_transaction->edit<nv::index::IColormap>(new_tag));

            // Need to set enabled state here because colormaps are handled differently than other
            // scene elements
            colormap->set_enabled(get_bool(dict.get("enabled", "true")));

            // Use colormap domain settings
            const mi::math::Vector<mi::Float32, 2> domain = get_vec_float32_2(dict.get("domain", "0 1"));
            colormap->set_domain(domain.x, domain.y);

            const std::string mode = dict.get("domain_boundary_mode", "clamp_to_edge");
            if (mode == "clamp_to_edge")
            {
                colormap->set_domain_boundary_mode(nv::index::IColormap::CLAMP_TO_EDGE);
            }
            else if (mode == "clamp_to_transparent")
            {
                colormap->set_domain_boundary_mode(nv::index::IColormap::CLAMP_TO_TRANSPARENT);
            }
            else
            {
                ERROR_LOG << "Invalid domain_boundary_mode '" << mode << "' for element " << elem_name;
            }
        }
    }
    else if (elem_type == "reservoir_colormap")
    {
        const std::string mode = dict.get("mode", "artificial");
        if (mode == "artificial")
        {
            nv::index::IColormap* colormap_element = m_scene->create_attribute<nv::index::IColormap>();
            std::vector<mi::math::Color_struct> colormap_entries;
            generate_colormap(colormap_entries);
            // Set the the colormap values.
            colormap_element->set_colormap(&(colormap_entries[0]), colormap_entries.size());
            mi::base::Handle<nv::index::IColormap> colormap(colormap_element);

            new_tag = m_dice_transaction->store_for_reference_counting(
                colormap.get(), mi::neuraylib::NULL_TAG, "reservoir_colormap");
        }
        else if (mode == "file")
        {
            // Use 39th ppm colormap from the GUI, this is the most colorful one we have
            // This is fixed simply because we can't change the colors on the fly anyway.
            mi::neuraylib::Tag colormap_tag = get_colormap_tag(mi::Uint32(39));
            new_tag = colormap_tag;
        }
        else
        {
            WARN_LOG << "Unsupported colormap mode '" << mode << "' for the reservoir";
        }
    }

    return new_tag;
}

mi::base::Handle<nv::index::IScene_element> Scene_element_importer::import_raster_attributes(
    const std::string& elem_name,
    const std::string& elem_type,
    const String_dict& dict)
{
    //
    // Label attributes
    //
    if (elem_type == "label_font")
    {
        mi::base::Handle<nv::index::IFont> font(
            m_scene->create_attribute<nv::index::IFont>());

        font->set_file_name(dict.get("filename", "").c_str());

        if (dict.is_defined("resolution"))
            font->set_font_resolution(get_float32(dict.get("resolution")));

        if (dict.is_defined("size"))
            WARN_LOG << "font::size is deprecated. Please use font::resolution.";

        return font;
    }
    else if (elem_type == "label_layout")
    {
        mi::base::Handle<nv::index::ILabel_layout> layout(
            m_scene->create_attribute<nv::index::ILabel_layout>());

        if (dict.is_defined("padding"))
            layout->set_padding(get_float32(dict.get("padding")));

        if (dict.is_defined("auto_flip"))
            layout->set_auto_flip(get_bool(dict.get("auto_flip")));

        layout->set_color(
            get_color(dict.get("foreground", "1 1 1 1")),
            get_color(dict.get("background", "0 0 0 1")));

        return layout;
    }
    //
    // Texture
    //
    else if (elem_type == "texture")
    {
        mi::base::Handle<nv::index::ITexture> tex(
            m_scene->create_attribute<nv::index::ITexture>());

        // load image file
        const std::string fname = dict.get("filename", "");
        if (fname.empty()) // create a test procedural texture
        {
            const mi::math::Vector<mi::Uint32, 2> res
                = get_vec_uint32_2(dict.get("resolution", "100 100"));
            tex->set_pixel_data(NULL, res.x, res.y, nv::index::ITexture::RGBA_UINT8);
        }
        else
        {
            std::vector< mi::math::Color_struct > ppm_color_st_buf;
            mi::Sint32 img_width  = -1;
            mi::Sint32 img_height = -1;
            std::string error_mes;
            if(!load_ppm(fname, ppm_color_st_buf, img_width, img_height, error_mes))
            {
                ERROR_LOG << "Texture failed to load [" << fname << "]: " << error_mes;
                return null_result();
            }
            tex->set_pixel_data(&ppm_color_st_buf[0], img_width, img_height, nv::index::ITexture::RGBA_FLOAT32);
        }

        return tex;
    }

    return null_result();
}
