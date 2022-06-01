/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief opengl drawing utility for Non-OpenGL (stubs)

#include "opengl_drawing_utilities.h"

#include <cassert>


//----------------------------------------------------------------------
void gl_initialize_opengl()
{
    // empty: stub
}

//----------------------------------------------------------------------
std::string gl_error()
{
    return "";
}

//----------------------------------------------------------------------
void gl_initialize_glew()
{
    // empty: stub
}

//----------------------------------------------------------------------
void gl_init_display_lists()
{
    // empty: stub
}

//----------------------------------------------------------------------
void gl_debug_draw_viewport(mi::Sint32 viewport_width, mi::Sint32 viewport_height)
{
    // empty: stub
}

//----------------------------------------------------------------------
void gl_debug_render_bbox(const mi::math::Bbox< mi::Float32, 3 >& bbox)
{
    // empty: stub
}

//----------------------------------------------------------------------
void gl_show_icamera_coordinate_system(
    bool show_coordinate_system,
    mi::math::Matrix<mi::Float32, 4, 4> const & mat)
{
    // empty: stub
}

//----------------------------------------------------------------------
// draw axis
void gl_draw_axis(
    const mi::neuraylib::Tag &        camera_tag,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    // empty: stub
}

//----------------------------------------------------------------------
void gl_display_colormap(
    const mi::neuraylib::Tag&          colormap_tag,
    mi::neuraylib::IDice_transaction * dice_transaction)
{
    // empty: stub
}

//----------------------------------------------------------------------
void gl_render_color_table(
    const mi::neuraylib::Tag&          colormap_tag,
    mi::neuraylib::IDice_transaction * dice_transaction)
{
    // empty: stub
}

//----------------------------------------------------------------------
void draw_host_workload(
    bool                      show_large,
    const std::vector<int>&   per_host_rendered_cubes,
    const std::vector<float>& per_host_time_rendering,
    const std::vector<float>& per_host_time_horizons,
    const std::vector<float>& per_host_time_volume,
    const std::vector<float>& per_host_time_volume_and_horizons,
    const std::vector<bool>&  per_host_using_cpu,
    const std::vector<bool>&  per_host_streaming)
{
    // empty: stub
}

//----------------------------------------------------------------------
void draw_rendering_workload(
    mi::Uint32                              nb_cluster_nodes,
    mi::Float32                             scale,
    nv::index::IPerformance_values*    perf_values
    )
{
    // empty: stub
}

//----------------------------------------------------------------------
void draw_compositing_workload(
    mi::Uint32                              nb_horizontal_spans,
    mi::Float32                             scale,
    nv::index::IPerformance_values*    perf_values
    )
{
    // empty: stub
}

//----------------------------------------------------------------------
void draw_performance(
    mi::Float32                             scale,
    nv::index::IPerformance_values*    perf_values
    )
{
    // empty: stub
}

//----------------------------------------------------------------------
void gl_initialize_viewport(
        mi::Sint32 width,
        mi::Sint32 height,
        const mi::math::Color_struct&   background_color,
        nv::index_common::String_dict const & opt)
{
    // empty: stub
}

//----------------------------------------------------------------------
void gl_setup_projection_matrix(
    const nv::index::ICamera* cam)
{
    // empty: stub
}

//----------------------------------------------------------------------
void gl_setup_lookat(mi::math::Vector<mi::Float32, 3> const & eye,
                     mi::math::Vector<mi::Float32, 3> const & target,
                     mi::math::Vector<mi::Float32, 3> const & up)
{
    // empty: stub
}

//----------------------------------------------------------------------
void gl_set_world_transform(
    mi::math::Matrix<mi::Float32, 4, 4> const & transform_mat)
{
    // empty: stub
}

//----------------------------------------------------------------------
void gl_visualize_roi(
    bool                                            is_show_roi,
    mi::math::Bbox_struct<mi::Float32, 3> const &   roi_bbox,
    const mi::math::Color&                          color,
    const mi::math::Matrix<mi::Float32, 4, 4>&      matrix)
{
    // empty: stub
}

//----------------------------------------------------------------------
void gl_visualize_bbox(
    bool                                            is_show_bbox,
    mi::math::Bbox_struct<mi::Float32, 3> const &   bbox,
    const mi::math::Color&                          color,
    const mi::math::Matrix<mi::Float32, 4, 4>&      matrix)
{
    // empty: stub
}

//----------------------------------------------------------------------
void gl_draw_slice(const mi::math::Vector<mi::Float32, 3>&    llf,
                   const mi::math::Vector<mi::Float32, 3>&    urb,
                   const mi::math::Matrix<mi::Float32, 4, 4>& transform)
{
    // empty: stub
}

//----------------------------------------------------------------------
void gl_draw_profile(mi::math::Vector_struct<mi::Float32, 2> const & v0,
                     mi::math::Vector_struct<mi::Float32, 2> const & v1,
                     mi::math::Bbox_struct<mi::Float32, 3>   const & ijk_space,
                     mi::math::Matrix<mi::Float32, 4, 4>     const & transform)
{
    // empty: stub
}

//----------------------------------------------------------------------
void gl_setup_2d_rendering(mi::Float32 width, mi::Float32 height)
{
    // empty: stub
}

//----------------------------------------------------------------------
void gl_visualize_horizontal_span(bool                             is_visualize,
                                  Span_renderer_IF*                p_span_buffer_renderer,
                                  nv::index::IPerformance_values*  perf_values)
{
    // empty: stub
}

//----------------------------------------------------------------------
void gl_render_color_table(
    bool is_show_color_table,
    const mi::neuraylib::Tag& colormap_tag,
    mi::neuraylib::IDice_transaction * dice_transaction)
{
    // empty: stub
}

//----------------------------------------------------------------------
void gl_push_matrix_state()
{
    // empty: stub
}

//----------------------------------------------------------------------
void gl_pop_matrix_state()
{
    // empty: stub
}

//----------------------------------------------------------------------
void gl_initialize_offscreen_context()
{
    // empty: stub
}

//----------------------------------------------------------------------
void gl_delete_offscreen_context(Offscreen_context*)
{
    // empty: stub
}

//----------------------------------------------------------------------
void gl_draw_heightfield_delete_polygon(
    const mi::math::Matrix<mi::Float32, 4, 4>& ijk_to_xyz_space)
{
    // empty: stub
}

//----------------------------------------------------------------------
