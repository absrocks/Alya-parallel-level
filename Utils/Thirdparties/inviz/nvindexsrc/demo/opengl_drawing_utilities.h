/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief opengl drawing utility

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_OPENGL_OPENGL_DRAWING_UTILITIES_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_OPENGL_OPENGL_DRAWING_UTILITIES_H

#include <mi/neuraylib/dice.h>
#include <mi/math/bbox.h>
#include <mi/math/color.h>
#include <mi/base/types.h>

#include <vector>
#include <string>

//
// Forward declarations
namespace nv {

namespace index_common {
class String_dict;
}

namespace index {
class ICamera;
class IPerformance_values;
class ITriangle_mesh_scene_element;
}

} // namespace nv

class Span_renderer_IF;
class Offscreen_context;

//----------------------------------------------------------------------
/// initialize OpenGL
void gl_initialize_opengl();

/// returns OpenGL error string or empty string
std::string gl_error();

/// clear GL buffer.
/// The color buffer is cleared by the current background color, depth
/// buffer is also cleared.
void gl_clear_buffer();

/// initialize glew
void gl_initialize_glew();

/// initialize display list
void gl_init_display_lists();

/// draw viewport for debug
void gl_debug_draw_viewport(mi::Sint32 viewport_width, mi::Sint32 viewport_height);

/// render bounding box without translation and scaling for debug
void gl_debug_render_bbox(const mi::math::Bbox< mi::Float32, 3 >& bbox);

/// show world coordinate system
void gl_show_icamera_coordinate_system(
    bool                                       show_coordinate_system,
    const mi::math::Matrix<mi::Float32, 4, 4>& mat);


/// draw axis
void gl_draw_axis(
    const mi::neuraylib::Tag &        camera_tag,
    mi::neuraylib::IDice_transaction* dice_transaction);

/// display colormap
void gl_display_colormap(
    const mi::neuraylib::Tag&          colormap_tag,
    mi::neuraylib::IDice_transaction * dice_transaction);

/// render color table
void gl_render_color_table(
    const mi::neuraylib::Tag&          colormap_tag,
    mi::neuraylib::IDice_transaction * dice_transaction);

/// draw host workload
void draw_host_workload(
    bool                      show_large,
    const std::vector<int>&   per_host_rendered_cubes,
    const std::vector<float>& per_host_time_rendering,
    const std::vector<float>& per_host_time_horizons,
    const std::vector<float>& per_host_time_seismic,
    const std::vector<float>& per_host_time_seismic_and_horizons,
    const std::vector<bool>&  per_host_using_cpu,
    const std::vector<bool>&  per_host_streaming);

void draw_rendering_workload(
    mi::Uint32                              nb_cluster_nodes,
    mi::Float32                             scale,
    nv::index::IPerformance_values*         perf_values);

void draw_compositing_workload(
    mi::Uint32                              nb_horizontal_spans,
    mi::Float32                             scale,
    nv::index::IPerformance_values*         perf_values);

void draw_performance(
    mi::Float32                             scale,
    nv::index::IPerformance_values*         perf_values);

//----------------------------------------------------------------------
/// initialize OpenGL viewport. Note viewport origin is always 0,0.
///
/// \param[in] width  viewport width
/// \param[in] height viewport height
/// \param[in] background_color background color
void gl_initialize_viewport(
    mi::Sint32                    width,
    mi::Sint32                    height,
    const mi::math::Color_struct& background_color);

//----------------------------------------------------------------------
/// setup OpenGL projection matrix
/// \param[in] cam     Perspective or orthographic camera
void gl_setup_projection_matrix(
    const nv::index::ICamera* cam);

//----------------------------------------------------------------------
/// setup OpenGL lookat camera matrix
///
/// \param[in] eye    eye position
/// \param[in] target target position
/// \param[in] up     up vector
void gl_setup_lookat(const mi::math::Vector<mi::Float32, 3>& eye,
                     const mi::math::Vector<mi::Float32, 3>& target,
                     const mi::math::Vector<mi::Float32, 3>& up);

//----------------------------------------------------------------------
/// set OpenGL world scaling and translation
/// \param[in] transform_mat transformation matrix for the scene
void gl_set_world_transform(
    const mi::math::Matrix<mi::Float32, 4, 4>& transform_mat);

//----------------------------------------------------------------------
/// visualize region of interest by OpenGL
///
/// \param[in] is_show_roi show region of interest when true
/// \param[in] roi_bbox    region of interest bounding box
/// \param[in] color       ROI bounding box color
/// \param[in] transform_mat transformation matrix applied before draw
void gl_visualize_roi(
    bool                                         is_show_roi,
    const mi::math::Bbox_struct<mi::Float32, 3>& roi_bbox,
    const mi::math::Color&                       color         = mi::math::Color(0.4f, 0.4f, 0.4f, 0.4f),
    const mi::math::Matrix<mi::Float32, 4, 4>&   transform_mat = mi::math::Matrix<mi::Float32, 4, 4>(1.f));

//----------------------------------------------------------------------
/// visualize volume bbox
///
/// \param[in] is_show_bbox show bounding box when true
/// \param[in] volume_bbox  volume bounding box
/// \param[in] color        ROI bounding box color
/// \param[in] transform_mat transformation matrix applied before draw
void gl_visualize_bbox(
    bool                                         is_show_bbox,
    const mi::math::Bbox_struct<mi::Float32, 3>& volume_bbox,
    const mi::math::Color&                       color         = mi::math::Color(0.4f, 0.4f, 0.4f, 0.4f),
    const mi::math::Matrix<mi::Float32, 4, 4>&   transform_mat = mi::math::Matrix<mi::Float32, 4, 4>(1.f));

//----------------------------------------------------------------------
/// draw slice by OpenGL
///
/// \param[in] llf lower left front
/// \param[in] urb up    right back
/// \param[in] transform_mat transformation matrix applied before draw
void gl_draw_slice(const mi::math::Vector<mi::Float32, 3>&    llf,
                   const mi::math::Vector<mi::Float32, 3>&    urb,
                   const mi::math::Matrix<mi::Float32, 4, 4>& transform_mat);

//----------------------------------------------------------------------
/// draw profile by OpenGL
///
/// \param[in] v0
/// \param[in] v1
/// \param[in] ijk_space
/// \param[in] transform_mat transformation matrix applied before draw
void gl_draw_profile(const mi::math::Vector_struct<mi::Float32, 2>& v0,
                     const mi::math::Vector_struct<mi::Float32, 2>& v1,
                     const mi::math::Bbox_struct<mi::Float32, 3>&   ijk_space,
                     const mi::math::Matrix<mi::Float32, 4, 4>&     transform_mat);

//----------------------------------------------------------------------
/// set up 2d rendering for OpenGL
///
/// \param[in] width  window width  (resolution x)
/// \param[in] height window height (resolution y)
void gl_setup_2d_rendering(mi::Float32 width, mi::Float32 height);

//----------------------------------------------------------------------
/// visualize horizontal span by OpenGL
///
/// \param[in] is_visualize visualize horizontal span when true
/// \param[in] p_span_buffer_renderer span buffer renderer
/// \param[in] perf_values   performance values for measurement
void gl_visualize_horizontal_span(
    bool                            is_visualize,
    Span_renderer_IF*               p_span_buffer_renderer,
    nv::index::IPerformance_values* perf_values);

//----------------------------------------------------------------------
/// render color table
///
/// \param[in] is_show_color_table show color table when true,
/// otherwise, show colormap
/// \param[in] colormap_tag colormap tag to show
/// \param[in] dice_transaction dice transaction
void gl_render_color_table(
    bool                               is_show_color_table,
    const mi::neuraylib::Tag&          colormap_tag,
    mi::neuraylib::IDice_transaction * dice_transaction);

//----------------------------------------------------------------------
/// push the matrix state
void gl_push_matrix_state();

//----------------------------------------------------------------------
/// pop the matrix state
void gl_pop_matrix_state();

//----------------------------------------------------------------------
/// initialize offscreen context
/// Note: implementation is in the opengl_offscreen_context.cpp
void gl_initialize_offscreen_context();


//----------------------------------------------------------------------
/// delete offscreen context
/// Note: implementation is in the opengl_offscreen_context.cpp
void gl_delete_offscreen_context(Offscreen_context* p_context);

//----------------------------------------------------------------------
/// get perspective parameters from OpenGL projection matrix
///
/// \param[out] mat_fovy_rad calculated fovy (in radian)
/// \param[out] mat_aspect   calculated aspect ratio
/// \param[out] mat_clip_min calculated clipping min plane distance
/// \param[out] mat_clip_max calculated clipping max plane distance
void gl_get_perspective_matrix_parameters(mi::Float32 & mat_fovy_rad,
                                          mi::Float32 & mat_aspect,
                                          mi::Float32 & mat_clip_min,
                                          mi::Float32 & mat_clip_max);

//----------------------------------------------------------------------
/// get modelview parameters from OpenGL projection matrix
///
/// \param[out] mat_eyepos   eye position
/// \param[out] mat_viewdir  viewing direction
/// \param[out] mat_updir    up direction vector
void gl_get_modelview_matrix_parameters(mi::math::Vector<mi::Float32, 3> & mat_eyepos,
                                        mi::math::Vector<mi::Float32, 3> & mat_viewdir,
                                        mi::math::Vector<mi::Float32, 3> & mat_updir);

//----------------------------------------------------------------------
/// draw heightfield delete polygon by OpenGL
///
/// \param[in] ijk_to_xyz_space transformation from IJK local space to XYZ global space
///
void gl_draw_heightfield_delete_polygon(
    const mi::math::Matrix<mi::Float32, 4, 4>& ijk_to_xyz_space);

//----------------------------------------------------------------------

#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_OPENGL_OPENGL_DRAWING_UTILITIES_H
