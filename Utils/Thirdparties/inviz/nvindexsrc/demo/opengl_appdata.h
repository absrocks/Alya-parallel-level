/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief opengl application side data. Independent from OpenGL library.

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_OPENGL_OPENGL_APPDATA_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_OPENGL_OPENGL_APPDATA_H

#include <mi/math/vector.h>

#include <vector>


//----------------------------------------------------------------------

// forward declaration
class Offscreen_context;
class Opengl_application_buffer;

/// OpenGL related application data
///
/// But this is not depends on OpenGL library. This keeps truck the
/// IndeX viewer's OpenGL related function state.
class OpenGL_AppData
{
public:
    /// Constructor
    OpenGL_AppData();
    /// Destructor
    virtual ~OpenGL_AppData();

    /// set show color table state
    /// \param[in] is_show when true show the color table
    void set_show_color_table(bool is_show){
        m_is_show_color_table = is_show;
    }
    /// get show color table state
    /// \return show the color table when true
    bool is_show_color_table() const {
        return m_is_show_color_table;
    }

    /// set color table position
    /// \param[in] table_pos color table position
    void set_color_table_position(mi::math::Vector< mi::Sint32, 2 > const & table_pos){
        m_color_table_position = table_pos;
    }
    /// get color table position
    /// \return color table position
    mi::math::Vector< mi::Sint32, 2 > get_color_table_position() const {
        return m_color_table_position;
    }

    /// set colormap scale
    /// \param[in] map_scale colormap scale
    void set_colormap_scale(mi::Float32 map_scale){
        m_colormap_scale = map_scale;
    }
    /// get colormap scale
    /// \return colormap scale
    mi::Float32 get_colormap_scale() const {
        return m_colormap_scale;
    }

    /// set axis display list
    /// \param[in] axis_list axis display list
    void set_axis_display_list(mi::Sint32 axis_list){
        m_axis_display_list = axis_list;
    }
    /// get axis display list
    /// \return display list for axis
    mi::Sint32 get_axis_display_list() const {
        return m_axis_display_list;
    }

    /// set quad display list
    /// \param[in] quad_list quad display list
    void set_quad_display_list(mi::Sint32 quad_list){
        m_quad_display_list = quad_list;
    }
    /// get quad display list
    /// \return display list for quad
    mi::Sint32 get_quad_display_list() const {
        return m_quad_display_list;
    }

    /// set cube display list
    /// \param[in] cube_list cube display list
    void set_cube_display_list(mi::Sint32 cube_list){
        m_cube_display_list = cube_list;
    }
    /// get cube display list
    /// \return display list for cube
    mi::Sint32 get_cube_display_list() const {
        return m_cube_display_list;
    }

    /// set offscreen context
    void set_offscreen_context(Offscreen_context * p_context);

public:
    //----------------------------------------------------------------------
    // experimental
    //----------------------------------------------------------------------

    /// get opengl application buffer pointer
    /// \return opengl application buffer. may 0.
    Opengl_application_buffer * get_opengl_application_buffer_ptr();

public:
    /// resize offscreen context
    ///
    /// Implementation depends on USE_OPENGL or not
    ///
    /// \param[in] width new width
    /// \param[in] height new height
    void resize_offscreen_context(mi::Sint32 width, mi::Sint32 height);

private:
    /// switch to show color table.
    bool        m_is_show_color_table;
    /// color table position
    mi::math::Vector< mi::Sint32, 2 > m_color_table_position;
    /// colormap scaling factor
    mi::Float32 m_colormap_scale;
    /// display list for axis
    mi::Sint32 m_axis_display_list;
    /// display list for quad
    mi::Sint32 m_quad_display_list;
    /// display list for cube
    mi::Sint32 m_cube_display_list;
    /// offscreen context
    Offscreen_context * m_p_context;

    //----------------------------------------------------------------------
    // experimental
    //----------------------------------------------------------------------
    /// opengl application buffer
    Opengl_application_buffer * m_p_opengl_application_buffer;

private:
    /// copy constructor. prohibit until proved useful.
    OpenGL_AppData(OpenGL_AppData const &);
    /// operator=. prohibit until proved useful.
    OpenGL_AppData const & operator=(OpenGL_AppData const &);
};

//----------------------------------------------------------------------
#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_OPENGL_OPENGL_APPDATA_H
