/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "opengl_drawing_utilities.h"

#include <cassert>
#include <cmath>

#include <mi/math/matrix.h>

#include <nv/index/icamera.h>
#include <nv/index/icolormap.h>
#include <nv/index/iperformance_values.h>
#include <nv/index/itriangle_mesh_scene_element.h>

#include <GL/glew.h>

#include "common/common_utility.h"
#include "common/forwarding_logger.h"
#include "common/string_dict.h"

#include "camera_utility.h"
#include "heightfield_workflow_functionality.h"
#include "nvindex_appdata.h"
#include "opengl_appdata.h"
#include "opengl_offscreen_context.h"
#include "span_renderer_if.h"

//----------------------------------------------------------------------
void gl_initialize_opengl()
{
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
}

//----------------------------------------------------------------------
std::string gl_error()
{
    GLenum code = glGetError();
    std::string error;
    if (code != GL_NO_ERROR)
    {
        error = reinterpret_cast<const char*>(gluErrorString(code));
    }
    return error;
}

//----------------------------------------------------------------------
void gl_clear_buffer()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

//----------------------------------------------------------------------
void gl_initialize_glew()
{
    GLenum err = glewInit();
    if (err) 
    {
        ERROR_LOG << "The OpenGL Extension Wrangler Library couldn't initialize.";
        exit(1);
    }
}

//----------------------------------------------------------------------
void gl_init_display_lists()
{
    OpenGL_AppData* ogl_appdata = Nvindex_AppData::instance()->get_opengl_appdata();
    assert(ogl_appdata != 0);
    
    if (ogl_appdata->get_quad_display_list() == 0) 
    {
        const GLfloat quad_extend = 0.2f;

        GLint quad_list = glGenLists(1);
        glNewList(quad_list, GL_COMPILE);

        glBegin(GL_QUADS);
        glVertex3f(-quad_extend, -quad_extend, 0.0f);
        glVertex3f( quad_extend, -quad_extend, 0.0f);
        glVertex3f( quad_extend,  quad_extend, 0.0f);
        glVertex3f(-quad_extend,  quad_extend, 0.0f);
        glEnd();

        glEndList();

        ogl_appdata->set_quad_display_list(quad_list);
    }

    if (ogl_appdata->get_axis_display_list() == 0)
    {
        const GLfloat axis_length = 1.0f;

        GLint axis_list = glGenLists(1);
        glNewList(axis_list, GL_COMPILE);

        glBegin(GL_LINES);
        glColor3f(0.8f, 0.2f, 0.2f);    // x axis starts - red
        glVertex3f(0.0f, 0.0f, 0.0f);
        glVertex3f(axis_length, 0.0f, 0.0f);
        glEnd();

        glBegin(GL_LINES);
        glColor3f(0.2f, 0.8f, 0.2f);    // y axis starts - green
        glVertex3f(0.0f, 0.0f, 0.0f);
        glVertex3f(0.0f, axis_length, 0.0f);
        glEnd();

        glBegin(GL_LINES);
        glColor3f(0.2f, 0.2f, 0.8f);    // z axis starts - blue
        glVertex3f(0.0f, 0.0f, 0.0f);
        glVertex3f(0.0f, 0.0f, axis_length);
        glEnd();

        glEndList();

        ogl_appdata->set_axis_display_list(axis_list);
    }

    if (ogl_appdata->get_cube_display_list() == 0)
    {
        const GLfloat cube_side = 1.0f;

        GLint cube_list = glGenLists(1);
        glNewList(cube_list, GL_COMPILE);

        glBegin(GL_LINE_LOOP);
        glVertex3f(-cube_side, -cube_side, -cube_side);
        glVertex3f( cube_side, -cube_side, -cube_side);
        glVertex3f( cube_side,  cube_side, -cube_side);
        glVertex3f(-cube_side,  cube_side, -cube_side);
        glEnd();

        glBegin(GL_LINE_LOOP);
        glVertex3f(-cube_side, -cube_side, cube_side);
        glVertex3f( cube_side, -cube_side, cube_side);
        glVertex3f( cube_side,  cube_side, cube_side);
        glVertex3f(-cube_side,  cube_side, cube_side);
        glEnd();

        glBegin(GL_LINES);
        glVertex3f(-cube_side, -cube_side, -cube_side);
        glVertex3f(-cube_side, -cube_side,  cube_side);
        glVertex3f( cube_side, -cube_side, -cube_side);
        glVertex3f( cube_side, -cube_side,  cube_side);
        glVertex3f( cube_side,  cube_side, -cube_side);
        glVertex3f( cube_side,  cube_side,  cube_side);
        glVertex3f(-cube_side,  cube_side, -cube_side);
        glVertex3f(-cube_side,  cube_side,  cube_side);
        glEnd();

        glEndList();

        ogl_appdata->set_cube_display_list(cube_list);
    }
}

//----------------------------------------------------------------------
void gl_debug_draw_viewport(mi::Sint32 viewport_width, mi::Sint32 viewport_height)
{
    // INFO_LOG << "DEBUG: resolution " << viewport_width << " " << viewport_height;
    glLineWidth(3);
    glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
    glBegin(GL_LINE_LOOP);
    {
        glVertex2i(1, 1);
        glVertex2i(viewport_width + 1, 1);
        glVertex2i(viewport_width + 1, viewport_height - 1);
        glVertex2i(1, viewport_height - 1);
    }
    glEnd();
}

//----------------------------------------------------------------------
/// bounding box draw (local function in this file)
template<class T>
void gl_render_bbox(const mi::math::Bbox_struct<T, 3>& bbox)
{
    const mi::math::Vector_struct<T, 3> min = bbox.min;
    const mi::math::Vector_struct<T, 3> max = bbox.max;
    const mi::Float32 i_range = mi::math::abs(mi::Float32(max.x)-mi::Float32(min.x));
    const mi::Float32 j_range = mi::math::abs(mi::Float32(max.y)-mi::Float32(min.y));
    const mi::Float32 k_range = mi::math::abs(mi::Float32(max.z)-mi::Float32(min.z));
    glPushMatrix();
    glTranslatef(mi::Float32(min.x), mi::Float32(min.y), mi::Float32(min.z));
    glScalef(i_range, j_range, k_range);
    glScalef(0.5f, 0.5f, 0.5f);
    glTranslatef(1.f, 1.f, 1.f);
    glCallList(Nvindex_AppData::instance()->get_opengl_appdata()->get_cube_display_list());
    glPopMatrix();
}

//----------------------------------------------------------------------
void gl_debug_render_bbox(const mi::math::Bbox< mi::Float32, 3 >& bbox)
{
    const mi::math::Vector< mi::Float32, 3 > min = bbox.min;
    const mi::math::Vector< mi::Float32, 3 > max = bbox.max;

    glBegin(GL_LINE_LOOP);
    {
        glVertex3f(min.x, min.y, min.z);
        glVertex3f(max.x, min.y, min.z);
        glVertex3f(max.x, min.y, max.z);
        glVertex3f(min.x, min.y, max.z);
    }
    glEnd();

    glBegin(GL_LINE_LOOP);
    {
        glVertex3f(min.x, max.y, min.z);
        glVertex3f(max.x, max.y, min.z);
        glVertex3f(max.x, max.y, max.z);
        glVertex3f(min.x, max.y, max.z);
    }
    glEnd();

    glBegin(GL_LINES);
    {
        glVertex3f(min.x, min.y, min.z);
        glVertex3f(min.x, max.y, min.z);

        glVertex3f(max.x, min.y, min.z);
        glVertex3f(max.x, max.y, min.z);

        glVertex3f(max.x, min.y, max.z);
        glVertex3f(max.x, max.y, max.z);

        glVertex3f(min.x, min.y, max.z);
        glVertex3f(min.x, max.y, max.z);
    }
    glEnd();
}

//----------------------------------------------------------------------
void gl_show_icamera_coordinate_system(
    bool                                       show_coordinate_system,
    const mi::math::Matrix<mi::Float32, 4, 4>& mat)
{
    if (show_coordinate_system)
    {
        mi::math::Matrix<mi::Float32, 4, 4> invmat = mat;
        invmat.invert();

        glLineWidth(2.0f);
        // undo the world coordinate change (assume the modelview is set)
        glMultMatrixf(&(invmat[0][0]));

        // render coord system in world space
        glCallList(Nvindex_AppData::instance()->get_opengl_appdata()->get_axis_display_list());

        // redo the world coordinate change
        glMultMatrixf(&(mat[0][0]));
    }
}

//----------------------------------------------------------------------
void gl_draw_axis(
    const mi::neuraylib::Tag&         camera_tag,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    OpenGL_AppData* ogl_appdata = Nvindex_AppData::instance()->get_opengl_appdata();
    assert(ogl_appdata != 0);

    mi::Sint32 viewport[4];
    glGetIntegerv(GL_VIEWPORT , viewport);

    glDisable(GL_DEPTH_TEST);

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(-0.25, (GLdouble)viewport[2]/100.0,
            -0.25, (GLdouble)viewport[3]/100.0, -1.0, 1.0);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColor4f(0.1f, 0.1f, 0.1f, 0.5f);
    glCallList(ogl_appdata->get_quad_display_list());
    glDisable(GL_BLEND);

    mi::base::Handle< const nv::index::ICamera > cam(
        dice_transaction->access< nv::index::ICamera >(camera_tag));
    assert(cam.is_valid_interface());

    mi::math::Vector<mi::Float32, 3> eye, center, up;
    Camera_tool::get_glu_lookat_vector(cam.get(), eye, center, up);

    // move the camera to the center to remove translation.
    center -= eye;
    eye     = mi::math::Vector<mi::Float32, 3>(0.0f, 0.0f, 0.0f); // => eye -= eye;
    gluLookAt(eye.x,    eye.y,    eye.z,
              center.x, center.y, center.z,
              up.x,     up.y,     up.z);

    glLineWidth(1.1f);
    glPushMatrix();
    glScalef(0.17f, 0.17f, 0.17f);
    glCallList(ogl_appdata->get_axis_display_list());
    glLineWidth(1.0f);
    glPopMatrix();

    glColor3f(0.5f, 0.5f, 0.5f);
    glPushMatrix();
    glScalef(0.12f, 0.12f, 0.12f);
    glCallList(ogl_appdata->get_cube_display_list());
    glPopMatrix();

    glEnable(GL_DEPTH_TEST);

    glPopMatrix();

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();

    glMatrixMode(GL_MODELVIEW);
}

//----------------------------------------------------------------------
void gl_display_colormap(
    const mi::neuraylib::Tag& colormap_tag,
    mi::neuraylib::IDice_transaction * dice_transaction)
{
    if (!colormap_tag.is_valid())
    {
        return;
    }

    assert(dice_transaction != 0);
    mi::base::Handle<const nv::index::IColormap> colormap(
        dice_transaction->access<nv::index::IColormap>(colormap_tag));
    if (!colormap)
    {
        return;
    }

    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    const mi::math::Vector< mi::Sint32, 2 > color_table_pos =
        Nvindex_AppData::instance()->get_opengl_appdata()->get_color_table_position();

    mi::Sint32 color_table_x = color_table_pos.x;
    mi::Sint32 color_table_y = color_table_pos.y - 15;
    mi::Sint32 color_table_height = 10;

    glColor4f(0.5f, 0.5f, 0.5f, 0.8f);
    glBegin(GL_LINE_LOOP);
    glVertex2i(color_table_x-1,         color_table_y-1);
    glVertex2i(color_table_x+1+255,     color_table_y-1);
    glVertex2i(color_table_x+1+255,     color_table_y+1+color_table_height);
    glVertex2i(color_table_x-1,         color_table_y+1+color_table_height);
    glEnd();

    const mi::Uint32 nb_colormap_entires = static_cast<mi::Uint32>(colormap->get_number_of_entries());
    glBegin(GL_LINES);
    std::vector<mi::math::Color_struct> color_table(nb_colormap_entires);
    for (mi::Uint32 i = 0; i < nb_colormap_entires; ++i) 
    {
        color_table[i] = colormap->get_color(i);
        const mi::math::Color& entry = color_table[i];

        glColor4f(entry[0], entry[1], entry[2], 1.0f);
        glVertex2i(color_table_x+i, color_table_y);
        glVertex2i(color_table_x+i, color_table_y+color_table_height);
    }
    glEnd();

    glPointSize(1);
    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
}

//----------------------------------------------------------------------
void gl_render_color_table(
    const mi::neuraylib::Tag&          colormap_tag,
    mi::neuraylib::IDice_transaction * dice_transaction)
{
    if (!colormap_tag.is_valid())
    {
        return;
    }

    assert(dice_transaction != 0);    
    mi::base::Handle<const nv::index::IColormap> colormap(
        dice_transaction->access<nv::index::IColormap>(colormap_tag));
    if (!colormap)
    {
        return;
    }

    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    const mi::math::Vector< mi::Sint32, 2 > color_table_pos =
        Nvindex_AppData::instance()->get_opengl_appdata()->get_color_table_position();

    const mi::Sint32 color_table_x = color_table_pos.x;
    const mi::Sint32 color_table_y = color_table_pos.y - 15;
    const mi::Sint32 color_table_height = 10;
    const mi::Sint32 alpha_table_y = color_table_pos.y + color_table_height + 255 + 10;
    const mi::Sint32 alpha_table_height = 30;

    glColor4f(0.5f, 0.5f, 0.5f, 0.8f);
    glBegin(GL_LINE_LOOP);
    glVertex2i(color_table_x-1,         color_table_y-1);
    glVertex2i(color_table_x+1+255,     color_table_y-1);
    glVertex2i(color_table_x+1+255,     color_table_y+1+color_table_height);
    glVertex2i(color_table_x-1,         color_table_y+1+color_table_height);
    glEnd();

    glBegin(GL_LINE_LOOP);
    glVertex2i(color_table_x-1,         alpha_table_y-1);
    glVertex2i(color_table_x+1+255,     alpha_table_y-1);
    glVertex2i(color_table_x+1+255,     alpha_table_y+1+alpha_table_height);
    glVertex2i(color_table_x-1,         alpha_table_y+1+alpha_table_height);
    glEnd();

    const mi::Uint32 nb_colormap_entires = static_cast<mi::Uint32>(colormap->get_number_of_entries());
    glBegin(GL_LINES);
    std::vector<mi::math::Color_struct> color_table(nb_colormap_entires);
    for (mi::Uint32 i=0; i<nb_colormap_entires; ++i)
    {
        color_table[i] = colormap->get_color(i);
        const mi::math::Color& entry = color_table[i];

        glColor4f(entry[0], entry[1], entry[2], 1.0f);
        glVertex2i(color_table_x+i, color_table_y);
        glVertex2i(color_table_x+i, color_table_y+color_table_height);

        glColor4f(entry[0], entry[1], entry[2], entry[3]);
        glVertex2i(color_table_x+i,  alpha_table_y);
        glVertex2i(color_table_x+i, alpha_table_y+alpha_table_height);

    }
    glEnd();

    /// Edit-area for alpha
    glColor4f(0.5f, 0.5f, 0.5f, 0.8f);
    glBegin(GL_LINE_LOOP);
    glVertex2i(color_table_x-1,     color_table_y-1);
    glVertex2i(color_table_x+255+1, color_table_y-1);
    glVertex2i(color_table_x+255+1, color_table_y+255+1);
    glVertex2i(color_table_x-1,     color_table_y+255+1);
    glEnd();

    glLineWidth(1);

    const mi::Float32 colormap_scale =
        Nvindex_AppData::instance()->get_opengl_appdata()->get_colormap_scale();

    glBegin(GL_LINES);
    for (mi::Uint32 i=0; i<nb_colormap_entires; ++i)
    {
        const mi::math::Color& entry = color_table[i];
        mi::Float32 alpha = 255*entry[3] * colormap_scale;
        alpha = mi::math::max(mi::math::min(alpha, 255.0f), 0.0f);
        glColor4f(entry[0], entry[1], entry[2], 0.9f);
        glVertex2i(static_cast<GLint>(color_table_x+i), static_cast<GLint>(color_table_y));
        glVertex2i(static_cast<GLint>(color_table_x+i), static_cast<GLint>(color_table_y+alpha));
    }
    glEnd();

    glPointSize(3);
    glBegin(GL_POINTS);
    for (mi::Uint32 i=0; i<nb_colormap_entires; ++i)
    {
        const mi::math::Color& entry = color_table[i];
        mi::Float32 alpha = 255*entry[3] * colormap_scale;
        alpha = mi::math::max(mi::math::min(alpha, 255.0f), 0.0f);
        glColor4f(0.3f, 0.3f, 0.3f, 0.5f);
        glVertex2f(static_cast<GLfloat>(color_table_x+i),
                   static_cast<GLfloat>(color_table_y+alpha));
    }
    glEnd();
    glPointSize(2);
    glBegin(GL_POINTS);
    for (mi::Uint32 i=0; i<nb_colormap_entires; ++i)
    {
        const mi::math::Color& entry = color_table[i];
        mi::Float32 alpha = 255*entry[3] * colormap_scale;
        alpha = mi::math::max(mi::math::min(alpha, 255.0f), 0.0f);
        glColor4f(0.5f, 0.5f, 0.5f, 0.75f);
        glVertex2i(static_cast<GLint>(color_table_x+i), static_cast<GLint>(color_table_y+alpha));
    }
    glEnd();
    glPointSize(1);
    glBegin(GL_POINTS);
    for (mi::Uint32 i=0; i<nb_colormap_entires; ++i)
    {
        const mi::math::Color& entry = color_table[i];
        mi::Float32 alpha = 255*entry[3] * colormap_scale;
        alpha = mi::math::max(mi::math::min(alpha, 255.0f), 0.0f);
        glColor4f(0.7f, 0.7f, 0.7f, 1.0f);
        glVertex2i(static_cast<GLint>(color_table_x+i), static_cast<GLint>(color_table_y+alpha));
    }
    glEnd();
    glPointSize(1);

    glDisable(GL_BLEND);

    glEnable(GL_DEPTH_TEST);
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
    const GLfloat WIDTH  = 1000;
    const GLfloat HEIGHT = 1000;
    const GLfloat ALPHA  = 0.5f; //(show_large ? 0.5f : 1.f);

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, WIDTH, 0, HEIGHT, 0, 1);

    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();

    if (!show_large)
    {
        glTranslatef(350.f, -10.f, 0.f);
        glScalef(0.6f, 0.15f, 1.f);
    }

    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    mi::Sint32 cubes_max = 0;
    GLfloat time_max = 0.f;
    for (size_t i=0; i < per_host_rendered_cubes.size(); ++i)
    {
        if (per_host_rendered_cubes[i] > cubes_max)
            cubes_max = per_host_rendered_cubes[i];

        if (per_host_time_rendering[i] > time_max)
            time_max = per_host_time_rendering[i];
        GLfloat accumulated = per_host_time_horizons[i] + per_host_time_volume[i] + per_host_time_volume_and_horizons[i];
        if (accumulated > time_max)
            time_max = accumulated;
    }

    GLfloat left = WIDTH * 0.1f;
    GLfloat skip = (WIDTH - WIDTH * 0.2f) / (per_host_rendered_cubes.size() * 2);

    GLfloat max_skip;
    if (show_large)
        max_skip = 100;
    else
        max_skip = 70;

    if (skip > max_skip)
        skip = max_skip;

    GLfloat border = HEIGHT * 0.1f;

    for (size_t i=0; i < per_host_rendered_cubes.size(); ++i)
    {
        GLfloat bar = skip * 0.5f;

        // Rendered cube count
        for (mi::Sint32 j=0; j < per_host_rendered_cubes[i]; ++j)
        {
            GLfloat bottom = (j / (GLfloat)cubes_max) * (HEIGHT - border * 2.f) + border;
            GLfloat top = ((j+1) / (GLfloat)cubes_max) * (HEIGHT - border * 2.f) + border;

            if (show_large) 
            {
                if (j % 2 == 0)
                    glColor4f(1.0f, 1.0f, 1.0f, 0.3f);
                else
                    glColor4f(0.0f, 0.0f, 0.0f, 0.3f);
            }
            else
            {
                if (j % 2 == 0)
                    glColor4f(0.85f, 0.85f, 0.85f, 1.0f);
                else
                    glColor4f(0.6f, 0.6f, 0.6f, 1.0f);
            }
            glBegin(GL_QUADS);
            glVertex2f(left, bottom);
            glVertex2f(left + bar, bottom);
            glVertex2f(left + bar, top);
            glVertex2f(left, top);
            glVertex2f(left, bottom);
            glEnd();

        }

        GLfloat height = (per_host_rendered_cubes[i] / (GLfloat)cubes_max) * (HEIGHT - border * 2.f) + border;
        glColor4f(0.85f, 0.85f, 0.85f, 1.0f);
        glBegin(GL_LINE_STRIP);
        glVertex2f(left, border);
        glVertex2f(left + bar, border);
        glVertex2f(left + bar, height);
        glVertex2f(left, height);
        glVertex2f(left, border);
        glEnd();

        left += bar;

        GLfloat bottom = border;

        // Volume time - orange
        glColor4f(0.95f, 0.53f, 0.0f, ALPHA);
        height = (per_host_time_volume[i] / time_max) * (HEIGHT - border * 2.f);
        glBegin(GL_QUADS);
        glVertex2f(left, bottom);
        glVertex2f(left + bar, bottom);
        glVertex2f(left + bar, bottom + height);
        glVertex2f(left, bottom + height);
        glVertex2f(left, bottom);
        glEnd();

        // Volume and heightfield time - yellow
        glColor4f(0.95f, 0.95f, 0.0f, ALPHA);
        bottom += height;
        height = (per_host_time_volume_and_horizons[i] / time_max) * (HEIGHT - border * 2.f);
        glBegin(GL_QUADS);
        glVertex2f(left, bottom);
        glVertex2f(left + bar, bottom);
        glVertex2f(left + bar, bottom + height);
        glVertex2f(left, bottom + height);
        glVertex2f(left, bottom);
        glEnd();

        // Heightfield time - greenish
        glColor4f(0.65f, 0.95f, 0.0f, ALPHA);
        bottom += height;
        height = (per_host_time_horizons[i] / time_max) * (HEIGHT - border * 2.f);
        glBegin(GL_QUADS);
        glVertex2f(left, bottom);
        glVertex2f(left + bar, bottom);
        glVertex2f(left + bar, bottom + height);
        glVertex2f(left, bottom + height);
        glVertex2f(left, bottom);
        glEnd();

        // Entire rendering time
        if (per_host_using_cpu[i])
            glColor4f(0.85f, 0.85f, 0.0f, 1.0f); // yellow
        else if (per_host_streaming[i])
            glColor4f(0.85f, 0.0f, 0.0f, 1.0f); // red
        else
            glColor4f(0.85f, 0.85f, 0.85f, 1.0f); // white

        height = (per_host_time_rendering[i] / time_max) * (HEIGHT - border * 2.f) + border;
        glBegin(GL_LINE_STRIP);
        glVertex2f(left, border);
        glVertex2f(left + bar, border);
        glVertex2f(left + bar, height);
        glVertex2f(left, height);
        glVertex2f(left, border);
        glEnd();


        left += skip * 2;
    }

    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);

    glPopMatrix();

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();

    glMatrixMode(GL_MODELVIEW);
}

//----------------------------------------------------------------------
static const mi::math::Color per_host_colors[] = {
    mi::math::Color(0.6f,    0.7f,    1.0f,    0.7f),
    mi::math::Color(1.0f,    0.7f,    1.0f,    0.7f),
    mi::math::Color(0.8f,    0.3f,    0.3f,    1.0f),
    mi::math::Color(0.7f,    0.7f,    0.7f,    1.0f),
    mi::math::Color(0.4f,    0.6f,    0.26f,   1.0f),
    mi::math::Color(0.48f,   0.12f,   0.37f,   1.0f),
    mi::math::Color(0.86f,   0.83f,   0.17f,   1.0f),
    mi::math::Color(0.48f,   0.25f,   0.57f,   1.0f),
    mi::math::Color(0.6565f, 0.233f,  0.2f,    1.0f),
    mi::math::Color(0.76f,   0.62f,   0.07f,   1.0f),
    mi::math::Color(0.75f,   0.5f,    0.5f,    1.0f),
    mi::math::Color(0.85f,   0.85f,   0.75f,   1.0f),
    mi::math::Color(0.5f,    0.95f,   0.95f,   1.0f),
    mi::math::Color(0.3f,    0.35f,   0.7f,    1.0f),
    mi::math::Color(0.5f,    0.5f,    0.35f,   1.0f),
    mi::math::Color(0.8f,    0.85f,   0.15f,   1.0f),
    mi::math::Color(0.25f,   0.25f,   0.35f,   1.0f)
};
static const mi::Uint32 MAX_PER_HOST_COLORS = 17;

//----------------------------------------------------------------------
void draw_performance(
    mi::Float32                             scale,
    nv::index::IPerformance_values*    perf_values
)
{
    static std::vector<mi::Float32> compositing_time_history;
    static std::vector<mi::Float32> rendering_time_history;
    static std::vector<mi::Float32> frame_time_history;

    const mi::Uint32 HISTORY_SIZE   = 100;
    const mi::Uint32 SAMPLE_SPACING = 5;
    const mi::Uint32 X_ORIGIN       = 30;
    const mi::Uint32 Y_ORIGIN       = 490;
    const mi::Uint32 HEIGHT         = 200;

    // Compositing time
    if (compositing_time_history.size()>HISTORY_SIZE)
        compositing_time_history.erase(compositing_time_history.begin());
    const mi::Float32 compositing_time = perf_values->get_time("time_total_compositing")-perf_values->get_time("time_total_final_compositing");
    if (compositing_time<1000)
        compositing_time_history.push_back(compositing_time);

    // Rendering time
    if (rendering_time_history.size()>HISTORY_SIZE)
        rendering_time_history.erase(rendering_time_history.begin());
    const mi::Float32 rendering_time = perf_values->get_time("time_total_rendering");
    if (rendering_time<1000)
        rendering_time_history.push_back(rendering_time);

    // Overall time
    if (frame_time_history.size()>HISTORY_SIZE)
        frame_time_history.erase(frame_time_history.begin());
    const mi::Float32 frame_time = perf_values->get_time("time_complete_frame");
    if (frame_time<1000)
        frame_time_history.push_back(frame_time);

    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // do it ....
    glColor4f(0.2f, 0.2f, 0.6f, 0.3f);
    glBegin(GL_QUADS);
    glVertex2i(X_ORIGIN-1,                               Y_ORIGIN);
    glVertex2i(X_ORIGIN-1,                               Y_ORIGIN+HEIGHT);
    glVertex2i(X_ORIGIN+1 + HISTORY_SIZE*SAMPLE_SPACING, Y_ORIGIN+HEIGHT);
    glVertex2i(X_ORIGIN+1 + HISTORY_SIZE*SAMPLE_SPACING, Y_ORIGIN);
    glVertex2i(X_ORIGIN-21,                            Y_ORIGIN);
    glVertex2i(X_ORIGIN-21,                            Y_ORIGIN+HEIGHT);
    glVertex2i(X_ORIGIN-4,                             Y_ORIGIN+HEIGHT);
    glVertex2i(X_ORIGIN-4,                             Y_ORIGIN);
    glEnd();
    glColor4f(0.4f, 0.4f, 0.4f, 1.0f);
    glBegin(GL_LINE_LOOP);
    glVertex2i(X_ORIGIN-1,                               Y_ORIGIN);
    glVertex2i(X_ORIGIN-1,                               Y_ORIGIN+HEIGHT);
    glVertex2i(X_ORIGIN+1 + HISTORY_SIZE*SAMPLE_SPACING, Y_ORIGIN+HEIGHT);
    glVertex2i(X_ORIGIN+1 + HISTORY_SIZE*SAMPLE_SPACING, Y_ORIGIN);
    glEnd();
    glBegin(GL_LINE_LOOP);
    glVertex2i(X_ORIGIN-21,                            Y_ORIGIN);
    glVertex2i(X_ORIGIN-21,                            Y_ORIGIN+HEIGHT);
    glVertex2i(X_ORIGIN-4,                             Y_ORIGIN+HEIGHT);
    glVertex2i(X_ORIGIN-4,                             Y_ORIGIN);
    glEnd();

    {
        glLineWidth(1);
        mi::Uint32 nb_entries = compositing_time_history.size();
        mi::Uint32 coming_from_right = X_ORIGIN + HISTORY_SIZE*SAMPLE_SPACING;
        glColor4f(.9f, .3f, .3f, 0.85f);
        if (nb_entries>1)
        {
            glBegin(GL_LINE_STRIP);
            for (mi::Uint32 i=nb_entries-1; i>0; --i)
            {
                glVertex2i(coming_from_right,
                           static_cast< GLint >(Y_ORIGIN+scale*compositing_time_history[i]));
                coming_from_right -= SAMPLE_SPACING;
            }
            glVertex2i(coming_from_right,
                       static_cast< GLint >(Y_ORIGIN+scale*compositing_time_history[0]));
            glEnd();
        }
        else if (nb_entries>0)
        {
            glBegin(GL_POINT);
            glVertex2i(coming_from_right,
                       static_cast< GLint >(Y_ORIGIN+scale*compositing_time_history[0]));
            glEnd();
        }

        nb_entries = rendering_time_history.size();
        coming_from_right = X_ORIGIN + HISTORY_SIZE*SAMPLE_SPACING;
        glColor4f(.3f, .3f, .9f, 0.85f);
        if (nb_entries>1)
        {
            glBegin(GL_LINE_STRIP);
            for (mi::Uint32 i=nb_entries-1; i>0; --i)
            {
                glVertex2i(coming_from_right,
                           static_cast< GLint >(Y_ORIGIN+scale*rendering_time_history[i]));
                coming_from_right -= SAMPLE_SPACING;
            }
            glVertex2i(coming_from_right,
                       static_cast< GLint >(Y_ORIGIN+scale*rendering_time_history[0]));
            glEnd();
        }
        else if (nb_entries>0)
        {
            glBegin(GL_POINT);
            glVertex2i(coming_from_right,
                       static_cast< GLint >(Y_ORIGIN+scale*frame_time_history[0]));
            glEnd();
        }

        nb_entries = frame_time_history.size();
        coming_from_right = X_ORIGIN + HISTORY_SIZE*SAMPLE_SPACING;
        glColor4f(.3f, .9f, .3f, 0.85f);
        if (nb_entries>1)
        {
            glBegin(GL_LINE_STRIP);
            for (mi::Uint32 i=nb_entries-1; i>0; --i)
            {
                glVertex2i(coming_from_right,
                           static_cast< GLint >(Y_ORIGIN+scale*frame_time_history[i]));
                coming_from_right -= SAMPLE_SPACING;
            }
            glVertex2i(coming_from_right, static_cast< GLint >(Y_ORIGIN+scale*frame_time_history[0]));
            glEnd();
        }
        else if (nb_entries>0)
        {
            glBegin(GL_POINT);
            glVertex2i(coming_from_right, static_cast< GLint >(Y_ORIGIN+scale*frame_time_history[0]));
            glEnd();
        }
    }

    glDisable(GL_BLEND);

    glLineWidth(1);

    glEnable(GL_DEPTH_TEST);
}

//----------------------------------------------------------------------
void draw_rendering_workload(
    mi::Uint32                              nb_cluster_nodes,
    mi::Float32                             scale,
    nv::index::IPerformance_values*    perf_values
)
{
    static std::map<mi::Uint32, std::vector<mi::Float32> > rendering_workload_history;

    const mi::Uint32 HISTORY_SIZE   = 100;
    const mi::Uint32 SAMPLE_SPACING = 5;
    const mi::Uint32 X_ORIGIN       = 30;
    const mi::Uint32 Y_ORIGIN       = 280;
    const mi::Uint32 HEIGHT         = 200;

    for (mi::Uint32 i=0; i<nb_cluster_nodes; ++i)
    {
        mi::Uint32 cluster_node_id = i+1;   // TODO: might not be a continuous list!!

        assert(perf_values != 0);
        mi::Float32 time = perf_values->get_time("time_rendering", cluster_node_id);

        std::map<mi::Uint32, std::vector<mi::Float32> >::iterator history_per_host = rendering_workload_history.find(cluster_node_id);
        if (history_per_host == rendering_workload_history.end())
        {
            std::vector<mi::Float32> new_history(1, time);
            rendering_workload_history[cluster_node_id] = new_history;
        }
        else
        {
            if (history_per_host->second.size()>HISTORY_SIZE)
                history_per_host->second.erase(history_per_host->second.begin());

            history_per_host->second.push_back(time);
        }
    }


    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // do it ....
    glColor4f(0.3f, 0.6f, 0.2f, 0.3f);
    glBegin(GL_QUADS);
    glVertex2i(X_ORIGIN-1,                               Y_ORIGIN);
    glVertex2i(X_ORIGIN-1,                               Y_ORIGIN+HEIGHT);
    glVertex2i(X_ORIGIN+1 + HISTORY_SIZE*SAMPLE_SPACING, Y_ORIGIN+HEIGHT);
    glVertex2i(X_ORIGIN+1 + HISTORY_SIZE*SAMPLE_SPACING, Y_ORIGIN);
    glVertex2i(X_ORIGIN-21,                            Y_ORIGIN);
    glVertex2i(X_ORIGIN-21,                            Y_ORIGIN+HEIGHT);
    glVertex2i(X_ORIGIN-4,                             Y_ORIGIN+HEIGHT);
    glVertex2i(X_ORIGIN-4,                             Y_ORIGIN);
    glEnd();
    glColor4f(0.4f, 0.4f, 0.4f, 1.0f);
    glBegin(GL_LINE_LOOP);
    glVertex2i(X_ORIGIN-1,                               Y_ORIGIN);
    glVertex2i(X_ORIGIN-1,                               Y_ORIGIN+HEIGHT);
    glVertex2i(X_ORIGIN+1 + HISTORY_SIZE*SAMPLE_SPACING, Y_ORIGIN+HEIGHT);
    glVertex2i(X_ORIGIN+1 + HISTORY_SIZE*SAMPLE_SPACING, Y_ORIGIN);
    glEnd();
    glBegin(GL_LINE_LOOP);
    glVertex2i(X_ORIGIN-21,                            Y_ORIGIN);
    glVertex2i(X_ORIGIN-21,                            Y_ORIGIN+HEIGHT);
    glVertex2i(X_ORIGIN-4,                             Y_ORIGIN+HEIGHT);
    glVertex2i(X_ORIGIN-4,                             Y_ORIGIN);
    glEnd();


    glDisable(GL_BLEND);

    std::map<mi::Uint32, std::vector<mi::Float32> >::iterator history_itr = rendering_workload_history.begin();
    for (; history_itr!=rendering_workload_history.end(); ++history_itr)
    {
        mi::Uint32 coming_from_right = X_ORIGIN + HISTORY_SIZE*SAMPLE_SPACING;
        mi::Uint32 cluster_node_id = history_itr->first;
        const std::vector<mi::Float32>& history = history_itr->second;

        const mi::Uint32 nb_entries = history.size();
        const mi::math::Color& color = per_host_colors[cluster_node_id-1];
        glColor4f(color.r, color.g, color.b, color.a);

        glLineWidth(2);
        glBegin(GL_LINES);
        glVertex2i(X_ORIGIN-20, Y_ORIGIN+3*cluster_node_id);
        glVertex2i(X_ORIGIN-5,  Y_ORIGIN+3*cluster_node_id);
        glEnd();

        glLineWidth(1);
        if (nb_entries>1)
        {
            glBegin(GL_LINE_STRIP);
            for (mi::Uint32 i=nb_entries-1; i>0; --i)
            {
                glVertex2i(coming_from_right, static_cast< GLint >(Y_ORIGIN+scale*history[i]));
                coming_from_right -= SAMPLE_SPACING;
            }
            glVertex2i(coming_from_right, static_cast< GLint >(Y_ORIGIN+scale*history[0]));
            glEnd();
        }
        else if (nb_entries>0)
        {
            glBegin(GL_POINT);
            glVertex2i(coming_from_right, static_cast< GLint >(Y_ORIGIN+scale*history[0]));
            glEnd();
        }
    }

    glLineWidth(1);

    glEnable(GL_DEPTH_TEST);
}

//----------------------------------------------------------------------
void draw_compositing_workload(
    mi::Uint32                              nb_horizontal_spans,
    mi::Float32                             scale,
    nv::index::IPerformance_values*    perf_values
)
{
    static std::map<mi::Uint32, std::vector<mi::Float32> > compositing_workload_history;

    const mi::Uint32 HISTORY_SIZE   = 100;
    const mi::Uint32 SAMPLE_SPACING = 5;
    const mi::Uint32 X_ORIGIN       = 30;
    const mi::Uint32 Y_ORIGIN       = 70;
    const mi::Uint32 HEIGHT         = 200;

    std::map<mi::Uint32, mi::Float32> compositing_workload;

    for (mi::Uint32 i=0; i<nb_horizontal_spans; ++i)
    {
        mi::base::Handle<nv::index::IPer_span_statistics> statistics(perf_values->get_per_span_statistics(i));
        if (!(statistics.is_valid_interface()))
            continue;

        mi::Uint32 cluster_node = statistics->get_cluster_node_id();

        // accumulate compositing workload per node
        compositing_workload[cluster_node] = mi::math::max(compositing_workload[cluster_node], statistics->get_time("time_compositing_stage"));
    }

    // Compositing time on each cluster node
    std::map<mi::Uint32, mi::Float32>::iterator itr = compositing_workload.begin();
    for (; itr!=compositing_workload.end(); ++itr)
    {
        if (itr->second>1000)
            continue;

        mi::Uint32 cluster_node_id = itr->first;
        std::map<mi::Uint32, std::vector<mi::Float32> >::iterator history_per_host = compositing_workload_history.find(cluster_node_id);
        if (history_per_host == compositing_workload_history.end())
        {
            std::vector<mi::Float32> new_history(1, itr->second);
            compositing_workload_history[cluster_node_id] = new_history;
        }
        else
        {
            if (history_per_host->second.size()>HISTORY_SIZE)
                history_per_host->second.erase(history_per_host->second.begin());

            history_per_host->second.push_back(itr->second);
        }
    }

    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // do it ....
    glColor4f(0.6f, 0.2f, 0.2f, 0.3f);
    glBegin(GL_QUADS);
    glVertex2i(X_ORIGIN-1,                               Y_ORIGIN);
    glVertex2i(X_ORIGIN-1,                               Y_ORIGIN+HEIGHT);
    glVertex2i(X_ORIGIN+1 + HISTORY_SIZE*SAMPLE_SPACING, Y_ORIGIN+HEIGHT);
    glVertex2i(X_ORIGIN+1 + HISTORY_SIZE*SAMPLE_SPACING, Y_ORIGIN);
    glVertex2i(X_ORIGIN-21,                            Y_ORIGIN);
    glVertex2i(X_ORIGIN-21,                            Y_ORIGIN+HEIGHT);
    glVertex2i(X_ORIGIN-4,                             Y_ORIGIN+HEIGHT);
    glVertex2i(X_ORIGIN-4,                             Y_ORIGIN);
    glEnd();
    glColor4f(0.4f, 0.4f, 0.4f, 1.0f);
    glBegin(GL_LINE_LOOP);
    glVertex2i(X_ORIGIN-1,                               Y_ORIGIN);
    glVertex2i(X_ORIGIN-1,                               Y_ORIGIN+HEIGHT);
    glVertex2i(X_ORIGIN+1 + HISTORY_SIZE*SAMPLE_SPACING, Y_ORIGIN+HEIGHT);
    glVertex2i(X_ORIGIN+1 + HISTORY_SIZE*SAMPLE_SPACING, Y_ORIGIN);
    glEnd();
    glBegin(GL_LINE_LOOP);
    glVertex2i(X_ORIGIN-21,                            Y_ORIGIN);
    glVertex2i(X_ORIGIN-21,                            Y_ORIGIN+HEIGHT);
    glVertex2i(X_ORIGIN-4,                             Y_ORIGIN+HEIGHT);
    glVertex2i(X_ORIGIN-4,                             Y_ORIGIN);
    glEnd();

    glDisable(GL_BLEND);

    std::map<mi::Uint32, std::vector<mi::Float32> >::iterator history_itr = compositing_workload_history.begin();
    for (; history_itr!=compositing_workload_history.end(); ++history_itr)
    {
        mi::Uint32 coming_from_right = X_ORIGIN + HISTORY_SIZE*SAMPLE_SPACING;
        mi::Uint32 cluster_node_id = history_itr->first;
        const std::vector<mi::Float32>& history = history_itr->second;

        const mi::Uint32 nb_entries = history.size();
        const mi::math::Color& color = per_host_colors[cluster_node_id-1];
        glColor4f(color.r, color.g, color.b, color.a);

        glLineWidth(2);
        glBegin(GL_LINES);
        glVertex2i(X_ORIGIN-20, Y_ORIGIN+3*cluster_node_id);
        glVertex2i(X_ORIGIN-5,  Y_ORIGIN+3*cluster_node_id);
        glEnd();

        glLineWidth(1);
        if (nb_entries>1)
        {
            glBegin(GL_LINE_STRIP);
            for (mi::Uint32 i=nb_entries-1; i>0; --i)
            {
                glVertex2i(coming_from_right, static_cast< GLint >(Y_ORIGIN+scale*history[i]));
                coming_from_right -= SAMPLE_SPACING;
            }
            glVertex2i(coming_from_right, static_cast< GLint >(Y_ORIGIN+scale*history[0]));
            glEnd();
        }
        else if (nb_entries>0)
        {
            glBegin(GL_POINT);
            glVertex2i(coming_from_right, static_cast< GLint >(Y_ORIGIN+scale*history[0]));
            glEnd();
        }
    }

    glLineWidth(1);

    glEnable(GL_DEPTH_TEST);
}


//----------------------------------------------------------------------
void gl_initialize_viewport(
    mi::Sint32 width, mi::Sint32 height,
    const mi::math::Color_struct& background_color)
{
    glViewport(0, 0, width, height);

    glClearColor(background_color.r, background_color.g, background_color.b, background_color.a);
    gl_clear_buffer();

    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

//----------------------------------------------------------------------
void gl_get_perspective_matrix_parameters(mi::Float32 & mat_fovy_rad,
                                          mi::Float32 & mat_aspect,
                                          mi::Float32 & mat_clip_min,
                                          mi::Float32 & mat_clip_max)
{
    GLfloat mat[16];
    glGetFloatv(GL_PROJECTION_MATRIX, mat);

    const GLfloat aa = mat[0];
    const GLfloat bb = mat[5];
    const GLfloat cc = mat[10];
    const GLfloat dd = mat[14];

    mat_aspect   = bb / aa;
    mat_fovy_rad = 2.0f * atan(1.0f / bb);

    const GLfloat kk = (cc - 1.0f) / (cc + 1.0f);
    mat_clip_min = (dd * (1.0f - kk)) / (2.0f * kk);
    mat_clip_max = kk * mat_clip_min;
}

//----------------------------------------------------------------------
void gl_get_modelview_matrix_parameters(mi::math::Vector<mi::Float32, 3> & mat_eyepos,
                                        mi::math::Vector<mi::Float32, 3> & mat_viewdir,
                                        mi::math::Vector<mi::Float32, 3> & mat_updir)
{
    GLfloat mat[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, mat);

    mi::math::Vector<mi::Float32, 3> xdir(0.0f, 0.0f, 0.0f);
    mi::math::Vector<mi::Float32, 3> ydir(0.0f, 0.0f, 0.0f);
    mi::math::Vector<mi::Float32, 3> zdir(0.0f, 0.0f, 0.0f);

    xdir[0] = mat[0];  ydir[0] = mat[1];  zdir[0] = mat[2];
    xdir[1] = mat[4];  ydir[1] = mat[5];  zdir[1] = mat[6];
    xdir[2] = mat[8];  ydir[2] = mat[9];  zdir[2] = mat[10];

    mi::math::Vector<mi::Float32, 3> bvec(-mat[12], -mat[13], -mat[14]);

    mi::math::Matrix<mi::Float32,3,3>
        basis_mat(xdir[0], xdir[1], xdir[2],
                  ydir[0], ydir[1], ydir[2],
                  zdir[0], zdir[1], zdir[2]);
    mi::math::Matrix<mi::Float32,3,3> basis_mat_inv(basis_mat);
    const bool is_inv = basis_mat_inv.invert();
    assert(is_inv); nv::index_common::no_unused_variable_warning_please(is_inv);

    mat_eyepos  = basis_mat_inv * bvec;
    mat_viewdir = -zdir;
    mat_updir   =  ydir;
}

//----------------------------------------------------------------------
void gl_setup_projection_matrix(const nv::index::ICamera* cam)
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    mi::base::Handle<const nv::index::IPerspective_camera>  perspective_camera(cam->get_interface<const nv::index::IPerspective_camera>());
    mi::base::Handle<const nv::index::IOrthographic_camera> ortho_camera(cam->get_interface<const nv::index::IOrthographic_camera>());

    if (ortho_camera)
    {
        const mi::Float64 film_width  = ortho_camera->get_aperture();
        const mi::Float64 film_height = film_width / ortho_camera->get_aspect();

        const mi::Float64 left    = -film_width / 2.0;
        const mi::Float64 right   =  film_width / 2.0;
        const mi::Float64 bottom  = -film_height / 2.0;
        const mi::Float64 top     =  film_height / 2.0;

        glOrtho(left, right, bottom, top, ortho_camera->get_clip_min(), ortho_camera->get_clip_max());
    }
    else if (perspective_camera)
    {
        const mi::Float64 fovy = perspective_camera->get_fov_y_rad() * (180.0 / MI_PI);
        gluPerspective(fovy, perspective_camera->get_aspect(), perspective_camera->get_clip_min(), perspective_camera->get_clip_max());
    }

    glMatrixMode(GL_MODELVIEW);
}

//----------------------------------------------------------------------
void gl_setup_lookat(const mi::math::Vector<mi::Float32, 3>& eye,
                     const mi::math::Vector<mi::Float32, 3>& target,
                     const mi::math::Vector<mi::Float32, 3>& up)
{
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(eye.x,    eye.y,    eye.z,
              target.x, target.y, target.z,
              up.x,     up.y,     up.z);
}

//----------------------------------------------------------------------
void gl_set_world_transform(
    const mi::math::Matrix<mi::Float32, 4, 4>& transform_mat)
{
    glMultMatrixf(&(transform_mat[0][0]));
}

//----------------------------------------------------------------------
void gl_visualize_roi(
    bool                                         is_show_roi,
    const mi::math::Bbox_struct<mi::Float32, 3>& roi_bbox,
    const mi::math::Color&                       color,
    const mi::math::Matrix<mi::Float32, 4, 4>&   transform_mat)
{
    if (!is_show_roi)
    {
        return;
    }

    glPushMatrix();
    glMultMatrixf(&transform_mat.xx);

    /// render bounding box of the entire scene
    mi::math::Bbox_struct<mi::Float32, 3> extent_bbox = roi_bbox;

    glLineWidth(1.0f);
    glColor4f(color.r, color.g, color.b, color.a);
    gl_render_bbox(extent_bbox);

    glPopMatrix();
}

//----------------------------------------------------------------------
void gl_visualize_bbox(
    bool                                         is_show_bbox,
    const mi::math::Bbox_struct<mi::Float32, 3>& bbox,
    const mi::math::Color&                       color,
    const mi::math::Matrix<mi::Float32, 4, 4>&   transform_mat)
{
    if (!is_show_bbox)
    {
        return;
    }

    glPushMatrix();
    glMultMatrixf(&transform_mat.xx);

    glLineWidth(1.0f);
    glColor4f(color.r, color.g, color.b, color.a);
    gl_render_bbox(bbox);

    glPopMatrix();
}

//----------------------------------------------------------------------
void gl_draw_slice(const mi::math::Vector<mi::Float32, 3>&    llf,
                   const mi::math::Vector<mi::Float32, 3>&    urb,
                   const mi::math::Matrix<mi::Float32, 4, 4>& transform)
{
    glPushMatrix();
    glMultMatrixf(&transform.xx);

    glBegin(GL_LINE_LOOP);
    if (llf.x == urb.x)
    {
        glVertex3f(llf.x, llf.y, llf.z);
        glVertex3f(llf.x, llf.y, urb.z);
        glVertex3f(urb.x, urb.y, urb.z);
        glVertex3f(llf.x, urb.y, llf.z);
    }
    else if (llf.y == urb.y)
    {
        glVertex3f(llf.x, llf.y, llf.z);
        glVertex3f(llf.x, llf.y, urb.z);
        glVertex3f(urb.x, urb.y, urb.z);
        glVertex3f(urb.x, urb.y, llf.z);
    }
    else if (llf.z == urb.z)
    {
        glVertex3f(llf.x, llf.y, llf.z);
        glVertex3f(llf.x, urb.y, urb.z);
        glVertex3f(urb.x, urb.y, urb.z);
        glVertex3f(urb.x, llf.y, llf.z);
    }
    glEnd();

    glPopMatrix();
}

//----------------------------------------------------------------------
void gl_draw_profile(const mi::math::Vector_struct<mi::Float32, 2>& v0,
                     const mi::math::Vector_struct<mi::Float32, 2>& v1,
                     const mi::math::Bbox_struct<mi::Float32, 3>&   ijk_space,
                     const mi::math::Matrix<mi::Float32, 4, 4>&     transform)
{
    glPushMatrix();
    glMultMatrixf(&transform.xx);

    glBegin(GL_LINE_LOOP);
    glVertex3f(v0.x, v0.y, ijk_space.min.z);
    glVertex3f(v0.x, v0.y, ijk_space.max.z);
    glVertex3f(v1.x, v1.y, ijk_space.max.z);
    glVertex3f(v1.x, v1.y, ijk_space.min.z);
    glEnd();

    glPopMatrix();
}

//----------------------------------------------------------------------
void gl_setup_2d_rendering(
    mi::Float32 width,
    mi::Float32 height)
{
    glLineWidth(1.0f);

    // reset all matrices, no more OpenGL 3D rendering
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    // Using window coords from here on ...
    glOrtho(0.f, (GLfloat)width,
            0.f, (GLfloat)height, 0.f, 1.f);

    glMatrixMode(GL_MODELVIEW);
}

//----------------------------------------------------------------------
void gl_visualize_horizontal_span(bool                             is_visualize,
                                  Span_renderer_IF*                p_span_buffer_renderer,
                                  nv::index::IPerformance_values*  perf_values)
{
    if (!is_visualize)
    {
        return;
    }

    // Visualize convex hull
    std::vector<mi::math::Bbox_struct<mi::Uint32, 2> > screen_space_subdivision;
    p_span_buffer_renderer->get_screen_space_subdivision(screen_space_subdivision);

    const mi::Uint32 nb_horizontal_spans = screen_space_subdivision.size();

    mi::Float32 max_compositing_time = 0;
    mi::Float32 min_compositing_time = 1000000;
    for (mi::Uint32 i=0; i<nb_horizontal_spans; ++i)
    {
        mi::base::Handle<nv::index::IPer_span_statistics> statistics(perf_values->get_per_span_statistics(i));
        if (!(statistics.is_valid_interface()))
        {
            continue;
        }
        const mi::Float32 time_image_comp = statistics->get_time("time_image_compositing");

        if (time_image_comp == 0.0f)
        {
            continue;
        }

        max_compositing_time = mi::math::max(max_compositing_time, time_image_comp);
        min_compositing_time = mi::math::min(min_compositing_time, time_image_comp);
    }


    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    for (mi::Uint32 i=0; i<nb_horizontal_spans; ++i)
    {
        mi::base::Handle<nv::index::IPer_span_statistics> statistics(perf_values->get_per_span_statistics(i));
        if (  !(statistics.is_valid_interface())
             || max_compositing_time <= min_compositing_time
             || statistics->get("time_image_compositing") == 0.f)
        {
            const mi::Float32 greyish = mi::Float32(i+1)/mi::Float32(nb_horizontal_spans);
            glColor4f(mi::math::sqrt(greyish), greyish, greyish, 0.2f);
        }
        else
        {
            const mi::Float32 compositing_time = statistics->get_time("time_image_compositing");
            const mi::Float32 blend_value = (compositing_time - min_compositing_time) / (max_compositing_time-min_compositing_time);
            glColor4f(0.2f+(0.4f*blend_value), 0.6f-(0.4f*blend_value), 0.1f, 0.5f);
        }

        const mi::math::Bbox<mi::Uint32,  2> & ss_sub_u = screen_space_subdivision[i];
        const mi::math::Bbox<mi::Float32, 2>   ss_sub_f(static_cast<GLfloat>(ss_sub_u.min.x),
                                                        static_cast<GLfloat>(ss_sub_u.min.y), 
                                                        static_cast<GLfloat>(ss_sub_u.max.x), 
                                                        static_cast<GLfloat>(ss_sub_u.max.y));

        glBegin(GL_QUADS);
        glVertex2f(ss_sub_f.min.x, ss_sub_f.min.y);
        glVertex2f(ss_sub_f.max.x, ss_sub_f.min.y);
        glVertex2f(ss_sub_f.max.x, ss_sub_f.max.y);
        glVertex2f(ss_sub_f.min.x, ss_sub_f.max.y);
        glEnd();

        glColor4f(0.6f, 0.6f, 0.1f, 0.7f);

        glBegin(GL_LINE_LOOP);
        glVertex2f(ss_sub_f.min.x, ss_sub_f.min.y);
        glVertex2f(ss_sub_f.max.x, ss_sub_f.min.y);
        glVertex2f(ss_sub_f.max.x, ss_sub_f.max.y);
        glVertex2f(ss_sub_f.min.x, ss_sub_f.max.y);
        glEnd();
    }
    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
}

//----------------------------------------------------------------------
void gl_render_color_table(
    bool is_show_color_table,
    const mi::neuraylib::Tag& colormap_tag,
    mi::neuraylib::IDice_transaction * dice_transaction)
{
    static bool is_shown_warning = false;
    if (!colormap_tag.is_valid())
    {
        if (!is_shown_warning)
        {
            WARN_LOG << "Invalid colormap_tag " << colormap_tag << ". Cannot show colormap. This message is only shown once.";
            is_shown_warning = true;
        }
        return;
    }

    if (is_show_color_table)
    {
        gl_render_color_table(colormap_tag, dice_transaction);
    }
    else
    {
        gl_display_colormap(colormap_tag, dice_transaction);
    }
}

//----------------------------------------------------------------------
void gl_push_matrix_state()
{
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
}

//----------------------------------------------------------------------
void gl_pop_matrix_state()
{
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
}

//----------------------------------------------------------------------
void gl_delete_offscreen_context(Offscreen_context* p_context)
{
    if (p_context != 0)
    {
        delete p_context;
    }
}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
void gl_draw_heightfield_delete_polygon(
    const mi::math::Matrix<mi::Float32, 4, 4>& ijk_to_xyz_space)
{
    // check the workflow is in the polygon delete state.
    assert(Heightfield_workflow_functionality::instance()->is_heightfield_workflow_delete_operation()
           || Heightfield_workflow_functionality::instance()->is_heightfield_workflow_elevation_change_operation()
           || Heightfield_workflow_functionality::instance()->is_heightfield_workflow_gridding_operation() );

    const std::vector<mi::math::Vector<mi::Float32, 3> >& polygon =
        Heightfield_workflow_functionality::instance()->get_bounding_polygon();
    const mi::Uint32 nb_vertices = polygon.size();

    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);

    glColor3f(0.8f, 0.2f, 0.2f);
    glLineWidth(2);
    glBegin(GL_LINE_LOOP);
    for (mi::Uint32 i=0; i<nb_vertices; ++i)
    {
        const mi::math::Vector<mi::Float32, 3>& point =
            mi::math::transform_point<mi::Float32, mi::Float32>(ijk_to_xyz_space, polygon[i]);
        glVertex3f(point.x, point.y, point.z);
    }
    glEnd();

    if (nb_vertices>2)
    {
        std::vector<mi::math::Vector<mi::Float32, 3> > convex_hull;
        Heightfield_workflow_functionality::instance()->get_bounding_polygon_convex_hull(convex_hull);
        const mi::Uint32 nb_convex_hull_vertices = convex_hull.size();

        glBegin(GL_LINE_LOOP);
        glColor3f(0.2f, 0.8f, 0.2f);
        for (mi::Uint32 i=0; i<nb_convex_hull_vertices; ++i)
        {
            const mi::math::Vector<mi::Float32, 3>& point =
                mi::math::transform_point<mi::Float32, mi::Float32>(ijk_to_xyz_space, convex_hull[i]);
            glVertex3f(point.x, point.y, point.z);
        }
        glEnd();
    }

    glDisable(GL_BLEND);
    glEnable(GL_DEPTH_TEST);
}

//----------------------------------------------------------------------
