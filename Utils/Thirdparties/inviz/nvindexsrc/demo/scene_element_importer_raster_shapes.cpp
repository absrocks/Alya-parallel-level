/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/

#include "scene_element_importer.h"

#include <nv/index/icircle.h>
#include <nv/index/iellipse.h>
#include <nv/index/iicon.h>
#include <nv/index/ilabel.h>
#include <nv/index/iline_set.h>
#include <nv/index/ipoint_set.h>
#include <nv/index/ipolygon.h>
#include <nv/index/ipath.h>

using namespace nv::index_common;

mi::base::Handle<nv::index::IScene_element> Scene_element_importer::import_raster_shapes(
    const std::string& elem_name,
    const std::string& elem_type,
    const String_dict& dict)
{
    //
    // Point set
    //
    if (elem_type == "point_set")
    {
        const std::string nb_vertices_str = dict.get("nb_vertices", "0");
        const mi::Uint32 nb_vertices = get_uint32(nb_vertices_str);

        // Default to number of vertices if these do not exist
        const mi::Uint32 nb_colors   = get_uint32(dict.get("nb_colors", nb_vertices_str));
        const mi::Uint32 nb_radii    = get_uint32(dict.get("nb_radii", nb_vertices_str));

        // Vertices and additional per-vertex colors and radii.
        std::vector<mi::math::Vector_struct<mi::Float32, 3> > vertices;
        std::vector<mi::math::Color_struct> per_vertex_color_values;
        std::vector<mi::Float32> per_vertex_radii;

        // Read the vertices ...
        for (mi::Uint32 i=0; i < nb_vertices; ++i)
        {
            std::ostringstream v;
            v << "v" << i;
            vertices.push_back(get_vec_float32_3(dict.get(v.str(), "0 0 0")));
        }

        // Colors
        const std::string color_default = "1 1 1 1";
        for (mi::Uint32 i=0; i < nb_colors; ++i)
        {
            std::ostringstream c;
            c << "c" << i;
            std::string color = dict.get(c.str());
            if (color.empty())
                color = color_default;
            per_vertex_color_values.push_back(get_color(color));
        }

        // Radii
        const std::string radius_default = dict.get("radius", "5");
        for (mi::Uint32 i=0; i < nb_radii; ++i)
        {
            std::ostringstream r;
            r << "r" << i;
            std::string radius = dict.get(r.str());
            if (radius.empty())
                radius = radius_default;
            per_vertex_radii.push_back(get_float32(radius));
        }

        // Create shape (point set) scene element using the scene's factory
        mi::base::Handle<nv::index::IPoint_set> shape(m_scene->create_shape<nv::index::IPoint_set>());
        // ... and the vertices with colors and radii:
        shape->set_vertices(&vertices[0], vertices.size());
        shape->set_colors(&per_vertex_color_values[0], per_vertex_color_values.size());
        shape->set_radii(&per_vertex_radii[0], per_vertex_radii.size());

        std::string style_str = dict.get("style", "flat_circle");
        nv::index::IPoint_set::Point_style style = nv::index::IPoint_set::SHADED_CIRCLE;
        if (style_str == "flat_circle")
            style = nv::index::IPoint_set::FLAT_CIRCLE;
        else if (style_str == "shaded_circle")
            style = nv::index::IPoint_set::SHADED_CIRCLE;
        else if (style_str == "shaded_flat_circle")
            style = nv::index::IPoint_set::SHADED_FLAT_CIRCLE;
        else if (style_str == "flat_square")
            style = nv::index::IPoint_set::FLAT_SQUARE;
        else if (style_str == "flat_triangle")
            style = nv::index::IPoint_set::FLAT_TRIANGLE;
        else
            ERROR_LOG << "Invalid point style for element '" << elem_name << "': " << style_str;

        // ... set the point styles (see IPoint_set):
        shape->set_point_style(style);

        return shape;
    }
    //
    // Line set
    //
    else if (elem_type == "line_set")
    {
        const mi::Uint32 demo_mode = get_uint32(dict.get("demo", "0"));

        // Create line set using the scene's factory
        mi::base::Handle<nv::index::ILine_set> shape(m_scene->create_shape<nv::index::ILine_set>());
        // ... set the line type. Currently only line segments are supported (see ILine_set)
        std::string line_type_str = dict.get("line_type", "segments");
        if (line_type_str == "segments")
            shape->set_line_type(nv::index::ILine_set::LINE_TYPE_SEGMENTS);
        else if(line_type_str == "path")
            shape->set_line_type(nv::index::ILine_set::LINE_TYPE_PATH);
        else if(line_type_str == "loop")
            shape->set_line_type(nv::index::ILine_set::LINE_TYPE_LOOP);

        // caps style
        std::string cap_style_str = dict.get("caps_style", "flat");
        if (cap_style_str == "flat")
            shape->set_cap_style(nv::index::ILine_set::CAP_STYLE_FLAT);
        else if (cap_style_str == "square")
            shape->set_cap_style(nv::index::ILine_set::CAP_STYLE_SQUARE);
        else if (cap_style_str == "round")
            shape->set_cap_style(nv::index::ILine_set::CAP_STYLE_ROUND);

        const std::string nb_vertices_str = dict.get("nb_vertices", "0");
        const mi::Uint32 nb_vertices = get_uint32(nb_vertices_str);

        // Default to number of vertices if these do not exist
        const mi::Uint32 nb_colors = get_uint32(dict.get("nb_colors", nb_vertices_str));
        const mi::Uint32 nb_widths = get_uint32(dict.get("nb_widths", nb_vertices_str));

        std::vector<mi::math::Vector_struct<mi::Float32, 3> > vertices;
        for (mi::Uint32 i=0; i < nb_vertices; ++i)
        {
            // vertex position
            mi::math::Vector_struct<mi::Float32, 3> vertex;
            if(demo_mode == 0)
            {
                std::ostringstream s;
                s << "v" << i;
                vertex = get_vec_float32_3(dict.get(s.str(), "0 0 0"));
            }
            else if(demo_mode == 1) // ring
            {
                mi::Float32 t = (mi::Float32)i/(mi::Float32)(nb_vertices);
                vertex.x = sinf(2.f*t*static_cast<mi::Float32>(M_PI))*100.f;
                vertex.y = sinf(2.f*t*static_cast<mi::Float32>(M_PI))*50.f;
                vertex.z = cosf(2.f*t*static_cast<mi::Float32>(M_PI))*100.f;
            }
            else if(demo_mode == 2) // lissajous
            {
                mi::Float32 t = (mi::Float32)i/(mi::Float32)(nb_vertices);
                vertex.x = sinf(2.f*float(M_PI)*t + float(M_PI)/2.f)*100.f;
                vertex.y = sinf(4.f*float(M_PI)*t)*100.f;
                vertex.z = sinf(2.f*float(M_PI)*t)*100.f;
            }

            vertices.push_back(vertex);
        }
        if (!vertices.empty())
            shape->set_lines(&vertices[0], vertices.size());

        // segment/vertex color
        const std::string color_default = dict.get("color", "0 0 1 1");
        std::vector<mi::math::Color_struct> per_vertex_color_values;
        for (mi::Uint32 i=0; i < nb_colors; ++i)
        {
            std::ostringstream c;
            c << "c" << i;
            std::string color = dict.get(c.str());
            if (color.empty())
                color = color_default;
            per_vertex_color_values.push_back(get_color(color));
        }
        if (!per_vertex_color_values.empty())
            shape->set_colors(&per_vertex_color_values[0], per_vertex_color_values.size());

        // segment/vertex width
        const std::string width_default = dict.get("width", "5");
        std::vector<mi::Float32> per_vertex_widths;
        for (mi::Uint32 i=0; i < nb_widths; ++i)
        {
            std::ostringstream w;
            w << "w" << i;
            std::string width = dict.get(w.str());
            if (width.empty())
                width = width_default;
            per_vertex_widths.push_back(get_float32(width));
        }
        if (!per_vertex_widths.empty())
            shape->set_widths(&per_vertex_widths[0], per_vertex_widths.size());

        std::string style_str = dict.get("style", "solid");
        nv::index::ILine_set::Line_style style = nv::index::ILine_set::LINE_STYLE_SOLID;
        if (style_str == "solid")
            style = nv::index::ILine_set::LINE_STYLE_SOLID;
        else if (style_str == "dashed")
            style = nv::index::ILine_set::LINE_STYLE_DASHED;
        else if (style_str == "dotted")
            style = nv::index::ILine_set::LINE_STYLE_DOTTED;
        else if (style_str == "center")
            style = nv::index::ILine_set::LINE_STYLE_CENTER;
        else if (style_str == "hidden")
            style = nv::index::ILine_set::LINE_STYLE_HIDDEN;
        else if (style_str == "phantom")
            style = nv::index::ILine_set::LINE_STYLE_PHANTOM;
        else if (style_str == "dashdot")
            style = nv::index::ILine_set::LINE_STYLE_DASHDOT;
        else if (style_str == "border")
            style = nv::index::ILine_set::LINE_STYLE_BORDER;
        else if (style_str == "divide")
            style = nv::index::ILine_set::LINE_STYLE_DIVIDE;
        else
            ERROR_LOG << "Invalid line style for element '" << elem_name << "': " << style_str;
        // ... set the line style such a dashed or dotted or solid linestyle (see ILine_set)
        shape->set_line_style(style);

        return shape;
    }
    //
    // Circle
    //
    else if (elem_type == "circle")
    {
        mi::base::Handle<nv::index::ICircle> circle(
            m_scene->create_shape<nv::index::ICircle>());

        mi::math::Vector<mi::Float32, 3> center = get_vec_float32_3(dict.get("center", "0 0 0"));
        mi::Float32 radius = get_float32(dict.get("radius", "10"));

        mi::math::Color_struct line_color = get_color(dict.get("line_color", "1 0 0 1"));
        mi::Float32 line_width = get_float32(dict.get("line_width", "1.5"));
        mi::math::Color_struct fill_color = get_color(dict.get("fill_color", "0 0 1 1"));

        circle->set_geometry(center, radius);
        circle->set_outline_style(line_color, line_width);
        if (dict.get("fill_style") == "empty")
            circle->set_fill_style(fill_color, nv::index::ICircle::FILL_EMPTY);
        else
            circle->set_fill_style(fill_color, nv::index::ICircle::FILL_SOLID);

        return circle;
    }
    //
    // Ellipses
    //
    else if (elem_type == "ellipse")
    {
        mi::base::Handle<nv::index::IEllipse> ellipse(
            m_scene->create_shape<nv::index::IEllipse>());

        mi::math::Vector<mi::Float32, 3> center = get_vec_float32_3(dict.get("center", "0 0 0"));
        mi::Float32 radius_x = get_float32(dict.get("radius_x", "10"));
        mi::Float32 radius_y = get_float32(dict.get("radius_y", "10"));
        mi::Float32 rotation = get_float32(dict.get("rotation", "0"));

        mi::math::Color_struct line_color = get_color(dict.get("line_color", "1 0 0 1"));
        mi::Float32 line_width = get_float32(dict.get("line_width", "1.5"));
        mi::math::Color_struct fill_color = get_color(dict.get("fill_color", "0 0 1 1"));

        ellipse->set_geometry(center, radius_x, radius_y, rotation);
        ellipse->set_outline_style(line_color, line_width);
        if (dict.get("fill_style") == "empty")
            ellipse->set_fill_style(fill_color, nv::index::IEllipse::FILL_EMPTY);
        else
            ellipse->set_fill_style(fill_color, nv::index::IEllipse::FILL_SOLID);

        return ellipse;
    }
    //
    // Polygon
    //
    else if (elem_type == "polygon")
    {
        mi::base::Handle<nv::index::IPolygon> polygon(
            m_scene->create_shape<nv::index::IPolygon>());

        mi::math::Vector<mi::Float32, 3> center = get_vec_float32_3(dict.get("center", "0 0 0"));
        const mi::Uint32 nb_vertices = get_uint32(dict.get("nb_vertices", "0"));

        // Vertices
        std::vector<mi::math::Vector_struct<mi::Float32, 2> > vertices;

        // Read the vertices ...
        for (mi::Uint32 i=0; i < nb_vertices; ++i)
        {
            std::ostringstream v;
            v << "v" << i;
            vertices.push_back(get_vec_float32_2(dict.get(v.str(), "0 0")));
        }

        if (polygon->set_geometry(vertices.data(), nb_vertices, center))
        {
            // line style
            //mi::math::Color_struct line_color = get_color(dict.get("line_color", "1 0 0 1"));
            //mi::Float32 line_width = get_float32(dict.get("line_width", "1.5"));
            // polygon->set_outline_style(line_color, line_width);

            // fill style
            mi::math::Color_struct fill_color = get_color(dict.get("fill_color", "0 0 1 1"));
            if (dict.get("fill_style") == "empty")
                polygon->set_fill_style(fill_color, nv::index::IPolygon::FILL_EMPTY);
            else
                polygon->set_fill_style(fill_color, nv::index::IPolygon::FILL_SOLID);

            return polygon;
        }
    }

    //
    // Label
    //
    else if (elem_type == "label_3d")
    {
        mi::base::Handle<nv::index::ILabel_3D> label(
            m_scene->create_shape<nv::index::ILabel_3D>());

        label->set_text(dict.get("text", "<no text specified>").c_str());

        mi::math::Vector<mi::Float32, 3> pos = get_vec_float32_3(dict.get("position", "0 0 0"));
        mi::math::Vector<mi::Float32, 3> right = get_vec_float32_3(dict.get("right", "1 0 0"));
        mi::math::Vector<mi::Float32, 3> up = get_vec_float32_3(dict.get("up", "0 1 0"));
        mi::Float32 height = get_float32(dict.get("height", "10"));
        mi::Float32 width = get_float32(dict.get("width", "0"));
        label->set_geometry(pos, right, up, height, width);

        label->set_pickable(get_bool(dict.get("pickable", "true")));

        return label;
    }
    else if (elem_type == "label_2d")
    {
        mi::base::Handle<nv::index::ILabel_2D> label(
            m_scene->create_shape<nv::index::ILabel_2D>());

        label->set_text(dict.get("text", "<no text specified>").c_str());

        mi::math::Vector<mi::Float32, 3> pos = get_vec_float32_3(dict.get("position", "0 0 0"));
        mi::math::Vector<mi::Float32, 2> right = get_vec_float32_2(dict.get("right", "1 0"));
        mi::math::Vector<mi::Float32, 2> up = get_vec_float32_2(dict.get("up", "0 1"));
        mi::Float32 height = get_float32(dict.get("height", "10"));
        mi::Float32 width = get_float32(dict.get("width", "0"));
        label->set_geometry(pos, right, up, height, width);

        label->set_pickable(get_bool(dict.get("pickable", "true")));

        return label;
    }
    //
    // Icon
    //
    else if (elem_type == "icon_3d")
    {
        mi::base::Handle<nv::index::IIcon_3D> icon(
            m_scene->create_shape<nv::index::IIcon_3D>());

        mi::math::Vector<mi::Float32, 3> pos = get_vec_float32_3(dict.get("position", "0 0 0"));
        mi::math::Vector<mi::Float32, 3> right = get_vec_float32_3(dict.get("right", "1 0 0"));
        mi::math::Vector<mi::Float32, 3> up = get_vec_float32_3(dict.get("up", "0 1 0"));
        mi::Float32 height = get_float32(dict.get("height", "10"));
        mi::Float32 width = get_float32(dict.get("width", "0"));
        icon->set_geometry(pos, right, up, height, width);

        icon->set_pickable(get_bool(dict.get("pickable", "true")));

        return icon;
    }
    else if (elem_type == "icon_2d")
    {
        mi::base::Handle<nv::index::IIcon_2D> icon(
            m_scene->create_shape<nv::index::IIcon_2D>());

        mi::math::Vector<mi::Float32, 3> pos = get_vec_float32_3(dict.get("position", "0 0 0"));
        mi::math::Vector<mi::Float32, 2> right = get_vec_float32_2(dict.get("right", "1 0"));
        mi::math::Vector<mi::Float32, 2> up = get_vec_float32_2(dict.get("up", "0 1"));
        mi::Float32 height = get_float32(dict.get("height", "10"));
        mi::Float32 width  = get_float32(dict.get("width", "0"));
        icon->set_geometry(pos, right, up, height, width);

        icon->set_pickable(get_bool(dict.get("pickable", "true")));

        return icon;
    }
    //
    // Line path
    //
    else if (elem_type == "line_path_3d")
    {
        const mi::Uint32 demo_mode = get_uint32(dict.get("demo", "0"));

        // points and additional colors and radii per point.
        std::vector<mi::math::Vector_struct<mi::Float32, 3> > points;
        std::vector<mi::Float32> radii;
        std::vector<mi::math::Color_struct> colors;
        std::vector<mi::Uint32> material_ids;
        std::vector<mi::Uint32> colormap_ids;

        const std::string nb_points_str = dict.get("nb_points", "0");
        const mi::Uint32 nb_points = get_uint32(nb_points_str);
        const std::string radius_default = dict.get("radius", "5");

        // Read the points ...
        if(demo_mode > 0)
        {
            mi::Float32 x=0.f, y=0.f;
            // in demo mode it create a line path following a exponential function
            // with some randomness
            const mi::Sint32 demo_seed = get_sint32(dict.get("demo_seed", "12357"));
            const mi::Float32 demo_randomness = get_float32(dict.get("demo_randomness", "0"));
            const mi::Float32 demo_scale = get_float32(dict.get("demo_scale", "100"));
            const std::string demo_color_str = dict.get("demo_color", "0");
            const mi::Uint32 demo_map_size = get_uint32(dict.get("demo_map_size", "0"));

            bool has_colors = (demo_color_str != "0");
            bool has_colormap = (demo_map_size > 0);

            srand(demo_seed);
            for (mi::Uint32 i=0; i < nb_points; i++)
            {
                x = 7.f*( static_cast<mi::Float32>(i)/static_cast<mi::Float32>(nb_points) );

                switch(demo_mode)
                {
                  case 1:
                      y = expf(-x) + sinf(x);
                      break;
                  case 2:
                      y = expf(-x) + cosf(2.f*x);
                      break;
                  case 3:
                      y = expf(-x) + sinf(2.f*x)/(1.f + 2.f*x);
                      break;
                }

                mi::math::Vector_struct<mi::Float32, 3> vertex;
                vertex.z = (x + (demo_randomness*rand())/(RAND_MAX + 1.f))*demo_scale;
                vertex.x = (y + (demo_randomness*rand())/(RAND_MAX + 1.f))*demo_scale;
                vertex.y = ((demo_randomness*rand())/(RAND_MAX + 1.f))*demo_scale;

                points.push_back(vertex);

                if(has_colors)
                {
                    mi::math::Color_struct vertex_color;

                    if(demo_color_str == "random")
                    {
                        vertex_color.r = rand()/(RAND_MAX + 1.f);
                        vertex_color.g = rand()/(RAND_MAX + 1.f);
                        vertex_color.b = rand()/(RAND_MAX + 1.f);

                        mi::Float32 max = vertex_color.r > vertex_color.g ? vertex_color.r : vertex_color.g;
                        max = vertex_color.b > max ? vertex_color.b : max;

                        vertex_color.r /= max;
                        vertex_color.g /= max;
                        vertex_color.b /= max;
                        vertex_color.a = 1.f;
                    }
                    else if(demo_color_str == "jetmap")
                    {
                        mi::Float32 t = (i+0.5f)/static_cast<mi::Float32>(nb_points);
                        vertex_color = jetmap(t, 0.f, 1.f);
                    }
                    else if(demo_color_str == "gradient")
                    {
                        mi::Float32 t = (i+0.5f)/static_cast<mi::Float32>(nb_points);
                        t = (int)(t*8.f)/8.f;

                        vertex_color.r = t;
                        vertex_color.g = 1.f-t;
                        vertex_color.b = 0.75;
                        vertex_color.a = 1.f;

                    }

                    colors.push_back(vertex_color);
                }

                if(has_colormap)
                {
                    mi::Float32 t = (i+0.5f)/static_cast<mi::Float32>(nb_points);
                    mi::Uint32 index = static_cast<mi::Uint32>(t*demo_map_size + 0.5f);

                    colormap_ids.push_back(index);
                }

            }
        }
        else
        {
            for (mi::Uint32 i=0; i < nb_points; i++)
            {
                std::ostringstream v;
                v << "v" << i;
                points.push_back(get_vec_float32_3(dict.get(v.str(), "0 0 0")));
            }

            const mi::Uint32 nb_radii = get_uint32(dict.get("nb_radii", nb_points_str));
            for (mi::Uint32 i=0; i < nb_radii; i++)
            {
                std::ostringstream r;
                r << "r" << i;
                if (dict.is_defined(r.str()))
                    radii.push_back(get_float32(dict.get(r.str())));
            }

            const mi::Uint32 nb_colors = get_uint32(dict.get("nb_colors", nb_points_str));
            for (mi::Uint32 i=0; i < nb_colors; i++)
            {
                std::ostringstream c;
                c << "c" << i;
                if (dict.is_defined(c.str()))
                    colors.push_back(get_color(dict.get(c.str())));
            }

            const mi::Uint32 nb_mat_id = get_uint32(dict.get("nb_mat_id", nb_points_str));
            for (mi::Uint32 i=0; i < nb_mat_id; i++)
            {
                std::ostringstream id;
                id << "mat_id" << i;
                if (dict.is_defined(id.str()))
                    material_ids.push_back(get_uint32(dict.get(id.str())));
            }

            const mi::Uint32 nb_map_id = get_uint32(dict.get("nb_map_id", nb_points_str));
            for (mi::Uint32 i=0; i < nb_map_id; i++)
            {
                std::ostringstream id;
                id << "map_id" << i;
                if (dict.is_defined(id.str()))
                    colormap_ids.push_back(get_uint32(dict.get(id.str())));
            }
        }

        // Create shape (point set) scene element using the scene's factory
        mi::base::Handle<nv::index::IPath_3D> shape(m_scene->create_shape<nv::index::IPath_3D>());
        shape->set_radius(get_float32(radius_default));

        if (!points.empty())
            shape->set_points(&points[0], points.size());
        if (!radii.empty())
            shape->set_radii(&radii[0], radii.size());
        if (!colors.empty())
            shape->set_colors(&colors[0], colors.size());
        if (!material_ids.empty())
            shape->set_material_ids(&material_ids[0], material_ids.size());
        if (!colormap_ids.empty())
            shape->set_color_map_indexes(&colormap_ids[0], colormap_ids.size());

        return shape;
    }
    else if (elem_type == "line_path_2d")
    {
        const mi::Uint32 demo_mode = get_uint32(dict.get("demo", "0"));
        const std::string radius_default = dict.get("radius", "5");

        // points and additional colors and radii per point.
        std::vector<mi::math::Vector_struct<mi::Float32, 3> > points;
        std::vector<mi::Float32> radii;
        std::vector<mi::math::Color_struct> colors;
        std::vector<mi::Uint32> colormap_ids;

        const std::string nb_points_str = dict.get("nb_points", "0");
        const mi::Uint32 nb_points = get_uint32(nb_points_str);

        // Read the points ...
        if(demo_mode > 0)
        {
            mi::Float32 x=0.f, y=0.f, z=0.f, t=0.f;
            // in demo mode it create a line path following a exponential function
            // with some randomness
            const mi::Sint32 demo_seed = get_sint32(dict.get("demo_seed", "12357"));
            const mi::Float32 demo_randomness = get_float32(dict.get("demo_randomness", "0"));
            const mi::Float32 demo_scale = get_float32(dict.get("demo_scale", "100"));
            const std::string demo_color_str = dict.get("demo_color", "0");
            const std::string demo_width_str = dict.get("demo_width", "0");
            const mi::Uint32 demo_map_size = get_uint32(dict.get("demo_map_size", "0"));

            bool has_colors = (demo_color_str != "0");
            bool has_radii = (demo_width_str != "0");
            bool has_colormap = (demo_map_size > 0);
            const mi::Float32 radius = get_float32(radius_default);

            srand(demo_seed);
            for (mi::Uint32 i=0; i < nb_points; i++)
            {

                if(demo_mode <= 3) //path demos
                {
                    x = 7.f*( (mi::Float32)i/(mi::Float32)nb_points );
                    z = 0.f;

                    switch(demo_mode)
                    {
                      case 1:
                          y = expf(-x) + sinf(x);
                          break;
                      case 2:
                          y = expf(-x) + cosf(2.f*x);
                          break;
                      case 3:
                          y = expf(-x) + sinf(2.f*x)/(1.f + 2.f*x);
                          break;
                    }
                }
                else // loop demos
                {
                    t = (mi::Float32)i/(mi::Float32)(nb_points-1);
                    switch(demo_mode)
                    {
                      case 4: //ring
                          x = sinf(t*2.f*float(M_PI));
                          y = cosf(t*2.f*float(M_PI));
                          z = 0.f;
                          break;

                      case 5: // frame
                          if(t < 0.25f)
                          {
                              t = t/0.25f - 0.5f;
                              x = t;
                              y = -0.5f;
                          }
                          else if(t < 0.5f)
                          {
                              t = (t-0.25f)/0.25f - 0.5f;
                              x = 0.5f;
                              y = t;
                          }
                          else if(t < 0.75f)
                          {
                              t = (t-0.5f)/0.25f - 0.5f;
                              x = -t;
                              y = 0.5f;
                          }
                          else // t < 1.f
                          {
                              t = (t-0.75f)/0.25f - 0.5f;
                              x = -0.5f;
                              y = -t;
                          }

                          z = 0.f;
                          break;

                      case 6: //infinite symbol
                          x = sinf(2.f*float(M_PI)*t + float(M_PI)/2.f);
                          y = sinf(4.f*float(M_PI)*t);
                          z = sinf(2.f*float(M_PI)*t);
                          break;

                    }
                }

                mi::math::Vector_struct<mi::Float32, 3> vertex;
                vertex.z = (x + (demo_randomness*rand())/(RAND_MAX + 1.f))*demo_scale;
                vertex.x = (y + (demo_randomness*rand())/(RAND_MAX + 1.f))*demo_scale;
                vertex.y = (z + (demo_randomness*rand())/(RAND_MAX + 1.f))*demo_scale;

                points.push_back(vertex);

                if(has_radii)
                {
                    t = (mi::Float32)i/(mi::Float32)(nb_points-1);
                    if(demo_width_str == "ramp")
                    {
                        x = 1.f + (radius - 1.f)*t;
                    }
                    else if(demo_width_str == "bump")
                    {
                        x = (0.25f + 0.75f*sinf(t*float(M_PI)))* radius;
                    }
                    else if(demo_width_str == "bowl")
                    {
                        x = (0.25f + 0.75f*(1.f - sinf(t*float(M_PI))))* radius;
                    }
                    else if(demo_width_str == "wave")
                    {
                        x = (0.25f + 0.75f*((cosf(4*t*float(M_PI)) + 1.f)*0.5f))* radius;
                    }
                    radii.push_back(x);

                }

                if(has_colors)
                {
                    mi::math::Color_struct vertex_color;

                    if(demo_color_str == "random")
                    {
                        vertex_color.r = rand()/(RAND_MAX + 1.f);
                        vertex_color.g = rand()/(RAND_MAX + 1.f);
                        vertex_color.b = rand()/(RAND_MAX + 1.f);

                        mi::Float32 max = vertex_color.r > vertex_color.g ? vertex_color.r : vertex_color.g;
                        max = vertex_color.b > max ? vertex_color.b : max;

                        vertex_color.r /= max;
                        vertex_color.g /= max;
                        vertex_color.b /= max;
                        vertex_color.a = 1.f;
                    }
                    else if(demo_color_str == "jetmap")
                    {
                        mi::Float32 t = (i+0.5f)/static_cast<mi::Float32>(nb_points);
                        vertex_color = jetmap(t, 0.f, 1.f);
                    }
                    else if(demo_color_str == "gradient")
                    {
                        mi::Float32 t = (i+0.5f)/static_cast<mi::Float32>(nb_points);
                        t = (int)(t*8.f)/8.f;

                        vertex_color.r = t;
                        vertex_color.g = 1.f-t;
                        vertex_color.b = 0.75;
                        vertex_color.a = 1.f;

                    }

                    colors.push_back(vertex_color);
                }

                if(has_colormap)
                {
                    mi::Float32 t = (i+0.5f)/static_cast<mi::Float32>(nb_points);
                    mi::Uint32 index = static_cast<mi::Uint32>(t*demo_map_size + 0.5f);

                    colormap_ids.push_back(index);
                }
            }
        }
        else
        {
            for (mi::Uint32 i=0; i < nb_points; i++)
            {
                std::ostringstream v;
                v << "v" << i;
                points.push_back(get_vec_float32_3(dict.get(v.str(), "0 0 0")));
            }

            const mi::Uint32 nb_radii = get_uint32(dict.get("nb_radii", nb_points_str));
            for (mi::Uint32 i=0; i < nb_radii; i++)
            {
                std::ostringstream r;
                r << "r" << i;
                if (dict.is_defined(r.str()))
                    radii.push_back(get_float32(dict.get(r.str())));
            }

            const mi::Uint32 nb_colors = get_uint32(dict.get("nb_colors", nb_points_str));
            for (mi::Uint32 i=0; i < nb_colors; i++)
            {
                std::ostringstream c;
                c << "c" << i;
                if (dict.is_defined(c.str()))
                    colors.push_back(get_color(dict.get(c.str())));
            }

            const mi::Uint32 nb_map_id = get_uint32(dict.get("nb_map_id", nb_points_str));
            for (mi::Uint32 i=0; i < nb_map_id; i++)
            {
                std::ostringstream id2;
                id2 << "map_id" << i;
                if (dict.is_defined(id2.str()))
                    colormap_ids.push_back(get_uint32(dict.get(id2.str())));
            }
        }

        // Create shape (point set) scene element using the scene's factory
        mi::base::Handle<nv::index::IPath_2D> shape(m_scene->create_shape<nv::index::IPath_2D>());
        shape->set_radius(get_float32(radius_default));

        if (!points.empty())
            shape->set_points(&points[0], points.size());
        if (!radii.empty())
            shape->set_radii(&radii[0], radii.size());
        if (!colors.empty())
            shape->set_colors(&colors[0], colors.size());
        if (!colormap_ids.empty())
            shape->set_color_map_indexes(&colormap_ids[0], colormap_ids.size());

        return shape;
    }
    //
    // Raster benchmark: Create a large amount of primities randomly for raster benchmark purposes
    //
    else if (elem_type == "raster_benchmark_points")
    {
        // Vertices and additional per-vertex colors and radii.
        std::vector<mi::math::Vector_struct<mi::Float32, 3> > vertices;
        std::vector<mi::math::Color_struct> per_vertex_color_values;
        std::vector<mi::Float32> per_vertex_radii;

        const mi::Uint32 nb_points = get_uint32(dict.get("nb_points", "100"));

        mi::math::Vector_struct<mi::Float32, 3> vmin;
        mi::math::Vector_struct<mi::Float32, 3> vmax;
        vmin.x = -500; vmin.y = -500; vmin.z = -500;
        vmax.x = 500; vmax.y = 500; vmax.z = 500;

        // random points
        for (mi::Uint32 i=0; i < nb_points; ++i)
        {
            mi::math::Vector_struct<mi::Float32, 3> vertex;
            vertex.x = random(vmin.x, vmax.x);
            vertex.y = random(vmin.y, vmax.y);
            vertex.z = random(vmin.z, vmax.z);
            vertices.push_back(vertex);

            mi::math::Color_struct color;
            color.r = random(0.1f, 1.f);
            color.g = random(0.1f, 1.f);
            color.b = random(0.1f, 1.f);

            mi::Float32 cmax = color.r;
            if(color.g > cmax) cmax = color.g;
            if(color.b > cmax) cmax = color.b;

            color.r /= cmax;
            color.g /= cmax;
            color.b /= cmax;

            color.a = random(0.5, 1.f);
            per_vertex_color_values.push_back(color);

            per_vertex_radii.push_back(random(2.f, 25.f));         // all have the same radius
        }

        // Create point set shape and add it to the scene graph
        // Create shape (point set) scene element using the scene's factory
        mi::base::Handle<nv::index::IPoint_set> shape(m_scene->create_shape<nv::index::IPoint_set>());
        // ... and the vertices with colors and radii:
        shape->set_vertices(&vertices[0], vertices.size());
        shape->set_colors(&per_vertex_color_values[0], per_vertex_color_values.size());
        shape->set_radii(&per_vertex_radii[0], per_vertex_radii.size());

        nv::index::IPoint_set::Point_style style = nv::index::IPoint_set::SHADED_CIRCLE;
        shape->set_point_style(style);
        return shape;
    }
    else if (elem_type == "raster_benchmark_lines")
    {
        mi::math::Color color(0.f, 0.f, 1.f, 1.f);

        std::vector<mi::math::Vector_struct<mi::Float32, 3> > vertices;
        std::vector<mi::math::Color_struct> per_vertex_color_values;
        std::vector<mi::Float32> per_vertex_widths;

        const mi::Uint32 nb_lines = get_uint32(dict.get("nb_lines", "100"));

        mi::math::Vector_struct<mi::Float32, 3> vmin;
        mi::math::Vector_struct<mi::Float32, 3> vmax;
        vmin.x = -500; vmin.y = -500; vmin.z = -500;
        vmax.x = 500; vmax.y = 500; vmax.z = 500;

        for (mi::Uint32 j=0; j < nb_lines; ++j)
        {
            color.r = random(0.1f, 1.f);
            color.g = random(0.1f, 1.f);
            color.b = random(0.1f, 1.f);

            mi::Float32 cmax = color.r;
            if(color.g > cmax) cmax = color.g;
            if(color.b > cmax) cmax = color.b;

            color.r /= cmax;
            color.g /= cmax;
            color.b /= cmax;

            color.a = random(0.5, 1.f);

            mi::Float32 width = random(1.f, 5.f);

            for (mi::Uint32 i=0; i < 2; ++i)
            {
                mi::math::Vector_struct<mi::Float32, 3> vertex;
                vertex.x = random(vmin.x, vmax.x);
                vertex.y = random(vmin.y, vmax.y);
                vertex.z = random(vmin.z, vmax.z);
                vertices.push_back(vertex);

                per_vertex_color_values.push_back(color);   // all have the same color
                per_vertex_widths.push_back(width);         // all have the same radius
            }
        }

        // Create line set using the scene's factory
        mi::base::Handle<nv::index::ILine_set> shape(m_scene->create_shape<nv::index::ILine_set>());
        // ... set the line type. Currently only line segments are supported (see ILine_set)
        shape->set_line_type(nv::index::ILine_set::LINE_TYPE_SEGMENTS);
        // ... and the lines with colors and widths all in line segment order.
        shape->set_lines(&vertices[0], vertices.size());
        shape->set_colors(&per_vertex_color_values[0], per_vertex_color_values.size());
        shape->set_widths(&per_vertex_widths[0], per_vertex_widths.size());

        nv::index::ILine_set::Line_style style = nv::index::ILine_set::LINE_STYLE_SOLID;
        shape->set_line_style(style);

        return shape;
    }

    return null_result();
}
