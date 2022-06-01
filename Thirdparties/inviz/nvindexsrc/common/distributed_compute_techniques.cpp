/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/

#include "distributed_compute_techniques.h"

#include <fstream>
#include <sstream>
#include <vector>
#include <limits>

#include <nv/index/iregular_heightfield_patch.h>

#include "forwarding_logger.h"
#include "ppm_io.h"

namespace nv {
namespace index_common {

// Distributed_compute_checkerboard_2d /////////////////////////////////////////////////////////////////////////////////
Distributed_compute_checkerboard_2d::Distributed_compute_checkerboard_2d()
  : m_buffer_format(nv::index::IDistributed_compute_destination_buffer_2d_texture::RGBA_FLOAT32),
    m_scale(32.0f),
    m_color_0(mi::math::Color(0.9f, 0.2f, 0.15f, 0.7f)),
    m_color_1(mi::math::Color(0.2f, 0.2f, 0.8f,  1.0f)),
    m_enabled(true)
{
    update_configuration();
}

Distributed_compute_checkerboard_2d::Distributed_compute_checkerboard_2d(
        const mi::math::Vector_struct<mi::Float32, 2>&                                     extent,
              nv::index::IDistributed_compute_destination_buffer_2d_texture::Buffer_format fmt)
  : m_extent(extent),
    m_buffer_format(fmt),
    m_scale(32.0f),
    m_color_0(mi::math::Color(0.9f, 0.2f, 0.15f, 0.7f)),
    m_color_1(mi::math::Color(0.2f, 0.2f, 0.8f,  1.0f)),
    m_enabled(true)
{
    update_configuration();
}

void Distributed_compute_checkerboard_2d::launch_compute(
        mi::neuraylib::IDice_transaction*                   dice_transaction,
        nv::index::IDistributed_compute_destination_buffer* dst_buffer) const
{
    using nv::index::IDistributed_compute_destination_buffer_2d_texture;

    // retrieve 2d texture buffer destination buffer interface
    const mi::base::Handle<IDistributed_compute_destination_buffer_2d_texture> dst_buffer_2d(
        dst_buffer->get_interface<IDistributed_compute_destination_buffer_2d_texture>());

    if (!dst_buffer_2d)
    {
        ERROR_LOG << "Distributed compute technique (checkerboard): "
                  << "Unable to retrieve valid 2d texture destination buffer interface.";
        return;
    }

    // define the texture buffer format
    // * this example allows to use: INTENSITY_UINT8, RGBA_UINT8, RGBA_FLOAT32
    const IDistributed_compute_destination_buffer_2d_texture::Buffer_format     texture_format = m_buffer_format;

    // define the texture buffer layout
    // * this example allows to use: ONE_DIMENSIONAL_ARRAY, TWO_DIMENSIONAL_ARRAY
    // * the screen-space area (in raster space) that the subregion covers enables the user-defined generation to
    //   choose, for instance, the most efficient texture generation technique (for example, based on geometry area
    //   or window, or on ray/geometry intersection).
    const IDistributed_compute_destination_buffer_2d_texture::Buffer_layout     texture_layout = IDistributed_compute_destination_buffer_2d_texture::TWO_DIMENSIONAL_ARRAY;

    // set the buffer config to the destination buffer, initializing the internal buffer storage to the required size.
    IDistributed_compute_destination_buffer_2d_texture::Buffer_config   texture_config;
    texture_config.format = texture_format;
    texture_config.layout = texture_layout;

    if (texture_layout == IDistributed_compute_destination_buffer_2d_texture::ONE_DIMENSIONAL_ARRAY)
    {
        const mi::base::Handle<nv::index::IDistributed_compute_intersection_points>
                intersection_points(dst_buffer_2d->get_intersection_points());

        mi::Uint32 nb_intersections = 0;
        const mi::math::Vector_struct<mi::Float32, 3>* intersections
            = intersection_points->generate_intersection_points(nb_intersections);

        texture_config.covered_area = mi::math::Bbox<mi::Float32, 2>(0.0f, 0.0f, 0.0f, 0.0f);
        texture_config.resolution   = mi::math::Vector<mi::Uint32, 2>(nb_intersections, 0u);

        if (!dst_buffer_2d->generate_buffer_storage(texture_config))
        {
            ERROR_LOG << "Distributed compute technique (checkerboard): "
                      << "unable to generate destination buffer storage.";
            return;
        }

        generate_texture_tile_1D(dice_transaction, dst_buffer_2d.get(), intersections, nb_intersections);
    }
    else
    {
        const mi::math::Bbox<mi::Float32, 2>   surface_area = dst_buffer_2d->get_surface_area();

        // Round plane area to the next full texel
        mi::math::Vector<mi::Sint32, 2> lower(
            static_cast< mi::Sint32 >(mi::math::floor(surface_area.min.x / m_scale)),
            static_cast< mi::Sint32 >(mi::math::floor(surface_area.min.y / m_scale)));
        mi::math::Vector<mi::Sint32, 2> upper(
            static_cast< mi::Sint32 >(mi::math::ceil(surface_area.max.x / m_scale)),
            static_cast< mi::Sint32 >(mi::math::ceil(surface_area.max.y / m_scale)));

        // Need to add additional boundary pixels depending on filtering mode (see ITexture_filter_mode)
        lower -= mi::math::Vector<mi::Sint32, 2>(1);
        upper += mi::math::Vector<mi::Sint32, 2>(1);

        // Compute the tile resolution
        texture_config.resolution.x = mi::math::max(upper.x - lower.x, 1);
        texture_config.resolution.y = mi::math::max(upper.y - lower.y, 1);

        // Compute the area that is covered by the tile, which is probably larger than the requested
        // surface_area due to rounding to full texels.
        texture_config.covered_area.min.x = lower.x * m_scale;
        texture_config.covered_area.min.y = lower.y * m_scale;
        texture_config.covered_area.max.x = upper.x * m_scale;
        texture_config.covered_area.max.y = upper.y * m_scale;

        if (!dst_buffer_2d->generate_buffer_storage(texture_config))
        {
            ERROR_LOG << "Distributed compute technique (checkerboard): "
                      << "unable to generate destination buffer storage.";
            return;
        }

        generate_texture_tile_2D(dice_transaction, dst_buffer_2d.get(), lower, upper);
        
    }

    VERBOSE_LOG << "Distributed Compute: Checkerboard tile for " << dst_buffer->get_subregion_bbox() << ":"
                << " resolution "           << texture_config.resolution
                << " covering "             << texture_config.covered_area;
}

void Distributed_compute_checkerboard_2d::generate_texture_tile_1D(
    mi::neuraylib::IDice_transaction*                                    dice_transaction,
    const nv::index::IDistributed_compute_destination_buffer_2d_texture* dst_buffer,
    const mi::math::Vector_struct<mi::Float32, 3>*                       intersections,
    const mi::Uint32                                                     nb_intersections) const
{
    if (dst_buffer->get_buffer_config().format == nv::index::IDistributed_compute_destination_buffer_2d_texture::RGBA_FLOAT32)
    {
        mi::math::Color_struct* texture_tile =
            reinterpret_cast<mi::math::Color_struct*>(dst_buffer->get_buffer_storage());

        for (mi::Uint32 i=0; i < nb_intersections; ++i)
        {
            const mi::math::Vector_struct<mi::Float32, 3> point = intersections[i];

            // Skip if there is no intersection
            if (point.x != mi::base::numeric_traits<mi::Float32>::max())
            {
                mi::Uint32 x = static_cast< mi::Uint32 >(point.x / m_scale);
                mi::Uint32 y = static_cast< mi::Uint32 >(point.y / m_scale);

                if (x % 2 == y % 2)
                    texture_tile[i] = m_color_0;
                else
                    texture_tile[i] = m_color_1;
            }
        }
    }
    else if (dst_buffer->get_buffer_config().format == nv::index::IDistributed_compute_destination_buffer_2d_texture::RGBA_UINT8)
    {
        const mi::math::Vector<mi::Uint8, 4> color_0(
            static_cast< mi::Uint8 >(m_color_0.r * 255.f),
            static_cast< mi::Uint8 >(m_color_0.g * 255.f),
            static_cast< mi::Uint8 >(m_color_0.b * 255.f),
            static_cast< mi::Uint8 >(m_color_0.a * 255.f));
        const mi::math::Vector<mi::Uint8, 4> color_1(
            static_cast< mi::Uint8 >(m_color_1.r * 255.f),
            static_cast< mi::Uint8 >(m_color_1.g * 255.f),
            static_cast< mi::Uint8 >(m_color_1.b * 255.f),
            static_cast< mi::Uint8 >(m_color_1.a * 255.f));

        // Allocate the texture buffer
        mi::math::Vector<mi::Uint8, 4>* texture_tile =
            reinterpret_cast<mi::math::Vector<mi::Uint8, 4>*>(dst_buffer->get_buffer_storage());

        for (mi::Uint32 i=0; i < nb_intersections; ++i)
        {
            const mi::math::Vector_struct<mi::Float32, 3> point = intersections[i];

            // Skip if there is no intersection
            if (point.x != mi::base::numeric_traits<mi::Float32>::max())
            {
                mi::Uint32 x = static_cast< mi::Uint32 >(point.x / m_scale);
                mi::Uint32 y = static_cast< mi::Uint32 >(point.y / m_scale);

                if (x % 2 == y % 2)
                    texture_tile[i] = color_0;
                else
                    texture_tile[i] = color_1;
            }
        }
    }
    else if (dst_buffer->get_buffer_config().format == nv::index::IDistributed_compute_destination_buffer_2d_texture::INTENSITY_UINT8)
    {
        unsigned char* texture_tile =
            reinterpret_cast<unsigned char*>(dst_buffer->get_buffer_storage());

        for (mi::Uint32 i=0; i < nb_intersections; ++i)
        {
            const mi::math::Vector_struct<mi::Float32, 3> point = intersections[i];

            // Skip if there is no intersection
            if (point.x != mi::base::numeric_traits<mi::Float32>::max())
            {
                mi::Uint32 x = static_cast< mi::Uint32 >(point.x / m_scale);
                mi::Uint32 y = static_cast< mi::Uint32 >(point.y / m_scale);

                if (x % 2 == y % 2)
                    texture_tile[i] = 0;
                else
                    texture_tile[i] = 255;
            }
        }
    }
}

void Distributed_compute_checkerboard_2d::generate_texture_tile_2D(
    mi::neuraylib::IDice_transaction*                                    dice_transaction,
    const nv::index::IDistributed_compute_destination_buffer_2d_texture* dst_buffer,
    const mi::math::Vector<mi::Sint32, 2>&                               ll_corner,
    const mi::math::Vector<mi::Sint32, 2>&                               ur_corner) const
{
    const mi::math::Vector_struct<mi::Uint32, 2>&   tile_resolution = 
        dst_buffer->get_buffer_config().resolution;

    if (dst_buffer->get_buffer_config().format == nv::index::IDistributed_compute_destination_buffer_2d_texture::RGBA_FLOAT32)
    {
        mi::math::Color_struct* texture_tile =
            reinterpret_cast<mi::math::Color_struct*>(dst_buffer->get_buffer_storage());

        // Paint the checkerboard pattern
        mi::Uint32 index = 0;
        for (mi::Uint32 y=0; y < tile_resolution.y; ++y)
        {
            for (mi::Uint32 x=0; x < tile_resolution.x; ++x)
            {
                if ((ll_corner.x + x) % 2 == (ll_corner.y + y) % 2)
                    texture_tile[index] = m_color_0;
                else
                    texture_tile[index] = m_color_1;

                ++index;
            }
        }
    }
    else if (dst_buffer->get_buffer_config().format == nv::index::IDistributed_compute_destination_buffer_2d_texture::RGBA_UINT8)
    {
        const mi::math::Vector<mi::Uint8, 4> color_0(
            static_cast< mi::Uint8 >(m_color_0.r * 255.f),
            static_cast< mi::Uint8 >(m_color_0.g * 255.f),
            static_cast< mi::Uint8 >(m_color_0.b * 255.f),
            static_cast< mi::Uint8 >(m_color_0.a * 255.f));
        const mi::math::Vector<mi::Uint8, 4> color_1(
            static_cast< mi::Uint8 >(m_color_1.r * 255.f),
            static_cast< mi::Uint8 >(m_color_1.g * 255.f),
            static_cast< mi::Uint8 >(m_color_1.b * 255.f),
            static_cast< mi::Uint8 >(m_color_1.a * 255.f));

        // Allocate the texture buffer
        mi::math::Vector<mi::Uint8, 4>* texture_tile =
            reinterpret_cast<mi::math::Vector<mi::Uint8, 4>*>(dst_buffer->get_buffer_storage());

        // Paint the checkerboard pattern
        mi::Uint32 index = 0;
        for (mi::Uint32 y=0; y < tile_resolution.y; ++y)
        {
            for (mi::Uint32 x=0; x < tile_resolution.x; ++x)
            {
                if ((ll_corner.x + x) % 2 == (ll_corner.y + y) % 2)
                    texture_tile[index] = color_0;
                else
                    texture_tile[index] = color_1;

                ++index;
            }
        }
    }
    else if (dst_buffer->get_buffer_config().format == nv::index::IDistributed_compute_destination_buffer_2d_texture::INTENSITY_UINT8)
    {
        // Allocate the texture buffer
        unsigned char* texture_tile =
            reinterpret_cast<unsigned char*>(dst_buffer->get_buffer_storage());

        // Paint the checkerboard pattern
        mi::Uint32 index = 0;
        for (mi::Uint32 y=0; y < tile_resolution.y; ++y)
        {
            for (mi::Uint32 x=0; x < tile_resolution.x; ++x)
            {
                if ((ll_corner.x + x) % 2 == (ll_corner.y + y) % 2)
                {
                    // Make use of the entire value range to show the effect of the color map
                    texture_tile[index] =
                        static_cast< mi::Uint8 >(((ll_corner.x + x) / (m_extent.x / m_scale)) * 255.f);
                }
                else
                {
                    // Many color maps set opacity to zero for this value
                    texture_tile[index] = 127;
                }

                ++index;
            }
        }
    }
}

const char* Distributed_compute_checkerboard_2d::get_configuration() const
{
    return m_configuration.c_str();
}

void Distributed_compute_checkerboard_2d::set_enabled(bool enable)
{
    m_enabled = enable;
}

bool Distributed_compute_checkerboard_2d::get_enabled() const
{
    return m_enabled;
}

void Distributed_compute_checkerboard_2d::set_color_0(const mi::math::Color& color_0)
{
    m_color_0 = color_0;
    update_configuration();
}

void Distributed_compute_checkerboard_2d::set_color_1(const mi::math::Color& color_1)
{
    m_color_1 = color_1;
    update_configuration();
}

mi::neuraylib::IElement* Distributed_compute_checkerboard_2d::copy() const
{
    Distributed_compute_checkerboard_2d* other = new Distributed_compute_checkerboard_2d();

    other->m_extent        = this->m_extent;
    other->m_buffer_format = this->m_buffer_format;
    other->m_color_0       = this->m_color_0;
    other->m_color_1       = this->m_color_1;
    other->m_scale         = this->m_scale;
    other->m_configuration = this->m_configuration;
    other->m_enabled       = this->m_enabled;

    return other;
}

const char* Distributed_compute_checkerboard_2d::get_class_name() const
{
    return "Distributed_compute_checkerboard_2d";
}

void Distributed_compute_checkerboard_2d::serialize(mi::neuraylib::ISerializer* serializer) const
{
    serializer->write(&m_extent.x, 2);

    const mi::Uint32 bf = m_buffer_format;
    serializer->write(&bf, 1);

    serializer->write(&m_color_0.r, 4);
    serializer->write(&m_color_1.r, 4);
    serializer->write(&m_scale, 1);

    mi::Uint32 nb_elements = mi::Uint32(m_configuration.size());
    serializer->write(&nb_elements, 1);
    serializer->write(reinterpret_cast<const mi::Uint8*>(m_configuration.c_str()), nb_elements);

    serializer->write(&m_enabled, 1);
}

void Distributed_compute_checkerboard_2d::deserialize(mi::neuraylib::IDeserializer* deserializer)
{
    deserializer->read(&m_extent.x, 2);
    mi::Uint32 bf = 0u;
    deserializer->read(&bf, 1);

    m_buffer_format = nv::index::IDistributed_compute_destination_buffer_2d_texture::Buffer_format(bf);

    deserializer->read(&m_color_0.r, 4);
    deserializer->read(&m_color_1.r, 4);
    deserializer->read(&m_scale, 1);

    mi::Uint32 nb_elements;
    deserializer->read(&nb_elements, 1);
    m_configuration.resize(nb_elements);
    deserializer->read(reinterpret_cast<mi::Uint8*>(&m_configuration[0]), nb_elements);

    deserializer->read(&m_enabled, 1);
}

void Distributed_compute_checkerboard_2d::get_references(
    mi::neuraylib::ITag_set* result) const
{
}

void Distributed_compute_checkerboard_2d::update_configuration()
{
    using nv::index::IDistributed_compute_destination_buffer_2d_texture;

    std::string fmt;
    switch (m_buffer_format)
    {
        case IDistributed_compute_destination_buffer_2d_texture::INTENSITY_UINT8:   fmt.assign("intensity_uint8"); break;
        case IDistributed_compute_destination_buffer_2d_texture::RGBA_UINT8:        fmt.assign("rgba_uint8"); break;
        case IDistributed_compute_destination_buffer_2d_texture::RGBA_FLOAT32:      // fall through default
        default:
            fmt.assign("rgba_float32"); break;
    }

    std::ostringstream s;
    s << "generator=checkerboard\n"
      << "resolution=" << m_extent.x << " " << m_extent.y << "\n"
      << "format=" << fmt << "\n"
      << "color0=" << m_color_0.r << " " << m_color_0.g << " " << m_color_0.b << " " << m_color_0.a << "\n"
      << "color1=" << m_color_1.r << " " << m_color_1.g << " " << m_color_1.b << " " << m_color_1.a << "\n";
    m_configuration = s.str();
}

// Distributed_compute_mandelbrot_2d ///////////////////////////////////////////////////////////////////////////////////
Distributed_compute_mandelbrot_2d::Distributed_compute_mandelbrot_2d()
  : m_enabled(true)
{
}

Distributed_compute_mandelbrot_2d::Distributed_compute_mandelbrot_2d(
        const mi::math::Vector_struct<mi::Float32, 2>& extent)
  : m_extent(extent),
    m_enabled(true)
{
    std::ostringstream s;
    s << "generator=mandelbrot\n"
      << "resolution=" << extent.x << " " << extent.y << "\n";
    m_configuration = s.str();
}

void Distributed_compute_mandelbrot_2d::launch_compute(
        mi::neuraylib::IDice_transaction*                   dice_transaction,
        nv::index::IDistributed_compute_destination_buffer* dst_buffer) const
{
    using nv::index::IDistributed_compute_destination_buffer_2d_texture;

    // retrieve 2d texture buffer destination buffer interface
    const mi::base::Handle<IDistributed_compute_destination_buffer_2d_texture> dst_buffer_2d(
        dst_buffer->get_interface<IDistributed_compute_destination_buffer_2d_texture>());

    if (!dst_buffer_2d)
    {
        ERROR_LOG << "Distributed compute technique (mandelbrot): "
                  << "Unable to retrieve valid 2d texture destination buffer interface.";
        return;
    }

    // define the texture buffer format
    // * this example allows to use: RGBA_FLOAT32
    const IDistributed_compute_destination_buffer_2d_texture::Buffer_format     texture_format = IDistributed_compute_destination_buffer_2d_texture::RGBA_FLOAT32;

    // define the texture buffer layout
    // * this example allows to use: ONE_DIMENSIONAL_ARRAY, TWO_DIMENSIONAL_ARRAY
    // * the screen-space area (in raster space) that the subregion covers enables the user-defined generation to
    //   choose, for instance, the most efficient texture generation technique (for example, based on geometry area
    //   or window, or on ray/geometry intersection).
    const IDistributed_compute_destination_buffer_2d_texture::Buffer_layout     texture_layout = IDistributed_compute_destination_buffer_2d_texture::TWO_DIMENSIONAL_ARRAY;

    // set the buffer config to the destination buffer, initializing the internal buffer storage to the required size.
    IDistributed_compute_destination_buffer_2d_texture::Buffer_config   texture_config;
    texture_config.format = texture_format;
    texture_config.layout = texture_layout;

    if (texture_layout == IDistributed_compute_destination_buffer_2d_texture::ONE_DIMENSIONAL_ARRAY)
    {
        const mi::base::Handle<nv::index::IDistributed_compute_intersection_points>
                intersection_points(dst_buffer_2d->get_intersection_points());

        mi::Uint32 nb_intersections = 0;
        const mi::math::Vector_struct<mi::Float32, 3>* intersections
            = intersection_points->generate_intersection_points(nb_intersections);

        texture_config.covered_area = mi::math::Bbox<mi::Float32, 2>(0.0f, 0.0f, 0.0f, 0.0f);
        texture_config.resolution   = mi::math::Vector<mi::Uint32, 2>(nb_intersections, 0u);

        if (!dst_buffer_2d->generate_buffer_storage(texture_config))
        {
            ERROR_LOG << "Distributed compute technique (mandelbrot): "
                      << "unable to generate destination buffer storage.";
            return;
        }

        generate_texture_tile_1D(dice_transaction, dst_buffer_2d.get(), intersections, nb_intersections);
    }
    else
    {
        const mi::math::Bbox<mi::Float32, 2>   surface_area = dst_buffer_2d->get_surface_area();

        // Round plane area to the next full texel
        mi::math::Vector<mi::Sint32, 2> lower(
            static_cast< mi::Sint32 >(mi::math::floor(surface_area.min.x)),
            static_cast< mi::Sint32 >(mi::math::floor(surface_area.min.y)));
        mi::math::Vector<mi::Sint32, 2> upper(
            static_cast< mi::Sint32 >(mi::math::ceil(surface_area.max.x)),
            static_cast< mi::Sint32 >(mi::math::ceil(surface_area.max.y)));

        // Need to add additional boundary pixels depending on filtering mode (see ITexture_filter_mode)
        lower -= mi::math::Vector<mi::Sint32, 2>(1);
        upper += mi::math::Vector<mi::Sint32, 2>(1);

        // Compute the tile resolution
        texture_config.resolution.x = mi::math::max(upper.x - lower.x, 1);
        texture_config.resolution.y = mi::math::max(upper.y - lower.y, 1);

        // Compute the area that is covered by the tile, which is probably larger than the requested
        // surface_area due to rounding to full texels.
        texture_config.covered_area.min.x = static_cast<mi::Float32>(lower.x);
        texture_config.covered_area.min.y = static_cast<mi::Float32>(lower.y);
        texture_config.covered_area.max.x = static_cast<mi::Float32>(upper.x);
        texture_config.covered_area.max.y = static_cast<mi::Float32>(upper.y);

        if (!dst_buffer_2d->generate_buffer_storage(texture_config))
        {
            ERROR_LOG << "Distributed compute technique (mandelbrot): "
                      << "unable to generate destination buffer storage.";
            return;
        }

        generate_texture_tile_2D(dice_transaction, dst_buffer_2d.get(), lower, upper);
        
    }

    VERBOSE_LOG << "Distributed Compute: Mandelbrot tile for " << dst_buffer->get_subregion_bbox() << ":"
                << " resolution "           << texture_config.resolution
                << " covering "             << texture_config.covered_area;
}

void Distributed_compute_mandelbrot_2d::generate_texture_tile_1D(
    mi::neuraylib::IDice_transaction*                                    dice_transaction,
    const nv::index::IDistributed_compute_destination_buffer_2d_texture* dst_buffer,
    const mi::math::Vector_struct<mi::Float32, 3>*                       intersections,
    const mi::Uint32                                                     nb_intersections) const
{
    mi::math::Color_struct* texture_tile =
        reinterpret_cast<mi::math::Color_struct*>(dst_buffer->get_buffer_storage());

    for (mi::Uint32 i=0; i < nb_intersections; ++i)
    {
        const mi::math::Vector_struct<mi::Float32, 3> point = intersections[i];

        // Skip if there is no intersection
        if (point.x != mi::base::numeric_traits<mi::Float32>::max())
        {
            // Convert point from [0.0, extent] to [0.0, 1.0]
            const mi::math::Vector<mi::Float32, 2> p(
                point.x / m_extent.x,
                point.y / m_extent.y);

            // Calculate color value
            texture_tile[i] = mandelbrot(p);
        }
    }
}

void Distributed_compute_mandelbrot_2d::generate_texture_tile_2D(
    mi::neuraylib::IDice_transaction*                                    dice_transaction,
    const nv::index::IDistributed_compute_destination_buffer_2d_texture* dst_buffer,
    const mi::math::Vector<mi::Sint32, 2>&                               ll_corner,
    const mi::math::Vector<mi::Sint32, 2>&                               ur_corner) const
{
    const mi::math::Vector_struct<mi::Uint32, 2>&   tile_resolution = 
        dst_buffer->get_buffer_config().resolution;

    mi::math::Color_struct* texture_tile =
        reinterpret_cast<mi::math::Color_struct*>(dst_buffer->get_buffer_storage());

    // Paint the checkerboard pattern
    mi::Uint32 index = 0;
    for (mi::Uint32 y=0; y < tile_resolution.y; ++y)
    {
        for (mi::Uint32 x=0; x < tile_resolution.x; ++x)
        {
            texture_tile[index] = mandelbrot(
                mi::math::Vector<mi::Float32, 2>(
                    static_cast<mi::Float32>(ll_corner.x + x) / m_extent.x,
                    static_cast<mi::Float32>(ll_corner.y + y) / m_extent.y));

            ++index;
        }
    }
}

const char* Distributed_compute_mandelbrot_2d::get_configuration() const
{
    return m_configuration.c_str();
}

void Distributed_compute_mandelbrot_2d::set_enabled(bool enable)
{
    m_enabled = enable;
}

bool Distributed_compute_mandelbrot_2d::get_enabled() const
{
    return m_enabled;
}

mi::neuraylib::IElement* Distributed_compute_mandelbrot_2d::copy() const
{
    Distributed_compute_mandelbrot_2d* other = new Distributed_compute_mandelbrot_2d();
    other->m_extent        = this->m_extent;
    other->m_configuration = this->m_configuration;
    other->m_enabled       = this->m_enabled;
    return other;
}

const char* Distributed_compute_mandelbrot_2d::get_class_name() const
{
    return "Distributed_compute_mandelbrot_2d";
}

void Distributed_compute_mandelbrot_2d::serialize(mi::neuraylib::ISerializer* serializer) const
{
    serializer->write(&m_extent.x, 2);

    mi::Uint32 nb_elements = mi::Uint32(m_configuration.size());
    serializer->write(&nb_elements, 1);
    serializer->write(reinterpret_cast<const mi::Uint8*>(m_configuration.c_str()), nb_elements);

    serializer->write(&m_enabled, 1);
}

void Distributed_compute_mandelbrot_2d::deserialize(mi::neuraylib::IDeserializer* deserializer)
{
    deserializer->read(&m_extent.x, 2);

    mi::Uint32 nb_elements;
    deserializer->read(&nb_elements, 1);
    m_configuration.resize(nb_elements);
    deserializer->read(reinterpret_cast<mi::Uint8*>(&m_configuration[0]), nb_elements);

    deserializer->read(&m_enabled, 1);
}

mi::math::Color Distributed_compute_mandelbrot_2d::mandelbrot(
    const mi::math::Vector<mi::Float32, 2>& tc) const
{
    mi::Float32 x0 = tc.x * 3.5f - 2.5f;
    mi::Float32 y0 = tc.y * 2.f - 1.f;
    mi::Float32 x = 0.f;
    mi::Float32 y = 0.f;

    mi::Sint32 iteration = 0;
    const mi::Sint32 MAX_ITERATION = 40;

    while (x*x + y*y < 4 && iteration < MAX_ITERATION)
    {
        mi::Float32 tmp = x*x - y*y + x0;
        y = 2 * x*y + y0;
        x = tmp;
        iteration++;
    }

    mi::math::Color color;
    if (iteration < 8)
        color = mi::math::Color(0.f, 0.f, (mi::Float32)iteration / 8, 1.f);
    else
        color = mi::math::Color(0.f, (mi::Float32)iteration / MAX_ITERATION, 0.f, 1.f);

    return color;
}


// Distributed_compute_bitmap_mapping_2d ///////////////////////////////////////////////////////////////////////////////
Distributed_compute_bitmap_mapping_2d::Distributed_compute_bitmap_mapping_2d()
  : m_enabled(true),
    m_use_cache(false)
{
}

Distributed_compute_bitmap_mapping_2d::Distributed_compute_bitmap_mapping_2d(
    const std::string&                              bm_filemane,
    bool                                            use_cache,
    const mi::math::Vector_struct<mi::Float32, 2>&  extent)
  : m_bm_filemane(bm_filemane),
    m_extent(extent),
    m_enabled(true),
    m_use_cache(use_cache)
{
    INFO_LOG << "BITMAP TEXTURE|filename : " << m_bm_filemane;
}

void Distributed_compute_bitmap_mapping_2d::launch_compute(
        mi::neuraylib::IDice_transaction*                   dice_transaction,
        nv::index::IDistributed_compute_destination_buffer* dst_buffer) const
{
    using nv::index::IDistributed_compute_destination_buffer_2d_texture;

    // retrieve 2d texture buffer destination buffer interface
    const mi::base::Handle<IDistributed_compute_destination_buffer_2d_texture> dst_buffer_2d(
        dst_buffer->get_interface<IDistributed_compute_destination_buffer_2d_texture>());

    if (!dst_buffer_2d)
    {
        ERROR_LOG << "Distributed compute technique (bitmap): "
                  << "Unable to retrieve valid 2d texture destination buffer interface.";
        return;
    }

    // define the texture buffer format
    // * this example allows to use: RGBA_FLOAT32
    const IDistributed_compute_destination_buffer_2d_texture::Buffer_format     texture_format = IDistributed_compute_destination_buffer_2d_texture::RGBA_FLOAT32;

    // define the texture buffer layout
    // * this example allows to use: ONE_DIMENSIONAL_ARRAY, TWO_DIMENSIONAL_ARRAY
    // * the screen-space area (in raster space) that the subregion covers enables the user-defined generation to
    //   choose, for instance, the most efficient texture generation technique (for example, based on geometry area
    //   or window, or on ray/geometry intersection).
    const IDistributed_compute_destination_buffer_2d_texture::Buffer_layout     texture_layout = IDistributed_compute_destination_buffer_2d_texture::TWO_DIMENSIONAL_ARRAY;

    // set the buffer config to the destination buffer, initializing the internal buffer storage to the required size.
    IDistributed_compute_destination_buffer_2d_texture::Buffer_config   texture_config;
    texture_config.format = texture_format;
    texture_config.layout = texture_layout;

    if (texture_layout == IDistributed_compute_destination_buffer_2d_texture::ONE_DIMENSIONAL_ARRAY)
    {
        const mi::base::Handle<nv::index::IDistributed_compute_intersection_points>
                intersection_points(dst_buffer_2d->get_intersection_points());

        mi::Uint32 nb_intersections = 0;
        const mi::math::Vector_struct<mi::Float32, 3>* intersections
            = intersection_points->generate_intersection_points(nb_intersections);

        texture_config.covered_area = mi::math::Bbox<mi::Float32, 2>(0.0f, 0.0f, 0.0f, 0.0f);
        texture_config.resolution   = mi::math::Vector<mi::Uint32, 2>(nb_intersections, 0u);

        if (!dst_buffer_2d->generate_buffer_storage(texture_config))
        {
            ERROR_LOG << "Distributed compute technique (bitmap): "
                      << "unable to generate destination buffer storage.";
            return;
        }

        generate_texture_tile_1D(dice_transaction, dst_buffer_2d.get(), intersections, nb_intersections);
    }
    else
    {
        const mi::math::Bbox<mi::Float32, 2>   surface_area = dst_buffer_2d->get_surface_area();

        // Round plane area to the next full texel
        mi::math::Vector<mi::Sint32, 2> lower(
            static_cast< mi::Sint32 >(mi::math::floor(surface_area.min.x)),
            static_cast< mi::Sint32 >(mi::math::floor(surface_area.min.y)));
        mi::math::Vector<mi::Sint32, 2> upper(
            static_cast< mi::Sint32 >(mi::math::ceil(surface_area.max.x)),
            static_cast< mi::Sint32 >(mi::math::ceil(surface_area.max.y)));

        // Need to add additional boundary pixels depending on filtering mode (see ITexture_filter_mode)
        lower -= mi::math::Vector<mi::Sint32, 2>(1);
        upper += mi::math::Vector<mi::Sint32, 2>(1);
        
        // Try importing texture from a cache file when available
        std::string cache_file_name;
        std::ifstream icache_file;

        if (m_use_cache)
        {
            std::ostringstream cache_file_name_str;
            cache_file_name_str << m_bm_filemane << "_"
                                << lower.x << "." << lower.y << "-"
                                << upper.x << "." << upper.y << "_texture_importer.tmp";
           
            cache_file_name = cache_file_name_str.str();
            icache_file.open(cache_file_name.c_str(), std::ios_base::in | std::ios_base::binary);
            
            if(!icache_file)
                VERBOSE_LOG << "Cache file now available yet. Loading bitmap from original source";
        }
        
        if(icache_file.is_open())
        {
            VERBOSE_LOG << "Loading bitmap cache file: " << cache_file_name;
            
            icache_file.read(reinterpret_cast<char*>(&texture_config.resolution.x), sizeof(mi::Uint32));
            icache_file.read(reinterpret_cast<char*>(&texture_config.resolution.y), sizeof(mi::Uint32));
            icache_file.read(reinterpret_cast<char*>(&texture_config.covered_area.min.x), sizeof(mi::Float32));
            icache_file.read(reinterpret_cast<char*>(&texture_config.covered_area.min.y), sizeof(mi::Float32));
            icache_file.read(reinterpret_cast<char*>(&texture_config.covered_area.max.x), sizeof(mi::Float32));
            icache_file.read(reinterpret_cast<char*>(&texture_config.covered_area.max.y), sizeof(mi::Float32));
            
            if (!dst_buffer_2d->generate_buffer_storage(texture_config))
            {
                ERROR_LOG << "Distributed compute technique (bitmap): "
                          << "unable to generate destination buffer storage.";
                return;
            }

            mi::math::Color_struct* texture_tile =
                reinterpret_cast<mi::math::Color_struct*>(dst_buffer_2d->get_buffer_storage());
            
            icache_file.read(reinterpret_cast<char*>(texture_tile), sizeof(mi::math::Color_struct)*texture_config.resolution.x*texture_config.resolution.y);
            icache_file.close();
            
        }
        else
        {
            mi::math::Vector<mi::Float32, 2> lower_fl(
                static_cast<mi::Float32>(lower.x)/m_extent.x,
                static_cast<mi::Float32>(lower.y)/m_extent.y);
            mi::math::Vector<mi::Float32, 2> upper_fl(
                static_cast<mi::Float32>(upper.x)/m_extent.x,
                static_cast<mi::Float32>(upper.y)/m_extent.y);
        
            // Render the bitmap
            std::vector< mi::math::Color_struct > ppm_color_st_buf;
            mi::Sint32 img_width  = -1;
            mi::Sint32 img_height = -1;
            std::string error_mes; 
            if(!nv::index_common::load_ppm(m_bm_filemane, ppm_color_st_buf, img_width, img_height, error_mes))
            {
                ERROR_LOG << "bitmap failed to load[" << m_bm_filemane << "].\n" << error_mes;
                return;
            }             
         
            // Compute the tile resolution
            mi::Sint32 min_i = static_cast<mi::Sint32>(mi::math::floor(lower_fl.x*img_width));
            if(min_i < 0) min_i = 0;
        
            mi::Sint32 max_i = static_cast<mi::Sint32>(mi::math::ceil(upper_fl.x*img_width));
            if(max_i > img_width) max_i = img_width;
        
            mi::Sint32 min_j = static_cast<mi::Sint32>(mi::math::floor(lower_fl.y*img_height));
            if(min_j < 0) min_j = 0;
        
            mi::Sint32 max_j = static_cast<mi::Sint32>(mi::math::ceil(upper_fl.y*img_height));
            if(max_j > img_height) max_j = img_height;
        
            // Compute the tile resolution
            texture_config.resolution.x = static_cast<mi::Uint32>(max_i - min_i);
            texture_config.resolution.y = static_cast<mi::Uint32>(max_j - min_j);

            INFO_LOG << "Image dimesions: " << img_width << ", " << img_height;
            INFO_LOG << "MIN_MAX[(" << min_i << ", " << min_j << ") (" << max_i << ", " << max_j << ")]";
            INFO_LOG << "Tile resolution: " << texture_config.resolution;

            // Compute the area that is covered by the tile, which is probably larger than the requested
            // surface_area due to rounding to full texels.
            texture_config.covered_area.min.x = static_cast<mi::Float32>(min_i)/static_cast<mi::Float32>(img_width)  * m_extent.x;
            texture_config.covered_area.min.y = static_cast<mi::Float32>(min_j)/static_cast<mi::Float32>(img_height) * m_extent.y;
            texture_config.covered_area.max.x = static_cast<mi::Float32>(max_i)/static_cast<mi::Float32>(img_width)  * m_extent.x;;
            texture_config.covered_area.max.y = static_cast<mi::Float32>(max_j)/static_cast<mi::Float32>(img_height) * m_extent.y;

            if (!dst_buffer_2d->generate_buffer_storage(texture_config))
            {
                ERROR_LOG << "Distributed compute technique (bitmap): "
                          << "unable to generate destination buffer storage.";
                return;
            }

            mi::math::Color_struct* texture_tile =
                reinterpret_cast<mi::math::Color_struct*>(dst_buffer_2d->get_buffer_storage());

            mi::Uint32 index=0;
            for (mi::Uint32 j=0; j < texture_config.resolution.y; ++j)
            {
                mi::Uint32 img_idx = (j+min_j)*img_width + min_i;
                for (mi::Uint32 i=0; i < texture_config.resolution.x; ++i, ++img_idx, ++index)
                {
                    texture_tile[index] = ppm_color_st_buf[img_idx];
                }
            }
            
            if(m_use_cache)
            {
                std::ofstream ocache_file;
                ocache_file.open(cache_file_name.c_str(), std::ios_base::out | std::ios_base::binary);
                if(!ocache_file)
                {
                    ERROR_LOG << "It could not open bitmap cache file: " << cache_file_name;
                }
                else
                {
                    ocache_file.write(reinterpret_cast<char*>(&texture_config.resolution.x), sizeof(mi::Uint32));
                    ocache_file.write(reinterpret_cast<char*>(&texture_config.resolution.y), sizeof(mi::Uint32));
                    ocache_file.write(reinterpret_cast<char*>(&texture_config.covered_area.min.x), sizeof(mi::Float32));
                    ocache_file.write(reinterpret_cast<char*>(&texture_config.covered_area.min.y), sizeof(mi::Float32));
                    ocache_file.write(reinterpret_cast<char*>(&texture_config.covered_area.max.x), sizeof(mi::Float32));
                    ocache_file.write(reinterpret_cast<char*>(&texture_config.covered_area.max.y), sizeof(mi::Float32));
                    ocache_file.write(reinterpret_cast<char*>(texture_tile), sizeof(mi::math::Color_struct)*texture_config.resolution.x*texture_config.resolution.y);
                    ocache_file.close();
                }
                
            }
        }
    }

    VERBOSE_LOG << "Distributed Compute: Bitmap tile for " << dst_buffer->get_subregion_bbox() << ":"
                << " resolution "     << texture_config.resolution
                << " covering "       << texture_config.covered_area;
}

void Distributed_compute_bitmap_mapping_2d::generate_texture_tile_1D(
    mi::neuraylib::IDice_transaction*                                    dice_transaction,
    const nv::index::IDistributed_compute_destination_buffer_2d_texture* dst_buffer,
    const mi::math::Vector_struct<mi::Float32, 3>*                       intersections,
    const mi::Uint32                                                     nb_intersections) const
{
    mi::math::Color_struct* texture_tile =
        reinterpret_cast<mi::math::Color_struct*>(dst_buffer->get_buffer_storage());

    // Render the bitmap
    std::vector< mi::math::Color_struct > ppm_color_st_buf;
    mi::Sint32 img_width  = -1;
    mi::Sint32 img_height = -1;
    std::string error_mes;

    if(!nv::index_common::load_ppm(m_bm_filemane, ppm_color_st_buf, img_width, img_height, error_mes))
    {
        ERROR_LOG << "texture failed to load[" << m_bm_filemane << "].\n" << error_mes;
        return;
    }             

    INFO_LOG << "BITMAP TEXTURE|width: " << img_width;
    INFO_LOG << "BITMAP TEXTURE|height: " << img_height;

    for (mi::Uint32 i=0; i < nb_intersections; ++i)
    {
        const mi::math::Vector_struct<mi::Float32, 3> point = intersections[i];

        // Skip if there is no intersection
        if (point.x != mi::base::numeric_traits<mi::Float32>::max())
        {
            mi::Uint32 index = (mi::Uint32)point.y*img_width + (mi::Uint32)point.x;

            // Calculate color value
            texture_tile[i] = ppm_color_st_buf[index];
        }
    }
}

void Distributed_compute_bitmap_mapping_2d::set_enabled(bool enable)
{
    m_enabled = enable;
}

bool Distributed_compute_bitmap_mapping_2d::get_enabled() const
{
    return m_enabled;
}

mi::neuraylib::IElement* Distributed_compute_bitmap_mapping_2d::copy() const
{
    Distributed_compute_bitmap_mapping_2d* other = new Distributed_compute_bitmap_mapping_2d();
    other->m_bm_filemane  = this->m_bm_filemane;
    other->m_extent       = this->m_extent;
    other->m_use_cache    = this->m_use_cache;
    other->m_enabled      = this->m_enabled;
    return other;
}

const char* Distributed_compute_bitmap_mapping_2d::get_class_name() const
{
    return "Distributed_compute_bitmap_mapping_2d";
}

void Distributed_compute_bitmap_mapping_2d::serialize(mi::neuraylib::ISerializer* serializer) const
{
    mi::Uint32 nb_elements = mi::Uint32(m_bm_filemane.size());
    serializer->write(&nb_elements, 1);
    serializer->write(reinterpret_cast<const mi::Uint8*>(&m_bm_filemane[0]), nb_elements); 
    
    serializer->write(&m_extent.x, 2);
    
    serializer->write(&m_use_cache, 1);
    serializer->write(&m_enabled, 1);
}

void Distributed_compute_bitmap_mapping_2d::deserialize(mi::neuraylib::IDeserializer* deserializer)
{
    mi::Uint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_bm_filemane.resize(nb_elements);
    deserializer->read(reinterpret_cast<mi::Uint8*>(&m_bm_filemane[0]), nb_elements); 
    
    deserializer->read(&m_extent.x, 2);

    deserializer->read(&m_use_cache, 1);
    deserializer->read(&m_enabled, 1);
}

// Distributed_compute_height_color_mapping_2d //////////////////////////////////////////////////////////////////////////
Distributed_compute_height_color_mapping_2d::Distributed_compute_height_color_mapping_2d()
  : m_enabled(true)
{
}

Distributed_compute_height_color_mapping_2d::Distributed_compute_height_color_mapping_2d(
        const mi::math::Vector<mi::Float32, 2>& height_range)
  : m_height_range(height_range),
    m_color0(0.f, 1.f, 0.f, 1.f), // green
    m_color1(1.f, 0.f, 0.f, 1.f), // red
    m_enabled(true)
{
}

void Distributed_compute_height_color_mapping_2d::launch_compute(
        mi::neuraylib::IDice_transaction*                   dice_transaction,
        nv::index::IDistributed_compute_destination_buffer* dst_buffer) const
{
    using nv::index::IDistributed_compute_destination_buffer_2d_texture;

    // retrieve 2d texture buffer destination buffer interface
    const mi::base::Handle<IDistributed_compute_destination_buffer_2d_texture> dst_buffer_2d(
        dst_buffer->get_interface<IDistributed_compute_destination_buffer_2d_texture>());

    if (!dst_buffer_2d)
    {
        ERROR_LOG << "Distributed compute technique (bitmap): "
                  << "Unable to retrieve valid 2d texture destination buffer interface.";
        return;
    }

    // define the texture buffer format
    // * this example allows to use: RGBA_FLOAT32
    const IDistributed_compute_destination_buffer_2d_texture::Buffer_format     texture_format = IDistributed_compute_destination_buffer_2d_texture::RGBA_FLOAT32;

    // define the texture buffer layout
    // * this example allows to use: ONE_DIMENSIONAL_ARRAY, TWO_DIMENSIONAL_ARRAY
    // * the screen-space area (in raster space) that the subregion covers enables the user-defined generation to
    //   choose, for instance, the most efficient texture generation technique (for example, based on geometry area
    //   or window, or on ray/geometry intersection).
    const IDistributed_compute_destination_buffer_2d_texture::Buffer_layout     texture_layout = IDistributed_compute_destination_buffer_2d_texture::TWO_DIMENSIONAL_ARRAY;

    // set the buffer config to the destination buffer, initializing the internal buffer storage to the required size.
    IDistributed_compute_destination_buffer_2d_texture::Buffer_config   texture_config;
    texture_config.format = texture_format;
    texture_config.layout = texture_layout;

    const mi::math::Bbox<mi::Float32, 2>   surface_area = dst_buffer_2d->get_surface_area();

    // Round plane area to the next full texel
    mi::math::Vector<mi::Sint32, 2> lower(
        static_cast< mi::Sint32 >(mi::math::floor(surface_area.min.x)),
        static_cast< mi::Sint32 >(mi::math::floor(surface_area.min.y)));
    mi::math::Vector<mi::Sint32, 2> upper(
        static_cast< mi::Sint32 >(mi::math::ceil(surface_area.max.x)),
        static_cast< mi::Sint32 >(mi::math::ceil(surface_area.max.y)));

    // Need to add additional boundary pixels depending on filtering mode (see ITexture_filter_mode)
    //lower -= mi::math::Vector<mi::Sint32, 2>(1);
    //upper += mi::math::Vector<mi::Sint32, 2>(1);

    // Compute the tile resolution
    texture_config.resolution.x = mi::math::max(upper.x - lower.x, 1);
    texture_config.resolution.y = mi::math::max(upper.y - lower.y, 1);

    texture_config.covered_area.min.x = static_cast<mi::Float32>(lower.x);
    texture_config.covered_area.min.y = static_cast<mi::Float32>(lower.y);
    texture_config.covered_area.max.x = static_cast<mi::Float32>(upper.x);
    texture_config.covered_area.max.y = static_cast<mi::Float32>(upper.y);

    if (!dst_buffer_2d->generate_buffer_storage(texture_config))
    {
        ERROR_LOG << "Distributed compute technique (height color-mapping): "
                  << "unable to generate destination buffer storage.";
        return;
    }

    mi::math::Color_struct* texture_tile_rgba =
        reinterpret_cast<mi::math::Color_struct*>(dst_buffer_2d->get_buffer_storage());

    // Retrieve patch data for the current subregion
    const mi::Float32* patch_array = 0;
    mi::Sint32 patch_array_size = 0;
    mi::math::Bbox<mi::Sint32, 2> patch_bbox;
    mi::math::Vector<mi::Sint32, 2> patch_size(0);

    if (dst_buffer_2d->get_subregion_geometry_data()) 
    {
        mi::base::Handle<const nv::index::IRegular_heightfield_patch> patch(
            dst_buffer_2d->get_subregion_geometry_data()->get_interface<
            const nv::index::IRegular_heightfield_patch>());

        if (patch)
        {
            assert(patch->get_nb_grid_points() < 
                   static_cast<mi::Size>(std::numeric_limits<mi::Sint32>::max()));
            patch_array      = patch->get_elevation_values();
            patch_array_size = static_cast<mi::Sint32>(patch->get_nb_grid_points());
            patch_bbox = patch->get_bounding_rectangle();
            patch_size = patch_bbox.max - patch_bbox.min;
        }
        else{
            ERROR_LOG << "No heightfield patch found.";
        }
    }

    // Compute the size of the border around the patch by looking at difference between tile and
    // patch size
    const mi::math::Bbox<mi::Uint32, 2> border(
        mi::math::max(lower.x - patch_bbox.min.x, 0),
        mi::math::max(lower.y - patch_bbox.min.y, 0),
        mi::math::max(patch_bbox.max.x - upper.x, 0),
        mi::math::max(patch_bbox.max.y - upper.y, 0));

    // Write color values into the tile
    mi::Uint32 index = 0;
    for (mi::Uint32 j=0; j < texture_config.resolution.y; ++j)
    {
        for (mi::Uint32 i=0; i < texture_config.resolution.x; ++i)
        {
            mi::math::Color color(1.f); // Default to white
            if (patch_array)
            {
                mi::Sint32 pos = (j + border.min.y) * patch_size.x + (i + border.min.x);
                if (pos >= 0 && pos < patch_array_size)
                {
                    // Read the height value corresponding to the current pixel from the patch
                    mi::Float32 height = patch_array[pos];

                    if (height >= 0.f)
                    {
                        // Normalize the height value according to the given range
                        mi::Float32 t = (height - m_height_range.x) / (m_height_range.y - m_height_range.x);

                        // Apply elementwise linear interpolation of colors
                        color = mi::math::lerp(m_color0, m_color1, mi::math::clamp(t, 0.f, 1.f));
                    }
                }
            }

            texture_tile_rgba[index] = color;
            ++index;
        }
    }

    VERBOSE_LOG << "Distributed Compute: Height color tile for " << dst_buffer->get_subregion_bbox() << ":"
                << " resolution "     << texture_config.resolution
                << " covering "       << texture_config.covered_area;
}

void Distributed_compute_height_color_mapping_2d::set_colors(
    const mi::math::Color& color0,
    const mi::math::Color& color1)
{
    m_color0 = color0;
    m_color1 = color1;
}

void Distributed_compute_height_color_mapping_2d::set_enabled(bool enable)
{
    m_enabled = enable;
}

bool Distributed_compute_height_color_mapping_2d::get_enabled() const
{
    return m_enabled;
}

const char* Distributed_compute_height_color_mapping_2d::get_configuration() const
{
    std::ostringstream s;
    s << "generator=height_color\n"
      << "color0=" << m_color0.r << " " << m_color0.g << " " << m_color0.b << " " << m_color0.a << "\n"
      << "color1=" << m_color1.r << " " << m_color1.g << " " << m_color1.b << " " << m_color1.a << "\n";
    m_configuration = s.str();

    return m_configuration.c_str();
}

const char* Distributed_compute_height_color_mapping_2d::get_class_name() const
{
    return "Distributed_compute_height_color_mapping_2d";
}

mi::neuraylib::IElement* Distributed_compute_height_color_mapping_2d::copy() const
{
    Distributed_compute_height_color_mapping_2d* other = new Distributed_compute_height_color_mapping_2d();
    other->m_height_range = this->m_height_range;
    other->m_color0       = this->m_color0;
    other->m_color1       = this->m_color1;
    other->m_enabled      = this->m_enabled;
    return other;
}

void Distributed_compute_height_color_mapping_2d::serialize(mi::neuraylib::ISerializer* serializer) const
{
    serializer->write(m_height_range.begin(), m_height_range.SIZE);
    serializer->write(m_color0.begin(), m_color0.SIZE);
    serializer->write(m_color1.begin(), m_color1.SIZE);
    serializer->write(&m_enabled, 1);
}

void Distributed_compute_height_color_mapping_2d::deserialize(mi::neuraylib::IDeserializer* deserializer)
{
    deserializer->read(m_height_range.begin(), m_height_range.SIZE);
    deserializer->read(m_color0.begin(), m_color0.SIZE);
    deserializer->read(m_color1.begin(), m_color1.SIZE);
    deserializer->read(&m_enabled, 1);
}

// Distributed_compute_checkerboard_3d /////////////////////////////////////////////////////////////////////////////////
Distributed_compute_checkerboard_3d::Distributed_compute_checkerboard_3d()
  : m_checker_size(10u)
  , m_enabled(true)
{
    std::ostringstream s;
    s << "generator=checkerboard_volume\n"
      << "checker_size=" << m_checker_size.x << " " << m_checker_size.y << " " << m_checker_size.z << "\n";
    m_configuration = s.str();
}

Distributed_compute_checkerboard_3d::Distributed_compute_checkerboard_3d(
        const mi::math::Vector<mi::Uint32, 3>& checker_size)
  : m_checker_size(checker_size)
  , m_enabled(true)
{
    std::ostringstream s;
    s << "generator=checkerboard_volume\n"
      << "checker_size=" << m_checker_size.x << " " << m_checker_size.y << " " << m_checker_size.z << "\n";
    m_configuration = s.str();
}

void Distributed_compute_checkerboard_3d::launch_compute(
        mi::neuraylib::IDice_transaction*                   dice_transaction,
        nv::index::IDistributed_compute_destination_buffer* dst_buffer) const
{
    using nv::index::IDistributed_compute_destination_buffer_3d_volume;
    using nv::index::IDistributed_compute_destination_buffer_3d_volume_uint8;

    // retrieve 2d texture buffer destination buffer interface
    const mi::base::Handle<IDistributed_compute_destination_buffer_3d_volume> dst_buffer_3d(
        dst_buffer->get_interface<IDistributed_compute_destination_buffer_3d_volume>());

    if (!dst_buffer_3d)
    {
        ERROR_LOG << "Distributed compute technique (checkerboard-volume): "
                  << "Unable to retrieve valid 3d volume destination buffer interface.";
        return;
    }

    const mi::base::Handle<IDistributed_compute_destination_buffer_3d_volume_uint8> dst_buffer_3d_uint8(
        dst_buffer_3d->get_interface<IDistributed_compute_destination_buffer_3d_volume_uint8>());

    if (!dst_buffer_3d_uint8)
    {
        ERROR_LOG << "Distributed compute technique (checkerboard-volume): "
                  << "Unable to retrieve valid uint8 3d volume destination buffer interface.";
        return;
    }

    const mi::math::Bbox<mi::Sint32, 3> brick_bbox_clipped = dst_buffer_3d->get_clipped_volume_brick_bbox();

    const mi::Size brick_width  = brick_bbox_clipped.max.x - brick_bbox_clipped.min.x;
    const mi::Size brick_height = brick_bbox_clipped.max.y - brick_bbox_clipped.min.y;
    const mi::Size brick_depth  = brick_bbox_clipped.max.z - brick_bbox_clipped.min.z;

    const mi::Size brick_size   = brick_width * brick_height * brick_depth * sizeof(mi::Uint8);

    mi::Uint8*     brick_array  = new mi::Uint8[brick_size];

    mi::Size i = 0ull;
    for (mi::Sint32 x = brick_bbox_clipped.min.x; x < brick_bbox_clipped.max.x; ++x)
    {
        for (mi::Sint32 y = brick_bbox_clipped.min.y; y < brick_bbox_clipped.max.y; ++y)
        {
            for (mi::Sint32 z = brick_bbox_clipped.min.z; z < brick_bbox_clipped.max.z; ++z)
            {
                const mi::Sint32 cx = x / static_cast<mi::Sint32>(m_checker_size.x);
                const mi::Sint32 cy = y / static_cast<mi::Sint32>(m_checker_size.y);
                const mi::Sint32 cz = z / static_cast<mi::Sint32>(m_checker_size.z);

                brick_array[i++] = (cx + cy + cz) % 2 == 0 ? 0u : 255u;
            }
        }
    }

    dst_buffer_3d_uint8->write(brick_bbox_clipped, brick_array, -1);

    delete [] brick_array;

    VERBOSE_LOG << "Distributed Compute: Checkerboard (volume) "
                << "tile for " << dst_buffer->get_subregion_bbox() << ": "
                << "brick bbox "           << brick_bbox_clipped;

}

const char* Distributed_compute_checkerboard_3d::get_configuration() const
{
    return m_configuration.c_str();
}

void Distributed_compute_checkerboard_3d::set_enabled(bool enable)
{
    m_enabled = enable;
}

bool Distributed_compute_checkerboard_3d::get_enabled() const
{
    return m_enabled;
}

mi::neuraylib::IElement* Distributed_compute_checkerboard_3d::copy() const
{
    Distributed_compute_checkerboard_3d* other = new Distributed_compute_checkerboard_3d();

    other->m_checker_size  = this->m_checker_size;
    other->m_configuration = this->m_configuration;
    other->m_enabled       = this->m_enabled;

    return other;
}

const char* Distributed_compute_checkerboard_3d::get_class_name() const
{
    return "Distributed_compute_checkerboard_3d";
}

void Distributed_compute_checkerboard_3d::serialize(mi::neuraylib::ISerializer* serializer) const
{
    serializer->write(&m_checker_size.x, 3);
    serializer->write(&m_enabled, 1);

    mi::Uint32 nb_elements = mi::Uint32(m_configuration.size());
    serializer->write(&nb_elements, 1);
    serializer->write(reinterpret_cast<const mi::Uint8*>(m_configuration.c_str()), nb_elements);

}

void Distributed_compute_checkerboard_3d::deserialize(mi::neuraylib::IDeserializer* deserializer)
{
    deserializer->read(&m_checker_size.x, 3);
    deserializer->read(&m_enabled, 1);

    mi::Uint32 nb_elements;
    deserializer->read(&nb_elements, 1);
    m_configuration.resize(nb_elements);
    deserializer->read(reinterpret_cast<mi::Uint8*>(&m_configuration[0]), nb_elements);

}

void Distributed_compute_checkerboard_3d::get_references(
    mi::neuraylib::ITag_set* result) const
{
}

} // namespace index_common
} // namespace nv
