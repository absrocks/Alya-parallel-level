/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Example texturing techniques that can be assigned to planes or heightfields.

#ifndef NVIDIA_INDEX_BIN_COMMON_DISTRIBUTED_COMPUTE_TECHNIQUES_H
#define NVIDIA_INDEX_BIN_COMMON_DISTRIBUTED_COMPUTE_TECHNIQUES_H

#include <nv/index/idistributed_compute_technique.h>
#include <nv/index/idistributed_compute_destination_buffer.h>

#include <mi/math/vector.h>

#include <string>

namespace nv {
namespace index_common {

/// Interface class for user-defined generation of tiles of a 2D texture on demand.
///
/// An implementation of the class can be assigned to a plane or height-field scene element.
/// When rendering the certain 2D area of the plane (defined by the renderer)
/// the texture creation technique is called and supposed to return a color
/// texture (a pointer to a color buffer representing a 2D array). The
/// renderer then takes care to map the texture onto the defined 2D plane area.
///
/// A typical use-case for the generation technique is to visualize data onto
/// the plane that result from a complex, possibly distributed computing process.
///
class Distributed_compute_checkerboard_2d :
    public nv::index::Distributed_compute_technique<0x57aeba48,0x3bc0,0x462e,0x9a,0xe4,0x6f,0xcd,0x1,0xc9,0x19,0xc6>
{
public:
    /// Empty default constructor for serialization
    Distributed_compute_checkerboard_2d();

    /// Constructor
    ///
    /// \param[in]  extent  The extent of the associated plane scene element.
    /// \param[in]  fmt     The buffer format
    Distributed_compute_checkerboard_2d(
        const mi::math::Vector_struct<mi::Float32, 2>&                               extent,
        nv::index::IDistributed_compute_destination_buffer_2d_texture::Buffer_format fmt);

    virtual void launch_compute(
        mi::neuraylib::IDice_transaction*                   dice_transaction,
        nv::index::IDistributed_compute_destination_buffer* dst_buffer) const;

    /// Returns configuration settings that may be used by the library for the session export
    /// mechanism provided by ISession::export_session().
    virtual const char* get_configuration() const;

    virtual void set_enabled(bool enable);
    virtual bool get_enabled() const;

    void set_color_0(const mi::math::Color& color_0);
    void set_color_1(const mi::math::Color& color_1);

    /// -------------------------------------------------------------------------------------------

    /// The copy needs to create and return a newly allocated object of the same
    /// class and has to copy all the fields which should be the same in the copy.
    ///
    /// \return     The new copy of the database element.
    ///
    virtual mi::neuraylib::IElement* copy() const;

    /// Get a human readable representation of the class name.
    ///
    virtual const char* get_class_name() const;

    /// Serialize the class to the given serializer.
    ///
    /// \param serializer       Write to this serializer.
    ///
    virtual void serialize(mi::neuraylib::ISerializer* serializer) const;

    /// Deserialize the class from the given deserializer.
    ///
    /// \param deserializer     Read from this deserializer.
    ///
    virtual void deserialize(mi::neuraylib::IDeserializer* deserializer);

    virtual void get_references(
        mi::neuraylib::ITag_set* result) const;

private:
    void update_configuration();

    /// Generates a texture tile in ONE_DIMENSIONAL_ARRAY mode.
    void generate_texture_tile_1D(
        mi::neuraylib::IDice_transaction*                                    dice_transaction,
        const nv::index::IDistributed_compute_destination_buffer_2d_texture* dst_buffer,
        const mi::math::Vector_struct<mi::Float32, 3>*                       intersections,
        const mi::Uint32                                                     nb_intersections) const;

    /// Generates a texture tile in TWO_DIMENSIONAL_ARRAY mode.
    void generate_texture_tile_2D(
        mi::neuraylib::IDice_transaction*                                    dice_transaction,
        const nv::index::IDistributed_compute_destination_buffer_2d_texture* dst_buffer,
        const mi::math::Vector<mi::Sint32, 2>&                               ll_corner,
        const mi::math::Vector<mi::Sint32, 2>&                               ur_corner) const;

    mi::math::Vector<mi::Float32, 2> m_extent;

    nv::index::IDistributed_compute_destination_buffer_2d_texture::Buffer_format m_buffer_format;

    mi::Float32                      m_scale;
    mi::math::Color                  m_color_0;
    mi::math::Color                  m_color_1;

    std::string                      m_configuration;
    bool                             m_enabled;
};

/// Another texture tile generator example which renders the Mandelbrot set.
class Distributed_compute_mandelbrot_2d :
    public nv::index::Distributed_compute_technique<0xab79b02e,0xcbdd,0x4319,0xbe,0x43,0x95,0xff,0x5,0x69,0xdf,0x23>
{
public:
    /// Empty default constructor for serialization
    Distributed_compute_mandelbrot_2d();

    /// \param[in]  extent  The extent of the associated plane scene element.
    Distributed_compute_mandelbrot_2d(const mi::math::Vector_struct<mi::Float32, 2>& extent);

    virtual void launch_compute(
        mi::neuraylib::IDice_transaction*                   dice_transaction,
        nv::index::IDistributed_compute_destination_buffer* dst_buffer) const;

    /// Returns configuration settings that may be used by the library for the session export
    /// mechanism provided by ISession::export_session().
    virtual const char* get_configuration() const;

    virtual void set_enabled(bool enable);
    virtual bool get_enabled() const;

    /// -----------------------------------------------------------
    virtual mi::neuraylib::IElement* copy() const;
    virtual const char* get_class_name() const;
    virtual void serialize(mi::neuraylib::ISerializer* serializer) const;
    virtual void deserialize(mi::neuraylib::IDeserializer* deserializer);

private:
    /// Generates a texture tile in ONE_DIMENSIONAL_ARRAY mode.
    void generate_texture_tile_1D(
        mi::neuraylib::IDice_transaction*                                    dice_transaction,
        const nv::index::IDistributed_compute_destination_buffer_2d_texture* dst_buffer,
        const mi::math::Vector_struct<mi::Float32, 3>*                       intersections,
        const mi::Uint32                                                     nb_intersections) const;

    /// Generates a texture tile in TWO_DIMENSIONAL_ARRAY mode.
    void generate_texture_tile_2D(
        mi::neuraylib::IDice_transaction*                                    dice_transaction,
        const nv::index::IDistributed_compute_destination_buffer_2d_texture* dst_buffer,
        const mi::math::Vector<mi::Sint32, 2>&                               ll_corner,
        const mi::math::Vector<mi::Sint32, 2>&                               ur_corner) const;

    mi::math::Color mandelbrot(
        const mi::math::Vector<mi::Float32, 2>& tc) const;

    mi::math::Vector<mi::Float32, 2> m_extent;
    std::string                      m_configuration;
    bool                             m_enabled;
};

/// Another texture tile generator that obtains the texture from a file
class Distributed_compute_bitmap_mapping_2d :
    public nv::index::Distributed_compute_technique<0xe577e314,0x2a03,0x4cf1,0x8f,0x6,0x1d,0x12,0x6d,0x44,0xdc,0xc8>
{
public:
    /// Empty default constructor for serialization
    Distributed_compute_bitmap_mapping_2d();

    /// \param[in]  bm_filemane The bitmap texture file name.
    /// \param[in]  use_cache   Enable or disable using a cache.
    /// \param[in]  extent      The extent of the associated plane scene element.
    Distributed_compute_bitmap_mapping_2d(
        const std::string&                              bm_filemane,
        bool                                            use_cache,
        const mi::math::Vector_struct<mi::Float32, 2>&  extent);
    
    virtual void launch_compute(
        mi::neuraylib::IDice_transaction*                   dice_transaction,
        nv::index::IDistributed_compute_destination_buffer* dst_buffer) const;

    virtual void set_enabled(bool enable);
    virtual bool get_enabled() const;

    /// -----------------------------------------------------------
    virtual mi::neuraylib::IElement* copy() const;
    virtual const char* get_class_name() const;
    virtual void serialize(mi::neuraylib::ISerializer* serializer) const;
    virtual void deserialize(mi::neuraylib::IDeserializer* deserializer);

private:
    /// Generates a texture tile in ONE_DIMENSIONAL_ARRAY mode.
    void generate_texture_tile_1D(
        mi::neuraylib::IDice_transaction*                                    dice_transaction,
        const nv::index::IDistributed_compute_destination_buffer_2d_texture* dst_buffer,
        const mi::math::Vector_struct<mi::Float32, 3>*                       intersections,
        const mi::Uint32                                                     nb_intersections) const;

    std::string                      m_bm_filemane;
    mi::math::Vector<mi::Float32, 2> m_extent;
    bool                             m_enabled;
    bool                             m_use_cache;
};

/// Example texturing technique for heightfields that generates colors based on height values.
class Distributed_compute_height_color_mapping_2d :
    public nv::index::Distributed_compute_technique<0x8cab541d,0xc61a,0x42f5,0xb6,0xe6,0xcf,0xdf,0x3b,0x2,0x31,0x2c>
{
public:
    /// Empty default constructor for serialization
    Distributed_compute_height_color_mapping_2d();

    Distributed_compute_height_color_mapping_2d(
        const mi::math::Vector<mi::Float32, 2>& height_range);


    virtual void launch_compute(
        mi::neuraylib::IDice_transaction*                   dice_transaction,
        nv::index::IDistributed_compute_destination_buffer* dst_buffer) const;

    /// Sets the two colors between which will be interpolated.
    ///
    /// \param[in] color0 Color representing minimum height values
    /// \param[in] color1 Color representing maximum height values
    void set_colors(
        const mi::math::Color& color0,
        const mi::math::Color& color1);

    virtual void set_enabled(bool enable);
    virtual bool get_enabled() const;

    virtual const char* get_configuration() const;

    virtual const char* get_class_name() const;
    virtual mi::neuraylib::IElement* copy() const;
    virtual void serialize(mi::neuraylib::ISerializer *serializer) const;
    virtual void deserialize(mi::neuraylib::IDeserializer* deserializer);

private:
    mi::math::Vector<mi::Float32, 2> m_height_range;
    mi::math::Color                  m_color0;
    mi::math::Color                  m_color1;
    bool                             m_enabled;
    mutable std::string              m_configuration;
};

/// Example volume compute technique
class Distributed_compute_checkerboard_3d :
    public nv::index::Distributed_compute_technique<0x3d8ee43e,0xcaac,0x4590,0x85,0xca,0xd,0xc1,0x2,0xd8,0x65,0x2f>
{
public:
    /// Empty default constructor for serialization
    Distributed_compute_checkerboard_3d();
    Distributed_compute_checkerboard_3d(const mi::math::Vector<mi::Uint32, 3>& checker_size);

    virtual void launch_compute(
        mi::neuraylib::IDice_transaction*                   dice_transaction,
        nv::index::IDistributed_compute_destination_buffer* dst_buffer) const;

    /// Returns configuration settings that may be used by the library for the session export
    /// mechanism provided by ISession::export_session().
    virtual const char* get_configuration() const;

    virtual void set_enabled(bool enable);
    virtual bool get_enabled() const;

    /// -------------------------------------------------------------------------------------------

    /// The copy needs to create and return a newly allocated object of the same
    /// class and has to copy all the fields which should be the same in the copy.
    ///
    /// \return     The new copy of the database element.
    ///
    virtual mi::neuraylib::IElement* copy() const;

    /// Get a human readable representation of the class name.
    ///
    virtual const char* get_class_name() const;

    /// Serialize the class to the given serializer.
    ///
    /// \param serializer       Write to this serializer.
    ///
    virtual void serialize(mi::neuraylib::ISerializer* serializer) const;

    /// Deserialize the class from the given deserializer.
    ///
    /// \param deserializer     Read from this deserializer.
    ///
    virtual void deserialize(mi::neuraylib::IDeserializer* deserializer);

    virtual void get_references(
        mi::neuraylib::ITag_set* result) const;

private:
    mi::math::Vector<mi::Uint32, 3>  m_checker_size;
    std::string                      m_configuration;
    bool                             m_enabled;

};

} // namespace index_common
} // namespace nv

#endif // NVIDIA_INDEX_BIN_COMMON_DISTRIBUTED_COMPUTE_TECHNIQUES_H
