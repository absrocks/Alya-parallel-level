/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Piecewise linear rgba colormap implementation

#ifndef GEOSPATIAL_BIN_COMMON_TRANSFERFUNCTION_COLORMAP_H
#define GEOSPATIAL_BIN_COMMON_TRANSFERFUNCTION_COLORMAP_H

#include <nv/index/icolormap.h>
#include <nv/index/iattribute.h>

#include <mi/math/color.h>

#include "spline.h"

#include <vector>
#include <sstream>

namespace nv {
namespace index_common {

/// A piecewise linear RGBA colormap implementation
///
/// - Initialize with global voxel value domain
/// - Insert/Remove the anchor points. An anchor point is a pair of {voxel_value, color}.
/// - Set/Get the anchor points voxel value and color. 
/// - Change the domain and domain resolution.
/// 
class Transferfunction_colormap : 
        public mi::neuraylib::Element<0xead24ac7,0xc0f4,0x441f,0xab,0xc9,0xee,0xa2,0x00,0x07,0x31,0x3b,
                                      nv::index::IColormap>
{
public:
    /// default constructor
    Transferfunction_colormap();

    /// destructor
    virtual ~Transferfunction_colormap();

    /// typedef of voxel_value, color pair
    struct Voxel_value_color_pair
    {
        /// default constructor
        Voxel_value_color_pair()
            :
            m_voxel_value(0.0f),
            m_color_value(0.0f, 0.0f, 0.0f, 1.0f)
        {
            // empty
        }
        /// constructor
        /// \param[in] voxel_value voxel value
        /// \param[in] color_value color value at the voxel value
        Voxel_value_color_pair(mi::Float32 voxel_value,
                               const mi::math::Color& color_value)
            :        
            m_voxel_value(voxel_value),
            m_color_value(color_value)
        {
            // empty
        }

        /// Get string representation of this object.
        /// \return string representation
        std::string get_string() const 
        {
            std::stringstream sstr;
            sstr << m_voxel_value << ", ";
            for (mi::Size i = 0; i < 4; ++i)
            {
                sstr << m_color_value[i] << " ";
            }
            return sstr.str();
        }

        /// voxel value
        mi::Float32     m_voxel_value;
        /// color value
        mi::math::Color m_color_value;
    };

    /// Voxel_value_color_pair comparison functor by voxel_value
    struct Voxel_value_comp 
    {
        /// comparison operator
        bool operator()(const Voxel_value_color_pair& lhs, 
                        const Voxel_value_color_pair& rhs)
        {
            return lhs.m_voxel_value < rhs.m_voxel_value;
        }
    };

    /// color space name 
    enum Color_space_e 
    {
        /// Color space RGB
        COLORSPACE_RGB,
        /// Color space HSV
        COLORSPACE_HSV,
        /// sentinel
        NB_COLORSPACE
    };

public:
    // Transferfunction colormap interface

    /// Set color space (default COLORSPACE_RGB)
    ///
    /// Note: this applies when update_colormap() called, since this
    /// colormap output is always in RBG space.
    ///
    /// \param[in] color_space color space.
    void set_color_space(Color_space_e color_space);
    
    /// Get current color space
    ///
    /// \return current color space.
    Color_space_e get_color_space() const;

    /// initialize the colormap. Remove all the anchor point and set
    /// two anchor points {min, max}.
    /// 
    /// \param[in] global_voxel_range     global [min, max] voxel value domain for the colormap
    /// \param[in] color_at_min           global color value at global minimal voxel value
    /// \param[in] color_at_max           global color value at global maximum voxel value
    /// \return true when succeeded
    bool initialize_colormap(const mi::math::Vector<mi::Float32, 2>& global_voxel_range, 
                             const mi::math::Color& color_at_min,
                             const mi::math::Color& color_at_max);

    /// Insert anchor point at the voxel_val
    ///
    /// Note: You can only insert an anchor at (0, size-1). This means
    /// not at 0, not at size-1.
    ///
    /// \param[in] anchor_point_idx anchor_point index to be inserted
    /// \param[in] voxel_val anchor point voxel value
    /// \param[in] color_val anchor point color value
    /// \return true when insert succeeded.
    bool insert_anchor_point(mi::Size anchor_point_idx, 
                             mi::Float32 voxel_val, 
                             const mi::math::Color& color_val);

    /// remove anchor point
    ///
    /// \param[in] anchor_point_idx anchor_point index to be removed
    /// \return true when succeeded
    bool remove_anchor_point(mi::Size anchor_point_idx);

    /// Get number anchor points
    ///
    /// \return number of anchor points
    mi::Size get_nb_anchor_point() const;

    /// Set anchor point color
    ///
    /// \param[in] anchor_point_idx anchor point index
    /// \param[in] color_val        color value of the anchor point
    /// \return true when succeeded
    bool set_anchor_color(mi::Size anchor_point_idx, const mi::math::Color& color_val);

    /// Set anchor point color
    ///
    /// \param[in]  anchor_point_idx anchor point index
    /// \param[out] color            color at anchor_point_idx
    /// \return true when succeeded
    bool get_anchor_color(mi::Size anchor_point_idx,
                          mi::math::Color& color) const;

    /// Set anchor point voxel value
    ///
    /// Note: Can not set the value over the next anchor point value,
    /// less value the previous anchor point value.
    ///
    /// \param[in] anchor_point_idx anchor point index
    /// \param[in] voxel_val  anchor point voxel value
    /// \return true when succeeded
    bool set_anchor_voxel_value(mi::Size anchor_point_idx, mi::Float32 voxel_val);
    
    /// Get anchor point voxel value
    ///
    /// \param[in]  anchor_point_idx anchor point index
    /// \param[out] voxel_val  anchor point voxel value
    /// \return true when succeeded
    bool get_anchor_voxel_value(mi::Size anchor_point_idx, mi::Float32& voxel_val) const;

    /// Set voxel domain resolution. This is effective when
    /// update_colormap() is called.
    ///
    /// \param[in] domain_resolution voxel domain resolution
    /// \return true when succeeded
    bool set_voxel_domain_resolution(mi::Size domain_resolution);

    /// Get voxel domain resolution
    ///
    /// \return domain resolution
    mi::Size get_voxel_domain_resolution() const;

    /// Set interpolation mode
    ///
    /// \param[in] is_liner true when linear interpolation mode
    void set_linear_interpolation_mode(bool is_liner);

    /// Get interpolation mode
    ///
    /// \return true when in the interpolation mode
    bool is_linear_interpolation_mode() const;

    /// Update the internal colormap entry
    void update_colormap();

public:
    // IColormap implementation
    virtual mi::math::Color_struct get_color(mi::Uint32 idx) const;
    virtual void set_color(
        mi::Uint32                    idx,
        const mi::math::Color_struct& color);
    virtual void set_colormap(
        const mi::math::Color_struct* colormap_values,
        mi::Uint32                    nb_entries);
    virtual const mi::math::Color_struct* get_colormap() const;
    virtual mi::Uint32 get_number_of_entries() const;
    virtual void get_non_transparent_range(
        mi::Uint32& min,
        mi::Uint32& max) const;

    virtual void set_domain(
        mi::Float32 value_first,
        mi::Float32 value_last);

    virtual void get_domain(
        mi::Float32& value_first,
        mi::Float32& value_last) const;

    virtual void set_domain_boundary_mode(Domain_boundary_mode mode);
    virtual Domain_boundary_mode get_domain_boundary_mode() const;

public:
    // implements IAttribute
    virtual mi::base::Uuid get_attribute_class() const;
    virtual bool multiple_active_instances() const;

public:
    // implements IScene_element
    virtual void set_enabled(bool enable);
    virtual bool get_enabled() const;
    virtual void set_meta_data(mi::neuraylib::Tag_struct tag);
    virtual mi::neuraylib::Tag_struct get_meta_data() const;

private:
    /// get current anchor_point_idx left and right voxel value
    ///
    /// \param[in]  idx        query index
    /// \param[out] left_right left (idx - 1), right (idx) value 
    /// \return false when no left_right value
    bool get_current_idx_left_right_voxel_value(
        mi::Size idx,
        mi::math::Vector<mi::Float32, 2>& left_right);

    /// check anchors are sorted or not
    ///
    /// \return true when the colormap entries are sorted by voxel value
    bool is_sorted_anchor_vec() const;

    /// Generate linear interpolated colormap
    void generate_linear_interpolated_colormap();

    /// Linear interpolation of the color between anchor points
    ///
    /// \param[in] voxel_val voxel value
    /// \return interpolated color
    mi::math::Color linear_interpolate_color(mi::Float32 voxel_val);

    /// Generate spline interpolated colormap
    void generate_spline_interpolated_colormap();

    /// Spline interpolation of the alpha between anchor points
    ///
    /// \param[in] voxel_val voxel value
    /// \return interpolated alpha
    mi::Float32 spline_interpolate_color(
        nv::index_common::Cubic_hermite_spline& chs, 
        mi::Float32 voxel_val);
    
    /// convert current each colormap value from HSV to RGB
    void convert_current_colormap_hsv_to_rgb();

public:    
    // IElement implementation
    virtual mi::neuraylib::IElement* copy() const;
    virtual const char* get_class_name() const
    {
        return "Transferfunction_colormap"; 
    }
    virtual void serialize(mi::neuraylib::ISerializer* serializer) const;
    virtual void deserialize(mi::neuraylib::IDeserializer* deserializer);

private:
    /// IScene_element: enable
    bool m_is_enabled;

    /// Colormap values: index -> Color
    std::vector<Voxel_value_color_pair> m_anchor_vec;
    /// global voxel value domain
    mi::math::Vector<mi::Float32, 2> m_global_voxel_value_domain;
    /// local  voxel value domain
    mi::math::Vector<mi::Float32, 2> m_local_voxel_value_domain;
    /// domain resolution (size)
    mi::Size m_domain_resolution;
    /// is linear interpolation mode
    bool m_is_linear_interpolation;
    /// current color space
    Color_space_e m_color_space;

    /// Colormap domain settings
    Domain_boundary_mode             m_domain_boundary_mode;

    /// discretized colormap (consider they are in the rgb space): not serialized
    std::vector<mi::math::Color_struct> m_color_vec;
};

} // namespace index_common
} // namespace nv

#endif // #ifndef GEOSPATIAL_BIN_COMMON_TRANSFERFUNCTION_COLORMAP_H
