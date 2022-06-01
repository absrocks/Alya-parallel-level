/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "transferfunction_colormap.h"

#include <mi/dice.h>

#include "forwarding_logger.h"
#include "color_space.h"

#include <cassert>

namespace nv {
namespace index_common {
//----------------------------------------------------------------------
//----------------------------------------------------------------------
Transferfunction_colormap::Transferfunction_colormap()
    :
    m_is_enabled(true),
    m_domain_resolution(255),
    m_is_linear_interpolation(true),
    m_color_space(COLORSPACE_RGB),
    m_domain_boundary_mode(CLAMP_TO_EDGE)
{
    // initialize_colormap(mi::math::Vector<mi::Float32, 2>(0.0f, 255.0f),
    //                     mi::math::Color(0.0f, 0.0f, 0.0f, 1.0f),
    //                     mi::math::Color(1.0f, 1.0f, 1.0f, 1.0f));
}

//----------------------------------------------------------------------
Transferfunction_colormap::~Transferfunction_colormap()
{
    m_anchor_vec.clear();
    m_color_vec.clear();
}

//----------------------------------------------------------------------
void Transferfunction_colormap::set_color_space(Transferfunction_colormap::Color_space_e color_space)
{
    if (static_cast<mi::Uint32>(color_space) >= static_cast<mi::Uint32>(NB_COLORSPACE))
    {
        ERROR_LOG << "Transferfunction_colormap::set_color_space: illegal color_space: "
                  << color_space << ". Ignored.";
        return;
    }

    m_color_space = color_space;
}

//----------------------------------------------------------------------
Transferfunction_colormap::Color_space_e Transferfunction_colormap::get_color_space() const
{
    return m_color_space;
}

//----------------------------------------------------------------------
bool Transferfunction_colormap::initialize_colormap(const mi::math::Vector<mi::Float32, 2>& global_voxel_domain,
                                                    const mi::math::Color& color_at_min,
                                                    const mi::math::Color& color_at_max)
{
    m_anchor_vec.clear();
    m_color_vec.clear();

    if (global_voxel_domain[0] >= global_voxel_domain[1])
    {
        ERROR_LOG << "Illegal voxel value range. (min >= max). Cannot initialize the colormap.";
        return false;
    }

    m_global_voxel_value_domain = global_voxel_domain;
    m_local_voxel_value_domain  = global_voxel_domain;

    m_anchor_vec.push_back(Voxel_value_color_pair(m_global_voxel_value_domain[0], color_at_min));
    m_anchor_vec.push_back(Voxel_value_color_pair(m_global_voxel_value_domain[1], color_at_max));

    // update_colormap();

    return true;
}

//----------------------------------------------------------------------
bool Transferfunction_colormap::insert_anchor_point(mi::Size    anchor_point_idx,
                                                    mi::Float32 voxel_val,
                                                    const mi::math::Color& color_val)
{
    assert(m_anchor_vec.size() >= 2);
    const std::string mn = "Transferfunction_colormap::insert_anchor_point: ";

    if (anchor_point_idx >= m_anchor_vec.size())
    {
        ERROR_LOG << mn << "Index is out of range, " << "anchor_point_idx = " << anchor_point_idx
                  << ", color_map_entries.size = " << m_anchor_vec.size();
        return false;
    }

    if (anchor_point_idx == 0)
    {
        ERROR_LOG << mn << "Cannot insert at min.";
        return false;
    }

    if (anchor_point_idx == m_anchor_vec.size())
    {
        ERROR_LOG << mn << "Cannot insert at max.";
        return false;
    }

    mi::math::Vector<mi::Float32, 2> left_right(0.0f, 0.0f);
    if (!get_current_idx_left_right_voxel_value(anchor_point_idx, left_right))
    {
        ERROR_LOG << mn << "Cannot get left and right value.";
        return false;
    }

    // exclude the same value for simplicity
    if ((left_right[0] >= voxel_val) || (voxel_val >= left_right[1]))
    {
        ERROR_LOG << mn << "voxel value (" << voxel_val << ") is out of range "
                  << left_right << ".";
        return false;
    }

    m_anchor_vec.insert(m_anchor_vec.begin() + anchor_point_idx,
                        Voxel_value_color_pair(voxel_val, color_val));

    return true;
}

//----------------------------------------------------------------------
bool Transferfunction_colormap::remove_anchor_point(mi::Size anchor_point_idx)
{
    if (anchor_point_idx == 0)
    {
        ERROR_LOG << "Transferfunction_colormap::remove_anchor_point: cannot remove the min anchor: "
                  << anchor_point_idx;
        return false;
    }

    if ((anchor_point_idx + 1) == m_anchor_vec.size())
    {
        ERROR_LOG << "Transferfunction_colormap::remove_anchor_point: cannot remove the max anchor: "
                  << anchor_point_idx;
        return false;
    }

    if ((anchor_point_idx + 1)> m_anchor_vec.size())
    {
        ERROR_LOG << "Transferfunction_colormap::remove_anchor_point: anchor index is out of range: "
                  << anchor_point_idx
                  << ", entry size: " << m_anchor_vec.size();
        return false;
    }

    m_anchor_vec.erase(m_anchor_vec.begin() + anchor_point_idx);

    return true;
}

//----------------------------------------------------------------------
mi::Size Transferfunction_colormap::get_nb_anchor_point() const
{
    return m_anchor_vec.size();
}

//----------------------------------------------------------------------
bool Transferfunction_colormap::set_anchor_color(mi::Size anchor_point_idx,
                                                 const mi::math::Color& color_val)
{
    if (anchor_point_idx >= m_anchor_vec.size())
    {
        ERROR_LOG << "Transferfunction_colormap::set_anchor_color: Index is out of range: "
                  << anchor_point_idx
                  << ", entry size: " << m_anchor_vec.size();
        return false;
    }

    m_anchor_vec[anchor_point_idx].m_color_value = color_val;

    return true;
}

//----------------------------------------------------------------------
bool Transferfunction_colormap::get_anchor_color(mi::Size anchor_point_idx,
                                                 mi::math::Color& color_val) const
{
    if (anchor_point_idx >= m_anchor_vec.size())
    {
        ERROR_LOG << "Transferfunction_colormap::set_anchor_color: Index is out of range: "
                  << anchor_point_idx
                  << ", entry size: " << m_anchor_vec.size();
        return false;
    }

    color_val = m_anchor_vec[anchor_point_idx].m_color_value;
    return true;
}

//----------------------------------------------------------------------
bool Transferfunction_colormap::set_anchor_voxel_value(mi::Size anchor_point_idx, mi::Float32 voxel_val)
{
    const std::string mn = "Transferfunction_colormap::set_anchor_voxel_value: ";
    if (anchor_point_idx == 0)
    {
        ERROR_LOG << mn << "Cannot move global min anchor: " << anchor_point_idx;
        return false;
    }

    if ((anchor_point_idx + 1) == m_anchor_vec.size())
    {
        ERROR_LOG << mn << "Cannot move global max anchor: " << anchor_point_idx;
        return false;
    }

    if (anchor_point_idx > m_anchor_vec.size())
    {
        ERROR_LOG << mn << "Index is out of range: " << anchor_point_idx
                  << ", entry size: " << m_anchor_vec.size();
        return false;
    }

    mi::math::Vector<mi::Float32, 2> left_range(0.0f, 0.0f);
    mi::math::Vector<mi::Float32, 2> right_range(0.0f, 0.0f);
    if (!get_current_idx_left_right_voxel_value(anchor_point_idx, left_range))
    {
        ERROR_LOG << mn << "Cannot get left range of the anchor: " << anchor_point_idx;
        return false;
    }

    if (!get_current_idx_left_right_voxel_value(anchor_point_idx + 1, right_range))
    {
        ERROR_LOG << mn << "Cannot get right range of the anchor: " << anchor_point_idx;
        return false;
    }

    const mi::Float32 left_min  = left_range[0];
    const mi::Float32 right_max = right_range[1];

    if ((left_min >= voxel_val) || (voxel_val >= right_max))
    {
        ERROR_LOG << mn << "voxel value (" << voxel_val << ") is out of range ["
                  << left_min << " " << right_max << "].";
        return false;
    }

    m_anchor_vec[anchor_point_idx].m_voxel_value = voxel_val;

    return true;
}

//----------------------------------------------------------------------
bool Transferfunction_colormap::get_anchor_voxel_value(mi::Size anchor_point_idx,
                                                       mi::Float32& voxel_val) const
{
    const std::string mn = "Transferfunction_colormap::get_anchor_voxel_value: ";
    if (anchor_point_idx >= m_anchor_vec.size())
    {
        ERROR_LOG << mn << "Index is out of range: " << anchor_point_idx
                  << ", entry size: " << m_anchor_vec.size();
        return false;
    }

    voxel_val = m_anchor_vec[anchor_point_idx].m_voxel_value;

    return true;
}

//----------------------------------------------------------------------
bool Transferfunction_colormap::set_voxel_domain_resolution(mi::Size domain_resolution)
{
    const std::string mn = "Transferfunction_colormap::set_voxel_domain_resolution: ";
    if (domain_resolution <= 1)
    {
        ERROR_LOG << mn << "illegal domain resolution. Should be > 1.";
        return false;
    }

    if (domain_resolution >= 65536)
    {
        ERROR_LOG << mn << "illegal domain resolution. Should be < 65536.";
        return false;
    }

    m_domain_resolution = domain_resolution;

    return true;
}

//----------------------------------------------------------------------
mi::Size Transferfunction_colormap::get_voxel_domain_resolution() const
{
    return m_domain_resolution;
}

//----------------------------------------------------------------------
void Transferfunction_colormap::set_linear_interpolation_mode(bool is_liner)
{
    m_is_linear_interpolation = is_liner;
}

//----------------------------------------------------------------------
bool Transferfunction_colormap::is_linear_interpolation_mode() const
{
    return m_is_linear_interpolation;
}

//----------------------------------------------------------------------
void Transferfunction_colormap::update_colormap()
{
    assert(m_domain_resolution > 2);
    assert(is_sorted_anchor_vec());

    if (m_is_linear_interpolation)
    {
        generate_linear_interpolated_colormap();
    }
    else
    {
        generate_spline_interpolated_colormap();
    }

    if (get_color_space() == COLORSPACE_RGB)
    {
        // no conversion
    }
    else if (get_color_space() == COLORSPACE_HSV)
    {
        convert_current_colormap_hsv_to_rgb();
    }
    else
    {
        assert(false);          // should not happen
    }
}

//----------------------------------------------------------------------
mi::math::Color_struct Transferfunction_colormap::get_color(mi::Uint32 idx) const
{
    const std::string mn = "Transferfunction_colormap::get_color: ";
    if (idx >= m_color_vec.size())
    {
        ERROR_LOG << mn << "idx is out of range.";
        return mi::math::Color(0.0f, 0.0f, 0.0f, 1.0f);
    }

    return m_color_vec[idx];
}

//----------------------------------------------------------------------
void Transferfunction_colormap::set_color(
    mi::Uint32                    idx,
    const mi::math::Color_struct& color)
{
    const std::string mn = "Transferfunction_colormap::set_color: ";
    WARN_LOG << mn << "invalid call. no effect. ignored.";
}

//----------------------------------------------------------------------
void Transferfunction_colormap::set_colormap(
    const mi::math::Color_struct* colormap_values,
    mi::Uint32                    nb_entries)
{
    const std::string mn = "Transferfunction_colormap::set_colormap: ";
    WARN_LOG << mn << "invalid call. cannot change the colormap. ignored.";
}

//----------------------------------------------------------------------
const mi::math::Color_struct* Transferfunction_colormap::get_colormap() const
{
    assert(!m_color_vec.empty()); // should not access if empty

    return &(m_color_vec[0]);
}

//----------------------------------------------------------------------
mi::Uint32 Transferfunction_colormap::get_number_of_entries() const
{
    return m_color_vec.size();
}

//----------------------------------------------------------------------
void Transferfunction_colormap::get_non_transparent_range(
    mi::Uint32& min,
    mi::Uint32& max) const
{
    min = 0;
    max = 255;
}

//----------------------------------------------------------------------
void Transferfunction_colormap::set_domain(
    mi::Float32 value_first,
    mi::Float32 value_last)
{
    // Note: In general and for full flexibility any value is allowed for the colormap domain, and
    // the following checks here are specific to the implementation of this class.

    const mi::math::Vector<mi::Float32, 2> local_voxel_range(value_first, value_last);
    const std::string mn = "Transferfunction_colormap::set_domain: ";
    if (local_voxel_range[0] >= local_voxel_range[1])
    {
        ERROR_LOG << mn << "illegal local voxel range: " << local_voxel_range;
        return;
    }

    if (local_voxel_range[0] < m_global_voxel_value_domain[0])
    {
        ERROR_LOG << mn << "illegal local voxel min: " << local_voxel_range[0]
                  << ", is out of global range [" << m_global_voxel_value_domain << "]";
        return;
    }

    if (local_voxel_range[1] > m_global_voxel_value_domain[1])
    {
        ERROR_LOG << mn << "illegal local voxel max: " << local_voxel_range[1]
                  << ", is out of global range [" << m_global_voxel_value_domain << "]";
        return;
    }

    m_local_voxel_value_domain = local_voxel_range;
}

//----------------------------------------------------------------------
void Transferfunction_colormap::get_domain(
        mi::Float32& value_first,
        mi::Float32& value_last) const
{
    value_first = m_local_voxel_value_domain[0];
    value_last  = m_local_voxel_value_domain[1];
}

//----------------------------------------------------------------------
void Transferfunction_colormap::set_domain_boundary_mode(Domain_boundary_mode mode)
{
    m_domain_boundary_mode = mode;
}

//----------------------------------------------------------------------
Transferfunction_colormap::Domain_boundary_mode Transferfunction_colormap::get_domain_boundary_mode() const
{
    return m_domain_boundary_mode;
}

//======================================================================
// IAttribute
//----------------------------------------------------------------------
mi::base::Uuid Transferfunction_colormap::get_attribute_class() const
{
    return IColormap::IID();
}

//----------------------------------------------------------------------
bool Transferfunction_colormap::multiple_active_instances() const
{
    return false;
}

//======================================================================
// IScene_element
//----------------------------------------------------------------------
void Transferfunction_colormap::set_enabled(bool enable)
{
    m_is_enabled = enable;
}

//----------------------------------------------------------------------
bool Transferfunction_colormap::get_enabled() const
{
    return m_is_enabled;
}

//----------------------------------------------------------------------
void Transferfunction_colormap::set_meta_data(mi::neuraylib::Tag_struct tag)
{
    // empty
}

//----------------------------------------------------------------------
mi::neuraylib::Tag_struct Transferfunction_colormap::get_meta_data() const
{
    return mi::neuraylib::NULL_TAG;
}

//======================================================================
// IElement
//----------------------------------------------------------------------
mi::neuraylib::IElement* Transferfunction_colormap::copy() const
{
    Transferfunction_colormap * other = new Transferfunction_colormap();
    other->m_is_enabled                = this->m_is_enabled;
    other->m_anchor_vec                = this->m_anchor_vec;
    other->m_global_voxel_value_domain = this->m_global_voxel_value_domain;
    other->m_local_voxel_value_domain  = this->m_local_voxel_value_domain;
    other->m_domain_resolution         = this->m_domain_resolution;
    other->m_is_linear_interpolation   = this->m_is_linear_interpolation;
    other->m_color_space               = this->m_color_space;
    other->m_color_vec                 = this->m_color_vec;
    other->m_domain_boundary_mode      = this->m_domain_boundary_mode;

    return other;
}

//----------------------------------------------------------------------
void Transferfunction_colormap::serialize(mi::neuraylib::ISerializer* serializer) const
{
    serializer->write(&m_is_enabled, 1);

    const mi::Size nb_anchor_point = m_anchor_vec.size();
    serializer->write(&nb_anchor_point, 1);
    for (mi::Size i = 0; i < nb_anchor_point; ++i)
    {
        const mi::Float32     voxel_value = m_anchor_vec[i].m_voxel_value;
        const mi::math::Color color_value = m_anchor_vec[i].m_color_value;
        serializer->write(&voxel_value, 1);
        serializer->write(&(color_value[0]), 4);
    }

    serializer->write(&(m_global_voxel_value_domain.x), 2);
    serializer->write(&(m_local_voxel_value_domain.x),  2);
    serializer->write(&m_domain_resolution, 1);
    serializer->write(&m_is_linear_interpolation, 1);
    const mi::Uint32 color_space = static_cast<mi::Uint32>(m_color_space);
    serializer->write(&color_space, 1);

    const mi::Uint32 mode = m_domain_boundary_mode;
    serializer->write(&mode, 1);

    // Do not serialized m_color_vec;
}

//----------------------------------------------------------------------
void Transferfunction_colormap::deserialize(mi::neuraylib::IDeserializer* deserializer)
{
    deserializer->read(&m_is_enabled, 1);

    m_anchor_vec.clear();
    mi::Size nb_anchor_point = 0;
    deserializer->read(&nb_anchor_point, 1);
    for (mi::Size i = 0; i < nb_anchor_point; ++i)
    {
        Voxel_value_color_pair vvcp;
        deserializer->read(&(vvcp.m_voxel_value),  1);
        deserializer->read(&(vvcp.m_color_value[0]), 4);

        m_anchor_vec.push_back(vvcp);
    }

    deserializer->read(&(m_global_voxel_value_domain.x), 2);
    deserializer->read(&(m_local_voxel_value_domain.x),  2);
    deserializer->read(&m_domain_resolution, 1);
    deserializer->read(&m_is_linear_interpolation, 1);
    mi::Uint32 color_space = 0;
    deserializer->read(&color_space, 1);
    m_color_space = static_cast<Color_space_e>(color_space);

    mi::Uint32 mode;
    deserializer->read(&mode, 1);
    m_domain_boundary_mode = static_cast<Domain_boundary_mode>(mode);

    // Do not deserialized m_color_vec;
    m_color_vec.clear();
}

//----------------------------------------------------------------------
bool Transferfunction_colormap::get_current_idx_left_right_voxel_value(
    mi::Size idx,
    mi::math::Vector<mi::Float32, 2>& left_right)
{
    if (idx == 0)
    {
        return false;
    }

    if (idx >= m_anchor_vec.size())
    {
        return false;
    }

    left_right[0] = m_anchor_vec[idx - 1].m_voxel_value;
    left_right[1] = m_anchor_vec[idx    ].m_voxel_value;

    return true;
}

//----------------------------------------------------------------------
bool Transferfunction_colormap::is_sorted_anchor_vec() const
{
    bool is_sorted = true;

    const mi::Size nb_entries = m_anchor_vec.size();
    for (mi::Size i = 1; i < nb_entries; ++i)
    {
        // exclude the same value for simplicity
        if (m_anchor_vec[i - 1].m_voxel_value >= m_anchor_vec[i].m_voxel_value)
        {
            ERROR_LOG << "DEBUG: is_sorted: i: " << i << ", nb_entries: " << nb_entries;
            is_sorted = false;
            break;
        }
    }

    return is_sorted;
}

//----------------------------------------------------------------------
void Transferfunction_colormap::generate_linear_interpolated_colormap()
{
    const mi::Size color_vec_size = m_domain_resolution;
    m_color_vec.resize(color_vec_size);

    const mi::Float32 domain_range = m_local_voxel_value_domain[1] - m_local_voxel_value_domain[0];
    const mi::Float32 domain_delta = domain_range/(static_cast<mi::Float32>(m_domain_resolution));
    assert(domain_delta > 0.0f);

    for (mi::Size i = 0; i < color_vec_size; ++i)
    {
        const mi::Float32 cur_voxel_val =
            m_local_voxel_value_domain[0] + static_cast<mi::Float32>(i) * domain_delta;
        // WARN_LOG << "update_colormap: idx: " << i << ", voxel_val: " << cur_voxel_val;

        m_color_vec[i] = linear_interpolate_color(cur_voxel_val);
    }
}

//----------------------------------------------------------------------
mi::math::Color Transferfunction_colormap::linear_interpolate_color(
    mi::Float32 voxel_val)
{
    // check the voxel_val is in range
    if ((voxel_val < m_local_voxel_value_domain[0]) ||
        (m_local_voxel_value_domain[1] < voxel_val))
    {
        WARN_LOG << "voxel_val out of range.";
        return mi::math::Color(0.0f, 0.0f, 0.0f, 1.0f);
    }

    // color value is dummy
    Voxel_value_color_pair voxel_only(voxel_val, mi::math::Color(0.0f, 0.0f, 0.0f, 1.0f));

    // find closest anchor (Lower bound is first not less than value)
    std::vector<Voxel_value_color_pair>::iterator left_i =
        std::lower_bound(m_anchor_vec.begin(), m_anchor_vec.end(), voxel_only, Voxel_value_comp());

    std::vector<Voxel_value_color_pair>::iterator right_i = left_i;
    if (voxel_val < left_i->m_voxel_value)
    {
        assert(left_i != m_anchor_vec.begin()); // should not happen.
        --left_i;                               // left is one step back.
    }
    else
    {
        assert(right_i != m_anchor_vec.end()); // should not happen.
        ++right_i;                             // right is one step further
    }

    const mi::Float32 left_val  = left_i ->m_voxel_value;
    const mi::Float32 right_val = right_i->m_voxel_value;
    assert(left_val < right_val);
    const mi::Float32 delta_val = right_val - left_val  ;
    const mi::math::Color delta_col = right_i->m_color_value - left_i->m_color_value;

    const mi::math::Color int_col = (delta_col / delta_val) * (voxel_val - left_val) + left_i->m_color_value;

    // WARN_LOG << "L: " << left_i->get_string() << ", R: " << right_i->get_string()
    //          << ", dc/dv: " << delta_col << "/" << delta_val << " * " << (voxel_val - left_val)
    //          << " + " << left_i->m_color_value << " = " << int_col;

    return int_col;
}

//----------------------------------------------------------------------
void Transferfunction_colormap::generate_spline_interpolated_colormap()
{
    // first: the color r, g, b are linearly interpolated (depends on
    // the application, color can be also interpolated.)

    generate_linear_interpolated_colormap();

    // get the anchor points as 2d points vector
    std::vector< mi::math::Vector<mi::Float32, 3> > point_vec;

    const mi::Size nb_anchor_point = m_anchor_vec.size();
    for (mi::Size i = 0; i < nb_anchor_point; ++i)
    {
        mi::math::Vector<mi::Float32, 3> point_2d(
            m_anchor_vec[i].m_voxel_value, m_anchor_vec[i].m_color_value.a, 0.0f);

        point_vec.push_back(point_2d);
    }

    Cubic_hermite_spline chs;
    const mi::Float32 bias    = 0.0f;
    const mi::Float32 tension = 0.0f;
    chs.build_curve(&(point_vec[0]), point_vec.size(), bias, tension);

    // interpolate alphas
    const mi::Size color_vec_size = m_domain_resolution;
    const mi::Float32 domain_range = m_local_voxel_value_domain[1] - m_local_voxel_value_domain[0];
    const mi::Float32 domain_delta = domain_range/(static_cast<mi::Float32>(m_domain_resolution));
    assert(domain_delta > 0.0f);

    for (mi::Size i = 0; i < color_vec_size; ++i)
    {
        const mi::Float32 cur_voxel_val =
            m_local_voxel_value_domain[0] + static_cast<mi::Float32>(i) * domain_delta;
        // WARN_LOG << "update_colormap: idx: " << i << ", voxel_val: " << cur_voxel_val;
        mi::Float32 interp_alpha = spline_interpolate_color(chs, cur_voxel_val);

        // clump alpha in [0, 1]
        if (interp_alpha < 0.0f)
        {
            interp_alpha = 0.0f;
        }
        else if (interp_alpha > 1.0f)
        {
            interp_alpha = 1.0f;
        }

        m_color_vec[i].a = interp_alpha;
    }
}

//----------------------------------------------------------------------
mi::Float32 Transferfunction_colormap::spline_interpolate_color(
    nv::index_common::Cubic_hermite_spline& chs,
    mi::Float32 voxel_val)
{
    // check the voxel_val is in range
    if ((voxel_val < m_local_voxel_value_domain[0]) ||
        (m_local_voxel_value_domain[1] < voxel_val))
    {
        WARN_LOG << "spline_interpolate_color: voxel_val out of range.";
        return 1.0f;
    }

    // find closest anchor (Lower bound is first not less than value)
    Voxel_value_color_pair voxel_only(voxel_val, mi::math::Color(0.0f, 0.0f, 0.0f, 1.0f)); // dummy
    std::vector<Voxel_value_color_pair>::iterator left_i =
        std::lower_bound(m_anchor_vec.begin(), m_anchor_vec.end(), voxel_only, Voxel_value_comp());

    std::vector<Voxel_value_color_pair>::iterator right_i = left_i;
    if (voxel_val < left_i->m_voxel_value)
    {
        assert(left_i != m_anchor_vec.begin()); // should not happen.
        --left_i;                               // left is one step back.
    }
    else
    {
        assert(right_i != m_anchor_vec.end()); // should not happen.
        ++right_i;                             // right is one step further
    }

    const mi::Float32 left_val  = left_i ->m_voxel_value;
    const mi::Float32 right_val = right_i->m_voxel_value;
    assert(left_val < right_val);
    const mi::Float32 delta_t   = voxel_val - left_val;
    const mi::Float32 delta_val = right_val - left_val;
    const mi::Float32 param_t   = delta_t / delta_val;

    mi::Size anchor_point_idx = left_i - m_anchor_vec.begin();

    mi::math::Vector<mi::Float32, 3> inter_pos = chs.get_value(anchor_point_idx, param_t);

    // WARN_LOG << "L: " << left_i->get_string() << ", R: " << right_i->get_string()
    //          << ", st/dv: " << delta_t << "/" << delta_val;

    return inter_pos.y;
}

//----------------------------------------------------------------------
void Transferfunction_colormap::convert_current_colormap_hsv_to_rgb()
{
    mi::Size nb_entries = m_color_vec.size();
    for (mi::Size i = 0; i < nb_entries; ++i)
    {
        mi::math::Color col = m_color_vec[i];
        m_color_vec[i] = hsv_to_rgb(col);
        // INFO_LOG << "convert[" << i << "]: " << col << " -> " << m_color_vec[i];
    }
}

//----------------------------------------------------------------------
} // namespace index_common
} // namespace nv
