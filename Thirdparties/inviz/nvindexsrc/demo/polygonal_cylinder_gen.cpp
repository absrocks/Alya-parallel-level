/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "polygonal_cylinder_gen.h"

#include <cassert>
#include <cstddef>
#include <cmath>
#include <sstream>

#include "common/forwarding_logger.h"

//----------------------------------------------------------------------
Polygonal_cylinder_gen::Polygonal_cylinder_gen()
    :
    m_poly_normal (0.0, 0.0, 1.0),
    m_poly_tangent(1.0, 0.0, 0.0),
    m_poly_binomal(0.0, 1.0, 0.0),
    // m_center_vec(),
    // m_radius_vec(),
    m_n_gon(0),
    m_is_gen_segment_tris(false)
    // m_vertex_vec(),
    // m_segment_face_vec(),
    // m_side_face_vec()
{
    // empty
}

//----------------------------------------------------------------------
Polygonal_cylinder_gen::~Polygonal_cylinder_gen()
{
    this->clear();
}

//----------------------------------------------------------------------
void Polygonal_cylinder_gen::clear()
{
    m_center_vec.clear();
    m_radius_vec.clear();

    m_vertex_vec      .clear();
    m_segment_face_vec.clear();
    m_side_face_vec   .clear();
}

//----------------------------------------------------------------------
bool Polygonal_cylinder_gen::append_center_point(
    mi::math::Vector< mi::Float32, 3 > const & center_point,
    mi::Float32 radius)
{
    assert(radius > 0.0);
    // sanity check. without this check still generates cylinders,
    // but, may be looked broken
    size_t list_sz = m_center_vec.size();
    if ((list_sz > 0) &&  (m_center_vec.at(list_sz - 1)[2] >= center_point[2])){
        ERROR_LOG << "z values must be acendent order. "
                  << (m_center_vec.at(list_sz - 1)) << ", " << center_point;
        return false;
    }

    m_center_vec.push_back(center_point);
    m_radius_vec.push_back(radius);

    return true;
}

//----------------------------------------------------------------------
void Polygonal_cylinder_gen::set_generate_segment_tris(bool is_gen)
{
    m_is_gen_segment_tris = is_gen;
}

//----------------------------------------------------------------------
void Polygonal_cylinder_gen::set_n_gon(mi::Sint32 n_gon)
{
    assert(n_gon >= 3);

    m_n_gon = n_gon;
}

//----------------------------------------------------------------------
void Polygonal_cylinder_gen::gen_cylinder()
{
    assert(this->is_able_to_gen());

    this->gen_vertex();
    if (m_is_gen_segment_tris){
        this->gen_segment_tris();
    }
    this->gen_side_polygon();
}

//----------------------------------------------------------------------
bool Polygonal_cylinder_gen::has_cylinder() const
{
    if(m_vertex_vec.empty() && m_side_face_vec.empty()){
        return false;
    }
    if((!m_vertex_vec.empty()) && (!m_side_face_vec.empty())){
        return true;
    }
    // Please add something like "FATAL_LOG" to declare errors as internal library errors.
    ERROR_LOG << "Invalid cylinder; please verify the input parameter ["
              << "vertex size: "      << m_vertex_vec.size()
              << ", side face size: " << m_side_face_vec.size()
              << "].";
    return false;
}

//----------------------------------------------------------------------
Polygonal_cylinder_gen::Float32_3_vec const & Polygonal_cylinder_gen::get_vertex_vec_ref() const
{
    return m_vertex_vec;
}

//----------------------------------------------------------------------
Polygonal_cylinder_gen::Sint32_3_vec  const & Polygonal_cylinder_gen::get_segment_face_ref() const
{
    return m_segment_face_vec;
}

//----------------------------------------------------------------------
Polygonal_cylinder_gen::Sint32_3_vec  const & Polygonal_cylinder_gen::get_side_face_ref() const
{
    return m_side_face_vec;
}

//----------------------------------------------------------------------
bool Polygonal_cylinder_gen::is_able_to_gen() const
{
    if (m_n_gon < 3){
        ERROR_LOG << "n-gon should be at least n >= 3.";
        return false;
    }
    if (m_center_vec.size() < 2) {
        ERROR_LOG << "The number of center points must be > 1.";
        return false;
    }
    return true;
}

//----------------------------------------------------------------------
void Polygonal_cylinder_gen::gen_vertex()
{
    m_vertex_vec.clear();
    mi::Float64 const step_angle = (2 * M_PI) / static_cast< mi::Float64 >(m_n_gon);

    for(size_t seg = 0; seg < m_center_vec.size(); ++seg){
        for(mi::Sint32 i = 0; i < m_n_gon; ++i){
            mi::Float64 const rad  = m_radius_vec[seg];
            mi::math::Vector< mi::Float32, 3 > const cp = m_center_vec[seg];
            mi::math::Vector< mi::Float32, 3 > vtx_point(0.0, 0.0, 0.0);
            vtx_point.x = static_cast<mi::Float32>(rad * cos(i * step_angle)) + cp[0];
            vtx_point.y = static_cast<mi::Float32>(rad * sin(i * step_angle)) + cp[1];
            vtx_point.z = cp[2];
            m_vertex_vec.push_back(vtx_point);
        }
    }
}

//----------------------------------------------------------------------
void Polygonal_cylinder_gen::gen_segment_tris()
{
    assert((m_vertex_vec.size()) == (m_center_vec.size() * m_n_gon));

    m_segment_face_vec.clear();
    mi::Sint32 const n = m_n_gon;
    for(size_t seg = 0; seg < m_center_vec.size(); ++seg){
        // base index of current processing segment triangles
        mi::Sint32 const n_seg = seg * n;
        for(mi::Sint32 i = 1; i < ((n + 1) - 2); ++i){
            mi::math::Vector< mi::Sint32, 3 > const fidx(0 + n_seg, i + n_seg, i + 1 + n_seg);
            m_segment_face_vec.push_back(fidx);
        }
    }
}

//----------------------------------------------------------------------
void Polygonal_cylinder_gen::gen_side_polygon()
{
    mi::Sint32 const seg_count = m_center_vec.size();
    assert(seg_count >= 2);

    m_side_face_vec.clear();
    mi::Sint32 const n = m_n_gon;

    for(mi::Sint32 seg = 0; seg < (seg_count - 1); ++seg){
        mi::Sint32 const bidx = seg * n;
        for(mi::Sint32 i = 0; i < (n - 1); ++i){
            mi::math::Vector< mi::Sint32, 3 > const f1(bidx + i, bidx + i + 1, bidx + n + i + 1);
            mi::math::Vector< mi::Sint32, 3 > const f2(bidx + i, bidx + n + i + 1, bidx + n + i);
            m_side_face_vec.push_back(f1);
            m_side_face_vec.push_back(f2);
        }
        // the last quad of this segment
        mi::math::Vector< mi::Sint32, 3 > const l_f1(bidx + n - 1, bidx + 0, bidx + n);
        mi::math::Vector< mi::Sint32, 3 > const l_f2(bidx + n - 1, bidx + n, bidx + (2 * n - 1));
        m_side_face_vec.push_back(l_f1);
        m_side_face_vec.push_back(l_f2);
    }
}

//----------------------------------------------------------------------
bool Polygonal_cylinder_gen::is_face_index_valid() const
{
    mi::Sint32 const vsize = m_vertex_vec.size();

    for(Sint32_3_vec::const_iterator vi = m_segment_face_vec.begin();
        vi != m_segment_face_vec.end(); ++vi)
    {
        for(mi::Sint32 i = 0; i < 3; ++i){
            mi::Sint32 const fidx = (*vi)[i];
            if((fidx < 0) || (fidx >= vsize)){
                ERROR_LOG << "segment face list has invalid face index.";
                return false;
            }
        }
    }

    for(Sint32_3_vec::const_iterator vi = m_side_face_vec.begin();
        vi != m_side_face_vec.end(); ++vi)
    {
        for(mi::Sint32 i = 0; i < 3; ++i){
            mi::Sint32 const fidx = (*vi)[i];
            if((fidx < 0) || (fidx >= vsize)){
                ERROR_LOG << "side face list has invalid face index.";
                return false;
            }
        }
    }

    return true;
}

//----------------------------------------------------------------------
std::string Polygonal_cylinder_gen::to_string() const
{
    std::stringstream sstr;
    sstr << "Polygonal_cylinder_gen: "
         << m_center_vec.size()       << " centers, "
         << m_n_gon                   << "-gon, "
         << m_vertex_vec.size()       << " vertices, "
         << m_segment_face_vec.size() << " seg tris, "
         << m_side_face_vec.size()    << " side faces";
    return sstr.str();
}

//----------------------------------------------------------------------
