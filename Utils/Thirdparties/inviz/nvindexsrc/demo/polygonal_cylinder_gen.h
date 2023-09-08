/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Generate a (polygonal) cylinder from a center point list

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_POLYGONAL_CYLINDER_GEN_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_POLYGONAL_CYLINDER_GEN_H

#include <mi/math/vector.h>
#include <vector>
#include <string>


/// Simple polygonal cylinder generator.
///
/// The normal of top and bottom polygons are always z+ (0,0,1)
class Polygonal_cylinder_gen
{
public:
    /// Float32_3 std vector
    typedef std::vector< mi::math::Vector< mi::Float32, 3 > > Float32_3_vec;
    /// Sint32_3 std vector
    typedef std::vector< mi::math::Vector< mi::Sint32, 3 >  > Sint32_3_vec;


public:
    /// constructor
    Polygonal_cylinder_gen();
    /// destructor
    virtual ~Polygonal_cylinder_gen();

    /// clear the input center points and generated data
    void clear();

    /// set polygon parameter
    /// \param[in] center_point polygon center point
    /// \param[in] radius       polygon radius
    /// \return true when append succeeded
    bool append_center_point(mi::math::Vector< mi::Float32, 3 > const & center_point,
                             mi::Float32 radius);

    /// set generate top and bottom triangles.
    /// \param[in] is_gen generate top and bottom triangles when True
    void set_generate_segment_tris(bool is_gen);

    /// set n of n-gon.
    /// \param[in] n_gon top and bottom polygon n-gon
    void set_n_gon(mi::Sint32 n_gon);

    /// generate a cylinder
    void gen_cylinder();

    /// Has the cylinder? (Has the cylinder been generated?)
    /// \return true when a cylinder exists.
    bool has_cylinder() const;

    /// get the vertices vector reference
    /// \return vertices vector reference
    Float32_3_vec const & get_vertex_vec_ref() const ;

    /// get generated segment face vector reference
    /// \return segment face vector reference. may empty when
    /// set_generate_segment_tris(false)
    Sint32_3_vec  const & get_segment_face_ref() const ;

    /// get generated side face vector reference
    /// \return side face vector reference
    Sint32_3_vec  const & get_side_face_ref() const ;

private:
    /// check we can generate a cylinder.
    // raise an exception if not.
    bool is_able_to_gen() const;

    /// generate vertices for the cylinder
    void gen_vertex();

    /// generate segment triangles (horizontal top/bottom of the cylinder).
    /// All triangles faces z+ direction.
    /// gen_vertex should be run before.
    void gen_segment_tris();

    /// generate cylinder side polygons
    void gen_side_polygon();

    /// check the face index's validity
    /// raise an exception when not valid
    bool is_face_index_valid() const;

    /// string representation
    /// \return string representation of this object
    std::string to_string() const;

private:
    /// segment polygon normal. basis n
    mi::math::Vector< mi::Float32, 3 > m_poly_normal;
    /// segment polygon tangent. basis t
    mi::math::Vector< mi::Float32, 3 > m_poly_tangent;
    /// segment polygon binomal. basis b
    mi::math::Vector< mi::Float32, 3 > m_poly_binomal;
    /// center point vector
    std::vector< mi::math::Vector< mi::Float32, 3 > > m_center_vec;
    /// radius vector
    std::vector< mi::Float32 >    m_radius_vec;
    /// n-gon. fixed for a cylinder
    mi::Sint32 m_n_gon;
    /// switch to generate segment polygon triangles
    bool m_is_gen_segment_tris;
    /// generated vertices vector
    Float32_3_vec m_vertex_vec;
    /// generated segment face vector
    Sint32_3_vec  m_segment_face_vec;
    /// generated side face vector
    Sint32_3_vec  m_side_face_vec;
};


#endif // #ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_POLYGONAL_CYLINDER_GEN_H
