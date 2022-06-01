/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "pipe.h"

namespace nv {
namespace index_common {

// constructor
Pipe::Pipe(const std::vector<mi::math::Vector_struct<mi::Float32, 3> > & points)
    : m_points(points) {}

// constructor
Pipe::Pipe(const std::vector<mi::math::Vector_struct<mi::Float32, 3> > & points,
	   const std::vector<mi::Float32> & radii)
    : m_points(points), m_radii(radii) {}

// constructor
Pipe::Pipe(const std::vector<mi::math::Vector_struct<mi::Float32, 3> > & points,
	   const std::vector<mi::Float32> & radii,
	   const std::vector<mi::Uint16> & materials)
    : m_points(points), m_radii(radii), m_materials(materials) {}

void Pipe::serialize(mi::neuraylib::ISerializer *serializer) const
{
    for (mi::Uint32 i=0; i < m_points.size(); i++)
	serializer->write(&m_points[i].x, 3);

    for (mi::Uint32 i=0; i < m_radii.size(); i++)
	serializer->write(&m_radii[i], 1);

    for (mi::Uint32 i=0; i < m_materials.size(); i++)
	serializer->write(&m_materials[i], 1);
}

void Pipe::deserialize(mi::neuraylib::IDeserializer* deserializer)
{
    mi::Uint32 n;
    deserializer->read(&n, 1);
    m_points.resize(n);
    for (mi::Uint32 i=0; i < n; i++)
	deserializer->read(&m_points[i].x, 3);

    deserializer->read(&n, 1);
    m_radii.resize(n);
    for (mi::Uint32 i=0; i < n; i++)
	deserializer->read(&m_radii[i], 1);

    deserializer->read(&n, 1);
    m_materials.resize(n);
    for (mi::Uint32 i=0; i < n; i++)
	deserializer->read(&m_materials[i], 1);
}

//----------------------------------------------------------------------
// Constructor
Pipe_set::Pipe_set(
    const std::vector<Pipe> & pipes,
    mi::Float32		radius,				// default radius for all pipes in set
    mi::math::Bbox_struct<mi::Float32, 3> & bbox)	// bbox of pipes with radii taken into account
{
    m_pipes = pipes;
    m_radius = radius;
    m_bbox = bbox;
    m_rendering_enabled = true;
    m_pickable = true;
    m_meta_data = mi::neuraylib::NULL_TAG;
}

mi::math::Bbox_struct<mi::Float32, 3> Pipe_set::get_bounding_box() const
{
    mi::math::Bbox_struct<mi::Float32, 3> bbox = m_bbox;
    return bbox;
}

const nv::index::IPipe* Pipe_set::get_pipe(mi::Uint32 index) const
{
    return &m_pipes[index];
}

mi::Size Pipe_set::get_nb_pipes() const
{
    return m_pipes.size();
}

mi::Float32 Pipe_set::get_radius() const
{
    return m_radius;
}

//----------------------------------------------------------------------
mi::neuraylib::IElement* Pipe_set::copy() const
{
    Pipe_set* other = new Pipe_set();
    other->m_pipes              = this->m_pipes;
    other->m_radius             = this->m_radius;
    other->m_bbox               = this->m_bbox;
    other->m_rendering_enabled  = this->m_rendering_enabled;
    other->m_pickable           = this->m_pickable;
    other->m_meta_data          = this->m_meta_data;
    return other;
}

//----------------------------------------------------------------------
const char* Pipe_set::get_class_name() const
{
    return "Pipe_set";
}

//----------------------------------------------------------------------
void Pipe_set::serialize(
    mi::neuraylib::ISerializer *serializer) const
{
    const mi::Size nb_elements = m_pipes.size();
    serializer->write(&nb_elements, 1);

    for (mi::Size i=0; i < nb_elements; i++)
	m_pipes[i].serialize(serializer);

    serializer->write(&m_radius, 1);

    serializer->write(&m_bbox.min.x, 3);
    serializer->write(&m_bbox.max.x, 3);

    serializer->write(&m_rendering_enabled, 1);
    serializer->write(&m_pickable, 1);
    serializer->write(&m_meta_data.id, 1);
}

//----------------------------------------------------------------------
void Pipe_set::deserialize(
    mi::neuraylib::IDeserializer* deserializer)
{
    mi::Uint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);

    m_pipes.resize(nb_elements);
    for (mi::Uint32 i=0; i < nb_elements; i++)
	m_pipes[i].deserialize(deserializer);

    deserializer->read(&m_radius, 1);

    deserializer->read(&m_bbox.min.x, 3);
    deserializer->read(&m_bbox.max.x, 3);

    deserializer->read(&m_rendering_enabled, 1);
    deserializer->read(&m_pickable, 1);
    deserializer->read(&m_meta_data.id, 1);
}

//----------------------------------------------------------------------
}} // namespace nv::index_common
