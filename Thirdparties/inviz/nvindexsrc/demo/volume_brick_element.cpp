/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "volume_brick_element.h"

#include <mi/dice.h>

#include "common/forwarding_logger.h"

#include <cassert>

//----------------------------------------------------------------------
Volume_brick_element::Volume_brick_element(
    const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
    const mi::Uint8*                                voxel_data)
    :
    m_bounding_box(bounding_box)
{
    mi::math::Bbox<mi::Sint32, 3> check_bounding_box(bounding_box);
    if(!check_bounding_box.is_volume()){
        ERROR_LOG << "cannot initialize."; // FIXME put this in the initialize after the fix.
        return;
    }

    const mi::Uint32 nb_values =
        (this->m_bounding_box.max.x - this->m_bounding_box.min.x) *
        (this->m_bounding_box.max.y - this->m_bounding_box.min.y) *
        (this->m_bounding_box.max.z - this->m_bounding_box.min.z);

    m_voxel_data = new mi::Uint8[nb_values];
    for(mi::Uint32 i=0; i<nb_values; ++i)
    {
        m_voxel_data[i] = voxel_data[i];
    }
}

//----------------------------------------------------------------------
Volume_brick_element::Volume_brick_element()
    :
    m_voxel_data(0)
{
    // empty
}

//----------------------------------------------------------------------
Volume_brick_element::~Volume_brick_element()
{
    free_local_memory();
}

//----------------------------------------------------------------------
const mi::Uint8* Volume_brick_element::get_voxel_data() const
{
    assert(m_voxel_data != 0);
    return &m_voxel_data[0];
}

//----------------------------------------------------------------------
//TODO: This may include the border, but do the functions using this know about that?
mi::math::Bbox_struct<mi::Sint32, 3> Volume_brick_element::get_bounding_box() const
{
    return m_bounding_box;
}

//----------------------------------------------------------------------
void Volume_brick_element::free_local_memory() const
{
    if (m_voxel_data)
    {
        delete[] m_voxel_data;
        m_voxel_data = 0;
    }
}

//----------------------------------------------------------------------
mi::neuraylib::IElement* Volume_brick_element::copy() const
{
    Volume_brick_element* other = new Volume_brick_element;

    other->m_bounding_box = this->m_bounding_box;

    const mi::Uint32 nb_values =
        (this->m_bounding_box.max.x - this->m_bounding_box.min.x) *
        (this->m_bounding_box.max.y - this->m_bounding_box.min.y) *
        (this->m_bounding_box.max.z - this->m_bounding_box.min.z);

    other->m_voxel_data = new mi::Uint8[nb_values];
    for(mi::Uint32 i=0; i<nb_values; ++i)
    {
        other->m_voxel_data[i] = this->m_voxel_data[i];
    }

    return other;
}

//----------------------------------------------------------------------
void Volume_brick_element::serialize(
    mi::neuraylib::ISerializer *serializer) const
{
    serializer->write(&m_bounding_box.min.x, 3);
    serializer->write(&m_bounding_box.max.x, 3);

    const mi::Uint32 nb_values =
        (m_bounding_box.max.x - m_bounding_box.min.x) *
        (m_bounding_box.max.y - m_bounding_box.min.y) *
        (m_bounding_box.max.z - m_bounding_box.min.z);

    serializer->write(&m_voxel_data[0], nb_values);

}

//----------------------------------------------------------------------
void Volume_brick_element::deserialize(
    mi::neuraylib::IDeserializer* deserializer)
{
    deserializer->read(&m_bounding_box.min.x, 3);
    deserializer->read(&m_bounding_box.max.x, 3);

    const mi::Uint32 nb_values =
        (m_bounding_box.max.x - m_bounding_box.min.x) *
        (m_bounding_box.max.y - m_bounding_box.min.y) *
        (m_bounding_box.max.z - m_bounding_box.min.z);

    m_voxel_data = new mi::Uint8[nb_values];
    deserializer->read(&m_voxel_data[0], nb_values);
}

//----------------------------------------------------------------------
