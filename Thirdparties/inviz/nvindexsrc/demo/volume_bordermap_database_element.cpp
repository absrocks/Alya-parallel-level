/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
#include <cassert>

#include <sstream>

#include "volume_bordermap_database_element.h"

#include "common/forwarding_logger.h"
#include "common/string_dict.h"

//----------------------------------------------------------------------
Volume_bordermap_database_element::Volume_bordermap_database_element()
    :
    m_bordertag_map()
{
    // empty
}

//----------------------------------------------------------------------
Volume_bordermap_database_element::~Volume_bordermap_database_element()
{
    // There is a case to have an reference. unless transaction is committed.
    // if(m_bordertag_map.size() > 0){
    //     ERROR_LOG << "Volume_bordermap_database_element: some volume data elements are not released.";
    // }
}

//----------------------------------------------------------------------
void Volume_bordermap_database_element::insert_bbox_volume_tag(
    mi::math::Bbox< mi::Sint64, 3 > const & bbox,
    const mi::neuraylib::Tag& volume_tag)
{
    std::string const key = nv::index_common::to_string(bbox);
    if(m_bordertag_map.find(key) != m_bordertag_map.end()){
        ERROR_LOG << "volume_tag[" << volume_tag.id << "] exists. ignored.";
        return;
    }
    m_bordertag_map[key] = volume_tag;
}

//----------------------------------------------------------------------
void Volume_bordermap_database_element::insert_bordermap(
    Volume_bordermap_database_element const * p_bmap)
{
    for(Border_map::const_iterator bi = p_bmap->get_bordermap_ref()->begin();
        bi != p_bmap->get_bordermap_ref()->end(); ++bi)
    {
        if(m_bordertag_map.find(bi->first) == m_bordertag_map.end()){
            m_bordertag_map[bi->first] = bi->second;
        }
        else{
            DEBUG_LOG << "duplication found, ignored.";
        }
    }
}

//----------------------------------------------------------------------
Volume_bordermap_database_element::Border_map const *
Volume_bordermap_database_element::get_bordermap_ref() const
{
    return & m_bordertag_map;
}

//----------------------------------------------------------------------
void Volume_bordermap_database_element::clear_bordermap(
    mi::neuraylib::IDice_transaction * dice_transaction)
{
    assert(dice_transaction != 0);

    // clean up all the volume tag element
    for(Border_map::const_iterator bi = m_bordertag_map.begin();
        bi != m_bordertag_map.end(); ++bi)
    {
        mi::neuraylib::Tag const volume_tag = bi->second;
        assert(volume_tag.is_valid());
        dice_transaction->remove(volume_tag);
        // DEBUG_LOG << "Volume_bordermap_database_element::clear_bordermap:: removed tag: "
        //           << volume_tag.id;
    }

    // clear tags
    m_bordertag_map.clear();
    assert(m_bordertag_map.empty());
}

//----------------------------------------------------------------------
mi::neuraylib::Tag Volume_bordermap_database_element::get_element_tag(
    mi::math::Bbox< mi::Sint64, 3 > const & bbox) const
{
    std::string const key = nv::index_common::to_string(bbox);
    Border_map::const_iterator bi = m_bordertag_map.find(key);
    if(bi == m_bordertag_map.end()){
        return mi::neuraylib::NULL_TAG;
    }

    return bi->second;
}

//----------------------------------------------------------------------
std::string Volume_bordermap_database_element::to_string() const
{
    std::stringstream sstr;
    sstr << "size = " << m_bordertag_map.size() << "\n";
    for(Border_map::const_iterator bi = m_bordertag_map.begin();
        bi != m_bordertag_map.end(); ++bi)
    {
        sstr << bi->first << " " << bi->second.id << "\n";
    }
    return sstr.str();
}

//----------------------------------------------------------------------
mi::neuraylib::IElement* Volume_bordermap_database_element::copy() const
{
    Volume_bordermap_database_element* other = new Volume_bordermap_database_element;
    other->m_bordertag_map = this->m_bordertag_map;

    return other;
}

//----------------------------------------------------------------------
void Volume_bordermap_database_element::serialize(
    mi::neuraylib::ISerializer *serializer) const
{
    mi::Sint32 const nb_elements = static_cast< mi::Sint32 >(m_bordertag_map.size());
    serializer->write(&nb_elements, 1);

    for(Border_map::const_iterator bi = m_bordertag_map.begin();
        bi != m_bordertag_map.end(); ++bi)
    {
        // key string
        mi::Uint32 const len = mi::Uint32(bi->first.size());
        serializer->write(&len, 1);
        serializer->write(reinterpret_cast<const mi::Uint8*>(&(bi->first[0])), len);

        // value tag
        serializer->write(&(bi->second.id), 1);
    }
}

//----------------------------------------------------------------------
void Volume_bordermap_database_element::deserialize(
    mi::neuraylib::IDeserializer* deserializer)
{
    m_bordertag_map.clear();

    mi::Sint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    for(mi::Sint32 i = 0; i < nb_elements; ++i){
        // key string
        mi::Uint32 len = 0;
        deserializer->read(&len, 1);
        std::string key;
        key.resize(len);
        deserializer->read(reinterpret_cast< mi::Uint8* >(&key[0]), len);

        // value tag
        mi::neuraylib::Tag tag;
        deserializer->read(&tag.id, 1);

        assert(m_bordertag_map.find(key) == m_bordertag_map.end());
        m_bordertag_map[key] = tag;
    }
}

//----------------------------------------------------------------------
