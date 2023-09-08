/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "colormap_manager.h"

#include <sstream>

#include "common/forwarding_logger.h"

//----------------------------------------------------------------------
Colormap_manager::Colormap_manager()
    :
    m_colormap_userdef_operation_type(0)
{
    // empty
}

//----------------------------------------------------------------------
Colormap_manager::~Colormap_manager()
{
    // actually not necessary, but cleaner.
    m_colormap_tag_map.clear();
}

//----------------------------------------------------------------------
mi::Uint32 Colormap_manager::get_number_of_colormap() const
{
    return static_cast< mi::Uint32 >(m_colormap_tag_map.size());
}

//----------------------------------------------------------------------
mi::neuraylib::Tag Colormap_manager::get_colormap_tag(mi::Uint32 id) const
{
    assert(this->is_valid_tag(id));

    id_tag_map::const_iterator mi = m_colormap_tag_map.find(id);
    return (mi->second);
}

//----------------------------------------------------------------------
bool Colormap_manager::set_tag(mi::Uint32 id,
                               mi::neuraylib::Tag tag)
{
    // check the id as tag id. should not be duplicated.
    if (this->is_valid_tag(id))
    {
        ERROR_LOG << "Colormap_manager::set_tag: duplicated tag. id = (" << id << ").";
        return false;
    }

    if (!tag.is_valid())
    {
        ERROR_LOG << "Colormap_manager::set_tag: cannot insert NULL_TAG. id = (" << id << ").";
        return false;
    }

    m_colormap_tag_map[id] = tag;
    return true;
}

//----------------------------------------------------------------------
bool Colormap_manager::is_valid_tag(mi::Uint32 id) const
{
    id_tag_map::const_iterator mi = m_colormap_tag_map.find(id);
    if ((mi != m_colormap_tag_map.end()) && (mi->second.is_valid()))
    {
        return true;
    }
    return false;
}

//----------------------------------------------------------------------
void Colormap_manager::set_colormap_userdef_operation_type(mi::Sint32 op_type)
{
    m_colormap_userdef_operation_type = op_type;
    DEBUG_LOG << "set op_type = " << this->get_colormap_userdef_operation_type();
}

//----------------------------------------------------------------------
mi::Sint32 Colormap_manager::get_colormap_userdef_operation_type() const
{
    return m_colormap_userdef_operation_type;
}

//----------------------------------------------------------------------
