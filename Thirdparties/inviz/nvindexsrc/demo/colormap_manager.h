/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief application defined colormap manager

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_COLORMAP_MANAGER_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_COLORMAP_MANAGER_H

#include <map>
#include <vector>
#include <cassert>

#include <mi/dice.h>


//----------------------------------------------------------------------

/// application colormap manager sample implementation.
///
/// This manager handles two maps:
///   - id -> colormap tag
class Colormap_manager
{
public:
    /// colormap id -> tag
    typedef std::map< mi::Uint32, mi::neuraylib::Tag > id_tag_map;

public:
    /// default constructor
    Colormap_manager();
    /// destructor
    virtual ~Colormap_manager();

public:
    //------------------------------------------------------------
    // {index -> tag, filename} mapping related methods
    //------------------------------------------------------------

    /// get number of colormap.
    /// some of them might not valid.
    /// \return number of colormaps according to the submitted tags.
    mi::Uint32 get_number_of_colormap() const;

    /// get colormap filename
    /// \param[in] id colormap id
    /// \return colormap tag by id. "" when invalid id.
    // std::string get_colormap_filename(mi::Uint32 id) const;

    /// get colormap tag
    /// \param[in] id colormap id
    /// \return colormap dice tag by id. NULL_TAG when invalid id.
    mi::neuraylib::Tag get_colormap_tag(mi::Uint32 id) const;

    /// set tag with id
    /// \param[in] id  colormap id
    /// \param[in] tag dice tag, should not be NULL_TAG
    /// \return true when the set succeeded.
    bool set_tag(mi::Uint32 id,
                 mi::neuraylib::Tag tag);

    /// is id's tag valid?
    /// \param[in] id colormap id
    /// \return true when id is valid for the tag
    bool is_valid_tag(mi::Uint32 id) const;

public:
    //------------------------------------------------------------
    // colormap userdefined operation related methods
    //------------------------------------------------------------

    /// set operation type
    ///
    /// \param[in] op_type operation type
    void set_colormap_userdef_operation_type(mi::Sint32 op_type);

    /// get operation type
    mi::Sint32  get_colormap_userdef_operation_type() const;

private:
    /// id -> tag
    id_tag_map    m_colormap_tag_map;

    /// user defined colormap operation type
    mi::Sint32    m_colormap_userdef_operation_type;

private:
    /// copy constructor. prohibit until proved useful.
    Colormap_manager(Colormap_manager const &);
    /// operator=. prohibit until proved useful.
    Colormap_manager const & operator=(Colormap_manager const &);
};
//----------------------------------------------------------------------


#endif // #ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_COLORMAP_MANAGER_H
