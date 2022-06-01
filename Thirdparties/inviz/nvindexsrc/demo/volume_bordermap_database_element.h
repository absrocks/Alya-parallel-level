/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief border tag database element for inset filter operation.

#ifndef GEOSPATIAL_BIN_GEOSPATIAL_STREAM_VIEWER_BORDERMAP_DATABASE_ELEMENT_H
#define GEOSPATIAL_BIN_GEOSPATIAL_STREAM_VIEWER_BORDERMAP_DATABASE_ELEMENT_H

#include <mi/dice.h>

#include <map>
#include <string>

/// volume border data set database element.
class Volume_bordermap_database_element :
        public mi::neuraylib::Element<0x7af30e68,0xd7fe,
                                      0x4d64,0x93,0x89,0x63,0xd6,0x29,0xf0,0x3e,0x4e,
                                      mi::neuraylib::IElement>
{
public:
    /// border map
    typedef std::map< std::string, mi::neuraylib::Tag > Border_map;

public:
    /// Default c'tor for serialization only
    Volume_bordermap_database_element();

    /// Destructor
    virtual ~Volume_bordermap_database_element();

    /// insert volume database element tag associated with the bounding box
    ///
    /// \param[in] bbox volume element bounding box
    /// \param[in] volume_tag volume data tag
    void insert_bbox_volume_tag(mi::math::Bbox< mi::Sint64, 3 > const & bbox,
                                const mi::neuraylib::Tag& volume_tag);

    /// append other Volume_bordermap_database_element. If duplication found,
    /// ignore it.
    ///
    /// \param[in] p_bmap  other bordermap_database_element
    void insert_bordermap(Volume_bordermap_database_element const * p_bmap);

    /// get borders tag set reference
    Border_map const * get_bordermap_ref() const;

    /// clear border tags and associated volume_scene_elements.
    void clear_bordermap(
        mi::neuraylib::IDice_transaction *dice_transaction);

    /// get the border element tag which corresponds to the bbox
    ///
    /// \param[in] bbox volume bounding box
    /// \return tag associated with bbox. NULL_TAG when not found.
    mi::neuraylib::Tag get_element_tag(mi::math::Bbox< mi::Sint64, 3 > const & bbox) const;

    /// get string representation of this object
    std::string to_string() const;

public:
    /// -------------------------------------------------------------------------------------------

    /// The copy needs to create and return a newly allocated object of the same class and has to
    /// copy all the fields which should be the same in the copy.
    /// \return The new copy
    virtual mi::neuraylib::IElement* copy() const;

    /// Give back a human readable representation of the class name
    virtual const char* get_class_name() const
    {
        return "Volume_bordermap_database_element";
    }

    /// Return the unique class id of this class
    /// \return The class id
    virtual mi::base::Uuid get_class_id() const
    {
        return IID();
    }

    /// Serialize the class to the given serializer.
    /// \param serializer Write to this serializer.
    virtual void serialize(mi::neuraylib::ISerializer *serializer) const;

    /// Deserialize the class from the given deserializer.
    /// \param deserializer Read from this deserializer.
    virtual void deserialize(mi::neuraylib::IDeserializer* deserializer);

    /// get reference
    virtual void get_references(mi::neuraylib::ITag_set* result) const {}

private:
    /// border tag map strage.
    Border_map m_bordertag_map;
};


#endif // #ifndef GEOSPATIAL_BIN_GEOSPATIAL_STREAM_VIEWER_BORDERMAP_DATABASE_ELEMENT_H
