/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief This class stores a subvolume as a volume brick database element.

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VOLUME_BRICK_ELEMENT_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VOLUME_BRICK_ELEMENT_H

#include <mi/dice.h>
#include <mi/math/bbox.h>

/// External implementation of the volume brick database element for
/// computing volume.
class Volume_brick_element :
    public mi::neuraylib::Element<0x20622925,0xbf71,0x40cb,0xbf,0x27,0xf0,0xd8,0x58,0xc9,0x8f,0xcb,
                                  mi::neuraylib::IElement>
{
public:
    /// Constructor of the subvolume database element
    /// voxel_range: I, J, K range inside the entire volume, definition range
    /// of the subvolume.
    ///
    /// Note: the contents of data (voxel_data) is copied.
    ///
    /// \param[in] bbox       volume data bounding box
    /// \param[in] voxel_data volume data, copied into the database element.
    Volume_brick_element(
        const mi::math::Bbox_struct<mi::Sint32, 3>& bbox,
        const mi::Uint8*                            voxel_data);

    /// Default c'tor for serialization only
    Volume_brick_element();

    /// Destructor
    virtual ~Volume_brick_element();

    /// Return pointer to volume  data
    virtual const mi::Uint8* get_voxel_data() const;

    /// Return volume bounding box
    /// \return bounding box of the volume
    virtual mi::math::Bbox_struct<mi::Sint32, 3> get_bounding_box() const;

    /// Free local memory allocated for this database element. This should only be called together
    /// with IDice_transaction::invalidate_job_results().
    virtual void free_local_memory() const;

public:
    // Implementation of IElement

    virtual mi::neuraylib::IElement* copy() const;

    virtual const char* get_class_name() const
    {
        return "Volume_brick_element";
    }

    virtual mi::base::Uuid get_class_id() const
    {
        return IID();
    }

    virtual void serialize(mi::neuraylib::ISerializer *serializer) const;

    virtual void deserialize(mi::neuraylib::IDeserializer* deserializer);

    virtual void get_references(mi::neuraylib::ITag_set* result) const 
    {
        // empty
    }

private:
    /// volume brick bounding box    
    mi::math::Bbox_struct<mi::Sint32, 3> m_bounding_box; 
    /// volume data
    mutable mi::Uint8*                   m_voxel_data;
};

#endif // #ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VOLUME_BRICK_ELEMENT_H
