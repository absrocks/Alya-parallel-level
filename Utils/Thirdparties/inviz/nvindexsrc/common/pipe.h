/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief implements Pipe, Pipe_set classes.

#ifndef NVIDIA_INDEX_PIPE_PRIMITIVE_H
#define NVIDIA_INDEX_PIPE_PRIMITIVE_H

#include <mi/dice.h>
#include <mi/base/interface_implement.h>
#include <mi/base/handle.h>
#include <mi/base/uuid.h>

#include <nv/index/ipipe.h>

#include <vector>

namespace nv {
namespace index_common {

/// Pipe shape implementation
class Pipe : public nv::index::IPipe
{
public:
    // constructor
    Pipe(const std::vector<mi::math::Vector_struct<mi::Float32, 3> > & points);
    Pipe(const std::vector<mi::math::Vector_struct<mi::Float32, 3> > & points,
	 const std::vector<mi::Float32> & radii);
    Pipe(const std::vector<mi::math::Vector_struct<mi::Float32, 3> > & points,
	 const std::vector<mi::Float32> & radii,
	 const std::vector<mi::Uint16> & materials);
    Pipe() {}

    // destructor
    ~Pipe() {}

    mi::Size get_nb_points() const { return m_points.size(); }
    const mi::math::Vector_struct<mi::Float32, 3> * get_points() const { return &m_points[0]; }

    mi::Size get_nb_radii() const { return m_radii.size(); }
    const mi::Float32 * get_radii() const { return m_radii.empty() ? 0 : &m_radii[0]; }

    mi::Size get_nb_materials() const { return m_materials.size(); }
    const mi::Uint16 * get_materials() const { return m_materials.empty() ? 0 : &m_materials[0]; }

    void serialize(mi::neuraylib::ISerializer *serializer) const;
    void deserialize(mi::neuraylib::IDeserializer* deserializer);

private:
    std::vector<mi::math::Vector_struct<mi::Float32, 3> > m_points;
    std::vector<mi::Float32> m_radii;
    std::vector<mi::Uint16> m_materials;
};

/// A set of Pipe.
class Pipe_set :
    public mi::neuraylib::Element<0x5d4e9c81,0x29f7,0x4936,0xa2,0x49,0x22,0xb0,0xef,0xe8,0xb4,0xb2,
                                  nv::index::IPipe_set>
{
public:
    /// Constructing the primitive
    Pipe_set(
        const std::vector<Pipe> &  pipes,
	mi::Float32		   radius,		// default radius if not given by pipe
	mi::math::Bbox_struct<mi::Float32, 3> & bbox);	// radii taken into account in bbox

    // default constructor
    Pipe_set()
        :
        // m_pipes(),
        m_radius(0.0f),
        m_bbox(),
        m_rendering_enabled(true),
        m_pickable(true)
    {
        // empty
    }

    virtual ~Pipe_set() {}

    virtual const nv::index::IPipe* get_pipe(mi::Uint32 index) const;
    virtual mi::Size get_nb_pipes() const;
    
    virtual mi::Float32 get_radius() const;

    // return object-space bounding box (implements IObject_space_shape::get_bounding_box)
    mi::math::Bbox_struct<mi::Float32, 3> get_bounding_box() const;

    virtual void set_enabled(bool enable)
    {
        m_rendering_enabled = enable;
    }

    virtual bool get_enabled() const
    {
        return m_rendering_enabled;
    }

    /// Each scene element can store additional user-defined meta data. Meta data,
    /// for instance, may include a string representing the scene element's name or
    /// domain specific attributes.
    /// A class that represents meta data has to be a database element and the scene
    /// element then referes to the database element by means of a tag.
    ///
    /// \param[in] tag  The tag that refers to the user-defined meta data associated
    ///                 with the scene element.
    ///
    virtual void set_meta_data(mi::neuraylib::Tag_struct tag) { m_meta_data = tag; }
    
    /// Retrieve the scene element's reference to the user-defined meta data.
    ///
    /// \return  Returns the tag that refers to the user-defined meta data
    ///          associated with the scene element.
    ///
    virtual mi::neuraylib::Tag_struct get_meta_data() const { return m_meta_data; }

    virtual bool get_pickable() const { return m_pickable; }
    virtual void set_pickable(bool pickable) { m_pickable = pickable; }

    /// -------------------------------------------------------------------------------------------
    /// \name Implemented IElement interface methods.
    ///@{
    /// The copy needs to create and return a newly allocated object of the same
    /// class and has to copy all the fields which should be the same in the copy.
    ///
    /// \return                 The new copy of the point set.
    ///
    virtual mi::neuraylib::IElement* copy() const;

    /// Get a human readable representation of the class name.
    ///
    /// \return                 The human readable string.
    ///
    virtual const char* get_class_name() const;

    /// Serialize the class to the given serializer.
    ///
    /// \param serializer       Write to this serializer.
    ///
    virtual void serialize(mi::neuraylib::ISerializer *serializer) const;

    /// Deserialize the class from the given deserializer.
    ///
    /// \param deserializer     Read from this deserializer.
    ///
    virtual void deserialize(mi::neuraylib::IDeserializer* deserializer);

    /// Get all tags referenced by the database elements.
    ///
    /// \param[in] result       The tags referenced by this database element.
    ///
    virtual void get_references(mi::neuraylib::ITag_set* result) const
    {
        if(m_meta_data.is_valid())
            result->add_tag(m_meta_data);
    }
    // ----------------------------------------------------------------------------------
    ///@}

private:
    std::vector<Pipe>                           m_pipes;
    mi::Float32                                 m_radius;
    mi::math::Bbox_struct<mi::Float32, 3>       m_bbox;
    bool                                        m_rendering_enabled;
    bool                                        m_pickable;
    mi::neuraylib::Tag                          m_meta_data;
};

}} // namespace nv::index_common

#endif // NVIDIA_INDEX_PIPE_PRIMITIVE_H
