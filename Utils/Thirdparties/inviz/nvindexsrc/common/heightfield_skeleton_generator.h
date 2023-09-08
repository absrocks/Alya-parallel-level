/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief heightfield skeleton generator

#ifndef NVIDIA_INDEX_HEIGHTFIELD_SKELETON_GENERATOR_H
#define NVIDIA_INDEX_HEIGHTFIELD_SKELETON_GENERATOR_H

#include <cassert>
#include <string>

#include <mi/dice.h>
#include <mi/base/interface_declare.h>
#include <mi/base/interface_implement.h>

#include <nv/index/idistributed_data_import_callback.h>
#include <nv/index/iregular_heightfield_patch.h>

#include "common_utility.h"
#include "heightfield_normal.h"

namespace nv
{
namespace index_common
{

/// The heightfield skeleton class. 
///
/// This heightfield generation technique generates distributed and
/// independent empty patch data. The heightfield skeleton is a place
/// holder which we can later fill in by the heightfield value via
/// heightfield compute task.
/// 
class Heightfield_skeleton_generator :
    public nv::index::Distributed_discrete_data_import_callback<0x22cbcdf8,0xf4b0,0x4d87,0x9f,0x80,0xde,0x1f,0x43,0x17,0xb2,0x88>
{
public:
    /// The default constructor is merely used for serialization/deserialization.
    Heightfield_skeleton_generator()
      : m_initial_height(nv::index::IRegular_heightfield_patch::get_hole_value())
    {
        // empty
    }

    /// The constructor requires the following parameter to generate the patches of an 
    /// synthetic heightfield:
    ///
    /// \param[in] init_height      initial height value of the heightfield placeholder
    ///
    /// \param[in] heightfield_size The size of the heightfield.
    ///
    Heightfield_skeleton_generator(
        const mi::Float32                             init_height,
        const mi::math::Vector_struct<mi::Uint32, 2>& heightfield_size)
        :
        m_initial_height(init_height),
        m_heightfield_size(heightfield_size)
    {
    }


    /// The destructor will be called by the NVIDIA IndeX library.
    virtual ~Heightfield_skeleton_generator()
    {
        // empty
    }

    // -------------------------------------------------------------------------
    // Implements the methods of the importer callback interface \c IDistributed_data_import_callback to
    // estimate memory consuptions and to create a volume brick.
    virtual mi::Size estimate(
        const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
        mi::neuraylib::IDice_transaction*               dice_transaction) const
    {
        // Size of the height field patch in memory
        const mi::Size range_i = bounding_box.max.x - bounding_box.min.x;
        const mi::Size range_j = bounding_box.max.y - bounding_box.min.y;
        const mi::Size patch_size = range_i*range_j;
        // ... elevation values + normal values per grid position:
        const mi::Size data_size = patch_size * (sizeof(mi::Float32) + sizeof(mi::Float32)*3);
        // note: color values not yet supported!
        return data_size;
    }


    virtual nv::index::IDistributed_data_subset* create(
        const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
        nv::index::IData_subset_factory*                factory,
        mi::neuraylib::IDice_transaction*               dice_transaction) const
    {
        // DEBUG_LOG << "Heightfield_skeleton_generator data: h:" << m_initial_height
        //           << ", sz: " << m_heightfield_size << ", bbox: " << bounding_box;

        // for more than 0.5G heightfield data entries, use 64bit.
        const mi::math::Vector<mi::Sint64, 2>  heightfield_dim(static_cast<mi::Sint64>(m_heightfield_size.x),
                                                               static_cast<mi::Sint64>(m_heightfield_size.y));
        mi::math::Bbox<mi::Sint64, 2> patch_raw_bbox(static_cast<mi::Sint64>(bounding_box.min.x),
                                                     static_cast<mi::Sint64>(bounding_box.min.y),
                                                     static_cast<mi::Sint64>(bounding_box.max.x),
                                                     static_cast<mi::Sint64>(bounding_box.max.y));

        patch_raw_bbox.min.x = mi::math::max(0ll, patch_raw_bbox.min.x);
        patch_raw_bbox.min.y = mi::math::max(0ll, patch_raw_bbox.min.y);
        patch_raw_bbox.max.x = mi::math::min(patch_raw_bbox.max.x, heightfield_dim.x);
        patch_raw_bbox.max.y = mi::math::min(patch_raw_bbox.max.y, heightfield_dim.y);

        // Size of the heightfield patch in memory
        const mi::math::Vector<mi::Sint64, 2> patch_raw_dim = patch_raw_bbox.max - patch_raw_bbox.min;
        assert((patch_raw_dim[0] > 0) && (patch_raw_dim[1] > 0));

        // create heightfield patch
        mi::base::Handle<nv::index::IRegular_heightfield_patch> heightfield_patch(
            factory->create<nv::index::IRegular_heightfield_patch>());
        if (!heightfield_patch.is_valid_interface())
        {
            ERROR_LOG << "Cannot create a heightfield_patch.";
            return 0;
        }
        const mi::math::Bbox_struct<mi::Sint32, 2> patch_rect_st = {
            bounding_box.min.x, bounding_box.min.y,
            bounding_box.max.x, bounding_box.max.y,
        };
        if (!heightfield_patch->initialize(patch_rect_st))
        {
            ERROR_LOG << "Fail to initialize heightfield patch";
            return 0;
        }

        mi::Float32* elevation_data = heightfield_patch->generate_elevation_storage();
        if (elevation_data == 0)
        {
            ERROR_LOG << "Cannot generate a evaluation storage.";
            return 0;
        }
        
        mi::math::Vector_struct<mi::Float32, 3>* normal_data = heightfield_patch->generate_normal_storage();
        if (normal_data == 0)
        {
            ERROR_LOG << "Cannot generate a normal storage.";
            return 0;
        }

        this->fill_height(patch_raw_bbox, elevation_data);
        this->fill_normal(patch_raw_bbox, normal_data);

        // determine extent of bounding box in k-dimension from height values
        const mi::Sint64 patch_raw_size = patch_raw_dim[0] * patch_raw_dim[1];
        mi::Float32 min_k = mi::base::numeric_traits<mi::Float32>::max();
        mi::Float32 max_k = mi::base::numeric_traits<mi::Float32>::negative_max();
        for(mi::Sint64 i = 0; i < patch_raw_size; ++i)
        {
            const mi::Float32 value = elevation_data[i];
            if (nv::index::IRegular_heightfield_patch::is_hole(value))
            {
                if (value < min_k)
                {
                    min_k = value;
                }

                if (value > max_k)
                {
                    max_k = value;
                }
            }
        }

        // k-min and k-max should be different to avoid missing planar heightfields
        if (max_k == min_k)
        {
            max_k = min_k + 1.f;
        }

        // store bounding box results. currently unused.
        // mi::math::Vector<mi::Float32, 2> height_range;
        // height_range.x = min_k;
        // height_range.y = max_k;

        heightfield_patch->retain(); // because this handle will be out of scope and the object must survive.
        return heightfield_patch.get();
    }

    virtual mi::base::Uuid subset_id() const
    {
        return nv::index::IRegular_heightfield_patch::IID();
    }

    /// Serialize the class to the given serializer.
    virtual void serialize(mi::neuraylib::ISerializer * serializer) const
    {
        serializer->write(&m_initial_height, 1);
        serializer->write(&m_heightfield_size.x, 2);
    }

    /// Deserialize the class from the given deserializer.
    virtual void deserialize(mi::neuraylib::IDeserializer * deserializer)
    {
        deserializer->read(&m_initial_height, 1);
        deserializer->read(&m_heightfield_size.x, 2);
    }

private:
    /// fill the height with the m_initial_height value
    ///
    /// \patam[in]  patch_raw_bbox bounding box of this local data
    /// \patam[out] height_data  raw height data (array) to be filled
    ///
    void fill_height(const mi::math::Bbox<mi::Sint64, 2>& patch_raw_bbox,
                     mi::Float32* const                   height_data) const
    {
        // Size of the heightfield patch in memory
        const mi::math::Vector<mi::Sint64, 2> patch_raw_dim = patch_raw_bbox.max - patch_raw_bbox.min;
        nv::index_common::no_unused_variable_warning_please(patch_raw_dim);
        assert((patch_raw_dim[0] > 0) && (patch_raw_dim[1] > 0));

        mi::math::Vector<mi::Sint64, 2> ij = patch_raw_bbox.min;
        for (ij[0] = patch_raw_bbox.min[0]; ij[0] < patch_raw_bbox.max[0]; ++ij[0])
        {
            for (ij[1] = patch_raw_bbox.min[1]; ij[1] < patch_raw_bbox.max[1]; ++ij[1])
            {
                mi::Sint64 const ij_idx = nv::index_common::get_heightfield_index(ij, patch_raw_bbox);
                assert(ij_idx >= 0);
                assert(ij_idx <  patch_raw_dim[0] * patch_raw_dim[1]);
                height_data[ij_idx] = m_initial_height;
            }
        }
    }

    /// fill the average normal
    ///
    /// \patam[in]  patch_raw_bbox bounding box of this local data
    /// \patam[in]  height_data  raw height data (array) for the normal calculation 
    /// \patam[out] normal_data  raw normal data (array) to be filled
    void fill_normal(
        const mi::math::Bbox<mi::Sint64, 2>&    patch_raw_bbox,
        mi::math::Vector_struct<mi::Float32, 3>* const normal_data) const
    {
        // The height_data is only m_initial_height, so no need to compute the normal.
        const mi::math::Vector<mi::Sint64, 2> patch_raw_dim = patch_raw_bbox.max - patch_raw_bbox.min;
        assert((patch_raw_dim[0] > 0) && (patch_raw_dim[1] > 0));
        mi::Sint64 const patch_raw_size = patch_raw_dim[0] * patch_raw_dim[1];

        const mi::math::Vector_struct<mi::Float32, 3> normal_vec = { 0.0f, 1.0f, 0.0f, };
        for (mi::Sint64 i = 0; i < patch_raw_size; i++)
        {
            normal_data[i] = normal_vec;
        }
    }

    mi::Float32                            m_initial_height;
    mi::math::Vector<mi::Uint32, 2>        m_heightfield_size;

private:
    /// copy constructor. prohibit until proved useful.
    Heightfield_skeleton_generator(Heightfield_skeleton_generator const &);
    /// operator=. prohibit until proved useful.
    Heightfield_skeleton_generator const & operator=(Heightfield_skeleton_generator const &);
};

}} // namespace nv::index_common
#endif // NVIDIA_INDEX_HEIGHTFIELD_SKELETON_GENERATOR_H
