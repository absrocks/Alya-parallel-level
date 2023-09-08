/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief image tile for distributed image store

#ifndef NVIDIA_INDEX_IMAGE_TILE_H
#define NVIDIA_INDEX_IMAGE_TILE_H


#include <mi/base/types.h>
#include <mi/dice.h>

#include <mi/math/color.h>

#include <vector>
#include <cassert>

#include "common/common_utility.h"
#include "common/encode_height.h"

#include "utilities.h"

/// image tile for distributed image store.
class Image_tile
{
public:
    /// default constructor
    Image_tile();
    /// destructor
    virtual ~Image_tile();
    /// copy constructor.
    Image_tile(Image_tile const & rhs);
    /// operator=.
    Image_tile & operator=(Image_tile const & rhs);

    /// resize buffer
    bool resize_buffer(mi::math::Bbox<mi::Sint64, 2 > const & ij_bbox);
    /// is valid
    bool is_valid_buffer() const;
    /// clear buffer
    /// \param[in] col color to fill
    void clear_buffer(mi::math::Color const & col);
    /// get this buffer size.
    /// \return (width, height)
    mi::math::Vector< mi::Sint64, 2 > get_buffer_size() const;
    /// get this bounding box.
    /// \return bounding box of this.
    mi::math::Bbox<mi::Sint64, 2 > get_bounding_box() const;
    /// put color
    /// \param[in] ij  ij position (global position)
    /// \param[in] col color to put
    void put_color(mi::math::Vector< mi::Sint64, 2 > const & ij, mi::math::Color const & col);
    /// get color
    /// \param[in] ij  ij position (global position)
    /// \return color value at ij position
    mi::math::Color get_color(mi::math::Vector< mi::Sint64, 2 > const & ij) const;

    /// partial image copy to this.
    ///
    /// \param[in]    src_image_tile source image tile
    /// \param[in]    copy_src_bbox  copy source image region. (global coords)
    /// \param[in]    dst_origin     destination min location position. (global coords)
    /// \return true when succeeded.
    bool copy_partial_image_tile(
        Image_tile const & src_image_tile,
        mi::math::Bbox<   mi::Sint64, 2 > const & copy_src_bbox,
        mi::math::Vector< mi::Sint64, 2 > const & dst_origin
        );

    /// save buffer.
    /// \param[in] fname filename to save
    bool save_buffer(std::string const & fname);
    /// load buffer.
    /// \param[in] fname filename to load
    bool load_buffer(std::string const & fname);

    /// serialize
    void serialize(mi::neuraylib::ISerializer* serializer) const;
    /// deserialize
    void deserialize(mi::neuraylib::IDeserializer* deserializer);

private:
    /// is a point inside the bbox?
    /// \param[in] ij_pos ij position
    /// \return true when ij_pos is inside the bounding box (m_ij_bbox)
    inline bool is_inside(mi::math::Vector< mi::Sint64, 2 > const & ij_pos) const
    {
        assert(this->is_valid_buffer());
        for(mi::Sint32 i = 0; i < 2; ++i){
            if((ij_pos[i] < m_ij_bbox.min[i]) || (ij_pos[i] >= m_ij_bbox.max[i])){
                return false;
            }
        }
        return true;
    }

    /// get index
    inline mi::Sint32 get_idx(mi::math::Vector< mi::Sint64, 2 > const & ij) const
    {
        assert(this->is_valid_buffer());
        assert(this->is_inside(ij));

        return static_cast<mi::Sint32>(nv::index_common::get_heightfield_index(ij, m_ij_bbox));
    }

private:
    /// image bbox
    mi::math::Bbox<mi::Sint64, 2 > m_ij_bbox;
    /// color pixel buffer
    std::vector< mi::math::Color > m_color_buffer;
};

//----------------------------------------------------------------------
#endif // NVIDIA_INDEX_IMAGE_TILE_H
